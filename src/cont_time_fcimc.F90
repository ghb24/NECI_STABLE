#include "macros.h"
module cont_time

    use CalcData, only: tContTimeFull, tOrthogonaliseReplicas, &
                        orthogonalise_iter, use_spawn_hash_table, &
                        tTruncInitiator, DiagSft, tau
    use fcimc_helper, only: rezero_iter_stats_each_iter, CalcParentFlag, &
                            create_particle, create_particle_with_hash_table
    use cont_time_rates, only: spawn_rate_full, cont_time_gen_excit_full
    use hash, only: remove_hash_table_entry, clear_hash_table
    use global_det_data, only: det_diagH, get_spawn_rate
    use orthogonalise, only: orthogonalise_replicas
    use dSFMT_interface, only: genrand_real2_dSFMT
    use AnnihilationMod, only: DirectAnnihilation
    use fcimc_iter_utils, only: update_iter_data
    use hphf_integrals, only: hphf_diag_helement
    use Determinants, only: get_helement
    use bit_reps, only: extract_bit_rep
    use LoggingData, only: FCIMCDebug
    use SystemData, only: nel, tHPHF
    use bit_reps, only: nullify_ilut
    use bit_rep_data, only: NIfTot
    use FciMCData
    use constants
    implicit none

contains

    ! TODO: Set the initial spawning rate of the first particle (In set global
    !       data).

    subroutine iterate_cont_time (iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'iterate_cont_time'

        real(dp) :: sgn(lenof_sign), rate, hdiag
        integer :: sgn_abs, iunused, flags, det(nel), j, p, TotWalkersNew
        integer :: part_type
        logical :: survives

        if (lenof_sign /= 1) then
            call stop_all(this_routine, "Complex and multi-run calculations &
                         &not supported using CONT-TIME")
        end if

        ! Reset free position list
        ValidSpawnedList = InitialSpawnedSlots
        FreeSlot(1:iEndFreeSlot) = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0
        if (use_spawn_hash_table) call clear_hash_table(spawn_ht)

        call rezero_iter_stats_each_iter(iter_data)
        call set_timer(walker_time)

        !
        ! Loop over the main list
        do j = 1, int(TotWalkers, sizeof_int)

            ! We are working in a non-contiguous list.
            fcimc_excit_gen_store%tFilled = .false.
            call extract_bit_rep(CurrentDets(:, j), det, sgn, flags, &
                                 fcimc_excit_gen_store)
            if (IsUnoccDet(sgn)) cycle

            ! Global stored data to make things efficient
            hdiag = det_diagH(j)
            if (tContTimeFull) then
                rate = get_spawn_rate(j)
            else
                call stop_all(this_routine, "Not Yet implemented")
            end if

            ! Calculate the flags that ought to be carried through
            if (tTruncInitiator) &
                call CalcParentFlag(j, iunused, hdiag)

            ! Loop over determinants, and the particles on the determinant
            do part_type = 1, lenof_sign

                ! n.b. only using integer particles atm.
                sgn_abs = int(abs(sgn(part_type)))
                do p = 1, sgn_abs

                    survives = process_part_cont_time( &
                                      CurrentDets(:,j), det, sgn(part_type), &
                                      part_type, rate, totimagtime, &
                                      totimagtime + tau, hdiag, iter_data)

                    if (.not. survives) then
                        iter_data%ndied = iter_data%ndied + 1
                        NoDied = NoDied + 1
                        sgn(part_type) = sgn(part_type) &
                                       - sign(1.0_dp, sgn(part_type))
                    end if
                end do
            end do

            ! If this particle has been completely destroyed, then remove it
            ! from the simulation
            if (all(sgn == 0)) then
                call remove_hash_table_entry(HashIndex, det, j)

                ! And add to the "freeslot" list
                iEndFreeSlot = iEndFreeSlot + 1
                FreeSlot(iEndFreeSlot) = j
                call nullify_ilut(CurrentDets(:, j))

            end if

        end do
        IFDEBUG(FCIMCDebug, 2) write(iout, '("Finnished loop over sites")')
        call halt_timer(walker_time)

        ! Send walkers to the correct nodes, and annihilate
        call set_timer(annihil_time)
        TotWalkersNew = int(TotWalkers, sizeof_int)
        call DirectAnnihilation(TotWalkersNew, iter_data, .false.)
        TotWalkers = TotWalkersNew
        call halt_timer(annihil_time)
        IFDEBUG(FCIMCDebug, 2) write(iout, '("Finished annihilation")')

        ! If we are orthogonalising the replica wavefunctions, to generate
        ! excited states, then do that here.
        if (tOrthogonaliseReplicas .and. iter > orthogonalise_iter) &
            call orthogonalise_replicas(iter_data)

        ! Update iteration counters
        call update_iter_data(iter_data)

    end subroutine

    recursive function process_part_cont_time (ilut, det, sgn, part_type, &
                                               rate, curr_time, annihil_time, &
                                               hdiag, iter_data) &
                                               result(survives)

        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer, intent(in) :: det(nel), part_type
        real(dp), intent(in) :: sgn, curr_time, annihil_time, rate, hdiag
        type(fcimc_iter_data), intent(inout) :: iter_data
        logical :: survives
        character(*), parameter :: this_routine = 'process_part_cont_time'

        real(dp) :: time, child(lenof_sign), dt, rate_adj, rate_spwn
        real(dp) :: hdiag_spwn, spwn_sgn
        integer :: nspawn, spawn_sgn, det_spwn(nel), ic, i
        integer(n_int) :: ilut_spwn(0:NIfTot)
        logical :: child_survives
        HElement_t :: hoffdiag

        ! A quick sanity check that we have calculated the spawning rate
        ! reasonably
        ASSERT(rate > 0)

        ! The starting time for the spawning
        time = curr_time

        survives = .true.
        do while (.true.)

            ! Calculate the _live_ spawning rate inside the particle loop so
            ! that the live rate gets updated as soon as the oversampling
            ! factors are tweaked
            if (tContTimeFull) then
                rate_adj = rate + abs(hdiag - DiagSft(part_type))
            else
                call stop_all(this_routine, "not yet implemented")
            end if

            ! Generate the next spawning event
            dt = - 1.0_dp * log(genrand_real2_dSFMT()) / rate_adj
            time = time + dt
            if (time > annihil_time) exit

            if (tContTimeFull) then
                nspawn = 1
                call cont_time_gen_excit_full (det, ilut, rate_adj, hdiag, &
                                               det_spwn, ilut_spwn, &
                                               hoffdiag, ic, part_type)
            else
                call stop_all(this_routine, "not yet implemented")
            end if

            if (.not. IsNullDet(det_spwn)) then

                ! If we generate the starting det, this is the same as particle
                ! death in normal FCIQMC. Simplify everything by just killing
                ! the current particle, rather than generating antiparticles.
                if (ic == 0) then
                    survives = .false.
                    exit
                end if

                if (tHPHF) then
                    hdiag_spwn = hphf_diag_helement (det_spwn, ilut_spwn)
                else
                    hdiag_spwn = get_helement (det_spwn, det_spwn, 0)
                end if
                iter_data%nborn = iter_data%nborn + nspawn
                NoBorn = NoBorn + nspawn
                if (tContTimeFull) then
                    rate_spwn = spawn_rate_full(det_spwn, ilut_spwn)
                else
                    call stop_all(this_routine, 'not yet implemented')
                end if
                spwn_sgn = sign(real(nspawn), hoffdiag)


                do i = 1, nspawn

                    child_survives = process_part_cont_time( &
                                        ilut_spwn, det_spwn, spwn_sgn, &
                                        part_type, rate_spwn, time, &
                                        annihil_time, hdiag_spwn, iter_data)

                    if (.not. child_survives) then
                        iter_data%ndied = iter_data%ndied + 1
                        NoDied = NoDied + 1
                        spwn_sgn = spwn_sgn - sign(1.0_dp, spwn_sgn)
                    end if
                end do

                ! And now create the generated particle (if it has survived
                ! this far)
                ! TODO: Ensure that the child iluts passed into this routine
                !       don't have the initiator flag set
                if (spwn_sgn /= 0) then
                    child = 0
                    child(part_type) = spwn_sgn
                    if (use_spawn_hash_table) then
                        call create_particle_with_hash_table( &
                                         det_spwn, ilut_spwn, child, &
                                         part_type, ilut)
                    else
                        call create_particle(det_spwn, ilut_spwn, child, &
                                             part_type, ilut)
                    end if
                end if

            end if ! Spawn not aborted

        end do ! Loop until annihil time


    end function

end module
