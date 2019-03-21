#include "macros.h"
module cont_time

    use CalcData, only: tContTimeFull, tOrthogonaliseReplicas, &
                        orthogonalise_iter, use_spawn_hash_table, &
                        tTruncInitiator, DiagSft, tau, tPairedReplicas
    use fcimc_helper, only: rezero_iter_stats_each_iter, CalcParentFlag, &
                            create_particle, create_particle_with_hash_table, &
                            SumEContrib, end_iter_stats
    use cont_time_rates, only: spawn_rate_full, cont_time_gen_excit_full, &
                               oversample_factors, cont_time_gen_excit, ostag,&
                               secondary_gen_store
    use hash, only: remove_hash_table_entry, clear_hash_table
    use DetBitOps, only: FindBitExcitLevel, count_open_orbs
    use global_det_data, only: det_diagH, get_spawn_rate
    use GenRandSymExcitNUMod, only: init_excit_gen_store
    use Determinants, only: get_helement, write_det
    use orthogonalise, only: orthogonalise_replicas
    use dSFMT_interface, only: genrand_real2_dSFMT
    use AnnihilationMod, only: DirectAnnihilation
    use bit_reps, only: nullify_ilut, encode_sign
    use fcimc_iter_utils, only: update_iter_data
    use hphf_integrals, only: hphf_diag_helement
    use SystemData, only: nel, tHPHF, LMS
    use bit_reps, only: extract_bit_rep
    use LoggingData, only: FCIMCDebug
    use bit_rep_data, only: NIfTot
    use rdm_data, only: rdm_definitions
    use FciMCData
    use constants
    use util_mod
    implicit none
    save

contains

    ! TODO: Set the initial spawning rate of the first particle (In set global
    !       data).

    subroutine iterate_cont_time (iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'iterate_cont_time'

        real(dp) :: sgn(lenof_sign), rate, hdiag
        integer :: sgn_abs, iunused, flags, det(nel), j, p, TotWalkersNew
        integer :: part_type, ic_hf, nopen
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

        call rezero_iter_stats_each_iter(iter_data, rdm_definitions)
        call set_timer(walker_time)

        !
        ! Loop over the main list
        do j = 1, int(TotWalkers, sizeof_int)

            fcimc_excit_gen_store%tFilled = .false.
            call extract_bit_rep(CurrentDets(:, j), det, sgn, flags, &
                                 fcimc_excit_gen_store)

            ! As the main list is not contiguous, skip (but store) empty sites.
            if (IsUnoccDet(sgn)) then
                iEndFreeSlot = iEndFreeSlot + 1
                FreeSlot(iEndFreeSlot) = j
                cycle
            end if

            IFDEBUG(FCIMCDebug, 3) then
                write(iout, "(A,I10,a)", advance='no') 'TW:', j, '['
                do part_type = 1, lenof_sign
                    write(iout, "(f10.5)", advance='no') sgn(part_type)
                end do
                write(iout, '(a)', advance='no') '] '
                call WriteBitDet(iout, CurrentDets(:,j), .true.)
                call neci_flush(iout) 
            end if

            ! Global stored data to make things efficient
            hdiag = det_diagH(j)
            if (tContTimeFull) then
                rate = get_spawn_rate(j)
                ASSERT(rate == spawn_rate_full(det, CurrentDets(:,j)))
            end if

            ! Calculate the flags that ought to be carried through
            if (tTruncInitiator) &
                call CalcParentFlag(j, iunused)

            ! Sum in the energy terms, yeah!
            ic_hf = FindBitExcitLevel(ilutRef(:,1), CurrentDets(:,j))
            call SumEContrib(det, ic_hf, sgn, CurrentDets(:,j), hdiag, 1.0_dp,&
                             tPairedReplicas, j)

            ! Needed for calculating oversample factors
            nopen = count_open_orbs(CurrentDets(:,j))

            ! Loop over determinants, and the particles on the determinant
            do part_type = 1, lenof_sign

                ! n.b. only using integer particles atm.
                sgn_abs = int(abs(sgn(part_type)))
                do p = 1, sgn_abs

                    survives = process_part_cont_time( &
                                      CurrentDets(:,j), det, sgn(part_type), &
                                      part_type, rate, totimagtime, &
                                      totimagtime + tau, hdiag, iter_data, &
                                      nopen, fcimc_excit_gen_store)

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
            else
                call encode_sign(CurrentDets(:,j), sgn)
            end if

        end do
        IFDEBUG(FCIMCDebug, 2) write(iout, '("Finnished loop over sites")')
        call halt_timer(walker_time)

        ! Update statistics. This is done before annihilation as the output
        ! statistics normally apply to the iteration _before_ the one that has
        ! just run (as they are accumulated during the normal loop)
        TotWalkersNew = int(TotWalkers, sizeof_int)
        call end_iter_stats(TotWalkersNew)

        ! Send walkers to the correct nodes, and annihilate
        call set_timer(annihil_time)
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
                                               hdiag, iter_data, nopen, &
                                               store) &
                                               result(survives)

        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer, intent(in) :: det(nel), part_type, nopen
        real(dp), intent(in) :: sgn, curr_time, annihil_time, rate, hdiag
        type(fcimc_iter_data), intent(inout) :: iter_data
        logical :: survives
        type(excit_gen_store_type), intent(inout) :: store
        character(*), parameter :: this_routine = 'process_part_cont_time'

        real(dp) :: time, child(lenof_sign), dt, rate_adj, rate_spwn
        real(dp) :: hdiag_spwn, spwn_sgn
        integer :: nspawn, spawn_sgn, det_spwn(nel), ic, i, y, nopen_spwn, err
        integer(n_int) :: ilut_spwn(0:NIfTot)
        logical :: child_survives
        HElement_t(dp) :: hoffdiag, htmp

        ! A quick sanity check that we have calculated the spawning rate
        ! reasonably
        ASSERT(rate > 0 .or. .not. tContTimeFull)

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
                rate_adj = sum(oversample_factors(:, nopen)) &
                         + abs(hdiag - DiagSft(part_type))
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
                call cont_time_gen_excit (det, ilut, rate_adj, hdiag, &
                                          det_spwn, ilut_spwn, hoffdiag, &
                                          ic, part_type, nopen, nspawn, store)
            end if

            if (.not. IsNullDet(det_spwn)) then

                ! If we generate the starting det, this is the same as particle
                ! death in normal FCIQMC. Simplify everything by just killing
                ! the current particle, rather than generating antiparticles.
                if (ic == 0) then
                    survives = .false.
                    IFDEBUG(FCIMCDebug, 3) write(iout, '("Particle died")')
                    exit
                end if

                ! We need the diagonal matrix element for calculating further
                ! spawning rates
                if (tHPHF) then
                    htmp = hphf_diag_helement (det_spwn, ilut_spwn)
                else
                    htmp = get_helement (det_spwn, det_spwn, 0)
                end if
                hdiag_spwn = real(htmp, dp) - Hii

                ! Calculate full rate if needed
                if (tContTimeFull) &
                    rate_spwn = spawn_rate_full(det_spwn, ilut_spwn)

                ! Particle accountancy
                iter_data%nborn = iter_data%nborn + nspawn
                NoBorn = NoBorn + nspawn

                ! This is the spawning coefficient. n.b. it is possible to
                ! spawn more than one particle whilst adjusting the overspawn
                ! factors
#ifndef __CMPLX
                spwn_sgn = - sign(1.0_dp, sgn) * sign(real(nspawn,dp),hoffdiag)
#else
                call stop_all(this_routine, "Not implemented")
#endif

                ! Need this for calculating further oversampling factors
                nopen_spwn = count_open_orbs(ilut_spwn)

                ! We can't store data in this structure recursively...
                secondary_gen_store%tFilled = .false.

                do i = 1, nspawn

                    child_survives = process_part_cont_time( &
                                        ilut_spwn, det_spwn, spwn_sgn, &
                                        part_type, rate_spwn, time, &
                                        annihil_time, hdiag_spwn, iter_data, &
                                        nopen_spwn, secondary_gen_store)

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

                    IFDEBUG(FCIMCDebug, 3) then
                        write(iout, '(a)', advance='no') 'SP: ['
                        do y = 1, lenof_sign
                            write(iout, '(f12.5)', advance='no') &
                                child(y)
                        end do
                        write(iout, '("] ")', advance='no')
                        call write_det(6, det_spwn, .true.)
                        call neci_flush(iout) 
                    end if

                    if (use_spawn_hash_table) then
                        call create_particle_with_hash_table(&
                                         det_spwn, ilut_spwn, child, &
                                         part_type, ilut, iter_data, err)
                    else
                        call create_particle(det_spwn, ilut_spwn, child, &
                                             part_type, err, ilut)
                    end if
                end if

            end if ! Spawn not aborted

        end do ! Loop until annihil time


    end function

end module
