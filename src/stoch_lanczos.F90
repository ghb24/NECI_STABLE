#include "macros.h"

module stoch_lanczos

    use stoch_lanczos_procs
    implicit none

contains

    subroutine perform_stochastic_lanczos(lanczos)

        use AnnihilationMod, only: DirectAnnihilation
        use bit_rep_data, only: NIfTot, NOffFlag, tUseFlags, test_flag
        use bit_reps, only: flag_deterministic, flag_determ_parent, set_flag
        use bit_reps, only: extract_bit_rep
        use CalcData, only: AvMCExcits, tSemiStochastic, tTruncInitiator, StepsSft
        use constants
        use DetBitOps, only: FindBitExcitLevel
        use FciMCData, only: fcimc_excit_gen_store, FreeSlot, iStartFreeSlot, iEndFreeSlot
        use FciMCData, only: TotWalkers, CurrentDets, CurrentH, iLutRef, max_calc_ex_level
        use FciMCData, only: iter_data_fciqmc, TotParts, NCurrH, exFlag, iter
        use FciMCData, only: indices_of_determ_states, partial_determ_vector
        use FciMCParMod, only: create_particle, CalcParentFlag, decide_num_to_spawn
        use FciMCParMod, only: calculate_new_shift_wrapper, walker_death, end_iter_stats
        use FciMCParMod, only: update_iter_data, CalcApproxpDoubles
        use LoggingData, only: tPopsFile
        use PopsFileMod, only: WriteToPopsFileParOneArr
        use procedure_pointers, only: generate_excitation, attempt_create, encode_child
        use procedure_pointers, only: new_child_stats, extract_bit_rep_avsign
        use semi_stoch_procs, only: is_core_state, check_determ_flag, deterministic_projection
        use soft_exit, only: ChangeVars
        use SystemData, only: nel
        
        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer :: iconfig, irepeat, ivec, iiter, idet, ireplica, ispawn
        integer :: nspawn, parent_flags, unused_flags, ex_level_to_ref
        integer :: TotWalkersNew, determ_index, ic, ex(2,2)
        integer :: nI_parent(nel), nI_child(nel)
        integer(n_int) :: ilut_child(0:NIfTot)
        integer(n_int), pointer :: ilut_parent(:)
        real(dp) :: prob, unused_rdm_real, AvMCExcits_Temp
        real(dp), dimension(lenof_sign) :: child_sign, parent_sign
        real(dp), dimension(lenof_sign) :: unused_sign1, unused_sign2
        logical :: tChildIsDeterm, tParentIsDeterm, tParentUnoccupied
        logical :: tParity, tSoftExitFound, tSingBiasChange, tWritePopsFound
        HElement_t :: HElGen

        integer(n_int) :: int_sign(lenof_sign*lanczos%nvecs)
        real(dp) :: test_sign(lenof_sign*lanczos%nvecs)

        call init_stoch_lanczos(lanczos)

        outer_loop: do iconfig = 1, lanczos%nconfigs

            do irepeat = 1, lanczos%nrepeats

                call init_stoch_lanczos_repeat(lanczos, irepeat)
                call WriteFCIMCStats()

                do ivec = 1, lanczos%nvecs

                    ! Copy the current state of CurrentDets to lanczos_vecs.
                    call store_lanczos_vec(ivec, lanczos%nvecs)

                    call calc_overlap_matrix_elems(lanczos, ivec)

                    do iiter = 1, lanczos%niters

                        if (iiter == 1) then
                            ! If the projected Hamiltonian will be calculated using spawned walkers.
                            AvMCExcits_Temp = AvMCExcits
                            AvMCExcits = lanczos%av_mc_excits_sl
                        end if

                        iter = iiter + (ivec-1)*lanczos%niters
                        call init_stoch_lanczos_iter(iter_data_fciqmc, determ_index)

                        !write(6,*) "CurrentDets:"
                        !do idet = 1, int(TotWalkers, sizeof_int)
                        !    call extract_bit_rep(CurrentDets(:, idet), nI_parent, parent_sign, unused_flags, &
                        !                          fcimc_excit_gen_store)
                        !    write(6,'(i3, i12, 4x, f18.7, 4x, f18.7, 4x, l1)') idet, CurrentDets(0, idet), parent_sign, &
                        !        test_flag(CurrentDets(:, idet), flag_deterministic)
                        !end do

                        do idet = 1, int(TotWalkers, sizeof_int)

                            ! The 'parent' determinant from which spawning is to be attempted.
                            ilut_parent => CurrentDets(:,idet)
                            parent_flags = 0

                            ! Indicate that the scratch storage used for excitation generation from the
                            ! same walker has not been filled (it is filled when we excite from the first
                            ! particle on a determinant).
                            fcimc_excit_gen_store%tFilled = .false.

                            call extract_bit_rep(ilut_parent, nI_parent, parent_sign, unused_flags, &
                                                  fcimc_excit_gen_store)

                            ex_level_to_ref = FindBitExcitLevel(iLutRef, ilut_parent, max_calc_ex_level)

                            tParentIsDeterm = check_determ_flag(ilut_parent)
                            tParentUnoccupied = IsUnoccDet(parent_sign)

                            ! If this determinant is in the deterministic space then store the relevant
                            ! data in arrays for later use.
                            if (tParentIsDeterm) then
                                ! Store the index of this state, for use in annihilation later.
                                indices_of_determ_states(determ_index) = idet

                                ! Add the amplitude to the deterministic vector.
                                partial_determ_vector(:,determ_index) = parent_sign

                                determ_index = determ_index + 1

                                ! The deterministic states are always kept in CurrentDets, even when
                                ! the amplitude is zero. Hence we must check if the amplitude is zero
                                ! and, if so, skip the state.
                                if (tParentUnoccupied) cycle
                            end if

                            ! If this slot is unoccupied (and also not a core determinant) then add it to
                            ! the list of free slots and cycle.
                            if (tParentUnoccupied) then
                                iEndFreeSlot = iEndFreeSlot + 1
                                FreeSlot(iEndFreeSlot) = idet
                                cycle
                            end if

                            if (tTruncInitiator) call CalcParentFlag(idet, idet, parent_flags)

                            do ireplica = 1, inum_runs

                                call decide_num_to_spawn(parent_sign(ireplica), AvMCExcits, nspawn)
                                
                                do ispawn = 1, nspawn

                                    ! Zero the bit representation, to ensure no extraneous data gets through.
                                    ilut_child = 0

                                    call generate_excitation (nI_parent, ilut_parent, nI_child, &
                                                        ilut_child, exFlag, ic, ex, tParity, prob, &
                                                        HElGen, fcimc_excit_gen_store)

                                    ! If a valid excitation.
                                    if (.not. IsNullDet(nI_child)) then

                                        call encode_child (ilut_parent, ilut_child, ic, ex)
                                        if (tUseFlags) ilut_child(nOffFlag) = 0

                                        if (tSemiStochastic) then
                                            tChildIsDeterm = is_core_state(ilut_child)

                                            ! Is the parent state in the core space?
                                            if (tParentIsDeterm) then
                                                ! If spawning is both from and to the core space, cancel it.
                                                if (tChildIsDeterm) cycle
                                                call set_flag(ilut_child, flag_determ_parent)
                                            else
                                                if (tChildIsDeterm) call set_flag(ilut_child, flag_deterministic)
                                            end if
                                        end if

                                        child_sign = attempt_create (nI_parent, ilut_parent, parent_sign, &
                                                            nI_child, ilut_child, prob, HElGen, ic, ex, tParity, &
                                                            ex_level_to_ref, ireplica, unused_sign2, unused_rdm_real)

                                    else
                                        child_sign = 0.0_dp
                                    end if

                                    ! If any (valid) children have been spawned.
                                    if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                                        call new_child_stats (iter_data_fciqmc, ilut_parent, &
                                                              nI_child, ilut_child, ic, ex_level_to_ref,&
                                                              child_sign, parent_flags, ireplica)

                                        call create_particle (nI_child, ilut_child, child_sign, parent_flags, &
                                                              ireplica, ilut_parent, parent_sign, &
                                                              ispawn, unused_rdm_real, nspawn)

                                    end if ! If a child was spawned.

                                end do ! Over mulitple particles on same determinant.

                            end do ! Over the replicas on the same determinant.

                            ! If this is a core-space determinant then the death step is done in
                            ! deterministic_projection.
                            if (.not. tParentIsDeterm) then
                                call walker_death (iter_data_fciqmc, nI_parent, ilut_parent, CurrentH(1,idet), &
                                                    parent_sign, unused_sign2, unused_sign1, idet, idet, &
                                                    ex_level_to_ref)
                            end if

                        end do ! Over all determinants.

                        if (tSemiStochastic) call deterministic_projection()

                        TotWalkersNew = int(TotWalkers, sizeof_int)
                        call end_iter_stats(TotWalkersNew)

                        call DirectAnnihilation (TotWalkersNew, iter_data_fciqmc, .false.)

                        TotWalkers = int(TotWalkersNew, int64)

                        if (iiter == 1) then
                            if ( .not. lanczos%exact_hamil) call calc_hamil_elems_direct(lanczos, ivec)
                            ! Reset AvMCExcits to its default value.
                            AvMCExcits = AvMCExcits_Temp
                        end if

                        call update_iter_data(iter_data_fciqmc)

                        if (mod(iter, StepsSft) == 0) then
                            call calculate_new_shift_wrapper(iter_data_fciqmc, TotParts)
                            call ChangeVars(tSingBiasChange, tSoftExitFound, tWritePopsFound)
                            if (tWritePopsFound) call WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
                            if (tSingBiasChange) call CalcApproxpDoubles()
                            if (tSoftExitFound) exit outer_loop
                        end if

                    end do ! Over all iterations between Lanczos vectors.

                end do ! Over all Lanczos vectors.

                if (lanczos%exact_hamil) call calc_hamil_exact(lanczos)

            ! Sum the overlap and projected Hamiltonian matrices from the various processors.
            call communicate_lanczos_matrices(lanczos)

            call output_lanczos_matrices(lanczos, iconfig, irepeat)

            end do ! Over all repeats for a given walker configuration.

        end do outer_loop ! Over all initial walker configurations.

        if (tPopsFile) call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)

    end subroutine perform_stochastic_lanczos

end module stoch_lanczos
