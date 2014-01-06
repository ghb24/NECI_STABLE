#include "macros.h"

module stoch_lanczos

    use bit_rep_data, only: test_flag, NIfTot, NOffFlag
    use bit_reps, only: flag_deterministic, flag_determ_parent, set_flag
    use CalcData, only: AvMcExcits, tSemiStochastic, tTruncInitiator
    use constants
    use DetBitOps, only: FindBitExcitLevel
    use FciMCData, only: fcimc_excit_gen_store, FreeSlot, iStartFreeSlot, iEndFreeSlot
    use FciMCData, only: TotWalkers, CurrentDets, CurrentH, iLutRef, max_calc_ex_level
    use FciMCData, only: iter_data_fciqmc, TotParts
    use FciMCData, only: NCurrH, exFlag, indices_of_determ_states, partial_determ_vector
    use FciMCParMod, only: create_particle, CalcParentFlag, find_num_to_spawn
    use FciMCParMod, only: calculate_new_shift_wrapper
    use procedure_pointers, only: generate_excitation, attempt_create, encode_child
    use procedure_pointers, only: new_child_stats, extract_bit_rep_avsign
    use semi_stoch_procs, only: is_core_state
    use stoch_lanczos_procs
    use SystemData, only: nel

    implicit none

contains

    subroutine perform_stoch_lanczos(lanczos)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer :: iconfig, irun, ivec, iiter, iwalker, ireplica, ispawn
        integer :: nspawn, parent_flags, unused_flags, walkExcitLevel
        integer :: determ_index, ic, ex(2,2)
        integer :: nIParent(nel), nIchild(nel)
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        real(dp) :: HDiagCurr, prob, unused_rdm_real
        real(dp), dimension(lenof_sign) :: child_sign, parent_sign
        real(dp), dimension(lenof_sign) :: unused_sign1, unused_sign2
        logical :: tInDetermSpace, tParity
        HElement_t :: HElGen

        call init_stoch_lanczos(lanczos)

        do iconfig = 1, lanczos%nconfigs

            do irun = 1, lanczos%nruns

                call create_initial_config(lanczos, irun)

                do ivec = 1, lanczos%nkrylov_vecs

                    do iiter = 1, lanczos%niters

                        do iwalker = 1, int(TotWalkers, sizeof_int)

                            ! Indicate that the scratch storage used for excitation generation
                            ! from the same walker has not been filled (it is filled when we
                            ! excite from the first particle on a determinant).
                            fcimc_excit_gen_store%tFilled = .false.

                            call extract_bit_rep_avsign(CurrentDets(:,iwalker), CurrentH(1:NCurrH,iwalker), &
                                                        nIParent, parent_sign, unused_flags, unused_sign1, &
                                                        unused_sign2, fcimc_excit_gen_store)

                            walkExcitLevel = FindBitExcitLevel(iLutRef, CurrentDets(:,iwalker), &
                                                                max_calc_ex_level)

                            ! If this state is in the deterministic space.
                            if (tSemiStochastic) then
                                if (test_flag(CurrentDets(:,iwalker), flag_deterministic)) then

                                    ! Store the index of this state, for use in annihilation later.
                                    indices_of_determ_states(determ_index) = iwalker

                                    ! Add this amplitude to the deterministic vector.
                                    partial_determ_vector(:,determ_index) = parent_sign

                                    determ_index = determ_index + 1

                                    ! The deterministic states are always kept in CurrentDets, even when
                                    ! the amplitude is zero. Hence we must check if the amplitude is zero,
                                    ! and if so, skip the state.
                                    if (IsUnoccDet(parent_sign)) cycle
                                 
                                end if
                            end if

                            ! If this slot is unoccupied, add it to the list of free slots and cycle.
                            if (IsUnoccDet(parent_sign)) then
                                iEndFreeSlot = iEndFreeSlot + 1
                                FreeSlot(iEndFreeSlot) = iwalker
                                cycle
                            end if

                            if (tTruncInitiator) call CalcParentFlag (iwalker, iwalker, parent_flags)
                            HDiagCurr = CurrentH(1,iwalker)

                            do ireplica = 1, inum_runs

                                call find_num_to_spawn(parent_sign(ireplica), AvMCExcits, nspawn)
                                
                                do ispawn = 1, nspawn

                                    ! Zero the bit representation, to ensure no extraneous data gets through.
                                    ilut_child = 0

                                    call generate_excitation (nIParent, CurrentDets(:,iwalker), nIChild, &
                                                        ilut_child, exFlag, ic, ex, tParity, prob, &
                                                        HElGen, fcimc_excit_gen_store)
                                    
                                    ! If a valid excitation.
                                    if (.not. IsNullDet(nIChild)) then

                                        if (tSemiStochastic) then
                                            call encode_child (CurrentDets(:,iwalker), ilut_child, ic, ex)

                                            ilut_child(nOffFlag) = 0
                                            
                                            ! Is the spawned state in the core space?
                                            tInDetermSpace = is_core_state(ilut_child)

                                            ! Is the parent state in the core space?
                                            if (test_flag(CurrentDets(:,iwalker), flag_deterministic)) then
                                                ! If spawning is from and to the core space, cancel it.
                                                if (tInDetermSpace) cycle
                                            else
                                                if (tInDetermSpace) call set_flag(ilut_child, flag_deterministic)
                                            end if

                                            ! If the walker being spawned is spawned from the deterministic space,
                                            ! then set the corresponding flag to specify this.
                                            if (test_flag(CurrentDets(:,iwalker), flag_deterministic)) &
                                                call set_flag(ilut_child, flag_determ_parent)

                                        end if

                                        child_sign = attempt_create (nIParent, CurrentDets(:,iwalker), parent_sign, &
                                                            nIChild, ilut_child, prob, HElGen, ic, ex, tParity, &
                                                            walkExcitLevel, ireplica, unused_sign2, unused_rdm_real)
                                    else
                                        child_sign = 0.0_dp
                                    end if

                                    ! If any (valid) children have been spawned.
                                    if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic < 2)) then

                                        ! Encode child if not done already.
                                        if (.not. (tSemiStochastic)) call encode_child (CurrentDets(:,iwalker), &
                                                                                        ilut_child, ic, ex)

                                        call new_child_stats (iter_data_fciqmc, CurrentDets(:,iwalker), &
                                                              nIChild, ilut_child, ic, walkExcitLevel,&
                                                              child_sign, parent_flags, ireplica)

                                        call create_particle (nIChild, ilut_child, child_sign, parent_flags, &
                                                              ireplica, CurrentDets(:,iwalker), parent_sign, &
                                                              ispawn, unused_rdm_real, nspawn)

                                    end if ! If a child was spawned.

                                end do ! Over mulitple particles on same determinant.

                            end do ! Over the two replicas on the same determinant.

                        end do ! Over all determinants.

                        call calculate_new_shift_wrapper(iter_data_fciqmc, TotParts)

                    end do ! Over all iterations between Lanczos vectors.

                end do ! Over all Lanczos vectors.

            end do ! Over all repeats for a given walker configuration.

        end do ! Over all initial walker configurations.

    end subroutine perform_stoch_lanczos

end module stoch_lanczos
