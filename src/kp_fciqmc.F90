! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

module kp_fciqmc

    use kp_fciqmc_hamil
    use kp_fciqmc_procs
    implicit none

contains

    subroutine perform_kp_fciqmc(kp)

        use AnnihilationMod, only: DirectAnnihilation
        use bit_rep_data, only: NIfTot, NOffFlag, tUseFlags, test_flag
        use bit_reps, only: flag_deterministic, flag_determ_parent, set_flag
        use bit_reps, only: extract_bit_rep, flag_is_initiator
        use CalcData, only: AvMCExcits, tSemiStochastic, tTruncInitiator, StepsSft
        use constants
        use DetBitOps, only: FindBitExcitLevel, return_ms
        use FciMCData, only: fcimc_excit_gen_store, FreeSlot, iEndFreeSlot
        use FciMCData, only: TotWalkers, CurrentDets, iLutRef, max_calc_ex_level
        use FciMCData, only: iter_data_fciqmc, TotParts, exFlag, iter
        use FciMCData, only: indices_of_determ_states, partial_determ_vector
        use fcimc_initialisation, only: CalcApproxpDoubles
        use fcimc_helper, only: SumEContrib, end_iter_stats, create_particle, &
                                CalcParentFlag, walker_death, &
                                decide_num_to_spawn
        use fcimc_output, only: end_iteration_print_warn
        use fcimc_iter_utils, only: calculate_new_shift_wrapper, &
                                    update_iter_data
        use global_det_data, only: det_diagH
        use LoggingData, only: tPopsFile
        use Parallel_neci, only: iProcIndex
        use ParallelHelper, only: root
        use PopsFileMod, only: WriteToPopsFileParOneArr
        use procedure_pointers, only: generate_excitation, attempt_create, encode_child
        use procedure_pointers, only: new_child_stats, extract_bit_rep_avsign
        use semi_stoch_procs, only: is_core_state, check_determ_flag, determ_projection
        use soft_exit, only: ChangeVars
        use SystemData, only: nel, lms, nbasis, tAllSymSectors, nOccAlpha, nOccBeta
        use timing_neci, only: set_timer, halt_timer

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: iiter, idet, ireplica, ispawn
        integer, target :: iconfig, irepeat, ivec
        integer :: nspawn, parent_flags, unused_flags, ex_level_to_ref
        integer :: TotWalkersNew, determ_ind, ic, ex(2,2), ms_parent
        integer :: nI_parent(nel), nI_child(nel), unused_vecslot
        integer(n_int) :: ilut_child(0:NIfTot)
        integer(n_int), pointer :: ilut_parent(:)
        real(dp) :: prob, unused_rdm_real, parent_hdiag
        real(dp), dimension(lenof_sign) :: child_sign, parent_sign
        real(dp), dimension(lenof_sign) :: unused_sign1, unused_sign2
        logical :: tChildIsDeterm, tParentIsDeterm, tParentUnoccupied
        logical :: tParity, tSoftExitFound, tSingBiasChange, tWritePopsFound
        HElement_t :: HElGen

        ! Variables to hold information output for the test suite.
        real(dp) :: s_sum, h_sum

        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: test_sign(lenof_sign_kp)
        type(ll_node), pointer :: temp_node

        kp%iconfig => iconfig
        kp%irepeat => irepeat
        kp%ivec => ivec
        call init_kp_fciqmc(kp)
        if (.not. tAllSymSectors) ms_parent = lms

        outer_loop: do iconfig = 1, kp%nconfigs

            do irepeat = 1, kp%nrepeats

                ! Point to the regions of memory where the projected Hamiltonian
                ! and overlap matrices for this repeat will be accumulated and stored.
                if (tStoreKPMatrices) then
                    kp%hamil_matrix => kp%hamil_matrices(:,:,irepeat)
                    kp%overlap_matrix => kp%overlap_matrices(:,:,irepeat)
                else
                    kp%hamil_matrix => kp%hamil_matrices(:,:,1)
                    kp%overlap_matrix => kp%overlap_matrices(:,:,1)
                end if

                call init_kp_fciqmc_repeat(kp)
                call WriteFCIMCStats()

                do ivec = 1, kp%nvecs

                    ! Copy the current state of CurrentDets to krylov_vecs.
                    call store_krylov_vec(kp)

                    call calc_overlap_matrix_elems(kp)

                    ! Calculate the overlap of the perturbed ground state vector
                    ! with the new Krylov vector, if requested.
                    if (tOverlapPert) call calc_perturbation_overlap(kp)

                    do iiter = 1, kp%niters(ivec)

                        call set_timer(walker_time)

                        iter = iter + 1
                        call init_kp_fciqmc_iter(iter_data_fciqmc, determ_ind)

                        !if (iter < 10) then
                        !    write(6,*) "CurrentDets before:"
                        !    do idet = 1, int(TotWalkers, sizeof_int)
                        !        call extract_bit_rep(CurrentDets(:, idet), nI_parent, parent_sign, unused_flags, &
                        !                              fcimc_excit_gen_store)
                        !        if (tUseFlags) then
                        !            write(6,'(i7, i12, 4x, f18.7, 4x, f18.7, 4x, l1)') idet, CurrentDets(0,idet), parent_sign, &
                        !                test_flag(CurrentDets(:,idet), flag_deterministic)
                        !        else
                        !            write(6,'(i7, i12, 4x, f18.7, 4x, f18.7)') idet, CurrentDets(0,idet), parent_sign
                        !        end if
                        !    end do

                        !    write(6,"(A)") "Hash Table: "
                        !    do idet = 1, nWalkerHashes
                        !        temp_node => HashIndex(idet)
                        !        if (temp_node%ind /= 0) then
                        !            write(6,'(i9)',advance='no') idet
                        !            do while (associated(temp_node))
                        !                write(6,'(i9)',advance='no') temp_node%ind
                        !                temp_node => temp_node%next
                        !            end do
                        !            write(6,'()',advance='yes')
                        !        end if
                        !    end do
                        !end if

                        do idet = 1, int(TotWalkers, sizeof_int)

                            ! The 'parent' determinant from which spawning is to be attempted.
                            ilut_parent => CurrentDets(:,idet)
                            parent_flags = 0_n_int

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
                                indices_of_determ_states(determ_ind) = idet
                                ! Add the amplitude to the deterministic vector.
                                partial_determ_vector(:,determ_ind) = parent_sign
                                determ_ind = determ_ind + 1

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

                            ! The current diagonal matrix element is stored persistently.
                            parent_hdiag = det_diagH(idet)

                            if (tTruncInitiator) call CalcParentFlag(idet, idet, parent_flags, parent_hdiag)

                            call SumEContrib (nI_parent, ex_level_to_ref, parent_sign, ilut_parent, &
                                               parent_hdiag, 1.0_dp, idet)

                            if (tAllSymSectors) then
                                ms_parent = return_ms(ilut_parent)
                                nOccAlpha = (nel+ms_parent)/2
                                nOccBeta = (nel-ms_parent)/2
                            end if

                            ! If this condition is not met (if all electrons have spin up or all have spin down)
                            ! then there will be no determinants to spawn to, so don't attempt spawning.
                            if (abs(ms_parent) /= nbasis/2) then

                                do ireplica = 1, lenof_sign

                                    call decide_num_to_spawn(parent_sign(ireplica), AvMCExcits, nspawn)
                                    
                                    do ispawn = 1, nspawn

                                        ! Zero the bit representation, to ensure no extraneous data gets through.
                                        ilut_child = 0_n_int

                                        call generate_excitation (nI_parent, ilut_parent, nI_child, &
                                                            ilut_child, exFlag, ic, ex, tParity, prob, &
                                                            HElGen, fcimc_excit_gen_store)

                                        ! If a valid excitation.
                                        if (.not. IsNullDet(nI_child)) then

                                            call encode_child (ilut_parent, ilut_child, ic, ex)
                                            if (tUseFlags) ilut_child(nOffFlag) = 0_n_int

                                            if (tSemiStochastic) then
                                                tChildIsDeterm = is_core_state(ilut_child, nI_child)

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

                            end if ! If connected determinants exist to spawn to.

                            ! If this is a core-space determinant then the death step is done in
                            ! determ_projection.
                            if (.not. tParentIsDeterm) then
                                call walker_death (iter_data_fciqmc, nI_parent, ilut_parent, parent_hdiag, &
                                                    parent_sign, unused_sign2, unused_sign1, unused_vecslot, idet, &
                                                    ex_level_to_ref)
                            end if

                        end do ! Over all determinants.

                        if (tSemiStochastic) call determ_projection()

                        TotWalkersNew = int(TotWalkers, sizeof_int)
                        call end_iter_stats(TotWalkersNew)
                        call end_iteration_print_warn(TotWalkersNew)

                        call halt_timer(walker_time)

                        call set_timer(annihil_time)

                        call DirectAnnihilation (TotWalkersNew, iter_data_fciqmc, .false.)

                        TotWalkers = int(TotWalkersNew, int64)

                        call halt_timer(annihil_time)

                        if (iiter == 1 .and. tHamilOnFly) call calc_hamil_on_fly(kp)

                        call update_iter_data(iter_data_fciqmc)

                        if (mod(iter, StepsSft) == 0) then
                            call calculate_new_shift_wrapper(iter_data_fciqmc, TotParts)
                            call ChangeVars(tSingBiasChange, tSoftExitFound, tWritePopsFound)
                            if (tWritePopsFound) call WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
                            if (tSingBiasChange) call CalcApproxpDoubles()
                            if (tSoftExitFound) exit outer_loop
                        end if

                    end do ! Over all iterations between Krylov vectors.

                end do ! Over all Krylov vectors.

                if (tExactHamil) then
                    call calc_hamil_exact(kp)
                else if (.not. tHamilOnFly) then
                    call calc_projected_hamil(kp)
                end if

                ! Sum the overlap and projected Hamiltonian matrices from the various processors.
                call communicate_kp_matrices(kp)

                call output_kp_matrices_wrapper(kp)

            end do ! Over all repeats for a fixed initial walker configuration.

            if (tOverlapPert) call average_and_comm_pert_overlaps(kp%nrepeats)

            if (iProcIndex == root .and. tStoreKPMatrices) then
                call average_kp_matrices_wrapper(kp)
                call find_and_output_lowdin_eigv(kp)
                call find_and_output_gs_eigv(kp)

                ! Calculate data for the testsuite.
                s_sum = sum(kp_overlap_mean)
                h_sum = sum(kp_hamil_mean)
            end if

        end do outer_loop ! Over all initial walker configurations.

        if (tPopsFile) call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)

        if (iProcIndex == root) call write_kpfciqmc_testsuite_data(s_sum, h_sum)

    end subroutine perform_kp_fciqmc

end module kp_fciqmc
