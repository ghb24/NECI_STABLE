#include "macros.h"

module kp_fciqmc_hamil

    use bit_rep_data, only: NIfTot, NIfDBO, NOffFlag
    use kp_fciqmc_procs
    use SystemData, only: nel
    implicit none

contains

    subroutine calc_projected_hamil(kp)

        use bit_reps, only: decode_bit_det
        use CalcData, only: tSemiStochastic
        use constants
        use DetBitOps, only: return_ms
        use FciMCData, only: fcimc_excit_gen_store, exFlag, partial_determ_vecs_kp
        use procedure_pointers, only: generate_excitation, encode_child, get_spawn_helement
        use semi_stoch_procs, only: is_core_state, check_determ_flag
        use semi_stoch_procs, only: deterministic_projection_kp_hamil
        use SystemData, only: nbasis, tAllSymSectors, nOccAlpha, nOccBeta
        use util_mod, only: stochastic_round
        
        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, ispawn, nspawn, i, j
        integer :: determ_ind, flag_ind, ic, ex(2,2), ms_parent
        integer :: nI_parent(nel), nI_child(nel)
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        real(dp) :: prob, tot_pop
        real(dp), dimension(lenof_sign_kp) :: child_sign, parent_sign
        integer(n_int) :: int_sign(lenof_sign_kp)
        logical :: tChildIsDeterm, tParentIsDeterm, tParentUnoccupied, tParity
        logical :: tNearlyFull, tFinished, tAllFinished
        HElement_t :: HElGen, HEl

        kp%hamil_matrix = 0.0_dp

        ilut_parent = 0_n_int
        if (tUseFlags) flag_ind = NIfDBO + lenof_sign_kp + 2
        ValidSpawnedList = InitialSpawnedSlots
        tNearlyFull = .false.
        tFinished = .false.
        determ_ind = 1

        do idet = 1, TotWalkersKp

            ! The 'parent' determinant from which spawning is to be attempted.
            ilut_parent(0:NIfDBO) = krylov_vecs(0:NIfDBO,idet)
            if (tUseFlags) ilut_parent(NOffFlag) = krylov_vecs(flag_ind,idet)

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call decode_bit_det(nI_parent, ilut_parent)
            int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, idet)
            parent_sign = transfer(int_sign, parent_sign)
            tot_pop = sum(abs(parent_sign))

            tParentIsDeterm = check_determ_flag(ilut_parent)
            tParentUnoccupied = IsUnoccDet(parent_sign)

            ! If this determinant is in the deterministic space then store the relevant
            ! data in arrays for later use.
            if (tParentIsDeterm) then
                ! Add the amplitude to the deterministic vector.
                partial_determ_vecs_kp(:,determ_ind) = parent_sign
                determ_ind = determ_ind + 1
            end if

            if (tParentUnoccupied) cycle

            if (tAllSymSectors) then
                ms_parent = return_ms(ilut_parent)
                ! If this condition is met (if all electrons have spin up or all have spin
                ! down) then there will be no determinants to spawn to, so don't attempt spawning.
                if (abs(ms_parent) == nbasis/2) cycle
                nOccAlpha = (nel+ms_parent)/2
                nOccBeta = (nel-ms_parent)/2
            end if

            nspawn = stochastic_round(av_mc_excits_kp*tot_pop)
            
            do ispawn = 1, nspawn

                ! Zero the bit representation, to ensure no extraneous data gets through.
                ilut_child = 0_n_int

                call generate_excitation (nI_parent, ilut_parent, nI_child, &
                                    ilut_child, exFlag, ic, ex, tParity, prob, &
                                    HElGen, fcimc_excit_gen_store)

                ! If a valid excitation.
                if (.not. IsNullDet(nI_child)) then

                    call encode_child (ilut_parent, ilut_child, ic, ex)

                    ! If spawning is both too and from the core space, abort it.
                    if (tSemiStochastic) then
                        tChildIsDeterm = is_core_state(ilut_child)
                        if (tParentIsDeterm .and. tChildIsDeterm) cycle
                    end if

                    HEl = get_spawn_helement(nI_parent, nI_child, ilut_parent, ilut_child, ic, ex, &
                                             tParity, HElGen)

                    child_sign = calc_amp_kp_hamil(parent_sign, prob, av_mc_excits_kp*tot_pop, HEl)
                else
                    child_sign = 0.0_dp
                end if

                ! If any (valid) children have been spawned.
                if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                    call create_particle_kp_hamil(nI_child, ilut_child, child_sign, tNearlyFull)

                    if (tNearlyFull) then
                        call add_in_hamil_contribs(kp, tFinished, tAllFinished)
                        tNearlyFull = .false.
                    end if

                end if ! If a child was spawned.

            end do ! Over mulitple spawns from the same determinant.

        end do ! Over all determinants.

        tFinished = .true.
        do
            call add_in_hamil_contribs(kp, tFinished, tAllFinished)
            if (tAllFinished) exit
        end do

        call calc_hamil_contribs_diag(kp)

        if (tSemiStochastic) then
            call deterministic_projection_kp_hamil()
            call calc_hamil_contribs_semistoch(kp)
        end if

        ! Symmetrise the projected Hamiltonian.
        do i = 1, kp%nvecs
            do j = 1, i-1
                kp%hamil_matrix(i,j) = kp%hamil_matrix(j,i)
            end do
        end do

    end subroutine calc_projected_hamil

    function calc_amp_kp_hamil(parent_sign, prob, av_nspawn, HEl) result(child_sign)

        real(dp), intent(in) :: parent_sign(lenof_sign_kp)
        real(dp), intent(in) :: prob
        real(dp), intent(in) :: av_nspawn
        HElement_t, intent(in) :: HEl
        real(dp) :: child_sign(lenof_sign_kp)
        real(dp) :: Hel_real, corrected_prob

        HEl_real = real(HEl, dp)
        corrected_prob = prob*av_nspawn
        child_sign = HEl*parent_sign/corrected_prob

    end function calc_amp_kp_hamil

    subroutine create_particle_kp_hamil (nI_child, ilut_child, child_sign, tNearlyFull)

        use hash, only: DetermineDetNode
        use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots

        integer, intent(in) :: nI_child(nel)
        integer(n_int), intent(in) :: ilut_child(0:NIfTot)
        real(dp), intent(in) :: child_sign(lenof_sign_kp)
        logical, intent(inout) :: tNearlyFull
        integer(n_int) :: int_sign(lenof_sign_kp)
        integer :: proc

        proc = DetermineDetNode(nI_child, 0)

        SpawnedPartsKP(0:NIfDBO, ValidSpawnedList(proc)) = ilut_child(0:NIfDBO)
        int_sign = transfer(child_sign, int_sign)
        SpawnedPartsKP(NOffSgn:NOffSgn+lenof_sign_kp-1, ValidSpawnedList(proc)) = int_sign

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1

        if (ValidSpawnedList(proc)-InitialSpawnedSlots(proc) > MaxSpawnedEachProc) tNearlyFull = .true.

    end subroutine create_particle_kp_hamil

    subroutine distribute_spawns_kp_hamil(nspawns_this_proc)

        use AnnihilationMod, only: SendProcNewParts
        use FciMCData, only: SpawnedParts, SpawnedParts2, SpawnedPartsKP, SpawnedPartsKP2

        integer, intent(out) :: nspawns_this_proc
        integer(n_int), pointer :: PointTemp(:,:), PointTemp2(:,:)

        ! We want SpawnedParts and SpawnedParts2 to point to the 'wider' spawning arrays which
        ! hold all the signs of *all* the Krylov vectors. This is because SendProcNewParts uses
        ! SpawnedParts and SpawnedParts2.
        PointTemp => SpawnedParts
        PointTemp2 => SpawnedParts2
        SpawnedParts => SpawnedPartsKP
        SpawnedParts2 => SpawnedPartsKP2

        call SendProcNewParts(nspawns_this_proc, .false.)

        ! Now we want SpawnedPartsKP to point to SpawnedParts2, which holds the output of
        ! the communication of the spawns.
        SpawnedPartsKP => SpawnedParts2
        SpawnedPartsKP2 => SpawnedParts

        ! Now let SpawnedParts1 and SpawnedParts2 to point to their original arrays.
        SpawnedParts => PointTemp
        SpawnedParts2 => PointTemp2

        nullify(PointTemp)
        nullify(PointTemp2)

    end subroutine distribute_spawns_kp_hamil

    subroutine add_in_hamil_contribs(kp, tFinished, tAllFinished)

        use Parallel_neci, only: MPIAllGather

        type(kp_fciqmc_data), intent(inout) :: kp
        logical, intent(in) :: tFinished
        logical, intent(out) :: tAllFinished
        logical :: tFinished_AllProcs(nProcessors)
        integer :: nspawns_this_proc, ierr

        call distribute_spawns_kp_hamil(nspawns_this_proc)

        call calc_hamil_contribs_spawn(kp, nspawns_this_proc)

        ValidSpawnedList = InitialSpawnedSlots

        call MPIAllGather(tFinished, 1, tFinished_AllProcs, 1, ierr)

        tAllFinished = all(tFinished_AllProcs)

    end subroutine add_in_hamil_contribs

    subroutine calc_hamil_contribs_spawn(kp, nspawns_this_proc)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer, intent(in) :: nspawns_this_proc
        integer :: idet, DetHash, det_ind, i, j
        integer(n_int) :: ilut_spawn(0:NIfTot)
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign_1(lenof_sign_kp), real_sign_2(lenof_sign_kp)
        type(ll_node), pointer :: temp_node
        logical :: tDetFound

        ilut_spawn = 0_n_int

        do idet = 1, nspawns_this_proc
            call neci_flush(6)
            ilut_spawn(0:NIfDBO) = SpawnedPartsKP(0:NIfDBO, idet)
            int_sign = SpawnedPartsKP(NOffSgn:NOffSgn+lenof_sign_kp-1, idet)
            real_sign_1 = transfer(int_sign, real_sign_1)
            call decode_bit_det(nI_spawn, ilut_spawn)
            DetHash = FindWalkerHash(nI_spawn, nhashes_kp)
            ! Point to the first node with this hash value in krylov_vecs.
            temp_node => krylov_vecs_ht(DetHash)
            if (temp_node%ind == 0) then
                ! If there are no determinants at all with this hash value in krylov_vecs.
                cycle
            else
                tDetFound = .false.
                do while (associated(temp_node))
                    if (DetBitEQ(ilut_spawn, krylov_vecs(:,temp_node%ind), NIfDBO)) then
                        ! If this CurrentDets determinant has been found in krylov_vecs.
                        det_ind = temp_node%ind
                        tDetFound = .true.
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    temp_node => temp_node%next
                end do
                if (tDetFound) then
                    int_sign = krylov_vecs(NOffSgn:NOffSgn+lenof_sign_kp-1, det_ind)
                    if (IsUnoccDet(int_sign)) cycle
                    real_sign_2 = transfer(int_sign, real_sign_1)

                    ! Finally, add in the contribution to the projected Hamiltonian for each pair of Krylov vectors.
                    do i = 1, kp%nvecs
                        do j = i, kp%nvecs
#ifdef __DOUBLERUN
                            kp%hamil_matrix(i,j) = kp%hamil_matrix(i,j) + &
                                (real_sign_1(2*i-1)*real_sign_2(2*j) + &
                                real_sign_1(2*i)*real_sign_2(2*j-1))/2.0_dp
#else
                            kp%hamil_matrix(i,j) = kp%hamil_matrix(i,j) + real_sign_1(i)*real_sign_2(j)
#endif
                        end do
                    end do
                end if
            end if
        end do

    end subroutine calc_hamil_contribs_spawn

    subroutine calc_hamil_contribs_diag(kp)
    
        use FciMCData, only: determ_proc_sizes
        use CalcData, only: tSemiStochastic

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, i, j, min_idet, hdiag_ind
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign(lenof_sign_kp)
        real(dp) :: h_diag

        ! In semi-stochastic calculations the diagonal elements of the core space are
        ! taken care of in the core Hamiltonian calculation, so don't skip them here.
        ! Core determinants are always kept at the top of the list, so they're simple
        ! to skip.
        if (tSemiStochastic) then
            min_idet = determ_proc_sizes(iProcIndex) + 1
        else
            min_idet = 1
        end if

        hdiag_ind = NIfDBO + lenof_sign_kp + 1

        do idet = min_idet, TotWalkersKp
            int_sign = krylov_vecs(NOffSgn:NOffSgn+lenof_sign_kp-1, idet)
            real_sign = transfer(int_sign, real_sign)
            h_diag = transfer(krylov_vecs(hdiag_ind, idet), h_diag) + Hii

            ! Finally, add in the contribution to the projected Hamiltonian for each pair of Krylov vectors.
            do i = 1, kp%nvecs
                do j = i, kp%nvecs
#ifdef __DOUBLERUN
                    kp%hamil_matrix(i,j) = kp%hamil_matrix(i,j) + &
                        h_diag*(real_sign(2*i-1)*real_sign(2*j) + &
                        real_sign(2*i)*real_sign(2*j-1))/2.0_dp
#else
                    kp%hamil_matrix(i,j) = kp%hamil_matrix(i,j) + h_diag*real_sign(i)*real_sign(j)
#endif
                end do
            end do

        end do

    end subroutine calc_hamil_contribs_diag

    subroutine calc_hamil_contribs_semistoch(kp)
    
        use FciMCData, only: partial_determ_vecs_kp

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, i, j
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign(lenof_sign_kp)

        do idet = 1, determ_proc_sizes(iProcIndex)
            int_sign = krylov_vecs(NOffSgn:NOffSgn+lenof_sign_kp-1, idet)
            real_sign = transfer(int_sign, real_sign)

            do i = 1, kp%nvecs
                do j = i, kp%nvecs
#ifdef __DOUBLERUN
                    kp%hamil_matrix(i,j) = kp%hamil_matrix(i,j) + &
                        (real_sign(2*i-1)*partial_determ_vecs_kp(2*j, idet) + &
                        real_sign(2*i)*partial_determ_vecs_kp(2*j-1, idet))/2.0_dp
#else
                    kp%hamil_matrix(i,j) = kp%hamil_matrix(i,j) + real_sign(i)*partial_determ_vecs_kp(j, idet)
#endif
                end do
            end do

        end do

    end subroutine calc_hamil_contribs_semistoch

end module kp_fciqmc_hamil
