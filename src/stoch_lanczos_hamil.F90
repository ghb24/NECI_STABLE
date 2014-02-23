#include "macros.h"

module stoch_lanczos_hamil

    use bit_rep_data, only: NIfTot, NIfDBO, NOffFlag
    use stoch_lanczos_procs
    use SystemData, only: nel
    implicit none

contains

    subroutine calc_projected_hamil(lanczos)

        use CalcData, only: tSemiStochastic, AvMCExcits
        use constants
        use FciMCData, only: fcimc_excit_gen_store, exFlag, partial_determ_vecs_sl
        use procedure_pointers, only: generate_excitation, encode_child, get_spawn_helement
        use semi_stoch_procs, only: is_core_state, check_determ_flag
        use semi_stoch_procs, only: deterministic_projection_sl_hamil
        
        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer :: idet, ispawn, nspawn, nspawns_this_proc, i, j
        integer :: determ_ind, flag_ind, ic, ex(2,2)
        integer :: nI_parent(nel), nI_child(nel)
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        real(dp) :: prob
        real(dp), dimension(lenof_sign_sl) :: child_sign, parent_sign
        integer(n_int) :: int_sign(lenof_sign_sl)
        logical :: tChildIsDeterm, tParentIsDeterm, tParity
        HElement_t :: HElGen, HEl

        lanczos%hamil_matrix = 0.0_dp

        ilut_parent = 0_n_int
        if (tUseFlags) flag_ind = NIfDBO + lenof_sign_sl + 2
        determ_ind = 1

        do idet = 1, TotWalkersLanczos

            ! The 'parent' determinant from which spawning is to be attempted.
            ilut_parent(0:NIfDBO) = lanczos_vecs(0:NIfDBO,idet)
            if (tUseFlags) ilut_parent(NOffFlag) = lanczos_vecs(flag_ind,idet)

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call decode_bit_det(nI_parent, ilut_parent)
            int_sign = lanczos_vecs(NIfDBO+1:NIfDBO+lenof_sign_sl, idet)
            parent_sign = transfer(int_sign, parent_sign)

            tParentIsDeterm = check_determ_flag(ilut_parent)

            ! If this determinant is in the deterministic space then store the relevant
            ! data in arrays for later use.
            if (tParentIsDeterm) then
                ! Add the amplitude to the deterministic vector.
                partial_determ_vecs_sl(:,determ_ind) = parent_sign
                determ_ind = determ_ind + 1
            end if

            nspawn = decide_num_to_spawn_sl_hamil(parent_sign, AvMCExcits)
            
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

                    child_sign = return_amp_sl_hamil(parent_sign, AvMCExcits, HEl)
                else
                    child_sign = 0.0_dp
                end if

                ! If any (valid) children have been spawned.
                if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                    call create_particle_sl_hamil (nI_child, ilut_child, child_sign)

                end if ! If a child was spawned.

            end do ! Over mulitple spawns from the same determinant.

        end do ! Over all determinants.

        call distribute_spawns_sl_hamil(nspawns_this_proc)

        call calc_hamil_contribs_spawn(lanczos, nspawns_this_proc)

        call calc_hamil_contribs_diag(lanczos)

        if (tSemiStochastic) then
            call deterministic_projection_sl_hamil()
            call calc_hamil_contribs_semistoch(lanczos)
        end if

        ! Symmetrise the projected Hamiltonian.
        do i = 1, lanczos%nvecs
            do j = 1, i-1
                lanczos%hamil_matrix(i,j) = lanczos%hamil_matrix(j,i)
            end do
        end do

    end subroutine calc_projected_hamil

    function decide_num_to_spawn_sl_hamil(parent_sign, av_spawns_per_walker) result(nspawn)

        real(dp), intent(in) :: parent_sign(lenof_sign_sl)
        real(dp), intent(in) :: av_spawns_per_walker
        integer :: nspawn
        real(dp) :: unrounded_nspawn, prob_extra_walker, r

        unrounded_nspawn = sum(abs(parent_sign*av_spawns_per_walker))
        nspawn = int(unrounded_nspawn)
        if (unrounded_nspawn - real(nspawn,dp) > 0) then
            prob_extra_walker = unrounded_nspawn - real(nspawn,dp)
            r = genrand_real2_dSFMT()
            if (prob_extra_walker > r) nspawn = nspawn + 1
        end if

    end function decide_num_to_spawn_sl_hamil

    function return_amp_sl_hamil(parent_sign, nattempts, HEl) result(child_sign)

        real(dp), intent(in) :: parent_sign(lenof_sign_sl)
        real(dp), intent(in) :: nattempts
        HElement_t, intent(in) :: HEl
        real(dp) :: child_sign(lenof_sign_sl)
        real(dp) :: Hel_real, prob

        HEl_real = real(HEl, dp)
        prob = prob*nattempts
        child_sign = HEl*parent_sign/prob

    end function return_amp_sl_hamil

    subroutine create_particle_sl_hamil (nI_child, ilut_child, child_sign)

        use hash, only: DetermineDetNode
        use FciMCData, only: ValidSpawnedList

        integer, intent(in) :: nI_child(nel)
        integer, intent(in) :: ilut_child(0:NIfTot)
        real(dp), intent(in) :: child_sign(lenof_sign_sl)
        integer :: proc
        integer :: int_sign(lenof_sign_sl)

        proc = DetermineDetNode(nI_child, 0)

        SpawnedPartsLanc(0:NIfDBO, ValidSpawnedList(proc)) = ilut_child(0:NIfDBO)
        int_sign = transfer(child_sign, int_sign)
        SpawnedPartsLanc(NOffSgn:NOffSgn+lenof_sign_sl-1, ValidSpawnedList(proc)) = int_sign

        ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1

    end subroutine create_particle_sl_hamil

    subroutine distribute_spawns_sl_hamil(nspawns_this_proc)

        use AnnihilationMod, only: SendProcNewParts
        use FciMCData, only: SpawnedParts, SpawnedParts2, SpawnedPartsLanc, SpawnedPartsLanc2

        integer, intent(out) :: nspawns_this_proc
        integer(n_int), pointer :: PointTemp(:,:), PointTemp2(:,:)

        ! We want SpawnedParts and SpawnedParts2 to point to the 'wider' spawning arrays which
        ! hold all the signs of *all* the Krylov vectors. This is because SendProcNewParts uses
        ! SpawnedParts and SpawnedParts2.
        PointTemp => SpawnedParts
        PointTemp2 => SpawnedParts2
        SpawnedParts => SpawnedPartsLanc
        SpawnedParts2 => SpawnedPartsLanc2

        call SendProcNewParts(nspawns_this_proc, .false.)

        ! Now we want SpawnedPartsLanc to point to SpawnedParts2, which holds the output of
        ! the communication of the spawns.
        SpawnedPartsLanc => SpawnedParts2
        SpawnedPartsLanc2 => SpawnedParts

        ! Now let SpawnedParts1 and SpawnedParts2 to point to their original arrays.
        SpawnedParts => PointTemp
        SpawnedParts2 => PointTemp2

        nullify(PointTemp)
        nullify(PointTemp2)

    end subroutine distribute_spawns_sl_hamil

    subroutine calc_hamil_contribs_spawn(lanczos, nspawns_this_proc)

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer, intent(in) :: nspawns_this_proc
        integer :: idet, DetHash, det_ind, i, j
        integer(n_int) :: ilut_spawn(0:NIfTot)
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_sign_sl)
        real(dp) :: real_sign_1(lenof_sign_sl), real_sign_2(lenof_sign_sl)
        type(ll_node), pointer :: temp_node
        logical :: tDetFound

        ilut_spawn = 0_n_int

        do idet = 1, nspawns_this_proc
            ilut_spawn(0:NIfDBO) = SpawnedParts(0:NIfDBO, idet)
            int_sign = SpawnedParts(NOffSgn:NOffSgn+lenof_sign_sl-1, idet)
            real_sign_1 = transfer(int_sign, real_sign_1)
            call decode_bit_det(nI_spawn, ilut_spawn)
            DetHash = FindWalkerHash(nI_spawn, nhashes_lanczos)
            ! Point to the first node with this hash value in lanczos_vecs.
            temp_node => lanczos_hash_table(DetHash)
            if (temp_node%ind == 0) then
                ! If there are no determinants at all with this hash value in lanczos_vecs.
                cycle
            else
                tDetFound = .false.
                do while (associated(temp_node))
                    if (DetBitEQ(ilut_spawn, lanczos_vecs(:,temp_node%ind), NIfDBO)) then
                        ! If this CurrentDets determinant has been found in lanczos_vecs.
                        det_ind = temp_node%ind
                        tDetFound = .true.
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    temp_node => temp_node%next
                end do
                if (tDetFound) then
                    int_sign = lanczos_vecs(NOffSgn:NOffSgn+lenof_sign_sl-1, det_ind)
                    if (IsUnoccDet(int_sign)) cycle
                    real_sign_2 = transfer(int_sign, real_sign_1)

                    ! Finally, add in the contribution to the projected Hamiltonian for each pair of Krylov vectors.
                    do i = 1, lanczos%nvecs
                        do j = i, lanczos%nvecs
                            lanczos%hamil_matrix(i,j) = lanczos%hamil_matrix(i,j) + &
                                (real_sign_1(2*i-1)*real_sign_2(2*j) + &
                                real_sign_1(2*i)*real_sign_2(2*j-1))/2
                        end do
                    end do

                end if
            end if
        end do

    end subroutine calc_hamil_contribs_spawn

    subroutine calc_hamil_contribs_diag(lanczos)
    
        use FciMCData, only: determ_proc_sizes
        use CalcData, only: tSemiStochastic

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer :: idet, i, j, min_idet, hdiag_ind
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_sign_sl)
        real(dp) :: real_sign(lenof_sign_sl)
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

        hdiag_ind = NIfDBO + lenof_sign_sl + 1

        do idet = min_idet, TotWalkersLanczos
            int_sign = lanczos_vecs(NOffSgn:NOffSgn+lenof_sign_sl-1, idet)
            real_sign = transfer(int_sign, real_sign)
            h_diag = transfer(lanczos_vecs(hdiag_ind, idet), h_diag) + Hii

            ! Finally, add in the contribution to the projected Hamiltonian for each pair of Krylov vectors.
            do i = 1, lanczos%nvecs
                do j = i, lanczos%nvecs
                    lanczos%hamil_matrix(i,j) = lanczos%hamil_matrix(i,j) + &
                        h_diag*(real_sign(2*i-1)*real_sign(2*j) + &
                        real_sign(2*i)*real_sign(2*j-1))/2
                end do
            end do

        end do

    end subroutine calc_hamil_contribs_diag

    subroutine calc_hamil_contribs_semistoch(lanczos)
    
        use FciMCData, only: partial_determ_vecs_sl

        type(stoch_lanczos_data), intent(inout) :: lanczos
        integer :: idet, i, j
        integer :: nI_spawn(nel)
        integer(n_int) :: int_sign(lenof_sign_sl)
        real(dp) :: real_sign(lenof_sign_sl)

        do idet = 1, determ_proc_sizes(iProcIndex)
            int_sign = lanczos_vecs(NOffSgn:NOffSgn+lenof_sign_sl-1, idet)
            real_sign = transfer(int_sign, real_sign)

            do i = 1, lanczos%nvecs
                do j = i, lanczos%nvecs
                    lanczos%hamil_matrix(i,j) = lanczos%hamil_matrix(i,j) + &
                        (real_sign(2*i-1)*partial_determ_vecs_sl(2*j, idet) + &
                        real_sign(2*i)*partial_determ_vecs_sl(2*j-1, idet))/2
                end do
            end do

        end do

    end subroutine calc_hamil_contribs_semistoch

end module stoch_lanczos_hamil
