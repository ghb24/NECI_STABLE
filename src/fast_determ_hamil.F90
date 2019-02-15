#include "macros.h"

module fast_determ_hamil

    use bit_rep_data, only: NIfTot, NIfDBO, NIfD
    use bit_reps, only: decode_bit_det
    use constants
    use DetBitOps, only: CountBits, TestClosedShellDet, CalcOpenOrbs
    use Determinants, only: get_helement
    use FciMCData, only: determ_space_size, determ_sizes, determ_displs, &
                         Hii, core_ham_diag, core_space
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_off_diag_helement_opt, hphf_off_diag_special_case
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier
    use sltcnd_mod, only: sltcnd_knowIC
    use sparse_arrays
    use SystemData, only: tHPHF, nel
    use timing_neci
    use util_mod, only: get_free_unit

    implicit none

    type auxiliary_array
        integer, allocatable :: pos(:)
    end type auxiliary_array

contains

    subroutine calc_determ_hamil_opt_hphf()

        use DetBitOps, only: MaskAlpha, MaskBeta
        use FciMCData, only: ll_node, s_first_ind, s_last_ind
        use hash, only: init_hash_table, clear_hash_table
        use hash, only: hash_table_lookup, add_hash_table_entry
        use sort_mod
        use SystemData, only: nOccAlpha, nOccBeta

        integer :: i, j, k, ierr
        integer :: IC, IC_beta
        integer :: ind_i, ind_j, ind_k, i_full
        integer :: ind, hash_val
        integer :: ind_beta, ind_alpha, hash_val_beta, hash_val_alpha, ind_alpha_conn
        integer :: ind_alpha_paired, ind_beta_paired
        logical :: tSuccess, tSuccess_a_paired, tSuccess_b_paired
        integer(n_int) :: tmp(0:NIfD)

        integer(n_int) :: alpha_ilut(0:NIfTot), alpha_m1_ilut(0:NIfTot), beta_ilut(0:NIfTot), beta_m1_ilut(0:NIfTot)
        integer(n_int) :: ilut_full(0:NIfTot), ilut_paired(0:NIfTot), beta_ilut_1(0:NIfD), beta_ilut_2(0:NIfD)
        integer(n_int), allocatable :: beta_list(:,:), alpha_list(:,:), alpha_m1_list(:,:), beta_m1_list(:,:)
        integer :: nI_alpha(nOccAlpha), nI_alpha_m1(nOccAlpha-1), nI_beta(nOccBeta), nI_beta_m1(nOccBeta-1)
        integer :: nI(nel), nI_paired(nel)
        integer :: nbeta, nalpha, nbeta_m1, nalpha_m1
        ! The number of beta and alpha strings in the 'unpaired' determinants
        ! in the HPHFs
        integer :: nbeta_unpaired, nalpha_unpaired

        integer :: nintersec
        integer, allocatable :: intersec_inds(:)

        type(ll_node), pointer :: beta_ht(:)
        type(ll_node), pointer :: alpha_ht(:)
        type(ll_node), pointer :: beta_m1_ht(:)
        type(ll_node), pointer :: alpha_m1_ht(:)

        integer, allocatable :: nbeta_dets(:), nalpha_dets(:), nbeta_m1_contribs(:), nalpha_m1_contribs(:)
        type(auxiliary_array), allocatable :: beta_dets(:)
        type(auxiliary_array), allocatable :: alpha_dets(:)
        type(auxiliary_array), allocatable :: beta_m1_contribs(:)
        type(auxiliary_array), allocatable :: alpha_m1_contribs(:)
        type(auxiliary_array), allocatable :: beta_with_alpha(:)

        integer, allocatable :: nalpha_alpha(:), nbeta_beta(:)
        type(auxiliary_array), allocatable :: alpha_alpha(:), beta_beta(:)

        HElement_t(dp) :: hel, hel_paired
        HElement_t(dp), allocatable :: hamil_row(:)
        integer, allocatable :: hamil_pos(:)
        integer, allocatable :: num_conns(:)

        ! HPHF arrays
        logical :: CS_I, CS_J, CS_K
        logical, allocatable :: cs(:)
        integer :: num_hphfs, hphf_ind
        integer :: OpenOrbsI, OpenOrbsJ, OpenOrbsK

        type(timer), save :: aux_time
        type(timer), save :: sort_aux_time
        type(timer), save :: ham_time
        type(timer), save :: sort_ham_time
        real(dp) :: total_time

        character(len=*), parameter :: t_r = "calc_determ_hamil_opt_hphf"

        allocate(beta_list(0:NIfD, 2*determ_space_size), stat=ierr)
        allocate(nbeta_dets(determ_space_size), stat=ierr)
        allocate(alpha_list(0:NIfD, 2*determ_space_size), stat=ierr)
        allocate(nalpha_dets(determ_space_size), stat=ierr)

        allocate(beta_ht( int(determ_space_size/10.0) ), stat=ierr)
        call init_hash_table(beta_ht)
        allocate(alpha_ht( int(determ_space_size/10.0) ), stat=ierr)
        call init_hash_table(alpha_ht)

        call set_timer(aux_time)

        ! --- Set up auxiliary arrays ------------------------------

        ! The number of different alpha and beta strings among the determinants
        nbeta = 0
        nalpha = 0
        ! The number of alpha strings with nel-1 electrons among the dets
        nbeta_m1 = 0
        nalpha_m1 = 0

        allocate(cs(determ_space_size), stat=ierr)
        do i = 1, determ_space_size
            cs(i) = TestClosedShellDet(core_space(:,i))
        end do

        ! First time through, just count the number entries in each of the
        ! lists to be set up. Then allocate those lists and fill them in.
        do hphf_ind = 1, 2

            if (hphf_ind == 2) then
                nbeta_unpaired = nbeta
                nalpha_unpaired = nalpha
            end if

            do i = 1, determ_space_size

                if (hphf_ind == 1) then
                    ilut_full(0:NIfD) = core_space(0:NIfD,i)
                else
                    ! For the second loop through, only looks at the determinants
                    ! 'paired' with the stored determinant in the HPHF.
                    if (cs(i)) then
                        ! For a closed shell determinant, there is no paired det
                        ! to consider.
                        cycle
                    else
                        call FindExcitBitDetSym(core_space(:,i), ilut_full)
                    end if
                end if

                ! --- Beta string -----------------------------
                beta_ilut(0:NIfD) = iand(ilut_full(0:NIfD), MaskBeta)
                call decode_bit_det(nI_beta, beta_ilut)
                call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)

                ! If this beta string has not been generated already, add it to
                ! the list of all beta strings, and the associated hash table.
                if (.not. tSuccess) then
                    nbeta = nbeta + 1
                    nbeta_dets(nbeta) = 0
                    beta_list(0:NIfD, nbeta) = beta_ilut(0:NIfD)
                    call add_hash_table_entry(beta_ht, nbeta, hash_val_beta)
                    ind_beta = nbeta
                end if

                ! Now add this determinant to the list of determinants with this
                ! beta string.
                ! Only store the 'unpaired' dets in this array
                if (hphf_ind == 1) nbeta_dets(ind_beta) = nbeta_dets(ind_beta) + 1

                ! --- Alpha string ---------------------------
                alpha_ilut(0:NIfD) = iand(ilut_full(0:NIfD), MaskAlpha)
                call decode_bit_det(nI_alpha, alpha_ilut)
                call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)

                ! If this alpha string has not been generated already, add it to
                ! the list of all alpha strings, and the associated hash table.
                if (.not. tSuccess) then
                    nalpha = nalpha + 1
                    nalpha_dets(nalpha) = 0
                    alpha_list(0:NIfD, nalpha) = alpha_ilut(0:NIfD)
                    call add_hash_table_entry(alpha_ht, nalpha, hash_val_alpha)
                    ind_alpha = nalpha
                end if

                ! Now add this determinant to the list of determinants with this
                ! alpha string.
                ! Only store the 'unpaired' dets in this array
                if (hphf_ind == 1) nalpha_dets(ind_alpha) = nalpha_dets(ind_alpha) + 1

            end do
        end do 

        allocate(beta_m1_list(0:NIfD, nbeta*(nOccBeta-1)), stat=ierr)
        allocate(nbeta_m1_contribs(nbeta*(nOccBeta-1)), stat=ierr)
        allocate(alpha_m1_list(0:NIfD, nalpha*(nOccAlpha-1)), stat=ierr)
        allocate(nalpha_m1_contribs(nalpha*(nOccAlpha-1)), stat=ierr)

        allocate(beta_m1_ht( int(nbeta*(nOccBeta-1)/10.0) ), stat=ierr)
        call init_hash_table(beta_m1_ht)
        allocate(alpha_m1_ht( int(nalpha*(nOccAlpha-1)/10.0) ), stat=ierr)
        call init_hash_table(alpha_m1_ht)

        ! --- Find the size of beta-1 arrays -----------------

        do i = 1, nbeta
            beta_ilut(0:NIfD) = beta_list(0:NIfD,i)
            call decode_bit_det(nI_beta, beta_ilut)

            do j = 1, nOccBeta
                ! Create the beta-1 orbitals and ilut
                if (j > 1) nI_beta_m1(1:j-1) = nI_beta(1:j-1)
                if (j < nOccBeta) nI_beta_m1(j:nOccBeta-1) = nI_beta(j+1:nOccBeta)
                beta_m1_ilut = beta_ilut
                clr_orb(beta_m1_ilut, nI_beta(j))

                call hash_table_lookup(nI_beta_m1, beta_m1_ilut, NIfD, beta_m1_ht, beta_m1_list, ind, hash_val, tSuccess)

                ! If this beta-1 string has not been generated already, add it to
                ! the list of all beta-1 strings, and the associated hash table.
                if (.not. tSuccess) then
                    nbeta_m1 = nbeta_m1 + 1
                    beta_m1_list(0:NIfD, nbeta_m1) = beta_m1_ilut(0:NIfD)
                    nbeta_m1_contribs(nbeta_m1) = 0
                    call add_hash_table_entry(beta_m1_ht, nbeta_m1, hash_val)
                    ind = nbeta_m1
                end if

                ! Now add this determinant to the list of determinants with this
                ! beta-1 string.
                nbeta_m1_contribs(ind) = nbeta_m1_contribs(ind) + 1
            end do
        end do

        ! --- Find the size of alpha-1 arrays -----------------

        do i = 1, nalpha
            alpha_ilut(0:NIfD) = alpha_list(0:NIfD,i)
            call decode_bit_det(nI_alpha, alpha_ilut)

            do j = 1, nOccAlpha
                ! Create the alpha-1 orbitals and ilut
                if (j > 1) nI_alpha_m1(1:j-1) = nI_alpha(1:j-1)
                if (j < nOccAlpha) nI_alpha_m1(j:nOccAlpha-1) = nI_alpha(j+1:nOccAlpha)
                alpha_m1_ilut = alpha_ilut
                clr_orb(alpha_m1_ilut, nI_alpha(j))

                call hash_table_lookup(nI_alpha_m1, alpha_m1_ilut, NIfD, alpha_m1_ht, alpha_m1_list, ind, hash_val, tSuccess)

                ! If this alpha-1 string has not been generated already, add it to
                ! the list of all alpha-1 strings, and the associated hash table.
                if (.not. tSuccess) then
                    nalpha_m1 = nalpha_m1 + 1
                    alpha_m1_list(0:NIfD, nalpha_m1) = alpha_m1_ilut(0:NIfD)
                    nalpha_m1_contribs(nalpha_m1) = 0
                    call add_hash_table_entry(alpha_m1_ht, nalpha_m1, hash_val)
                    ind = nalpha_m1
                end if

                ! Now add this determinant to the list of determinants with this
                ! alpha-1 string.
                nalpha_m1_contribs(ind) = nalpha_m1_contribs(ind) + 1
            end do
        end do

        ! Now we know the size of the auxiliary arrays

        ! --- Allocate the auxiliary arrays ---------------
        allocate(beta_dets(nbeta_unpaired), stat=ierr)
        allocate(alpha_dets(nalpha_unpaired), stat=ierr)
        allocate(beta_with_alpha(nalpha_unpaired), stat=ierr)

        allocate(beta_m1_contribs(nbeta_m1), stat=ierr)
        allocate(alpha_m1_contribs(nalpha_m1), stat=ierr)

        do i = 1, nbeta_unpaired
            allocate(beta_dets(i)%pos(nbeta_dets(i)), stat=ierr)
        end do
        do i = 1, nalpha_unpaired
            allocate(alpha_dets(i)%pos(nalpha_dets(i)), stat=ierr)
        end do
        do i = 1, nalpha_unpaired
            allocate(beta_with_alpha(i)%pos(nalpha_dets(i)), stat=ierr)
        end do
        do i = 1, nbeta_m1
            allocate(beta_m1_contribs(i)%pos(nbeta_m1_contribs(i)), stat=ierr)
        end do
        do i = 1, nalpha_m1
            allocate(alpha_m1_contribs(i)%pos(nalpha_m1_contribs(i)), stat=ierr)
        end do

        ! --- Now fill the auxiliary arrays ----------------
        nbeta_dets(1:nbeta_unpaired) = 0
        nalpha_dets(1:nalpha_unpaired) = 0
        nbeta_m1_contribs(1:nbeta_m1) = 0
        nalpha_m1_contribs(1:nalpha_m1) = 0

        do i = 1, determ_space_size
            ! If this is a closed shell HPHF then only one determinant contributes.
            ! Otherwise, two do, and they both must be considered.
            num_hphfs = 1

            do hphf_ind = 1, num_hphfs

                ilut_full(0:NIfD) = core_space(0:NIfD,i)

                ! --- Beta string -----------------------------
                beta_ilut(0:NIfD) = iand(ilut_full(0:NIfD), MaskBeta)
                call decode_bit_det(nI_beta, beta_ilut)
                call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)

                ! --- Alpha string -----------------------------
                alpha_ilut(0:NIfD) = iand(ilut_full(0:NIfD), MaskAlpha)
                call decode_bit_det(nI_alpha, alpha_ilut)
                call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)

                ! Now add this determinant to the list of determinants with this
                ! beta string.
                nbeta_dets(ind_beta) = nbeta_dets(ind_beta) + 1
                beta_dets(ind_beta)%pos(nbeta_dets(ind_beta)) = i

                ! Now add this determinant to the list of determinants with this
                ! alpha string.
                nalpha_dets(ind_alpha) = nalpha_dets(ind_alpha) + 1
                alpha_dets(ind_alpha)%pos(nalpha_dets(ind_alpha)) = i
                beta_with_alpha(ind_alpha)%pos(nalpha_dets(ind_alpha)) = ind_beta
            end do
        end do

        do i = 1, nbeta
            beta_ilut(0:NIfD) = beta_list(0:NIfD,i)
            call decode_bit_det(nI_beta, beta_ilut)

            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)
            ! --- Beta-1 string ----------------------------
            do j = 1, nOccBeta
                ! Create the beta-1 orbitals and ilut
                if (j > 1) nI_beta_m1(1:j-1) = nI_beta(1:j-1)
                if (j < nOccBeta) nI_beta_m1(j:nOccBeta-1) = nI_beta(j+1:nOccBeta)
                beta_m1_ilut = beta_ilut
                clr_orb(beta_m1_ilut, nI_beta(j))

                call hash_table_lookup(nI_beta_m1, beta_m1_ilut, NIfD, beta_m1_ht, beta_m1_list, ind, hash_val, tSuccess)

                if (.not. tSuccess) call stop_all(t_r, "beta-1 string not found.")

                ! Now add this determinant to the list of determinants with this
                ! beta-1 string.
                nbeta_m1_contribs(ind) = nbeta_m1_contribs(ind) + 1
                beta_m1_contribs(ind)%pos(nbeta_m1_contribs(ind)) = ind_beta
            end do
        end do

        do i = 1, nalpha
            alpha_ilut(0:NIfD) = alpha_list(0:NIfD,i)
            call decode_bit_det(nI_alpha, alpha_ilut)

            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)
            ! --- Alpha-1 string ---------------------------
            do j = 1, nOccAlpha
                ! Create the alpha-1 orbitals and ilut
                if (j > 1) nI_alpha_m1(1:j-1) = nI_alpha(1:j-1)
                if (j < nOccAlpha) nI_alpha_m1(j:nOccAlpha-1) = nI_alpha(j+1:nOccAlpha)
                alpha_m1_ilut = alpha_ilut
                clr_orb(alpha_m1_ilut, nI_alpha(j))

                call hash_table_lookup(nI_alpha_m1, alpha_m1_ilut, NIfD, alpha_m1_ht, alpha_m1_list, ind, hash_val, tSuccess)

                if (.not. tSuccess) call stop_all(t_r, "alpha-1 string not found.")

                ! Now add this determinant to the list of determinants with this
                ! alpha-1 string.
                nalpha_m1_contribs(ind) = nalpha_m1_contribs(ind) + 1
                alpha_m1_contribs(ind)%pos(nalpha_m1_contribs(ind)) = ind_alpha
            end do
        end do

        deallocate(beta_m1_list)
        deallocate(alpha_m1_list)

        call clear_hash_table(beta_m1_ht)
        deallocate(beta_m1_ht, stat=ierr)
        nullify(beta_m1_ht)

        call clear_hash_table(alpha_m1_ht)
        deallocate(alpha_m1_ht, stat=ierr)
        nullify(alpha_m1_ht)

        allocate(nbeta_beta(nbeta), stat=ierr)
        allocate(nalpha_alpha(nalpha), stat=ierr)
        allocate(beta_beta(nbeta), stat=ierr)
        allocate(alpha_alpha(nalpha), stat=ierr)

        nbeta_beta = 0
        nalpha_alpha = 0

        ! Find the size of the beta_beta array to be created
        do i = 1, nbeta_m1
            do j = 1, nbeta_m1_contribs(i)
                !nbeta_beta(beta_m1_contribs(i)%pos(j)) = nbeta_beta(beta_m1_contribs(i)%pos(j)) + nbeta_m1_contribs(i) - 1
                do k = 1, nbeta_m1_contribs(i)
                    if (j == k) cycle
                    ! Don't need to consider beta values that aren't in the unpaired detereminants
                    if (beta_m1_contribs(i)%pos(k) > nbeta_unpaired) cycle
                    nbeta_beta(beta_m1_contribs(i)%pos(j)) = nbeta_beta(beta_m1_contribs(i)%pos(j)) + 1
                end do
            end do
        end do

        ! Allocate the beta_beta array...
        do i = 1, nbeta
            allocate(beta_beta(i)%pos(nbeta_beta(i)), stat=ierr)
        end do

        ! Rezero this so that we can use it as a counter for the following
        nbeta_beta = 0

        ! ...and finally fill the beta_beta array.
        do i = 1, nbeta_m1
            do j = 1, nbeta_m1_contribs(i)
                do k = 1, nbeta_m1_contribs(i)
                    if (j == k) cycle
                    ! Don't need to consider beta values that aren't in the unpaired detereminants
                    if (beta_m1_contribs(i)%pos(k) > nbeta_unpaired) cycle
                    nbeta_beta(beta_m1_contribs(i)%pos(j)) = nbeta_beta(beta_m1_contribs(i)%pos(j)) + 1
                    beta_beta( beta_m1_contribs(i)%pos(j) )%pos(nbeta_beta( beta_m1_contribs(i)%pos(j) )) = beta_m1_contribs(i)%pos(k)
                end do
            end do
        end do

        ! Find the size of the alpha_alpha array to be created
        do i = 1, nalpha_m1
            do j = 1, nalpha_m1_contribs(i)
                !nalpha_alpha(alpha_m1_contribs(i)%pos(j)) = nalpha_alpha(alpha_m1_contribs(i)%pos(j)) + nalpha_m1_contribs(i) - 1
                do k = 1, nalpha_m1_contribs(i)
                    if (j == k) cycle
                    ! Don't need to consider alpha values that aren't in the unpaired detereminants
                    if (alpha_m1_contribs(i)%pos(k) > nalpha_unpaired) cycle
                    nalpha_alpha(alpha_m1_contribs(i)%pos(j)) = nalpha_alpha(alpha_m1_contribs(i)%pos(j)) + 1
                end do
            end do
        end do

        ! Allocate the alpha_alpha array...
        do i = 1, nalpha
            allocate(alpha_alpha(i)%pos(nalpha_alpha(i)), stat=ierr)
        end do

        ! Rezero this so that we can use it as a counter for the following
        nalpha_alpha = 0

        ! ...and finally fill the alpha_alpha array.
        do i = 1, nalpha_m1
            do j = 1, nalpha_m1_contribs(i)
                do k = 1, nalpha_m1_contribs(i)
                    if (j == k) cycle
                    ! Don't need to consider alpha values that aren't in the unpaired detereminants
                    if (alpha_m1_contribs(i)%pos(k) > nalpha_unpaired) cycle
                    nalpha_alpha(alpha_m1_contribs(i)%pos(j)) = nalpha_alpha(alpha_m1_contribs(i)%pos(j)) + 1
                    alpha_alpha( alpha_m1_contribs(i)%pos(j) )%pos(nalpha_alpha( alpha_m1_contribs(i)%pos(j) )) = alpha_m1_contribs(i)%pos(k)

                end do
            end do
        end do

        call halt_timer(aux_time)

        write(6,'("Time to create auxiliary arrays:", f9.3)') get_total_time(aux_time); call neci_flush(6)

        ! Sort auxiliary arrays into the required order
        call set_timer(sort_aux_time)

        do i = 1, nalpha_unpaired
            call sort(beta_with_alpha(i)%pos, alpha_dets(i)%pos)
        end do

        do i = 1, nbeta
            call sort(beta_beta(i)%pos)
        end do

        do i = 1, nalpha
            call sort(alpha_alpha(i)%pos)
        end do

        call halt_timer(sort_aux_time)

        ! Deallocate the arrays which we don't need for the following.

        do i = nbeta_m1, 1, -1
            deallocate(beta_m1_contribs(i)%pos, stat=ierr)
        end do
        do i = nalpha_m1, 1, -1
            deallocate(alpha_m1_contribs(i)%pos, stat=ierr)
        end do

        deallocate(nbeta_m1_contribs)
        deallocate(nalpha_m1_contribs)
        deallocate(beta_m1_contribs)
        deallocate(alpha_m1_contribs)

        write(6,'("Time to sort auxiliary arrays:", f9.3)') get_total_time(sort_aux_time); call neci_flush(6)

        ! Actually create the Hamiltonian
        call set_timer(ham_time)

        allocate(intersec_inds(nbeta), stat=ierr)

        allocate(hamil_row(determ_space_size), stat=ierr)
        allocate(hamil_pos(determ_space_size), stat=ierr)

        allocate(sparse_core_ham(determ_sizes(iProcIndex)), stat=ierr)
        allocate(core_ham_diag(determ_sizes(iProcIndex)), stat=ierr)

        allocate(num_conns(determ_space_size), stat=ierr)
        num_conns = 0

        ! Loop over the determinants on this process
        do i = 1, determ_sizes(iProcIndex)
            ! Find the index in the *full* list of determinants
            i_full = i + determ_displs(iProcIndex)
            CS_I = cs(i_full)

            call decode_bit_det(nI, core_space(:,i_full))

            if (.not. CS_I) then
                call CalcOpenOrbs(core_space(:,i_full), OpenOrbsI)
                call FindExcitBitDetSym(core_space(:,i_full), ilut_paired)
                call decode_bit_det(nI_paired, ilut_paired)
            end if

            ! --- Beta string -----------------------------
            beta_ilut(0:NIfD) = iand(core_space(0:NIfD,i_full), MaskBeta)
            call decode_bit_det(nI_beta, beta_ilut)
            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)

            do j = 1, nbeta_dets(ind_beta)
                ind_j = beta_dets(ind_beta)%pos(j)
                if (i_full == ind_j) cycle
                CS_J = cs(ind_j)
                tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_j))
                tmp = iand(core_space(0:NIfD,i_full), tmp)
                IC = CountBits(tmp, NIfD)

                if (IC <= 2) then
                    hel = hphf_off_diag_helement_opt(nI, core_space(:,i_full), &
                                                     core_space(:,ind_j), IC, CS_I, CS_J)

                    if (abs(hel) > 0.0_dp) then
                        num_conns(i) = num_conns(i) + 1
                        hamil_pos(num_conns(i)) = ind_j
                        hamil_row(num_conns(i)) = hel
                    end if
                end if
            end do

            ! Now consider any open-open connections that would not have been found above
            if (.not. CS_I) then
                beta_ilut(0:NIfD) = iand(ilut_paired(0:NIfD), MaskBeta)
                call decode_bit_det(nI_beta, beta_ilut)
                call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, &
                                       ind_beta_paired, hash_val_beta, tSuccess_b_paired)

                ! If the beta string of the paired determinant is in one of the
                ! unpaired determinants - it might not be!
                !if (tSuccess_b_paired) then
                if (ind_beta_paired <= nbeta_unpaired) then
                    do j = 1, nbeta_dets(ind_beta_paired)
                        ind_j = beta_dets(ind_beta_paired)%pos(j)
                        ! This is a diagonal elements, included later
                        if (i_full == ind_j) cycle
                        CS_J = cs(ind_j)
                        ! We're only looking for open-open connections here
                        if (CS_J) cycle

                        ! First we need to make sure that the 'unpaired' determinant isn't
                        ! connected to this determinant. If it is then we will have included
                        ! it already, so don't want to include it again:
                        tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_j))
                        tmp = iand(core_space(0:NIfD,i_full), tmp)
                        IC = CountBits(tmp, NIfD)
                        if (IC <= 2) cycle

                        ! Finally, check if this paired determiannt is connected to this one
                        ! that we've just generated
                        tmp = ieor(ilut_paired, core_space(0:NIfD,ind_j))
                        tmp = iand(ilut_paired, tmp)
                        IC = CountBits(tmp, NIfD)
                        if (IC <= 2) then
                            hel = hphf_off_diag_special_case(nI_paired, ilut_paired, core_space(:,ind_j), IC, OpenOrbsI)
                            if (abs(hel) > 0.0_dp) then
                                num_conns(i) = num_conns(i) + 1
                                hamil_pos(num_conns(i)) = ind_j
                                hamil_row(num_conns(i)) = hel
                            end if
                        end if

                    end do
                end if
            end if

            ! --- Alpha string -----------------------------
            alpha_ilut(0:NIfD) = iand(core_space(0:NIfD,i_full), MaskAlpha)
            call decode_bit_det(nI_alpha, alpha_ilut)
            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)

            do j = 1, nalpha_dets(ind_alpha)
                ind_j = alpha_dets(ind_alpha)%pos(j)
                if (i_full == ind_j) cycle
                CS_J = cs(ind_j)
                tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_j))
                tmp = iand(core_space(0:NIfD,i_full), tmp)
                IC = CountBits(tmp, NIfD)

                if (IC <= 2) then
                    hel = hphf_off_diag_helement_opt(nI, core_space(:,i_full), &
                                                     core_space(:,ind_j), IC, CS_I, CS_J)
                    if (abs(hel) > 0.0_dp) then
                        num_conns(i) = num_conns(i) + 1
                        hamil_pos(num_conns(i)) = ind_j
                        hamil_row(num_conns(i)) = hel
                    end if
                end if

            end do

            ! Now consider any open-open connections that would not have been found above
            if (.not. CS_I) then
                alpha_ilut(0:NIfD) = iand(ilut_paired(0:NIfD), MaskAlpha)
                call decode_bit_det(nI_alpha, alpha_ilut)
                call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, &
                                       ind_alpha_paired, hash_val_alpha, tSuccess_a_paired)

                ! If the alpha string of the paired determinant is in one of the
                ! unpaired determinants - it might not be!
                !if (tSuccess_a_paired) then
                if (ind_alpha_paired <= nalpha_unpaired) then
                    do j = 1, nalpha_dets(ind_alpha_paired)
                        ind_j = alpha_dets(ind_alpha_paired)%pos(j)
                        ! This is a diagonal elements, included later
                        if (i_full == ind_j) cycle
                        CS_J = cs(ind_j)
                        ! We're only looking for open-open connections here
                        if (CS_J) cycle

                        ! First we need to make sure that the 'unpaired' determinant isn't
                        ! connected to this determinant. If it is then we will have included
                        ! it already, so don't want to include it again:
                        tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_j))
                        tmp = iand(core_space(0:NIfD,i_full), tmp)
                        IC = CountBits(tmp, NIfD)
                        if (IC <= 2) cycle

                        ! Finally, check if this paired determiannt is connected to this one
                        ! that we've just generated
                        tmp = ieor(ilut_paired, core_space(0:NIfD,ind_j))
                        tmp = iand(ilut_paired, tmp)
                        IC = CountBits(tmp, NIfD)
                        if (IC <= 2) then
                            hel = hphf_off_diag_special_case(nI_paired, ilut_paired, core_space(:,ind_j), IC, OpenOrbsI)
                            if (abs(hel) > 0.0_dp) then
                                num_conns(i) = num_conns(i) + 1
                                hamil_pos(num_conns(i)) = ind_j
                                hamil_row(num_conns(i)) = hel
                            end if
                        end if

                    end do
                end if
            end if

            ! Loop through alpha strings connected to this ind_alpha by a single excitation
            do j = 1, nalpha_alpha(ind_alpha)
                ! This is the index of the connected alpha string
                ind_alpha_conn = alpha_alpha(ind_alpha)%pos(j)

                call find_intersec( nalpha_dets(ind_alpha_conn), nbeta_beta(ind_beta), &
                               beta_with_alpha(ind_alpha_conn)%pos, beta_beta(ind_beta)%pos, &
                               intersec_inds, nintersec )

                do k = 1, nintersec
                    ind_k = alpha_dets(ind_alpha_conn)%pos( intersec_inds(k) )
                    CS_K = cs(ind_k)

                    hel = hphf_off_diag_helement_opt(nI, core_space(:,i_full), core_space(:,ind_k), 2, CS_I, CS_K)
                    if (abs(hel) > 0.0_dp) then
                        num_conns(i) = num_conns(i) + 1
                        hamil_pos(num_conns(i)) = ind_k
                        hamil_row(num_conns(i)) = hel
                    end if
                end do
            end do

            ! Now consider any open-open connections that would not have been found above
            if (.not. CS_I) then
                !if (tSuccess_a_paired .and. tSuccess_b_paired) then
                    do j = 1, nalpha_alpha(ind_alpha_paired)
                        ! This is the index of the connected alpha string
                        ind_alpha_conn = alpha_alpha(ind_alpha_paired)%pos(j)

                        call find_intersec( nalpha_dets(ind_alpha_conn), nbeta_beta(ind_beta_paired), &
                                       beta_with_alpha(ind_alpha_conn)%pos, beta_beta(ind_beta_paired)%pos, &
                                       intersec_inds, nintersec )

                        do k = 1, nintersec
                            ind_k = alpha_dets(ind_alpha_conn)%pos( intersec_inds(k) )
                            CS_K = cs(ind_k)
                            ! We're only looking for open-open connections here
                            if (CS_K) cycle

                            ! First we need to make sure that the 'unpaired' determinant isn't
                            ! connected to this determinant. If it is then we will have included
                            ! it already, so don't want to include it again:
                            tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_k))
                            tmp = iand(core_space(0:NIfD,i_full), tmp)
                            IC = CountBits(tmp, NIfD)
                            if (IC <= 2) cycle

                            hel = hphf_off_diag_special_case(nI_paired, ilut_paired, core_space(:,ind_k), 2, OpenOrbsI)
                            if (abs(hel) > 0.0_dp) then
                                num_conns(i) = num_conns(i) + 1
                                hamil_pos(num_conns(i)) = ind_k
                                hamil_row(num_conns(i)) = hel
                            end if

                        end do
                    end do
                !end if
            end if

            ! Calculate and add the diagonal element
            num_conns(i) = num_conns(i) + 1
            hel = hphf_diag_helement(nI, core_space(:,i_full)) - Hii
            hamil_pos(num_conns(i)) = i_full
            hamil_row(num_conns(i)) = hel

            ! Now finally allocate and fill in the actual deterministic
            ! Hamiltonian row for this determinant
            allocate(sparse_core_ham(i)%elements(num_conns(i)), stat=ierr)
            allocate(sparse_core_ham(i)%positions(num_conns(i)), stat=ierr)

            sparse_core_ham(i)%num_elements = num_conns(i)
            sparse_core_ham(i)%elements(1:num_conns(i)) = hamil_row(1:num_conns(i))
            sparse_core_ham(i)%positions(1:num_conns(i)) = hamil_pos(1:num_conns(i))

            ! Fill the array of diagonal elements
            core_ham_diag(i) = hel
        end do

        call halt_timer(ham_time)

        write(6,'("Time to create the Hamiltonian:", f9.3)') get_total_time(ham_time); call neci_flush(6)

        ! Optional: sort the Hamiltonian? This could speed up subsequent
        ! multiplications, as we don't jump about in memory so much
        call set_timer(sort_ham_time)
        !do i = 1, determ_sizes(iProcIndex)
        !    call sort(sparse_core_ham(i)%positions, sparse_core_ham(i)%elements)
        !    write(6,*) "i:", i
        !    do j = 1, sparse_core_ham(i)%num_elements
        !        write(6,*) "pos:", sparse_core_ham(i)%positions(j), "elem:", sparse_core_ham(i)%elements(j)
        !    end do
        !end do
        call halt_timer(sort_ham_time)
        !write(6,'("Time to sort the Hamiltonian:", f9.3)') get_total_time(sort_ham_time); call neci_flush(6)

        total_time = get_total_time(aux_time) + get_total_time(sort_aux_time) + &
                      get_total_time(ham_time) + get_total_time(sort_ham_time)

        total_time = get_total_time(aux_time) + get_total_time(sort_aux_time) + get_total_time(ham_time)
        write(6,'("total_time:", f9.3)') total_time; call neci_flush(6)

        ! --- Deallocate all auxiliary arrays -------------

        do i = nbeta_unpaired, 1, -1
            deallocate(beta_dets(i)%pos, stat=ierr)
        end do
        do i = nalpha_unpaired, 1, -1
            deallocate(alpha_dets(i)%pos, stat=ierr)
        end do
        do i = nalpha_unpaired, 1, -1
            deallocate(beta_with_alpha(i)%pos, stat=ierr)
        end do

        deallocate(nbeta_dets)
        deallocate(nalpha_dets)

        deallocate(beta_dets)
        deallocate(alpha_dets)
        deallocate(beta_with_alpha)

        deallocate(num_conns)

        call clear_hash_table(beta_ht)
        deallocate(beta_ht, stat=ierr)
        nullify(beta_ht)

        call clear_hash_table(alpha_ht)
        deallocate(alpha_ht, stat=ierr)
        nullify(alpha_ht)

        deallocate(beta_list)
        deallocate(alpha_list)

        deallocate(intersec_inds)
        deallocate(hamil_row)
        deallocate(hamil_pos)

        call MPIBarrier(ierr, tTimeIn=.false.)

    end subroutine calc_determ_hamil_opt_hphf

    subroutine calc_determ_hamil_opt()

        use DetBitOps, only: MaskAlpha, MaskBeta
        use FciMCData, only: ll_node, s_first_ind, s_last_ind
        use hash, only: init_hash_table, clear_hash_table
        use hash, only: hash_table_lookup, add_hash_table_entry
        use sort_mod
        use SystemData, only: nOccAlpha, nOccBeta

        integer :: i, j, k, ierr
        integer :: IC, IC_beta
        integer :: ind_i, ind_j, ind_k, i_full
        integer :: ind, hash_val
        integer :: ind_beta, ind_alpha, hash_val_beta, hash_val_alpha, ind_alpha_conn
        logical :: tSuccess
        integer(n_int) :: tmp(0:NIfD)

        integer(n_int) :: alpha_ilut(0:NIfTot), alpha_m1_ilut(0:NIfTot), beta_ilut(0:NIfTot), beta_m1_ilut(0:NIfTot)
        integer(n_int) :: beta_ilut_1(0:NIfD), beta_ilut_2(0:NIfD)
        integer(n_int), allocatable :: beta_list(:,:), alpha_list(:,:), alpha_m1_list(:,:), beta_m1_list(:,:)
        integer :: nI_alpha(nOccAlpha), nI_alpha_m1(nOccAlpha-1), nI_beta(nOccBeta), nI_beta_m1(nOccBeta-1)
        integer :: nbeta, nalpha, nbeta_m1, nalpha_m1

        integer :: nintersec
        integer, allocatable :: intersec_inds(:)

        type(ll_node), pointer :: beta_ht(:)
        type(ll_node), pointer :: alpha_ht(:)
        type(ll_node), pointer :: beta_m1_ht(:)
        type(ll_node), pointer :: alpha_m1_ht(:)

        integer, allocatable :: nbeta_dets(:), nalpha_dets(:), nbeta_m1_contribs(:), nalpha_m1_contribs(:)
        type(auxiliary_array), allocatable :: beta_dets(:)
        type(auxiliary_array), allocatable :: alpha_dets(:)
        type(auxiliary_array), allocatable :: beta_m1_contribs(:)
        type(auxiliary_array), allocatable :: alpha_m1_contribs(:)
        type(auxiliary_array), allocatable :: beta_with_alpha(:)

        integer, allocatable :: nalpha_alpha(:), nbeta_beta(:)
        type(auxiliary_array), allocatable :: alpha_alpha(:), beta_beta(:)

        integer :: nI(nel), nJ(nel), nK(nel)

        HElement_t(dp) :: hel
        HElement_t(dp), allocatable :: hamil_row(:)
        integer, allocatable :: hamil_pos(:)
        integer, allocatable :: num_conns(:)

        type(timer), save :: aux_time
        type(timer), save :: sort_aux_time
        type(timer), save :: ham_time
        type(timer), save :: sort_ham_time
        real(dp) :: total_time

        character(len=*), parameter :: t_r = "calc_determ_hamil_opt"

        allocate(beta_list(0:NIfD, determ_space_size), stat=ierr)
        allocate(nbeta_dets(determ_space_size), stat=ierr)
        allocate(alpha_list(0:NIfD, determ_space_size), stat=ierr)
        allocate(nalpha_dets(determ_space_size), stat=ierr)

        allocate(beta_ht( int(determ_space_size/10.0) ), stat=ierr)
        call init_hash_table(beta_ht)
        allocate(alpha_ht( int(determ_space_size/10.0) ), stat=ierr)
        call init_hash_table(alpha_ht)

        call set_timer(aux_time)

        ! --- Set up auxiliary arrays ------------------------------

        ! The number of different alpha and beta strings among the determinants
        nbeta = 0
        nalpha = 0
        ! The number of alpha strings with nel-1 electrons among the dets
        nbeta_m1 = 0
        nalpha_m1 = 0

        ! First time through, just count the number entries in each of the
        ! lists to be set up. Then allocate those lists and fill them in.
        do i = 1, determ_space_size
            ! --- Beta string -----------------------------
            beta_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskBeta)
            call decode_bit_det(nI_beta, beta_ilut)
            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)

            ! If this beta string has not been generated already, add it to
            ! the list of all beta strings, and the associated hash table.
            if (.not. tSuccess) then
                nbeta = nbeta + 1
                nbeta_dets(nbeta) = 0
                beta_list(0:NIfD, nbeta) = beta_ilut(0:NIfD)
                call add_hash_table_entry(beta_ht, nbeta, hash_val_beta)
                ind_beta = nbeta
            end if

            ! Now add this determinant to the list of determinants with this
            ! beta string.
            nbeta_dets(ind_beta) = nbeta_dets(ind_beta) + 1

            ! --- Alpha string ---------------------------
            alpha_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskAlpha)
            call decode_bit_det(nI_alpha, alpha_ilut)
            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)

            ! If this alpha string has not been generated already, add it to
            ! the list of all alpha strings, and the associated hash table.
            if (.not. tSuccess) then
                nalpha = nalpha + 1
                nalpha_dets(nalpha) = 0
                alpha_list(0:NIfD, nalpha) = alpha_ilut(0:NIfD)
                call add_hash_table_entry(alpha_ht, nalpha, hash_val_alpha)
                ind_alpha = nalpha
            end if

            ! Now add this determinant to the list of determinants with this
            ! alpha string.
            nalpha_dets(ind_alpha) = nalpha_dets(ind_alpha) + 1
        end do

        ! Allocate beta-1 and alpha-1 arrays
        allocate(beta_m1_list(0:NIfD, nbeta*(nOccBeta-1)), stat=ierr)
        allocate(nbeta_m1_contribs(nbeta*(nOccBeta-1)), stat=ierr)
        allocate(alpha_m1_list(0:NIfD, nalpha*(nOccAlpha-1)), stat=ierr)
        allocate(nalpha_m1_contribs(nalpha*(nOccAlpha-1)), stat=ierr)

        allocate(beta_m1_ht( int(nbeta*(nOccBeta-1)/10.0) ), stat=ierr)
        call init_hash_table(beta_m1_ht)
        allocate(alpha_m1_ht( int(nalpha*(nOccAlpha-1)/10.0) ), stat=ierr)
        call init_hash_table(alpha_m1_ht)

        ! --- Find the size of beta-1 arrays -----------------

        do i = 1, nbeta
            beta_ilut(0:NIfD) = beta_list(0:NIfD,i)
            call decode_bit_det(nI_beta, beta_ilut)

            do j = 1, nOccBeta
                ! Create the beta-1 orbitals and ilut
                if (j > 1) nI_beta_m1(1:j-1) = nI_beta(1:j-1)
                if (j < nOccBeta) nI_beta_m1(j:nOccBeta-1) = nI_beta(j+1:nOccBeta)
                beta_m1_ilut = beta_ilut
                clr_orb(beta_m1_ilut, nI_beta(j))

                call hash_table_lookup(nI_beta_m1, beta_m1_ilut, NIfD, beta_m1_ht, beta_m1_list, ind, hash_val, tSuccess)

                ! If this beta-1 string has not been generated already, add it to
                ! the list of all beta-1 strings, and the associated hash table.
                if (.not. tSuccess) then
                    nbeta_m1 = nbeta_m1 + 1
                    beta_m1_list(0:NIfD, nbeta_m1) = beta_m1_ilut(0:NIfD)
                    nbeta_m1_contribs(nbeta_m1) = 0
                    call add_hash_table_entry(beta_m1_ht, nbeta_m1, hash_val)
                    ind = nbeta_m1
                end if

                ! Now add this determinant to the list of determinants with this
                ! beta-1 string.
                nbeta_m1_contribs(ind) = nbeta_m1_contribs(ind) + 1
            end do
        end do

        ! --- Find the size of alpha-1 arrays -----------------

        do i = 1, nalpha
            alpha_ilut(0:NIfD) = alpha_list(0:NIfD,i)
            call decode_bit_det(nI_alpha, alpha_ilut)

            do j = 1, nOccAlpha
                ! Create the alpha-1 orbitals and ilut
                if (j > 1) nI_alpha_m1(1:j-1) = nI_alpha(1:j-1)
                if (j < nOccAlpha) nI_alpha_m1(j:nOccAlpha-1) = nI_alpha(j+1:nOccAlpha)
                alpha_m1_ilut = alpha_ilut
                clr_orb(alpha_m1_ilut, nI_alpha(j))

                call hash_table_lookup(nI_alpha_m1, alpha_m1_ilut, NIfD, alpha_m1_ht, alpha_m1_list, ind, hash_val, tSuccess)

                ! If this alpha-1 string has not been generated already, add it to
                ! the list of all alpha-1 strings, and the associated hash table.
                if (.not. tSuccess) then
                    nalpha_m1 = nalpha_m1 + 1
                    alpha_m1_list(0:NIfD, nalpha_m1) = alpha_m1_ilut(0:NIfD)
                    nalpha_m1_contribs(nalpha_m1) = 0
                    call add_hash_table_entry(alpha_m1_ht, nalpha_m1, hash_val)
                    ind = nalpha_m1
                end if

                ! Now add this determinant to the list of determinants with this
                ! alpha-1 string.
                nalpha_m1_contribs(ind) = nalpha_m1_contribs(ind) + 1
            end do
        end do

        ! Now we know the size of the auxiliary arrays

        ! --- Allocate the auxiliary arrays ---------------
        allocate(beta_dets(nbeta), stat=ierr)
        allocate(alpha_dets(nalpha), stat=ierr)
        allocate(beta_m1_contribs(nbeta_m1), stat=ierr)
        allocate(alpha_m1_contribs(nalpha_m1), stat=ierr)
        allocate(beta_with_alpha(nalpha), stat=ierr)

        do i = 1, nbeta
            allocate(beta_dets(i)%pos(nbeta_dets(i)), stat=ierr)
        end do
        do i = 1, nalpha
            allocate(alpha_dets(i)%pos(nalpha_dets(i)), stat=ierr)
        end do
        do i = 1, nbeta_m1
            allocate(beta_m1_contribs(i)%pos(nbeta_m1_contribs(i)), stat=ierr)
        end do
        do i = 1, nalpha_m1
            allocate(alpha_m1_contribs(i)%pos(nalpha_m1_contribs(i)), stat=ierr)
        end do
        do i = 1, nalpha
            allocate(beta_with_alpha(i)%pos(nalpha_dets(i)), stat=ierr)
        end do

        ! --- Now fill the auxiliary arrays ----------------
        nbeta_dets(1:nbeta) = 0
        nalpha_dets(1:nalpha) = 0
        nbeta_m1_contribs(1:nbeta_m1) = 0
        nalpha_m1_contribs(1:nalpha_m1) = 0

        do i = 1, determ_space_size
            ! --- Beta string -----------------------------
            beta_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskBeta)
            call decode_bit_det(nI_beta, beta_ilut)
            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)

            ! --- Alpha string -----------------------------
            alpha_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskAlpha)
            call decode_bit_det(nI_alpha, alpha_ilut)
            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)

            ! Now add this determinant to the list of determinants with this
            ! beta string.
            nbeta_dets(ind_beta) = nbeta_dets(ind_beta) + 1
            beta_dets(ind_beta)%pos(nbeta_dets(ind_beta)) = i

            ! Now add this determinant to the list of determinants with this
            ! alpha string.
            nalpha_dets(ind_alpha) = nalpha_dets(ind_alpha) + 1
            alpha_dets(ind_alpha)%pos(nalpha_dets(ind_alpha)) = i
            beta_with_alpha(ind_alpha)%pos(nalpha_dets(ind_alpha)) = ind_beta
        end do

        do i = 1, nbeta
            beta_ilut(0:NIfD) = beta_list(0:NIfD,i)
            call decode_bit_det(nI_beta, beta_ilut)

            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)
            ! --- Beta-1 string ----------------------------
            do j = 1, nOccBeta
                ! Create the beta-1 orbitals and ilut
                if (j > 1) nI_beta_m1(1:j-1) = nI_beta(1:j-1)
                if (j < nOccBeta) nI_beta_m1(j:nOccBeta-1) = nI_beta(j+1:nOccBeta)
                beta_m1_ilut = beta_ilut
                clr_orb(beta_m1_ilut, nI_beta(j))

                call hash_table_lookup(nI_beta_m1, beta_m1_ilut, NIfD, beta_m1_ht, beta_m1_list, ind, hash_val, tSuccess)

                if (.not. tSuccess) call stop_all(t_r, "beta-1 string not found.")

                ! Now add this determinant to the list of determinants with this
                ! beta-1 string.
                nbeta_m1_contribs(ind) = nbeta_m1_contribs(ind) + 1
                beta_m1_contribs(ind)%pos(nbeta_m1_contribs(ind)) = ind_beta
            end do
        end do

        do i = 1, nalpha
            alpha_ilut(0:NIfD) = alpha_list(0:NIfD,i)
            call decode_bit_det(nI_alpha, alpha_ilut)

            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)
            ! --- Alpha-1 string ---------------------------
            do j = 1, nOccAlpha
                ! Create the alpha-1 orbitals and ilut
                if (j > 1) nI_alpha_m1(1:j-1) = nI_alpha(1:j-1)
                if (j < nOccAlpha) nI_alpha_m1(j:nOccAlpha-1) = nI_alpha(j+1:nOccAlpha)
                alpha_m1_ilut = alpha_ilut
                clr_orb(alpha_m1_ilut, nI_alpha(j))

                call hash_table_lookup(nI_alpha_m1, alpha_m1_ilut, NIfD, alpha_m1_ht, alpha_m1_list, ind, hash_val, tSuccess)

                if (.not. tSuccess) call stop_all(t_r, "alpha-1 string not found.")

                ! Now add this determinant to the list of determinants with this
                ! alpha-1 string.
                nalpha_m1_contribs(ind) = nalpha_m1_contribs(ind) + 1
                alpha_m1_contribs(ind)%pos(nalpha_m1_contribs(ind)) = ind_alpha
            end do
        end do

        deallocate(beta_m1_list)
        deallocate(alpha_m1_list)

        call clear_hash_table(beta_m1_ht)
        deallocate(beta_m1_ht, stat=ierr)
        nullify(beta_m1_ht)

        call clear_hash_table(alpha_m1_ht)
        deallocate(alpha_m1_ht, stat=ierr)
        nullify(alpha_m1_ht)

        allocate(nbeta_beta(nbeta), stat=ierr)
        allocate(nalpha_alpha(nalpha), stat=ierr)
        allocate(beta_beta(nbeta), stat=ierr)
        allocate(alpha_alpha(nalpha), stat=ierr)

        nbeta_beta = 0
        nalpha_alpha = 0

        ! Find the size of the beta_beta array to be created
        do i = 1, nbeta_m1
            do j = 1, nbeta_m1_contribs(i)
                nbeta_beta(beta_m1_contribs(i)%pos(j)) = nbeta_beta(beta_m1_contribs(i)%pos(j)) + nbeta_m1_contribs(i) - 1
            end do
        end do

        ! Allocate the beta_beta array...
        do i = 1, nbeta
            allocate(beta_beta(i)%pos(nbeta_beta(i)), stat=ierr)
        end do

        ! Rezero this so that we can use it as a counter for the following
        nbeta_beta = 0

        ! ...and finally fill the beta_beta array.
        do i = 1, nbeta_m1
            do j = 1, nbeta_m1_contribs(i)
                do k = j+1, nbeta_m1_contribs(i)
                    nbeta_beta( beta_m1_contribs(i)%pos(j) ) = nbeta_beta( beta_m1_contribs(i)%pos(j) ) + 1
                    beta_beta( beta_m1_contribs(i)%pos(j) )%pos(nbeta_beta( beta_m1_contribs(i)%pos(j) )) = beta_m1_contribs(i)%pos(k)

                    nbeta_beta( beta_m1_contribs(i)%pos(k) ) = nbeta_beta( beta_m1_contribs(i)%pos(k) ) + 1
                    beta_beta( beta_m1_contribs(i)%pos(k) )%pos(nbeta_beta( beta_m1_contribs(i)%pos(k) )) = beta_m1_contribs(i)%pos(j)
                end do
            end do
        end do

        ! Find the size of the alpha_alpha array to be created
        do i = 1, nalpha_m1
            do j = 1, nalpha_m1_contribs(i)
                nalpha_alpha(alpha_m1_contribs(i)%pos(j)) = nalpha_alpha(alpha_m1_contribs(i)%pos(j)) + nalpha_m1_contribs(i) - 1
            end do
        end do

        ! Allocate the alpha_alpha array...
        do i = 1, nalpha
            allocate(alpha_alpha(i)%pos(nalpha_alpha(i)), stat=ierr)
        end do

        ! Rezero this so that we can use it as a counter for the following
        nalpha_alpha = 0

        ! ...and finally fill the alpha_alpha array.
        do i = 1, nalpha_m1
            do j = 1, nalpha_m1_contribs(i)
                do k = j+1, nalpha_m1_contribs(i)
                    nalpha_alpha( alpha_m1_contribs(i)%pos(j) ) = nalpha_alpha( alpha_m1_contribs(i)%pos(j) ) + 1
                    alpha_alpha( alpha_m1_contribs(i)%pos(j) )%pos(nalpha_alpha( alpha_m1_contribs(i)%pos(j) )) = alpha_m1_contribs(i)%pos(k)

                    nalpha_alpha( alpha_m1_contribs(i)%pos(k) ) = nalpha_alpha( alpha_m1_contribs(i)%pos(k) ) + 1
                    alpha_alpha( alpha_m1_contribs(i)%pos(k) )%pos(nalpha_alpha( alpha_m1_contribs(i)%pos(k) )) = alpha_m1_contribs(i)%pos(j)
                end do
            end do
        end do

        call halt_timer(aux_time)

        write(6,'("Time to create auxiliary arrays:", f9.3)') get_total_time(aux_time); call neci_flush(6)

        ! Sort auxiliary arrays into the required order
        call set_timer(sort_aux_time)

        do i = 1, nbeta
            call sort(beta_beta(i)%pos)
        end do

        do i = 1, nalpha
            call sort(alpha_alpha(i)%pos)
        end do

        do i = 1, nalpha
            call sort(beta_with_alpha(i)%pos, alpha_dets(i)%pos)
        end do

        call halt_timer(sort_aux_time)

        ! Deallocate the arrays which we don't need for the following.

        do i = nbeta_m1, 1, -1
            deallocate(beta_m1_contribs(i)%pos, stat=ierr)
        end do
        do i = nalpha_m1, 1, -1
            deallocate(alpha_m1_contribs(i)%pos, stat=ierr)
        end do

        deallocate(nbeta_m1_contribs)
        deallocate(nalpha_m1_contribs)
        deallocate(beta_m1_contribs)
        deallocate(alpha_m1_contribs)

        write(6,'("Time to sort auxiliary arrays:", f9.3)') get_total_time(sort_aux_time); call neci_flush(6)

        ! Actually create the Hamiltonian
        call set_timer(ham_time)

        allocate(intersec_inds(nbeta), stat=ierr)

        allocate(hamil_row(determ_space_size), stat=ierr)
        allocate(hamil_pos(determ_space_size), stat=ierr)

        allocate(sparse_core_ham(determ_sizes(iProcIndex)), stat=ierr)
        allocate(core_ham_diag(determ_sizes(iProcIndex)), stat=ierr)

        allocate(num_conns(determ_space_size), stat=ierr)
        num_conns = 0

        ! Loop over the determinants on this process
        do i = 1, determ_sizes(iProcIndex)
            ! Find the index in the *full* list of determinants
            i_full = i + determ_displs(iProcIndex)

            call decode_bit_det(nI, core_space(:,i_full))

            ! --- Beta string -----------------------------
            beta_ilut(0:NIfD) = iand(core_space(0:NIfD,i_full), MaskBeta)
            call decode_bit_det(nI_beta, beta_ilut)
            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind_beta, hash_val_beta, tSuccess)

            do j = 1, nbeta_dets(ind_beta)
                ind_j = beta_dets(ind_beta)%pos(j)
                if (i_full == ind_j) cycle
                tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_j))
                tmp = iand(core_space(0:NIfD,i_full), tmp)
                IC = CountBits(tmp, NIfD)
                if (IC <= 2) then
                    call decode_bit_det(nJ, core_space(:,ind_j))
                    hel = get_helement(nI, nJ, IC, core_space(:,i_full), core_space(:,ind_j))
                    if (abs(hel) > 0.0_dp) then
                        num_conns(i) = num_conns(i) + 1
                        hamil_pos(num_conns(i)) = ind_j
                        hamil_row(num_conns(i)) = hel
                    end if
                end if
            end do

            ! --- Alpha string -----------------------------
            alpha_ilut(0:NIfD) = iand(core_space(0:NIfD,i_full), MaskAlpha)
            call decode_bit_det(nI_alpha, alpha_ilut)
            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind_alpha, hash_val_alpha, tSuccess)

            do j = 1, nalpha_dets(ind_alpha)
                ind_j = alpha_dets(ind_alpha)%pos(j)
                if (i_full == ind_j) cycle
                tmp = ieor(core_space(0:NIfD,i_full), core_space(0:NIfD,ind_j))
                tmp = iand(core_space(0:NIfD,i_full), tmp)
                IC = CountBits(tmp, NIfD)
                if (IC <= 2) then
                    call decode_bit_det(nJ, core_space(:,ind_j))
                    hel = get_helement(nI, nJ, IC, core_space(:,i_full), core_space(:,ind_j))
                    if (abs(hel) > 0.0_dp) then
                        num_conns(i) = num_conns(i) + 1
                        hamil_pos(num_conns(i)) = ind_j
                        hamil_row(num_conns(i)) = hel
                    end if
                end if
            end do

            ! Loop through alpha strings connected to this ind_alpha by a single excitation
            do j = 1, nalpha_alpha(ind_alpha)
                ! This is the index of the connected alpha string
                ind_alpha_conn = alpha_alpha(ind_alpha)%pos(j)

                call find_intersec( nalpha_dets(ind_alpha_conn), nbeta_beta(ind_beta), &
                               beta_with_alpha(ind_alpha_conn)%pos, beta_beta(ind_beta)%pos, &
                               intersec_inds, nintersec )

                do k = 1, nintersec
                    ind_k = alpha_dets(ind_alpha_conn)%pos( intersec_inds(k) )
                    call decode_bit_det(nK, core_space(:,ind_k))
                    hel = get_helement(nI, nK, 2, core_space(:,i_full), core_space(:,ind_k))
                    if (abs(hel) > 0.0_dp) then
                        num_conns(i) = num_conns(i) + 1
                        hamil_pos(num_conns(i)) = ind_k
                        hamil_row(num_conns(i)) = hel
                    end if
                end do
            end do

            ! Calculate and add the diagonal element
            hel = get_helement(nI, nI, 0) - Hii
            num_conns(i) = num_conns(i) + 1
            hamil_pos(num_conns(i)) = i_full
            hamil_row(num_conns(i)) = hel

            ! Now finally allocate and fill in the actual deterministic
            ! Hamiltonian row for this determinant
            allocate(sparse_core_ham(i)%elements(num_conns(i)), stat=ierr)
            allocate(sparse_core_ham(i)%positions(num_conns(i)), stat=ierr)
            sparse_core_ham(i)%num_elements = num_conns(i)
            sparse_core_ham(i)%elements(1:num_conns(i)) = hamil_row(1:num_conns(i))
            sparse_core_ham(i)%positions(1:num_conns(i)) = hamil_pos(1:num_conns(i))

            ! Fill the array of diagonal elements
            core_ham_diag(i) = hel
        end do

        call halt_timer(ham_time)

        write(6,'("Time to create the Hamiltonian:", f9.3)') get_total_time(ham_time); call neci_flush(6)

        ! Optional: sort the Hamiltonian? This could speed up subsequent
        ! multiplications, as we don't jump about in memory so much
        !call set_timer(sort_ham_time)
        !do i = 1, determ_sizes(iProcIndex)
        !    call sort(sparse_core_ham(i)%positions, sparse_core_ham(i)%elements)
        !end do
        !call halt_timer(sort_ham_time)
        !write(6,'("Time to sort the Hamiltonian:", f9.3)') get_total_time(sort_ham_time); call neci_flush(6)

        !total_time = get_total_time(aux_time) + get_total_time(sort_aux_time) + &
        !              get_total_time(ham_time) + get_total_time(sort_ham_time)

        total_time = get_total_time(aux_time) + get_total_time(sort_aux_time) + get_total_time(ham_time)
        write(6,'("total_time:", f9.3)') total_time; call neci_flush(6)

        ! --- Deallocate all auxiliary arrays -------------

        do i = nbeta, 1, -1
            deallocate(beta_dets(i)%pos, stat=ierr)
        end do
        do i = nalpha, 1, -1
            deallocate(alpha_dets(i)%pos, stat=ierr)
        end do
        do i = nalpha, 1, -1
            deallocate(beta_with_alpha(i)%pos, stat=ierr)
        end do

        deallocate(nbeta_dets)
        deallocate(nalpha_dets)

        deallocate(beta_dets)
        deallocate(alpha_dets)
        deallocate(beta_with_alpha)

        deallocate(num_conns)

        call clear_hash_table(beta_ht)
        deallocate(beta_ht, stat=ierr)
        nullify(beta_ht)

        call clear_hash_table(alpha_ht)
        deallocate(alpha_ht, stat=ierr)
        nullify(alpha_ht)

        deallocate(beta_list)
        deallocate(alpha_list)

        deallocate(intersec_inds)
        deallocate(hamil_row)
        deallocate(hamil_pos)

        call MPIBarrier(ierr, tTimeIn=.false.)

    end subroutine calc_determ_hamil_opt

    pure subroutine find_intersec( nelem_in_1, nelem_in_2, arr_1, arr_2, intersec, nelem_out )

        integer, intent(in) :: nelem_in_1, nelem_in_2
        integer, intent(in) :: arr_1(1:), arr_2(1:)
        integer, intent(inout) :: intersec(1:)
        integer, intent(out) :: nelem_out

        integer :: i, j

        nelem_out = 0
        i = 1
        j = 1

        do while (i <= nelem_in_1 .and. j <= nelem_in_2)
            if (arr_1(i) < arr_2(j)) then
                i = i + 1
            else if (arr_1(i) > arr_2(j)) then
                j = j + 1
            else
                nelem_out = nelem_out + 1
                !intersec(nelem_out) = arr_1(i)
                intersec(nelem_out) = i
                i = i + 1
                j = j + 1
            end if
        end do

    end subroutine find_intersec

    subroutine calc_determ_hamil_opt_old()

        use DetBitOps, only: MaskAlpha, MaskBeta
        use FciMCData, only: ll_node, s_first_ind, s_last_ind
        use hash, only: init_hash_table, clear_hash_table
        use hash, only: hash_table_lookup, add_hash_table_entry
        use sort_mod
        use SystemData, only: nOccAlpha, nOccBeta

        integer :: i, j, k, ierr
        integer(n_int) :: alpha_ilut(0:NIfTot), alpha_m1_ilut(0:NIfTot), beta_ilut(0:NIfTot)
        integer(n_int) :: beta_ilut_1(0:NIfD), beta_ilut_2(0:NIfD)
        integer(n_int), allocatable :: beta_list(:,:), alpha_list(:,:), alpha_m1_list(:,:)
        integer :: nI_alpha(nOccAlpha), nI_alpha_m1(nOccAlpha-1), nI_beta(nOccBeta)
        integer :: ind, hash_val, nbeta, nalpha, nalpha_m1
        logical :: tSuccess

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC, IC_beta, ind_i, ind_j, i_new, j_new, true_ind

        type(ll_node), pointer :: beta_ht(:)
        type(ll_node), pointer :: alpha_ht(:)
        type(ll_node), pointer :: alpha_m1_ht(:)

        integer, allocatable :: nbeta_dets(:), nalpha_dets(:), nalpha_m1_dets(:)
        type(auxiliary_array), allocatable :: beta_dets(:)
        type(auxiliary_array), allocatable :: alpha_dets(:)
        type(auxiliary_array), allocatable :: alpha_m1_dets(:)

        integer, allocatable :: num_conns(:)
        integer, allocatable :: nI_list(:,:)
        HElement_t(dp) :: hel

        character(len=*), parameter :: t_r = "calc_determ_hamil_opt"

        type(timer), save :: a_time
        type(timer), save :: b_time
        type(timer), save :: c_time
        type(timer), save :: d_time
        type(timer), save :: e_time
        type(timer), save :: f_time
        type(timer), save :: g_time
        type(timer), save :: sort_time
        real(dp) :: total_time

        call set_timer(a_time)

        allocate(num_conns(determ_space_size), stat=ierr)

        allocate(beta_list(0:NIfD, determ_space_size), stat=ierr)
        allocate(nbeta_dets(determ_space_size), stat=ierr)
        allocate(alpha_list(0:NIfD, determ_space_size), stat=ierr)
        allocate(nalpha_dets(determ_space_size), stat=ierr)
        allocate(alpha_m1_list(0:NIfD, determ_space_size*(nOccAlpha-1)), stat=ierr)
        allocate(nalpha_m1_dets(determ_space_size*(nOccAlpha-1)), stat=ierr)

        allocate(beta_ht(determ_space_size), stat=ierr)
        call init_hash_table(beta_ht)
        allocate(alpha_ht(determ_space_size), stat=ierr)
        call init_hash_table(alpha_ht)

        allocate(alpha_m1_ht(determ_space_size*(nOccAlpha-1)), stat=ierr)
        call init_hash_table(alpha_m1_ht)

        ! --- Set up auxiliary arrays ------------------------------

        nbeta = 0
        nalpha = 0
        nalpha_m1 = 0

        ! First time through, just count the number entries in each of the
        ! lists to be set up. Then allocate those lists and fill them in.
        do i = 1, determ_space_size
            ! --- Beta string -----------------------------
            beta_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskBeta)
            call decode_bit_det(nI_beta, beta_ilut)
            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind, hash_val, tSuccess)

            ! If this beta string has not been generated already, add it to
            ! the list of all beta strings, and the associated hash table.
            if (.not. tSuccess) then
                nbeta = nbeta + 1
                nbeta_dets(nbeta) = 0
                beta_list(0:NIfD, nbeta) = beta_ilut(0:NIfD)
                call add_hash_table_entry(beta_ht, nbeta, hash_val)
                ind = nbeta
            end if

            ! Now add this determinant to the list of determinants with this
            ! beta string.
            nbeta_dets(ind) = nbeta_dets(ind) + 1

            ! --- Alpha string ---------------------------
            alpha_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskAlpha)
            call decode_bit_det(nI_alpha, alpha_ilut)
            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind, hash_val, tSuccess)

            ! If this alpha string has not been generated already, add it to
            ! the list of all alpha strings, and the associated hash table.
            if (.not. tSuccess) then
                nalpha = nalpha + 1
                nalpha_dets(nalpha) = 0
                alpha_list(0:NIfD, nalpha) = alpha_ilut(0:NIfD)
                call add_hash_table_entry(alpha_ht, nalpha, hash_val)
                ind = nalpha
            end if

            ! Now add this determinant to the list of determinants with this
            ! alpha string.
            nalpha_dets(ind) = nalpha_dets(ind) + 1

            ! --- Set up alpha-1 arrays -----------------

            do j = 1, nOccAlpha
                ! Create the alpha-1 orbitals and ilut
                if (j > 1) nI_alpha_m1(1:j-1) = nI_alpha(1:j-1)
                if (j < nOccAlpha) nI_alpha_m1(j:nOccAlpha-1) = nI_alpha(j+1:nOccAlpha)
                alpha_m1_ilut = alpha_ilut
                clr_orb(alpha_m1_ilut, nI_alpha(j))

                call hash_table_lookup(nI_alpha_m1, alpha_m1_ilut, NIfD, alpha_m1_ht, alpha_m1_list, ind, hash_val, tSuccess)

                ! If this alpha-1 string has not been generated already, add it to
                ! the list of all alpha-1 strings, and the associated hash table.
                if (.not. tSuccess) then
                    nalpha_m1 = nalpha_m1 + 1
                    alpha_m1_list(0:NIfD, nalpha_m1) = alpha_m1_ilut(0:NIfD)
                    nalpha_m1_dets(nalpha_m1) = 0
                    call add_hash_table_entry(alpha_m1_ht, nalpha_m1, hash_val)
                    ind = nalpha_m1
                end if

                ! Now add this determinant to the list of determinants with this
                ! alpha-1 string.
                nalpha_m1_dets(ind) = nalpha_m1_dets(ind) + 1
            end do
        end do

        ! Now we know the size of the auxiliary arrays

        ! Allocate the auxiliary arrays
        allocate(beta_dets(nbeta), stat=ierr)
        allocate(alpha_dets(nalpha), stat=ierr)
        allocate(alpha_m1_dets(nalpha_m1), stat=ierr)

        do i = 1, nbeta
            allocate(beta_dets(i)%pos(nbeta_dets(i)), stat=ierr)
        end do
        do i = 1, nalpha
            allocate(alpha_dets(i)%pos(nalpha_dets(i)), stat=ierr)
        end do
        do i = 1, nalpha_m1
            allocate(alpha_m1_dets(i)%pos(nalpha_m1_dets(i)), stat=ierr)
        end do

        ! --- Now fille the auxiliary arrays --------------
        nbeta_dets(1:nbeta) = 0
        nalpha_dets(1:nalpha) = 0
        nalpha_m1_dets(1:nalpha_m1) = 0

        do i = 1, determ_space_size
            ! --- Beta string -----------------------------
            beta_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskBeta)
            call decode_bit_det(nI_beta, beta_ilut)
            call hash_table_lookup(nI_beta, beta_ilut, NIfD, beta_ht, beta_list, ind, hash_val, tSuccess)

            ! Now add this determinant to the list of determinants with this
            ! beta string.
            nbeta_dets(ind) = nbeta_dets(ind) + 1
            beta_dets(ind)%pos(nbeta_dets(ind)) = i

            ! --- Alpha string -----------------------------
            alpha_ilut(0:NIfD) = iand(core_space(0:NIfD,i), MaskAlpha)
            call decode_bit_det(nI_alpha, alpha_ilut)
            call hash_table_lookup(nI_alpha, alpha_ilut, NIfD, alpha_ht, alpha_list, ind, hash_val, tSuccess)

            ! Now add this determinant to the list of determinants with this
            ! alpha string.
            nalpha_dets(ind) = nalpha_dets(ind) + 1
            alpha_dets(ind)%pos(nalpha_dets(ind)) = i

            ! --- Alpha-1 string ---------------------------
            do j = 1, nOccAlpha
                ! Create the alpha-1 orbitals and ilut
                if (j > 1) nI_alpha_m1(1:j-1) = nI_alpha(1:j-1)
                if (j < nOccAlpha) nI_alpha_m1(j:nOccAlpha-1) = nI_alpha(j+1:nOccAlpha)
                alpha_m1_ilut = alpha_ilut
                clr_orb(alpha_m1_ilut, nI_alpha(j))

                call hash_table_lookup(nI_alpha_m1, alpha_m1_ilut, NIfD, alpha_m1_ht, alpha_m1_list, ind, hash_val, tSuccess)

                if (.not. tSuccess) call stop_all(t_r, "alpha-1 string not found.")

                ! Now add this determinant to the list of determinants with this
                ! alpha-1 string.
                nalpha_m1_dets(ind) = nalpha_m1_dets(ind) + 1
                alpha_m1_dets(ind)%pos(nalpha_m1_dets(ind)) = i
            end do
        end do

        call clear_hash_table(beta_ht)
        deallocate(beta_ht, stat=ierr)
        nullify(beta_ht)

        call clear_hash_table(alpha_ht)
        deallocate(alpha_ht, stat=ierr)
        nullify(alpha_ht)

        call clear_hash_table(alpha_m1_ht)
        deallocate(alpha_m1_ht, stat=ierr)
        nullify(alpha_m1_ht)

        deallocate(beta_list)
        deallocate(alpha_list)
        deallocate(alpha_m1_list)

        allocate(nI_list(nel, determ_space_size), stat=ierr)
        do i = 1, determ_space_size
            call decode_bit_det(nI_list(:,i), core_space(:,i))
        end do

        call halt_timer(a_time)

        call set_timer(sort_time)

        do k = 1, nbeta
            call sort(beta_dets(k)%pos)
        end do
        do k = 1, nalpha
            call sort(alpha_dets(k)%pos)
        end do
        do k = 1, nalpha_m1
            call sort(alpha_m1_dets(k)%pos)
        end do

        call halt_timer(sort_time)

        call set_timer(b_time)

        num_conns = 0

        do k = 1, nbeta
            do j = 1, nbeta_dets(k)

                ind_j = beta_dets(k)%pos(j)
                if (ind_j < s_first_ind .or. ind_j > s_last_ind) cycle
                j_new = ind_j - determ_displs(iProcIndex)

                !do i = j+1, nbeta_dets(k)
                do i = 1, nbeta_dets(k)
                    if (i == j) cycle

                    ind_i = beta_dets(k)%pos(i)

                    tmp = ieor(core_space(0:NIfD,ind_j), core_space(0:NIfD,ind_i))
                    tmp = iand(core_space(0:NIfD,ind_j), tmp)
                    IC = CountBits(tmp, NIfD)

                    if (IC <= 2) then
                        num_conns(j_new) = num_conns(j_new) + 1
                        !if (ind_i >= s_first_ind .and. ind_i <= s_last_ind) then
                        !    i_new = ind_i - determ_displs(iProcIndex)
                        !    num_conns(i_new) = num_conns(i_new) + 1
                        !end if
                    end if

                end do
            end do
        end do

        do k = 1, nalpha
            do j = 1, nalpha_dets(k)

                ind_j = alpha_dets(k)%pos(j)
                if (ind_j < s_first_ind .or. ind_j > s_last_ind) cycle
                j_new = ind_j - determ_displs(iProcIndex)

                !do i = j+1, nalpha_dets(k)
                do i = 1, nalpha_dets(k)
                    if (i == j) cycle

                    ind_i = alpha_dets(k)%pos(i)

                    tmp = ieor(core_space(0:NIfD,ind_j), core_space(0:NIfD,ind_i))
                    tmp = iand(core_space(0:NIfD,ind_j), tmp)
                    IC = CountBits(tmp, NIfD)

                    if (IC <= 2) then
                        num_conns(j_new) = num_conns(j_new) + 1
                        !if (ind_i >= s_first_ind .and. ind_i <= s_last_ind) then
                        !    i_new = ind_i - determ_displs(iProcIndex)
                        !    num_conns(i_new) = num_conns(i_new) + 1
                        !end if
                    end if

                end do
            end do
        end do

        call halt_timer(b_time)

        call set_timer(c_time)

        do k = 1, nalpha_m1
            do j = 1, nalpha_m1_dets(k)

                ind_j = alpha_m1_dets(k)%pos(j)
                if (ind_j < s_first_ind .or. ind_j > s_last_ind) cycle
                j_new = ind_j - determ_displs(iProcIndex)

                beta_ilut_1(0:NIfD) = iand(core_space(0:NIfD,ind_j), MaskBeta)

                !do i = j+1, nalpha_m1_dets(k)
                do i = 1, nalpha_m1_dets(k)
                    if (i == j) cycle

                    ind_i = alpha_m1_dets(k)%pos(i)

                    beta_ilut_2(0:NIfD) = iand(core_space(0:NIfD,ind_i), MaskBeta)

                    tmp = ieor(beta_ilut_1(0:NIfD), beta_ilut_2(0:NIfD))
                    tmp = iand(beta_ilut_1(0:NIfD), tmp)
                    IC_beta = CountBits(tmp, NIfD)

                    if (IC_beta == 1) then

                        tmp = ieor(core_space(0:NIfD,ind_j), core_space(0:NIfD,ind_i))
                        tmp = iand(core_space(0:NIfD,ind_j), tmp)
                        IC = CountBits(tmp, NIfD)

                        if (IC == 2) then

                            num_conns(j_new) = num_conns(j_new) + 1
                            !if (ind_i >= s_first_ind .and. ind_i <= s_last_ind) then
                            !    i_new = ind_i - determ_displs(iProcIndex)
                            !    num_conns(i_new) = num_conns(i_new) + 1
                            !end if
                        end if
                    end if

                end do
            end do
        end do

        call halt_timer(c_time)

        call set_timer(d_time)

        allocate(sparse_core_ham(determ_sizes(iProcIndex)), stat=ierr)
        allocate(core_ham_diag(determ_sizes(iProcIndex)), stat=ierr)

        do i = 1, determ_sizes(iProcIndex)
            ! Add 1 for diagonal elements
            num_conns(i) = num_conns(i) + 1
            !write(6,*) "i: ", i, "num_elements: ", sparse_core_ham(i)%num_elements

            allocate(sparse_core_ham(i)%elements(num_conns(i)), stat=ierr)
            allocate(sparse_core_ham(i)%positions(num_conns(i)), stat=ierr)
            sparse_core_ham(i)%num_elements = num_conns(i)
        end do

        num_conns = 0

        call halt_timer(d_time)

        call set_timer(e_time)

        do k = 1, nbeta
            do j = 1, nbeta_dets(k)

                ind_j = beta_dets(k)%pos(j)
                if (ind_j < s_first_ind .or. ind_j > s_last_ind) cycle
                j_new = ind_j - determ_displs(iProcIndex)

                !do i = j+1, nbeta_dets(k)
                do i = 1, nbeta_dets(k)
                    if (i == j) cycle

                    ind_i = beta_dets(k)%pos(i)

                    tmp = ieor(core_space(0:NIfD,ind_j), core_space(0:NIfD,ind_i))
                    tmp = iand(core_space(0:NIfD,ind_j), tmp)
                    IC = CountBits(tmp, NIfD)

                    if (IC <= 2) then
                        hel = get_helement(nI_list(:,ind_i), nI_list(:,ind_j), IC, core_space(:,ind_i), core_space(:,ind_j))
                        num_conns(j_new) = num_conns(j_new) + 1
                        sparse_core_ham(j_new)%positions(num_conns(j_new)) = ind_i
                        sparse_core_ham(j_new)%elements(num_conns(j_new)) = hel

                        !if (ind_i >= s_first_ind .and. ind_i <= s_last_ind) then
                        !    i_new = ind_i - determ_displs(iProcIndex)
                        !    num_conns(i_new) = num_conns(i_new) + 1
                        !    sparse_core_ham(i_new)%positions(num_conns(i_new)) = ind_j
                        !    sparse_core_ham(i_new)%elements(num_conns(i_new)) = hel
                        !end if
                    end if

                end do
            end do
        end do

        do k = 1, nalpha
            do j = 1, nalpha_dets(k)

                ind_j = alpha_dets(k)%pos(j)
                if (ind_j < s_first_ind .or. ind_j > s_last_ind) cycle
                j_new = ind_j - determ_displs(iProcIndex)

                !do i = j+1, nalpha_dets(k)
                do i = 1, nalpha_dets(k)
                    if (i == j) cycle

                    ind_i = alpha_dets(k)%pos(i)

                    tmp = ieor(core_space(0:NIfD,ind_j), core_space(0:NIfD,ind_i))
                    tmp = iand(core_space(0:NIfD,ind_j), tmp)
                    IC = CountBits(tmp, NIfD)

                    if (IC <= 2) then
                        hel = get_helement(nI_list(:,ind_i), nI_list(:,ind_j), IC, core_space(:,ind_i), core_space(:,ind_j))
                        num_conns(j_new) = num_conns(j_new) + 1
                        sparse_core_ham(j_new)%positions(num_conns(j_new)) = ind_i
                        sparse_core_ham(j_new)%elements(num_conns(j_new)) = hel

                        !if (ind_i >= s_first_ind .and. ind_i <= s_last_ind) then
                        !    i_new = ind_i - determ_displs(iProcIndex)
                        !    num_conns(i_new) = num_conns(i_new) + 1
                        !    sparse_core_ham(i_new)%positions(num_conns(i_new)) = ind_j
                        !    sparse_core_ham(i_new)%elements(num_conns(i_new)) = hel
                        !end if
                    end if

                end do
            end do
        end do

        call halt_timer(e_time)

        call set_timer(f_time)

        do k = 1, nalpha_m1
            do j = 1, nalpha_m1_dets(k)

                ind_j = alpha_m1_dets(k)%pos(j)
                if (ind_j < s_first_ind .or. ind_j > s_last_ind) cycle
                j_new = ind_j - determ_displs(iProcIndex)

                beta_ilut_1(0:NIfD) = iand(core_space(0:NIfD,ind_j), MaskBeta)

                !do i = j+1, nalpha_m1_dets(k)
                do i = 1, nalpha_m1_dets(k)
                    if (i == j) cycle

                    ind_i = alpha_m1_dets(k)%pos(i)

                    beta_ilut_2(0:NIfD) = iand(core_space(0:NIfD,ind_i), MaskBeta)

                    tmp = ieor(beta_ilut_1(0:NIfD), beta_ilut_2(0:NIfD))
                    tmp = iand(beta_ilut_1(0:NIfD), tmp)
                    IC_beta = CountBits(tmp, NIfD)

                    if (IC_beta == 1) then

                        tmp = ieor(core_space(0:NIfD,ind_j), core_space(0:NIfD,ind_i))
                        tmp = iand(core_space(0:NIfD,ind_j), tmp)
                        IC = CountBits(tmp, NIfD)

                        if (IC == 2) then

                            hel = get_helement(nI_list(:,ind_i), nI_list(:,ind_j), IC, core_space(:,ind_i), core_space(:,ind_j))
                            num_conns(j_new) = num_conns(j_new) + 1
                            sparse_core_ham(j_new)%positions(num_conns(j_new)) = ind_i
                            sparse_core_ham(j_new)%elements(num_conns(j_new)) = hel

                            !if (ind_i >= s_first_ind .and. ind_i <= s_last_ind) then
                            !    i_new = ind_i - determ_displs(iProcIndex)
                            !    num_conns(i_new) = num_conns(i_new) + 1
                            !    sparse_core_ham(i_new)%positions(num_conns(i_new)) = ind_j
                            !    sparse_core_ham(i_new)%elements(num_conns(i_new)) = hel
                            !end if
                        end if
                    end if

                end do
            end do
        end do

        call halt_timer(f_time)

        call set_timer(g_time)

        ! Add diagonal elements in
        do i = 1, determ_sizes(iProcIndex)
            ! Add 1 for diagonal elements
            true_ind = i + determ_displs(iProcIndex)
            num_conns(i) = num_conns(i) + 1
            sparse_core_ham(i)%positions(num_conns(i)) = true_ind
            sparse_core_ham(i)%elements(num_conns(i)) = get_helement(nI_list(:,true_ind), nI_list(:,true_ind), 0) - Hii
            core_ham_diag(i) = sparse_core_ham(i)%elements(num_conns(i))
        end do

        call halt_timer(g_time)

        write(6,'("a_time:", f9.3)') get_total_time(a_time)
        write(6,'("b_time:", f9.3)') get_total_time(b_time)
        write(6,'("c_time:", f9.3)') get_total_time(c_time)
        write(6,'("d_time:", f9.3)') get_total_time(d_time)
        write(6,'("e_time:", f9.3)') get_total_time(e_time)
        write(6,'("f_time:", f9.3)') get_total_time(f_time)
        write(6,'("g_time:", f9.3)') get_total_time(g_time)
        write(6,'("sort_time:", f9.3)') get_total_time(sort_time)

        total_time = get_total_time(a_time) + get_total_time(b_time) + get_total_time(c_time) + &
          get_total_time(d_time) + get_total_time(e_time) + get_total_time(f_time) + get_total_time(g_time)
        write(6,'("total_time:", f9.3)') total_time

        do i = nbeta, 1, -1
            deallocate(beta_dets(i)%pos, stat=ierr)
        end do
        do i = nalpha_m1, 1, -1
            deallocate(alpha_m1_dets(i)%pos, stat=ierr)
        end do

        deallocate(nbeta_dets)
        deallocate(nalpha_m1_dets)
        deallocate(beta_dets)
        deallocate(alpha_m1_dets)

        call MPIBarrier(ierr, tTimeIn=.false.)

    end subroutine calc_determ_hamil_opt_old

end module fast_determ_hamil
