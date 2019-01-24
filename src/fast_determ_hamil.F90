#include "macros.h"

module fast_determ_hamil

    use bit_rep_data, only: NIfTot, NIfDBO, NIfD
    use bit_reps, only: decode_bit_det
    use constants
    use DetBitOps, only: CountBits, TestClosedShellDet
    use Determinants, only: get_helement
    use FciMCData, only: determ_space_size, determ_sizes, determ_displs, &
                         Hii, core_ham_diag, core_space
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_off_diag_helement_opt
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier
    use sparse_arrays
    use SystemData, only: tHPHF, nel
    use timing_neci
    use util_mod, only: get_free_unit

    implicit none

    type auxiliary_array
        integer, allocatable :: pos(:)
    end type auxiliary_array

contains

    subroutine calc_determ_hamil_opt()

        use DetBitOps, only: MaskAlpha, MaskBeta
        use FciMCData, only: ll_node, s_first_ind, s_last_ind
        use hash, only: init_hash_table, clear_hash_table
        use hash, only: hash_table_lookup, add_hash_table_entry
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

        a_time%timer_name = 'a_time'

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

    end subroutine calc_determ_hamil_opt

end module fast_determ_hamil
