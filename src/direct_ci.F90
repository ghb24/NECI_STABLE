#include "macros.h"

module direct_ci

    use DetBitOps, only: get_single_parity
    use enumerate_excitations, only: simple_excit_store
    use Integrals_neci, only: get_umat_el
    use IntegralsData, only: ptr_getumatel
    use OneEInts, only: GetTMatEl
    use ras

    implicit none

contains

    subroutine perform_multiplication(ras, classes, vec_in, vec_out)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(inout), allocatable, dimension(:) :: classes
        type(ras_vector), intent(in) :: vec_in(size(classes), size(classes), 0:7)
        type(ras_vector), intent(inout) :: vec_out(size(classes), size(classes), 0:7)
        type(ras_vector) :: c(size(classes), size(classes), 0:7)
        integer :: num_classes, class_i, class_j, class_k, class_m, par_1, par_2
        integer :: i, j, k, l, m, ind_i, ind_j, ind_k, temp_1, temp_2, nras1, nras3
        integer :: sym_i, sym_j, sym_k, sym_m, sym_prod1, sym_prod2
        integer :: ex1(2), ex2(2)
        real(dp), allocatable, dimension(:) :: factor
        integer :: string_i(tot_nelec), string_j(tot_nelec), string_k(tot_nelec)
        integer(n_int) :: ilut_i(0:NIfTot), ilut_j(0:NIfTot), ilut_k(0:NIfTot)
        type(ras_factors) :: factors(size(classes), 0:7), r(size(classes), 0:7), &
                             v(size(classes), 0:7)
        type(simple_excit_store), target :: gen_store
        logical :: none_left

        num_classes = size(classes)

        do class_i = 1, num_classes
            do sym_i = 0, 7
                allocate(factors(class_i,sym_i)%elements(1:classes(class_i)%num_sym(sym_i)))
                allocate(r(class_i,sym_i)%elements(1:classes(class_i)%num_sym(sym_i)))
                allocate(v(class_i,sym_i)%elements(1:classes(class_i)%num_sym(sym_i)))
            end do
        end do

        do class_i = 1, num_classes
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)
                do sym_i = 0, 7
                    sym_j = ieor(HFSym, sym_i)

                    if (classes(class_i)%num_sym(sym_i) == 0 .or. classes(class_j)%num_sym(sym_j) == 0) cycle
                    allocate(c(class_i,class_j,sym_i)%elements(1:classes(class_i)%num_sym(sym_i), &
                                                               1:classes(class_j)%num_sym(sym_j)))
                    vec_out(class_i,class_j,sym_i)%elements(:,:) = 0.0_dp
                end do
            end do
        end do

        write(6,*) "Sigma 1 and 2:"
        call neci_flush(6)

        ! Loop over all classes.
        do i = 1, num_classes

            write(6,*) "Class:", i
            call neci_flush(6)

            call generate_first_full_string(string_i, ras, classes(i))

            ! Loop over all strings in class i.
            do
                call zero_factors_array(num_classes, factors)

                call encode_string(string_i, ilut_i)

                sym_i = get_abelian_sym(string_i)
                ! Shift the addresses so that each block of symmetries start with address 1.
                ind_i = classes(i)%address_map(get_address(classes(i), string_i))
                write(6,*) "ind:", ind_i
                call neci_flush(6)
                ind_i = ind_i - sum(classes(i)%num_sym(0:sym_i-1))

                ! To initialise the excitation generator in the first cycle.
                ilut_k(0) = -1

                ! Loop over all single excitations from string_i.
                do
                    nras1 = classes(i)%nelec_1
                    nras3 = classes(i)%nelec_3
                    
                    call gen_next_single_ex(string_i, ilut_i, string_k, ilut_k, ex1, nras1, nras3, ras, gen_store)

                    if (ilut_k(0) == -1) exit

                    class_k = ras%class_label(nras1, nras3)
                    sym_k = get_abelian_sym(string_k)
                    par_1 = get_single_parity(ilut_i, ex1(1), ex1(2)) 

                    ind_k = classes(k)%address_map(get_address(classes(k), string_k))
                    ind_k = ind_k - sum(classes(k)%num_sym(0:sym_k-1))

                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) + &
                            par_1*GetTMatEl(BRR(2*ex1(i)),BRR(2*ex1(j)))

                    do j = 1, tot_norbs
                        factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) + &
                                0.5_dp*par_1*get_umat_el(ptr_getumatel, ex1(2), j, j, ex1(1))
                    end do

                    ! To initialise the excitation generator in the first cycle.
                    ilut_j(0) = -1

                    ! Loop over all single excitations from string_k.
                    do
                        nras1 = classes(class_k)%nelec_1
                        nras3 = classes(class_k)%nelec_3
                    
                        call gen_next_single_ex(string_k, ilut_k, string_j, ilut_j, ex2, nras1, nras3, ras, gen_store)

                        if (ilut_j(0) == -1) exit
                        ! Only need to consider excitations where (ij) >= (jk).
                        if ( (ex2(2)-1)*tot_norbs + ex2(1) < (ex1(2)-1)*tot_norbs + ex1(1)) cycle

                        class_j = ras%class_label(nras1, nras3)
                        sym_j = get_abelian_sym(string_j)
                        par_2 = get_single_parity(ilut_j, ex2(1), ex2(2)) 

                        ind_j = classes(j)%address_map(get_address(classes(j), string_j))
                        ind_j = ind_j - sum(classes(j)%num_sym(0:sym_j-1))

                        ! Avoid overcounting for the case that the indices are the same.
                        if (ex1(1) == ex2(1) .and. ex1(2) == ex2(2)) then
                            factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                              0.5_dp*par_1*par_2*get_umat_el(ptr_getumatel, ex2(2), ex2(1), ex1(2), ex1(1))
                        else
                            factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                              par_1*par_2*get_umat_el(ptr_getumatel, ex2(2), ex2(1), ex1(2), ex1(1))
                        end if

                    end do

                end do

                ! The factors array has now been fully generated for string_i. Now we just have to add in the
                ! contirbution to vec_out from this state.

                ! Loop over all classes connected to the class of string_i.
                do j = 1, classes(i)%num_comb
                    class_j = classes(i)%allowed_combns(j)

                    ! The *required* symmetry,
                    sym_j = ieor(HFSym, sym_i)
                    ! If there are no states in this class with the required symmetry.
                    if (classes(class_j)%num_sym(sym_j) == 0) cycle

                        ! Loop over all classes connected to string_j.
                        do k = 1, classes(class_j)%num_comb
                            class_k = classes(class_j)%allowed_combns(k)
                            ! The *required* symmetry, sym_k = ieor(HFSym, sym_j) = sym_i.
                            sym_k = sym_i

                            ! If there are no states in this class with the required symmetry.
                            if (classes(class_k)%num_sym(sym_k) == 0) cycle

                            ! Finally, update the output vector.
                            ! Add in sigma_2.
                            vec_out(class_j, i, sym_j)%elements(:, ind_i) = &
                                vec_out(class_j, i, sym_j)%elements(:, ind_i) + &
                                matmul(vec_in(class_j, class_k, sym_j)%elements(:,:), factors(class_k, sym_k)%elements(:))

                            ! Add in sigma_1.
                            vec_out(i, class_j, sym_i)%elements(ind_i, :) = &
                                vec_out(i, class_j, sym_i)%elements(ind_i, :) + &
                                matmul(factors(class_k, sym_k)%elements(:), vec_in(class_k, class_j, sym_k)%elements(:,:))

                        end do ! Over all classes connected to string_j.
                    
                end do ! Over all classes connected to string_i.


                call generate_next_string(string_i, ras, classes(i), none_left)
                ! If no strings left in this class, go to the next class.
                if (none_left) exit
            end do ! Over all strings in this class.
        end do ! Over all classes.


        ! Next, calculate contirbution from the alpha-beta term, sigma_3.

        write(6,*) "sigma 3:"
        call neci_flush(6)
        
        ! Loop over all combinations of two spatial indices, (kl).
        do k = 1, tot_norbs
            do l = 1, tot_norbs

                write(6,*) "k, l:", k, l
                call neci_flush(6)

                ex1(1) = k
                ex1(2) = l

                do i = 1, num_classes
                    write(6,*) "class:", i
                    call neci_flush(6)
                    call generate_first_full_string(string_i, ras, classes(i))

                    nras1 = classes(i)%nelec_1
                    nras3 = classes(i)%nelec_3

                    ! Loop over all strings in class i.
                    do
                        sym_i = get_abelian_sym(string_i)
                        ind_i = classes(i)%address_map(get_address(classes(i), string_i))
                        write(6,*) "ind:", ind_i
                        call neci_flush(6)
                        ind_i = ind_i - sum(classes(i)%num_sym(0:sym_i-1))
                        call encode_string(string_i, ilut_i)

                        if ((.not. IsOcc(ilut_i,BRR(2*k))) .or. &
                                (IsOcc(ilut_i,BRR(2*l)) .and. k /= l)) then
                            r(i, sym_i)%elements(ind_i) = 0.0_dp
                            do m = 1, classes(i)%num_comb
                                class_m = classes(i)%allowed_combns(m)
                                sym_m = ieor(HFSym, sym_i)
                                if (classes(class_m)%num_sym(sym_m) == 0) cycle

                                c(i,class_m,sym_i)%elements(ind_i,:) = 0.0_dp
                            end do
                        else
                            call get_excit_details(string_i, ilut_i, ex1, ras, nras1, nras3, string_j, sym_j, class_j)
                            ind_j = classes(class_j)%address_map(get_address(classes(class_j), string_j))
                            ind_j = ind_j - sum(classes(class_j)%num_sym(0:sym_j-1))

                            r(i, sym_i)%elements(ind_i) = real(get_single_parity(ilut_i, k, l), dp)

                            do m = 1, classes(i)%num_comb
                                class_m = classes(i)%allowed_combns(m)
                                sym_m = ieor(HFSym, sym_i)
                                if (.not. class_comb_allowed(ras, classes(class_j), classes(class_m))) cycle
                                if (classes(class_m)%num_sym(sym_m) == 0) cycle 

                                write(6,*) "i, class_m, sym_i", i, class_m, sym_i
                                write(6,*) "size1:", size(c(i,class_m,sym_i)%elements, 2)
                                write(6,*) "size2:", size(vec_in(class_j,class_m,sym_j)%elements, 2)
                                call neci_flush(6)
                                c(i,class_m,sym_i)%elements(ind_i,:) = &
                                    vec_in(class_j,class_m,sym_j)%elements(ind_j,:)*r(i, sym_i)%elements(ind_i)
                            end do
                        end if

                        call generate_next_string(string_i, ras, classes(i), none_left)
                        ! If no strings left in this class, go to the next class.
                        if (none_left) exit
                    end do
                end do

                write(6,*) "init done."
                call neci_flush(6)

                do i = 1, num_classes

                    call generate_first_full_string(string_i, ras, classes(i))

                    ! Loop over all strings in class i.
                    do
                        call zero_factors_array(num_classes, factors)

                        call encode_string(string_i, ilut_i)

                        sym_i = get_abelian_sym(string_i)
                        ind_i = classes(i)%address_map(get_address(classes(i), string_i))
                        write(6,*) "ind:", ind_i
                        call neci_flush(6)
                        ind_i = ind_i - sum(classes(i)%num_sym(0:sym_i-1))

                        ! To initialise the excitation generator in the first cycle.
                        ilut_j(0) = -1

                        ! Loop over all single excitations from string_i.
                        do
                            nras1 = classes(i)%nelec_1
                            nras3 = classes(i)%nelec_3
                            
                            call gen_next_single_ex(string_i, ilut_i, string_j, ilut_j, ex1, nras1, nras3, ras, gen_store)

                            if (ilut_j(0) == -1) exit

                            class_j = ras%class_label(nras1, nras3)
                            sym_j = get_abelian_sym(string_k)
                            par_1 = get_single_parity(ilut_i, ex1(1), ex1(2)) 

                            ind_j = classes(j)%address_map(get_address(classes(j), string_j))
                            ind_j = ind_j - sum(classes(j)%num_sym(0:sym_j-1))

                            factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                              par_1*get_umat_el(ptr_getumatel, ex1(2), ex1(1), k, l)

                        end do

                        do j = 1, classes(i)%num_comb

                            class_j = classes(i)%allowed_combns(j)
                            sym_j = ieor(HFSym, sym_i)

                            do ind_j = 1, classes(class_j)%num_sym(sym_j)

                                if (.not. abs(r(class_j, sym_j)%elements(ind_j)) > 0.0_dp) cycle

                                do m = 1, classes(class_j)%num_comb

                                    class_m = classes(class_j)%allowed_combns(m)
                                    sym_m = sym_i

                                    if (classes(class_m)%num_sym(sym_m) == 0) cycle

                                    v(class_j, sym_j)%elements(ind_j) = dot_product(factors(class_m,sym_m)%elements(:), &
                                                                     c(class_j, class_m, sym_j)%elements(ind_j, :))

                                end do ! Over all classes connected to class_j.

                                vec_out(class_j, i, sym_j)%elements(ind_j, ind_i) = &
                                    vec_out(class_j, i, sym_j)%elements(ind_j, ind_i) + v(class_j, sym_j)%elements(ind_j)

                            end do ! Over all states in class_j with symmetry sym_j.

                        end do ! Over all classes.


                        call generate_next_string(string_i, ras, classes(i), none_left)
                        ! If no strings left in this class, go to the next class.
                        if (none_left) exit

                    end do ! Over all states in a class.
                
                end do ! Over all classes.

            end do
        end do

    end subroutine perform_multiplication

    subroutine encode_string(string, ilut)

        integer, intent(in) :: string(tot_nelec)
        integer(n_int), intent(out) :: ilut(0:NIfTot)
        integer :: i, pos

        ilut = 0

        do i = 1, tot_nelec
            pos = (BRR(string(i)*2) - 1) / bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(BRR(string(i)*2)-1, bits_n_int))
        end do
                
    end subroutine encode_string

    subroutine get_excit_details(string_i, ilut_i, ex, ras, nras1, nras3, string_j, sym_j, class_j)

        integer, intent(in) :: string_i(tot_nelec)
        integer(n_int), intent(in) :: ilut_i(0:NIfTot)
        integer, intent(in) :: ex(2)
        type(ras_parameters), intent(in) :: ras
        integer, intent(in) :: nras1, nras3
        integer, intent(out) :: string_j(tot_nelec)
        integer, intent(out) :: sym_j, class_j
        integer :: i, new_1, new_3

        new_1 = nras1
        new_3 = nras3

        if (ex(1) <= ras%size_1) then
            new_1 = new_1 - 1
        else if (ex(1) > ras%size_1 + ras%size_2) then
            new_3 = new_3 - 1
        end if

        if (ex(2) <= ras%size_1) then
            new_1 = new_1 + 1
        else if (ex(2) > ras%size_1 + ras%size_2) then
            new_3 = new_3 + 1
        end if

        class_j = ras%class_label(new_1, new_3)

        string_j = string_i
        do i = 1, tot_nelec
            if (string_j(i) == ex(1))then
                string_j(i) = ex(2)
                exit
            end if
        end do
        call sort(string_j)

        sym_j = get_abelian_sym(string_j)

    end subroutine get_excit_details

    subroutine zero_factors_array(num_classes, factors)

        integer, intent(in) :: num_classes
        type(ras_factors), intent(inout) :: factors(num_classes, 0:7)
        integer :: class_i, sym_i

        do class_i = 1, num_classes
            do sym_i = 0, 7
                if (allocated(factors(class_i,sym_i)%elements)) &
                    factors(class_i,sym_i)%elements(:) = 0.0_dp
            end do
        end do

    end subroutine zero_factors_array

    subroutine gen_next_single_ex(string_i, ilut_i, string_k, ilut_k, &
                                  ex, nras1, nras3, ras, gen_store)

        integer, intent(in) :: string_i(tot_nelec)
        integer(n_int), intent(in) :: ilut_i(0:NIfTot)
        integer, intent(out) :: string_k(tot_nelec)
        integer(n_int), intent(out) :: ilut_k(0:NIfTot)
        integer, intent(out) :: ex(2)
        integer, intent(inout) :: nras1, nras3
        type(ras_parameters), intent(in) :: ras
        type(simple_excit_store), intent(inout), target :: gen_store

        integer, pointer :: i, j
        integer :: orb1, orb2, temp_1, temp_3

        ! Map the local variables onto the store.
        i => gen_store%i;
        j => gen_store%j

        ! Initialise the generator.
        if (ilut_k(0) == -1) then
            i = 1
            j = 0
        end if

        ! Find the next possible single excitation. Loop over electrons and vacant orbitals.
        ! Interrupt loop when we find what we need.
        do i = i, tot_nelec

            if (j == 0) then
                orb1 = BRR(string_i(i)*2)
            end if

            j = j + 1
            do j = j, tot_norbs

                orb2 = BRR(2*j)

                ! Cannot excite to an occupied orbital, unless it is the orbital that we are
                ! exciting from.
                if (IsOcc(ilut_i, orb2) .and. (string_i(i) /= j)) cycle

                temp_1 = nras1
                temp_3 = nras3

                ! Store the values of nras1 and nras 3 for the new string.
                if (ex(1) <= ras%size_1) then
                    temp_1 = temp_1 - 1
                else if (ex(1) > ras%size_1 + ras%size_2) then
                    temp_3 = temp_3 - 1
                end if

                if (ex(2) <= ras%size_1) then
                    temp_1 = temp_1 + 1
                else if (ex(2) > ras%size_1 + ras%size_2) then
                    temp_3 = temp_3 + 1
                end if

                ! We don't have to consider excitations to outside the ras space.
                if (.not. class_allowed(ras, temp_1, temp_3)) cycle

                ! If the class is allowed then the final test is passed, so update nras1 and nras3
                ! ready to be output.
                nras1 = temp_1
                nras3 = temp_3

                ! Generate the determinant and interrupt the loop.
                ilut_k = ilut_i
                clr_orb(ilut_k, orb1)
                set_orb(ilut_k, orb2)

                ex(1) = string_i(i)
                ex(2) = j

                ! Encode new string.
                string_k = string_i
                string_k(i) = j
                call sort(string_k)

                return
            end do
            j = 0 ! Reset loop.
        end do

        ilut_k(0) = -1

    end subroutine gen_next_single_ex

end module direct_ci
