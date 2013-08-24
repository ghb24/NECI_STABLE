#include "macros.h"

module direct_ci

    use bit_rep_data, only: NIfD
    use DetBitOps, only: get_single_parity
    use enumerate_excitations, only: simple_excit_store
    use Integrals_neci, only: get_umat_el
    use IntegralsData, only: ptr_getumatel
    use OneEInts, only: GetTMatEl
    use ras

    implicit none

    type direct_ci_excit
        integer(sp), allocatable, dimension(:) :: excit_ind
        integer(sp), allocatable, dimension(:) :: par
        integer(sp), allocatable, dimension(:,:) :: orbs
        integer(sp) :: nexcit
    end type direct_ci_excit

contains

    subroutine perform_multiplication(ras, classes, ras_strings, ras_iluts, ras_excit, vec_in, vec_out)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        integer(sp), intent(in) :: ras_strings(-1:tot_nelec, core_ras%num_strings)
        integer(n_int), intent(in) :: ras_iluts(0:NIfD, core_ras%num_strings)
        type(direct_ci_excit), intent(in) :: ras_excit(core_ras%num_strings)
        type(ras_vector), intent(in) :: vec_in(ras%num_classes, ras%num_classes, 0:7)
        type(ras_vector), intent(inout) :: vec_out(ras%num_classes, ras%num_classes, 0:7)

        integer(sp) :: class_i, class_j, class_k, class_m, par_1, par_2
        integer(sp) :: i, j, k, l, m
        integer(sp) :: excit_j, excit_k
        integer(sp) :: ind_i, ind_j, ind_k, full_ind_j, full_ind_k
        integer(sp) :: nras1, nras3, min_ind, max_ind
        integer(sp) :: sym_i, sym_j, sym_k, sym_m, HFSym_sp
        integer(sp) :: ex1(2), ex2(2)

        type(ras_vector) :: c(ras%num_classes, ras%num_classes)
        real(dp), allocatable, dimension(:) :: factor
        integer(sp) :: string_i(tot_nelec), string_j(tot_nelec), string_k(tot_nelec)
        integer(n_int) :: ilut_i(0:NIfD), ilut_j(0:NIfD), ilut_k(0:NIfD)
        type(ras_factors) :: factors(ras%num_classes, 0:7), v(ras%num_classes, 0:7)
        type(ras_factors) :: r(ras%num_classes)
        type(simple_excit_store), target :: gen_store_1, gen_store_2
        logical :: none_left, in_ras_space

        do class_i = 1, ras%num_classes
            do sym_i = 0, 7
                allocate(factors(class_i,sym_i)%elements(1:classes(class_i)%num_sym(sym_i)))
                allocate(v(class_i,sym_i)%elements(1:classes(class_i)%num_sym(sym_i)))
            end do
            allocate(r(class_i)%elements(1:classes(class_i)%class_size))
        end do

        do class_i = 1, ras%num_classes
            do class_j = 1, ras%num_classes
                allocate(c(class_i,class_j)%elements(1:classes(class_i)%class_size, &
                           1:classes(class_j)%class_size))
                c(class_i,class_j)%elements = 0.0_dp
            end do
        end do

        do class_i = 1, ras%num_classes
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)
                do sym_i = 0, 7
                    sym_j = ieor(HFSym_sp, sym_i)

                    if (classes(class_i)%num_sym(sym_i) == 0 .or. classes(class_j)%num_sym(sym_j) == 0) cycle
                    vec_out(class_i,class_j,sym_i)%elements(:,:) = 0.0_dp
                end do
            end do
        end do

        HFSym_sp = int(HFSym,sp)

        ! Loop over all strings.
        do i = 1, ras%num_strings

            call zero_factors_array(ras%num_classes, factors)

            ! Get info for string_i.
            sym_i = ras_strings(-1,i)
            class_i = ras_strings(0,i)
            ind_i = i - sum(classes(1:class_i-1)%class_size) - sum(classes(class_i)%num_sym(0:sym_i-1))

            ! Loop over all single excitations from string_i (-> string_k).
            do excit_k = 1, ras_excit(i)%nexcit

                full_ind_k = ras_excit(i)%excit_ind(excit_k)
                
                sym_k = ras_strings(-1,full_ind_k)
                class_k = ras_strings(0,full_ind_k)
                par_1 = ras_excit(i)%par(excit_k)
                ex1 = ras_excit(i)%orbs(:,excit_k)

                ind_k = full_ind_k - sum(classes(class_k)%num_sym(0:sym_k-1))

                factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) + &
                        par_1*GetTMatEl(BRR(2*ex1(2)),BRR(2*ex1(1)))

                do j = 1, ex1(2)-1
                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) - &
                            par_1*get_umat_el(ptr_getumatel, BRR(2*ex1(2))/2, BRR(2*j)/2, &
                                                             BRR(2*j)/2, BRR(2*ex1(1))/2)
                end do

                if (ex1(2) > ex1(1)) then
                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) - &
                            par_1*get_umat_el(ptr_getumatel, BRR(2*ex1(2))/2, BRR(2*ex1(2))/2, &
                                                             BRR(2*ex1(2))/2, BRR(2*ex1(1))/2)
                else if (ex1(2) == ex1(1)) then
                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) - &
                            0.5_dp*par_1*get_umat_el(ptr_getumatel, BRR(2*ex1(2))/2, BRR(2*ex1(2))/2, &
                                                                    BRR(2*ex1(2))/2, BRR(2*ex1(1))/2)
                end if

                ! Loop over all single excitations from string_k (-> string_j).
                do excit_j = 1, ras_excit(full_ind_k)%nexcit

                    ex2 = ras_excit(full_ind_k)%orbs(:,excit_j)
                
                    ! Only need to consider excitations where (ij) >= (kl).
                    if ( (ex2(2)-1)*tot_norbs + ex2(1) < (ex1(2)-1)*tot_norbs + ex1(1)) cycle

                    full_ind_j = ras_excit(full_ind_k)%excit_ind(excit_j)
                    
                    sym_k = ras_strings(-1,full_ind_j)
                    class_k = ras_strings(0,full_ind_j)
                    par_1 = ras_excit(full_ind_k)%par(excit_j)

                    ind_j = full_ind_j - sum(classes(class_j)%num_sym(0:sym_j-1))

                    ! Avoid overcounting for the case that the indices are the same.
                    if (ex1(1) == ex2(1) .and. ex1(2) == ex2(2)) then
                        factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                          0.5_dp*par_1*par_2*get_umat_el(ptr_getumatel, BRR(2*ex2(2))/2, BRR(2*ex1(2))/2, &
                                                                        BRR(2*ex2(1))/2, BRR(2*ex1(1))/2)
                    else
                        factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                          par_1*par_2*get_umat_el(ptr_getumatel, BRR(2*ex2(2))/2, BRR(2*ex1(2))/2, &
                                                                 BRR(2*ex2(1))/2, BRR(2*ex1(1))/2)
                    end if

                end do

            end do

            ! The factors array has now been fully generated for string_i. Now we just have to
            ! add in the contribution to vec_out from this state.

            ! Loop over all classes connected to the class of string_i.
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)

                ! The *required* symmetry,
                sym_j = ieor(HFSym_sp, sym_i)
                ! If there are no states in this class with the required symmetry.
                if (classes(class_j)%num_sym(sym_j) == 0) cycle

                ! Loop over all classes connected to string_j.
                do k = 1, classes(class_j)%num_comb
                    class_k = classes(class_j)%allowed_combns(k)
                    ! The *required* symmetry, sym_k = ieor(HFSym_sp, sym_j) = sym_i.
                    sym_k = sym_i

                    ! If there are no states in this class with the required symmetry.
                    if (classes(class_k)%num_sym(sym_k) == 0) cycle

                    ! Finally, update the output vector.
                    ! Add in sigma_2.
                    vec_out(class_j, class_i, sym_j)%elements(:, ind_i) = &
                        vec_out(class_j, class_i, sym_j)%elements(:, ind_i) + &
                        matmul(vec_in(class_j, class_k, sym_j)%elements(:,:), factors(class_k, sym_k)%elements(:))

                    ! Add in sigma_1.
                    vec_out(class_i, class_j, sym_i)%elements(ind_i, :) = &
                        vec_out(class_i, class_j, sym_i)%elements(ind_i, :) + &
                        matmul(factors(class_k, sym_k)%elements(:), vec_in(class_k, class_j, sym_k)%elements(:,:))

                end do ! Over all classes connected to string_j.
                
            end do ! Over all classes connected to string_i.

        end do ! Over all strings.

        ! Next, calculate contirbution from the alpha-beta term, sigma_3.

        ! Loop over all combinations of two spatial indices, (kl).
        do k = 1, tot_norbs
            do l = 1, tot_norbs

                ! Zero arrays.
                do class_i = 1, ras%num_classes
                    r(class_i)%elements = 0.0_dp
                    do class_j = 1, ras%num_classes
                        c(class_i,class_j)%elements = 0.0_dp
                    end do
                end do

                ex1(1) = l
                ex1(2) = k

                ! Loop over all strings.
                do i = 1, ras%num_strings

                    ! Get info for string_i.
                    sym_i = ras_strings(-1,i)
                    class_i = ras_strings(0,i)
                    nras1 = classes(class_i)%nelec_1
                    nras3 = classes(class_i)%nelec_3
                    string_i = ras_strings(1:tot_nelec,i)
                    ilut_i = ras_iluts(:,i)
                    ind_i = i - sum(classes(1:class_i-1)%class_size) - sum(classes(class_i)%num_sym(0:sym_i-1))

                    if ( IsOcc(ilut_i,l) .and. &
                            (.not. (IsOcc(ilut_i,k) .and. k /= l)) ) then
                        ! Temporarily set this for get_excit_details to use.
                        sym_j = sym_i
                        call get_excit_details(string_i, ex1, ras, nras1, nras3, string_j, sym_j, class_j, in_ras_space)
                        if (in_ras_space) then
                            ind_j = classes(class_j)%address_map(get_address(classes(class_j), string_j))
                            ind_j = ind_j - sum(classes(class_j)%num_sym(0:sym_j-1))

                            r(i)%elements(ind_i) = real(get_single_parity(ilut_i, int(l,sizeof_int), int(k,sizeof_int)), dp)

                            do m = 1, classes(class_j)%num_comb
                                class_m = classes(class_j)%allowed_combns(m)
                                sym_m = ieor(HFSym_sp, sym_j)
                                if (classes(class_m)%num_sym(sym_m) /= 0) then
                                    min_ind = sum(classes(class_m)%num_sym(0:sym_m-1)) + 1
                                    max_ind = sum(classes(class_m)%num_sym(0:sym_m))
                                    c(i,class_m)%elements(ind_i,min_ind:max_ind) = &
                                        c(i,class_m)%elements(ind_i,min_ind:max_ind) + &
                                        vec_in(class_j,class_m,sym_j)%elements(ind_j,:)*r(i)%elements(ind_i)
                                end if
                            end do
                        end if
                    end if

                end do

                ! Loop over all strings.
                do i = 1, ras%num_strings

                    call zero_factors_array(ras%num_classes, factors)

                    ! Get info for string_i.
                    sym_i = ras_strings(-1,i)
                    class_i = ras_strings(0,i)
                    ind_i = i - sum(classes(1:class_i-1)%class_size) - sum(classes(class_i)%num_sym(0:sym_i-1))

                    ! Loop over all single excitations from string_i (-> string_j).
                    do excit_j = 1, ras_excit(i)%nexcit

                        full_ind_j = ras_excit(i)%excit_ind(excit_j)
                        
                        sym_j = ras_strings(-1,full_ind_j)
                        class_j = ras_strings(0,full_ind_j)
                        par_1 = ras_excit(i)%par(excit_j)
                        ex1 = ras_excit(i)%orbs(:,excit_j)

                        ind_j = full_ind_j - sum(classes(class_j)%num_sym(0:sym_j-1))

                        factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                          par_1*get_umat_el(ptr_getumatel, BRR(2*ex1(2))/2, BRR(2*k)/2, BRR(2*ex1(1))/2, BRR(2*l)/2)

                    end do

                    ! Loop over all classes allowed with class_i.
                    do j = 1, classes(class_i)%num_comb

                        class_j = classes(class_i)%allowed_combns(j)
                        sym_j = ieor(HFSym_sp, sym_i)

                        do ind_j = 1, classes(class_j)%num_sym(sym_j)

                            full_ind_j = ind_j + sum(classes(class_j)%num_sym(0:sym_j-1))

                            if (.not. abs(r(class_j)%elements(full_ind_j)) > 0.0_dp) cycle

                            call zero_factors_array(ras%num_classes, v)

                            do class_m = 1, ras%num_classes

                                do sym_m = 0, 7

                                    min_ind = sum(classes(class_m)%num_sym(0:sym_m-1)) + 1
                                    max_ind = sum(classes(class_m)%num_sym(0:sym_m))

                                    v(class_j, sym_j)%elements(ind_j) = v(class_j, sym_j)%elements(ind_j) + &
                                        dot_product(factors(class_m,sym_m)%elements(:), &
                                        c(class_j, class_m)%elements(full_ind_j, min_ind:max_ind))

                                end do

                            end do ! Over all classes connected to class_j.

                            vec_out(class_j, class_i, sym_j)%elements(ind_j, ind_i) = &
                                vec_out(class_j, class_i, sym_j)%elements(ind_j, ind_i) + &
                                v(class_j, sym_j)%elements(ind_j)

                        end do ! Over all states in class_j with symmetry sym_j.

                    end do ! Over all classes allowed with class_i.

                end do ! Over all strings.

            end do
        end do

    end subroutine perform_multiplication

    subroutine encode_string(string, ilut)

        integer(sp), intent(in) :: string(tot_nelec)
        integer(n_int), intent(out) :: ilut(0:NIfD)
        integer(sp) :: i, pos

        ilut = 0

        do i = 1, tot_nelec
            pos = (string(i)-1)/bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(string(i)-1, bits_n_int))
        end do
                
    end subroutine encode_string

    subroutine get_excit_details(string_i, ex, ras, nras1, nras3, string_j, sym_j, class_j, in_ras_space)

        integer(sp), intent(in) :: string_i(tot_nelec)
        integer(sp), intent(in) :: ex(2)
        type(ras_parameters), intent(in) :: ras
        integer(sp), intent(in) :: nras1, nras3
        integer(sp), intent(out) :: string_j(tot_nelec)
        integer(sp), intent(inout) :: sym_j
        integer(sp), intent(out) :: class_j
        logical, intent(out) :: in_ras_space
        integer(sp) :: i, new1, new3
        integer(sp) :: sym_prod

        new1 = nras1
        new3 = nras3

        if (ex(1) <= ras%size_1) then
            new1 = new1 - 1
        else if (ex(1) > ras%size_1 + ras%size_2) then
            new3 = new3 - 1
        end if

        if (ex(2) <= ras%size_1) then
            new1 = new1 + 1
        else if (ex(2) > ras%size_1 + ras%size_2) then
            new3 = new3 + 1
        end if

        if (class_allowed(ras, new1, new3)) then
            in_ras_space = .true.
        else
            in_ras_space = .false.
            return
        end if

        class_j = ras%class_label(new1, new3)

        string_j = string_i
        do i = 1, tot_nelec
            if (string_j(i) == ex(1)) then
                string_j(i) = ex(2)
                exit
            end if
        end do
        call sort(string_j)

        sym_prod = int(ieor(G1(BRR(ex(1)*2))%Sym%S, G1(BRR(ex(1)*2))%Sym%S),sp)
        sym_j = ieor(sym_j,sym_prod)

    end subroutine get_excit_details

    subroutine zero_factors_array(num_classes, factors)

        integer(sp), intent(in) :: num_classes
        type(ras_factors), intent(inout) :: factors(num_classes, 0:7)
        integer(sp) :: class_i, sym_i

        do class_i = 1, num_classes
            do sym_i = 0, 7
                if (allocated(factors(class_i,sym_i)%elements)) &
                    factors(class_i,sym_i)%elements(:) = 0.0_dp
            end do
        end do

    end subroutine zero_factors_array

    subroutine gen_next_single_ex(string_i, ilut_i, nras1, nras3, ind, par, ex, ras, classes, gen_store, tgen, tcount)

        integer(sp), intent(in) :: string_i(tot_nelec)
        integer(n_int), intent(in) :: ilut_i(0:NIfD)
        integer(sp), intent(in) :: nras1, nras3
        integer(sp), intent(out) :: ind, par
        integer(sp), intent(out) :: ex(2)
        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        type(simple_excit_store), intent(inout), target :: gen_store
        logical, intent(inout) :: tgen
        logical, intent(in) :: tcount

        integer, pointer :: i, j
        integer(sp) :: orb1, orb2, temp1, temp3, class_k
        integer(sp) :: string_k(tot_nelec)

        ! Map the local variables onto the store.
        i => gen_store%i;
        j => gen_store%j

        ! Initialise the generator.
        if (tgen) then
            i = 1
            j = 0
            tgen = .false.
        end if

        ! Find the next possible single excitation. Loop over electrons and vacant orbitals.
        ! Interrupt loop when we find what we need.
        do i = i, tot_nelec

            orb1 = string_i(i)

            j = j + 1
            do j = j, tot_norbs

                orb2 = j

                ! Cannot excite to an occupied orbital, unless it is the orbital that we are
                ! exciting from.
                if (IsOcc(ilut_i, orb2) .and. (string_i(i) /= j)) cycle

                ex(1) = string_i(i)
                ex(2) = j

                temp1 = nras1
                temp3 = nras3

                ! Store the values of nras1 and nras 3 for the new string.
                if (ex(1) <= ras%size_1) then
                    temp1 = temp1 - 1
                else if (ex(1) > ras%size_1 + ras%size_2) then
                    temp3 = temp3 - 1
                end if

                if (ex(2) <= ras%size_1) then
                    temp1 = temp1 + 1
                else if (ex(2) > ras%size_1 + ras%size_2) then
                    temp3 = temp3 + 1
                end if

                ! We don't have to consider excitations to outside the ras space.
                if (.not. class_allowed(ras, temp1, temp3)) cycle

                ! If we get to this point then (orb1, orb2) is a valid excitation. As an
                ! optimisation, if we are only counting the excitations then return now.
                if (tcount) return

                ! Encode new string.
                string_k = string_i
                string_k(i) = j
                call sort(string_k)

                class_k = ras%class_label(temp1, temp3)
                ind = classes(class_k)%address_map(get_address(classes(class_k), string_k))
                ind = ind + sum(classes(1:class_k-1)%class_size)
                par = get_single_parity(ilut_i, int(ex(1),sizeof_int), int(ex(2),sizeof_int)) 

                return
            end do
            j = 0 ! Reset loop.
        end do

        tgen = .true.

    end subroutine gen_next_single_ex

    subroutine create_direct_ci_arrays(ras, classes, ras_strings, ras_iluts, ras_excit)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(inout), allocatable, dimension(:) :: classes
        integer(sp), intent(out) :: ras_strings(-1:tot_nelec, core_ras%num_strings)
        integer(n_int), intent(out) :: ras_iluts(0:NIfD, core_ras%num_strings)
        type(direct_ci_excit), intent(out) :: ras_excit(core_ras%num_strings)
        integer(n_int) :: ilut_i(0:NIfD)
        integer(sp) :: class_i, ind, new_ind, counter
        integer(sp) :: nras1, nras3, par
        integer(sp) :: ex(2)
        integer(sp) :: string_i(tot_nelec)
        logical :: none_left, tgen
        type(simple_excit_store), target :: gen_store_1

        ilut_i = 0

        ! Loop over all classes.
        do class_i = 1, ras%num_classes
            nras1 = classes(class_i)%nelec_1
            nras3 = classes(class_i)%nelec_3
            call generate_first_full_string(string_i, ras, classes(class_i))
            do
                ind = classes(class_i)%address_map(get_address(classes(class_i), string_i))
                ind = ind + sum(classes(1:class_i-1)%class_size)

                ras_strings(1:tot_nelec,ind) = string_i
                ras_strings(-1,ind) = get_abelian_sym(string_i)
                ras_strings(0,ind) = class_i
                call encode_string(string_i, ilut_i)
                ras_iluts(:,ind) = ilut_i

                ! Now count the number of excitations.
                counter = 0
                ! Initialise the excitation generator.
                tgen = .true.
                do
                    call gen_next_single_ex(string_i, ilut_i, nras1, nras3, new_ind, par, &
                            ex, ras, classes, gen_store_1, tgen, .true.)
                    if (tgen) exit
                    counter = counter + 1
                end do

                ras_excit(ind)%nexcit = counter

                ! Now that they have been counted, allocate the excitation array and store them.
                allocate(ras_excit(ind)%excit_ind(counter))
                allocate(ras_excit(ind)%par(counter))
                allocate(ras_excit(ind)%orbs(2,counter))

                counter = 0
                tgen = .true.
                do
                    call gen_next_single_ex(string_i, ilut_i, nras1, nras3, new_ind, par, &
                            ex, ras, classes, gen_store_1, tgen, .false.)
                    if (tgen) exit
                    counter = counter + 1
                    ras_excit(ind)%excit_ind(counter) = new_ind
                    ras_excit(ind)%par(counter) = par
                    ras_excit(ind)%orbs(:,counter) = ex
                end do

                call generate_next_string(string_i, ras, classes(class_i), none_left)
                ! If no strings left in this class, go to the next class.
                if (none_left) exit
            end do
        end do

    end subroutine create_direct_ci_arrays

end module direct_ci
