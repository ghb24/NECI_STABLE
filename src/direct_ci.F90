#include "macros.h"

! This module contains routines to perform the multiplication of the Hamiltonian by a vector
! using the direct CI approach. Specifically, it takes the approach of Olsen, Roos, Jorgensen
! and Jensen in J. Chem. Phys. 89, 2185 (1988). The code follows the description they give,
! and ocassionally refers to it in comments. This approach is also described in Molecular
! Electronic-Structure Theory (Helgaker, Jorgensen, Olsen) in section 11.8.4. The
! multiplication is performed for the Hamiltonian in a RAS space, and abelian point-group
! symmetry is used.

! *Note*, only ms=0 subspaces are implemented currently (it is assumed that the possible alpha
! and beta strings are the same). It also shouldn't be used in UHF calculations.

! Also note that the diagonal elements of H in this approach don't contain ecore, so this
! must be added to the result separately.

! Because alpha and beta strings are used, the ordering of orbitals is different to the rest
! of the code. Therefore, some vector components will differ by a minus sign.

! The vector components are stored in matrices rather than a vector, due to using alpha and
! beta strings. Each matrix will refer to a pair of RAS classes and symmetries for the
! corresponding alpha and beta strings. To convert an input vector to this block matrix form,
! use the transfer_from_block_form and transfer_to_block_form routines in this module. These
! can then be fed into the perform_multiplcation routine. However, note that the differing
! signs due to the altered ordering (see above) are not corrected for by these routines.

! See comments in the perform_multiplication routine and the ras module for which routines
! to call to create the required arrays.

! TODO: - The calls to obtain the integrals are very slow, probably due to having to go through
! the function UMatInd. A lot of checks performed here are unnecessary in this case, so a
! separate UMatInd function should be written for direct CI.
!       - The sort in get_excit_details is probably slow and definitely unnecessary.

module direct_ci

    use bit_rep_data, only: NIfD
    use DetBitOps, only: get_single_parity
    use enumerate_excitations, only: simple_excit_store
    use procedure_pointers, only: get_umat_el
    use OneEInts, only: GetTMatEl
    use ras
    use ras_data

    implicit none

contains

    subroutine perform_multiplication(ras, classes, ras_strings, ras_iluts, ras_excit, vec_in, vec_out)

        ! In: ras - Details of the RAS space being used.
        ! In: classes - Details of the RAS classes for the above RAS space. This array (and also the ras
        !     variable) should be created by the the initialise_ras_space routine in the ras module.
        ! In: ras_strings - All of the alpha (and beta) strings for the RAS space.
        ! In: ras_iluts - All of the alpha (and beta) iluts for the RAS space.
        ! In: ras_excit - For each RAS string, store the addresses of strings which can be excited to
        !     by a single excitation. This array, and the two above, are created by the
        !     create_direct_ci_arrays routine in this module.
        ! In: vec_in - The input vector to be multiplied, in block form. A standard vector can be bought
        !     to this form by using the transfer_to_block_form routine in this module.
        ! Inout: vec_out - The resulting output vector after multiplication.

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        integer, intent(in) :: ras_strings(-1:tot_nelec, ras%num_strings)
        integer(n_int), intent(in) :: ras_iluts(0:NIfD, ras%num_strings)
        type(direct_ci_excit), intent(in) :: ras_excit(ras%num_strings)
        type(ras_vector), intent(in) :: vec_in(ras%num_classes, ras%num_classes, 0:7)
        type(ras_vector), intent(inout) :: vec_out(ras%num_classes, ras%num_classes, 0:7)

        type(ras_factors) :: factors(ras%num_classes, 0:7)
        real(dp), allocatable, dimension(:) :: alpha_beta_fac
        integer, allocatable, dimension(:) :: r
        real(dp), allocatable, dimension(:,:) :: c
        integer, allocatable, dimension(:,:) :: excit_info
        integer :: string_i(tot_nelec), string_j(tot_nelec)
        integer(n_int) :: ilut_i(0:NIfD)
        logical :: in_ras_space
        real(dp) :: v

        integer :: class_i, class_j, class_k, class_m
        integer :: excited_class, excited_sym
        integer :: par_1, par_2
        integer :: i, j, k, l, m
        integer :: excit_j, excit_k
        integer :: ind_i, ind_j, ind_k, full_ind_j, full_ind_k
        integer :: nras1, nras3, min_ind, max_ind
        integer :: sym_i, sym_j, sym_k, sym_m
        integer :: ex1(2), ex2(2)
        integer :: BRR_ex1(2), BRR_ex2(2)
        integer :: ierr

        ! Initialisation - allocate arrays.

        do class_i = 1, ras%num_classes
            do sym_i = 0, 7
                allocate(factors(class_i,sym_i)%elements(1:classes(class_i)%num_sym(sym_i)), &
                     stat = ierr)
            end do
        end do

        do class_i = 1, ras%num_classes
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)
                do sym_i = 0, 7
                    sym_j = ieor(HFSym_ras, sym_i)
                    if (classes(class_i)%num_sym(sym_i) == 0 .or. &
                        classes(class_j)%num_sym(sym_j) == 0) cycle
                    vec_out(class_i,class_j,sym_i)%elements(:,:) = 0.0_dp
                end do
            end do
        end do

        allocate(alpha_beta_fac(ras%num_strings))
        allocate(c(ras%num_strings,ras%num_strings))
        allocate(r(ras%num_strings))
        alpha_beta_fac = 0.0_dp
        c = 0.0_dp
        r = 0

        allocate(excit_info(2, ras%num_strings))
        excit_info = 0


        ! Beta and beta-beta (sigma_1) and alpha and alpha-alpha (sigma_2) contributions.
        ! See Eq. 20 for algorithm.
        ! Loop over all strings.
        do i = 1, ras%num_strings

            ! Set factors array equal to zero.
            call zero_factors_array(ras%num_classes, factors)

            ! Get info for string_i.
            sym_i = ras_strings(-1,i)
            class_i = ras_strings(0,i)
            ind_i = i - ras%cum_classes(class_i) - classes(class_i)%cum_sym(sym_i)

            ! Loop over all single excitations from string_i (-> string_k).
            do excit_k = 1, ras_excit(i)%nexcit

                ! The full address of string k.
                full_ind_k = ras_excit(i)%excit_ind(excit_k)
                
                sym_k = ras_strings(-1,full_ind_k)
                class_k = ras_strings(0,full_ind_k)
                par_1 = ras_excit(i)%par(excit_k)
                ex1 = ras_excit(i)%orbs(:,excit_k)
                BRR_ex1(1) = BRR(2*ex1(1))/2
                BRR_ex1(2) = BRR(2*ex1(2))/2

                ! The shifted address (shifted so that the first string with this symmetry and
                ! in this RAS class has address 1).
                ind_k = full_ind_k - ras%cum_classes(class_k) - classes(class_k)%cum_sym(sym_k)

                ! In Eq. 20, this term is sgn(kl)*h_{kl}.
                factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) + &
                        par_1*GetTMatEl(BRR_ex1(2)*2,BRR_ex1(1)*2)

                ! The following lines add in g_{kl} from Eq. 29.
                do j = 1, ex1(2)-1
                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) - &
                            par_1*get_umat_el(BRR_ex1(2), BRR(2*j)/2, BRR(2*j)/2, BRR_ex1(1))
                end do

                if (ex1(2) > ex1(1)) then
                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) - &
                            par_1*get_umat_el(BRR_ex1(2), BRR_ex1(2), BRR_ex1(2), BRR_ex1(1))
                else if (ex1(2) == ex1(1)) then
                    factors(class_k, sym_k)%elements(ind_k) = factors(class_k, sym_k)%elements(ind_k) - &
                            0.5_dp*par_1*get_umat_el(BRR_ex1(2), BRR_ex1(2), BRR_ex1(2), BRR_ex1(1))
                end if

                ! Loop over all single excitations from string_k (-> string_j).
                do excit_j = 1, ras_excit(full_ind_k)%nexcit

                    ex2 = ras_excit(full_ind_k)%orbs(:,excit_j)

                    ! Only need to consider excitations where (ij) >= (kl) (see Eq. 28).
                    if ((ex2(2)-1)*tot_norbs + ex2(1) < (ex1(2)-1)*tot_norbs + ex1(1)) cycle

                    BRR_ex2(1) = BRR(2*ex2(1))/2
                    BRR_ex2(2) = BRR(2*ex2(2))/2

                    ! The full address for string j.
                    full_ind_j = ras_excit(full_ind_k)%excit_ind(excit_j)
                    
                    sym_j = ras_strings(-1,full_ind_j)
                    class_j = ras_strings(0,full_ind_j)
                    par_2 = ras_excit(full_ind_k)%par(excit_j)

                    ! The shifted address (shifted so that the first string with this symmetry and
                    ! in this RAS class has address 1).
                    ind_j = full_ind_j - ras%cum_classes(class_j) - classes(class_j)%cum_sym(sym_j)

                    ! Avoid overcounting for the case that the indices are the same.
                    if (ex1(1) == ex2(1) .and. ex1(2) == ex2(2)) then
                        factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                          0.5_dp*par_1*par_2*get_umat_el(BRR_ex2(2), BRR_ex1(2), BRR_ex2(1), BRR_ex1(1))
                    else
                        factors(class_j, sym_j)%elements(ind_j) = factors(class_j, sym_j)%elements(ind_j) + &
                          par_1*par_2*get_umat_el(BRR_ex2(2), BRR_ex1(2), BRR_ex2(1), BRR_ex1(1))
                    end if

                end do

            end do

            ! The factors array has now been fully generated for string_i. Now we just have to
            ! add in the contribution to vec_out from this string.
            ! This performs the final line in Eq. 20.

            ! Loop over all classes connected to the class of string_i.
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)

                ! The *required* symmetry,
                sym_j = ieor(HFSym_ras, sym_i)
                ! If there are no states in this class with the required symmetry.
                if (classes(class_j)%num_sym(sym_j) == 0) cycle

                ! Loop over all classes connected to string_j.
                do k = 1, classes(class_j)%num_comb
                    class_k = classes(class_j)%allowed_combns(k)
                    ! The *required* symmetry, sym_k = ieor(HFSym_ras, sym_j) = sym_i.
                    sym_k = sym_i

                    ! If there are no states in this class with the required symmetry.
                    if (classes(class_k)%num_sym(sym_k) == 0) cycle

                    ! Add in sigma_2.
                    vec_out(class_j, class_i, sym_j)%elements(:, ind_i) = &
                        vec_out(class_j, class_i, sym_j)%elements(:, ind_i) + &
                        matmul(vec_in(class_j, class_k, sym_j)%elements(:,:), factors(class_k, sym_k)%elements)

                    ! Add in sigma_1.
                    vec_out(class_i, class_j, sym_i)%elements(ind_i, :) = &
                        vec_out(class_i, class_j, sym_i)%elements(ind_i, :) + &
                        matmul(factors(class_k, sym_k)%elements, vec_in(class_k, class_j, sym_k)%elements(:,:))

                end do ! Over all classes connected to string_j.
                
            end do ! Over all classes connected to string_i.

        end do ! Over all strings.
        ! Next, calculate the contribution from the alpha-beta term (sigma_3) (Eq. 23).
        ! Loop over all combinations of two spatial indices, (kl).
        do k = 1, tot_norbs
            do l = 1, tot_norbs

                ! Zero arrays.
                r = 0
                ! The array c here is C' in Eq. 23.
                c = 0.0_dp

                ex1(1) = l
                ex1(2) = k
                ! Store these, for quick access later.
                BRR_ex1(1) = BRR(2*l)/2
                BRR_ex1(2) = BRR(2*k)/2

                ! Loop over all strings (equivalent to J_{beta} in Eq. 23).
                do i = 1, ras%num_strings

                    ! Get info for string_i.
                    sym_i = ras_strings(-1,i)
                    class_i = ras_strings(0,i)
                    nras1 = classes(class_i)%nelec_1
                    nras3 = classes(class_i)%nelec_3
                    string_i = ras_strings(1:tot_nelec,i)
                    ilut_i = ras_iluts(:,i)

                    ! This is the condition for (kl) to be an allowed excitation from the current string. 
                    if ( IsOcc(ilut_i,l) .and. &
                        (.not. (IsOcc(ilut_i,k) .and. k /= l)) ) then
                        ! Temporarily set these for get_excit_details to use.
                        sym_j = sym_i
                        class_j = class_i
                        string_j = string_i
                        call get_excit_details(ex1, ras, nras1, nras3, string_j, sym_j, class_j, in_ras_space)
                        ! If the excitation is to outside the RAS space, then we don't need to consider it.
                        if (in_ras_space) then
                            ! Store the class and symmetry of the excited string for later, to save some speed.
                            excit_info(1,i) = class_j
                            excit_info(2,i) = sym_j

                            ! The full address for string j.
                            ind_j = classes(class_j)%address_map(get_address(classes(class_j), string_j))
                            ! The shifted address, so that the first string in this class will have address 1.
                            ind_j = ind_j - classes(class_j)%cum_sym(sym_j)

                            r(i) = get_single_parity(ilut_i, l, k)

                            ! Construct array C' in Eq 23. To do so, loop over all connected strings.
                            do m = 1, classes(class_j)%num_comb
                                class_m = classes(class_j)%allowed_combns(m)
                                sym_m = ieor(HFSym_ras, sym_j)
                                if (classes(class_m)%num_sym(sym_m) /= 0) then
                                    min_ind = ras%cum_classes(class_m) + classes(class_m)%cum_sym(sym_m) + 1
                                    max_ind = min_ind + classes(class_m)%num_sym(sym_m) - 1
                                    c(min_ind:max_ind, i) = real(r(i),dp)*vec_in(class_j,class_m,sym_j)%elements(ind_j,:)
                                end if
                            end do
                         end if
                      endif

                end do

                ! Loop over all strings.
                do i = 1, ras%num_strings

                    ! This is equiavlent to F(J_{beta}) in Eq. 23 (except that we don't need to store
                    ! the value for all strings at once).
                    alpha_beta_fac = 0.0_dp

                    ! Get info for string_i.
                    sym_i = ras_strings(-1,i)
                    class_i = ras_strings(0,i)
                    ind_i = i - ras%cum_classes(class_i) - classes(class_i)%cum_sym(sym_i)

                    ! Loop over all single excitations from string_i (-> string_j).
                    do excit_j = 1, ras_excit(i)%nexcit

                        full_ind_j = ras_excit(i)%excit_ind(excit_j)
                        
                        class_j = ras_strings(0,full_ind_j)
                        par_1 = ras_excit(i)%par(excit_j)
                        ex2 = ras_excit(i)%orbs(:,excit_j)
                        BRR_ex2(1) = BRR(2*ex2(1))/2
                        BRR_ex2(2) = BRR(2*ex2(2))/2

                        alpha_beta_fac(full_ind_j) = alpha_beta_fac(full_ind_j) + &
                          par_1*get_umat_el(BRR_ex2(2), BRR_ex1(2), BRR_ex2(1), BRR_ex1(1))

                    end do

                    ! Loop over all classes allowed with class_i.
                    do j = 1, classes(class_i)%num_comb

                        class_j = classes(class_i)%allowed_combns(j)
                        sym_j = ieor(HFSym_ras, sym_i)

                        do ind_j = 1, classes(class_j)%num_sym(sym_j)

                            full_ind_j = ind_j + ras%cum_classes(class_j) + classes(class_j)%cum_sym(sym_j)

                            ! If this is true then (kl) isn't a valid excitation from string j.
                            if (r(full_ind_j) == 0) cycle

                            v = 0.0_dp
                            ! excited_class and excited_sym give the class and symmetry of the string
                            ! which we get by exciting string j with (kl).
                            excited_class = excit_info(1,full_ind_j)
                            excited_sym = excit_info(2,full_ind_j)

                            ! Loop over all classes connected to excited_class.
                            do m = 1, classes(excited_class)%num_comb
                                ! Get the class and required symmetry.
                                class_m = classes(excited_class)%allowed_combns(m)
                                sym_m = ieor(HFSym_ras, excited_sym)
                                if (classes(class_m)%num_sym(sym_m) /= 0) then
                                    min_ind = ras%cum_classes(class_m) + classes(class_m)%cum_sym(sym_m) + 1
                                    max_ind = min_ind + classes(class_m)%num_sym(sym_m) - 1
                                    v = v + dot_product(alpha_beta_fac(min_ind:max_ind), c(min_ind:max_ind, full_ind_j))
                                end if
                            end do

                            vec_out(class_j, class_i, sym_j)%elements(ind_j, ind_i) = &
                                vec_out(class_j, class_i, sym_j)%elements(ind_j, ind_i) + v

                        end do ! Over all strings in class_j with symmetry sym_j.

                    end do ! Over all classes allowed with class_i.

                end do ! Over all strings.

            end do ! Over all orbitals, l.
        end do ! Over all orbitals, k.
        ! Deallocate arrays.
        do class_i = 1, ras%num_classes
            do sym_i = 0, 7
               deallocate(factors(class_i,sym_i)%elements)
            end do
        end do
        deallocate(alpha_beta_fac)
        deallocate(c)
        deallocate(r)
    end subroutine perform_multiplication

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

    subroutine get_excit_details(ex, ras, nras1, nras3, string_j, sym_j, class_j, in_ras_space)

        integer, intent(in) :: ex(2)
        type(ras_parameters), intent(in) :: ras
        integer, intent(in) :: nras1, nras3
        integer, intent(inout) :: string_j(tot_nelec)
        integer, intent(inout) :: sym_j
        integer, intent(inout) :: class_j
        logical, intent(out) :: in_ras_space
        integer :: i, new1, new3
        integer :: sym_prod

        in_ras_space = .true.
        if (ex(1) == ex(2)) return

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

        if (.not. class_allowed(ras, new1, new3)) then
            in_ras_space = .false.
            return
        end if

        class_j = ras%class_label(new1, new3)

        do i = 1, tot_nelec
            if (string_j(i) == ex(1)) then
                string_j(i) = ex(2)
                exit
            end if
        end do
        call sort(string_j)

        sym_prod = int(ieor(G1(BRR(ex(1)*2))%Sym%S, G1(BRR(ex(2)*2))%Sym%S))
        sym_j = ieor(sym_j,sym_prod)

    end subroutine get_excit_details

    subroutine create_direct_ci_arrays(ras, classes, ras_strings, ras_iluts, ras_excit)

        ! This routine should be used before perform_multiplication is called. It will
        ! create some arrays which will make this multiplication quicker.

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        integer, intent(out) :: ras_strings(-1:tot_nelec, ras%num_strings)
        integer(n_int), intent(out) :: ras_iluts(0:NIfD, ras%num_strings)
        type(direct_ci_excit), intent(out) :: ras_excit(ras%num_strings)
        integer(n_int) :: ilut_i(0:NIfD)
        integer :: class_i, ind, new_ind, counter
        integer :: nras1, nras3, par, k
        integer :: ex(2)
        integer :: string_i(tot_nelec)
        logical :: none_left, tgen
        type(simple_excit_store), target :: gen_store_1

        ilut_i = 0_n_int

        ! Loop over all classes.
        do class_i = 1, ras%num_classes
            nras1 = classes(class_i)%nelec_1
            nras3 = classes(class_i)%nelec_3
            call generate_first_full_string(string_i, ras, classes(class_i))
            ! Loop over all strings.
            do
                ind = classes(class_i)%address_map(get_address(classes(class_i), string_i))
                do k = 1, class_i-1
                   ind = ind + classes(k)%class_size
                end do

                ! Store the string, along with the symmetry, RAS class and ilut.
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

                ! Store the number of excitations.
                ras_excit(ind)%nexcit = counter

                ! Now we know how many excitations there are, allocate the excitation array and store them.
                allocate(ras_excit(ind)%excit_ind(ras_excit(ind)%nexcit))
                allocate(ras_excit(ind)%par(ras_excit(ind)%nexcit))
                allocate(ras_excit(ind)%orbs(2,ras_excit(ind)%nexcit))

                counter = 0
                tgen = .true.
                do
                    call gen_next_single_ex(string_i, ilut_i, nras1, nras3, new_ind, par, &
                            ex, ras, classes, gen_store_1, tgen, .false.)
                    if (tgen) exit
                    counter = counter + 1
                    ! Store the index...
                    ras_excit(ind)%excit_ind(counter) = new_ind
                    ! ...and the parity of the excitation...
                    ras_excit(ind)%par(counter) = par
                    ! ...and the two orbitals involved in the excitation.
                    ras_excit(ind)%orbs(:,counter) = ex
                end do

                call generate_next_string(string_i, ras, classes(class_i), none_left)
                ! If no strings left in this class, then go to the next class.
                if (none_left) exit
            end do
        end do

    end subroutine create_direct_ci_arrays

    subroutine encode_string(string, ilut)

        ! Encode an alpha or beta string as an ilut.

        integer, intent(in) :: string(tot_nelec)
        integer(n_int), intent(out) :: ilut(0:NIfD)
        integer :: i, pos

        ilut = 0

        do i = 1, tot_nelec
            pos = (string(i)-1)/bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(string(i)-1, bits_n_int))
        end do
                
    end subroutine encode_string

    subroutine gen_next_single_ex(string_i, ilut_i, nras1, nras3, ind, par, ex, ras, classes, gen_store, tgen, tcount)

        ! Generate the next single excitation (within the RAS space - we don't need to consider
        ! excitations to outside it) and also return the associated parity and the class and
        ! address of the created string. If tcount is true then the routine is only being used
        ! to count the excitation, so these are not returned in this case.

        integer, intent(in) :: string_i(tot_nelec)
        integer(n_int), intent(in) :: ilut_i(0:NIfD)
        integer, intent(in) :: nras1, nras3
        integer, intent(out) :: ind, par
        integer, intent(out) :: ex(2)
        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        type(simple_excit_store), intent(inout), target :: gen_store
        logical, intent(inout) :: tgen
        logical, intent(in) :: tcount

        integer, pointer :: i, j
        integer :: orb1, orb2, temp1, temp3, class_k, m
        integer :: string_k(tot_nelec)

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

                ! Store the values of nras1 and nras3 for the new string.
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
                do m = 1, class_k-1
                   ind = ind  + classes(m)%class_size
                end do
                par = get_single_parity(ilut_i, ex(1), ex(2))

                return
            end do
            j = 0 ! Reset loop.
        end do

        tgen = .true.

    end subroutine gen_next_single_ex

    subroutine transfer_to_block_form(ras, classes, full_vec, ras_vec)

        ! Take a standard vector (as a 1d array) and transform it to a RAS vector, stored
        ! in matrix form. This can then be input into the perform_multiplication routine.

        ! Important: Because orbitals are ordered differently when using alpha and beta
        ! strings compared to the rest of NECI, some components of the vector will differ
        ! by a factor or -1. However, this routine does *not* apply these minus signs,
        ! so one should not take a vector from somewhere in NECI, use this routine and
        ! then the perform_multiplication routine. The extra minuses should be applied
        ! separately (although there is no routine to do this currently...)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        real(dp), intent(in) :: full_vec(:)
        type(ras_vector), intent(inout) :: ras_vec(ras%num_classes, ras%num_classes, 0:7)
        integer :: class_i, class_j, j, sym_i, sym_j, ind_i, ind_j
        integer :: counter

        counter = 0

        ! Over all RAS classes.
        do class_i = 1, ras%num_classes
            ! Over all classes allowed with class_i.
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)
                ! Over all symmetry labels.
                do sym_i = 0, 7
                    ! The symmetry label of string_j is fixed by that of string_i.
                    sym_j = ieor(HFSym_ras, sym_i)
                    ! If there are no strings with the correct symmetries in these classes.
                    if (classes(class_i)%num_sym(sym_i) == 0) cycle
                    if (classes(class_j)%num_sym(sym_j) == 0) cycle
                    ! Over all strings with symmetry sym_i and class class_i.
                    do ind_i = 1, classes(class_i)%num_sym(sym_i)
                        ! Over all strings with symmetry sym_j and class class_j.
                        do ind_j = 1, classes(class_j)%num_sym(sym_j)
                            counter = counter + 1
                            ! Copy the amplitude across.
                            ras_vec(class_i, class_j, sym_i)%elements(ind_i, ind_j) = full_vec(counter)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine transfer_to_block_form

    subroutine transfer_from_block_form(ras, classes, full_vec, ras_vec)

        ! See comments for transfer_to_block_form. This routine does the opposite transformation.

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        real(dp), intent(out) :: full_vec(:)
        type(ras_vector), intent(inout) :: ras_vec(ras%num_classes, ras%num_classes, 0:7)
        integer :: class_i, class_j, j, sym_i, sym_j, ind_i, ind_j
        integer :: counter

        counter = 0

        do class_i = 1, ras%num_classes
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)
                do sym_i = 0, 7
                    sym_j = ieor(HFSym_ras, sym_i)
                    if (classes(class_i)%num_sym(sym_i) == 0) cycle
                    if (classes(class_j)%num_sym(sym_j) == 0) cycle
                    do ind_i = 1, classes(class_i)%num_sym(sym_i)
                        do ind_j = 1, classes(class_j)%num_sym(sym_j)
                            counter = counter + 1
                            full_vec(counter) = ras_vec(class_i, class_j, sym_i)%elements(ind_i, ind_j)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine transfer_from_block_form

    subroutine create_ham_diag_direct_ci(ras, classes, ras_strings, ham_diag)

        ! The Davidson routine needs the Hamiltonian diagonal for the preconditioner.
        ! This routine creates and stores this array in the RAS space being used.

        use Determinants, only: get_helement
        use SystemData, only: nel

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(ras%num_classes)
        integer, intent(in) :: ras_strings(-1:tot_nelec, ras%num_strings)
        real(dp), intent(inout) :: ham_diag(:)
        integer :: class_i_ind, class_j_ind, sym_i_ind, sym_j_ind
        integer :: class_i, class_j, j, sym_i, sym_j
        integer :: ind_i, ind_j, full_ind_i, full_ind_j
        integer :: counter, k
        integer :: string_i(tot_nelec), string_j(tot_nelec), nI(nel)

        counter = 0

        ! Over all classes.
        do class_i = 1, ras%num_classes
            class_i_ind = ras%cum_classes(class_i)
            ! Over all classes connected to class_i.
            do j = 1, classes(class_i)%num_comb
                class_j = classes(class_i)%allowed_combns(j)
                class_j_ind = ras%cum_classes(class_j)
                ! Over all symmetry labels.
                do sym_i = 0, 7
                    ! The symmetry label for string_j is fixed by the label for string_i.
                    sym_j = ieor(HFSym_ras, sym_i)
                    ! The index of the first string with this symmetry and class.
                    sym_i_ind = class_i_ind + classes(class_i)%cum_sym(sym_i)
                    sym_j_ind = class_j_ind + classes(class_j)%cum_sym(sym_j)
                    ! If there are no strings with this symmetry and class.
                    if (classes(class_i)%num_sym(sym_i) == 0) cycle
                    if (classes(class_j)%num_sym(sym_j) == 0) cycle
                    ! Over all strings with this symmetry and class, for string_i.
                    do ind_i = 1, classes(class_i)%num_sym(sym_i)
                        full_ind_i = sym_i_ind + ind_i
                        ! Lookup the string corresponding to this address.
                        string_i = ras_strings(1:tot_nelec,full_ind_i)
                        ! Over all strings with this symmetry and class, for string_j.
                        do ind_j = 1, classes(class_j)%num_sym(sym_j)
                            counter = counter + 1
                            full_ind_j = sym_j_ind + ind_j
                            string_j = ras_strings(1:tot_nelec,full_ind_j)
                            ! Beta string.
                            nI(1:tot_nelec) = string_i*2-1
                            ! Alpha string.
                            nI(tot_nelec+1:nel) = string_j*2
                            ! Replace all orbital numbers, orb, with the true orbital
                            ! numbers, BRR(orb). Also, sort this list.
                            do k = 1, nel
                                nI(k) = BRR(nI(k))
                            end do
                            call sort(nI)
                            ! Finally calculate and store the corresponding element.
                            ham_diag(counter) = get_helement(nI, nI, 0)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine create_ham_diag_direct_ci

end module direct_ci
