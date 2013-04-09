#include "macros.h"

! A RAS space is defined by splitting all spatial orbitals into
! three subspaces, RAS1, RAS2 and RAS3. A minimum number of
! electrons (both alpha and beta) is specified to occupy the RAS1
! orbitals and a maximum number of electrons to occupy the RAS3
! orbitals. There are no restrictions imposed on the occupation
! of RAS2. Examples of RAS spaces include the spaces of all
! excitations up to any level (e.g. all singles and doubles), all
! RAS spaces, and all single and double excitations from a RAS space.

! To use a RAS space, create a variable of type ras_parameters and
! an allocatable 1d array of type ras_class_data. Define the five
! parameters in the ras_parameters variable. Then input these two
! variables into the subroutine initialise_ras_space. This will allocate
! and fill in the ras_class_data array.

! Note 1: All orbitals numbers in this module refer to energy order of
! an orbital, rather than the actual orbital number. Hence for an orbital,
! orb, in this module, use BRR(orb) to get the true orbital number.

! Note 2: Currently it is assumed that we are in an Ms=0 subspace so that
! the set of alpha strings is the same as the set of beta strings.

module ras

    use bit_reps, only: NIfTot
    use constants
    use DetBitOps, only: EncodeBitDet
    use FciMCData, only: HFDet
    use sort_mod, only: sort
    use sym_mod, only: getsym
    use SystemData, only: G1, nbasismax, nel, basisfn, BRR
    use util_mod, only: find_next_comb

    implicit none

    ! The five parameters defining a RAS space.
    type ras_parameters
        ! The number of (spatial) orbitals in each of the RAS spaces.
        integer :: size_1, size_2, size_3
        ! The minimum number of electrons (both alpha and beta) in RAS1,
        ! and the maximum number in RAS3.
        integer :: min_1, max_3
    end type

    ! A RAS class refers to a collection of strings with an fixed number
    ! of electrons in RAS1, RAS2 and RAS3. Thus, the collection of all
    ! RAS strings will, in general, be formed from many RAS classes.
    type ras_class_data
        ! The number of electrons in each of RAS1, RAS2 and RAS3.
        integer :: nelec_1, nelec_2, nelec_3
        ! The total number of strings in this class.
        integer :: class_size
        ! The vertex weights used in the addressing scheme.
        integer, allocatable, dimension(:,:) :: vertex_weights
        ! The number of classes which can be combined with this one in
        ! a full determinant to give the correct overall RAS parameters.
        integer :: num_comb
        ! The labels of the classes which can be combined with this one.
        integer, allocatable, dimension(:) :: allowed_combns
        ! A one-to-one map from the address of a string, as obtained by the
        ! function get_address, to one where states are sorted by their symmetry.
        integer, allocatable, dimension(:) :: address_map ! (class_size)
        ! The number of strings in this class with each of the symmetry labels.
        integer :: num_sym(0:7)
        ! The number of symmetry labels which have atleast one string in this
        ! class (i.e. the number of non-zero elements in num_sym).
        integer :: non_zero_sym
    end type

    type sym_combinations
        ! The number of symmetry values which can be combined with another
        ! symmetry to get the correct total symmetry for the full determinant.
        integer :: num_comb
        ! The symmetry labels which can be combined correctly. This has 8
        ! elements, but not all will be used.
        integer :: allowed_combns(8)
    end type

    type(sym_combinations) :: sym_combs(0:7)

    ! The number of electrons occupying one alpha or beta string. As only
    ! Ms=0 is implemented, this is just nOccAlpha=nOccBeta=nEl/2
    integer :: tot_nelec
    ! The number of spatial orbitals.
    integer :: tot_norbs
    ! The total symmetry of the Hartree-Fock state (always 0 for the only case
    ! that can be treated so far...)
    integer :: HFSym

contains

    subroutine initialise_ras_space(ras, ras_classes)

        ! Given the five defining parameters for a RAS space, initialise the ras
        ! space by finding the ras classes corresponding to the ras space.

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(inout), allocatable, dimension(:) :: ras_classes
        integer :: i, j, k, counter
        integer :: num_ras_classes, class_size
        type(basisfn) :: hfbasisfn
        integer :: string(tot_nelec)
        integer :: space_size

        ! First we need to find the different classes. A class is defined by the number
        ! of electrons in RAS1 and RAS3. Thus, we need to find all possible allowed
        ! combinations.

        ! First count the number of RAS classes...
        counter = 0
        do i = max(0, ras%min_1-ras%size_1), min(tot_nelec, ras%size_1)
            do j = max(0, tot_nelec-i-2*ras%size_2) , min(tot_nelec-i, ras%max_3)
                counter = counter + 1
            end do
        end do

        num_ras_classes = counter
        counter = 0

        ! ...then fill the classes in.
        allocate(ras_classes(num_ras_classes))
        do i = max(0, ras%min_1-ras%size_1), min(tot_nelec, ras%size_1)
            do j = max(0, tot_nelec-i-ras%size_2) , min(tot_nelec-i, ras%max_3)
                counter = counter + 1
                ras_classes(counter)%nelec_1 = i
                ras_classes(counter)%nelec_3 = j
                ras_classes(counter)%nelec_2 = tot_nelec-i-j
                if (ras_classes(counter)%nelec_2 < 0) then
                    write(6,*) "nelec_1, nelec_2, nelec_3:", ras_classes(counter)%nelec_1, &
                        ras_classes(counter)%nelec_2, ras_classes(counter)%nelec_3
                    call stop_all("initialise_ras_space", &
                                  "Current RAS combination is not possible.")
                end if
                allocate(ras_classes(counter)%vertex_weights(0:tot_norbs, 0:tot_nelec))
                ras_classes(counter)%vertex_weights = 0
            end do
        end do

        ! Form the vertex weights for each class.
        do i = 1, num_ras_classes
            ras_classes(i)%vertex_weights(0, 0) = 1
            do j = 1, tot_norbs
                do k = 0, tot_nelec
                    ! If vertex not allowed then leave the corresponding vertex weight as 0.
                    if (vertex_not_allowed(ras_classes(i)%nelec_1, &
                                           ras_classes(i)%nelec_3, j, k, ras)) cycle
                    ! (Eq. 11.8.2)
                    if (k == 0) then
                        ! The first columns will always contain 1's.
                        ras_classes(i)%vertex_weights(j,k) = 1
                    else
                        ras_classes(i)%vertex_weights(j,k) = ras_classes(i)%vertex_weights(j-1,k) + &
                                                              ras_classes(i)%vertex_weights(j-1,k-1)
                    end if
                end do
            end do
        end do

        do i = 1, num_ras_classes
            call print_matrix(ras_classes(i)%vertex_weights, tot_norbs, tot_nelec)
        end do

        ! Find the allowed combinations of classes in the full RAS space.
        do i = 1, num_ras_classes
            counter = 0
            allocate(ras_classes(i)%allowed_combns(num_ras_classes))
            ras_classes(i)%allowed_combns = 0
            do j = 1, num_ras_classes
                ! If the total number of electrons in the RAS spaces are correct with
                ! this combination.
                if (ras_classes(i)%nelec_1+ras_classes(j)%nelec_1 >= ras%min_1 .and. &
                    ras_classes(i)%nelec_3+ras_classes(j)%nelec_3 <= ras%max_3) then
                    counter = counter + 1
                    ras_classes(i)%allowed_combns(counter) = j
                end if
            end do
            ras_classes(i)%num_comb = counter
        end do

        do i = 1, num_ras_classes
            call setup_ras_class(ras, ras_classes(i))
        end do

        call getsym(HFDet, nel, G1, nbasismax, hfbasisfn)
        HFSym = hfbasisfn%Sym%S

    end subroutine initialise_ras_space

    function vertex_not_allowed(n_elec_1, n_elec_3, orb, elec, ras) result(not_allowed)

        integer, intent(in) :: n_elec_1, n_elec_3
        integer, intent(in) :: orb, elec
        type(ras_parameters), intent(in) :: ras
        integer :: n_elec_2
        logical :: not_allowed

        not_allowed = .true.

        n_elec_2 = tot_nelec - n_elec_1 - n_elec_3

        ! For the current (orb, elec) vertex to be allowed, it must lie within one of three
        ! trapeziums, defined by the RAS parameters. The first corresponds to RAS1 orbitals, 
        ! the second to RAS2 orbitals and the third to RAS3 orbitals. The three if statements
        ! below correspond to these three cases, and the if statements within to the
        ! trapeziums which they must lie in.

        if (orb <= ras%size_1) then
            ! If a RAS1 orbital.
            ! Condition for (orb,elec) combination to be allowed.
            if (orb-elec <= ras%size_1-n_elec_1 .and. orb-elec >= 0 &
                .and. elec <= n_elec_1) not_allowed = .false.

        else if (orb > ras%size_1+ras%size_2) then
            ! If a RAS3 orbital.
            if (orb+(tot_nelec-elec) <= tot_norbs .and. &
                orb+(tot_nelec-elec) >= tot_norbs-(ras%size_3-n_elec_3) .and. &
                elec >= (n_elec_1+n_elec_2)) not_allowed = .false.
        
        else
            ! If a RAS2 orbital.
            if (orb-(elec-n_elec_1) <= ras%size_1+(ras%size_2-n_elec_2) .and. &
                orb-(elec-n_elec_1) >= ras%size_1 .and. &
                elec >= n_elec_1 .and. elec <= (n_elec_1+n_elec_2)) not_allowed = .false.

        end if

    end function vertex_not_allowed

    function get_address(ras_class, string) result(address)

        type(ras_class_data), intent(in) :: ras_class
        integer, intent(in) :: string(:)
        integer :: elec
        integer :: address

        ! (Eq. 11.8.3)
        address = 1
        do elec = 1, size(string)
            address = address + ras_class%vertex_weights(string(elec), elec) - &
                                ras_class%vertex_weights(string(elec)-1, elec-1)
        end do

    end function get_address

    subroutine setup_ras_class(ras, ras_class)

        ! Create a one-to-one map from the original address of a string, as found by
        ! the function get_address, to a new address where states are sorted in order of
        ! their symmetry. Also count the number of strings with each symmetry label.

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(inout) :: ras_class
        integer :: string(ras_class%nelec_1+ras_class%nelec_2+ras_class%nelec_3)
        integer, allocatable, dimension(:) :: symmetries, addresses
        integer :: i, j, counter
        logical :: none_left

        ! Find the number of strings in this class, the maximum address. To find this, put
        ! all electrons in their maximum possible orbtials.
        counter = 0
        do j = 1, ras_class%nelec_1
            counter = counter + 1
            string(counter) = ras%size_1 - ras_class%nelec_1 + j
        end do
        do j = 1, ras_class%nelec_2
            counter = counter + 1
            string(counter) = ras%size_1+ras%size_2 - ras_class%nelec_2 + j
        end do
        do j = 1, ras_class%nelec_3
            counter = counter + 1
            string(counter) = ras%size_1+ras%size_2+ras%size_3 - ras_class%nelec_3 + j
        end do

        ras_class%class_size = get_address(ras_class, string)

        ! Now that we have the class size.
        allocate(symmetries(ras_class%class_size))
        allocate(ras_class%address_map(ras_class%class_size))
        allocate(addresses(ras_class%class_size))

        ! Now generate the string with the lowest address.
        call generate_first_full_string(string, ras, ras_class)
        symmetries(1) = get_abelian_sym(string)
        addresses(1) = 1
        
        ! Loop over all possible strings in this class, and store the symmetry of each.
        do i = 2, ras_class%class_size
            call generate_next_string(string, ras, ras_class, none_left)

            addresses(i) = get_address(ras_class, string)

            symmetries(i) = get_abelian_sym(string)
        end do

        ! Check that the last string in the above loop really was the last string (as it should
        ! be) by calling the generating routine again.
        call generate_next_string(string, ras, ras_class, none_left)
        if (.not. none_left) call stop_all("find_symmetries", "Incorrect number of states found.")

        ! Now, sort the symmetries array from smallest to largest.
        call sort(symmetries, addresses)

        ! The addresses array was sorted with symmetries. If we have the final position of a state
        ! *after* the sort and we put that position into addresses, we get out the address of this
        ! state (as found by get_address). We want the inverse of this, a map that will take
        ! in an address from get_address and give out the position after the sort.
        do i = 1, ras_class%class_size
            ! g(f(x)) = x => g=f^(-1), so the following creates the inverse map.
            ras_class%address_map(addresses(i)) = i
        end do

        ! Finally, count the number of states with each symmetry label.
        ras_class%num_sym = 0
        do i = 1, ras_class%class_size
            ras_class%num_sym(symmetries(i)) = ras_class%num_sym(symmetries(i)) + 1
        end do

        deallocate(symmetries)
        deallocate(addresses)

    end subroutine setup_ras_class

    subroutine generate_first_subspace_string(string, n_elec)

        ! Generate the first string (lowest orbitals all occupied) in a RAS subspace.
        ! For RAS2 and RAS3 the orbital numbers should have been shifted so that the
        ! first orbital in these subspaces is 1, *not* the actualt orbital number.

        integer, intent(in) :: n_elec
        integer, intent(out) :: string(n_elec)
        integer :: i

        do i = 1, n_elec
            string(i) = i
        end do

    end subroutine generate_first_subspace_string

    subroutine generate_first_full_string(string, ras, ras_class)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: ras_class
        integer, intent(inout) :: string(ras_class%nelec_1+ras_class%nelec_2+&
                                         ras_class%nelec_3)
        integer :: i, counter

        ! In each space subspace (RAS1, RAS2, RAS3) put each electron in the lowest orbitals.
        counter = 1
        do i = 1, ras_class%nelec_1
            string(i) = counter
            counter = counter + 1
        end do
        counter = 1
        do i = ras_class%nelec_1+1, ras_class%nelec_1+ras_class%nelec_2
            string(i) = ras%size_1 + counter
            counter = counter + 1
        end do
        counter = 1
        do i = ras_class%nelec_1+ras_class%nelec_2+1, tot_nelec
            string(i) = ras%size_1+ras%size_2 + counter
            counter = counter + 1
        end do

    end subroutine generate_first_full_string

    subroutine generate_next_string(string, ras, ras_class, none_left)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: ras_class
        integer, intent(inout) :: string(ras_class%nelec_1+ras_class%nelec_2+&
                                         ras_class%nelec_3)
        logical, intent(out) :: none_left
        integer :: string_1(ras_class%nelec_1)
        integer :: string_2(ras_class%nelec_2)
        integer :: string_3(ras_class%nelec_3)

        ! The strings in the 3 RAS spaces, shifted so that the first orbital in each space would
        ! be labelled 1 (this is needed for the routine that generates all combinations).
        string_1 = string(1:ras_class%nelec_1)
        string_2 = string(ras_class%nelec_1+1:ras_class%nelec_1+ras_class%nelec_2) - ras%size_1
        string_3 = string(ras_class%nelec_1+ras_class%nelec_2+1:tot_nelec) - &
                   (ras%size_1+ras%size_2)

        loop1: do ! Over RAS1.

            loop2: do ! Over RAS2.

                loop3: do ! Over RAS3.

                    call find_next_comb(string_3, ras_class%nelec_3, ras%size_3, none_left)
                    if (.not. none_left) exit loop1
                    ! If we had the last string last time, go back to the first one.
                    call generate_first_subspace_string(string_3, ras_class%nelec_3)
                    exit loop3

                end do loop3

                call find_next_comb(string_2, ras_class%nelec_2, ras%size_2, none_left)
                if (.not. none_left) exit loop1
                call generate_first_subspace_string(string_2, ras_class%nelec_2)
                exit loop2

            end do loop2

            call find_next_comb(string_1, ras_class%nelec_1, ras%size_1, none_left)
            exit loop1

        end do loop1

        ! If the last string in RAS1, RAS2 and RAS3 was found last time, then last_time will be
        ! equal to .true. at this point and this value will be returned, signalling the end.

        string(1:ras_class%nelec_1) = string_1
        string(ras_class%nelec_1+1:ras_class%nelec_1+ras_class%nelec_2) = string_2 + ras%size_1
        string(ras_class%nelec_1+ras_class%nelec_2+1:tot_nelec) = string_3 + (ras%size_1+ras%size_2)

    end subroutine generate_next_string

    subroutine generate_entire_ras_space(ras, ras_classes, space_size, ilut_list)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: ras_classes(:)
        integer, intent(in) :: space_size
        integer(n_int), intent(out) :: ilut_list(0:NIfTot, space_size)
        integer :: nI(nel)
        integer :: string(tot_nelec)
        integer, allocatable, dimension(:,:) :: string_list
        integer, allocatable, dimension(:,:) :: min_indices
        integer :: num_ras_classes, num_strings, temp_class
        integer :: string_address, block_address
        integer :: i, j, k, l, m, n, o, counter
        logical :: none_left

        ! Calculate the total number of alpha/beta strings.
        num_ras_classes = size(ras_classes)
        num_strings = 0
        do i = 1, num_ras_classes
            num_strings = num_strings + ras_classes(i)%class_size
        end do

        allocate(string_list(tot_nelec, num_strings))

        ! min_indices(i,j) stores the index of the first state in class i with symmetry label j.
        allocate(min_indices(num_ras_classes, 0:7))

        ! Loop over all classes, and for each, generate all states in that class. Store all these
        ! states in blocks, one after another, in a 1d array (sting_list). Within each block
        ! find_address and ras_classes(i)%address_map are used to to find an address *relative to
        ! the first state* in that class.
        block_address = 0
        do i = 1, num_ras_classes
            ! Address of the last string in the last class.
            if (i > 1) block_address = block_address + ras_classes(i-1)%class_size

            do j = 0, 7
                min_indices(i,j) = block_address + sum(ras_classes(i)%num_sym(0:j-1)) + 1
            end do

            call generate_first_full_string(string, ras, ras_classes(i))
            string_address = ras_classes(i)%address_map(get_address(ras_classes(i), string)) + &
                             block_address
            string_list(:, string_address) = string

            do j = 2, ras_classes(i)%class_size
                call generate_next_string(string, ras, ras_classes(i), none_left)
                string_address = ras_classes(i)%address_map(get_address(ras_classes(i), string)) + &
                                 block_address
                string_list(:, string_address) = string
            end do
        end do

        ! Now combine the alpha and beta strings together to create the full list of states.

        counter = 0
        ilut_list = 0
        ! Loop over all classes, call it class_i.
        do i = 1, num_ras_classes
            ! For this class, loop over all classes which can be combined with this one, call
            ! it class_j.
            do j = 1, ras_classes(i)%num_comb
                temp_class = ras_classes(i)%allowed_combns(j)
                ! Loop over all symmetries.
                do k = 0, 7
                    ! Loop over all symmetries.
                    do l = 0, 7
                        ! If these two symmetries give the correct total symmetry then add all
                        ! combinations of states in these symmetry blocks to ilut_list.
                        if (ieor(k, l) == HFSym) then
                            do m = min_indices(i,k), min_indices(i,k)+&
                                    ras_classes(i)%num_sym(k)-1
                                do n = min_indices(temp_class,l), min_indices(temp_class,l)+&
                                        ras_classes(temp_class)%num_sym(l)-1

                                    ! Beta string.
                                    nI(1:tot_nelec) = string_list(:, m)*2-1
                                    ! Alpha string.
                                    nI(tot_nelec+1:nel) = string_list(:, n)*2

                                    ! Replace all orbital numbers, orb, with the true orbital
                                    ! numbers, BRR(orb). Also, sort this list.
                                    do o = 1, nel
                                        nI(o) = BRR(nI(o))
                                    end do
                                    call sort(nI)

                                    ! Find bitstring representation.
                                    counter = counter+1
                                    call EncodeBitDet(nI, ilut_list(:,counter))

                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do

        if (counter /= space_size) call stop_all("generate_entire_ras_space", "Wrong number of &
                                                 &states found in generating loop.")

        deallocate(string_list)
        deallocate(min_indices)

    end subroutine generate_entire_ras_space

    function get_abelian_sym(string) result(sym)

        ! Note, the orbital numbers in the input string refer to the spatial orbitals, so we
        ! mnultiply these by 2 when used in G1.

        integer, intent(in) :: string(:)
        integer :: sym
        integer :: i

        sym = G1(BRR(string(1)*2))%Sym%S

        do i = 2, size(string)
            sym = ieor(sym, G1(BRR(string(i)*2))%Sym%S)
        end do

    end function get_abelian_sym

    subroutine find_ras_size(ras, ras_classes, space_size)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: ras_classes(:)
        integer, intent(out) :: space_size
        integer :: i, j, k, l
        integer :: num_ras_classes, temp_class

        space_size = 0
        num_ras_classes = size(ras_classes)

        ! Loop over all classes, call it class_i.
        do i = 1, num_ras_classes
            ! For this class, loop over all classes which can be combined with this one,
            ! call it class_j.
            do j = 1, ras_classes(i)%num_comb
                temp_class = ras_classes(i)%allowed_combns(j)
                ! Loop over all symmetries for class_i.
                do k = 0, 7
                    ! Loop over all symmetries for class_j.
                    do l = 0, 7
                        ! Finally, if these two symmetries combine to give the correct
                        ! overall symmetry then add the total number of combinations of
                        ! strings from the two classes with these symmetry labels.
                        if (ieor(k, l) == HFSym) space_size = space_size + &
                            ras_classes(i)%num_sym(k)*ras_classes(temp_class)%num_sym(l)
                    end do
                end do
            end do
        end do

    end subroutine find_ras_size

    subroutine find_sym_combinations()

        ! Fill in suym_combs variable. sym_combs(i) stores the symmetry labels which can
        ! be combined with label i to give the same overall symmetry as the HF state.

        integer :: i, j, counter

        do i = 0, 7
            counter = 0
            do j = 0, 7
                if (ieor(i,j) == HFSym) then
                    counter = counter + 1
                    sym_combs(i)%allowed_combns(counter) = j
                end if
            end do
            sym_combs(i)%num_comb = counter
        end do

    end subroutine find_sym_combinations

    subroutine print_matrix(matrix, N, M)

        integer :: I, J
        integer, intent(in) :: N, M
        integer, intent(in) :: matrix(0:N, 0:M)

        write(6,*)

        do I = 0, N
            do J = 0, M
                write(6, fmt = '(I4,A1)', advance = 'no') matrix(I,J), ' '
            end do
            print *,
        end do

        write(6,*)

    end subroutine print_matrix

end module ras
