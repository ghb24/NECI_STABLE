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
    use FciMCData, only: HFSym
    use ras_data
    use sort_mod, only: sort
    use sym_mod, only: getsym
    use SystemData, only: G1, nbasismax, nel, nbasis, basisfn, BRR, tHub
    use util_mod, only: find_next_comb

    implicit none

contains

    subroutine initialise_ras_space(ras, classes)

        ! Given the five defining parameters for a RAS space, initialise the ras
        ! space by finding the ras classes corresponding to the ras space.

        type(ras_parameters), intent(inout) :: ras
        type(ras_class_data), intent(inout), allocatable, dimension(:) :: classes
        integer :: i, j, k, counter
        integer :: lower_ras3, upper_ras3

        tot_nelec = nel/2
        tot_norbs = nbasis/2

        ! Check that the RAS parameters are possible.
        if (ras%size_1+ras%size_2+ras%size_3 /= tot_norbs .or. &
            ras%min_1 > ras%size_1*2 .or. ras%max_3 > ras%size_3*2) &
            call stop_all("generate_ras", "RAS parameters are not possible.")
        if (mod(nel, 2) /= 0) call stop_all("generate_ras", "RAS-core only implmented for &
                                            & closed shell molecules.")

        ! First we need to find the different classes. A class is defined by the number
        ! of electrons in RAS1 and RAS3. Thus, we need to find all possible allowed
        ! combinations.

        ras%lower_ras1 = max(0, ras%min_1-ras%size_1)
        ras%upper_ras1 = min(tot_nelec, ras%size_1)

        allocate(ras%class_label(ras%lower_ras1:ras%upper_ras1, 0:tot_nelec))
        ras%class_label = 0

        ! First count the number of RAS classes...
        counter = 0
        do i = ras%lower_ras1, ras%upper_ras1
            lower_ras3 = max(0, tot_nelec-i-ras%size_2)
            upper_ras3 = min(tot_nelec-i, ras%max_3)
            do j = lower_ras3, upper_ras3
                counter = counter + 1
                ras%class_label(i,j) = counter
            end do
        end do

        ras%num_classes = counter
        counter = 0

        ! ...then fill the classes in.
        allocate(classes(ras%num_classes))
        do i = ras%lower_ras1, ras%upper_ras1
            lower_ras3 = max(0, tot_nelec-i-2*ras%size_2)
            upper_ras3 = min(tot_nelec-i, ras%max_3)
            do j = lower_ras3, upper_ras3
                counter = counter + 1
                classes(counter)%nelec_1 = i
                classes(counter)%nelec_3 = j
                classes(counter)%nelec_2 = tot_nelec-i-j
                if (classes(counter)%nelec_2 < 0) then
                    call stop_all("initialise_ras_space", &
                                  "Current RAS combination is not possible.")
                end if
                allocate(classes(counter)%vertex_weights(0:tot_norbs, 0:tot_nelec))
                classes(counter)%vertex_weights = 0
            end do
        end do

        ! Form the vertex weights for each class.
        do i = 1, ras%num_classes
            classes(i)%vertex_weights(0, 0) = 1
            do j = 1, tot_norbs
                do k = 0, tot_nelec
                    ! If vertex not allowed then leave the corresponding vertex weight as 0.
                    if (vertex_not_allowed(classes(i)%nelec_1, &
                                           classes(i)%nelec_3, j, k, ras)) cycle
                    ! (Eq. 11.8.2)
                    if (k == 0) then
                        ! The first columns will always contain 1's.
                        classes(i)%vertex_weights(j,k) = 1
                    else
                        classes(i)%vertex_weights(j,k) = classes(i)%vertex_weights(j-1,k) + &
                                                              classes(i)%vertex_weights(j-1,k-1)
                    end if
                end do
            end do
        end do

        ! Find the allowed combinations of classes in the full RAS space.
        do i = 1, ras%num_classes
            counter = 0
            allocate(classes(i)%allowed_combns(ras%num_classes))
            classes(i)%allowed_combns = 0
            do j = 1, ras%num_classes
                ! If the total number of electrons in the RAS spaces are correct with
                ! this combination.
                if (class_comb_allowed(ras, classes(i), classes(j))) then
                    counter = counter + 1
                    classes(i)%allowed_combns(counter) = j
                end if
            end do
            classes(i)%num_comb = counter
        end do

        allocate(ras%cum_classes(ras%num_classes))
        ras%cum_classes(1) = 0

        ras%num_strings = 0
        do i = 1, ras%num_classes
            call setup_ras_class(ras, classes(i))
            ras%num_strings = ras%num_strings + classes(i)%class_size
            if (i > 1) ras%cum_classes(i) = ras%cum_classes(i-1) + classes(i-1)%class_size
        end do

        HFSym_ras = int(HFSym%Sym%S)

    end subroutine initialise_ras_space

    pure function class_allowed(ras, n_elec_1, n_elec_3) result (allowed)

        ! This function assumes that the total number of electrons is equal to tot_nelec, so
        ! that we don't have to check the number of electrons in RAS2. It also assumes
        ! obvious things like the numbers of electrons not being negative.

        type(ras_parameters), intent(in) :: ras
        integer, intent(in) :: n_elec_1, n_elec_3
        integer :: lower_ras3, upper_ras3
        logical :: allowed

        allowed = .false.

        if (n_elec_1 >= ras%lower_ras1 .and. n_elec_1 <= ras%upper_ras1) then
            lower_ras3 = max(0, tot_nelec-n_elec_1-ras%size_2)
            upper_ras3 = min(tot_nelec-n_elec_1, ras%max_3)
            if (n_elec_3 >= lower_ras3 .and. n_elec_3 <= upper_ras3) then
                allowed = .true.
            end if
        end if

    end function class_allowed

    pure function class_comb_allowed(ras, class_1, class_2) result (allowed)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: class_1, class_2
        logical :: allowed

        allowed = (class_1%nelec_1+class_2%nelec_1 >= ras%min_1 .and. &
                    class_1%nelec_3+class_2%nelec_3 <= ras%max_3)

    end function class_comb_allowed

    pure function vertex_not_allowed(n_elec_1, n_elec_3, orb, elec, ras) result(not_allowed)

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

    pure function get_address(ras_class, string) result(address)

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

        ! Count the number of states with each symmetry label.
        ras_class%num_sym = 0
        do i = 1, ras_class%class_size
            ras_class%num_sym(symmetries(i)) = ras_class%num_sym(symmetries(i)) + 1
        end do

        do i = 0, 7
            ras_class%cum_sym(i) = sum(ras_class%num_sym(0:i-1))
        end do

        deallocate(symmetries)
        deallocate(addresses)

    end subroutine setup_ras_class

    pure subroutine generate_first_subspace_string(string, n_elec)

        ! Generate the first string (lowest orbitals all occupied) in a RAS subspace.
        ! For RAS2 and RAS3 the orbital numbers should have been shifted so that the
        ! first orbital in these subspaces is 1, *not* the actual orbital number.

        integer, intent(in) :: n_elec
        integer, intent(out) :: string(n_elec)
        integer :: i

        do i = 1, n_elec
            string(i) = i
        end do

    end subroutine generate_first_subspace_string

    pure subroutine generate_first_full_string(string, ras, ras_class)

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

    pure subroutine generate_next_string(string, ras, ras_class, none_left)

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

    subroutine mysort(vec, parity)
          integer, dimension(:), intent(inout) :: vec
          integer, intent(out) :: parity 
          integer :: temp, bubble, lsup, j
          lsup = size(vec)
          parity = 1
          do while (lsup > 1)
            bubble = 0 !bubble in the greatest element out of order
            do j = 1, (lsup-1)
              if (vec(j) > vec(j+1)) then
                temp = vec(j)
                vec(j) = vec(j+1)
                vec(j+1) = temp
                bubble = j
                parity=-1*parity
              endif 
            enddo
            lsup = bubble   
          enddo  
    end subroutine
    subroutine generate_entire_ras_space(ras, classes, space_size, ilut_list, parities)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(:)
        integer, intent(in) :: space_size
        integer(n_int), intent(out) :: ilut_list(0:NIfTot, space_size)
        integer, intent(out), optional :: parities(space_size)
        integer :: nI(nel)
        integer :: string(tot_nelec)
        integer, allocatable, dimension(:,:) :: string_list
        integer, allocatable, dimension(:,:) :: min_indices
        integer :: temp_class
        integer :: string_address, block_address
        integer :: i, j, k, l, m, n, o, counter
        logical :: none_left
        integer :: parity

        allocate(string_list(tot_nelec, ras%num_strings))

        ! min_indices(i,j) stores the index of the first state in class i with symmetry label j.
        allocate(min_indices(ras%num_classes, 0:7))

        ! Loop over all classes, and for each, generate all states in that class. Store all these
        ! states in blocks, one after another, in a 1d array (sting_list). Within each block
        ! find_address and classes(i)%address_map are used to to find an address *relative to
        ! the first state* in that class.
        block_address = 0
        do i = 1, ras%num_classes
            ! Address of the last string in the last class.
            if (i > 1) block_address = block_address + classes(i-1)%class_size

            do j = 0, 7
                min_indices(i,j) = block_address + sum(classes(i)%num_sym(0:j-1)) + 1
            end do

            call generate_first_full_string(string, ras, classes(i))
            string_address = classes(i)%address_map(get_address(classes(i), string)) + &
                             block_address
            string_list(:, string_address) = string

            do j = 2, classes(i)%class_size
                call generate_next_string(string, ras, classes(i), none_left)
                string_address = classes(i)%address_map(get_address(classes(i), string)) + &
                                 block_address
                string_list(:, string_address) = string
            end do
        end do

        ! Now combine the alpha and beta strings together to create the full list of states.

        counter = 0
        ilut_list = 0
        ! Loop over all classes, call it class_i.
        do i = 1, ras%num_classes
            ! For class_i, loop over all classes which can be combined with this one.
            do j = 1, classes(i)%num_comb
                temp_class = classes(i)%allowed_combns(j)
                ! Loop over all symmetries.
                do k = 0, 7
                    l = ieor(HFSym_ras, k)
                    ! Add all combinations of states in these symmetry blocks to ilut_list.
                    do m = min_indices(i,k), min_indices(i,k)+&
                            classes(i)%num_sym(k)-1
                        do n = min_indices(temp_class,l), min_indices(temp_class,l)+&
                                classes(temp_class)%num_sym(l)-1

                            ! Beta string.
                            nI(1:tot_nelec) = string_list(:, m)*2-1
                            !nI(1:nel-1:2) = string_list(:, m)*2-1
                            ! Alpha string.
                            nI(tot_nelec+1:nel) = string_list(:, n)*2
                            !nI(2:nel:2) = string_list(:, n)*2

                            ! Replace all orbital numbers, orb, with the true orbital
                            ! numbers, BRR(orb). Also, sort this list.
                            do o = 1, nel
                                nI(o) = BRR(nI(o))
                            end do
                            !call sort(nI, par=parity)
                            call mysort(nI, parity)

                            ! Find bitstring representation.
                            counter = counter+1
                            call EncodeBitDet(int(nI,sizeof_int), ilut_list(:,counter))
                            if(present(parities))then
                                parities(counter) = parity
                            end if
                        end do
                    end do
                end do
            end do
        end do

        if (counter /= space_size) call stop_all("generate_entire_ras_space", "Wrong number of &
                                                 &states found in generating loop.")

        deallocate(string_list)
        deallocate(min_indices)

    end subroutine generate_entire_ras_space

    pure function get_abelian_sym(string) result(sym)

        ! Note, the orbital numbers in the input string refer to the spatial orbitals, so we
        ! mnultiply these by 2 when used in G1.

        integer, intent(in) :: string(:)
        integer :: sym
        integer(int64) :: temp_sym
        integer :: i

        if(tHub) then
            !Since RAS is originally developed for molucules, it cannot handle kpoint symmetries.
            !As a quick fix, we ignore symmetry labels of the Hubbard model.
            sym = 0
        else
            temp_sym = G1(BRR(string(1)*2))%Sym%S

            do i = 2, size(string)
                temp_sym = ieor(temp_sym, G1(BRR(string(i)*2))%Sym%S)
            end do

            sym = int(temp_sym, sizeof_int)
        end if

    end function get_abelian_sym

    pure subroutine find_ras_size(ras, classes, space_size)

        type(ras_parameters), intent(in) :: ras
        type(ras_class_data), intent(in) :: classes(:)
        integer, intent(out) :: space_size
        integer :: i, j, k, l
        integer :: temp_class

        space_size = 0

        ! Loop over all classes, call it class_i.
        do i = 1, ras%num_classes
            ! For this class, loop over all classes which can be combined with this one,
            ! call it class_j.
            do j = 1, classes(i)%num_comb
                temp_class = classes(i)%allowed_combns(j)
                ! Loop over all symmetries for class_i.
                do k = 0, 7
                    ! Required symmetry for class_j.
                    l = ieor(HFSym_ras, k)
                    ! Finally, add the total number of combinations of strings from the two
                    ! classes with these symmetry labels.
                    space_size = space_size + &
                        int(classes(i)%num_sym(k)*classes(temp_class)%num_sym(l),sizeof_int)
                end do
            end do
        end do

    end subroutine find_ras_size

end module ras
