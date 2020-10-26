#include "macros.h"
#:include "macros.fpph"

module gasci
    use SystemData, only: nBasis
    use util_mod, only: cumsum, stop_all, operator(.div.)
    implicit none

    private
    public :: operator(==), operator(/=), possible_GAS_exc_gen, &
        GAS_exc_gen, GAS_specification, GASSpec_t, &
        user_input_GAS_exc_gen, get_name, construct_GASSpec_t



    type :: GAS_exc_gen_t
        integer :: val
    end type

    type :: possible_GAS_exc_gen_t
        type(GAS_exc_gen_t) :: &
            DISCONNECTED = GAS_exc_gen_t(1), &
            GENERAL = GAS_exc_gen_t(2), &
            DISCARDING = GAS_exc_gen_t(3), &
            GENERAL_PCHB = GAS_exc_gen_t(4)
    end type

    type(possible_GAS_exc_gen_t), parameter :: possible_GAS_exc_gen = possible_GAS_exc_gen_t()

    type(GAS_exc_gen_t) :: GAS_exc_gen = possible_GAS_exc_gen%GENERAL
    type(GAS_exc_gen_t), allocatable :: user_input_GAS_exc_gen

    interface operator(==)
        module procedure eq_GAS_exc_gen_t
    end interface

    interface operator(/=)
        module procedure neq_GAS_exc_gen_t
    end interface

    ! NOTE: At the current state of implementation `GASSpec_t` is a completely immutable
    ! datastructure to outside code after the constructor has been called.
    ! If possible keep it like this!

    !> Speficies the GAS spaces.
    type :: GASSpec_t
        !> The indices are:
        !>  cn_min(1 : nGAS), cn_max(1 : nGAS)
        !> cn_min(iGAS) specifies the **cumulated** minimum particle number per GAS space.
        !> cn_max(iGAS) specifies the **cumulated** maximum particle number per GAS space.
        private
        integer, allocatable :: cn_min(:), cn_max(:)
        !> GAS_table(i) returns the GAS space for the i-th spin orbital
        integer, allocatable :: GAS_table(:)
        !> The number of spin orbitals per GAS space
        integer, allocatable  :: GAS_sizes(:)
        !> maxval(GAS_sizes)
        integer :: largest_GAS_size
        !> splitted_orbitals orbitals is the preimage of GAS_specification%GAS_table.
        !> An array that contains the spin orbitals per GAS space.
        !> splitted_orbitals(1 : maxval(GAS_sizes), 1 : nGAS)
        !> only splitted_orbitals(i, j), 1 <= i <= GAS_sizes(j)
        !> is defined.
        integer, allocatable :: splitted_orbitals(:, :)
        !> These lookup variables stay valid, because the data structure is
        !>  immutable
        logical :: lookup_is_connected
    contains
        ! All member functions should be public.
        procedure :: contains_det
        procedure :: contains_supergroup
        procedure :: is_connected => get_is_connected
        procedure :: is_valid
        procedure :: nGAS => get_nGAS
        procedure :: nEl => get_nEl
        procedure :: n_spin_orbs => get_nOrbs
        procedure :: max_GAS_size => get_max_GAS_size
        procedure :: GAS_size => get_GAS_size
        procedure :: get_iGAS
        procedure :: get_orb_idx
        procedure :: cumulated_min
        procedure :: cumulated_max
        procedure :: split_per_GAS
        procedure :: count_per_GAS
        procedure :: write_to
    end type

    interface GASSpec_t
        module procedure construct_GASSpec_t
    end interface

    type(GASSpec_t), allocatable :: GAS_specification

contains

    logical elemental function eq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        eq_GAS_exc_gen_t = lhs%val == rhs%val
    end function

    logical elemental function neq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        neq_GAS_exc_gen_t = lhs%val /= rhs%val
    end function

    !> @brief
    !> Returns the total number of GAS spaces.
    integer pure function get_nGAS(self)
        class(GASSpec_t), intent(in) :: self
        get_nGAS = size(self%GAS_sizes)
    end function

    !> @brief
    !> Returns the size of the largest GAS space.
    integer pure function get_max_GAS_size(self)
        class(GASSpec_t), intent(in) :: self
        get_max_GAS_size = self%largest_GAS_size
    end function

    !> @brief
    !>  Returns the size of the i-th GAS space.
    integer elemental function get_GAS_size(self, iGAS)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: iGAS
        get_GAS_size = self%GAS_sizes(iGAS)
    end function

    !> @brief
    !> Returns the GAS space for a given spin orbital index.
    integer elemental function get_iGAS(self, spin_orb_idx)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: spin_orb_idx
        get_iGAS = self%GAS_table(spin_orb_idx)
    end function

    !> @brief
    !> Returns the cumulated minimum particle number for a given GAS space.
    integer elemental function cumulated_min(self, iGAS)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: iGAS
        cumulated_min = self%cn_min(iGAS)
    end function

    !> @brief
    !> Returns the cumulated maximum particle number for a given GAS space.
    integer elemental function cumulated_max(self, iGAS)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: iGAS
        cumulated_max = self%cn_max(iGAS)
    end function

    integer elemental function get_nEl(self)
        class(GASSpec_t), intent(in) :: self
        get_nEl = self%cn_min(size(self%cn_min))
    end function

    integer elemental function get_nOrbs(self)
        class(GASSpec_t), intent(in) :: self
        get_nOrbs = size(self%GAS_table)
    end function

    !> @brief
    !> Returns the i-th spin orbital in the iGAS GAS space.
    !>
    !> @details
    !>  Can be seen as the preimage of get_iGAS (which is usually not injective).
    integer elemental function get_orb_idx(self, i, iGAS)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: i, iGAS
        character(*), parameter :: this_routine = 'get_orb_idx'

        @:pure_ASSERT(1 <= i .and. i <= self%GAS_size(iGAS))
        @:pure_ASSERT(1 <= iGAS .and. iGAS <= self%nGAS())

        get_orb_idx = self%splitted_orbitals(i, iGAS)
    end function


    !>  @brief
    !>      Constructor of GASSpec_t
    !>
    !>  @details
    !>
    !>  @param[in] n_min, Cumulative minimum particle number.
    !>  @param[in] n_max, Cumulative maximum particle number
    !>  @param[in] spat_GAS_orbs, GAS space for the i-th **spatial** orbital.
    pure function construct_GASSpec_t(n_min, n_max, spat_GAS_orbs) result(GAS_spec)
        integer, intent(in) :: n_min(:), n_max(:)
        integer, intent(in) :: spat_GAS_orbs(:)

        type(GASSpec_t) :: GAS_spec
        character(*), parameter :: this_routine = 'construct_GASSpec_t'

        integer :: n_spin_orbs, max_GAS_size
        integer, allocatable :: splitted_orbitals(:, :), GAS_table(:), GAS_sizes(:)
        integer :: i, iel, iGAS, nGAS

        nGAS = maxval(spat_GAS_orbs)
        GAS_sizes = 2 * frequency(spat_GAS_orbs)

        max_GAS_size = maxval(GAS_sizes)
        n_spin_orbs = sum(GAS_sizes)

        allocate(GAS_table(n_spin_orbs))
        GAS_table(1::2) = spat_GAS_orbs
        GAS_table(2::2) = spat_GAS_orbs

        allocate(splitted_orbitals(max_GAS_size, nGAS))
        block
            integer :: counter(nGAS), all_orbs(n_spin_orbs)
            integer :: splitted_sizes(nGAS)
            all_orbs = [(i, i = 1, n_spin_orbs)]
            splitted_sizes = 0
            do iel = 1, size(all_orbs)
                iGAS = GAS_table(iel)
                splitted_sizes(iGAS) = splitted_sizes(iGAS) + 1
                splitted_orbitals(splitted_sizes(iGAS), iGAS) = all_orbs(iel)
            end do
            @:pure_ASSERT(all(GAS_sizes == splitted_sizes))
        end block

        GAS_spec = GASSpec_t(&
                n_min, n_max, GAS_table, &
                GAS_sizes, max_GAS_size, splitted_orbitals, &
                any(n_min(:) /= n_max(:)))

        @:pure_ASSERT(GAS_spec%is_valid())

        contains

        pure function frequency(N) result(res)
            integer, intent(in) :: N(:)
            integer, allocatable :: res(:)
            integer :: i
            @:pure_ASSERT(minval(N) == 1)
            allocate(res(maxval(N)), source=0)
            do i = 1, size(N)
                res(N(i)) = res(N(i)) + 1
            end do
        end function
    end function


    !>  @brief
    !>      Query if there are connected GAS spaces under the GAS specification.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    logical pure function get_is_connected(self)
        class(GASSpec_t), intent(in) :: self
        get_is_connected = self%lookup_is_connected
    end function


    !>  @brief
    !>      Query wether a determinant is contained in the GAS space.
    !>
    !>  @details
    !>  It is **assumed** that the determinant is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  The return value is not defined, if that is not the case!
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] occupied, An index of occupied spin orbitals.
    pure function contains_det(self, occupied) result(res)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: occupied(:)

        logical :: res

        res = self%contains_supergroup(self%count_per_GAS(occupied))
    end function

    !>  @brief
    !>      Query wether a supergroup is contained in the GAS space.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] supergroup, A supergroup.
    pure function contains_supergroup(self, supergroup) result(res)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: supergroup(:)

        logical :: res

        !> Cumulated number of particles per iGAS
        integer :: cn_particle(size(supergroup))

        cn_particle = cumsum(supergroup)

        res = all(self%cn_min(:) <= cn_particle(:) &
            .and. cn_particle(:) <= self%cn_max(:))
    end function


    !>  @brief
    !>      Check if the GAS specification is valid
    !>
    !>  @details
    !>   If the number of particles or the number of spin orbitals
    !>   is provided, then the consistency with these numbers
    !>   is checked as well.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] n_particles, Optional.
    !>  @param[in] n_basis, Optional. The number of spin orbitals.
    logical pure function is_valid(self, n_particles, n_basis)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in), optional :: n_particles, n_basis

        logical :: shapes_match, nEl_correct, pauli_principle, monotonic, &
            n_orbs_correct
        integer :: iGAS, i

        associate(GAS_sizes => self%GAS_sizes, n_min => self%cn_min, &
                  n_max => self%cn_max, nGAS => self%nGAS())

            shapes_match = &
                all([size(GAS_sizes), size(n_min), size(n_max)] == nGAS) &
                .and. maxval(self%GAS_table) == nGAS

            pauli_principle = all(n_min(:) <= cumsum(GAS_sizes))

            if (nGAS >= 2) then
                monotonic = all([all(n_max(2:) >= n_max(: nGAS - 1)), &
                                 all(n_min(2:) >= n_min(: nGAS - 1))])
            else
                monotonic = .true.
            end if

            if (present(n_particles)) then
                nEl_correct = all([n_min(nGAS), n_max(nGAS)] == n_particles)
            else
                nEl_correct = n_min(nGAS) == n_max(nGAS)
            end if

            if (present(n_basis)) then
                n_orbs_correct = sum(GAS_sizes) == n_basis
            else
                n_orbs_correct = .true.
            end if
        end associate

        is_valid = all([shapes_match, nEl_correct, pauli_principle, &
                        monotonic, n_orbs_correct])
    end function

    pure subroutine split_per_GAS(self, occupied, splitted, splitted_sizes)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: occupied(:)
        integer, intent(out) :: &
            splitted(get_max_GAS_size(self), get_nGAS(self)), &
            splitted_sizes(get_nGAS(self))

        integer :: iel, iGAS

        splitted_sizes = 0
        do iel = 1, size(occupied)
            iGAS = self%get_iGAS(occupied(iel))
            splitted_sizes(iGAS) = splitted_sizes(iGAS) + 1
            splitted(splitted_sizes(iGAS), iGAS) = occupied(iel)
        end do
    end subroutine

    pure function count_per_GAS(self, occupied) result(splitted_sizes)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: occupied(:)

        integer :: splitted_sizes(get_nGAS(self))

        integer :: iel, iGAS

        splitted_sizes = 0
        do iel = 1, size(occupied)
            iGAS = self%get_iGAS(occupied(iel))
            splitted_sizes(iGAS) = splitted_sizes(iGAS) + 1
        end do
    end function

    pure function get_name(impl) result(res)
        type(GAS_exc_gen_t), intent(in) :: impl
        character(len=:), allocatable :: res
        if (impl == possible_GAS_exc_gen%DISCONNECTED) then
            res = 'Heat-bath on-the-fly GAS implementation for disconnected spaces'
        else if (impl == possible_GAS_exc_gen%GENERAL) then
            res = 'Heat-bath on-the-fly GAS implementation'
        else if (impl == possible_GAS_exc_gen%DISCARDING) then
            res = 'Discarding GAS implementation'
        else if (impl == possible_GAS_exc_gen%GENERAL_PCHB) then
            res = 'PCHB GAS implementation'
        end if
    end function

    !> @brief
    !> Write a string representation of this GAS specification to iunit
    subroutine write_to(self, iunit)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: iunit
        integer :: iGAS, iorb

        write(iunit, '(A)') 'n_i: number of spatial orbitals per i-th GAS space'
        write(iunit, '(A)') 'cn_min_i: cumulative minimum number of particles per i-th GAS space'
        write(iunit, '(A)') 'cn_max_i: cumulative maximum number of particles per i-th GAS space'
        write(iunit, '(A10, 1x, A10, 1x, A10)') 'n_i', 'cn_min_i', 'cn_max_i'
        write(iunit, '(A)') '--------------------------------'
        do iGAS = 1, self%nGAS()
            write(iunit, '(I10, 1x, I10, 1x, I10)') self%GAS_size(iGAS) .div. 2, self%cumulated_min(iGAS), self%cumulated_max(iGAS)
        end do
        write(iunit, '(A)') '--------------------------------'
        write(iunit, '(A)') 'The distribution of spatial orbitals to GAS spaces is given by:'
        do iorb = 1, self%n_spin_orbs(), 2
            write(iunit, '(I0, 1x)', advance='no') self%get_iGAS(iorb)
        end do
        write(iunit, *)
    end subroutine

end module gasci
