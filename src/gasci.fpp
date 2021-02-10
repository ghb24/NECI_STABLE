#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci
    use SystemData, only: nBasis
    use util_mod, only: cumsum, stop_all, operator(.div.)
    use excitation_types, only: SingleExc_t, DoubleExc_t
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==)
    use util_mod, only: EnumBase_t
    implicit none

    private
    public :: possible_GAS_exc_gen, &
        GAS_exc_gen, GAS_specification, GASSpec_t, &
        user_input_GAS_exc_gen, get_name, construct_GASSpec_t

    type, extends(EnumBase_t) :: GAS_exc_gen_t
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
        integer, allocatable :: min(:), max(:)
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
        procedure :: n_spin_orbs => get_nOrbs
        procedure :: max_GAS_size => get_max_GAS_size
        procedure :: GAS_size => get_GAS_size
        procedure :: get_iGAS
        procedure :: get_orb_idx
        procedure :: split_per_GAS
        procedure :: count_per_GAS
        procedure :: write_to
        generic :: is_allowed => is_allowed_single
        generic :: is_allowed => is_allowed_double
        procedure, private :: is_allowed_single
        procedure, private :: is_allowed_double
        generic :: get_min => get_min_i
        generic :: get_min => get_min_all
        procedure, private :: get_min_i
        procedure, private :: get_min_all
        generic :: get_max => get_max_i
        generic :: get_max => get_max_all
        procedure, private :: get_max_i
        procedure, private :: get_max_all
    end type

    interface GASSpec_t
        module procedure construct_GASSpec_t
    end interface

    type(GASSpec_t), allocatable :: GAS_specification

contains

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
    !> Returns the minimum particle number for a given GAS space.
    integer elemental function get_min_i(self, iGAS)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: iGAS
        get_min_i = self%min(iGAS)
    end function

    !> @brief
    !> Returns the minimum particle number for a given GAS space.
    pure function get_min_all(self) result(res)
        class(GASSpec_t), intent(in) :: self
        integer, allocatable :: res(:)
        res = self%min(:)
    end function

    !> @brief
    !> Returns the maximum particle number for a given GAS space.
    integer elemental function get_max_i(self, iGAS)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in) :: iGAS
        get_max_i = self%max(iGAS)
    end function

    !> @brief
    !> Returns the maximum particle number for a given GAS space.
    pure function get_max_all(self) result(res)
        class(GASSpec_t), intent(in) :: self
        integer, allocatable :: res(:)
        res = self%max(:)
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
            integer :: all_orbs(n_spin_orbs), splitted_sizes(nGAS)
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
        res = all(self%min(:) <= supergroup .and. supergroup <= self%max(:))
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
    logical pure function is_valid(self, n_basis)
        class(GASSpec_t), intent(in) :: self
        integer, intent(in), optional :: n_basis

        logical :: shapes_match, pauli_principle, n_orbs_correct

        associate(GAS_sizes => self%GAS_sizes, n_min => self%min, &
                  n_max => self%max, nGAS => self%nGAS())

            shapes_match = &
                all([size(GAS_sizes), size(n_min), size(n_max)] == nGAS) &
                .and. maxval(self%GAS_table) == nGAS

            pauli_principle = all(n_min(:) <= GAS_sizes)

            if (present(n_basis)) then
                n_orbs_correct = sum(GAS_sizes) == n_basis
            else
                n_orbs_correct = .true.
            end if
        end associate

        is_valid = all([shapes_match, pauli_principle, n_orbs_correct])
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
        write(iunit, '(A)') 'n_min_i: minimum number of particles per i-th GAS space'
        write(iunit, '(A)') 'n_max_i: maximum number of particles per i-th GAS space'
        write(iunit, '(A10, 1x, A10, 1x, A10)') 'n_i', 'cn_min_i', 'cn_max_i'
        write(iunit, '(A)') '--------------------------------'
        do iGAS = 1, self%nGAS()
            write(iunit, '(I10, 1x, I10, 1x, I10)') self%GAS_size(iGAS) .div. 2, self%get_min(iGAS), self%get_max(iGAS)
        end do
        write(iunit, '(A)') '--------------------------------'
        write(iunit, '(A)') 'The distribution of spatial orbitals to GAS spaces is given by:'
        do iorb = 1, self%n_spin_orbs(), 2
            write(iunit, '(I0, 1x)', advance='no') self%get_iGAS(iorb)
        end do
        write(iunit, *)

        if (any(self%GAS_sizes < self%max)) then
            write(iunit, '(A)') 'In at least one GAS space, the maximum allowed particle number by GAS constraints'
            write(iunit, '(A)') '   is larger than the particle number allowed by the Pauli principle.'
            write(iunit, '(A)') '   Was this intended when preparing your input?'
        end if
    end subroutine

    !> @brief
    !> Check if a single excitation is allowed.
    !>
    !> @details
    !> Is called once at initialization, so it does not have to be super fast.
    logical pure function is_allowed_single(this, exc, supergroup)
        class(GASSpec_t), intent(in) :: this
        type(SingleExc_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)

        integer :: excited_supergroup(size(supergroup))
        integer :: src_space, tgt_space

        src_space = this%get_iGAS(exc%val(1))
        tgt_space = this%get_iGAS(exc%val(2))

        if (src_space == tgt_space) then
            ! All electrons come from the same space and there are no restrictions
            ! regarding recoupling or GAS.
            is_allowed_single = .true.
        else
            ! Ensure that GAS specifications contain supergroup **after** excitation.
            excited_supergroup = supergroup
            excited_supergroup(src_space) = excited_supergroup(src_space) - 1
            excited_supergroup(tgt_space) = excited_supergroup(tgt_space) + 1

            is_allowed_single = this%contains_supergroup(excited_supergroup)
        end if
    end function


    !> @brief
    !> Check if a double excitation is allowed.
    !>
    !> @details
    !> Is called once at initialization, so it does not have to be super fast.
    !> `recoupling` allows recoupling excitations that change the spin projection
    !> of individual GAS spaces.
    logical pure function is_allowed_double(this, exc, supergroup, recoupling)
        class(GASSpec_t), intent(in) :: this
        type(DoubleExc_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)
        logical, intent(in) :: recoupling

        integer :: excited_supergroup(size(supergroup))
        integer :: src_spaces(2), tgt_spaces(2)

        src_spaces = this%get_iGAS(exc%val(1, :))
        tgt_spaces = this%get_iGAS(exc%val(2, :))

        if (all(src_spaces == tgt_spaces) .and. src_spaces(1) == src_spaces(2)) then
            ! All electrons come from the same space and there are no restrictions
            ! regarding recoupling or GAS.
            is_allowed_double = .true.
        else
            ! Ensure that GAS specifications contain supergroup **after** excitation.
            excited_supergroup = supergroup
            excited_supergroup(src_spaces) = excited_supergroup(src_spaces) - 1
            excited_supergroup(tgt_spaces) = excited_supergroup(tgt_spaces) + 1

            is_allowed_double = this%contains_supergroup(excited_supergroup)

            if (is_allowed_double .and. .not. recoupling) then
                block
                    type(SpinProj_t) :: src_spins(2), tgt_spins(2)
                    #:set spin_swap = functools.partial(swap, 'SpinProj_t', "", 0)

                    src_spins = calc_spin_raw(exc%val(1, :))
                    tgt_spins = calc_spin_raw(exc%val(2, :))

                    if (src_spaces(1) > src_spaces(2)) then
                        @:spin_swap(src_spins(1), src_spins(2))
                    end if

                    if (tgt_spaces(1) > tgt_spaces(2)) then
                        @:spin_swap(tgt_spins(1), tgt_spins(2))
                    end if

                    is_allowed_double = all(src_spins == tgt_spins)
                end block
            end if
        end if
    end function
end module gasci
