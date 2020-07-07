#include "macros.h"
#:include "macros.fpph"

module gasci
    use SystemData, only: nBasis
    use util_mod, only: cumsum
    implicit none

    private
    public :: operator(==), operator(/=), possible_GAS_exc_gen, &
        GAS_exc_gen, GAS_specification, GASSpec_t, user_input_GAS_exc_gen

    type :: GAS_exc_gen_t
        integer :: val
    end type

    type :: possible_GAS_exc_gen_t
        type(GAS_exc_gen_t) :: &
            DISCONNECTED = GAS_exc_gen_t(1), &
            GENERAL = GAS_exc_gen_t(2)
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

    !> Speficies the GAS spaces.
    type :: GASSpec_t
        !> The indices are:
        !>  cn_min(1 : nGAS), cn_max(1 : nGAS), GAS_table(1 : nBasis)
        !> cn_min(iGAS) specifies the **cumulated** minimum particle number per GAS space.
        !> cn_max(iGAS) specifies the **cumulated** maximum particle number per GAS space.
        !> GAS_table(i) returns the GAS space for the i-th spin orbital
        integer, allocatable :: cn_min(:), cn_max(:), GAS_table(:)
        !> The number of GAS spaces
        integer :: nGAS
        !> The number of spin orbitals per GAS space
        integer, allocatable  :: GAS_sizes(:)
        !> maxval(GAS_sizes)
        integer :: max_GAS_size
        !> splitted_orbitals orbitals is the preimage of GAS_specification%GAS_table.
        !> An array that contains the spin orbitals per GAS space.
        !> splitted_orbitals(1 : maxval(GAS_sizes), 1 : nGAS)
        !> only splitted_orbitals(i, j), 1 <= i <= GAS_sizes(j)
        !> is defined.
        integer, allocatable :: splitted_orbitals(:, :)
    contains
        procedure :: contains => contains_det
        procedure :: is_connected
        procedure :: is_valid
        procedure :: split_per_GAS
        procedure :: count_per_GAS
    end type

    interface GASSpec_t
        module procedure construct_GASSpec_t
    end interface


    type(GASSpec_t) :: GAS_specification


contains

    logical pure function eq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        eq_GAS_exc_gen_t = lhs%val == rhs%val
    end function

    logical pure function neq_GAS_exc_gen_t(lhs, rhs)
        type(GAS_exc_gen_t), intent(in) :: lhs, rhs
        neq_GAS_exc_gen_t = lhs%val /= rhs%val
    end function


    !>  @brief
    !>      Constructor of GASSpec_t
    !>
    !>  @details
    !>
    !>  @param[in] n_min, Cumulative minimum particle number.
    !>  @param[in] n_max, Cumulative maximum particle number
    !>  @param[in] spat_GAS_orbs, GAS space for the i-th **spatial** orbital.
    function construct_GASSpec_t(n_min, n_max, spat_GAS_orbs) result(GAS_spec)
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
            @:ASSERT(all(GAS_sizes == splitted_sizes))
        end block

        GAS_spec = GASSpec_t(&
                n_min, n_max, GAS_table, nGAS, &
                GAS_sizes, max_GAS_size, splitted_orbitals)
        @:ASSERT(GAS_spec%is_valid())

        contains

        DEBUG_IMPURE function frequency(N) result(res)
            integer, intent(in) :: N(:)
            integer, allocatable :: res(:)
            integer :: i
            @:ASSERT(minval(N) == 1, N)
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
    logical pure function is_connected(GAS_spec)
        class(GASSpec_t), intent(in) :: GAS_spec
        is_connected = any(GAS_spec%cn_min(:) /= GAS_spec%cn_max(:))
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
    function contains_det(GAS_spec, occupied) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: occupied(:)

        logical :: res

        !> Cumulated number of particles per iGAS
        integer :: cum_n_particle(GAS_spec%nGAS), i

        cum_n_particle = cumsum(GAS_spec%count_per_GAS(occupied))

        res = all(GAS_spec%cn_min(:) <= cum_n_particle(:) &
            .and. cum_n_particle(:) <= GAS_spec%cn_max(:))
    end function


    logical pure function is_valid(GAS_spec, n_particles, n_basis)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in), optional :: n_particles, n_basis

        logical :: shapes_match, nEl_correct, pauli_principle, monotonic, &
            n_orbs_correct
        integer :: nGAS, iGAS, i

        associate(GAS_sizes => GAS_spec%GAS_sizes, n_min => GAS_spec%cn_min, &
                  n_max => GAS_spec%cn_max)

            nGAS = GAS_spec%nGAS

            shapes_match = &
                all([size(GAS_sizes), size(n_min), size(n_max)] == nGAS) &
                .and. maxval(GAS_spec%GAS_table) == nGAS

            if (present(n_particles)) then
                nEl_correct = all([n_min(nGAS), n_max(nGAS)] == n_particles)
            else
                nEl_correct = n_min(nGAS) == n_max(nGAS)
            end if

            pauli_principle = all(n_min(:) <= cumsum(GAS_sizes))

            if (nGAS >= 2) then
                monotonic = all([all(n_max(2:) >= n_max(: nGAS - 1)), &
                                 all(n_min(2:) >= n_min(: nGAS - 1))])
            else
                monotonic = .true.
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


    subroutine split_per_GAS(GAS_spec, occupied, splitted, splitted_sizes)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: occupied(:)
        integer, intent(out) :: &
            splitted(GAS_spec%max_GAS_size, GAS_spec%nGAS), &
            splitted_sizes(GAS_spec%nGAS)

        integer :: iel, iGAS

        splitted_sizes = 0
        do iel = 1, size(occupied)
            iGAS = GAS_spec%GAS_table(occupied(iel))
            splitted_sizes(iGAS) = splitted_sizes(iGAS) + 1
            splitted(splitted_sizes(iGAS), iGAS) = occupied(iel)
        end do
    end subroutine

    function count_per_GAS(GAS_spec, occupied) result(splitted_sizes)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: occupied(:)

        integer :: splitted_sizes(GAS_spec%nGAS)

        integer :: iel, iGAS

        splitted_sizes = 0
        do iel = 1, size(occupied)
            iGAS = GAS_spec%GAS_table(occupied(iel))
            splitted_sizes(iGAS) = splitted_sizes(iGAS) + 1
        end do
    end function
end module gasci
