#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']


module gasci
    use constants, only: n_int, dp
    use SystemData, only: nBasis
    use util_mod, only: cumsum, stop_all, operator(.div.)
    use excitation_types, only: SingleExc_t, DoubleExc_t
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==)
    use util_mod, only: lex_leq, cumsum, operator(.div.), near_zero, binary_search_first_ge, &
        operator(.isclose.), custom_findloc, EnumBase_t
    use sets_mod, only: disjoint, operator(.U.), is_sorted, operator(.complement.)
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==), operator(/=), operator(-), sum, &
        alpha, beta
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite, get_last_tgt, set_last_tgt, UNKNOWN
    use sltcnd_mod, only: sltcnd_excit
    use dSFMT_interface, only: genrand_real2_dSFMT
    use bit_rep_data, only: nIfTot
    use bit_reps, only: decode_bit_det
    use growing_buffers, only: buffer_int_2D_t, buffer_int_1D_t

    implicit none
    private
    public :: possible_GAS_exc_gen, &
        GAS_exc_gen, GAS_specification, GASSpec_t, &
        user_input_GAS_exc_gen, get_name, LocalGASSpec_t, CumulGASSpec_t


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

    type, abstract :: GASSpec_t
        !! Speficies the GAS spaces.
        private
        integer, allocatable :: GAS_table(:)
            !! GAS_table(i) returns the GAS space for the i-th spin orbital
        integer, allocatable  :: GAS_sizes(:)
            !! The number of spin orbitals per GAS space
        integer :: largest_GAS_size
            !! maxval(GAS_sizes)
        integer, allocatable :: splitted_orbitals(:, :)
            !! This is the preimage of `GASSpec_t%GAS_table.`
            !!
            !! An array that contains the spin orbitals per GAS space of dimension
            !! `splitted_orbitals(1 : maxval(GAS_sizes), 1 : nGAS)`.
            !! Only `splitted_orbitals(i, j), 1 <= i <= GAS_sizes(j)`
            !! is defined.
        logical :: lookup_is_connected
            !! These lookup variables stay valid, because the data structure is
            !!  immutable
        logical :: exchange_recoupling
    contains
        ! Nearly all member functions should be public.
        procedure(contains_supergroup_t), deferred :: contains_supergroup
        procedure(is_valid_t), deferred :: is_valid
        procedure(write_to_t), deferred :: write_to
        procedure(get_possible_spaces_t), deferred :: get_possible_spaces
        procedure :: get_possible_holes
        procedure :: contains_conf
        procedure :: contains_ilut
        procedure :: is_connected => get_is_connected
        procedure :: nGAS => get_nGAS
        procedure :: n_spin_orbs => get_nOrbs
        procedure :: max_GAS_size => get_max_GAS_size
        procedure :: recoupling

        generic :: GAS_size => get_GAS_size_i
        generic :: GAS_size => get_GAS_size_idx
        generic :: GAS_size => get_GAS_size_all
        procedure, private :: get_GAS_size_i
        procedure, private :: get_GAS_size_idx
        procedure, private :: get_GAS_size_all

        procedure :: get_iGAS
        procedure :: get_orb_idx
        procedure :: count_per_GAS
        generic :: is_allowed => is_allowed_single
        generic :: is_allowed => is_allowed_double
        procedure, private :: is_allowed_single
        procedure, private :: is_allowed_double
    end type

    abstract interface
        logical pure function contains_supergroup_t(this, supergroup)
            import :: GASSpec_t
            implicit none
            class(GASSpec_t), intent(in) :: this
            integer, intent(in) :: supergroup(:)
        end function

        logical pure function is_valid_t(this, n_basis)
            import :: GASSpec_t
            implicit none
            class(GASSpec_t), intent(in) :: this
            integer, intent(in), optional :: n_basis
        end function

        subroutine write_to_t(this, iunit)
            !! Write a string representation of this GAS specification to iunit
            import :: GASSpec_t
            implicit none
            class(GASSpec_t), intent(in) :: this
            integer, intent(in) :: iunit
        end subroutine

        pure function get_possible_spaces_t(this, supergroup, add_holes, add_particles, n_total) result(spaces)
            !!  Return the GAS spaces, where one particle can be created.
            !!
            !!  The returned array can be empty (allocated, but size == 0).
            import :: GASSpec_t
            implicit none
            class(GASSpec_t), intent(in) :: this
                !!  Specification of GAS spaces.
            integer, intent(in) :: supergroup(size(this%GAS_sizes))
                !!  The particles per GAS space.
            integer, intent(in), optional :: add_holes(:)
                !!  An index of orbitals where particles should be deleted
                !!  before creating the new particle.
            integer, intent(in), optional :: add_particles(:)
                !!  Index of orbitals where particles should be created
                !! before creating the new particle.
            integer, intent(in), optional :: n_total
                !! The total number of particles that will be created.
                !! Defaults to one.  (Relevant for double excitations)
            integer, allocatable :: spaces(:)
        end function
    end interface

    type, extends(GASSpec_t) :: LocalGASSpec_t
        private
        !> The indices are:
        !> `min(1 : nGAS), max(1 : nGAS)`.
        !> `min(iGAS)` specifies the minimum particle number per GAS space.
        !> `max(iGAS)` specifies the maximum particle number per GAS space.
        integer, allocatable :: min(:), max(:)
    contains
        procedure :: contains_supergroup => Local_contains_supergroup
        procedure :: is_valid => Local_is_valid
        procedure :: write_to => Local_write_to
        procedure :: get_possible_spaces => Local_get_possible_spaces

        generic :: get_min => get_min_i
        generic :: get_min => get_min_all
        procedure, private :: get_min_i
        procedure, private :: get_min_all
        generic :: get_max => get_max_i
        generic :: get_max => get_max_all
        procedure, private :: get_max_i
        procedure, private :: get_max_all
    end type

    interface LocalGASSpec_t
        module procedure construct_LocalGASSpec_t
    end interface

    type, extends(GASSpec_t) :: CumulGASSpec_t
        private
        !> The indices are:
        !> `c_min(1 : nGAS), c_max(1 : nGAS)`.
        !> `c_min(iGAS)` specifies the cumulated minimum particle number per GAS space.
        !> `c_max(iGAS)` specifies the cumulated maximum particle number per GAS space.
        integer, allocatable :: c_min(:), c_max(:)
    contains
        procedure :: contains_supergroup => Cumul_contains_supergroup
        procedure :: is_valid => Cumul_is_valid
        procedure :: write_to => Cumul_write_to
        procedure :: get_possible_spaces => Cumul_get_possible_spaces

        generic :: get_cmin => get_cmin_i
        generic :: get_cmin => get_cmin_all
        procedure, private :: get_cmin_i
        procedure, private :: get_cmin_all
        generic :: get_cmax => get_cmax_i
        generic :: get_cmax => get_cmax_all
        procedure, private :: get_cmax_i
        procedure, private :: get_cmax_all
    end type

    interface CumulGASSpec_t
        module procedure construct_CumulGASSpec_t
    end interface


    class(GASSpec_t), allocatable :: GAS_specification

contains

    !> Returns the total number of GAS spaces.
    integer pure function get_nGAS(this)
        class(GASSpec_t), intent(in) :: this
        get_nGAS = size(this%GAS_sizes)
    end function

    !> Returns the size of the largest GAS space.
    integer pure function get_max_GAS_size(this)
        class(GASSpec_t), intent(in) :: this
        get_max_GAS_size = this%largest_GAS_size
    end function

    !>  Returns the size of the i-th GAS space in number of spin orbitals.
    integer pure function get_GAS_size_i(this, iGAS)
        class(GASSpec_t), intent(in) :: this
        integer, intent(in) :: iGAS
        get_GAS_size_i = this%GAS_sizes(iGAS)
    end function

    !>  Returns the sizes for GAS spaces specified in idx.
    pure function get_GAS_size_idx(this, idx) result(res)
        class(GASSpec_t), intent(in) :: this
        integer, intent(in) :: idx(:)
        integer :: res(size(idx))
        res(:) = this%GAS_sizes(idx)
    end function

    !>  Returns the sizes for all GAS spaces.
    pure function get_GAS_size_all(this) result(res)
        class(GASSpec_t), intent(in) :: this
        integer :: res(size(this%GAS_sizes))
        res(:) = this%GAS_sizes(:)
    end function

    !> Returns the GAS space for a given spin orbital index.
    integer elemental function get_iGAS(this, spin_orb_idx)
        class(GASSpec_t), intent(in) :: this
        integer, intent(in) :: spin_orb_idx
        get_iGAS = this%GAS_table(spin_orb_idx)
    end function

    integer elemental function get_nOrbs(this)
        class(GASSpec_t), intent(in) :: this
        get_nOrbs = size(this%GAS_table)
    end function

    !> Returns the i-th spin orbital in the iGAS GAS space.
    !>
    !> Can be seen as the preimage of get_iGAS (which is usually not injective).
    integer elemental function get_orb_idx(this, i, iGAS)
        class(GASSpec_t), intent(in) :: this
        integer, intent(in) :: i, iGAS
        character(*), parameter :: this_routine = 'get_orb_idx'

        @:pure_ASSERT(1 <= i .and. i <= this%GAS_size(iGAS))
        @:pure_ASSERT(1 <= iGAS .and. iGAS <= this%nGAS())

        get_orb_idx = this%splitted_orbitals(i, iGAS)
    end function


    logical pure function get_is_connected(this)
        !! Query if there are connected GAS spaces under the GAS specification.
        class(GASSpec_t), intent(in) :: this
            !!  Specification of GAS spaces.
        get_is_connected = this%lookup_is_connected
    end function


    !>  Query wether a determinant or CSF is contained in the GAS space.
    !>
    !>  It is **assumed** that the configuration is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  The return value is not defined, if that is not the case!
    pure function contains_conf(this, nI) result(res)
        !> Specification of GAS spaces.
        class(GASSpec_t), intent(in) :: this
        !> An index of occupied spin orbitals.
        integer, intent(in) :: nI(:)
        logical :: res
        res = this%contains_supergroup(this%count_per_GAS(nI))
    end function

    logical elemental function recoupling(this)
        class(GASSpec_t), intent(in) :: this
        recoupling = this%exchange_recoupling
    end function

    !>  Query wether a determinant in bitmask format is contained in the GAS space.
    !>
    !>  The function in nI-format is faster!
    !>  It is **assumed** that the determinant is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  The return value is not defined, if that is not the case!
    pure function contains_ilut(this, ilut) result(res)
        !>  Specification of GAS spaces.
        class(GASSpec_t), intent(in) :: this
        !>  An index of occupied spin orbitals.
        integer(n_int), intent(in) :: ilut(0 : nIfTot)
        logical :: res
        integer :: nI(sum(popcnt(ilut)))
        call decode_bit_det(nI, ilut)
        res = this%contains_conf(nI)
    end function

    !> Count the particles per GAS space. i.e. return the supergroup.
    pure function count_per_GAS(this, occupied) result(supergroup)
        class(GASSpec_t), intent(in) :: this
        integer, intent(in) :: occupied(:)

        integer :: supergroup(get_nGAS(this))

        integer :: iel, iGAS

        supergroup = 0
        do iel = 1, size(occupied)
            iGAS = this%get_iGAS(occupied(iel))
            supergroup(iGAS) = supergroup(iGAS) + 1
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

    !> Check if a single excitation is allowed.
    !>
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


    !> Check if a double excitation is allowed.
    !>
    !> Is called once at initialization, so it does not have to be super fast.
    !> `recoupling` allows recoupling excitations that change the spin projection
    !> of individual GAS spaces.
    logical pure function is_allowed_double(this, exc, supergroup)
        class(GASSpec_t), intent(in) :: this
        type(DoubleExc_t), intent(in) :: exc
        integer, intent(in) :: supergroup(:)

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

            if (is_allowed_double .and. .not. this%exchange_recoupling) then
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

    !>  Return the possible holes where a particle can be created under GAS constraints.
    !>
    !>  This function uses `get_possible_spaces` to find possible GAS spaces
    !>  where a particle can be created and returns only unoccupied
    !>  sites of correct spin.
    !>
    !>  "Trivial" excitations are avoided. That means, that a site is only counted
    !>  as unoccupied if it was unoccupied in nI from the beginning on.
    !>  (A double excitation where a particle is deleted, but immediately
    !>  recreated would be such a trivial excitations.)
    pure function get_possible_holes(this, det_I, add_holes, add_particles, n_total, excess) result(possible_holes)
        !>  Specification of GAS spaces.
        class(GASSpec_t), intent(in) :: this
        !>  The starting determinant
        integer, intent(in) :: det_I(:)
        !>  An index of orbitals
        !>      where particles should be deleted before creating the new particle.
        integer, intent(in), optional :: add_holes(:)
        !>  An index of orbitals
        !>      where particles should be created before creating the new particle.
        integer, intent(in), optional :: add_particles(:)
        !>  The total number of particles
        !>      that will be created. Defaults to one.
        integer, intent(in), optional :: n_total
        !>  The current excess of spin projections.
        !>      If a beta electron was deleted, the excess is \(1 \cdot \alpha)\.
        type(SpinProj_t), intent(in), optional :: excess
        character(*), parameter :: this_routine = 'get_possible_holes'

        integer, allocatable :: possible_holes(:)

        integer, allocatable :: spaces(:)
        integer :: n_total_

        @:def_default(n_total_, n_total, 1)
        @:pure_ASSERT(1 <= n_total_)
        if (present(excess)) then
            @:pure_ASSERT(abs(excess%val) <= n_total_)
        end if


        ! Note, that non-present optional arguments can be passed
        ! into optional arguments without checking!
        spaces = this%get_possible_spaces(&
             this%count_per_GAS(det_I), add_holes=add_holes, &
             add_particles=add_particles, n_total=n_total_)

        if (size(spaces) == 0) then
            possible_holes = [integer::]
            return
        end if

        block
            integer :: i, iGAS, incr, curr_value, iGAS_min_val, idx_space
            integer :: L(size(spaces)), counter(size(spaces))
            integer, allocatable :: possible_values(:)
            type(SpinProj_t) :: m_s

            m_s = SpinProj_t(0)
            if (present(excess)) then
                if (abs(excess%val) == n_total_) m_s = -SpinProj_t(sign(1, excess%val))
            end if


            L(:) = this%GAS_size(spaces)
            if (m_s == beta) then
                allocate(possible_values(sum(L) .div. 2))
                counter(:) = 1
                incr = 2
            else if (m_s == alpha) then
                allocate(possible_values(sum(L) .div. 2))
                counter(:) = 2
                incr = 2
            else
                allocate(possible_values(sum(L)))
                counter(:) = 1
                incr = 1
            end if

            ! Here we merge the values from splitted_orbitals sortedly into
            ! possible_values
            i = 1
            do while (any(counter <= L))
                curr_value = huge(curr_value)
                ! foreach iGAS in spaces
                do idx_space = 1, size(spaces)
                    iGAS = spaces(idx_space)
                    if (counter(idx_space) <= L(idx_space)) then
                        if (this%get_orb_idx(counter(idx_space), iGAS) < curr_value) then
                            curr_value = this%get_orb_idx(counter(idx_space), iGAS)
                            iGAS_min_val = idx_space
                        end if
                    end if
                end do

                counter(iGAS_min_val) = counter(iGAS_min_val) + incr
                possible_values(i) = curr_value
                i = i + 1
            end do

            if (present(add_particles)) then
                possible_holes = possible_values .complement. (det_I .U. add_particles)
            else
                possible_holes = possible_values .complement. det_I
            end if
        end block
    end function

    !>  Constructor of LocalGASSpec_t
    pure function construct_LocalGASSpec_t(n_min, n_max, spat_GAS_orbs, recoupling) result(GAS_spec)
        !> Minimum particle number per GAS space.
        integer, intent(in) :: n_min(:)
        !> Maximum particle number per GAS space.
        integer, intent(in) :: n_max(:)
        !> GAS space for the i-th **spatial** orbital.
        integer, intent(in) :: spat_GAS_orbs(:)
        !> Exchange double excitations that recouple the spin are allowed
        logical, intent(in), optional :: recoupling
        logical :: recoupling_

        type(LocalGASSpec_t) :: GAS_spec
        character(*), parameter :: this_routine = 'construct_LocalGASSpec_t'

        integer :: n_spin_orbs, max_GAS_size
        integer, allocatable :: splitted_orbitals(:, :), GAS_table(:), GAS_sizes(:)
        integer :: i, iel, iGAS, nGAS

        @:def_default(recoupling_, recoupling, .true.)

        nGAS = maxval(spat_GAS_orbs)
        GAS_sizes = 2 * frequency(spat_GAS_orbs)

        if (any(min(n_max, GAS_sizes) < n_min)) then
            call stop_all(this_routine, 'any(min(n_max, GAS_sizes) < n_min) violated')
        end if


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


        ! This block and the additional declaration of new_n_max
        ! is just necessary because of shitty intel compilers (Ifort 18).
        ! If you can remove it, I am happy.
        block
            integer :: new_n_max(size(n_max))
            new_n_max = min(n_max, GAS_sizes)
            GAS_spec = LocalGASSpec_t(&
                    min=n_min, max=new_n_max, GAS_table=GAS_table, &
                    GAS_sizes=GAS_sizes, largest_GAS_size=max_GAS_size, &
                    splitted_orbitals=splitted_orbitals, &
                    lookup_is_connected=any(n_min(:) /= n_max(:)), &
                    exchange_recoupling=recoupling_)
        end block

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

    !>  Query wether a supergroup is contained in the GAS space.
    pure function Local_contains_supergroup(this, supergroup) result(res)
        !> Specification of GAS spaces.
        class(LocalGASSpec_t), intent(in) :: this
        !> A supergroup.
        integer, intent(in) :: supergroup(:)
        logical :: res
        res = all(this%min(:) <= supergroup &
                  .and. supergroup <= min(this%max(:), this%GAS_sizes(:)))
    end function

    !> Check if the GAS specification is valid
    !>
    !> If the number of particles or the number of spin orbitals
    !> is provided, then the consistency with these numbers
    !> is checked as well.
    logical pure function Local_is_valid(this, n_basis)
        !>  Specification of GAS spaces with local constraints.
        class(LocalGASSpec_t), intent(in) :: this
        integer, intent(in), optional :: n_basis

        logical :: shapes_match, pauli_principle, n_orbs_correct

        associate(GAS_sizes => this%GAS_sizes, n_min => this%min, &
                  n_max => this%max, nGAS => this%nGAS())

            shapes_match = &
                all([size(GAS_sizes), size(n_min), size(n_max)] == nGAS) &
                .and. maxval(this%GAS_table) == nGAS

            pauli_principle = all(n_min(:) <= GAS_sizes)

            if (present(n_basis)) then
                n_orbs_correct = sum(GAS_sizes) == n_basis
            else
                n_orbs_correct = .true.
            end if
        end associate

        Local_is_valid = all([shapes_match, pauli_principle, n_orbs_correct])
    end function

    subroutine Local_write_to(this, iunit)
        class(LocalGASSpec_t), intent(in) :: this
        integer, intent(in) :: iunit
        integer :: iGAS, iorb

        write(iunit, '(A)') 'Local constraints'
        write(iunit, '(A)') 'n_i: number of spatial orbitals per i-th GAS space'
        write(iunit, '(A)') 'n_min_i: minimum of particle number per i-th GAS space'
        write(iunit, '(A)') 'n_max_i: maximum of particle number per i-th GAS space'
        write(iunit, '(A10, 1x, A10, 1x, A10)') 'n_i', 'n_min_i', 'n_max_i'
        write(iunit, '(A)') '--------------------------------'
        do iGAS = 1, this%nGAS()
            write(iunit, '(I10, 1x, I10, 1x, I10)') this%GAS_size(iGAS) .div. 2, this%get_min(iGAS), this%get_max(iGAS)
        end do
        write(iunit, '(A)') '--------------------------------'
        write(iunit, '(A)') 'The distribution of spatial orbitals to GAS spaces is given by:'
        do iorb = 1, this%n_spin_orbs(), 2
            write(iunit, '(I0, 1x)', advance='no') this%get_iGAS(iorb)
        end do
        write(iunit, *)

        if (any(this%GAS_sizes < this%max)) then
            write(iunit, '(A)') 'In at least one GAS space, the maximum allowed particle number by GAS constraints'
            write(iunit, '(A)') '   is larger than the particle number allowed by the Pauli principle.'
            write(iunit, '(A)') '   Was this intended when preparing your input?'
        end if

        if (.not. this%recoupling()) then
            write(iunit, '(A)') 'Double excitations with exchange are forbidden.'
        end if
    end subroutine


    !> Returns the minimum particle number for a given GAS space.
    integer elemental function get_min_i(this, iGAS)
        class(LocalGASSpec_t), intent(in) :: this
        integer, intent(in) :: iGAS
        get_min_i = this%min(iGAS)
    end function

    !> Returns the minimum particle number for all GAS spaces.
    pure function get_min_all(this) result(res)
        class(LocalGASSpec_t), intent(in) :: this
        integer, allocatable :: res(:)
        res = this%min(:)
    end function

    !> Returns the maximum particle number for a given GAS space.
    integer elemental function get_max_i(this, iGAS)
        class(LocalGASSpec_t), intent(in) :: this
        integer, intent(in) :: iGAS
        get_max_i = this%max(iGAS)
    end function

    !> Returns the maximum particle number for all GAS spaces.
    pure function get_max_all(this) result(res)
        class(LocalGASSpec_t), intent(in) :: this
        integer, allocatable :: res(:)
        res = this%max(:)
    end function

    pure function Local_get_possible_spaces(this, supergroup, add_holes, add_particles, n_total) result(spaces)
        class(LocalGASSpec_t), intent(in) :: this
        integer, intent(in) :: supergroup(size(this%GAS_sizes))
        integer, intent(in), optional :: add_holes(:), add_particles(:), n_total
        integer, allocatable :: spaces(:)

        integer :: n_total_, iGAS
        integer :: &
        !> pumber of particles per iGAS
            n_particle(this%nGAS()), &
        !> deficit per iGAS
            deficit(this%nGAS()), &
        !> vacant orbitals per iGAS
            vacant(this%nGAS())
        type(buffer_int_1D_t) :: space_buffer

        @:def_default(n_total_, n_total, 1)

        block
            integer :: B(this%nGAS()), C(this%nGAS())
            if (present(add_holes)) then
                B = this%count_per_GAS(add_holes)
            else
                B = 0
            end if
            if (present(add_particles)) then
                C = this%count_per_GAS(add_particles)
            else
                C = 0
            end if
            n_particle = supergroup - B + C
        end block

        deficit(:) = this%get_min() - n_particle(:)
        vacant(:) = this%get_max() - n_particle(:)

        if (n_total_ < sum(deficit) .or. sum(vacant) < n_total_) then
            spaces = [integer::]
        else if (any(deficit == n_total_)) then
            spaces = [custom_findloc(deficit == n_total_, val=.true.)]
        else
            call space_buffer%init()
            do iGAS = 1, this%nGAS()
                if (vacant(iGAS) > 0) then
                    call space_buffer%push_back(iGAS)
                end if
            end do
            call space_buffer%dump_reset(spaces)
        end if
    end function


    !>  Constructor of CumulGASSpec_t
    pure function construct_CumulGASSpec_t(cn_min, cn_max, spat_GAS_orbs, recoupling) result(GAS_spec)
        !>  Cumulative minimum particle number.
        integer, intent(in) :: cn_min(:)
        !>  Cumulative maximum particle number.
        integer, intent(in) :: cn_max(:)
        !>  GAS space for the i-th **spatial** orbital.
        integer, intent(in) :: spat_GAS_orbs(:)
        !> Exchange double excitations that recouple the spin are allowed
        logical, intent(in), optional :: recoupling
        logical :: recoupling_

        type(CumulGASSpec_t) :: GAS_spec
        character(*), parameter :: this_routine = 'construct_CumulGASSpec_t'

        integer :: n_spin_orbs, max_GAS_size
        integer, allocatable :: splitted_orbitals(:, :), GAS_table(:), GAS_sizes(:)
        integer :: i, iel, iGAS, nGAS

        @:def_default(recoupling_, recoupling, .true.)

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

        GAS_spec = CumulGASSpec_t(&
                c_min=cn_min, c_max=cn_max, GAS_table=GAS_table, &
                GAS_sizes=GAS_sizes, largest_GAS_size=max_GAS_size, &
                splitted_orbitals=splitted_orbitals, &
                lookup_is_connected=any(cn_min(:) /= cn_max(:)), &
                exchange_recoupling=recoupling_)

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

    !> Query wether a supergroup is contained in the GAS space.
    pure function Cumul_contains_supergroup(this, supergroup) result(res)
        class(CumulGASSpec_t), intent(in) :: this
        integer, intent(in) :: supergroup(:)
        logical :: res
        associate(cumulated => cumsum(supergroup))
            res = all(this%c_min(:) <= cumulated &
                      .and. cumulated <= this%c_max(:) &
                      .and. supergroup <= this%GAS_sizes(:))
        end associate
    end function

    !> Check if the GAS specification is valid
    !>
    !> If the number of particles or the number of spin orbitals
    !> is provided, then the consistency with these numbers
    !> is checked as well.
    pure function Cumul_is_valid(this, n_basis) result(is_valid)
        !>  Specification of GAS spaces.
        class(CumulGASSpec_t), intent(in) :: this
        !>  The number of spin orbitals.
        integer, intent(in), optional :: n_basis
        logical :: is_valid

        logical :: shapes_match, pauli_principle, monotonic, &
            n_orbs_correct

        associate(GAS_sizes => this%GAS_size(), n_min => this%get_cmin(), &
                  n_max => this%get_cmax(), nGAS => this%nGAS())

            shapes_match = &
                all([size(GAS_sizes), size(n_min), size(n_max)] == nGAS) &
                .and. maxval(this%GAS_table) == nGAS

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

        is_valid = all([shapes_match, pauli_principle, monotonic, n_orbs_correct])
    end function

    !> Write a string representation of this GAS specification to iunit
    subroutine Cumul_write_to(this, iunit)
        class(CumulGASSpec_t), intent(in) :: this
        integer, intent(in) :: iunit
        integer :: iGAS, iorb

        write(iunit, '(A)') 'Cumulative constraints'
        write(iunit, '(A)') 'n_i: number of spatial orbitals per i-th GAS space'
        write(iunit, '(A)') 'cn_min_i: minimum of cumulative particle number per i-th GAS space'
        write(iunit, '(A)') 'cn_max_i: maximum of cumulative particle number per i-th GAS space'
        write(iunit, '(A10, 1x, A10, 1x, A10)') 'n_i', 'cn_min_i', 'cn_max_i'
        write(iunit, '(A)') '--------------------------------'
        do iGAS = 1, this%nGAS()
            write(iunit, '(I10, 1x, I10, 1x, I10)') &
                this%GAS_size(iGAS) .div. 2, this%get_cmin(iGAS), this%get_cmax(iGAS)
        end do
        write(iunit, '(A)') '--------------------------------'
        write(iunit, '(A)') 'The distribution of spatial orbitals to GAS spaces is given by:'
        do iorb = 1, this%n_spin_orbs(), 2
            write(iunit, '(I0, 1x)', advance='no') this%get_iGAS(iorb)
        end do
        write(iunit, *)

        if (any(this%GAS_size() < (this%get_cmax() - eoshift(this%get_cmin(), -1)))) then
            write(iunit, '(A)') 'In at least one GAS space, the maximum allowed particle number by GAS constraints'
            write(iunit, '(A)') '   is larger than the particle number allowed by the Pauli principle.'
            write(iunit, '(A)') '   Was this intended when preparing your input?'
        end if

        if (.not. this%recoupling()) then
            write(iunit, '(A)') 'Double excitations with exchange are forbidden.'
        end if
    end subroutine


    !> Returns the minimum particle number for a given GAS space.
    integer elemental function get_cmin_i(this, iGAS)
        class(CumulGASSpec_t), intent(in) :: this
        integer, intent(in) :: iGAS
        get_cmin_i = this%c_min(iGAS)
    end function

    !> Returns the minimum particle number for all GAS spaces.
    pure function get_cmin_all(this) result(res)
        class(CumulGASSpec_t), intent(in) :: this
        integer, allocatable :: res(:)
        res = this%c_min(:)
    end function

    !> Returns the maximum particle number for a given GAS space.
    integer elemental function get_cmax_i(this, iGAS)
        class(CumulGASSpec_t), intent(in) :: this
        integer, intent(in) :: iGAS
        get_cmax_i = this%c_max(iGAS)
    end function

    !> Returns the maximum particle number for all GAS spaces.
    pure function get_cmax_all(this) result(res)
        class(CumulGASSpec_t), intent(in) :: this
        integer, allocatable :: res(:)
        res = this%c_max(:)
    end function

    pure function Cumul_get_possible_spaces(this, supergroup, add_holes, add_particles, n_total) result(spaces)
        class(CumulGASSpec_t), intent(in) :: this
        integer, intent(in) :: supergroup(size(this%GAS_sizes))
        integer, intent(in), optional :: add_holes(:), add_particles(:), n_total
        !> Lower and upper bound for spaces where a particle can be created.
        !> If no particle can be created, then spaces == 0 .
        integer, allocatable :: spaces(:)

        integer :: n_total_, iGAS, lower_bound, upper_bound
        type(buffer_int_1D_t) :: space_buffer

        integer :: &
        !> Cumulated number of particles per iGAS
            cum_n_particle(this%nGAS()), &
        !> Cumulated deficit per iGAS
            deficit(this%nGAS()), &
        !> Cumulated vacant orbitals per iGAS
            vacant(this%nGAS())

        @:def_default(n_total_, n_total, 1)

        block
            integer :: B(this%nGAS()), C(this%nGAS())
            if (present(add_holes)) then
                B = this%count_per_GAS(add_holes)
            else
                B = 0
            end if
            if (present(add_particles)) then
                C = this%count_per_GAS(add_particles)
            else
                C = 0
            end if
            cum_n_particle = cumsum(supergroup - B + C)
        end block

        deficit(:) = this%get_cmin() - cum_n_particle(:)
        vacant(:) = this%get_cmax() - cum_n_particle(:)

        if (any(n_total_ < deficit) .or. all(vacant < n_total_)) then
            spaces = [integer :: ]
            return
        end if

        ! We assume that it is possible to create a particle at least in
        ! the last GAS space.
        ! Search from behind the first occurence where it is not possible
        ! anymore to create a particle.
        ! The lower bound is one GAS index above.
        lower_bound = custom_findloc(vacant <= 0, .true., back=.true.) + 1
        ! Find the first index, where a particle has to be created.
        upper_bound = custom_findloc(deficit, n_total_)

        if (lower_bound > upper_bound .or. lower_bound > this%nGAS()) then
            spaces = [integer :: ]
        else
        block
            integer :: new_deficit(size(deficit)), new_vacant(size(vacant)), jGAS
            logical :: allowed
            call space_buffer%init()
            do iGAS = lower_bound, upper_bound
                new_deficit = deficit
                new_deficit(iGAS :) = new_deficit(iGAS :) - 1
                new_vacant = vacant
                new_vacant(iGAS :) = new_vacant(iGAS :) - 1
                allowed = .true.
                do jGAS = 1, size(new_deficit) - 1
                    if (any(new_deficit(jGAS) > new_vacant(jGAS + 1 : ))) then
                        allowed = .false.
                        exit
                    end if
                end do
                if (allowed) call space_buffer%push_back(iGAS)
            end do
            call space_buffer%dump_reset(spaces)
        end block
        end if
    end function
end module gasci
