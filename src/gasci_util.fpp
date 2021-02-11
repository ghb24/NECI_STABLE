#:include "macros.fpph"
#:include "algorithms.fpph"


#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

!> This module contains functions for GAS that are not bound
!>  to a specific GAS excitation generator.
module gasci_util
    use constants, only: n_int, dp
    use SystemData, only: nEl
    use gasci, only: GASSpec_t, LocalGASSpec_t, CumulGASSpec_t, GAS_specification
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==), operator(/=), operator(-), sum, &
        alpha, beta
    use sort_mod, only: sort
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite, get_last_tgt, set_last_tgt, UNKNOWN
    use util_mod, only: lex_leq, cumsum, operator(.div.), near_zero, binary_search_first_ge, &
        operator(.isclose.), custom_findloc
    use dSFMT_interface, only: genrand_real2_dSFMT
    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet
    use bit_rep_data, only: NIfTot, NIfD
    use sltcnd_mod, only: sltcnd_excit
    use sets_mod, only: disjoint, union, complement, is_sorted
    use growing_buffers, only: buffer_int_2D_t, buffer_int_1D_t
    use symexcit3, only: FCI_gen_all_excits => gen_excits
    implicit none
    private
    public :: get_available_singles, get_available_doubles, &
        get_possible_spaces, get_possible_holes, gen_all_excits, gen_all_excits_wrapper

    public :: get_cumulative_list, draw_from_cum_list

    interface get_cumulative_list
        #:for Excitation_t in ExcitationTypes
            module procedure get_cumulative_list_${Excitation_t}$
        #:endfor
    end interface

contains

    !>   @brief
    !>   Return all configurations that are connected to nI as
    !>   array of iluts (det_list(0:niftot, n_excits)).
    !>
    !>  @details
    !>  The routine is not efficient!
    !>  Improve performance if used in tight code.
    !>  Triple excitations are not supported.
    !>
    !>  @param[in] nI, The configuration from which to excite.
    !>  @param[out] det_list, The connected configurations in ilut format.
    !>                  (det_list(0:niftot, n_excits))
    !>  @param[in] ex_flag, The requested excitations. (1 = singles, 2 = doubles)
    !>          If ommited all excitations will be generated.
    subroutine gen_excits(GAS_spec, nI, det_list, ic)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(out), allocatable :: det_list(:, :)
        integer, optional, intent(in) :: ic
        character(*), parameter :: this_routine = "gen_excits"

        integer :: N
        integer :: FCI_n_excits
        integer(n_int), allocatable :: FCI_det_list(:, :), tmp_list(:, :)
        integer :: i

        @:ASSERT(GAS_spec%contains_det(nI))
        call FCI_gen_all_excits(nI, FCI_n_excits, FCI_det_list, ic)
        allocate(tmp_list, mold=FCI_det_list)

        N = 0
        do i = 1, FCI_n_excits
            if (GAS_spec%contains_ilut(FCI_det_list(:, i))) then
                N = N + 1
                tmp_list(:, N) = FCI_det_list(:, i)
            end if
        end do

        allocate(det_list(0 : nIfTot, N))
        det_list(:, :) = tmp_list(:, : N)
    end subroutine


    !>  @brief
    !>  Get all single excitated determinants from det_I that are allowed under GAS constraints.
    function get_available_singles(GAS_spec, det_I) result(singles_exc_list)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: singles_exc_list(:, :)
        character(*), parameter :: this_routine = 'get_available_singles'

        integer(n_int), allocatable :: det_list(:, :)
        integer :: i

        @:pure_ASSERT(GAS_spec%contains_det(det_I))
        call gen_excits(GAS_spec, det_I, det_list, ic=1)
        allocate(singles_exc_list(sum(popcnt(det_list(:, 1))), size(det_list, 2)))
        do i = 1, size(det_list, 2)
            call decode_bit_det(singles_exc_list(:, i), det_list(:, i))
        end do
        @:sort(integer, singles_exc_list, rank=2, along=2, comp=lex_leq)
    end function

    !>  @brief
    !>  Get all double excitated determinants from det_I that are allowed under GAS constraints.
    function get_available_doubles(GAS_spec, det_I) result(doubles_exc_list)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: doubles_exc_list(:, :)
        character(*), parameter :: this_routine = 'get_available_doubles'

        integer(n_int), allocatable :: det_list(:, :)
        integer :: i

        @:pure_ASSERT(GAS_spec%contains_det(det_I))
        call gen_excits(GAS_spec, det_I, det_list, ic=2)
        allocate(doubles_exc_list(sum(popcnt(det_list(:, 1))), size(det_list, 2)))
        do i = 1, size(det_list, 2)
            call decode_bit_det(doubles_exc_list(:, i), det_list(:, i))
        end do
        @:sort(integer, doubles_exc_list, rank=2, along=2, comp=lex_leq)
    end function


    !>  @brief
    !>  Get all excitated determinants from det_I that are allowed under GAS constraints.
    subroutine gen_all_excits(GAS_spec, nI, n_excits, det_list)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: nI(:)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)

        integer, allocatable :: singles(:, :), doubles(:, :)
        integer :: i, j

        singles = get_available_singles(GAS_spec, nI)
        doubles = get_available_doubles(GAS_spec, nI)

        n_excits = size(singles, 2) + size(doubles, 2)
        allocate(det_list(0:niftot, n_excits))
        j = 1
        do i = 1, size(singles, 2)
            call EncodeBitDet(singles(:, i), det_list(:, j))
            j = j + 1
        end do

        do i = 1, size(doubles, 2)
            call EncodeBitDet(doubles(:, i), det_list(:, j))
            j = j + 1
        end do

        call sort(det_list, ilut_lt, ilut_gt)
    end subroutine gen_all_excits

    subroutine gen_all_excits_wrapper(nI, n_excits, det_list)
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)

        call gen_all_excits(GAS_specification, nI, n_excits, det_list)
    end subroutine gen_all_excits_wrapper



    !>  @brief
    !>      Return the GAS spaces, where one particle can be created.
    !>
    !>  @details
    !>  It can be proven, that the spaces where a particle can
    !>  be created are contigous, so a two element integer array
    !>  [lower_bound, upper_bound] is returned.
    !>  As long as lower_bound <= iGAS <= upper_bound is true,
    !>  the created particle will lead to a valid Slater-Determinant.
    !>
    !>  It is **assumed** that the input determinant is contained in the
    !>  Full CI space and obeys e.g. the Pauli principle.
    !>  Checks are only performed in DEBUG compilation mode and
    !>  the return value is undefined, if this is not the case!
    !>
    !>  It is possible to delete additional particles with the
    !>  optional argument `add_holes` before checking the
    !>  validity of particle creation.
    !>  On the other hand it is possible to create particles before
    !>  checking the validity of particle creation.
    !>
    !>  If more than one particle should be created, the optional argument
    !>  n_total (default 1) should be used.
    !>
    !>  If no creation is allowed by the GAS constraints,
    !>  the bounds will be returned as integer constant EMPTY_BOUNDS.
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] particles_per_GAS, The particles per GAS space.
    !>  @param[in] add_holes, optional, An index of orbitals
    !>      where particles should be deleted before creating the new particle.
    !>  @param[in] add_particles, optional, An index of orbitals
    !>      where particles should be created before creating the new particle.
    !>  @param[in] n_total, optional, The total number of particles
    !>      that will be created. Defaults to one (integer).
    pure function get_possible_spaces(GAS_spec, particles_per_GAS, add_holes, add_particles, n_total) result(spaces)
        type(LocalGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: particles_per_GAS(GAS_spec%nGAS())
        integer, intent(in), optional :: add_holes(:), add_particles(:), n_total
        integer, allocatable :: spaces(:)

        integer :: n_total_, iGAS
        integer :: &
        !> pumber of particles per iGAS
            n_particle(GAS_spec%nGAS()), &
        !> deficit per iGAS
            deficit(GAS_spec%nGAS()), &
        !> vacant orbitals per iGAS
            vacant(GAS_spec%nGAS())
        type(buffer_int_1D_t) :: space_buffer

        @:def_default(n_total_, n_total, 1)

        block
            integer :: B(GAS_spec%nGAS()), C(GAS_spec%nGAS())
            if (present(add_holes)) then
                B = GAS_spec%count_per_GAS(add_holes)
            else
                B = 0
            end if
            if (present(add_particles)) then
                C = GAS_spec%count_per_GAS(add_particles)
            else
                C = 0
            end if
            n_particle = particles_per_GAS - B + C
        end block

        deficit(:) = GAS_spec%get_min() - n_particle(:)
        vacant(:) = GAS_spec%get_max() - n_particle(:)

        if (n_total_ < sum(deficit) .or. sum(vacant) < n_total_) then
            spaces = [integer::]
        else if (any(deficit == n_total_)) then
            spaces = [custom_findloc(deficit == n_total_, val=.true.)]
        else
            call space_buffer%init()
            do iGAS = 1, GAS_spec%nGAS()
                if (vacant(iGAS) > 0) then
                    call space_buffer%push_back(iGAS)
                end if
            end do
            call space_buffer%dump_reset(spaces)
        end if
    end function


    !>  @brief
    !>  Return the possible holes where a particle can be created under GAS constraints.
    !>
    !>  @details
    !>  This function uses `get_possible_spaces` to find possible GAS spaces
    !>  where a particle can be created and returns only unoccupied
    !>  sites of correct spin.
    !>
    !>  "Trivial" excitations are avoided. That means, that a site is only counted
    !>  as unoccupied if it was unoccupied in nI from the beginning on.
    !>  (A double excitation where a particle is deleted, but immediately
    !>  recreated would be such a trivial excitations.)
    !>
    !>  @param[in] GAS_spec, Specification of GAS spaces (GASSpec_t).
    !>  @param[in] particles_per_GAS, The particles per GAS space.
    !>  @param[in] add_holes, optional, An index of orbitals
    !>      where particles should be deleted before creating the new particle.
    !>  @param[in] add_particles, optional, An index of orbitals
    !>      where particles should be created before creating the new particle.
    !>  @param[in] n_total, optional, The total number of particles
    !>      that will be created. Defaults to one (integer).
    !>  @param[in] excess, optional, The current excess of spin projections.
    !>      If a beta electron was deleted, the excess is 1 * alpha.
    pure function get_possible_holes(GAS_spec, det_I, add_holes, add_particles, n_total, excess) result(possible_holes)
        type(LocalGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, intent(in), optional :: add_holes(:)
        integer, intent(in), optional :: add_particles(:)
        integer, intent(in), optional :: n_total
        type(SpinProj_t), intent(in), optional :: excess
        character(*), parameter :: this_routine = 'get_possible_holes'

        integer, allocatable :: possible_holes(:)

        integer :: &
            splitted(GAS_spec%max_GAS_size(), GAS_spec%nGAS()), &
            splitted_sizes(GAS_spec%nGAS())

        integer, allocatable :: spaces(:)
        integer :: n_total_

        @:def_default(n_total_, n_total, 1)
        @:pure_ASSERT(1 <= n_total_)
        if (present(excess)) then
            @:pure_ASSERT(abs(excess%val) <= n_total_)
        end if

        call GAS_spec%split_per_GAS(det_I, splitted, splitted_sizes)

        ! Note, that non-present optional arguments can be passed
        ! into optional arguments without checking!
        spaces = get_possible_spaces(&
             GAS_spec, splitted_sizes, add_holes=add_holes, &
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


            L(:) = GAS_spec%GAS_size(spaces)
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
                        if (GAS_spec%get_orb_idx(counter(idx_space), iGAS) < curr_value) then
                            curr_value = GAS_spec%get_orb_idx(counter(idx_space), iGAS)
                            iGAS_min_val = idx_space
                        end if
                    end if
                end do

                counter(iGAS_min_val) = counter(iGAS_min_val) + incr
                possible_values(i) = curr_value
                i = i + 1
            end do

            if (present(add_particles)) then
                possible_holes = complement(possible_values, union(det_I, add_particles))
            else
                possible_holes = complement(possible_values, det_I)
            end if
        end block
    end function


#:for excitation_t in ExcitationTypes
    !>  @brief
    !>  Build up a cumulative list of matrix elements.
    !>
    !>  @details
    !>  Calculate the matrix elements for the possible excitations from det_I
    !>  to the possible holes using the incomplete defined excitation.
    !>
    !>  @param[in] det_I, Reference determinant in "nI-format".
    !>  @param[in] incomplete_exc, An excitation where the last target is unknown.
    !>  @param[in] possible_holes, Possible holes for the last target.
    function get_cumulative_list_${excitation_t}$(det_I, incomplete_exc, possible_holes) result(cSum)
        integer, intent(in) :: det_I(:)
        type(${excitation_t}$), intent(in) :: incomplete_exc
        integer, intent(in) :: possible_holes(:)
        real(dp) :: cSum(size(possible_holes))
        character(*), parameter :: this_routine = 'get_cumulative_list_${excitation_t}$'

        real(dp) :: previous
        type(${excitation_t}$) :: exc
        integer :: i

        @:ASSERT(get_last_tgt(exc) == UNKNOWN)
        exc = incomplete_exc

        ! build the cumulative list of matrix elements <src|H|tgt>
        previous = 0.0_dp
        do i = 1, size(possible_holes)
            call set_last_tgt(exc, possible_holes(i))
            cSum(i) = abs(sltcnd_excit(det_I, exc, .false.)) + previous
            previous = cSum(i)
        end do

        ! Normalize
        if (near_zero(cSum(size(cSum)))) then
            cSum(:) = 0.0_dp
        else
            cSum(:) = cSum(:) / cSum(size(cSum))
        end if
    end function get_cumulative_list_${excitation_t}$
#:endfor


    !>  @brief
    !>  Draw from a cumulative list.
    subroutine draw_from_cum_list(c_sum, idx, pgen)
        real(dp), intent(in) :: c_sum(:)
        integer, intent(out) :: idx
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'draw_from_cum_list'

        @:pure_ASSERT((c_sum(size(c_sum)) .isclose. 0.0_dp) &
            .or. (c_sum(size(c_sum)) .isclose. 1.0_dp))
        @:pure_ASSERT(is_sorted(c_sum))

        ! there might not be such an excitation
        if (c_sum(size(c_sum)) > 0) then
            ! find the index of the target hole in the cumulative list
            idx = binary_search_first_ge(c_sum, genrand_real2_dSFMT())
            @:pure_ASSERT(1 <= idx .and. idx <= size(c_sum))

            ! adjust pgen with the probability for picking tgt from the cumulative list
            if (idx == 1) then
                pgen = c_sum(1)
            else
                pgen = c_sum(idx) - c_sum(idx - 1)
            end if
        else
            idx = 0
            pgen = 0.0_dp
        end if
    end subroutine

end module
