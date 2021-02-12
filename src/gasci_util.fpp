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
    use bit_reps, only: decode_bit_det
    implicit none
    private
    public :: get_available_singles, get_available_doubles, gen_all_excits, gen_all_excits_wrapper

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
        integer, intent(in) :: nI(:)
        integer(n_int), intent(out), allocatable :: det_list(:, :)
        integer, optional, intent(in) :: ic
        character(*), parameter :: this_routine = "gen_excits"

        integer :: N
        integer :: FCI_n_excits
        integer(n_int), allocatable :: FCI_det_list(:, :), tmp_list(:, :)
        integer :: i

        @:ASSERT(size(nI) == nEL)
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

        @:ASSERT(GAS_spec%contains_det(det_I))
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
