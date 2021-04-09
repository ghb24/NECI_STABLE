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
    use bit_reps, only: decode_bit_det
    implicit none
    private
    public :: get_available_singles, get_available_doubles, gen_all_excits

    public :: get_cumulative_list, draw_from_cum_list

    interface get_cumulative_list
        #:for Excitation_t in ExcitationTypes
            module procedure get_cumulative_list_${Excitation_t}$
        #:endfor
    end interface

contains


    !>  @brief
    !>  Get all single excitated determinants from det_I that are allowed under GAS constraints.
    pure function get_available_singles(GAS_spec, det_I) result(singles_exc_list)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: singles_exc_list(:, :)
        character(*), parameter :: this_routine = 'get_available_singles'

        integer, allocatable :: possible_holes(:)
        integer :: i, j, src1, tgt1
        type(SpinProj_t) :: m_s_1
        type(buffer_int_2D_t) :: buffer


        @:pure_ASSERT(GAS_spec%contains_det(det_I))

        call buffer%init(size(det_I))

        do i = 1, size(det_I)
            src1 = det_I(i)
            m_s_1 = calc_spin_raw(src1)
            possible_holes = GAS_spec%get_possible_holes(&
                        det_I, add_holes=det_I(i:i), &
                        excess=-m_s_1)
            do j = 1, size(possible_holes)
                tgt1 = possible_holes(j)
                call buffer%push_back(excite(det_I, SingleExc_t(src1, tgt1)))
            end do
        end do
        call buffer%dump_reset(singles_exc_list)

        @:sort(integer, singles_exc_list, rank=2, along=2, comp=lex_leq)
    end function

    !>  @brief
    !>  Get all double excitated determinants from det_I that are allowed under GAS constraints.
    pure function get_available_doubles(GAS_spec, det_I) result(doubles_exc_list)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: doubles_exc_list(:, :)
        character(*), parameter :: this_routine = 'get_available_doubles'

        integer, allocatable :: first_pick_possible_holes(:), second_pick_possible_holes(:), deleted(:)
        integer :: i, j, k, l, src1, src2, tgt1, tgt2
        type(SpinProj_t) :: m_s_1
        type(buffer_int_2D_t) :: buffer

        @:pure_ASSERT(GAS_spec%contains_det(det_I))

        call buffer%init(size(det_I))

        do i = 1, size(det_I)
            do j = i + 1, size(det_I)
                src1 = det_I(i)
                src2 = det_I(j)
                deleted = det_I([i, j])
                first_pick_possible_holes = GAS_spec%get_possible_holes(det_I, &
                                        add_holes=deleted, excess=-sum(calc_spin_raw(deleted)), &
                                        n_total=2)
                @:pure_ASSERT(disjoint(first_pick_possible_holes, det_I))
                do k = 1, size(first_pick_possible_holes)
                    tgt1 = first_pick_possible_holes(k)
                    m_s_1 = calc_spin_raw(tgt1)
                    @:pure_ASSERT(any(m_s_1 == [alpha, beta]))

                    second_pick_possible_holes = GAS_spec%get_possible_holes(&
                            det_I, add_holes=deleted, &
                            add_particles=[tgt1], &
                            n_total=1, excess=m_s_1 - sum(calc_spin_raw(deleted)))

                    @:pure_ASSERT(disjoint(second_pick_possible_holes, [tgt1]))
                    @:pure_ASSERT(disjoint(second_pick_possible_holes, det_I))

                    do l = 1, size(second_pick_possible_holes)
                        tgt2 = second_pick_possible_holes(l)
                        call buffer%push_back(excite(det_I, DoubleExc_t(src1, tgt1, src2, tgt2)))
                    end do
                end do
            end do
        end do

        call buffer%dump_reset(doubles_exc_list)

        @:sort(integer, doubles_exc_list, rank=2, along=2, comp=lex_leq)

        remove_double_appearances : block
            integer, allocatable :: tmp_buffer(:, :)
            allocate(tmp_buffer(size(doubles_exc_list, 1), size(doubles_exc_list, 2)))
            j = 1
            tmp_buffer(:, j) = doubles_exc_list(:, 1)
            do i = 2, size(doubles_exc_list, 2)
                if (any(doubles_exc_list(:, i - 1) /= doubles_exc_list(:, i))) then
                    j = j + 1
                    tmp_buffer(:, j) = doubles_exc_list(:, i)
                end if
            end do
            doubles_exc_list = tmp_buffer(:, : j)
        end block remove_double_appearances
    end function


    !>  @brief
    !>  Get all excitated determinants from
    !>  det_I that are allowed under GAS constraints.
    !>
    !>  @param[in] GAS_spec GAS specification
    !>  @param[in] nI Starting determinant
    !>  @param[out] n_excits Number of determinants
    !>  @param[out] det_list Allocatable array of determinants in ilut format
    !>  @param[in] ic Optional input for excitation level (ic=1 => singles, ic=2 => doubles)
    !>      If ommited generate all.
    subroutine gen_all_excits(GAS_spec, nI, n_excits, det_list, ic)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: nI(:)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        integer, intent(in), optional :: ic

        integer, allocatable :: singles(:, :), doubles(:, :)
        integer :: i, j

        if (present(ic)) then
            select case(ic)
            case(1)
                singles = get_available_singles(GAS_spec, nI)
                allocate(doubles(0, 0))
            case(2)
                allocate(singles(0, 0))
                doubles = get_available_doubles(GAS_spec, nI)
            end select
        else
            singles = get_available_singles(GAS_spec, nI)
            doubles = get_available_doubles(GAS_spec, nI)
        end if

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
