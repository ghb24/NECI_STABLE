#:include "macros.fpph"
#:include "algorithms.fpph"


#:set ExcitationTypes = ['Excite_1_t', 'Excite_2_t']

!> This module contains functions for GAS that are not bound
!>  to a specific GAS excitation generator.
module gasci_util
    use constants, only: n_int, dp, int64
    use SystemData, only: nEl, nBasis
    use gasci, only: GASSpec_t, LocalGASSpec_t, CumulGASSpec_t, GAS_specification
    use gasci_supergroup_index, only: get_supergroups
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, sum, alpha, beta
    use sort_mod, only: sort
    use excitation_types, only: Excite_1_t, Excite_2_t, excite, get_last_tgt, set_last_tgt, UNKNOWN
    use util_mod, only: lex_leq, cumsum, operator(.div.), near_zero, binary_search_first_ge, &
        operator(.isclose.), custom_findloc, choose_i64
    use dSFMT_interface, only: genrand_real2_dSFMT
    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet
    use bit_rep_data, only: NIfTot, NIfD
    use sltcnd_mod, only: sltcnd_excit
    use sets_mod, only: disjoint, is_sorted
    use growing_buffers, only: buffer_int_2D_t, buffer_int_1D_t
    use bit_reps, only: decode_bit_det
    implicit none
    private
    public :: get_available_singles, get_available_doubles, &
        gen_all_excits, get_n_SDs


    public :: get_cumulative_list, draw_from_cum_list

    public :: write_GAS_info, t_output_GAS_sizes

    interface get_cumulative_list
        #:for Excitation_t in ExcitationTypes
            module procedure get_cumulative_list_${Excitation_t}$
        #:endfor
    end interface

    logical :: t_output_GAS_sizes = .false.

contains


    !>  @brief
    !>  Get all single excitated determinants from det_I that are allowed under GAS constraints.
    pure function get_available_singles(GAS_spec, det_I) result(singles_exc_list)
        class(GASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, allocatable :: singles_exc_list(:, :)
            !! Dimension is (nEl, n_configurations)
        character(*), parameter :: this_routine = 'get_available_singles'

        integer, allocatable :: possible_holes(:)
        integer :: i, j, src1, tgt1
        type(SpinProj_t) :: m_s_1
        type(buffer_int_2D_t) :: buffer


        @:pure_ASSERT(GAS_spec%contains_conf(det_I))

        call buffer%init(size(det_I))

        do i = 1, size(det_I)
            src1 = det_I(i)
            m_s_1 = calc_spin_raw(src1)
            possible_holes = GAS_spec%get_possible_holes(&
                        det_I, add_holes=det_I(i:i), &
                        excess=-m_s_1)
            do j = 1, size(possible_holes)
                tgt1 = possible_holes(j)
                call buffer%push_back(excite(det_I, Excite_1_t(src1, tgt1)))
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

        @:pure_ASSERT(GAS_spec%contains_conf(det_I))

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
                        call buffer%push_back(excite(det_I, Excite_2_t(src1, tgt1, src2, tgt2)))
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

    pure function get_alpha_supergroups(sg, n_orbs, S_z) result(sg_alpha)
        !! Return the possible supergroups/distributions for alpha electrons.
        !!
        !! If `sg_alpha` and `sg_beta` are the distributions of alpha/beta electrons among the
        !! GAS spaces. Then we have `sg(:) = sg_alpha(:) + sg_beta(:)`.
        !! We want to generate all possible `sg_alpha` such that
        !! The GAS constraints are still fullfilled and the total
        !! spin projection is maintained. (`sum(sg_alpha) + sum(sg_beta) == 2*S_z`)
        integer, intent(in) :: sg(:)
            !! The overall supergroup
        integer, intent(in) :: n_orbs(size(sg))
            !! The number of spatial orbitals per GAS space
        type(SpinProj_t), intent(in) :: S_z
            !! The Spin projection
        integer, allocatable :: sg_alpha(:, :)
            !! All possible distributions of alpha electrons among the
            !!      GAS spaces.
        integer :: N, i, j

        N = sum(sg)
        if (mod(S_z%val + N, 2) == 0) then
        block
            type(LocalGASSpec_t) :: sg_alpha_constraint
            integer :: N_alpha
            N_alpha = (N + S_z%val) .div. 2
            sg_alpha_constraint = LocalGASSpec_t(&
                                n_min=max(sg - n_orbs, 0), &
                                n_max=min(sg, n_orbs), &
                                spat_GAS_orbs=[((j, i = 1, n_orbs(j)), j = 1, size(n_orbs))])
            sg_alpha = get_supergroups(sg_alpha_constraint, N_alpha)
        end block
        else
            allocate(sg_alpha(0, 0))
        end if
    end function

    elemental function get_n_SDs(GAS_spec, N, S_z) result(n_SDs)
        !! Return the number of Slater-determinants.
        class(GASSpec_t), intent(in) :: GAS_spec
            !! GAS specification.
        integer, intent(in) :: N
            !! The number of particles
        type(SpinProj_t), intent(in) :: S_z
            !! Spin projection
        integer(int64) :: n_SDs
        integer, allocatable :: supergroups(:, :), alpha_supergroups(:, :)
        integer :: i, j, iGAS

        supergroups = get_supergroups(GAS_spec, N)

        block
            integer :: N_alpha(GAS_spec%nGAS()), N_beta(GAS_spec%nGAS()), n_spat_orbs(GAS_spec%nGAS())
                !! These are the number of alpha/beta electrons and the number of
                !!  spatial orbitals per GAS space
            n_spat_orbs = GAS_spec%GAS_size() .div. 2

            n_SDs = 0_int64
            do j = 1, size(supergroups, 2)
                alpha_supergroups = get_alpha_supergroups(supergroups(:, j), n_spat_orbs, S_z)
                do i = 1, size(alpha_supergroups, 2)
                    N_alpha = alpha_supergroups(:, i)
                    N_beta = supergroups(:, j) - N_alpha(:)
                    ! For a given supergroup and S_z in each GAS space
                    !   (coming from the distribution of alpha electrons)
                    !   the number of possible configurations in each GAS space is calculated.
                    ! The overall number is just the product over the GAS spaces.
                    n_SDs = n_SDs + product([(choose_i64(n_spat_orbs(iGAS), N_alpha(iGAS)) &
                                            * choose_i64(n_spat_orbs(iGAS), N_beta(iGAS)), iGAS = 1, GAS_spec%nGAS())])
                end do
            end do
        end block
    end function

    subroutine write_GAS_info(GAS_spec, N, S_z, iunit)
        !! Write info about the GAS constraints to `iunit`
        !!
        !! The routine especially compares the CAS and GAS Hilbert space sizes.
        class(GASSpec_t), intent(in) :: GAS_spec
            !! GAS constraints.
        integer, intent(in) :: N
            !! The particle number.
        type(SpinProj_t), intent(in) :: S_z
            !! The total spin projection.
        integer, intent(in) :: iunit

        call GAS_spec%write_to(iunit)

        @:unused_var(N, S_z)

        if (t_output_GAS_sizes) then
        block
            integer :: n_alpha, n_beta, n_spat_orbs
            integer(int64) :: size_CAS, size_GAS
            N_alpha = (N + S_z%val) .div. 2
            N_beta = N - N_alpha
            n_spat_orbs = GAS_spec%n_spin_orbs() .div. 2
            size_CAS = choose_i64(n_spat_orbs, N_alpha) * choose_i64(n_spat_orbs, N_beta)
            write(iunit, '(A, 1x, I0)') 'The size of the CAS space is:', size_CAS
            size_GAS = get_n_SDs(GAS_spec, nEl, S_z)
            write(iunit, '(A, 1x, I0)') 'The size of the GAS space is:', size_GAS
            write(iunit, '(A, 1x, E10.5)') 'The fraction of the GAS space is:', real(size_GAS, dp) / real(size_CAS, dp)
        end block
        end if
    end subroutine
end module
