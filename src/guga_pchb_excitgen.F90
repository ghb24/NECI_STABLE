#include "macros.h"
module guga_pchb_excitgen

    use aliasSampling, only: AliasSampler_1D_t
    use constants, only: n_int, dp, maxExcit, int64, iout, int_rdm
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel, G1, current_stepvector, t_pchb_weighted_singles, &
                          nBasis, nSpatOrbs, ElecPairs, currentOcc_int, &
                          t_analyze_pchb, t_old_pchb, t_exchange_pchb
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles, MaxTau
    use guga_data, only: tNewDet, ExcitationInformation_t, gen_type, excit_type
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: near_zero, fuseIndex, intswap, binary_search_first_ge, &
                        get_free_unit
    use CalcData, only: t_matele_cutoff, matele_cutoff, frq_ratio_cutoff, &
                        max_frequency_bound, n_frequency_bins, &
                        t_hist_tau_search, t_truncate_spawns
    use sym_general_mod, only: ClassCountInd
    use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel
    use UMatCache, only: gtID
    use guga_excitations, only: assign_excitinfo_values_single, &
                                createStochasticExcitation_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, get_guga_integral_contrib_spat, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo
    use guga_procedure_pointers, only: gen_single_excit_guga, gen_double_excit_guga
    use guga_bitrepops, only: identify_excitation, encode_excit_info, extract_excit_info, &
                              contract_2_rdm_ind
    use bit_reps, only: decode_bit_det
    use shared_array, only: shared_array_int64_t, shared_array_real_t
    use MPI_wrapper, only: iProcIndex_intra, iprocindex
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use procedure_pointers, only: get_umat_el


    implicit none

    private

    public :: pick_orbitals_double_pchb, pick_orbitals_pure_uniform_singles, &
              calc_orbital_pgen_contr_pchb, calc_orbital_pgen_contr_start_pchb, &
              calc_orbital_pgen_contr_end_pchb, init_guga_pchb_excitgen, &
              calc_pgen_guga_pchb, pick_uniform_spatial_hole, &
              calc_orb_pgen_uniform_singles, setup_pchb_sampler_conditional, &
              tgtOrbs, guga_pchb_sampler, finalize_pchb_excitgen_guga, &
              store_pchb_analysis

    type :: GugaAliasSampler_t

        private

        type(AliasSampler_1D_t) :: alias_sampler

        ! i don't need to make two, but lets keep it similar to Kais
        ! implementation
        type(shared_array_int64_t) :: all_info_table
        type(shared_array_int64_t), allocatable :: info_tables(:)


        ! for analysis also keep track: (will be removed after optimization)
        type(shared_array_int64_t) :: all_counts
        type(shared_array_int64_t), allocatable :: count_tables(:)

        type(shared_array_int64_t) :: all_invalids
        type(shared_array_int64_t), allocatable :: invalid_tables(:)

        type(shared_array_real_t) :: all_sums
        type(shared_array_real_t), allocatable :: sums_tables(:)

        type(shared_array_real_t) :: all_worst_orb
        type(shared_array_real_t), allocatable :: worst_orb_table(:)

        type(shared_array_real_t) :: all_high_pgen
        type(shared_array_real_t), allocatable :: high_pgen_table(:)

        type(shared_array_real_t) :: all_pgen
        type(shared_array_real_t), allocatable :: pgen_table(:)

        type(shared_array_real_t) :: all_low_pgen
        type(shared_array_real_t), allocatable :: low_pgen_table(:)

    contains

        private

        procedure :: setup_info_table
        procedure :: setup_entry_info
        procedure :: info_table_destructor
        procedure :: get_info

        ! for analysis (will be removed after optimization)

        procedure :: setup_count_table
        procedure :: setup_entry_count_vec
        procedure :: setup_entry_count_scalar
        procedure :: count_table_destructor
        procedure :: get_count

        procedure :: setup_invalid_table
        procedure :: setup_entry_invalid_vec
        procedure :: setup_entry_invalid_scalar
        procedure :: invalid_table_destructor
        procedure :: get_invalid

        procedure :: setup_sum_table
        procedure :: setup_entry_sum_vec
        procedure :: setup_entry_sum_scalar
        procedure :: sum_table_destructor
        procedure :: get_sum

        procedure :: setup_worst_orb_table
        procedure :: setup_entry_worst_orb_vec
        procedure :: setup_entry_worst_orb_scalar
        procedure :: worst_orb_table_destructor
        procedure :: get_worst_orb

        procedure :: setup_pgen_table
        procedure :: setup_entry_pgen_vec
        procedure :: setup_entry_pgen_scalar
        procedure :: pgen_table_destructor
        procedure :: get_pgen

        procedure :: setup_high_pgen_table
        procedure :: setup_entry_pgen_high_vec
        procedure :: setup_entry_pgen_high_scalar
        procedure :: high_pgen_table_destructor
        procedure :: get_high_pgen

        procedure :: setup_low_pgen_table
        procedure :: setup_entry_pgen_low_vec
        procedure :: setup_entry_pgen_low_scalar
        procedure :: low_pgen_table_destructor
        procedure :: get_low_pgen



    end type GugaAliasSampler_t


    ! start with one sampler for now!
    type(GugaAliasSampler_t) :: guga_pchb_sampler
    integer, allocatable :: tgtOrbs(:,:)
    interface calc_orb_pgen_uniform_singles
        module procedure calc_orb_pgen_uniform_singles_exmat
        module procedure calc_orb_pgen_uniform_singles_excitInfo
    end interface calc_orb_pgen_uniform_singles

    interface calc_orb_pgen_guga_pchb_double
        ! module procedure calc_orb_pgen_guga_pchb_double_exmat
        module procedure calc_orb_pgen_guga_pchb_double_excitInfo
    end interface calc_orb_pgen_guga_pchb_double

contains

    function get_pchb_integral_contrib(i, j, a, b, typ) result(integral)
        ! specialized function to obtain the guga-integral contrib for
        ! the pchb weights
        integer, intent(in) :: a, i, b, j, typ
        real(dp) :: integral
        debug_function_name("get_pchb_integral_contrib")
        logical :: flag_
        real(dp) :: cpt1, cpt2, cpt3, cpt4

        ASSERT(a > 0 .and. a <= nSpatOrbs)
        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(b > 0 .and. b <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)


        if (t_old_pchb) then
            integral = get_guga_integral_contrib_spat([i,j],a,b)
            return
        end if

        flag_ = t_exchange_pchb

        select case (typ)

        ! i need to get the correct indices for the integral contributions!
        case (excit_type%single_overlap_L_to_R)

            integral = abs(get_umat_el(a, a, i, j) + get_umat_el(a, a, j, i)) / 2.0_dp

        case (excit_type%single_overlap_R_to_L)

            integral = abs(get_umat_el(a, b, i, i) + get_umat_el(b, a, i, i)) / 2.0_dp

        case (excit_type%double_lowering, excit_type%double_raising)

            cpt1 = (get_umat_el(a, b, i, j) + get_umat_el(b, a, j, i))
            cpt2 = (get_umat_el(b, a, i, j) + get_umat_el(a, b, j, i))
            cpt3 = abs(cpt1 + cpt2) / 2.0_dp
            cpt4 = abs(cpt1 - cpt2) / 2.0_dp

            cpt1 = abs(cpt1)
            cpt2 = abs(cpt2)

            if (flag_) then
                integral = cpt4
            else
                integral = maxval([cpt1, cpt2, cpt3, cpt4])
            end if

        case (excit_type%double_L_to_R_to_L, excit_type%double_R_to_L_to_R, &
              excit_type%double_L_to_R, excit_type%double_R_to_L)

            cpt1 = (get_umat_el(a, b, j, i) + get_umat_el(b, a, i, j))
            cpt2 = (get_umat_el(b, a, j, i) + get_umat_el(a, b, i, j))

            cpt3 = abs(cpt2 - cpt1)
            ! cpt4 = abs(cpt1 + cpt2)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1)/2.0_dp

            if (flag_) then
                integral = abs(cpt2)
            else
                integral = maxval([abs(cpt1),abs(cpt2)/2.0_dp, cpt3, cpt4])
            end if


        case (excit_type%fullstop_lowering, excit_type%fullstart_raising)

            integral = abs(get_umat_el(a, a, i, j) + get_umat_el(a, a, j, i))/2.0_dp

        case (excit_type%fullstop_raising, excit_type%fullstart_lowering)

            integral = abs(get_umat_el(a, b, i, i) + get_umat_el(b, a, i, i))/2.0_dp

        case (excit_type%fullstop_R_to_L, excit_type%fullstop_L_to_R)

            integral = abs(get_umat_el(b, a, j, b) + get_umat_el(a, b, b, j))/2.0_dp

        case (excit_type%fullstart_L_to_R, excit_type%fullstart_R_to_L)

            integral = abs(get_umat_el(a, b, i, a) + get_umat_el(b, a, a, i))/2.0_dp

        case (excit_type%fullstart_stop_alike)

            integral = abs(get_umat_el(a, a, i, i))/2.0_dp

        case (excit_type%fullstart_stop_mixed)

            integral = abs(get_umat_el(a, i, i, a) + get_umat_el(i, a, a, i)) / 2.0_dp

#ifdef DEBUG_
        case default
            call stop_all(this_routine, "wrong excit-type")
#endif

        end select

    end function get_pchb_integral_contrib


    subroutine store_pchb_analysis(h_element, pgen, excitInfo, not_valid)
        real(dp), intent(in) :: h_element, pgen
        type(ExcitationInformation_t), intent(in) :: excitInfo
        logical, intent(in) :: not_valid
        integer :: ij, ab


        associate(a => excitInfo%i, i => excitInfo%j, b => excitInfo%k, &
                  j => excitInfo%l)


            select case (excitInfo%typ)

            ! all the cases without recalculating
            case (excit_type%single_overlap_L_to_R, &
                  excit_type%single_overlap_R_to_L, &
                  excit_type%double_lowering, &
                  excit_type%double_raising, &
                  excit_type%double_L_to_R_to_L, &
                  excit_type%double_R_to_L_to_R, &
                  excit_type%double_L_to_R, &
                  excit_type%double_R_to_L, &
                  excit_type%fullstop_lowering, &
                  excit_type%fullstop_raising, &
                  excit_type%fullstart_lowering, &
                  excit_type%fullstart_raising, &
                  excit_type%fullstart_stop_alike, &
                  excit_type%fullstop_R_to_L, &
                  excit_type%fullstop_L_to_R, &
                  excit_type%fullstart_L_to_R, &
                  excit_type%fullstart_R_to_L, &
                  excit_type%fullstart_stop_mixed)

                ij = fuseIndex(i,j)
                ab = fuseIndex(a,b)

                if (not_valid) then
                    call guga_pchb_sampler%setup_entry_invalid_scalar(ij,ab)
                else
                    call guga_pchb_sampler%setup_entry_count_scalar(ij,ab)
                    call guga_pchb_sampler%setup_entry_sum_scalar(ij,ab,abs(h_element))
                    call guga_pchb_sampler%setup_entry_worst_orb_scalar(ij,ab, &
                        abs(h_element) / guga_pchb_sampler%alias_sampler%get_prob(ij,ab))
                    call guga_pchb_sampler%setup_entry_pgen_scalar(ij,ab,pgen)
                    call guga_pchb_sampler%setup_entry_pgen_high_scalar(ij,ab, &
                        abs(h_element) / pgen)
                    call guga_pchb_sampler%setup_entry_pgen_low_scalar(ij,ab, &
                        abs(h_element) / pgen)
                end if

               ! in the other cases no contribution to PCHB
            end select

        end associate

    end subroutine store_pchb_analysis

! **************** analysis functions (to be removed after optimization) ******

    subroutine setup_invalid_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%invalid_tables(nEntries))

        total_size = nEntries * entrySize

        call this%all_invalids%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%invalid_tables(iEntry)%ptr => this%all_invalids%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_invalid_table

    subroutine setup_count_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%count_tables(nEntries))

        total_size = nEntries * entrySize

        call this%all_counts%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%count_tables(iEntry)%ptr => this%all_counts%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_count_table

    subroutine setup_sum_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%sums_tables(nEntries))

        total_size = nEntries * entrySize

        call this%all_sums%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%sums_tables(iEntry)%ptr => this%all_sums%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_sum_table

    subroutine setup_pgen_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%pgen_table(nEntries))

        total_size = nEntries * entrySize

        call this%all_pgen%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%pgen_table(iEntry)%ptr => this%all_pgen%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_pgen_table


    subroutine setup_worst_orb_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%worst_orb_table(nEntries))

        total_size = nEntries * entrySize

        call this%all_worst_orb%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%worst_orb_table(iEntry)%ptr => this%all_worst_orb%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_worst_orb_table

    subroutine setup_high_pgen_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%high_pgen_table(nEntries))

        total_size = nEntries * entrySize

        call this%all_high_pgen%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%high_pgen_table(iEntry)%ptr => this%all_high_pgen%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_high_pgen_table

    subroutine setup_low_pgen_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%low_pgen_table(nEntries))

        total_size = nEntries * entrySize

        call this%all_low_pgen%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%low_pgen_table(iEntry)%ptr => this%all_low_pgen%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_low_pgen_table

    subroutine setup_entry_invalid_vec(this, iEntry, invalids)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        integer(int64), intent(in) :: invalids(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%invalid_tables(iEntry)%ptr = invalids
        end if

        ! then sync:
        call this%all_invalids%sync()

    end subroutine setup_entry_invalid_vec

    subroutine setup_entry_invalid_scalar(this, iEntry, tgt)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%invalid_tables(iEntry)%ptr(tgt) = this%get_invalid(iEntry,tgt) &
                + 1_int64
        end if

        ! then sync:
        call this%all_invalids%sync()

    end subroutine setup_entry_invalid_scalar



    subroutine setup_entry_count_vec(this, iEntry, counts)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        integer(int64), intent(in) :: counts(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%count_tables(iEntry)%ptr = counts
        end if

        ! then sync:
        call this%all_counts%sync()

    end subroutine setup_entry_count_vec

    subroutine setup_entry_count_scalar(this, iEntry, tgt)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%count_tables(iEntry)%ptr(tgt) = this%get_count(iEntry,tgt) + 1_int64
        end if

        ! then sync:
        call this%all_counts%sync()

    end subroutine setup_entry_count_scalar

    subroutine setup_entry_sum_vec(this, iEntry, sums)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        real(dp), intent(in) :: sums(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%sums_tables(iEntry)%ptr = sums
        end if

        ! then sync:
        call this%all_sums%sync()

    end subroutine setup_entry_sum_vec

    subroutine setup_entry_pgen_vec(this, iEntry, pgens)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        real(dp), intent(in) :: pgens(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%pgen_table(iEntry)%ptr = pgens
        end if

        ! then sync:
        call this%all_pgen%sync()

    end subroutine setup_entry_pgen_vec


    subroutine setup_entry_sum_scalar(this, iEntry, tgt, sigma)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp), intent(in) :: sigma

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%sums_tables(iEntry)%ptr(tgt) = this%get_sum(iEntry,tgt) + sigma
        end if

        ! then sync:
        call this%all_sums%sync()

    end subroutine setup_entry_sum_scalar

    subroutine setup_entry_pgen_scalar(this, iEntry, tgt, pgen)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp), intent(in) :: pgen

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%pgen_table(iEntry)%ptr(tgt) = this%get_pgen(iEntry,tgt) + pgen
        end if

        ! then sync:
        call this%all_pgen%sync()

    end subroutine setup_entry_pgen_scalar

    subroutine setup_entry_worst_orb_vec(this, iEntry, worst)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        real(dp), intent(in) :: worst(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%worst_orb_table(iEntry)%ptr = worst
        end if

        ! then sync:
        call this%all_worst_orb%sync()
    end subroutine setup_entry_worst_orb_vec


    subroutine setup_entry_pgen_low_vec(this, iEntry, worst)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        real(dp), intent(in) :: worst(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%low_pgen_table(iEntry)%ptr = worst
        end if

        ! then sync:
        call this%all_low_pgen%sync()
    end subroutine setup_entry_pgen_low_vec



    subroutine setup_entry_pgen_high_vec(this, iEntry, worst)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        real(dp), intent(in) :: worst(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%high_pgen_table(iEntry)%ptr = worst
        end if

        ! then sync:
        call this%all_high_pgen%sync()
    end subroutine setup_entry_pgen_high_vec

    subroutine setup_entry_worst_orb_scalar(this, iEntry, tgt, worst)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp), intent(in) :: worst

        ! i think this is everyhing...
        ! do i need to change this iProcIndex_intra if? i guess so..
        if (iProcIndex_intra == 0) then
            this%worst_orb_table(iEntry)%ptr(tgt) = max(worst, this%get_worst_orb(iEntry,tgt))
        end if

        ! then sync:
        call this%all_worst_orb%sync()
    end subroutine setup_entry_worst_orb_scalar

    subroutine setup_entry_pgen_low_scalar(this, iEntry, tgt, worst)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp), intent(in) :: worst

        real(dp) :: old

        if (near_zero(worst)) return

        ! i think this is everyhing...
        ! do i need to change this iProcIndex_intra if? i guess so..
        if (iProcIndex_intra == 0) then
            old = this%get_low_pgen(iEntry,tgt)
            if (near_zero(old)) then
                this%low_pgen_table(iEntry)%ptr(tgt) = worst
            else
                this%low_pgen_table(iEntry)%ptr(tgt) = min(worst, old)
            end if
        end if

        ! then sync:
        call this%all_low_pgen%sync()

    end subroutine setup_entry_pgen_low_scalar


    subroutine setup_entry_pgen_high_scalar(this, iEntry, tgt, worst)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp), intent(in) :: worst

        ! i think this is everyhing...
        ! do i need to change this iProcIndex_intra if? i guess so..
        if (iProcIndex_intra == 0) then
            this%high_pgen_table(iEntry)%ptr(tgt) = max(worst, this%get_high_pgen(iEntry,tgt))
        end if

        ! then sync:
        call this%all_high_pgen%sync()

    end subroutine setup_entry_pgen_high_scalar


    subroutine invalid_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_invalids%shared_dealloc()

    end subroutine invalid_table_destructor


    subroutine count_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_counts%shared_dealloc()

    end subroutine count_table_destructor

    subroutine sum_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_sums%shared_dealloc()

    end subroutine sum_table_destructor

    subroutine pgen_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_pgen%shared_dealloc()

    end subroutine pgen_table_destructor


    subroutine worst_orb_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_worst_orb%shared_dealloc()

    end subroutine worst_orb_table_destructor

    subroutine high_pgen_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_high_pgen%shared_dealloc()

    end subroutine high_pgen_table_destructor

    subroutine low_pgen_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_low_pgen%shared_dealloc()

    end subroutine low_pgen_table_destructor



    function get_invalid(this, iEntry, tgt) result(cnt)
        debug_function_name("get_invalid")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: cnt

        ASSERT(associated(this%invalid_tables(iEntry)%ptr))

        cnt = this%invalid_tables(iEntry)%ptr(tgt)

    end function get_invalid


    function get_count(this, iEntry, tgt) result(cnt)
        debug_function_name("get_count")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: cnt

        ASSERT(associated(this%count_tables(iEntry)%ptr))

        cnt = this%count_tables(iEntry)%ptr(tgt)

    end function get_count

    function get_sum(this, iEntry, tgt) result(sigma)
        debug_function_name("get_sum")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp) :: sigma

        ASSERT(associated(this%sums_tables(iEntry)%ptr))

        sigma = this%sums_tables(iEntry)%ptr(tgt)

    end function get_sum

    function get_pgen(this, iEntry, tgt) result(pgen)
        debug_function_name("get_pgen")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp) :: pgen

        ASSERT(associated(this%pgen_table(iEntry)%ptr))

        pgen = this%pgen_table(iEntry)%ptr(tgt)

    end function get_pgen

    function get_worst_orb(this, iEntry, tgt) result(worst)
        debug_function_name("get_worst_orb")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp) :: worst

        ASSERT(associated(this%worst_orb_table(iEntry)%ptr))

        worst = this%worst_orb_table(iEntry)%ptr(tgt)

    end function get_worst_orb

    function get_high_pgen(this, iEntry, tgt) result(worst)
        debug_function_name("get_high_pgen")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp) :: worst

        ASSERT(associated(this%high_pgen_table(iEntry)%ptr))

        worst = this%high_pgen_table(iEntry)%ptr(tgt)

    end function get_high_pgen

    function get_low_pgen(this, iEntry, tgt) result(worst)
        debug_function_name("get_low_pgen")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        real(dp) :: worst

        ASSERT(associated(this%low_pgen_table(iEntry)%ptr))

        worst = this%low_pgen_table(iEntry)%ptr(tgt)

    end function get_low_pgen




    subroutine print_pchb_statistics
        ! routine to print out the accumulated PCHB excit-gen statistics
        integer :: i, j, a, b, iunit, ij, ab, dist, overlap
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: weight, sums, worst_orb, ratio_orb, pgen_sum, high_pgen, &
                    ratio_pgen, low_pgen
        integer(int64) :: counts, invalids
        integer(int_rdm) :: ijkl

        ! can i do this on the root? do i have to accumulate over all nodes??
        ! i guess so.. for now test on 1 node..
        if (iProcIndex_intra == 0) then

            iunit = get_free_unit()
            open(iunit, file = 'pchb-stats', status = 'unknown')

            write(iunit, '(a)') &
                '# E_{a,i} E_{b,j}, rdm-ind, dist, overlap, weight, counts, sums, worst_orb-case, &
                &sums/(counts*weights), pgen-sum, high-pgen, low-pgen, sums/pgens, typ'

            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    ij = fuseIndex(i,j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a,b)

                            weight = guga_pchb_sampler%alias_sampler%get_prob(ij,ab)

                            if (.not. near_zero(weight)) then
                                call extract_excit_info(guga_pchb_sampler%get_info(ij, ab), &
                                    excitInfo)

                                counts = guga_pchb_sampler%get_count(ij,ab)
                                invalids = guga_pchb_sampler%get_invalid(ij,ab)
                                if (counts > 0 .or. invalids > 0) then
                                    sums = guga_pchb_sampler%get_sum(ij,ab)
                                    worst_orb = guga_pchb_sampler%get_worst_orb(ij,ab)
                                    pgen_sum = guga_pchb_sampler%get_pgen(ij,ab)
                                    high_pgen = guga_pchb_sampler%get_high_pgen(ij,ab)
                                    low_pgen = guga_pchb_sampler%get_low_pgen(ij,ab)
                                    dist = excitInfo%fullEnd - excitInfo%fullstart
                                    overlap = excitInfo%firstEnd - excitInfo%secondStart
                                    ijkl = contract_2_rdm_ind(excitInfo%i, &
                                        excitInfo%j, excitInfo%k, excitInfo%l)


                                    if (counts > 0_int64) then
                                        ratio_orb = sums / (real(counts,dp) * weight)
                                    else
                                        ratio_orb = 0.0_dp
                                    end if

                                    if (.not. near_zero(pgen_sum)) then
                                        ratio_pgen = sums / pgen_sum
                                    else
                                        ratio_pgen = 0.0_dp
                                    end if

                                    write(iunit, '(4I3,I12,2I4,E15.8,I12,7E15.8,I3,I12)') &
                                        excitInfo%i, excitInfo%j, &
                                        excitInfo%k, excitInfo%l, ijkl, dist, overlap, weight, counts, &
                                        sums, worst_orb, ratio_orb, pgen_sum, &
                                        high_pgen, low_pgen, ratio_pgen, excitInfo%typ, &
                                        invalids
                                end if
                            end if
                        end do
                    end do
                end do
            end do

            close(iunit)
        end if

    end subroutine print_pchb_statistics

! *************END  analysis functions (to be removed after optimization) ******

    subroutine setup_info_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%info_tables(nEntries))

        total_size = nEntries * entrySize

        call this%all_info_table%shared_alloc(total_size)

        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%info_tables(iEntry)%ptr => this%all_info_table%ptr(windowStart:windowEnd)
        end do

    end subroutine setup_info_table

    function get_info(this, iEntry, tgt) result(info)
        debug_function_name("get_info")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: info

        ASSERT(associated(this%info_tables(iEntry)%ptr))

        info = this%info_tables(iEntry)%ptr(tgt)

    end function get_info

    subroutine info_table_destructor(this)
        class(GugaAliasSampler_t), intent(inout) :: this

        call this%all_info_table%shared_dealloc()

    end subroutine info_table_destructor

    subroutine setup_entry_info(this, iEntry, infos)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        integer(int64), intent(in) :: infos(:)

        ! i think this is everyhing...
        if (iProcIndex_intra == 0) then
            this%info_tables(iEntry)%ptr = infos
        end if

        ! then sync:
        call this%all_info_table%sync()

    end subroutine setup_entry_info

    subroutine init_guga_pchb_excitgen
        integer :: a, b
        integer(int64) :: memCost, ijMax, abMax, ab

        ! also set some more strict defaults for the PCHB implo:
        root_print "Setting reasonable defaults for GUGA-PCHB:"
        if (near_zero(MaxTau) .or. MaxTau > 1e-3) then
            root_print "max-tau zero or > 1e-3. setting it to: 1e-3"
            MaxTau = 1e-3
        end if

        if (t_hist_tau_search) then
            if (frq_ratio_cutoff < 0.999999) then
                root_print "setting frequency cutoff to 0.999999"
                frq_ratio_cutoff = 0.999999
            end if
            if (max_frequency_bound < 1e5) then
                root_print "setting  max_frequency_bound to 1e5"
                max_frequency_bound = 1e5
            end if
            if (n_frequency_bins < 1e5) then
                root_print "setting n_frequency_bins to 1e5"
                n_frequency_bins = 1e5
            end if
        end if

        if (.not. t_truncate_spawns) then
            root_print "'truncate-spawns' not activated! consider turning it &
                &on if too many blooms happen!"
        end if

        write(iout,*) "Allocating GUGA PCHB excitation generator objects"
        ! initialize the mapping ab -> (a,b)
        abMax = fuseIndex(nSpatOrbs,nSpatOrbs)
        allocate(tgtOrbs(2,0:abMax), source = 0)
        do a = 1, nSpatOrbs
            do b = a, nSpatOrbs
                ab = fuseIndex(a,b)
                tgtOrbs(1,ab) = a
                tgtOrbs(2,ab) = b
            end do
        end do

        ! enable catching exceptions
        tgtOrbs(:,0) = 0

        ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        memCost = abMax*ijMax*24*2
        write(iout,*) "Excitation generator requires", &
            real(memCost,dp)/2.0_dp**30, "GB of memory"
        write(iout,*) "Generating samplers for PCHB excitation generator"

        call setup_pchb_sampler_conditional()

        write(iout,*) "Finished GUGA PCHB excitation generator initialization"

    end subroutine init_guga_pchb_excitgen

    subroutine setup_pchb_sampler_conditional()
        integer :: i, j, ij, a, b, ab
        integer(int64) :: ijMax, abMax
        integer(int64), allocatable :: excit_info(:), counts(:)
        real(dp), allocatable :: w(:), x(:)

        ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        abMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        ! just to be save also code up a setup function with the
        ! necessary if statements, so i can check if my optimized
        ! sampler is correct (and if it is even faster..)
        ! weights per pair
        allocate(w(abMax), source = 0.0_dp)
        ! excitation information per pair:
        allocate(excit_info(abMax))

        ! allocate: all samplers have the same size
        call guga_pchb_sampler%alias_sampler%shared_alloc(int(ijMax), int(abMax), "GUGA-pchb")
        ! todo: do the same for the excit_info array!
        call guga_pchb_sampler%setup_info_table(ijMax, abMax)

        if (t_analyze_pchb) then
            ! also for the analysis for now..
            call guga_pchb_sampler%setup_count_table(ijMax, abMax)
            call guga_pchb_sampler%setup_invalid_table(ijMax, abMax)
            call guga_pchb_sampler%setup_sum_table(ijMax, abMax)
            call guga_pchb_sampler%setup_worst_orb_table(ijMax, abMax)
            call guga_pchb_sampler%setup_high_pgen_table(ijMax, abMax)
            call guga_pchb_sampler%setup_low_pgen_table(ijMax, abMax)
            call guga_pchb_sampler%setup_pgen_table(ijMax, abMax)

            allocate(x(abMax), source = 0.0_dp)
            allocate(counts(abMax), source = 0_int64)
        end if

        ! the encode_excit_info function encodes as: E_{ai}E_{bj)
        ! but for some index combinations we have to change the input so
        ! actually E_{aj}E_{bi} is encoded! convention to use only
        ! certain type of excit-types
        do i = 1, nSpatOrbs
            do j = i, nSpatOrbs
                w = 0.0_dp
                ij = fuseIndex(i,j)
                excit_info = 0_int64
                do a = 1, nSpatOrbs
                    do b = a, nSpatOrbs
                        ! shoud i point group symmetry restrctions here?
                        ! this would avoid unnecessary if statements..

                        if (RandExcitSymLabelProd(SpinOrbSymLabel(2*i), &
                            SpinOrbSymLabel(2*j)) /= &
                            RandExcitSymLabelProd(SpinOrbSymLabel(2*a), &
                            SpinOrbSymLabel(2*b))) cycle

                        ab = fuseIndex(a,b)
                        if (i == j) then
                            if (a == b) then
                                ! here we only have a contribution if
                                ! a != i
                                if (a < i) then
                                    ! _RR_(a) -> ^RR^(i)
                                    w(ab) = get_pchb_integral_contrib(i, j, a, b, &
                                        typ=excit_type%fullstart_stop_alike)
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%fullstart_stop_alike, &
                                        a=a, i=i, b=b, j=j)
                                else if (a > i) then
                                    ! _LL_(i) > ^LL^(a)
                                    w(ab) = get_pchb_integral_contrib(i, j, a, b, &
                                        typ = excit_type%fullstart_stop_alike)
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%fullstart_stop_alike, &
                                        a=a, i=i, b=b, j=j)
                                end if
                            elseif (a /= b) then
                                ! here we have to determine where (a) and
                                ! (b) are relative to (i=j) and a == i or
                                ! b == i are NOT allowed!
                                if (a < i) then
                                    if (b < i) then
                                        ! _R(a) -> _RR(b) -> ^RR^(i)
                                        w(ab) = get_pchb_integral_contrib(i, j, a, b, &
                                            typ=excit_type%fullstop_raising)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%fullstop_raising, &
                                            a=a, i=i, b=b, j=j)
                                    else if (b > i) then
                                        ! _R(a) -> ^RL_(i) -> ^L(b)
                                        w(ab) = get_pchb_integral_contrib(i, j, a, b, &
                                            excit_type%single_overlap_R_to_L)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%single_overlap_R_to_L, &
                                            a=a, i=i, b=b, j=j)
                                    end if
                                else if (a > i) then
                                    ! since b > a ensured only:
                                    ! _LL_(i) -> ^LL(a) -> ^L(b)
                                    w(ab) = get_pchb_integral_contrib(i, j, a, b, &
                                        excit_type%fullstart_lowering)
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%fullstart_lowering, &
                                        a=a, i=i, b=b, j=j)
                                end if
                            end if
                        else if ( i /= j) then
                            if (a == b) then
                                ! a == i or a == j NOT allowed!
                                if (a < i) then
                                    ! _RR_(a) -> ^RR(i) -> ^R(j)
                                    w(ab) = get_pchb_integral_contrib(i, j, a, b,&
                                        typ = excit_type%fullstart_raising)
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%fullstart_raising, &
                                        a=a, i=i, b=b, j=j)
                                else if (a > i .and. a < j) then
                                    ! _L(i) -> ^LR_(a) -> ^R(j)
                                    w(ab) = get_pchb_integral_contrib(i, j, a, b,&
                                        excit_type%single_overlap_L_to_R)
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%single_overlap_L_to_R, &
                                        a=a, i=i, b=b, j=j)
                                else if (a > j) then
                                    ! _L(i) -> _LL(j) -> ^LL^(a)
                                    w(ab) = get_pchb_integral_contrib(i, j, a, b,&
                                        excit_type%fullstop_lowering)
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%fullstop_lowering, &
                                        a=a, i=i, b=b, j=j)
                                end if
                            else if (a /= b) then
                                ! this is the most general case. a lot of IFs
                                if (a < i) then
                                    if (b < i) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=b, b=a, &
                                            typ=excit_type%double_raising)
                                        ! in E_{ai}E_{bi} form this would be
                                        ! _R(a) -> RR_(b) -> ^RR(i) -> ^R(j)
                                        ! which would have the opposite
                                        ! sign conventions for the x1
                                        ! elements as in the Shavitt 81
                                        ! paper.
                                        ! (technically it would not matter
                                        !  here since both are flipped,
                                        !  but lets stay consistent!)
                                        ! for E_{bj}E_{ai}
                                        ! _R(a) -> _RR(b) -> RR^(i) -> ^R(j)
                                        ! it would fit! so:
                                        excit_info(ab) = encode_excit_info(&
                                            typ = excit_type%double_raising, &
                                            a=b, i=j, b=a, j=i)
                                    ! else if (b == i) then
                                        ! b == i also NOT allowed here,
                                        ! since this would correspond to
                                        ! a single!
                                    else if (b > i .and. b < j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%double_R_to_L_to_R)
                                        ! for E_{ai}E_{bj} this would
                                        ! correspond to a non-overlap:
                                        ! _R(a) -> ^R(i) + _R(b) > ^R(j)
                                        ! which are not directly sampled
                                        ! However for E_{aj}E_{bj} this is
                                        ! _R{a} -> _LR(i) -> ^LR(b) -> ^R(j)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%double_R_to_L_to_R, &
                                            a=a, i=j, b=b, j=i)
                                    else if (b == j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%fullstop_R_to_L)
                                        ! here we also have to switch
                                        ! indices to E_{aj}E_{bi} to get:
                                        ! _R(i) -> _LR(a) -> ^LR^(j)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%fullstop_R_to_L, &
                                            a=a, i=j, b=b, j=i)

                                    else if (b > j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%double_R_to_L)
                                        ! here we have to switch to
                                        ! E_{aj}E_{bi} to get:
                                        ! _R(a) > _LR(i) -> LR^(j) -> ^L(b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%double_R_to_L, &
                                            a=a, i=j, b=b, j=i)
                                    end if
                                else if (a == i) then
                                    ! b > i is ensured here since b > a in here!
                                    if (b < j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%fullstart_L_to_R)
                                        ! here we have to switch again:
                                        ! E_{aj}E_{bi}:
                                        ! _RL_(i) -> ^LR(b) -> ^R(j)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%fullstart_L_to_R, &
                                            a=a, i=j, b=b, j=i)
                                    else if (b == j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%fullstart_stop_mixed)
                                        ! switch again: E_{aj}E_{bi}
                                        ! _RL_(i) -> _RL_(j)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%fullstart_stop_mixed, &
                                            a=a, i=j, b=b, j=i)
                                    else if (b > j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%fullstart_R_to_L)
                                        ! switch again: E_{aj}E_{bi}
                                        ! _RL_(i) -> ^RL(j) -> ^L(b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%fullstart_R_to_L, &
                                            a=a, i=j, b=b, j=i)
                                    end if
                                else if (a > i .and. a < j) then
                                    ! b > a still ensured!
                                    if (b < j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%double_L_to_R)
                                        ! switch: E_{aj}E_{bi}:
                                        ! _L(j) -> _RL(a) -> ^LR(b) -> ^R(j)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%double_L_to_R, &
                                            a=a, i=j, b=b, j=i)

                                    else if (b == j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%fullstop_L_to_R)
                                        ! switch: E_{aj}E_{bi}
                                        ! _L(i) -> _RL(a) -> ^RL^(j)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%fullstop_L_to_R, &
                                            a=a, i=j, b=b, j=i)
                                    else if (b > j) then
                                        w(ab) = get_pchb_integral_contrib(&
                                            i=j, j=i, a=a, b=b, &
                                            typ=excit_type%double_L_to_R_to_L)
                                        ! switch: E_{aj}E_{bi}
                                        ! _L(i) -> _RL(a) - > ^RL(j) -> ^L(b)
                                        excit_info(ab) = encode_excit_info(&
                                            typ=excit_type%double_L_to_R_to_L, &
                                            a=a, i=j, b=b, j=i)
                                    end if
                                ! else if (a == j) then
                                    ! a == j also NOT allowed here!
                                else if (a > j) then
                                    ! b > a > j implied here!
                                    w(ab) = get_pchb_integral_contrib(&
                                        i=i, j=j, a=a, b=b, &
                                        typ=excit_type%double_lowering)
                                    ! E_{ai}E_{bj} would lead to:
                                    ! _L(i) -> LL_(j) -> ^LL(a) -> ^L(b)
                                    ! which has the correct sign convention
                                    ! as in the Shavitt 81 paper
                                    excit_info(ab) = encode_excit_info(&
                                        typ=excit_type%double_lowering, &
                                        a=a, i=i, b=b, j=j)
                                end if
                            end if
                        end if
                    end do
                end do
                call guga_pchb_sampler%alias_sampler%setup_entry(ij,w)
                call guga_pchb_sampler%setup_entry_info(ij,excit_info)

                if (t_analyze_pchb) then
                    ! for analysis
                    call guga_pchb_sampler%setup_entry_sum_vec(ij,x)
                    call guga_pchb_sampler%setup_entry_count_vec(ij,counts)
                    call guga_pchb_sampler%setup_entry_invalid_vec(ij,counts)
                    call guga_pchb_sampler%setup_entry_pgen_vec(ij,x)
                    call guga_pchb_sampler%setup_entry_worst_orb_vec(ij,x)
                    call guga_pchb_sampler%setup_entry_pgen_high_vec(ij,x)
                    call guga_pchb_sampler%setup_entry_pgen_low_vec(ij,x)
                end if

                ! todo: do the same for the excit_info array!
            end do
        end do

    end subroutine setup_pchb_sampler_conditional

    subroutine finalize_pchb_excitgen_guga()

        if (t_analyze_pchb) then
            ! for analysis:
            call print_pchb_statistics()
            call guga_pchb_sampler%count_table_destructor()
            call guga_pchb_sampler%invalid_table_destructor()
            call guga_pchb_sampler%sum_table_destructor()
            call guga_pchb_sampler%worst_orb_table_destructor()
            call guga_pchb_sampler%high_pgen_table_destructor()
            call guga_pchb_sampler%low_pgen_table_destructor()
            call guga_pchb_sampler%pgen_table_destructor()
        end if

        call guga_pchb_sampler%alias_sampler%finalize()
        call guga_pchb_sampler%info_table_destructor()
        deallocate(tgtOrbs)


    end subroutine finalize_pchb_excitgen_guga

    subroutine pick_orbitals_double_pchb(ilut, nI, excitInfo, pgen)
        debug_function_name("pick_orbitals_double_pchb")
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), sym_prod, sum_ml, i, j, a, b, orbs(2), ij, ab
        real(dp) :: pgen_elec, pgen_orbs

        unused_var(ilut)

        ! maybe I will also produce a weighted electron pickin in the
        ! GUGA formalism.. but for now pick them uniformly:
        call pick_elec_pair_uniform_guga(nI, src, sym_prod, sum_ml, pgen_elec)
        ! make a different picker, which does not bias towards doubly
        ! occupied orbitals here too!

        ASSERT( src(1) < src(2) )

        i = gtID(src(1))
        j = gtID(src(2))

        if (i /= j) then
            if (current_stepvector(i) == 3) pgen_elec = 2.0_dp * pgen_elec
            if (current_stepvector(j) == 3) pgen_elec = 2.0_dp * pgen_elec
        end if

        ! use the sampler for this electron pair -> order of src electrons
        ! does not matter
        ij = fuseIndex(i, j)

        ! get a pair of orbitals using the precomputed weights
        call guga_pchb_sampler%alias_sampler%sample(ij,ab,pgen_orbs)

        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if(near_zero(pgen_orbs)) then
            excitInfo%valid = .false.
            pgen = 0.0_dp
            ! Yes, print. Those events are signficant enough to be always noted in the output
            print *, "WARNING: Generated excitation with probability of 0"
            return
        endif


        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = tgtOrbs(:,ab)

        a = orbs(1)
        b = orbs(2)

        ! check if the picked orbs are a valid choice - if they are the same, match one
        ! occupied orbital or are zero (maybe because there are no allowed picks for
        ! the given source) abort
        ! these are the 'easy' checks for GUGA.. more checks need to be done
        ! to see if it is actually a valid combination..
        if (any(orbs == 0) .or. (current_stepvector(a) == 3) .or. &
            (current_stepvector(b) == 3) .or. (a == b .and. &
            current_stepvector(b) /= 0)) then
            excitInfo%valid = .false.
            return
        end if

        ! setup a getInfo functionality in the sampler!
        call extract_excit_info(guga_pchb_sampler%get_info(ij, ab), excitInfo)

        pGen = pgen_elec * pgen_orbs

    end subroutine pick_orbitals_double_pchb

    subroutine pick_orbitals_pure_uniform_singles(ilut, nI, excitInfo, pgen)
        debug_function_name("pick_orbitals_pure_uniform_singles")
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: elec, orb, nOcc
        integer :: s_elec(nEl)

        unused_var(ilut)

        ASSERT(isProperCSF_ilut(ilut, .true.))

        ! init for safety and so we can return on abort
        pgen = 0.0_dp

        ! pick random electron
        ! have to modify pgen for doubly occupied orbs! since double the
        elec = 1 + int(genrand_real2_dSFMT() * nel)

        ! chance!
        ! i have this fake bias towards double occupied spatial orbitals
        ! this way...
        ! try it with a 'real pure' uniform picking here..
        s_elec = gtID(nI)

        ! r = genrand_real2_dSFMT() * s_elec(nel)
        ! elec = binary_search_first_ge(real(s_elec,dp), r)
        elec = s_elec(elec)
        nOcc = count(currentOcc_int /= 0)

        call pick_uniform_spatial_hole(elec, orb, pgen)

        if (near_zero(pgen) .or. orb == 0) then
            pgen = 0.0_dp
            orb = 0
            return
        end if

        ! pgen = pgen / real(nOcc, dp)
        pgen = pgen / real(nel, dp)
        if (current_stepvector(elec) == 3) pgen = 2.0_dp * pgen

        if (orb < elec) then
            excitInfo = assign_excitinfo_values_single(gen_type%R, orb, elec, &
                orb, elec)
        else
            excitInfo = assign_excitinfo_values_single(gen_type%L, orb, elec, &
                elec, orb)
        end if

    end subroutine pick_orbitals_pure_uniform_singles

    subroutine pick_uniform_spatial_hole(s_elec, s_orb, pgen)
        integer, intent(in) :: s_elec
        integer, intent(out) :: s_orb
        real(dp), intent(out) :: pgen

        integer :: cc_i, nOrb, sym_index, orb, nValid
        integer, allocatable :: sym_orbs(:), valid_orbs(:)
        logical, allocatable :: mask(:)

        pgen = 0.0_dp
        s_orb = 0

        ! get the symmetry index for this electron
        cc_i = ClassCountInd(1, SpinOrbSymLabel(2*s_elec), G1(2*s_elec)%ml)

        ! and the number of orbitals
        nOrb = OrbClassCount(cc_i)

        if (nOrb == 0) return

        ! get the symmetry index for later use
        sym_index = SymLabelCounts2(1, cc_i)

        ! now keep drawing from the symmetry orbitals until we pick an
        ! 'empty' (not doubly occupied!) one

        ! or better determine the non-double occupied orbitals in this
        ! symmetry sector
        ! and also is not the electron index!
        allocate(sym_orbs(nOrb), source = sym_label_list_spat(sym_index:sym_index+nOrb-1))
        allocate(mask(nOrb), &
            source = ((current_stepvector(sym_orbs) /= 3) .and. &
                     (sym_orbs /= s_elec)))

        nValid = count(mask)

        allocate(valid_orbs(nValid), source = pack(sym_orbs, mask))

        if (nValid == 0) return

        orb = 1 + floor(genrand_real2_dSFMT() * nValid)

        s_orb = valid_orbs(orb)

        ! do i need nOrb now or the actual number of unoccupied in the
        ! CSF? i think the second..
        pgen = 1.0_dp / real(nValid, dp)


    end subroutine pick_uniform_spatial_hole

    function calc_orb_pgen_uniform_singles_exmat(ex) result(pgen)
        integer, intent(in) :: ex(2,2)
        real(dp) :: pgen
        debug_function_name("calc_orb_pgen_uniform_singles_exmat")

        integer :: src(2), so_elec, cc_i, nOrb, nOcc, sym_index, nUnocc

        ASSERT(all(ex(:,2) == 0))
        ASSERT(all(ex(:,1) > 0))
        ASSERT(all(ex(:,1) <= nBasis))

        pgen = 0.0_dp

        if (gtID(ex(1,1)) == gtID(ex(2,1))) return
        if (current_stepvector(gtID(ex(1,1))) == 0) return
        if (current_stepvector(gtID(ex(2,1))) == 3) return

        src = get_src(ex)

        ! determine the number of occupied spatial orbs:
        nOcc = count(currentOcc_int /= 0)

        ! find the number of non-double occupied symmetry allowed orbitals /= s_elec
        so_elec = src(1)
        cc_i = ClassCountInd(1, SpinOrbSymLabel(so_elec), G1(so_elec)%ml)
        nOrb = OrbClassCount(cc_i)
        ! get the symmetry index for later use
        sym_index = SymLabelCounts2(1, cc_i)

        nUnocc = count(&
            (current_stepvector(sym_label_list_spat(sym_index:sym_index+nOrb-1)) /= 3) &
            .and. (sym_label_list_spat(sym_index:sym_index+nOrb-1) /= gtID(so_elec)))


        pgen = 1.0_dp / real(nOcc * nUnocc, dp)

    end function calc_orb_pgen_uniform_singles_exmat

    function calc_orb_pgen_uniform_singles_excitInfo(excitInfo) result(pgen)
        debug_function_name("calc_orb_pgen_uniform_singles_excitInfo")
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: nOrb, so_elec, cc_i, nOcc, sym_index, nUnocc

        ASSERT(excitInfo%i > 0 .and. excitInfo%i <= nSpatOrbs)
        ASSERT(excitInfo%j > 0 .and. excitInfo%j <= nSpatOrbs)

        pgen = 0.0_dp

        if (excitInfo%i == excitInfo%j) return
        if (current_stepvector(excitInfo%i) == 3) return
        if (current_stepvector(excitInfo%j) == 0) return

        nOcc = count(currentOcc_int /= 0)

        so_elec = 2*excitInfo%j

        cc_i = ClassCountInd(1, SpinOrbSymLabel(so_elec), G1(so_elec)%ml)
        nOrb = OrbClassCount(cc_i)
        ! get the symmetry index for later use
        sym_index = SymLabelCounts2(1, cc_i)


        nUnocc = count(&
            (current_stepvector(sym_label_list_spat(sym_index:sym_index+nOrb-1)) /= 3) &
            .and. (sym_label_list_spat(sym_index:sym_index+nOrb-1) /= excitInfo%j))

        pgen = 1.0_dp / real(nOcc * nUnocc, dp)


    end function calc_orb_pgen_uniform_singles_excitInfo

    function calc_pgen_guga_pchb(ilutI, ilutJ, excitInfo_in) result(pgen)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(GugaBits%len_tot)
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_in
        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: pgen
        integer :: ic
        integer :: nI(nel), nJ(nel)

        if (present(excitInfo_in)) then
            excitInfo = excitInfo_in
        else
            excitInfo = identify_excitation(ilutI, ilutJ)
        end if

        ic = get_excit_level_from_excitInfo(excitInfo)

        if (ic == 1) then
            if (t_pchb_weighted_singles) then
                call decode_bit_det(nI, ilutI)
                call decode_bit_det(nJ, ilutJ)
                pgen = calc_pgen_mol_guga_single(ilutI, nI, ilutJ, &
                    nJ, excitInfo)
            else
                pgen = calc_orb_pgen_uniform_singles(excitInfo)
            end if

            pgen = pgen * pSingles

        else if (ic == 2) then

            pgen = pDoubles * calc_orb_pgen_guga_pchb_double(excitInfo)

        else
            pgen = 0.0_dp
        end if

    end function calc_pgen_guga_pchb

    function calc_orb_pgen_guga_pchb_double_excitInfo(excitInfo) result(pgen)
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: ij, ab
        real(dp) :: p_elec

        ! i, j, k and l entries have the info about the electrons and
        ! holes of the excitation: E_{ij}E_{kl} is the convention ...
        ij = fuseIndex(excitInfo%j, excitInfo%l)

        p_elec = 1.0_dp / real(ElecPairs, dp)

        ab = fuseIndex(excitInfo%i, excitInfo%k)
        pgen = p_elec * guga_pchb_sampler%alias_sampler%get_prob(ij, ab)

    end function calc_orb_pgen_guga_pchb_double_excitInfo

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    subroutine calc_orbital_pgen_contr_pchb(ilut, occ_orbs, cpt_a, cpt_b)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b

        integer :: ij
        unused_var(ilut)

        ! this function can in theory be called with both i < j and i > j..
        ! to take the correct values here!
        ! although for the fuseIndex function the order is irrelevant

        ! i think i have to consider both the above and below contribution..
        ! but i am not so sure how.. in this PCHB case..
        ! maybe in PCHB those 2 are just the same.. i am confused
        ij = fuseIndex(gtID(occ_orbs(1)), gtID(occ_orbs(2)))

        ! and both I and J are electron and hole indices here
        cpt_a = guga_pchb_sampler%alias_sampler%get_prob(ij, ij) / 2.0_dp
        ! but since they will get added in later routines i actually have
        ! do divide by 2 here..

        ! and in the PCHB there is no difference between the 2!
        cpt_b = cpt_a

    end subroutine calc_orbital_pgen_contr_pchb


    ! i think it would be better if i 'just' reimplement:
    subroutine calc_orbital_pgen_contr_start_pchb(occ_orbs, a, orb_pgen)
        debug_function_name("calc_orbital_pgen_contr_start_pchb")
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen

        integer :: i, j, ij, ab

        ! depending on type (R->L / L->R) a can be > j or < j, but always > i
        !
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ASSERT( i < a )

        ! here i is both electron and hole index!

        ij = fuseIndex(i, j)
        ab = fuseIndex(i, a)

        orb_pgen = guga_pchb_sampler%alias_sampler%get_prob(ij, ab)

    end subroutine calc_orbital_pgen_contr_start_pchb

    subroutine calc_orbital_pgen_contr_end_pchb(occ_orbs, a, orb_pgen)
        debug_function_name("calc_orbital_pgen_contr_end_pchb")
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen

        integer :: i, j, ij, ab

        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! and j is at the same time electron and hole index!

        ! depending on L->R or R->L type a can be > i o r < i, but always < j!
        ASSERT( a < j )

        ! j here is both elec and hole ind!
        ij = fuseIndex(i,j)

        ab = fuseIndex(a,j)

        orb_pgen = guga_pchb_sampler%alias_sampler%get_prob(ij,ab)

    end subroutine calc_orbital_pgen_contr_end_pchb

end module guga_pchb_excitgen

