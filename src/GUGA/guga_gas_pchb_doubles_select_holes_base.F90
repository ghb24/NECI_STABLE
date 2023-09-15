#include "macros.h"
module guga_gas_pchb_doubles_select_holes_base

    use aliasSampling, only: AliasSampler_2D_t
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm, bits_n_int
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel, G1, &
                          nBasis, nSpatOrbs, ElecPairs, &
                          t_analyze_pchb, t_old_pchb, t_exchange_pchb
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type, ExcitationType_t, &
        DistinctDoubleType_t, distinct_doubles, DistinctDouble_t, DistinctDoubleTypeVals_t, ijab_Index_t
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut, CSF_Info_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: near_zero, fuseIndex, binary_search_first_ge, &
                        get_free_unit, stop_all, operator(.isclose.), operator(.div.), &
                        EnumBase_t
    use CalcData, only: t_matele_cutoff, matele_cutoff, t_truncate_spawns
    use tau_main, only: max_tau, tau_search_method, possible_tau_search_methods
    use tau_search_hist, only: frq_ratio_cutoff, max_frequency_bound, n_frequency_bins
    use fortran_strings, only: to_upper
    use sym_general_mod, only: ClassCountInd
    use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel
    use excitation_types, only: Excite_1_t, Excite_2_t
    use UMatCache, only: gtID
    use guga_excitations, only: assign_excitinfo_values_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo, &
                                checkCompatibility_single
    use guga_bitrepops, only: identify_excitation, encode_excit_info, extract_excit_info, &
                              contract_2_rdm_ind
    use guga_base_class, only: DoublesGUGABase_t
    use bit_reps, only: decode_bit_det
    use shared_array, only: shared_array_int64_t, shared_array_real_t
    use MPI_wrapper, only: iProcIndex_intra, iprocindex, root
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use procedure_pointers, only: get_umat_el

    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use sets_mod, only: is_sorted, is_set

    use guga_gas_pchb_doubles_select_particles, only: ParticleSelector_t, &
        UniformParticles_t, PCHB_particle_selection_vals, &
        PCHB_ParticleSelection_t, particle_allocate_and_init => allocate_and_init


    better_implicit_none

    private
    public :: PCHB_HoleAlgorithm_t, PCHB_HoleAlgorithm_vals_t, &
        WeightingChoice_t, WeightingChoice_vals_t, GAS_PCHB_Doubles_t

    abstract interface
        real(dp) pure function get_weight_t(exc)
            import :: DistinctDouble_t, dp
            type(DistinctDouble_t), intent(in) :: exc
        end function
    end interface

    type, extends(EnumBase_t) :: WeightingChoice_t
        character(20) :: str
    end type

    type :: WeightingChoice_vals_t
        type(WeightingChoice_t) :: &
            INITIAL_WERNER = WeightingChoice_t(1, 'INITIAL-WERNER')
    contains
        procedure, nopass :: from_str => WeightingChoice_from_str
    end type

    type(WeightingChoice_vals_t), parameter :: weight_choice_vals = WeightingChoice_vals_t()

    type, extends(EnumBase_t) :: PCHB_HoleAlgorithm_t
        private
        character(9) :: str
    end type

    type :: PCHB_HoleAlgorithm_vals_t
        type(PCHB_HoleAlgorithm_t) :: &
            FAST_FAST = PCHB_HoleAlgorithm_t(1, 'FAST-FAST')
    contains
        procedure, nopass :: from_str => HoleAlgorithm_from_str
    end type

    type(PCHB_HoleAlgorithm_vals_t), parameter :: algorithm_vals = PCHB_HoleAlgorithm_vals_t()

    type, extends(DoublesGUGABase_t), abstract :: GAS_PCHB_Doubles_t
        private

        !> Use a lookup for the supergroup index in global_det_data
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.

        class(ParticleSelector_t), allocatable :: particle_selector

        type(shared_array_int64_t) :: all_info_table
        type(shared_array_int64_t), allocatable :: info_tables(:)

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(GASSpec_t), allocatable :: GAS_spec
        procedure(get_weight_t), pointer, nopass :: get_weight => null()

    contains
        private
        procedure, public :: init => init_GugaAliasSampler_t
        procedure, public :: finalize => finalize_GugaAliasSampler_t


        procedure :: new_info_table
        procedure :: set_info_entry
        procedure :: get_info_entry

        procedure(calculate_pchb_sampler_t), deferred :: setup_pchb_sampler_conditional
        procedure(additional_finalize_t), deferred :: finalize_pchb_sampler
    end type

    abstract interface
        subroutine calculate_pchb_sampler_t(this, ij_weights)
            import :: GAS_PCHB_Doubles_t, dp
            class(GAS_PCHB_Doubles_t), intent(inout) :: this
            real(dp), allocatable, intent(out) :: ij_weights(:, :, :)
        end subroutine

        subroutine additional_finalize_t(this)
            import :: GAS_PCHB_Doubles_t
            implicit none
            class(GAS_PCHB_Doubles_t), intent(inout) :: this
        end subroutine
    end interface



contains

    subroutine new_info_table(this, nEntries, entrySize)
        class(GAS_PCHB_Doubles_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%info_tables(nEntries))
        total_size = nEntries * entrySize
        call this%all_info_table%shared_alloc(total_size)
        do iEntry = 1, nEntries
            windowStart = (iEntry - 1) * entrySize + 1
            windowEnd = windowStart + entrySize - 1

            this%info_tables(iEntry)%ptr => this%all_info_table%ptr(windowStart : windowEnd)
        end do
    end subroutine

    function get_info_entry(this, iEntry, tgt) result(info)
        debug_function_name("get_info_entry")
        class(GAS_PCHB_Doubles_t) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: info

        ASSERT(associated(this%info_tables(iEntry)%ptr))
        info = this%info_tables(iEntry)%ptr(tgt)
    end function

    subroutine set_info_entry(this, ij, ab, info)
        class(GAS_PCHB_Doubles_t) :: this
        integer, intent(in) :: ij, ab
        integer(int64), intent(in) :: info

        if (iProcIndex_intra == root) then
            this%info_tables(ij)%ptr(ab) = info
        end if
        call this%all_info_table%sync()
    end subroutine

    subroutine init_GugaAliasSampler_t(this, GAS_spec, particle_selection, weighting_selection, use_lookup, create_lookup)
        class(GAS_PCHB_Doubles_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        type(PCHB_ParticleSelection_t), intent(in) :: particle_selection
        type(WeightingChoice_t), intent(in) :: weighting_selection
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = "init_GugaAliasSampler_t"

        real(dp), allocatable :: ij_weights(:, :, :)

        ! also set some more strict defaults for the PCHB implo:
        root_print "Setting reasonable defaults for GUGA-PCHB:"
        if (tau_search_method == possible_tau_search_methods%HISTOGRAMMING) then
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

        ASSERT(GAS_spec%recoupling())

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'GUGA GAS PCHB doubles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'GAS PCHB doubles is using the supergroup lookup'


        block
            integer :: i, j, a, b, ij, ab
            ! The info for double excitations is independent of the super group
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    ij = fuseIndex(i, j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a, b)
                            call this%set_info_entry(ij, ab, get_info(ijab_Index_t([i, j, a, b])))
                        end do
                    end do
                end do
            end do
        end block

        call set_get_weight_pointer(weighting_selection, this%get_weight)
        call this%setup_pchb_sampler_conditional(ij_weights)
        call particle_allocate_and_init(particle_selection, this%particle_selector)
        write(stdout, *) "Finished excitation generator initialization"
    end subroutine

    subroutine finalize_GugaAliasSampler_t(this)
        class(GAS_PCHB_Doubles_t), intent(inout) :: this
        call this%all_info_table%shared_dealloc()
        deallocate(this%info_tables)
        call this%finalize_pchb_sampler()
    end subroutine


    pure function initial_werner_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable

        integer :: i, j, k, l
        HElement_t(dp) :: cpt1, cpt2, cpt3, cpt4

        call exc%idx%extract_ijkl(i, j, k, l)

        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val)
            weight = abs(g(i, j, k, l) + g(i, l, k, j)) / 2.0_dp
        case(dd%ijij%val, dd%jiji%val)
            weight = abs(g(i, j, k, l)) / 2.0_dp
        case(dd%jikk%val, dd%ijkk%val)
            weight = abs(g(i, k, k, j))
        case(dd%iikj%val, dd%iijk%val, dd%iijj%val)
            weight = abs(g(i, l, k, i))
        case(dd%kilj%val, dd%ikjl%val)
            cpt1 = 2 * g(i, j, k, l)
            cpt2 = 2 * g(i, l, k, j)

            cpt3 = abs(cpt1 + cpt2) / 2.0_dp
            cpt4 = abs(cpt1 - cpt2) / 2.0_dp

            if (t_exchange_pchb) then
                weight = abs(cpt4)
            else
                weight = maxval(abs([cpt1, cpt2, cpt3, cpt4]))
            end if
        case(dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            cpt1 = 2 * g(i, j, k, l)
            cpt2 = 2 * g(i, l, k, j)

            cpt3 = abs(cpt2 - cpt1)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1) / 2.0_dp

            if (t_exchange_pchb) then
                weight = abs(cpt2)
            else
                weight = maxval(abs([cpt1, cpt2 / 2.0_dp, cpt3, cpt4]))
            end if
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select

    contains

        elemental function g(i, j, k, l) result(res)
            !! The two-electron integral in chemist's notation
            !!
            !! Convention in the purple book[@purple_book] or **the**
            !! Shavit paper.[@Hinze1981]
            integer, intent(in) :: i, j, k, l
            HElement_t(dp) :: res
            res = get_umat_el(i, k, j, l)
        end function
    end function

    elemental function get_info(idx) result(info)
        type(ijab_Index_t), intent(in) :: idx
        integer(int64) :: info
        type(DistinctDouble_t) :: distinct_double

        info = -1_int64

        distinct_double = DistinctDouble_t(idx)
        if (distinct_double%typ /= distinct_doubles%invalid) then
            block
                integer :: i, j, a, b
                type(ijab_Index_t) :: new_idx

                associate(tmp => distinct_double%get_repr_index())
                    new_idx = tmp%convert()
                end associate

                call new_idx%extract_ijab(i, j, a, b)
                info = encode_excit_info( &
                       distinct_double%typ%repr, &
                       a=a, i=i, b=b, j=j)
            end block
        end if
    end function


    pure subroutine set_get_weight_pointer(weight_choice, ptr)
        type(WeightingChoice_t), intent(in) :: weight_choice
        procedure(get_weight_t), pointer, intent(out) :: ptr
        routine_name("set_get_weight_pointer")

        if (weight_choice == weight_choice_vals%INITIAL_WERNER) then
            ptr => initial_werner_get_weight
        else
            call stop_all(this_routine, "Invalid weight choice: "//weight_choice%str)
        end if
    end subroutine

    pure function WeightingChoice_from_str(str) result(res)
        character(*), intent(in) :: str
        type(WeightingChoice_t) :: res
        routine_name("WeightingChoice_from_str")

        select case(to_upper(str))
        case(weight_choice_vals%INITIAL_WERNER%str)
            res = weight_choice_vals%INITIAL_WERNER
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

    pure function HoleAlgorithm_from_str(str) result(res)
        character(*), intent(in) :: str
        type(PCHB_HoleAlgorithm_t) :: res
        routine_name("WeightingChoice_from_str")

        select case(to_upper(str))
        case(algorithm_vals%FAST_FAST%str)
            res = algorithm_vals%FAST_FAST
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

    logical elemental function symmetry_allowed(exc)
        type(DistinctDouble_t), intent(in) :: exc
        integer :: i, j, a, b
        call exc%idx%extract_ijkl(a, i, b, j)
        symmetry_allowed = &
            RandExcitSymLabelProd(SpinOrbSymLabel(2 * i), SpinOrbSymLabel(2 * j)) &
            == RandExcitSymLabelProd(SpinOrbSymLabel(2 * a), SpinOrbSymLabel(2 * b))
    end function

end module
