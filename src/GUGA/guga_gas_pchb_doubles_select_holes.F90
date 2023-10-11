#include "macros.h"
module guga_gas_pchb_doubles_select_holes

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
    public :: GAS_PCHB_FastFast_t, PCHB_HoleAlgorithm_t, PCHB_HoleAlgorithm_vals_t, &
        WeightingChoice_t, WeightingChoice_vals_t

    abstract interface
        real(dp) function get_weight_t(exc)
            import :: DistinctDouble_t, dp
            type(DistinctDouble_t), intent(in) :: exc
        end function
    end interface

    type, extends(EnumBase_t) :: WeightingChoice_t
        character(20) :: str
    end type

    type :: WeightingChoice_vals_t
        type(WeightingChoice_t) :: &
            INITIAL_WERNER = WeightingChoice_t(1, 'INITIAL-WERNER'), &
            UNIFORM = WeightingChoice_t(2, 'UNIFORM'), &
            ABS_INTEGRAL = WeightingChoice_t(3, 'ABS-INTEGRAL'), &
            N4_TETRAHEDRON = WeightingChoice_t(4, 'N4-TETRAHEDRON'), &
            RANDOM = WeightingChoice_t(5, 'RANDOM'), &
            TWEAKED_WERNER = WeightingChoice_t(6, 'TWEAKED-WERNER'), &
            NEW_WERNER = WeightingChoice_t(7, 'NEW-WERNER')
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
        procedure, public :: init => init_GAS_PCHB_Doubles_t
        procedure, public :: finalize => finalize_GugaAliasSampler_t


        procedure :: new_info_table
        procedure :: set_info_entry
        procedure :: get_info_entry

        procedure(calculate_pchb_sampler_t), deferred :: setup_pchb_sampler_conditional
        procedure(additional_finalize_t), deferred :: finalize_pchb_sampler
    end type

    abstract interface
        subroutine calculate_pchb_sampler_t(this, rank_with_info, ij_weights)
            import :: GAS_PCHB_Doubles_t, dp
            class(GAS_PCHB_Doubles_t), intent(inout) :: this
            integer, intent(in) :: rank_with_info
            real(dp), allocatable, intent(out) :: ij_weights(:, :, :)
        end subroutine

        subroutine additional_finalize_t(this)
            import :: GAS_PCHB_Doubles_t
            implicit none
            class(GAS_PCHB_Doubles_t), intent(inout) :: this
        end subroutine
    end interface



    type, extends(GAS_PCHB_Doubles_t) :: GAS_PCHB_FastFast_t
        private

        !> The shape is (fused_number_of_double_excitations, n_supergroup)
        type(AliasSampler_2D_t) :: alias_sampler
        integer, allocatable :: tgtOrbs(:, :)

    contains
        private
        procedure, public :: pickOrbitals
        procedure, public :: calc_pgen

        procedure, public :: calc_orbital_pgen_contr
        procedure, public :: calc_orbital_pgen_contr_start
        procedure, public :: calc_orbital_pgen_contr_end

        procedure :: setup_pchb_sampler_conditional
        procedure :: finalize_pchb_sampler
    end type

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

    subroutine init_GAS_PCHB_Doubles_t(this, GAS_spec, particle_selection, weighting_selection, use_lookup, create_lookup)
        class(GAS_PCHB_Doubles_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        type(PCHB_ParticleSelection_t), intent(in) :: particle_selection
        type(WeightingChoice_t), intent(in) :: weighting_selection
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = "init_GAS_PCHB_Doubles_t"

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

        call set_get_weight_pointer(weighting_selection, this%get_weight)

        call this%setup_pchb_sampler_conditional(root, ij_weights)
        call particle_allocate_and_init(particle_selection, this%GAS_spec, ij_weights, root, use_lookup, this%particle_selector)
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
    end function

    pure function new_werner_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable

        integer :: i, j, k, l

        call exc%idx%extract_ijkl(i, j, k, l)

        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val)
            weight = abs(g(i, j, k, l)) + abs(g(i, l, k, j))
        case(dd%ijij%val, dd%jiji%val)
            weight = abs(g(i, j, k, l)) + 1e-1_dp
        case(dd%jikk%val, dd%ijkk%val)
            weight = abs(g(i, k, k, j))
        case(dd%iikj%val, dd%iijk%val, dd%iijj%val)
            weight = abs(g(i, l, k, i))
        case(dd%kilj%val, dd%ikjl%val, dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = 2 * (abs(g(i, j, k, l)) + abs(g(i, l, k, j)))
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    pure function tweaked_werner_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable

        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = initial_werner_get_weight(exc) + 1e-4_dp
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    elemental function g(i, j, k, l) result(res)
        !! The two-electron integral in chemist's notation
        !!
        !! Convention in the purple book[@purple_book] or **the**
        !! Shavit paper.[@Hinze1981]
        integer, intent(in) :: i, j, k, l
        HElement_t(dp) :: res
        res = get_umat_el(i, k, j, l)
    end function

    function uniform_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("uniform_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = 1._dp
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    function random_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("uniform_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = 1._dp + get_r() * 100
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select

    contains
        ! We cannot use the own random number generator of NECI,
        ! it is not yet initialized
        real(dp) function get_r()
            call random_number(get_r)
        end function
    end function

    function abs_integral_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("abs_integral_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        integer :: i, j, k, l

        call exc%idx%extract_ijkl(i, j, k, l)
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = abs(g(i, j, k, l)) + abs(g(i, l, k, j))
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
    end function

    function N4_tetrahedron_get_weight(exc) result(weight)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: weight
        routine_name("N4_tetrahedron_get_weight")
        type(DistinctDoubleTypeVals_t), parameter :: dd = distinct_doubles
            !! Just a shorter name for an often used variable
        integer :: i, j, k, l

        if (nSpatOrbs /= 12) then
            call stop_all(this_routine, "N4 tetrahedron required")
        end if

        call exc%idx%extract_ijkl(i, j, k, l)
        select case (exc%typ%val)
        case(dd%jijk%val, dd%ijkj%val, dd%ikjk%val, dd%kikj%val, dd%ijik%val, dd%jiki%val, &
                dd%ijij%val, dd%jiji%val, dd%jikk%val, dd%ijkk%val, &
                dd%iikj%val, dd%iijk%val, dd%iijj%val, dd%kilj%val, dd%ikjl%val, &
                dd%ijkl%val, dd%jilk%val, dd%ijlk%val, dd%jikl%val)
            weight = merge(3, 1, mod(i - 1, 3) == mod(j - 1, 3)) * merge(3, 1, mod(k - 1, 3) == mod(l - 1, 3))
        case(dd%invalid%val)
            weight = 0._dp
        case default
            call stop_all(this_routine, 'forgotten case')
        end select
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


    subroutine set_get_weight_pointer(weight_choice, ptr)
        type(WeightingChoice_t), intent(in) :: weight_choice
        procedure(get_weight_t), pointer, intent(out) :: ptr
        routine_name("set_get_weight_pointer")

        if (weight_choice == weight_choice_vals%INITIAL_WERNER) then
            ptr => initial_werner_get_weight
        else if (weight_choice == weight_choice_vals%TWEAKED_WERNER) then
            ptr => tweaked_werner_get_weight
        else if (weight_choice == weight_choice_vals%NEW_WERNER) then
            ptr => new_werner_get_weight
        else if (weight_choice == weight_choice_vals%UNIFORM) then
            ptr => uniform_get_weight
        else if (weight_choice == weight_choice_vals%ABS_INTEGRAL) then
            ptr => abs_integral_get_weight
        else if (weight_choice == weight_choice_vals%N4_TETRAHEDRON) then
            ptr => N4_tetrahedron_get_weight
        else if (weight_choice == weight_choice_vals%RANDOM) then
            ptr => random_get_weight
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
        case(weight_choice_vals%TWEAKED_WERNER%str)
            res = weight_choice_vals%TWEAKED_WERNER
        case(weight_choice_vals%NEW_WERNER%str)
            res = weight_choice_vals%NEW_WERNER
        case(weight_choice_vals%UNIFORM%str)
            res = weight_choice_vals%UNIFORM
        case(weight_choice_vals%ABS_INTEGRAL%str)
            res = weight_choice_vals%ABS_INTEGRAL
        case(weight_choice_vals%N4_TETRAHEDRON%str)
            res = weight_choice_vals%N4_TETRAHEDRON
        case(weight_choice_vals%RANDOM%str)
            res = weight_choice_vals%RANDOM
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


    subroutine setup_pchb_sampler_conditional(this, rank_with_info, ij_weights)
        class(GAS_PCHB_FastFast_t), intent(inout) :: this
        integer, intent(in) :: rank_with_info
        real(dp), allocatable, intent(out) :: ij_weights(:, :, :)
        integer :: i, j, ij, a, b, ab, i_sg
        integer(int64) :: ijMax, abMax, memCost
        integer, allocatable :: supergroups(:, :)
        real(dp), allocatable :: w(:)

        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! initialize the mapping ab -> (a, b)
        abMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        ijMax = abMax
        allocate(this%tgtOrbs(2, 0:abMax), source=0)
        do a = 1, nSpatOrbs
            do b = 1, a
                ab = fuseIndex(a, b)
                this%tgtOrbs(1, ab) = b
                this%tgtOrbs(2, ab) = a
            end do
        end do


        calculate_memory_demand: block
            integer(int64) :: number_of_fused_indices, bytes_per_sampler, n_supergroup
            number_of_fused_indices = int(ijMax, int64)
            bytes_per_sampler = (int(abMax, int64) * 3_int64 * 8_int64)
            n_supergroup = size(supergroups, 2, kind=int64)
            memCost = number_of_fused_indices * bytes_per_sampler * n_supergroup
        endblock calculate_memory_demand

        write(stdout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(stdout, *) "The number of supergroups is", size(supergroups, 2)
        write(stdout, *) "Generating samplers for PCHB excitation generator"
        write(stdout, *) "Depending on the number of supergroups this can take up to 10min."
        call this%alias_sampler%shared_alloc([int(ijMax), size(supergroups, 2)], int(abMax), 'PCHB')
        call this%new_info_table(int(ijMax, int64), int(abMax, int64))


        if (iProcIndex_intra == rank_with_info) then
            allocate(w(abMax), source=0._dp)
            allocate(ij_weights(nSpatOrbs, nSpatOrbs, size(supergroups, 2)), source=0._dp)
        else
            allocate(w(0))
            allocate(ij_weights(0, 0, 0))
        end if

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

        do i_sg = 1, size(supergroups, 2)
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w = 0.0_dp
                    ij = fuseIndex(i, j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a, b)
                            if (iProcIndex_intra == rank_with_info) then
                                associate(exc => DistinctDouble_t(ijab_Index_t([i, j, a, b])))
                                if (this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) .and. symmetry_allowed(exc)) then
                                    w(ab) = this%get_weight(exc)
                                end if
                                end associate
                            end if
                        end do
                    end do
                    if (iProcIndex_intra == root) then
                        ij_weights(i, j, i_sg) = sum(w)
                        ij_weights(j, i, i_sg) = ij_weights(i, j, i_sg)
                    end if
                    call this%alias_sampler%setup_entry(ij, i_sg, rank_with_info, w)
                end do
            end do
        end do
    end subroutine setup_pchb_sampler_conditional

    subroutine finalize_pchb_sampler(this)
        class(GAS_PCHB_FastFast_t), intent(inout) :: this
        call this%alias_sampler%finalize()
        deallocate(this%tgtOrbs)
    end subroutine

    subroutine pickOrbitals(this, ilut, nI, csf_i, excitInfo, pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), a, b, orbs(2), ij, ab, i_sg
        real(dp) :: pgen_elec, pgen_orbs
        routine_name("pickOrbitals")

        unused_var(ilut)
        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        call this%particle_selector%draw(nI, ilut, csf_i, i_sg, src, pgen_elec)
        ASSERT(src(1) <= src(2))
        if (src(1) == 0) then
            excitInfo%valid = .false.
            pgen = 0.0_dp
            return
        end if

        ! use the sampler for this electron pair -> order of src electrons
        ! does not matter
        ij = fuseIndex(src(1), src(2))

        ! get a pair of orbitals using the precomputed weights
        call this%alias_sampler%sample(ij, i_sg, ab, pgen_orbs)

        ! can happen for GAS and is ok
        if (ab == 0) then
            excitInfo%valid = .false.
            pgen = 0.0_dp
            return
        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        else if (near_zero(pgen_orbs)) then
            excitInfo%valid = .false.
            pgen = 0.0_dp
            ! Yes, print. Those events are signficant enough to be always noted in the output
            print *, "WARNING: Generated excitation with probability of 0"
            return
        end if

        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = this%tgtOrbs(:, ab)

        a = orbs(1); b = orbs(2)

        ! check if the picked orbs are a valid choice - if they are the same, match one
        ! occupied orbital or are zero (maybe because there are no allowed picks for
        ! the given source) abort
        ! these are the 'easy' checks for GUGA.. more checks need to be done
        ! to see if it is actually a valid combination..
        if (any(orbs == 0) &
               .or. (csf_i%stepvector(a) == 3) &
               .or. (csf_i%stepvector(b) == 3) &
               .or. (a == b .and. csf_i%stepvector(b) /= 0)) then
            excitInfo%valid = .false.
            return
        end if

        call extract_excit_info(this%get_info_entry(ij, ab), excitInfo)
        excitInfo%i_sg_start = i_sg

        pGen = pgen_elec * pgen_orbs
    end subroutine pickOrbitals

    function calc_pgen(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: ij, ab
        real(dp) :: p_elec
        unused_var(nI); unused_var(ilutI); unused_var(csf_i)

        ! i, j, k and l entries have the info about the electrons and
        ! holes of the excitation: E_{ij}E_{kl} is the convention ...
        ij = fuseIndex(excitInfo%j, excitInfo%l)
        ab = fuseIndex(excitInfo%i, excitInfo%k)

        p_elec = this%particle_selector%get_pgen(nI, csf_i, excitInfo%i_sg_start, excitInfo%j, excitInfo%l)
        pgen = pDoubles * p_elec * this%alias_sampler%get_prob(ij, excitInfo%i_sg_start, ab)
    end function calc_pgen

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    pure subroutine calc_orbital_pgen_contr(this, csf_i, occ_orbs, above_cpt, below_cpt)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt

        integer :: ij, i_sg

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if
        ! this function can in theory be called with both i < j and i > j..
        ! to take the correct values here!
        ! although for the fuseIndex function the order is irrelevant

        ! i think i have to consider both the above and below contribution..
        ! but i am not so sure how.. in this PCHB case..
        ! maybe in PCHB those 2 are just the same.. i am confused
        ij = fuseIndex(gtID(occ_orbs(1)), gtID(occ_orbs(2)))

        ! and both I and J are electron and hole indices here
        above_cpt = this%alias_sampler%get_prob(ij, i_sg, ij) / 2.0_dp
        ! but since they will get added in later routines i actually have
        ! do divide by 2 here..

        ! and in the PCHB there is no difference between the 2!
        below_cpt = above_cpt

    end subroutine calc_orbital_pgen_contr

    pure subroutine calc_orbital_pgen_contr_start(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = 'calc_orbital_pgen_contr_start'

        integer :: i, j, ij, ab, i_sg

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        ! depending on type (R->L / L->R) a can be > j or < j, but always > i
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))
        ASSERT(i < orb_a)

        ! here i is both electron and hole index!
        ij = fuseIndex(i, j)
        ab = fuseIndex(i, orb_a)

        orb_pgen = this%alias_sampler%get_prob(ij, i_sg, ab)
    end subroutine calc_orbital_pgen_contr_start

    subroutine calc_orbital_pgen_contr_end(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = "calc_orbital_pgen_contr_end"

        integer :: i, j, ij, ab, i_sg

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! and j is at the same time electron and hole index!

        ! depending on L->R or R->L type a can be > i o r < i, but always < j!
        ASSERT(orb_a < j)

        ! j here is both elec and hole ind!
        ij = fuseIndex(i, j)
        ab = fuseIndex(orb_a, j)

        orb_pgen = this%alias_sampler%get_prob(ij, i_sg, ab)
    end subroutine calc_orbital_pgen_contr_end
end module
