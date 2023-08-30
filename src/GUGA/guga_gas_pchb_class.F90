#include "macros.h"
module guga_gas_pchb_mod

    use aliasSampling, only: AliasSampler_2D_t
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm, bits_n_int
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel, G1, t_guga_pchb_weighted_singles, &
                          nBasis, nSpatOrbs, ElecPairs, &
                          t_analyze_pchb, t_old_pchb, t_exchange_pchb
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type, ExcitationType_t, &
        DistinctDoubleType_t, distinct_doubles, DistinctDouble_t, ijab_Index_t
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut, CSF_Info_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: near_zero, fuseIndex, binary_search_first_ge, &
                        get_free_unit, stop_all, operator(.isclose.), operator(.div.)
    use CalcData, only: t_matele_cutoff, matele_cutoff, t_truncate_spawns
    use tau_main, only: max_tau, tau_search_method, possible_tau_search_methods
    use tau_search_hist, only: frq_ratio_cutoff, max_frequency_bound, n_frequency_bins
    use sym_general_mod, only: ClassCountInd
    use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel
    use excitation_types, only: Excite_1_t, Excite_2_t
    use UMatCache, only: gtID
    use guga_excitations, only: assign_excitinfo_values_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, get_guga_integral_contrib_spat, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo, &
                                checkCompatibility_single
    use guga_bitrepops, only: identify_excitation, encode_excit_info, extract_excit_info, &
                              contract_2_rdm_ind
    use guga_base_class, only: GUGABaseClass_t
    use bit_reps, only: decode_bit_det
    use shared_array, only: shared_array_int64_t, shared_array_real_t
    use MPI_wrapper, only: iProcIndex_intra, iprocindex, root
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use procedure_pointers, only: get_umat_el

    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use sets_mod, only: is_sorted, is_set

    better_implicit_none

    private
    public :: GAS_PCHB_t

    type, extends(GUGABaseClass_t) :: GAS_PCHB_t
        private
        !> Use a lookup for the supergroup index in global_det_data
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.

        !> The shape is (fused_number_of_double_excitations, n_supergroup)
        type(AliasSampler_2D_t) :: alias_sampler
        integer, allocatable :: tgtOrbs(:, :)

        type(shared_array_int64_t) :: all_info_table
        type(shared_array_int64_t), allocatable :: info_tables(:)

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(GASSpec_t), allocatable :: GAS_spec

        integer(n_int), allocatable :: allowed_holes(:, :, :)
            !! Bitmask for the allowed holes for a single excitation
        integer :: L_spat_bits = -1
            !! The number of  integer(n_int) numbers to
            !! store information about spatial orbitals as bitmask
    contains
        private
        procedure, public :: init => init_GugaAliasSampler_t
        procedure, public :: finalize => finalize_GugaAliasSampler_t

        procedure, public :: pickOrbitals_double
        procedure, public :: pickOrbitals_single
        procedure, public :: get_valid_orbs
        procedure :: singles_uniform_init
        procedure, public :: calc_pgen

        procedure :: calc_orb_pgen_uniform_singles

        procedure, public :: calc_orbital_pgen_contr
        procedure, public :: calc_orbital_pgen_contr_start
        procedure, public :: calc_orbital_pgen_contr_end

        procedure :: new_info_table
        procedure :: set_info_entry
        procedure :: get_info_entry

        procedure :: calc_orb_pgen_guga_pchb_double

        procedure :: setup_pchb_sampler_conditional

    endtype GAS_PCHB_t

contains

    subroutine new_info_table(this, nEntries, entrySize)
        class(GAS_PCHB_t) :: this
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
        class(GAS_PCHB_t) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: info

        ASSERT(associated(this%info_tables(iEntry)%ptr))
        info = this%info_tables(iEntry)%ptr(tgt)
    end function

    subroutine set_info_entry(this, ij, ab, info)
        class(GAS_PCHB_t) :: this
        integer, intent(in) :: ij, ab
        integer(int64), intent(in) :: info

        if (iProcIndex_intra == 0) then
            this%info_tables(ij)%ptr(ab) = info
        end if
        call this%all_info_table%sync()
    end subroutine

    subroutine init_GugaAliasSampler_t(this, GAS_spec, use_lookup, create_lookup)
        class(GAS_PCHB_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = "init_GugaAliasSampler_t"

        integer :: a, b
        integer(int64) :: abMax, ab

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

        ! initialize the mapping ab -> (a, b)
        abMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        allocate(this%tgtOrbs(2, 0:abMax), source=0)
        do a = 1, nSpatOrbs
            do b = 1, a
                ab = fuseIndex(a, b)
                this%tgtOrbs(1, ab) = b
                this%tgtOrbs(2, ab) = a
            end do
        end do

        call this%setup_pchb_sampler_conditional()

        call this%singles_uniform_init()

        write(stdout, *) "Finished excitation generator initialization"
    end subroutine init_GugaAliasSampler_t

    subroutine setup_pchb_sampler_conditional(this)
        class(GAS_PCHB_t), intent(inout) :: this
        integer :: i, j, ij, a, b, ab, i_sg
        integer(int64) :: ijMax, abMax, memCost
        integer, allocatable :: supergroups(:, :)
        real(dp), allocatable :: w(:)

        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        abMax = ijMax

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

        ! the encode_excit_info function encodes as: E_{ai}E_{bj)
        ! but for some index combinations we have to change the input so
        ! actually E_{aj}E_{bi} is encoded! convention to use only
        ! certain type of excit-types

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

        allocate(w(abMax), source=0.0_dp)
        do i_sg = 1, size(supergroups, 2)
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w = 0.0_dp
                    ij = fuseIndex(i, j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a, b)
                            if (this%GAS_spec%is_allowed(Excite_2_t(2 * i, 2 * a, 2 * j, 2 * b), supergroups(:, i_sg)) &
                                    .and. symmetry_allowed(i, j, a, b)) then
                                w(ab) = get_weight(DistinctDouble_t(ijab_Index_t([i, j, a, b])))
                            end if
                        end do
                    end do
                    call this%alias_sampler%setup_entry(ij, i_sg, root, w)
                end do
            end do
        end do
    contains
        logical function symmetry_allowed(i, j, a, b)
            integer, intent(in) :: i, j, a, b
            ! shoud i point group symmetry restrctions here?
            ! this would avoid unnecessary if statements..
            symmetry_allowed = &
                RandExcitSymLabelProd(SpinOrbSymLabel(2 * i), SpinOrbSymLabel(2 * j)) &
                == RandExcitSymLabelProd(SpinOrbSymLabel(2 * a), SpinOrbSymLabel(2 * b))
        end function
    end subroutine setup_pchb_sampler_conditional

    subroutine finalize_GugaAliasSampler_t(this)
        class(GAS_PCHB_t), intent(inout) :: this
        call this%alias_sampler%finalize()
        call this%all_info_table%shared_dealloc()
        deallocate(this%tgtOrbs, this%info_tables, this%allowed_holes)
    end subroutine

    subroutine pickOrbitals_double(this, ilut, nI, csf_i, excitInfo, pgen)
        class(GAS_PCHB_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pickOrbitals_double"

        integer :: src(2), sym_prod, sum_ml, i, j, a, b, orbs(2), ij, ab, i_sg
        real(dp) :: pgen_elec, pgen_orbs

        unused_var(ilut)
        if (this%use_lookup) then
            ! TODO implement
            call stop_all(this_routine, 'has to be implemented')
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        ! maybe I will also produce a weighted electron pickin in the
        ! GUGA formalism.. but for now pick them uniformly:
        call pick_elec_pair_uniform_guga(nI, src, sym_prod, sum_ml, pgen_elec)
        ! make a different picker, which does not bias towards doubly
        ! occupied orbitals here too!

        ASSERT(src(1) < src(2))

        i = gtID(src(1))
        j = gtID(src(2))

        if (i /= j) then
            if (csf_i%stepvector(i) == 3) pgen_elec = 2.0_dp * pgen_elec
            if (csf_i%stepvector(j) == 3) pgen_elec = 2.0_dp * pgen_elec
        end if

        ! use the sampler for this electron pair -> order of src electrons
        ! does not matter
        ij = fuseIndex(i, j)

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

        a = orbs(1)
        b = orbs(2)

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

        ! setup a getInfo functionality in the sampler!
        ! TODO: include supergroup into bitmask
        call extract_excit_info(this%get_info_entry(ij, ab), excitInfo)
        excitInfo%i_sg_start = i_sg

        pGen = pgen_elec * pgen_orbs
    end subroutine pickOrbitals_double

    subroutine singles_uniform_init(this)
        class(GAS_PCHB_t), intent(inout) :: this

        integer :: tgt, src, i_sg
        integer, allocatable :: supergroups(:, :)

        this%L_spat_bits = nSpatOrbs .div. bits_n_int

        supergroups = this%indexer%get_supergroups()

        allocate(this%allowed_holes(0:this%L_spat_bits, nSpatOrbs, size(supergroups, 2)), &
                 source=0_n_int)


        do i_sg = 1, size(supergroups, 2)
            do src = 1, nSpatOrbs
                do tgt = 1, nSpatOrbs
                block
                    type(Excite_1_t) :: exc
                    exc = Excite_1_t(2 * src, 2 * tgt)
                    if (src /= tgt &
                            .and. this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) &
                            .and. symmetry_allowed(exc)) then
                        call my_set_orb(this%allowed_holes(:, src, i_sg), tgt)
                    end if
                end block
                end do
            end do
        end do

    contains

        ! Ugly Fortran syntax rules forces us to not use the macro.
        ! Even associate would not help here
        ! https://stackoverflow.com/questions/65734764/non-one-indexed-array-from-associate
        pure subroutine my_set_orb(ilut, orb)
            integer(n_int), intent(inout) :: ilut(0:this%L_spat_bits)
            integer, intent(in) :: orb
            set_orb(ilut, orb)
        end subroutine

        logical pure function symmetry_allowed(exc)
            use SymExcitDataMod, only: SpinOrbSymLabel
            type(Excite_1_t), intent(in) :: exc
            symmetry_allowed = SpinOrbSymLabel(exc%val(1, 1)) == SpinOrbSymLabel(exc%val(2, 1))
        end function
    end subroutine

    subroutine pickOrbitals_single(this, ilut, nI, csf_i, excitInfo, pgen)
        debug_function_name("pickOrbitals_single")
        class(GAS_PCHB_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src_orb, tgt_orb, i
        integer, allocatable :: occ_spat_orbs(:)
            !! The occupied spatial orbitals
        integer, allocatable :: valid_orbs(:)
        unused_var(ilut); unused_var(nI)

        occ_spat_orbs = pack([(i, i=1, size(csf_i%Occ_int))], csf_i%Occ_int /= 0)
        src_orb = occ_spat_orbs(1 + int(genrand_real2_dSFMT() * size(occ_spat_orbs)))

        valid_orbs = this%get_valid_orbs(csf_i, src_orb)

        if (size(valid_orbs) == 0) then
            pgen = 0.0_dp
            tgt_orb = 0
            excitInfo%valid = .false.
        else
            pgen = 1._dp / real(size(occ_spat_orbs) * size(valid_orbs), dp)
            tgt_orb = valid_orbs(1 + floor(genrand_real2_dSFMT() * size(valid_orbs)))


            if (tgt_orb < src_orb) then
                excitInfo = assign_excitinfo_values_single( &
                            gen_type%R, tgt_orb, src_orb, tgt_orb, src_orb)
            else
                excitInfo = assign_excitinfo_values_single( &
                            gen_type%L, tgt_orb, src_orb, src_orb, tgt_orb)
            end if
            excitInfo%valid = checkCompatibility_single(csf_i, excitInfo)
            if (excitInfo%valid) then
                ASSERT(pgen .isclose. this%calc_orb_pgen_uniform_singles(csf_i, excitInfo))
            else
                pgen = 0._dp
            end if
        end if
    end subroutine pickOrbitals_single

    function get_valid_orbs(this, csf_i, src_orb) result(valid_orbs)
        class(GAS_PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: src_orb
        integer, allocatable :: valid_orbs(:)
        character(*), parameter :: this_routine = 'get_valid_orbs'

        integer(n_int) :: &
            ilut_doubly_occupied(0:this%L_spat_bits), ilut_valid(0:this%L_spat_bits)
        integer :: i_sg

        if (this%use_lookup) then
            ! TODO implement
            call stop_all(this_routine, 'has to be implemented')
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        block
            integer, allocatable :: doubly_occupied(:)
            integer :: i
            doubly_occupied = pack([(i, i=1, size(csf_i%stepvector))], csf_i%stepvector == 3)
            ilut_doubly_occupied = 0_n_int
            do i = 1, size(doubly_occupied)
                set_orb(ilut_doubly_occupied, doubly_occupied(i))
            end do
        endblock

        ilut_valid = iand(this%allowed_holes(:, src_orb, i_sg), not(ilut_doubly_occupied))
        allocate(valid_orbs(sum(popcnt(ilut_valid))))
        if (size(valid_orbs) /= 0) then
            call decode_bit_det(valid_orbs, ilut_valid)
        end if
    end function get_valid_orbs

    function calc_orb_pgen_uniform_singles(this, csf_i, excitInfo) result(pgen)
        debug_function_name("calc_orb_pgen_uniform_singles")
        class(GAS_PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: nOcc, nUnocc
        integer, allocatable :: valid_orbs(:)

        ASSERT(1 <= excitInfo%i .and. excitInfo%i <= nSpatOrbs)
        ASSERT(1 <= excitInfo%j .and. excitInfo%j <= nSpatOrbs)

        if (excitInfo%i == excitInfo%j &
               .or. csf_i%stepvector(excitInfo%i) == 3 &
               .or. csf_i%stepvector(excitInfo%j) == 0) then
            pgen = 0.0_dp
        else
            nOcc = count(csf_i%Occ_int /= 0)
            valid_orbs = this%get_valid_orbs(csf_i, excitInfo%j)
            nUnocc = size(valid_orbs)
            pgen = 1.0_dp / real(nOcc * nUnocc, dp)
        end if
    end function calc_orb_pgen_uniform_singles

    function calc_pgen(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(GAS_PCHB_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        integer :: ic

        ic = get_excit_level_from_excitInfo(excitInfo)

        if (ic == 1) then
            if (t_guga_pchb_weighted_singles) then
                pgen = calc_pgen_mol_guga_single(nI, ilutI, csf_i, excitInfo)
            else
                pgen = this%calc_orb_pgen_uniform_singles(csf_i, excitInfo)
            end if

            pgen = pgen * pSingles

        else if (ic == 2) then

            pgen = pDoubles * this%calc_orb_pgen_guga_pchb_double(excitInfo)

        else
            pgen = 0.0_dp
        end if
    end function calc_pgen

    pure function calc_orb_pgen_guga_pchb_double(this, excitInfo) result(pgen)
        class(GAS_PCHB_t), intent(in) :: this
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: ij, ab
        real(dp) :: p_elec

        ! i, j, k and l entries have the info about the electrons and
        ! holes of the excitation: E_{ij}E_{kl} is the convention ...
        ij = fuseIndex(excitInfo%j, excitInfo%l)

        p_elec = 1.0_dp / real(ElecPairs, dp)

        ab = fuseIndex(excitInfo%i, excitInfo%k)
        pgen = p_elec * this%alias_sampler%get_prob(ij, excitInfo%i_sg_start, ab)

    end function

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    pure subroutine calc_orbital_pgen_contr(this, csf_i, occ_orbs, above_cpt, below_cpt)
        class(GAS_PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt
        character(*), parameter :: this_routine = "calc_orbital_pgen_contr"

        integer :: ij, i_sg

        if (this%use_lookup) then
            ! TODO implement
            call stop_all(this_routine, 'has to be implemented')
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
        class(GAS_PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = 'calc_orbital_pgen_contr_start'

        integer :: i, j, ij, ab, i_sg

        if (this%use_lookup) then
            ! TODO implement
            call stop_all(this_routine, 'has to be implemented')
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

    pure subroutine calc_orbital_pgen_contr_end(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(GAS_PCHB_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = "calc_orbital_pgen_contr_end"

        integer :: i, j, ij, ab, i_sg

        if (this%use_lookup) then
            ! TODO implement
            call stop_all(this_routine, 'has to be implemented')
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

    pure function get_weight(exc) result(integral)
        type(DistinctDouble_t), intent(in) :: exc
        real(dp) :: integral

        integer :: i, j, a, b
        logical :: flag_
        HElement_t(dp) :: cpt1, cpt2, cpt3, cpt4

        flag_ = t_exchange_pchb
        associate(ijab_idx => exc%idx%convert())
            call ijab_idx%extract(i, j, a, b)
        end associate
        select case (exc%typ%val)
        case(distinct_doubles%jijk%val)
            integral = abs(get_umat_el(a, a, i, j) + get_umat_el(a, a, j, i)) / 2.0_dp
        case(distinct_doubles%ijkj%val)
            integral = abs(get_umat_el(a, b, i, i) + get_umat_el(b, a, i, i)) / 2.0_dp
        case(distinct_doubles%ikjk%val)
            integral = abs(get_umat_el(a, b, i, i) + get_umat_el(b, a, i, i)) / 2.0_dp
        case(distinct_doubles%kikj%val)
            integral = abs(get_umat_el(a, a, i, j) + get_umat_el(a, a, j, i)) / 2.0_dp
        case(distinct_doubles%jikk%val)
            integral = abs(get_umat_el(b, a, i, b) + get_umat_el(a, b, b, i)) / 2.0_dp
        case(distinct_doubles%ijkk%val)
            integral = abs(get_umat_el(j, a, i, b) + get_umat_el(a, b, j, i)) / 2.0_dp
        case(distinct_doubles%ijik%val)
            integral = abs(get_umat_el(a, a, i, j) + get_umat_el(a, a, j, i)) / 2.0_dp
        case(distinct_doubles%jiki%val)
            integral = abs(get_umat_el(a, b, i, i) + get_umat_el(b, a, i, i)) / 2.0_dp
        case(distinct_doubles%iikj%val)
            integral = abs(get_umat_el(a, b, j, a) + get_umat_el(b, a, a, j)) / 2.0_dp
        case(distinct_doubles%iijk%val)
            integral = abs(get_umat_el(a, b, j, a) + get_umat_el(b, a, a, j)) / 2.0_dp
        case(distinct_doubles%ijij%val)
            integral = abs(get_umat_el(a, a, i, i)) / 2.0_dp
        case(distinct_doubles%jiji%val)
            integral = abs(get_umat_el(a, a, i, i)) / 2.0_dp
        case(distinct_doubles%iijj%val)
            integral = abs(get_umat_el(a, j, j, a) + get_umat_el(j, a, a, j)) / 2.0_dp
        case(distinct_doubles%ikjl%val)
            cpt1 = (get_umat_el(b, a, j, i) + get_umat_el(a, b, i, j))
            cpt2 = (get_umat_el(a, b, j, i) + get_umat_el(b, a, i, j))
            cpt3 = abs(cpt1 + cpt2) / 2.0_dp
            cpt4 = abs(cpt1 - cpt2) / 2.0_dp

            cpt1 = abs(cpt1)
            cpt2 = abs(cpt2)

            if (flag_) then
                integral = abs(cpt4)
            else
                integral = maxval(abs([cpt1, cpt2, cpt3, cpt4]))
            end if
        case(distinct_doubles%kilj%val)
            cpt1 = (get_umat_el(a, b, i, j) + get_umat_el(b, a, j, i))
            cpt2 = (get_umat_el(b, a, i, j) + get_umat_el(a, b, j, i))
            cpt3 = abs(cpt1 + cpt2) / 2.0_dp
            cpt4 = abs(cpt1 - cpt2) / 2.0_dp

            cpt1 = abs(cpt1)
            cpt2 = abs(cpt2)

            if (flag_) then
                integral = abs(cpt4)
            else
                integral = maxval(abs([cpt1, cpt2, cpt3, cpt4]))
            end if
        case(distinct_doubles%ijkl%val)
            cpt1 = (get_umat_el(a, b, i, j) + get_umat_el(b, a, j, i))
            cpt2 = (get_umat_el(b, a, i, j) + get_umat_el(a, b, j, i))

            cpt3 = abs(cpt2 - cpt1)
            ! cpt4 = abs(cpt1 + cpt2)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1) / 2.0_dp

            if (flag_) then
                integral = abs(cpt2)
            else
                integral = maxval(abs([cpt1, cpt2 / 2.0_dp, cpt3, cpt4]))
            end if
        case(distinct_doubles%jilk%val)
            cpt1 = (get_umat_el(a, b, i, j) + get_umat_el(b, a, j, i))
            cpt2 = (get_umat_el(b, a, i, j) + get_umat_el(a, b, j, i))

            cpt3 = abs(cpt2 - cpt1)
            ! cpt4 = abs(cpt1 + cpt2)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1) / 2.0_dp

            if (flag_) then
                integral = abs(cpt2)
            else
                integral = maxval(abs([cpt1, cpt2 / 2.0_dp, cpt3, cpt4]))
            end if
        case(distinct_doubles%ijlk%val)
            cpt1 = (get_umat_el(a, b, i, j) + get_umat_el(b, a, j, i))
            cpt2 = (get_umat_el(b, a, i, j) + get_umat_el(a, b, j, i))

            cpt3 = abs(cpt2 - cpt1)
            ! cpt4 = abs(cpt1 + cpt2)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1) / 2.0_dp

            if (flag_) then
                integral = abs(cpt2)
            else
                integral = maxval(abs([cpt1, cpt2 / 2.0_dp, cpt3, cpt4]))
            end if
        case(distinct_doubles%jikl%val)
            cpt1 = (get_umat_el(a, b, i, j) + get_umat_el(b, a, j, i))
            cpt2 = (get_umat_el(b, a, i, j) + get_umat_el(a, b, j, i))

            cpt3 = abs(cpt2 - cpt1)
            ! cpt4 = abs(cpt1 + cpt2)
            cpt4 = abs(cpt2 - 2.0_dp * cpt1) / 2.0_dp

            if (flag_) then
                integral = abs(cpt2)
            else
                integral = maxval(abs([cpt1, cpt2 / 2.0_dp, cpt3, cpt4]))
            end if
        case(distinct_doubles%invalid%val)
            integral = 0._dp
        end select
    end function

    function get_info(idx) result(info)
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

                call new_idx%extract(i, j, a, b)
                info = encode_excit_info( &
                       distinct_double%typ%repr, &
                       a=a, i=i, b=b, j=j)
            end block
        end if
    end function

end module guga_gas_pchb_mod
