#include "macros.h"
module guga_gas_pchb_doubles_select_holes

    use iso_c_binding, only: c_sizeof
    use aliasSampling, only: AliasSampler_2D_t, AliasSampler_3D_t, do_direct_calculation
    use constants, only: n_int, dp, int64, stdout, int_rdm, bits_n_int
    use bit_rep_data, only: GugaBits, nIfD
    use bit_reps, only: decode_bit_det, UnoccupiedGetter_t
    use SystemData, only: nel, nSpatOrbs, nBasis
    use guga_data, only: ExcitationInformation_t, &
        distinct_doubles, DistinctDouble_t, DistinctDoubleTypeVals_t, ijab_Index_t, &
        ijkl_Index_t, fused_to_DistinctDouble
    use guga_bitrepops, only: CSF_Info_t
    use util_mod, only: near_zero, stop_all, operator(.isclose.), EnumBase_t, &
        swap, ByteSize_t
    use sets_mod, only: operator(.in.)
    use indexing_mod, only: fuse_symm_idx
    use CalcData, only: t_truncate_spawns
    use fortran_strings, only: to_upper
    use SymExcitDataMod, only: SpinOrbSymLabel
    use UMatCache, only: gtID
    use guga_bitrepops, only: encode_excit_info, extract_excit_info
    use guga_base_class, only: DoublesGUGABase_t
    use shared_array, only: shared_array_int64_t
    use MPI_wrapper, only: iProcIndex_intra, iProcIndex, root
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use sort_mod, only: sort

    use guga_pchb_doubles_weights_mod, only: WeightingChoice_t, get_doubles_weight_t, &
        WeightingChoice_vals_t, set_doubles_weight_pointer

    use growing_buffers, only: buffer_int64_1D_t

    use gasci, only: GASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer

    use guga_gas_pchb_doubles_select_particles, only: ParticleSelector_t, &
        PCHB_ParticleSelection_t, particle_allocate_and_init => allocate_and_init


    better_implicit_none

    private
    public :: GAS_PCHB_FastFast_t, PCHB_HoleAlgorithm_t, PCHB_HoleAlgorithm_vals_t, &
        WeightingChoice_t, WeightingChoice_vals_t, &
        GAS_PCHB_FullFull_t

    type, extends(EnumBase_t) :: PCHB_HoleAlgorithm_t
        private
        character(9) :: str
    contains
        procedure :: to_str => HoleAlgorithm_to_str
    end type

    type :: PCHB_HoleAlgorithm_vals_t
        type(PCHB_HoleAlgorithm_t) :: &
            FAST_FAST = PCHB_HoleAlgorithm_t(1, 'FAST-FAST'), &
            FULL_FULL = PCHB_HoleAlgorithm_t(2, 'FULL-FULL')
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
        procedure(get_doubles_weight_t), pointer, nopass :: get_weight => null()

    contains
        private
        procedure, public :: init => init_GAS_PCHB_Doubles_t
        procedure, public :: finalize => finalize_GugaAliasSampler_t

        procedure, public :: gen_all_distinct

        procedure :: new_info_table
        procedure :: set_info_entry
        procedure :: get_info_entry

        procedure, public :: calc_pgen => calc_pgen_Base

        procedure, public :: calc_orbital_pgen_contr => calc_orbital_pgen_contr_Base
        procedure, public :: calc_orbital_pgen_contr_start => calc_orbital_pgen_contr_start_Base
        procedure, public :: calc_orbital_pgen_contr_end => calc_orbital_pgen_contr_end_Base

        procedure(calc_hole_pgen_t), deferred :: calc_hole_pgen

        procedure(calculate_pchb_sampler_t), deferred :: setup_pchb_sampler_conditional
        procedure(additional_finalize_t), deferred :: finalize_pchb_sampler
        procedure(GetMemoryDemand_t), public, deferred :: get_memory_demand
    end type

    abstract interface
        subroutine calculate_pchb_sampler_t(this, rank_with_info, ij_weights)
            import :: GAS_PCHB_Doubles_t, dp
            implicit none
            class(GAS_PCHB_Doubles_t), intent(inout) :: this
            integer, intent(in) :: rank_with_info
            real(dp), allocatable, intent(out) :: ij_weights(:, :, :)
        end subroutine

        subroutine additional_finalize_t(this)
            import :: GAS_PCHB_Doubles_t
            implicit none
            class(GAS_PCHB_Doubles_t), intent(inout) :: this
        end subroutine

        function calc_hole_pgen_t(this, nI, ilutI, csf_i, ij, i_sg, a, b) result(pgen)
            import :: GAS_PCHB_Doubles_t, dp, CSF_Info_t, nEl, n_int, GugaBits
            implicit none
            class(GAS_PCHB_Doubles_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
            integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: ij, i_sg, a, b
            real(dp) :: pgen
        end function

        type(ByteSize_t) elemental function GetMemoryDemand_t(this)
            import :: GAS_PCHB_Doubles_t, ByteSize_t
            implicit none
            class(GAS_PCHB_Doubles_t), intent(in) :: this
        end function
    end interface

    type, extends(GAS_PCHB_Doubles_t) :: GAS_PCHB_FastFast_t
        private
        !> The shape is (fused_number_of_particles i <= j, n_supergroup)
        type(AliasSampler_2D_t) :: alias_sampler
        integer, allocatable :: tgtOrbs(:, :)
    contains
        private
        procedure, public :: pickOrbitals => pickOrbitals_FastFast
        procedure :: calc_hole_pgen => calc_hole_pgen_FastFast
        procedure :: setup_pchb_sampler_conditional => setup_pchb_sampler_conditional_FastFast
        procedure :: finalize_pchb_sampler => finalize_pchb_sampler_FastFast
        procedure, public :: get_memory_demand => get_memory_demand_FastFast
    end type

    type, extends(GAS_PCHB_Doubles_t) :: GAS_PCHB_FullFull_t
        private

        !> The shape is (fused_number_of_double_excitations, n_supergroup) -> orbitals
        type(AliasSampler_2D_t) :: A_sampler
        !> The shape is (orbitals, fused_number_of_double_excitations, n_supergroup) -> orbitals
        type(AliasSampler_3D_t) :: B_sampler

        type(UnoccupiedGetter_t) :: unoccupied
    contains
        private
        procedure, public :: pickOrbitals => pickOrbitals_FullFull
        procedure :: calc_hole_pgen => calc_hole_pgen_FullFull
        procedure :: setup_pchb_sampler_conditional => setup_pchb_sampler_conditional_FullFull
        procedure :: finalize_pchb_sampler => finalize_pchb_sampler_FullFull
        procedure, public :: get_memory_demand => get_memory_demand_FullFull
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

        call set_doubles_weight_pointer(weighting_selection, this%get_weight)

        call this%setup_pchb_sampler_conditional(root, ij_weights)
        call particle_allocate_and_init(particle_selection, this%GAS_spec, ij_weights, root, use_lookup, this%particle_selector)

        write(stdout, *) "Finished excitation generator initialization"
    end subroutine

    subroutine finalize_GugaAliasSampler_t(this)
        class(GAS_PCHB_Doubles_t), intent(inout) :: this
        if (allocated(this%info_tables)) then
            ! Yes, we assume, that either all or none are allocated
            call this%all_info_table%shared_dealloc()
            call this%finalize_pchb_sampler()
            call this%particle_selector%finalize()
            deallocate(this%particle_selector, this%info_tables, this%GAS_spec, this%indexer)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
        end if
    end subroutine

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    subroutine calc_orbital_pgen_contr_Base(this, nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
        class(GAS_PCHB_Doubles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt

        integer :: i, j, ij, i_sg
        real(dp) :: p_gen

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))
        ij = fuse_symm_idx(i, j)
        p_gen = this%particle_selector%get_pgen(nI, csf_i, i_sg, i, j) &
                    * this%calc_hole_pgen(nI, ilut, csf_i, ij, i_sg, i, j)

        above_cpt =  p_gen / 2.0_dp
        below_cpt = above_cpt
    end subroutine


    subroutine calc_orbital_pgen_contr_start_Base(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
        class(GAS_PCHB_Doubles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = 'calc_orbital_pgen_contr_start'

        integer :: i, j, ij, i_sg

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
        ij = fuse_symm_idx(i, j)

        orb_pgen = this%particle_selector%get_pgen(nI, csf_i, i_sg, i, j) &
                        * this%calc_hole_pgen(nI, ilut, csf_i, ij, i_sg, i, orb_a)
    end subroutine

    subroutine calc_orbital_pgen_contr_end_Base(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
        class(GAS_PCHB_Doubles_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = "calc_orbital_pgen_contr_end"

        integer :: i, j, ij, i_sg

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
        ij = fuse_symm_idx(i, j)

        orb_pgen = this%particle_selector%get_pgen(nI, csf_i, i_sg, i, j) &
                        * this%calc_hole_pgen(nI, ilut, csf_i, ij, i_sg, orb_a, j)
    end subroutine

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

    pure function HoleAlgorithm_from_str(str) result(res)
        character(*), intent(in) :: str
        type(PCHB_HoleAlgorithm_t) :: res
        routine_name("WeightingChoice_from_str")

        select case(to_upper(str))
        case(algorithm_vals%FAST_FAST%str)
            res = algorithm_vals%FAST_FAST
        case(algorithm_vals%FULL_FULL%str)
            res = algorithm_vals%FULL_FULL
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


    subroutine setup_pchb_sampler_conditional_FastFast(this, rank_with_info, ij_weights)
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
        abMax = fuse_symm_idx(nSpatOrbs, nSpatOrbs)
        ijMax = abMax
        allocate(this%tgtOrbs(2, 0:abMax), source=0)
        do a = 1, nSpatOrbs
            do b = 1, a
                ab = fuse_symm_idx(a, b)
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

        block
            type(ByteSize_t) :: memory_demand
            memory_demand = this%get_memory_demand()
            write(stdout, *) "The number of supergroups is", size(supergroups, 2)
            write(stdout, *) "The excitation generator requires roughly: ", &
                memory_demand%to_human_str(binary=.true.), ' or ', memory_demand%to_human_str()
        end block


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
                ij = fuse_symm_idx(i, j)
                do a = 1, nSpatOrbs
                    do b = a, nSpatOrbs
                        ab = fuse_symm_idx(a, b)
                        call this%set_info_entry(ij, ab, get_info(ijab_Index_t([i, j, a, b])))
                    end do
                end do
            end do
        end do

        do i_sg = 1, size(supergroups, 2)
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w = 0.0_dp
                    ij = fuse_symm_idx(i, j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuse_symm_idx(a, b)
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
    end subroutine setup_pchb_sampler_conditional_FastFast

    elemental function get_memory_demand_FastFast(this) result(res)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(ByteSize_t) :: res
        res = this%alias_sampler%get_memory_demand() &
              + ByteSize_t(size(this%tgtOrbs) * int(c_sizeof(this%tgtOrbs(1, 1)), int64)) &
              + this%all_info_table%get_memory_demand()
    end function

    subroutine finalize_pchb_sampler_FastFast(this)
        class(GAS_PCHB_FastFast_t), intent(inout) :: this
        if (allocated(this%tgtOrbs)) then
            deallocate(this%tgtOrbs)
            call this%alias_sampler%finalize()
        end if
    end subroutine

    subroutine pickOrbitals_FastFast(this, nI, ilut, csf_i, excitInfo, pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), a, b, orbs(2), ij, ab, i_sg
        real(dp) :: pgen_elec, pgen_orbs
        routine_name("pickOrbitals_FastFast")

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
        ij = fuse_symm_idx(src(1), src(2))

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
        ASSERT(pgen .isclose. this%calc_pgen(nI, ilut, csf_i, excitInfo))
    end subroutine pickOrbitals_FastFast


    function calc_hole_pgen_FastFast(this, nI, ilutI, csf_i, ij, i_sg, a, b) result(pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: ij, i_sg, a, b
        real(dp) :: pgen
        unused_var(nI); unused_var(ilutI); unused_var(csf_i)
        pgen = this%alias_sampler%get_prob(ij, i_sg, fuse_symm_idx(a, b))
    end function


    subroutine setup_pchb_sampler_conditional_FullFull(this, rank_with_info, ij_weights)
        class(GAS_PCHB_FullFull_t), intent(inout) :: this
        integer, intent(in) :: rank_with_info
        real(dp), allocatable, intent(out) :: ij_weights(:, :, :)
        integer :: i, j, ij, a, b, ab, i_sg
        integer(int64) :: ijMax, abMax
        real(dp), allocatable :: w_A(:), w_B(:), w_B_spat(:)
        integer, allocatable :: supergroups(:, :)

        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! initialize the mapping ab -> (a, b)
        abMax = fuse_symm_idx(nSpatOrbs, nSpatOrbs)
        ijMax = abMax

        this%unoccupied = UnoccupiedGetter_t(nBasis)

        call this%a_sampler%shared_alloc([int(ijMax), size(supergroups, 2)], nBasis, 'PCHB')
        call this%b_sampler%shared_alloc([nBasis, int(ijMax), size(supergroups, 2)], nBasis, 'PCHB')
        call this%new_info_table(int(ijMax, int64), int(abMax, int64))

        block
            type(ByteSize_t) :: memory_demand
            memory_demand = this%get_memory_demand()
            write(stdout, *) "The number of supergroups is", size(supergroups, 2)
            write(stdout, *) "The excitation generator requires roughly: ", &
                memory_demand%to_human_str(binary=.true.), ' or ', memory_demand%to_human_str()
        end block

        if (iProcIndex_intra == rank_with_info) then
            allocate(w_A(nBasis), w_B(nBasis), w_B_spat(nSpatOrbs), source=0._dp)
            allocate(ij_weights(nSpatOrbs, nSpatOrbs, size(supergroups, 2)), source=0._dp)
        else
            allocate(w_A(0), w_B(0), w_B_spat(0))
            allocate(ij_weights(0, 0, 0))
        end if

        ! The info for double excitations is independent of the super group
        do i = 1, nSpatOrbs
            do j = i, nSpatOrbs
                ij = fuse_symm_idx(i, j)
                do a = 1, nSpatOrbs
                    do b = a, nSpatOrbs
                        ab = fuse_symm_idx(a, b)
                        call this%set_info_entry(ij, ab, get_info(ijab_Index_t([i, j, a, b])))
                    end do
                end do
            end do
        end do

        do i_sg = 1, size(supergroups, 2)
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w_A(:) = 0.0_dp
                    ij = fuse_symm_idx(i, j)
                    do a = 1, nSpatOrbs
                        w_B_spat(:) = 0.0_dp
                        do b = 1, nSpatOrbs
                            if (iProcIndex_intra == rank_with_info) then
                                associate(idx => ijab_Index_t([i, j, a, b]))
                                associate(exc => DistinctDouble_t(idx%normalize()))
                                if (this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) &
                                        .and. symmetry_allowed(exc)) then
                                    w_B_spat(b) = this%get_weight(exc)
                                end if
                                end associate
                                end associate
                            end if
                        end do
                        if (iProcIndex_intra == root) then
                            w_B(1::2) = w_B_spat(:)
                            w_B(2::2) = w_B_spat(:)
                            w_B(2 * a - 1) = 0._dp
                        end if
                        call this%B_sampler%setup_entry(2 * a - 1, IJ, i_sg, root, w_B)
                        if (iProcIndex_intra == root) then
                            call swap(w_B(2 * a - 1), w_B(2 * a))
                        end if
                        call this%B_sampler%setup_entry(2 * a, IJ, i_sg, root, w_B)
                        if (iProcIndex_intra == root) then
                            w_A(2 * A - 1) = sum(w_B_spat)
                            w_A(2 * A) = w_A(2 * A - 1)
                        end if
                    end do
                    call this%A_sampler%setup_entry(IJ, i_sg, root, w_A)
                    if (iProcIndex_intra == root) then
                        ij_weights(i, j, i_sg) = sum(w_A)
                        ij_weights(j, i, i_sg) = ij_weights(i, j, i_sg)
                    end if
                end do
            end do
        end do
    end subroutine setup_pchb_sampler_conditional_FullFull

    elemental function get_memory_demand_FullFull(this) result(res)
        class(GAS_PCHB_FullFull_t), intent(in) :: this
        type(ByteSize_t) :: res
        res = this%A_sampler%get_memory_demand() &
                + this%B_sampler%get_memory_demand() &
                + this%all_info_table%get_memory_demand()
    end function

    subroutine finalize_pchb_sampler_FullFull(this)
        class(GAS_PCHB_FullFull_t), intent(inout) :: this
        call this%A_sampler%finalize()
        call this%B_sampler%finalize()
    end subroutine

    subroutine pickOrbitals_FullFull(this, nI, ilut, csf_i, excitInfo, pgen)
        class(GAS_PCHB_FullFull_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), ij, i_sg
        integer :: unoccupied(nBasis - nEl), spin_tgts(2), tgts(2)
            !! `unoccupied` lists the unoccupied "spin"-orbitals in nI-format.
            !! `spin_tgts` lists the tgts as fake "spin"-orbitals.
            !! `tgts` lists the tgts with spatial orbitals.
        real(dp) :: pgen_elec, pgen_hole
        real(dp) :: renorm_first, p_first(2)
        real(dp) :: renorm_second(2), p_second(2)
        integer(n_int) :: ilut_unoccupied(0 : nIfD)
        routine_name("pickOrbitals_FullFull")

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
        ij = fuse_symm_idx(src(1), src(2))

        picking: block
            integer :: dummy
            call this%unoccupied%get(ilut(0 : nIfD), ilut_unoccupied, unoccupied)
            renorm_first = 1._dp - sum(this%A_sampler%get_prob(ij, i_sg, nI))
            if (do_direct_calculation(renorm_first)) then
                renorm_first = sum(this%A_sampler%get_prob(ij, i_sg, unoccupied))
            end if

            call this%A_sampler%constrained_sample(&
                ij, i_sg, unoccupied, ilut_unoccupied, renorm_first, &
                dummy, spin_tgts(1), p_first(1))
            if (spin_tgts(1) == 0) then
                call make_invalid()
                return
            end if

            renorm_second(1) = 1._dp - sum(this%B_sampler%get_prob(spin_tgts(1), ij, i_sg, nI))
            if (do_direct_calculation(renorm_second(1))) then
                renorm_second(1) = sum(this%B_sampler%get_prob(spin_tgts(1), ij, i_sg, unoccupied))
            end if
            ! Note that p(A | A, IJ) is automatically zero and cannot be drawn
            call this%B_sampler%constrained_sample(&
                spin_tgts(1), ij, i_sg, unoccupied, ilut_unoccupied, renorm_second(1), &
                dummy, spin_tgts(2), p_second(1))
            if (spin_tgts(2) == 0) then
                call make_invalid()
                return
            end if
        end block picking

        tgts(:) = gtID(spin_tgts(:))
        if (tgts(1) == tgts(2)) then
            ASSERT(csf_i%occ_int(tgts(1)) == 0)
            pgen = pgen_elec * 2._dp * p_first(1) * p_second(1)
        else
            ! From now on we assume tgts(1) /= tgts(2)!

            ! We could have picked them the other way round.
            ! Account for that.
            ! The renormalization for the first hole is the same
            p_first(2) = this%A_sampler%constrained_getProb(&
                ij, i_sg, unoccupied, renorm_first, spin_tgts(2))
            renorm_second(2) = 1._dp - sum(this%B_sampler%get_prob(spin_tgts(2), ij, i_sg, nI))
            if (do_direct_calculation(renorm_second(2))) then
                renorm_second(2) = sum(this%B_sampler%get_prob(spin_tgts(2), ij, i_sg, unoccupied))
            end if
            p_second(2) = this%B_sampler%constrained_getProb(&
                spin_tgts(2), ij, i_sg, unoccupied, renorm_second(2), spin_tgts(1))

            associate(holes => (2 - csf_i%occ_int(tgts(:))))
                pgen_hole = sum(p_first(:) * p_second(:) * holes(:) * holes(2:1:-1))
            end associate
            pgen = pgen_elec * pgen_hole
        end if
        if (tgts(1) > tgts(2)) call swap(tgts(1), tgts(2))

        block
            integer :: ab
            ab = fuse_symm_idx(tgts(1), tgts(2))
            call extract_excit_info(this%get_info_entry(ij, ab), excitInfo)
            excitInfo%i_sg_start = i_sg
        end block

        ASSERT(pgen .isclose. this%calc_pgen(nI, ilut, csf_i, excitInfo))
        contains
            subroutine make_invalid()
                excitInfo%valid = .false.
                pgen = 0.0_dp
            end subroutine
    end subroutine pickOrbitals_FullFull

    function calc_hole_pgen_FullFull(this, nI, ilutI, csf_i, ij, i_sg, a, b) result(pgen)
        class(GAS_PCHB_FullFull_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: ij, i_sg, a, b
        real(dp) :: pgen
        debug_function_name("calc_hole_pgen_FullFull")

        integer :: unoccupied(nBasis - nEl)
            !! `unoccupied` lists the unoccupied "spin"-orbitals in nI-format.
        integer(n_int) :: ilut_unoccupied(0 : nIfD)

        unused_var(nI)
        ASSERT(all(csf_i%occ_int([a, b]) < 2))
        call this%unoccupied%get(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)
        if (a == b) then
            if (csf_i%occ_int(a) /= 0) then
                pgen = 0._dp
            else
            associate(a_u => 2 * a - 1, a_d => 2 * a)
                pgen = 2 &
                    * this%a_sampler%constrained_getProb(&
                        ij, i_sg, unoccupied, sum(this%a_sampler%get_prob(ij, i_sg, unoccupied)), a_u) &
                    * this%b_sampler%constrained_getProb(&
                        a_u, ij, i_sg, unoccupied, sum(this%b_sampler%get_prob(a_u, ij, i_sg, unoccupied)), a_d)
            end associate
            end if
        else
        block
            ! We have to account for occupation number,
            ! holes from completely empty orbitals are twice as likely to be picked.
            ! This is achieved by using "fake" spin orbitals.
            ! We also have to account for that we can pick `A, B` or `B, A`,
            ! i.e. reverse the order.
            integer :: spin_tgts(2)
            real(dp) :: p_first(2), p_second(2)
            real(dp) :: renorm_first, renorm_second(2)

            ! If the stepvector value is 0, i.e. fully empty
            ! then the "fake" spin orb index could take both values, u or d.
            ! We decide on d, because it allows a simple ternary logical operation.
            spin_tgts(1) = merge(2 * a, 2 * a - 1, csf_i%stepvector(a) <= 1)
            spin_tgts(2) = merge(2 * b, 2 * b - 1, csf_i%stepvector(b) <= 1)

            renorm_first = sum(this%A_sampler%get_prob(ij, i_sg, unoccupied))

            p_first(1) = this%A_sampler%constrained_getProb(&
                ij, i_sg, unoccupied, renorm_first, spin_tgts(1))
            p_first(2) = this%A_sampler%constrained_getProb(&
                ij, i_sg, unoccupied, renorm_first, spin_tgts(2))

            renorm_second(1) = sum(this%B_sampler%get_prob(spin_tgts(1), ij, i_sg, unoccupied))
            renorm_second(2) = sum(this%B_sampler%get_prob(spin_tgts(2), ij, i_sg, unoccupied))

            p_second(1) = this%B_sampler%constrained_getProb(&
                spin_tgts(1), ij, i_sg, unoccupied, renorm_second(1), spin_tgts(2))
            p_second(2) = this%B_sampler%constrained_getProb(&
                spin_tgts(2), ij, i_sg, unoccupied, renorm_second(2), spin_tgts(1))

            pgen = sum((2 - csf_i%occ_int([a, b])) * p_first(:) * (2 - csf_i%occ_int([b, a])) * p_second(:))
        end block
        end if
    end function


    function calc_pgen_Base(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(GAS_PCHB_Doubles_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: ij, i_sg

        ! i, j, k and l entries have the info about the electrons and
        ! holes of the excitation: E_{ij}E_{kl} is the convention ...
        i_sg = excitInfo%i_sg_start
        ij = fuse_symm_idx(excitInfo%j, excitInfo%l)

        pgen = this%particle_selector%get_pgen(nI, csf_i, i_sg, excitInfo%j, excitInfo%l) &
               * this%calc_hole_pgen(nI, ilutI, csf_i, ij, i_sg, excitInfo%i, excitInfo%k)
    end function



    pure function gen_all_distinct(this, csf_i, sort_fused_key) result(res)
        !! This function returns all possible excitations for the orbital picker
        !!
        !! This means it considers only occupation, not spin recoupling.
        class(GAS_PCHB_Doubles_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        logical, intent(in), optional :: sort_fused_key
            !! The fused index of the excitations should be sorted.
        type(DistinctDouble_t), allocatable :: res(:)
        res = gen_all_distinct_GAS(this%GAS_spec, csf_i, sort_fused_key)
    end function

    pure function gen_all_distinct_GAS(GAS_spec, csf_i, sort_fused_key) result(res)
        !! This function returns all possible excitations for the orbital picker
        !!
        !! This means it considers only occupation, not spin recoupling.
        class(GASSpec_t), intent(in) :: GAS_spec
        type(CSF_Info_t), intent(in) :: csf_i
        logical, intent(in), optional :: sort_fused_key
            !! The fused index of the excitations should be sorted.
        type(DistinctDouble_t), allocatable :: res(:)

        type(buffer_int64_1D_t) :: buffer
        type(DistinctDouble_t) :: exc
        type(ijkl_Index_t) :: idx
        logical :: sort_fused_key_
        integer :: i, j, a, b

        def_default(sort_fused_key_, sort_fused_key, .false.)

        call buffer%init()

        do i = 1, nSpatOrbs
            if (csf_i%occ_int(i) == 0) cycle
            do j = i, nSpatOrbs
                if (csf_i%occ_int(j) == 0) cycle
                if (i == j .and. csf_i%occ_int(i) /= 2) cycle
                do a = 1, nSpatOrbs
                    do b = a, nSpatOrbs
                        if (any(csf_i%stepvector([a, b]) == 3)) cycle
                        if (a == b .and. csf_i%stepvector(a) /= 0) cycle
                        idx = ijkl_Index_t([a, i, b, j])
                        idx = idx%normalize()
                        exc = DistinctDouble_t(idx)
                        if (exc%typ == distinct_doubles%invalid &
                            .or. .not. GAS_spec%is_allowed(exc, GAS_spec%count_per_GAS(csf_i))) then
                            cycle
                        end if
                        call buffer%push_back(exc%fuse())
                    end do
                end do
            end do
        end do

        if (sort_fused_key_) then
            call sort(buffer%buf(:buffer%size()))
            res = fused_to_DistinctDouble(buffer%buf(:buffer%size()))
        else
            res = fused_to_DistinctDouble(buffer%buf(:buffer%size()))
        end if
    end function

    pure function HoleAlgorithm_to_str(this) result(res)
        class(PCHB_HoleAlgorithm_t), intent(in) :: this
        character(:), allocatable :: res
        res = trim(this%str)
    end function

end module
