#include "macros.h"

module guga_gas_pchb_excitgen
    use constants, only: n_int, dp, int64, maxExcit, stdout, stderr, bits_n_int, int32
    use bit_rep_data, only: IlutBits, GugaBits, nIfTot
    use util_mod, only: fuseIndex, near_zero, stop_all
    use aliasSampling, only: AliasSampler_2D_t
    use UMatCache, only: gtID, numBasisIndices
    use SymExcitDataMod, only: pDoubNew, ScratchSize, OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite

    use parallel_neci, only: iProcIndex_intra

    use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel

    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles, MaxTau

    use SystemData, only: nel, G1, current_stepvector, t_pchb_weighted_singles, &
                          nBasis, nSpatOrbs, ElecPairs, currentOcc_int, &
                          t_analyze_pchb, t_old_pchb, t_exchange_pchb

    use procedure_pointers, only: get_umat_el

    use guga_data, only: tNewDet, ExcitationInformation_t, gen_type, excit_type
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut
    use guga_excitations, only: assign_excitinfo_values_single, &
                                createStochasticExcitation_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, get_guga_integral_contrib_spat, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo
    use guga_procedure_pointers, only: gen_single_excit_guga, gen_double_excit_guga
    use guga_bitrepops, only: identify_excitation, encode_excit_info, extract_excit_info, &
                              contract_2_rdm_ind
    use shared_array, only: shared_array_int64_t

    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer

    use excitation_generators, only: &
            ExcitationGenerator_t, SingleExcitationGenerator_t, &
            DoubleExcitationGenerator_t, gen_exc_sd, get_pgen_sd, gen_all_excits_sd
    implicit none

    private
    public :: GAS_doubles_PCHB_ExcGenerator_t

    !> The GAS PCHB excitation generator for doubles
    type :: GAS_doubles_PCHB_ExcGenerator_t
        private
        !> Use a lookup for the supergroup index in global_det_data
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.

        !> The shape is (fused_number_of_double_excitations, n_supergroup)
        type(AliasSampler_2D_t) :: alias_sampler
        type(shared_array_int64_t) :: all_info_table
        type(shared_array_int64_t), allocatable :: info_tables(:)

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(GASSpec_t), allocatable :: GAS_spec
        integer, allocatable :: tgtOrbs(:, :)
    contains
        private
        procedure, public :: init => GAS_doubles_PCHB_init
        procedure, public :: finalize => GAS_doubles_PCHB_finalize
        procedure, public :: pick_orbitals_double_pchb

        procedure, public :: calc_orb_pgen_guga_pchb_double_excitInfo
        procedure, public :: calc_orbital_pgen_contr_pchb
        procedure, public :: calc_orbital_pgen_contr_start_pchb
        procedure, public :: calc_orbital_pgen_contr_end_pchb

        procedure :: compute_samplers => GAS_doubles_PCHB_compute_samplers
        procedure :: create_info_table
        procedure :: get_info
        procedure :: set_info_entry
    end type

contains

    !>  @brief
    !>  Initialize the pchb excitation generator
    !>
    !>  @details
    !>  This does two things:
    !>  1. setup the lookup table for the mapping ab -> (a,b)
    !>  2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
    subroutine GAS_doubles_PCHB_init(this, GAS_spec, use_lookup, create_lookup)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_init'

        integer :: ab, a, b, abMax

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
        allocate(this%tgtOrbs(2, 0 : abMax), source=0)
        do a = 1, nSpatOrbs
            do b = 1, a
                ab = fuseIndex(a, b)
                this%tgtOrbs(1, ab) = b
                this%tgtOrbs(2, ab) = a
            end do
        end do

        ! setup the alias table
        call this%compute_samplers()

        write(stdout, *) "Finished excitation generator initialization"
    end subroutine GAS_doubles_PCHB_init

    subroutine create_info_table(this, n_entries, entry_size)
        class(GAS_doubles_PCHB_ExcGenerator_t) :: this
        integer(int64), intent(in) :: n_entries, entry_size

        integer(int64) :: total_size
        integer(int64) :: i, left, right

        allocate(this%info_tables(n_entries))
        total_size = n_entries * entry_size
        call this%all_info_table%shared_alloc(total_size)

        left = 1_int64
        do i = 1_int64, n_entries
            right = left + entry_size - 1_int64
            this%info_tables(i)%ptr => this%all_info_table%ptr(left : right)
            left = left + entry_size
        end do
    end subroutine

    pure function get_info(this, iEntry, tgt) result(info)
        debug_function_name("get_info")
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: info

        ASSERT(associated(this%info_tables(iEntry)%ptr))
        info = this%info_tables(iEntry)%ptr(tgt)
    end function get_info

    subroutine set_info_entry(this, iEntry, infos)
        class(GAS_doubles_PCHB_ExcGenerator_t) :: this
        integer, intent(in) :: iEntry
        integer(int64), intent(in) :: infos(:)
        if (iProcIndex_intra == 0) then
            this%info_tables(iEntry)%ptr = infos
        end if
        call this%all_info_table%sync()
    end subroutine set_info_entry


    !>  @brief
    !>  Deallocate the sampler and the mapping ab -> (a,b)
    subroutine GAS_doubles_PCHB_finalize(this)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this

        call this%alias_sampler%finalize()
        call this%all_info_table%shared_dealloc()
        deallocate(this%info_tables)
        deallocate(this%tgtOrbs)

        deallocate(this%indexer)

        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine


    subroutine pick_orbitals_double_pchb(this, ilut, nI, store, excitInfo, pgen)
        debug_function_name("pick_orbitals_double_pchb")
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(excit_gen_store_type), intent(in) :: store
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), sym_prod, sum_ml, i, j, a, b, orbs(2), ij, ab
        integer :: i_sg
        real(dp) :: pgen_elec, pgen_orbs
        unused_var(ilut)

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
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
            if (current_stepvector(i) == 3) pgen_elec = 2.0_dp * pgen_elec
            if (current_stepvector(j) == 3) pgen_elec = 2.0_dp * pgen_elec
        end if

        ! use the sampler for this electron pair -> order of src electrons
        ! does not matter
        ij = fuseIndex(i, j)

        ! get a pair of orbitals using the precomputed weights
        call this%alias_sampler%sample(ij, i_sg, ab, pgen_orbs)

        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if (near_zero(pgen_orbs)) then
            excitInfo%valid = .false.
            pgen = 0.0_dp
            ! Yes, print. Those events are signficant enough to be always noted in the output
            write(stderr, '(A)') "WARNING: Generated excitation with probability of 0"
            return
        endif

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
                .or. (current_stepvector(a) == 3) &
                .or. (current_stepvector(b) == 3) &
                .or. (a == b .and. current_stepvector(b) /= 0)) then
            excitInfo%valid = .false.
            return
        end if

        ! setup a getInfo functionality in the sampler!
        call extract_excit_info(this%get_info(ij, ab), excitInfo)

        pGen = pgen_elec * pgen_orbs
    end subroutine pick_orbitals_double_pchb


    subroutine GAS_doubles_PCHB_compute_samplers(this)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        integer :: i, j, ij, ijMax
        integer :: a, b, ab, abMax
        integer :: ex(2, 2)
        integer(int64) :: memCost
        real(dp), allocatable :: w(:)
        integer, allocatable :: supergroups(:, :)
        integer(int64), allocatable :: excit_info(:)
        integer :: i_sg
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! number of possible source orbital pairs
        ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        abMax = ijMax

        ! number_of_fused_indices * (bytes_per_sampler) * n_supergroup
        memCost = int(ijMax, int64) * (int(abMax, int64) * 3_int64 * 8_int64) &
                  * size(supergroups, 2, kind=int64)

        write(stdout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(stdout, *) "The number of supergroups is", size(supergroups, 2)
        write(stdout, *) "Generating samplers for PCHB excitation generator"
        write(stdout, *) "Depending on the number of supergroups this can take up to 10min."
        call this%alias_sampler%shared_alloc([ijMax, size(supergroups, 2)], abMax, 'PCHB')
        call this%create_info_table(int(ijMax, int64), int(abMax, int64))
        allocate(w(abMax))

        do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(stdout, *) 'Still generating the samplers'
            ! allocate: all samplers have the same size
            do i = 1, nSpatOrbs
                ! TODO: The comment is not correct
                ! as we order a,b, we can assume j <= i
                do j = i, nSpatOrbs
                    w(:) = 0.0_dp
                    excit_info = 0_int64
                    ij = fuseIndex(i,j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            if (RandExcitSymLabelProd(&
                                        SpinOrbSymLabel(2 * i), SpinOrbSymLabel(2 * j)) &
                                /= RandExcitSymLabelProd(&
                                        SpinOrbSymLabel(2 * a), SpinOrbSymLabel(2 * b))) then
                                cycle
                            end if

                            ab = fuseIndex(a, b)
                            if (this%GAS_spec%is_allowed(DoubleExc_t(2*i, 2*a, 2*j, 2*b), supergroups(:, i_sg))) then
                                call get_weight_and_info(i, j, a, b, w(ab), excit_info(ab))
                            end if
                        end do
                    end do
                    call this%alias_sampler%setup_entry(ij, i_sg, w)
                    call this%set_info_entry(ij, excit_info)
                end do
            end do
        end do
    end subroutine

    pure function calc_orb_pgen_guga_pchb_double_excitInfo(this, excitInfo, i_sg) result(pgen)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer, intent(in) :: i_sg
            !! The supergroup index
        real(dp) :: pgen

        integer :: ij, ab
        real(dp) :: p_elec

        ! i, j, k and l entries have the info about the electrons and
        ! holes of the excitation: E_{ij}E_{kl} is the convention ...
        ij = fuseIndex(excitInfo%j, excitInfo%l)
        p_elec = 1.0_dp / real(ElecPairs, dp)
        ab = fuseIndex(excitInfo%i, excitInfo%k)
        pgen = p_elec * this%alias_sampler%get_prob(ij, i_sg, ab)
    end function calc_orb_pgen_guga_pchb_double_excitInfo

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    pure subroutine calc_orbital_pgen_contr_pchb(this, ilut, occ_orbs, i_sg, cpt_a, cpt_b)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: occ_orbs(2)
        integer, intent(in) :: i_sg
            !! The supergroup index
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
        cpt_a = this%alias_sampler%get_prob(ij, i_sg, ij) / 2.0_dp
        ! but since they will get added in later routines i actually have
        ! to divide by 2 here..

        ! and in the PCHB there is no difference between the 2!
        cpt_b = cpt_a
    end subroutine calc_orbital_pgen_contr_pchb

    ! i think it would be better if i 'just' reimplement:
    pure function calc_orbital_pgen_contr_start_pchb(this, occ_orbs, a, i_sg) result(orb_pgen)
        debug_function_name("calc_orbital_pgen_contr_start_pchb")
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: occ_orbs(2), a
        integer, intent(in) :: i_sg
            !! The supergroup index
        real(dp) :: orb_pgen

        integer :: i, j, ij, ab

        ! depending on type (R->L / L->R) a can be > j or < j, but always > i
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ASSERT( i < a )

        ! here i is both electron and hole index!
        ij = fuseIndex(i, j)
        ab = fuseIndex(i, a)
        orb_pgen = this%alias_sampler%get_prob(ij, i_sg, ab)
    end function calc_orbital_pgen_contr_start_pchb

    pure function calc_orbital_pgen_contr_end_pchb(this, occ_orbs, a, i_sg) result(orb_pgen)
        debug_function_name("calc_orbital_pgen_contr_end_pchb")
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: occ_orbs(2), a
        integer, intent(in) :: i_sg
            !! The supergroup index
        real(dp) :: orb_pgen

        integer :: i, j, ij, ab

        i = gtID(occ_orbs(1)); j = gtID(occ_orbs(2))
        ! j is at the same time electron and hole index!

        ! depending on L->R or R->L type a can be > i o r < i, but always < j!
        ASSERT( a < j )

        ! j here is both elec and hole ind!
        ij = fuseIndex(i, j)
        ab = fuseIndex(a, j)
        orb_pgen = this%alias_sampler%get_prob(ij, i_sg, ab)
    end function calc_orbital_pgen_contr_end_pchb


    subroutine get_weight_and_info(i, j, a, b, w, info)
        integer, intent(in) :: i, j, a, b
        real(dp), intent(out) :: w
        integer(int64), intent(out) :: info
        if (i == j) then
            if (a == b) then
                ! here we only have a contribution if
                ! a != i
                if (a < i) then
                    ! _RR_(a) -> ^RR^(i)
                    w = get_pchb_integral_contrib(i, j, a, b, &
                        typ=excit_type%fullstart_stop_alike)
                    info = encode_excit_info(&
                        typ=excit_type%fullstart_stop_alike, &
                        a=a, i=i, b=b, j=j)
                else if (a > i) then
                    ! _LL_(i) > ^LL^(a)
                    w = get_pchb_integral_contrib(i, j, a, b, &
                        typ = excit_type%fullstart_stop_alike)
                    info = encode_excit_info(&
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
                        w = get_pchb_integral_contrib(i, j, a, b, &
                            typ=excit_type%fullstop_raising)
                        info = encode_excit_info(&
                            typ=excit_type%fullstop_raising, &
                            a=a, i=i, b=b, j=j)
                    else if (b > i) then
                        ! _R(a) -> ^RL_(i) -> ^L(b)
                        w = get_pchb_integral_contrib(i, j, a, b, &
                            excit_type%single_overlap_R_to_L)
                        info = encode_excit_info(&
                            typ=excit_type%single_overlap_R_to_L, &
                            a=a, i=i, b=b, j=j)
                    end if
                else if (a > i) then
                    ! since b > a ensured only:
                    ! _LL_(i) -> ^LL(a) -> ^L(b)
                    w = get_pchb_integral_contrib(i, j, a, b, &
                        excit_type%fullstart_lowering)
                    info = encode_excit_info(&
                        typ=excit_type%fullstart_lowering, &
                        a=a, i=i, b=b, j=j)
                end if
            end if
        else if ( i /= j) then
            if (a == b) then
                ! a == i or a == j NOT allowed!
                if (a < i) then
                    ! _RR_(a) -> ^RR(i) -> ^R(j)
                    w = get_pchb_integral_contrib(i, j, a, b,&
                        typ = excit_type%fullstart_raising)
                    info = encode_excit_info(&
                        typ=excit_type%fullstart_raising, &
                        a=a, i=i, b=b, j=j)
                else if (a > i .and. a < j) then
                    ! _L(i) -> ^LR_(a) -> ^R(j)
                    w = get_pchb_integral_contrib(i, j, a, b,&
                        excit_type%single_overlap_L_to_R)
                    info = encode_excit_info(&
                        typ=excit_type%single_overlap_L_to_R, &
                        a=a, i=i, b=b, j=j)
                else if (a > j) then
                    ! _L(i) -> _LL(j) -> ^LL^(a)
                    w = get_pchb_integral_contrib(i, j, a, b,&
                        excit_type%fullstop_lowering)
                    info = encode_excit_info(&
                        typ=excit_type%fullstop_lowering, &
                        a=a, i=i, b=b, j=j)
                end if
            else if (a /= b) then
                ! this is the most general case. a lot of IFs
                if (a < i) then
                    if (b < i) then
                        w = get_pchb_integral_contrib(&
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
                        info = encode_excit_info(&
                            typ = excit_type%double_raising, &
                            a=b, i=j, b=a, j=i)
                    ! else if (b == i) then
                        ! b == i also NOT allowed here,
                        ! since this would correspond to
                        ! a single!
                    else if (b > i .and. b < j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%double_R_to_L_to_R)
                        ! for E_{ai}E_{bj} this would
                        ! correspond to a non-overlap:
                        ! _R(a) -> ^R(i) + _R(b) > ^R(j)
                        ! which are not directly sampled
                        ! However for E_{aj}E_{bj} this is
                        ! _R{a} -> _LR(i) -> ^LR(b) -> ^R(j)
                        info = encode_excit_info(&
                            typ=excit_type%double_R_to_L_to_R, &
                            a=a, i=j, b=b, j=i)
                    else if (b == j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%fullstop_R_to_L)
                        ! here we also have to switch
                        ! indices to E_{aj}E_{bi} to get:
                        ! _R(i) -> _LR(a) -> ^LR^(j)
                        info = encode_excit_info(&
                            typ=excit_type%fullstop_R_to_L, &
                            a=a, i=j, b=b, j=i)

                    else if (b > j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%double_R_to_L)
                        ! here we have to switch to
                        ! E_{aj}E_{bi} to get:
                        ! _R(a) > _LR(i) -> LR^(j) -> ^L(b)
                        info = encode_excit_info(&
                            typ=excit_type%double_R_to_L, &
                            a=a, i=j, b=b, j=i)
                    end if
                else if (a == i) then
                    ! b > i is ensured here since b > a in here!
                    if (b < j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%fullstart_L_to_R)
                        ! here we have to switch again:
                        ! E_{aj}E_{bi}:
                        ! _RL_(i) -> ^LR(b) -> ^R(j)
                        info = encode_excit_info(&
                            typ=excit_type%fullstart_L_to_R, &
                            a=a, i=j, b=b, j=i)
                    else if (b == j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%fullstart_stop_mixed)
                        ! switch again: E_{aj}E_{bi}
                        ! _RL_(i) -> _RL_(j)
                        info = encode_excit_info(&
                            typ=excit_type%fullstart_stop_mixed, &
                            a=a, i=j, b=b, j=i)
                    else if (b > j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%fullstart_R_to_L)
                        ! switch again: E_{aj}E_{bi}
                        ! _RL_(i) -> ^RL(j) -> ^L(b)
                        info = encode_excit_info(&
                            typ=excit_type%fullstart_R_to_L, &
                            a=a, i=j, b=b, j=i)
                    end if
                else if (a > i .and. a < j) then
                    ! b > a still ensured!
                    if (b < j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%double_L_to_R)
                        ! switch: E_{aj}E_{bi}:
                        ! _L(j) -> _RL(a) -> ^LR(b) -> ^R(j)
                        info = encode_excit_info(&
                            typ=excit_type%double_L_to_R, &
                            a=a, i=j, b=b, j=i)

                    else if (b == j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%fullstop_L_to_R)
                        ! switch: E_{aj}E_{bi}
                        ! _L(i) -> _RL(a) -> ^RL^(j)
                        info = encode_excit_info(&
                            typ=excit_type%fullstop_L_to_R, &
                            a=a, i=j, b=b, j=i)
                    else if (b > j) then
                        w = get_pchb_integral_contrib(&
                            i=j, j=i, a=a, b=b, &
                            typ=excit_type%double_L_to_R_to_L)
                        ! switch: E_{aj}E_{bi}
                        ! _L(i) -> _RL(a) - > ^RL(j) -> ^L(b)
                        info = encode_excit_info(&
                            typ=excit_type%double_L_to_R_to_L, &
                            a=a, i=j, b=b, j=i)
                    end if
                ! else if (a == j) then
                    ! a == j also NOT allowed here!
                else if (a > j) then
                    ! b > a > j implied here!
                    w = get_pchb_integral_contrib(&
                        i=i, j=j, a=a, b=b, &
                        typ=excit_type%double_lowering)
                    ! E_{ai}E_{bj} would lead to:
                    ! _L(i) -> LL_(j) -> ^LL(a) -> ^L(b)
                    ! which has the correct sign convention
                    ! as in the Shavitt 81 paper
                    info = encode_excit_info(&
                        typ=excit_type%double_lowering, &
                        a=a, i=i, b=b, j=j)
                end if
            end if
        end if
    end subroutine

    function get_pchb_integral_contrib(i, j, a, b, typ) result(integral)
        ! specialized function to obtain the guga-integral contrib for
        ! the pchb weights
        integer, intent(in) :: a, i, b, j, typ
        real(dp) :: integral
        debug_function_name("get_pchb_integral_contrib")
        logical :: flag_
        real(dp) :: cpt1, cpt2, cpt3, cpt4

        ASSERT(0 < a .and. a <= nSpatOrbs)
        ASSERT(0 < i .and. i <= nSpatOrbs)
        ASSERT(0 < b .and. b <= nSpatOrbs)
        ASSERT(0 < j .and. j <= nSpatOrbs)


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

        case default

            call stop_all(this_routine, "wrong excit-type")

        end select
    end function get_pchb_integral_contrib
end module
