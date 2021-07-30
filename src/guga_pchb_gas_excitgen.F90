module guga_pchb_gas_excitgen
    use constants, only: n_int, dp, int64, maxExcit, iout, bits_n_int, int32
    ! use util_mod, only: fuseIndex, getSpinIndex, near_zero, intswap, operator(.div.), operator(.implies.), EnumBase_t
    ! use dSFMT_interface, only: genrand_real2_dSFMT
    ! use get_excit, only: make_double, exciteIlut
    ! use SymExcitDataMod, only: pDoubNew, ScratchSize
    ! use excitation_types, only: SingleExc_t, DoubleExc_t, excite
    ! use sltcnd_mod, only: sltcnd_excit
    ! use procedure_pointers, only: generate_single_excit_t
    ! use aliasSampling, only: AliasSampler_2D_t
    ! use UMatCache, only: gtID, numBasisIndices
    ! use FciMCData, only: pSingles, excit_gen_store_type, pParallel, projEDet
    ! use excit_gens_int_weighted, only: pick_biased_elecs
    ! use shared_ragged_array, only: shared_ragged_array_int32_t
    ! use growing_buffers, only: buffer_int32_1D_t
    ! use parallel_neci, only: iProcIndex_intra
    ! use sets_mod, only: complement, operator(.complement.)
    ! use get_excit, only: make_single
    ! use growing_buffers, only: buffer_int_2D_t
    ! use timing_neci, only: timer, set_timer, halt_timer
    !
    ! use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    ! use bit_rep_data, only: NIfTot, nIfD
    ! use bit_reps, only: decode_bit_det
    ! use sort_mod, only: sort
    ! use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt

    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer

    use excitation_generators, only: &
            ExcitationGenerator_t, SingleExcitationGenerator_t, &
            DoubleExcitationGenerator_t, gen_exc_sd, get_pgen_sd, gen_all_excits_sd
    implicit none

    !> The GAS PCHB excitation generator for doubles
    type, extends(DoubleExcitationGenerator_t) :: GAS_doubles_PCHB_ExcGenerator_t
        private
        !> Use a lookup for the supergroup index in global_det_data
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.

        !> The shape is (fused_number_of_double_excitations, n_supergroup)
        type(AliasSampler_2D_t) :: pchb_samplers
        type(shared_array_int64_t) :: all_info_table
        type(shared_array_int64_t), allocatable :: info_tables(:)

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(GASSpec_t), allocatable :: GAS_spec
        integer, allocatable :: tgtOrbs(:, :)
    contains
        private
        procedure, public :: init => GAS_doubles_PCHB_init
        procedure, public :: finalize => GAS_doubles_PCHB_finalize
        procedure, public :: gen_exc => GAS_doubles_PCHB_gen_exc
        procedure, public :: get_pgen => GAS_doubles_PCHB_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_PCHB_gen_all_excits

        procedure :: compute_samplers => GAS_doubles_PCHB_compute_samplers

        procedure :: setup_info_table
        procedure :: get_info
        procedure :: setup_entry_info
    end type

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
        integer :: nBI

        ASSERT(GAS_spec%recoupling())

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(iout, *) 'GUGA GAS PCHB doubles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(iout, *) 'GAS PCHB doubles is using the supergroup lookup'

        write(iout, *) "Allocating PCHB excitation generator objects"
        ! number of spatial orbs
        nBI = numBasisIndices(this%GAS_spec%n_spin_orbs() .div. 2)
        ! initialize the mapping ab -> (a, b)
        abMax = fuseIndex(nBI, nBI)
        allocate(this%tgtOrbs(2, 0 : abMax), source=0)
        do a = 1, nBI
            do b = 1, a
                ab = fuseIndex(a, b)
                this%tgtOrbs(1, ab) = b
                this%tgtOrbs(2, ab) = a
            end do
        end do

        ! setup the alias table
        call this%compute_samplers(nBI)

        write(iout, *) "Finished excitation generator initialization"
    end subroutine GAS_doubles_PCHB_init

    subroutine setup_info_table(this, nEntries, entrySize)
        class(GugaAliasSampler_t) :: this
        integer(int64), intent(in) :: nEntries, entrySize

        integer(int64) :: total_size
        integer(int64) :: iEntry, windowStart, windowEnd

        allocate(this%info_tables(nEntries))
        total_size = nEntries * entrySize
        call this%all_info_table%shared_alloc(total_size)

        windowStart = 1_int64
        do iEntry = 1, nEntries
            windowEnd = windowStart + entrySize - 1

            this%info_tables(iEntry)%ptr => this%all_info_table%ptr(windowStart : windowEnd)
            windowStart = windowStart + entrySize
        end do
    end subroutine setup_info_table

    pure function get_info(this, iEntry, tgt) result(info)
        debug_function_name("get_info")
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry, tgt
        integer(int64) :: info

        ASSERT(associated(this%info_tables(iEntry)%ptr))
        info = this%info_tables(iEntry)%ptr(tgt)
    end function get_info

    subroutine setup_entry_info(this, iEntry, infos)
        class(GugaAliasSampler_t) :: this
        integer, intent(in) :: iEntry
        integer(int64), intent(in) :: infos(:)

        if (iProcIndex_intra == 0) then
            this%info_tables(iEntry)%ptr = infos
        end if

        call this%all_info_table%sync()
    end subroutine setup_entry_info


    !>  @brief
    !>  Deallocate the sampler and the mapping ab -> (a,b)
    subroutine GAS_doubles_PCHB_finalize(this)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this

        call this%pchb_samplers%finalize()
        call this%all_info_table%shared_dealloc()
        deallocate(this%info_tables)
        deallocate(this%tgtOrbs)

        deallocate(this%indexer)

        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine


    !>  @brief
    !>  Given the initial determinant (both as nI and ilut), create a random double
    !>  excitation using the hamiltonian matrix elements as weights
    subroutine GAS_doubles_PCHB_gen_exc(&
                    this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                    ex, tParity, pGen, hel, store, part_type)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        integer :: elecs(2), src(2), sym_prod, ispn, sum_ml, ij
        integer :: orbs(2), ab
        real(dp) :: pGenHoles
        logical :: invalid
        integer :: spin(2), samplerIndex
        integer :: i_sg

        @:unused_var(exFlag, part_type)
        ic = 2
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0.0_dp)
#endif

        ! first, pick two random elecs
        call pick_elec_pair_uniform_guga(nI, src, sym_prod, sum_ml, pgen_elec)
        ASSERT(src(1) < src(2))

        i = gtID(src(1)); j = gtID(src(2))

        ! TODO: put current_stepvector out here
        ! if (i /= j) then
        !     if (current_stepvector(i) == 3) pgen_elec = 2.0_dp * pgen_elec
        !     if (current_stepvector(j) == 3) pgen_elec = 2.0_dp * pgen_elec
        ! end if

        ij = fuseIndex(i, j)

        call guga_pchb_sampler%alias_sampler%sample(ij, ab, pgen_orbs)

        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if(near_zero(pgen_orbs)) then
            excitInfo%valid = .false.
            pgen = 0.0_dp
            ! Yes, print. Those events are signficant enough to be always noted in the output
            write(iout, *) "WARNING: Generated excitation with probability of 0"
            return
        end if

        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = tgtOrbs(:, ab)
        a = orbs(1); b = orbs(2)

        ! check if the picked orbs are a valid choice - if they are the same, match one
        ! occupied orbital or are zero (maybe because there are no allowed picks for
        ! the given source) abort
        ! these are the 'easy' checks for GUGA.. more checks need to be done
        ! to see if it is actually a valid combination..
        if (any(orbs == 0) &
            .or. (current_stepvector(a) == 3)
            .or. (current_stepvector(b) == 3) &
            .or. (a == b .and. current_stepvector(b) /= 0)) then
            excitInfo%valid = .false.
            return
        end if

        ! setup a getInfo functionality in the sampler!
        call extract_excit_info(this%get_info(ij, ab), excitInfo)

        pGen = pgen_elec * pgen_orbs
    end subroutine


    !>  @brief
    !>  Calculate the probability of drawing a given double excitation ex
    !>
    !>  @param[in] ex  2x2 excitation matrix
    !>
    !>  @return pgen  probability of generating this double with the pchb double excitgen
    function GAS_doubles_PCHB_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_get_pgen'

        integer :: ab, ij, nex(2, 2), samplerIndex, i_sg

    end function


    subroutine GAS_doubles_PCHB_compute_samplers(this)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        integer :: i, j, ij, ijMax
        integer :: a, b, ab, abMax
        integer :: ex(2, 2)
        integer(int64) :: memCost
        real(dp), allocatable :: w(:)
        integer, allocatable :: supergroups(:, :)
        integer(int64), allocatable :: excit_info(:), counts(:)
        integer :: i_sg, i_exch
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! number of possible source orbital pairs
        ijMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        abMax = ijMax

        !> n_supergroup * number_of_fused_indices * (bytes_per_sampler)
        memCost = size(supergroups, 2, kind=int64) &
                    * int(ijMax, int64) &
                    * (int(abMax, int64) * 3_int64 * 8_int64)

        write(iout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(iout, *) "The number of supergroups is", size(supergroups, 2)
        write(iout, *) "Generating samplers for PCHB excitation generator"
        write(iout, *) "Depending on the number of supergroups this can take up to 10min."
        call this%pchb_samplers%shared_alloc([ijMax, size(supergroups, 2)], abMax, 'PCHB')
        call this%setup_info_table(ijMax, abMax)
        allocate(w(abMax))

        do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(iout, *) 'Still generating the samplers'
            ! allocate: all samplers have the same size
            do i = 1, nSpatOrbs
                ! TODO: The comment is not correct
                ! as we order a,b, we can assume j <= i
                do j = i, nSpatOrbs
                    w(:) = 0.0_dp
                    ij = fuseIndex(i,j)
                    excit_info = 0_int64
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            if (RandExcitSymLabelProd(&
                                        SpinOrbSymLabel(2 * i), SpinOrbSymLabel(2 * j))
                                /= RandExcitSymLabelProd(&
                                        SpinOrbSymLabel(2 * a), SpinOrbSymLabel(2 * b))) then
                                cycle
                            end if

                            ab = fuseIndex(a, b)
                            if (this%GAS_spec%is_allowed(DoubleExc_t(ex), supergroups(:, i_sg))) then
                                call get_weight_and_info(i, j, a, b, w(ab), excit_info(ab))
                            else
                                w(ab) = 0._dp
                            end if
                        end do
                    end do
                    call this%alias_sampler%setup_entry(ij, i_sg, w)
                    call this%setup_entry_info(ij, excit_info)
                end do
            end do
        end do
    end subroutine


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



    subroutine GAS_doubles_PCHB_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=2)
    end subroutine

    elemental subroutine get_weight_and_info(i, j, a, b, w, info)
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
end module
