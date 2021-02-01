#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

! The main idea of the precomputed Heat bath sampling (PCHB) is taken from
!    J. Li, M. Otten, A. A. Holmes, S. Sharma, and C. J. Umrigar, J. Comput. Phys. 149, 214110 (2018).
! and described there.
! The main "ingredient" are precomputed probability distributions p(ab | ij) to draw a, b holes
! when i, j electrons were chosen.
! This requires #{i, j | i < j} probability distributions.
!
! The improved version to use spatial orbital indices to save memory is described in
!    Guther K. et al., J. Chem. Phys. 153, 034107 (2020);
! and described there.
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx.
! This gives #{i, j | i < j} * 3 probability distributions.
!
! The generalization to GAS spaces is not yet published.
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and GAS constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in gasci_supergroup_index.f90

module gasci_pchb
    use constants, only: n_int, dp, int64, maxExcit, iout, bits_n_int, int32
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==), operator(/=), alpha, beta
    use util_mod, only: fuseIndex, getSpinIndex, near_zero, intswap, operator(.div.), operator(.implies.), EnumBase_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use SymExcitDataMod, only: pDoubNew, ScratchSize
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite
    use sltcnd_mod, only: sltcnd_excit
    use procedure_pointers, only: generate_single_excit_t
    use aliasSampling, only: AliasSampler_3D_t
    use UMatCache, only: gtID, numBasisIndices
    use FciMCData, only: pSingles, excit_gen_store_type, pParallel, projEDet
    use excit_gens_int_weighted, only: pick_biased_elecs
    use shared_ragged_array, only: shared_ragged_array_int32_t
    use growing_buffers, only: buffer_int32_1D_t
    use parallel_neci, only: iProcIndex_intra
    use sets_mod, only: complement, operator(.complement.)
    use get_excit, only: make_single
    use growing_buffers, only: buffer_int_2D_t
    use timing_neci, only: timer, set_timer, halt_timer

    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs, tGASSpinRecoupling
    use bit_rep_data, only: NIfTot, nIfD
    use bit_reps, only: decode_bit_det
    use sort_mod, only: sort
    use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt

    use gasci, only: GASSpec_t
    use gasci_general, only: GAS_singles_heat_bath_ExcGen_t
    use gasci_util, only: get_available_singles, get_available_doubles
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use exc_gen_class_wrappers, only: UniformSingles_t

    use excitation_generators, only: &
            ExcitationGenerator_t, SingleExcitationGenerator_t, &
            DoubleExcitationGenerator_t, gen_exc_sd, get_pgen_sd, gen_all_excits_sd
    implicit none

    private
    public :: GAS_PCHB_ExcGenerator_t, use_supergroup_lookup, GAS_doubles_PCHB_ExcGenerator_t, &
        possible_GAS_singles, GAS_PCHB_singles_generator

    logical :: use_supergroup_lookup = .true.

    ! there are three pchb_samplers for each supergroup:
    ! 1 - same-spin case
    ! 2 - opp spin case without exchange
    ! 3 - opp spin case with exchange
    integer, parameter :: SAME_SPIN = 1, OPP_SPIN_NO_EXCH = 2, OPP_SPIN_EXCH = 3

    !> The precomputed GAS uniform excitation generator
    type, extends(SingleExcitationGenerator_t) :: GAS_singles_PC_uniform_ExcGenerator_t
        private
        ! allowed_holes(:, src, i_sg)
        ! is a bitmask that returns for a given supergroup `i_sg` and `src`
        ! the GAS allowed holes.
        integer(n_int), allocatable :: allowed_holes(:, :, :)
        type(GASSpec_t) :: GAS_spec
        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        !> Use a lookup for the supergroup index in global_det_data.
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.
    contains
        private
        procedure, public :: init => GAS_singles_uniform_init
        procedure, public :: finalize => GAS_singles_uniform_finalize

        !> Get the GAS allowed holes for a given determinant and a chosen particle.
        procedure, public :: get_possible_holes => GAS_singles_uniform_possible_holes
        procedure, public :: gen_exc => GAS_singles_uniform_gen_exc
        procedure, public :: get_pgen => GAS_singles_uniform_get_pgen
        procedure, public :: gen_all_excits => GAS_singles_uniform_gen_all_excits
    end type


    type, extends(SingleExcitationGenerator_t) :: GAS_singles_DiscardingGenerator_t
        private
        type(UniformSingles_t) :: FCI_singles_generator
        type(GASSpec_t) :: GAS_spec
    contains
        private
        procedure, public :: finalize => GAS_discarding_singles_finalize
        procedure, public :: gen_exc => GAS_discarding_singles_gen_exc
        procedure, public :: get_pgen => GAS_discarding_singles_get_pgen
        procedure, public :: gen_all_excits => GAS_discarding_singles_gen_all_excits
    end type


    type, extends(EnumBase_t) :: GAS_used_singles_t
    end type

    type :: possible_GAS_singles_t
        type(GAS_used_singles_t) :: &
            ON_FLY_HEAT_BATH = GAS_used_singles_t(1), &
            DISCARDING_UNIFORM = GAS_used_singles_t(2), &
            PC_UNIFORM = GAS_used_singles_t(3)
    end type

    type(possible_GAS_singles_t), parameter :: possible_GAS_singles = possible_GAS_singles_t()

    type(GAS_used_singles_t) :: GAS_PCHB_singles_generator = possible_GAS_singles%PC_UNIFORM

    !> The GAS PCHB excitation generator for doubleles
    type, extends(DoubleExcitationGenerator_t) :: GAS_doubles_PCHB_ExcGenerator_t
        private
        !> Use a lookup for the supergroup index in global_det_data
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.

        !> The shape is (fused_number_of_double_excitations, 3, n_supergroup)
        type(AliasSampler_3D_t) :: pchb_samplers


        type(SuperGroupIndexer_t), pointer :: indexer => null()
        type(GASSpec_t) :: GAS_spec
        real(dp), allocatable :: pExch(:, :)
        integer, allocatable :: tgtOrbs(:, :)
    contains
        private
        procedure, public :: init => GAS_doubles_PCHB_init
        procedure, public :: finalize => GAS_doubles_PCHB_finalize
        procedure, public :: gen_exc => GAS_doubles_PCHB_gen_exc
        procedure, public :: get_pgen => GAS_doubles_PCHB_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_PCHB_gen_all_excits

        procedure :: compute_samplers => GAS_doubles_PCHB_compute_samplers
    end type

    type, extends(ExcitationGenerator_t) :: GAS_PCHB_ExcGenerator_t
        private
        type(GAS_doubles_PCHB_ExcGenerator_t) :: doubles_generator
        ! NOTE: Change into class(SingleExcitationGenerator_t), allocatable
        !   if you want to change singles_generators at runtime.
        class(SingleExcitationGenerator_t), allocatable :: singles_generator
    contains
        private
        procedure, public :: init => GAS_PCHB_init
        procedure, public :: finalize => GAS_PCHB_finalize
        procedure, public :: gen_exc => GAS_PCHB_gen_exc
        procedure, public :: get_pgen => GAS_PCHB_get_pgen
        procedure, public :: gen_all_excits => GAS_PCHB_gen_all_excits
    end type

    interface GAS_singles_DiscardingGenerator_t
        module procedure construct_GAS_singles_DiscardingGenerator_t
    end interface

contains

    subroutine GAS_singles_uniform_init(this, GAS_spec, use_lookup, create_lookup)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this
        type(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        integer, allocatable :: supergroups(:, :)
        character(*), parameter :: this_routine = 'GAS_singles_uniform_init'

        integer :: i_sg, src, tgt

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec))
        this%create_lookup = create_lookup
        if (create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(iout, *) 'GAS singles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        this%use_lookup = use_lookup
        if (use_lookup) write(iout, *) 'GAS singles is using the supergroup lookup'
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        allocate(this%allowed_holes(0 : nIfD, this%GAS_spec%n_spin_orbs(), size(supergroups, 2)), source=0_n_int)

        ! Find for each supergroup (i_sg)
        ! the allowed holes `tgt` for a given `src`.
        do i_sg = 1, size(supergroups, 2)
            ! Note that the loop cannot take a symmetric shape.
            ! It may be, that (1 -> 5) is allowed, but (5 -> 1) is not.
            do src = 1, this%GAS_spec%n_spin_orbs()
                do tgt = 1, this%GAS_spec%n_spin_orbs()
                    if (src /= tgt &
                            .and. calc_spin_raw(src) == calc_spin_raw(tgt) &
                            .and. this%GAS_spec%is_allowed(SingleExc_t(src, tgt), supergroups(:, i_sg))) then
                        call my_set_orb(this%allowed_holes(:, src, i_sg), tgt)
                    end if
                end do
            end do
        end do

        contains

            ! Ugly Fortran syntax rules forces us to not use the macro.
            ! Even associate would not help here
            ! https://stackoverflow.com/questions/65734764/non-one-indexed-array-from-associate
            pure subroutine my_set_orb(ilut, orb)
                integer(n_int), intent(inout) :: ilut(0 : nIfD)
                integer, intent(in) :: orb
                set_orb(ilut, orb)
            end subroutine
    end subroutine


    subroutine GAS_singles_uniform_finalize(this)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this

        deallocate(this%allowed_holes)
        deallocate(this%indexer)

        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine

    !> @brief
    !> For a determinant nI and a spin orbital src return
    !>  the GAS allowed orbitals with the same spin as src which are not occupied in nI.
    function GAS_singles_uniform_possible_holes(this, nI, ilutI, src, use_lookup, store) result(unoccupied)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel), src
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
        logical, intent(in) :: use_lookup
        type(excit_gen_store_type), optional, intent(in) :: store
        integer, allocatable :: unoccupied(:)
        integer(n_int) :: ilut_unoccupied(0 : nIfD)
        integer :: i_sg
        character(*), parameter :: this_routine = 'GAS_PC_possible_holes'

        @:ASSERT(use_lookup .implies. present(store))
        if (use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        ilut_unoccupied = iand(this%allowed_holes(:, src, i_sg), not(ilutI))
        allocate(unoccupied(sum(popcnt(ilut_unoccupied))))
        call decode_bit_det(unoccupied, ilut_unoccupied)
    end function


    !> @brief
    !> This is the uniform singles excitation generator which uses precomputed indices
    !> to generate only GAS allowed excitations.
    subroutine GAS_singles_uniform_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        integer :: elec, src, tgt
        integer, allocatable :: unoccupied(:)

        @:unused_var(exFlag, part_type)
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif
        ic = 1

        elec = int(genrand_real2_dSFMT() * nel) + 1
        src = nI(elec)

        unoccupied = this%get_possible_holes(nI, ilutI, src, this%use_lookup, store)

        ! NOTE: this is actually possible for some systems.
        if (size(unoccupied) == 0) then
            pgen = 1._dp / nEl
            nJ(1) = 0
            ilutJ = 0_n_int
            return
        end if

        ! NOTE: The `tgt` could be drawn from `unoccupied` with weights according
        !  to the matrix elements < nI | H |  E_{src}^{tgt} nI >.
        !  Probably not worth it and I am too lazy to implement it now.
        tgt = unoccupied(int(genrand_real2_dSFMT() * size(unoccupied)) + 1)
        pgen = 1._dp / (nEl * size(unoccupied))

        call make_single(nI, nJ, elec, tgt, ex, tParity)
        ilutJ = ilutI
        clr_orb(ilutJ, src)
        set_orb(ilutJ, tgt)
    end subroutine


    function GAS_singles_uniform_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_PC_get_pgen'

        integer :: src
        integer, allocatable :: unoccupied(:)

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 1)
        src = ex(1, 1)
        unoccupied = this%get_possible_holes(nI, ilutI, src, use_lookup=.false.)
        @:ASSERT(size(unoccupied) > 0)
        pgen = 1._dp / (nEl * size(unoccupied))
    end function


    subroutine GAS_singles_uniform_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_singles_PC_uniform_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        integer, allocatable :: singles(:, :)
        integer :: i

        singles = get_available_singles(this%GAS_spec, nI)

        n_excits = size(singles, 2)
        allocate(det_list(0:niftot, n_excits))
        do i = 1, size(singles, 2)
            call EncodeBitDet(singles(:, i), det_list(:, i))
        end do

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine

    pure function construct_GAS_singles_DiscardingGenerator_t(GAS_spec) result(res)
        type(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_singles_DiscardingGenerator_t) :: res
        res%GAS_spec = GAS_spec
        res%FCI_singles_generator = UniformSingles_t()
    end function

    subroutine GAS_discarding_singles_finalize(this)
        class(GAS_singles_DiscardingGenerator_t), intent(inout) :: this
        call this%FCI_singles_generator%finalize()
    end subroutine

    subroutine GAS_discarding_singles_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ex, tParity, pGen, hel, store, part_type)
        class(GAS_singles_DiscardingGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'GAS_discarding_singles_gen_exc'

        integer :: src_copy(maxExcit)

        ASSERT(this%GAS_spec%contains_det(nI))

        call this%FCI_singles_generator%gen_exc(&
                    nI, ilutI, nJ, ilutJ, exFlag, ic, &
                    ex, tParity, pGen, hel, store, part_type)
        if (nJ(1) /= 0) then
            if (.not. this%GAS_spec%contains_det(nJ)) then
                src_copy(:ic) = ex(1, :ic)
                call sort(src_copy)
                ex(1, :ic) = src_copy(:ic)
                ex(2, :ic) = ex(2, :ic)
                nJ(1) = 0
                ilutJ = 0_n_int
            end if
        end if
    end subroutine

    function GAS_discarding_singles_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_singles_DiscardingGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        pgen = this%FCI_singles_generator%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
    end function


    subroutine GAS_discarding_singles_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_singles_DiscardingGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        integer, allocatable :: singles(:, :)
        integer :: i

        singles = get_available_singles(this%GAS_spec, nI)

        n_excits = size(singles, 2)
        allocate(det_list(0:niftot, n_excits))
        do i = 1, size(singles, 2)
            call EncodeBitDet(singles(:, i), det_list(:, i))
        end do

        call sort(det_list, ilut_lt, ilut_gt)
    end subroutine


    !>  @brief
    !>  Initialize the pchb excitation generator
    !>
    !>  @details
    !>  This does two things:
    !>  1. setup the lookup table for the mapping ab -> (a,b)
    !>  2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
    subroutine GAS_doubles_PCHB_init(this, GAS_spec, use_lookup, create_lookup, recoupling)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        type(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup, recoupling
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_init'

        integer :: ab, a, b, abMax
        integer :: nBI

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec) )
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(iout, *) 'GAS PCHB doubles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(iout, *) 'GAS PCHB doubles is using the supergroup lookup'

        write(iout, *) "Allocating PCHB excitation generator objects"
        ! number of spatial orbs
        nBI = numBasisIndices(this%GAS_spec%n_spin_orbs())
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
        call this%compute_samplers(nBI, recoupling)

        write(iout, *) "Finished excitation generator initialization"

        ! this is some bias used internally by CreateSingleExcit - not used here
        pDoubNew = 0.0
    end subroutine GAS_doubles_PCHB_init


    !>  @brief
    !>  Deallocate the sampler and the mapping ab -> (a,b)
    subroutine GAS_doubles_PCHB_finalize(this)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this

        call this%pchb_samplers%finalize()
        deallocate(this%tgtOrbs)
        deallocate(this%pExch)

        deallocate(this%indexer)

        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine


    !>  @brief
    !>  Given the initial determinant (both as nI and ilut), create a random double
    !>  excitation using the hamiltonian matrix elements as weights
    !>
    !> @param[in] nI  determinant to excite from
    !> @param[in] elec_map  map to translate electron picks to orbitals
    !> @param[in] ilut  determinant to excite from in ilut format
    !> @param[out] nJ  on return, excited determinant
    !> @param[out] excitMat  on return, excitation matrix nI -> nJ
    !> @param[out] tParity  on return, the parity of the excitation nI -> nJ
    !> @param[out] pGen  on return, the probability of generating the excitation nI -> nJ
    !> @param[in] idet Optional index of determinant in the CurrentDets array.
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
        call pick_biased_elecs(nI, elecs, src, sym_prod, ispn, sum_ml, pGen)
        if (src(1) > src(2)) call intswap(src(1), src(2))

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        invalid = .false.
        ! use the sampler for this electron pair -> order of src electrons does not matter
        ij = fuseIndex(gtID(src(1)), gtID(src(2)))
        ! the spin of the electrons: 0 - alpha, 1 - beta
        spin = getSpinIndex(src)
        ! determine type of spin-excitation: same-spin, opp spin w exchange, opp spin w/o exchange
        if (spin(1) == spin(2)) then
            ! same spin
            samplerIndex = SAME_SPIN
        else
            ! else, pick exchange with...some ij-spin bias
            if (genrand_real2_dSFMT() < this%pExch(ij, i_sg)) then
                samplerIndex = OPP_SPIN_EXCH
                ! adjust pgen
                pGen = pGen * this%pExch(ij, i_sg)
                ! the spins of the target are the opposite of the source spins
                call intswap(spin(1), spin(2))
            else
                samplerIndex = OPP_SPIN_NO_EXCH
                ! adjust pgen
                pGen = pGen * (1.0_dp - this%pExch(ij, i_sg))
            end if
        end if
        ! get a pair of orbitals using the precomputed weights
        call this%pchb_samplers%sample(ij, samplerIndex, int(i_sg), ab, pGenHoles)
        ! split the index ab (using a table containing mapping ab -> (a,b))
        orbs = this%tgtOrbs(:, ab)
        ! convert orbs to spin-orbs with the same spin
        orbs = 2 * orbs - spin

        ! check if the picked orbs are a valid choice - if they are the same, match one
        ! occupied orbital or are zero (maybe because there are no allowed picks for
        ! the given source) abort
        invalid = (any(orbs == 0) .or. any(orbs(1) == nI) &
                   .or. any(orbs(2) == nI)) .or. (orbs(1) == orbs(2))
        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if (near_zero(pGenHoles)) then
            invalid = .true.
            ! Yes, print. Those events are signficant enough to be always noted in the output
            print *, "WARNING: Generated excitation with probability of 0"
        end if

        pGen = pGen * pGenHoles
        if (invalid) then
            ! if 0 is returned, there are no excitations for the chosen elecs
            ! -> return nulldet
            nJ = 0
            ilutJ = 0_n_int
            ex(2, 1 : 2) = orbs
            ex(1, 1 : 2) = src
        else
            ! else, construct the det from the chosen orbs/elecs

            call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

            ilutJ = exciteIlut(ilutI, src, orbs)
        end if
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

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 2)

        i_sg = this%indexer%idx_nI(nI)
        ! spatial orbitals of the excitation
        nex = gtID(ex(:, : 2))
        ij = fuseIndex(nex(1, 1), nex(1, 2))
        ! the probability of picking the two electrons: they are chosen uniformly
        ! check which sampler was used
        if (is_beta(ex(1, 1)) .eqv. is_beta(ex(1, 2))) then
            pgen = pParallel / par_elec_pairs
            ! same-spin case
            samplerIndex = SAME_SPIN
        else
            pgen = (1.0_dp - pParallel) / AB_elec_pairs
            ! excitations without spin-exchange OR to the same spatial orb
            if ((is_beta(ex(1, 1)) .eqv. is_beta(ex(2, 1))) .or. (nex(2, 1) == nex(2, 2))) then
                ! opp spin case without exchange
                samplerIndex = OPP_SPIN_NO_EXCH
                pgen = pgen * (1 - this%pExch(ij, i_sg))
            else
                ! opp spin case with exchange
                samplerIndex = OPP_SPIN_EXCH
                pgen = pgen * this%pExch(ij, i_sg)
            end if
        end if

        ! look up the probability for this excitation in the sampler
        ab = fuseIndex(nex(2, 1), nex(2, 2))
        pgen = pgen * this%pchb_samplers%get_prob(ij, samplerIndex, i_sg, ab)
    end function


    subroutine GAS_doubles_PCHB_compute_samplers(this, nBI, recoupling)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nBI
        logical, intent(in) :: recoupling
        integer :: i, j, ij, ijMax
        integer :: a, b, ab, abMax
        integer :: ex(2, 2)
        integer(int64) :: memCost
        real(dp), allocatable :: w(:), pNoExch(:)
        integer, allocatable :: supergroups(:, :)
        integer :: i_sg, i_exch
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! number of possible source orbital pairs
        ijMax = fuseIndex(nBI, nBI)
        abMax = ijMax
        ! allocate the bias for picking an exchange excitation
        allocate(this%pExch(ijMax, size(supergroups, 2)), source=0.0_dp)
        ! temporary storage for the unnormalized prob of not picking an exchange excitation

        !> n_supergroup * number_of_fused_indices * 3 * (bytes_per_sampler)
        memCost = size(supergroups, 2) * ijMax * 3 * (abMax * 3 * 8)

        call this%pchb_samplers%shared_alloc([ijMax, 3, size(supergroups, 2)], abMax)
        write(iout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(iout, *) "The number of supergroups is", size(supergroups, 2)
        write(iout, *) "Generating samplers for PCHB excitation generator"
        write(iout, *) "Depending on the number of supergroups this can take up to 10min."
        ! weights per pair
        allocate(w(abMax))
        ! initialize the three samplers
        do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(iout, *) 'Still generating the samplers'
            pNoExch = 1.0_dp - this%pExch(:, i_sg)
            do i_exch = 1, 3
                ! allocate: all samplers have the same size
                do i = 1, nBI
                    ! map i to alpha spin (arbitrary choice)
                    ex(1, 1) = 2 * i
                    ! as we order a,b, we can assume j <= i
                    do j = 1, i
                        w(:) = 0.0_dp
                        ! for samplerIndex == 1, j is alpha, else, j is beta
                        ex(1, 2) = map_orb(j, [SAME_SPIN])
                        ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
                        ! weights to prepare the sampler
                        do a = 1, nBI
                            ! a is alpha for same-spin (1) and opp spin w/o exchange (2)
                            ex(2, 2) = map_orb(a, [SAME_SPIN, OPP_SPIN_NO_EXCH])
                            do b = 1, a
                                ! exception: for sampler 3, a!=b
                                if (i_exch == OPP_SPIN_EXCH .and. a == b) cycle
                                ab = fuseIndex(a, b)
                                ! ex(2,:) is in ascending order
                                ! b is alpha for sampe-spin (1) and opp spin w exchange (3)
                                ex(2, 1) = map_orb(b, [SAME_SPIN, OPP_SPIN_EXCH])
                                ! use the actual matrix elements as weights
                                if (this%GAS_spec%is_allowed(DoubleExc_t(ex), supergroups(:, i_sg), recoupling)) then
                                    w(ab) = abs(sltcnd_excit(projEDet(:, 1), DoubleExc_t(ex), .false.))
                                else
                                    w(ab) = 0._dp
                                end if
                            end do
                        end do
                        ij = fuseIndex(i, j)
                        call this%pchb_samplers%setup_entry(ij, i_exch, i_sg, w)
                        if (i_exch == OPP_SPIN_EXCH) this%pExch(ij, i_sg) = sum(w)
                        if (i_exch == OPP_SPIN_NO_EXCH) pNoExch(ij) = sum(w)
                    end do
                end do
            end do
            ! normalize the exchange bias (where normalizable)
            where (near_zero(this%pExch(:, i_sg) + pNoExch))
                this%pExch(:, i_sg) = 0._dp
            else where
                this%pExch(:, i_sg) = this%pExch(:, i_sg) / (this%pExch(:, i_sg) + pNoExch)
            end where
        end do
    contains
        function map_orb(orb, alphaSamplers) result(sorb)
            ! map spatial orbital to the spin orbital matching the current samplerIndex
            ! Input: orb - spatial orbital to be mapped
            !        alphaSamplers - list of samplerIndex values for which the mapping shall be to alpha
            ! Output: sorb - corresponding spin orbital
            integer, intent(in) :: orb
            integer, intent(in) :: alphaSamplers(:)
            integer :: sorb

            if (any(i_exch == alphaSamplers)) then
                sorb = 2 * orb
            else
                sorb = 2 * orb - 1
            end if
        end function map_orb
    end subroutine


    subroutine GAS_doubles_PCHB_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_doubles_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        integer, allocatable :: doubles(:, :)
        integer :: i

        doubles = get_available_doubles(this%GAS_spec, nI)
        n_excits = size(doubles, 2)
        allocate(det_list(0:niftot, n_excits))
        do i = 1, size(doubles, 2)
            call EncodeBitDet(doubles(:, i), det_list(:, i))
        end do
        call sort(det_list, ilut_lt, ilut_gt)
    end subroutine


    !> @brief
    !> Initialize the PCHB excitation generator. The doubles generator is implemented
    !>  and fixed. A singles generator is required.
    !>
    !>  @param[in] GAS_spec The GAS specifications for the excitation generator.
    !>  @param[in] use_lookup Use a lookup for the supergroup indexing.
    !>  @param[in] recoupling Allow double excitations that change the
    !>                  spin projection per GAS space.
    !>  @param[in] singles_generator A **fully** initialized singles_generator
    !>                  whose cleanup happens outside. Has to be a target.
    subroutine GAS_PCHB_init(this, GAS_spec, use_lookup, create_lookup, recoupling, used_singles_generator)
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
        type(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup, recoupling
        type(GAS_used_singles_t), intent(in) :: used_singles_generator

        if (used_singles_generator == possible_GAS_singles%DISCARDING_UNIFORM) then
            write(iout, *) 'GAS discarding singles activated'
            allocate(this%singles_generator, source=GAS_singles_DiscardingGenerator_t(GAS_spec))
        else if (used_singles_generator == possible_GAS_singles%PC_UNIFORM) then
            write(iout, *) 'GAS precomputed singles activated'
            allocate(GAS_singles_PC_uniform_ExcGenerator_t :: this%singles_generator)
            select type(generator => this%singles_generator)
            type is(GAS_singles_PC_uniform_ExcGenerator_t)
                ! NOTE: only one of the excitation generators should manage the
                !   supergroup lookup!
                call generator%init(GAS_spec, use_lookup, create_lookup=.false.)
            end select
        else if (used_singles_generator == possible_GAS_singles%ON_FLY_HEAT_BATH) then
            write(iout, *) 'GAS heat bath on the fly singles activated'
            allocate(this%singles_generator, source=GAS_singles_heat_bath_ExcGen_t(GAS_spec))
        end if

        call this%doubles_generator%init(GAS_spec, use_lookup, create_lookup, recoupling)
    end subroutine


    subroutine GAS_PCHB_finalize(this)
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
        call this%doubles_generator%finalize()
        call this%singles_generator%finalize()
        deallocate(this%singles_generator)
    end subroutine


    subroutine GAS_PCHB_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                ex, tParity, pGen, hel, store, part_type)
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        call gen_exc_sd(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                        ex, tParity, pGen, hel, store, part_type, &
                        this%singles_generator, this%doubles_generator)
    end subroutine


    real(dp) function GAS_PCHB_get_pgen(&
            this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
        class(GAS_PCHB_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        GAS_PCHB_get_pgen = get_pgen_sd(&
                        nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2, &
                        this%singles_generator, this%doubles_generator)
    end function


    subroutine GAS_PCHB_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_PCHB_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)
        call gen_all_excits_sd(nI, n_excits, det_list, &
                               this%singles_generator, this%doubles_generator)
    end subroutine
end module gasci_pchb
