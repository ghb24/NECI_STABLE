#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_doubles_spatorb_fastweighted
    !! precomputed heat bath implementation for GASCI using spatial orbitals
    use constants, only: n_int, dp, int64, maxExcit, stdout
    use util_mod, only: fuseIndex, getSpinIndex, near_zero, swap, &
        operator(.implies.), operator(.isclose.), stop_all
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use SymExcitDataMod, only: pDoubNew, ScratchSize
    use excitation_types, only: DoubleExc_t, excite
    use sltcnd_mod, only: sltcnd_excit
    use aliasSampling, only: AliasSampler_3D_t
    use UMatCache, only: gtID, numBasisIndices
    use FciMCData, only: excit_gen_store_type, projEDet
    use SystemData, only: nEl
    use bit_rep_data, only: NIfTot
    use MPI_wrapper, only: root
    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pchb_doubles_select_particles, only: &
        allocate_and_init, ParticleSelector_t, PCHB_ParticleSelection_t
    use excitation_generators, only: DoubleExcitationGenerator_t
    better_implicit_none

    private
    public :: GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t

    ! there are three pchb_samplers for each supergroup:
    ! 1 - same-spin case
    ! 2 - opp spin case without exchange
    ! 3 - opp spin case with exchange
    integer, parameter :: SAME_SPIN = 1, OPP_SPIN_NO_EXCH = 2, OPP_SPIN_EXCH = 3

    !> The GAS PCHB excitation generator for doubles
    type, extends(DoubleExcitationGenerator_t) :: GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t
        private
        !> Use a lookup for the supergroup index in global_det_data
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.

        !> The shape is (fused_number_of_double_excitations, 3, n_supergroup)
        type(AliasSampler_3D_t) :: pchb_samplers

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(ParticleSelector_t), allocatable :: particle_selector
        class(GASSpec_t), allocatable :: GAS_spec
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
    end type GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t

contains


    !>  @brief
    !>  Initialize the pchb excitation generator
    !>
    !>  @details
    !>  This does two things:
    !>  1. setup the lookup table for the mapping ab -> (a,b)
    !>  2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
    subroutine GAS_doubles_PCHB_init(this, GAS_spec, &
            use_lookup, create_lookup, PCHB_particle_selection)
        class(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_init'

        integer :: ab, a, b, abMax
        integer :: nBI

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'GAS PCHB (RHF) doubles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'GAS PCHB doubles is using the supergroup lookup'

        write(stdout, *) "Allocating PCHB excitation generator objects"
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
        call this%compute_samplers(nBI, PCHB_particle_selection)

        write(stdout, *) "Finished excitation generator initialization"

        ! this is some bias used internally by CreateSingleExcit - not used here
        pDoubNew = 0.0
    end subroutine GAS_doubles_PCHB_init


    !>  @brief
    !>  Deallocate the sampler and the mapping ab -> (a,b)
    subroutine GAS_doubles_PCHB_finalize(this)
        class(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t), intent(inout) :: this

        if (allocated(this%particle_selector)) then
            call this%pchb_samplers%finalize()
            call this%particle_selector%finalize()
            ! Yes we assume that either all or none are allocated.
            deallocate(this%particle_selector, this%tgtOrbs, this%pExch, this%indexer)
            if (this%create_lookup) then
                nullify(lookup_supergroup_indexer)
            end if
        end if
    end subroutine GAS_doubles_PCHB_finalize


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
        class(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_gen_exc'

        integer :: elecs(2), src(2), ij
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

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if


        ! first, pick two random elecs
        call this%particle_selector%draw(nI, ilutI, i_sg, elecs, src, pGen)
        if (src(1) == 0) then
            call invalidate()
            return
        end if


        @:ASSERT(pGen .isclose. this%particle_selector%get_pgen(nI, i_sg, src(1), src(2)))

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
                call swap(spin(1), spin(2))
            else
                samplerIndex = OPP_SPIN_NO_EXCH
                ! adjust pgen
                pGen = pGen * (1.0_dp - this%pExch(ij, i_sg))
            end if
        end if
        ! get a pair of orbitals using the precomputed weights
        call this%pchb_samplers%sample(ij, samplerIndex, i_sg, ab, pGenHoles)
        @:ASSERT(ab /= 0 .implies. (pGenHoles .isclose.  this%pchb_samplers%get_prob(ij, samplerIndex, i_sg, ab)))


        if (ab == 0) then
            invalid = .true.
            orbs = 0
        else
            ! split the index ab (using a table containing mapping ab -> (a,b))
            orbs = this%tgtOrbs(:, ab)
            ! convert orbs to spin-orbs with the same spin
            orbs = 2 * orbs - spin
            @:ASSERT(all(orbs /= 0) .implies. orbs(1) /= orbs(2))
            ! note: nI is an array, not a scalar, so we need two `any` checks below
            invalid = any(orbs == 0) .or. any(orbs(1) == nI) .or. any(orbs(2) == nI)
        end if

        ! unfortunately, there is a super-rare case when, due to floating point error,
        ! an excitation with pGen=0 is created. Those are invalid, too
        if (.not. invalid .and. near_zero(pGenHoles)) then
            invalid = .true.
            ! Yes, print. Those events are signficant enough to be always noted in the output
            write(stdout, *) "WARNING: Generated excitation with probability of 0"
        end if

        if (invalid) then
            ! if 0 is returned, there are no excitations for the chosen elecs
            ! -> return nulldet
            call invalidate()
        else
            ! else, construct the det from the chosen orbs/elecs

            call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

            ilutJ = exciteIlut(ilutI, src, orbs)

            pGen = pGen * pGenHoles

            block
                integer :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
                @:ASSERT(pgen .isclose. this%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2))
            end block

        end if

    contains

        subroutine invalidate()
            nJ = 0
            ilutJ = 0_n_int
            ex(1, 1 : 2) = src
            ex(2, 1 : 2) = orbs
        end subroutine invalidate
    end subroutine GAS_doubles_PCHB_gen_exc


    !>  @brief
    !>  Calculate the probability of drawing a given double excitation ex
    !>
    !>  @param[in] ex  2x2 excitation matrix
    !>
    !>  @return pgen  probability of generating this double with the pchb double excitgen
    function GAS_doubles_PCHB_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_get_pgen'

        integer :: ab, ij, nex(2, 2), samplerIndex, i_sg

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 2)

        nex = gtID(ex(:, : ic))
        ij = fuseIndex(nex(1, 1), nex(1, 2))
        ab = fuseIndex(nex(2, 1), nex(2, 2))
        i_sg = this%indexer%idx_nI(nI)

        pgen = this%particle_selector%get_pgen(nI, i_sg, ex(1, 1), ex(1, 2))

        ! the probability of picking the two electrons: they are chosen uniformly
        ! check which sampler was used
        if (is_beta(ex(1, 1)) .eqv. is_beta(ex(1, 2))) then
            ! same-spin case
            samplerIndex = SAME_SPIN
        else
            ! excitations without spin-exchange OR to the same spatial orb
            if ((is_beta(ex(1, 1)) .eqv. is_beta(ex(2, 1))) .or. (nex(2, 1) == nex(2, 2))) then
                ! opp spin case without exchange
                samplerIndex = OPP_SPIN_NO_EXCH
                pGen = pGen * (1.0_dp - this%pExch(ij, i_sg))
            else
                ! opp spin case with exchange
                samplerIndex = OPP_SPIN_EXCH
                pGen = pGen * this%pExch(ij, i_sg)
            end if
        end if

        pgen = pgen * this%pchb_samplers%get_prob(ij, samplerIndex, i_sg, ab)
    end function GAS_doubles_PCHB_get_pgen


    subroutine GAS_doubles_PCHB_compute_samplers(this, nBI, PCHB_particle_selection)
        !! computes and stores values for the alias sampling table
        class(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nBI
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        integer :: i, j, ij, ijMax
        integer :: a, b, ab, abMax
        integer :: ex(2, 2)
        integer(int64) :: memCost
        real(dp), allocatable :: w(:), pNoExch(:), IJ_weights(:, :, :)
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
        memCost = size(supergroups, 2, kind=int64) &
                    * int(ijMax, int64) &
                    * 3_int64 &
                    * (int(abMax, int64) * 3_int64 * 8_int64)

        write(stdout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(stdout, *) "The number of supergroups is", size(supergroups, 2)
        write(stdout, *) "Generating samplers for PCHB excitation generator"
        write(stdout, *) "Depending on the number of supergroups this can take up to 10min."
        call this%pchb_samplers%shared_alloc([ijMax, 3, size(supergroups, 2)], abMax, 'PCHB_RHF')
        ! weights per pair
        allocate(w(abMax))
        allocate(IJ_weights(nBI * 2, nBI * 2, size(supergroups, 2)), source=0._dp)
        ! initialize the three samplers
        do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(stdout, *) 'Still generating the samplers'
            pNoExch = 1.0_dp - this%pExch(:, i_sg)
            do i_exch = 1, 3
                ! allocate: all samplers have the same size
                do i = 1, nBI
                    ! map i to alpha spin (arbitrary choice)
                    ex(1, 1) = to_spin_orb(i, is_alpha=.true.)
                    ! as we order a,b, we can assume j <= i
                    do j = 1, i
                        w(:) = 0.0_dp
                        ! for samplerIndex == 1, j is alpha, else, j is beta
                        ex(1, 2) = to_spin_orb(j, i_exch == SAME_SPIN)
                        ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
                        ! weights to prepare the sampler
                        do a = 1, nBI
                            ! a is alpha for same-spin (1) and opp spin w/o exchange (2)
                            ex(2, 2) = to_spin_orb(a, any(i_exch == [SAME_SPIN, OPP_SPIN_NO_EXCH]))
                            do b = 1, a
                                ab = fuseIndex(a, b)
                                ! ex(2,:) is in ascending order
                                ! b is alpha for sampe-spin (1) and opp spin w exchange (3)
                                ex(2, 1) = to_spin_orb(b, any(i_exch == [SAME_SPIN, OPP_SPIN_EXCH]))

                                ! exception: for sampler 3, a!=b
                                if (i_exch == OPP_SPIN_EXCH .and. a == b &
                                        .or. any(ex(1, 1) == ex(2, :)) .or. any(ex(1, 2) == ex(2, :)) &
                                        .or. .not. this%GAS_spec%is_allowed(DoubleExc_t(ex), supergroups(:, i_sg))) then
                                    w(ab) = 0._dp
                                else
                                    w(ab) = abs(sltcnd_excit(projEDet(:, 1), DoubleExc_t(ex), .false.))
                                end if
                            end do
                        end do
                        ij = fuseIndex(i, j)

                        call this%pchb_samplers%setup_entry(ij, i_exch, i_sg, root, w)
                        if (i_exch == OPP_SPIN_EXCH) this%pExch(ij, i_sg) = sum(w)
                        if (i_exch == OPP_SPIN_NO_EXCH) pNoExch(ij) = sum(w)

                        associate(I => ex(1, 1), J => ex(1, 2))
                            IJ_weights(I, J, i_sg) = IJ_weights(I, J, i_sg) + sum(w)
                            IJ_weights(J, I, i_sg) = IJ_weights(J, I, i_sg) + sum(w)
                        end associate
                        if (i /= j) then
                            ! sum over alpha and beta of the same orbital
                            if (i_exch == SAME_SPIN) then
                                associate(I => ex(1, 1) - 1, J => ex(1, 2) - 1)
                                    IJ_weights(I, J, i_sg) = IJ_weights(I, J, i_sg) + sum(w)
                                    IJ_weights(J, I, i_sg) = IJ_weights(J, I, i_sg) + sum(w)
                                end associate
                            else
                                associate(I => ex(1, 1) - 1, J => ex(1, 2) + 1)
                                    IJ_weights(I, J, i_sg) = IJ_weights(I, J, i_sg) + sum(w)
                                    IJ_weights(J, I, i_sg) = IJ_weights(J, I, i_sg) + sum(w)
                                end associate
                            end if
                        end if

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


        call allocate_and_init(PCHB_particle_selection, this%GAS_spec, IJ_weights, this%use_lookup, this%particle_selector)

    contains
        elemental function to_spin_orb(orb, is_alpha) result(sorb)
            ! map spatial orbital to the spin orbital matching the current samplerIndex
            ! Input: orb - spatial orbital to be mapped
            ! Output: sorb - corresponding spin orbital
            integer, intent(in) :: orb
            logical, intent(in) :: is_alpha
            integer :: sorb

            sorb = merge(2 * orb, 2 * orb - 1, is_alpha)
        end function to_spin_orb
    end subroutine GAS_doubles_PCHB_compute_samplers


    subroutine GAS_doubles_PCHB_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_PCHB_DoublesSpatOrbFastWeightedExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=2)
    end subroutine GAS_doubles_PCHB_gen_all_excits

end module gasci_pchb_doubles_spatorb_fastweighted
