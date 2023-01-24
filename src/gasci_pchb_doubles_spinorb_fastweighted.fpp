#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_doubles_spinorb_fastweighted
    !! spin-orbital-resolved GASCI-PCHB using the "fast weighted" scheme

    ! @jph reformulate as notes instead of self-ramblings when completed and I
    !    know what's going on
    !==========================================================================!
    ! key point is I can get rid of this line in gasci_pchb
    ! integer, parameter :: SAME_SPIN = 1, OPP_SPIN_NO_EXCH = 2, OPP_SPIN_EXCH = 3
    !   and any loops over `samplerIndex`
    ! which simplifies the code, purportedly...
    ! heavy lifting will be in this file, but I will need a pchb_uhf_excitgen too
    ! not sure where I will put it but I will need to also specify to use this
    ! excit gen when both PCHB and UHF are selected as input
    ! also very important: unit tests
    ! unit_tests/gasci/test_gasci_general_pchb.F90
    ! unit_tests/pcpp_excitgen/test_pchb_excitgen.F90
    ! unit_tests/pcpp_excitgen/test_aliasTables.F90 (don't need to change)
    !==========================================================================!
    use constants, only: dp, n_int, maxExcit, stdout, int64
    use util_mod, only: operator(.isclose.), near_zero, fuseIndex, stop_all, operator(.implies.)
    use get_excit, only: make_double, exciteIlut
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: excit_gen_store_type
    use UMatCache, only: GTID, numBasisIndices
    use SystemData, only: nEl, nBasis
    use SymExcitDataMod, only: ScratchSize
    use sltcnd_mod, only: sltcnd_excit
    use bit_rep_data, only: nIfTot
    use excitation_generators, only: doubleExcitationGenerator_t
    use FciMCData, only: ProjEDet, excit_gen_store_type
    use excitation_types, only: DoubleExc_t
    use aliasSampling, only: AliasSampler_2D_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pchb_doubles_select_particles, only: ParticleSelector_t, PCHB_ParticleSelection_t, &
                                  PCHB_particle_selection_vals, PC_FullyWeightedParticles_t, &
                                  PC_FastWeightedParticles_t, UniformParticles_t
    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    better_implicit_none

    private
    public :: GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t

    type, extends(DoubleExcitationGenerator_t) :: GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t
        !! GAS PCHB excitation generator for doubles using spin orbitals, with
        !! "fast weighting"
        !! This means we choose holes via \(p(AB|IJ)\)
        private
        logical, public :: use_lookup = .false.
            !! use a lookup for the supergroup index
        logical, public :: create_lookup = .false.
            !! create and manage the supergroup index

        class(ParticleSelector_t), allocatable :: particle_selector
            !! particle selector for getting I, J
        type(AliasSampler_2D_t) :: AB_sampler
            !! sampler for holes, dubbed \(AB\) (fused index)
            !! It yields \( p(AB | IJ, i_{\text{sg}}) \)
            !! where IJ is a fused index I < J and `i_sg` is the supergroup.

        type(SuperGroupIndexer_t), pointer :: indexer => null()
        class(GASSpec_t), allocatable :: GAS_spec
        integer, allocatable :: tgtOrbs(:, :)
    contains
        private
        procedure, public :: init => GAS_doubles_PCHB_UHF_init
        procedure, public :: finalize => GAS_doubles_PCHB_UHF_finalize
        procedure, public :: gen_exc => GAS_doubles_PCHB_uhf_gen_exc
        procedure, public :: get_pgen => GAS_doubles_PCHB_uhf_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_PCHB_uhf_gen_all_excits

        procedure :: compute_samplers => GAS_doubles_PCHB_uhf_compute_samplers
        ! @jph do we need a get_unoccupied?
    end type GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t

contains

! @jph TODO implementation for doubles

    subroutine GAS_doubles_PCHB_UHF_init(this, GAS_spec, use_lookup, &
                            create_lookup, PCHB_particle_selection)
        !! initalises the spinorb-resolved PCHB doubles excitation generator
        !!
        !! more specifically, sets up a lookup table for ab -> (a,b) and
        !! sets up the alias table for picking ab given ij with prob ~ Hijab
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_UHF_init'

        integer :: ab, a, b, abMax, nBI

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'GAS PCHB (UHF) doubles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'GAS PCHB doubles is using the supergroup lookup'

        write(stdout, *) "Allocating PCHB excitation generator objects"
        ! number of *spin* orbs
        nBI = numBasisIndices(this%GAS_spec%n_spin_orbs())
        ! initialize the mapping ab -> (a, b)
        abMax = fuseIndex(nBI, nBI)
        allocate(this%tgtOrbs(2, 0:abMax), source=0)
        do a = 1, nBI
            do b = 1, a
                ab = fuseIndex(a, b)
                this%tgtOrbs(1, ab) = b
                this%tgtOrbs(2, ab) = a
            end do
        end do

        ! set up the alias table
        call this%compute_samplers(nBI, PCHB_particle_selection)

        write(stdout, *) "Finished excitation generator initialization"

    end subroutine GAS_doubles_PCHB_UHF_init

    subroutine GAS_doubles_PCHB_UHF_finalize(this)
        !! deallocates the sampler and mapper
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this

        if (allocated(this%particle_selector)) then
            call this%AB_sampler%finalize()
            call this%particle_selector%finalize()
            ! Yes, we assume, that either all or none are allocated
            deallocate(this%particle_selector, this%tgtOrbs, this%indexer, this%GAS_spec)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
        end if
    end subroutine GAS_doubles_PCHB_UHF_finalize

    subroutine GAS_doubles_PCHB_uhf_gen_exc(&
                this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                ex, tParity, pGen, hel, store, part_type)
                ! @jph
        !! given the initial determinant (both as nI and ilut), create a random
        !! doubles excitation using the Hamiltonian matrix elements as weights
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
            !! the exctitation generator
        integer, intent(in) :: nI(nel), exFlag
            !! determinant to excite from
            !! unused in this generator
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
            !! determint from which to excite, ilut format
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
            !! the excited determinant upon return
            !! excitation order (for doubles generator, always == 2)
            !! excitation matrix nI -> nJ
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
            !! excited determinant, ilut format
        real(dp), intent(out) :: pGen
            !! probability of generating the excitation
        logical, intent(out) :: tParity
            !! the parity of the excitation nI -> nJ
        HElement_t(dp), intent(out) :: hel
            !! matrix element Hijab
        type(excit_gen_store_type), intent(inout), target :: store
            !!
        integer, intent(in), optional :: part_type
            !! unused in this generator
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_uhf_gen_exc'

        integer :: i_sg ! supergroup index
        integer :: src(2) ! particles (I, J)
        integer :: elecs(2) ! particle indices in nI
        integer :: tgt(2)

        integer :: ij
        logical :: invalid

        integer :: ab
        real(dp) :: pGenHoles


        @:unused_var(exFlag, part_type)
        ic = 2
#ifdef WARNING_WORKAROUND_
        hel = h_cast(0._dp) ! macro
#endif
        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(store%idx_curr_dets, nI)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        ! pick two random electrons
        call this%particle_selector%draw(nI, ilutI, i_sg, elecs, src, pGen)
        if (src(1) == 0) then
            call invalidate()
            return
        end if

        @:ASSERT(pgen .isclose. this%particle_selector%get_pgen(nI, i_sg, src(1), src(2)))

        invalid = .false.
        ij = fuseIndex(gtId(src(1)), gtId(src(2)))

        ! get a pair of orbitals using the precomputed weights
        call this%AB_sampler%sample(ij, i_sg, ab, pGenHoles)
        @:ASSERT(ab /= 0 .implies. (pGenHoles .isclose. this%AB_sampler%get_prob(ij, i_sg, ab)))

        if (ab == 0) then
            invalid = .true.
            tgt = 0
        else
            ! split ab -> a,b
            tgt = this%tgtOrbs(:, ab)
            @:ASSERT(all(tgt /= 0) .implies. tgt(1) /= tgt(2))
            invalid = any(tgt == 0) .or. any(tgt == nI)
        end if

        ! as in the spatially-resolved case, there is a very rare case where (due
        ! to floating point error) an excitation with pgen=0 is created. Invalidate.
        if (.not. invalid .and. near_zero(pGenHoles)) then
            invalid = .true.
            ! Yes, print. Those events are signficant enough to be always noted in the output
            write(stdout, *) "WARNING: Generated excitation with probability of 0"
        end if

        if (invalid) then
            ! if 0 returned, no excitations -> nulldet
            call invalidate()
        else
            ! construct the determinant given the selected orbitals and electrons
            call make_double(nI, nJ, elecs(1), elecs(2), tgt(1), tgt(2), ex, tParity)

            ilutJ = exciteIlut(ilutI, src, tgt)

            pgen = pgen * pgenholes

            block
                integer :: classCount2(ScratchSize), classCountUnocc2(ScratchSize)
                @:ASSERT(pGen .isclose. this%get_pgen(nI, ilutI, ex, ic, classCount2, classCountUnocc2))
            end block
        end if

        ! @jph review above

    contains

        subroutine invalidate()
            nJ = 0
            ilutJ = 0_n_int
            ex(1, 1 : 2) = src
            ex(2, 1 : 2) = tgt
        end subroutine invalidate

    end subroutine GAS_doubles_PCHB_uhf_gen_exc

    real(dp) function GAS_doubles_PCHB_uhf_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        !! calculates the probability of drawing a given double excitation
        !! parametrised by the excitation matrix ex
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
            !! excitation matrix
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        integer :: i_sg, IJ
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_uhf_get_pgen'

        ! real(dp) :: p_first(2), p_second(2), pgen_particle
        ! integer :: IJ, i_sg
        ! integer :: unoccupied(nBasis - nEl)

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        ! double excitation, so ic==2
        @:ASSERT(ic == 2)

        IJ = fuseIndex(ex(1, 1), ex(1, 2))
        i_sg = this%indexer%idx_nI(nI)



        ! pgen = this%particle_selector%get_pgen(nI, i_sg, ex(1, 1), ex(1, 2))

        ! @jph stub
        ! ! block
        ! !     integer(n_int) :: ilut_unoccupied(0 : nIfD)
        ! !     call this%get_unoccupied(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)
        ! ! end block

        ! associate(A => ex(2, 1), B => ex(2, 2)) ! hole indices
        !     p_first(1) = this%A_sampler%getProb(IJ, i_sg, unoccupied, renorm_first, A)
        !     p_first(2) = this%A_sampler%getProb(IJ, i_sg, unoccupied, renorm_first, B)

        !     renorm_second(1) = 1._dp - sum(this%B_sampler%get_prob(A, IJ, i_sg, nI))
        !     if (do_direct_calculation(renorm_second(1))) then
        !         renorm_second(1) = sum(this%B_sampler%get_prob(A, IJ, i_sg, unoccupied))
        !     end if
        !     p_second(1) = this%B_sampler%getProb(&
        !         A, IJ, i_sg, unoccupied, renorm_second(1), B)

        !     renorm_second(2) = 1._dp - sum(this%B_sampler%get_prob(B, IJ, i_sg, nI))
        !     if (do_direct_calculation(renorm_second(2))) then
        !         renorm_second(2) = sum(this%B_sampler%get_prob(B, IJ, i_sg, unoccupied))
        !     end if
        !     p_second(2) = this%B_sampler%getProb(&
        !         B, IJ, i_sg, unoccupied, renorm_second(2), A)
        !     p_second(2) = this%B_sampler%
        ! end associate

        ! pgen = pgen_particle * sum(p_first * p_second)

        ! @jph review above (first stab)

    end function GAS_doubles_PCHB_uhf_get_pgen

    subroutine GAS_doubles_PCHB_uhf_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=2)

    end subroutine GAS_doubles_PCHB_uhf_gen_all_excits

    subroutine GAS_doubles_PCHB_uhf_compute_samplers(this, nBI, PCHB_particle_selection)
        !! computes and stores values for the alias (spin-independent) sampling table
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nBI
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        integer :: i, j, ij, ijMax
        integer :: a, b, ab, abMax
        integer :: ex(2, 2)
        integer(int64) :: memCost
            !! n_supergroup * num_fused_indices * bytes_per_sampler
        real(dp), allocatable :: w(:), IJ_weights(:, :, :)
            ! @jph don't think I understand what exactly these weights are
        integer, allocatable :: supergroups(:, :)
        integer :: i_sg
        character(*), parameter :: this_routine = "GAS_doubles_PCHB_uhf_compute_samplers"
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! number of possible spin orbital pairs
        ijMax = fuseIndex(nBI, nBI)
        abMax = ijMax

        memCost = size(supergroups, 2, kind=int64) &
                  * int(ijMax, int64) &
                  * int(abMax, int64) * 8_int64

        write(stdout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(stdout, *) "The number of supergroups is", size(supergroups, 2)
        write(stdout, *) "Generating samplers for PCHB excitation generator"
        write(stdout, *) "Depending on the number of supergroups this can take up to 10min."

        ! @jph
        call this%AB_sampler%shared_alloc([ijMax, size(supergroups, 2)], abMax, 'PCHB_UHF')
        allocate(w(abMax))
        allocate(IJ_weights(nBI, nBI, size(supergroups, 2)), source=0._dp)

        do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(stdout, *) 'Still generating the samplers'
            do i = 1, nBI
                ex(1, 1) = i ! already a spin orbital
                ! a,b are ordered, so we can assume j <= i
                do j = 1, i
                    ex(1, 2) = j
                    w(:) = 0._dp
                    ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
                    ! weights to prepare the sampler
                    do a = 1, nBI
                        ex(2, 2) = a
                        do b = 1, a
                            ex(2, 1) = b
                            ab = fuseIndex(a, b)
                            w(ab) = abs(sltcnd_excit(projEDet(:, 1), DoubleExc_t(ex), .false.))
                        end do ! b
                    end do ! a
                    ij = fuseIndex(i, j)
                    ! @jph
                    call this%AB_sampler%setup_entry(ij, i_sg, w)

                    associate(I => ex(1, 1), J => ex(1, 2))
                        IJ_weights(I, J, i_sg) = IJ_weights(I, J, i_sg) + sum(w)
                        IJ_weights(J, I, i_sg) = IJ_weights(J, I, i_sg) + sum(w)
                    end associate
                end do ! j
            end do ! i
        end do ! i_sg


        if (PCHB_particle_selection == PCHB_particle_selection_vals%FULLY_WEIGHTED) then
            allocate(PC_FullyWeightedParticles_t :: this%particle_selector)
            select type(particle_selector => this%particle_selector)
            type is(PC_FullyWeightedParticles_t)
                call particle_selector%init(this%GAS_spec, IJ_weights, this%use_lookup, .false.)
            end select
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%FAST_WEIGHTED) then
            allocate(PC_FastWeightedParticles_t :: this%particle_selector)
            select type(particle_selector => this%particle_selector)
            type is(PC_FastWeightedParticles_t)
                call particle_selector%init(this%GAS_spec, IJ_weights, this%use_lookup, .false.)
            end select
        else if (PCHB_particle_selection == PCHB_particle_selection_vals%UNIFORM) then
            allocate(UniformParticles_t :: this%particle_selector)
        else
            call stop_all(this_routine, 'not yet implemented')
        end if

        ! @jph review above and compare with paper -- not sure I understand
    end subroutine GAS_doubles_PCHB_uhf_compute_samplers


end module gasci_pchb_doubles_spinorb_fastweighted
