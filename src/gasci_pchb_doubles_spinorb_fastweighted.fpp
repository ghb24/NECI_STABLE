#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_doubles_spinorb_fastweighted
    !! spin-orbital-resolved GASCI-PCHB using the "fast weighted" scheme
    use constants, only: dp, n_int, maxExcit, stdout, int64
    use util_mod, only: operator(.isclose.), near_zero, fuseIndex, stop_all, operator(.implies.)
    use get_excit, only: make_double, exciteIlut
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: excit_gen_store_type
    use UMatCache, only: numBasisIndices
    use SystemData, only: nEl, nBasis
    use SymExcitDataMod, only: ScratchSize
    use sltcnd_mod, only: nI_invariant_sltcnd_excit
    use bit_rep_data, only: nIfTot
    use excitation_generators, only: doubleExcitationGenerator_t
    use FciMCData, only: ProjEDet, excit_gen_store_type
    use MPI_wrapper, only: root
    use excitation_types, only: Excite_2_t, canonicalize
    use aliasSampling, only: AliasSampler_2D_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pchb_doubles_select_particles, only: ParticleSelector_t, PCHB_ParticleSelection_t, &
                                  PCHB_particle_selection_vals, PC_FullyWeightedParticles_t, &
                                  PC_FastWeightedParticles_t, UniformParticles_t, allocate_and_init, &
                                  get_PCHB_weight
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
        procedure, public :: init => GAS_doubles_PCHB_spinorb_init
        procedure, public :: finalize => GAS_doubles_PCHB_spinorb_finalize
        procedure, public :: gen_exc => GAS_doubles_PCHB_spinorb_gen_exc
        procedure, public :: get_pgen => GAS_doubles_PCHB_spinorb_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_PCHB_spinorb_gen_all_excits

        procedure :: compute_samplers => GAS_doubles_PCHB_spinorb_compute_samplers
    end type GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t

contains

    subroutine GAS_doubles_PCHB_spinorb_init(this, GAS_spec, use_lookup, &
                            create_lookup, PCHB_particle_selection)
        !! initalises the spinorb-resolved PCHB doubles excitation generator
        !!
        !! more specifically, sets up a lookup table for ab -> (a,b) and
        !! sets up the alias table for picking ab given ij with prob ~ Hijab
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_spinorb_init'

        integer :: AB, A, B, abMax, nBI

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'GAS PCHB (spinorb) doubles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'GAS PCHB doubles is using the supergroup lookup'

        write(stdout, *) "Allocating PCHB excitation generator objects"
        ! number of *spin* orbs
        nBI = this%GAS_spec%n_spin_orbs()
        ! initialize the mapping ab -> (a, b)
        abMax = fuseIndex(nBI, nBI)
        allocate(this%tgtOrbs(2, 0:abMax), source=0)
        do A = 1, nBI
            do B = 1, A
                AB = fuseIndex(A, B)
                ! b comes first as this is an ordered list
                this%tgtOrbs(1, AB) = B
                this%tgtOrbs(2, AB) = A
            end do
        end do

        ! set up the alias table
        call this%compute_samplers(nBI, PCHB_particle_selection)

        write(stdout, *) "Finished excitation generator initialization"

    end subroutine GAS_doubles_PCHB_spinorb_init

    subroutine GAS_doubles_PCHB_spinorb_finalize(this)
        !! deallocates the sampler and mapper
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this

        if (allocated(this%particle_selector)) then
            call this%AB_sampler%finalize()
            call this%particle_selector%finalize()
            ! Yes, we assume, that either all or none are allocated
            deallocate(this%particle_selector, this%tgtOrbs, this%indexer, this%GAS_spec)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
        end if
    end subroutine GAS_doubles_PCHB_spinorb_finalize

    subroutine GAS_doubles_PCHB_spinorb_gen_exc(&
                this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                ex, tParity, pGen, hel, store, part_type)
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
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_spinorb_gen_exc'

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
        ij = fuseIndex(src(1), src(2))

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
            invalid = any(tgt == 0) .or. any(tgt(1) == nI) .or. any(tgt(2) == nI)
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

    contains

        subroutine invalidate()
            nJ = 0
            ilutJ = 0_n_int
            ex(1, 1 : 2) = src
            ex(2, 1 : 2) = tgt
        end subroutine invalidate

    end subroutine GAS_doubles_PCHB_spinorb_gen_exc

    real(dp) function GAS_doubles_PCHB_spinorb_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        !! calculates the probability of drawing a given double excitation
        !! parametrised by the excitation matrix ex
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
            !! excitation matrix
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        integer :: i_sg, IJ, AB
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_spinorb_get_pgen'

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 2)

        IJ = fuseIndex(ex(1, 1), ex(1, 2))
        AB = fuseIndex(ex(2, 1), ex(2, 2))
        i_sg = this%indexer%idx_nI(nI)

        pgen = this%particle_selector%get_pgen(nI, i_sg, ex(1, 1), ex(1, 2))
        pgen = pgen * this%AB_sampler%get_prob(IJ, i_sg, AB)

    end function GAS_doubles_PCHB_spinorb_get_pgen

    subroutine GAS_doubles_PCHB_spinorb_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=2)

    end subroutine GAS_doubles_PCHB_spinorb_gen_all_excits

    subroutine GAS_doubles_PCHB_spinorb_compute_samplers(this, nBI, PCHB_particle_selection)
        !! computes and stores values for the alias (spin-independent) sampling table
        class(GAS_PCHB_DoublesSpinOrbFastWeightedExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nBI
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        integer :: I, J, IJ, IJMax
        integer :: A, B, AB, ABMax
        integer :: ex(2, 2)
        integer(int64) :: memCost
            !! n_supergroup * num_fused_indices * bytes_per_sampler
        real(dp), allocatable :: w(:), IJ_weights(:, :, :)
        integer, allocatable :: supergroups(:, :)
        integer :: i_sg
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

        call this%AB_sampler%shared_alloc([ijMax, size(supergroups, 2)], abMax, 'AB_PCHB_spinorb_sampler')

        ! One could allocate only on the intra-node-root here, if memory
        ! at initialization ever becomes an issue.
        ! Look at `gasci_pchb_doubles_spin_fulllyweighted.fpp` for inspiration.
        allocate(w(abMax))
        allocate(IJ_weights(nBI, nBI, size(supergroups, 2)), source=0._dp)

        supergroup: do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(stdout, *) 'Still generating the samplers'
            first_particle: do I = 1, nBI
                ex(1, 1) = I ! already a spin orbital
                ! A,B are ordered, so we can assume J < I
                second_particle: do J = 1, I - 1
                    ex(1, 2) = J
                    w(:) = 0._dp
                    ! for each (I,J), get all matrix elements <IJ|H|AB> and use them as
                    ! weights to prepare the sampler
                    first_hole: do A = 1, nBI
                        if (any(A == [I, J])) cycle
                        ex(2, 1) = A
                        second_hole: do B = 1, nBI
                            if (A == B .or. any(B == [I, J])) cycle
                            ex(2, 2) = B
                            AB = fuseIndex(A, B)
                            associate(exc => canonicalize(Excite_2_t(ex)))
                                if (this%GAS_spec%is_allowed(exc, supergroups(:, i_sg))) then
                                    w(AB) = get_PCHB_weight(exc)
                                else
                                    w(AB) = 0._dp
                                end if
                            end associate
                        end do second_hole
                    end do first_hole
                    IJ = fuseIndex(I, J)
                    call this%AB_sampler%setup_entry(IJ, i_sg, root, w)

                    IJ_weights(I, J, i_sg) = sum(w)
                    IJ_weights(J, I, i_sg) = IJ_weights(I, J, i_sg)
                end do second_particle
            end do first_particle
        end do supergroup

        call allocate_and_init(PCHB_particle_selection, this%GAS_spec, &
            IJ_weights, root, this%use_lookup, this%particle_selector)

    end subroutine GAS_doubles_PCHB_spinorb_compute_samplers

end module gasci_pchb_doubles_spinorb_fastweighted
