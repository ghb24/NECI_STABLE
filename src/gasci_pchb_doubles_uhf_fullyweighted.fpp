#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_pchb_doubles_UHF_fullyweighted
    !! precomputed heat bath implementation for GASCI using spin orbitals and fully weighting
    use constants, only: n_int, dp, int64, maxExcit, stdout, bits_n_int
    use util_mod, only: fuseIndex, getSpinIndex, near_zero, swap, &
        operator(.implies.), operator(.isclose.), isclose, swap, stop_all
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use bit_reps, only: decode_bit_det
    use bit_rep_data, only: nIfD
    use SymExcitDataMod, only: pDoubNew, ScratchSize
    use excitation_types, only: DoubleExc_t, excite
    use sltcnd_mod, only: sltcnd_excit
    use CDF_sampling_mod, only: CDF_Sampler_t
    use aliasSampling, only: AliasSampler_2D_t, AliasSampler_3D_t
    use UMatCache, only: gtID, numBasisIndices
    use FciMCData, only: excit_gen_store_type, projEDet
    use excit_gens_int_weighted, only: pick_biased_elecs
    use SystemData, only: nEl, nBasis
    use sets_mod, only: set, operator(.cap.)
    use bit_rep_data, only: NIfTot
    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use gasci_pchb_doubles_select_particles, only: &
        ParticleSelector_t, PC_FullyWeightedParticles_t, &
        PC_FastWeightedParticles_t, UniformParticles_t, &
        PCHB_ParticleSelection_t, PCHB_particle_selection_vals, &
        allocate_and_init
    use excitation_generators, only: DoubleExcitationGenerator_t
    better_implicit_none

    private
    public :: GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t

    type, extends(DoubleExcitationGenerator_t) :: GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t
        !! The GAS PCHB excitation generator for doubles using spin orbitals
        !! and doing full weighting.
        !! This means that first a hole is chosen via \( p( A | I J) |_{A \notin D_i} )\
        !! then a second hole is chosen via \( p( B | I J A) |_{B \notin D_i} )\.
        !! This means it is guaranteed that only unoccupied sites are sampled.
        private
        logical, public :: use_lookup = .false.
            !! Use a lookup for the supergroup index in global_det_data
        logical, public :: create_lookup = .false.
            !! Create **and** manage! the supergroup index lookup in global_det_data.

        class(ParticleSelector_t), allocatable :: particle_selector
            !! The particle selector for I, J
        type(AliasSampler_2D_t) :: A_sampler
            !! The sampler for the first hole `A`.
            !! It yields \( p(A | IJ, i_{\text{sg}}) \)
            !! where IJ is a fused index I < J and `i_sg` is the supergroup.
        type(AliasSampler_3D_t) :: B_sampler
            !! The sampler for the second hole `B`.
            !! It yields \( p(B | A, IJ, i_{\text{sg}}) \)
            !! where IJ is a fused index I < J, A is the first hole and `i_sg` is the supergroup.
        class(GASSpec_t), allocatable :: GAS_spec
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        logical :: hole_excess
            !! If there are more holes than particles we can
            !! speed up the calculation by normalizing via the complement.
            !! \( \sum_{A \notin D_i} p(A) = 1 - \sum_{I \in D_i} p(I) \)
        integer(n_int), private :: last_possible_occupied
            !! The last element of the ilut array has some elements
            !! which are not used, if the number of spinorbitals is not a multiple of
            !! bitsize_n_int.
            !! To correctly zero them this bitmask is 1 wherever a determinant
            !! could be occupied in the last element, and 0 otherwise.
    contains
        private
        procedure, public :: init => GAS_doubles_PCHB_init
        procedure, public :: finalize => GAS_doubles_PCHB_finalize
        procedure, public :: gen_exc => GAS_doubles_PCHB_gen_exc
        procedure, public :: get_pgen => GAS_doubles_PCHB_get_pgen
        procedure, public :: gen_all_excits => GAS_doubles_PCHB_gen_all_excits

        procedure :: compute_samplers => GAS_doubles_PCHB_compute_samplers
        procedure :: get_unoccupied
    end type GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t

    real(dp), parameter :: direct_calculation = 1e-5_dp
        !! For weights of constrained subsets that are below `direct_calculation`
        !! we calculate the weight not indirectly via the complement, but directly.
        !! Theoretically we have
        !! \Sum{ p_i}{i \in D_i} = 1 - \Sum{ p_i }{ i \notin D_i }
        !! which for small LHS is not true for floating point numbers.

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
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_init'

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

        this%hole_excess = nBasis > nEl * 2

        this%last_possible_occupied = 0_n_int
        block
            integer :: i
            do i = 0, ilut_off(nBasis)
                this%last_possible_occupied = ibset(this%last_possible_occupied, i)
            end do
        end block

        call this%compute_samplers(PCHB_particle_selection)

        write(stdout, *) "Finished excitation generator initialization"
    end subroutine GAS_doubles_PCHB_init


    subroutine GAS_doubles_PCHB_finalize(this)
        !! Finalize everything
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(inout) :: this

        if (allocated(this%particle_selector)) then
            ! Yes, we assume, that either all or none are allocated
            call this%A_sampler%finalize()
            call this%B_sampler%finalize()
            call this%particle_selector%finalize()
            deallocate(this%particle_selector, this%GAS_spec, this%indexer)
            if (this%create_lookup) nullify(lookup_supergroup_indexer)
        end if

    end subroutine GAS_doubles_PCHB_finalize


    subroutine GAS_doubles_PCHB_gen_exc(&
                    this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                    ex, tParity, pGen, hel, store, part_type)
        !>  Given the initial determinant (both as nI and ilut), create a random double
        !>  excitation using the hamiltonian matrix elements as weights
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(inout) :: this
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

        integer :: elecs(2), src(2), tgt(2), IJ
        real(dp) :: p_IJ, p_first(2), p_second(2), renorm_first, renorm_second(2)
        integer :: unoccupied(nBasis - nEl)
        integer(n_int) :: ilut_unoccupied(0 : nIfD)
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
        call this%particle_selector%draw(nI, ilutI, i_sg, elecs, src, p_IJ)
        if (src(1) == 0) then
            call invalidate()
            return
        end if

        @:ASSERT(p_IJ .isclose. this%particle_selector%get_pgen(nI, i_sg, src(1), src(2)))

        call this%get_unoccupied(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)
        IJ = fuseIndex(src(1), src(2))

        renorm_first = 1._dp - sum(this%A_sampler%get_prob(IJ, i_sg, nI))
        block
            integer :: dummy
            call this%A_sampler%constrained_sample(&
                IJ, i_sg, unoccupied, ilut_unoccupied, renorm_first, dummy, tgt(1), p_first(1))
        end block
        if (tgt(1) == 0) then
            call invalidate()
            return
        end if

        renorm_second(1) = 1._dp - sum(this%B_sampler%get_prob(tgt(1), IJ, i_sg, nI))
        if (renorm_second(1) >= 1._dp .or. renorm_second(1) < direct_calculation) then
            renorm_second(1) = sum(this%B_sampler%get_prob(tgt(1), IJ, i_sg, unoccupied))
        end if
        if (near_zero(renorm_second(1))) then
            call invalidate()
            return
        end if

        block
            integer :: dummy
            call this%B_sampler%constrained_sample(&
                tgt(1), IJ, i_sg, unoccupied, ilut_unoccupied, &
                renorm_second(1), dummy, tgt(2), p_second(1))
        end block
        if (tgt(2) == 0) then
            call invalidate()
            return
        end if

        ! We could have picked them the other way round.
        ! Account for that.
        ! The renormalization for the first hole is the same
        p_first(2) = this%A_sampler%constrained_getProb(&
            IJ, i_sg, unoccupied, renorm_first, tgt(2))
        renorm_second(2) = 1._dp - sum(this%B_sampler%get_prob(tgt(2), IJ, i_sg, nI))
        if (renorm_second(2) >= 1._dp .or. renorm_second(2) < direct_calculation) then
            renorm_second(2) = sum(this%B_sampler%get_prob(tgt(2), IJ, i_sg, unoccupied))
        end if
        p_second(2) = this%B_sampler%constrained_getProb(&
            tgt(2), IJ, i_sg, unoccupied, renorm_second(2), tgt(1))

        pGen = p_IJ * sum(p_first * p_second)

        call make_double(nI, nJ, elecs(1), elecs(2), tgt(1), tgt(2), ex, tParity)
        ilutJ = exciteIlut(ilutI, src, tgt)
        block
            integer :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
            @:ASSERT(isclose(pgen, this%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2), atol=1e-10_dp), &
                        pgen, &
                        this%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2), &
                        abs(pgen - this%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)) &
                    )
        end block

    contains

        subroutine invalidate()
            nJ = 0
            ilutJ = 0_n_int
            ex(1, 1 : 2) = src
            ex(2, 1 : 2) = tgt
        end subroutine invalidate
    end subroutine GAS_doubles_PCHB_gen_exc


    !>  @brief
    !>  Calculate the probability of drawing a given double excitation ex
    !>
    !>  @param[in] ex  2x2 excitation matrix
    !>
    !>  @return pgen  probability of generating this double with the pchb double excitgen
    function GAS_doubles_PCHB_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        real(dp) :: pgen_particle
        character(*), parameter :: this_routine = 'GAS_doubles_PCHB_get_pgen'

        integer :: IJ, i_sg
        integer :: unoccupied(nBasis - nEl)
        real(dp) :: renorm_first, p_first(2), renorm_second(2), p_second(2)

        @:unused_var(ilutI, ClassCount2, ClassCountUnocc2)
        @:ASSERT(ic == 2)

        IJ = fuseIndex(ex(1, 1), ex(1, 2))
        i_sg = this%indexer%idx_nI(nI)

        pgen_particle = this%particle_selector%get_pgen(nI, i_sg, ex(1, 1), ex(1, 2))

        block
            integer(n_int) :: ilut_unoccupied(0 : nIfD)
            call this%get_unoccupied(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)
        end block
        ! The renormalization for the first hole is the same, regardless of order.
        renorm_first = 1._dp - sum(this%A_sampler%get_prob(IJ, i_sg, nI))

        associate(A => ex(2, 1), B => ex(2, 2))
            p_first(1) = this%A_sampler%constrained_getProb(IJ, i_sg, unoccupied, renorm_first, A)
            p_first(2) = this%A_sampler%constrained_getProb(IJ, i_sg, unoccupied, renorm_first, B)

            renorm_second(1) = 1._dp - sum(this%B_sampler%get_prob(A, IJ, i_sg, nI))
            if (renorm_second(1) >= 1._dp .or. renorm_second(1) < direct_calculation) then
                ! In this rare occasion it might be that the whole distribution
                ! has zero probability, then renorm_second == 1 would be wrong.
                ! We have to calculate exactly:
                renorm_second(1) = sum(this%B_sampler%get_prob(A, IJ, i_sg, unoccupied))
            end if
            p_second(1) = this%B_sampler%constrained_getProb(&
                A, IJ, i_sg, unoccupied, renorm_second(1), B)

            renorm_second(2) = 1._dp - sum(this%B_sampler%get_prob(B, IJ, i_sg, nI))
            if (renorm_second(2) >= 1._dp .or. renorm_second(2) < direct_calculation) then
                ! In this rare occasion it might be that the whole distribution
                ! has zero probability, then renorm_second == 1 would be wrong.
                ! We have to calculate exactly:
                renorm_second(2) = sum(this%B_sampler%get_prob(B, IJ, i_sg, unoccupied))
            end if
            p_second(2) = this%B_sampler%constrained_getProb(&
                B, IJ, i_sg, unoccupied, renorm_second(2), A)
        end associate

        pgen = pgen_particle * sum(p_first * p_second)

    end function GAS_doubles_PCHB_get_pgen


    subroutine GAS_doubles_PCHB_compute_samplers(this, PCHB_particle_selection)
        !! computes and stores values for the alias sampling table
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(inout) :: this
        type(PCHB_ParticleSelection_t), intent(in) :: PCHB_particle_selection
        integer :: I, J, IJ, IJ_max, A, B ! Uppercase because they are indexing spin orbitals
        integer :: ex(2, 2)
        integer(int64) :: memCost
        real(dp), allocatable :: w_A(:), w_B(:), IJ_weights(:, :, :)
        integer, allocatable :: supergroups(:, :)
        integer :: i_sg
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! number of possible source orbital pairs
        IJ_max = fuseIndex(nBasis, nBasis)

        !> n_supergroup * number_of_fused_indices * 3 * (bytes_per_sampler)
        memCost = size(supergroups, 2, kind=int64) &
                    * (int(IJ_max, int64) * (int(nBasis, int64) * 3_int64 * 8_int64) & ! p(A | IJ)
                        + int(IJ_max, int64) * (int(nBasis, int64)**2 * 3_int64 * 8_int64) & ! P (B | IJ A)
                    )

        write(stdout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(stdout, *) "The number of supergroups is", size(supergroups, 2)
        write(stdout, *) "Generating samplers for PCHB excitation generator"
        write(stdout, *) "Depending on the number of supergroups this can take up to 10min."
        call this%A_sampler%shared_alloc([IJ_max, size(supergroups, 2)], nBasis, 'A_sampler_PCHB_UHF_FullyWeighted')
        call this%B_sampler%shared_alloc([nBasis, IJ_max, size(supergroups, 2)], nBasis, 'A_sampler_PCHB_UHF_FullyWeighted')
        allocate(w_A(nBasis), w_B(nBasis))

        allocate(IJ_weights(nBasis, nBasis, size(supergroups, 2)), source=0._dp)
        ! initialize the three samplers
        do i_sg = 1, size(supergroups, 2)
            if (mod(i_sg, 100) == 0) write(stdout, *) 'Still generating the samplers'
            first_particle: do I = 1, nBasis
                ex(1, 1) = I
                second_particle: do J = 1, I - 1
                    ex(1, 2) = J
                    IJ = fuseIndex(I, J)
                    w_A(:) = 0.0_dp
                    first_hole: do A = 1, nBasis
                        if (any(A == [I, J])) cycle
                        ex(2, 1) = A
                        w_B(:) = 0.0_dp
                        second_hole: do B = 1, nBasis
                            if (A == B .or. any(B == [I, J])) cycle
                            ex(2, 2) = B
                            if (this%GAS_spec%is_allowed(DoubleExc_t(ex), supergroups(:, i_sg))) then
                                w_B(B) = abs(sltcnd_excit(projEDet(:, 1), DoubleExc_t(ex), .false.))
                            else
                                w_B(B) = 0._dp
                            end if
                        end do second_hole
                        call this%B_sampler%setup_entry(A, IJ, i_sg, w_B(1 : nBasis))
                        w_A(A) = sum(w_B)
                    end do first_hole
                    call this%A_sampler%setup_entry(IJ, i_sg, w_A(1 : nBasis))

                    IJ_weights(I, J, i_sg) = sum(w_A)
                    IJ_weights(J, I, i_sg) = sum(w_A)
                end do second_particle
            end do first_particle
        end do

        call allocate_and_init(PCHB_particle_selection, this%GAS_spec, IJ_weights, this%use_lookup, this%particle_selector)

    end subroutine GAS_doubles_PCHB_compute_samplers


    subroutine GAS_doubles_PCHB_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), allocatable, intent(out) :: det_list(:,:)

        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list, ic=2)
    end subroutine GAS_doubles_PCHB_gen_all_excits


    pure subroutine get_unoccupied(this, ilutI, ilut_unoccupied, unoccupied)
        !! Return a bitmask and enumeration of the unoccupied spin orbitals.
        class(GAS_PCHB_Doubles_UHF_FullyWeighted_ExcGenerator_t), intent(in) :: this
        integer(n_int), intent(in) :: ilutI(0 : nIfD)
        integer(n_int), intent(out) :: ilut_unoccupied(0 : nIfD)
        integer, intent(out) :: unoccupied(nBasis - nEl)
        ilut_unoccupied = not(ilutI)
        ilut_unoccupied(nIfd) = iand(ilut_unoccupied(nIfd), this%last_possible_occupied)
        call decode_bit_det(unoccupied, ilut_unoccupied)
    end subroutine

end module gasci_pchb_doubles_UHF_fullyweighted
