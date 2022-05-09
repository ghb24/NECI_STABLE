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
!    Guther K. et al., J. Chem. Phys. 153, 034107 (2020).
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx.
! This gives #{i, j | i < j} * 3 probability distributions.
!
! The generalization to GAS spaces is available in a preprint (should be available in JCTC soon as well)
!    https://chemrxiv.org/engage/chemrxiv/article-details/61447e60b1d4a605d589af2e
! The main "ingredient" are precomputed probability distributions p(ab | ij, s_idx, i_sg) to draw a, b holes
! when i, j electrons were chosen for three distinc spin cases given by s_idx and a supergroup index i_sg
! This gives #{i, j | i < j} * 3 * n_supergroup probability distributions.
! Depending on the supergroup and GAS constraints certain excitations can be forbidden by setting p to zero.
!
! The details of calculating i_sg can be found in gasci_supergroup_index.f90

module gasci_pchb
    use constants, only: n_int, dp, int64, maxExcit, stdout, bits_n_int, int32
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw, operator(==), operator(/=), alpha, beta
    use util_mod, only: fuseIndex, getSpinIndex, near_zero, swap, operator(.div.), operator(.implies.), EnumBase_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use get_excit, only: make_double, exciteIlut
    use SymExcitDataMod, only: pDoubNew, ScratchSize
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite
    use sltcnd_mod, only: sltcnd_excit
    use procedure_pointers, only: generate_single_excit_t
    use aliasSampling, only: AliasSampler_2D_t
    use UMatCache, only: gtID, numBasisIndices
    use FciMCData, only: pSingles, excit_gen_store_type, pParallel, projEDet
    use excit_gens_int_weighted, only: pick_biased_elecs
    use shared_ragged_array, only: shared_ragged_array_int32_t
    use growing_buffers, only: buffer_int32_1D_t
    use parallel_neci, only: iProcIndex_intra
    use get_excit, only: make_single
    use growing_buffers, only: buffer_int_2D_t
    use timing_neci, only: timer, set_timer, halt_timer

    use SystemData, only: nEl, AB_elec_pairs, par_elec_pairs
    use bit_rep_data, only: NIfTot, nIfD
    use bit_reps, only: decode_bit_det
    use sort_mod, only: sort
    use DetBitOps, only: EncodeBitDet, ilut_lt, ilut_gt

    use gasci, only: GASSpec_t
    use gasci_general, only: GAS_singles_heat_bath_ExcGen_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use exc_gen_class_wrappers, only: UniformSingles_t

    use excitation_generators, only: &
            ExcitationGenerator_t, SingleExcitationGenerator_t, &
            DoubleExcitationGenerator_t, gen_exc_sd, get_pgen_sd, gen_all_excits_sd
    implicit none

    private
    public :: GAS_PCHB_ExcGenerator_t, use_supergroup_lookup, GAS_doubles_PCHB_ExcGenerator_t, &
        possible_GAS_singles, GAS_PCHB_singles_generator

    !> The precomputed GAS uniform excitation generator
    type :: PC_WeightedParticles_t
        private
        !> The shape is (number_of_spin_orbs, n_supergroup)
        type(AliasSampler_2D_t) :: samplers

        class(GASSpec_t), allocatable :: GAS_spec
        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        !> Use a lookup for the supergroup index in global_det_data.
        logical, public :: use_lookup = .false.
        !> Create **and** manage! the supergroup index lookup in global_det_data.
        logical, public :: create_lookup = .false.
    contains
        private
        procedure, public :: init
        procedure, public :: finalize

        procedure, public :: draw
        procedure, public :: get_pgen
    end type
contains

    subroutine init(this, GAS_spec)
        class(PC_weighted_ExcGenerator_t), intent(in) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        character(*), parameter :: this_routine = 'init'

        integer(int64) :: memCost
        real(dp), allocatable :: w(:)
        integer, allocatable :: supergroups(:, :)
        integer :: nBI
        integer :: i_sg

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        this%use_lookup = use_lookup

        if (this%create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'PC weighted singles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        if (this%use_lookup) write(stdout, *) 'PC weighted singles is using the supergroup lookup'

        ASSERT(this%GAS_spec%n_spin_orbs() == numBasisIndices(this%GAS_spec%n_spin_orbs()))
        nBI = this%GAS_spec%n_spin_orbs()

        supergroups = this%indexer%get_supergroups()
        !> n_supergroup * number_of_spin_orbs * (bytes_per_sampler)
        memCost = size(supergroups, 2, kind=int64) &
                    * nBI
                    * (int(abMax, int64) * 3_int64 * 8_int64)


        write(stdout, *) "Excitation generator requires", real(memCost, dp) / 2.0_dp**30, "GB of memory"
        write(stdout, *) "The number of supergroups is", size(supergroups, 2)
        write(stdout, *) "Generating samplers for weighted singles"
        write(stdout, *) "Depending on the number of supergroups this can take up to 10min."
        call this%pchb_samplers%shared_alloc([nBI, size(supergroups, 2)], nBI, 'PC_singles')
        allocate(w(nBI))

        do i_sg = 1, size(supergroups, 2)
            do i = 1, nBI
                ex(1, 1) = i
                w(:) = 0.0_dp
                do j = 1, nBI
                    ex(1, 2) = j
                    if (i == j) cycle
                    do a = 1, nBI
                        ex(2, 2) = a
                        do b = 1, nBI
                            ! ex(2, :) is in ascending order
                            ex(2, 1) = b
                            if (a /= b &
                                    .and. all(a /= ex(1, :)) &
                                    .and. all(b /= ex(1, :)) &
                                    .and. this%GAS_spec%is_allowed(DoubleExc_t(ex), supergroups(:, i_sg))) then
                                w(j) = w(j) + abs(sltcnd_excit(projEDet(:, 1), DoubleExc_t(ex), .false.))
                            end if
                        end do
                    end do
                end do
                call this%pchb_samplers%setup_entry(i, i_sg, w)
            end do
        end do

    end subroutine


    subroutine finalize(this)
        class(PC_weighted_ExcGenerator_t), intent(inout) :: this

        call this%pchb_samplers%finalize()
        deallocate(this%indexer)
        if (this%create_lookup) then
            nullify(lookup_supergroup_indexer)
        end if
    end subroutine
end module
