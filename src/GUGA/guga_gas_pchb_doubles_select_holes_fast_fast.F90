#include "macros.h"
module guga_gas_pchb_doubles_select_holes_fast_fast

    use aliasSampling, only: AliasSampler_2D_t
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm, bits_n_int
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel, G1, &
                          nBasis, nSpatOrbs, ElecPairs, &
                          t_analyze_pchb, t_old_pchb, t_exchange_pchb
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type, ExcitationType_t, &
        DistinctDoubleType_t, distinct_doubles, DistinctDouble_t, DistinctDoubleTypeVals_t, ijab_Index_t
    use guga_bitrepops, only: convert_ilut_toGUGA, isProperCSF_ilut, CSF_Info_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: near_zero, fuseIndex, binary_search_first_ge, &
                        get_free_unit, stop_all, operator(.isclose.), operator(.div.), &
                        EnumBase_t
    use CalcData, only: t_matele_cutoff, matele_cutoff, t_truncate_spawns
    use tau_main, only: max_tau, tau_search_method, possible_tau_search_methods
    use tau_search_hist, only: frq_ratio_cutoff, max_frequency_bound, n_frequency_bins
    use fortran_strings, only: to_upper
    use sym_general_mod, only: ClassCountInd
    use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat, SpinOrbSymLabel
    use excitation_types, only: Excite_1_t, Excite_2_t
    use UMatCache, only: gtID
    use guga_excitations, only: assign_excitinfo_values_single, &
                                pick_elec_pair_uniform_guga, &
                                excitationIdentifier_double, &
                                calc_pgen_mol_guga_single, get_excit_level_from_excitInfo, &
                                checkCompatibility_single
    use guga_bitrepops, only: identify_excitation, encode_excit_info, extract_excit_info, &
                              contract_2_rdm_ind
    use guga_base_class, only: DoublesGUGABase_t
    use bit_reps, only: decode_bit_det
    use shared_array, only: shared_array_int64_t, shared_array_real_t
    use MPI_wrapper, only: iProcIndex_intra, iprocindex, root
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use procedure_pointers, only: get_umat_el

    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use sets_mod, only: is_sorted, is_set

    use guga_gas_pchb_doubles_select_holes_base, only: GAS_PCHB_Doubles_t

    use guga_gas_pchb_doubles_select_particles, only: ParticleSelector_t, &
        UniformParticles_t, PCHB_particle_selection_vals, &
        PCHB_ParticleSelection_t, particle_allocate_and_init => allocate_and_init


    better_implicit_none

    private
    public :: GAS_PCHB_FastFast_t


    type, extends(GAS_PCHB_Doubles_t) :: GAS_PCHB_FastFast_t
        private

        !> The shape is (fused_number_of_double_excitations, n_supergroup)
        type(AliasSampler_2D_t) :: alias_sampler
        integer, allocatable :: tgtOrbs(:, :)

    contains
        private
        procedure, public :: pickOrbitals
        procedure, public :: calc_pgen

        procedure, public :: calc_orbital_pgen_contr
        procedure, public :: calc_orbital_pgen_contr_start
        procedure, public :: calc_orbital_pgen_contr_end

        procedure :: setup_pchb_sampler_conditional
        procedure :: finalize_pchb_sampler
    end type

contains

    subroutine setup_pchb_sampler_conditional(this, ij_weights)
        class(GAS_PCHB_FastFast_t), intent(inout) :: this
        real(dp), allocatable, intent(out) :: ij_weights(:, :, :)
        integer :: i, j, ij, a, b, ab, i_sg
        integer(int64) :: ijMax, abMax, memCost
        integer, allocatable :: supergroups(:, :)
        real(dp), allocatable :: w(:)

        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        ! initialize the mapping ab -> (a, b)
        abMax = fuseIndex(nSpatOrbs, nSpatOrbs)
        ijMax = abMax
        allocate(this%tgtOrbs(2, 0:abMax), source=0)
        do a = 1, nSpatOrbs
            do b = 1, a
                ab = fuseIndex(a, b)
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


        if (iProcIndex_intra == root) then
            allocate(w(abMax), source=0._dp)
            allocate(ij_weights(nSpatOrbs, nSpatOrbs, size(supergroups, 2)), source=0._dp)
        else
            allocate(w(0))
            allocate(ij_weights(0, 0, 0))
        end if

        do i_sg = 1, size(supergroups, 2)
            do i = 1, nSpatOrbs
                do j = i, nSpatOrbs
                    w = 0.0_dp
                    ij = fuseIndex(i, j)
                    do a = 1, nSpatOrbs
                        do b = a, nSpatOrbs
                            ab = fuseIndex(a, b)
                            if (iProcIndex_intra == root) then
                                associate(exc => DistinctDouble_t(ijab_Index_t([i, j, a, b])))
                                if (this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) .and. symmetry_allowed(exc)) then
                                    w(ab) = this%get_weight(exc)
                                end if
                                end associate
                            end if
                        end do
                    end do
                    call this%alias_sampler%setup_entry(ij, i_sg, root, w)
                end do
            end do
        end do
    end subroutine setup_pchb_sampler_conditional

    subroutine finalize_pchb_sampler(this)
        class(GAS_PCHB_FastFast_t), intent(inout) :: this
        call this%alias_sampler%finalize()
        deallocate(this%tgtOrbs)
    end subroutine

    subroutine pickOrbitals(this, ilut, nI, csf_i, excitInfo, pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src(2), a, b, orbs(2), ij, ab, i_sg
        real(dp) :: pgen_elec, pgen_orbs

        unused_var(ilut)
        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        call this%particle_selector%draw(nI, ilut, csf_i, src, pgen_elec)

        ! use the sampler for this electron pair -> order of src electrons
        ! does not matter
        ij = fuseIndex(src(1), src(2))

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
    end subroutine pickOrbitals

    function calc_pgen(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: ij, ab
        real(dp) :: p_elec
        unused_var(nI); unused_var(ilutI); unused_var(csf_i)

        ! i, j, k and l entries have the info about the electrons and
        ! holes of the excitation: E_{ij}E_{kl} is the convention ...
        ij = fuseIndex(excitInfo%j, excitInfo%l)
        ab = fuseIndex(excitInfo%i, excitInfo%k)

        p_elec = this%particle_selector%get_pgen(nI, csf_i, excitInfo%j, excitInfo%l)
        pgen = pDoubles * p_elec * this%alias_sampler%get_prob(ij, excitInfo%i_sg_start, ab)
    end function calc_pgen

    ! I need the pgen-recalculation routines for exchange type excitations
    ! also for the PCHB excit-gen
    pure subroutine calc_orbital_pgen_contr(this, csf_i, occ_orbs, above_cpt, below_cpt)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt

        integer :: ij, i_sg

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_supergroup_idx(csf_I)
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
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = 'calc_orbital_pgen_contr_start'

        integer :: i, j, ij, ab, i_sg

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
        ij = fuseIndex(i, j)
        ab = fuseIndex(i, orb_a)

        orb_pgen = this%alias_sampler%get_prob(ij, i_sg, ab)
    end subroutine calc_orbital_pgen_contr_start

    subroutine calc_orbital_pgen_contr_end(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(GAS_PCHB_FastFast_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = "calc_orbital_pgen_contr_end"

        integer :: i, j, ij, ab, i_sg

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
        ij = fuseIndex(i, j)
        ab = fuseIndex(orb_a, j)

        orb_pgen = this%alias_sampler%get_prob(ij, i_sg, ab)
    end subroutine calc_orbital_pgen_contr_end

end module
