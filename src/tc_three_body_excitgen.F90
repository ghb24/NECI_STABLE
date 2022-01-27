#include "macros.h"
module tc_three_body_excitgen
    use constants
    use SystemData, only: nel, nOccAlpha, nOccBeta, nBasis, G1, t_ueg_3_body, &
                          tContact, t_exclude_3_body_excits, symmetry, tNoSymGenRandExcits
    use lattice_mod, only: sort_unique
    use bit_rep_data, only: NIfTot
    use k_space_hubbard, only: make_triple
    use tc_three_body_data
    use FciMCData, only: excit_gen_store_type, pDoubles, pSingles
    use dSFMT_interface, only: genrand_real2_dSFMT
    use lattice_models_utils, only: make_ilutJ
    use util_mod, only: choose, intswap
    use excit_gens_int_weighted, only: pick_biased_elecs, pick_oppspin_elecs
    use GenRandSymExcitNUMod, only: calc_pgen_symrandexcit2, ScratchSize, &
                                    createSingleExcit, createDoubExcit, construct_class_counts, &
                                    gen_rand_excit
    use procedure_pointers, only: generate_two_body_excitation
    use sym_general_mod, only: ClassCountInd
    use sym_mod, only: symprod, symconj

    implicit none
    ! Factors accounting for permutation of electrons
    real(dp), parameter :: same_spin_perm = 6.0_dp
    real(dp), parameter :: opp_spin_perm = 2.0_dp

contains

    subroutine gen_excit_mol_tc(nI, ilut, nJ, ilutJ, exFlag, ic, ExcitMat, &
                                tParity, pGen, HelGen, store, part_type)
        use ueg_excit_gens, only: gen_ueg_excit
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2, maxExcit)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_mol_tc'

        real(dp) :: r

        r = genrand_real2_dSFMT()
        ! select if a triple shall be generate
        if (r < pTriples) then
            call generate_triple_excit(nI, ilut, nJ, ilutJ, ExcitMat, tParity, pGen, &
                                       HelGen, store)
            pGen = pGen * pTriples
            IC = 3
        else
            call generate_two_body_excitation(nI, ilut, nJ, ilutJ, exFlag, ic, ExcitMat, &
                                              tParity, pGen, HelGen, store, part_type)
            pGen = pGen * (1.0 - pTriples)
        end if
    end subroutine gen_excit_mol_tc

!------------------------------------------------------------------------------------------!

    function calc_pgen_triple(nI, ex) result(pgen)
        ! get the probability to get excitation `ex` from a determinant `nI`
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: ex(2, 3)
        real(dp) :: pgen
        integer :: ms, i, tgt_spin
        character(*), parameter :: t_r = "calc_pgen_triple"

        ! get the spin
        ms = 0
        ! sum up the spin of the single orbitals
        do i = 1, 3
            ms = ms + G1(ex(2, i))%ms
        end do

        ! start with pTriples
        pgen = pTriples
        ! then, add the spin bias and electron picking probs
        if (ms == -3) then
            pgen = pgen * p0A * pgen3B
        else if (ms == -1) then
            pgen = pgen * p2B * pgen2B
        else if (ms == 1) then
            pgen = pgen * (1 - p2B - p0A - p0B) * pgen1B
        else if (ms == 3) then
            pgen = pgen * p0B * pgen0B
        else
            call stop_all(t_r, "Invalid spin")
        end if

        ! Add the probability of picking these three target orbitals
        if (tNoSymGenRandExcits) then
            call calc_pgen_triple_target_nosym(ms, pgen)
        else
            call calc_pgen_triple_target_sym(nI, ex, ms, pgen)
        end if
    end function calc_pgen_triple

    !------------------------------------------------------------------------------------------!

    !> Calculates the probability of picking three orbitals with total spin ms without symmetry
    !> @param[in] ms  total spin of the picked orbitals
    !> @param[inout] pgen  on call, the probability of picking the electrons, on return, the total
    !! probability
    pure subroutine calc_pgen_triple_target_nosym(ms, pgen)
        integer, intent(in) :: ms
        real(dp), intent(inout) :: pgen

        if (ms == -3) then
            pgen = pgen * same_spin_perm &
                   / (nUnoccBeta * (nUnoccBeta - 1) * (nUnoccBeta - 2))
        else if (ms == -1) then
            pgen = pgen * opp_spin_perm &
                   / (nUnoccBeta * (nUnoccBeta - 1) * nUnoccAlpha)
        else if (ms == 1) then
            pgen = pgen * opp_spin_perm &
                   / (nUnoccBeta * nUnoccAlpha * (nUnoccAlpha - 1))
        else
            pgen = pgen * same_spin_perm &
                   / (nUnoccAlpha * (nUnoccAlpha - 1) * (nUnoccAlpha - 2))
        end if
    end subroutine calc_pgen_triple_target_nosym

    !------------------------------------------------------------------------------------------!

    !> Calculates the probability of picking three orbitals with total spin ms with symmetry
    !> @param[in] nI  determinant the excitation was made from
    !> @param[in] ex  the excitation matrix (2x3 array)
    !> @param[in] ms  total spin of the picked orbitals
    !> @param[inout] pgen  on call, the probability of picking the electrons, on return, the total
    !! probability
    subroutine calc_pgen_triple_target_sym(nI, ex, ms, pgen)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: ex(2, 3)
        integer, intent(in) :: ms
        real(dp), intent(inout) :: pgen

        ! Temporary: Store the probability of picking the canonical order
        real(dp) :: pgen_pick
        integer :: pool_sizes(3), tgt_spin, tmp_ex(2, 3)
        integer :: cc_unocc(ScratchSize), cc_occ(ScratchSize), cc_ind

        ! manually mimick pass-by-value
        tmp_ex = ex
        ! get the number of available orbs per symmetry sector
        call construct_class_counts(nI, cc_occ, cc_unocc)

        ! Now, get the prob for picking these three target orbitals

        ! The first two are chosen uniformly, and have the majority spin
        if (ms > 0) then
            pool_sizes(1:2) = nUnoccAlpha
            ! tgt_spin is the spin of the third, assign it to majority spin
            tgt_spin = 1
        else
            pool_sizes(1:2) = nUnoccBeta
            tgt_spin = 2
        end if

        if (abs(ms) /= 3) then
            ! Now, we have to be careful: The pick-algorithm has a convention on the order
            ! Alpha/Beta are picked, but this order is not preserved => Sort

            ! The convention is: the first two orbitals have the same spin
            if (G1(tmp_ex(2, 1))%MS /= G1(tmp_ex(2, 2))%MS) call intswap(tmp_ex(2, 2), tmp_ex(2, 3))
            ! Also, the third electron has minority spin now, so swap tgt_spin
            ! (map 1 -> 2 and 2 -> 1
            tgt_spin = 3 - tgt_spin
        end if
        ! Get the index of the symmetry class of the third (symmetry-restricted) orb
        cc_ind = ClassCountInd(tgt_spin, G1(tmp_ex(2, 3))%Sym%S, 0)
        ! And the from that number of available orbs
        pool_sizes(3) = cc_unocc(cc_ind)

        ! The pool for the second is one smaller, because the first one is not available anymore
        pool_sizes(2) = pool_sizes(2) - 1
        pgen_pick = pgen / product(pool_sizes)

        ! now, account for permutations
        call add_permutations_to_pgen(pgen, pgen_pick, pool_sizes, ms, tmp_ex(2, :), cc_unocc)

    end subroutine calc_pgen_triple_target_sym

!------------------------------------------------------------------------------------------!

    subroutine setup_mol_tc_excitgen()
        ! initialize the biases and auxiliary variables for the molecular
        ! transcorrelated 3-body excitation generator

        call init_mol_tc_biases()
        call precompute_pgen()
    end subroutine setup_mol_tc_excitgen

!------------------------------------------------------------------------------------------!

    subroutine precompute_pgen()

        ! set the number of unoccupied alpha/beta
        nUnoccAlpha = nBasis / 2 - nOccAlpha
        nUnoccBeta = nBasis / 2 - nOccBeta

        ! the number of valid triple excitations is just given by the binomial coefficients
        pgen3B = nOccBeta * (nOccBeta - 1) * (nOccBeta - 2)
        pgen3B = scaleInvert(same_spin_perm, pgen3B)

        pgen2B = nOccBeta * (nOccBeta - 1) * nOccAlpha
        pgen2B = scaleInvert(opp_spin_perm, pgen2B)

        pgen1B = nOccBeta * nOccAlpha * (nOccAlpha - 1)
        pgen1B = scaleInvert(opp_spin_perm, pgen1B)

        pgen0B = nOccAlpha * (nOccAlpha - 1) * (nOccAlpha - 2)
        pgen0B = scaleInvert(same_spin_perm, pgen0B)

    contains
        pure function scaleInvert(scl, p) result(sp)
            real(dp), intent(in) :: scl, p
            real(dp) :: sp

            if (p > eps) then
                sp = scl / p
            else
                sp = 0.0_dp
            end if
        end function scaleInvert
    end subroutine precompute_pgen

!------------------------------------------------------------------------------------------!

    subroutine init_mol_tc_biases()
        use SystemData, only: tSmallBasisForThreeBody
        ! reference determinant for initializing the biases
        real(dp) :: normalization
        ! if we read in a value, use that one
        if (abs(pTriples) < eps) then
            pTriples = 0.1
        end if
        ! for 2 electrons, there are obviously no
        ! triple excitations
        if (nel == 2 .or. .not. tSmallBasisForThreeBody) t_exclude_3_body_excits = .true.
        ! For contact interaction we also exclude the too small basis-sets

        if (t_exclude_3_body_excits) then
            pTriples = 0.0_dp
            return
        end if
        write(stdout, *) "pTriples set to ", pTriples
        ! pSingles and pDoubles add to 1, and pTriples is an additional bias not to do
        ! a two-body excitation
        if (tContact) then
!       We do not have those kind of excitations for triples, where all the
!       fermions are up or down.
            p0A = 0.0_dp
            p0B = 0.0_dp
!       We determine the rate uniformly between all the possible exciations
            p2B = choose(nOccBeta, 2) * nOccAlpha
            normalization = p2B + choose(nOccAlpha, 2) * nOccBeta
            p2B = p2B / normalization
        else
            ! scale the probabilities with the number of possible picks
            normalization = choose(nel, 3)
            p0A = choose(nOccBeta, 3) / normalization
            p0B = choose(noccAlpha, 3) / normalization
            p2B = choose(nOccBeta, 2) * nOccAlpha / normalization
        end if
        p1B = 1.0_dp - p0A - p0B - p2B
    end subroutine init_mol_tc_biases

!------------------------------------------------------------------------------------------!

    subroutine generate_triple_excit(nI, ilutI, nJ, ilutJ, ExcitMat, tParity, pGen, &
                                     HelGen, store)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ExcitMat(2, maxExcit)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)

        integer :: sym_prod, src(3), tgt(3), elecs(3)
        integer :: ms
        unused_var(store)
        HElGen = 0.0_dp

        ! first, pick three electrons at random
        call pick_three_elecs(nI, elecs, src, sym_prod, pgen, ms)

        ! if three electrons can be picked
        if (src(3) /= 0) then
            ! get three unoccupied orbitals with the same ms

            if (t_ueg_3_body) then
                call pick_three_orbs_ueg(nI, src, tgt, pgen, ms)
            else if (tNoSymGenRandExcits) then
                ! picking three orbs regardless of symmetry is easier
                call pick_three_orbs_nosym(nI, tgt, pgen, ms)
            else
                call pick_three_orbs_sym(nI, src, tgt, pgen, ms)
            end if

            if (tgt(3) == 0) then
                ! the excitation is invalid
                nJ = 0
                ilutJ = 0_n_int
                ExcitMat = 0
                tParity = .false.
                return
            end if

            ! and create a triple excitation
            call make_triple(nI, nJ, elecs, tgt, ExcitMat, tParity)

            ilutJ = make_ilutJ(ilutI, ExcitMat, 3)
        else
            ! else, the excitation is invalid
            nJ = 0
            ilutJ = 0_n_int
            ExcitMat = 0
            tParity = .false.
        end if

    end subroutine generate_triple_excit

!------------------------------------------------------------------------------------------!

    subroutine pick_three_elecs(nI, elecs, src, sym_prod, pgen, ms)
        ! picks three random electrons from nI, biased towards 0, 1, 2 or 3 beta electrons
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(3), src(3), sym_prod, ms
        real(dp), intent(out) :: pgen
        integer :: ispn, sum_ml
        integer :: nOcc
        real(dp) :: r
        logical :: pickAlpha

        if (tContact) then
            call pick_oppspin_elecs(nI, elecs(1:2), src(1:2), sym_prod, ispn, sum_ml, pgen)
        else
            call pick_biased_elecs(nI, elecs(1:2), src(1:2), sym_prod, ispn, sum_ml, &
                                   pgen, p0B + p0A, p0B)
        end if

        if (ispn == 3 .or. ispn == 1) then
            ! all elecs have the same spin
            pickAlpha = ispn == 3
            if (pickAlpha) then
                nOcc = nOccAlpha
                ms = 3
            else
                nOcc = nOccBeta
                ms = -3
            end if
            call get_missing_elec(nI, elecs, nOcc, 2, pickAlpha, pgen)
            pgen = pgen * 3.0_dp
        else
            ! first picked one alpha and the beta
            r = genrand_real2_dSFMT()
            if (r < p2B) then
                ! we have 2 beta elecs + 1 alpha
                ! then pick the second beta
                call get_missing_elec(nI, elecs, nOccBeta, 1, .false., pgen)
                ms = -1
                pgen = pgen * p2B / (1.0_dp - (p0B + p0A))
            else
                ! we have 2 alpha elecs + 1 beta
                ! then pick the second alpha
                call get_missing_elec(nI, elecs, nOccAlpha, 1, .true., pgen)
                ms = 1
                pgen = pgen * p1B / (1.0_dp - (p0B + p0A))
            end if
            pgen = pgen * 2.0_dp
        end if

        ! sort the generated electrons
        elecs = sort_unique(elecs)
        ! check for invalid excitation
        ! after sorting, invalid electrons are definitly at position 1
        if (elecs(1) == 0) then
            src = 0
        else
            src = nI(elecs)
        end if

    end subroutine pick_three_elecs

!------------------------------------------------------------------------------------------!

    subroutine get_missing_elec(nI, elecs, nOcc, nPicked, tAlpha, pgen)
        ! after picking two electrons using the symrandexcit routines,
        ! get a third one with the missing spin
        integer, intent(in) :: nI(nel)
        integer, intent(inout) :: elecs(3) ! picked electrons, elecs(1:2) has to be assigned on entry
        integer, intent(in) :: nOcc    ! number of occupied orbs of the type of the missing one
        integer, intent(in) :: nPicked ! number of already picked elecs
        logical, intent(in) :: tAlpha  ! if the missing electron is alpha
        real(dp), intent(inout) :: pgen ! probability of generation
        real(dp) :: r
        integer :: index, i, count

        elecs(3) = 0
        if (nOcc <= nPicked) then
            ! there are not enought elecs available - invalid excitation
            return
        end if

        r = genrand_real2_dSFMT()
        ! pick a random index of an electron
        index = int((nOcc - nPicked) * r) + 1

        ! get the corresponding electron
        count = 0
        do i = 1, nel
            if (tAlpha .neqv. is_beta(nI(i))) then
                count = count + 1
                if (count == index) elecs(3) = i
            end if
        end do
        ! the picked electrons are not counted, so skip their indices

        ! if we need, skip the first index
        call skipElec(1)
        ! if we need, skip the second index
        call skipElec(2)
        ! note that elecs(2) > elecs(1), so we cannot end up on elecs(1) at the end of
        ! skipping

        ! uniformly chosen
        pgen = pgen / (nOcc - nPicked)

    contains

        subroutine skipElec(ind)
            integer, intent(in) :: ind

            ! if we need to skip an index
            if ((tAlpha .neqv. is_beta(nI(elecs(ind))))) then
                ! if we are above the index, we need to add 1 more, because we did not
                ! take elecs(ind) into account when picking elecs(3)
                if (elecs(3) >= elecs(ind)) then
                    ! jump to the next electron with the right spin
                    elecs(3) = elecs(3) + 1
                    do while (.not. (tAlpha .neqv. is_beta(nI(elecs(3)))))
                        elecs(3) = elecs(3) + 1
                    end do
                end if
            end if
        end subroutine skipElec
    end subroutine get_missing_elec

    !------------------------------------------------------------------------------------------!

    !> picks three random unoccupied orbitals, given the occupied orbitals, ignoring symmetry
    !! This is a more efficient version of pick_three_orbs_sym for the case that point-group symmetry
    !! is not used
    subroutine pick_three_orbs_nosym(nI, tgt, pgen, ms)
        integer, intent(in) :: nI(nel), ms
        integer, intent(out) :: tgt(3)
        real(dp), intent(inout) :: pgen

        integer :: i, msCur, msOrb

        msCur = ms
        do i = 0, 2
            ! get the ms of this orb
            ! we take ms = 1 until the total leftover ms is negative
            if (msCur > 0) then
                msOrb = 1
            else
                ! then we start taking ms = -1
                msOrb = -1
            end if
            call get_rand_orb(nI, tgt, msOrb, i, pgen)
            ! the remaining ms
            msCur = msCur - msOrb
        end do

        ! adjust the probability by taking permutations into account
        pgen = pgen * 2 * abs(ms)

    end subroutine pick_three_orbs_nosym

    !------------------------------------------------------------------------------------------!

    subroutine pick_three_orbs_sym(nI, src, tgt, pgen, ms)
        ! picks three random unoccupied orbitals, given the occupied orbitals
        integer, intent(in) :: nI(nel), ms, src(3)
        integer, intent(out) :: tgt(3)
        real(dp), intent(inout) :: pgen

        integer :: i, msCur, msOrb
        integer :: cc_occ(ScratchSize), cc_unocc(ScratchSize), tgt_sym
        ! pool to draw each orbital from
        integer, allocatable :: pool(:)
        ! size of the pools
        integer :: pool_sizes(3)
        real(dp) :: pgen_pick

        msCur = ms
        pgen_pick = pgen
        ! First, pick two arbitrary orbitals
        do i = 0, 1
            ! get the ms of this orb
            ! we take ms = 1 until the total leftover ms is negative
            ! This ensures the first two electrons have the same spin
            if (msCur > 0) then
                msOrb = 1
            else
                ! then we start taking ms = -1
                msOrb = -1
            end if
            call create_full_pool(nI, tgt, msOrb, pool, i)
            call get_orb_from_pool(tgt, pool, i, pgen_pick)
            ! the remaining ms
            msCur = msCur - msOrb
            ! keep the size of the pool in memory
            pool_sizes(i + 1) = size(pool)
        end do

        ! Get the number of orbitals per symmetry/occupation
        call construct_class_counts(nI, cc_occ, cc_unocc)
        ! Now, pick the last one according to symmetry
        tgt_sym = get_tgt_sym(tgt, src)
        call create_sym_pool(nI, tgt, msCur, pool, i, tgt_sym, cc_unocc)
        call get_orb_from_pool(tgt, pool, i, pgen_pick)
        pool_sizes(i + 1) = size(pool)

        ! There might be no symmetry-allowed picks for the last orbital - in that case
        ! abort the excitation
        if (any(tgt == 0)) then
            tgt = 0
            return
        end if

        call add_permutations_to_pgen(pgen, pgen_pick, pool_sizes, ms, tgt, cc_unocc)

        ! sort the target orbitals for further usage
        tgt = sort_unique(tgt)
    end subroutine pick_three_orbs_sym

    !------------------------------------------------------------------------------------------!

    subroutine add_permutations_to_pgen(pgen, pgen_pick, pool_sizes, ms, tgt, cc_unocc)
        real(dp), intent(inout) :: pgen
        real(dp), intent(inout) :: pgen_pick
        integer, intent(in) :: pool_sizes(3), ms, tgt(3), cc_unocc(ScratchSize)

        integer :: swap_pool_size, spin_ind, cc_ind, i
        integer :: irreps(3)

        do i = 1, 3
            irreps(i) = int(G1(tgt(i))%Sym%S)
        end do
        ! adjust the probability by taking permutations into account
        ! these depend on the spin
        if (abs(ms) == 1) then
            ! We can only swap the first two electrons => factor of two
            pgen_pick = pgen_pick * 2
        else
            ! The first two could have been picked in any order
            pgen_pick = pgen_pick * 2
            ! get the spin index. Convention: 1 == alpha, 2 == beta
            if (ms > 0) then
                spin_ind = 1
            else
                spin_ind = 2
            end if

            ! Either the first or the second might have been chosen last as well
            do i = 1, 2
                ! The first might have been the last orb
                cc_ind = ClassCountInd(spin_ind, irreps(i), 0)
                swap_pool_size = cc_unocc(cc_ind)
                ! Check if the first irrep appears multiple times (shrinks pool)
                if (irreps(i) == irreps(mod(i, 3) + 1)) swap_pool_size = swap_pool_size - 1
                if (irreps(i) == irreps(mod(i + 1, 3) + 1)) swap_pool_size = swap_pool_size - 1
                ! The first two could have been picked in any order (same pool_size as above)
                pgen_pick = pgen_pick + 2.0_dp * pgen / (pool_sizes(1)) * 1.0_dp / (pool_sizes(2)) &
                            * 1.0_dp / (swap_pool_size)
            end do
        end if
        pgen = pgen_pick

    end subroutine add_permutations_to_pgen

    !------------------------------------------------------------------------------------------!

    subroutine pick_three_orbs_ueg(nI, src, tgt, pgen, ms)
        use SymExcitDataMod, only: kPointToBasisFn
        use SystemData, only: nBasisMax, tOrbECutoff, OrbECutoff, nmaxx, nmaxy, nmaxz
        use sym_mod, only: mompbcsym

        ! picks three random unoccupied orbitals, given the occupied orbitals
        integer, intent(in) :: nI(nel), ms, src(3)
        integer, intent(out) :: tgt(3)
        real(dp), intent(inout) :: pgen

        integer :: i, msCur, msOrb, k1(3), k2(3), k3(3), p1(3), p2(3), p3(3), j
        integer :: k, l
        real(dp) :: testE
        logical :: i_occup, is_allowed

        msCur = ms
        do i = 0, 1
            ! get the ms of this orb
            ! we take ms = 1 until the total leftover ms is negative
            if (msCur > 0) then
                msOrb = 1
            else
                ! then we start taking ms = -1
                msOrb = -1
            end if
            call get_rand_orb(nI, tgt, msOrb, i, pgen)
            ! the remaining ms
            msCur = msCur - msOrb
        end do

        if (msCur > 0) then
            msOrb = 1
        else
            ! then we start taking ms = -1
            msOrb = -1
        end if

        do i = 1, 3
            k1(i) = G1(src(1))%k(i)
            k2(i) = G1(src(2))%k(i)
            k3(i) = G1(src(3))%k(i)
        end do

        do i = 1, 3
            p1(i) = G1(tgt(1))%k(i)
            p2(i) = G1(tgt(2))%k(i)
        end do

        p3 = k1 + k2 + k3 - p1 - p2

        ! Is p3 allowed by the size of the space?
        testE = sum(p3**2)
        if (abs(p3(1)) <= nmaxx .and. abs(p3(2)) <= nmaxy .and. &
            abs(p3(3)) <= nmaxz .and. &
            (.not. (tOrbECutoff .and. (testE > (OrbECutoff + 1.d-12))))) then
            is_allowed = .true.
        else
            is_allowed = .false.
        end if

        if (.not. is_allowed) then
            tgt(3) = 0
            pgen = 0.0_dp
            return
        end if

        i = kPointToBasisFn(p3(1), p3(2), p3(3), (msOrb + 1) / 2 + 1)

        if (i > nBasis .or. i < 1) then
            print *, 'bug kPointToBasisFn', p3, msOrb, i
            do i = -1, 1
            do j = -1, 1
            do k = -1, 1
            do l = -1, 1
                print *, i, j, k, l, kPointToBasisFn(i, j, k, (l + 1) / 2 + 1)
            end do
            end do
            end do
            end do
            call stop_all('pick_three_orbs_ueg', "kPointToBasisFn")
        end if

        ! i occupied? check:
        i_occup = .false.
        if (i == tgt(1) .or. i == tgt(2)) then
            i_occup = .true.
        else
            do j = 1, nel
                if (i == nI(j)) i_occup = .true.
            end do
        end if

        if (i_occup) then
            tgt(3) = 0
            pgen = 0.0_dp
            return
        else
            tgt(3) = i
            tgt = sort_unique(tgt)
        end if

        ! adjust the probability by taking permutations into account
        pgen = pgen * 2 * abs(ms)

    end subroutine pick_three_orbs_ueg

!------------------------------------------------------------------------------------------!

    subroutine get_rand_orb(nI, tgt, ms, nPicked, pgen)
        integer, intent(inout) :: tgt(3)
        integer, intent(in) :: ms, nPicked, nI(nel)
        real(dp), intent(inout) :: pgen

        integer :: pool ! available orbitals
        integer :: iOrb, i
        real(dp) :: r

        ! get the number of possible orbitals
        if (ms > 0) then
            pool = nUnoccAlpha
        else
            pool = nUnoccBeta
        end if
        ! we need to see how many same spin orbs have been picked so far
        do i = 1, nPicked
            if ((ms > 0) .neqv. is_beta(tgt(i))) pool = pool - 1
        end do

        ! pick a random index
        r = genrand_real2_dSFMT()
        iOrb = spinOrb(int(r * pool) + 1)

        ! check if the orb is already targeted
        call skipPicked()
        ! assign the orbital
        tgt(nPicked + 1) = iOrb

        ! adjust the probability
        pgen = pgen / pool

        ! we need to sort tgt (see above)
        tgt(1:(nPicked + 1)) = sort_unique(tgt(1:(nPicked + 1)))

    contains
        pure function spinOrb(orb) result(sorb)
            integer, intent(in) :: orb
            integer :: sorb

            sorb = 2 * orb + (ms - 1) / 2
        end function spinOrb

        subroutine skipPicked()
            integer :: i
            integer :: invalidOrbs(nel + nPicked)

            invalidOrbs(1:nel) = nI
            invalidOrbs((nel + 1):(nel + nPicked)) = tgt(1:nPicked)

            invalidOrbs = sort_unique(invalidOrbs)

            do i = 1, (nPicked + nel)
                ! check if the orb is already targeted
                ! assumes tgt is sorted
                if (invalidOrbs(i) <= iOrb .and. G1(invalidOrbs(i))%ms == ms) then
                    iOrb = iOrb + 2
                end if
            end do
        end subroutine skipPicked
    end subroutine get_rand_orb

    !------------------------------------------------------------------------------------------!

    !> Randomly pick an orbital from a pre-arranged pool of possible orbitals
    !> @param[inout] tgt  array of size 3, contains the target orbitals
    !> @param[in] pool  array containing the available orbitals to pick from
    !> @param[in] nPicked  number of already picked orbitals
    !> @param[inout] pgen  probabaility of choosing this orbital
    subroutine get_orb_from_pool(tgt, pool, nPicked, pgen)
        integer, intent(inout) :: tgt(3)
        integer, intent(in) :: pool(:), nPicked
        real(dp), intent(inout) :: pgen

        integer :: iOrb, i, pool_size
        real(dp) :: r

        ! pick a random index
        pool_size = size(pool)
        if (pool_size > 0) then
            r = genrand_real2_dSFMT()
            iOrb = pool(int(r * pool_size) + 1)

            ! assign the orbital
            tgt(nPicked + 1) = iOrb

            ! adjust the probability
            pgen = pgen / pool_size
        else
            tgt(nPicked + 1) = 0
        end if
    end subroutine get_orb_from_pool

    !------------------------------------------------------------------------------------------!

    subroutine create_full_pool(nI, tgt, ms, pool, nPicked)
        integer, intent(inout) :: tgt(3)
        integer, intent(in) :: ms, nI(nel), nPicked
        integer, allocatable, intent(out) :: pool(:)

        integer :: i, pool_size, k

        ! get the number of possible orbitals
        if (ms > 0) then
            pool_size = nUnoccAlpha
        else
            pool_size = nUnoccBeta
        end if

        ! we need to see how many same spin orbs have been picked so far
        do i = 1, nPicked
            if ((ms > 0) .neqv. is_beta(tgt(i))) pool_size = pool_size - 1
        end do

        allocate(pool(pool_size))

        k = 0

        do i = 1, nBasis
            ! Check if this has the right spin
            if (ms > 0 .neqv. is_beta(i)) then
                ! Check if the orb is both unocc and
                if (all(tgt(1:nPicked) /= i) .and. all(nI /= i)) then
                    k = k + 1
                    pool(k) = i
                end if
            end if
        end do

        if (k /= pool_size) then
            write(stdout, *) "Error: wrong number of targets", k, "/=", pool_size
            call stop_all('create_full_pool', 'size mismatch')
        end if
    end subroutine create_full_pool

    !------------------------------------------------------------------------------------------!

    !> create a pool of unoccupied orbitals with a given symmetry
    !> @param[in] nI  determinant to excite from
    !> @param[in] tgt  array of size 3, contains the already picked target orbitals
    !> @param[in] ms  spin of the pool. ms>0 is alpha, ms<0 beta
    !> @param[out] pool  on return, a list of the available orbitals to excite to
    !> @param[in] nPicked  number of already picked target orbitals
    !> @param[in] tgt_sym  symmetry of the pool
    !> @param[in] cc_unocc  array containing the number of unoccupied orbitals per irrep
    subroutine create_sym_pool(nI, tgt, ms, pool, nPicked, tgt_sym, cc_unocc)
        integer, intent(inout) :: tgt(3)
        integer, intent(in) :: ms, nI(nel), nPicked
        integer, allocatable, intent(out) :: pool(:)
        integer, intent(in) :: tgt_sym
        integer, intent(in) :: cc_unocc(ScratchSize)

        integer :: i, pool_size, k, cc_ind, spin_ind

        if (ms > 0) then
            spin_ind = 1
        else
            spin_ind = 2
        end if
        ! get the number of possible orbitals
        cc_ind = ClassCountInd(spin_ind, tgt_sym, 0)
        pool_size = cc_unocc(cc_ind)

        ! we need to see how many same spin orbs have been picked so far
        do i = 1, nPicked
            if ((sign(1, ms) == G1(tgt(i))%Ms) .and. &
                G1(tgt(i))%Sym%S == tgt_sym) pool_size = pool_size - 1
        end do

        allocate(pool(pool_size))

        k = 0

        do i = 1, nBasis
            ! Check if this has the right spin
            if ((ms > 0 .neqv. is_beta(i)) .and. (G1(i)%Sym%S == tgt_sym)) then
                ! Check if the orb is both unocc and
                if (all(tgt(1:nPicked) /= i) .and. all(nI /= i)) then
                    k = k + 1
                    pool(k) = i
                end if
            end if
        end do

        if (k /= pool_size) then
            write(stdout, *) "Error: wrong number of targets", k, "/=", pool_size
            call stop_all('create_sym_pool', 'size mismatch')
        end if
    end subroutine create_sym_pool

    !------------------------------------------------------------------------------------------!

    !> Determine the symmetry of the third orbital
    !> @param[in] tgt  array of size 3, the first two entries are two orbitals to excite to
    !> @param[in] src  array of size 3, the three orbitals excited from
    !> @return sym  symmetry of the last orbital to excite to
    function get_tgt_sym(tgt, src) result(sym)
        integer, intent(in) :: tgt(3)
        integer, intent(in) :: src(3)
        integer :: sym

        type(symmetry) :: s_tmp

        ! Get the symmetry of the target orb
        s_tmp = symprod(G1(src(1))%Sym, G1(src(2))%Sym)
        s_tmp = symprod(G1(src(3))%Sym, s_tmp)
        s_tmp = symprod(symconj(G1(tgt(1))%Sym), s_tmp)
        s_tmp = symprod(symconj(G1(tgt(2))%Sym), s_tmp)
        s_tmp = symconj(s_tmp)
        sym = int(s_tmp%S)
    end function get_tgt_sym

!------------------------------------------------------------------------------------------!

end module tc_three_body_excitgen
