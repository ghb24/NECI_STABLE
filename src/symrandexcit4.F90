#include "macros.h"

module excit_gens_int_weighted

    use SystemData, only: nel, nbasis, nOccAlpha, nOccBeta, G1, nOccAlpha, &
                          nOccBeta, tExch, AA_elec_pairs, BB_elec_pairs, &
                          AB_elec_pairs, par_elec_pairs, AA_hole_pairs, &
                          par_hole_pairs, AB_hole_pairs
    use SymExcit3, only: CountExcitations3, GenExcitations3
    use SymExcitDataMod, only: SymLabelList2, SymLabelCounts2, OrbClassCount, &
                               pDoubNew, ScratchSize
    use sym_general_mod, only: ClassCountInd, ClassCountInv, class_count_ms
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles, &
                         rand_excit_opp_bias
    use dSFMT_interface, only: genrand_real2_dSFMT
    use Determinants, only: get_helement, write_det
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, ilut_lt, ilut_gt
    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: decode_bit_det
    use symdata, only: nSymLabels
    use procedure_pointers, only: get_umat_el
    use UMatCache, only: gtid
    use OneEInts, only: GetTMATEl
    use sltcnd_mod, only: sltcnd_1
    use GenRandSymExcitNUMod, only: PickElecPair, gen_rand_excit, &
                                    init_excit_gen_store, &
                                    construct_class_counts, &
                                    RandExcitSymLabelProd
    use constants
    use get_excit, only: make_double, make_single
    use sort_mod
    use util_mod
    implicit none
    save

    integer :: ncnt, nsel

contains

    subroutine gen_excit_hel_weighted (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                       ExcitMat, tParity, pGen, HelGen, store)

        ! A really laborious, slow, explicit and brute force method to
        ! generating all excitations in proportion to their connection
        ! strength. This demonstrates the maximum possible value of tau that
        ! can be used.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_hel_weighted'

        integer :: nsing, ndoub, nexcit

        ! Count how many singles and doubles there are!
        call CountExcitations3 (nI, 3, nsing, ndoub)
        nexcit = nsing + ndoub

        call gen_excit_hel_local (nI, ilutI, nJ, ilutJ, exFlag, ic, ExcitMat, &
                                  tParity, pGen, HElGen, store, nexcit)

    end subroutine


    subroutine gen_excit_hel_local (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                    ExcitMat, tParity, pGen, HElGen, store, &
                                    nexcit)

        ! A really laborious, slow, explicit and brute force method to
        ! generating all excitations in proportion to their connection
        ! strength. This demonstrates the maximum possible value of tau that
        ! can be used.

        integer, intent(in) :: nI(nel), exFlag, nexcit
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_hel_weighted'

        integer(n_int) :: iluts(0:NIfTot, nexcit)
        real(dp) :: hels(nexcit), hel_sum, hel_cum
        integer :: excit_count, ex(2,2), i, flag
        logical :: found_all, par

        ! Generate two lists. One with all of the available excitations, and
        ! one with their HElement values
        excit_count = 0
        found_all = .false.
        hel_sum = 0
        nJ = 0
        flag = 3
        ex = 0
        call GenExcitations3(nI, ilutI, nJ, flag, ex, par, found_all, &
                             .false.)
        do while (.not. found_all)
            excit_count = excit_count + 1
            call EncodeBitDet(nJ, iluts(:, excit_count))
            hels(excit_count) = abs(get_helement(nI, nJ))
            hel_sum = hel_sum + hels(excit_count)
            call GenExcitations3(nI, ilutI, nJ, flag, ex, par, found_all, &
                                 .false.)
        end do

        if (excit_count /= nexcit) &
            call stop_all(this_routine,"Incorrect number of excitations found")

        ! Sort the lists!!!
        call sort(hels, iluts)

        ! Pick a random cumulative helement value
        hel_cum = hel_sum * genrand_real2_dSFMT()

        ! Work from the end, and stop when we have got where we need to be!
        do i = nexcit, 1, -1
            hel_cum = hel_cum - hels(i)
            if (hel_cum <= 0) exit
        end do

        ! Just in case we get shafted by rounding errors
        if (i < 1) i = 1

        ! Now we know what we are returning
        ilutJ = iluts(:, i)
        call decode_bit_det (nJ, ilutJ)
        ExcitMat(1,1) = 2
        call GetBitExcitation (ilutI, ilutJ, ExcitMat, tParity)
        ic = FindBitExcitLevel(ilutI, ilutJ, 2)
        pgen = hels(i) / hel_sum

    end subroutine


    !
    !
    ! ---------------------------------------------------------------------
    ! Now we look at Ali's new biased excitation scheme
    ! ---------------------------------------------------------------------
    !
    !
    subroutine init_4ind_bias ()

        character(*), parameter :: this_routine = 'init_4ind_bias'
        integer :: i, j

        ncnt = 0
        nsel = 0

        rand_excit_opp_bias = 1.0_dp

    end subroutine


    subroutine gen_excit_4ind_weighted (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ExcitMat, tParity, pGen, HelGen, store)

        ! TODO: description
        !
        ! n.b. (ij|kl) <= sqrt( (ij|ij) * (kl|kl) )
        !      This provides quite a good description of the large elements

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_4ind_weighted'
        integer :: orb

        ncnt = ncnt + 1

        ! We now use the class counts to do the construction. This is an
        ! O[N] opearation, and gives the number of occupied/unoccupied
        ! spin-orbitals with each given spin/symmetry.
        if (.not. store%tFilled) then
            call construct_class_counts (nI, store%ClassCountOcc, &
                                         store%ClassCountUnocc)
            store%tFilled = .true.
        end if

        ! Choose if we want to do a single or a double excitation
        ! TODO: We can (in principle) re-use this random number by subdivision
        if (genrand_real2_dSFMT() < pSingles) then

            ic = 1
            call gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, &
                                     tParity, pGen)
            pgen = pgen * pSingles

        else

            ! OK, we want to do a double excitation
            ic = 2
            call gen_double_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, tParity, &
                                     pGen, store)
            pgen = pgen * pDoubles

        end if

        if (nJ(1) /= 0) nsel = nsel + 1

    end subroutine


    function get_paired_cc_ind (cc_ind, sym_product, iSpn) result(cc_ret)

        ! Get the paired Class Count index, given a sym product and iSpn
        ! n.b. Ignoring momentum
        !
        ! Returns index == -1 if no pair exists

        integer, intent(in) :: cc_ind, sym_product, iSpn
        integer :: cc_ret
        integer :: sym_i, sym_j, spn_i, spn_j, mom

        ! Get the relevant details about the class count index
        call ClassCountInv (cc_ind, sym_i, spn_i, mom)

        ! Get the paired symmetry
        !sym_j = ieor(sym_i, sym_product)
        sym_j = RandExcitSymLabelProd (sym_i, sym_product)

        ! Consider the paired symmetries permitted
        spn_j = spn_i
        if (iSpn == 2) & ! i.e. if alpha/beta
            spn_j = 3 - spn_i

        ! If we are restricted to alpha/alpha and beta/beta, and the
        ! specified Class Count index doesn't match, then reject this pair
        ! by returning -1
        if ((spn_i == 1 .and. iSpn == 1) .or. &
            (spn_i == 2 .and. iSpn == 3)) then
            cc_ret = -1
            return
        end if

        ! Get the new index
        cc_ret = ClassCountInd (spn_j, sym_j, mom)

    end function


    subroutine gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ex, par, pgen)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NifTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        integer, intent(out) :: ex(2,2)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elec, src, tgt, cc_index

        ! In this version of the excitation generator, we pick an electron
        ! at random. Then we construct a list of connection
        ! strengths to each of the available orbitals in the correct symmetry.

        ! We could pick the electron based on the number of orbitals available.
        ! Currently, it is just picked uniformly.
        elec = 1 + floor(genrand_real2_dSFMT() * nel)
        src = nI(elec)

        ! What is the symmetry category?
        cc_index = ClassCountInd (get_spin(src), G1(src)%Sym%S, 0)

        ! Select the target orbital by approximate connection strength
        tgt = select_orb_sing (nI, ilutI, src, cc_index, pgen)
        if (tgt == 0) then
            nJ(1) = 0
            return
        end if

        ! Construct the new determinant, excitation matrix and parity
        call make_single (nI, nJ, elec, tgt, ex, par)

        ! Update the target ilut
        ilutJ = ilutI
        clr_orb (ilutJ, src)
        set_orb (ilutJ, tgt)

        ! And the generation probability
        pgen = pgen / real(nel, dp)

    end subroutine


    function select_orb_sing (nI, ilut, src, cc_index, pgen) result(orb)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: cc_index, src
        real(dp), intent(out) :: pgen

        real(dp) :: cum_sum, cumulative_arr(OrbClassCount(cc_index)), r
        real(dp) :: cpt_arr(OrbClassCount(cc_index))
        integer :: orb, norb, label_index, orb_index, i, j
        integer :: n_id(nel), id_src, id
        HElement_t :: hel
        
        ! How many orbitals of the correct symmetry are there?
        norb = OrbClassCount(cc_index)
        label_index = SymLabelCounts2(1, cc_index)

        ! Spatial orbital IDs
        n_id = gtID(nI)
        id_src = gtID(src)
        ASSERT(tExch)

        ! Construct the cumulative list of strengths
        cum_sum = 0
        do i = 1, norb

            orb = SymLabelList2(label_index + i - 1)
            hel = 0
            if (IsNotOcc(ilut, orb)) then
                ASSERT(G1(orb)%Ms == G1(src)%Ms)

                ! For now, we want to assume that we can just use the
                ! one-electron integrals as an approximation. If this is
                ! rubbish, then we will need to do more work...
                !cum_sum = cum_sum &
                !        + abs_l1(GetTMATEl(src, orb))

                ! This is based on an extract from sltcnd_1. We can assume
                ! tExch, and SymLabelList2 ensures the spins are equal
                id = gtID(orb)
                do j = 1, nel
                    if (nI(j) == src) cycle
                    hel = hel + get_umat_el (id_src, n_id(j), id, n_id(j))
                    if (is_beta(src) .eqv. is_beta(nI(j))) &
                        hel = hel - get_umat_el (id_src, n_id(j), n_id(j), id)
                end do
                hel = hel + GetTMATEl(src, orb)

            end if

            ! And store the values for later searching
            cpt_arr(i) = abs_l1(hel)
            cum_sum = cum_sum + cpt_arr(i)
            cumulative_arr(i) = cum_sum

        end do

        ! Select a particular orbital to use, or abort.
        if (cum_sum == 0) then
            orb = 0
        else
            r = genrand_real2_dSFMT() * cum_sum
            orb_index = binary_search_first_ge(cumulative_arr, r)
            orb = SymLabelList2(label_index + orb_index - 1)

            ! And the impact on the generation probability
            pgen = cpt_arr(orb_index) / cum_sum
        end if

    end function



    
    subroutine gen_double_4ind_ex (nI, ilutI, nJ, ilutJ, ex, par, pgen, store)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        type(excit_gen_store_type), intent(in) :: store
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elecs(2), orbs(2), src(2)
        real(dp) :: int_cpt(2), cum_sum(2), cc_cum(ScratchSize), cc_tot, r
        integer :: sym_product, ispn, sum_ml
        integer :: cc_i, cc_j

        ! Initially, select a pair of electrons
        call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)

        ! Construct the list of possible symmetry pairs listed according to
        ! a unique choice of symmetry A (we enforce that sym(A) < sym(B)).
        cc_tot = 0
        do cc_i = 1, ScratchSize

            cc_j = get_paired_cc_ind (cc_i, sym_product, iSpn)

            ! As cc_i > 0, the following test also excludes any rejected
            ! pairings (where cc_j == -1).
            if (cc_j >= cc_i) then
                cc_tot = cc_tot + store%ClassCountUnocc(cc_i) &
                                * store%ClassCountUnocc(cc_j)
                !cc_tot = cc_tot + OrbClassCount(cc_i) * OrbClassCount(cc_j)
            end if
            cc_cum(cc_i) = cc_tot
        end do
    
        ! Pick a specific pairing of spin-symmetries
        r = genrand_real2_dSFMT() * cc_tot
        cc_i = binary_search_first_ge (cc_cum, r)
        cc_j = get_paired_cc_ind (cc_i, sym_product, iSpn)
        pgen = pgen * store%ClassCountUnocc(cc_i) &
                    * store%ClassCountUnocc(cc_j) / cc_tot
!        pgen = pgen * OrbClassCount(cc_i) * OrbClassCount(cc_j) / cc_tot

        ! Select the two orbitals
        orbs(1) = select_orb (ilutI, src, cc_i, -1, int_cpt(1), cum_sum(1))
        if (orbs(1) /= 0) &
            orbs(2) = select_orb (ilutI, src, cc_j, orbs(1), int_cpt(2), &
                                  cum_sum(2))

        ! Just as a check, abort any excitation where we have excited two
        ! electrons into the same orbital
        if (any(orbs == 0)) then
            nJ(1) = 0
            return
        end if

        ! Adjust the probabilities. If we are selecting two orbitals from
        ! the same list, then we could have selected them either way around.
        if (cc_i == cc_j) then
            pgen = pgen * ( &
                   (int_cpt(1) / cum_sum(1) * int_cpt(2) / cum_sum(2)) &
                 + (int_cpt(2) / cum_sum(1) * &
                    int_cpt(1) / (cum_sum(1) - int_cpt(2))))
        else
            pgen = pgen * (int_cpt(1) / cum_sum(1)) * (int_cpt(2) / cum_sum(2))
        end if

        ! And generate the actual excitation.
        call make_double (nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                          ex, par)
        ilutJ = ilutI
        clr_orb (ilutJ, src(1))
        clr_orb (ilutJ, src(2))
        set_orb (ilutJ, orbs(1))
        set_orb (ilutJ, orbs(2))

    end subroutine


    subroutine pick_biased_elecs (nI, elecs, src, sym_prod, ispn, sum_ml, pgen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), src(2), sym_prod, ispn, sum_ml
        real(dp), intent(out) :: pgen

        real(dp) :: ntot, r
        integer :: al_req, be_req, al_num(2), be_num(2), elecs_found, i, idx
        integer :: al_count, be_count

        ntot = AA_elec_pairs + BB_elec_pairs &
              + AB_elec_pairs * rand_excit_opp_bias

        ! The overall generation probability for this section is remarkably
        ! simple!
        pgen = 1.0_dp / ntot

        ! We want to have the n'th alpha, or beta electrons in the determinant
        ! Select them according to the availability of pairs (and the
        ! weighting of opposite-spin pairs relative to same-spin ones).
        r = genrand_real2_dSFMT() * ntot
        if (r < AA_elec_pairs) then
            al_req = 2
            be_req = 0
            idx = floor(r)
            al_num(1) = ceiling((1 + sqrt(9 + 8*real(idx, dp))) / 2)
            al_num(2) = idx + 1 - ((al_num(1) - 1) * (al_num(1) - 2)) / 2
            iSpn = 3
        else if (r < par_elec_pairs) then
            al_req = 0
            be_req = 2
            idx = floor(r - AA_elec_pairs)
            be_num(1) = ceiling((1 + sqrt(9 + 8*real(idx, dp))) / 2)
            be_num(2) = idx + 1 - ((be_num(1) - 1) * (be_num(1) - 2)) / 2
            iSpn = 1
        else
            al_req = 1
            be_req = 1
            pgen = pgen * rand_excit_opp_bias
            idx = floor((r - par_elec_pairs) / rand_excit_opp_bias)
            al_num(1) = 1 + mod(idx, nOccAlpha)
            be_num(1) = 1 + floor(idx / real(nOccAlpha,dp))
            iSpn = 2
        end if

        ! Loop through the determiant, and select the relevant orbitals from
        ! it.
        ! This is not as clean as the implementation in PickElecPair, which
        ! doesn't consider the spins. Would work MUCH better in an
        ! environment where the alpha/beta spins were stored seperately...
        al_count = 0
        be_count = 0
        elecs_found = 0
        do i = 1, nel
            if (is_alpha(nI(i))) then
                al_count = al_count + 1
                if (al_req > 0) then
                    if (al_count == al_num(al_req)) then
                        elecs_found = elecs_found + 1
                        elecs(elecs_found) = i
                        al_req = al_req - 1
                    end if
                end if
            else
                be_count = be_count + 1
                if (be_req > 0) then
                    if (be_count == be_num(be_req)) then
                        elecs_found = elecs_found + 1
                        elecs(elecs_found) = i
                        be_req = be_req - 1
                    end if
                end if
            end if
            if (al_req == 0 .and. be_req == 0) exit
        end do

        ! Generate the orbitals that are being considered
        src = nI(elecs)

        ! The Ml value is obtained from the orbitals
        sum_ml = sum(G1(src)%Ml)

        ! Get the symmetries
        sym_prod = RandExcitSymLabelProd (int(G1(src(1))%Sym%S), &
                                          int(G1(src(2))%Sym%S))

    end subroutine


    function select_orb (ilut, src, cc_index, orb_pair, cpt, cum_sum) &
            result(orb)

        ! For an excitation from electrons 1,2, consider all of the pairs
        ! (e1 a|e1 a) == <e1 e1 | a a> and (e2 a | e2 a) == <e2 e2 | a a>
        ! as appropriate to bias selection of the electron

        integer, intent(in) :: src(2), orb_pair, cc_index
        real(dp), intent(out) :: cpt, cum_sum
        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer :: orb

        integer :: label_index, orb_index, norb, i, orbid, srcid(2), ms
        integer :: src_orb, src_id

        ! Our biasing arrays must consider all of the possible orbitals with
        ! the correct symmetry.
        real(dp) :: cumulative_arr(OrbClassCount(cc_index)), r

        ! How many orbitals are there with the given symmetry?
        !cc_index = ClassCountInd(spin, sym, 0)
        label_index = SymLabelCounts2(1, cc_index)
        norb = OrbClassCount(cc_index)

        ! Here we construct a cumulative list of the absolute values of the
        ! integrals between the source electrons (orbitals), and the available
        ! orbitals in the specified symmetry class.
        ! --> These can then be binary searched to pick an orbital
        ! --> The L1-norm is used when complex integrals are being used, as
        !     the cumulative spawning rate is related to the sum of the values
        !     rather than the norm of the complex number.

        if (G1(src(1))%Ms == G1(src(2))%Ms) then

            ! Both of the electrons have the same spin. Therefore we need to
            ! include both electron-hole interactions.
            cum_sum = 0
            srcid = gtID(src)
            do i = 1, norb
            
                orb = SymLabelList2(label_index + i - 1)
                if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                    orbid = gtID(orb)
                    cum_sum = cum_sum &
                    + sqrt(abs_l1(get_umat_el(srcid(1), srcid(1), orbid, orbid)))&
                    + sqrt(abs_l1(get_umat_el(srcid(2), srcid(2), orbid, orbid)))
                end if
                cumulative_arr(i) = cum_sum

            end do

        else
            
            ! The two electrons have differing spin. Therefore, only the 
            ! electron-hole interaction with the same spin is required.
            ms = class_count_ms(cc_index)
            if (ms == G1(src(1))%Ms) then
                src_orb = src(1)
            else
                src_orb = src(2)
            end if

            cum_sum = 0
            src_id = gtID(src_orb)
            do i = 1, norb

                orb = SymLabelList2(label_index + i - 1)
                if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                    orbid = gtID(orb)
                    cum_sum = cum_sum &
                    + sqrt(abs_l1(get_umat_el(src_id, src_id, orbid, orbid)))
                end if
                cumulative_arr(i) = cum_sum

            end do

        end if


        ! If there are no available orbitals to pair with, we need to abort
        if (cum_sum == 0) then
            orb = 0
            return
        end if

        ! Binary search within this list to choose a value.
        r = genrand_real2_dSFMT() * cum_sum
        orb_index = binary_search_first_ge(cumulative_arr(1:norb), r)

        ! And return the relevant value.
        ! TODO: We don't need to call get_umat_ell. Already known.
        orb = SymLabelList2(label_index + orb_index - 1)
        orbid = gtID(orb)
        if (G1(src(1))%Ms == G1(src(2))%Ms) then
            cpt = sqrt(abs_l1(get_umat_el(srcid(1), srcid(1), orbid, orbid))) + &
                  sqrt(abs_l1(get_umat_el(srcid(2), srcid(2), orbid, orbid)))
        else
            cpt = sqrt(abs_l1(get_umat_el(src_id, src_id, orbid, orbid)))
        end if

    end function


    subroutine gen_excit_4ind_reverse (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                       ExcitMat, tParity, pGen, HelGen, store)

        ! TODO: description
        !
        ! n.b. (ij|kl) <= sqrt( (ij|ij) * (kl|kl) )
        !      This provides quite a good description of the large elements

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_4ind_reverse'

        integer :: orb

        if (genrand_real2_dSFMT() < pSingles) then

            ! Do we want to use the forward, or reverse, singles generator?
            ic = 1
            call gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, &
                                      tParity, pGen)
            pgen = pgen * pSingles

        else

            ! OK, we want to do a double excitation
            ic = 2
            call gen_double_4ind_rev (nI, ilutI, nJ, ilutJ, ExcitMat, tParity,&
                                      pgen)
            pgen = pgen * pDoubles

        end if

    end subroutine


    subroutine gen_single_4ind_rev (nI, ilutI, nJ, ilutJ, ex, par, pgen)

        ! Select a vacant hole, and then select an electron that is connected
        ! to it.

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NifTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        integer, intent(out) :: ex(2,2)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: nholes, src, elec, tgt

        ! In this version of the excitation generator, we pick a hole at
        ! random. Then we construct a list of connection strengths to each of
        ! the available electrons.

        ! Select a hole (uniformly) at random.
        tgt = pick_hole (ilutI, -1)
        nholes = nbasis - nel

        ! Select the source electron by connection strength
        elec = select_elec_sing (nI, tgt, src, pgen)
        if (elec == 0) then
            nJ(1) = 0
            return
        end if

        ! Construct the new determinant, electron matrix, and parity
        call make_single (nI, nJ, elec, tgt, ex, par)

        ! Update the target ilut
        ilutJ = ilutI
        clr_orb (ilutJ, src)
        set_orb (ilutJ, tgt)

        ! And the generation probability
        pgen = pgen / real(nholes, dp)

    end subroutine


    function select_elec_sing (nI, tgt, src, pgen) result(elec)

        integer, intent(in) :: nI(nel), tgt
        integer, intent(out) :: src
        real(dp), intent(out) :: pgen
        integer :: elec

        real(dp) :: cum_sum, cum_arr(nel), cpt_arr(nel), r
        integer :: cc_index, n_id(nel), id_tgt, id_src, cc_src, j, src_elec
        HElement_t :: hel

        ! What is the symmetry category we are considering?
        cc_index = ClassCountInd (get_spin(tgt), G1(tgt)%Sym%S, 0)

        ! Spatial orbital IDs
        n_id = gtID(nI)
        id_tgt = gtID(tgt)
        ASSERT(tExch)

        ! Construct the cumulative list of strengths
        cum_sum = 0
        do src_elec = 1, nel

            ! What is the symmetry/spin class of the source electron
            src = nI(src_elec)
            cc_src = ClassCountInd (get_spin(src), G1(src)%Sym%S, 0)

            ! If the symmetry/spin are not the same, then the matrix element
            ! will be zero --> don't calculate it!
            hel = 0
            if (cc_src == cc_index) then

                ! Construct the hamiltonian matrix element (ignoring overall
                ! sign). We can assume tExch, and we have just enforced that
                ! the spins and symmetries are equal.
                id_src = n_id(src_elec)
                do j = 1, nel
                    if (j == src_elec) cycle
                    hel = hel + get_umat_el (id_src, n_id(j), id_tgt, n_id(j))
                    if (is_beta(src) .eqv. is_beta(nI(j))) &
                        hel = hel - get_umat_el (id_src, n_id(j), n_id(j), id_tgt)
                end do
                hel = hel + GetTMATEl(src, tgt)

            end if

            ! And store the values for later searching
            cpt_arr(src_elec) = abs_l1(hel)
            cum_sum = cum_sum + cpt_arr(src_elec)
            cum_arr(src_elec) = cum_sum

        end do
        
        ! Select a particulor electron, or abort
        if (cum_sum == 0) then
            elec = 0
        else
            r = genrand_real2_dSFMT() * cum_sum
            elec = binary_search_first_ge(cum_arr, r)
            src = nI(elec)

            ! And the impact on the generation probability
            pgen = cpt_arr(elec) / cum_sum
        end if

    end function


    subroutine gen_double_4ind_rev (nI, ilutI, nJ, ilutJ, ex, par, pgen)

        ! Here we do the 'opposite' of gen_double_4ind_ex. We select two
        ! holes, and then pick the electrons to excite based on the strengths
        ! of the connecting matrix elements. Due to the substantially reduced
        ! number of terms, we can actually calculate the EXACT matrix elements
        ! for this, which is quite nice!

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: src(2), tgt(2), elecs(2), id_tgt(2), id(2)
        integer :: sym_prod, ispn, e_ispn, e_sym_prod, idx, i, j
        HElement_t :: hel
        real(dp) :: cum_arr(nel * (nel - 1) / 2)
        real(dp) :: val_arr(nel * (nel - 1) / 2)
        real(dp) :: val, cum_val, r

        ! Select a pair of holes (vacant orbitals)
        call pick_hole_pair_biased (ilutI, tgt, pgen, sym_prod, ispn)
        id_tgt = gtID(tgt)

        ! Now we want to consider all of the possible pairs of electrons
        idx = 0
        cum_val = 0
        do i = 2, nel
            src(1) = nI(i)
            do j = 1, i-1

                ! Get the symmetries
                src(2) = nI(j)
                e_ispn = get_ispn(src(1), src(2))
                e_sym_prod = RandExcitSymLabelProd(int(G1(src(1))%Sym%s), &
                                                   int(G1(src(2))%Sym%s))

                ! Get the weight (HElement) associated with the elecs/holes
                if (e_ispn == iSpn .and. e_sym_prod == sym_prod) then
                    hel = 0
                    id = gtID(src)
                    if (G1(src(1))%Ms == G1(tgt(1))%Ms .and.&
                        G1(src(2))%Ms == G1(tgt(2))%Ms) &
                        hel = hel + get_umat_el (id(1), id(2), id_tgt(1), &
                                                 id_tgt(2))
                    if (G1(src(1))%Ms == G1(tgt(2))%Ms .and.&
                        G1(src(2))%Ms == G1(tgt(1))%Ms) &
                        hel = hel - get_umat_el (id(1), id(2), id_tgt(2), &
                                                 id_tgt(1))
                    val = abs_l1(hel)
                else
                    val = 0
                end if

                ! Store the terms
                idx = idx + 1
                cum_val = cum_val + val
                val_arr(idx) = val
                cum_arr(idx) = cum_val

            end do
        end do

        ! If there are no choices, then we have to abort...
        if (cum_val == 0) then
            nJ(1) = 0
            return
        end if

        ! Pick the electron pair we are going to use.
        r = genrand_real2_dSFMT() * cum_val
        idx = binary_search_first_ge (cum_arr, r)
        elecs(2) = ceiling((1 + sqrt(1 + 8*real(idx, dp))) / 2)
        elecs(1) = idx - ((elecs(2) - 1) * (elecs(2) - 2)) / 2

        ! Get the generation probability
        pgen = pgen * val_arr(idx) / cum_val

!>>>!        write(6,*) 'ELECS', elecs
!>>>!        write(6,*) 'SRC', nI(elecs)
!>>>!        write(6,*) 'TGT', tgt

        ! Generate the new determinant and ilut
        call make_double (nI, nJ, elecs(1), elecs(2), tgt(1), tgt(2), ex, par)
        ilutJ = ilutI
        clr_orb (ilutJ, nI(elecs(1)))
        clr_orb (ilutJ, nI(elecs(2)))
        set_orb (ilutJ, tgt(1))
        set_orb (ilutJ, tgt(2))


    end subroutine

    function get_ispn (orbi, orbj) result(ispn)

        ! ispn == 1 --> beta/beta
        ! ispn == 2 --> alpha/beta
        ! ispn == 3 --> alpha/alpha

        integer, intent(in) :: orbi, orbj
        integer :: ispn

        if (is_alpha(orbi) .eqv. is_alpha(orbj)) then
            if (is_alpha(orbi)) then
                ispn = 3
            else
                ispn = 1
            end if
        else
            ispn = 2
        end if

    end function


    subroutine pick_hole_pair_biased (ilut, orbs, pgen, sym_prod, ispn)
        
        ! This is a biased version of pick_hole_pair below. Whereas the
        ! other function picks hole pairs entirely uniformly throughout the
        ! space, this function biases that selection towards picking opposite
        ! spin pairs (according to the biasing factor rand_excit_opp_bias).

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: orbs(2), sym_prod, ispn
        real(dp), intent(out) :: pgen

        real(dp) :: ntot, r
        integer :: spn(2)

        ntot = par_hole_pairs + AB_hole_pairs * rand_excit_opp_bias

        ! The overall generation probability is remarkably simple!
        pgen = 1.0_dp / ntot

        ! Pick what sort of hole pair we want. AA, BB, or AB.
        r = genrand_real2_dSFMT() * ntot
        if (r < AA_hole_pairs) then
            ! alpha/alpha
            iSpn = 3
            spn = (/1, 1/)
        else if (r < par_hole_pairs) then
            ! beta/beta
            iSpn = 1
            spn = (/2, 2/)
        else
            ! alpha/beta
            iSpn = 2
            spn = (/1, 2/)

            ! We need to adjust the generation probability to account for the
            ! selection bias
            pgen = pgen * rand_excit_opp_bias
        end if

        ! Select orbitals at random with the given spins
        orbs(1) = pick_hole_spn (ilut, spn(1), -1)
        orbs(2) = pick_hole_spn (ilut, spn(2), orbs(1))

        ! What are the symmetry/spin properties of this pick?
        sym_prod = RandExcitSymLabelProd(int(G1(orbs(1))%Sym%S), &
                                         int(G1(orbs(2))%Sym%S))
        ASSERT(iSpn == get_ispn(orbs(1), orbs(2)))

    end subroutine

    function pick_hole_spn (ilut, spn, exclude_orb) result(orb)

        ! Pick an unoccupied orbital at random, uniformly, from the available
        ! basis set (excluding a selected orbital if desired), with the given
        ! spin

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: exclude_orb, spn
        integer :: orb, attempts, norb
        integer, parameter :: max_attempts = 1000
        character(*), parameter :: t_r = 'pick_hole'

        ASSERT(.not. btest(nbasis, 0))
        norb = nbasis / 2

        ! Draw orbitals randomly until we find the first orbital
        attempts = 0
        do while(.true.)

            ! Select an orbital randomly with the correct spin
            if (spn == 1) then
                ! alpha
                orb = 2 + 2 * int(genrand_real2_dSFMT() * norb)
            else
                ! beta
                orb = 1 + 2 * int(genrand_real2_dSFMT() * norb)
            end if

            ! Accept this selection if it is a hole.
            if (IsNotOcc(ilut, orb) .and. orb /= exclude_orb) exit

            ! Just in case, add a check
            if (attempts > max_attempts) then
                write(6,*) 'Unable to find unoccupied orbital'
                call writebitdet (6, ilut, .true.)
                call stop_all(t_r, 'Out of attempts')
            end if
        end do

    end function

    
    subroutine pick_hole_pair (ilut, orbs, pgen, sym_prod, ispn)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: orbs(2), sym_prod, ispn
        real(dp), intent(out) :: pgen
        integer :: attempts, nvac, nchoose

        ! Get two unoccupied orbitals
        orbs(1) = pick_hole(ilut, -1)
        orbs(2) = pick_hole(ilut, orbs(1))

        ! What is the generation probability of picking this selection?
        nvac = nbasis - nel
        nchoose = nvac * (nvac - 1) / 2
        pgen = 1.0_dp / real(nchoose, dp)

        ! What are the symmetry/spin properties of this pick?
        sym_prod = RandExcitSymLabelProd (int(G1(orbs(1))%Sym%S), &
                                          int(G1(orbs(2))%Sym%S))
        ispn = get_ispn(orbs(1), orbs(2))

    end subroutine


    function pick_hole (ilut, exclude_orb) result(orb)

        ! Pick an unoccupied orbital at random, uniformly, from the available
        ! basis set (excluding a selected orbital if desired).

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: exclude_orb
        integer :: orb, attempts
        integer, parameter :: max_attempts = 1000
        character(*), parameter :: t_r = 'pick_hole'

        ! Draw orbitals randomly until we find the first orbital
        attempts = 0
        do while(.true.)

            ! Select an orbital randomly. Accept it if it is a hole.
            orb = 1 + int(genrand_real2_dSFMT() * nbasis)
            if (IsNotOcc(ilut, orb) .and. orb /= exclude_orb) exit

            ! Just in case, add a check
            if (attempts > max_attempts) then
                write(6,*) 'Unable to find unoccupied orbital'
                call writebitdet (6, ilut, .true.)
                call stop_all(t_r, 'Out of attempts')
            end if
        end do

    end function


    subroutine test_excit_gen_4ind (ilut, iterations)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: iterations
        character(*), parameter :: this_routine = 'test_excit_gen_4ind'

        integer :: src_det(nel), det(nel), nsing, ndoub, nexcit, ndet, ex(2,2)
        integer :: flag, ngen, pos, iunit, i, ic
        type(excit_gen_store_type) :: store
        integer(n_int) :: tgt_ilut(0:NifTot)
        integer(n_int), allocatable :: det_list(:,:)
        real(dp), allocatable :: contrib_list(:)
        logical, allocatable :: generated_list(:)
        logical :: found_all, par
        real(dp) :: contrib, pgen
        HElement_t :: helgen

        ! Decode the determiant
        call decode_bit_det (src_det, ilut)

        ! Initialise
        call init_excit_gen_store (store)

        ! How many connected determinants are we expecting?
        call CountExcitations3 (src_det, 3, nsing, ndoub)
        nexcit = nsing + ndoub
        allocate(det_list(0:NIfTot, nexcit))

        ! Loop through all of the possible excitations
        ndet = 0
        found_all = .false.
        ex = 0
        flag = 3
        write(6,'("*****************************************")')
        write(6,'("Enumerating excitations")')
        write(6,'("Starting from: ")', advance='no')
        call write_det (6, src_det, .true.)
        write(6,*) 'Expecting ', nexcit, "excitations"
        call GenExcitations3 (src_det, ilut, det, flag, ex, par, found_all, &
                              .false.)
        do while (.not. found_all)
            ndet = ndet + 1
            call EncodeBitDet (det, det_list(:,ndet))

            call GenExcitations3 (src_det, ilut, det, flag, ex, par, &
                                  found_all, .false.)
        end do
        if (ndet /= nexcit) &
            call stop_all(this_routine,"Incorrect number of excitations found")

        ! Sort the dets, so they are easy to find by binary searching
        call sort(det_list, ilut_lt, ilut_gt)

        ! Lists to keep track of things
        allocate(generated_list(nexcit))
        allocate(contrib_list(nexcit))
        generated_list = .false.
        contrib_list = 0

        ! Repeated generation, and summing-in loop
        ngen = 0
        contrib = 0
        do i = 1, iterations
            if (mod(i, 10000) == 0) &
                write(6,*) i, '/', iterations, ' - ', contrib / real(ndet*i)

            call gen_excit_4ind_weighted (src_det, ilut, det, tgt_ilut, 3, &
                                          ic, ex, par, pgen, helgen, store)
            if (det(1) == 0) cycle

            call EncodeBitDet (det, tgt_ilut)
            pos = binary_search(det_list, tgt_ilut, NIfD+1)
            if (pos < 0) then
                write(6,*) det
                write(6,'(b64)') tgt_ilut(0)
                write(6,*) 'FAILED DET', tgt_ilut
                call writebitdet(6, tgt_ilut, .true.)
                call stop_all(this_routine, 'Unexpected determinant generated')
            else
                generated_list(pos) = .true.

                ! Count this det, and sum in its contribution.
                ngen = ngen + 1
                contrib = contrib + 1.0_dp / pgen
                contrib_list(pos) = contrib_list(pos) + 1.0_dp / pgen
            end if
        end do

        ! How many of the iterations generated a good det?
        write(6,*) ngen, " dets generated in ", iterations, " iterations."
        write(6,*) 100_dp * (iterations - ngen) / real(iterations), &
                   '% abortion rate'
        ! Contribution averages
        write(6, '("Averaged contribution: ", f15.10)') &
                contrib / real(ndet * iterations)

        ! Output the determinant specific contributions
        iunit = get_free_unit()
        open(iunit, file="contribs_4ind", status='unknown')
        do i = 1, ndet
            call writebitdet(iunit, det_list(:,i), .false.)
            write(iunit, *) contrib_list(i) / real(iterations, dp)
        end do
        close(iunit)

        ! Check that all of the determinants were generated!!!
        if (.not. all(generated_list)) then
            write(6,*) count(.not.generated_list), '/', size(generated_list), &
                       'not generated'
            do i = 1, ndet
                if (.not. generated_list(i)) &
                    call writebitdet(6, det_list(:,i), .true.)
            end do
            call stop_all(this_routine, "Determinant not generated")
        end if
        if (any(abs(contrib_list / iterations - 1.0) > 0.01)) &
            call stop_all(this_routine, "Insufficiently uniform generation")

        ! Clean up
        deallocate(det_list)
        deallocate(contrib_list)
        deallocate(generated_list)

    end subroutine

end module
