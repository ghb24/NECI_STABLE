#include "macros.h"

module excit_gens_int_weighted

    use SystemData, only: nel, nbasis, nOccAlpha, nOccBeta, G1, nOccAlpha, &
                          nOccBeta, tExch, AA_elec_pairs, BB_elec_pairs, &
                          AB_elec_pairs, par_elec_pairs, AA_hole_pairs, &
                          par_hole_pairs, AB_hole_pairs, iMaxLz, &
                          tGen_4ind_part_exact, tGen_4ind_lin_exact, &
                          tGen_4ind_unbound, t_iiaa, t_ratio, UMatEps, tGUGA
    use CalcData, only: matele_cutoff, t_matele_cutoff
    use SymExcit3, only: CountExcitations3, GenExcitations3
    use SymExcitDataMod, only: SymLabelList2, SymLabelCounts2, OrbClassCount, &
                               pDoubNew, ScratchSize, SpinOrbSymLabel, &
                               SymInvLabel
    use sym_general_mod, only: ClassCountInd, ClassCountInv, class_count_ms, &
                               class_count_ml
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles, pParallel
    use dSFMT_interface, only: genrand_real2_dSFMT
    use Determinants, only: get_helement, write_det
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, ilut_lt, ilut_gt
    use bit_rep_data, only: NIfTot, NIfD, test_flag
    use bit_reps, only: decode_bit_det, get_initiator_flag
    use symdata, only: nSymLabels
    use procedure_pointers, only: get_umat_el
    use UMatCache, only: gtid, UMat2d
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
    use LoggingData, only: t_log_ija, ija_bins_para, ija_bins_anti, ija_thresh, &
                           ija_orbs_para, ija_orbs_anti, ija_bins_sing, ija_orbs_sing

    implicit none
    save

contains

    subroutine gen_excit_hel_weighted (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                       ExcitMat, tParity, pGen, HelGen, store, part_type)

        ! A really laborious, slow, explicit and brute force method to
        ! generating all excitations in proportion to their connection
        ! strength. This demonstrates the maximum possible value of tau that
        ! can be used.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

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
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_hel_weighted'

        integer(n_int) :: iluts(0:NIfTot, nexcit)
        real(dp) :: hels(nexcit), hel_sum, hel_cum
        integer :: excit_count, ex(2,2), i, flag
        logical :: found_all, par

        ! Generate two lists. One with all of the available excitations, and
        ! one with their HElement values
        HElGen = HEl_zero
        excit_count = 0
        found_all = .false.
        hel_sum = 0
        nJ = 0
        flag = 3
        ex = 0
        call GenExcitations3(nI, ilutI, nJ, flag, ex, par, found_all, &
                             .false.)
        if (tGUGA) then 
            call stop_all(this_routine, "modify get_helement for GUGA")
        end if
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

    subroutine gen_excit_4ind_weighted (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ExcitMat, tParity, pGen, HelGen, store, part_type)

        ! TODO: description
        !
        ! n.b. (ij|kl) <= sqrt( (ij|ij) * (kl|kl) )
        !      This provides quite a good description of the large elements

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(in), optional :: part_type

        character(*), parameter :: this_routine = 'gen_excit_4ind_weighted'
        integer :: orb
        real(dp) :: pgen2

        HElGen = HEl_zero

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

!             call write_det(6, nJ, .true.)
!             print *, "pgen singles: ", pgen

        else

            ! OK, we want to do a double excitation
            ic = 2
            call gen_double_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, tParity, &
                                     pGen, store)
            pgen = pgen * pDoubles

!             call write_det(6, nJ, .true.)
!             print *, "pgen doubles: ", pgen

        end if

        ! And a careful check!
#ifdef __DEBUG
        if (.not. IsNullDet(nJ)) then
            pgen2 = calc_pgen_4ind_weighted(nI, ilutI, ExcitMat, ic, &
                                            store%ClassCountUnocc)
            if (abs(pgen - pgen2) > 1.0e-6_dp) then
                write(6,*) 'Calculated and actual pgens differ.'
                write(6,*) 'This will break HPHF calculations'
                call write_det(6, nI, .false.)
                write(6, '(" --> ")', advance='no')
                call write_det(6, nJ, .true.)
                write(6,*) 'Excitation matrix: ', ExcitMat(1,1:ic), '-->', &
                           ExcitMat(2,1:ic)
                write(6,*) 'Generated pGen:  ', pgen
                write(6,*) 'Calculated pGen: ', pgen2
                call stop_all(this_routine, "Invalid pGen")
            end if
        end if
#endif

    end subroutine



    function calc_pgen_4ind_weighted (nI, ilutI, ex, ic, ClassCountUnocc) &
            result(pgen)

        ! What is the probability of the excitation _from_ determinant nI
        ! described by the excitation matrix ex, and the excitation level ic,
        ! being generated according to the 4ind_weighted excitaiton generator?

        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer, intent(in) :: ClassCountUnocc(ScratchSize)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'calc_pgen_4ind_weighted'

        integer :: cc_index, tgt(2), sum_ml, iSpn, sym_product
        integer :: cc_i, cc_j, cc_i_final, cc_j_final, cc_a, cc_b
        real(dp) :: cpt, cpt_tgt, cum_sum, cum_sums(2), int_cpt(2), ntot
        real(dp) :: cpt_pair(2), sum_pair(2)
        HElement_t(dp) :: hel


        if (ic == 1) then

            ! Singles are calculated in the same way for the _ex and the
            ! _reverse excitation generators
            pgen = pSingles * pgen_single_4ind (nI, ilutI, ex(1,1), ex(2,1))

        else if (ic == 2) then

            ! Obviously bias by pDoubles
            pgen = pDoubles

            ! We want to select a pair of electrons in a way which is biased
            ! towards pairs with opposing spins. This takes a simple form.
            if (is_alpha(ex(1,1)) .eqv. is_alpha(ex(1,2))) then
                pgen = pgen * pParallel / par_elec_pairs
                if (is_alpha(ex(1,1))) then
                    iSpn = 3
                else
                    iSpn = 1
                end if
            else
                pgen = pgen * (1.0_dp - pParallel) / AB_elec_pairs
                iSpn = 2
            end if

            ! What is the likelihood of picking the given symmetry?
            sym_product = RandExcitSymLabelProd(SpinOrbSymLabel(ex(1,1)), &
                                                SpinOrbSymLabel(ex(1,2)))
            sum_ml = sum(G1(ex(1,:))%Ml)
            cc_a = ClassCountInd(get_spin(ex(2,1)), SpinOrbSymLabel(ex(2,1)), &
                                 G1(ex(2,1))%Ml)
            cc_b = ClassCountInd(get_spin(ex(2,2)), SpinOrbSymLabel(ex(2,2)), &
                                 G1(ex(2,2))%Ml)
            cc_i_final = min(cc_a, cc_b)
            cc_j_final = max(cc_a, cc_b)
            cum_sum = 0
            do cc_i = 1, ScratchSize
                cc_j = get_paired_cc_ind (cc_i, sym_product, sum_ml, iSpn)

                ! We restrict cc_i > cc_j, and rejected pairings where.
                if (cc_j >= cc_i) then
                    cpt = ClassCountUnocc(cc_i) * ClassCountUnocc(cc_j)
                    cum_sum = cum_sum + cpt
                    if (cc_i == cc_i_final) then
                        ASSERT(cc_j == cc_j_final)
                        cpt_tgt = cpt

                        ! ENSURE that the tgt orbitals match cc_i, cc_j as
                        ! the choice of orbitals following is not symmetric.
                        if (cc_i == cc_a) then
                            tgt = ex(2,:)
                        else
                            tgt = [ex(2,2), ex(2,1)]
                        end if
                    end if
                end if
            end do

            ! Andjust the probability for this symmetry stuff
            if (cum_sum < EPS) then
                pgen = 0
                return
            else
                pgen = pgen * cpt_tgt / cum_sum
            end if

            ! What is the likelihood of selecting the first orbital? And the
            ! second, given the first?
            call pgen_select_orb (ilutI, ex(1,:), -1, tgt(1), int_cpt(1), &
                                  cum_sums(1))
            call pgen_select_orb (ilutI, ex(1,:), tgt(1), tgt(2), int_cpt(2), &
                                  cum_sums(2))

            ! Deal with cases when there are no available excitations
            ! with the given pathway.
            if (any(cum_sums < EPS)) then
                cum_sums = 1.0
                int_cpt = 0.0
            end if

            if (cc_i_final == cc_j_final) then
                if (tGen_4ind_lin_exact .or. tGen_4ind_part_exact) then
                    call pgen_select_orb(ilutI, ex(1,:), -1, tgt(2), &
                                         cpt_pair(1), sum_pair(1))
                    call pgen_select_orb(ilutI, ex(1,:), tgt(2), tgt(1), &
                                         cpt_pair(2), sum_pair(2))
                else
                    cpt_pair(1) = int_cpt(2)
                    cpt_pair(2) = int_cpt(1)
                    sum_pair(1) = cum_sums(1)
                    sum_pair(2) = cum_sums(1) - int_cpt(2)
                end if
                
                ! Deal with cases when there are no available excitations
                ! with the given pathway.
                If (any(sum_pair < EPS)) then
                    sum_pair = 1.0
                    cpt_pair = 0.0
                end if

                pgen = pgen * (product(int_cpt) / product(cum_sums) + &
                               product(cpt_pair) / product(sum_pair))
            else
                pgen = pgen * (int_cpt(1) / cum_sums(1)) &
                            * (int_cpt(2) / cum_sums(2))
            end if

        else
            ! IC /= 1, 2 --> not connected by the excitation generator.
            ! 
            ! If we are hitting here, then earlier checks have failed. We can
            ! return the correct (zero) value at runtime, but really this
            ! should be fixed elsewhere
            ASSERT(.false.)
            pgen = 0.0_dp
        end if

    end function



    function get_paired_cc_ind (cc_ind, sym_product, sum_ml, iSpn) &
            result(cc_ret)

        ! Get the paired Class Count index, given a sym product and iSpn
        ! n.b. Ignoring momentum
        !
        ! Returns index == -1 if no pair exists

        integer, intent(in) :: cc_ind, sym_product, iSpn, sum_ml
        integer :: cc_ret
        integer :: sym_i, sym_j, spn_i, spn_j, mom_i, mom_j
        character(*), parameter :: this_routine = 'get_paired_cc_ind'

        ! Get the relevant details about the class count index
        call ClassCountInv (cc_ind, sym_i, spn_i, mom_i)

        ! Get the paired symmetry
        !sym_j = ieor(sym_i, sym_product)
        sym_j = RandExcitSymLabelProd (SymInvLabel(sym_i), sym_product)

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

        ! The sum of the momentum terms must equal sum_ml. If there isn't a
        ! paired orbital that gives the relevant momentum, then fail.
        mom_j = sum_ml - mom_i
        if (abs(mom_j) > iMaxLz) then
            cc_ret = -1
            return
        end if

        ! Get the new index
        cc_ret = ClassCountInd (spn_j, sym_j, mom_j)

#ifdef __DEBUG
        if (cc_ret /= -1) then
            if (class_count_ml(cc_ind) /= mom_i) &
                call stop_all(this_routine, 'wrong_mom_i')
            if (class_count_ml(cc_ret) /= mom_j) &
                call stop_all(this_routine, 'wrong_mom_j')
        end if
#endif

    end function


    subroutine gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ex, par, pgen)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NifTot)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        integer, intent(out) :: ex(2,2)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_single_4ind_ex"

        integer :: elec, src, tgt, cc_index
        real(dp) :: pgen_elec
        integer :: loc, dummy_src(2)

        ! In this version of the excitation generator, we pick an electron
        ! at random. Then we construct a list of connection
        ! strengths to each of the available orbitals in the correct symmetry.

        ! We could pick the electron based on the number of orbitals available.
        ! Currently, it is just picked uniformly.
        elec = 1 + floor(genrand_real2_dSFMT() * nel)

        src = nI(elec)

        ! What is the symmetry category?
        cc_index = ClassCountInd (get_spin(src), SpinOrbSymLabel(src), &
                                  G1(src)%Ml)

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

        pgen = pgen / real(nel, dp)

!         print *, "pgen single: ", pgen

    end subroutine


    function pgen_single_4ind (nI, ilutI, src, tgt) result(pgen)

        integer, intent(in) :: nI(nel), src, tgt
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "pgen_single_4ind"

        integer :: cc_index, label_index, norb, n_id(nel), id_src, id_tgt
        integer :: i, j, orb
        real(dp) :: cum_sum, cpt, cpt_tgt
        HElement_t(dp) :: hel

        ! The electron to excite is picked uniformly at random
        pgen = 1.0_dp / real(nel, dp)

        ! The class count index of the target orbital is known
        cc_index = ClassCountInd (get_spin(tgt), SpinOrbSymLabel(tgt), &
                                  G1(tgt)%Ml)
         
        ! How many orbitals of the correct symmetry are there?
        norb = OrbClassCount(cc_index)
        label_index = SymLabelCounts2(1, cc_index)

        ! Some ids for utility
        id_src = gtID(src)
        n_id = gtID(nI)

        ! Generate the cumulative sum, as used in the excitation generator,
        ! and store the relevant term for generating the excitation.
        cum_sum = 0
        do i = 1, norb
            orb = SymLabelList2(label_index + i - 1)
            if (IsNotOcc(ilutI, orb)) then
                hel = 0
                id_tgt = gtID(orb)
                do j = 1, nel
                    if (nI(j) == src) cycle
                    hel = hel + get_umat_el (id_src, n_id(j), id_tgt, &
                                             n_id(j))
                    if (is_beta(src) .eqv. is_beta(nI(j))) then
                        hel = hel - get_umat_el (id_src, n_id(j), n_id(j),&
                                                  id_tgt)
                    end if
                end do
                hel = hel + GetTMATEl(src, orb)
                cpt = abs_l1(hel)

                if (t_matele_cutoff) then 
                    if (cpt < matele_cutoff) then 
                        cpt = 0.0_dp
                    end if
                end if
                cum_sum = cum_sum + cpt
                if (orb == tgt) cpt_tgt = cpt
            end if
        end do

        ! Adjust the generation probability for the relevant values.
        if (cum_sum < EPS) then
            pgen = 0.0_dp
        else
            pgen = pgen * cpt_tgt / cum_sum
        end if


    end function


    function select_orb_sing (nI, ilut, src, cc_index, pgen) result(orb)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: cc_index, src
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'select_orb_sing'

        real(dp) :: cum_sum, cumulative_arr(OrbClassCount(cc_index)), r
        real(dp) :: cpt_arr(OrbClassCount(cc_index))
        integer :: orb, norb, label_index, orb_index, i, j
        integer :: n_id(nel), id_src, id
        HElement_t(dp) :: hel
        real(dp) :: cpt
        
        ! How many orbitals of the correct symmetry are there?
        norb = OrbClassCount(cc_index)
        label_index = SymLabelCounts2(1, cc_index)

        ! Spatial orbital IDs
        n_id = gtID(nI)
        id_src = gtID(src)
        ASSERT(tExch)

        ! Construct the cumulative list of strengths
        cum_sum = 0.0_dp
        do i = 1, norb

            orb = SymLabelList2(label_index + i - 1)
            hel = 0.0_dp
            if (IsNotOcc(ilut, orb)) then
                ASSERT(G1(orb)%Ms == G1(src)%Ms)
                ASSERT(G1(orb)%Ml == G1(src)%Ml)
                ASSERT(SpinOrbSymLabel(orb) == SpinOrbSymLabel(src))

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
            cpt = abs_l1(hel) 

            if (t_matele_cutoff) then
                if (cpt < matele_cutoff) cpt = 0.0_dp
            end if

            cpt_arr(i) = cpt
            cum_sum = cum_sum + cpt_arr(i)
            cumulative_arr(i) = cum_sum

        end do

        !  for testing purposes: 
!         if (cum_sum < UMatEps) then 
!             print *, "===================================="
!             print *, "single excitation: "
!             print *, "norb: ", norb 
!             print *, "cum_sum: ", cum_sum
!             print *, " cumulative_arr: ", cumulative_arr
!             print *, "nI: ", nI 
!             print *, "src: ", src
!         end if

        ! Select a particular orbital to use, or abort.
        ! ok i really think we have to be consistent with this matrix element 
        ! cutoff.. because i think by ignoring some, we allow other excitations 
        ! which should have 0 matrix element to slip through and cause major 
        ! headache..
        if (cum_sum < EPS) then
            orb = 0
            pgen = 0.0_dp
            return
        end if

        if (t_log_ija) then 
            if (cum_sum < ija_thresh) then 
                ija_bins_sing(id_src) = ija_bins_sing(id_src) + 1
                ija_orbs_sing(id_src) = norb
            end if
        end if

        r = genrand_real2_dSFMT() * cum_sum
        orb_index = binary_search_first_ge(cumulative_arr, r)
        orb = SymLabelList2(label_index + orb_index - 1)

        ! And the impact on the generation probability
        pgen = cpt_arr(orb_index) / cum_sum

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
        real(dp) :: cpt_pair(2), sum_pair(2)
        integer :: sym_product, ispn, sum_ml
        integer :: cc_i, cc_j

        ! Initially, select a pair of electrons
        call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)
!        call pick_weighted_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)

        ! Construct the list of possible symmetry pairs listed according to
        ! a unique choice of symmetry A (we enforce that sym(A) < sym(B)).
        cc_tot = 0
        do cc_i = 1, ScratchSize

            cc_j = get_paired_cc_ind (cc_i, sym_product, sum_ml, iSpn)

            ! As cc_i > 0, the following test also excludes any rejected
            ! pairings (where cc_j == -1).
            if (cc_j >= cc_i) then
                cc_tot = cc_tot + store%ClassCountUnocc(cc_i) &
                                * store%ClassCountUnocc(cc_j)
                !cc_tot = cc_tot + OrbClassCount(cc_i) * OrbClassCount(cc_j)
            end if
            cc_cum(cc_i) = cc_tot
        end do

        ! If there are no excitation from this electron pair, then we should
        ! go no further.
        if (cc_tot == 0) then
            nJ(1) = 0
            return
        end if

        ! Pick a specific pairing of spin-symmetries
        r = genrand_real2_dSFMT() * cc_tot
        cc_i = binary_search_first_ge (cc_cum, r)
        cc_j = get_paired_cc_ind (cc_i, sym_product, sum_ml, iSpn)
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
            if (tGen_4ind_lin_exact .or. tGen_4ind_part_exact) then
                call pgen_select_orb(ilutI, src, -1, orbs(2), cpt_pair(1), &
                                     sum_pair(1))
                call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                     cpt_pair(2), sum_pair(2))
            else
                cpt_pair(1) = int_cpt(2)
                cpt_pair(2) = int_cpt(1)
                sum_pair(1) = cum_sum(1)
                sum_pair(2) = cum_sum(1) - int_cpt(2)
            end if
            pgen = pgen * (product(int_cpt) / product(cum_sum) + &
                           product(cpt_pair) / product(sum_pair))
        else
            pgen = pgen * (product(int_cpt) / product(cum_sum))
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

        ! We want to have the n'th alpha, or beta electrons in the determinant
        ! Select them according to the availability of pairs (and the
        ! weighting of opposite-spin pairs relative to same-spin ones).
        r = genrand_real2_dSFMT()
        if (r < pParallel) then
            ! Same spin case
            pgen = pParallel / real(par_elec_pairs, dp)
            r = (r / pParallel) * par_elec_pairs
            idx = floor(r)
            if (idx < AA_elec_pairs) then
                al_req = 2
                be_req = 0
                iSpn = 3
                al_num(1) = ceiling((1 + sqrt(9 + 8*real(idx, dp))) / 2)
                al_num(2) = idx + 1 - ((al_num(1) - 1) * (al_num(1) - 2)) / 2
            else
                al_req = 0
                be_req = 2
                iSpn = 1
                idx = idx - AA_elec_pairs
                be_num(1) = ceiling((1 + sqrt(9 + 8*real(idx, dp))) / 2)
                be_num(2) = idx + 1 - ((be_num(1) - 1) * (be_num(1) - 2)) / 2
            end if
        else
            ! Opposite spin case
            iSpn = 2
            al_req = 1
            be_req = 1
            pgen = (1.0_dp - pParallel) / real(AB_elec_pairs, dp)
            r = ((r - pParallel) / (1.0_dp - pParallel)) * AB_elec_pairs
            idx = floor(r)
            al_num(1) = 1 + mod(idx, nOccAlpha)
            be_num(1) = 1 + floor(idx / real(nOccAlpha, dp))
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
        sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                         SpinOrbSymLabel(src(2)))

    end subroutine


    subroutine pgen_select_orb (ilut, src, orb_pair, tgt, cpt, cum_sum)

        integer, intent(in) :: src(2), orb_pair, tgt
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp) :: cpt, cum_sum

        integer :: norb, label_index, orb, srcid(2), i, src_id
        integer :: cc_index, ms
        real(dp) :: tmp

        ! How many orbitals are available with the given symmetry?
        cc_index = ClassCountInd(get_spin(tgt), SpinOrbSymLabel(tgt), &
                                 G1(tgt)%Ml)
        label_index = SymLabelCounts2(1, cc_index)
        norb = OrbClassCount(cc_index)

        ! We perform different sums depending on the relative spins :-(
        cum_sum = 0.0_dp
        if (is_beta(src(1)) .eqv. is_beta(src(2))) then
            ! Both electrons have the same spin. So we need to include both
            ! electron-hole interactions.
            srcid = gtID(src)
            do i = 1, norb
                orb = SymLabelList2(label_index + i - 1)
                if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                    tmp = same_spin_pair_contrib(srcid(1), srcid(2), &
                                                 orb, orb_pair)
                    cum_sum = cum_sum + tmp
                    if (orb == tgt) cpt = tmp
                end if
            end do
        else
            ! The two electrons have differing spin. Therefore, only the
            ! electron-hole interaction with the same spin is required
            ms = class_count_ms(cc_index)
            if (ms == G1(src(1))%Ms) then
                srcid = gtID(src)
            else
                srcid(1) = gtID(src(2))
                srcid(2) = gtID(src(1))
            end if

            do i = 1, norb
                orb = SymLabelList2(label_index + i - 1)
                if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                    tmp = opp_spin_pair_contrib(srcid(1), srcid(2), &
                                                orb, orb_pair)
                    cum_sum = cum_sum + tmp
                    if (orb == tgt) cpt = tmp
                end if
            end do
        end if

    end subroutine


    function opp_spin_pair_contrib(indi, indj, orba, orbb) result(contrib)

        ! What term should we be summing into the list of contributions?
        !
        ! n.b. Because this routine is only passed the spatial index of i,j
        !      it does _no_ checking that the appropriate spin terms have been
        !      selected

        integer, intent(in) :: indi, indj, orba, orbb
        real(dp) :: contrib
        integer :: inda, indb

        if (tGen_4ind_part_exact .and. orbb > 0) then
            ! Include a contribution of: sqrt(abs(<ij|ab>))
            ! n.b. This can only be used for the case <ij|ba> == 0.
            inda = gtID(orba)
            indb = gtID(orbb)
            if (tGen_4ind_unbound) then
                contrib = abs(get_umat_el(indi, indj, inda, indb))
            else
                contrib = max(sqrt(abs(get_umat_el(indi, indj, inda, indb))), 0.0001_dp)
!                 contrib = sqrt(abs(get_umat_el(indi, indj, inda, indb)))
            end if
        else if (tGen_4ind_lin_exact) then
            if (orbb > 0) then
                ! Include a contribution of abs(<ij|ab>)
                ! n.b. This can only be used for the case <ij|ba> == 0.
                inda = gtID(orba)
                indb = gtID(orbb)
                contrib = abs(get_umat_el(indi, indj, inda, indb))
            else
                ! Select first orbital linearly
                contrib = 1.0_dp
            end if
        else
            ! Include the contribution of this term sqrt(<ia|ia>)
            inda = gtID(orba)
!             
!             if (t_iiaa .and. t_ratio) then 
!                 ! ok.. maybe i have to talk to ali about that, what he 
!                 ! meant with this splitting of p(a|ij) = p(j)*p(a|i) 
!                 ! because i am not sure about that ..
!                 ! althoug i should be carefull if we do not divide by 
!                 ! 0 here.. 
!                 ! NOTE: by testing it was seen that the ratio approach causes the 
!                 ! pgens to be much too low and thus the H_ij/pgen ratios to 
!                 ! explode
! 
!                 contrib = sqrt(abs(get_umat_el(indi, inda, indi, inda) / & 
!                            max(abs(get_umat_el(indj, inda, indj, inda)), 0.0001_dp))) &
!                         + sqrt(abs(get_umat_el(indj, inda, indi, inda) / & 
!                            max(abs(get_umat_el(indi, inda, indj, inda)), 0.0001_dp)))

            if (t_iiaa) then 
                
                contrib = sqrt(abs(get_umat_el(indi, inda, indi, inda)))

!             else if (t_ratio) then 
!                 ! also here i have to check if i actually should take care 
!                 ! of the indj influence.. 
! 
!                 contrib = sqrt(abs(UMat2D(max(indi, inda), min(indi, inda))) / & 
!                            max(abs(UMat2D(max(indj, inda), min(indj, inda))), 0.0001_dp)) &
!                         + sqrt(abs(UMat2D(max(indj, inda), min(indj, inda))) / & 
!                            max(abs(UMat2D(max(indi, inda), min(indi, inda))), 0.0001_dp))

            else 

                contrib = sqrt(abs_l1(UMat2D(max(indi, inda), min(indi, inda))))
                
            end if

        end if

        if (t_matele_cutoff) then
            if (contrib < matele_cutoff) contrib = 0.0_dp
        end if

    end function


    function same_spin_pair_contrib(indi, indj, orba, orbb) result(contrib)

        ! What term should be be summing into the list of the contributions
        ! from the ij term

        integer, intent(in) :: indi, indj, orba, orbb
        real(dp) :: contrib
        integer :: inda, indb

        if (tGen_4ind_part_exact .and. orbb > 0) then
            ! Include a contribution of:
            ! sqrt(abs(<ij|ab> - <ij|ba>))
            inda = gtID(orba)
            indb = gtID(orbb)
            if (tGen_4ind_unbound) then
                contrib = abs(get_umat_el(indi, indj, inda, indb) &
                                - get_umat_el(indi, indj, indb, inda))
            else
                ! finally get rid of this arbitrary thresholds..
                contrib = max(sqrt(abs(get_umat_el(indi, indj, inda, indb) &
                                - get_umat_el(indi, indj, indb, inda))), 0.00001_dp)
!                 contrib = sqrt(abs(get_umat_el(indi, indj, inda, indb) &
!                                 - get_umat_el(indi, indj, indb, inda)))
            end if
        else if (tGen_4ind_lin_exact) then
            if (orbb > 0) then
                ! Include a contribution of:
                ! abs(<ij|ab> - <ij|ba>)
                inda = gtID(orba)
                indb = gtID(orbb)
                contrib = abs(get_umat_el(indi, indj, inda, indb) &
                            - get_umat_el(indi, indj, indb, inda))
            else
                ! Select first orbital linearly.
                contrib = 1.0_dp
            end if
        else
            ! Include a contribution of (orb can be a or b):
            ! sqrt((ii|aa) + (jj|aa))
            inda = gtID(orba)

            if (t_iiaa .and. t_ratio) then 
                ! ok.. maybe i have to talk to ali about that, what he 
                ! meant with this splitting of p(a|ij) = p(j)*p(a|i) 
                ! because i am not sure about that ..
                ! althoug i should be carefull if we do not divide by 
                ! 0 here.. 

                contrib = sqrt(abs(get_umat_el(indi, inda, indi, inda) / & 
                           max(abs(get_umat_el(indj, inda, indj, inda)), 0.0001_dp))) &
                        + sqrt(abs(get_umat_el(indj, inda, indj, inda) / & 
                           max(abs(get_umat_el(indi, inda, indi, inda)), 0.0001_dp)))

            else if (t_iiaa) then 
                
                contrib = sqrt(abs(get_umat_el(indi, inda, indi, inda))) & 
                        + sqrt(abs(get_umat_el(indj, inda, indj, inda)))

            else if (t_ratio) then 
                ! also here i have to check if i actually should take care 
                ! of the indj influence.. 

                contrib = sqrt(abs(UMat2D(max(indi, inda), min(indi, inda))) / & 
                           max(abs(UMat2D(max(indj, inda), min(indj, inda))), 0.0001_dp)) &
                        + sqrt(abs(UMat2D(max(indj, inda), min(indj, inda))) / & 
                           max(abs(UMat2D(max(indi, inda), min(indi, inda))), 0.0001_dp))

            else 


                contrib = sqrt(abs_l1(UMat2D(max(indi, inda), min(indi, inda))))&
                        + sqrt(abs_l1(UMat2D(max(indj, inda), min(indj, inda))))
            end if
            !sqrt(abs_l1(get_umat_el(srcid(1), srcid(1), inda, inda))) + &
            !sqrt(abs_l1(get_umat_el(srcid(2), srcid(2), inda, inda)))
        end if

        if (t_matele_cutoff) then
            if (contrib < matele_cutoff) contrib = 0.0_dp
        end if

    end function


    function select_orb (ilut, src, cc_index, orb_pair, cpt, cum_sum) &
            result(orb)

        ! For an excitation from electrons 1,2, consider all of the pairs
        ! (e1 a|e1 a) == <e1 e1 | a a> and (e2 a | e2 a) == <e2 e2 | a a>
        ! as appropriate to bias selection of the electron

        integer, intent(in) :: src(2), orb_pair, cc_index
        real(dp), intent(out) :: cpt, cum_sum
        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer :: orb
        character(*), parameter :: this_routine = "select_orb"

        integer :: label_index, orb_index, norb, i, srcid(2), ms
        integer :: src_id

        ! Our biasing arrays must consider all of the possible orbitals with
        ! the correct symmetry.
        real(dp) :: cumulative_arr(OrbClassCount(cc_index)), r

        logical :: t_par
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

            t_par = .true.
            ! Both of the electrons have the same spin. Therefore we need to
            ! include both electron-hole interactions.
            cum_sum = 0.0_dp
            srcid = gtID(src)
            do i = 1, norb
            
                orb = SymLabelList2(label_index + i - 1)
                if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                    cum_sum = cum_sum + same_spin_pair_contrib(&
                                             srcid(1), srcid(2), orb, orb_pair)
                end if
                cumulative_arr(i) = cum_sum

            end do

        else
            
            t_par = .false.
            ! The two electrons have differing spin. Therefore, only the 
            ! electron-hole interaction with the same spin is required.
            ms = class_count_ms(cc_index)
            if (ms == G1(src(1))%Ms) then
                srcid = gtID(src)
            else
                srcid(1) = gtID(src(2))
                srcid(2) = gtID(src(1))
            end if

            cum_sum = 0.0_dp
            do i = 1, norb

                orb = SymLabelList2(label_index + i - 1)
                if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                    cum_sum = cum_sum + opp_spin_pair_contrib(&
                                            srcid(1), srcid(2), orb, orb_pair)
                end if
                cumulative_arr(i) = cum_sum

            end do

        end if

        ! also check here ig the problem is overall low pgens for certain 
        ! excitations 


        ! If there are no available orbitals to pair with, we need to abort
        if (cum_sum < EPS) then
            orb = 0
            return
        end if

        ! do i need the matele cutoff here too? i shouldnt.. 
        ! check in debug mode!
#ifdef __DEBUG
        if (t_matele_cutoff) then 
            if (cum_sum < matele_cutoff) then 
                call stop_all(this_routine, &
                    "although matrix-cutoff something slipped through..")
            end if
        end if
#endif

        ! the dead end we want to log are here actually.. 
        if (t_log_ija .and. cum_sum < ija_thresh) then 
            ! here source is not yet sorted! 
            ! but only take the unique (ij) combinations!
            ! and do i want to have more information? 
            ! maybe i want to know how many symmetry allowed orbitals there 
            ! are for this kind of excitation... yes!
            ! and maybe i only want to store the spatial orbitals and 
            ! the info if it is a parallel spin excitation or an opposite 
            ! spin excitation.. this would reduce the output amount even 
            ! farther yes! 
            if (t_par) then 
                ija_bins_para(minval(srcid),maxval(srcid),gtID(orb_pair)) = &
                    ija_bins_para(minval(srcid),maxval(srcid),gtID(orb_pair)) + 1

                ija_orbs_para(minval(srcid),maxval(srcid),gtID(orb_pair)) = norb
            else
                ija_bins_anti(minval(srcid),maxval(srcid),gtID(orb_pair)) = &
                    ija_bins_anti(minval(srcid),maxval(srcid),gtID(orb_pair)) + 1

                ija_orbs_anti(minval(srcid),maxval(srcid),gtID(orb_pair)) = norb

            end if
        end if
        
!         if (cum_sum < 1.0e-4_dp) then 
!             print *, "========================================"
!             if (t_par) then 
!                 print *, "parallel double excitation: "
!             else 
!                 print *, "opposite double excitation: "
!             end if
!             print *, "norb: ", norb
!             print *, "cum_sum: ", cum_sum
!             print *, "cumulative_arr: ", cumulative_arr
!             print *, "(i,j): ", src
!             print *, "(a): ", orb_pair
!         end if

        ! Binary search within this list to choose a value.
        r = genrand_real2_dSFMT() * cum_sum
        orb_index = binary_search_first_ge(cumulative_arr(1:norb), r)

        ! And return the relevant value.
        orb = SymLabelList2(label_index + orb_index - 1)
        if (orb_index == 1) then
            cpt = cumulative_arr(1)
        else
            cpt = cumulative_arr(orb_index) - cumulative_arr(orb_index - 1)
        end if

    end function


    subroutine gen_excit_4ind_reverse (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                       ExcitMat, tParity, pGen, HelGen, store, part_type)

        ! TODO: description
        !
        ! n.b. (ij|kl) <= sqrt( (ij|ij) * (kl|kl) )
        !      This provides quite a good description of the large elements

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(in), optional :: part_type

        character(*), parameter :: this_routine = 'gen_excit_4ind_reverse'

        integer :: orb
        real(dp) :: pgen2

        HElGen = HEl_zero
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

        ! And a careful check!
#ifdef __DEBUG
        if (.not. IsNullDet(nJ)) then
            pgen2 = calc_pgen_4ind_reverse (nI, ilutI, ExcitMat, ic)
            if (abs(pgen - pgen2) > 1.0e-6_dp) then
                write(6,*) 'Calculated and actual pgens differ.'
                write(6,*) 'This will break HPHF calculations'
                call write_det(6, nI, .false.)
                write(6, '(" --> ")', advance='no')
                call write_det(6, nJ, .true.)
                write(6,*) 'Excitation matrix: ', ExcitMat(1,1:ic), '-->', &
                           ExcitMat(2,1:ic)
                write(6,*) 'Generated pGen:  ', pgen
                write(6,*) 'Calculated pGen: ', pgen2
                call stop_all(this_routine, "Invalid pGen")
            end if
        end if
#endif

    end subroutine


    function calc_pgen_4ind_reverse (nI, ilutI, ex, ic) result(pgen)

        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer(n_int), intent(in) :: ilutI(0:NifTot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'calc_pgen_4ind_reverse'

        integer :: iSpn, sym_product, i, j, src(2), e_ispn, e_sym_prod, id(2)
        integer :: tgt(2), id_tgt(2), sum_ml, e_sum_ml
        real(dp) :: ntot, cum_sum, cpt, cpt_tgt
        HElement_t(dp) :: hel

        if (ic == 1) then

            ! Singles are calculated in the same way for the _ex and the
            ! _reverse excitation generators
            pgen = pSingles * pgen_single_4ind (nI, ilutI, ex(1,1), ex(2,1))

        else if (ic == 2) then

            ! Obviously bias by pDoubles...
            pgen = pDoubles 

            ! We want to select a pair of orbitals in a way which is biased
            ! towards pairs with opposing spins. This takes a simple form.
            if (is_alpha(ex(2,1)) .eqv. is_alpha(ex(2,2))) then
                pgen = pgen * pParallel / real(par_hole_pairs, dp)
                if (is_alpha(ex(2,1))) then
                    iSpn = 3
                else
                    iSpn = 1
                end if
            else
                pgen = pgen * (1.0_dp - pParallel) / real(AB_hole_pairs, dp)
                iSpn = 2
            end if

            ! Now consider all of the possible pairs of electrons
            tgt = ex(2,1:2)
            id_tgt = gtID(tgt)
            sym_product = RandExcitSymLabelProd(SpinOrbSymLabel(tgt(1)), &
                                                SpinOrbSymLabel(tgt(2)))
            sum_ml = sum(G1(tgt)%Ml)
            cum_sum = 0
            do i = 2, nel
                src(1) = nI(i)
                do j = 1, i-1

                    ! Get the symmetries
                    src(2) = nI(j)
                    e_ispn = get_ispn(src(1), src(2))
                    e_sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                                       SpinOrbSymLabel(src(2)))
                    e_sum_ml = sum(G1(src)%Ml)
                    
                    ! Calculate the cross HElements
                    if (e_ispn == iSpn .and. e_sym_prod == sym_product &
                        .and. e_sum_ml == sum_ml) then
                        hel = 0
                        id = gtID(src)
                        if ((is_beta(src(1)) .eqv. is_beta(tgt(1))) .and. &
                            (is_beta(src(2)) .eqv. is_beta(tgt(2)))) &
                            hel = hel + get_umat_el (id(1), id(2), id_tgt(1), &
                                                     id_tgt(2))
                        if ((is_beta(src(1)) .eqv. is_beta(tgt(2))) .and. &
                            (is_beta(src(2)) .eqv. is_beta(tgt(1)))) &
                            hel = hel - get_umat_el (id(1), id(2), id_tgt(2), &
                                                     id_tgt(1))
                        cpt = abs_l1(hel)
                    else
                        cpt = 0
                    end if
                    cum_sum = cum_sum + cpt

                    ! And if this is the excitation, store the value.
                    if (minval(src) == minval(ex(1,1:2)) .and. &
                        maxval(src) == maxval(ex(1,1:2))) then
                        cpt_tgt = cpt
                    end if

                end do
            end do

            ! And account for the case where this is not a connected excitation
            ! actually this comparison with 0 should be removed..
            if (cum_sum == 0) then
!             if (cum_sum < EPS) then
                pgen = 0
            else
                pgen = pgen * cpt_tgt / cum_sum
            end if

        else
            ! IC /= 1, 2 --> not connected by the excitation generator.
            ! 
            ! If we are hitting here, then earlier checks have failed. We can
            ! return the correct (zero) value at runtime, but really this
            ! should be fixed elsewhere
            ASSERT(.false.)
            pgen = 0.0_dp
        end if

    end function



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
        character(*), parameter :: this_routine = 'select_elec_sing'

        real(dp) :: cum_sum, cum_arr(nel), cpt_arr(nel), r
        integer :: n_id(nel), id_tgt, id_src, cc_src, j, src_elec
        integer :: tgt_sym, tgt_ml
        logical :: tgt_beta
        HElement_t(dp) :: hel

        ! Don't consider symmetry categories, do the components separately.
        ! --> Slightly quicker
        tgt_beta = is_beta(tgt)
        tgt_sym = SpinOrbSymLabel(tgt)
        tgt_ml = G1(tgt)%Ml

        ! Spatial orbital IDs
        n_id = gtID(nI)
        id_tgt = gtID(tgt)
        ASSERT(tExch)

        ! Construct the cumulative list of strengths
        cum_sum = 0
        do src_elec = 1, nel

            ! What is the symmetry/spin class of the source electron
            src = nI(src_elec)

            ! If the symmetry/spin are not the same, then the matrix element
            ! will be zero --> don't calculate it!
            hel = 0
            if ((is_beta(src) .eqv. tgt_beta) .and. &
                SpinOrbSymLabel(src) == tgt_sym .and. &
                G1(src)%Ml == tgt_ml) then

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
!         if (cum_sum < EPS) then
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

        integer :: src(2), tgt(2), elecs(2), id_tgt(2), id(2), e_sum_ml
        integer :: sym_prod, ispn, e_ispn, e_sym_prod, idx, i, j, sum_ml
        HElement_t(dp) :: hel
        real(dp) :: cum_arr(nel * (nel - 1) / 2)
        real(dp) :: val_arr(nel * (nel - 1) / 2)
        real(dp) :: val, cum_val, r

        ! Select a pair of holes (vacant orbitals)
        call pick_hole_pair_biased (ilutI, tgt, pgen, sym_prod, sum_ml, ispn)
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
                e_sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                                   SpinOrbSymLabel(src(2)))
                e_sum_ml = sum(G1(src)%Ml)

                ! Get the weight (HElement) associated with the elecs/holes
                if (e_ispn == iSpn .and. e_sym_prod == sym_prod .and. &
                    e_sum_ml == sum_ml) then
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


    subroutine pick_hole_pair_biased (ilut, orbs, pgen, sym_prod, sum_ml, ispn)
        
        ! This is a biased version of pick_hole_pair below. Whereas the
        ! other function picks hole pairs entirely uniformly throughout the
        ! space, this function biases that selection towards picking opposite
        ! spin pairs (according to the biasing factor rand_excit_par_bias).

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: orbs(2), sym_prod, ispn, sum_ml
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'pick_hole_pair_biased'

        real(dp) :: ntot, r
        integer :: spn(2)

        ! Pick what sort of hole pair we want. AA, BB, or AB
        r = genrand_real2_dSFMT()
        if (r < pParallel) then
            r = (r / pParallel) * par_hole_pairs
            pgen = pParallel / real(par_hole_pairs, dp)
            if (r < AA_hole_pairs) then
                ! alpha/alpha
                iSpn = 3
                spn = (/1, 1/)
            else
                ! beta/beta
                iSpn = 1
                spn = (/2, 2/)
            end if

        else
            ! Opposite spin
            iSpn = 2
            spn = (/1, 2/)
            pgen = (1.0_dp - pParallel) / real(AB_hole_pairs, dp)
        end if

        ! Select orbitals at random with the given spins
        orbs(1) = pick_hole_spn (ilut, spn(1), -1)
        orbs(2) = pick_hole_spn (ilut, spn(2), orbs(1))

        ! What is the cumulative Ml value?
        sum_ml = sum(G1(orbs)%Ml)

        ! What are the symmetry/spin properties of this pick?
        sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(orbs(1)), &
                                         SpinOrbSymLabel(orbs(2)))
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
        character(*), parameter :: this_routine = 'pick_hole'

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
                call stop_all(this_routine, 'Out of attempts')
            end if
        end do

    end function

    
    subroutine pick_hole_pair (ilut, orbs, pgen, sym_prod, sum_ml, ispn)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: orbs(2), sym_prod, ispn, sum_ml
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
        sym_prod = RandExcitSymLabelProd (SpinOrbSymLabel(orbs(1)), &
                                          SpinOrbSymLabel(orbs(2)))
        ispn = get_ispn(orbs(1), orbs(2))
        sum_ml = sum(G1(orbs)%Ml)

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
        integer(n_int), pointer :: det_list(:,:)
        real(dp), allocatable :: contrib_list(:), pgen_list(:), matEle_list(:)
        logical, allocatable :: generated_list(:)
        logical :: found_all, par
        real(dp) :: contrib, pgen, sum_pgens, sum_helement
        HElement_t(dp) :: helgen, hel
        character(255) :: filename

        ! Decode the determiant
        call decode_bit_det (src_det, ilut)

        ! could just use my new calc_all_excitations
        call calc_all_excitations(ilut, det_list, nexcit)

        ! Initialise
        call init_excit_gen_store (store)
        store%tFilled = .false.

!         ! How many connected determinants are we expecting?
!         call CountExcitations3 (src_det, 3, nsing, ndoub)
!         nexcit = nsing + ndoub
!         allocate(det_list(0:NIfTot, nexcit))
! 
!         ! Loop through all of the possible excitations
!         ndet = 0
!         found_all = .false.
!         ex = 0
!         flag = 3
        write(6,'("*****************************************")')
        write(6,'("Enumerating excitations")')
        write(6,'("Starting from: ")', advance='no')
        call write_det (6, src_det, .true.)
        write(6,*) 'Expecting ', nexcit, "excitations"
!         call GenExcitations3 (src_det, ilut, det, flag, ex, par, found_all, &
!                               .false.)
!         do while (.not. found_all)
!             ndet = ndet + 1
!             call EncodeBitDet (det, det_list(:,ndet))
! 
!             call GenExcitations3 (src_det, ilut, det, flag, ex, par, &
!                                   found_all, .false.)
!         end do
!         if (ndet /= nexcit) &
!             call stop_all(this_routine,"Incorrect number of excitations found")
! 
!         ! Sort the dets, so they are easy to find by binary searching
!         call sort(det_list, ilut_lt, ilut_gt)

        ! Lists to keep track of things
        allocate(generated_list(nexcit))
        allocate(contrib_list(nexcit))
        allocate(pgen_list(nexcit))
        allocate(matEle_list(nexcit))
        generated_list = .false.
        contrib_list = 0.0_dp
        pgen_list = 0.0_dp
        matEle_list = 0.0_dp

        ! Repeated generation, and summing-in loop
        ngen = 0
        contrib = 0.0_dp
        do i = 1, iterations
            if (mod(i, 10000) == 0) &
                write(6,*) i, '/', iterations, ' - ', contrib / (real(nexcit,dp)*i)

            !pSingles = 0.0
            !pDoubles = 1.0

            call gen_excit_4ind_weighted (src_det, ilut, det, tgt_ilut, 3, &
                                          ic, ex, par, pgen, helgen, store)
            if (det(1) == 0) cycle

            call EncodeBitDet (det, tgt_ilut)

            helgen = get_helement(src_det, det, ic, ex, par) 

            if (abs(helgen) < 1.0e-6_dp) cycle

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
                matEle_list(pos) = helgen
                pgen_list(pos) = pgen
            end if
        end do

        ! How many of the iterations generated a good det?
        write(6,*) ngen, " dets generated in ", iterations, " iterations."
        write(6,*) 100_dp * (iterations - ngen) / real(iterations,dp), &
                   '% abortion rate'
        ! Contribution averages
        write(6, '("Averaged contribution: ", f15.10)') &
                contrib / (real(nexcit,dp) * iterations)

        ! Output the determinant specific contributions
        iunit = get_free_unit()
        open(iunit, file="contribs_4ind", status='unknown')
        do i = 1, nexcit
            call writebitdet(iunit, det_list(:,i), .false.)
            write(iunit, *) contrib_list(i) / real(iterations, dp)
        end do
        close(iunit)

        ! Check that all of the determinants were generated!!!
        if (.not. all(generated_list)) then
            write(6,*) count(.not.generated_list), '/', size(generated_list), &
                       'not generated'
            found_all = .true.
            do i = 1, nexcit
                if (.not. generated_list(i)) then
                    call decode_bit_det(det, det_list(:,i))
                    hel = get_helement(src_det, det, ilut, det_list(:,i))
                    if (abs(hel) > 1.0e-6) then
                        found_all = .false.
                        call writebitdet(6, det_list(:,i), .false.)
                        write(6,*) hel
                    end if
                end if
            end do
            if (.not. found_all) &
                call stop_all(this_routine, "Determinant not generated")
        end if
        if (any(abs(contrib_list / iterations - 1.0_dp) > 0.01_dp)) then
            do i = 1, nexcit
                call writebitdet(6, det_list(:,i), .false.)
                write(6,*) contrib_list(i) / (iterations - 1.0_dp)
            end do
!             call stop_all(this_routine, "Insufficiently uniform generation")
        end if

        sum_pgens = sum(pgen_list)
        sum_helement = sum(abs(matEle_list))

        iunit = get_free_unit()
        call get_unique_filename("pgen_vs_matrixElements",.true.,.true.,1,filename)
        open(iunit, file=filename,status='unknown')
        write(iunit,*) "pgens and matrix elements for CSF:"
        call write_det(6, src_det, .true.)
        do i = 1, nExcit
            write(iunit, "(f16.7)", advance = 'no') pgen_list(i) !/sum_pgens
            write(iunit, "(f16.7)", advance = 'no') matEle_list(i) !/sum_helement
            write(iunit, "(f16.7)") contrib_list(i) / real(iterations,dp)
        end do
        close(iunit)

        ! Clean up
        deallocate(det_list)
        deallocate(contrib_list)
        deallocate(generated_list)
        deallocate(pgen_list)
        deallocate(matEle_list)

    end subroutine


    subroutine pick_weighted_elecs(nI, elecs, src, sym_prod, ispn, sum_ml, &
                                   pgen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), src(2), sym_prod, ispn, sum_ml
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = 'pick_weighted_elecs'

        logical, parameter :: tlinear = .true.
        logical, parameter :: tlinear2 = .false.

        real(dp) :: cum_list(nel), cum_sum, cpt, final_cpt, r
        real(dp) :: cum_sum_ii, cum_list_ii(nel), cpt_ii, cpt_jj
        integer :: i, ind1, ind2, inds(nel)

        ! Pick the first electron uniformly, or according to <ii|ii>,
        ! and then pick the second electron according to <ji|ji>

        inds = gtID(nI)
        if (tlinear) then
            ! Pick the first electron uniformly
            elecs(1) = 1 + floor(genrand_real2_dSFMT() * nel)
            cpt_ii = 1.0_dp
            cum_sum_ii = nel
        else
            ! Pick the first electron according to <ii|ii>
            cum_sum_ii = 0
            do i = 1, nel
                cum_sum_ii = cum_sum_ii &
                           + abs(get_umat_el(inds(i), inds(i), &
                                             inds(i), inds(i)))
                cum_list_ii(i) = cum_sum_ii
            end do

            r = genrand_real2_dSFMT() * cum_sum_ii
            elecs(1) = binary_search_first_ge(cum_list_ii, r)
            if (elecs(1) == 1) then
                cpt_ii = cum_list_ii(1)
            else
                cpt_ii = cum_list_ii(elecs(1)) - cum_list_ii(elecs(1) - 1)
            end if
        end if

        ! Construct the weighted list <ji|ji>
        cum_sum = 0
        ind1 = inds(elecs(1))
        do i = 1, nel
            if (i == elecs(1)) then
                cpt = 0
            else
                if (tlinear2) then
                    cpt = 1.0
                else
                    cpt = abs(get_umat_el(ind1, inds(i), ind1, inds(i)))
                end if
                if (is_beta(nI(i)) .eqv. is_beta(nI(elecs(1)))) then
                    cpt = cpt * pParallel
                else
                    cpt = cpt * (1.0_dp - pParallel)
                end if
            end if
            cum_sum = cum_sum + cpt
            cum_list(i) = cum_sum
        end do

        ! Get the appropriate index
        r = genrand_real2_dSFMT() * cum_sum
        elecs(2) = binary_search_first_ge(cum_list, r)

        ! Calculate the (partial) probability
        if (elecs(2) == 1) then
            cpt = cum_list(1)
        else
            cpt = cum_list(elecs(2)) - cum_list(elecs(2) - 1)
        end if
        pgen = (cpt_ii / cum_sum_ii) * cpt / cum_sum

        ! To calculate the generation probability, we need to consider the
        ! possibility of having picked the electrons in the order j, i.
        ! Construct the reverse list
        cum_sum = 0
        ind2 = inds(elecs(2))
        do i = 1, nel
            if (i == elecs(2)) then
                cpt = 0
            else
                if (tlinear2) then
                    cpt = 1.0
                else
                    cpt = abs(get_umat_el(ind2, inds(i), ind2, inds(i)))
                end if
                if (is_beta(nI(i)) .eqv. is_beta(nI(elecs(2)))) then
                    cpt = cpt * pParallel
                else
                    cpt = cpt * (1.0_dp - pParallel)
                end if
            endif
            if (i == elecs(1)) final_cpt = cpt
            cum_sum = cum_sum + cpt
        end do

        ! And get the component for the second electron in the initial list
        if (tlinear) then
            cpt_jj = 1.0_dp
        else
            if (elecs(2) == 1) then
                cpt_jj = cum_list_ii(1)
            else
                cpt_jj = cum_list_ii(elecs(2)) - cum_list_ii(elecs(2) - 1)
            end if
        end if

        ! Adjust the probability for the j,i choice and then the 1/N one
        pgen = pgen + ((cpt_jj / cum_sum_ii) * (final_cpt / cum_sum))

        ! Generate the orbitals under consideration
        src = nI(elecs)

        if (is_beta(src(1)) .eqv. is_beta(src(2))) then
            if (is_beta(src(1))) then
                iSpn = 1
            else
                iSpn = 3
            end if
        else
            iSpn = 2
        end if

        ! The Ml value is obtained from the orbitals
        sum_ml = sum(G1(src)%Ml)

        ! And the spatial symmetries
        sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                         SpinOrbSymLabel(src(2)))

    end subroutine pick_weighted_elecs


    function pgen_weighted_elecs(nI, orbs) result(prob)

        ! Return the probability of having picked the electrons, and the

        integer, intent(in) :: nI(nel), orbs(2)
        real(dp) :: prob
        logical, parameter :: tlinear = .true.
        logical, parameter :: tlinear2 = .false.

        real(dp) :: cpt1(2), cpt2(2), cum_sum1, cum_sum2(2), cpt, cpts(2)
        integer :: i, inds(nel), ind1(2)

        ! Enumerate possibilities for the first electron, and select the
        ! relevant orbitals.
        inds = gtID(nI)
        ind1 = gtID(orbs)
        if (tlinear) then
            cum_sum1 = 0
            cum_sum2 = 0
            do i = 1, nel

                ! Contributions for selection of electron 1
                if (.not. tlinear) then
                    cpt = abs(get_umat_el(inds(i), inds(i), inds(i), inds(i)))
                    cum_sum1 = cum_sum1 + cpt
                end if

                ! Extract the necessary components for calculating first e-
                ! probability, and grab list components for the second one.
                ! Contributions for selection of electron 2
                if (nI(i) == orbs(1)) then
                    cpt1(1) = cpt
                    cpts(1) = 0.0
                else
                    if (tlinear2) then
                        cpts(1) = 1.0
                    else
                        cpts(1) = abs(get_umat_el(ind1(1), inds(i), ind1(1), &
                                      inds(i)))
                    end if
                    if (is_beta(nI(i)) .eqv. is_beta(orbs(1))) then
                        cpts(1) = cpts(1) * pParallel
                    else
                        cpts(1) = cpts(1) * (1.0_dp - pParallel)
                    end if
                end if
                if (nI(i) == orbs(2)) then
                    cpt1(2) = cpt
                    cpts(2) = 0.0
                else
                    if (tlinear2) then
                        cpts(2) = 1.0
                    else
                        cpts(2) = abs(get_umat_el(ind1(2), inds(i), ind1(2), &
                                      inds(i)))
                    end if
                    if (is_beta(nI(i)) .eqv. is_beta(orbs(2))) then
                        cpts(2) = cpts(2) * pParallel
                    else
                        cpts(2) = cpts(2) * (1.0_dp - pParallel)
                    end if
                end if

                ! And extract the correct second electron components
                ! n.b. this choice is for the second electron, hence 2,1
                if (nI(i) == orbs(2)) cpt2(1) = cpts(1)
                if (nI(i) == orbs(1)) cpt2(2) = cpts(2)
                cum_sum2 = cum_sum2 + cpts
            end do
        end if

        ! If we are selecting the first electron in a linear manner, then
        ! we know the values exactly.
        if (tlinear) then
            cpt1 = 1.0
            cum_sum1 = nel
        end if

        ! And extract the combined probability
        prob = ((cpt1(1) / cum_sum1) * (cpt2(1) / cum_sum2(1))) &
             + ((cpt1(2) / cum_sum1) * (cpt2(2) / cum_sum2(2)))

    end function

    subroutine calc_all_excitations(ilut, non_zero_list, n_non_zero)
        ! routine which calculates all the excitations for a given 
        ! determinant.. have to move that to a different file later on.. 
        ! does not fit in here!
        use bit_reps, only: decode_bit_det
        use GenRandSymExcitNUMod, only: init_excit_gen_store
        use SymExcit3, only: CountExcitations3, GenExcitations3
        use Determinants, only: write_det
        use DetBitOps, only: EncodeBitDet
        use sort_mod, only: sort
        
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer(n_int), intent(out), pointer :: non_zero_list(:,:)
        integer, intent(out) :: n_non_zero
        character(*), parameter :: this_routine = "calc_all_excitations"

        integer :: src_det(nel), nsing, ndoub, ndet, ex(2,2), flag, det(nel), &
                   nexcit, i, cnt
        type(excit_gen_store_type) :: store
        logical :: found_all, par
        real(dp), allocatable :: hel_list(:)
        integer(n_int), pointer :: det_list(:,:)
        logical, allocatable :: t_non_zero(:)

        ! Decode the determiant
        call decode_bit_det (src_det, ilut)

        ! Initialise
        call init_excit_gen_store (store)

        ! How many connected determinants are we expecting?
        call CountExcitations3 (src_det, 3, nsing, ndoub)
        nexcit = nsing + ndoub
        allocate(det_list(0:NIfTot, nexcit))
        allocate(hel_list(nexcit))
        hel_list = 0.0_dp

        ! Loop through all of the possible excitations
        ndet = 0
        found_all = .false.
        ex = 0
        flag = 3
!         write(6,'("*****************************************")')
!         write(6,'("Enumerating excitations")')
!         write(6,'("Starting from: ")', advance='no')
!         call write_det (6, src_det, .true.)
!         write(6,*) 'Expecting ', nexcit, "excitations"
        call GenExcitations3 (src_det, ilut, det, flag, ex, par, found_all, &
                              .false.)

        if (tGUGA) then 
            call stop_all(this_routine, "modify get_helement for GUGA")
        end if

        do while (.not. found_all)
            ndet = ndet + 1
            call EncodeBitDet (det, det_list(:,ndet))

            hel_list(ndet) = get_helement(src_det, det, ilut, det_list(:,ndet))

            call GenExcitations3 (src_det, ilut, det, flag, ex, par, &
                                  found_all, .false.)

        end do

        if (ndet /= nexcit) &
            call stop_all(this_routine,"Incorrect number of excitations found")

        ! now i should loop over the hel_list and check the actual number 
        ! of non-zero excitations 
        allocate(t_non_zero(nexcit))
        t_non_zero = .false.

        n_non_zero = 0

        do i = 1, ndet 
            if (abs(hel_list(i)) > 1.0e-6) then
                n_non_zero = n_non_zero + 1 
                t_non_zero(i) = .true. 
            end if 
        end do

        allocate(non_zero_list(0:niftot,n_non_zero))

        cnt = 0

        do i = 1, ndet 
            if (t_non_zero(i)) then
                cnt = cnt + 1
                non_zero_list(:,cnt) = det_list(:,i)
            end if
        end do

        deallocate(det_list)
        deallocate(t_non_zero)
        deallocate(hel_list)

        ! now i have a list of all non-zero excitations.. 

        ! for now this creates all excitation, independent if the matrix 
        ! element is actually zero.. 
        ! so also check matrix element! 

        ! Sort the dets, so they are easy to find by binary searching
        call sort(non_zero_list, ilut_lt, ilut_gt)

    end subroutine calc_all_excitations


end module
