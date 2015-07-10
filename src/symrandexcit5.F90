#include "macros.h"

module excit_gen_5

    use excit_gens_int_weighted, only: gen_single_4ind_ex, pgen_single_4ind, &
                                       get_paired_cc_ind, select_orb, &
                                       opp_spin_pair_contrib, &
                                       same_spin_pair_contrib
    use SymExcitDataMod, only: SpinOrbSymLabel, SymInvLabel, ScratchSize
    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles
    use SystemData, only: G1, tUHF, tStoreSpinOrbs, nbasis, nel
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, &
                                    construct_class_counts
    use dSFMT_interface, only: genrand_real2_dSFMT
    use procedure_pointers, only: get_umat_el
    use sym_general_mod, only: ClassCountInd
    use get_excit, only: make_double
    use bit_rep_data, only: NIfTot
    use UMatCache, only: gtid
    use constants
    use util_mod
    implicit none

contains

    subroutine gen_excit_4ind_weighted2 (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                     ExcitMat, tParity, pGen, HelGen, store)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_4ind_weighted'

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

        else

            ! OK, we want to do a double excitation
            ic = 2
            call gen_double_4ind_ex2 (nI, ilutI, nJ, ilutJ, ExcitMat, tParity, &
                                      pGen, store)
            pgen = pgen * pDoubles

        end if

        ! And a careful check!
#ifdef __DEBUG
        if (.not. IsNullDet(nJ)) then
            pgen2 = calc_pgen_4ind_weighted2(nI, ilutI, ExcitMat, ic, &
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


    function calc_pgen_4ind_weighted2 (nI, ilutI, ex, ic, ClassCountUnocc) &
            result(pgen)

        ! What is the probability of the excitation _from_ determinant nI
        ! described by the excitation matrix ex, and the excitation level ic,
        ! being generated according to the 4ind_weighted excitaiton generator?

        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer, intent(in) :: ClassCountUnocc(ScratchSize)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'calc_pgen_4ind_weighted'

        integer :: cc_index, src, tgt, id_src, id_tgt, n_id(nel)
        integer :: norb, label_index, orb, i, j, iSpn, sum_ml
        integer :: cc_i, cc_j, cc_i_final, cc_j_final, sym_product
        real(dp) :: cpt, cpt_tgt, cum_sum, cum_sums(2), int_cpt(2), ntot
        real(dp) :: cpt_pair(2), sum_pair(2)
        HElement_t :: hel


        if (ic == 1) then

            ! Singles are calculated in the same way for the _ex and the
            ! _reverse excitation generators
            pgen = pSingles * pgen_single_4ind (nI, ilutI, ex(1,1), ex(2,1))

        else if (ic == 2) then

            pgen = 0.0_dp

        end if

    end function


    subroutine gen_double_4ind_ex2 (nI, ilutI, nJ, ilutJ, ex, par, pgen, store)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        type(excit_gen_store_type), intent(in) :: store
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elecs(2), src(2), orbs(2), ispn, sum_ml, cc_a, cc_b
        integer :: sym_product
        real(dp) :: int_cpt(2), cum_sum(2)

        ! Pick the electrons in a weighted fashion
        call pick_weighted_elecs(nI, elecs, src, sym_product, ispn, sum_ml, &
                                 pgen)

        ! Select the A orbital _excluding_ knowledge of symmetry information.
        ! Only exclude terms that have no coupling elements.
        ! TODO: The pick_a_orb selection should only work if there are avail
        !       b orbs of the appropriate symmetry!
        orbs(1) = pick_a_orb(ilutI, src, int_cpt(1), cum_sum(1))
        if (orbs(1) == 0) then
            nJ(1) = 0
            return
        end if

        ! Determine the symmtry of the desired B orbital
        cc_a = ClassCountInd(orbs(1))
        cc_b = get_paired_cc_ind(cc_a, sym_product, sum_ml, iSpn)

        ! Select the B orbital, in the same way as before!!
        orbs(2) = select_orb(ilutI, src, cc_b, orbs(1), int_cpt(2), cum_sum(2))
        if (orbs(2) == 0) then
            nJ(1) = 0
            return
        end if

        ! Calculate the pgens. For now just discard the paired excitations
        ! until we have things working!!!
        if (orbs(2) < orbs(1)) then
            nJ(1) = 0
            return
        end if
        pgen = pgen * product(int_cpt) / product(cum_sum)

        ! And generate the actual excitation
        call make_double (nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                          ex, par)
        ilutJ = ilutI
        clr_orb (ilutJ, src(1))
        clr_orb (ilutJ, src(2))
        set_orb (ilutJ, orbs(1))
        set_orb (ilutJ, orbs(2))

    end subroutine gen_double_4ind_ex2


    subroutine pick_weighted_elecs(nI, elecs, src, sym_prod, ispn, sum_ml, &
                                   pgen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), src(2), sym_prod, ispn, sum_ml
        real(dp), intent(out) :: pgen

        logical, parameter :: tlinear = .true.

        real(dp) :: cum_list(nel), cum_sum, cpt, final_cpt, r
        integer :: i, ind1, ind2

        ! Pick the first electron uniformly, or according to <ii|ii>,
        ! and then pick the second electron according to <ji|ji>

        if (tlinear) then
            elecs(1) = 1 + floor(genrand_real2_dSFMT() * nel)
        else
            ! Still TODO: <ii|ii>
            ! Will impact the probability below too.
            ASSERT(.false.)
        end if

        ! Construct the weighted list <ji|ji>
        cum_sum = 0
        ind1 = gtID(nI(elecs(1)))
        do i = 1, nel
            if (i == elecs(1)) then
                cpt = 0
            else
                ! TODO: Can we do the biasing of parallel/opposite here?
                ind2 = gtID(nI(i))
                cpt = abs(get_umat_el(ind1, ind2, ind1, ind2))
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
        pgen = cpt / cum_sum

        ! To calculate the generation probability, we need to consider the
        ! possibility of having picked the electrons in the order j, i.
        ! Construct the reverse list
        cum_sum = 0
        ind2 = gtID(nI(elecs(2)))
        do i = 1, nel
            if (i == elecs(2)) then
                cpt = 0
            else
                ind1 = gtID(nI(i))
                cpt = abs(get_umat_el(ind2, ind1, ind2, ind1))
            endif
            if (i == elecs(1)) final_cpt = cpt
            cum_sum = cum_sum + cpt
        end do

        ! Adjust the probability for the j,i choice and then the 1/N one
        pgen = (pgen + (final_cpt / cum_sum)) / nel

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


    function pick_a_orb(ilut, src, cpt, cum_sum) result(orb)

        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer, intent(in) :: src(2)
        real(dp), intent(out) :: cpt, cum_sum
        integer :: orb
        character(*), parameter :: this_routine = 'ASSUMES RHF orbitals'

        real(dp) :: cum_arr(nbasis), r
        integer :: start_ind, nused, srcid(2)
        logical :: occa, occb

        ! Just in case. eeep.
        ! Note that we scale spin/spatial orbitals here with factors of 2
        ! which need more care if using spin orbitals
        if (tUHF .or. tStoreSpinOrbs) &
            call stop_all(this_routine, "ASSUMES RHF orbitals")

        if (is_beta(src(1)) .eqv. is_beta(src(2))) then

            ! If the two spins are aligned, then we only have to consider half
            ! of the available orbitals.

            ! Beta orbitals come before alpha orbs
            start_ind = 2
            if (is_beta(src(1))) start_ind = 1

            nused = nbasis / 2
            srcid = gtID(src)
            cum_sum = 0
            do orb = start_ind, nbasis, 2
                if (IsNotOcc(ilut, orb)) then
                    cpt = same_spin_pair_contrib(srcid(1), srcid(2), orb, -1)
                    cum_sum = cum_sum + cpt
                endif
                ! Note the gtID in this line...
                cum_arr(gtID(orb)) = cum_sum
            end do

        else

            ! If the spins are not parallel, then we need to consider all
            ! alpha/beta orbitals. On the plus side, we only need to consider
            ! one element at a time...
            nused = nbasis
            srcid = gtID(src)
            cum_sum = 0
            do orb = 1, nbasis - 1, 2

                occa = IsNotOcc(ilut, orb+1)
                occb = IsNotOcc(ilut, orb)
                if (.not. (occa .and. occb)) then
                    cpt = opp_spin_pair_contrib(srcid(1), srcid(2), orb, -1)
                end if
                if (.not. occb) then
                    cum_sum = cum_sum + cpt
                end if
                cum_arr(orb) = cum_sum
                if (.not. occa) then
                    cum_sum = cum_sum + cpt
                end if
                cum_arr(orb + 1) = cum_sum

            end do

        end if

        ! And exit if invalid
        if (cum_sum == 0) then
            orb = 0
            return
        end if

        ! Ensure that we get the correct spin when searching
        r = genrand_real2_dSFMT() * cum_sum
        orb = binary_search_first_ge(cum_arr(1:nused), r)
        if (orb == 1) then
            cpt = cum_arr(1)
        else
            cpt = cum_arr(orb) - cum_arr(orb - 1)
        end if

        ! If the spins are equal, the above selection will select in spatial
        ! orbitals. Convert to spin orbs
        if (is_beta(src(1)) .eqv. is_beta(src(2))) then
            orb = 2 * orb
            if (is_beta(src(1))) then
                orb = get_beta(orb)
            else
                orb = get_alpha(orb)
            end if
        end if

    end function pick_a_orb


end module
