#include "macros.h"

module excit_gens

    use SystemData, only: nel, nbasis
    use SymExcit3, only: CountExcitations3, GenExcitations3
    use SymExcitDataMod, only: SymLabelList2, SymLabelCounts2, OrbClassCount
    use sym_general_mod, only: ClassCountInd
    use FciMCData, only: excit_gen_store_type, pSingles
    use dSFMT_interface, only: genrand_real2_dSFMT
    use Determinants, only: get_helement, write_det
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use bit_rep_data, only: NIfTot
    use bit_reps, only: decode_bit_det
    use symdata, only: nSymLabels
    use procedure_pointers, only: get_umat_el
    use UMatCache, only: gtid
    use GenRandSymExcitNUMod, only: PickElecPair, gen_rand_excit
    use constants
    use get_excit, only: make_double
    use sort_mod
    use util_mod
    implicit none
    save

    ! Data for the 4ind-integral biasing scheme
    real(dp), allocatable :: epair_4ind_ints(:)

contains

    subroutine gen_excit_hel_weighted (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                       ExcitMat, tParity, pGen, HelGen, store)

        ! A really laborious, slow, explicit and brute force method to
        ! generating all excitations in proportion to their connection
        ! strength. This demonstrates the maximum possible value of tau that
        ! can be used.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in), target :: ilutI(0:NIfTot)
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
        integer(n_int), intent(in), target :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_hel_weighted'

        integer(n_int) :: iluts(0:NIfTot, nexcit)
        HElement_t :: hels(nexcit), hel_sum, hel_cum
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

        ! Just a quick sanity check
        if (allocated(epair_4ind_ints)) &
            call stop_all(this_routine, "Electron pair data array already &
                                        &allocated")

!        write(6,*) 'NBASIS', nbasispairs
!        do i = 1, nbasis
!        end do

    end subroutine


    subroutine gen_excit_4ind_weighted (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ExcitMat, tParity, pGen, HelGen, store)

        ! TODO: description
        !
        ! n.b. (ij|kl) <= sqrt( (ij|ij) * (kl|kl) )
        !      This provides quite a good description of the large elements

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in), target :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_4ind_weighted'


        ! Choose if we want to do a single or a double excitation
        ! TODO: We can (in principle) re-use this random number by subdivision
        if (genrand_real2_dSFMT() < pSingles) then

            ! Currently, just use the normal single excitation generator
            call gen_rand_excit (nI, ilutI, nJ, ilutJ, 1, ic, ExcitMat, &
                                 tParity, pGen, HElGen, store)

        else

            ! OK, we want to do a double excitation
            ic = 2
            call gen_double_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, tParity, &
                                     pGen)

        end if

    end subroutine

    
    subroutine gen_double_4ind_ex (nI, ilutI, nJ, ilutJ, ex, par, pgen)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elecs(2), syms(2), spn(2), orbs(2), src(2)
        real(dp) :: int_cpt(2), cum_sum(2)
        integer :: sym_product, ispn, sum_ml

        ! Initially, select a pair of electrons
        ! (Use routine from symrandexcit2, as it is known to work!)
        call PickElecPair (nI, elecs(1), elecs(2), sym_product, ispn, sum_ml, &
                           -1)
        src = nI(elecs)
        pgen = 2.0_dp / real(nel * (1 + nel), dp)

        ! Pick a a symmetry A. This selects symmetry B
        ! n.b. There is only one way to pick each case where symA == symB, and
        !      in all other cases we could have picked symB first --> double
        !      the generation probability.
        syms(1) = int(nSymLabels * genrand_real2_dSFMT())
        syms(2) = ieor(syms(1), sym_product)
        if (syms(1) /= syms(2)) &
            pgen = 2.0_dp * pgen

        ! What combinations of spins are permitted
        if (iSpn == 1) then
            spn = 2 ! All electrons are beta
        else if (iSpn == 3) then
            spn = 1 ! All electrons are alpha
        else if (genrand_real2_dSFMT() > 0.5) then
            ! We need to make a choice
            spn(1:2) = (/1, 2/)
        else
            spn(1:2) = (/2, 1/)
        end if

        ! Select the two orbitals
        orbs(1) = select_orb (ilutI, src, syms(1), spn(1), -1, int_cpt(1), & 
                              cum_sum(1))
        if (orbs(1) /= 0) &
            orbs(2) = select_orb (ilutI, src, syms(2), spn(2), orbs(1), &
                                  int_cpt(2), cum_sum(2))

        ! Just as a check, abort any excitation where we have excited two
        ! electrons into the same orbital
        if (any(orbs == 0)) then
            nJ(1) = 0
            return
        end if

        ! Adjust the probabilities. If we are selecting two orbitals from
        ! the same list, then we could have selected them either way around.
        if (syms(1) == syms(2) .and. spn(1) == spn(2)) then
            pgen = pgen * ( &
                   (int_cpt(1) / cum_sum(1) * int_cpt(2) / cum_sum(2)) &
                 + (int_cpt(2) / cum_sum(1) * &
                    int_cpt(1) / (cum_sum(1) - int_cpt(2))))
        else
            pgen = pgen * int_cpt(1) / cum_sum(1) * int_cpt(2) / cum_sum(2)
        end if

        ! And generate the actual excitation.
        call make_double (nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                          ex, par)

    end subroutine


    function select_orb (ilut, src, sym, spin, orb_pair, cpt, cum_sum) &
            result(orb)

        ! For an excitation from electrons 1,2, consider all of the pairs
        ! (e1 a|e1 a) == <e1 e1 | a a> and (e2 a | e2 a) == <e2 e2 | a a>
        ! as appropriate to bias selection of the electron

        integer, intent(in) :: sym, src(2), orb_pair, spin
        real(dp), intent(out) :: cpt, cum_sum
        integer(n_int), intent(in) :: ilut(0:NifTot)
        integer :: orb

        integer :: cc_index, label_index, orb_index, norb, i, orbid, srcid(2)

        ! Our biasing arrays must consider all of the possible orbitals with
        ! the correct symmetry.
        real(dp) :: cumulative_arr(nbasis), r

        ! How many orbitals are there with the given symmetry?
        cc_index = ClassCountInd(spin, sym, 0)
        label_index = SymLabelCounts2(1, cc_index)
        norb = OrbClassCount(cc_index)


        ! Construct a list of orbitals to excite to, and 
        cum_sum = 0
        srcid = gtID(src)
        do i = 1, norb
            
            orb = SymLabelList2(label_index + i - 1)
            if (IsNotOcc(ilut, orb) .and. orb /= orb_pair) then
                orbid = gtID(orb)
                cum_sum = cum_sum &
                        + abs(get_umat_el(srcid(1), srcid(1), orbid, orbid)) &
                        + abs(get_umat_el(srcid(2), srcid(2), orbid, orbid))
            end if
            cumulative_arr(i) = cum_sum

        end do

        ! If there are no available orbitals to pair with, we need to abort
        if (cum_sum == 0) then
            orb = 0
            return
        end if

        ! Binary search within this list to choose a value.
        r = genrand_real2_dSFMT() * cum_sum
        orb_index = binary_search_first_ge(cumulative_arr, r)

        ! And return the relevant value.
        orb = SymLabelList2(label_index + orb_index - 1)
        orbid = gtID(orb)
        cpt = abs(get_umat_el(srcid(1), srcid(1), orbid, orbid)) + &
              abs(get_umat_el(srcid(2), srcid(2), orbid, orbid))

    end function

end module
