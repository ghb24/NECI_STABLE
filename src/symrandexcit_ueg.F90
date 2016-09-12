#include "macros.h"

module ueg_excit_gens

    use SystemData, only: nel, nbasis, tOrbECutoff, ElecPairs, OrbECutoff, &
                          nmaxx, nmaxy, nmaxz, G1
    use dSFMT_interface, only: genrand_real2_dSFMT
    use SymExcitDataMod, only: KPointToBasisFn
    use FciMCData, only: excit_gen_store_type
    use DeterminantData, only: write_det
    use get_excit, only: make_double
    use bit_rep_data, only: NIfTot
    use sltcnd_mod, only: sltcnd_2
    use constants
    use util_mod
    implicit none

contains

    subroutine gen_ueg_excit (nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tPar, &
                              pgen, HelGen, store)

        ! This is a new excitation generator, modelled on the lines of the
        ! 4ind-weighted excitation generator used for determinants
        !
        ! N.B. in the UEG, only double excitations exist, and given a
        !      a choice of electron A, the choice of electron B is
        !      already made.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,2)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HelGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)

        real(dp) :: cum_sum, cum_arr(nbasis), pelec, porb, elem, r, testE
        integer :: eleci, elecj, spnb, ind, iSpn, iUnused
        integer :: orbi, orbj, orba, orbb
        integer :: ki(3), kj(3), ka(3), kb(3)

        ! Mitigate warnings
        HelGen = 0.0_dp; iUnused=exFlag; iUnused=store%nopen

        ! Pick a pair of electrons (i,j) to generate from.
        ! This uses a triangular mapping to pick them uniformly.
        ind = 1 + int(ElecPairs * genrand_real2_dSFMT())
        eleci = ceiling((1 + sqrt(1 + 8*real(ind, dp))) / 2)
        elecj = ind - ((eleci - 1) * (eleci - 2)) / 2
        pelec = 1.0_dp / real(ElecPairs, dp)

        ! Obtain the orbitals and their momentum vectors for the given elecs.
        orbi = nI(eleci)
        orbj = nI(elecj)
        ki = G1(orbi)%k
        kj = G1(orbj)%k
        ex(1, 1) = orbi
        ex(1, 2) = orbj

        ! And spin properties
        if (is_beta(orbi) .eqv. is_beta(orbj)) then
            if (is_beta(orbi)) then
                iSpn = 1 ! Beta-Beta
            else
                iSpn = 3 ! Alpha-Alpha
            end if
        else
            iSpn = 2 ! Alpha-Beta
        end if

        ! Loop through all available orbitals, A. If it is unoccupied, then
        ! find the orbital, B, that will complete the excitation given i,j.
        ! If this is also unoccupied, then contribute to the cumulative list
        ! for making selections
        cum_sum = 0.0_dp
        do orba = 1, nbasis

            ! TODO: Symmetry restrictions on A (if parallel, can't pick opp)
            elem = 0.0_dp
            if (IsNotOcc(ilutI, orba) .and. &
                (.not. ((iSpn == 1 .and. .not. is_beta(orba)) .or. &
                       (iSpn == 3 .and. is_beta(orba))))) then

                ! Obtain the new momentum vectors
                ka = G1(orba)%k
                kb = ki + kj - ka

                ! Is kb allowed by the size of the space?
                testE = sum(kb**2)
                if (abs(kb(1)) <= nmaxx .and. abs(kb(2)) <= nmaxy .and. &
                    abs(kb(3)) <= nmaxz .and. &
                    (.not. (tOrbECutoff .and. (testE > OrbECutoff)))) then

                    ! What would the spin of orbital b be? (2=al, 1=be)
                    if (iSpn == 1) then
                        spnb = 1
                    elseif (iSpn == 2) then
                        ! w.d: bug found by peter jeszenski and confirmed by 
                        ! simon! 
                        ! messed up alpa and beta spin here.. 
                        spnb = (-G1(orba)%Ms + 1)/2 + 1
!                         spnb = (G1(orba)%Ms)/2 + 1
                    elseif(iSpn == 3) then
                        spnb = 2
                    end if

                    ! Obtain the b-orbital
                    orbb = KPointToBasisFn(kb(1), kb(2), kb(3), spnb)

                    ! n.b. we enforce strict selection a-b, not b-a
                    if (orbb > orba .and. IsNotOcc(ilutI, orbb)) then

                        ! We don't need to worry about which a,b is which, as
                        ! we don't care about the overall sign.
                        ex(2, 1) = orba
                        ex(2, 2) = orbb
                        elem = abs(sltcnd_2(ex, .false.))
                    end if
                end if
            end if

            ! Increment the cumulative sum
            cum_sum = cum_sum + elem
            cum_arr(orba) = cum_sum

        end do

        ! If there are no available excitations, then we need to reject this
        ! excitation
        if (abs(cum_sum) < 1.0e-12_dp) then
            nJ(1) = 0
            return
        end if

        ! Pick a pair of orbitals A,B according to the cumulative probabilites
        ! already generated.
        r = genrand_real2_dSFMT() * cum_sum
        orba = binary_search_first_ge(cum_arr, r)

        ! Given A, calculate B in the same way as before
        ka = G1(orba)%k
        kb = ki + kj - ka
        if (iSpn == 1) then
            spnb = 1
        elseif (iSpn == 2) then
        ! w.d: bug found by peter jeszenski and confirmed by 
        ! simon! 
        ! messed up alpa and beta spin here.. 
        spnb = (-G1(orba)%Ms + 1)/2 + 1
!       spnb = (G1(orba)%Ms)/2 + 1
        elseif(iSpn == 3) then
            spnb = 2
        end if
        orbb = KPointToBasisFn(kb(1), kb(2), kb(3), spnb)

        ! Calculate the orbital selection probability
        elem = cum_arr(orba)
        if (orba > 1) elem = elem - cum_arr(orba - 1)
        porb = elem / cum_sum

        ! Construct and return the determinant
        call make_double(nI, nJ, eleci, elecj, orba, orbb, ex, tPar)
        ilutJ = ilutI
        clr_orb(ilutJ, orbi)
        clr_orb(ilutJ, orbj)
        set_orb(ilutj, orba)
        set_orb(ilutj, orbb)
        ic = 2
        pgen = pelec * porb

    end subroutine

end module
