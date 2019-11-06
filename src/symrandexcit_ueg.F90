#include "macros.h"

module ueg_excit_gens

    use SystemData, only: nel, nbasis, tOrbECutoff, ElecPairs, OrbECutoff, &
                          nmaxx, nmaxy, nmaxz, G1, TContact, tTrcorrExgen
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: excit_gen_store_type
    use DeterminantData, only: write_det
    use get_excit, only: make_double
    use bit_rep_data, only: NIfTot
    use sltcnd_mod, only: sltcnd_2_kernel, sltcnd_2
    use UMatCache, only: gtID
    use constants
    use util_mod
    use back_spawn, only: is_allowed_ueg_k_vector, get_orb_from_kpoints, get_ispn
    implicit none

contains

    subroutine gen_ueg_excit (nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tPar, &
                              pgen, HelGen, store, part_type)

        ! This is a new excitation generator, modelled on the lines of the
        ! 4ind-weighted excitation generator used for determinants
        !
        ! N.B. in the UEG, only double excitations exist, and given a
        !      a choice of electron A, the choice of electron B is
        !      already made.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,maxExcit)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HelGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        integer, intent(in), optional :: part_type

        unused_var(exFlag); unused_var(store); unused_var(part_type);
        HelGen = 0.0_dp

        ! W.D:
        ! split this functionality to allow back-spawning to reuse code
        call gen_double_ueg(nI, ilutI, nJ, ilutJ, tPar,ex, pgen)
        ic = 2
    end subroutine gen_ueg_excit

    subroutine gen_double_ueg(nI, ilutI, nJ, ilutJ, tPar, ex, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2,maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen

        integer :: elecs(2), eleci, elecj, ispn, orba, orbb, orbi, orbj
        real(dp) :: pelec, cum_sum, elem, porb, r, cum_arr(nBasis)

        ! Pick a pair of electrons (i,j) to generate from.
        ! This uses a triangular mapping to pick them uniformly.
        call pick_uniform_elecs(elecs, pelec)

        eleci = elecs(1)
        elecj = elecs(2)

        ! Obtain the orbitals and their momentum vectors for the given elecs.
        orbi = nI(eleci)
        orbj = nI(elecj)
        ex(1, 1) = orbi
        ex(1, 2) = orbj

        ! Loop through all available orbitals, A. If it is unoccupied, then
        ! find the orbital, B, that will complete the excitation given i,j.
        ! If this is also unoccupied, then contribute to the cumulative list
        ! for making selections
        if (TContact) then
                call create_ab_list_ua(nI,ilutI, [orbi,orbj], cum_arr, cum_sum)
        else
                call create_ab_list_ueg(ilutI, [orbi,orbj], cum_arr, cum_sum)
        endif

        ! If there are no available excitations, then we need to reject this
        ! excitation
        if (cum_sum < EPS) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! Pick a pair of orbitals A,B according to the cumulative probabilites
        ! already generated.
        r = genrand_real2_dSFMT() * cum_sum
        orba = binary_search_first_ge(cum_arr, r)

        ! write a cleaner implementation
        orbb = get_orb_from_kpoints(orbi, orbj, orba)

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
        pgen = pelec * porb

    end subroutine gen_double_ueg

    subroutine pick_uniform_elecs(elecs, pelec)
        integer, intent(out) :: elecs(2)
        real(dp), intent(out) :: pelec

        integer :: ind, eleci, elecj

        ind = 1 + int(ElecPairs * genrand_real2_dSFMT())
        elecs(1) = ceiling((1 + sqrt(1 + 8*real(ind, dp))) / 2)
        elecs(2) = ind - ((elecs(1) - 1) * (elecs(1) - 2)) / 2
        pelec = 1.0_dp / real(ElecPairs, dp)

    end subroutine pick_uniform_elecs

    subroutine create_ab_list_ueg(ilutI, src, cum_arr, cum_sum)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: src(2)
        real(dp), intent(out) :: cum_arr(nbasis), cum_sum

        integer :: ex(2,2), orba, orbb, ispn
        real(dp) :: elem, testE

        ex(1,:) = src

        ispn = get_ispn(src)

        cum_sum = 0.0_dp
        do orba = 1, nbasis

            ! TODO: Symmetry restrictions on A (if parallel, can't pick opp)
            elem = 0.0_dp
            if (IsNotOcc(ilutI, orba) .and. &
                (.not. ((iSpn == 1 .and. .not. is_beta(orba)) .or. &
                       (iSpn == 3 .and. is_beta(orba))))) then

                if (is_allowed_ueg_k_vector(src(1), src(2), orba)) then

                    orbb = get_orb_from_kpoints(src(1), src(2), orba)

                   ! n.b. we enforce strict selection a-b, not b-a
                    if (orbb > orba .and. IsNotOcc(ilutI, orbb)) then

                        ! We don't need to worry about which a,b is which, as
                        ! we don't care about the overall sign.
                        ex(2, 1) = orba
                        ex(2, 2) = orbb
                        elem = abs(sltcnd_2_kernel(ex))
                    end if
                end if
            end if

            ! Increment the cumulative sum
            cum_sum = cum_sum + elem
            cum_arr(orba) = cum_sum

        end do

    end subroutine create_ab_list_ueg

    subroutine create_ab_list_ua(nI,ilutI, src, cum_arr, cum_sum)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: src(2)
        real(dp), intent(out) :: cum_arr(nbasis), cum_sum

        integer :: ex(2,2), orba, orbb, ispn
        real(dp) :: elem, testE

        ex(1,:) = src

        ispn = get_ispn(src)

        cum_sum = 0.0_dp
        do orba = 1, nbasis

            ! TODO: Symmetry restrictions on A (if parallel, can't pick opp)
            elem = 0.0_dp
            if (IsNotOcc(ilutI, orba) .and. &
                (.not. ((iSpn == 1 .and. .not. is_beta(orba)) .or. &
                       (iSpn == 3 .and. is_beta(orba))))) then

                if (is_allowed_ueg_k_vector(src(1), src(2), orba)) then

                    orbb = get_orb_from_kpoints(src(1), src(2), orba)

                   ! n.b. we enforce strict selection a-b, not b-a
                    if (orbb > orba .and. IsNotOcc(ilutI, orbb)) then

                        ! We don't need to worry about which a,b is which, as
                        ! we don't care about the overall sign.
                        ex(2, 1) = orba
                        ex(2, 2) = orbb
                        elem = abs(sltcnd_2(nI,ex,.false.))
!                       elem = 1.0_dp
                    end if
                end if
            end if

            ! Increment the cumulative sum
            cum_sum = cum_sum + elem
            cum_arr(orba) = cum_sum

        end do

    end subroutine create_ab_list_ua

     function calc_pgen_ueg(ilutI, ex, ic) result(pgen)
        ! i also have to write a pgen recalculator for the pgens with this
        ! new UEG excitation generator.. i am a bit confused why this has
        ! not been done yet i have to admit..
        ! and i need this function if i want to use it with HPHF..
        integer, intent(in) :: ex(2,2), ic
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen

        real(dp) :: p_elec, p_orb, cum_arr(nBasis), cum_sum
        integer :: src(2), orb_a, tgt(2)

        if (ic /= 2) then
            pgen = 0.0_dp
            return
        end if

        src = get_src(ex)
        tgt = get_tgt(ex)

        ! i should return 0 if b < a.. since those excitations
        ! are never created.. but are the actually called ever? hm..
        ! the question is how does the ex() determination work actually?
        ! is it ever possible in the HPHF case eg. that b > a when
        ! comparing two dets? hm.. check!
!         if (tgt(1) > tgt(2)) then
!             pgen = 0.0_dp
!             return
!         end if

        ! now the real part..
        p_elec = 1.0_dp /real(ElecPairs, dp)

        orb_a = tgt(1)

        ! i need to recalct the cum_arr
        call create_ab_list_ueg(ilutI, src, cum_arr, cum_sum)

        if (cum_sum < EPS) then
            pgen = 0.0_dp
            return
        end if
        ! and now i have to check orba probability
        if (orb_a == 1) then
            p_orb = cum_arr(1) / cum_sum
        else
            p_orb = (cum_arr(orb_a) - cum_arr(orb_a - 1)) / cum_sum
        end if

        pgen = p_orb * p_elec

    end function calc_pgen_ueg

end module
