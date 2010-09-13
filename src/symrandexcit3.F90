#include "macros.h"

module symrandexcit3

    ! This is another version of the excitation generators. It creates 
    ! random excitations with a calculable, but non-uniform, probability.
    !
    ! Motivation:
    !     i) Use fewer random numbers than symrandexcit2
    !    ii) Generate (a bit) more uniform generation probabilities.

    use SystemData, only: nel, tFixLz, G1, ElecPairs
    use SymExcitDataMod, only: ScratchSize, SpinOrbSymLabel
    use dSFMT_interface, only: genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd
    use FciMCData, only: pDoubles
    use bit_reps, only: niftot
    use constants, only: dp, n_int
    implicit none

contains

    subroutine gen_rand_excit3 (nI, ilutI, nJ, ilutJ, exFlag, IC, ExcitMat, &
                                tParity, pGen, HElGen, tFilled, scratch1, &
                                scratch2, scratch3)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        integer, intent(inout) :: scratch1(ScratchSize)
        integer, intent(inout) :: scratch2(ScratchSize)
        integer, intent(inout) :: scratch3(ScratchSize)
        logical, intent(inout) :: tFilled
        HElement_t, intent(out) :: HElGen

        real(dp) :: r
        character(*), parameter :: this_routine = 'gen_rand_excit3'

        ! Just in case
        ilutJ(0) = -1

        ! If exFlag is 3, select singles or doubles randomly, according
        ! to the value in pDoubles. Otherwise exFlag = 1 gives a single,
        ! and exFlag = 2 gives a double.
        ASSERT(IC <= 3 .and. IC >= 1)
        IC = exFlag
        if (exFlag > 2) then
            r = genrand_real2_dSFMT()
            if (r < pDoubles) then
                IC = 2
            else
                IC = 1
            endif
        endif

        if (IC == 2) then
            pGen = gen_double (nI, nJ, iLutI)
        else
           ! call gen_single ()
        endif

    end subroutine


    function gen_double (nI, nJ, iLutI) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(in) :: iLutI(0:niftot)
        real(dp) :: pGen

        real(dp) :: pElecs
        integer :: elecs(2), ispn, sumMl, sym_prod

        ! Pick and unbiased, distinct, electron pair.
        pElecs = pick_elec_pair (nI, elecs, sym_prod, ispn, sumMl)

        ! Return the final probability
        pGen = pElecs
    end function

    function pick_elec_pair (nI, elecs, sym_prod, ispn, sumMl) result(pElec)

        ! Use a triangular indexing system.
        !  --> Only need one random number to pick two distinct electrons
        !      from N(N-1)/2 pairs.
        !
        ! i.e.  54  53  52  51         10  9   8   7
        !           43  42  41   ==>       6   5   4
        !               32  31                 3   2
        !                   21                     1
        !
        ! We can obtain the first index, A, by considering the largest
        ! integer, i, which can give an element on that row. For an integer
        ! 1 <= i <= npair:
        ! 
        !   --> A = ceil((1 + sqrt(1 + 8*i)) / 2)
        !   --> B = i - (A-1)(A-2)/2

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), sym_prod, ispn, sumMl
        real(dp) :: pElec

        integer :: ind, orbs(2)

        ! Generate a random integer 1 <= i <= nel(nel-1)/2 (ElecPairs)
        ind = ceiling(ElecPairs * genrand_real2_dSFMT())

        ! Generate the two indices, and obtain the associated orbitals
        elecs(1) = ceiling((1 + sqrt(1 + 8*real(ind,dp))) / 2)
        elecs(2) = ind - ((elecs(1) - 1) * (elecs(1) - 2)) / 2
        orbs = nI(elecs)

        ! Obtain the symmetry product label
        sym_prod = RandExcitSymLabelProd (SpinOrbSymLabel(orbs(1)), &
                                          SpinOrbSymLabel(orbs(2)))

        ! Obtain spin product
        if (product(G1(orbs)%Ms) == -1) then
            ! alpha / beta pair)
            ispn = 2
        elseif (G1(orbs(1))%Ms == 1) then
            ! alpha / alpha pair
            ispn = 3
        else
            ! beta / beta pair
            ispn = 1
        endif

        ! Store the Lz value to preserve
        if (tFixLz) sumMl = product(G1(orbs)%Ml)

        ! Return the generation probability
        pElec = 2 / real(nel*(nel-1), dp)

    end function


end module

