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
        integer, intent(inout) :: CCOcc(ScratchSize)
        integer, intent(inout) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: scratch3(ScratchSize)
        logical, intent(inout) :: tFilled
        HElement_t, intent(out) :: HElGen

        real(dp) :: r
        character(*), parameter :: this_routine = 'gen_rand_excit3'

        ! Just in case
        ilutJ(0) = -1

        ! UEG and Hubbard interjection for now
        ! TODO: This should be made into its own fn-pointered case.
        if ((tUEG .and. tLatticeGens) .or. (tHub .and. tLatticeGens)) then
            call CreateExcitLattice (nI, iLut, nJ, tParity, ExcitMat, pGen)
            IC = 2
            return
        endif

        ! Count occupied/unoccupied orbitals in each symmetry class. This
        ! is an O[nel] operation. For efficiency, store these arrays
        ! between invocations of the excitation generator.
        if (.not. tFilled) then
            call ConstructClassCounts (nI, CCOcc, CCUnocc)
            tFilled = .true.
        endif

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

        ! Call the actual single/double excitation generators.
        if (IC == 2) then
            pGen = gen_double (nI, nJ, iLutI, CCUnocc)
        else
           ! call gen_single ()
        endif

    end subroutine


    function gen_double (nI, nJ, iLutI, CCUnocc) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer(n_int), intent(in) :: iLutI(0:niftot)
        real(dp) :: pGen

        real(dp) :: pElecs
        integer :: elecs(2), spn(2), sumMl, sym_prod, syms(2)

        ! Pick and unbiased, distinct, electron pair.
        pElecs = pick_elec_pair (nI, elecs, sym_prod, ispn, sumMl)

        call pick_syms (syms, elecs, sym_prod, spn, sumMl, CCUnocc)





        ! Return the final probability
        pGen = pElecs
    end function

    function pick_elec_pair (nI, elecs, sym_prod, spn, sumMl) result(pElec)

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
        integer, intent(out) :: elecs(2), sym_prod, spn(2), sumMl
        real(dp) :: pElec

        integer :: ind, orbs(2)

        ! Generate a random integer 1 <= i <= nel(nel-1)/2 (ElecPairs)
        ind = 1 + int(ElecPairs * genrand_real2_dSFMT())

        ! Generate the two indices, and obtain the associated orbitals
        elecs(1) = ceiling((1 + sqrt(1 + 8*real(ind,dp))) / 2)
        elecs(2) = ind - ((elecs(1) - 1) * (elecs(1) - 2)) / 2
        orbs = nI(elecs)

        ! Obtain the symmetry product label
        sym_prod = RandExcitSymLabelProd (SpinOrbSymLabel(orbs(1)), &
                                          SpinOrbSymLabel(orbs(2)))

        ! Obtain spins
        spn = get_spin(orbs(1))

        ! Store the Lz value to preserve
        if (tFixLz) sumMl = product(G1(orbs)%Ml)

        ! Return the generation probability
        pElec = 2 / real(nel*(nel-1), dp)

    end function
    
    subroutine pick_syms (syms, elecs, sym_prod, spn, sumMl, CCUnocc)

        integer, intent(out) :: syms(2)
        integer, intent(in) :: elecs(2), sym_prod, spn(2), sumMl
        integer, intent(in) :: CCUnocc(ScratchSize)

        integer :: num_pairs(0:nSymLabels)

        integer :: symA, symB, indA, indB, rint, tmp_tot

        ! TODO: Are we going to be able to store this one in ScratchSize
        !       arrays as well?
        num_pairs  = 0
        do symA = 0, nSymLabels - 1

            ! Get the required B symmetry to complement the chosen A sym.
            symB = RandExcitSymLabelProd(SymInvLabel(symA), sym_prod)

            ! We must avoid double counting A/B pairs
            if (symB < symA) cycle

            ASSERT(.not. tFixLz)
            indA = ClassCountInd(spn(1), symA, -1)
            indB = ClassCountInd(spn(2), symB, -1)
            if (symA == symB) then
                if (spn(1) /= spn(2)) then
                    num_pairs(symA) = (CCUnocc(indA) - 1) &
                                    * (CCUnocc(indB) - 1)
                else
                    num_pairs(symA) = (CCUnocc(indA) - 2)                                                         * (CCUnocc(indB) - 3) / 2
                endif
            else
                num_pairs(symA) = (CCUnocc(indA) - 1) * (CCUnocc(indB) - 1)
                if (spn(1) /= spn(2)) then
                    indA = ClassCountInd(spn(2), symA, -1)
                    indB = ClassCountInd(spn(1), symB, -1)
                    num_pairs(symA) = num_pairs(symA) + &
                                      (CCUnocc(indA) * CCUnocc(indB)) 
                endif
            endif
        enddo

        ! Count the total number of pairs.
        tot_pairs = sum(num_pairs)

        ! Select a symA/symB pair biased by the number of possible excitations
        ! which can be made into them.
        rint = 1 + int(genrand_real2_dSFMT() * tot_pairs)
        tmp_tot = 0
        do symA = 0, nSymLabels - 1
            tmp_tot = tmp_tot + num_pairs(symA)
            if (tmp_tot >= rint) then
                syms(1) = symA
                syms(2) = RandExcitsymLabelProd(SymInvLabel(symA), sym_prod)
                exit
            endif
        enddo

    end subroutine



end module

