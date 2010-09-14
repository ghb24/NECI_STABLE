#include "macros.h"

module symrandexcit3

    ! This is another version of the excitation generators. It creates 
    ! random excitations with a calculable, but non-uniform, probability.
    !
    ! Motivation:
    !     i) Use fewer random numbers than symrandexcit2
    !    ii) Generate (a bit) more uniform generation probabilities.

    use SystemData, only: nel, tFixLz, G1, ElecPairs, tUEG, tHub, tLatticeGens
    use SymExcitDataMod, only: ScratchSize, SpinOrbSymLabel, SymInvLabel, &
                               SymLabelList2, SymLabelCounts2
    use SymData, only: nSymLabels
    use dSFMT_interface, only: genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, ClassCountInd, &
                                    construct_class_counts,CreateSingleExcit,&
                                    CreateExcitLattice
    use FciMCData, only: pDoubles
    use bit_reps, only: niftot
    use constants, only: dp, n_int, bits_n_int
    implicit none

contains

    subroutine gen_rand_excit3 (nI, ilutI, nJ, ilutJ, exFlag, IC, ExcitMat, &
                                tParity, pGen, HElGen, tFilled, CCOcc, &
                                CCUnocc, pair_list)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        integer, intent(inout) :: CCOcc(ScratchSize)
        integer, intent(inout) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: pair_list(ScratchSize)
        logical, intent(inout) :: tFilled
        HElement_t, intent(out) :: HElGen

        real(dp) :: r
        character(*), parameter :: this_routine = 'gen_rand_excit3'

        ! Just in case
        ilutJ(0) = -1

        ! UEG and Hubbard interjection for now
        ! TODO: This should be made into its own fn-pointered case.
        if ((tUEG .and. tLatticeGens) .or. (tHub .and. tLatticeGens)) then
            call CreateExcitLattice (nI, iLutI, nJ, tParity, ExcitMat, pGen)
            IC = 2
            return
        endif

        ! Count occupied/unoccupied orbitals in each symmetry class. This
        ! is an O[nel] operation. For efficiency, store these arrays
        ! between invocations of the excitation generator.
        if (.not. tFilled) then
            call construct_class_counts (nI, CCOcc, CCUnocc)
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
            pGen = gen_double (nI, nJ, iLutI, ExcitMat, tParity, CCUnocc, &
                               pair_list)
        else
            call CreateSingleExcit (nI, nJ, CCOcc, CCUnocc, &
                                    ilutI, ExcitMat, tParity, pGen)
           ! call gen_single ()
        endif

    end subroutine


    function gen_double (nI, nJ, iLutI, ExcitMat, tParity, CCUnocc, &
                         pair_list) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: pair_list(ScratchSize)
        integer(n_int), intent(in) :: iLutI(0:niftot)
        integer, intent(out) :: ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp) :: pGen

        real(dp) :: pElecs
        integer :: elecs(2), spn(2), orbs(2), sym_inds(2)
        integer :: sumMl, sym_prod, rint, tot_pairs

        ! Pick and unbiased, distinct, electron pair.
        pElecs = pick_elec_pair (nI, elecs, sym_prod, spn, sumMl)

        ! Pick a pair of symmetries, such that 
        tot_pairs = count_orb_pairs (sym_prod, spn, sumMl, CCUnocc, pair_list)

        ! If there are no possible excitations for the electron pair picked, 
        ! then we abort the excitation
        if (tot_pairs == 0) then
            nJ(1) = 0
            return
        endif

        ! Given a random number, the remainder of the generation is entirely
        ! deterministic
        rint = 1 + int(genrand_real2_dSFMT() * tot_pairs)

        ! Select a pair of symmetries to choose from
        call select_syms(rint, sym_inds, sym_prod, spn, sumMl, CCUnocc, &
                         pair_list)

        ! Select a pair of orbitals from the symmetries above.
        call select_orb_pair (rint, sym_inds, ilutI, orbs, CCUnocc)

        ! Generate the final determinant.
        call create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)
       
        ! Return the final probability
        pGen = pDoubles * pElecs / real(tot_pairs, dp)

    end function


    function pick_elec_pair (nI, elecs, sym_prod, spn, sumMl) result(pElec)

        ! Use a triangular indexing system.
        !  --> Only need one random number to pick two distinct electrons
        !      from N(N-1)/2 pairs.
        !
        ! i.e.              21                     1
        !               32  31   ==>           3   2
        !           43  42  41             6   5   4
        !       54  53  52  51         10  9   8   7
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
    
    function count_orb_pairs (sym_prod, spn, sumMl, CCUnocc, num_pairs) &
                              result(tot_pairs)

        integer, intent(in) :: sym_prod, spn(2), sumMl
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: num_pairs(0:nSymLabels-1)
        integer :: tot_pairs


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
                num_pairs(symA+1) = (CCUnocc(indA) - 1) * (CCUnocc(indB) - 1)
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

    end function

    subroutine select_syms (rint, sym_inds, sym_prod, spn, sumMl, CCUnocc, &
                            num_pairs)

        integer, intent(inout) :: rint, spn(2)
        integer, intent(out) :: sym_inds(2)
        integer, intent(in) :: sumMl, sym_prod
        integer, intent(in) :: CCUnocc(ScratchSize), num_pairs(0:nSymLabels-1)

        integer :: tmp_tot, syms(2), npairs, inds(2), tmp, symA

        ! Select a symA/symB pair biased by the number of possible 
        ! excitations which can be made into them.
        tmp_tot = 0
        do symA = 0, nSymLabels - 1
            if (tmp_tot + num_pairs(symA) >= rint) then
                syms(1) = symA
                syms(2) = RandExcitsymLabelProd(SymInvLabel(symA), sym_prod)
                exit
            endif
            tmp_tot = tmp_tot + num_pairs(symA)
        enddo

        ! Modify rint such that it now specifies which of the orbital pairs
        ! within the selected symmetry categories is desired.
        rint = rint - tmp_tot

        ! If there are more than one symmetry index corresponding to these
        !
        ! n.b. need to look in SymLabelList2(SymLabelCounts2(1,Ind)+i)
        !      which is a list of length SymLabelCounts2(2,ind)
        if (spn(1) /= spn(2) .and. syms(1) /= syms(2)) then
            inds = ClassCountInd(spn, syms, -1)
            npairs = (CCUnocc(inds(1)) - 1) * (CCUnocc(inds(2)) -1)
            if (rint > npairs) then
                ! Swap the spins around --> the other possibility, and adjust
                ! the random number to select a pair in this set.
                tmp = spn(1)
                spn(1) = spn(2)
                spn(2) = tmp
                rint = rint - npairs
            endif
        endif

        ! Return the symmetry indices, rather than the symmetry labels
        ! as that is what we will need for the selections.        
        sym_inds = ClassCountInd(spn, syms, -1)

    end subroutine


    subroutine select_orb_pair (rint, sym_inds, ilutI, orbs, CCUnocc)

        integer, intent(in) :: rint, sym_inds(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer :: orbs(2)
        character(*), parameter :: this_routine = 'select_orb_pair'

        integer :: a, b, i, nvac, pos, orb, offset, norbs

        if (sym_inds(1) == sym_inds(2)) then
            ! We are picking two orbitals from the same category
            ! --> Use the triangular scheme previously for selecting
            !     electrons.

            ! Select the positions of the two orbitals in the vacant list.
            a = ceiling((1 + sqrt(1 + 8*real(rint,dp))) / 2)
            b = i - ((a - 1) * (a - 2) / 2)
            orbs(1) = min(a, b)
            orbs(2) = max(a, b)

            nvac = 0
            pos = 1
            offset = SymLabelCounts2(1, sym_inds(1))
            norbs = SymLabelCounts2(2, sym_inds(1))
            do i = 0, norbs
                orb = SymLabelList2(offset + i)
                if (IsNotOcc(ilutI, orb)) then
                    nvac = nvac + 1
                    if (nvac == orbs(pos)) then
                        orbs(pos) = orb
                        pos = pos + 1
                        if (pos > 2) exit
                    endif
                endif
            enddo

            if (i == norbs) &
                call stop_all (this_routine, 'Unable to find enough vacant &
                                             &orbitals in SymLabelList2')
        else
            ! We are picking two orbitals from different categories
            ! --> use a 'rectangular', mapping scheme.
            orbs(1) = mod(rint - 1, CCUnocc(sym_inds(1))) + 1
            orbs(2) = floor((real(rint,dp) - 1) / CCUnocc(sym_inds(1))) + 1

            do pos = 1, 2
                nvac = 0
                offset = SymLabelCounts2(1, sym_inds(pos))
                norbs = SymLabelCounts2(2, sym_inds(pos))
                do i = 0, norbs
                    orb = SymLabelList2(offset + 1)
                    if (IsNotOcc(ilutI, orb)) then
                        nvac = nvac + 1
                        if (nvac == orbs(pos)) then
                            orbs(pos) = orb
                            exit
                        endif
                    endif
                enddo

                if (i == norbs) &
                    call stop_all (this_routine, 'Unable to find enough &
                                          &vacant orbitals in SymLabelList2')
            enddo

        endif

    end subroutine


    subroutine create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)

        integer, intent(in) :: nI(nel), elecs(2), orbs(2)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2)
        logical, intent(out) :: tParity

        ! TODO: Update the bit representations here as well?
        ExcitMat(1,:) = elecs
        ExcitMat(2,:) = orbs
        nJ = nI

        ! TODO FindExcitDet (excit.F) is not modularised/interfaced yet
        call FindExcitDet(ExcitMat, nJ, 2, tParity)

    end subroutine

end module

