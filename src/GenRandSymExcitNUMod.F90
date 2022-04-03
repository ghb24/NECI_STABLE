#include "macros.h"

MODULE GenRandSymExcitNUMod
!This is a new random excitation generator. It still creates excitations with a normalised and calculable probability,
!but these probabilities are not uniform. They should also be quicker to generate, and have small or no excitation generators
!to store. nI is the root determinant and nJ is the excitation.
!This currently is only available for abelian symmetry. This is because I label the symmetries as integers, but this
!algorithm could theoretically be applied relatively easily to non-abelian symmetry. The symmetries would have to be stored
!in a better way. For molecular systems, symmetry irrep is just symmetry_label-1. Symlabels will need to be used for more
!complicated symmetry tables.

!THEORY

    !  P[i,j,a,b] = P_doub x ( P[i,j] x P[a,b|i,j] ) x (M-N)/(M-N-Q)
    !  P[a,b|i,j] = P[b|i,j] x P[a|b,i,j] + P[a|i,j] x P[b|a,i,j]
    !  P[i,j] = 2/N(N-1) always - others are more difficult to calculate
    !  Difficulties arise if there is ever no allowed final unoccupied orbital to pick.

    !  Spin and spatial symmetry needs to be considered.
    !  The fourth orbital chosen (i.e. one of the unoccupied ones) will have its
    !  symmetry class specified by the classes of the other ones.

    !  i and j are picked randomly. They can be alpha,beta | beta,alpha | alpha,alpha | beta,beta
    !  a is then picked from all orbitals with the allowed spin. This means that a beta,beta i,j pair
    !  will force a to also be beta. There are also a number of forbidden a orbitals which are counted.
    !  These are forbidden since they have no possible b orbital which will give rise to a symmetry and
    !  spin allowed unoccupied a,b pair. The number of these orbitals, Q, is needed to calculate the
    !  normalised probability of generating the excitation.
    use SystemData, only: ALAT, iSpinSkip, tFixLz, iMaxLz, tUEG, tNoFailAb, &
                          tLatticeGens, tHub, nEl, G1, nBasis, nBasisMax, &
                          tNoSymGenRandExcits, Arr, nMax, tCycleOrbs, &
                          nOccAlpha, nOccBeta, ElecPairs, MaxABPairs, &
                          tKPntSym, lzTot, tNoBrillouin, tUseBrillouin, &
                          tUEG2, kvec, tAllSymSectors, NMAXX, NMAXY, NMAXZ, &
                          tOrbECutoff, OrbECutoff, tGUGA, t_pcpp_excitgen, &
                          t_new_real_space_hubbard
    use CalcData, only: tSpinProject
    use FciMCData, only: pDoubles, psingles, iter, excit_gen_store_type, iluthf
    use Parallel_neci
    use IntegralsData, only: UMat
    use Determinants, only: get_helement, write_det
    use SymData, only: nSymLabels, TwoCycleSymGens, SymLabelList, &
                       SymLabelCounts
    use dSFMT_interface, only: genrand_real2_dSFMT
    use SymExcitDataMod
    use DetBitOps, only: FindExcitBitDet, EncodeBitDet, detbiteq
    use excitation_types, only: SingleExc_t
    use sltcnd_mod, only: sltcnd_excit
    use constants, only: dp, n_int, bits_n_int, maxExcit, stdout
    use bit_reps, only: NIfTot
    use sym_mod, only: mompbcsym, GetLz
    use timing_neci
    use sym_general_mod
    use get_excit, only: make_single, make_double
    use procedure_pointers, only: get_umat_el
    use back_spawn, only: pick_virtual_electrons_double_hubbard
    use back_spawn, only: pick_occupied_orbital_hubbard, check_electron_location, get_ispn
    use neci_intfce
    use bit_rep_data, only: test_flag

    IMPLICIT NONE

contains

    subroutine gen_rand_excit(nI, ilut, nJ, ilutnJ, exFlag, IC, ExcitMat, &
                              tParity, pGen, HElGen, store, part_type)

        ! This routine is the same as GenRandSymexcitNu, but you can pass in
        ! the class count arrays so that they do not have to be recalculated
        ! each time for the same excitation. If tFilled is false, it will
        ! assume they are unfilled and calculate them, returning the arrays
        ! and tFilled=T.
        ! The two arrays want to be integers, both of size (1, 1:nSymLabels)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: iLut(0:niftot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2, maxExcit)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        type(excit_gen_store_type), intent(inout), target :: store

        ! Not used
        integer(n_int), intent(out) :: ilutnJ(0:niftot)
        HElement_t(dp), intent(out) :: HElGen

        integer, intent(in), optional :: part_type

        real(dp) :: r
        character(*), parameter :: this_routine = 'gen_rand_excit'

        ! Just in case
        ilutnJ(0) = -1
        HElGen = 0.0_dp

        IF ((tUEG .and. tLatticeGens) .or. (tHub .and. tLatticeGens)) THEN
            call CreateExcitLattice(nI, iLut, nJ, tParity, ExcitMat, pGen, part_type)
            IC = 2
            RETURN
        end if

        !TODO: Not quite sure what conditions we need to check for now...
        if (.not. store%tFilled) then
!First, we need to do an O[N] operation to find the number of occupied alpha electrons, number of occupied beta electrons
!and number of occupied electrons of each symmetry class and spin. This is similar to the ClassCount array.
!This has the format (Spn,sym), where Spin=1,2 corresponding to alpha and beta.
!For molecular systems, sym runs from 0 to 7.
!This is stored to save doing this multiple times, but shouldn't be too costly an operation.
            CALL construct_class_counts(nI, store%ClassCountOcc, &
                                        store%ClassCountUnocc)
            store%tFilled = .true.
        end if

!ExFlag is 1 for singles, 2 for just doubles, and 3 for both.
        IF (ExFlag == 3) THEN
!Choose whether to generate a double or single excitation. Prob of generating a double is given by pDoub.
            pDoubNew = pDoubles

            r = genrand_real2_dSFMT()
            IF (r < pDoubNew) THEN
!A double excitation has been chosen to be created.
                IC = 2
            ELSE
                IC = 1
            end if
        else if (ExFlag == 2) THEN
            IC = 2
            pDoubNew = 1.0_dp
        else if (ExFlag == 1) THEN
            IC = 1
            pDoubNew = 0.0_dp
        ELSE
            CALL Stop_All(this_routine, "Error in choosing excitations to create.")
        end if

        IF (IC == 2) THEN
            CALL CreateDoubExcit(nI, nJ, store%ClassCountUnocc, ILUT, &
                                 ExcitMat, tParity, pGen)
        ELSE
            CALL CreateSingleExcit(nI, nJ, store%ClassCountOcc, &
                                   store%ClassCountUnocc, ILUT, ExcitMat, &
                                   tParity, pGen)


        end if

    end subroutine

    SUBROUTINE CreateDoubExcit(nI, nJ, ClassCountUnocc2, ILUT, ExcitMat, tParity, pGen)
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel), ExcitMat(2, maxExcit)
        integer, intent(in) :: ClassCountUnocc2(ScratchSize)
        integer(n_int), intent(in) :: iLut(0:NIfTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity

        integer :: nExcitOtherWay, orbB, nExcitB, SpinOrbA, OrbA, SymA, SymB
        integer :: nExcitA, sumMl, mlA, mlB, iSpn, Elec1Ind, Elec2Ind
        integer :: SymProduct, ForbiddenOrbs
        logical :: tAOrbFail

!First, we need to pick an unbiased distinct electron pair.
!These have symmetry product SymProduct, and spin pair iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha
        CALL PickElecPair(nI, Elec1Ind, Elec2Ind, SymProduct, iSpn, SumMl, -1)


!This routine finds the number of orbitals which are allowed by spin,
! but not part of any spatial symmetry allowed unoccupied pairs.
!This number is needed for the correct normalisation of the probability of drawing
! any given A orbital since these can be chucked and redrawn.
        IF (tNoSymGenRandExcits) THEN
            CALL FindNumForbiddenOrbsNoSym(ForbiddenOrbs, ClassCountUnocc2, iSpn)
        ELSE
            CALL FindNumForbiddenOrbs(ForbiddenOrbs, ClassCountUnocc2, SymProduct, iSpn, SumMl)
        end if

!Now we have to pick the first unoccupied orbital. If an orbital is not present in any allowed pairs,
! it is chucked and a new one drawn.
!The number NExcit is the number of unoccupied orbitals that the orbital was chosen from
! (including symmetry-forbidden orbital pairs) Arguments:
        CALL PickAOrb(nI,iSpn,ILUT,ClassCountUnocc2,NExcitA,Elec1Ind,Elec2Ind,SpinOrbA,OrbA,SymA,SymB, &
      &             SymProduct,SumMl,MlA,MlB,ForbiddenOrbs,tAOrbFail)
        IF(tAOrbFail) THEN
            nJ(1)=0
            pGen=HUGE(0.0_dp)
            RETURN
        end if

!This routine will pick an unoccupied orbital at random from a specified spin and symmetry class.
!There should definitely be a possible spinorbital, since A was chosen so that there was one.
!We have to make sure with alpha/alpha or beta/beta pairs and when SymProduct=0,
! that we don't choose the same unoccupied orbital.
!If we do this, then we should chuck and redraw,
! since there should definitely be another allowed spinorbital in the class.
!We return the number of allowed B's for the A we picked in NExcitB,
! however we also need to know the number of allowed A's if we
!had picked B first. This will be returned in NExcitOtherWay.
        CALL PickBOrb(nI, iSpn, ILUT, ClassCountUnocc2, SpinOrbA, OrbA, SymA, OrbB, SymB, NExcitB, MlA, MlB, NExcitOtherWay)

        call make_double(nI, nJ, elec1ind, elec2ind, orbA, orbB, &
                         ExcitMat, tParity)

        CALL FindDoubleProb(ForbiddenOrbs, NExcitA, NExcitB, NExcitOtherWay, pGen)

    END SUBROUTINE CreateDoubExcit

!This routine finds the probability of creating the excitation.
! See the header of the file for more information on how this works.
    SUBROUTINE FindDoubleProb(ForbiddenOrbs, NExcitA, NExcitB, NExcitOtherWay, pGen)
        INTEGER, INTENT(IN) :: ForbiddenOrbs, NExcitA, NExcitB, NExcitOtherWay
        real(dp), INTENT(OUT) :: pGen!,PabGivenij

pGen = pDoubles * ((1.0_dp / real(NExcitB, dp)) + (1.0_dp / real(NExcitOtherWay, dp))) / (REAL((ElecPairs * (NExcitA - ForbiddenOrbs)), dp))

    END SUBROUTINE FindDoubleProb

    SUBROUTINE PickBOrb(nI, iSpn, ILUT, ClassCountUnocc2, SpinOrbA, OrbA, SymA, OrbB, SymB, NExcit, MlA, MlB, NExcitOtherWay)
        integer, intent(in) :: nI(nel), iSpn, SpinOrbA, OrbA, SymA, SymB
        integer, intent(in) :: MlA, MlB
        integer, intent(in) :: ClassCountUnocc2(ScratchSize)
        integer, intent(out) :: nExcitOtherWay, nExcit, OrbB
        integer(n_int), intent(in) :: iLut(0:NIfTot)
        integer :: norbs, i, z, ind, ChosenUnocc, attempts, SpinOrbB
        real(dp) :: r

!We want to calculate the number of possible B's given the symmetry
! and spin it has to be since we have already picked A.
!We have calculated in NExcit the number of orbitals available for B given A,
! but we also need to know the number of orbitals to choose from for A IF
!we had picked B first.
        ind = 0
        IF (iSpn == 2) THEN
!If iSpn=2, then we want to find a spinorbital of the opposite spin of SpinOrbA
            IF (SpinOrbA == -1) THEN
!We have already picked a beta orbital, so now we want to pick an alpha orbital.
! Find out how many of these there are.
                NExcit=ClassCountUnocc2(ClassCountInd(1,SymB,MlB))
                NExcitOtherWay=ClassCountUnocc2(ClassCountInd(2,SymA,MlA))
                SpinOrbB=1  !This is defined differently to SpinOrbA. 1=Alpha, 2=Beta.
            ELSE
!Want to pick an beta orbital.
                NExcit = ClassCountUnocc2(ClassCountInd(2, SymB, MlB))
                NExcitOtherWay = ClassCountUnocc2(ClassCountInd(1, SymA, MlA))
                SpinOrbB = 2
            end if
        else if (iSpn == 1) THEN
!Definitely want a beta orbital
            NExcit = ClassCountUnocc2(ClassCountInd(2, SymB, MlB))
            NExcitOtherWay = ClassCountUnocc2(ClassCountInd(2, SymA, MlA))
            SpinOrbB = 2
        ELSE
!Definitely want an alpha orbital
            NExcit = ClassCountUnocc2(ClassCountInd(1, SymB, MlB))
            NExcitOtherWay = ClassCountUnocc2(ClassCountInd(1, SymA, MlA))
            SpinOrbB = 1
        end if

        IF ((iSpn /= 2) .and. (SymA == SymB) .and. (MlA == MlB)) THEN
!In this case, we need to check that we do not pick the same orbital as OrbA. If we do this, then we need to redraw.
!Only when SymProduct=0 will the classes of a and b be the same, and the spins will be different if iSpn=2,
!  so this is the only possibility of a clash.
            NExcit = NExcit - 1     !Subtract 1 from the number of possible orbitals since we cannot choose orbital A.
            NExcitOtherWay = NExcitOtherWay - 1     !The same goes for the probabilities the other way round.
        end if

!All orbitals with the specified symmetry and spin should be allowed unless it is OrbA.
!There will be NExcit of these. Pick one at random.
!Check that orbital is not in ILUT and is not = OrbA (Although this can only happen
!in the circumstance indicated earlier).
!Now we need to choose the final unoccupied orbital.
!There are two ways to do this.
!  We can either choose the orbital we want out of the NExcit possible unoccupied orbitals.
!It would then be necessary to cycle through all orbitals of that symmetry and spin,
! only counting the unoccupied ones to find the correct determinant.
!Alternatively, we could pick orbitals at random and redraw until we find an allowed one.
! This would probably be preferable for larger systems.
        IF (tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals in symmetry class with needed spin to find allowed one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
            r = genrand_real2_dSFMT()
            ChosenUnocc = INT(NExcit * r) + 1

!Now run through all allowed orbitals until we find this one.
            IF (tNoSymGenRandExcits) THEN
                nOrbs = nBasis / 2      !No symmetry, therefore all orbitals of allowed spin possible to generate.
            ELSE
                Ind = ClassCountInd(SpinOrbB, SymB, MlB)
                nOrbs = SymLabelCounts2(2, Ind)
            end if
            z = 0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
            do i = 0, nOrbs - 1
                IF (tNoSymGenRandExcits) THEN
                    OrbB = (2 * (i + 1)) - (SpinOrbB - 1)
                ELSE
!Find the spin orbital index. SymLabelCounts has the index of the state for the given symmetry.
                    OrbB = SymLabelList2(SymLabelCounts2(1, Ind) + i)
                ENDIF

!Find out if the orbital is in the determinant, or is the other unocc picked
                IF ((.not. (BTEST(ILUT((OrbB - 1) / bits_n_int), MOD((OrbB - 1), bits_n_int)))) .and. (OrbB /= OrbA)) THEN
!The orbital is not found in the original determinant - increment z
                    z = z + 1
                    IF (z == ChosenUnocc) THEN
!We have got to the determinant that we want to pick.
                        EXIT
                    end if
                end if

            end do

!We now have our final orbitals.
            IF (z /= ChosenUnocc) THEN
                write(stdout, *) "This is a problem, since there should definitely be an allowed beta orbital once alpha is chosen..."
                CALL Stop_All("PickBOrb", "Could not find allowed unoccupied orbital to excite to.")
            end if

        ELSE
! METHOD 2 (Keep drawing orbitals from the desired symmetry and spin until we find one unoccupied)
! =========================

            IF (tNoSymGenRandExcits) THEN
                nOrbs = nBasis / 2
            ELSE
                Ind = ClassCountInd(SpinOrbB, SymB, MlB)
                nOrbs = SymLabelCounts2(2, Ind)
            end if

            Attempts = 0
            do while (.true.)

!Draw randomly from the set of orbitals
                r = genrand_real2_dSFMT()
                ChosenUnocc = INT(nOrbs * r)
                IF (tNoSymGenRandExcits) THEN
                    OrbB = (2 * (ChosenUnocc + 1)) - (SpinOrbB - 1)
                ELSE
                    OrbB = SymLabelList2(SymLabelCounts2(1, Ind) + ChosenUnocc)
                ENDIF

!Find out if orbital is in nI or not. Accept if it isn't in it.
                IF ((.not. (BTEST(ILUT((OrbB - 1) / bits_n_int), MOD((OrbB - 1), bits_n_int)))) .and. (OrbB /= OrbA)) THEN
!Orbital not in nI. Accept.
                    EXIT
                end if

                IF (Attempts > 1000) THEN
                    write(stdout, *) "Cannot find double excitation unoccupied orbital after 1000 attempts..."
                    write(stdout, *) "This is a problem, since there should definitly be an allowed beta orbital once alpha is chosen..."
                    write(stdout, *) "nI: "
                    call write_det(6, nI, .TRUE.)
                    write(stdout, *) "iSpn: ", iSpn
                    write(stdout, *) "ClassCountUnocc2: ", ClassCountUnocc2(:)
                    write(stdout, *) "NExcit", NExcit
                    CALL Stop_All("PickBOrb", "Cannot find double excitation unoccupied orbital after 250 attempts...")
                end if
                Attempts = Attempts + 1

            end do

        end if

    END SUBROUTINE PickBOrb

!This routine does the same as the FindNumForbiddenOrbs routine,
! but is optimised for when there are no spatial symmetry considerations.
    SUBROUTINE FindNumForbiddenOrbsNoSym(ForbiddenOrbs, ClassCountUnocc2, iSpn)
        integer, intent(in) :: ClassCountUnocc2(ScratchSize), iSpn
        integer, intent(out) :: ForbiddenOrbs

!We know that all orbitals are totally symmetric, and that the symproduct=0

        ForbiddenOrbs = 0
        IF (iSpn == 2) THEN
            IF (ClassCountUnocc2(1) == 0) THEN
!There are no unoccupied alpha orbitals - are there any beta spins which are now forbidden?
                ForbiddenOrbs = ClassCountUnocc2(2)
            end if
            IF (ClassCountUnocc2(2) == 0) THEN
!There are no unoccupied beta orbitals - are there any alpha spins which are now forbidden?
                ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(1)
            end if

        else if (iSpn == 1) THEN
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
            IF (ClassCountUnocc2(2) == 1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                ForbiddenOrbs = 1
            end if
        else if (iSpn == 3) THEN
            IF (ClassCountUnocc2(1) == 1) THEN
                ForbiddenOrbs = 1
            end if
        end if

    END SUBROUTINE FindNumForbiddenOrbsNoSym

!This routine finds the number of orbitals which are allowed by spin,
! but not part of any spatial symmetry allowed unoccupied pairs.
!This number is needed for the correct normalisation of the probability of
! drawing any given A orbital since these can be chucked and redrawn.
!For Lz symmetry, it is generally quicker to count the allowed orbitals, and subtract from all possible ones,
!  rather than directly counting
!the forbidden ones.
    SUBROUTINE FindNumForbiddenOrbs(ForbiddenOrbs, ClassCountUnocc2, SymProduct, iSpn, SumMl)
        integer, intent(in) :: ClassCountUnocc2(ScratchSize)
        integer, intent(in) :: SumMl, iSpn, SymProduct
        integer, intent(out) :: ForbiddenOrbs
        integer :: i, k, Ind, AllowedOrbs, SymInd, OrbAMl, SymOrbs, ConjSym

        ForbiddenOrbs = 0

        IF (tFixLz) THEN

            AllowedOrbs = 0
            IF (iSpn == 2) THEN
                Ind = 1
                do k = -iMaxLz, iMaxLz
                    OrbAMl = SumMl - k
                    IF ((k <= OrbAMl) .and. (abs(OrbAMl) <= iMaxLz)) THEN
                        do i = 0, nSymLabels - 1
                            SymInd = ClassCountInd(1, IEOR(SymProduct, i), OrbAMl)  !Alpha of the corresponding a orbital
                            IF (ClassCountUnocc2(Ind) /= 0) THEN
                                !Check the beta conjugate orbital
                                SymOrbs = ClassCountUnocc2(SymInd + 1)
                                IF (SymOrbs /= 0) THEN
                                    IF (k == OrbAMl) THEN
                                        !The Ml symmetries are the same! Don't double count.
                                        AllowedOrbs = AllowedOrbs + SymOrbs
                                    ELSE
                                        AllowedOrbs = AllowedOrbs + SymOrbs + ClassCountUnocc2(Ind)
                                    end if
                                end if
                            end if
                            IF (ClassCountUnocc2(Ind + 1) /= 0) THEN
                                SymOrbs = ClassCountUnocc2(SymInd)
                                IF (SymOrbs /= 0) THEN
                                    IF (k == OrbAMl) THEN
                                        !The Ml symmetries are the same! Don't double count.
                                        AllowedOrbs = AllowedOrbs + SymOrbs
                                    ELSE
                                        AllowedOrbs = AllowedOrbs + SymOrbs + ClassCountUnocc2(Ind + 1)
                                    end if
                                end if
                            end if
                            Ind = Ind + 2
                        end do
                    ELSE
                        !Move onto the next k-block of B orbitals.
                        Ind = Ind + nSymLabels * 2
                    end if
                end do
                ForbiddenOrbs = nBasis - NEl - AllowedOrbs
            else if (iSpn == 3) THEN  !alpha/alpha - run through all alpha orbitals
                Ind = 1
                do k = -iMaxLz, iMaxLz
                    OrbAMl = SumMl - k
                    IF ((k <= OrbAMl) .and. (abs(OrbAMl) <= iMaxLz)) THEN
                        do i = 0, nSymLabels - 1
                            IF (ClassCountUnocc2(Ind) /= 0) THEN
                                IF ((SymProduct == 0) .and. (OrbAMl == k)) THEN
                                    IF (ClassCountUnocc2(Ind) > 1) THEN
                                        AllowedOrbs = AllowedOrbs + ClassCountUnocc2(Ind)
                                    end if
                                ELSE
                                    SymOrbs = ClassCountUnocc2(ClassCountInd(1, IEOR(SymProduct, i), OrbAMl))
                                    IF (SymOrbs /= 0) THEN
                                        IF (k == OrbAMl) THEN
                                            !The Ml symmetries are the same! Don't double count.
                                            AllowedOrbs = AllowedOrbs + SymOrbs
                                        ELSE
                                            AllowedOrbs = AllowedOrbs + SymOrbs + ClassCountUnocc2(Ind)
                                        end if
                                    end if
                                end if
                            end if
                            Ind = Ind + 2
                        end do
                    ELSE
                        Ind = Ind + nSymLabels * 2
                    end if
                end do
                ForbiddenOrbs = (nBasis / 2) - nOccAlpha - AllowedOrbs
            ELSE
                Ind = 2
                do k = -iMaxLz, iMaxLz
                    OrbAMl = SumMl - k
                    IF ((k <= OrbAMl) .and. (abs(OrbAMl) <= iMaxLz)) THEN
                        do i = 0, nSymLabels - 1
                            IF (ClassCountUnocc2(Ind) /= 0) THEN
                                IF ((SymProduct == 0) .and. (OrbAMl == k)) THEN
                                    IF (ClassCountUnocc2(Ind) > 1) THEN
                                        AllowedOrbs = AllowedOrbs + ClassCountUnocc2(Ind)
                                    end if
                                ELSE
                                    SymOrbs = ClassCountUnocc2(ClassCountInd(2, IEOR(SymProduct, i), OrbAMl))
                                    IF (SymOrbs /= 0) THEN
                                        IF (k == OrbAMl) THEN
                                            !The Ml symmetries are the same! Don't double count.
                                            AllowedOrbs = AllowedOrbs + SymOrbs
                                        ELSE
                                            AllowedOrbs = AllowedOrbs + SymOrbs + ClassCountUnocc2(Ind)
                                        end if
                                    end if
                                end if
                            end if
                            Ind = Ind + 2
                        end do
                    ELSE
                        Ind = Ind + nSymLabels * 2
                    end if
                end do
                ForbiddenOrbs = (nBasis / 2) - nOccBeta - AllowedOrbs
            end if

        ELSE
            !Not Lz symmetry...
            IF (iSpn == 2) THEN
!                write(stdout,*) "Alpha/Beta"
!i,j are an alpha/beta pair. The number of forbidden orbitals includes all alphas and betas.

                Ind = 1

                do i = 0, nSymLabels - 1
!Run though all symmetries of possible "a" orbital syms. If there aren't any, then we
!know the corresponding "b" orbitals are excluded.
                    IF (ClassCountUnocc2(Ind) == 0) THEN
!This symmetry has no unoccupied alpha orbitals - does its symmetry conjugate have any
!unoccupied beta orbitals which are now forbidden?
!If there are no unoccupied orbitals in this conjugate symmetry, then it won't increase
!the forbidden orbital number, since it can never be chosen.
!                        ConjSym=IEOR(SymProduct,i)
                        ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(ClassCountInd(2,   &
                &           RandExcitSymLabelProd(SymProduct, SymInvLabel(i)), 0))
                        !No unocc alphas in i, therefore all betas in ConjSym are forbidden
!                        write(stdout,*) ClassCountUnocc2(2,ConjSym),i,ConjSym
                    end if
                    IF (ClassCountUnocc2(Ind + 1) == 0) THEN
!This symmetry has no unoccupied beta orbitals - does its symmetry conjugate have any
!unoccupied alpha orbitals which are now forbidden?
!If there are no unoccupied orbitals in this conjugate symmetry, then it won't increase
!the forbidden orbital number, since it can never be chosen.
                        ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(ClassCountInd(1,   &
                &           RandExcitSymLabelProd(SymProduct, SymInvLabel(i)), 0))
                    end if
                    Ind = Ind + 2
                end do

            else if (iSpn == 1) THEN
                Ind = 2
                IF (.not. tKPntSym) THEN
                    !With molecular systems, the irreps are their own inverses, so it is a little simpler to do the
                    !two cases seperately.
                    IF (SymProduct /= 0) THEN
!i,j are a beta/beta pair. The number of forbidden orbitals is just betas
                        do i = 0, nSymLabels - 1
                            IF (ClassCountUnocc2(Ind) == 0) THEN
                                !                            ConjSym=IEOR(SymProduct,i)
                                ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(ClassCountInd(2, RandExcitSymLabelProd(SymProduct, i), 0))
                            end if
                            Ind = Ind + 2
                        end do
                    ELSE
!There is a subtle point here, which could change the probabilities.
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
                        do i = 0, nSymLabels - 1
                            IF (ClassCountUnocc2(Ind) == 1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                                ForbiddenOrbs = ForbiddenOrbs + 1
                            end if
                            Ind = Ind + 2
                        end do
                    end if
                ELSE
                    !With KPntSym, we have to work out if we are in the case that sym_a^* = sym_b
                    !Unfortunately, I don't think you can tell from SymProduct when this case is going to be satisfied.
                    do i = 0, nSymLabels - 1        !Run over symmetries of the orbitals
                        IF (ClassCountUnocc2(Ind) <= 1) THEN
                            ConjSym = RandExcitSymLabelProd(SymProduct, SymInvLabel(i))
                            IF (ConjSym == i) THEN
                                !A and B come from the same symmetry, so we must have more than one
                                !orbitals available from there..
                                IF (ClassCountUnocc2(Ind) == 1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                                    ForbiddenOrbs = ForbiddenOrbs + 1
                                end if
                            ELSE
                                IF (ClassCountUnocc2(Ind) == 0) THEN
                                    ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(ClassCountInd(2, ConjSym, 0))
                                end if
                            end if
                        end if

                        Ind = Ind + 2
                    end do
                end if

            else if (iSpn == 3) THEN
                Ind = 1
                IF (.not. tKPntSym) THEN
                    IF (SymProduct /= 0) THEN
!i,j are a alpha/alpha pair. The number of forbidden orbitals is just alphas
                        do i = 0, nSymLabels - 1
                            IF (ClassCountUnocc2(Ind) == 0) THEN
                                ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(ClassCountInd(1, IEOR(SymProduct, i), 0))
                            end if
                            Ind = Ind + 2
                        end do
                    ELSE
!There is a subtle point here, which could change the probabilities.
!If the symmetry product of the occupied orbitals is 0, then the a,b pair want to be taken from the same class.
!This means that if there is only one spin-allowed orbital in that class, it has no symmetry-allowed pairs, and so is forbidden.
                        do i = 0, nSymLabels - 1
                            IF (ClassCountUnocc2(Ind) == 1) THEN
!The one alpha orbital in this class is forbidden, since it cannot form a pair.
                                ForbiddenOrbs = ForbiddenOrbs + 1
                            end if
                            Ind = Ind + 2
                        end do
                    end if
                ELSE
                    !With KPntSym, we have to work out if we are in the case that sym_a^* = sym_b
                    !Unfortunately, I don't think you can tell from SymProduct when this case is going to be satisfied.
                    do i = 0, nSymLabels - 1        !Run over symmetries of the orbitals
                        IF (ClassCountUnocc2(Ind) <= 1) THEN
                            ConjSym = RandExcitSymLabelProd(SymProduct, SymInvLabel(i))
                            IF (ConjSym == i) THEN
                                !A and B come from the same symmetry, so we must have more than one
                                !orbitals available from there..
                                IF (ClassCountUnocc2(Ind) == 1) THEN
!The one beta orbital in this class is forbidden, since it cannot form a pair.
                                    ForbiddenOrbs = ForbiddenOrbs + 1
                                end if
                            ELSE
                                IF (ClassCountUnocc2(Ind) == 0) THEN
                                    !There aren't any a's in this pair of syms, so the b's are forbidden
                                    ForbiddenOrbs = ForbiddenOrbs + ClassCountUnocc2(ClassCountInd(1, ConjSym, 0))
                                end if
                            end if
                        end if
                        Ind = Ind + 2
                    end do
                end if
            end if
        end if

    END SUBROUTINE FindNumForbiddenOrbs

!Choose the first unoccupied orbital to excite to. Spatial symmetry does not need to be considered, but spin symmetry does.
!After one has been chosen, we have to make sure that an allowed symmetry pair can be formed, otherwise the orbital is
!chucked and a new one drawn.
!Arguments:     NExcit = Number of possible spin-allowed unoccupied spinorbitals, including forbidden orbs (these will be chucked)
!               SpinOrbA = Spin of the chosen spin-orbital. 1 is alpha, -1 is beta.
!               OrbA = Index of the chosen spinorbital
!               SymB = Symmetry required of the second unoccupied spinorbital, so that sym(A) x sym(B) = SymProduct
!               SymProduct = (intent in), the symmetry of electron pair chosen
!               SumMl = Sum of the Ml values for the two picked electrons
!               MlA = Ml of the a orbital chosen
!               MlB = Required Ml of the b orbital still to be chosen
    SUBROUTINE PickAOrb(nI, iSpn, ILUT, ClassCountUnocc2, NExcit, Elec1Ind, Elec2Ind, SpinOrbA, OrbA, SymA, SymB, &
    &   SymProduct, SumMl, MlA, MlB, ForbiddenOrbs, tAOrbFail)
        INTEGER, INTENT(IN) :: nI(NEl), iSpn, Elec1Ind, Elec2Ind, ForbiddenOrbs, SymProduct, SumMl, ClassCountUnocc2(ScratchSize)
        INTEGER, INTENT(OUT) :: SpinOrbA, SymA, MlA, MlB, NExcit, SymB, OrbA
        INTEGER(KIND=n_int), INTENT(IN) :: ILUT(0:NIfTot)
        LOGICAL, INTENT(OUT) :: tAOrbFail
        INTEGER :: AttemptsOverall, ChosenUnocc, z, i, Attempts
        real(dp) :: r

        IF (iSpn == 2) THEN
!There is no restriction on whether we choose an alpha or beta spin, so there are nBasis-NEl possible spinorbitals to choose from.
!Therefore the spinOrbA variable has to be set after we have chosen one. For the other iSpn types, we can set it now.
            NExcit = nBasis - NEl
        else if (iSpn == 1) THEN
!SpinOrbA indicates the spin of the chosen A orb. This is only really useful for iSpn=2
            SpinOrbA = -1  !Going to pick a beta orb
            NExcit = (nBasis / 2) - nOccBeta  !This is the number of unoccupied beta spinorbitals there are.
        ELSE    !iSpn is 3 - alpha/alpha pair.
            SpinOrbA = 1 !Going to pick an alpha orb
            NExcit = (nBasis / 2) - nOccAlpha !This is the number of unoccupied alpha spinorbitals there are.
        end if

        IF (NExcit == ForbiddenOrbs) THEN
            tAOrbFail = .true.
            RETURN
        ELSE
            tAOrbFail = .false.
        end if

        AttemptsOverall = 0
        do while (.true.)
!Keep drawing unoccupied orbitals, until we find one which has an allowed partner to form a symmetry-allowed unoccupied pair.

            IF (iSpn == 2) THEN
!Electrons chosen were an alpha/beta pair, therefore first randomly chosen orbital can be an alpha OR beta orbital - no restriction.

                IF (tCycleOrbs .or. ((NExcit - ForbiddenOrbs) <= 3)) THEN
! METHOD 1 (Run though all orbitals to find desired one out of NExcit-ForbiddenOrbs...now
!only one random number needed, and guarenteed to find excitation.)
! ==========================

!Choose the unoccupied orbital to excite to
                    r = genrand_real2_dSFMT()
                    ChosenUnocc = INT((NExcit - ForbiddenOrbs) * r) + 1

                    z = 1
                    do i = 0, nBasis - 1

                        !We need to find if allowed
                        IF (.not. BTEST(ILUT(i / bits_n_int), MOD(i, bits_n_int))) THEN
                            !Is not in the original determinant - allow if sym allowed
                            !Now check that its symmetry allowed
                            SpinOrbA = G1(i + 1)%Ms
                            IF (IsAOrbSymAllowed(iSpn, i + 1, SpinOrbA, SymProduct, SumMl, SymA, SymB, MlA, MlB, ClassCountUnocc2)) THEN
                                IF (z == ChosenUnocc) THEN
                                    !We have found our allowed orbital
                                    OrbA = i + 1
                                    RETURN
                                end if
                                z = z + 1
                            end if
                        end if

                    end do

                    CALL Stop_All("PickAOrb", "Should not get here - have not found allowed A Orb")

                ELSE
! METHOD 2 (Keep drawing orbitals randomly until we find one unoccupied). This should
!be more efficient, unless we have v. small basis sets.
! =========================

                    Attempts = 0
                    do while (.true.)

!Draw randomly from the set of orbitals
                        r = genrand_real2_dSFMT()
                        ChosenUnocc = INT(nBasis * r)

!Find out if orbital is in nI or not. Accept if it isn't in it.
                        IF (.not. (BTEST(ILUT(ChosenUnocc / bits_n_int), MOD(ChosenUnocc, bits_n_int)))) THEN
!Orbital not in nI. Accept.
                            EXIT
                        end if

                        IF (Attempts > 250) THEN
                            write(stdout, *) "Cannot find A unoccupied orbital after 250 attempts..."
                            call write_det(6, nI, .TRUE.)
                            CALL Stop_All("PickAOrb", "Cannot find A unoccupied orbital after 250 attempts...")
                        end if
                        Attempts = Attempts + 1

                    end do

                    OrbA = ChosenUnocc + 1  !This is the allowed orbital

                end if

                SpinOrbA = G1(OrbA)%Ms

            ELSE
!We are either constrained to choose a beta orbital, or an alpha orbital - we know that there are NExcit of these.

                IF (tCycleOrbs .or. ((NExcit - ForbiddenOrbs) <= 3)) THEN
! METHOD 1 (Run though all orbitals to find desired one out of NExcit-ForbiddenOrbs)
! ==========================

                    r = genrand_real2_dSFMT()
                    ChosenUnocc = INT((NExcit - ForbiddenOrbs) * r) + 1
                    SpinOrbA = 0  !The spin of the chosen A is not important, since we have already defined it from iSpn

                    z = 1
                    do i = 1, nBasis / 2
                        IF (iSpn == 1) THEN
!We want to run through all beta orbitals (odd numbered basis function)
                            OrbA = (2 * i) - 1
                        ELSE
!We want to run through all alpha orbitals (odd numbered basis functions)
                            OrbA = (2 * i)
                        end if

                        !We need to find if allowed
                        IF (.not. BTEST(ILUT((OrbA - 1) / bits_n_int), MOD((OrbA - 1), bits_n_int))) THEN
                            !Is not in the original determinant - allow
                            !Now check that its symmetry allowed
                            IF (IsAOrbSymAllowed(iSpn, OrbA, SpinOrbA, SymProduct, SumMl, SymA, SymB, MlA, MlB, ClassCountUnocc2)) THEN
                                IF (z == ChosenUnocc) THEN
                                    !We have found our allowed orbital
                                    !OrbA is the allowed orbital
                                    RETURN
                                end if
                                z = z + 1
                            end if
                        end if

                    end do

                    CALL Stop_All("PickAOrb", "Should not get here - have not found allowed A Orb")

                ELSE
! METHOD 2 (Keep drawing orbitals randomly until we find one unoccupied). This should be more
!efficient, unless we have small basis sets.
! =========================

                    Attempts = 0
                    do while (.true.)

!Draw randomly from the set of orbitals
                        r = genrand_real2_dSFMT()
                        ChosenUnocc = INT((nBasis / 2) * r) + 1
                        IF (iSpn == 1) THEN
!We only want to choose beta orbitals(i.e. odd).
                            ChosenUnocc = (2 * ChosenUnocc) - 1
                        ELSE
!We only want to choose alpha orbitals(i.e. even).
                            ChosenUnocc = 2 * ChosenUnocc
                        end if

!Find out if orbital is in nI or not. Accept if it isn't in it.
                        IF (.not. (BTEST(ILUT((ChosenUnocc - 1) / bits_n_int), MOD((ChosenUnocc - 1), bits_n_int)))) THEN
!Orbital not in nI. Accept.
                            EXIT
                        end if

                        IF (Attempts > 250) THEN
                            write(stdout, *) "Cannot find A unoccupied orbital after 250 attempts..."
                            call write_det(6, nI, .TRUE.)
                            CALL Stop_All("PickAOrb", "Cannot find A unoccupied orbital after 250 attempts...")
                        end if
                        Attempts = Attempts + 1

                    end do

                    OrbA = ChosenUnocc  !This is the allowed orbital

                end if

            end if

!We now need to test whether this orbital has any symmetry-allowed unoccupied orbitals to form a pair with.
!To test this, we need to find the needed symmetry of B, in order that Sym(A) x Sym(B) = SymProduct
            IF (IsAOrbSymAllowed(iSpn, OrbA, SpinOrbA, SymProduct, SumMl, SymA, SymB, MlA, MlB, ClassCountUnocc2)) THEN
                !OrbA picked is allowed, exit from loop
                EXIT
            end if

            IF (AttemptsOverall > (nBasis * 20)) THEN
                write(stdout, *) "Cannot find first allowed unoccupied orbital for given i,j pair after ", nBasis * 20, " attempts."
                write(stdout, *) "It may be that there are no possible excitations from this i,j pair, in which case "
                write(stdout, *) "the given algorithm is inadequate to describe excitations from such a small space."
                write(stdout, *) "Try reverting to old excitation generators."
                call write_det(6, nI, .TRUE.)
                write(stdout, *) "***", NExcit, ForbiddenOrbs
                write(stdout, *) "I,J pair; sym_i, sym_j: ", nI(Elec1Ind), nI(Elec2Ind), G1(nI(Elec1Ind))%Sym%S, G1(nI(Elec2Ind))%Sym%S
                CALL Stop_All("PickAOrb", "Cannot find first allowed unocc orb for double excitation")
            end if
            AttemptsOverall = AttemptsOverall + 1

        end do

    END SUBROUTINE PickAOrb

!This routine will look at an orbital (OrbA) and check whether it is an allowed A orbital to
!pick, i.e. it has allowed B orbitals, given the i,js.
    LOGICAL FUNCTION IsAOrbSymAllowed(iSpn, OrbA, SpinOrbA, SymProduct, SumMl, SymA, SymB, MlA, MlB, ClassCountUnocc2)
        INTEGER, INTENT(IN) :: iSpn, OrbA, SpinOrbA, SymProduct, SumMl, ClassCountUnocc2(ScratchSize)
        INTEGER, INTENT(OUT) :: SymA, SymB, MlA, MlB

        IsAOrbSymAllowed = .false.
        IF (tNoSymGenRandExcits) THEN
            SymA = 0
            SymB = 0
            MlA = 0
            MlB = 0
        else if (tKPntSym) THEN
            SymA = SpinOrbSymLabel(OrbA)
            SymB = RandExcitSymLabelProd(SymInvLabel(SymA), SymProduct)
            MlB = 0
            MlA = 0
        ELSE
            SymA = INT(G1(OrbA)%Sym%S, 4)
            SymB = IEOR(SymA, SymProduct)
            MlB = 0
            MlA = 0
            IF (tFixLz) THEN
                MlA = G1(OrbA)%Ml
                MlB = SumMl - MlA
            end if
        end if

        IF (abs(MlB) <= iMaxLz) THEN
!Make sure that the B orbital that we would need to pick to conserve momentum is actually in the available range of Ml values.
            IF (iSpn == 2) THEN
!We want an alpha/beta unocc pair.
                IF (SpinOrbA == 1) THEN
!We have picked an alpha orbital - check to see if there are allowed beta unoccupied orbitals from the SymB Class.
                    IF (ClassCountUnocc2(ClassCountInd(2, SymB, MlB)) /= 0) THEN
!Success! We have found an allowed A orbital!
                        IsAOrbSymAllowed = .true.
                    end if
                ELSE
!We have picked a beta orbital - check to see if there are allowed alpha unoccupied orbitals from the SymB Class.
                    IF (ClassCountUnocc2(ClassCountInd(1, SymB, MlB)) /= 0) THEN
!Success! We have found an allowed A orbital!
                        IsAOrbSymAllowed = .true.
                    end if
                end if
            else if (iSpn == 1) THEN
!We want a beta/beta pair.
                IF ((SymB /= SymA) .or. (MlA /= MlB)) THEN
!Check to see if there are any unoccupied beta orbitals in the SymB Class.
                    IF (ClassCountUnocc2(ClassCountInd(2, SymB, MlB)) /= 0) THEN
!Success! We have found an allowed A orbital!
                        IsAOrbSymAllowed = .true.
                    end if
                ELSE
!We want an orbital from the same class. Check that this isn't the only unoccupied beta orbital in the class.
                    IF (ClassCountUnocc2(ClassCountInd(2, SymB, MlB)) /= 1) THEN
!Success! We have found an allowed A orbital!
                        IsAOrbSymAllowed = .true.
                    end if
                end if
            ELSE
!We want an alpha/alpha pair.
                IF ((SymA /= SymB) .or. (MlA /= MlB)) THEN
!Check to see if there are any unoccupied alpha orbitals in the SymB Class.
                    IF (ClassCountUnocc2(ClassCountInd(1, SymB, MlB)) /= 0) THEN
!Success! We have found an allowed A orbital!
                        IsAOrbSymAllowed = .true.
                    end if
                ELSE
!We want an orbital from the same class. Check that this isn't the only unoccupied alpha orbital in the class.
                    IF (ClassCountUnocc2(ClassCountInd(1, SymB, MlB)) /= 1) THEN
!Success! We have found an allowed A orbital!
                        IsAOrbSymAllowed = .true.
                    end if
                end if
            end if
        end if

    END FUNCTION IsAOrbSymAllowed

!This routine takes determinant nI and returns two randomly chosen electrons, whose index in nI is Elec1Ind and Elec2Ind.
!These electrons have symmetry product SymProduct and spin pairing iSpn, where iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha.
!If IndInp = -1, the pair is picked randomly with prob = 1/ElecPairs. Otherwise, it will choose electron pair given by index IndInp.
    SUBROUTINE PickElecPair(nI, Elec1Ind, Elec2Ind, SymProduct, iSpn, SumMl, IndInp)
        INTEGER, INTENT(IN) :: nI(NEl), IndInp
        INTEGER, INTENT(OUT) :: Elec1Ind, Elec2Ind, SymProduct, iSpn, SumMl
        INTEGER :: Ind, X, K, Orb1, Orb2
        real(dp) :: r
!Triangular indexing system.
!This is used for picking two distinct electrons out of all N(N-1)/2 pairs.
!
!   12  13  14  15          1   2   3   4
!       23  24  25              5   6   7
!           34  35      =>          8   9
!               45                      10

!For a given index, ind, there are [N(N-1)/2 - ind] elements at positions larger
!than ind. Call this number X.
!We want to find out how many rows there are after the row containing the element
!ind. t rows has t(t+1)/2 elements in it.
!Therefore, to find out the number of rows after ind, we want to find the largest K, such that K(K+1)/2 =< X
!The solution to this is that K =< (SQRT(8*X+1)-1)/2, therefore we can find K (the
!largest integer which satisfies this).
!We then know the number of rows after the element ind. Therefore, since there are N-1
!rows in total, we know we are on row N-1-K.
!This gives us the index of the first electron.
!To find the second electron (i.e. the column), we know that out of the X elements in
!positions larger than ind, K(K+1)/2 are in the next rows.
!This means that X - K(K+1)/2 are in the same row. There are N-(N-1-K) = 1+K elements
!in the row chosen, and so the number of elements into the
!row it is, is given by (1+K) - (X-K(K+1)/2). However, in row z, the first column
!index is z+1. Therefore the index of the second electron is
!(1+K) - (X-K(K+1)/2) + N-K-1 = N-X+(K(K+1)/2).

!        ElecPairs=(NEl*(NEl-1))/2

! If we want to find an index randomly, IndInp will be -1.
        IF (IndInp == -1) THEN
!Find an index randomly.
            r = genrand_real2_dSFMT()
            Ind = INT(ElecPairs * r) + 1
        ELSE
!Otherwise we are looking for a specific electron pair specified by IndInp
            Ind = IndInp
        end if

!X is number of elements at positions larger than ind
        X = ElecPairs - Ind
!K is the number of complete rows after the element ind
        K = INT((SQRT(8.0_dp * REAL(X, dp) + 1.0_dp) - 1.0_dp) / 2.0_dp)
        Elec1Ind = NEl - 1 - K
        Elec2Ind = NEl - X + ((K * (K + 1)) / 2)

        Orb1 = nI(Elec1Ind)
        Orb2 = nI(Elec2Ind)

!We now want to find the symmetry product label of the two electrons, and the spin product of the two electrons.
        SymProduct = RandExcitSymLabelProd(SpinOrbSymLabel(Orb1), SpinOrbSymLabel(Orb2))

        IF ((G1(Orb1)%Ms) * (G1(Orb2)%Ms) == -1) THEN
!We have an alpha beta pair of electrons.
            iSpn = 2
        ELSE
            IF (G1(Orb1)%Ms == 1) THEN
!We have an alpha alpha pair of electrons.
                iSpn = 3
            ELSE
!We have a beta beta pair of electrons.
                iSpn = 1
            end if
        end if

        SumMl = G1(Orb1)%Ml + G1(Orb2)%Ml

    END SUBROUTINE PickElecPair

    SUBROUTINE CheckIfSingleExcits(ElecsWNoExcits, ClassCount2, ClassCountUnocc2, nI)
        INTEGER, intent(out) :: ElecsWNoExcits
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        INTEGER, intent(in) :: nI(NEl)
        integer :: i

!First, we need to find out if there are any electrons which have no possible excitations.
!This is because these will need to be redrawn and so will affect the probabilities.
        ElecsWNoExcits = 0

        IF (tFixLz .or. tKPntSym) THEN
!Here, we also have to check that the electron is momentum allowed.
!Since there are many more irreps, it will be quicker here to check all electrons, rather than all the symmetries.

            do i = 1, NEl

                IF (G1(nI(i))%Ms == 1) THEN
                    IF (ClassCountUnocc2(ClassCountInd(1, SpinOrbSymLabel(nI(i)), G1(nI(i))%Ml)) == 0) THEN
                        ElecsWNoExcits = ElecsWNoExcits + 1
                    end if
                ELSE
                    IF (ClassCountUnocc2(ClassCountInd(2, SpinOrbSymLabel(nI(i)), G1(nI(i))%Ml)) == 0) THEN
                        ElecsWNoExcits = ElecsWNoExcits + 1
                    end if
                end if
            end do

!            do i=1,ScratchSize
!!Run through all labels
!                IF((ClassCount2(i).ne.0).and.(ClassCountUnocc2(i).eq.0)) THEN
!!If there are electrons in this class with no possible unoccupied orbitals in the
!same class, these electrons have no single excitations.
!                    ElecsWNoExcits=ElecsWNoExcits+ClassCount2(i)
!                end if
!            end do

        ELSE

            do i = 1, ScratchSize
!Run through all labels
                IF ((ClassCount2(i) /= 0) .and. (ClassCountUnocc2(i) == 0)) THEN
!If there are electrons in this class with no possible unoccupied orbitals in the same
!class, these electrons have no single excitations.
                    ElecsWNoExcits = ElecsWNoExcits + ClassCount2(i)
                end if
            end do
        end if

!Rather than choosing a double now if there are no singles, just return a null det.
!        IF(ElecsWNoExcits.eq.NEl) THEN
!!There are no single excitations from this determinant at all. This means the
!probability to create a double excitation = 1
!!Then we will create a double excitation instead.
!            pDoubNew=1.0_dp
!            RETURN
!        end if
!
    END SUBROUTINE CheckIfSingleExcits

    !> Wrapper function for creating a uniform single excitation
    subroutine uniform_single_excit_wrapper(nI, ilutI, nJ, ilutJ, ex, tpar, store, pgen)
        implicit none
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pGen
        type(excit_gen_store_type), intent(inout), target :: store

        ! First, count the number of available orbitals per irrep
        if (.not. store%tFilled) then
            call construct_class_counts(nI, store%ClassCountOcc, &
                                        store%ClassCountUnocc)
            store%tFilled = .true.
        end if
        pDoubNew = 1.0 - pSingles
        ! Then, chose uniformly from them
        call createSingleExcit(nI, nJ, store%ClassCountOcc, store%classCountUnocc, ilutI, &
                               ex, tPar, pGen)

        if (nJ(1) /= 0) then
            ! set the output ilut
            ilutJ = ilutI
            clr_orb(ilutJ, ex(1, 1))
            set_orb(ilutJ, ex(2, 1))
        else
            ilutJ = 0_n_int
        end if
    end subroutine uniform_single_excit_wrapper

    SUBROUTINE CreateSingleExcit(nI, nJ, ClassCount2, ClassCountUnocc2, ILUT, ExcitMat, tParity, pGen)
        INTEGER :: ElecsWNoExcits, i, Attempts, nOrbs, z, Orb
        INTEGER :: Eleci, ElecSym, nI(NEl), nJ(NEl), NExcit, iSpn, ChosenUnocc
        INTEGER :: ExcitMat(2, maxExcit)
        INTEGER :: ClassCount2(ScratchSize)
        INTEGER :: ClassCountUnocc2(ScratchSize), ElecK, SymIndex
        INTEGER(KIND=n_int) :: ILUT(0:NIfTot)
        real(dp) :: r, pGen
        LOGICAL :: tParity
        character(*), parameter :: t_r = 'CreateSingleExcit'

        CALL CheckIfSingleExcits(ElecsWNoExcits, ClassCount2, ClassCountUnocc2, nI)
        IF (ElecsWNoExcits == NEl) THEN
!There are no single excitations from this determinant - return a null excitation
            nJ(1) = 0
            pgen = 0.0_dp
            RETURN
        end if

        Attempts = 0
        elecK = 0
        orb = 0
        do while (.true.)

!Choose an electron randomly...
            r = genrand_real2_dSFMT()
            Eleci = INT(NEl * r) + 1

!Find symmetry of chosen electron
            IF (tNoSymGenRandExcits) THEN
                ElecSym = 0
            ELSE
!For abelian symmetry, the irrep of i and a must be the same.
!For solids, this means that the excitation must be within the same k-point
                ElecSym = SpinOrbSymLabel(nI(Eleci))
                ElecK = G1(nI(Eleci))%Ml
            end if

            IF (G1(nI(Eleci))%Ms == 1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                iSpn = 1
            ELSE
!Beta orbital
                iSpn = 2
            end if

            SymIndex = ClassCountInd(iSpn, ElecSym, ElecK)
            NExcit = ClassCountUnocc2(SymIndex)

            IF (NExcit /= 0) EXIT    !Have found electron with allowed excitations

            IF (Attempts > 250) THEN
                write(stdout, *) "Cannot find single excitation from electrons after 250 attempts..."
                call write_det(6, nI, .true.)
                write(stdout, *) "***"
                call stop_all(t_r, "Cannot find single excitation from &
                                   &electrons after 250 attempts...")
            end if
            Attempts = Attempts + 1

        end do

!Now we need to choose the unoccupied orbital for the chosen electron.
!There are two ways to do this. We can either choose the orbital we want out of the
!NExcit possible unoccupied orbitals.
!It would then be necessary to cycle through all orbitals of that symmetry and spin,
!only counting the unoccupied ones to find the correct determinant.
!Alternatively, we could pick orbitals at random and redraw until we find an allowed
!one. This would probably be preferable for larger systems.
        IF (tCycleOrbs) THEN
! METHOD 1 (Run though all orbitals in symmetry class with needed spin to find allowed one out of NExcit)
! ==========================

!Choose the unoccupied orbital to exite to
            r = genrand_real2_dSFMT()
            ChosenUnocc = INT(NExcit * r) + 1

!Now run through all allowed orbitals until we find this one.
            IF (tNoSymGenRandExcits) THEN
                nOrbs = nBasis / 2
            ELSE
                nOrbs = OrbClassCount(SymIndex)
            end if

            z = 0     !z is the counter for the number of allowed unoccupied orbitals we have gone through
            do i = 0, nOrbs - 1
!Find the spin orbital index. SymLabelCounts has the index of the state for the given symmetry.
                IF (tNoSymGenRandExcits) THEN
                    Orb = (2 * (i + 1)) - (iSpn - 1)
                ELSE
                    Orb = SymLabelList2(SymLabelCounts2(1, SymIndex) + i)
                end if

!Find out if the orbital is in the determinant.
                IF (.not. (BTEST(ILUT((Orb - 1) / bits_n_int), MOD(Orb - 1, bits_n_int)))) THEN
!The orbital is not found in the original determinant - increment z
                    z = z + 1
                    IF (z == ChosenUnocc) THEN
!We have got to the determinant that we want to pick.
                        EXIT
                    end if
                end if

            end do

!We now have our final orbitals. i=nI(Eleci). a=Orb.
            IF (z /= ChosenUnocc) THEN
                call stop_all(t_r, "Could not find allowed unoccupied orbital &
                                   &to excite to.")
            end if

        ELSE
! METHOD 2 (Keep drawing orbitals from the desired symmetry and spin until we find one unoccupied)
! =========================

            IF (tNoSymGenRandExcits) THEN
                nOrbs = nBasis / 2
            ELSE
                nOrbs = OrbClassCount(SymIndex)
            end if
            Attempts = 0
            do while (.true.)

!Draw randomly from the set of orbitals
                r = genrand_real2_dSFMT()
                ChosenUnocc = INT(nOrbs * r)
                IF (tNoSymGenRandExcits) THEN
                    Orb = (2 * (ChosenUnocc + 1)) - (iSpn - 1)
                ELSE
                    Orb = SymLabelList2(SymLabelCounts2(1, SymIndex) + ChosenUnocc)
                end if

!Find out if orbital is in nI or not. Accept if it isn't in it.
                IF (.not. (BTEST(ILUT((Orb - 1) / bits_n_int), MOD(Orb - 1, bits_n_int)))) THEN
!Orbital not in nI. Accept.
                    EXIT
                end if

                IF (Attempts > 250) THEN
                    write(stdout, *) "Cannot find single excitation unoccupied orbital after 250 attempts..."
                    write(stdout, *) "Desired symmetry of unoccupied orbital = ", ElecSym
                    write(stdout, *) "Number of orbitals (of correct spin) in symmetry = ", nOrbs
                    write(stdout, *) "Number of orbitals to legitimatly pick = ", NExcit
                    call write_det(6, nI, .true.)
                    call stop_all(t_r, "Cannot find single excitation &
                                   &unoccupied orbital after 250 attempts...")
                end if
                Attempts = Attempts + 1

            end do

        end if

        ! Construct the new determinant, excitation matrix and parity
        call make_single(nI, nJ, eleci, orb, ExcitMat, tParity)

#ifdef DEBUG_
        ! These are useful (but O[N]) operations to test the determinant
        ! generated. If there are any problems with then excitations, I
        ! recommend uncommenting these tests to check the results.
        if (.not. SymAllowedExcit(nI, nJ, 1, ExcitMat)) &
            call stop_all(t_r, 'Generated excitation invalid')
#endif

!Now we need to find the probability of creating this excitation.
!This is: P_single x P(i) x P(a|i) x N/(N-ElecsWNoExcits)
        pGen = (1 - pDoubNew) / (REAL((NExcit * (NEl - ElecsWNoExcits)), dp))

    END SUBROUTINE CreateSingleExcit

    subroutine construct_class_counts(nI, CCOcc, CCUnocc)

        ! Return two arrays of length ScratchSize, containing information on
        ! the number of orbitals - occupied and unoccupied - in each symmetry.
        !
        ! The arrays are indexed via the indices returned by ClassCountInd
        ! n.b. this is O[nel], so we should store this if we can.

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: CCOcc(ScratchSize), CCUnocc(ScratchSize)

        integer :: ind_alpha, ind_beta, i, ind

        CCOcc = 0
        CCUnocc = OrbClassCount

        if (tNoSymGenRandExcits .or. t_new_real_space_hubbard) then
            ind_alpha = ClassCountInd(1, 0, 0)
            ind_beta = ClassCountInd(2, 0, 0)
            CCOcc(ind_alpha) = nOccAlpha
            CCOcc(ind_beta) = nOccBeta
            CCUnocc(ind_alpha) = CCUnocc(ind_alpha) - nOccAlpha
            CCUnocc(ind_beta) = CCUnocc(ind_beta) - nOccBeta
        else
            do i = 1, nel
                ind = ClasscountInd(nI(i))

                CCOcc(ind) = CCOcc(ind) + 1
                CCUnocc(ind) = CCUnocc(ind) - 1
            end do
        end if

    end subroutine

    subroutine calc_pgen_symrandexcit2(nI, ex, ic, ClassCount2, &
                                       ClassCountUnocc2, pDoub, pGen)

        ! This routine will calculate the PGen between two connected
        ! determinants, nI and nJ which are IC excitations of each other, using
        ! the unbiased scheme.
        !
        ! Only the excitation matrix is needed (1,*) are the i,j orbs, and
        ! (2,*) are the a,b orbs. This is the prob of generating nJ FROM nI,
        ! not the other way round. Passed in is also the ClassCount2 arrays for
        ! nI, and the probability of picking a double.
        !
        ! A word of warning: The routine does not check that the determinants
        ! are indeed connected, and may well return a non-zero probability even
        ! if they arent. Therefore, make sure that they are at most double
        ! excitations of each other.
        !
        ! nI is the determinant from which the excitation comes from.

        integer, intent(in) :: nI(nel), ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize)
        integer, intent(in) :: ClassCountUnocc2(ScratchSize)
        real(dp), intent(in) :: pDoub
        real(dp), intent(out) :: pGen

        integer :: ForbiddenOrbs, symA, symB, sumMl, MlA, MlB, Elec1Ml
        integer :: ElecsWNoExcits, NExcitOtherWay, OrbI, OrbJ, iSpn, NExcitA
        integer :: NExcitB, ElecSym, orbA, orbB

        IF (IC == 1) THEN

!First, we need to find out if there are any electrons which have no possible excitations.
!This is because these will need to be redrawn and so
!will affect the probabilities.
            CALL CheckIfSingleExcits(ElecsWNoExcits, ClassCount2, ClassCountUnocc2, nI)

            IF (tNoSymGenRandExcits) THEN
!Find symmetry of chosen electron
                ElecSym = 0
            ELSE
                ElecSym = SpinOrbSymLabel((Ex(1, 1)))
                IF (tFixLz) THEN
                    Elec1Ml = G1(Ex(1, 1))%Ml
                end if
            end if

            IF (G1(Ex(1, 1))%Ms == 1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                NExcitA = ClassCountUnocc2(ClassCountInd(1, ElecSym, Elec1Ml))
            ELSE
!Beta orbital
                NExcitA = ClassCountUnocc2(ClassCountInd(2, ElecSym, Elec1Ml))
            end if

!Now we need to find the probability of creating this excitation.
!This is: P_single x P(i) x P(a|i) x N/(N-ElecsWNoExcits)
!Prob of generating a single is 1-pDoub
!            pGen=(1.0_dp-pDoub)*(1.0_dp/real(NEl,dp))*(1.0_dp/real(NExcitA,dp))*((real(NEl,dp))/(real((NEl-ElecsWNoExcits),dp)))
            pGen = (1 - pDoub) / (REAL((NExcitA * (NEl - ElecsWNoExcits)), dp))

        ELSE
!Prob of generating a double excitation.
!Find the I and J orbitals
            OrbI = Ex(1, 1)
            OrbJ = Ex(1, 2)
            OrbA = Ex(2, 1)
            OrbB = Ex(2, 2)
!            write(stdout,*) OrbI,OrbJ,OrbA,OrbB

!Find the spin-product of the occupied pair
            IF ((G1(OrbI)%Ms) * (G1(OrbJ)%Ms) == -1) THEN
!We have an alpha beta pair of electrons.
                iSpn = 2
!NExcit is the number of allowed unoccupied orbitals to choose
                NExcitA = nBasis - NEl
            ELSE
                IF (G1(OrbI)%Ms == 1) THEN
!We have an alpha alpha pair of electrons.
                    iSpn = 3
                    NExcitA = (nBasis / 2) - nOccAlpha     !This is the number of unocc alpha spinorbs
                ELSE
!We have a beta beta pair of electrons.
                    iSpn = 1
                    NExcitA = (nBasis / 2) - nOccBeta      !This is the number of unocc beta spinorbs
                end if
            end if

            IF (tFixLz) THEN
                SumMl = G1(OrbI)%Ml + G1(OrbJ)%Ml
                MlA = G1(OrbA)%Ml
                MlB = SumMl - MlA
            ELSE
                SumMl = 0
                MlA = 0
                MlB = 0
            end if
            IF (tNoSymGenRandExcits) THEN
                ElecSym = 0
                CALL FindNumForbiddenOrbsNoSym(ForbiddenOrbs, ClassCountUnocc2, iSpn)
                SymA = 0
                SymB = 0
            ELSE
!Calculate the symmetry product of the occupied orbital pair
                ElecSym = RandExcitSymLabelProd(SpinOrbSymLabel(OrbI), SpinOrbSymLabel(OrbJ))
!                ElecSym=INT(IEOR(G1(OrbI)%Sym%S,G1(OrbJ)%Sym%S),4)
!This will calculate the A orbitals which will have no B pair
                CALL FindNumForbiddenOrbs(ForbiddenOrbs, ClassCountUnocc2, ElecSym, iSpn, SumMl)
!Need to find the symmetries of the unoccupied A and B orbitals.
                SymA = SpinOrbSymLabel(OrbA)
                SymB = RandExcitSymLabelProd(SymInvLabel(SymA), ElecSym)
!                SymA=INT(G1(OrbA)%Sym%S,4)
!                SymB=IEOR(SymA,ElecSym)
            end if
!            write(stdout,*) "Check: ",ForbiddenOrbs,ElecSym,iSpn,SymA,SymB,OrbA,OrbB

!We want to calculate the number of possible B's given the symmetry and spin it has to be since we have already picked A.
!We have calculated in NExcit the number of orbitals available for B given A,
!but we also need to know the number of orbitals to choose from for A IF
!we had picked B first.
            IF (iSpn == 2) THEN
!If iSpn=2, then we want to find a spinorbital of the opposite spin of SpinOrbA
                IF ((G1(OrbA)%Ms) == -1) THEN
!We have already picked a beta orbital, so now we want to pick an alpha orbital. Find out how many of these there are.
                    NExcitB = ClassCountUnocc2(ClassCountInd(1, SymB, MlB))
                    NExcitOtherWay = ClassCountUnocc2(ClassCountInd(2, SymA, MlA))
                ELSE
!Want to pick an beta orbital.
                    NExcitB = ClassCountUnocc2(ClassCountInd(2, SymB, MlB))
                    NExcitOtherWay = ClassCountUnocc2(ClassCountInd(1, SymA, MlA))
                end if
            else if (iSpn == 1) THEN
!Definitely want a beta orbital
                NExcitB = ClassCountUnocc2(ClassCountInd(2, SymB, MlB))
                NExcitOtherWay = ClassCountUnocc2(ClassCountInd(2, SymA, MlA))
            ELSE
!Definitely want an alpha orbital
                NExcitB = ClassCountUnocc2(ClassCountInd(1, SymB, MlB))
                NExcitOtherWay = ClassCountUnocc2(ClassCountInd(1, SymA, MlA))
            end if

            IF ((iSpn /= 2) .and. (SymA == SymB) .and. (MlB == MlA)) THEN
!In this case, we need to check that we do not pick the same orbital as OrbA.
!If we do this, then we need to redraw.
!Only when ElecSym=0 will the classes of a and b be the same, and the spins will
!be different if iSpn=2, so this is the only possibility of a clash.
                NExcitB = NExcitB - 1     !Subtract 1 from the number of possible orbitals since we cannot choose orbital A.
                NExcitOtherWay = NExcitOtherWay - 1     !The same goes for the probabilities the other way round.
            end if

           pGen = pDoub * ((1.0_dp / real(NExcitB, dp)) + &
               (1.0_dp / real(NExcitOtherWay, dp))) / (REAL((ElecPairs * (NExcitA - ForbiddenOrbs)), dp))

        end if

    end subroutine

    !***********************  BIASED EXCITATION GENERATION ROUTINES *************************!

    SUBROUTINE CreateSingleExcitBiased(nI, nJ, iLut, ExcitMat, tParity, ElecsWNoExcits, nParts, WSign, Tau, iCreate)
        INTEGER :: ClassCount2(ScratchSize), i, Attempts, OrbA
        INTEGER :: ClassCountUnocc2(ScratchSize), Ind
        INTEGER :: ElecsWNoExcits, nParts, iCreate, nI(NEl), nJ(NEl)
        REAL(dp) :: WSign
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: ExcitMat(2, maxExcit), SpawnOrb(nBasis), Eleci, ElecSym, NExcit, VecInd, ispn, EndSymState, j
        real(dp) :: Tau, SpawnProb(nBasis), NormProb, r, rat
        LOGICAL :: tParity
        HElement_t(dp) :: rh
        debug_function_name("CreateSingleExcitBiased")

!First, we need to do an O[N] operation to find the number of occupied alpha electrons, number of occupied beta electrons
!and number of occupied electrons of each symmetry class and spin. This is similar to the ClassCount array.
!This has the format (Spn,sym), where Spin=1,2 corresponding to alpha and beta.
!For molecular systems, sym runs from 0 to 7. This is NOT general and should be made so using SymLabels.
!This could be stored to save doing this multiple times, but shouldn't be too costly an operation.
        CALL construct_class_counts(nI, ClassCount2, ClassCountUnocc2)

!We need to find out if there are any electrons which have no possible excitations.
!This is because these will need to be redrawn and so will affect the probabilities.
        ElecsWNoExcits = 0

!Need to look for forbidden electrons through all the irreps.
        do i = 0, nSymLabels - 1
!Run through all labels
            IF ((ClassCount2(ClassCountInd(1, i, 0)) /= 0) .and. (ClassCountUnocc2(ClassCountInd(1, i, 0)) == 0)) THEN
!If there are alpha electrons in this class with no possible unoccupied alpha orbitals
!in the same class, these alpha electrons have no single excitations.
                ElecsWNoExcits = ElecsWNoExcits + ClassCount2(ClassCountInd(1, i, 0))
            end if
            IF ((ClassCount2(ClassCountInd(2, i, 0)) /= 0) .and. (ClassCountUnocc2(ClassCountInd(2, i, 0)) == 0)) THEN
                ElecsWNoExcits = ElecsWNoExcits + ClassCount2(ClassCountInd(2, i, 0))
            end if
        end do

        IF (ElecsWNoExcits == NEl) THEN
!There are no single excitations from this determinant at all.
!Then we will create a double excitation instead.
            RETURN
        end if

!We want to pick an occupied electron - Prob = 1/(N-ElecsWNoExcits)
        Attempts = 0
        do while (.true.)

            Attempts = Attempts + 1

!Choose an electron randomly...
            r = genrand_real2_dSFMT()
            Eleci = INT(NEl * r) + 1

!Find symmetry of chosen electron
            ElecSym = INT((G1(nI(Eleci))%Sym%S), 4)

            IF (G1(nI(Eleci))%Ms == 1) THEN
!Alpha orbital - see how many single excitations there are from this electron...
                NExcit = ClassCountUnocc2(ClassCountInd(1, ElecSym, 0))
                ispn = 1  !Alpha ispn=1, Beta ispn=2
            ELSE
!Beta orbital
                NExcit = ClassCountUnocc2(ClassCountInd(2, ElecSym, 0))
                ispn = 2
            end if

            IF (NExcit /= 0) EXIT    !Have found electron with allowed excitations

            IF (Attempts > 250) THEN
                write(stdout, *) "Cannot find single excitation from electrons after 250 attempts..."
                call write_det(6, nI, .true.)
                write(stdout, *) "***"
                CALL Stop_All("CreateSingleExcit", "Cannot find single excitation from electrons after 250 attempts...")
            end if

        end do
        ExcitMat(1, 1) = nI(Eleci)

!Now we want to run through all sym+spin allowed excitations of the chosen electron, and determine their matrix elements
!To run just through the states of the required symmetry we want to use SymLabelCounts.

!We also want to take into account spin. We want the spin of the chosen unoccupied
!orbital to be the same as the chosen occupied orbital.
!Run over all possible a orbitals
        Ind = ClassCountInd(ispn, ElecSym, 0)
        EndSymState = SymLabelCounts2(1, Ind) + SymLabelCounts2(2, Ind) - 1

        VecInd = 1
        NormProb = 0.0_dp

        do j = SymLabelCounts2(1, Ind), EndSymState

!We want to look through all beta orbitals
            OrbA = SymLabelList2(j)     !This is the spin orbital chosen for a

            IF (BTEST(ILUT((OrbA - 1) / bits_n_int), MOD((OrbA - 1), bits_n_int))) THEN
!Orbital is in nI...not an unoccupied orbital
                CYCLE
            end if

!Now we want to find the information about this excitation
            ExcitMat(2, 1) = OrbA
            rh = sltcnd_excit(nI, SingleExc_t(ExcitMat(:, 1)), .false.)

            SpawnProb(VecInd) = abs(REAL(rh, dp))
            SpawnOrb(VecInd) = OrbA
            NormProb = NormProb + SpawnProb(VecInd)
            VecInd = VecInd + 1
            IF (VecInd >= nBasis) THEN
                CALL Stop_All("CreateSingleExcitBiased", "Finding too many virtual pairs...")
            end if

        end do
        IF ((VecInd - 1) /= NExcit) THEN
            CALL Stop_All("CreateSingleExcitBiased", "Wrong number of excitations found from chosen electron")
        end if
        IF (VecInd == 1) THEN
            CALL Stop_All("CreateSingleExcitBiased", "No excitations found from electron, but some should exist")
        end if

!We now have to find out how many children to spawn, based on the value of normprob.
        rat = Tau * NormProb * REAL((NEl - ElecsWNoExcits) * nParts, dp) / (1.0_dp - PDoubNew)
        iCreate = INT(rat)
        rat = rat - REAL(iCreate, dp)
        r = genrand_real2_dSFMT()
        IF (rat > r) THEN
!Child is created
            iCreate = iCreate + 1
        end if

        if (tGUGA) then
            call stop_all("CreateSingleExcitBiased", &
                          "modify get_helement for GUGA")
        end if
        IF (iCreate > 0) THEN
!We want to spawn particles. This only question now is where. Run
!through the ab pairs again and choose based on the SpawnProb element.
            r = genrand_real2_dSFMT()
            r = r * NormProb

            i = 0
            do while (r > 0.0_dp)
                i = i + 1
                r = r - SpawnProb(i)
            end do
            IF (i > VecInd - 1) THEN
                CALL Stop_All("CreateSingleExcitBiased", "Chosen virtual does not correspond to allowed orbital")
            end if

            OrbA = SpawnOrb(i)

            ! Construct the new determinant, excitation matrix and parity
            call make_single(nI, nJ, eleci, orbA, ExcitMat, tParity)

!These are useful (but O[N]) operations to test the determinant generated. If there are any problems with then
!excitations, I recommend uncommenting these tests to check the results.
         ASSERT(SymAllowedExcit(nI, nJ, 1, ExcitMat))

!Once we have the definitive determinant, we also want to find out what sign the particles we want to create are.
!iCreate is initially positive, so its sign can change depending on the sign of the connection and of the parent particle(s)
            rh = get_helement(nI, nJ, 1, ExcitMat, tParity)

            IF (WSign > 0.0_dp) THEN
                !Parent particle is positive
                IF (real(rh, dp) > 0.0_dp) THEN
                    iCreate = -iCreate     !-ve walker created
                end if
            ELSE
                IF (real(rh, dp) < 0.0_dp) THEN
                    iCreate = -iCreate    !-ve walkers created
                end if
            end if

        end if

    END SUBROUTINE CreateSingleExcitBiased

    SUBROUTINE CreateDoubExcitBiased(nI, nJ, iLut, ExcitMat, tParity, nParts, WSign, Tau, iCreate)
        INTEGER :: nI(NEl), nJ(NEl), ExcitMat(2, maxExcit), iCreate, iSpn, OrbA, OrbB, SymProduct
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: Elec1Ind, Elec2Ind, nParts, SumMl
        REAL(dp) :: WSign
        HElement_t(dp) :: rh
        LOGICAL :: tParity
        real(dp) :: Tau

!First, we need to pick an unbiased distinct electron pair.
!These have symmetry product SymProduct, and spin pair iSpn = 1=beta/beta; 2=alpha/beta; 3=alpha/alpha
!The probability for doing this is 1/ElecPairs.
        CALL PickElecPair(nI, Elec1Ind, Elec2Ind, SymProduct, iSpn, SumMl, -1)

!This routine runs through all distinct ab pairs for the chosen ij and stochastically chooses how many particles to create.
!If spawning wants to occur, then it runs through the list again and chooses a pair, which it returns.
        CALL CalcAllab(nI, iLut, Elec1Ind, Elec2Ind, SymProduct, iSpn, OrbA, OrbB, nParts, iCreate, Tau)

        if (tGUGA) then
            call stop_all("CreateDoubExcitBiased", &
                          "modify get_helement for GUGA")
        end if

!We now know that we want to create iCreate particles, from orbitals nI(Elec1/2Ind) -> OrbA + OrbB.
        IF (iCreate > 0) THEN
            call make_double(nI, nJ, elec1ind, elec2ind, orbA, orbB, &
                             ExcitMat, tParity)

!Once we have the definitive determinant, we also want to find out what sign the particles we want to create are.
!iCreate is initially positive, so its sign can change depending on the sign of the connection and of the parent particle(s)
            rh = get_helement(nI, nJ, 2, ExcitMat, tParity)

            IF (WSign > 0) THEN
                !Parent particle is positive
                IF (real(rh, dp) > 0.0_dp) THEN
                    iCreate = -iCreate     !-ve walker created
                end if
            ELSE
                IF (real(rh, dp) < 0.0_dp) THEN
                    iCreate = -iCreate    !-ve walkers created
                end if
            end if

        end if

    END SUBROUTINE CreateDoubExcitBiased

    SUBROUTINE CalcAllab(nI, ILUT, Elec1Ind, Elec2Ind, SymProduct, iSpn, OrbA, OrbB, nParts, iCreate, Tau)
        INTEGER :: nI(NEl), Elec1Ind, Elec2Ind, SymProduct, iSpn, OrbA, OrbB, iCreate
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: SpatOrbi, SpatOrbj, Spini, Spinj, i, aspn, bspn, SymA, SymB, SpatOrba, EndSymState, VecInd
        real(dp) :: Tau, SpawnProb(MaxABPairs), NormProb, rat, r
        INTEGER :: SpawnOrbs(2, MaxABPairs), j, nParts, SpinIndex, Ind
        HElement_t(dp) :: HEl

!We want the spatial orbital number for the ij pair (Elec1Ind is the index in nI).
!Later, we'll have to use GTID for UHF.
        SpatOrbi = ((nI(Elec1Ind) - 1) / 2) + 1
        SpatOrbj = ((nI(Elec2Ind) - 1) / 2) + 1
        Spini = G1(nI(Elec1Ind))%Ms
        Spinj = G1(nI(Elec2Ind))%Ms
        VecInd = 1
        NormProb = 0.0_dp

        do i = 1, nBasis
!Run through all a orbitals
            IF (mod(i, 2) == 0) THEN
!We have an alpha spin...
                aspn = 1
                IF (iSpn == 1) THEN
!ij is an beta/beta pair, so we only want to pick beta a orbitals.
                    CYCLE
                end if
            ELSE
!We have a beta spin...
                aspn = -1
                IF (iSpn == 3) THEN
!ij is a alpha/alpha pair, so we only want to pick alpha a orbitals.
                    CYCLE
                end if
            end if

!We also want to check that the a orbital we have picked is not already in the determinant we are exciting from.
            IF (BTEST(ILUT((i - 1) / bits_n_int), MOD((i - 1), bits_n_int))) THEN
!Orbital is in nI...do not make a pair from this.
                CYCLE
            end if

!Now we have to run over all b's which can go with a. We only want unique pairs, so we constrain b to be at least a+1.
!We only want to run over symmetry and spin allowed b's though.
!Find the required symmetry of b.
            SymA = INT(G1(i)%Sym%S, 4)
            SymB = IEOR(SymA, SymProduct)
            SpatOrba = ((i - 1) / 2) + 1

!We also want to take into account spin.
            IF (ispn == 1) THEN
                bspn = -1  !Want beta spin b orbitals
                SpinIndex = 2
            else if (ispn == 3) THEN
                bspn = 1  !Want alpha spin b orbitals
                SpinIndex = 1
            ELSE
!ij pair is an alpha/beta spin pair, therefore b wants to be of opposite spin to a.
                IF (aspn == -1) THEN
!a orbital is a beta orbital, therefore we want b to be an alpha orbital.
                    bspn = 1
                    SpinIndex = 1
                ELSE
                    bspn = -1
                    SpinIndex = 2
                end if
            end if

!To run just through the states of the required symmetry we want to use SymLabelCounts.
!            StartSymState=SymLabelCounts(1,SymB+1)
            Ind = ClassCountInd(SpinIndex, SymB, 0)
            EndSymState = SymLabelCounts2(1, Ind) + SymLabelCounts2(2, Ind) - 1

!Run over all possible b orbitals
            do j = SymLabelCounts2(1, Ind), EndSymState

!                IF(bspn.eq.-1) THEN
!                    OrbB=(2*SymLabelList(j))-1     !This is the spin orbital chosen for b
!                ELSE
!                    OrbB=(2*SymLabelList(j))
!                end if

                OrbB = SymLabelList2(j)     !This is the spin orbital chosen for b

                IF (OrbB <= i) THEN
!Since we only want unique ab pairs, ensure that b > a.
                    CYCLE
                end if

                IF (BTEST(ILUT((OrbB - 1) / bits_n_int), MOD((OrbB - 1), bits_n_int))) THEN
!Orbital is in nI...do not make a pair from this.
                    CYCLE
                end if

!We have now found an allowed ab pair to go with the ij pair chosen previously - record its stats.
                IF (Spini == aspn .and. Spinj == bspn) THEN
                    Hel = get_umat_el(SpatOrbi, SpatOrbj, SpatOrba, j)
                ELSE
                    Hel = (0.0_dp)
                end if
                IF (Spini == bspn .and. Spinj == aspn) THEN
                    Hel = Hel - get_umat_el(SpatOrbi, SpatOrbj, j, SpatOrba)
                end if

                SpawnProb(VecInd) = abs(REAL(Hel, dp))
                SpawnOrbs(1, VecInd) = i
                SpawnOrbs(2, VecInd) = OrbB
                NormProb = NormProb + SpawnProb(VecInd)
                VecInd = VecInd + 1
                IF (VecInd >= MaxABPairs) THEN
                    CALL Stop_All("CalcAllab", "Finding too many ab pairs...")
                end if

            end do
        end do

        VecInd = VecInd - 1     !This now indicates the total number of ab pairs we have found for the chosen ij.

        IF (VecInd == 0) THEN
            CALL Stop_All("CalcAllab", "No ab pairs found for the chosen ij")
        end if

!We now have to find out how many children to spawn, based on the value of normprob.
        rat = Tau * NormProb * REAL(ElecPairs * nParts, dp) / PDoubNew
        iCreate = INT(rat)
        rat = rat - REAL(iCreate, dp)
        r = genrand_real2_dSFMT()
        IF (rat > r) THEN
!Child is created
            iCreate = iCreate + 1
        end if

        IF (iCreate > 0) THEN
!We want to spawn particles. This only question now is where. Run through the
!ab pairs again and choose based on the SpawnProb element.
            r = genrand_real2_dSFMT()
            r = r * NormProb

            i = 0
            do while (r > 0.0_dp)
                i = i + 1
                r = r - SpawnProb(i)
            end do
            IF (i > VecInd) THEN
                CALL Stop_All("CalcAllab", "Chosen ab pair does not correspond to allowed pair")
            end if

            OrbA = SpawnOrbs(1, i)
            OrbB = SpawnOrbs(2, i)

        end if

    END SUBROUTINE CalcAllab

! For a given ij pair in the UEG or Hubbard model, this generates ab as a double excitation efficiently.
! This takes into account the momentum conservation rule, i.e. that kb=ki+ki-ka(+G).
    SUBROUTINE CreateDoubExcitLattice(nI, iLutnI, nJ, tParity, ExcitMat, pGen, Elec1Ind, Elec2Ind, iSpn, part_type)

        INTEGER :: i, nI(NEl), nJ(NEl), Elec1Ind, Elec2Ind, iSpn, kb_ms
        INTEGER(KIND=n_int) :: iLutnI(0:NIfTot)
        INTEGER :: ChosenUnocc, Hole1BasisNum, Hole2BasisNum, ki(3), kj(3), ka(3), kb(3), ExcitMat(2, maxExcit), iSpinIndex, TestEnergyB
        LOGICAL :: tAllowedExcit, tParity
        real(dp) :: r, pGen, pAIJ
        integer, intent(in), optional :: part_type

        integer :: elecs(2), loc, temp_part_type

        ! [W.D.] just a workaround for now, after i rehaul all the lattice
        ! excitation generators, this will be more optimized
        if (present(part_type)) then
            temp_part_type = part_type
        else
            temp_part_type = 1
        end if

        ! This chooses an a of the correct spin, excluding occupied orbitals
        ! This currently allows b orbitals to be created that are disallowed
        hole2basisnum = 0
        DO
            r = genrand_real2_dSFMT()
            ! Choose a
            IF (iSpn == 2) THEN ! alpha/beta combination
                ChosenUnocc = INT(nBasis * r) + 1
            ELSE ! alpha/alpha, beta/beta
                ChosenUnocc = 2 * (INT(nBasis / 2 * r) + 1) & ! 2*(a number between 1 and nbasis/2) gives the alpha spin
                & - (1 - (iSpn / 3)) ! alpha numbered even, iSpn/3 returns 1 for alpha/alpha, 0 for beta/beta
            end if
            ! Check a isn't occupied
            IF (.not. (BTEST(iLutnI((ChosenUnocc - 1) / bits_n_int), MOD(ChosenUnocc - 1, bits_n_int)))) THEN
                !Orbital not in nI. Accept.
                EXIT
            end if
        end do

        Hole1BasisNum = ChosenUnocc
        IF ((Hole1BasisNum <= 0) .or. (Hole1BasisNum > nBasis)) THEN
            CALL Stop_All("CreateDoubExcitLattice", "Incorrect basis function generated")
        end if

        !=============================================
        if (tUEG2) then

            ! kb is now uniquely defined
            ki = kvec(nI(Elec1Ind), 1:3)
            kj = kvec(nI(Elec2Ind), 1:3)
            ka = kvec(Hole1BasisNum, 1:3)
            kb = ki + kj - ka

            ! Find the spin of b
            IF (iSpn == 2) THEN ! alpha/beta required, therefore b has to be opposite spin to a
                kb_ms = 1 - ((G1(Hole1BasisNum)%Ms + 1) / 2) * 2
            ELSE ! b is the same spin as a
                kb_ms = G1(Hole1BasisNum)%Ms
            end if

            ! Is kb allowed by the size of the space?
            tAllowedExcit = .true.
            TestEnergyB = kb(1)**2 + kb(2)**2 + kb(3)**2
            IF (tOrbECutoff .and. (TestEnergyB > OrbECutoff)) tAllowedExcit = .false.
            IF (.not. tAllowedExcit) THEN
                nJ(1) = 0
                RETURN
            end if
            iSpinIndex = (kb_ms + 1) / 2 + 1
            Hole2BasisNum = kPointToBasisFn(kb(1), kb(2), kb(3), iSpinIndex)

            IF (Hole2BasisNum == -1 .or. Hole1BasisNum == Hole2BasisNum) THEN
                nJ(1) = 0
                RETURN
            end if

            ! Is b occupied?
            IF (BTEST(iLutnI((Hole2BasisNum - 1) / bits_n_int), MOD(Hole2BasisNum - 1, bits_n_int))) THEN
                !Orbital is in nI. Reject.
                tAllowedExcit = .false.
            end if

            IF (.not. tAllowedExcit) THEN
                nJ(1) = 0
                RETURN
            end if

            ! Check that the correct kb has been found -- can be commented out later
            DO i = 1, 3
                IF ((kvec(nI(Elec2Ind), i) + kvec(nI(Elec1Ind), i) &
                     - kvec(Hole1BasisNum, i) - kvec(Hole2BasisNum, i)) /= 0) THEN
                    write(stdout, *) "Tried to excite "
                    write(stdout, *) "ki ", ki
                    write(stdout, *) "kj ", kj
                    write(stdout, *) "ka ", ka
                    write(stdout, *) "kb should be ", kb
                    write(stdout, *) "but found as ", kvec(Hole2BasisNum, 1:3)
                    CALL Stop_All("CreateDoubExcitLattice", "Wrong b found")
                end if
            end do

            ! Find the new determinant
            call make_double(nI, nJ, elec1ind, elec2ind, Hole1BasisNum, &
                             Hole2BasisNum, ExcitMat, tParity)

            !Calculate generation probabilities
            IF (iSpn == 2) THEN
                pAIJ = 1.0_dp / (nBasis - Nel)
            else if (iSpn == 1) THEN
                pAIJ = 1.0_dp / (nBasis / 2 - nOccBeta)
            ELSE
                !iSpn = 3
                pAIJ = 1.0_dp / (nBasis / 2 - nOccAlpha)
            end if
            ! Note, p(b|ij)=p(a|ij) for this system

            pGen = 2.0_dp / (NEl * (NEl - 1)) * 2.0_dp * pAIJ

            return
        end if
        !=============================================

        ! kb is now uniquely defined
        ki = G1(nI(Elec1Ind))%k
        kj = G1(nI(Elec2Ind))%k
        ka = G1(Hole1BasisNum)%k
        kb = ki + kj - ka

        ! Find the spin of b
        IF (iSpn == 2) THEN ! alpha/beta required, therefore b has to be opposite spin to a
            kb_ms = 1 - ((G1(Hole1BasisNum)%Ms + 1) / 2) * 2
        ELSE ! b is the same spin as a
            kb_ms = G1(Hole1BasisNum)%Ms
        end if

        IF (tHub) THEN

            ! Re-map back into cell
            CALL MomPbcSym(kb, nBasisMax)

            ! This is the look-up table method of finding the kb orbital
            iSpinIndex = (kb_ms + 1) / 2 + 1
            ! [W.D.]
            ! make the access to kpoint same in hubbard and ueg
!             Hole2BasisNum=kPointToBasisFn(kb(1),kb(2),1,iSpinIndex)
            Hole2BasisNum = kPointToBasisFn(kb(1), kb(2), kb(3), iSpinIndex)
!Hole2BasisNum will be -1 if that orb is frozen

            ! Is b occupied?
            IF (Hole2BasisNum == -1 .or. BTEST(iLutnI((Hole2BasisNum - 1) / bits_n_int), MOD(Hole2BasisNum - 1, bits_n_int))) THEN
                !Orbital is in nI. Reject.
                nJ(1) = 0
                RETURN
            end if

            ! Reject if a=b.
            IF (Hole1BasisNum == Hole2BasisNum) THEN
                nJ(1) = 0
                RETURN
            end if

        end if

        IF (tUEG) THEN

            ! Is kb allowed by the size of the space?
            ! Currently only applies when NMAXX etc. are set by the CELL keyword
            ! Not sure what happens when an energy cutoff is set
            tAllowedExcit = .true.
            IF (ABS(kb(1)) > NMAXX) tAllowedExcit = .false.
            IF (ABS(kb(2)) > NMAXY) tAllowedExcit = .false.
            IF (ABS(kb(3)) > NMAXZ) tAllowedExcit = .false.
            TestEnergyB = kb(1)**2 + kb(2)**2 + kb(3)**2
            IF (tOrbECutoff .and. (TestEnergyB > OrbECutoff)) tAllowedExcit = .false.
            IF (.not. tAllowedExcit) THEN
                nJ(1) = 0
                RETURN
            end if

            ! Which orbital has momentum kb?
            !IF (iSpn.eq.2) THEN
            !    Hole2BasisNum=2*((NMAXZ*2+1)*(NMAXY*2+1)*(kb(1)+NMAXX)+(NMAXZ*2+1)*
            !(kb(2)+NMAXY)+(kb(3)+NMAXZ)+1)-(1+G1(Hole1BasisNum)%Ms)/2
            !ELSE
            !    Hole2BasisNum=2*((NMAXZ*2+1)*(NMAXY*2+1)*(kb(1)+NMAXX)+(NMAXZ*2+1)*
            !(kb(2)+NMAXY)+(kb(3)+NMAXZ)+1)-1+iSpn/3
            ! end if

            iSpinIndex = (kb_ms + 1) / 2 + 1
            Hole2BasisNum = kPointToBasisFn(kb(1), kb(2), kb(3), iSpinIndex)

            IF (Hole2BasisNum == -1 .or. Hole1BasisNum == Hole2BasisNum) THEN
                nJ(1) = 0
                RETURN
            end if

            ! Is b occupied?
            IF (BTEST(iLutnI((Hole2BasisNum - 1) / bits_n_int), MOD(Hole2BasisNum - 1, bits_n_int))) THEN
                !Orbital is in nI. Reject.
                tAllowedExcit = .false.
            end if

            IF (.not. tAllowedExcit) THEN
                nJ(1) = 0
                RETURN
            end if

            ! Check that the correct kb has been found -- can be commented out later
            DO i = 1, 3
                IF ((G1(nI(Elec2Ind))%k(i) + G1(nI(Elec1Ind))%k(i) - G1(Hole1BasisNum)%k(i) - G1(Hole2BasisNum)%k(i)) /= 0) THEN
                    write(stdout, *) "Tried to excite "
                    write(stdout, *) "ki ", ki
                    write(stdout, *) "kj ", kj
                    write(stdout, *) "ka ", ka
                    write(stdout, *) "kb should be ", kb
                    write(stdout, *) "but found as ", G1(Hole2BasisNum)%k
                    CALL Stop_All("CreateDoubExcitLattice", "Wrong b found")
                end if
            end do

        end if

        ! Find the new determinant
        call make_double(nI, nJ, elec1ind, elec2ind, Hole1BasisNum, &
                         Hole2BasisNum, ExcitMat, tParity)

        IF (tHub) THEN
            ! Debug to test the resultant determinant
            ! i think i also have to check here if kpoint symmetry is turned
            ! on... or is it generally turned on in the hubbard case??
            IF (.not. (IsMomentumAllowed(nJ)) .and. tKPntSym) THEN
                CALL Stop_All("CreateDoubExcitLattice", "Incorrect kb generated -- momentum not conserved")
            end if
        end if

        !Calculate generation probabilities
        IF (iSpn == 2) THEN
            ! otherwise take the pAIJ set above
            pAIJ = 1.0_dp / (nBasis - Nel)
        else if (iSpn == 1) THEN
            pAIJ = 1.0_dp / (nBasis / 2 - nOccBeta)
        ELSE
            !iSpn = 3
            pAIJ = 1.0_dp / (nBasis / 2 - nOccAlpha)
        end if
        ! Note, p(b|ij)=p(a|ij) for this system
        if (tUEG) then
            pGen = 2.0_dp / (NEl * (NEl - 1)) * 2.0_dp * pAIJ
        else ! i.e. if hubbard model, use modified probabilities
            ! hubbard model can't spawn alpha/alpha and beta/beta type excitations
            pGen = 1.0_dp / (nOccAlpha * nOccBeta) * 2.0_dp * pAIJ
        end if

    END SUBROUTINE CreateDoubExcitLattice

    ! This creates a general excitation for the UEG and Hubbard model. Each ij pair is
    !picked with equal probability, and then
    ! each ab is picked for that ij pair with equal probability.
    ! For UEG there is a more sophisticated algorithm that allows more ijab choices to be
    !rejected before going back to the main
    ! code.
    SUBROUTINE CreateExcitLattice(nI, iLutnI, nJ, tParity, ExcitMat, pGen, part_type)

        INTEGER :: i, j ! Loop variables
        INTEGER :: Elec1, Elec2
        INTEGER :: nI(NEl), nJ(NEl), Elec1Ind, Elec2Ind, ExcitMat(2, maxExcit), iSpn
        INTEGER(KIND=n_int) :: iLutnI(0:NIfTot)
        INTEGER :: ki(3), kj(3), kTrial(3), iElecInExcitRange, iExcludedKFromElec1, iAllowedExcites
        INTEGER :: KaXLowerLimit, KaXUpperLimit, KaXRange, KaYLowerLimit, KaYUpperLimit, KaYRange, KaZLowerLimit, KaZUpperLimit, KaZRange
        LOGICAL :: tParity, tDoubleCount, tExtraPoint
        real(dp) :: r(2), pGen, pAIJ
        integer, intent(in), optional :: part_type
        INTEGER, ALLOCATABLE :: Excludedk(:, :)
        character(*), parameter :: this_routine = "CreateExcitLattice"
        integer :: elecs(2), src(2), sym_prod, sum_ml
        real(dp) :: pgen_back

        integer :: temp_part_type
!        CALL PickElecPair(nI,Elec1Ind,Elec2Ind,SymProduct,iSpn,SumMl,-1)

        ! Completely random ordering of electrons is important when considering ij->ab ij/->ba.
        !This affects pgens for alpha/beta pairs.
        kaxrange = 0
        ielecinexcitrange = 0
        kazrange = 0
        kayrange = 0

        if (present(part_type)) then
            temp_part_type = part_type
        else
            temp_part_type = 1
        end if

        DO
            r(1) = genrand_real2_dSFMT()
            Elec1 = INT(r(1) * NEl + 1)
            DO
                r(2) = genrand_real2_dSFMT()
                Elec2 = INT(r(2) * NEl + 1)
                IF (Elec2 /= Elec1) EXIT
            end do
            Elec1Ind = Elec1
            Elec2Ind = Elec2
            IF ((G1(nI(Elec1Ind))%Ms == -1) .and. (G1(nI(Elec2Ind))%Ms == -1)) THEN
                iSpn = 1
            ELSE
                IF ((G1(nI(Elec1Ind))%Ms == 1) .and. (G1(nI(Elec2Ind))%Ms == 1)) THEN
                    iSpn = 3
                ELSE
                    iSpn = 2
                end if
            end if
            IF ((tHub .and. iSpn == 2) .or. (tUEG)) EXIT ! alpha/beta pairs are the only pairs generated for the hubbard model
        end do

        IF (tNoFailAb) THEN ! pGen is calculated first because there might be no excitations available for this ij pair

            IF (tHub) CALL Stop_All("CreateExcitLattice", "This doesn't work with Hubbard Model")

            ! Find the upper and lower ranges of kx, ky and kz from the point of view of electron i
            ki = G1(nI(Elec1Ind))%k
            kj = G1(nI(Elec2Ind))%k
            !================================================
            if (tUEG2) then
                ki = kvec(nI(Elec1Ind), 1:3)
                kj = kvec(nI(Elec2Ind), 1:3)
            end if
            !================================================
            KaXLowerLimit = MAX(-NMAXX, ki(1) - (NMAXX - kj(1)))
            KaXUpperLimit = MIN(NMAXX, ki(1) + (NMAXX + kj(1)))
            KaXRange = KaXUpperLimit - KaXLowerLimit + 1
            KaYLowerLimit = MAX(-NMAXY, ki(2) - (NMAXY - kj(2)))
            KaYUpperLimit = MIN(NMAXY, ki(2) + (NMAXY + kj(2)))
            KaYRange = KaYUpperLimit - KaYLowerLimit + 1
            KaZLowerLimit = MAX(-NMAXZ, ki(3) - (NMAXZ - kj(3)))
            KaZUpperLimit = MIN(NMAXZ, ki(3) + (NMAXZ + kj(3)))
            KaZRange = KaZUpperLimit - KaZLowerLimit + 1

            ! How many disallowed excitations are there due to electrons 'blocking' i->a or j->b excitations
            iElecInExcitRange = 0
            tExtraPoint = .false.
            allocate(Excludedk(NEl, 3))
            ! Is a = b allowed by momentum and disallowed by spin symmetry?
            IF ((G1(nI(Elec1Ind))%Ms == G1(nI(Elec2Ind))%Ms) .and. &
            &       (MOD(ki(1) - kj(1), 2) == 0) .and. &
            &       (MOD(ki(2) - kj(2), 2) == 0) .and. &
            &       (MOD(ki(3) - kj(3), 2) == 0)) THEN ! This is the disallowed double-excitation to the same orbital
                iElecInExcitRange = iElecInExcitRange + 1
                Excludedk(iElecInExcitRange, :) = (ki + kj) / 2 ! Integer division okay because we're already checked for mod2=0
                tExtraPoint = .true.
            end if
            ! How many disallowed i->a are there?
            DO i = 1, NEl
                IF (G1(nI(i))%Ms /= G1(nI(Elec1Ind))%Ms) CYCLE
                kTrial = G1(nI(i))%k
                !================================================
                if (tUEG2) then
                    kTrial = kvec(nI(i), 1:3)
                end if
                !================================================
                IF (kTrial(1) < KaXLowerLimit) CYCLE
                IF (kTrial(1) > KaXUpperLimit) CYCLE
                IF (kTrial(2) < KaYLowerLimit) CYCLE
                IF (kTrial(2) > KaYUpperLimit) CYCLE
                IF (kTrial(3) < KaZLowerLimit) CYCLE
                IF (kTrial(3) > KaZUpperLimit) CYCLE
                IF (tExtraPoint) THEN
                    IF ((kTrial(1) == Excludedk(1, 1))    &
                    & .and. (kTrial(2) == Excludedk(1, 2))   &
                    & .and. (kTrial(3) == Excludedk(1, 3))) CYCLE
                end if
                iElecInExcitRange = iElecInExcitRange + 1
                Excludedk(iElecInExcitRange, :) = kTrial
            end do
            iExcludedKFromElec1 = iElecInExcitRange
            ! How many disallowed j->b are there, given that some k-points have already been elimated by i->a being disallowed?
            DO i = 1, NEl
                IF (G1(nI(i))%Ms /= G1(nI(Elec2Ind))%Ms) CYCLE
                kTrial = ki + kj - G1(nI(i))%k
                !================================================
                if (tUEG2) then
                    kTrial = ki + kj - kvec(nI(i), 1:3)
                end if
                !================================================
                IF (kTrial(1) < KaXLowerLimit) CYCLE
                IF (kTrial(1) > KaXUpperLimit) CYCLE
                IF (kTrial(2) < KaYLowerLimit) CYCLE
                IF (kTrial(2) > KaYUpperLimit) CYCLE
                IF (kTrial(3) < KaZLowerLimit) CYCLE
                IF (kTrial(3) > KaZUpperLimit) CYCLE
                ! Need to check for this k-point already having been eliminated
                ! by the previous loop over electrons
                tDoubleCount = .false.
                DO j = 1, iExcludedKFromElec1
                    IF ((Excludedk(j, 1) == kTrial(1)) .and. (Excludedk(j, 2) == kTrial(2)) .and. Excludedk(j, 3) == kTrial(3)) then
                        tDoubleCount = .true.
                    end if
                end do
                IF (.not. tDoubleCount) iElecInExcitRange = iElecInExcitRange + 1
            end do

            DEallocate(Excludedk)

            ! If there are no excitations for this ij pair then this subroutine must exit
            iAllowedExcites = KaXRange * KaYRange * KaZRange - iElecInExcitRange
            IF (iAllowedExcites == 0) THEN
                !    write(stdout,*) "No allowed excitations from this ij pair"
                nJ(1) = 0
                pGen = 1.0_dp
                RETURN
            end if
        end if

        DO i = 1, 10000
            IF (i == 10000) THEN ! Arbitrary termination of the loop to prevent hanging
                ! due to no excitations being allowed from a certain ij pair
                write(stdout, *) "nI:", nI
                write(stdout, *) "i & j", nI(Elec1Ind), nI(Elec2Ind)
                write(stdout, *) "Allowed Excitations", iAllowedExcites
                CALL Stop_All("CreateExcitLattice", "Failure to generate a valid excitation from an electron pair combination")
            end if
            CALL CreateDoubExcitLattice(nI, iLutnI, nJ, tParity, ExcitMat, pGen, Elec1Ind, Elec2Ind, iSpn, temp_part_type)
            IF (.not. tNoFailAb) then
                ! For an invalid excitation we still need to have a valid value for pGen
                if (nJ(1) == 0) pGen = 1.0_dp
                RETURN
            end if
            IF (nJ(1) /= 0) EXIT ! i.e. if we are using the NoFail algorithm only exit on successful nJ(1)!=0
        end do

        ! ***This part of the code is only used if tNoFailAb is ON***
        ! Else the pgen used is from CreateDoubExcitLattice

        ! Now calculate pgen
        pAIJ = 1.0_dp / (KaXRange * KaYRange * KaZRange - iElecInExcitRange)
        ! pBIJ is zero for this kind of excitation generator for antiparallel spins
        ! but is equal to pAIJ for parallel spins.
        IF (G1(nI(Elec1Ind))%Ms /= G1(nI(Elec2Ind))%Ms) THEN
            pGen = 2.0_dp / (NEl * (NEl - 1)) * pAIJ ! Spins not equal
        ELSE
            pGen = 2.0_dp / (NEl * (NEl - 1)) * 2.0 * pAIJ ! Spins equal
        end if

        IF (pAIJ <= 0.0_dp) CALL Stop_All("CreateExcitLattice", "pAIJ is less than 0")

    END SUBROUTINE CreateExcitLattice

    !As with CalcNonUniPgens:
    !This routine will calculate the PGen between two connected determinants, nI and nJ which
    !are IC excitations of each other, using the unbiased scheme.
    !Only the excitation matrix is needed (1,*) are the i,j orbs, and (2,*) are the a,b orbs
    !This routine does it for the lattice models: UEG and hubbard model
    SUBROUTINE CalcPGenLattice(Ex, pGen)
        use back_spawn, only: get_ispn

        INTEGER :: Ex(2, 2), iSpin, jSpin
        real(dp) :: pGen, pAIJ
        character(*), parameter :: this_routine = "CalcPGenLattice"

        ASSERT(.not. tNoFailAb)

        ! can i make this easier?
        iSpin = get_ispn(get_src(ex))
        if (iSpin == 1) then
            pAIJ = 1.0_dp / (nBasis / 2 - nOccBeta)
        else if (iSpin == 2) then
            pAIJ = 1.0_dp / (nBasis - Nel)
        else if (iSpin == 3) then
            pAIJ = 1.0_dp / (nBasis / 2 - nOccAlpha)
        end if

        ! Note, p(b|ij)=p(a|ij) for this system
        ! there is a small difference here and in the excitation generator
        ! part. is is only multiplied by 2 if it is a parallel excitation
        ! in the excitation generator.. but not here anymore..
        if (tUEG) then
            pGen = 2.0_dp / (NEl * (NEl - 1)) * 2.0_dp * pAIJ
        else ! i.e. if hubbard model, use modified probabilities
            ! hubbard model can't spawn alpha/alpha and beta/beta type excitations
            pGen = 1.0_dp / (nOccAlpha * nOccBeta) * 2.0_dp * pAIJ
        end if

    END SUBROUTINE CalcPGenLattice

    FUNCTION IsMomentumAllowed(nJ)

        use sym_mod, only: mompbcsym
        use SystemData, only: kvec, tUEG2, tAllSymSectors

        LOGICAL :: IsMomentumAllowed ! Returns whether the determinant is momentum allowed for
        ! UEG and Hubbard models
        ! Compares the total k from a determinant nI with kTotal
        INTEGER :: nJ(NEl), kx, ky, kz, ktrial(3), i

        if (tAllSymSectors) then
            IsMomentumAllowed = .true.
            return
        end if

        IsMomentumAllowed = .false.

        !====================================================
        if (tUEG2) then
            ! The momentum constraint for UEG: every determinant must have a total momentum
            ! which is equal.
            kx = 0
            ky = 0
            kz = 0
            do i = 1, NEl
                kx = kx + kvec(nJ(i), 1)
                ky = ky + kvec(nJ(i), 2)
                kz = kz + kvec(nJ(i), 3)
            end do
            IF (kx == kTotal(1) .and. ky == kTotal(2) .and. kz == kTotal(3)) THEN
                IsMomentumAllowed = .true.
            end if
            return
        end if
        !====================================================
        ! The momentum constraint for UEG: every determinant must have a total momentum
        ! which is equal.
        IF (tUEG) THEN
            kx = 0
            ky = 0
            kz = 0
            do i = 1, NEl
                kx = kx + G1(nJ(i))%k(1)
                ky = ky + G1(nJ(i))%k(2)
                kz = kz + G1(nJ(i))%k(3)
            end do
            IF (kx == kTotal(1) .and. ky == kTotal(2) .and. kz == kTotal(3)) THEN
                IsMomentumAllowed = .true.
            end if
        end if

        ! The momentum constraint from Hubbard model: every determinant must have a total momentum
        ! which is equal to within a reciprocal lattice vector.
        IF (tHub) THEN
            kx = 0
            ky = 0
            kz = 0
            do i = 1, NEl
                kx = kx + G1(nJ(i))%k(1)
                ky = ky + G1(nJ(i))%k(2)
                kz = kz + G1(nJ(i))%k(3)
            end do
            ktrial = (/kx, ky, 0/)
            CALL MomPbcSym(ktrial, nBasisMax) ! This re-maps the total momentum under PBCs: equivalent to this being equal to
            ! a value to within a reciproval lattice vector.
            IF (ktrial(1) == kTotal(1) .and. ktrial(2) == kTotal(2)) THEN
                IsMomentumAllowed = .true.
            else
                print *, "ERRONEOUS MOMENTUM: ", ktrial, ktotal
            end if
        end if

    END FUNCTION IsMomentumAllowed

    ! Initialise or clean the excitation generation storage.
    subroutine init_excit_gen_store(store)

        type(excit_gen_store_type), intent(inout) :: store

        call clean_excit_gen_store(store)

        store%tFilled = .false.

        allocate(store%ClassCountOcc(ScratchSize1))
        allocate(store%ClassCountUnocc(ScratchSize2))
        allocate(store%scratch3(ScratchSize3))

        if (tBuildOccVirtList) then
            allocate(store%occ_list(nel, ScratchSize1))
            allocate(store%virt_list(maxval(OrbClassCount), Scratchsize1))
        end if

        if (tBuildSpinSepLists) then
            allocate(store%nI_alpha(nel))
            allocate(store%nI_beta(nel))
            allocate(store%nI_alpha_inds(nel))
            allocate(store%nI_beta_inds(nel))
            !store%nel_alpha = 0
        end if

        if (tSpinProject) then
            allocate(store%dorder_i(nel))
            allocate(store%dorder_j(nel))
        end if

        if (t_pcpp_excitgen) allocate(store%elec_map(nel))

    end subroutine

    subroutine clean_excit_gen_store(store)

        type(excit_gen_store_type), intent(inout) :: store

        if (associated(store%ClassCountOcc)) then
            deallocate(store%ClassCountOcc)
            nullify (store%ClassCountOcc)
        end if
        if (associated(store%ClassCountUnocc)) then
            deallocate(store%ClassCountUnocc)
            nullify (store%ClassCountUnocc)
        end if
        if (associated(store%scratch3)) then
            deallocate(store%scratch3)
            nullify (store%scratch3)
        end if
        if (associated(store%occ_list)) then
            deallocate(store%occ_list)
            nullify (store%occ_list)
        end if
        if (associated(store%virt_list)) then
            deallocate(store%virt_list)
            nullify (store%virt_list)
        end if
        if (associated(store%dorder_i)) then
            deallocate(store%dorder_i)
            nullify (store%dorder_i)
        end if
        if (associated(store%dorder_j)) then
            deallocate(store%dorder_j)
            nullify (store%dorder_j)
        end if
        if (associated(store%elec_map)) then
            deallocate(store%elec_map)
            nullify (store%elec_map)
        end if

    end subroutine

END MODULE GenRandSymExcitNUMod

!Sometimes (especially UHF orbitals), the symmetry routines will not set up the orbitals correctly.
!Therefore, this routine will set up symlabellist and symlabelcounts
!to be cast in terms of spin orbitals, and the symrandexcit2 routines will use these arrays.
SUBROUTINE SpinOrbSymSetup()
    use SymExcitDataMod, only: ScratchSize, ScratchSize1, ScratchSize2, &
                               SymLabelList2, SymLabelCounts2, kTotal, &
                               OrbClassCount, kPointToBasisFn, sym_label_list_spat
    use SymExcitDataMod, only: SpinOrbSymLabel, SymInvLabel, &
                               SymTableLabels, KPntInvSymOrb
    use GenRandSymExcitNUMod, only: ClassCountInd
    use SymData, only: nSymLabels, SymClasses, SymLabels
    use SystemData, only: NMAXZ, tFixLz, iMaxLz, nBasis, tUEG, tKPntSym, G1, tHub, nBasisMax, NEl
    use SystemData, only: Symmetry, tHPHF, tSpn, tISKFuncs, Arr, tNoSymGenRandExcits, Elecpairs
    use SystemData, only: MaxABPairs, tUEG2, kvec, t_k_space_hubbard
    use Determinants, only: FDet
    use umatcache, only: gtID
    use sym_mod, only: mompbcsym, SymProd, writesym
    use constants, only: dp, stdout
    use lattice_mod, only: lat

    IMPLICIT NONE

    INTEGER :: i, j, k, SymInd, Lab
    INTEGER :: Spin, ierr, OrbSym, InvSym
    real(dp) :: OrbEnergy
    INTEGER, ALLOCATABLE :: Temp(:)
    ! These are for the hubbard and UEG model look-up table
    INTEGER :: kmaxX, kmaxY, kminX, kminY, kminZ, kmaxz, iSpinIndex, ktrial(3)
    type(Symmetry) :: SymProduct, SymI, SymJ
    character(len=*), parameter :: this_routine = 'SpinOrbSymSetup'
    integer :: sym0

    ElecPairs = (NEl * (NEl - 1)) / 2
    MaxABPairs = (nBasis * (nBasis - 1) / 2)
!    write(stdout,*) "SETTING UP SYMMETRY!!",nBasis,elecpairs

    IF (tFixLz) THEN
!Calculate the upper array bound for the ClassCount2 arrays. This will be dependant on the number of symmetries needed.
        ScratchSize = 2 * nSymLabels * (2 * iMaxLz + 1)
    ELSE
        ScratchSize = 2 * nSymLabels    !For k-point sym, this can be large
    end if

    IF (tNoSymGenRandExcits .or. tUEG) THEN
        ScratchSize = 2
    end if

    ! We only need the first 2 scratch arrays
    ScratchSize1 = ScratchSize
    ScratchSize2 = ScratchSize

    !Create SpinOrbSymLabel array.
    !This array will return a number between 0 and nSymLabels-1.
    !For molecular systems, this IS the character of the irrep
    !For k-point systems, this is an arbitrary label, and is equal to the standard label - 1.
    !This is chosen so that the indexing works with the rest of the excitation generators.
    if (allocated(SpinOrbSymLabel)) deallocate(SpinOrbSymLabel)
    allocate(SpinOrbSymLabel(nBasis))
    do i = 1, nBasis
        if (tNoSymGenRandExcits .or. tUEG) then
            SpinOrbSymLabel(i) = 0
        else if (tKPntSym .or. tHub) then
            SpinOrbSymLabel(i) = SymClasses(((i + 1) / 2)) - 1        !This ensures that the symmetry labels go from 0 -> nSymLabels-1
        else
            SpinOrbSymLabel(i) = INT(G1(i)%Sym%S, 4)
        end if
    end do
#ifdef DEBUG_
    write(stdout, *) "SpinOrbSymLabel: "
    do i = 1, nBasis
        write(stdout, *) i, SpinOrbSymLabel(i)
    end do
#endif

    if (tKPntSym) then
        allocate(SymTableLabels(0:nSymLabels - 1, 0:nSymLabels - 1))
        SymTableLabels(:, :) = -9000    !To make it easier to track bugs
        do i = 0, nSymLabels - 1
            do j = 0, i
                SymI = SymLabels(i + 1)        !Convert to the other symlabel convention to use SymLabels -
                !TODO: I will fix this to make them consistent when working (ghb24)!
                SymJ = SymLabels(j + 1)

                SymProduct = SymProd(SymI, SymJ)
                !Now, we need to find the label according to this symmetry!
                !Run through all symmetries to make working (this could be far more efficient, but its only once, so sod it...
                do Lab = 1, nSymLabels
                    if (SymLabels(Lab)%S == SymProduct%S) then
                        EXIT
                    end if
                end do
                if (Lab == nSymLabels + 1) then
                    call stop_all("SpinOrbSymSetup", "Cannot find symmetry label")
                end if
                SymTableLabels(i, j) = Lab - 1
                SymTableLabels(j, i) = Lab - 1
            end do
        end do
#ifdef DEBUG_
        write(stdout, *) "SymTable:"
        do i = 0, nSymLabels - 1
            do j = 0, nSymLabels - 1
                write(stdout, "(I6)", advance='no') SymTableLabels(i, j)
            end do
            write(stdout, *) ""
        end do
#endif

    end if
!SymInvLabel takes the label (0 -> nSymLabels-1) of a spin orbital, and returns the inverse symmetry label, suitable for
!use in ClassCountInd.
    if (allocated(SymInvLabel)) deallocate(SymInvLabel)
    allocate(SymInvLabel(0:nSymLabels - 1))
    SymInvLabel = -999

    ! Dongxia changes the gamma point away from center.
    ! SDS: Provide a default sym0 for cases where this doesn't apply
    sym0 = 0
    do i = 1, nsymlabels
        if (symlabels(i)%s == 0) sym0 = i - 1
    end do

    do i = 0, nSymLabels - 1
        if (tKPntSym) then
            ! Change the sym label back to the representation used by the rest
            ! of the code, use SymConjTab, then change back to other rep of
            ! labels SymConjTab only works when all irreps are self-inverse.
            ! Therefore, instead, we will calculate the inverses by just
            ! finding the symmetry which will give A1.
            do j = 0, nSymLabels - 1
                ! Run through all labels to find what gives totally symmetric
                ! rep
                if (SymTableLabels(i, j) == sym0) then
                    if (SymInvLabel(i) /= -999) then
                        write(stdout, *) "SymLabel: ", i
                        call stop_all(this_routine, &
                                      "Multiple inverse irreps found - error")
                    end if
                    ! This is the inverse
                    SymInvLabel(i) = j
                end if
            end do
            if (SymInvLabel(i) == -999) then
                write(stdout, *) "SymLabel: ", i
                call stop_all(this_routine, "No inverse symmetry found - error")
            end if
        else
            ! If not using k-point sym, then we are self-inverse
            SymInvLabel(i) = i
        end if
    end do
#ifdef DEBUG_
    write(stdout, *) "SymInvLabel: "
    do i = 0, nSymLabels - 1
        write(stdout, *) i, SymInvLabel(i)
    end do
#endif

    if (tISKFuncs) then
        write(stdout, *) "Setting up inverse orbital lookup for use with ISK functions..."
        if (.not. tKPntSym) call stop_all(this_routine, "Cannot use ISK funcs without KPoint symmetry")
        if (tSpn) call stop_all(this_routine, "Cannot use ISK on open shell systems")
        if (tHPHF) call stop_all(this_routine, "HPHF is not yet working with ISK")
        allocate(KPntInvSymOrb(1:nBasis), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Allocation error")
        do i = 1, nBasis
            !Find k-inverse orbital for each spin orbital.
            !If the orbital is a self inverse, then just lookup itself.
            OrbSym = SpinOrbSymLabel(i)
            InvSym = SymInvLabel(OrbSym)
            if (InvSym == OrbSym) then
                !Orbital is self-inverse
                KPntInvSymOrb(i) = i
                cycle
            end if
            !Run through all orbitals looking for inverse
            OrbEnergy = Arr(i, 2)  !Fock energy of orbital
            !Ensure that we get a same-spin orbital
            do j = 1, nBasis
                if ((SpinOrbSymLabel(j) == InvSym) .and. (mod(i, 2) == mod(j, 2))) then
                    !This orbital is the right symmetry - is it the inverse orbital? Check Energy.
                    if ((abs(OrbEnergy - Arr(j, 2))) < 1.0e-7_dp) then
                        !Assume that this is the inverse orbital.
                        KPntInvSymOrb(i) = j
                        exit
                    end if
                end if
            end do
            if (j > nBasis) then
                write(stdout, *) "Orbital: ", i
                write(stdout, *) "Symmetry label: ", OrbSym
                write(stdout, *) "Inverse Symmetry label: ", InvSym
                write(stdout, *) "Fock energy: ", Arr(i, 2)
                write(stdout, *) "SpinOrbSymLabel: "
                do k = 1, nBasis
                    write(stdout, *) k, SpinOrbSymLabel(k)
                end do
                write(stdout, *) "SymInvLavel: "
                do k = 0, nSymLabels - 1
                    write(stdout, *) k, SymInvLabel(k)
                end do
                call stop_all(this_routine, "Could not find inverse orbital pair for ISK setup.")
            end if
        end do
        do i = 1, nBasis
            if (KPntInvSymOrb(i) /= i) exit
        end do
        if (i > nBasis) then
            write(stdout, *) "!! All kpoints are self-inverse, i.e. at the Gamma point or BZ boundary !!"
            write(stdout, *) "This means that ISK functions cannot be constructed."
            write(stdout, *) "However, through correct rotation of orbitals, all orbitals should be " &
            & //"made real, and the code run in real mode (with tRotatedOrbsReal set)."
            write(stdout, *) "Alternatively, run again in complex mode without ISK functions."
            write(stdout, *) "If ISK functions are desired, the kpoint lattice must be shifted from this position."
!            call stop_all(this_routine,"All kpoints are self-inverse")
        end if
        write(stdout, *) "All inverse kpoint orbitals correctly assigned."
        write(stdout, *) "Orbital     Inverse Orbital"
        do i = 1, nBasis
            write(stdout, *) i, KPntInvSymOrb(i)
        end do
    end if

!SymLabelList2 and SymLabelCounts2 are now organised differently, so that it is more efficient, and easier to add new symmetries.
!SymLabelCounts is of size (2,ScratchSize), where 1,x gives the index in SymlabelList2 where the orbitals of symmetry x start.
!SymLabelCounts(2,x) tells us the number of orbitals of spin & symmetry x there are.

!Therefore, if you want to run over all orbitals of a specific symmetry, you want to run over
!SymLabelList from SymLabelCounts(1,sym) to SymLabelCounts(1,sym)+SymLabelCounts(2,sym)-1

    allocate(SymLabelList2(nBasis))

    if (allocated(sym_label_list_spat)) deallocate(sym_label_list_spat)

    allocate(sym_label_list_spat(nBasis))
    allocate(SymLabelCounts2(2, ScratchSize))
    SymLabelList2(:) = 0          !Indices:   spin-orbital number
    SymLabelCounts2(:, :) = 0      !Indices:   index/Number , symmetry(inc. spin)
    allocate(Temp(ScratchSize))

    do j = 1, nBasis
        IF (G1(j)%Ms == 1) THEN
            Spin = 1
        ELSE
            Spin = 2
        end if
!        write(stdout,*) "BASIS FN ",j,G1(j)%Sym,SymClasses((j+1)/2)
        SymInd = ClassCountInd(Spin, SpinOrbSymLabel(j), G1(j)%Ml)
        SymLabelCounts2(2, SymInd) = SymLabelCounts2(2, SymInd) + 1
    end do
    SymLabelCounts2(1, 1) = 1
    do j = 2, ScratchSize
        SymLabelCounts2(1, j) = SymLabelCounts2(1, j - 1) + SymLabelCounts2(2, j - 1)
    end do
    Temp(:) = SymLabelCounts2(1, :)
    do j = 1, nBasis
        IF (G1(j)%Ms == 1) THEN
            Spin = 1
        ELSE
            Spin = 2
        end if
        SymInd = ClassCountInd(Spin, SpinOrbSymLabel(j), G1(j)%Ml)
        SymLabelList2(Temp(SymInd)) = j
        Temp(SymInd) = Temp(SymInd) + 1
    end do
    ! get also only the spatial orbitals
    sym_label_list_spat = gtID(SymLabelList2)

!    write(stdout,*) "SymLabelCounts2: ",SymLabelCounts2(1,:)
!    write(stdout,*) "SymLabelCounts2: ",SymLabelCounts2(2,:)
    Deallocate(Temp)

    allocate(OrbClassCount(ScratchSize))
    OrbClassCount(:) = 0

    IF (tNoSymGenRandExcits .or. tUEG) THEN
!All orbitals are in irrep 0
        OrbClassCount(ClassCountInd(1, 0, 0)) = (nBasis / 2)
        OrbClassCount(ClassCountInd(2, 0, 0)) = (nBasis / 2)
    ELSE
        do i = 1, nBasis
            IF (G1(i)%Ms == 1) THEN
!                write(stdout,*) "Index: ",ClassCountInd(1,SpinOrbSymLabel(i),G1(i)%Ml)
!                write(stdout,*) i,"SpinOrbSymLabel: ",SpinOrbSymLabel(i)
                OrbClassCount(ClassCountInd(1, SpinOrbSymLabel(i), G1(i)%Ml)) = &
                & OrbClassCount(ClassCountInd(1, SpinOrbSymLabel(i), G1(i)%Ml)) + 1
            ELSE
!                write(stdout,*) "Index: ",ClassCountInd(1,SpinOrbSymLabel(i),G1(i)%Ml)
!                write(stdout,*) i,"SpinOrbSymLabel: ",SpinOrbSymLabel(i)
                OrbClassCount(ClassCountInd(2, SpinOrbSymLabel(i), G1(i)%Ml)) = &
                & OrbClassCount(ClassCountInd(2, SpinOrbSymLabel(i), G1(i)%Ml)) + 1
            end if
        end do
    end if

!    write(stdout,*) "*******",OrbClassCount(:)

!        ELSE
!!SymLabelCounts(2,1:nSymLabels) gives the number of *states* in each symmetry class.
!!There are therefore equal number of alpha and beta orbitals in each state from which to calculate the unoccupied classcount.
!            do i=1,nSymLabels
!                OrbClassCount(ClassCountInd(1,i-1,0))=SymLabelCounts2(1,2,i)
!                OrbClassCount(ClassCountInd(2,i-1,0))=SymLabelCounts2(2,2,i)
!            end do
!        end if

    ! This makes a 3D lookup table kPointToBasisFn(kx,ky,1,ms_index) which gives the orbital number for a given kx, ky and ms_index
    IF (tHub .and. .not. (NMAXZ /= 0 .and. NMAXZ /= 1)) THEN
!        IF(NMAXZ.ne.0.and.NMAXZ.ne.1) CALL Stop_All("SpinOrbSymSetup","This routine doesn't work with non-2D Hubbard model")
        kmaxX = 0
        kminX = 0
        kmaxY = 0
        kminY = 0
        ! [W.D:] can we make this the same as in the UEG:
        kminZ = 0
        kmaxZ = 0
        do i = 1, nBasis ! In the hubbard model with tilted lattice boundary conditions, it's unobvious what the maximum values of
            ! kx and ky are, so this should be found
            IF (G1(i)%k(1) > kmaxX) kmaxX = G1(i)%k(1)
            IF (G1(i)%k(1) < kminX) kminX = G1(i)%k(1)
            IF (G1(i)%k(2) > kmaxY) kmaxY = G1(i)%k(2)
            IF (G1(i)%k(2) < kminY) kminY = G1(i)%k(2)
            if (G1(i)%k(3) > kmaxz) kmaxz = g1(i)%k(3)
            if (G1(i)%k(3) < kminz) kminz = g1(i)%k(3)
        end do
!         allocate(kPointToBasisFn(kminX:kmaxX,kminY:kmaxY,1,2))
        allocate(kPointToBasisFn(kminX:kmaxX, kminY:kmaxY, kminz:kmaxz, 2))
        kPointToBasisFn = -1 !Init to invalid
        do i = 1, nBasis
            iSpinIndex = (G1(i)%Ms + 1) / 2 + 1 ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
            kPointToBasisFn(G1(i)%k(1), G1(i)%k(2), G1(i)%k(3), iSpinIndex) = i
        end do
    end if

    !======================================================
    if (tUEG2) then
        kmaxX = 0
        kminX = 0
        kmaxY = 0
        kminY = 0
        kminZ = 0
        kmaxZ = 0

        do i = 1, nBasis
            IF (kvec(i, 1) > kmaxX) kmaxX = kvec(i, 1)
            IF (kvec(i, 1) < kminX) kminX = kvec(i, 1)
            IF (kvec(i, 2) > kmaxY) kmaxY = kvec(i, 2)
            IF (kvec(i, 2) < kminY) kminY = kvec(i, 2)
            IF (kvec(i, 3) > kmaxZ) kmaxZ = kvec(i, 3)
            IF (kvec(i, 3) < kminZ) kminZ = kvec(i, 3)
        end do

        allocate(kPointToBasisFn(kminX:kmaxX, kminY:kmaxY, kminZ:kmaxZ, 2))
        kPointToBasisFn = -1 !Init to invalid
        do i = 1, nBasis
            iSpinIndex = (G1(i)%Ms + 1) / 2 + 1 ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
            kPointToBasisFn(kvec(i, 1), kvec(i, 2), kvec(i, 3), iSpinIndex) = i
        end do

        kTotal(1) = 0
        kTotal(2) = 0
        kTotal(3) = 0
        do j = 1, NEl
            kTotal(1) = kTotal(1) + kvec(FDet(j), 1)
            kTotal(2) = kTotal(2) + kvec(FDet(j), 2)
            kTotal(3) = kTotal(3) + kvec(FDet(j), 3)
        end do
        write(stdout, *) "Total momentum is", kTotal

        return
    end if
    !======================================================

    ! This makes a 3D lookup table kPointToBasisFn(kx,ky,kz,ms_index) which
    !gives the orbital number for a given kx, ky, kz and ms_index
    IF (tUEG) THEN
        kmaxX = 0
        kminX = 0
        kmaxY = 0
        kminY = 0
        kminZ = 0
        kmaxZ = 0

        do i = 1, nBasis
            IF (G1(i)%k(1) > kmaxX) kmaxX = G1(i)%k(1)
            IF (G1(i)%k(1) < kminX) kminX = G1(i)%k(1)
            IF (G1(i)%k(2) > kmaxY) kmaxY = G1(i)%k(2)
            IF (G1(i)%k(2) < kminY) kminY = G1(i)%k(2)
            IF (G1(i)%k(3) > kmaxZ) kmaxZ = G1(i)%k(3)
            IF (G1(i)%k(3) < kminZ) kminZ = G1(i)%k(3)
        end do
        allocate(kPointToBasisFn(kminX:kmaxX, kminY:kmaxY, kminZ:kmaxZ, 2))
        kPointToBasisFn = -1 !Init to invalid
        do i = 1, nBasis
            iSpinIndex = (G1(i)%Ms + 1) / 2 + 1 ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
            kPointToBasisFn(G1(i)%k(1), G1(i)%k(2), G1(i)%k(3), iSpinIndex) = i
        end do
    end if

    IF (tUEG .or. tHUB) THEN
        kTotal = 0
        do j = 1, NEl
            if (t_k_space_hubbard) then
                kTotal = lat%add_k_vec(kTotal, G1(FDet(j))%k)
            else
                kTotal = kTotal + G1(FDet(j))%k
            end if
            ! just to be sure..
            ! not necessary anymore with new add_k_vec
!             if (t_k_space_hubbard) then
!                 ktotal = lat%map_k_vec(ktotal)
!             end if
        end do
        if (tHub .and. .not. t_k_space_hubbard) then
            ! is this turned off correctly?! check:
            ktrial = (/kTotal(1), kTotal(2), 0/)
            CALL MomPbcSym(ktrial, nBasisMax)
            kTotal(1) = ktrial(1)
            kTotal(2) = ktrial(2)
        end if
        write(stdout, *) "Total momentum is", kTotal
    end if

END SUBROUTINE SpinOrbSymSetup

logical function IsMomAllowedDet(nJ)
    use SystemData, only: nEl
    use FciMCData, only: HFSym
    implicit none
    logical :: IsMomAllowedDetAnyParent
    integer :: nJ(nEl)
    IsMomAllowedDet = IsMomAllowedDetAnyParent(nJ, HFSym%Sym)
end function IsMomAllowedDet

LOGICAL FUNCTION IsMomAllowedDetAnyParent(nJ, parentSym)
    use sym_mod
    use SystemData, only: G1, NEl, Symmetry, nBasisMax, BasisFN
!    use GenRandSymExcitNUMod , only : Counter
    IMPLICIT NONE
    Type(Symmetry) :: SYM1
    Type(Symmetry) :: parentSym
    Type(BasisFN) :: iSym
    INTEGER :: i, nJ(NEl), KPnt(3)

    IsMomAllowedDetAnyParent = .false.

    SYM1%S = 0
    do i = 1, NEl
        SYM1 = SYMPROD(SYM1, G1(nJ(i))%Sym)
    end do

    IF (SYM1%S /= parentSym%S) THEN
        write(stdout, *) "nJ: ", nJ(:)
        write(stdout, *) "parentSym,SYM1: ", parentSym%S, SYM1%S
!        write(stdout,*) "Counter: ",Counter
        CALL DecomposeAbelianSym(SYM1%S, KPnt)
        write(stdout, "(A,3I5)") "KPnt for nJ: ", KPnt(1), KPnt(2), KPnt(3)
        CALL Stop_All("IsMomAllowedDet", "Momentum forbidden excitation created1.")
    ELSE
        IsMomAllowedDetAnyParent = .true.
    end if

    CALL GETSYM(nJ, NEl, G1, nBasisMax, iSym)

    IF (iSym%Sym%S /= parentSym%S) THEN
        write(stdout, *) "nJ: ", nJ(:)
        write(stdout, *) "parentSym,SYM1: ", parentSym%S, iSym%Sym%S
!        write(stdout,*) "Counter: ",Counter
        CALL DecomposeAbelianSym(iSym%Sym%S, KPnt)
        write(stdout, "(A,3I5)") "KPnt for nJ: ", KPnt(1), KPnt(2), KPnt(3)
        CALL Stop_All("IsMomAllowedDet", "Momentum forbidden excitation created2.")
    ELSE
        IsMomAllowedDetAnyParent = .true.
    end if

END FUNCTION IsMomAllowedDetAnyParent

