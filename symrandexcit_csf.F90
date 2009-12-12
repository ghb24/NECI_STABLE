#include "macros.h"
! This is a random excitation generator for use with csfs.
! Generate using a normalised and calculable probability.
module GenRandSymExcitCSF
    use Systemdata, only: nel, NIftot, tNoSymGenRandExcits, G1, LMS, nbasis
    use SystemData, only: nbasismax, lztot, tFixLz, iMaxLz
    use SymExcitDataMod
    use SymData, only: TwoCycleSymGens, nSymLabels
    use csf, only: csf_orbital_mask, csf_test_bit, csf_apply_random_yama
    use csf, only: get_num_csfs, csf_apply_yama, csf_get_yamas, write_yama
    use mt95, only: genrand_real2
    use GenRandSymExcitNUMod, only: ClassCountInd
    use DetBitOps, only: EncodeBitDet, DecodeBitDet, is_canonical_ms_order
    use DetBitOps, only: shift_det_bit_singles_to_beta, count_open_orbs
    use Parallel
    implicit none

    ! Non-modularised functions (sigh)
    interface
        real*8 pure function choose(N,R)
            integer, intent(in) :: N,R
        end function
        logical function int_arr_eq (a, b, len)
            integer, intent(in), dimension(:) :: a, b
            integer, intent(in), optional :: len
        end function
    end interface
    
contains
    subroutine GenRandSymCSFExcit (nI, iLut, nJ, pSingle, pDouble, IC, &
                                   ExcitMat, exFlag, pGen, CCDbl, CCSgl, &
                                   CCUn, tFilled)
        integer, intent(in)    :: nI(nel), iLut(0:NIfTot), exFlag
        integer, intent(out)   :: nJ(nel), IC, ExcitMat(2,2)
        integer, intent(inout) :: CCDbl(ScratchSize) ! ClassCountDoubleOcc2
        integer, intent(inout) :: CCSgl(ScratchSize) ! ClassCountSingleOcc2
        integer, intent(inout) :: CCUn(ScratchSize)  ! ClassCountUnocc2
        integer :: CCDblS(ScratchSize/2), CCSglS(ScratchSize/2)
        integer :: CCUnS(ScratchSize/2) ! Only consider spacial terms
        logical, intent(inout) :: tFilled
        real*8,  intent(in)    :: pSingle, pDouble
        real*8,  intent(out)   :: pGen
        character(*), parameter   :: this_routine = 'GenRandSymExcitCSF'
        integer :: Attempts, nopen, ncsf
        logical :: bSingle

        ! Count the open shell electrons
        nopen = count_open_orbs(iLut) 

        ! If the array is not already populated, perform an O[N] operation to
        ! find the number of occupied alpha/beta electrons, and number of
        ! occupied e- of each symmetry class and spin.
        if (.not. tFilled) then
            if ((.not.TwoCycleSymgens) .and. (.not.tNoSymGenRandExcits)) then
                write(6,'("GenRandSymExcitCSF can only be used for molecular&
                          & systems")')
                write(6,'("This is because of difficulties with other &
                          &symmetries setup.")')
                write(6,'("If you want to use these excitation generators &
                          &then add NOSYMGEN to the input to ignor symmetry &
                          &while generating excitations.")')
                call flush(6)
                call stop_all(this_routine,"GenRandsymExcitCSF can only be &
                                  &used for molecular systems using symmetry")
            endif

            call ConstructClassCounts(nI, nel-nopen, CCDbl, CCSgl, CCUn)
            tFilled = .true.
        endif

        ! TODO: make this covered by tFilled as well...
        call ConstructClassCountsSpacial(nI, nel-nopen, CCDblS, CCSglS, CCUnS)

        ! Select type of excitation depending on ExcitFlag.
        select case (ExFlag)
            case (1)
                IC = 0
            case (2)
                IC = 1
            case (4)
                IC = 2
            case default
                ! TODO: Pick IC randomly, depending on pSingle/pDouble.
                !       Options depend on btest of exflag
                call stop_all (this_routine, "Unsupported excitation mode")
        end select

        ! Do what we need to do
        select case (IC)
        case (0)
            ! Change only the Yamanouchi symbol. If it is not possible to
            ! change it (i.e. if ncsf = 0,1) then return 0 determinant.
            nJ = nI
            call csf_apply_random_yama (nJ, nopen, real(LMS,8)/2, ncsf,.true.)
            if (ncsf < 2) then
                nJ(1) = 0
            else
                pGen = (1-pSingle-pDouble) / (ncsf - 1)
            endif
        case (1) ! Create single excitation
            call CSFCreateSingleExcit (nI, nJ, CCDbl, CCSgl, CCUn, iLut, &
                                    ExcitMat, nopen, pSingle, pGen)
        case (2) ! Create double excitation
            call CSFCreateDoubleExcit (nI, nJ, CCDblS, CCSglS, CCUnS, iLut, &
                                       ExcitMat, nopen, pDouble, pGen)
        endselect
    end subroutine

    subroutine CSFCreateDoubleExcit (nI, nJ, CCDblS, CCSglS, CCUnS, iLut, &
                                     ExcitMat, nopen, pDouble, pGen)
        use GenRandSymExcitNUMod, only: CreateDoubExcit
        use csf, only: iscsf
        use SymExcitDataMod, only: pDoubNew
        integer, intent(in)  :: nI(nel), iLut(0:NIfTot), nopen
        integer, intent(out) :: nJ(nel)
        integer, intent(inout) :: ExcitMat(2,2)
        integer, intent(in)  :: CCSglS (ScratchSize/2) ! ClassCountSingleOcc2
        integer, intent(in)  :: CCDblS (ScratchSize/2) ! ClassCountDoubleOcc2
        integer, intent(in)  :: CCUnS (ScratchSize/2)  ! ClassCountUnocc2
        integer :: CCSglDelta (ScratchSize/2)
        real*8,  intent(in)  :: pDouble
        real*8, intent(out) :: pGen
        character(*), parameter :: this_routine = 'CSFCreateDoubleExcit'
        integer :: nexcits, elecs(2), orbs(2,2), symProd, sumMl, sym(2), Ml(2)
        integer :: CCOcc (ScratchSize), iLutTmp(0:NIfTot), lnopen, nup, nEquiv
        integer :: ncsf, ListOpen(nopen), nSing, nDoub, nVac, nopen2
        integer :: orbsWNoPair
        integer, save :: excitcount = 0
        real*8 :: tmpR

        !do nup=1,46
        !    print*, nup, G1(nup)%Sym%S, G1(nup)%Ms, btest(G1(nup)%Sym%S, 0)
        !enddo
        !call stop_all (this_routine, 'listed')

        ! General working values.
        lnopen = nopen
        pGen = pDouble
        nSing = sum(CCSglS)
        nDoub = sum(CCDblS)
        nVac = sum(CCUnS)

        ! Pick an electron pair at random, and modulate pGen
        call CSFPickElecPair (nI, iLut, elecs, orbs(1,:), symProd, sumMl, &
                              nSing, nDoub, CCSglDelta, lnopen, pGen)

        ! Count the number of orbitals we can't excite to due to there being
        ! no available pair to make the symmetry product == symProd.
        orbsWNoPair = CSFCalcOrbsWNoSymmetryPair (CCDblS, CCSglS, CCUnS, &
                                                  CCSglDelta, symProd)

        ! If there are no orbitals we can excite to, then fail.
        ! TODO: Is this really < 2 (ie is 1 impossible, as pair valid?)
        if ((nSing + nVac) == orbsWNoPair) then
            nJ(1) = 0
            return
        endif

        ! Pick the first orbital, as one which has an allowed other excitation
        call CSFPickAOrb (nI, iLut, CCDblS, CCSglS, CCUnS, CCSglDelta, &
                          orbs, symProd, sumMl, orbs(2,1), sym, Ml, &
                          lnopen, nSing, nVac)

        ! Pick a second orbital, with symmetry such that the product gives
        ! symProd (conserves symmetry of original electron pair).
        call CSFPickBOrb (nI, iLut, CCDblS, CCSglS, CCUnS, CCSglDelta, &
                          orbs(1,:), symProd, sumMl, orbs(2,1), sym, &
                          Ml, orbs(2,2), lnopen)

        ! Set this up for the find_excit_det_general routine
        ExcitMat (1,1) = min(elecs(1), elecs(2))
        ExcitMat (2,1) = min(orbs(2,1), orbs(2,2))
        ExcitMat (1,2) = max(elecs(1), elecs(2))
        ExcitMat (2,2) = max(orbs(2,1), orbs(2,2))
        !print*, 'ExcitMat', iand(nI(ExcitMat(1,1)),csf_orbital_mask), &
        !iand(nI(ExcitMat(1,2)),csf_orbital_mask), Excitmat(2,1), ExcitMat(2,2)
        !print*, 'a,b', orbs(2,1), orbs(2,2)

        ! Generate the excited determinant. ncsf specifies that we are
        ! only generating one determinant, and returns the degeneracy of the
        ! Yamanouchi symbol.
        nJ = nI
        ncsf = 1
        call csf_find_excit_det_general (ExcitMat, nJ, iLut, nopen, 2, ncsf, &
                                         nopen2)

        ! Adjust the probability for the selection of orbitals and for the
        ! degeneracy of the Yamanouchi symbol
        pGen = pGen * CSFPickOrbsProb (nI, iLut, CCDblS, CCSglS, CCUnS, &
                                       CCSglDelta, orbsWNoPair, nSing, nVac,&
                                       orbs(2,:), sym, Ml)
        pGen = pGen / real(ncsf)

        !excitcount = excitcount + 1
        !if (excitcount == 6) call stop_all (this_routine,  "DONE")
        !print*, 'count', excitcount
        if (iand(NJ(nel), csf_orbital_mask) == 0)&
            call stop_all (this_routine,"alert")
    end subroutine

    integer function CSFCalcOrbsWNoSymmetryPair (CCDblS, CCSglS, CCUnS, &
                                  CCSglDelta, symProd)
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCDblS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: CCSglDelta (ScratchSize/2)
        integer, intent(in) :: symProd
        character(*), parameter :: this_routine = 'CSFPickAOrb'
        integer :: numB, ind, i, orbsWNoPair

        ! Loop over all beta symmetry classes and exclude orbitals with no
        ! pair which would give symProd as the resultant.
        orbsWNoPair = 0
        if (tFixLz) then
            call stop_all (this_routine, 'Not implemented yet')
            ! TODO: Understand tFixLz here
            !do i=-iMaxLz,iMaxLz
            !    mlA = sumMl - i
            !    ! TODO: Why i < mlA?
            !    if ((i < mlA) .and. (abs(mlA) < iMaxLz)) then
            !        do i=1,nSymLabels-1 ! Loop over the betas (spacial only)
            !            ind = CCindS (ibclr(ieor(i,symProd)), mlA)

            !        enddo
            !    endif
            !enddo
        else
            ! Loop over all sym indices for alpha elecs. If there are none, 
            ! then orbital A cannot be from the associated sym.
            do i=1,2*nSymlabels-1,2
                numB = CCUnS(i) + CCSglS(i) + CCSglDelta(i)
                ! If no B orbitals, then A not allowed.
                if (numB == 0) then
                    ! sym index for required A orbitals for this B
                    ind = CCindS (ieor(i,symProd), 0)
                    orbsWNoPair = orbsWNoPair + CCUnS(ind) + CCSglS(ind) + &
                                  CCSglDelta(ind)
                ! Cannot excite two electrons to one orbital.
                else if (numB == 1) then
                    ind = CCindS (ieor(i,symProd), 0)
                    if (ind == i) orbsWNoPair = orbsWNoPair + 1
                endif
            enddo
        endif

        CSFCalcOrbsWNoSymmetryPair = orbsWNoPair
    end function

    subroutine CSFPickAOrb (nI, iLut, CCDblS, CCSglS, CCUnS, CCSglDelta, &
                            orbs, symProd, sumMl, orbA, sym, Ml, &
                            nopen, nSing, nVac)
        integer, intent(in) :: nI(nel), iLut(0:NIfTot), symProd
        integer, intent(in) :: sumMl, orbs(2), nSing, nVac
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCDblS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: CCSglDelta (ScratchSize/2)
        integer, intent(out) :: orbA, sym(2), Ml(2)
        integer, intent(inout) :: nopen
        character(*), parameter :: this_routine = 'CSFPickAOrb'
        real*8 :: r
        logical :: bSingle
        integer :: orb, orb2, ind, i, orbsWNoPair, numB

        ! Draw orbitals randomly until we find one unoccupied
        do i=1,250
            call genrand_real2 (r)
            orb = 2 * int(r*nBasis/2) ! 0 based, beta

            ! If exciting to self, or to the other source, this is either
            ! invalid, or is a single (also disallowed).
            if ((orb+1 == orbs(1)) .or. (orb+2 == orbs(1)) .or. &
                (orb+1 == orbs(2)) .or. (orb+2 == orbs(2))) cycle

            ! Is this occupied
            bSingle = .false.
            if (btest(ilut(orb/32), mod(orb,32))) then
                ! Is it a single, otherwise a double, so try again.
                if (.not.(btest(iLut((orb+1)/32), mod(orb+1,32)))) then
                    orb = orb + 1
                    nopen = nopen - 1 ! Creating a double
                else
                    cycle
                endif
            else
                bSingle = .true. ! Exciting into vacant -> create single
                nopen = nopen + 1
            endif
            orb = orb + 1 ! Now 1 based index

            ! Check that the Ml of another chosen orbital would be valid
            Ml = 0
            if (tFixLz) then
                Ml(1) = G1(orb)%Ml
                Ml(2) = SumMl - Ml(1)
                if (abs(Ml(2)) > iMaxLz) cycle
            endif

            ! Is there a symmetry allowed B orbital to match this A?
            sym(1) = int(G1(orb)%Sym%S,4)
            sym(2) = ieor(sym(1), symProd)
            ind = CCIndS(sym(2), Ml(2))

            if (CCUnS(ind) > 0) exit
            if (CCSglS(ind)+CCSglDelta(ind) > 0) exit
        enddo

        if (i > 250) then
            write(6,'("Cannot find an unoccupied orbital for a double")')
            write(6,'("excitation after 250 attempts.")')
            write(6,'("Desired symmetry of orbital pair =",i3)') &
                symProd
            call writedet(6,nI,nel,.true.)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a double excitation after 250 attempts.")
        endif

        ! Set output
        orbA = orb
    end subroutine

    subroutine CSFPickBOrb  (nI, iLut, CCDblS, CCSglS, CCUnS, CCSglDelta, &
                             orbs, symProd, sumMl, orbA, sym, Ml, &
                             orbB, nopen)
        integer, intent(in) :: nI(nel), iLut(0:NIfTot), symProd
        integer, intent(in) :: sumMl, orbs(2), orbA, sym(2), Ml(2)
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCDblS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: CCSglDelta (ScratchSize/2)
        integer, intent(out) :: orbB
        integer, intent(inout) :: nopen
        character(*), parameter :: this_routine = 'CSFPickBOrb'
        real*8 :: r
        integer :: norbs, ind, orb, i, sym_ind, full_ind, numB
        logical :: bSingle

        ! Draw orbitals randomly until we find one unoccupied
        sym_ind = CCIndS(sym(2), Ml(2))
        full_ind = ClassCountInd(2, sym(2), Ml(2))
        norbs = OrbClassCount(full_ind)
        do i=1,250
            call genrand_real2 (r)
            orb = int(r * norbs)
            ind = SymLabelCounts2(1,full_ind) + orb
            orb = SymLabelList2(ind)

            ! Debugging check
            if (G1(orb)%Ms /= -1) then
                call stop_all (this_routine, 'Invalid spin generated')
            endif

            ! If exciting to self, or to the other source, this is either
            ! invalid, or is a single (also disallowed).
            if ((orb == orbs(1)) .or. (orb+1 == orbs(1)) .or. &
                (orb == orbs(2)) .or. (orb+1 == orbs(2))) cycle

            ! Is this occupied
            bSingle = .false.
            if (IsOcc(iLut, orb) .or. (orb == orbA)) then
                ! Is it a single, otherwise a double, so try again.
                if (IsNotOcc (iLut, orb+1) .and. (orb+1 /= orbA)) then
                    orb = orb + 1
                    nopen = nopen - 1 ! Creating a double
                    exit
                endif
            else
                bSingle = .true. ! Exciting into vacant -> create single
                nopen = nopen + 1
                exit
            endif
        enddo

        if (i > 250) then
            write(6,'("Cannot find an unoccupied orbital for a double")')
            write(6,'("excitation after 250 attempts.")')
            write(6,'("Desired symmetry of orbital pair =",i3)') &
                symProd
            call writedet(6,nI,nel,.true.)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a double excitation after 250 attempts.")
        endif

        ! Set the output
        orbB = orb
    end subroutine

    real*8 function CSFPickOrbsProb (nI, iLut, CCDblS, CCSglS, CCUnS, &
                                     CCSglDelta, orbsWNoPair, nSing, nVac, &
                                     orbs, sym, Ml)
        integer, intent(in) :: nI(nel), ilut(0:NIfTot)
        integer, intent(in) :: orbs(2), sym(2), Ml(2), orbsWNoPair
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCDblS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: CCSglDelta (ScratchSize/2)
        integer, intent(in) :: nSing, nVac
        real*8 :: p
        integer :: numB, permutations, orbA, orbB, symB, MlB, sym_ind, i

        ! Consider the possibility of generating this either way around
        ! (The alternative was to restrict orbB > orbA, and throw away
        !  some of the random numbers)
        permutations = 2
        if (is_in_pair(orbs(1), orbs(2))) permutations = 1
        !print*, 'sym', sym(1), sym(2)
        !print*, 'ml', ml(1), ml(2)
        !print*, 'orbs', orbs(1), orbs(2)
        !print*, 'nsing, vac, nopair', nSing, nVac, orbsWNoPair
        !print*, 'permutations', permutations

        !print*, 'ccsgls', ccsgls
        !print*, 'ccuns', ccuns

        CSFPickOrbsProb = 0
        do i=1,permutations
            orbA = orbs(i)
            orbB = orbs(3-i)
            symB = sym(3-i)
            MlB = Ml(3-i)

            ! Probability of picking first orbital
            p = 1 / real (nSing + nVac - orbsWNoPair)

            ! Number of available B orbitals given A
            ! TODO: if exciting from single, reduce numA and numB!!!
            !       (numA should be taken care of w/ orbsWNoExcit)
            sym_ind = CCIndS(symB, MlB)
            numB = CCUnS(sym_ind) + CCSglS(sym_ind)! + CCSglDelta(sym_ind)
            if (sym(2) == sym(1) .and. IsOcc (ilut, ibclr(orbA-1,0)+1)) then
                numB = numB - 1
            endif
            !print*, 'perm', i, 'numb', numb

            ! Total probability of picking A-B
            p = p * (1 / real(numB))

            ! Add contribution for this permutation.
            CSFPickOrbsProb = CSFPickOrbsProb + p
        enddo
    end function

    ! Pick an electron pair randomly, which is allowed for a double excit.
    subroutine CSFPickElecPair (nI, iLut, elecs, orbs, symProd, sumMl, &
                                nSing, nDoub, CCSglDelta, nopen, pGen)
        integer, intent(in) :: nI(nel), iLut(0:NIfTot)
        integer, intent(out) :: elecs(2), orbs(2), symProd, sumMl
        integer, intent(in) :: nSing, nDoub
        integer, intent(out) :: CCSglDelta (ScratchSize/2)
        integer, intent(inout) :: nopen
        real*8, intent(inout) :: pGen
        character(*), parameter :: this_routine = 'CSFPickElecPair'
        integer ::i, elec, orb, orb2, ind, found
        real*8 :: r, pElec, pElec2
        logical :: bSecond, bSingle

        ! Pick two electrons randomly.
        ! TODO: Remove the calculation of the effect on nopen here (do it in 
        !       csf_get_excit_general).
        found = 0
        elec = 0
        CCSglDelta = 0
        pElec = 0
        do while (found < 2)
            bSingle = .false.
            elecs(1) = elec
            orbs(1) = orb
            do i=1,250
                call genrand_real2 (r)
                elec = int(nel*r) + 1
                orb = iand(nI(elec), csf_orbital_mask)

                ! Cannot pick the same electron twice
                if ((found == 1) .and. (elec == elecs(1))) cycle

                ! Is this in a doubly occupied orbital?
                ! Singly occupied must give beta elec
                ! Only pick alpha elec from double, except that we can pick
                ! the beta electron from an already chosen double.
                orb2 = ieor((orb-1),1)
                if (btest(iLut(orb2/32), mod(orb2,32))) then
                    if (G1(orb)%Ms == 1) then
                        ind = CCIndS(int(G1(orb)%Sym%S), G1(orb)%Ml)
                        CCSglDelta(ind) = CCSglDelta(ind) + 1
                        nopen = nopen + 1
                        exit
                    else if ((found == 1) .and. (orb2+1 == orbs(1))) then
                        ind = CCIndS(int(G1(orb)%Sym%S), G1(orb)%Ml)
                        CCSglDelta(ind) = CCSglDelta(ind) - 1
                        nopen = nopen - 1
                        exit
                    endif
                else
                    if (G1(orb)%Ms == 1) call stop_all (this_routine, &
                                                        "Invalid spin")
                    ind = CCIndS(int(G1(orb)%Sym%S), G1(orb)%Ml)
                    CCSglDelta(ind) = CCSglDelta(ind) - 1
                    nopen = nopen - 1
                    bSingle = .true.
                    exit
                endif
            enddo
            found = found + 1

            ! Generate the probability of picking an electron and then
            ! another. If the second electron is not the pair of the first
            ! then it is also possible to generate the other way around, so
            ! include that too.
            if ((found == 1) .or. &
                ((found == 2) .and. (ieor(elecs(1)-1,1)+1 /= elec))) then
                if (bSingle) then
                    pElec = pElec + 1/real((nSing + nDoub - 1)*(nSing+nDoub))
                else
                    pElec = pElec + 1/real((nSing + nDoub)**2)
                endif
            endif

            if (i > 250) then
                write(6,'("Cannot find excitable electron after 250 &
                          &attempts")')
                call stop_all(this_routine, "Cannot find excitable electron &
                                            &after 250 attempts")
            endif
        enddo
        elecs(2) = elec
        orbs(2) = orb
        pGen = pGen * pElec

        ! Generate the symmetry product of these two electrons
        symProd = ieor(G1(orbs(1))%Sym%S,G1(orbs(2))%Sym%S)

        ! Sum the Mls if required.
        if (tFixLz) sumMl = G1(orbs(1))%Ml + G1(orbs(2))%Ml
    end subroutine

    ! Generate a single (random) excitation on an arbitrary CSF.
    ! This routine returns the generation probability (pGen) and the
    ! excitation matrix
    subroutine CSFCreateSingleExcit (nI, nJ, CCDbl, CCSgl, CCUn, iLut, &
                                     ExcitMat, nopen, pSingle, pGen)
        integer, intent(in)  :: nI(nel), iLut(0:NIfTot), nopen
        integer, intent(out) :: nJ(nel)
        integer, intent(inout) :: ExcitMat(2,2)
        integer, intent(in)  :: CCSgl (ScratchSize) ! ClassCountSingleOcc2
        integer, intent(in)  :: CCDbl (ScratchSize) ! ClassCountDoubleOcc2
        integer, intent(in)  :: CCUn (ScratchSize)  ! ClassCountUnocc2
        real*8,  intent(in)  :: pSingle
        real*8,  intent(out) :: pGen
        character(*), parameter :: this_routine = 'CreateSingleExcit'
        integer :: elecsWNoExcits, elec, orb, orb2, spn, ind, norbs, symEx
        integer :: lnopen, ncsf, i, sym_ind, nexcit
        real*8 :: r, S
        logical :: bSingle

        ! TODO: Check that this condition is not necessary!!
        !if (tNoSingsPossible .or. tNoSymGenRandExcits) then
        !    print*, 'nosingspossible', tNoSingsPossible
        !    print*, 'nosymgenrandexcits', tNoSymGenRandExcits
        !    call stop_all (this_routine, "This condition should never be &
        !                   &broken? See symrandexcit2.F90")
        !endif

        lnopen = nopen

        ! Loop over e- pairs to count e- with no excitations
        elecsWNoExcits = 0
        do i=2,ScratchSize,2
            ! Disallowed from beta e- in doubles
            if (CCDbl(i) /= 0) then
                elecsWNoExcits = elecsWNoExcits + CCDbl(i)
            endif

            ! Only from singles if more singles or vacant with same sym
            if ((CCSgl(i)/=0) .and. (CCSgl(i)<2) .and. (CCUn(i)==0)) then
                elecsWNoExcits = elecsWNoExcits + CCSgl(i)
            endif

            ! Alpha e- from doubles to singles or vacancies
            if ((CCDbl(i-1)/=0) .and. (CCSgl(i)==0) .and. (CCUn(i)==0)) then
                elecsWNoExcits = elecsWNoExcits + CCDbl(i-1)
            endif
        enddo

        ! 250 attempts to pick randomly
        do i=1,250
            ! Pick an electron at random, and extract its orbital
            call genrand_real2(r)
            elec = int(nel*r) + 1
            orb = iand(nI(elec), csf_orbital_mask)

            ! Obtain the symmetry index and Ms indicator (1=alpha)
            spn = 2
            if (G1(orb)%Ms == 1) spn = 1
            sym_ind = ClassCountInd(spn, int(G1(orb)%Sym%S,4), G1(orb)%Ml)
            
            ! Is this electron in a singly or doubly occupied orbital?
            ! Test if there are any allowed excitations
            orb2 = ieor((orb-1), 1)
            if (btest(iLut(orb2/32), mod(orb2,32))) then
                bSingle = .false.
                if (spn == 2) then
                    nexcit = 0
                else
                    nexcit = CCSgl(sym_ind+1) + CCUn(sym_ind+1)
                endif
            else
                if (spn == 1) call stop_all(this_routine, "Invalid Spin")
                bSingle = .true.
                nexcit = CCSgl(sym_ind) + (CCUn(sym_ind)-1)
            endif
            if (nexcit /= 0) exit
        enddo

        if (i > 250) then
            write(6,'("Cannot find single excitation after 250 attempts")')
            call stop_all(this_routine, "Cannot find single excitation after &
                                        &250 attempts")
        endif

        ! If we are exciting a singly occupied e-, then flip its spin.
        symEx = sym_ind
        if (bSingle) symEx = symEx - 1

        ! Choose an (allowed) unoccupied orbital to excite to. Draw orbitals
        ! from the desired symmetry and spin until we find one unoccupied.
        ! This is method 2 from symrandexcit2.
        norbs = OrbClassCount(symEx)
        do i=1,250
            call genrand_real2(r)
            orb2 = int(norbs*r)
            ind = SymLabelCounts2(1,symEx) + orb2
            orb2 = SymLabelList2(ind)

            ! Cannot excite a single to itself
            if (bSingle .and. (orb2 == orb+1)) cycle

            ! If target available, then select it. If exciting to a vacant
            ! orbital, then select the beta version.
            if (.not.(btest(iLut((orb2-1)/32), mod(orb2-1,32)))) then
                if (.not.btest(iLut((orb2-2)/32),mod(orb2-2,32))) then
                    orb2 = orb2-1
                    if (.not.bSingle) lnopen = lnopen + 2
                else if (bSingle) then
                    lnopen = lnopen - 2
                endif

                ! If this orbital is excluded (source of 1st excit), try again
                ! Must be after above tests, as orb2 conditionally decremented
                if (orb2 == ExcitMat(1,1)) then
                    lnopen = nopen
                    cycle
                endif

                ! Found --> leave the loop.
                exit
            endif
        enddo

        if (i > 250) then
            write(6,'("Cannot find an unoccupied orbital for a single")')
            write(6,'("excitation after 250 attempts.")')
            write(6,'("Desired symmetry of unoccupied orbital =",i3)') &
                int(G1(orb)%Sym%S, 4)
            write(6,'("Num. orbitals (of correct spin) in symmetry =",i4)') &
                norbs
            write(6,'("Number of orbitals to legitimately pick =",i4)') nexcit
            call writedet(6,nI,nel,.true.)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a single excitation after 250 attempts.")
        endif

        nJ = iand(nI, csf_orbital_mask)
        ! ExcitMat is the index of the orbital to excite from, and the actual
        ! orbital to excite to
        ExcitMat(1,1) = elec
        ExcitMat(2,1) = orb2

        ! Call csf_find_excit_det to generate the determinants and return the
        ! number of possible csfs associated with the spacial configuration.
        ncsf = 1
        call csf_find_excit_det (ExcitMat(:,1), nJ, iLut, nopen, &
                                 lnopen, ncsf, .true.)

        ! Generation probability
        pGen = pSingle / real(nexcit * (nel - elecsWNoExcits) * ncsf) 
        !print*, 'orbs', orb, orb2
        !print*, 'sym init, final', G1(orb)%Sym%S, G1(orb2)%Sym%S
    end subroutine

    ! ClassCountIndex for the spacial arrays
    integer function CCIndS (sym, mom)
        integer, intent(in) :: sym
        integer, intent(in) :: mom

        CCIndS =  ((ClassCountInd(1,sym,mom)-1)/2) + 1
        return
    end function

    ! Generate three arrays indicating the number of spacial orbitals of
    ! each possible symmetry which are doubly, singly and un-occupied
    ! nb. CCDbl is the number of _pairs_ of spacial orbitals
    subroutine ConstructClassCountsSpacial (nI, nclosed, CCDblS, CCSglS,CCUnS)
        integer, intent(in) :: nI(nel), nclosed
        integer, intent(out) :: CCDblS (ScratchSize/2)
        integer, intent(out) :: CCSglS (ScratchSize/2)
        integer, intent(out) :: CCUnS (ScratchSize/2)
        character(*), parameter :: this_routine = 'ConstructClassCounts'
        integer :: i, orb, ind

        ! nb. Unoccupied array is produced from overall orbital array minus
        !     the occupied electrons
        CCDblS = 0
        CCSglS = 0
        forall (i=1:ScratchSize/2) CCUnS(i) = OrbClassCount(2*(i-1)+1)
        if (tNoSymGenRandExcits) then
            ! TODO: Implement ConstructClassCounts for tNoGenRandExcits
            call stop_all (this_routine, 'Unimplemented')
        else
            ! First loop over the closed shell electrons
            do i = 1,nclosed-1, 2
                ! Place e- into ClassCountDoubleOcc, and remove from Unocc.
                ! ind(beta) = ind(alpha) + 1 --> Can do both in one step.
                orb = iand (nI(i), csf_orbital_mask)
                ind = CCIndS (int(G1(orb)%Sym%S,4), G1(orb)%Ml)
                CCDblS(ind) = CCDblS(ind) + 1
                CCUnS(ind) = CCUnS(ind) - 1
            enddo

            ! Now loop over the open shell electrons
            do i = nclosed+1, nel
                orb = iand (nI(i), csf_orbital_mask)
                ind = CCIndS (int(G1(orb)%Sym%S,4), G1(orb)%Ml)
                CCSglS(ind) = CCSglS(ind) + 1
                CCUnS(ind) = CCUnS(ind) - 1
            enddo
        endif
    end subroutine

    ! Generate three arrays indicating the number of orbitals of each
    ! possible symmetry which are doubly-, singly- and un-occupied.
    subroutine ConstructClassCounts (nI, nclosed, CCDbl, CCSgl, CCUn)
        integer, intent(in) :: nI(nel), nclosed
        integer, intent(out) :: CCDbl (ScratchSize) ! ClassCountDoubleOcc2
        integer, intent(out) :: CCSgl (ScratchSize) ! ClassCountSingleOcc2
        integer, intent(out) :: CCUn (ScratchSize)  ! ClassCountUnocc2
        character(*), parameter :: this_routine = 'ConstructClassCounts'
        integer :: i, orb, ind

        ! nb. Unoccupied array is produced from overall orbital array minus
        !     the occupied electrons
        CCDbl = 0
        CCSgl = 0
        CCUn = OrbClasscount
        if (tNoSymGenRandExcits) then
            ! TODO: Implement ConstructClassCounts for tNoGenRandExcits
            call stop_all (this_routine, 'Unimplemented')
        else
            ! First loop over the closed shell electrons
            do i = 1,nclosed-1, 2
                ! Place e- into ClassCountDoubleOcc, and remove from Unocc.
                ! ind(beta) = ind(alpha) + 1 --> Can do both in one step.
                orb = iand (nI(i), csf_orbital_mask)
                ind = ClassCountInd(1, int(G1(orb)%Sym%S,4), G1(orb)%Ml)
                CCDbl(ind:ind+1) = CCDbl(ind:ind+1) + 1
                CCUn(ind:ind+1) = CCUn(ind:ind+1) - 1
            enddo

            ! Now loop over the open shell electrons
            do i = nclosed+1, nel
                orb = iand (nI(i), csf_orbital_mask)
                ind = ClassCountInd(2, int(G1(orb)%Sym%S,4), G1(orb)%Ml)
                CCSgl(ind) = CCSgl(ind) + 1
                CCUn(ind) = CCUn(ind) - 1
            enddo
        endif
    end subroutine

    ! Generate a general excitation
    ! ExcitMat (1,:) are the electron indices to excite from
    ! ExcitMat (2,:) are the target orbitals.
    ! ExcitMat must be sorted, ExcitMat(1,1) < ExcitMat(1,2) and for (2,...)
    ! If array yamas is provided, then nopen2 must be correct. Otherwise
    ! returns the new value of nopen.
    subroutine csf_find_excit_det_general (ExcitMat, nJ, iLut, nopen, IC, &
                                           ncsf, nopen2, yamas)
        integer, intent(in) :: nopen, iLut(0:NIfTot), IC
        integer, intent(inout) :: nJ(ncsf,nel), ExcitMat(2,1:IC)
        integer, intent(inout) :: ncsf
        integer, intent(inout) :: nopen2
        integer, intent(in), optional :: yamas(ncsf, nopen2)
        character(*), parameter :: this_routine = 'csf_find_excit_det_general'
        integer, dimension(1:IC) :: src, src_pair, srcB, exA, exB
        integer :: nI(nel), lnopen
        integer :: singNew (4), singOld (4)
        integer :: newS, oldS
        integer :: pos, dpos, spos, epos, srcpos, oldpos, newpos, nclosed, i
        logical :: bExcitIsDoub, bExcitFromDoub, bRetry

        ! Initial values
        newS = 0
        oldS = 0
        lnopen = nopen
        nclosed = nel - nopen

        ! Strip csf data from the determinant
        nI = iand(nJ(1,:), csf_orbital_mask)
        !nJ = 0

        ! Obtain the orbital and their pairs for each of the source orbs
        src = nI(ExcitMat(1,:))
        src_pair = ieor(src-1,1) + 1
        srcB = ibclr(src-1,0) + 1

        ! The alpha and beta orbitals for the specified spacial excitation.
        exB = ibclr(ExcitMat(2,:)-1,0) + 1
        exA = ibset(Excitmat(2,:)-1,0) + 1

        ! Are the target orbitals a new double?
        bExcitIsDoub = .false.
        if ((IC == 2) .and. (ExcitMat(2,1)==(ieor(ExcitMat(2,2)-1,1)+1))) then
            bExcitIsDoub = .true.
        endif

        ! Are the source orbitals a double?
        bExcitFromDoub = .false.
        if ((IC == 2) .and. (src(1) == src_pair(2))) then
            bExcitFromDoub = .true.
        endif

        ! Fill up all of the doubles correctly.
        dpos = 1
        epos = 1
        srcpos = 1
        do i=1,nel,2
            ! Have we exhausted all the possible new doubles?
            if ((epos > IC) .and. (dpos > nclosed)) exit

            bRetry = .true.
            do while (bRetry .and. ((epos<=IC) .or. (dpos<=nclosed)))
                ! Are we inserting an entirely new excitation?
                if (bExcitIsDoub .and. (epos == 1)) then
                    if ((dpos>nclosed) .or. (ExcitMat(2,1) < nI(dpos))) then
                        nJ(1,i:i+1) = ExcitMat(2,1:2)
                        epos = 3
                        bRetry = .false.
                        cycle
                    endif
                endif

                ! Is the next item an excitation of a single?
                if ((epos<=IC) .and. (.not.bExcitIsDoub) .and. &
                    ((dpos>nclosed) .or. (ExcitMat(2,epos) < nI(dpos)))) then
                    
                    ! TODO: If we assume that all is correct, we can just test
                    !       if excitmat --> alpha or beta!!!
                    if (IsOcc(ilut, exB(epos))) then
                        if (ExcitMat(2,epos) /= exA(epos)) then
                            print*, 'exB occ, excitmat/=exA', ExcitMat(2,epos),&
                            exA(epos)
                            call stop_all (this_routine, "Excit. out of order")
                        endif
                        nJ(1,i) = exB(epos)
                        nJ(1,i+1) = exA(epos)
                        oldS = oldS + 1
                        singOld(oldS) = exA(epos)
                        lnopen = lnopen - 1
                        bRetry = .false.
                    else
                        if (ExcitMat(2,epos) /= exB(epos)) then
                            print*, 'exB unocc, excitmat/=exB', ExcitMat(2,epos),&
                            exB(epos)
                            call stop_all (this_routine, "Excit. out of order")
                        endif
                        newS = newS + 1
                        singNew(newS) = ExcitMat(2,epos)
                    endif
                    epos = epos + 1
                else if (dpos < nclosed) then
                    if ((srcpos <= IC) .and. (nI(dpos)==srcB(srcpos))) then
                        if (.not. bExcitFromDoub) then
                            newS = newS + 1
                            singNew(newS) = nI(dpos)
                        endif
                        srcpos = srcpos + 1
                    else
                        nJ(1,i:i+1) = nI(dpos:dpos+1)
                        bRetry = .false.
                    endif
                    dpos = dpos + 2
                endif
            enddo
                if (bRetry) exit
        enddo

        ! Now we can fill up the singles!
        spos = nclosed + 1
        newpos = 1
        oldpos = 1
        srcpos = 1
        do i=i,nel
            if ((newpos <= newS) .and. &
               &((spos > nel) .or. (singNew(newpos) < nI(spos)))) then
                nJ(1,i) = singNew(newpos)
                newpos = newpos + 1
                lnopen = lnopen + 1
            else if (spos <= nel) then
                if ((oldpos<=oldS) .and. (nI(spos) == singOld(oldpos))) then
                    oldpos = oldpos + 1
                else if ((srcpos<=IC) .and. (nI(spos) == src(srcpos))) then
                    srcpos = srcpos + 1
                else
                    nJ(1,i) = nI(spos)
                endif
                spos = spos + 1
            endif

            if ((newpos > newS) .and. (spos > nel)) exit
        enddo

        ! Put the orbital numbers rather than electron indices into ExcitMat
        Excitmat(1,:) = src

        ! Make this into a csf that iscsf would recognise
        nJ(1,:) = ibset(nJ(1,:), csf_test_bit)

        ! Apply a random Yamanouchi symbol.
        if (lnopen > 0) then
            if (present(yamas)) then
                if (nopen2 /= lnopen) then
                    call stop_all (this_routine, "Incorrect value of nopen2")
                endif

                forall (i=2:ncsf) nJ(i,:) = nJ(1,:)
                do i=1,ncsf
                    call csf_apply_yama (nJ(i,:), yamas(i,:))
                enddo
            else
                call csf_apply_random_yama (nJ, lnopen, real(LMS/2,8), ncsf, &
                                            .false.)
            endif
        endif

        if (.not. present(yamas)) nopen2 = lnopen
    end subroutine

    ! Generate a determinant for the excitation specified in ExcitMat
    ! Note that this ASSUMES that you have got the allowed csf excitations
    ! correctly (a='alpha', b='beta', really just the occupation of spacial
    ! orbitals, always filling 'beta' first. Not really spin symmetry here)
    !:
    ! __SINGLES__:
    ! double a -> single a,   double a -> vacant b
    ! single b -> single a,   single b -> vacant b
    !
    ! This ONLY performs a single. Don't worry about doubles, as we do those
    ! as 2x a single
    subroutine csf_find_excit_det (ExcitMat, nJ, iLut, nopen, nopen_new,&
                                   ncsf, bApplyYama, yama)
        integer, intent(in) :: nopen, nopen_new, iLut(0:nIfTot)
        integer, intent(in), optional :: yama(ncsf, nopen_new)
        integer, intent(inout) :: ncsf, nJ(ncsf,nel), ExcitMat(2)
        integer :: i, pos, exbeta, exalpha, sralpha, srbeta
        integer :: ins(4), nclosed, src
        logical, intent(in) :: bApplyYama
        character(*), parameter :: this_routine = "csf_find_excit_det"

        ! Strip csf data from the determinant.
        nJ(1,:) = iand(nJ(1,:), csf_orbital_mask)

        ! Useful integer values (purely for readability)
        exbeta = ibclr(ExcitMat(2)-1,0)+1
        exalpha = ibset(ExcitMat(2)-1,0)+1
        srbeta = ibclr(nJ(1,ExcitMat(1))-1,0)+1
        sralpha = ibset(nJ(1,ExcitMat(1))-1,0)+1
        src = nJ(1, Excitmat(1))
        nclosed = nel - nopen

        !>>>! print*, 'attempting single excitation'
        !>>>! print*, ExcitMat(1), nel-nopen
        ! Are we exciting from a doubly occupied orbital?
        if (ExcitMat(1) <= nel-nopen) then
            !>>>! print*, 'exciting from a doubly occupied orbital'
            ! Are we exciting to a singly occupied, or a vacant orbital
            if (btest(iLut((exbeta-1)/32),mod(exbeta-1,32))) then
                !>>>! print*, 'exciting to single', sralpha, exalpha
                ! Find the index of the beta e- in the spacial orbital
                ! we are exciting to.
                do i=nclosed+1,nel
                    if (nJ(1,i) == exbeta) exit
                enddo
                if (i > nel) call stop_all (this_routine, &
                                            "Could not find orbital")

                ! Remove this single from the list, and place the beta
                ! e- from the original double into it.
                nJ(1,i:nel-1) = nJ(1,i+1:nel)
                ins(1) = srbeta
                call int_list_merge (nJ(1,nclosed+1:nel),ins(1:1),nopen,1)

                ! Add the new double 
                nJ(1,ExcitMat(1)-1:nclosed-2) = &
                        nJ(1,ExcitMat(1)+1:nclosed)
                ins(1) = exbeta
                ins(2) = exalpha
                call int_list_merge (nJ(1,1:nclosed),ins(1:2),nclosed-2,2)
            else
                !>>>! print*, 'exciting to vacant', sralpha, exbeta
                ! Create two singles, both of them are betas.
                ! One is the remaining e- from the original double.
                ins(1) = min(srbeta,exbeta)
                ins(2) = max(srbeta,exbeta)
                nJ(1,ExcitMat(1)-1:nel-2) = nJ(1,ExcitMat(1)+1:nel)
                call int_list_merge (nJ(1,nclosed-1:nel), ins(1:2), nopen, 2)
            endif
        else
            !>>>! print*,'exciting from a single'
            ! Are we exciting to a singly occupied, or a vacant orbital
            if (btest(iLut((exbeta-1)/32),mod(exbeta-1,32))) then
                ! Test exciting a beta -> alpha
                ! Remove both orbitals from singles
                pos = nel
                do i=nel,nclosed+1,-1
                    if (i == ExcitMat(1)) cycle
                    if (nJ(1,i) == ExcitMat(2)) cycle
                    nJ(1,pos) = nJ(1,i)
                    pos = pos - 1
                enddo
                ! Add to the list of doubles.
                ins(1) = exbeta
                ins(2) = ExcitMat(2)
                call int_list_merge (nJ(1,1:nclosed+2),ins(1:2),nclosed,2)
            else
                ! Test exciting a beta -> beta
                ! Remove orbital from singles, and insert a new one.
                nJ(1,ExcitMat(1):nel-1) = nJ(1,ExcitMat(1)+1:nel)
                inS(1) = ExcitMat(2)
                call int_list_merge (nJ(1,nclosed+1:nel), ins(1:1), nopen-1,1)                    
            endif
        endif

        ! Replace the orbital index with the orbital number in ExcitMat
        ExcitMat(1) = src

        ! Make this into a csf that iscsf would recognise
        nJ(1,:) = ibset(nJ(1,:), csf_test_bit)

        ! If we have specified a yamanouchi symbol(s) apply it. Otherwise
        ! we must pick a random one.
        if (bApplyYama) then
            if (present(yama)) then
                forall (i=2:ncsf) nJ(i,:) = nJ(1,:)
                do i=1,ncsf
                    call csf_apply_yama (nJ(i,:), yama(i,:))
                enddo
            else
                call csf_apply_random_yama (nJ, nopen_new, real(LMS/2,8), &
                                            ncsf, .false.)
            endif
        endif
    end subroutine

    subroutine csf_gen_elec_pair (nI, iLut, nopen, elecA, elecB, delta_nopen)
        integer, intent(in) :: nI(nel), iLut(0:NIfTot), nopen
        integer, intent(inout) :: elecA, elecB
        integer, intent(out) :: delta_nopen
        character(*), parameter :: this_routine = 'csf_gen_elec_pair'
        integer :: nclosed, orbA, orbB

        ! Get the first pair if elecA is set to -1
        nclosed = nel - nopen
        if (elecA == -1) then
            elecA = 1
            elecB = 2
        else
            ! Useful values
            orbA = iand(nI(elecA), csf_orbital_mask)
            orbB = iand(nI(elecB), csf_orbital_mask)

            if ((elecA <= nclosed) .and. (elecB <= nclosed)) then
                if (is_beta(orbA)) then
                    if (elecB /= elecA + 1) then
                        call stop_all (this_routine, "Invalid orbital specified")
                    endif
                    elecA = elecB
                endif

                if (elecB < nclosed - 1) then
                    elecB = elecB + 2
                else
                    if (elecB /= nclosed) then
                        call stop_all (this_routine, "Invalid orbital specified")
                    endif

                    if (nclosed < nel) then
                        elecB = elecB + 1
                    else
                        if (elecA < nclosed - 1) then
                            elecA = elecA + 1
                            elecB = elecA + 1
                        else
                            elecA = -1
                            elecB = -1
                        endif
                    endif
                endif
            else if (elecA < nopen) then
                if (elecB < nel) then
                    elecB = elecB + 1
                else
                    if ((elecA < nclosed) .or. (nopen > 1)) then
                        elecA = elecA + 1
                        elecB = elecA + 1
                    else
                        elecA = -1
                        elecB = -1
                    endif
                endif
            else
                if (elecB < nel) then
                    elecB = elecB + 1
                else
                    if (elecA < nel - 1) then
                        elecA = elecA + 1
                        elecB = elecA + 1
                    else
                        elecA = -1
                        elecB = -1
                    endif
                endif
            endif
        endif

        ! Calculate the change of nopen.
        if (elecA /= -1) then
            if (elecA > nclosed) then ! Both A, B > nclosed (as elecB > elecA)
                delta_nopen = -2
            else if (elecB > nclosed) then ! Just one > nclosed
                delta_nopen = 0
            else ! Both from doubles
                if (is_beta(iand(nI(elecA), csf_orbital_mask))) then
                    delta_nopen = 0
                else
                    delta_nopen = 2
                endif
            endif
        endif
    end subroutine


    ! TODO: Change from LMS --> STOT?
    ! TODO: Combine the counting and generation into one step - either by
    !       specifing a maximum number to generate, generating iteratively (ie
    !       generate next excit given current), or by looping twice.
    subroutine csf_gen_excits (nI, iLut, nopen, exFlag, CCDbl, &
                               CCSgl, CCUn, CCDblS, CCSglS, CCUnS, nexcit, nJ)
        use symexcit3, only: GenExcitations3
        integer, intent(in) :: nI(nel), ilut(0:NIfTot), nopen
        integer, intent(in) :: CCDbl(ScratchSize), CCSgl(ScratchSize)
        integer, intent(in) :: CCUn(ScratchSize)
        integer, intent(in) :: CCSglS (ScratchSize/2), CCDblS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: exFlag
        logical :: bYama, bDouble, bSingle
        integer, intent(out) :: nexcit
        integer, intent(out), dimension(:,:), allocatable, optional :: nJ
        character(*), parameter :: this_routine = 'csf_gen_excits'
        integer :: i, j, k, l, nclosed, elecA, elecB, symA, symB, MlB
        integer :: ind, sym_ind, spn, orb, orb2, numcsfs(-2:2), excit
        integer :: paircount, orbi, orbj, orbs(2), orb_sym(2), symProd
        integer :: ierr, orb3, ExcitMat(2,2), sumMl, syms(2)
        integer :: ncsf_S, ncsf_V, numS, numV, ncsf
        integer :: delta_nopen, lnopen, numB, lnopen2, ndets
        integer :: tmp_yama (nopen)
        integer, allocatable, dimension(:,:) :: csf0, csfp, csfm, csfpp, csfmm
        real*8 :: S
        ! Eurgh
        integer nK(nel), nL(nel), excitmat3(2,2)
        logical tFoundAll, tParity

        ! Interpret exFlag
        bYama = .false.
        bDouble = .false.
        bSingle = .false.
        if (btest(exFlag, 0)) bYama = .true.
        if (btest(exFlag, 1)) bSingle = .true.
        if (btest(exFlag, 2)) bDouble = .true.

        ! Calculate number of different Yamanouchi symbols given S
        ! and the possible values of nopen
        S = real(LMS) / 2
        nclosed = nel - nopen
        numcsfs(0) = get_num_csfs (nopen, S)
        if (nopen<nel-1) numcsfs(1) = get_num_csfs (nopen+2, S)
        if (nopen>1) numcsfs(-1) = get_num_csfs (nopen-2, S)
        if (nopen<nel-3) numcsfs(2) = get_num_csfs (nopen+4, S)
        if (nopen>3) numcsfs(-2) = get_num_csfs (nopen-4, S)

        nexcit = 0
        if (bYama) nexcit = nexcit + numcsfs(0) - 1

        if (bSingle) then
            ! Iterate over all the electrons. Select those with allowed
            ! transitions and sum the possible transitions
            do i=1,nel
                ! Obtain the orbital and its Ms/symmetry values
                orb = iand(nI(i), csf_orbital_mask)
                spn = (3 - G1(orb)%Ms) / 2 ! alpha=1, beta=2
                sym_ind = ClassCountInd(spn, int(G1(orb)%Sym%S,4), G1(orb)%Ml)

                ! Is it doubly or singly occupied
                orb2 = ieor((orb-1), 1)
                if (btest(iLut(orb2/32), mod(orb2,32))) then
                    ! Only allow transitions from doubly occupied alpha
                    if (spn == 1) then
                        nexcit = nexcit + (numcsfs(0)*CCSgl(sym_ind+1))
                        nexcit = nexcit + (numcsfs(1)*CCUn(sym_ind+1))
                    endif
                else
                    ! Only beta electrons allowed for singly occupied
                    if (spn == 1) call stop_all(this_routine, "Invalid spin")
                    nexcit = nexcit + numcsfs(-1)*(CCSgl(sym_ind)-1)
                    nexcit = nexcit + (numcsfs(0)*CCUn(sym_ind))
                endif
            enddo
        endif

        if (bDouble) then
            ! We need to iterate through all of the electron _pairs_
            elecA = -1
            elecB = -1
            paircount = 0
            call csf_gen_elec_pair(nI, iLut, nopen, elecA, elecB, delta_nopen)
            do while (elecA /= -1)
                orbs(1) = iand(nI(elecA), csf_orbital_mask)
                orbs(2) = iand(nI(elecB), csf_orbital_mask)
                syms = G1(orbs)%Sym%S
                symProd = ieor(syms(1), syms(2))
                sumMl = G1(orbs(1))%Ml + G1(orbs(2))%Ml

                ! Loop through all of the possible A orbitals.
                do i=1,nbasis-1,2
                    ! If this is doubly occupied, we cannot excite to it.
                    if (IsDoub(ilut, i)) cycle

                    ! If we are exciting from it, we cannot excite to it.
                    if (is_in_pair(orbs(1), i) .or. is_in_pair(orbs(2), i)) &
                        cycle

                    ! If exciting to a single, creates double (decrease nopen)
                    ! If to a vacant orbital, increases nopen.
                    lnopen = nopen + delta_nopen
                    orb = i
                    if (IsOcc(ilut, i)) then
                        lnopen = lnopen - 1
                        orb = orb + 1 ! Need to excite to alpha
                    else
                        lnopen = lnopen + 1
                    endif

                    ! What symmetry must the B orbital be to give the correct
                    ! symmetry product?
                    symA = int(G1(orb)%Sym%S)
                    symB = ieor(symA, symProd)

                    ! What Ml must we have?
                    MlB = 0
                    if (tFixLz) MlB = sumMl - G1(orb)%Ml

                    ! Loop over all possible B orbitals.
                    sym_ind = ClassCountInd(2, symB, MlB) ! beta
                    ind = SymLabelCounts2(1, sym_ind)
                    do j=1,OrbClassCount(sym_ind)
                        orb2 = SymLabelList2(ind+j-1)

                        ! Ensure that we only count excitations once
                        if (orb2 < orb) cycle

                        ! If this is doubly occupied, we cannot excite to it
                        if (IsDoub(ilut, orb2)) cycle

                        ! If we are exciting from it, we cannot excite to it.
                        if (is_in_pair(orbs(1), i).or.is_in_pair(orbs(2), i))&
                            cycle

                        ! If exciting to a single, create double (decrease
                        ! nopen). If to a vacant orbital, increases nopen.
                        lnopen2 = lnopen
                        if (IsOcc(ilut, orb2)) then
                            orb2 = orb2 + 1
                            ! Are we already exciting to here?
                            if (orb == orb2) cycle
                            lnopen2 = lnopen2 - 1
                        else
                            if (orb == orb2) then
                                orb2 = orb2 + 1
                                lnopen2 = lnopen2 - 1
                            else
                                lnopen2 = lnopen2 + 1
                            endif
                        endif

                        nexcit = nexcit + &
                                 num_csf_dets(numcsfs((lnopen2-nopen)/2))
                    enddo
                enddo

                ! Generate the next electon pair
                call csf_gen_elec_pair (nI, iLut, nopen, elecA, elecB, &
                                        delta_nopen)
            enddo
        endif

        !nexcit = paircount
        if (present(nJ)) then
            ! Allocate the required memory
            allocate(nJ(nexcit,nel), csf0(numcsfs(0),nopen), stat=ierr)
            if (bSingle .or. bDouble) then
                if ((ierr == 0) .and. (nopen < nel-1)) &
                    allocate(csfp(numcsfs(1), nopen+2), stat=ierr)
                if ((ierr == 0) .and. (nopen > 1)) & 
                    allocate(csfm(numcsfs(-1), nopen-2), stat=ierr)
                if (ierr /= 0) call stop_all(this_routine,"Allocation failed")
                forall (i=1:nexcit) nJ(i,:) = nI

                ! Get all the required csfs required for singles and doubles.
                call csf_get_yamas (nopen, S, csf0, numcsfs(0))
                if (nopen<nel-1) call csf_get_yamas (nopen+2, S, csfp, &
                                                     numcsfs(1))
                if (nopen>1) call csf_get_yamas (nopen-2, S, csfm, &
                                                 numcsfs(-1))
            endif
            if (ierr /= 0) call stop_all(this_routine,"Allocation failed")

            ! Generate all allowed changes of Yamanouchi symbol.
            excit = 1
            if (bYama .and. (nopen /= 0)) then
                call get_csf_yama (nI, tmp_yama)
                do i=1,numcsfs(0)
                    if (.not. int_arr_eq (tmp_yama, csf0(i,:)) excit
                        call csf_apply_yama (nJ(excit,:), csf0(i,:))
                        excit = excit + 1
                    endif
                enddo
            endif

            ! Generate all the allowed singles
            if (bSingle) then
                !print*, 'generating singles'
                do i=1,nel
                    if (excit > nexcit) &
                        call stop_all(this_routine, "Generated too many csfs")

                    ! Obtain the orbital/symmetry to excite from
                    orb = iand(nI(i), csf_orbital_mask)
                    spn = (3 - G1(orb)%Ms) / 2 ! alpha=1, beta=2
                    sym_ind = ClassCountInd(spn, int(G1(orb)%Sym%S,4), &
                                            G1(orb)%Ml)
                                            
                    ! Is the source orbital doubly occupied?
                    orb2 = ieor((orb-1), 1) ! Spacial pair (zero based)
                    if (btest(iLut(orb2/32), mod(orb2,32))) then
                        ! Only promote alpha e- from doubly occupied orbitals
                        if (spn /= 1) cycle

                        ! Loop through all symmetry related orbitals
                        ind = SymLabelCounts2(1,sym_ind)
                        !print*, 'orb, sym_ind, ind, count', orb, sym_ind, ind, OrbClassCount(sym_ind)
                        do j=1,OrbClassCount(sym_ind)
                            ExcitMat(1,1) = i
                            ! If the target orbital is filled, skip it
                            orb2 = SymLabelList2(ind+j-1)
                            if (btest(iLut((orb2-1)/32),mod(orb2-1,32))) cycle

                            ! Is this a vacant spacial orbital, or a single
                            orb3 = ieor((orb2-1),1) ! zero based
                            if (.not.btest(iLut(orb3/32),mod(orb3,32))) then
                                ! Excite into beta orbital of vacant pair
                                !print*, 'beta of vac', orb3+1
                                ExcitMat(2,1) = orb3+1
                                call csf_find_excit_det (ExcitMat(:,1), &
                                      nJ(excit:excit+numcsfs(1)-1,:), iLut, &
                                      nopen, nopen+2, numcsfs(1),.true.,csfp)
                                excit = excit + numcsfs(1)
                            else
                                ! Excite into alpha orbital of single
                                !print*, 'alpha of sing, numcsfs', orb2, numcsfs(0)
                                ExcitMat(2,1) = orb2
                                call csf_find_excit_det (ExcitMat(:,1), &
                                      nJ(excit:excit+numcsfs(0)-1,:), iLut, &
                                      nopen, nopen, numcsfs(0),.true.,csf0)
                                excit = excit + numcsfs(0)
                            endif
                        enddo
                    else ! Now consider excitations from singles.
                         ! Loop through all symmetry related orbitals
                         ! nb. spn == 2 (beta)
                         !print*, 'excit from single'
                         ind = SymLabelCounts2(1,sym_ind)
                         do j=1,OrbClassCount(sym_ind)
                            ExcitMat(1,1) = i
                             ! Cannot excite to self
                             orb2 = SymLabelList2(ind+j-1)
                             if (orb2 == orb) cycle

                             ! Is 'beta' orbital occipied?
                             if (btest(iLut((orb2-1)/32),mod(orb2-1,32))) then
                                 ! Check if 'alpha' is vacant (single->single)
                                 orb3 = ieor((orb2-1),1) ! zero based
                                 if (.not.btest(iLut(orb3/32),mod(orb3,32))) then
                                     ExcitMat(2,1) = orb3+1
                                     call csf_find_excit_det (ExcitMat(:,1), &
                                          nJ(excit:excit+numcsfs(-1)-1,:), iLut,&
                                          nopen, nopen-2, numcsfs(-1),.true.,csfm)
                                     excit = excit + numcsfs(-1)
                                 endif
                             else ! Exciting to vacant spacial pair (stay beta)
                                 ExcitMat(2,1) = orb2
                                 call csf_find_excit_det (ExcitMat(:,1), &
                                      nJ(excit:excit+numcsfs(0)-1,:), iLut, &
                                      nopen, nopen, numcsfs(0), .true.,  csf0)
                                 excit = excit + numcsfs(0)
                             endif
                         enddo
                    endif
                enddo
            endif

            if (bDouble) then
                ! Allocate csfs only needed for doubles
                if (nopen > 3) allocate(csfmm(numcsfs(-2),nopen-4), stat=ierr)
                if ((ierr == 0) .and. (nopen < nel-3)) &
                    allocate (csfpp(numcsfs(2), nopen+4), stat=ierr)
                if (ierr /= 0) call stop_all(this_routine,"Allocation failed")

                ! Acquire csfs
                if (nopen > 3) call csf_get_yamas (nopen-4, S, csfmm, &
                                                   numcsfs(-2))
                if (nopen < nel-3) call csf_get_yamas (nopen+4, S, csfpp, &
                                                       numcsfs(2))

                ! Loop through all electron pairs
                elecA = -1
                elecB = -1
                call csf_gen_elec_pair(nI, iLut, nopen, elecA, elecB, &
                                       delta_nopen)
                do while (elecA /= -1)
                    if (excit > nexcit) &
                        call stop_all(this_routine, "Generated too many csfs")
                    orbs(1) = iand(nI(elecA), csf_orbital_mask)
                    orbs(2) = iand(nI(elecB), csf_orbital_mask)
                    syms = G1(orbs)%Sym%S
                    symProd = ieor(syms(1), syms(2))
                    sumMl = G1(orbs(1))%Ml + G1(orbs(2))%Ml

                    ! Loop through all the orbitals for orbital A
                    do i=1,nbasis-1,2
                        ! If this is doubly occupied, we cannot excite to it.
                        if (IsDoub(ilut, i)) cycle

                        ! If we are exciting from it, we cannot excite to it.
                        if (is_in_pair(orbs(1), i).or.is_in_pair(orbs(2), i))&
                            cycle

                        ! If exciting to a single, create double (decrease
                        ! nopen). If to a vacant orbital, increases nopen
                        lnopen = nopen + delta_nopen
                        orb = i
                        if (IsOcc(ilut, orb)) then
                            lnopen = lnopen - 1
                            orb = orb + 1 ! Need to excite to alpha
                        else
                            lnopen = lnopen + 1
                        endif
                        
                        ! What symmetry must the B orbital be to give the 
                        ! correct symmetry product?
                        symA = int(G1(orb)%Sym%S)
                        symB = ieor(symA, symProd)

                        ! What Ml must we have?
                        MlB = 0
                        if (tFixLz) MlB = sumMl - G1(orb)%Ml

                        ! Loop over all possible B orbitals.
                        sym_ind = ClassCountInd(2, symB, MlB)
                        ind = SymLabelCounts2(1, sym_ind)
                        do j=1,OrbClassCount(sym_ind)
                            orb2 = SymLabelList2(ind+j-1)

                            ! Ensure that we only count excitations once
                            if (orb2 < orb) cycle

                            ! If this is doubly occupied, we cannot excite to
                            ! it.
                            if (IsDoub(ilut, orb2)) cycle

                            ! If we are exciting from it, we cannot excite 
                            ! to it.
                            if (is_in_pair(orbs(1), i) .or. &
                                is_in_pair(orbs(2), i)) &
                                cycle

                            ! If exciting to a single, create double (decrease
                            ! nopen). If to a vacant orbital, increases nopen.
                            lnopen2 = lnopen
                            if (IsOcc(ilut, orb2)) then
                                orb2 = orb2 + 1
                                ! Are we already exciting to here?
                                if (orb == orb2) cycle
                                lnopen2 = lnopen2 - 1
                            else
                                if (orb == orb2) then
                                    orb2 = orb2 + 1
                                    lnopen2 = lnopen2 - 1
                                else
                                    lnopen2 = lnopen2 + 1
                                endif
                            endif

                            ExcitMat(1,1) = elecA
                            ExcitMat(1,2) = elecB
                            ExcitMat(2,1) = min(orb, orb2)
                            ExcitMat(2,2) = max(orb, orb2)
                            ndets = num_csf_dets(numcsfs((lnopen2-nopen)/2))
                            ! TODO: make it nice and easy to access csf0/p/m
                            select case (lnopen2)
                            case (-4)
                                call csf_find_excit_det_general (ExcitMat, &
                                      nJ(excit:excit+ndets-1, :), ilut, &
                                      nopen, 2, ndets, lnopen2, csfmm)
                            case (-2)
                                call csf_find_excit_det_general (ExcitMat, &
                                      nJ(excit:excit+ndets-1, :), ilut, &
                                      nopen, 2, ndets, lnopen2, csfm)
                            case (0)
                                call csf_find_excit_det_general (ExcitMat, &
                                      nJ(excit:excit+ndets-1, :), ilut, &
                                      nopen, 2, ndets, lnopen2, csf0)
                            case (2)
                                call csf_find_excit_det_general (ExcitMat, &
                                      nJ(excit:excit+ndets-1, :), ilut, &
                                      nopen, 2, ndets, lnopen2, csfp)
                            case(4)
                                call csf_find_excit_det_general (ExcitMat, &
                                      nJ(excit:excit+ndets-1, :), ilut, &
                                      nopen, 2, ndets, lnopen2, csfpp)
                            case default
                                print*, 'invalid lnopen2', lnopen2
                                call stop_all (this_routine,"invalid lnopen2")
                            end select

                            excit = excit + ndets
                        enddo
                    enddo

                    ! Generate the next electon pair
                    call csf_gen_elec_pair (nI, iLut, nopen, elecA, elecB, &
                                            delta_nopen)
                enddo
            endif
            !print*, 'total excits generated', excit-1

            ! Clear up
            if (allocated(csf0)) deallocate (csf0)
            if (allocated(csfp)) deallocate (csfp)
            if (allocated(csfm)) deallocate (csfm)
            if (allocated(csfpp)) deallocate (csfpp)
            if (allocated(csfmm)) deallocate (csfmm)
        endif
    end subroutine

    integer function num_csf_dets (ncsf)
        integer, intent(in) :: ncsf
        if (ncsf == 0) then
            num_csf_dets = 1
        else
            num_csf_dets = ncsf
        endif
    end function

    subroutine TestCSF123 (nI)
        integer, intent(in) :: nI(nel)
        integer :: iLut(0:NIfTot), nopen
        integer :: CCDbl(ScratchSize), CCSgl(ScratchSize), CCUn(ScratchSize)
        integer :: CCDblS(ScratchSize/2), CCSglS(ScratchSize/2)
        integer :: CCUnS(ScratchSize/2)
        integer :: ierr, nexcit, i
        integer, allocatable, dimension(:,:) :: nK
        character(*), parameter :: this_routine = 'TestGenRandSymCSFExcit'

        ! call TestGenRandSymCSFExcit (nI, 1000000, 1.0, 1, 10000)

        ! Generate bit representation, and count open shell electrons
        call EncodeBitDet (nI, iLut)
        nopen = count_open_orbs (iLut)

        print*, 'Starting determinant:'
        call writedet(6, nI, nel, .true.)

        ! Obtain the orbital symmetries for the following steps
        call ConstructClassCounts(nI, nel-nopen, CCDbl, CCSgl, CCUn)
        call ConstructClasscountsSpacial(nI, nel-nopen, CCDblS, CCSglS, CCUnS)

        ! Enumerate all possible excitations
        call csf_gen_excits (nI, iLut, nopen, 2, CCDbl, CCSgl,&
                             CCUn, CCDblS, CCSglS, CCUnS, nexcit, nK)
        print*, 'Excitations'
        do i=1,nexcit
            call writedet(6, nK(i,:), nel, .true.)
            call TestGenRandSymCSFExcit (nK(i,:), 1000000, 1.d0, 0.d0, 2 &
                                         ,10000)
        enddo
        deallocate(nK)
    end subroutine

    ! A test routine for the CSF excitation generators. Initially generate
    ! (and count) all of the excited csfs. Then generate excititans randomly
    ! and histogram the generation probabilities.
    ! TODO: Only changing the Yamanouchi symbol
    subroutine TestGenRandSymCSFExcit (nI, iterations, pSingle, pDouble, &
                                       exFlag, writeInterval)
        integer, intent(in) :: nI(nel), iterations, exFlag, writeInterval
        real*8,  intent(in) :: pSingle, pDouble
        integer :: iLut(0:NIfTot), nJ(nel), ExcitMat(2,2), IC, nopen
        integer :: CCDbl(ScratchSize), CCSgl(ScratchSize), CCUn(ScratchSize)
        integer :: CCDblS(ScratchSize/2), CCSglS(ScratchSize/2)
        integer :: CCUnS(ScratchSize/2)
        integer :: i, j, k, l, ierr, nexcit, ind(4)
        logical :: tFilled, bTestList
        real*8  :: pGen, pYama, avContrib, avContribAll
        integer, allocatable, dimension(:,:) :: nK
        logical, allocatable :: ex_list(:)
        real*8,  allocatable, dimension(:,:) :: SinglesHist, AllSinglesHist
        real*8,  allocatable, dimension(:,:,:,:) :: DoublesHist,AllDoublesHist
        character(*), parameter :: this_routine = 'TestGenRandSymCSFExcit'

        ! Generate bit representation, and count open shell electrons
        call EncodeBitDet (nI, iLut)
        nopen = count_open_orbs (iLut)

        ! Obtain the orbital symmetries for the following steps
        call ConstructClassCounts(nI, nel-nopen, CCDbl, CCSgl, CCUn)
        call ConstructClasscountsSpacial(nI, nel-nopen, CCDblS, CCSglS, CCUnS)

        ! Probability of changing only the Yamanouchi symbol
        pYama = 1 - pSingle - pDouble

        ! Enumerate all possible excitations
        bTestList = .true.
        if (bTestList) then
            call csf_gen_excits (nI, iLut, nopen, exFlag, CCDbl, &
                            CCSgl, CCUn, CCDblS, CCSglS, CCUnS, nexcit, nK)
            write(6,*), 'Excitations'
            allocate(ex_list(nexcit))
            ex_list = .false.
            do i=1,nexcit
                call writedet(6, nK(i,:), nel, .true.)
            enddo
            write(6,'("------------------")')
        else
            call csf_gen_excits (nI, iLut, nopen, exFlag, CCDbl, &
                                CCSgl, CCUn, CCDblS, CCSglS, CCUnS, nexcit)
        endif

        ! Allocate memory for the histograms
        allocate (SinglesHist(nBasis,nBasis), &
                  AllSinglesHist(nBasis,nBasis), &
                  DoublesHist(nBasis,nBasis,nBasis,nBasis), &
                  AllDoublesHist(nBasis,nBasis,nBasis,nBasis), stat=ierr)
        if (ierr /= 0) call stop_all (this_routine,"Memory allocation failed")

        avContrib = 0
        avContribAll = 0
        SinglesHist = 0
        AllSinglesHist = 0
        DoublesHist = 0
        AllDoublesHist = 0
        tFilled = .true.
        open(9, file='AvContrib', status='unknown', position='append')
        do i=1,iterations
            ! Generate a random excitation
            call GenRandSymCSFExcit (nI, iLut, nJ, pSingle, pDouble, IC, &
                                     ExcitMat, exFlag, pGen, CCDbl, CCSgl, &
                                     CCUn,tFilled)

            ! Only average etc. for an allowed transition
            if (nJ(1) /= 0) then
                avContrib = avContrib + 1/pGen

                ! Test if this CSF was enumerated as an excitation
                if (bTestList) then
                    do j=1,nexcit
                        if (int_arr_eq(nJ,nK(j,:))) exit
                    enddo
                    ! If we haven't found it, abort.
                    if (j <= nexcit) then
                        ex_list(j) = .true.
                    else
                        call writedet (6, nJ, nel, .true.)
                        print*, 'excit', excitmat(1,1), excitmat(2,1), excitmat
                        print*, 'sym', G1(excitmat(1,1))%Sym%S, G1(excitmat(2,1))%Sym%S
                        ierr = ClassCountInd(2, int(G1(excitmat(1,1))%Sym%S,4), G1(excitmat(2,1))%Ml)
                        ic = SymLabelCounts2(1, ierr)
                        print*, 'ind (sym/list)', ierr, ic
                        print*, 'orb list', SymLabelList2(ic:ic+OrbClassCount(ierr)-1)
                        call stop_all (this_routine, "CSF not in list of &
                                      &enumerated CSFs")
                    endif
                endif

                ! Histogram the pGens
                select case (IC)
                case (1)
                    SinglesHist(ExcitMat(1,1),ExcitMat(2,1)) = &
                        SinglesHist(ExcitMat(1,1),ExcitMat(2,1)) + (1/pGen)
                case (2)
                    ! Ensure that we store this in a canonical order
                    ind(1) = min(ExcitMat(1,1), ExcitMat(1,2))
                    ind(2) = max(Excitmat(1,1), ExcitMat(1,2))
                    ind(3) = min(Excitmat(2,1), Excitmat(2,2))
                    ind(4) = max(Excitmat(2,1), Excitmat(2,2))
                    DoublesHist(ind(1),ind(2),ind(3),ind(4)) = &
                           DoublesHist(ind(1),ind(2),ind(3),ind(4)) + (1/pGen)
                case default
                    call stop_all (this_routine, "Invalid excitation.")
                endselect
            endif

            ! Take the average contribution over all processors, and write 
            ! out on node 0
            if (mod(i,writeInterval) == 0) then
                avContribAll = 0
#ifdef PARALLEL
                call MPI_Reduce (avContrib, avContribAll, 1, &
                                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                                 MPI_COMM_WORLD, ierr)
#else
                avContriball = avContrib
#endif
                if (iProcIndex == 0) then
                    print*, i, avcontribAll/real(i*nexcit*nProcessors)
                    write(9,*) i, avContribAll/real(i*nexcit*nProcessors)
                endif
            endif
        enddo
        close(9)


#ifdef PARALLEL
        ! Sum the histograms over all processors
        call MPI_Reduce (SinglesHist, AllSinglesHist, nBasis**2, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, &
                         ierr)
        call MPI_Reduce (DoublesHist, AllDoublesHist, nBasis**4, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, &
                         ierr)
#else
        AllSinglesHist = SinglesHist
        AllDoublesHist = DoublesHist
#endif

        ! Normalise the histograms and output in a readable form.
        ! These should tend to 0 or ncsf (an integer) for the excited csf.
        open (9,file="SinglesHist",status='unknown', position='append')
        do i=1,nbasis
            do j=1,nbasis
                if (AllSinglesHist(i,j) > 0) then
                    write(9,*)AllSinglesHist(i,j)/real(iterations*nProcessors)
                    pGen = AllSinglesHist(i,j)/real(iterations*nProcessors)
                endif
            enddo
        enddo
        close(9)

        ! Similarly for the doubles histograms
        open (9,file="DoublesHist",status='unknown',position='append')
        do i=1,nbasis !-1
            do j=i+1,nbasis
                do k=1,nbasis
                    do l=k+1,nbasis
                        if (AllDoublesHist(i,j,k,l) > 0) then
                            write(9,*) AllDoublesHist(i,j,k,l) / &
                                           real(iterations*nprocessors)
                        endif
                    enddo
                enddo
            enddo
        enddo
 
        if (bTestList) then
            do i=1,nexcit
                if (ex_list(i) .eqv. .false.) then
                    call stop_all (this_routine, "Excitation in list not &
                                  &generated stochastically")
                endif
            enddo
        endif

        ! Clean up
        deallocate (SinglesHist, AllSinglesHist, DoublesHist, AllDoublesHist)
        if (allocated(nK)) deallocate(nK)
        if (allocated(ex_list)) deallocate(ex_list)
    end subroutine
end module
