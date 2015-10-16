! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"
! This is a random excitation generator for use with csfs.
! Generate using a normalised and calculable probability.
module GenRandSymExcitCSF
    use Systemdata, only: nel, tNoSymGenRandExcits, G1, nbasis, STOT, LMS, &
                          nbasismax, lztot, tFixLz, iMaxLz, tTruncateCSF, &
                          csf_trunc_level, LMS
    use SymExcitDataMod
    use FciMCData, only: pSingles, pDoubles, excit_gen_store_type
    use SymData, only: TwoCycleSymGens, nSymLabels
    use csf, only: csf_orbital_mask, csf_test_bit, csf_apply_random_yama, &
                   get_num_csfs, csf_apply_yama, csf_get_yamas, write_yama, &
                   get_csf_yama, num_csf_dets, csf_get_random_det, iscsf, &
                   det_to_random_csf, get_csf_bit_yama
    use dSFMT_interface, only: genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: ClassCountInd, gen_rand_excit, &
                                    init_excit_gen_store,clean_excit_gen_store
    use DetBitOps, only: EncodeBitDet, is_canonical_ms_order, &
                         shift_det_bit_singles_to_beta, count_open_orbs
    use Determinants, only: write_det
    use Parallel_neci
    use constants, only: n_int, bits_n_int
    use bit_reps, only: NIfTot,NIfD
    use sym_general_mod, only: CCIndS
    implicit none

contains

    subroutine gen_csf_excit (nI, iLut, nJ, ilutnJ, exFlag, IC, &
                              excitMat, tParity, pGen, HElGen, store)

        ! Generate an excitation from a CSF at random, as specified by exFlag,
        ! and return the Excitation matrix and the probability of generating
        ! that excitation.
        !
        ! Exflag: bit 0 -> Enable Yamanouchi only
        !         bit 1 -> Enable Single Excitations
        !         bit 2 -> Enable Double Excitations
        !
        ! In:  nI        - CSF/determinant to excite from
        !      iLut      - Bit representation of nI
        !      tFilled   - Are the utility (scratch) arrays CCDbl... already
        !                  filled for this case
        !      CCDblS... - Arrays to fill with counts of spatial orbitals
        !                  which are singles, doubles or unoccupied

        integer, intent(in)    :: nI(nel), exFlag
        integer(kind=n_int), intent(in)    :: iLut(0:NIfTot)
        integer, intent(out)   :: nJ(nel), IC, ExcitMat(2,2)
        logical, intent(out)   :: tParity
        real(dp),  intent(out)   :: pGen
        type(excit_gen_store_type), intent(inout), target :: store

        ! Unused:
        integer(kind=n_int), intent(out) :: iLutnJ(0:niftot)
        HElement_t(dp), intent(out) :: HElGen

        ! We only need the spatial terms for the CSF stuff. However, keep the
        ! full 1-ScratchSize array, so we can pass it through to the normal
        ! excitation routines for truncated mode.
        integer, pointer :: CCDblS(:), CCSglS(:)
        integer, pointer :: CCUnS(:)

        character(*), parameter   :: this_routine = 'GenRandSymExcitCSF'
        integer :: nopen, ncsf, exTmp
        real(dp) :: r

        ! Point to the correct bits
        CCDblS => store%ClassCountOcc
        CCSglS => store%ClassCountUnocc
        CCUnS => store%scratch3

        ! Count the open shell electrons
        nopen = count_open_orbs(iLut) 

        ! If we are above the truncation level, then generate a normal,
        ! determinential, excitation rather than using CSF specific routines.
        if (tTruncateCSF .and. csf_trunc_level /= 0 .and. &
            nopen > csf_trunc_level .and. .not. iscsf(nI)) then

            select case (exFlag)
                case (2)
                    exTmp = 1
                case (4)
                    exTmp = 2
                case default
                    exTmp = 3
            end select
            call gen_rand_excit (nI, iLut, nJ, iLutnJ, exTmp, IC, ExcitMat, &
                                 tParity, pGen, HElGen, store)

            ! If we have fallen back below the truncation level, then
            ! regenerate a CSF (pick Yamanouchi symbol at random).
            call EncodeBitDet(nJ, iLutnJ)
            nopen = count_open_orbs(iLutnJ)
            if (nopen <= csf_trunc_level) then
                ncsf = det_to_random_csf (nJ)

                ! This is probably not the most efficient way to do this in
                ! light of the above, but it does get the correct Yama symbol,
                ! and the singles having been shifted to betas.
                call EncodeBitdet(nJ, iLutnJ)
                
                ! All of the cases where nopen will FALL below csf_trunc_level
                ! require nopen to decrease. All of the possibilities for this
                ! have no degenerate excitations giving the same spatial
                ! configuration for the given determinant
                ! --> Don't need to adjust pGen before dividing by ncsf.
                ! --> w00t!
                pGen = pgen / ncsf
            endif
            return
        endif
        tParity = .false.

        ! If the array is not already populated, perform an O[N] operation to
        ! find the number of occupied alpha/beta electrons, and number of
        ! occupied e- of each symmetry class and spin.
        if (.not. store%tFilled) then
            if ((.not.TwoCycleSymgens) .and. (.not.tNoSymGenRandExcits)) then
                write(6,'("GenRandSymExcitCSF can only be used for molecular&
                          & systems")')
                write(6,'("This is because of difficulties with other &
                          &symmetries setup.")')
                write(6,'("If you want to use these excitation generators &
                          &then add NOSYMGEN to the input to ignore symmetry &
                          &while generating excitations.")')
                call neci_flush(6)
                call stop_all(this_routine,"GenRandsymExcitCSF can only be &
                                  &used for molecular systems using symmetry")
            endif

            call ConstructClassCountsSpatial(nI, nel-nopen, CCDblS, CCSglS, &
                                             CCUnS)
            store%tFilled = .true.
        endif

        ! Select type of excitation depending on ExcitFlag.
        select case (ExFlag)
            case (1)
                IC = 0
            case (2)
                IC = 1
            case (4)
                IC = 2
            case default
                r = genrand_real2_dSFMT()
                if ((r < pSingles) .and. btest(exFlag, 1)) then
                    IC = 1
                else if ((r < pSingles+pDoubles) .and. btest(exFlag,2)) then
                    IC = 2
                else if (btest(exFlag, 0)) then
                    IC = 0
                else
                    call stop_all (this_routine, "Mismatch between pSingles, &
                                  &pDoubles and exFlag")
                endif
        end select

        ! Do what we need to do
        select case (IC)
        case (0)
            ! Change only the Yamanouchi symbol. If it is not possible to
            ! change it (i.e. if ncsf = 0,1) then return 0 determinant.
            nJ = nI
            if (nopen == 0 .or. nopen == 1) then
                nJ(1) = 0
            else if (tTruncateCSF .and. nopen > csf_trunc_level) then
                nJ(1) = 0
                call stop_all (this_routine, "We should be exciting this as &
                                             &a determinant.")
            else
                call csf_apply_random_yama (nJ, nopen, STOT, ncsf, .true.)
                if (ncsf < 2) then
                    nJ(1) = 0
                else
                    pGen = (1-pSingles-pDoubles) / (ncsf - 1)
                endif
            endif
        case (1) ! Create single excitation
            call CSFCreateSingleExcit (nI, nJ, CCDblS, CCSglS, CCUnS, iLut, &
                                    ExcitMat, nopen, pSingles, pGen)
        case (2) ! Create double excitation
            call CSFCreateDoubleExcit (nI, nJ, CCDblS, CCSglS, CCUnS, iLut, &
                                       ExcitMat, nopen, pDoubles, pGen)
        endselect

        if (.not. IsNullDet(nJ)) call EncodeBitDet (nJ, iLutnJ)
    end subroutine

    subroutine CSFCreateDoubleExcit (nI, nJ, CCDblS, CCSglS, CCUnS, iLut, &
                                     ExcitMat, nopen, pDouble, pGen)
        use csf, only: iscsf
        use SymExcitDataMod, only: pDoubNew
        integer, intent(in)  :: nI(nel), nopen 
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        integer, intent(inout) :: ExcitMat(2,2)
        integer, intent(in)  :: CCSglS (ScratchSize/2) ! ClassCountSingleOcc2
        integer, intent(in)  :: CCDblS (ScratchSize/2) ! ClassCountDoubleOcc2
        integer, intent(in)  :: CCUnS (ScratchSize/2)  ! ClassCountUnocc2
        integer :: CCSglDelta (ScratchSize/2)
        real(dp),  intent(in)  :: pDouble
        real(dp), intent(out) :: pGen
        character(*), parameter :: this_routine = 'CSFCreateDoubleExcit'
        integer :: elecs(2), orbs(2,2), symProd, sumMl, sym(2), Ml(2)
        integer :: lnopen
        integer :: ncsf, nSing, nDoub, nVac, nopen2
        integer :: orbsWNoPair

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
        orbsWNoPair = CSFCalcOrbsWNoSymmetryPair (CCSglS, CCUnS, CCSglDelta, &
                                                  orbs(1,:), symProd)

        ! If there are no orbitals we can excite to, then fail.
        ! TODO: Is this really < 2 (ie is 1 impossible, as pair valid?)
        if ((nSing + nVac) == orbsWNoPair) then
            nJ(1) = 0
            return
        endif

        ! Pick the first orbital, as one which has an allowed other excitation
        call CSFPickAOrb (nI, iLut, CCSglS, CCUnS, CCSglDelta, orbs(1,:), &
                          symProd, sumMl, orbs(2,1), sym, Ml, lnopen)

        ! Pick a second orbital, with symmetry such that the product gives
        ! symProd (conserves symmetry of original electron pair).
        call CSFPickBOrb (nI, iLut, orbs(1,:), symProd, orbs(2,1), sym, Ml, &
                          orbs(2,2), lnopen)

        ! Set this up for the find_excit_det_general routine
        ExcitMat (1,1) = min(elecs(1), elecs(2))
        ExcitMat (2,1) = min(orbs(2,1), orbs(2,2))
        ExcitMat (1,2) = max(elecs(1), elecs(2))
        ExcitMat (2,2) = max(orbs(2,1), orbs(2,2))

        ! Generate the excited determinant. ncsf specifies that we are
        ! only generating one determinant, and returns the degeneracy of the
        ! Yamanouchi symbol.
        nJ = nI
        ncsf = 1
        call csf_find_excit_det_general (ExcitMat, nJ, iLut, nopen, 2, ncsf, &
                                         nopen2)

        ! Adjust the probability for the selection of orbitals and for the
        ! degeneracy of the Yamanouchi symbol
        pGen = pGen * CSFPickOrbsProb (iLut, CCSglS, CCUnS, CCSglDelta, &
                                       orbsWNoPair, nSing, nVac, orbs(2,:), &
                                       sym, Ml)
        pGen = pGen / real(ncsf,dp)

        if (iand(NJ(nel), csf_orbital_mask) == 0)&
            call stop_all (this_routine,"alert")
    end subroutine

    integer function CSFCalcOrbsWNoSymmetryPair (CCSglS, CCUnS, CCSglDelta, &
                                                 orbs, symProd)
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: CCSglDelta (ScratchSize/2)
        integer, intent(in) :: orbs(2), symProd
        character(*), parameter :: this_routine = 'CSFCalcOrbsWNoSymmetryPair'
        integer :: numB, symA, ind, ind2, i, orbsWNoPair

        ! Do we need to correct orbsWNoPair for excitations from singles
        ! (ie we cannot excite to there)
        orbsWNoPair = 0
        if (.not. is_in_pair(orbs(1), orbs(2))) then
            do i=1,2
                if (is_beta(orbs(i))) then 
                    orbsWNoPair = orbsWNoPair + 1
                endif
            enddo
        endif

        ! Loop over all beta symmetry classes and exclude orbitals with no
        ! pair which would give symProd as the resultant.
        if (tFixLz) then
            call stop_all (this_routine, 'Not implemented yet')
            ! TODO: Understand tFixLz here
        else
            ! Loop over all sym indices for alpha elecs. If there are none, 
            ! then we cannot excite to any orbitals (A) with the paired 
            ! symmetry (ind2)
            do i=0,nSymLabels-1
                ind = CCIndS(i, 0)
                numB = CCUnS(ind) + CCSglS(ind) + CCSglDelta(ind)
                symA = ieor(i, symProd)
                ind2 = CCIndS (symA, 0)

                ! If no B orbitals, then A not allowed. Also cannot excite two
                ! electrons to the same orbital (ie. if same sym, numB >= 2)
                if ((numB == 0) .or. &
                    (numB==1 .and. ind2 == ind .and. CCUnS(ind) == 0)) then
                    orbsWNoPair = orbsWNoPair + CCUnS(ind2) + CCSglS(Ind2) + &
                                  CCSglDelta(ind2)
                endif
            enddo
        endif

        CSFCalcOrbsWNoSymmetryPair = orbsWNoPair
    end function

    subroutine CSFPickAOrb (nI, iLut, CCSglS, CCUnS, CCSglDelta, &
                            orbs, symProd, sumMl, orbA, sym, Ml, nopen)
        integer, intent(in) :: nI(nel),symProd
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot) 
        integer, intent(in) :: sumMl, orbs(2)
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(inout) :: CCSglDelta (ScratchSize/2)
        integer, intent(out) :: orbA, sym(2), Ml(2)
        integer, intent(inout) :: nopen
        character(*), parameter :: this_routine = 'CSFPickAOrb'
        real(dp) :: r
        logical :: bSingle
        integer :: orb, ind, i, locDelta, lnopen

        ! Draw orbitals randomly until we find one unoccupied
        lnopen = 0
        do i=1,250
            r = genrand_real2_dSFMT()
            orb = 2 * int(r*nBasis/2) + 1 ! beta

            ! If exciting to self, or to the other source, this is either
            ! invalid, or is a single (also disallowed).
            if (is_in_pair(orb, orbs(1)) .or. is_in_pair(orb, orbs(2))) &
                cycle

            ! Is this occupied
            bSingle = .false.
            lnopen = nopen
            if (IsOcc(iLut, orb)) then
                ! Is it a single, otherwise a double, so try again.
                if (IsNotOcc(iLut, orb+1)) then
                    orb = orb + 1
                    lnopen = lnopen - 1
                else
                    cycle
                endif
            else
                bSingle = .true. ! Exciting into vacant -> create single
                lnopen = lnopen + 1
            endif

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

            ! If we are picking a single of the correct symmetry, it is not
            ! available to pick as a B orbital --> small correction.
            locDelta = 0
            if ((ind == CCIndS(sym(1), Ml(1))) .and. (.not. bSingle)) then
                locDelta = -1
            endif

            if (CCUnS(ind) > 0) exit
            if (CCSglS(ind)+CCSglDelta(ind)+locDelta > 0) exit
        enddo

        if (i > 250) then
            write(6,'("Cannot find an unoccupied orbital for a double")')
            write(6,'("excitation after 250 attempts.")')
            write(6,'("Desired symmetry of orbital pair =",i3)') &
                symProd
            write(6,*) 'src', orbs(1), orbs(2)
            call write_det (6, nI, .true.)
            call neci_flush(6)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a double excitation after 250 attempts.")
        endif

        ! Set output
        orbA = orb
        nopen = lnopen
    end subroutine

    subroutine CSFPickBOrb  (nI, iLut, orbs, symProd, orbA, sym, Ml, orbB, &
                             nopen)
        integer, intent(in) :: nI(nel), symProd
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot) 
        integer, intent(in) :: orbs(2), orbA, sym(2), Ml(2)
        integer, intent(out) :: orbB
        integer, intent(inout) :: nopen
        character(*), parameter :: this_routine = 'CSFPickBOrb'
        real(dp) :: r
        integer :: norbs, ind, orb, i, sym_ind, full_ind, lnopen
        logical :: bSingle

        ! Draw orbitals randomly until we find one unoccupied
        sym_ind = CCIndS(sym(2), Ml(2))
        full_ind = ClassCountInd(2, sym(2), Ml(2))
        norbs = OrbClassCount(full_ind)
        lnopen = 0
        do i=1,250
            r = genrand_real2_dSFMT()
            orb = int(r * norbs)
            ind = SymLabelCounts2(1,full_ind) + orb
            orb = SymLabelList2(ind)

            ! Debugging check
            if (G1(orb)%Ms /= -1) then
                call stop_all (this_routine, 'Invalid spin generated')
            endif

            ! If exciting to self, or to the other source, this is either
            ! invalid, or is a single (also disallowed).
            if (is_in_pair(orbs(1), orb) .or. is_in_pair(orbs(2), orb)) &
                cycle

            ! Is this occupied
            bSingle = .false.
            lnopen = nopen
            if (IsOcc(iLut, orb) .or. (orb == orbA)) then
                ! Is it a single, otherwise a double, so try again.
                if (IsNotOcc (iLut, orb+1) .and. (orb+1 /= orbA)) then
                    orb = orb + 1
                    lnopen = lnopen - 1 ! Creating a double
                    exit
                endif
            else
                bSingle = .true. ! Exciting into vacant -> create single
                lnopen = lnopen + 1
                exit
            endif
        enddo

        if (i > 250) then
            write(6,'("Cannot find an unoccupied orbital for a double")')
            write(6,'("excitation after 250 attempts.")')
            write(6,'("Desired symmetry of orbital pair =",i3)') &
                symProd
            call write_det (6, nI, .true.)
            write(6,*) 'src', orbs(1), orbs(2)
            write(6,*) 'tgt', orbA
            write(6,*) 'num orbs', nbasis
            call neci_flush(6)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a double excitation after 250 attempts.")
        endif

        ! Set the output
        orbB = orb
        nopen = lnopen
    end subroutine

    real(dp) function CSFPickOrbsProb (iLut, CCSglS, CCUnS, CCSglDelta, &
                                     orbsWNoPair, nSing, nVac, orbs, sym, Ml)
        integer(kind=n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: orbs(2), sym(2), Ml(2), orbsWNoPair
        integer, intent(in) :: CCSglS (ScratchSize/2)
        integer, intent(in) :: CCUnS (ScratchSize/2)
        integer, intent(in) :: CCSglDelta (ScratchSize/2)
        integer, intent(in) :: nSing, nVac
        real(dp) :: p
        integer :: numB, permutations, orbA, symB, MlB, sym_ind, i, borb

        ! Consider the possibility of generating this either way around
        ! (The alternative was to restrict orbB > orbA, and throw away
        !  some of the random numbers)
        permutations = 2
        if (is_in_pair(orbs(1), orbs(2))) permutations = 1

        CSFPickOrbsProb = 0
        do i=1,permutations
            orbA = orbs(i)
            symB = sym(3-i)
            MlB = Ml(3-i)

            ! Number of available B orbitals given A
            sym_ind = CCIndS(symB, MlB)
            numB = CCUnS(sym_ind) + CCSglS(sym_ind) + CCSglDelta(sym_ind)
            borb=get_beta(orbA)
            if (sym(2) == sym(1) .and. IsOcc (ilut, borb)) then
                numB = numB - 1
            endif

            ! Total probability of picking A-B
            p = 1 / real(numB,dp)

            ! Add contribution for this permutation.
            CSFPickOrbsProb = CSFPickOrbsProb + p
        enddo

        ! Probability of picking first orbital
        CSFPickOrbsProb = CSFPickOrbsProb  / real(nSing + nVac - orbsWNoPair,dp)
    end function

    ! Pick an electron pair randomly, which is allowed for a double excit.
    subroutine CSFPickElecPair (nI, iLut, elecs, orbs, symProd, sumMl, &
                                nSing, nDoub, CCSglDelta, nopen, pGen)
        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        integer, intent(out) :: elecs(2), orbs(2), symProd, sumMl
        integer, intent(in) :: nSing, nDoub
        integer, intent(out) :: CCSglDelta (ScratchSize/2)
        integer, intent(inout) :: nopen
        real(dp), intent(inout) :: pGen
        character(*), parameter :: this_routine = 'CSFPickElecPair'
        integer ::i, elec, orb, orb2, ind, found
        real(dp) :: r, pElec
        logical :: bSingle

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
                r = genrand_real2_dSFMT()
                elec = int(nel*r) + 1
                orb = iand(nI(elec), csf_orbital_mask)

                ! Cannot pick the same electron twice
                if ((found == 1) .and. (elec == elecs(1))) cycle

                ! Is this in a doubly occupied orbital?
                ! Singly occupied must give beta elec
                ! Only pick alpha elec from double, except that we can pick
                ! the beta electron from an already chosen double.
                orb2 = ab_pair(orb)
                if (IsOcc(ilut, orb2)) then
                    if (G1(orb)%Ms == 1) then
                        nopen = nopen + 1
                        exit
                    else if ((found == 1) .and. (orb2 == orbs(1))) then
                        nopen = nopen - 1
                        exit
                    endif
                else
                    if (G1(orb)%Ms == 1) then
                        call write_det (6, nI, .true.)
                        write (6,'("elec, orb: ",2i4)') elec, orb
                        call stop_all (this_routine, "Invalid spin")
                    endif
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
                ((found == 2) .and. (.not. is_in_pair(orbs(1),orb)))) then
                if (bSingle) then
                    pElec = pElec + 1/real((nSing + nDoub - 1)*(nSing+nDoub),dp)
                else
                    pElec = pElec + 1/real((nSing + nDoub)**2,dp)
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
        symProd = int(ieor(G1(orbs(1))%Sym%S,G1(orbs(2))%Sym%S),sizeof_int)

        ! Sum the Mls if required.
        if (tFixLz) sumMl = G1(orbs(1))%Ml + G1(orbs(2))%Ml
    end subroutine

    ! Generate a single (random) excitation on an arbitrary CSF.
    ! This routine returns the generation probability (pGen) and the
    ! excitation matrix
    subroutine CSFCreateSingleExcit (nI, nJ, CCDblS, CCSglS, CCUnS, iLut, &
                                     ExcitMat, nopen, pSingle, pGen)
        integer, intent(in)  :: nI(nel), nopen 
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot) 
        integer, intent(out) :: nJ(nel)
        integer, intent(inout) :: ExcitMat(2,2)
        integer, intent(in)  :: CCUnS(ScratchSize/2), CCDblS(ScratchSize/2)
        integer, intent(in)  :: CCSglS(ScratchSize/2)
        real(dp),  intent(in)  :: pSingle
        real(dp),  intent(out) :: pGen
        character(*), parameter :: this_routine = 'CSFCreateSingleExcit'
        integer :: elecsWNoExcits, elec, orb, orb2, spn, ind, norbs, symEx
        integer :: lnopen, ncsf, i, sym_ind, nexcit
        real(dp) :: r
        logical :: bSingle

        ! TODO: Check that this condition is not necessary!!
        !if (tNoSingsPossible .or. tNoSymGenRandExcits) then
        !    print*, 'nosingspossible', tNoSingsPossible
        !    print*, 'nosymgenrandexcits', tNoSymGenRandExcits
        !    call stop_all (this_routine, "This condition should never be &
        !                   &broken? See symrandexcit2.F90")
        !endif

        ! Keep track of the number of open shell electrons
        lnopen = nopen

        ! Loop over e- pairs to count e- with no excitations
        elecsWNoExcits = 0
        do i=1,ScratchSize/2
            ! Disallowed from beta e- in doubles
            if (CCDblS(i) /= 0) then
                elecsWNoExcits = elecsWNoExcits + CCDblS(i)
            endif

            ! Only from singles if more singles or vacant with same sym
            if ((CCSglS(i)/=0) .and. (CCSglS(i)<2) .and. (CCUnS(i)==0)) then
                elecsWNoExcits = elecsWNoExcits + CCSglS(i)
            endif

            ! Alpha e- from doubles to singles or vacancies
            if ((CCDblS(i)/=0) .and. (CCSglS(i)==0) .and. (CCUnS(i)==0)) then
                elecsWNoExcits = elecsWNoExcits + CCDblS(i)
            endif
        enddo

        ! If there are no singles from this CSF, return failure
        if (elecsWNoExcits == nel) then
            nJ(1) = 0
            return
        endif

        ! 250 attempts to pick randomly
        do i=1,250
            ! Pick an electron at random, and extract its orbital
            r = genrand_real2_dSFMT()
            elec = int(nel*r) + 1
            orb = iand(nI(elec), csf_orbital_mask)

            ! Obtain the symmetry index and Ms indicator (1=alpha)
            spn = 2
            if (G1(orb)%Ms == 1) spn = 1
            sym_ind = CCIndS(G1(orb)%Sym%S, G1(orb)%Ml)
            
            ! Is this electron in a singly or doubly occupied orbital?
            ! Test if there are any allowed excitations
            orb2 = ieor((orb-1), 1)
            if (btest(iLut(orb2/bits_n_int), mod(orb2,bits_n_int))) then
                bSingle = .false.
                if (spn == 2) then
                    nexcit = 0
                else
                    nexcit = CCSglS(sym_ind) + CCUnS(sym_ind)
                endif
            else
                if (spn == 1) call stop_all(this_routine, "Invalid Spin")
                bSingle = .true.
                nexcit = CCSglS(sym_ind) + (CCUnS(sym_ind)-1)
            endif
            if (nexcit /= 0) exit
        enddo

        if (i > 250) then
            write(6,'("Cannot find single excitation after 250 attempts")')
            call write_det (6, nI, .true.)
            call neci_flush(6)
            call stop_all(this_routine, "Cannot find single excitation after &
                                        &250 attempts")
        endif

        ! If we are exciting a singly occupied e-, then flip its spin.
        symEx = ClassCountInd(spn,G1(orb)%Sym%S, G1(orb)%Ml)
        if (bSingle) symEx = symEx - 1

        ! Choose an (allowed) unoccupied orbital to excite to. Draw orbitals
        ! from the desired symmetry and spin until we find one unoccupied.
        ! This is method 2 from symrandexcit2.
        norbs = OrbClassCount(symEx)
        do i=1,250
         r = genrand_real2_dSFMT()
            orb2 = int(norbs*r)
            ind = SymLabelCounts2(1,symEx) + orb2
            orb2 = SymLabelList2(ind)
            
            ! Cannot excite a single to itself
            if (bSingle .and. (orb2 == orb+1)) cycle

            ! If target available, then select it. If exciting to a vacant
            ! orbital, then select the beta version.
            if (.not.(btest(iLut((orb2-1)/bits_n_int), mod(orb2-1,bits_n_int)))) then
                if (.not.btest(iLut((orb2-2)/bits_n_int),mod(orb2-2,bits_n_int))) then
                    orb2 = orb2-1
                    if (.not.bSingle) lnopen = lnopen + 2
                else if (bSingle) then
                    lnopen = lnopen - 2
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
            write(6,'("Exciting from: ",i4)') orb
            print*, G1(orb)%Sym%S, symEx
            print*, SymLabelList2(SymLabelCounts2(1,symEx):SymLabelCounts2(1,symEx)+OrbClassCount(symEx)-1)
            write(6,'("Number of orbitals to legitimately pick =",i4)') nexcit
            call write_det (6, nI, .true.)
            call neci_flush(6)
            call stop_all(this_routine, "Cannot find an unoccupied orbital &
                         &for a single excitation after 250 attempts.")
        endif

        nJ = iand(nI, csf_orbital_mask)

        ! ExcitMat is the index of the orbital to excite from, and the actual
        ! orbital to excite to
        ExcitMat(1,1) = elec
        ExcitMat(2,1) = orb2

        ! Call csf_find_excit_det to generate the determinants and return the
        ! number of possible csfs associated with the spatial configuration.
        ncsf = 1
        call csf_find_excit_det_general (Excitmat(:,1), nJ, iLut, nopen, 1, &
                                         ncsf, lnopen)
        ! TODO: decide which version to use here.
        !call csf_find_excit_det (ExcitMat(:,1), nJ, iLut, nopen, &
        !                         lnopen, ncsf, .true.)

        ! Generation probability
        pGen = pSingle / real(nexcit * (nel - elecsWNoExcits) * ncsf,dp) 
    end subroutine

    ! Generate three arrays indicating the number of spatial orbitals of
    ! each possible symmetry which are doubly, singly and un-occupied
    ! nb. CCDbl is the number of _pairs_ of spatial orbitals
    subroutine ConstructClassCountsSpatial (nI, nclosed, CCDblS, CCSglS,CCUnS)
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
                ind = CCIndS (G1(orb)%Sym%S, G1(orb)%Ml)
                CCDblS(ind) = CCDblS(ind) + 1
                CCUnS(ind) = CCUnS(ind) - 1
            enddo

            ! Now loop over the open shell electrons
            do i = nclosed+1, nel
                orb = iand (nI(i), csf_orbital_mask)
                ind = CCIndS (G1(orb)%Sym%S, G1(orb)%Ml)
                CCSglS(ind) = CCSglS(ind) + 1
                CCUnS(ind) = CCUnS(ind) - 1
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
        integer, intent(inout) :: ncsf
        integer, intent(in) :: nopen, IC
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        integer, intent(inout) :: nJ(ncsf,nel), ExcitMat(2,1:IC)
        integer, intent(inout) :: nopen2
        integer, intent(in), optional :: yamas(ncsf, nopen2)
        character(*), parameter :: this_routine = 'csf_find_excit_det_general'
        integer, dimension(1:IC) :: src, src_pair, srcA, srcB, exA, exB
        integer :: nI(nel), lnopen
        integer :: singNew (4), singOld (4)
        integer :: newS, oldS, tmp
        integer :: dpos, spos, epos, srcpos, oldpos, newpos, nclosed, i
        logical :: bExcitIsDoub, bExcitFromDoub, bRetry, bContinue

        ! Initial values
        newS = 0
        oldS = 0
        lnopen = nopen
        nclosed = nel - nopen

        ! Strip csf data from the determinant
        nI = iand(nJ(1,:), csf_orbital_mask)
        !nJ = 0

        ! Obtain the orbital and their pairs for each of the source orbs
        ! Ensure canonical ordering for src (not closed before open).
        src = nI(ExcitMat(1,:))
        if (IC == 2) then
            if (src(1) > src(2)) then
                tmp = src(1)
                src(1) = src(2)
                src(2) = tmp
            endif
        endif
        src_pair = ieor(src-1,1) + 1
        srcB = ibclr(src-1,0) + 1
        srcA = ibset(src-1,0) + 1

        ! The alpha and beta orbitals for the specified spatial excitation.
        exB = ibclr(ExcitMat(2,:)-1,0) + 1
        exA = ibset(Excitmat(2,:)-1,0) + 1

        bExcitIsDoub = .false.
        bExcitFromDoub = .false.
        if (IC == 2) then
            ! Are the target orbitals a new double?
            if (ExcitMat(2,1)==(ieor(ExcitMat(2,2)-1,1)+1)) then
                bExcitIsDoub = .true.
            endif

            ! Are the source orbitals a double?
            if (src(1) == src_pair(2)) bExcitFromDoub = .true.
        endif

        ! Fill up all of the doubles correctly, whilst creating a list of
        ! singles to be destroyed (which have become doubles), and new singles
        ! to be created (where a double has been removed)
        dpos = 1   ! Position in doubles list
        epos = 1   ! Position in excitmat
        srcpos = 1 ! Position in list of src orbitals

        do i=1,nel,2
            ! Have we exhausted all the possible new doubles?
            if ((epos > IC) .and. (dpos > nclosed)) exit

            ! bRetry is set to false when we have filled nJ(i,i+1), so need to
            ! increment i.
            bRetry = .true.
            do while (bRetry .and. ((epos<=IC) .or. (dpos<=nclosed)))
                ! bContinue avoids tests which break array bounds.
                bContinue = .true.

                ! Are we inserting an entirely new excitation?
                if (bExcitIsDoub .and. (epos == 1)) then
                    if (dpos > nclosed) then
                        bContinue = .false.
                    else if (ExcitMat(2,1) < nI(dpos)) then
                        bContinue = .false.
                    endif
                    if (.not. bContinue) then
                        nJ(1,i:i+1) = ExcitMat(2,1:2)
                        epos = 3
                        bRetry = .false.
                    endif
                endif

                ! Is the next item an excitation of a single?
                if (bContinue .and. (epos<=IC) .and. (.not.bExcitIsDoub)) then
                    if (dpos > nclosed) then
                        bContinue = .false.
                    else if (ExcitMat(2,epos) < nI(dpos)) then
                        bContinue = .false.
                    endif
                    if (.not. bContinue) then
                        ! TODO: If we assume that all is correct, we can just
                        !       test if excitmat --> alpha or beta!!!
                        if (IsOcc(ilut, exB(epos))) then
                            if (ExcitMat(2,epos) /= exA(epos)) then
                                write(6,*) 'Excitation failed not alpha'
                                call write_det (6, nI, .true.)
                                write(6,*) 'excitmat: ', excitmat
                                call stop_all(this_routine, &
                                              "Excit. out of order")
                            endif

                            nJ(1,i) = exB(epos)
                            nJ(1,i+1) = exA(epos)
                            oldS = oldS + 1
                            singOld(oldS) = exB(epos)
                            lnopen = lnopen - 1
                            bRetry = .false.
                        else
                            if (ExcitMat(2,epos) /= exB(epos)) then
                                write(6,*) 'Excitation failed not beta'
                                call write_det (6, nI, .true.)
                                write(6,*) 'excitmat: ', excitmat
                                call stop_all(this_routine, &
                                              "Excit. out of order")
                            endif

                            newS = newS + 1
                            singNew(newS) = ExcitMat(2,epos)
                        endif
                        epos = epos + 1
                    endif
                endif

                ! No excitation to insert. Consider copying the next nI pair.
                if (bContinue .and. dpos < nclosed) then
                    ! Keep the source position up to date.
                    do while (srcpos <= IC)
                        if (srcB(srcpos) >= nI(dpos)) exit
                        srcpos = srcpos + 1
                    enddo

                    ! If the current pair is being excited from, don't copy it
                    ! in the doubles regime. If needed, create a new single.
                    if (srcpos <= IC) then
                        if (nI(dpos)==srcB(srcpos)) then
                            bContinue = .false.
                            if (.not. bExcitFromDoub) then
                                newS = newS + 1
                                singNew(newS) = nI(dpos)
                                srcpos = srcpos + 1
                            else
                                srcpos = srcpos + 2
                            endif
                        endif
                    endif
                    
                    ! If not disallowed above, then add the current double.
                    if (bContinue) then
                        nJ(1,i:i+1) = nI(dpos:dpos+1)
                        bRetry = .false.
                    endif
                    dpos = dpos + 2
                endif
            enddo
            if (bRetry) exit
        enddo

        ! Now we can fill up the singles.
        spos = nclosed + 1 ! Source for singles in nI
        newpos = 1         ! For newly created singles
        oldpos = 1         ! Singles to destroy
        srcpos = 1         ! Position in array src
        do i=i,nel
            ! bRetry set to false when a new single is added to nJ --> inc i
            bRetry = .true.
            do while (bRetry)
                ! bContinue avoids tests which break array bounds.
                bContinue = .true.

                ! Are there any new singles to add. Only do so if in order.
                if (newpos <= newS) then
                    if (spos > nel) then
                        bContinue = .false.
                    else if (singNew(newpos) < nI(spos)) then
                        bContinue = .false.
                    endif
                    if (.not. bContinue) then
                        nJ(1,i) = singNew(newpos)
                        bRetry = .false.
                        newpos = newpos + 1
                        lnopen = lnopen + 1
                    endif
                endif
                
                ! Consider adding existing singles from nI.
                if (bContinue .and. spos <= nel) then
                    ! Is the next single in the list of singles to remove?
                    if (oldpos<=oldS) then
                        if (nI(spos) == singOld(oldpos)) then
                            bContinue = .false.
                            oldpos = oldpos + 1
                        endif
                    endif

                    ! Is the next single being excited from
                    if (bContinue) then
                        ! Keep the source position up to date
                        do while (srcpos<=IC)
                            if (src(srcpos) >= nI(spos)) exit
                            srcpos = srcpos + 1
                        enddo

                        if (srcpos<=IC) then
                            if (nI(spos) == src(srcpos)) then
                                bContinue = .false.
                                srcpos = srcpos + 1
                            endif
                        endif
                    endif

                    ! If not disallowed above, then add the current single.
                    if (bContinue) then
                        nJ(1,i) = nI(spos)
                        bRetry = .false.
                    endif
                    spos = spos + 1
                endif

                ! If nothing left then exit
                if ((newpos > newS) .and. (spos > nel)) bRetry = .false.
            enddo

            ! If nothing left, then exit.
            if ((newpos > newS) .and. (spos > nel)) exit
        enddo

        do i=1,IC
           if (IsOcc(ilut, srcB(i)) &
           .and. IsNotOcc(ilut, srcA(i))) &
                lnopen = lnopen - 1
        enddo

        ! Put the orbital numbers rather than electron indices into ExcitMat
        Excitmat(1,:) = src

        ! Make this into a csf that iscsf would recognise
        nJ(1,:) = ibset(nJ(1,:), csf_test_bit)

        ! Apply a random Yamanouchi symbol.
        if (lnopen > 0) then
            if (present(yamas)) then
                if (nopen2 /= lnopen) then
                    call write_det (6, nJ(1,:), .true.)
                    call stop_all (this_routine, "Incorrect value of nopen2")
                endif

                forall (i=2:ncsf) nJ(i,:) = nJ(1,:)
                do i=1,ncsf
                    call csf_apply_yama (nJ(i,:), yamas(i,:))
                enddo
            else if (tTruncateCSF .and. (lnopen > csf_trunc_level)) then
                ! Use ncsf here to return number of determinants possible
                ! rather than 
                ncsf = csf_get_random_det (nJ, lnopen, LMS)
            else
                call csf_apply_random_yama (nJ, lnopen, STOT, ncsf, .false.)
            endif
        endif

        if (.not. present(yamas)) nopen2 = lnopen
    end subroutine

    ! Generate a determinant for the excitation specified in ExcitMat
    ! Note that this ASSUMES that you have got the allowed csf excitations
    ! correctly (a='alpha', b='beta', really just the occupation of spatial
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
        integer, intent(in) :: nopen, nopen_new 
        integer, intent(inout) :: ncsf
        integer(kind=n_int), intent(in) :: iLut(0:nIfTot)
        integer, intent(in), optional :: yama(ncsf, nopen_new)
        integer, intent(inout) :: nJ(ncsf,nel), ExcitMat(2)
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

        ! Are we exciting from a doubly occupied orbital?
        if (ExcitMat(1) <= nel-nopen) then
            ! Are we exciting to a singly occupied orbital?
            if (btest(iLut((exbeta-1)/bits_n_int),mod(exbeta-1,bits_n_int))) then
                ! Find the index of the beta e- in the spatial orbital
                ! we are exciting to.
                do i=nclosed+1,nel
                    if (nJ(1,i) == exbeta) exit
                enddo
                if (i > nel) call stop_all (this_routine, &
                                            "Could not find orbital")

                ! Remove this single from the list, and place the beta
                ! e- from the original double into it.
                if (i < nel) nJ(1,i:nel-1) = nJ(1,i+1:nel)
                ins(1) = srbeta
                call int_list_merge (nJ(1,nclosed+1:nel),ins(1:1),nopen-1,1)

                ! Add the new double 
                nJ(1,ExcitMat(1)-1:nclosed-2) = &
                        nJ(1,ExcitMat(1)+1:nclosed)
                ins(1) = exbeta
                ins(2) = exalpha
                call int_list_merge (nJ(1,1:nclosed),ins(1:2),nclosed-2,2)
                
            ! Exciting to vacant orbital
            else
                !>>>! print*, 'exciting to vacant', sralpha, exbeta
                ! Create two singles, both of them are betas.
                ! One is the remaining e- from the original double.
                ins(1) = min(srbeta,exbeta)
                ins(2) = max(srbeta,exbeta)
                nJ(1,ExcitMat(1)-1:nel-2) = nJ(1,ExcitMat(1)+1:nel)
                call int_list_merge (nJ(1,nclosed-1:nel), ins(1:2), nopen, 2)
            endif
        ! Exciting from a singly occupied orbital
        else
            ! Are we exciting to a singly occupied, or a vacant orbital
            if (btest(iLut((exbeta-1)/bits_n_int),mod(exbeta-1,bits_n_int))) then
                ! Test exciting a beta -> alpha
                ! Remove both orbitals from singles
                pos = nel
                do i=nel,nclosed+1,-1
                    if (i == ExcitMat(1)) cycle
                    if (nJ(1,i) == exbeta) cycle
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
        if (bApplyYama .and. nopen_new > 0) then
            if (present(yama)) then
                forall (i=2:ncsf) nJ(i,:) = nJ(1,:)
                do i=1,ncsf
                    call csf_apply_yama (nJ(i,:), yama(i,:))
                enddo
            else
                call csf_apply_random_yama (nJ, nopen_new, STOT, ncsf,.false.)
            endif
        endif
    end subroutine

    ! Given a pair of electrons, or an initialisation value, generate the next
    ! pair of electrons for a given system. Selecting always the alpha elec
    ! from a double unless getting an AB pair.
    subroutine csf_gen_elec_pair (nI, nopen, elecA, elecB, delta_nopen)
        integer, intent(in) :: nI(nel), nopen
        integer, intent(inout) :: elecA, elecB
        integer, intent(out) :: delta_nopen
        character(*), parameter :: this_routine = 'csf_gen_elec_pair'
        integer :: nclosed, orbA, orbB

        nclosed = nel - nopen

        ! If there are fewer than two electrons, we cannot generate a pair
        if (nel < 2) then
            elecA = -1
            elecB = -1
        ! Get the first pair if elecA is set to -1
        else if (elecA == -1) then
            elecA = 1
            elecB = 2
        ! Given a valid pair, generate the next one.
        else
            ! The orbitals associated with the electrons
            orbA = iand(nI(elecA), csf_orbital_mask)
            orbB = iand(nI(elecB), csf_orbital_mask)

            ! Both electrons in the closed shell region
            if ((elecA <= nclosed) .and. (elecB <= nclosed)) then
                ! If AB pair, move elecA to alpha for remainder
                if (is_beta(orbA)) then
                    if (elecB /= elecA + 1) call stop_all (this_routine, &
                                                  "Invalid orbital specified")
                    elecA = elecB
                endif

                ! Increment (alpha) elec B to next closed pair if possible.
                ! Otherwise move elecB into open region. If no open region,
                ! increment to next AB pair. If on last pair, signal end.
                if (elecB < nclosed - 1) then
                    elecB = elecB + 2
                else
                    if (elecB /= nclosed) call stop_all (this_routine, &
                                                  "Invalid orbital specified")

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
            ! elecB in open region, elecA in open or closed region (elecA 
            ! needs to be incremented by one in either case, ether to new open
            ! position, or to the beta part of an AB pair.
            else
                ! Increment elecB. If at end, increment elecA and place elecB
                ! after it. If no more space, then signal end.
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

        ! Calculate the change of nopen if we remove the two selected elecs
        ! from the CSF nI
        if (elecA /= -1) then
            ! Removing two open shell elecs
            if (elecA > nclosed) then ! Both A, B > nclosed (as elecB > elecA)
                delta_nopen = -2
            ! Remove one open shell, and one closed elec (cancels)
            else if (elecB > nclosed) then ! Just one > nclosed
                delta_nopen = 0
            ! Two closed shell elecs. If from same pair, then no effect on
            ! nopen, if different pairs, creates two new singles.
            else
                if (is_beta(iand(nI(elecA), csf_orbital_mask))) then
                    delta_nopen = 0
                else
                    delta_nopen = 2
                endif
            endif
        endif
    end subroutine

    ! TODO: Combine the counting and generation into one step - either by
    !       specifing a maximum number to generate, generating iteratively (ie
    !       generate next excit given current), or by looping twice.
    subroutine csf_gen_excits (nI, iLut, nopen, exFlag, CCSglS, CCUnS, &
                               nexcit, nJ, return_excits)
        use symexcit3, only: GenExcitations3
        integer, intent(in) :: nI(nel), nopen
        integer(kind=n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: CCSglS(ScratchSize/2)
        integer, intent(in) :: CCUnS(ScratchSize/2)
        integer, intent(in) :: exFlag
        integer, intent(out) :: nexcit
        character(*), parameter :: this_routine = 'csf_gen_excits'
        integer :: i, j, ierr, excit, ndets

        ! The outputted excitations if required.
        integer, intent(out), dimension(:,:), allocatable :: nJ
        logical :: return_excits

        ! What are we intending to generate
        logical :: bYama, bDouble, bSingle

        ! Orbitals and electrons. orbs(2) are the source orbitals for doub.
        ! Excitation matrix indicates src elecs/tgt orbitals
        integer :: orbs(2), orb, orb2, orb3, elecA, elecB, ExcitMat(2,2)

        ! Indices, spin, Ml and symmetry values.
        integer :: ind, sym_ind, spn, MlB, sumMl
        integer :: syms(2), symA, symB, symProd

        ! Temporary values of nopen.
        integer :: delta_nopen, lnopen, lnopen2

        ! Arrays of Yamanouchi symbols, and their counts. Index of numcsfs is
        ! delta_nopen/2 [(lnopen - nopen)/2].
        integer :: numcsfs(-2:2), tmp_yama(nopen)
        integer, allocatable, dimension(:,:) :: csf0, csfp, csfm, csfpp, csfmm
        
        ! Debugging parameters
        ! integer paircount, acount
        ! logical tAllowed

        ! Interpret exFlag
        bYama = .false.
        bDouble = .false.
        bSingle = .false.
        if (btest(exFlag, 0)) bYama = .true.
        if (btest(exFlag, 1)) bSingle = .true.
        if (btest(exFlag, 2)) bDouble = .true.

        ! Calculate number of different Yamanouchi symbols given S
        ! and the possible values of nopen
        numcsfs(0) = get_num_csfs (nopen, STOT)
        if (nopen<nel-1) numcsfs(1) = get_num_csfs (nopen+2, STOT)
        if (nopen>1) numcsfs(-1) = get_num_csfs (nopen-2, STOT)
        if (nopen<nel-3) numcsfs(2) = get_num_csfs (nopen+4, STOT)
        if (nopen>3) numcsfs(-2) = get_num_csfs (nopen-4, STOT)

        nexcit = 0
        if (bYama) nexcit = nexcit + numcsfs(0) - 1

        if (bSingle) then
            ! Iterate over all the electrons. Select those with allowed
            ! transitions and sum the possible transitions
            do i=1,nel
                ! Obtain the orbital and its Ms/symmetry values
                orb = iand(nI(i), csf_orbital_mask)
                spn = (3 - G1(orb)%Ms) / 2 ! alpha=1, beta=2
                sym_ind = CCIndS (G1(orb)%Sym%S, G1(orb)%Ml)

                ! Is it doubly or singly occupied
                orb2 = ieor((orb-1), 1)
                if (btest(iLut(orb2/bits_n_int), mod(orb2,bits_n_int))) then
                    ! Only allow transitions from doubly occupied alpha
                    if (spn == 1) then
                        nexcit = nexcit + (numcsfs(0)*CCSglS(sym_ind))
                        nexcit = nexcit + (numcsfs(1)*CCUnS(sym_ind))
                    endif
                else
                    ! Only beta electrons allowed for singly occupied
                    if (spn == 1) call stop_all(this_routine, "Invalid spin")
                    nexcit = nexcit + numcsfs(-1)*(CCSglS(sym_ind)-1)
                    nexcit = nexcit + (numcsfs(0)*CCUnS(sym_ind))
                endif
            enddo
        endif

        if (bDouble) then
            ! We need to iterate through all of the electron _pairs_
            elecA = -1
            elecB = -1
            !paircount = 0
            !acount = 0
            call csf_gen_elec_pair(nI, nopen, elecA, elecB, delta_nopen)
            do while (elecA /= -1)
                orbs(1) = iand(nI(elecA), csf_orbital_mask)
                orbs(2) = iand(nI(elecB), csf_orbital_mask)
                syms(1) = int(G1(orbs(1))%Sym%S,sizeof_int)
                syms(1) = int(G1(orbs(1))%Sym%S,sizeof_int)
                !syms = int(G1(orbs)%Sym%S,sizeof_int)
                symProd = ieor(syms(1), syms(2))
                sumMl = G1(orbs(1))%Ml + G1(orbs(2))%Ml
                !paircount = paircount + 1

                ! Loop through all of the possible A orbitals.
                do i=1,nbasis-1,2
                    ! If this is doubly occupied, we cannot excite to it.
                    if(IsOcc(ilut,i).and. &
                    IsOcc(ilut,ab_pair(i))) cycle

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
                    !tAllowed = .false.
                    do j=1,OrbClassCount(sym_ind)
                        orb2 = SymLabelList2(ind+j-1)

                        ! Ensure that we only count excitations once
                        ! For valid acount (debugging), comment out this line.
                        if (orb2 < orb) cycle

                        ! If this is doubly occupied, we cannot excite to it
                        if(IsOcc(ilut,orb2).and. &
                         IsOcc(ilut,ab_pair(orb2))) cycle

                        ! If we are exciting from it, we cannot excite to it.
                        if (is_in_pair(orbs(1), orb2) .or. &
                            is_in_pair(orbs(2), orb2))&
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
                        !tallowed = .true.
                    enddo
                    !if (tallowed) then
                    !    acount = acount + 1
                    !endif
                enddo

                ! Generate the next electon pair
                call csf_gen_elec_pair (nI, nopen, elecA, elecB, delta_nopen)
            enddo
        endif

        if (return_excits .and. nexcit /= 0) then
            ! Allocate the required memory, init. and get Yamanouchi symbols
            allocate(nJ(nexcit,nel), csf0(numcsfs(0),nopen), stat=ierr)
            forall (i=1:nexcit) nJ(i,:) = nI
            call csf_get_yamas (nopen, STOT, csf0, numcsfs(0))

            if (bSingle .or. bDouble) then
                if ((ierr == 0) .and. (nopen < nel-1)) &
                    allocate(csfp(numcsfs(1), nopen+2), stat=ierr)
                if ((ierr == 0) .and. (nopen > 1)) & 
                    allocate(csfm(numcsfs(-1), nopen-2), stat=ierr)
                if (ierr /= 0) call stop_all(this_routine,"Allocation failed")

                ! Get all the required csfs required for singles.
                if (nopen<nel-1) call csf_get_yamas (nopen+2, STOT, csfp, &
                                                     numcsfs(1))
                if (nopen>1) call csf_get_yamas (nopen-2, STOT, csfm, &
                                                 numcsfs(-1))
            endif
            if (ierr /= 0) call stop_all(this_routine,"Allocation failed")

            ! Generate all allowed changes of Yamanouchi symbol.
            excit = 1
            if (bYama .and. (nopen > 2)) then
                call get_csf_yama (nI, tmp_yama, nopen)
                do i=1,numcsfs(0)
                    if (any(tmp_yama /= csf0(i,:))) then
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
                    sym_ind = ClassCountInd(spn, G1(orb)%Sym%S, &
                                            G1(orb)%Ml)
                                            
                    ! Is the source orbital doubly occupied?
                    orb2 = ieor((orb-1), 1) ! Spatial pair (zero based)
                    if (btest(iLut(orb2/bits_n_int), mod(orb2,bits_n_int))) then
                        ! Only promote alpha e- from doubly occupied orbitals
                        if (spn /= 1) cycle

                        ! Loop through all symmetry related orbitals
                        ind = SymLabelCounts2(1,sym_ind)
                        !print*, 'orb, sym_ind, ind, count', orb, sym_ind, ind, OrbClassCount(sym_ind)
                        do j=1,OrbClassCount(sym_ind)
                            ExcitMat(1,1) = i
                            ! If the target orbital is filled, skip it
                            orb2 = SymLabelList2(ind+j-1)
                            if (btest(iLut((orb2-1)/bits_n_int),mod(orb2-1,bits_n_int))) cycle

                            ! Is this a vacant spatial orbital, or a single
                            orb3 = ieor((orb2-1),1) ! zero based
                            if (.not.btest(iLut(orb3/bits_n_int),mod(orb3,bits_n_int))) then
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
                             if (btest(iLut((orb2-1)/bits_n_int),mod(orb2-1,bits_n_int))) then
                                 ! Check if 'alpha' is vacant (single->single)
                                 orb3 = ieor((orb2-1),1) ! zero based
                                 if (.not.btest(iLut(orb3/bits_n_int),mod(orb3,bits_n_int))) then
                                     ExcitMat(2,1) = orb3+1
                                     call csf_find_excit_det (ExcitMat(:,1), &
                                          nJ(excit:excit+numcsfs(-1)-1,:), iLut,&
                                          nopen, nopen-2, numcsfs(-1),.true.,csfm)
                                     excit = excit + numcsfs(-1)
                                 endif
                             else ! Exciting to vacant spatial pair (stay beta)
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
                if (nopen > 3) call csf_get_yamas (nopen-4, STOT, csfmm, &
                                                   numcsfs(-2))
                if (nopen < nel-3) call csf_get_yamas (nopen+4, STOT, csfpp, &
                                                       numcsfs(2))

                ! Loop through all electron pairs
                elecA = -1
                elecB = -1
                call csf_gen_elec_pair(nI, nopen, elecA, elecB, delta_nopen)
                do while (elecA /= -1)
                    if (excit > nexcit) &
                        call stop_all(this_routine, "Generated too many csfs")
                    orbs(1) = iand(nI(elecA), csf_orbital_mask)
                    orbs(2) = iand(nI(elecB), csf_orbital_mask)
                    syms(1) = int(G1(orbs(1))%Sym%S, sizeof_int)
                    syms(2) = int(G1(orbs(2))%Sym%S, sizeof_int)
                    !syms = int(G1(orbs)%Sym%S,sizeof_int)
                    symProd = ieor(syms(1), syms(2))
                    sumMl = G1(orbs(1))%Ml + G1(orbs(2))%Ml

                    ! Loop through all the orbitals for orbital A
                    do i=1,nbasis-1,2
                        ! If this is doubly occupied, we cannot excite to it.
                        if(IsOcc(ilut,i).and. &
                         IsOcc(ilut,ab_pair(i))) cycle

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
                            if(IsOcc(ilut,orb2).and. &
                             IsOcc(ilut,ab_pair(orb2))) cycle

                            ! If we are exciting from it, we cannot excite 
                            ! to it.
                            if (is_in_pair(orbs(1), orb2) .or. &
                                is_in_pair(orbs(2), orb2)) &
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
                            select case (lnopen2-nopen)
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
                                call write_det(6, nJ(excit,:), .true.)
                                print*, excitmat(1,:), excitmat(2,:)
                                print*, 'invalid lnopen2', lnopen2, nopeN
                                print*, 'case', int((lnopen2 - nopen)/2)
                                call stop_all (this_routine,"invalid lnopen2")
                            end select

                            excit = excit + ndets
                        enddo
                    enddo

                    ! Generate the next electon pair
                    call csf_gen_elec_pair (nI, nopen, elecA, elecB, &
                                            delta_nopen)
                enddo
            endif

            ! Final test
            if (excit /= nexcit + 1) &
                call stop_all (this_routine, &
                               'Incorrect number of CSFs generated')

            ! Clear up
            if (allocated(csf0)) deallocate (csf0)
            if (allocated(csfp)) deallocate (csfp)
            if (allocated(csfm)) deallocate (csfm)
            if (allocated(csfpp)) deallocate (csfpp)
            if (allocated(csfmm)) deallocate (csfmm)
        endif

        ! For debugging probabilities.
        !nexcit = acount
        !nexcit = paircount
    end subroutine

    subroutine TestCSF123 (nI)
        integer, intent(in) :: nI(nel)
        integer(kind=n_int) :: iLut(0:NIfTot)
        integer :: CCDblS(ScratchSize/2), CCSglS(ScratchSize/2), nopen
        integer :: CCUnS(ScratchSize/2)
        integer :: nexcit, i, nFreeze
        integer, allocatable, dimension(:,:) :: nK
        !integer, dimension(10) :: nJ=(/-2147483631,-2147483630,-2147483615,-2147483614,
        !-1073741823,-2147483645,-1073741805,-1073741803,-2147483625,-2147483623/)
        !integer, dimension(10) :: nJ=(/1,2,3,4,17,18,21,22,33,34/)
        !integer, dimension(10) :: nJ=(/3,4,17,18,23,24,31,32,33,34/)
        character(*), parameter :: this_routine = 'TestGenRandSymCSFExcit'

        ! Generate bit representation, and count open shell electrons
        call EncodeBitDet (nI, iLut)
        nopen = count_open_orbs (iLut)

        write (6, '("Starting determinant: ")', advance='no')
        call write_det (6, nI, .true.)

        ! Obtain the orbital symmetries for the following steps
        call ConstructClassCountsSpatial(nI, nel-nopen, CCDblS, CCSglS, CCUnS)

        ! Enumerate all possible excitations
        call csf_gen_excits (nI, iLut, nopen, 7, CCSglS, CCUnS, &
                             nexcit, nK, .true.)

        ! Run the testing routine or all of the excitatitons of the starting
        ! CSF. Currently counts the number of excitations frozen out if 
        ! Freeze(2,0) is enabled.
        nFreeze = 0
        write (6, '(i5, " excitations found:")') nexcit
        do i=1,nexcit
            if (iand(nK(i,1), csf_orbital_mask) /= 1 .or. &
                iand(nK(i,2), csf_orbital_mask) /= 2) then
                nFreeze = nFreeze + 1
            endif

            write (6, '(i6,": ")', advance='no') i
            call write_det (6, nK(i,:), .true.)
            call TestGenRandSymCSFExcit (nK(i,:), 4000000, 0.2_dp, 0.75_dp, 7, &
                                         10000)
        enddo
        deallocate(nK)

        write (6, '("Assuming freezing two core electrons")')
        write (6, '("Num frozen: ", i5)') nFreeze
        write (6, '("Num remaining: ", i5)') nexcit - nFreeze

        ! Test the starting determinants excitations
        call TestGenRandSymCSFExcit (nI, 4000000, 0.0_dp, 1.0_dp, 4, 10000)
        call stop_all(this_routine, 'end of test')
    end subroutine

    ! A test routine for the CSF excitation generators. Initially generate
    ! (and count) all of the excited csfs. Then generate excititans randomly
    ! and histogram the generation probabilities.
    subroutine TestGenRandSymCSFExcit (nI, iterations, pSingle, pDouble, &
                                       exFlag, writeInterval)
        integer, intent(in) :: nI(nel), iterations, exFlag, writeInterval
        real(dp),  intent(in) :: pSingle, pDouble
        character(*), parameter :: this_routine = 'TestGenRandSymCSFExcit'
        integer(kind=n_int) :: iLut(0:NIfTot)
        integer :: nJ(nel), ExcitMat(2,2), IC, nopen
        integer :: i, j, k, l, ierr, nexcit, ind(4)
        logical :: bTestList, tParity
        real(dp)  :: pGen, avContrib, avContribAll
        ! Store the generated excitations and if they have been generated.
        integer, allocatable, dimension(:,:) :: nK
        logical, allocatable, dimension(:)   :: ex_list
        ! Histogram the generation probabilities.
        real(dp),  allocatable, dimension(:,:) :: SinglesHist, AllSinglesHist
        real(dp),  allocatable, dimension(:,:,:,:) :: DoublesHist,AllDoublesHist
        ! Unused
        integer(kind=n_int) :: iLutnJ(0:niftot)
        type(excit_gen_store_type), target :: store
        integer, pointer :: CCDblS(:), CCSglS(:), CCUnS(:)
        HElement_t(dp) :: HElGen

        ! Generate bit representation, and count open shell electrons
        call EncodeBitDet (nI, iLut)
        nopen = count_open_orbs (iLut)

        ! Generate excitation store
        call init_excit_gen_store (store)

        ! Obtain the orbital symmetries for the following steps
        CCDblS => store%ClassCountOcc
        CCSglS => store%ClassCountUnocc
        CCUnS => store%scratch3
        call ConstructClassCountsSpatial(nI, nel-nopen, CCDblS, CCSglS, CCUnS)
        store%tFilled = .true.

        ! Set the global variables appropriately.
        pSingles = pSingle
        pDoubles = pDouble

        ! Count the excitations
        bTestList = .true.
        if (bTestList) then
            call csf_gen_excits (nI, iLut, nopen, exFlag, CCSglS, CCUnS, &
                                 nexcit, nK, .true.)
            if (nexcit > 0) then
                allocate(ex_list(nexcit), stat=ierr)
                if (ierr /= 0) call stop_all (this_routine, &
                                              "Memory allocation failed")
                ex_list = .false.

                ! Output a list of all generated CSFs
                open (9, file='genCSF', status='unknown', position='append')
                write(9,'(i5, "Excitations from ")', advance='no') nexcit
                call write_det (9, nI, .true.)
                do i=1,nexcit
                    call write_det (9, nK(i,:), .true.)
                enddo
                write(9,'("==================")')
                close(9)
            endif
        else
            call csf_gen_excits (nI, iLut, nopen, exFlag, CCSglS, CCUnS, &
                                 nexcit, nK, .false.)
        endif

        ! If there are no possible excitations, don't mess around.
        ! (nothing will be allocated if nexcit == 0)
        if (nexcit == 0) then
            write (6,*) 'No excitations from this CSF'
            return
        endif

        ! Allocate memory for the histograms
        allocate (SinglesHist(nBasis,nBasis), &
                  AllSinglesHist(nBasis,nBasis), &
                  DoublesHist(nBasis,nBasis,nBasis,nBasis), &
                  AllDoublesHist(nBasis,nBasis,nBasis,nBasis), stat=ierr)
        if (ierr /= 0) call stop_all (this_routine,"Memory allocation failed")

        ! Initialise averages and histograms
        avContrib = 0
        avContribAll = 0
        SinglesHist = 0
        AllSinglesHist = 0
        DoublesHist = 0
        AllDoublesHist = 0

        ! Main testing loop
        open(9, file='AvContrib', status='unknown', position='append')
        do i=1,iterations
            ! Generate a random excitation
            call gen_csf_excit (nI, iLut, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                                tParity, pGen, HElGen, store)

            ! Only average etc. for an allowed transition
            if (nJ(1) /= 0) then
                avContrib = avContrib + 1/pGen

                ! Test if this CSF was enumerated as an excitation
                if (bTestList) then
                    do j=1,nexcit
                        if (all(nJ == nK(j,:))) exit
                    enddo
                    ! If we haven't found it, abort.
                    if (j <= nexcit) then
                        ex_list(j) = .true.
                    else
                        call write_det (6, nJ, .true.)
                        call stop_all (this_routine, "CSF not in list of &
                                      &enumerated CSFs")
                    endif
                endif

                ! Histogram the pGens
                select case (IC)
                case (0)
                    continue ! Not doing anything here yet...
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

                call MPIReduce(avContrib,MPI_SUM,avContribAll)

                if (iProcIndex == 0) then
                    !print*, i, avcontribAll/real(i*nexcit*nProcessors)
                    write(9,*) i, avContribAll/real(i*nexcit*nProcessors,dp)
                endif
            endif
        enddo
        close(9)
        !print*, avcontribAll/real(i*nexcit*nProcessors)

        call MPIReduce(SinglesHist,MPI_SUM,AllSinglesHist)
        call MPIReduce(DoublesHist,MPI_SUM,AllDoublesHist)

        ! Normalise the histograms and output in a readable form.
        ! These should tend to 0 or ncsf (an integer) for the excited csf.
        open (9,file="SinglesHist",status='unknown', position='append')
        do i=1,nbasis
            do j=1,nbasis
                if (AllSinglesHist(i,j) > 0) then
                    write(9,*) AllSinglesHist(i,j) / &
                                   real(iterations*nProcessors,dp)
                endif
            enddo
        enddo
        close(9)

        ! Similarly for the doubles histograms
        open (9,file="DoublesHist",status='unknown',position='append')
        call write_det (9, ni, .true.)
        do i=1,nbasis !-1
            do j=i+1,nbasis
                do k=1,nbasis
                    do l=k+1,nbasis
                        if (AllDoublesHist(i,j,k,l) > 0) then
                            write(9,*) AllDoublesHist(i,j,k,l) / &
                                           real(iterations*nprocessors,dp)
                        endif
                    enddo
                enddo
            enddo
        enddo
 
        ! Test if all of the possible excitations have been generated by the
        ! random routine.
        if (bTestList) then
            do i=1,nexcit
                if (ex_list(i) .eqv. .false.) then
                    call write_det (6, nK(i,:), .true.)
                    call stop_all (this_routine, "Excitation in list not &
                                  &generated stochastically")
                endif
            enddo
        endif

        ! Clean up
        deallocate (SinglesHist, AllSinglesHist, DoublesHist, AllDoublesHist)
        if (allocated(nK)) deallocate(nK)
        if (allocated(ex_list)) deallocate(ex_list)
        call clean_excit_gen_store (store)
    end subroutine

    ! Insert the (already sorted) items of n2 into the sorted array
    ! n1 which only contains items up to position n1, but is large
    ! enough to contain list2
    subroutine int_list_merge (list1, list2, n1, n2)
        integer, intent(in) :: n1, n2, list2(1:n2)
        integer, intent(inout) :: list1(1:n1+n2)
        integer i, j, pos

        ! Work backwards through the lists, therefore never have to move
        ! items out of the way.
        ! j contains the reading position from list1
        ! pos contains the writing position into list1
        pos = n1+n2+1
        j = n1
        do i=n2,1,-1
            do j=j,1,-1
                pos = pos - 1
                if (list2(i) > list1(j)) then
                    list1(pos) = list2(i)
                    exit
                else
                    list1(pos) = list1(j)
                endif
            enddo
            if (j == 0) exit
        enddo

        list1(1:i) = list2(1:i)
    end subroutine

end module


! N.B. This is outside the module *sigh*
subroutine csf_sym_setup ()

    use SymExcitDataMod, only: ScratchSize, ScratchSize3
    implicit none

    ! We use the third scratch array to store the third class count.

    call SpinOrbSymSetup ()
    ScratchSize3 = ScratchSize / 2

end subroutine
    
