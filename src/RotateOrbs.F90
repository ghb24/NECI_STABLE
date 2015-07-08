module RotateOrbsMod

    use Global_utilities
    use Parallel_neci 
    use IntegralsData, only: UMAT, nFrozen, ChemPot
    use UMatCache, only: UMatInd
    use constants, only: dp, PI
    use SystemData, only: ConvergedForce, TimeStep, tLagrange, tShake, tShakeApprox, ShakeConverged
    use SystemData, only: tROIteration, ROIterMax, tShakeIter, ShakeIterMax, OrbEnMaxAlpha
    use SystemData, only: G1, ARR, NEl, nBasis, LMS, ECore, tSeparateOccVirt, Brr, nBasisMax, OrbOrder
    use SystemData, only: lNoSymmetry, tRotatedOrbs, tERLocalization, tRotateOccOnly
    use SystemData, only: tOffDiagMin, DiagWeight, OffDiagWeight, tRotateVirtOnly, tOffDiagSqrdMax
    use SystemData, only: tOffDiagSqrdMin, tOffDiagMax, tDoubExcMin, tOneElIntMax, tOnePartOrbEnMax
    use SystemData, only: tShakeDelay, ShakeStart, tVirtCoulombMax, tVirtExchangeMin, MaxMinFac
    use SystemData, only: tMaxHLGap, tHijSqrdMin, OneElWeight, DiagMaxMinFac, OneElMaxMinFac
    use SystemData, only: tDiagonalizehij, tHFSingDoubExcMax, tSpinOrbs, tReadInCoeff, tUseMP2VarDenMat
    use SystemData, only: tStoreSpinOrbs, tROHF, tFindCINatOrbs, tUseHFOrbs, tUEG
    use LoggingData, only: tROHistogramAll, tROFciDump, tROHistER, tROHistOffDiag, tROHistDoubExc, tPrintRODump
    use LoggingData, only: tROHistSingExc, tROHistOnePartOrbEn, tROHistOneElInts, tROHistVirtCoulomb
    use LoggingData, only: tPrintInts, tTruncRODump, NoTruncOrbs, NoDumpTruncs, tTruncDumpbyVal, TruncEvalues, tWriteTransMat
    use OneEInts, only: TMAT2D
    use SymData, only: TwoCycleSymGens, SymLabelList, SymLabelCounts
    use Timing_neci, only: end_timing, print_timing_report
    use Soft_exit, only: test_SOFTEXIT
    use RotateOrbsData 
    use sort_mod
    use util_mod, only: get_free_unit

    implicit none

    integer, allocatable :: Lab(:,:), LabVirtOrbs(:), LabOccOrbs(:), SymLabelList3_rotInv(:)
    real(dp), allocatable :: CoeffCorT2(:,:), CoeffUncorT2(:,:)
    real(dp), allocatable :: Lambdas(:,:), ArrNew(:,:), ArrDiagNew(:), TMAT2DTemp(:,:), TMAT2DRot(:,:), TMAT2DPartRot01(:,:)
    real(dp), allocatable :: TMAT2DPartRot02(:,:)
    real(dp), allocatable :: DerivCoeff(:,:), UMATTemp01(:,:,:,:), UMATTemp02(:,:,:,:)
    real(dp), allocatable :: DerivLambda(:,:), ForceCorrect(:,:), Correction(:,:), ShakeLambdaNew(:), ConstraintCor(:)
    real(dp), allocatable :: Constraint(:), ShakeLambda(:), DerivConstrT1(:,:,:), DerivConstrT2(:,:,:), DerivConstrT1T2(:,:)
    real(dp), allocatable :: DerivConstrT1T2Diag(:), FourIndInts(:,:,:,:)
    real(dp), allocatable :: TwoIndInts01(:,:,:,:), TwoIndInts02(:,:,:,:), ThreeIndInts01(:,:,:,:), FourIndInts02(:,:,:,:)
    real(dp), allocatable :: ThreeIndInts02(:,:,:,:), ThreeIndInts03(:,:,:,:), ThreeIndInts04(:,:,:,:)  
    real(dp), allocatable :: DiagTMAT2Dfull(:), TMAT2DNew(:,:) 
    real(dp), allocatable :: TwoIndIntsER(:,:,:), ThreeIndInts01ER(:,:), ThreeIndInts02ER(:,:), FourIndIntsER(:)
    integer(TagIntType) :: TwoIndIntsERTag, ThreeIndInts01ERTag, ThreeIndInts02ERTag, FourIndIntsERTag
    integer(TagIntType) :: TwoIndInts01Tag, TwoIndInts02Tag, ThreeIndInts01Tag, ThreeIndInts02Tag, ThreeIndInts03Tag
    integer(TagIntType) :: FourIndInts02Tag, ThreeIndInts04Tag, UMATTemp02Tag
    integer(TagIntType) :: TMAT2DTempTag, TMAT2DRotTag, TMAT2DPartRot01Tag, TMAT2DPartRot02Tag
    integer(TagIntType) :: LabTag, ForceCorrectTag, CorrectionTag, FourIndIntsTag, ArrDiagNewTag, ArrNewTag, UMATTemp01Tag
    integer :: ShakeIterInput, NoOcc, LowBound02, HighBound02, Iteration, TotNoConstraints
    integer(TagIntType) :: CoeffCorT2Tag, CoeffUncorT2Tag, LambdasTag, DerivCoeffTag, DerivLambdaTag
    integer(TagIntType) :: ShakeLambdaNewTag
    integer(TagIntType) :: ShakeLambdaTag, ConstraintTag, ConstraintCorTag, DerivConstrT1Tag, DerivConstrT2Tag, DerivConstrT1T2Tag
    integer(TagIntType) :: DerivConstrT1T2DiagTag
    integer(TagIntType) :: LabVirtOrbsTag, LabOccOrbsTag
    integer :: MinOccVirt, MaxOccVirt, MinMZ, MaxMZ, error, LowBound, HighBound
    integer :: NoInts01, NoInts02, NoInts03, NoInts04, NoInts05, NoInts06
    integer(TagIntType) :: DiagTMAT2DfullTag, TMAT2DNewTag, SymLabelList3_rotInvTag
    logical :: tNotConverged, tInitIntValues
    real(dp) :: OrthoNorm, ERPotEnergy, HijSqrdPotEnergy, OffDiagPotEnergy, CoulPotEnergy, PotEnergy, Force, TwoEInts, DistCs
    real(dp) :: OrthoForce, DistLs, LambdaMag, PEInts, PEOrtho
    real(dp) :: ForceInts, TotCorrectedForce
    real(dp) :: ijOccVirtPotEnergy, EpsilonMin, MaxTerm
    real(dp) :: DiagOneElPotInit, ERPotInit, ijVirtOneElPotInit, ijVirtCoulPotInit, ijVirtExchPotInit
    real(dp) :: singCoulijVirtInit, singExchijVirtInit, singCoulconHFInit, singExchconHFInit, ijklPotInit, ijklantisymPotInit
    real(dp) :: ijOccVirtOneElPotInit, ijOccVirtCoulPotInit, ijOccVirtExchPotInit
    real(dp) :: OrthoFac = 1.0_dp, ROHistSing(2, 4002), ROHistOffDiag(2, 4002), ROHistDoubExc(2, 4002), ROHistER(2, 4002)
    real(dp) :: ROHistHijVirt(2, 4002), ROHistHijOccVirt(2, 4002), ROHistHii(2, 4002)
    real(dp) :: ROHistOnePartOrbEn(2, 4002), ROHistDCijOcklVir(2, 4002), ROHistDEijOcklVir(2, 4002), ROHistDCijklVir(2, 4002)
    real(dp) :: ROHistDEijklVir(2, 4002)
    real(dp) :: ROHistSCikOcjVir(2, 4002), ROHistSEikOcjVir(2, 4002), ROHistSCkOcijVir(2, 4002), ROHistSEkOcijVir(2, 4002)
    real(dp) :: ROHistSCijkVir(2, 4002), ROHistSEijkVir(2, 4002)
    real(dp) :: ROHistSASikOcjVir(2, 4002), ROHistSASkOcijVir(2, 4002), ROHistSASijkVir(2, 4002), ROHistASijklVir(2, 4002)
    real(dp) :: ROHistASijOcklVir(2, 4002)
    type(timer), save :: Rotation_Time, FullShake_Time, Shake_Time, Findtheforce_Time, Transform2ElInts_Time
    type(timer), save :: findandusetheforce_time, CalcDerivConstr_Time, TestOrthoConver_Time
    type(timer), save :: RefillUMAT_Time, PrintROFCIDUMP_Time
! In this routine, alpha (a), beta (b), gamma (g) and delta (d) refer to the unrotated (HF) orbitals where 
!possible such that < a b | g d > is an unrotated four index integral.   
! For the rotated orbitals, the letter i, j, k and l are generally used, i.e. < i j | k l > refers to 
!a transformed four index integral.
! Differentiation of the potential energy (to find the force) is done with respect to coefficient 
!c(z,m) (or c(a,m)), where zeta (z) or a refers to the HF index, and m to the rotated.
    
contains

    subroutine RotateOrbs()

        if (iProcIndex == Root) then

! If we are reading in our own transformation matrix (coeffT1) don't need a lot of the initialisation stuff.
            if (tReadInCoeff.or.tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs) then

                tNotConverged = .false.
                call FindNatOrbitals()

            else
! Need to actually find the coefficient matrix and then use it.

                tNotConverged = .true.
                call InitLocalOrbs()        ! Set defaults, allocate arrays, write out headings 
                                            ! for OUTPUT, set integarals to HF values.

                if (tDiagonalizehij) then

                    call Diagonalizehij()
                    tNotConverged = .false.
                    Iteration = 2

                elseif (tMaxHLGap) then
                    call EquateDiagFock()
                    tNotConverged = .false.

                else

                    tNotConverged = .true.

                    call WriteStats()           ! write out the original stats before any rotation.
                   
                    call set_timer(Rotation_Time, 30)

                    do while(tNotConverged)     ! rotate the orbitals until the sum of the four index 
                                                ! integral falls below a chose convergence value.

                        Iteration = Iteration+1
                        
                        call FindNewOrbs()      ! bulk of the calculation.
                                                ! do the actual transformations, moving the coefficients by 
                                                !a timestep according to the calculated force. 

                        call WriteStats()       ! write out the stats for this iteration.

                    end do           

                    call halt_timer(Rotation_Time)
                    
                    write(6,*) "Convergence criterion met. Finalizing new orbitals..."

                end if

! Make symmetry, orbitals, one/two-electron integrals consistent with rest of NECI
                call FinalizeNewOrbs()

                call writeBASIS(6, G1, nBasis, ARR, BRR)

                call DeallocateMem()

            end if            

            call neci_flush(6)
            call neci_flush(transform_unit)
        end if


    end subroutine RotateOrbs

    subroutine FindNatOrbitals()

! This routine simply takes a transformation matrix and rotates the integrals to produce a new FCIDUMP file.
! In one case the transformation matrix is read in from a file TRANSFORMMAT.
! In the other, the transformation matrix is calculated from the MP2 variational density matrix.

! MP2VDM = D2_ab = sum_ijc [ t_ij^ac ( 2 t_ij^bc - t_ji^bc ) ]
! Where :  t_ij^ac = - < ab | ij > / ( E_a - E_i + E_b - Ej )
! Ref : J. Chem. Phys. 131, 034113 (2009) - note: in Eqn 1, the cb indices are the wrong way round (should be bc).

        use NatOrbsMod, only: SetUpNatOrbLabels, FindNatOrbs, FillCoeffT1, DeallocateNatOrbs, PrintOccTable
        integer :: i, a, ierr, MinReadIn, MaxReadIn, iunit
        character(len=*), parameter :: this_routine = 'FindNatOrbitals'

        if (tUseMP2VarDenMat) write(6,*) '*** Transforming the HF orbitals into the MP2 approximate natural orbitals. ***'
        if (tFindCINatOrbs) then
            write(6,*) '*** Transforming the HF orbitals into approximate natural orbitals'
            write(6,*) 'based on the one-electron density matrix found from the wavefunction calculated above. ***'
        end if

        if (tSpinOrbs) then
            if (.not.tStoreSpinOrbs) then
                write(6,*) "We want to use spin orbitals - turning on tStoreSpinOrbs."
                tStoreSpinOrbs = .true.
            end if
        end if

        if (tROHF .and. tStoreSpinOrbs) call Stop_All(this_routine,"Cannot compress open shell systems into spatial " &
            & //"orbitals when rotating, turn off ROHF.")

        if (tTruncRODump .and. (.not.tTruncDumpbyVal)) then 
            NoFrozenVirt = NoTruncOrbs(1)
        elseif (tTruncRODump) then
            ! If the 'number of frozen orbitals' is given as a cutoff - take NoFrozenVirt to be 0 
            !for all the allocation purposes - will set this later when
            ! we have the eigenvalues and know how many orbitals lie below it.
            NoFrozenVirt = 0
            TruncEval = TruncEvalues(1)
        else
            NoFrozenVirt = 0
        end if

        SpatOrbs = nBasis/2
        if (tStoreSpinOrbs) then
            NoOrbs = nBasis
            NoOcc = NEl
            MinReadIn = 1
            MaxReadIn = nBasis
            if (tRotateVirtOnly) MinReadIn = NEl+1
            if (tRotateOccOnly) MaxReadIn = NEl
            ! If tStoreSpinOrbs ARR(:,2) is not filled, but we want to use it later, so just fill it here.
            do i = 1, NoOrbs
                ARR(BRR(i), 2) = ARR(i,1)
            end do
            allocate(SymLabelCounts2_rot(2, 32), stat = ierr)
            call LogMemAlloc('SymLabelCounts2_rot', 2*32, 4, this_routine, SymLabelCounts2_rotTag, ierr)
            SymLabelCounts2_rot(:,:) = 0
            ! first 8 refer to the occupied, and the second to the virtual beta spin.
            ! third and fourth to the occupied and virtual alpha spin.
 
        else
            NoOrbs = SpatOrbs
            NoOcc = NEl/2
            MinReadIn = 1
            MaxReadIn = SpatOrbs
            if (tRotateVirtOnly) MinReadIn = (NEl/2)+1
            if (tRotateOccOnly) MaxReadIn = NEl/2
            allocate(SymLabelCounts2_rot(2,16), stat = ierr)
            call LogMemAlloc('SymLabelCounts2_rot', 2*16, 4, this_routine, SymLabelCounts2_rotTag, ierr)
            SymLabelCounts2_rot(:,:) = 0
            ! first 8 refer to the occupied, and the second to the virtual.

        end if
        NoRotOrbs = NoOrbs

        call ApproxMemReq()

        allocate(SymLabelList2_rot(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelList2_rot', NoOrbs, 4, this_routine, SymLabelList2_rotTag, ierr)
        SymLabelList2_rot(:) = 0                     
        allocate(SymLabelList3_rot(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelList3_rot', NoOrbs, 4, this_routine, SymLabelList3_rotTag, ierr)
        SymLabelList3_rot(:) = 0                     
 
        allocate(SymLabelListInv_rot(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelListInv_rot', NoOrbs, 4, this_routine, SymLabelListInv_rotTag, ierr)
        SymLabelListInv_rot(:) = 0                     


        if (tReadInCoeff.or.tUseHFOrbs) then
! No symmetry, so no reordering of the orbitals - symlabellist just goes from 1-NoOrbs.        
! When we are just reading in the coefficients and transforming, it does not matter about the ordering of the orbitals.
            do i = 1, NoOrbs
                SymLabelList2_rot(i) = i
                SymLabelListInv_rot(i) = i
            end do

        elseif (tFindCINatOrbs.or.tUseMP2VarDenMat) then

            call SetupNatOrbLabels() 

        end if

! Yet another labelling system, SymLabelList3_rot is created here.
! This indicates the label of the transformed orbital.
! In the case where we are truncating the space, the transformed orbitals are ordered according 
!to the size of the eigenvalues of the MP2VDM
! matrix when it is diagonalised.  We wish to keep them in this order when transforming 
!the integrals etc, so that when we truncate the last 
! NoFrozenVirt orbitals, we are removing those with the smallest MP2VDM eigenvalues (occupation numbers).
! In the case where no truncation is made however, SymLabelList3_rot is the same as SymLabelList2_rot, 
!so that the indexes remain the same as previously.
! This allows for the option of going straight into a spawning calc from the rotation, which is not 
!possible when a truncation is performed 
! because of the messed up indices.
        if (tTruncRODump) then
            if (MOD(NoFrozenVirt, 2) /= 0) call Stop_All(this_routine,"Must freeze virtual spin orbitals in pairs of 2.")
            if (tStoreSpinOrbs) then
                NoRotOrbs = NoOrbs-NoFrozenVirt
            else
                NoFrozenVirt = NoFrozenVirt/2
                NoRotOrbs = NoOrbs-NoFrozenVirt
            end if            
            do i = 1, NoOrbs
                SymLabelList3_rot(i) = i
            end do
        else
            do i = 1, NoOrbs
                SymLabelList3_rot(i) = SymLabelList2_rot(i)
            end do
        end if

        allocate(CoeffT1(NoOrbs, NoRotOrbs), stat = ierr)
        call LogMemAlloc(this_routine, NoRotOrbs*NoOrbs, 8, this_routine, CoeffT1Tag, ierr)
        CoeffT1(:,:) = 0.0_dp
        if (tSeparateOccVirt) then
            do i = 1, NoRotOrbs
                CoeffT1(i,i) = 1.0_dp
            end do
        end if

        if (tUEG) then

            call FindNatOrbs()

            call FillCoeffT1()

        else
            if (tReadInCoeff) then

                write(6,'(A)') " Reading in the transformation matrix from TRANSFORMMAT, and using this to rotate the HF orbitals."

                iunit  =  get_free_unit()
                open(iunit, file='TRANSFORMMAT', status='old')
                do i = 1, NoOrbs
                    do a = 1, NoOrbs
                        READ(iunit,*) CoeffT1(a,i)
                    end do
                end do
                close(iunit)
          
            elseif (tFindCINatOrbs.or.tUseMP2VarDenMat.or.tUseHFOrbs) then
                
                if (.not.tUseHFOrbs) call FindNatOrbs()
                
                if (tUseHFOrbs) then
                    call PrintOccTable()
                else
                    call FillCoeffT1()
                end if

            end if

            if (tPrintRODump) then
                allocate(FourIndInts(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
                call LogMemAlloc('FourIndInts',(NoOrbs**4), 8, this_routine, FourIndIntsTag, ierr)

                ! Then, transform2ElInts
                write(6,*) 'Transforming the four index integrals'
                call Transform2ElIntsMemSave()

                write(6,*) 'Re-calculating the fock matrix'
                call CalcFOCKMatrix()

                write(6,*) 'Refilling the UMAT and TMAT2D'
                ! The ROFCIDUMP is also printed out in here.        
                call RefillUMATandTMAT2D()        

                call neci_flush(6)

                if ((tFindCINatOrbs.or.tUseMP2VarDenMat) .and. (NoDumpTruncs > 1)) call ReTruncROFciDump()

                if ((.not.tUseHFOrbs) .and. (.not.tReadInCoeff)) call DeallocateNatOrbs()
            end if

            if (tWriteTransMat) call WriteTransformMat()
     
! If a truncation is being made, the new basis will not be in the correct energetic ordering - this does not matter, as we
! never go straight into a spawning and they will be reordered when the ROFCIDUMP file is read in again. 
            call writeBASIS(6, G1, nBasis, ARR, BRR)

            deallocate(CoeffT1)
            call LogMemDeAlloc(this_routine, CoeffT1Tag)
            deallocate(SymLabelList2_rot)
            call LogMemDeAlloc(this_routine, SymLabelList2_rotTag)
            deallocate(SymLabelListInv_rot)
            call LogMemDeAlloc(this_routine, SymLabelListInv_rotTag)
            if (tPrintRODump) then
                deallocate(FourIndInts)
                call LogMemDeAlloc(this_routine, FourIndIntsTag)
            end if
        end if

    end subroutine FindNatOrbitals 

    subroutine ReTruncROFciDump()

        use NatOrbsMod, only: FillCoeffT1

        integer :: i, j, ierr
        character(len=*), parameter :: this_routine = 'ReTruncROFciDump'

        do i = 2, NoDumpTruncs

            deallocate(ArrDiagNew)
            call LogMemDeAlloc(this_routine, ArrDiagNewTag)
            deallocate(CoeffT1)
            call LogMemDeAlloc(this_routine, CoeffT1Tag)
            deallocate(FourIndInts)
            call LogMemDeAlloc(this_routine, FourIndIntsTag)
            deallocate(SymOrbs_rot)
            call LogMemDeAlloc(this_routine, SymOrbs_rotTag)
            deallocate(TMAT2DNew)
            call LogMemDeAlloc(this_routine, TMAT2DNewTag)
            deallocate(EvaluesTrunc)
            call LogMemDeAlloc(this_routine, EvaluesTruncTag)


            if (tTruncDumpbyVal) then
                NoFrozenVirt = 0
                TruncEval = TruncEvalues(i)
            else
                if (tStoreSpinOrbs) then
                    NoFrozenVirt = NoTruncOrbs(i)
                else
                    NoFrozenVirt = NoTruncOrbs(i)/2
                end if            
            end if
            NoRotOrbs = NoOrbs-NoFrozenVirt
 
            if (MOD(NoFrozenVirt, 2) /= 0) call Stop_All(this_routine,"Must freeze virtual spin orbitals in pairs of 2.")

            allocate(CoeffT1(NoOrbs, NoRotOrbs), stat = ierr)
            call LogMemAlloc(this_routine, NoRotOrbs*NoOrbs, 8, this_routine, CoeffT1Tag, ierr)
            CoeffT1(:,:) = 0.0_dp
            if (tSeparateOccVirt) then
                do j = 1, NoRotOrbs
                    CoeffT1(i,i) = 1.0_dp
                end do
            end if

            call FillCoeffT1()


            allocate(FourIndInts(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('FourIndInts',(NoOrbs**4), 8, this_routine, FourIndIntsTag, ierr)

            ! Then, transform2ElInts.
            write(6,*) 'Transforming the four index integrals.'
            call Transform2ElIntsMemSave()

            write(6,*) 'Re-calculating the fock matrix.'
            call CalcFOCKMatrix()

            write(6,*) 'Refilling the UMAT and TMAT2D.'
            ! The ROFCIDUMP is also printed out in here.        
            call RefillUMATandTMAT2D()        

            call neci_flush(6)

        end do

    end subroutine ReTruncROFciDump

    subroutine ApproxMemReq()

        ! This routine makes a quick sum of the memory that will be require to
        ! transform the integrals from the HF to the new basis.

        ! Main arrays required are:
        MemAllocRot = 0
        
        ! Symmetry/Labelling:
        !   - SymLabelLists(NoOrbs) x 3 
        !   - SymLabelCounts(32/16 - Spin/Spat)
        MemAllocRot = MemAllocRot+(3*NoOrbs*4)
        MemAllocRot = MemAllocRot+(32*4)
        
        ! Finding transformation matrices
        !   - NatOrbsMat(NoOrbs, NoOrbs)
        !   - Evalues(NoOrbs) x 2
        MemAllocRot = MemAllocRot+((NoOrbs**2)*8)
        MemAllocRot = MemAllocRot+(2*NoOrbs*8)

        ! Transformation of integrals
        !   - CoeffT1(NoOrbs, NoRotOrbs) 
        !   - FourIndInts(NoRotOrbs, NoRotOrbs, NoOrbs, NoOrbs)
        !   - Temp4indints(NoRotOrbs, NoOrbs)
        if (tPrintRODump) then
            MemAllocRot = MemAllocRot+(NoOrbs*NoRotOrbs*8*2)
            MemAllocRot = MemAllocRot+((NoRotOrbs**2)*(NoOrbs**2)*8)

                ! Transform fock
                !   - ArrNew(NoOrbs) - reduce this?
                MemAllocRot = MemAllocRot+(NoOrbs*8)

                ! RefillTMAT2D
                !   - TMAT2D(nBasis, nBasis) 
                MemAllocRot = MemAllocRot+((nBasis**2)*8)
            end if

            write(6,'(A72, F20.10, A15)') "Rough estimate of the memory required for the orbital transformation  =  ", &
               real(MemAllocRot, dp)/1048576.0_dp," Mb/Processor"

    end subroutine ApproxMemReq

    subroutine WriteTransformMat()

        integer :: w, x, i, a, b, iunit

! This file is printed to be used to produce cube files from QChem.
! Line 1 is the coefficients of HF spatial orbitals 1 2 3 ... which form transformed orbital 1 etc.

        iunit  =  get_free_unit()
        open(iunit, file='MOTRANSFORM', FORM = 'UNFORMATTED', access = 'direct', recl = 8)
        ! Need to put this back into the original order. 

        x  =  0
        if (tStoreSpinOrbs) then
            do i = 1, NoOrbs-1, 2
                ! SymLabelList2_rot(i) gives the orbital label (from Dalton or QChem) corresponding to our
                ! label i.
                ! SymLabelListInv_rot(j) therefore gives the label used in CoeffT1 corresponding to the
                ! Qchem/Dalton label j.
                    
                do a = 1, NoOrbs-1, 2
                    b = SymLabelListInv_rot(a)
                    write(iunit, rec = x) CoeffT1(b,i)
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                end do
            end do
            do i = 2, NoOrbs, 2
                do a = 2, NoOrbs, 2
                    b = SymLabelListInv_rot(a)
                    write(iunit, rec = x) CoeffT1(b,i)
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                end do
            end do
        else
            w = 1
            x = 1   !keep a counter of record number
            do while (w <= 2)
                do i = 1, SpatOrbs
                    ! SymLabelList2_rot(i) gives the orbital label (from Dalton or QChem) corresponding to our
                    ! label i.
                    ! SymLabelListInv_rot(j) therefore gives the label used in CoeffT1 corresponding to the
                    ! Qchem/Dalton label j.
                    do a = 1, SpatOrbs
                        b = SymLabelListInv_rot(a)
                        write(iunit, rec = x) CoeffT1(b,i)
                        x = x+1
                        ! a/b are the original (HF) orbitals, and i/j the transformed
                    end do
                end do
                w = w+1
                ! print the whole matrix twice, once for alpha spin, once for beta.
            end do
        end if
        close(iunit)
 
        open(iunit, file='MOTRANSFORM02')
        ! Need to put this back into the original order. 
        w = 1
        x = 1   !keep a counter of record number
        do while (w <= 2)
            do i = 1, SpatOrbs
                ! SymLabelList2_rot(i) gives the orbital label (from Dalton or QChem) corresponding to our
                ! label i.
                ! SymLabelListInv_rot(j) therefore gives the label used in CoeffT1 corresponding to the
                ! Qchem/Dalton label j.
                do a = 1, SpatOrbs
                    b = SymLabelListInv_rot(a)
                    write(iunit,'(F20.10)', advance='no') CoeffT1(b,i)
                    x = x+1
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                end do
                write(iunit,*) ''
            end do
            w = w+1
            ! print the whole matrix twice, once for alpha spin, once for beta.
        end do
        close(iunit)

        open(iunit, file='TRANSFORMMAT', status='unknown')
        do i = 1, NoOrbs
            do a = 1, NoOrbs
                b = SymLabelListInv_rot(a)
                write(iunit,*) CoeffT1(b,i)
            end do
        end do
        call neci_flush(iunit)
        close(iunit)

    end subroutine WriteTransformMat

    subroutine InitLocalOrbs()

        character(len=*), parameter :: this_routine = 'InitLocalOrbs'
        integer :: ierr

        ! Writing to output which PE is being maximised/minimised.
        write(6,*) '*****'
        if (tERLocalization) then
            write(6,*) "Calculating new molecular orbitals based on Edmiston-Reudenberg localisation,"
            write(6,*) "i.e. maximisation of the <ii|ii> integrals..."
            write(6,*) "*****"
        end if
        if (tVirtCoulombMax) then
            write(6,*) "Calculating new molecular orbitals based on maximisation of the sum of the"
            write(6,*) "<ij|ij> integrals, where i and j are both virtuals..."
            write(6,*) "*****"
        end if
        if (tOffDiagSqrdMin) then
            write(6,*) "Calculating new molecular orbitals based on mimimisation "
            write(6,*) "of <ij|kl>^2 integrals..."
            write(6,*) "*****"
        end if
        if (tOffDiagMin) then
            write(6,*) "Calculating new molecular orbitals based on mimimisation "
            write(6,*) "of <ij|kl> integrals..."
            write(6,*) "*****"
        end if
        if (tDoubExcMin) then
            write(6,*) "Calculating new molecular orbitals based on mimimisation "
            write(6,*) "of the double excitation hamiltonian elements."
            write(6,*) "*****"
        end if
        if (tOnePartOrbEnMax) then
            write(6,*) "Calculating new molecular orbitals based on maximisation "
            write(6,*) "of the virtual one particle orbital energies."
            write(6,*) "*****"
        elseif (tMaxHLGap) then
            ! This will transform all the orbitals within a particlar group to
            ! have the same diagonal fock matrix element.
            write(6,*) "Transforming orbitals based on equating their diagonal fock matrix elements."
            write(6,*) "*****"
        end if

        ! Writing out which orthonormalisation method is being used...
        if (tLagrange) then
            if (tShake) then
                call neci_flush(6)
                call Stop_All(this_routine,"ERROR. Both LAGRANGE and SHAKE keywords present in the input. &
                & These two orthonormalisation methods clash.")
            end if
            write(6,*) "Using a Lagrange multiplier to attempt to rotate orbitals in a way to maintain orthonormality"
        elseif (tShake) then
            write(6,*) "Using the shake algorithm to iteratively find lambdas which maintain "
            write(6,*) "orthonormalisation with rotation"
        else
            write(6,*) "Explicity reorthonormalizing orbitals after each rotation."
        end if
        
        ! Check for a few possible errors.
        if (.not.TwoCycleSymGens) then
            call neci_flush(6)
            call Stop_All(this_routine,"ERROR. TwoCycleSymGens is false.  Symmetry is not abelian.") 
        end if
        if ((tRotateOccOnly.or.tRotateVirtOnly) .and. (.not.tSeparateOccVirt)) then
            tSeparateOccVirt = .true.
            write(6,*) "NOTE. Cannot rotate only occupied or virtual without first separating them."
            write(6,*) "SEPARATEOCCVIRT keyword is being turned on."
        end if        
        if ((tOffDiagSqrdMax .and. tOffDiagSqrdMin).or.(tOffDiagMax .and. tOffDiagMin)) then
            call neci_flush(6)
            call Stop_All(this_routine,"ERROR. Cannot both maximise and minimise off diagonal elements simultaneously")
        end if
        if (tOnePartOrbEnMax .and. (.not.tSeparateOccVirt)) then
            call neci_flush(6)
            call Stop_All(this_routine, &
            "ERROR. Cannot currently maximise the one particle orbital energies without separating occupied and virtual.") 
        end if
        write(6,*) "*****"

        ! Zero values.
        OrthoNorm = 0.0_dp
        ERPotEnergy = 0.0_dp
        PotEnergy = 0.0_dp
        Force = 0.0_dp
        TwoEInts = 0.0_dp
        PEInts = 0.0_dp
        PEOrtho = 0.0_dp
        ForceInts = 0.0_dp
        DistCs = 0.0_dp
        DistLs = 0.0_dp
        LambdaMag = 0.0_dp
        SpatOrbs = nBasis/2
        if (tStoreSpinOrbs) then
            NoOrbs = nBasis
            NoOcc = NEl
        else
            NoOrbs = SpatOrbs
            NoOcc = NEl/2
        end if
        NoRotOrbs = NoOrbs
        Iteration = 0
        OrthoForce = 0.0_dp
        ShakeIterInput = ShakeIterMax
        TotNoConstraints = (NoOrbs*(NoOrbs+1))/2

        ! When maximising the one particle orbital energies, choose the zero
        ! value (Epsilon min).
        if (tRotateVirtOnly .and. tOnePartOrbEnMax) then
            EpsilonMin = ARR(NEl+1, 1)
            write(6,*) 'Taking EpsilonMin to be the LUMO of the HF orbitals...'
            write(6,*) 'EpsilonMin  =  ', EpsilonMin
        elseif (tOnePartOrbEnMax) then
            EpsilonMin = ChemPot
            write(6,*) 'Taking EpsilonMin to be the chemical potential (midway between HF HOMO and LUMO)...'
            write(6,*) 'therefore EpsilonMin  =  ', EpsilonMin
        end if

        ! Set timed routine names.
        Rotation_Time%timer_name = 'RotateTime'
        Shake_Time%timer_name = 'ShakeTime'
        FullShake_Time%timer_name = 'FullShakeTime'
        Findtheforce_Time%timer_name = 'FindtheForceTime'
        Transform2ElInts_Time%timer_name = 'Transform2ElIntsTime'
        findandusetheforce_time%timer_name = 'Findandusetheforce'
        CalcDerivConstr_Time%timer_name = 'CalcDerivConstr'
        TestOrthoConver_Time%timer_name = 'TestOrthoConver'
       
        ! Allocate memory.

        allocate(CoeffT1(NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('CoeffT1', NoOrbs**2, 8, this_routine, CoeffT1Tag, ierr)
        allocate(CoeffCorT2(NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('CoeffCorT2', NoOrbs**2, 8, this_routine, CoeffCorT2Tag, ierr)
        allocate(CoeffUncorT2(NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('CoeffUncT2', NoOrbs**2, 8, this_routine, CoeffUncorT2Tag, ierr)
        CoeffUncorT2(:,:) = 0.0_dp
         
        allocate(DerivCoeff(NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('DerivCoeff', NoOrbs**2, 8, this_routine, DerivCoeffTag, ierr)
  
        allocate(DiagTMAT2Dfull(NoOrbs-(NoOcc)), stat = ierr)
        call LogMemAlloc('DiagTMAT2Dfull',(NoOrbs-(NoOcc)), 8, this_routine, DiagTMAT2DfullTag, ierr)
        allocate(UMATTemp01(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('UMATTemp01', NoOrbs**4, 8, this_routine, UMATTemp01Tag, ierr)
        allocate(TwoIndInts01(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('TwoIndInts01', NoOrbs**4, 8, this_routine, TwoIndInts01Tag, ierr)
        allocate(ThreeIndInts02(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('ThreeIndInts02', NoOrbs**4, 8, this_routine, ThreeIndInts02Tag, ierr)
        allocate(FourIndInts(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('FourIndInts', NoOrbs**4, 8, this_routine, FourIndIntsTag, ierr)
        allocate(FourIndInts02(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('FourIndInts02', NoOrbs**4, 8, this_routine, FourIndInts02Tag, ierr)

        ! Partially transformed temporary arrays.
        if (tERLocalization .and. (.not.tStoreSpinOrbs)) then
            allocate(TwoIndIntsER(NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('TwoIndIntsER', NoOrbs**3, 8, this_routine, TwoIndIntsERTag, ierr)
            allocate(ThreeIndInts01ER(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ThreeIndInts01ER', NoOrbs**2, 8, this_routine, ThreeIndInts01ERTag, ierr)
            allocate(ThreeIndInts02ER(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ThreeIndInts02ER', NoOrbs**2, 8, this_routine, ThreeIndInts02ERTag, ierr)
            allocate(FourIndIntsER(NoOrbs), stat = ierr)
            call LogMemAlloc('FourIndIntsER', NoOrbs, 8, this_routine, FourIndIntsERTag, ierr)
        else
            allocate(TMAT2DTemp(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('TMAT2DTemp', NoOrbs**2, 8, this_routine, TMAT2DTempTag, ierr)
            allocate(TMAT2DPartRot01(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('TMAT2DPartRot01', NoOrbs**2, 8, this_routine, TMAT2DPartRot01Tag, ierr)
            allocate(TMAT2DPartRot02(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('TMAT2DPartRot02', NoOrbs**2, 8, this_routine, TMAT2DPartRot02Tag, ierr)
            allocate(TMAT2DRot(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('TMAT2DRot', NoOrbs**2, 8, this_routine, TMAT2DRotTag, ierr)
            allocate(UMATTemp02(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('UMATTemp02', NoOrbs**4, 8, this_routine, UMATTemp02Tag, ierr)
 
            ! Partially transformed combined arrays.
            allocate(TwoIndInts02(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('TwoIndInts02', NoOrbs**4, 8, this_routine, TwoIndInts02Tag, ierr)
            allocate(ThreeIndInts01(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ThreeIndInts01', NoOrbs**4, 8, this_routine, ThreeIndInts01Tag, ierr)
            allocate(ThreeIndInts03(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ThreeIndInts03', NoOrbs**4, 8, this_routine, ThreeIndInts03Tag, ierr)
            allocate(ThreeIndInts04(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ThreeIndInts04', NoOrbs**4, 8, this_routine, ThreeIndInts04Tag, ierr)
        end if

        ! Allocate according to orthonormalisation method being used.
        if (tLagrange) then
            allocate(Lambdas(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('Lambdas', NoOrbs**2, 8, this_routine, LambdasTag, ierr)
            allocate(DerivLambda(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('DerivLambda', NoOrbs**2, 8, this_routine, DerivLambdaTag, ierr)
            Lambdas(:,:) = 0.0_dp
            DerivLambda(:,:) = 0.0_dp
        end if

        if (tShake) then
            allocate(ShakeLambda(TotNoConstraints), stat = ierr)
            call LogMemAlloc('ShakeLambda', TotNoConstraints, 8, this_routine, ShakeLambdaTag, ierr)
            ShakeLambda(:) = 0.0_dp                     
            allocate(ShakeLambdaNew(TotNoConstraints), stat = ierr)
            call LogMemAlloc('ShakeLambdaNew', TotNoConstraints, 8, this_routine, ShakeLambdaNewTag, ierr)
            ShakeLambdaNew(:) = 0.0_dp                     
            allocate(Constraint(TotNoConstraints), stat = ierr)
            call LogMemAlloc('Constraint', TotNoConstraints, 8, this_routine, ConstraintTag, ierr)
            allocate(ConstraintCor(TotNoConstraints), stat = ierr)
            call LogMemAlloc('ConstraintCor', TotNoConstraints, 8, this_routine, ConstraintCorTag, ierr)
            allocate(DerivConstrT1(NoOrbs, NoOrbs, TotNoConstraints), stat = ierr)
            call LogMemAlloc('DerivConstrT1', NoOrbs*TotNoConstraints*NoOrbs, 8, this_routine, DerivConstrT1Tag, ierr)
            DerivConstrT1(:,:,:) = 0.0_dp
            allocate(DerivConstrT2(NoOrbs, NoOrbs, TotNoConstraints), stat = ierr)
            call LogMemAlloc('DerivConstrT2', NoOrbs*TotNoConstraints*NoOrbs, 8, this_routine, DerivConstrT2Tag, ierr)
            allocate(ForceCorrect(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ForceCorrect', NoOrbs**2, 8, this_routine, ForceCorrectTag, ierr)
            allocate(Correction(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('Correction', NoOrbs**2, 8, this_routine, CorrectionTag, ierr)
            if (tShakeApprox) then
                allocate(DerivConstrT1T2Diag(TotNoConstraints), stat = ierr)
                call LogMemAlloc('DerivConstrT1T2Diag', TotNoConstraints, 8, this_routine, DerivConstrT1T2DiagTag, ierr)
                DerivConstrT1T2Diag(:) = 0.0_dp
            else
                allocate(DerivConstrT1T2(TotNoConstraints, TotNoConstraints), stat = ierr)
                call LogMemAlloc('DerivConstrT1T2', TotNoConstraints**2, 8, this_routine, DerivConstrT1T2Tag, ierr)
            end if
        end if               

        ! Indexing arrays.
        allocate(SymLabelList2_rot(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelList2_rot', NoOrbs, 4, this_routine, SymLabelList2_rotTag, ierr)
        SymLabelList2_rot(:) = 0                     
        allocate(SymLabelList3_rot(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelList3_rot', NoOrbs, 4, this_routine, SymLabelList3_rotTag, ierr)
        SymLabelList3_rot(:) = 0                     
 
        allocate(SymLabelListInv_rot(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelListInv_rot', NoOrbs, 4, this_routine, SymLabelListInv_rotTag, ierr)
        SymLabelListInv_rot(:) = 0                     
   
        allocate(Lab(2, TotNoConstraints), stat = ierr)
        call LogMemAlloc('Lab', 2*TotNoConstraints, 4, this_routine, LabTag, ierr)
        Lab(:,:) = 0                     

        ! Do any initial calculations, and set up starting values for arrays
        ! used in rotation.
        call InitRotCalc()

        ! Write out the headings for the results file.
        transform_unit  =  get_free_unit()
        open(transform_unit, file='Transform', status='unknown')
        if (tLagrange) then
            write(transform_unit,"(A12, 11A18)") "# Iteration","2.PotEnergy","3.PEInts","4.PEOrtho","5.Force","6.ForceInts", &
            "7.OrthoForce","8.Sum<ij|kl>^2",&
                        &"9.OrthoNormCondition","10.DistMovedbyCs","11.DistMovedByLs","12.LambdaMag"
            write(6,"(A12, 11A19)") "Iteration","2.PotEnergy","3.PEInts","4.PEOrtho","5.Force","6.ForceInts","7.OrthoForce", &
                    "8.Sum<ij|kl>^2",&
                    &"9.OrthoNormCondition","10.DistMovedbyCs","11.DistMovedbyLs","12.LambdaMag"
        elseif (tERLocalization .and. tHijSqrdMin) then
            write(transform_unit,"(A12, 7A24)") "# Iteration","2.ERPotEnergy","3.HijSqrdPotEnergy","4.PotEnergy","5.Force", &
                "6.Totalcorrforce","7.OrthoNormCondition","8.DistMovedbyCs"
            write(6,"(A12, 7A24)") "# Iteration","2.ERPotEnergy","3.HijSqrdPotEnergy","4.PotEnergy","5.Force",              &
                "6.Totalcorrforce","7.OrthoNormCondition","8.DistMovedbyCs"
        elseif (tERLocalization) then
            write(transform_unit,"(A12, 5A24)") "# Iteration","2.Sum_i<ii|ii>","3.Force","4.TotCorrForce",   &
                "5.OrthoNormCondition","6.DistMovedbyCs"
            write(6,"(A12, 5A24)") "Iteration","2.Sum_i<ii|ii>","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        else
            write(transform_unit,"(A12, 5A24)") "# Iteration","2.PotEnergy","3.Force","4.Totalcorrforce",    &
                "5.OrthoNormCondition","6.DistMovedbyCs"
            write(6,"(A12, 5A24)") "Iteration","2.PotEnergy","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        end if

    end subroutine InitLocalOrbs

    subroutine InitRotCalc()

        ! Sets up the initial arrays to be used in the orbital rotation.    

        character(len=*), parameter :: this_routine = 'InitRotCalc'
        real(dp) :: RAN2
        integer :: i, j, Const, iseed = -8, MinRot, MaxRot

        call InitSymmArrays()
! Creates an indexing system for each of the cases with symmetry on/off, and mixing all orbitals or separating
! the occuppied from virtual.
! The arrays used in this routine are labelled with a 2 (SymLabelList2_rot and SymLabelCount2), so as to not
! mess up the spawing/FCI calcs.
        do i = 1, NoOrbs
            SymLabelList3_rot(i) = SymLabelList2_rot(i)
        end do

        ! Set up constraint labels.  Constraint l is the dot product of i.j.
        Const = 0
        do i = 1, NoOrbs
            do j = i, NoOrbs
                Const = Const+1
                Lab(1, Const) = i
                Lab(2, Const) = j
            end do
        end do

        ! Just a check that the number of constraints labeled is the same as
        ! that calculated above.
        write(6,*) 'Total number of constraints  =  ', TotNoConstraints
        if (Const /= TotNoConstraints) then
            call Stop_all(this_routine,'ERROR in the number of constraints calculated.  lmax does not equal TotNoConstraints')
        end if
 
! Zero/initialise the arrays
! In the case where symmetry is kept, the starting transformation matrix is just the identity.  Starting with a symmetric system
! means the symmetry is never broken.
! When we are breaking the symmetry, the starting transformation matrix is completely random, (and then orthonormalised).
! The ordering of the orbitals in CoeffT1 follow the ordering in SymLabelList2_rot.

        CoeffT1(:,:) = 0.0_dp
        if (tRotateOccOnly) then
            MinRot = 1
            MaxRot = NoOcc
        elseif (tRotateVirtOnly) then
            MinRot = NoOcc+1
            MaxRot = NoOrbs
        else
            MinRot = 1
            MaxRot = NoOrbs
        end if
        do i = 1, NoOrbs
            CoeffT1(i,i) = 1.0_dp
        end do
        ! If the symmetry is kept on, start with the symmetric identity matrix of coefficients, and it will be maintained.

! When only the occupied or virtual orbitals are rotated, the non rotated orbitals are still included in the transformation matrix,
! but start as the identity, and remain that way throughout.
        if (lNoSymmetry) then
            do i = MinRot, MaxRot
                do j = MinRot, MaxRot
                    CoeffT1(j,i) = RAN2(iseed)*(1E-02_dp)
                end do
            end do
        end if

        ! Ensures transformation matrix elements between the occupied and
        ! virtual orbitals are 0 (should be the case anyway though).
        if (tSeparateOccVirt) call ZeroOccVirtElements(CoeffT1)

        ! Orthonormalise starting matrix.        
        call GRAMSCHMIDT(CoeffT1, NoOrbs)

! A UMATTemp is created (from the UMAT read in from the FCIDUMP) using the rotate orbs indexing system.
! i.e. in UMatTemp(1,1,1,1), the orbital involved is the first in SymLabelList2_rot.
! Doing this now, rather than using UMatInd in each transform2elint routine proved a lot faster.
        DerivCoeff(:,:) = 0.0_dp
        UMATTemp01(:,:,:,:) = 0.0_dp
        if (((.not.tERLocalization) .and. (.not.tReadInCoeff) .and. (.not.tUseMP2VarDenMat) &
             .and. (.not.tFindCINatOrbs) .and. (.not.tUseHFOrbs))&
            &.or.(tERLocalization .and. tStoreSpinOrbs)) UMATTemp02(:,:,:,:) = 0.0_dp

        call CopyAcrossUMAT()

        call TestOrthonormality()
    
        ! With UMAT with the correct indexing and the starting coefficient, 
        ! find the partially transformed four index integrals (and hence the
        ! initial potential energy), and then the initial force.

        if (tERLocalization .and. (.not.tStoreSpinOrbs)) then
            call Transform2ElIntsERlocal()
        else
            call Transform2ElInts()
        end if
        call FindTheForce()

        if (tPrintInts) then
            ! This sets the initial values for the integral sums being printed.
            ! Values printed are then relative to these initial sums, per
            ! integral.
            tInitIntValues  =  .true.
            call PrintIntegrals()
            tInitIntValues  =  .false.
        end if

    end subroutine InitRotCalc

    subroutine CopyAcrossUMAT()

        integer :: a, b, g, d, i, j, k, l
        real(dp) :: s, t

        if (((.not.tERLocalization) .and. (.not.tReadInCoeff) .and. (.not.tUseMP2VarDenMat) .and. (.not.tFindCINatOrbs))&
        &.or.(tERLocalization .and. tStoreSpinOrbs)) TMAT2DTemp(:,:)  =  0.0_dp

        ! These loops can be sped up with spatial symmetry and pairwise
        ! permutation symmetry if needed.
        do a = 1, NoOrbs
            i = SymLabelList2_rot(a) ! The spin orbital we are looking for.
            do g = 1, a
                j = SymLabelList2_rot(g)
                if (((.not.tERLocalization) .and. (.not.tReadInCoeff) .and. (.not.tUseMP2VarDenMat) .and. (.not.tFindCINatOrbs))&
                &.or.(tERLocalization .and. tStoreSpinOrbs)) then
                    if (tStoreSpinOrbs) then
                        s = real(TMAT2D(i,j), dp)
                        TMAT2DTemp(a,g) = s
                        TMAT2DTemp(g,a) = s
                    else
                        s = real(TMAT2D(2*i, 2*j), dp)
                        TMAT2DTemp(a,g) = s
                        TMAT2DTemp(g,a) = s
                    end if
                end if

                do b = 1, NoOrbs
                    k = SymLabelList2_rot(b)
                    do d = 1, b
                        l = SymLabelList2_rot(d)
                        t = real(UMAT(UMatInd(i, k, j, l, 0, 0)), dp)
                        UMATTemp01(a,g,b,d) = t    ! a, g, d, b chosen to make 'transform2elint' steps more efficient.
                        UMATTemp01(g,a,b,d) = t
                        UMATTemp01(a,g,d,b) = t
                        UMATTemp01(g,a,d,b) = t
                        if (((.not.tERLocalization) .and. (.not.tReadInCoeff) .and. &
                            (.not.tUseMP2VarDenMat) .and. (.not.tFindCINatOrbs))&
                        &.or.(tERLocalization .and. tStoreSpinOrbs)) then
                            UMATTemp02(d,b,a,g) = t    ! d, b, a, g order also chosen to speed up the transformation.
                            UMATTemp02(d,b,g,a) = t
                            UMATTemp02(b,d,a,g) = t
                            UMATTemp02(b,d,g,a) = t
                        end if
                    end do
                end do
            end do
        end do

    end subroutine CopyAcrossUMAT

    subroutine WriteStats()

        if (tLagrange) then
            write(6,"(I12, 11F18.10)") Iteration, PotEnergy, PEInts, PEOrtho, Force, ForceInts, OrthoForce, TwoEInts, &
                OrthoNorm, DistCs, DistLs, LambdaMag
            write(transform_unit,"(I12, 11F18.10)") Iteration, PotEnergy, PEInts, PEOrtho, Force, ForceInts, OrthoForce, &
                TwoEInts, OrthoNorm, DistCs, DistLs, LambdaMag
        elseif (tERLocalization .and. tHijSqrdMin) then
            if (Mod(Iteration, 10) == 0) then
                write(6,"(I12, 7F24.10)") Iteration, ERPotEnergy, HijSqrdPotEnergy, PotEnergy, Force, TotCorrectedForce, &
                                           OrthoNorm, DistCs
                write(transform_unit,"(I12, 7F24.10)") Iteration, ERPotEnergy, HijSqrdPotEnergy, PotEnergy, Force, &
                                                        TotCorrectedForce, OrthoNorm, DistCs
            end if
        else
            if (Mod(Iteration, 10) == 0) then
                write(6,"(I12, 5F24.10)") Iteration, PotEnergy, Force, TotCorrectedForce, OrthoNorm, DistCs
                write(transform_unit,"(I12, 5F24.10)") Iteration, PotEnergy, Force, TotCorrectedForce, OrthoNorm, DistCs
            end if
        end if
        call neci_flush(6)
        call neci_flush(transform_unit)

        ! After writing out stats, test for SOFTEXIT.
        if (test_SOFTEXIT()) then
            write(6,*) 'SOFTEXIT detected, finalizing new orbitals.'
            tNotConverged = .false.
        end if

    end subroutine WriteStats

    subroutine InitSymmArrays()

! This routine creates indexing arrays for the cases with symmetry on/off, and either mixing all orbitals or 
! separating the occupied and virtuals.
! The arrays used specific to the orbital rotation are named with a 2. 

! The arrays produced are as follows...
! SymLabelList2_rot(NoOrbs) contains the spatial orbitals, ordered in groups of increasing symmetry label.
! - when the orbitals are being separated, the first NoOcc of SymLabelList2_rot are the occupied, and the rest are virtual.
! - essentially this array relates the orbital labelling used in the orbital rotation (1, 2, 3 according to the order 
! - in SymLabelList2_rot) to the labels used in arrays being fed in/out of this routine (UMAT etc).

! SymLabelCounts2_rot(1:Sym) is the index in SymLabelList where the symmetry block S starts
! SymLabelCounts2_rot(2:Sym) is the number of orbitals in symmetry block S.
! E.g. if symmetry S starts at index 2 and has 3 orbitals.
! SymLabelList2_rot(2)->SymLabelList2_rot(4) will give the indexes of these orbitals.

        use sym_mod, only: GenSymStatePairs
        integer :: j, i, ierr
        character(len=*), parameter :: this_routine = 'InitSymmArrays'

        if (.not.tSeparateOccVirt) then
            SymLabelCounts(:,:) = 0
            SymLabelList(:) = 0
            if (tStoreSpinOrbs) call Stop_All(this_routine, &
                "There may be a problem with GENSymStatePairs when using spin orbitals.")
            call GENSymStatePairs(SpatOrbs,.false.)
        end if
! Sets up the SymLabelList and SymLabelCounts arrays used in the spawing etc. (When the rotate
! orbs routine is called, this has not been done yet).
! If the symmetry is on, and all orbitals are being mixed, this will end up being the same as SymLabelList2_rot.


        if (tSeparateOccVirt) then
            MinOccVirt = 1
            MaxOccVirt = 2
            if (tRotateOccOnly) then
                MaxOccVirt = 1
            elseif (tRotateVirtOnly) then
                MinOccVirt = 2
            end if
            call InitOrbitalSeparation()
            ! rewrite all the symmetry lists to account for the separation and have simple option if
            ! symmetry is off.
        else
            MinOccVirt = 1
            MaxOccVirt = 1
            allocate(SymLabelCounts2_rot(2,8), stat = ierr)
            call LogMemAlloc('SymLabelCounts2_rot', 2*8, 4, this_routine, SymLabelCounts2_rotTag, ierr)
            SymLabelCounts2_rot(:,:) = 0                     
            do i = 1, SpatOrbs   
                if (tStoreSpinOrbs) then
                    SymLabelList2_rot(2*i) = 2*SymLabelList(i)
                    SymLabelList2_rot(2*i-1) = (2*SymLabelList(i))-1
                else
                    SymLabelList2_rot(i) = SymLabelList(i)
                end if
            end do
            if (lNoSymmetry) then
                SymLabelCounts2_rot(1,1) = 1
                SymLabelCounts2_rot(2,1) = NoOrbs
            else
                do j = 1, 8
                    if (tStoreSpinOrbs) then
                        SymLabelCounts2_rot(1,j) = (2*SymLabelCounts(1,j))-1
                        SymLabelCounts2_rot(2,j) = 2*SymLabelCounts(2,j)
                    else
                        do i = 1, 2
                            SymLabelCounts2_rot(i,j) = SymLabelCounts(i,j)
                        end do
                    end if
                end do
            end if
        end if

        do i = 1, NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i)) = i
        end do
        
    end subroutine InitSymmArrays

    subroutine EquateDiagFock()

        integer :: irr, NumInSym, Orbi, Orbj, w, i, j, k, ConjInd, OrbjConj
        real(dp) :: Angle, AngleConj, Check, Norm

        CoeffT1(:,:) = 0.0_dp

        ! Do virtual and occupied orbitals seperately.
        do w = MinOccVirt, MaxOccVirt

            ! Loop over irreps.
            do irr = 1, 8

                NumInSym = SymLabelCounts2_rot(2,(w-1)*8+irr)

                ! Loop over the j-orthogonal vectors to create in this
                ! symmetry block.
                do j = 1, NumInSym

                    Orbj = SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+j)

                    ! See if this vector has already been done.
                    Check = 0.0_dp
                    do i = 1, NoOrbs
                        Check = Check+CoeffT1(i, Orbj)
                    end do
                    if (Check /= 0.0_dp) then
                        ! This vector is a conjugate pair of another vector and
                        ! has already been worked out...
                        cycle
                    end if

                    ! Find out if we this vector will be complex. It will be
                    ! real if j = N or j = N/2
                    if (j == NumInSym) then
                    !i The vector will be the normalized 1, 1, 1 vector.

                        do i = 1, NumInSym
                            Orbi = SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)
                            CoeffT1(Orbi, Orbj) = 1/SQRT(real(NumInSym, dp))
                        end do

                    elseif ((mod(NumInSym, 2) == 0) .and. (j == (NumInSym/2))) then

                        do i = 1, NumInSym
                            Orbi = SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)
                            if (mod(i,2) == 1) then
                                CoeffT1(Orbi, Orbj) = -1/SQRT(real(NumInSym, dp))
                            else
                                CoeffT1(Orbi, Orbj) = 1/SQRT(real(NumInSym, dp))
                            end if
                        end do

                    else
                        ! Vector is complex - find its conjugate vector - do
                        ! these at the same time.
                        ConjInd = NumInSym-j
                        OrbjConj = SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+ConjInd)
                        
                        do i = 1, NumInSym
                            
                            Orbi = SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)

                            Angle = real(i*j*2, dp)*PI/real(NumInSym, dp)
                            AngleConj = real(i*ConjInd*2, dp)*PI/real(NumInSym, dp)

                            CoeffT1(Orbi, Orbj) = (1/SQRT(real(2*NumInSym, dp)))*(COS(Angle)+COS(AngleConj))
                            CoeffT1(Orbi, OrbjConj) = (1/SQRT(real(2*NumInSym, dp)))*(SIN(Angle)-SIN(AngleConj))

                        end do

                    end if

                end do

            end do
        end do

        do j = 1, NoOrbs
            Norm = 0.0_dp
            do i = 1, NoOrbs
                Norm = Norm+(CoeffT1(i,j)**2)
            end do
            if (Norm == 0.0_dp) then
                CoeffT1(j,j) = 1.0_dp
            end if
        end do

        do j = 1, NoOrbs
            do i = 1, NoOrbs
                write(6,"(G13.5)", advance='no') CoeffT1(j,i)
            end do
            write(6,*) ""
        end do

        !Check normalization.
        do j = 1, NoOrbs
            Norm = 0.0_dp
            do i = 1, NoOrbs
                Norm = Norm+(CoeffT1(i,j)**2)
            end do
            if (abs(Norm-1.0_dp) > 1.0e-7_dp) then
                call Stop_All("EquateDiagFock","Rotation Coefficients not normalized")
            end if
        end do

        ! Check orthogonality.
        do j = 1, NoOrbs
            do i = 1, NoOrbs
                if (i == j) cycle
                Norm = 0.0_dp
                do k = 1, NoOrbs
                    Norm = Norm+(CoeffT1(k,j)*CoeffT1(k,i))
                end do
                if (abs(Norm) > 1.0e-7_dp) then
                    write(6,*) "COLUMNS: ", j, i
                    call Stop_All("EquateDiagFock","RotationCoefficients not orthogonal")
                end if
            end do
        end do

    end subroutine EquateDiagFock

    subroutine InitOrbitalSeparation()

! This subroutine is called if the SEPARATEOCCVIRT keyword is present in the input, it sets up SymLabelList2_rot so that the first 
! NoOcc orbitals are the HF occupied, and the rest the virtual.  Within this separation, orbitals are ordered in symmetry 
! groups. 
! This means that two iterations of the rotate orbs routine will be performed, the first treats the occupied orbitals and the second
! the virtual.

        integer :: i, j, ierr, SymCurr, Symi
        integer(TagIntType) :: SymVirtOrbsTag, SymOccOrbsTag
        integer :: lo, hi
        integer, allocatable :: SymVirtOrbs(:), SymOccOrbs(:)
        character(len=*), parameter :: this_routine = 'InitOrbitalSeparation'

        allocate(SymLabelCounts2_rot(2, 16), stat = ierr)
        call LogMemAlloc('SymLabelCounts2_rot', 2*16, 4, this_routine, SymLabelCounts2_rotTag, ierr)
        SymLabelCounts2_rot(:,:) = 0
        ! first 8 refer to the occupied, and the second to the virtual.

        allocate(LabVirtOrbs(NoOrbs-NoOcc), stat = ierr)
        call LogMemAlloc('LabVirtOrbs',(NoOrbs-NoOcc), 4, this_routine, LabVirtOrbsTag, ierr)
        LabVirtOrbs(:) = 0
        allocate(LabOccOrbs(NoOcc), stat = ierr)
        call LogMemAlloc('LabOccOrbs',(NoOcc), 4, this_routine, LabOccOrbsTag, ierr)
        LabOccOrbs(:) = 0
        allocate(SymVirtOrbs(NoOrbs-NoOcc), stat = ierr)
        call LogMemAlloc('SymVirtOrbs',(NoOrbs-NoOcc), 4, this_routine, SymVirtOrbsTag, ierr)
        SymVirtOrbs(:) = 0
        allocate(SymOccOrbs(NoOcc), stat = ierr)
        call LogMemAlloc('SymOccOrbs',(NoOcc), 4, this_routine, SymOccOrbsTag, ierr)
        SymOccOrbs(:) = 0


! First fill SymLabelList2_rot.

! This picks out the NoOcc lowest energy orbitals from BRR as these will be the occupied.
! these are then ordered according to symmetry, and the same done to the virtual.
        do i = 1, NoOcc
            if (tStoreSpinOrbs) then
                LabOccOrbs(i) = BRR(i)
                SymOccOrbs(i) = int(G1(LabOccOrbs(i))%sym%S)
            else
                LabOccOrbs(i) = (BRR(2*i))/2
                SymOccOrbs(i) = int(G1(LabOccOrbs(i)*2)%sym%S)
            end if
        end do
        
        call sort (SymOccOrbs, LabOccOrbs)
        ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in terms of symmetry). 

        do i = 1, NoOrbs-NoOcc
            if (tStoreSpinOrbs) then
                LabVirtOrbs(i) = BRR(i+NEl)
                SymVirtOrbs(i) = int(G1(LabVirtOrbs(i))%sym%S)
            else
                LabVirtOrbs(i) = (BRR((2*i)+NEl))/2
                SymVirtOrbs(i) = int(G1(LabVirtOrbs(i)*2)%sym%S)
            end if
        end do
        
        call sort (SymVirtOrbs, LabVirtOrbs)

! SymLabelList2_rot is then filled with the symmetry ordered occupied then virtual arrays.        
        do i = 1, NoOcc
            SymLabelList2_rot(i) = LabOccOrbs(i)
        end do
        j = 0
        do i = NoOcc+1, NoOrbs
            j = j+1
            SymLabelList2_rot(i) = LabVirtOrbs(j)
        end do

!************
! Second fill SymLabelCounts2_rot.
! - the first 8 places of SymLabelCounts2_rot(1,:) and SymLabelCounts2_rot(2,:) refer to the occupied orbitals 
! - and the second 8 to the virtuals.

        if (lNoSymmetry) then
            ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
            SymLabelCounts2_rot(1,1) = 1
            SymLabelCounts2_rot(1,9) = NoOcc+1
            SymLabelCounts2_rot(2,1) = NoOcc
            SymLabelCounts2_rot(2,9) = NoOrbs-NoOcc
        else 
            ! otherwise we run through the occupied orbitals, counting the number with each symmetry
            ! and noting where in SymLabelList2_rot each symmetry block starts.
            SymCurr = 0
            SymLabelCounts2_rot(1,1) = 1
            do i = 1, NoOcc
                if (tStoreSpinOrbs) then
                    Symi = int(G1(SymLabelList2_rot(i))%sym%S)
                else
                    Symi = int(G1(SymLabelList2_rot(i)*2)%sym%S)
                end if
                SymLabelCounts2_rot(2,(Symi+1)) = SymLabelCounts2_rot(2,(Symi+1))+1
                if (Symi > SymCurr) then
                    SymLabelCounts2_rot(1,(Symi+1)) = i
                    SymCurr = Symi
                end if
            end do
            ! the same is then done for the virtuals.
            SymCurr = 0
            SymLabelCounts2_rot(1,9) = NoOcc+1
            do i = NoOcc+1, NoOrbs
                if (tStoreSpinOrbs) then
                    Symi = int(G1(SymLabelList2_rot(i))%sym%S)
                else
                    Symi = int(G1(SymLabelList2_rot(i)*2)%sym%S)
                end if
                SymLabelCounts2_rot(2,(Symi+9)) = SymLabelCounts2_rot(2,(Symi+9))+1
                if (Symi > SymCurr) then
                    SymLabelCounts2_rot(1,(Symi+9)) = i
                    SymCurr = Symi
                end if
            end do
        end if

        ! Go through each symmetry group, making sure the orbital pairs are ordered lowest to highest.
        do i = 1, 16
            if (SymLabelCounts2_rot(2,i) /= 0) then
                lo  =  SymLabelCounts2_rot(1,i)
                hi  =  lo + SymLabelCounts2_rot(2,i) - 1
                call sort (SymLabelList2_rot(lo:hi))
            end if
        end do

! Deallocate the arrays just used in this routine.
        deallocate(LabOccOrbs)
        call LogMemDealloc(this_routine, LabOccOrbsTag)
        deallocate(LabVirtOrbs)
        call LogMemDealloc(this_routine, LabVirtOrbsTag)
        deallocate(SymOccOrbs)
        call LogMemDealloc(this_routine, SymOccOrbsTag)
        deallocate(SymVirtOrbs)
        call LogMemDealloc(this_routine, SymVirtOrbsTag)

    end subroutine InitOrbitalSeparation

    subroutine Diagonalizehij()

! This routine takes the original <i|h|j> matrix and diagonalises it.  The resulting coefficients from this process 
! are then the rotation coefficients to be applied to the four index integrals etc.
! This eliminates the <i|h|j> elements from the single excitations, and leaves only coulomb and exchange terms.
! In order to maintain the same HF energy, only the virtual elements are diagonalised, within symmetry blocks.

        integer :: i, j, Sym, ierr, NoSymBlock, WorkSize, WorkCheck, SymStartInd
        integer(TagIntType) WorkTag, DiagTMAT2DBlockTag, TMAT2DSymBlockTag
        real(dp), allocatable :: TMAT2DSymBlock(:,:), DiagTMAT2DBlock(:), Work(:)
        character(len=*), parameter :: this_routine = 'Diagonalizehij'

        write(6,*) 'The original coefficient matrix'
        do i = 1, NoOrbs
            do j = 1, NoOrbs
                write(6,'(F20.10)', advance='no') CoeffT1(j,i)
            end do
            write(6,*) ''
        end do

        write(6,*) 'The original TMAT2D matrix'
        do i = 1, NoOrbs
            do j = 1, NoOrbs
                write(6,'(F20.10)', advance='no') TMAT2DTemp(j,i)
            end do
            write(6,*) ''
        end do
        TMAT2DRot(:,:) = 0.0_dp
        DiagTMAT2Dfull(:) = 0.0_dp

! Now need to pick out symmetry blocks, from the virtual orbitals and diagonalize them.

! Take first symmetry, (0) and find the number of virtual orbitals with this symmetry.  If this is greater than 1, 
! take the block, diagonlize it, and put it into TMAT2DRot.
        
        Sym = 0
        WorkSize = -1
        do while (Sym <= 7)

            NoSymBlock = SymLabelCounts2_rot(2, Sym+9)

            SymStartInd = SymLabelCounts2_rot(1, Sym+9)-1
            ! This is one less than the index that the symmetry starts, so that when we run through i = 1,..., we can
            ! start at SymStartInd+i.

            if (NoSymBlock > 1) then
                allocate(TMAT2DSymBlock(NoSymBlock, NoSymBlock), stat = ierr)
                call LogMemAlloc('TMAT2DSymBlock', NoSymBlock**2, 8, this_routine, TMAT2DSymBlockTag, ierr)
                allocate(DiagTMAT2DBlock(NoSymBlock), stat = ierr)
                call LogMemAlloc('DiagTMAT2DBlock', NoSymBlock, 8, this_routine, DiagTMAT2DBlockTag, ierr)

                WorkCheck = 3*NoSymBlock+1
                WorkSize = WorkCheck
                allocate(Work(WorkSize), stat = ierr)
                call LogMemAlloc('Work', WorkSize, 8, this_routine, WorkTag, ierr)

                do j = 1, NoSymBlock
                    do i = 1, NoSymBlock
                        TMAT2DSymBlock(i,j) = TMAT2DTemp(SymStartInd+i, SymStartInd+j)
                    end do
                end do

                write(6,*) '*****'
                write(6,*) 'Symmetry ', Sym,' has ', NoSymBlock,' orbitals .'
                write(6,*) 'The TMAT2D for this symmetry block is '
                do i = 1, NoSymBlock
                    do j = 1, NoSymBlock
                        write(6,'(F20.10)', advance='no') TMAT2DSymBlock(j,i)
                    end do
                    write(6,*) ''
                end do

                call DSYEV('V','U', NoSymBlock, TMAT2DSymBlock, NoSymBlock, DiagTMAT2Dblock, Work, WorkSize, ierr)
                ! TMAT2DSymBlock goes in as the original TMAT2DSymBlock, comes out as the eigenvectors (Coefficients).
                ! TMAT2DBlock comes out as the eigenvalues in ascending order.
                if (ierr /= 0) then
                    write(6,*) 'Problem with symmetry, ', Sym,' of TMAT2D'
                    call neci_flush(6)
                    call Stop_All(this_routine,"Diagonalization of TMAT2DSymBlock failed...")
                end if

                write(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are,'
                do i = 1, NoSymBlock
                    write(6,'(F20.10)', advance='no') DiagTMAT2Dblock(i)
                end do
                write(6,*) ''
                write(6,*) 'These go from orbital,', SymStartInd+1,' to ', SymStartInd+NoSymBlock
               
                do i = 1, NoSymBlock
                    DiagTMAT2Dfull(SymStartInd+i-NoOcc) = DiagTMAT2DBlock(i)
                end do

                ! CAREFUL if eigenvalues are put in ascending order, this may not be correct, with the labelling system.
                ! may be better to just take coefficients and transform TMAT2DRot in transform2elints.
                ! a check that comes out as diagonal is a check of this routine anyway.

                write(6,*) 'The eigenvectors (coefficients) for symmtry block ', Sym
                do i = 1, NoSymBlock
                    do j = 1, NoSymBlock
                        write(6,'(F20.10)', advance='no') TMAT2DSymBlock(j,i)
                    end do
                    write(6,*) ''
                end do

             
                ! Directly fill the coefficient matrix with the eigenvectors from the diagonalization.
                do j = 1, NoSymBlock
                    do i = 1, NoSymBlock
                        CoeffT1(SymStartInd+i, SymStartInd+j) = TMAT2DSymBlock(i,j)
                    end do
                end do

                deallocate(Work)
                call LogMemDealloc(this_routine, WorkTag)

                deallocate(DiagTMAT2DBlock)
                call LogMemDealloc(this_routine, DiagTMAT2DBlockTag)

                deallocate(TMAT2DSymBlock)
                call LogMemDealloc(this_routine, TMAT2DSymBlockTag)
            elseif (NoSymBlock == 1) then
                DiagTMAT2Dfull(SymStartInd+1-NoOcc) = TMAT2DTemp(SymStartInd+1, SymStartInd+1)
                write(6,*) '*****'
                write(6,*) 'Symmetry ', Sym,' has only one orbital.'
                write(6,*) 'Copying diagonal element,', SymStartInd+1,'to DiagTMAT2Dfull'
            end if

            Sym = Sym+1
        end do
 
        write(6,*) '*****'
        write(6,*) 'The final coefficient matrix'
        do i = 1, NoOrbs
            do j = 1, NoOrbs
                write(6,'(F20.10)', advance='no') CoeffT1(j,i)
            end do
            write(6,*) ''
        end do

        write(6,*) '*****'
        write(6,*) 'The diagonal elements of TMAT2D'
        do i = 1,(NoOrbs-NoOcc)
            write(6,*) DiagTMAT2Dfull(i)
        end do

    end subroutine Diagonalizehij

    subroutine ZeroOccVirtElements(Coeff)

! This routine sets all the elements of the coefficient matrix that connect occupied and virtual orbitals to 0.
! This ensures that only occupied mix with occupied and virtual mix with virtual.

        real(dp) :: Coeff(NoOrbs, NoOrbs)
        integer :: i, j

        do i = 1, NoOcc
            do j = NoOcc+1, NoOrbs
                Coeff(i,j) = 0.0_dp
                Coeff(j,i) = 0.0_dp
            end do
        end do

    end subroutine ZeroOccVirtElements

    subroutine FindNewOrbs()
           
        if (tERLocalization .and. (.not.tStoreSpinOrbs)) then
            call Transform2ElIntsERlocal()
        else
! Find the partially (and completely) transformed 4 index integrals to be used in further calcs.
            call Transform2ElInts()     
        end if


!Find derivatives of the c and lambda matrices and print the sum of off-diagonal matrix elements.
        call FindTheForce()
        ! This finds the unconstrained force (unless the lagrange keyword is present).
      
!Update coefficents by moving them in direction of force. Print sum of squared changes in coefficients. 
        if (tShake) then
            call ShakeConstraints()
            ! Find the force that moves the coefficients while keeping them orthonormal, and use it 
            ! to get these new coefficients.
        else
            call UseTheForce()
            ! This can be either completely unconstrained, or have the lagrange constraints imposed.
        end if
!The coefficients coefft1(a,m) are now those that have been shifted by the time step.

!Test these for orthonomaility and then convergence.
!If they do not meet the convergence criteria, they go back into the previous step to produce another set of coefficients.

        call set_timer(testorthoconver_time, 30)        

        call TestOrthonormality()
!Force should go to zero as we end in minimum - test for this

        call TestForConvergence()

        call halt_timer(testorthoconver_time)

    end subroutine FindNewOrbs
    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.

    subroutine Transform2ElInts()

        integer :: i, j, k, l, a, b, g, d
        real(dp) :: t, Temp4indints(NoRotOrbs, NoOrbs)
        real(dp) :: Temp4indints02(NoRotOrbs, NoRotOrbs)  

        call set_timer(Transform2ElInts_time, 30)

!Zero arrays from previous transform

        TwoIndInts01(:,:,:,:) = 0.0_dp
        FourIndInts(:,:,:,:) = 0.0_dp

        if (tNotConverged) then
            TwoIndInts02(:,:,:,:) = 0.0_dp
            ThreeIndInts01(:,:,:,:) = 0.0_dp
            ThreeIndInts02(:,:,:,:) = 0.0_dp
            ThreeIndInts03(:,:,:,:) = 0.0_dp
            ThreeIndInts04(:,:,:,:) = 0.0_dp
            FourIndInts02(:,:,:,:) = 0.0_dp
        end if

! ************
!Transform the 1 electron, 2 index integrals (<i|h|j>).
        if (tNotConverged) then
            TMAT2DRot(:,:) = 0.0_dp
            TMAT2DPartRot01(:,:) = 0.0_dp
            TMAT2DPartRot02(:,:) = 0.0_dp

            call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs,    &
                TMAT2DTemp(:,:), NoOrbs, 0.0_dp, TMAT2DPartRot01(:,:), NoOrbs)

            call dgemm('T','T', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs,    &
                TMAT2DTemp(:,:), NoOrbs, 0.0_dp, TMAT2DPartRot02(:,:), NoOrbs)
     
            call dgemm('T','T', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs,    &
                TMAT2DPartRot01(:,:), NoOrbs, 0.0_dp, TMAT2DRot(:,:), NoOrbs)

        end if

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i, j, k, l, 0, 0)

        do b = 1, NoOrbs
            do d = 1, b
                Temp4indints(:,:) = 0.0_dp
                call dgemm('T','N', NoRotOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, UMatTemp01(:,:, d, b), NoOrbs, &
                    0.0_dp, Temp4indints(:,:), NoRotOrbs)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, we need the transpose of this.

                Temp4indints02(:,:) = 0.0_dp
                call dgemm('T','T', NoRotOrbs, NoRotOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, Temp4indints(:,:), NoRotOrbs, &
                    0.0_dp, Temp4indints02(:,:), NoRotOrbs)
                ! Get Temp4indits02(i,k)

                do i = 1, NoRotOrbs
                    do k = 1, i
                        TwoIndInts01(d,b,k,i) = Temp4indints02(k,i)
                        TwoIndInts01(b,d,k,i) = Temp4indints02(k,i)
                        TwoIndInts01(d,b,i,k) = Temp4indints02(k,i)
                        TwoIndInts01(b,d,i,k) = Temp4indints02(k,i)

                    end do
                end do
            end do
        end do

! These calculations are unnecessary when this routine is calculated to finalize the new orbs.
        if (tNotConverged) then
            do g = 1, NoOrbs                
                do a = 1, g
                    Temp4indints(:,:) = 0.0_dp
                    call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, UMatTemp02(:,:, a, g), NoOrbs, &
                        0.0_dp, Temp4indints(:,:), NoOrbs)

                    Temp4indints02(:,:) = 0.0_dp
                    call dgemm('T','T', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, Temp4indints(:,:), NoOrbs, &
                        0.0_dp, Temp4indints02(:,:), NoOrbs)

                    do l = 1, NoOrbs
                        do j = 1, l
                            TwoIndInts02(g,a,j,l) = Temp4indints02(j,l)
                            TwoIndInts02(a,g,j,l) = Temp4indints02(j,l)
                            TwoIndInts02(g,a,l,j) = Temp4indints02(j,l)
                            TwoIndInts02(a,g,l,j) = Temp4indints02(j,l)
                        end do
                    end do
                end do
            end do
        end if

! Calculating the 3 transformed, 4 index integrals. 01 = a untransformed, 02 = b, 03 = g, 04 = d

        do i = 1, NoRotOrbs
            do k = 1, i
                Temp4indints(:,:) = 0.0_dp
                call dgemm('T','N', NoRotOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, TwoIndInts01(:,:, k, i), NoOrbs, &
                    0.0_dp, Temp4indints(:,:), NoRotOrbs)

                if (tNotConverged) then
                    do b = 1, NoOrbs
                        do l = 1, NoOrbs
                            ThreeIndInts02(i,k,l,b) = Temp4indints(l,b)
                            ThreeIndInts02(k,i,l,b) = Temp4indints(l,b)
                        end do
                    end do
                    Temp4indints02(:,:) = 0.0_dp
                    call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, TwoIndInts01(:,:, k, i), NoOrbs, &
                        0.0_dp, Temp4indints02(:,:), NoRotOrbs)
                    do d = 1, NoOrbs
                        do j = 1, NoOrbs
                            ThreeIndInts04(k,i,j,d) = Temp4indints02(j,d)
                            ThreeIndInts04(i,k,j,d) = Temp4indints02(j,d)
                        end do
                    end do
                end if
                Temp4indints02(:,:) = 0.0_dp
                call dgemm('T','T', NoRotOrbs, NoRotOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, Temp4indints(:,:), NoRotOrbs, &
                    0.0_dp, Temp4indints02(:,:), NoRotOrbs)
                do l = 1, NoRotOrbs
                    do j = 1, l
                        FourIndInts(i,j,k,l) = Temp4indints02(j,l)
                        FourIndInts(i,l,k,j) = Temp4indints02(j,l)
                        FourIndInts(k,j,i,l) = Temp4indints02(j,l)
                        FourIndInts(k,l,i,j) = Temp4indints02(j,l)

                        if (tNotConverged) then
                            FourIndInts02(j,k,l,i) = Temp4indints02(j,l)
                            FourIndInts02(j,i,l,k) = Temp4indints02(j,l)
                            FourIndInts02(l,k,j,i) = Temp4indints02(j,l)
                            FourIndInts02(l,i,j,k) = Temp4indints02(j,l)
                        end if
                    end do
                end do
            end do
        end do

        if (tNotConverged) then
            do l = 1, NoOrbs
                do j = 1, l
                    Temp4indints(:,:) = 0.0_dp
                    call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, TwoIndInts02(:,:, j, l), NoOrbs, &
                        0.0_dp, Temp4indints(:,:), NoOrbs)
                    do a = 1, NoOrbs
                        do k = 1, NoOrbs
                            ThreeIndInts01(k,j,l,a) = Temp4indints(k,a)
                            ThreeIndInts01(k,l,j,a) = Temp4indints(k,a)
                        end do
                    end do
                    Temp4indints(:,:) = 0.0_dp
                    call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, TwoIndInts02(:,:, j, l), NoOrbs, &
                        0.0_dp, Temp4indints(:,:), NoOrbs)
                    do g = 1, NoOrbs
                        do i = 1, NoOrbs
                            ThreeIndInts03(i,l,j,g) = Temp4indints(i,g)
                            ThreeIndInts03(i,j,l,g) = Temp4indints(i,g)
                        end do
                    end do
                end do
            end do
        end if

! ***************************
! Calc the potential energies for this iteration (with these transformed integrals).        

! This can be sped up by merging the calculations of the potentials with the transformations, but while 
! we are playing around with different potentials, it is simpler to keep these separate.
    
        if ((.not.tReadInCoeff) .and. (.not.tUseMP2VarDenMat) .and. (.not.tFindCINatOrbs) .and. (.not.tUseHFOrbs)) then
            
            PotEnergy = 0.0_dp
            TwoEInts = 0.0_dp
            PEInts = 0.0_dp
            call CalcPotentials()

            if (tPrintInts) call PrintIntegrals()
            if ((Iteration == 0).or.((.not.tNotConverged) .and. (Iteration > 1))) call WriteDoubHisttofile()
            if (tROHistSingExc .and. (Iteration == 0)) call WriteSingHisttofile()

! If doing Lagrange orthormalisations, find the change of the potential energy due to the orthonormality 
! of the orbitals...
            if (tLagrange) then
                PEOrtho = 0.0_dp
                do i = 1, NoOrbs
                    do j = 1, NoOrbs
                        t = 0.0_dp
                        do a = 1, NoOrbs
                            t = CoeffT1(a,i)*CoeffT1(a,j)
                        end do
                        if (i == j) t = t-1.0_dp
                        PEOrtho = PEOrtho-Lambdas(i,j)*t
                        PotEnergy = PotEnergy-Lambdas(i,j)*t
                    end do
                end do
            end if
        end if

        call halt_timer(Transform2ElInts_Time)


    end subroutine Transform2ElInts

! This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
! This is v memory inefficient and currently does not use any spatial symmetry information.

    subroutine Transform2ElIntsMemSave()

        integer :: i, j, k, l, a, b, g, d, ierr, a2, b2, g2, d2
        integer(TagIntType) Temp4indintsTag
        real(dp), allocatable :: Temp4indints(:,:)

#ifdef __CMPLX
        call stop_all('Transform2ElIntsMemSave', 'Rotating orbitals not implemented for complex orbitals.')
#endif
        
        Transform2ElInts_Time%timer_name = 'Transform2ElIntsTime'
        call set_timer(Transform2ElInts_time, 30)

        ! Zero arrays from previous transform.
 
        allocate(Temp4indints(NoRotOrbs, NoOrbs), stat = ierr)
        call LogMemAlloc('Temp4indints', NoRotOrbs*NoOrbs, 8,'Transform2ElIntsMemSave', Temp4indintsTag, ierr)
        if (ierr /= 0) call Stop_All('Transform2ElIntsMemSave','Problem allocating memory to Temp4indints.')
 
        FourIndInts(:,:,:,:) = 0.0_dp

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i, j, k, l, 0, 0)

        do b = 1, NoOrbs
            if (tTurnStoreSpinOff) then
                b2 = CEILING(real(SymLabelList2_rot(b), dp)/2.0_dp)
            else
                b2 = SymLabelList2_rot(b)
            end if
            do d = 1, b
                if (tTurnStoreSpinOff) then
                    d2 = CEILING(real(SymLabelList2_rot(d), dp)/2.0_dp)
                else
                    d2 = SymLabelList2_rot(d)
                end if
                do a = 1, NoOrbs
                    if (tTurnStoreSpinOff) then
                        a2 = CEILING(real(SymLabelList2_rot(a), dp)/2.0_dp)
                    else
                        a2 = SymLabelList2_rot(a)
                    end if
                    do g = 1, a
                        if (tTurnStoreSpinOff) then
                            g2 = CEILING(real(SymLabelList2_rot(g), dp)/2.0_dp)
                        else
                            g2 = SymLabelList2_rot(g)
                        end if
                        FourIndInts(a,g,b,d) = real(UMAT(UMatInd(a2, b2, g2, d2, 0, 0)), dp)
                        FourIndInts(g,a,b,d) = real(UMAT(UMatInd(a2, b2, g2, d2, 0, 0)), dp)
                        FourIndInts(a,g,d,b) = real(UMAT(UMatInd(a2, b2, g2, d2, 0, 0)), dp)
                        FourIndInts(g,a,d,b) = real(UMAT(UMatInd(a2, b2, g2, d2, 0, 0)), dp)
                    end do
                end do
                Temp4indints(:,:) = 0.0_dp
                call dgemm('T','N', NoRotOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, FourIndInts(1:NoOrbs,1:NoOrbs,b,d),&
                    NoOrbs, 0.0_dp, Temp4indints(1:NoRotOrbs, 1:NoOrbs), NoRotOrbs)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, we need the transpose of this.

                call dgemm('T','T', NoRotOrbs, NoRotOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, Temp4indints(1:NoRotOrbs,1:NoOrbs),&
                    NoRotOrbs, 0.0_dp, FourIndInts(1:NoRotOrbs, 1:NoRotOrbs, b, d), NoRotOrbs)
                ! Get Temp4indits02(i,k)

                do i = 1, NoRotOrbs
                    do k = 1, i
                        FourIndInts(i,k,d,b) = FourIndInts(i,k,b,d)
                        FourIndInts(k,i,d,b) = FourIndInts(i,k,b,d)
                        FourIndInts(i,k,b,d) = FourIndInts(i,k,b,d)
                        FourIndInts(k,i,b,d) = FourIndInts(i,k,b,d)
                    end do
                end do
            end do
        end do
        
! Calculating the 3 transformed, 4 index integrals. 01 = a untransformed, 02 = b, 03 = g, 04 = d
        do i = 1, NoRotOrbs
            do k = 1, i

                Temp4indints(:,:) = 0.0_dp
                call dgemm('T','N', NoRotOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, FourIndInts(i,k,1:NoOrbs,1:NoOrbs),&
                    NoOrbs, 0.0_dp, Temp4indints(1:NoRotOrbs, 1:NoOrbs), NoRotOrbs)

                call dgemm('T','T', NoRotOrbs, NoRotOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, Temp4indints(1:NoRotOrbs,1:NoOrbs),&
                    NoRotOrbs, 0.0_dp, FourIndInts(i, k, 1:NoRotOrbs, 1:NoRotOrbs), NoRotOrbs)
                do l = 1, NoRotOrbs
                    do j = 1, l
                        FourIndInts(k,i,j,l) = FourIndInts(i,k,j,l)
                        FourIndInts(k,i,l,j) = FourIndInts(i,k,j,l)
                        FourIndInts(i,k,j,l) = FourIndInts(i,k,j,l)
                        FourIndInts(i,k,l,j) = FourIndInts(i,k,j,l)
                    end do
                end do
            end do
        end do
 
        deallocate(Temp4indints)
        call LogMemDeAlloc('Transform2ElIntsMemSave', Temp4indintsTag)
 
        call halt_timer(Transform2ElInts_Time)

    end subroutine Transform2ElIntsMemSave
   
! This is a transformation of the four index integrals for the ERlocalisation, in this only the <ii|ii> integrals are needed 
! therefore the process may be much simpler.

    subroutine Transform2ElIntsERlocal()

        integer :: i, j, a, b, g, d, m
        real(dp) :: t, Temp4indints(NoOrbs, NoOrbs)
        real(dp) :: Temp4indints02(NoOrbs)  

        call set_timer(Transform2ElInts_time, 30)

        ! Zero arrays from previous transform.

        TwoIndIntsER(:,:,:) = 0.0_dp

        ThreeIndInts01ER(:,:) = 0.0_dp
        ThreeIndInts02ER(:,:) = 0.0_dp

        FourIndIntsER(:) = 0.0_dp

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i, j, k, l, 0, 0)

        do d = 1, NoOrbs
            do b = 1, d
                Temp4indints(:,:) = 0.0_dp
                Temp4indints02(:) = 0.0_dp
                call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, UMATTemp01(:,:, b, d), &
                    NoOrbs, 0.0_dp, Temp4indints(:,:), NoOrbs)
                ! a -> m. Temp4indints(m,g) comes out of here.
                ! Want to transform g to m as well.

                do m = 1, NoOrbs
                    do g = 1, NoOrbs
                        Temp4indints02(m) = Temp4indints02(m)+(Temp4indints(m,g)*CoeffT1(g,m))
                    end do
                end do
                ! Now have Temp4indints(m,m) for each b and d.

                do m = 1, NoOrbs
                    TwoIndIntsER(b, d, m) = Temp4indints02(m)
                    TwoIndIntsER(d, b, m) = Temp4indints02(m)
                end do
            end do
        end do
        
! Now want to transform g to get one of the 3-transformed 4-index integrals <a m | m m>.
! These can be stored in 2-D arrays, as they can be specified by only m and z.

        do m = 1, NoOrbs
            do b = 1, NoOrbs
                do d = 1, NoOrbs
                    ThreeIndInts01ER(b,m) = ThreeIndInts01ER(b,m)+(TwoIndIntsER(b, d, m)*CoeffT1(d,m))
                end do
            end do
        end do
        ! ThreeIndInts01ER(z,m) is where z is alpha (a).


        TwoIndIntsER(:,:,:) = 0.0_dp
        do d = 1, NoOrbs
            do b = 1, NoOrbs
                Temp4indints(:,:) = 0.0_dp
                Temp4indints02(:) = 0.0_dp
                call dgemm('T','N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, CoeffT1(:,:), NoOrbs, UMATTemp01(:,:, b, d), &
                    NoOrbs, 0.0_dp, Temp4indints(:,:), NoOrbs)
                ! a -> m. Temp4indints(m,g) comes out of here.
                ! Want to transform g to m as well.

                do m = 1, NoOrbs
                    do g = 1, NoOrbs
                        Temp4indints02(m) = Temp4indints02(m)+(Temp4indints(m,g)*CoeffT1(g,m))
                    end do
                end do
                ! Now have Temp4indints(m,m) for each a and g.

                do m = 1, NoOrbs
                    TwoIndIntsER(b, d, m) = Temp4indints02(m)
                    TwoIndIntsER(d, b, m) = Temp4indints02(m)
                end do
            end do
        end do

! Now want to transform g to get one of the 3-transformed 4-index integrals <a m | m m>.
! These can be stored in 2-D arrays, as they can be specified by only m and z.

        do m = 1, NoOrbs
            do b = 1, NoOrbs
                do d = 1, NoOrbs
                    ThreeIndInts02ER(b,m) = ThreeIndInts02ER(b,m)+(TwoIndIntsER(b, d, m)*CoeffT1(d,m))
                end do
            end do
        end do
        ! ThreeIndInts02ER(z,m) is where z is beta (b).

! Find the <ii|ii> integrals, to calculate the potential energy.

        do m = 1, NoOrbs
            do a = 1, NoOrbs
                FourIndIntsER(m) = FourIndIntsER(m)+(ThreeIndInts01ER(a,m)*CoeffT1(a,m))
            end do
        end do

! ***************************
! Calc the potential energies for this iteration (with these transformed integrals).        

! This can be sped up by merging the calculations of the potentials with the transformations, but while 
! we are playing around with different potentials, it is simpler to keep these separate.
        
        PotEnergy = 0.0_dp
        TwoEInts = 0.0_dp
        PEInts = 0.0_dp
        call CalcPotentials()

        if (tPrintInts) call PrintIntegrals()
        if ((Iteration == 0).or.((.not.tNotConverged) .and. (Iteration > 1))) call WriteDoubHisttofile()
        if (tROHistSingExc .and. (Iteration == 0)) call WriteSingHisttofile()


! If doing Lagrange orthormalisations, find the change of the potential energy due to the orthonormality 
! of the orbitals...
        if (tLagrange) then
            PEOrtho = 0.0_dp
            do i = 1, NoOrbs
                do j = 1, NoOrbs
                    t = 0.0_dp
                    do a = 1, NoOrbs
                        t = CoeffT1(a,i)*CoeffT1(a,j)
                    end do
                    if (i == j) t = t-1.0_dp
                    PEOrtho = PEOrtho-Lambdas(i,j)*t
                    PotEnergy = PotEnergy-Lambdas(i,j)*t
                end do
            end do
        end if

        call halt_timer(Transform2ElInts_Time)

    end subroutine Transform2ElIntsERlocal

    subroutine CalcPotentials()

        ! Only temporarily like this, can tidy it up majorly

        integer :: i, j, k, l, Starti, Finishi
        real(dp) :: MaxTerm

        l  =  0
        if (tERLocalization .and. (.not.tStoreSpinOrbs)) then
            ERPotEnergy = 0.0_dp
            if (tRotateVirtOnly) then 
                Starti = NoOcc+1
                Finishi = NoOrbs
            elseif (tRotateOccOnly) then
                Starti = 1
                Finishi = NoOcc
            else
                Starti = 1
                Finishi = NoOrbs
            end if
            CoulPotEnergy = 0.0_dp
            OffDiagPotEnergy = 0.0_dp
            do i = Starti, Finishi
                    ERPotEnergy = ERPotEnergy+FourIndIntsER(i)
                    if (FourIndIntsER(i) < 0) then
                        call neci_flush(6)
                        call Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    end if
                    PotEnergy = PotEnergy+FourIndIntsER(i)
                    TwoEInts = TwoEInts+FourIndIntsER(i)
                    PEInts = PEInts+FourIndIntsER(i)
            end do
        elseif (tERLocalization) then
            ERPotEnergy = 0.0_dp
            PotEnergy = 0.0_dp
            if (tRotateVirtOnly) then 
                Starti = NoOcc+1
                Finishi = NoOrbs
            elseif (tRotateOccOnly) then
                Starti = 1
                Finishi = NoOcc
            else
                Starti = 1
                Finishi = NoOrbs
            end if
            do i = Starti, Finishi
                if (tStoreSpinOrbs) then
                    if (MOD(i,2) == 0) then
                        j = i-1
                    else
                        j = i+1
                    end if
                    ERPotEnergy = ERPotEnergy+FourIndInts(i,j,i,j)
                    if ((FourIndInts(i,j,i,j) < 0).or.(FourIndInts(j,i,j,i) < 0)) then
                        call neci_flush(6)
                        call Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    end if
                    PotEnergy = PotEnergy+FourIndInts(i,j,i,j)
                    TwoEInts = TwoEInts+FourIndInts(i,j,i,j)
                    PEInts = PEInts+FourIndInts(i,j,i,j)
                else
                    ERPotEnergy = ERPotEnergy+FourIndInts(i,i,i,i)
                    if ((FourIndInts(i,i,i,i) < 0)) then
                        call neci_flush(6)
                        call Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    end if
                    PotEnergy = PotEnergy+FourIndInts(i,i,i,i)
                    TwoEInts = TwoEInts+FourIndInts(i,i,i,i)
                    PEInts = PEInts+FourIndInts(i,i,i,i)
                end if
            end do
        end if

        if (tOffDiagSqrdMin.or.tOffDiagSqrdMax.or.tOffDiagMin.or.tOffdiagMax) then
            do l = 1, NoOrbs
                do j = 1, l-1
                    do k = 1, j-1
                        do i = 1, k-1
                            if (tOffDiagSqrdMin.or.tOffDiagSqrdMax) then
                                if (((i /= j) .and. (j /= l)) .and. ((i /= k).or.(j /= l))) then
                                    PotEnergy = PotEnergy+(FourIndInts(i,j,k,l)**2)
                                    TwoEInts = TwoEInts+(FourIndInts(i,j,k,l)**2)
                                    PEInts = PEInts+(FourIndInts(i,j,k,l)**2)
                                end if
                            end if
                            if (tOffDiagMin.or.tOffDiagMax) then
                                if (.not.((k == i).or.(j == l))) then
                                    PotEnergy = PotEnergy+FourIndInts(i,j,k,l)
                                    TwoEInts = TwoEInts+FourIndInts(i,j,k,l)
                                    PEInts = PEInts+FourIndInts(i,j,k,l)
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end if
      
        if (tDoubExcMin) then
            do i = 1, NoOrbs
                do j = 1, NoOrbs
                    do k = 1, i-1
                        if ((k == l) .and. (k == i)) cycle
                        do l = 1, j-1
                            if ((j == k) .and. (j == l)) cycle
                            if ((j == k) .and. (j == i)) cycle
                            if ((j == l) .and. (j == i)) cycle
                            PotEnergy = PotEnergy+(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            TwoEInts = TwoEInts+(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            PEInts = PEInts+(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                        end do
                    end do
                end do
            end do
        end if
 
        if (tOnePartOrbEnMax.or.tOneElIntMax) then
            do i = NoOcc+1, NoOrbs
                MaxTerm = 0.0_dp
                MaxTerm = TMAT2DRot(i,i)
                if (tOnePartOrbEnMax) then
                    do j = 1, NoOcc
                        MaxTerm = MaxTerm+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                    end do
                    MaxTerm = MaxTerm-EpsilonMin
                    MaxTerm = MaxTerm**OrbEnMaxAlpha
                end if
                PotEnergy = PotEnergy+MaxTerm
            end do
        end if

        if (tHijSqrdMin) then
            HijSqrdPotEnergy = 0.0_dp
            do i = NoOcc+1, NoOrbs
                do j = NoOcc+1, NoOrbs
                    if (j > i) then
                        PotEnergy = PotEnergy+(TMAT2DRot(i,j)**2)
                        HijSqrdPotEnergy = HijSqrdPotEnergy+(TMAT2DRot(i,j)**2)
                    end if
                end do
            end do
        end if

        if (tVirtCoulombMax) then
            ERPotEnergy = 0.0_dp
            ijOccVirtPotEnergy = 0.0_dp
            do i = 1, NoOrbs
                if (i <= NoOcc) then
                    do j = NoOcc+1, NoOrbs
                        ijOccVirtPotEnergy = ijOccVirtPotEnergy+FourIndInts(i,j,i,j)
                    end do
                end if
                if (i > NoOcc) then
                    ERPotEnergy = ERPotEnergy+FourIndInts(i,i,i,i)
                    do j = NoOcc+1, NoOrbs
                        if (j <= i) cycle
                        PotEnergy = PotEnergy+FourIndInts(i,j,i,j)
                        TwoEInts = TwoEInts+FourIndInts(i,j,i,j)
                    end do
                end if
            end do
        end if

        if (tHFSingDoubExcMax) then
            do i = 1, NoOcc
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do l = NoOcc+1, NoOrbs
                            PotEnergy = PotEnergy+(FourIndInts(i,j,k,l)**2)
                        end do
                        
                        ! Sing excitations <ij|ik> where i and j are occ, k virt.
                        PotEnergy = PotEnergy+(FourIndInts(i,j,i,k)**2)
                    end do
                end do
            end do
        end if

    end subroutine CalcPotentials

    subroutine FindTheForce()

        integer :: m, z, i, j, k, l, a, Symm, w, x, y, SymMin
        real(dp) :: OffDiagForcemz, DiagForcemz, OneElForcemz, LambdaTerm1, LambdaTerm2
        real(dp) :: NonDerivTerm, DerivPot
        logical :: leqm, jeqm, keqm
      
        ! Running over m and z, covers all matrix elements of the force matrix (derivative 
        ! of equation we are minimising, with respect to each translation coefficient) filling 
        ! them in as it goes.
        call set_timer(FindtheForce_time, 30)
        
        DerivCoeff(:,:) = 0.0_dp
        Force = 0.0_dp
        ForceInts = 0.0_dp
        OrthoForce = 0.0_dp
        OffDiagForceMZ  =  0
        SymMin  =  0
        i  =  0

        ! If the orbitals are being separated, do this whole loop twice, once for occupied and once for virtual
        ! i.e w  =  1, 2. Otherwise do them all at once.
        do w = MinOccVirt, MaxOccVirt
            if (w == 1) then
                SymMin = 1
                MinMZ = 1
                if (tSeparateOccVirt) then
                    MaxMZ = NoOcc
                else
                    MaxMZ = NoOrbs
                end if
            else
                SymMin = 9
                MinMZ = NoOcc+1
                MaxMZ = NoOrbs
            end if
! If we are localising the occupied and virtual orbitals separately, the above block ensures that we loop over
! first the occupied then the virtual.  If we are not separating the orbitals we just run over all orbitals.

            do m = MinMZ, MaxMZ
                if (tStoreSpinOrbs) then
                    SymM = int(G1(SymLabelList2_rot(m))%sym%S)
                else
                    SymM = int(G1(SymLabelList2_rot(m)*2)%sym%S)
                end if
                do z = SymLabelCounts2_rot(1, SymM+SymMin), &
                        (SymLabelCounts2_rot(1, SymM+SymMin) + &
                            SymLabelCounts2_rot(2, SymM+SymMin)-1)
 
                    ! Find the force on a coefficient c(m,z). 
                    OffDiagForcemz = 0.0_dp
                    ! OffDiagForce is from any of the OffDiagMin/Max (Sqrd or not), or the double/single excitation
                    ! max/min, as only one of these terms may be used at once.

                    DiagForcemz = 0.0_dp
                    ! DiagForce includes ER localisation, and the coulomb terms <ij|ij>.

                    OneElForcemz = 0.0_dp
                    ! OneElForce includes that from the one electron integrals <i|h|j> and the one particle orbital 
                    ! energies.

                    ! DIAG TERMS
                    ! Maximise <ii|ii>, self interaction terms. 
                    if (tERLocalization .and. (.not.tStoreSpinOrbs)) then
                        DiagForcemz = DiagForcemz+(2*ThreeIndInts01ER(z,m))+(2*ThreeIndInts02ER(z,m))
                        ! Derivative of <ii|ii> only non-zero when i = m.
                        ! each of the four terms then correspond to zeta  =  a, b, g, then d in the unrotated basis.
                    elseif (tERLocalization) then
                        ! Looking at <ij|ij> terms where j = i+1 (i.e. i is alpha of spin orbital and j is beta - or vice versa).
                        if (tStoreSpinOrbs) then
                            if (MOD(m,2) == 0) then      ! m  =  j
                                i = m-1
                                DiagForcemz = DiagForcemz+ThreeIndInts01(m,i,i,z)+ThreeIndInts01(z,i,i,m)+ &
                                    ThreeIndInts01(i,m,z,i)+ThreeIndInts01(i,z,m,i)
                            else
                                j = m+1
                                DiagForcemz = DiagForcemz+ThreeIndInts01(m,j,j,z)+ThreeIndInts01(z,j,j,m)+ &
                                    ThreeIndInts01(j,m,z,j)+ThreeIndInts01(j,z,m,j)
                            end if
                        else
                            DiagForcemz = DiagForcemz+ThreeIndInts01(m,m,m,z)+ThreeIndInts02(m,m,m,z)+ &
                                ThreeIndInts03(m,m,m,z)+ThreeIndInts04(m,m,m,z)
                        end if
                        ! First term when m = i and z = a, second when m = i and z = g.
                    end if

                    ! Maximise <ij|ij>, coulomb terms, where i<j, i occ or virt, j virt only.
                    if (tVirtCoulombMax) then
                        do i = 1, NoOrbs
                            if (i == m) then
                                do j = NoOcc+1, NoOrbs
                                    if (j <= i) cycle        ! i<j.
                                    DiagForcemz = DiagForcemz+ThreeIndInts01(m,j,j,z)+ThreeIndInts03(m,j,j,z)
                                    ! First term for when m = i and z = a, second when m = i and z = g.
                                end do
                            end if
                            if ((m > NoOcc) .and. (m > i)) DiagForcemz = &
                                                            DiagForcemz+ThreeIndInts02(i,i,m,z)+ThreeIndInts04(i,i,m,z)
                            ! This only contributes when j = m (no point in running over all j.
                            ! First term when m = j and z = b, second when m = j and z = d.
                        end do
                    end if

                    ! ONE ELECTRON TERMS
                    ! Minimise |<i|h|j>|^2 where either one or bot of i and j are virtual, but i<j.
                    if (tHijSqrdMin) then
                        do j = NoOcc+1, NoOrbs
                            if (m /= j) OneElForcemz = OneElForcemz+(2*TMAT2DRot(m,j)*TMAT2DPartRot02(z,j))
                            ! m = i and z = a.
                        end do
                        do i = NoOcc+1, NoOrbs
                            if (m /= i) OneElForcemz = OneElForcemz+(2*TMAT2DRot(i,m)*TMAT2DPartRot01(i,z))
                            ! m = j and z = b
                        end do
                    end if

                    ! OnePartOrbEnMax ; Maximisie sum_i [E_i - E_min]^Alpha
                    ! where E_i  =  <i|h|i> + sum_j <ij||ij> and E_min is either E_LUMO (rotating virtual only) or the chemical 
                    ! potential (midway between LUMO and HOMO, when rotating all), Alpha specified in input.
                    ! The derivative of the one part orb energies is then Alpha * NonDerivTerm * DerivPot^(Alpha-1)  

                    ! OneElIntMax ; Maximise <i|h|i>
                    if (tOnePartOrbEnMax.or.tOneElIntMax) then
                        do i = NoOcc+1, NoOrbs
                            DerivPot = 0.0_dp
                            DerivPot = DerivPot+TMAT2DPartRot02(z,m)+TMAT2DPartRot01(m,z)
                            ! First term when m = i and z = a, second when m = i and z = b.
                            ! This is all that is needed for OneElIntMax

                            if (tOnePartOrbEnMax) then
                                NonDerivTerm = 0.0_dp
                                if (OrbEnMaxAlpha /= 1.0_dp) then 
                                    ! The non-derived term in the chain rule, <i|h|i> + sum_j <ij||ij> - E_min.
                                    NonDerivTerm = NonDerivTerm+TMAT2DRot(i,i)-EpsilonMin
                                    do j = 1, NoOcc
                                        NonDerivTerm = NonDerivTerm+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                                    end do
                                    NonDerivTerm = OrbEnMaxAlpha*(NonDerivTerm**(OrbEnMaxAlpha-1))
                                else
                                    ! If Alpha  =  1, the NonDerivTerm will be raised to the power of 0, thus always 1.
                                    NonDerivTerm = 1.0
                                end if
                                if (i == m) then
                                    do j = 1, NoOcc
                                        DerivPot = DerivPot+(2*ThreeIndInts01(m,j,j,z))-ThreeIndInts01(j,j,m,z)+ &
                                            (2*ThreeIndInts03(m,j,j,z))-ThreeIndInts03(m,j,z,j)
                                        ! First part is for when m = i and z = a, the second is for when m = i and z = g
                                    end do
                                end if
                                ! When m = j, for a particular i.
                                ! m and z run only over virtual, and j is over occupied. m will never  =  j.
                            else
                                NonDerivTerm = 1.0_dp
                            end if

                            OneElForcemz = OneElForcemz+(NonDerivTerm*DerivPot)
                        end do
                    end if

                    ! OFFDIAGTERMS
                    ! Maximises the square of the single and double excitation integrals connected to the HF.
                    ! I.e maximises <ij|kl> where i, j are occupied and k, l are virtual (doubles), except k may be occuppied if
                    ! equal to i (<ij|il> singles).
                    ! Currently this is only used for rotating virtual only, so m can only equal k or l.
                    if (tHFSingDoubExcMax) then
                        do i = 1, NoOcc
                            do j = 1, NoOcc
                                do k = NoOcc+1, NoOrbs
                                    if (k == m) then
                                        do l = NoOcc+1, NoOrbs
                                            OffDiagForcemz = OffDiagForcemz+(2*FourIndInts(i,j,m,l)*ThreeIndInts03(i,j,l,z))
                                            ! m = k and z = g.
                                        end do
                                    end if
                                    
                                    OffDiagForcemz = OffDiagForcemz+(2*FourIndInts(i,j,k,m)*ThreeIndInts04(i,k,j,z))
                                    ! m = l and z = d. 
                                    
                                    ! Sing excitations <ij|il> where i and j are occ, l virt.
                                    OffDiagForcemz = OffDiagForcemz+(2*FourIndInts(i,j,i,m)*ThreeIndInts04(i,i,j,z))
                                    ! m = l
                                end do
                            end do
                        end do
                    end if

                    ! OffDiag Sqrd/notSqrd Min/Max treats the elements <ij|kl>
                    ! i<k and j<l.
                    if (tOffDiagSqrdMin.or.tOffDiagSqrdMax.or.tOffDiagMin.or.tOffDiagMax.or.tDoubExcMin) then
                        do l = 1, NoOrbs
                            if (l == m) then
                                leqm = .true.
                            else
                                leqm = .false.
                            end if
                            do j = 1, l-1                        
                                if (j == m) then
                                    jeqm = .true.
                                else
                                    jeqm = .false.
                                end if
                                do k = 1, j-1                 
                                    if (k == l) cycle
                                    if (k == m) then
                                        keqm = .true.
                                    else
                                        keqm = .false.
                                    end if
                                    ! only i with symmetry equal to j x k x l will have integrals with overall
                                    ! symmetry A1 and therefore be non-zero.


                                    ! Running across i, ThreeIndInts01 only contributes 
                                    if ((m <= k-1) .and. (m /= j) .and. ((i /= k).or.(j /= l))) then
                                        if (tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz = OffDiagForcemz+2* &
                                            (FourIndInts02(j,k,l,m)*ThreeIndInts01(k,j,l,z))
                                        if (tOffDiagMin.or.tOffDiagMax) OffDiagForcemz = OffDiagForcemz+ThreeIndInts01(k,j,l,z)
                                        if (tDoubExcMin) OffDiagForcemz = OffDiagForcemz+ThreeIndInts01(k,j,l,z)- &
                                            ThreeIndInts01(l,j,k,z)
                                    end if

                                    if (jeqm) then
                                        do i = 1, k-1
                                            if ((i /= j) .and. ((i /= k).or.(j /= l))) then
                                                if (tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz = OffDiagForcemz+2* &
                                                    (FourIndInts(i,j,k,l)*ThreeIndInts02(i,k,l,z))
                                                if (tOffDiagMin.or.tOffDiagMax) &
                                                    OffDiagForcemz = OffDiagForcemz+ThreeIndInts02(i,k,l,z)
                                                if (tDoubExcMin) OffDiagForcemz = OffDiagForcemz+(ThreeIndInts02(i,k,l,z))- &
                                                    ThreeIndInts02(i,l,k,z)
                                            end if
                                        end do
                                    end if

                                    if (keqm) then
                                        do i = 1, k-1
                                            if ((i /= j) .and. ((i /= k).or.(j /= l))) then
                                                if (tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz = OffDiagForcemz+2* &
                                                    (FourIndInts(i,j,k,l)*ThreeIndInts03(i,j,l,z))
                                                if (tOffDiagMin.or.tOffDiagSqrdMax) OffDiagForcemz = OffDiagForcemz+ &
                                                    ThreeIndInts03(i,j,l,z)
                                                if (tDoubExcMin) OffDiagForcemz = OffDiagForcemz+(ThreeIndInts03(i,j,l,z))- &
                                                    ThreeIndInts03(i,j,z,l)
                                            end if
                                        end do
                                    end if

                                    if (leqm) then
                                        do i = 1, k-1
                                            if ((i /= j) .and. ((i /= k).or.(j /= l))) then
                                                if (tOffDiagSqrdMin.or.tOffDiagSqrdMin) OffDiagForcemz = OffDiagForcemz+2* &
                                                    (FourIndInts(i,j,k,l)*ThreeIndInts04(i,k,j,z))
                                                if (tOffDiagMin.or.tOffDiagMax) &
                                                    OffDiagForcemz = OffDiagForcemz+ThreeIndInts04(i,k,j,z)
                                                if (tDoubExcMin) OffDiagForcemz = OffDiagForcemz+(ThreeIndInts04(i,k,j,z))- &
                                                    ThreeIndInts04(i,z,j,k)
                                            end if
                                        end do
                                    end if
                                end do
                            end do
                        end do
                    end if
                    DerivCoeff(z,m) = (MaxMinFac*OffDiagWeight*OffDiagForcemz)+(DiagMaxMinFac*DiagWeight*DiagForcemz)+ &
                        (OneElMaxMinFac*OneElWeight*OneElForcemz)
                    Force = Force+ABS(DerivCoeff(z,m))
                end do
            end do
        end do

        Force = Force/real(NoOrbs**2, dp)

! Calculate the derivatives of orthogonalisation condition.
! Have taken this out of the m and z loop to make the shake faster, but can put it back in if start using it a lot.
        if (tLagrange) then
            do x = MinMZ, MaxMZ
                m = SymLabelList2_rot(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM = int(G1(m*2)%sym%S)
                do y = SymLabelCounts2_rot(1, SymM+SymMin), &
                        (SymLabelCounts2_rot(1, SymM+SymMin) + &
                            SymLabelCounts2_rot(2, SymM+SymMin)-1)
                    z = SymLabelList2_rot(y)
      
                    LambdaTerm1 = 0.0_dp
                    LambdaTerm2 = 0.0_dp
                    
                    do j = 1, NoOrbs
                        LambdaTerm1 = LambdaTerm1+(Lambdas(m,j)*CoeffT1(z,j))
                        LambdaTerm2 = LambdaTerm2+(Lambdas(j,m)*CoeffT1(z,j))
                    end do

! DerivCoeff is 'the force'.  I.e. the derivative of |<ij|kl>|^2 with 
! respect to each transformation coefficient.  It is the values of this matrix that will tend to 0 as
! we minimise the sum of the |<ij|kl>|^2 values.
! With the Lagrange keyword this includes orthonormality conditions, otherwise it is simply the unconstrained force.
                    DerivCoeff(z,m) = (2*OffDiagForcemz)-LambdaTerm1-LambdaTerm2
                    OrthoForce = OrthoForce-LambdaTerm1-LambdaTerm2
                end do
            end do
 
! If doing a Lagrange calc we also need to find the force on the lambdas to ensure orthonormality...
            OrthoForce = OrthoForce/real(NoOrbs**2, dp)
            DerivLambda(:,:) = 0.0_dp
            do i = 1, NoOrbs
                do j = 1, i
                    do a = 1, NoOrbs
                        DerivLambda(i,j) = DerivLambda(i,j)+CoeffT1(a,i)*CoeffT1(a,j)
                    end do
                    DerivLambda(j,i) = DerivLambda(i,j)
                end do
            end do
            do i = 1, NoOrbs
                DerivLambda(i,i) = DerivLambda(i,i)-1.0_dp
            end do
        end if

        call halt_timer(FindtheForce_Time)

    end subroutine FindTheForce
    
    subroutine UseTheForce()

        ! This routine takes the old translation coefficients and Lambdas and moves them by a timestep in the direction 
        ! of the calculated force.

        integer :: m, w, z, i, j, Symm, SymMin
        real(dp) :: NewCoeff, NewLambda

        DistCs = 0.0_dp 
    
        do w = MinOccVirt, MaxOccVirt
            if (w == 1) then
                SymMin = 1
                MinMZ = 1
                if (tSeparateOccVirt) then
                    MaxMZ = NoOcc
                else
                    MaxMZ = NoOrbs
                end if
            else
                SymMin = 9
                MinMZ = NoOcc+1
                MaxMZ = NoOrbs
            end if

            do m = MinMZ, MaxMZ
                if (tStoreSpinOrbs) then
                    SymM = int(G1(SymLabelList2_rot(m))%sym%S)
                else
                    SymM = int(G1(SymLabelList2_rot(m)*2)%sym%S)
                end if

                ! Symmetry requirement that z must be from the same irrep as m.
                do z = SymLabelCounts2_rot(1, SymM+SymMin), &
                        (SymLabelCounts2_rot(1, SymM+SymMin) + &
                            SymLabelCounts2_rot(2, SymM+SymMin)-1)

                    ! Only coeffs with sym of m and z the same have non-zero coeffs.    
                    NewCoeff = 0.0_dp
                    NewCoeff = CoeffT1(z,m)-(TimeStep*DerivCoeff(z,m))
                    DistCs = DistCs+abs(TimeStep*DerivCoeff(z,m))
                    CoeffT1(z,m) = NewCoeff
                end do
            end do
        end do

        DistCs = DistCs/(real(NoOrbs**2, dp))

        if (tLagrange) then

            DistLs = 0.0_dp
            LambdaMag = 0.0_dp
            do i = 1, NoOrbs
                do j = 1, NoOrbs
                    NewLambda = 0.0_dp
                    NewLambda = Lambdas(i,j)-(TimeStep*DerivLambda(i,j))  ! Timestep must be specified in the input file.
                    DistLs = DistLs+abs(TimeStep*DerivLambda(i,j))
                    Lambdas(i,j) = NewLambda
                    LambdaMag = LambdaMag+abs(NewLambda)
                end do
            end do
            DistLs = DistLs/(real(NoOrbs**2, dp))
            LambdaMag = LambdaMag/(real(NoOrbs**2, dp))

        end if

    end subroutine UseTheForce
   
    subroutine TestOrthonormality()
        integer :: i, j
        real(dp) :: OrthoNormDP

        OrthoNorm = 0.0_dp
        do i = 1, NoOrbs
            do j = 1, i
                OrthoNormDP = 0.0_dp
                OrthoNormDP = Dot_Product(CoeffT1(:,i), CoeffT1(:,j))
                OrthoNorm = OrthoNorm+ABS(OrthoNormDP)
            end do
        end do
        OrthoNorm = OrthoNorm-real(NoOrbs, dp)
        OrthoNorm = (OrthoNorm*2.0_dp)/real((NoOrbs*(NoOrbs+1.0_dp)), dp)

    end subroutine TestOrthonormality

    subroutine TestForConvergence()

        ! This just tests the convergence on the grounds that the force is
        ! smaller that the input parameter: ConvergedForce

        if (tLagrange) then
            if ((abs(Force) < ConvergedForce) .and. (abs(OrthoForce) < ConvergedForce)) then
                tNotConverged = .false.
            end if
        elseif (tROIteration) then
            if (Iteration == ROIterMax) then
                tNotConverged = .false.
            end if
        elseif (abs(TotCorrectedForce) < ConvergedForce) then
            tNotConverged = .false.
        end if
! if an ROIteration value is specified, use this to specify the end of the orbital rotation, otherwise use the 
! conversion limit (ConvergedForce).

    end subroutine TestForConvergence

    subroutine ShakeConstraints()

        ! DerivCoeff(k,a) is the unconstrained force on the original coefficients (CoeffT1(a,k)). 

        integer :: w, l, a, m, ShakeIteration, ConvergeCount, SymM, SymMin
        real(dp) :: TotCorConstraints, TotConstraints, TotLambdas
        real(dp) :: TotUncorForce, TotDiffUncorCoeffs, TotDiffCorCoeffs
        logical :: tShakeNotConverged
        integer, save :: shake_io

        if (Iteration == 1) then
            shake_io  =  get_free_unit()
            open(shake_io, file='SHAKEstats', status='unknown')
            write(shake_io,'(A20, 4A35, A20)') 'Shake Iteration','Sum Lambdas','Total of corrected forces', &
                & 'Sum unconstrained constraints',&
                                        &'Sum corrected constraints','Converge count'
        end if
        if (Mod(Iteration, 10) == 0) write(shake_io,*) 'Orbital rotation iteration  =  ', Iteration


        ShakeIteration = 0
        tShakeNotConverged = .true.

! Before we start iterating, take the current coefficients and find the derivative of the constraints with respect to them.

        call CalcDerivConstr(CoeffT1, DerivConstrT1)

! Then find the coefficients at time t2, when moved by the completely unconstrained force and the values of the each 
! constraint at these positions.

        Correction(:,:) = 0.0_dp
        call FindandUsetheForce(TotUncorForce, TotDiffUncorCoeffs, CoeffUncorT2)


        call CalcConstraints(CoeffUncorT2, Constraint, TotConstraints)

        call set_timer(Shake_Time, 30)

        if (tShakeDelay) then
            if (Iteration < ShakeStart) then
                ShakeIterMax = 1
            else
                ShakeIterMax = ShakeIterInput
            end if
        end if

        ! Actually starting the calculation.
        do while (tShakeNotConverged)

            ShakeIteration = ShakeIteration+1
            
            ForceCorrect(:,:) = 0.0_dp          ! Zeroing terms that are re-calculated each iteration.
            CoeffCorT2(:,:) = 0.0_dp
            ConstraintCor(:) = 0.0_dp
            DerivConstrT2(:,:,:) = 0.0_dp
            TotLambdas = 0.0_dp
            TotCorrectedForce = 0.0_dp
            TotDiffCorCoeffs = 0.0_dp

            if (ShakeIteration /= 1) then
                call UpdateLambdas()
            end if
            
            ShakeLambdaNew(:) = 0.0_dp

            ! For a particular set of coefficients cm:
            ! Force(corrected) = Force(uncorrected)-Lambdas.DerivConstrT1
            ! Use these derivatives, and the current lambdas to find the trial corrected force.
            ! Then use this to get the (trial) shifted coefficients.

            ! Use the lambdas of this iteration to calculate the correction to the force due to the constraints.
            Correction(:,:) = 0.0_dp
                
            do w = MinOccVirt, MaxOccVirt
                if (w == 1) then
                    SymMin = 1
                    MinMZ = 1
                    if (tSeparateOccVirt) then
                        MaxMZ = NoOcc
                    else
                        MaxMZ = NoOrbs
                    end if
                else
                    SymMin = 9
                    MinMZ = NoOcc+1
                    MaxMZ = NoOrbs
                end if

                do m = MinMZ, MaxMZ
                    if (tStoreSpinOrbs) then
                        SymM = int(G1(SymLabelList2_rot(m))%sym%S)
                    else
                        SymM = int(G1(SymLabelList2_rot(m)*2)%sym%S)
                    end if
                    do a = SymLabelCounts2_rot(1, SymM+SymMin), &
                            (SymLabelCounts2_rot(1, SymM+SymMin) + &
                                SymLabelCounts2_rot(2, SymM+SymMin)-1)
                        do l = 1, TotNoConstraints
                            Correction(a,m) = Correction(a,m)+(ShakeLambda(l)*DerivConstrT1(a, m, l)) 
                        end do
                    end do
                end do
            end do

            call FindandUsetheForce(TotCorrectedForce, TotDiffCorCoeffs, CoeffCorT2)

            ! Use these new shifted coefficients to calculate the derivative of the constraints 
            ! (at time t2).
           
            call CalcDerivConstr(CoeffCorT2, DerivConstrT2) 
            
            
! Test for convergence, if convergence is reached, make the new coefficients the original ones to start the whole process again.
! Then exit out of this do loop and hence the subroutine.
            call TestShakeConvergence(ConvergeCount, TotCorConstraints, ShakeIteration, tShakeNotConverged)

! If the convergence criteria is met, exit out of this subroutine, a rotation has been made which keeps the coefficients 
! orthogonal.

! and to SHAKEstats file:
            call neci_flush(6)
            call neci_flush(shake_io)
            if (Mod(Iteration, 10) == 0) then
                write(shake_io,'(I20, 4F35.20, I20)') ShakeIteration, TotLambdas, TotCorrectedForce, TotConstraints, &
                    TotCorConstraints, ConvergeCount 
            end if

! If the convergence criteria is not met, use either the full matrix inversion method to 
!find a new set of lambdas, or the shake algorithm 
! (in which case SHAKEAPPROX is required in the system block of the input).

            if (tShakeApprox .and. tShakeNotConverged) then
                call ShakeApproximation()
            elseif (tShakeNotConverged) then
                call FullShake()
            else
                DistCs = TotDiffCorCoeffs
            end if
   
    end do

    call halt_timer(Shake_Time)

    end subroutine ShakeConstraints

    subroutine CalcDerivConstr(CurrCoeff, DerivConstr)

        ! This calculates the derivative of each of the orthonormalisation
        ! constraints, l, with respect to each set of coefficients cm.

        integer :: l, i, j, a
        real(dp) :: CurrCoeff(NoOrbs, NoOrbs)
        real(dp) :: DerivConstr(NoOrbs, NoOrbs, TotNoConstraints)

        call set_timer(CalcDerivConstr_Time, 30)

            DerivConstr(:,:,:) = 0.0_dp
            do l = 1, TotNoConstraints
                i = lab(1,l)
                j = lab(2,l)
                if (i == j) then
                    do a = 1, NoOrbs
                        DerivConstr(a, i, l) = CurrCoeff(a,i)*2
                    end do
                else
                    do a = 1, NoOrbs
                        DerivConstr(a, j, l) = CurrCoeff(a,i) 
                    end do
                    do a = 1, NoOrbs
                        DerivConstr(a, i, l) = CurrCoeff(a,j)
                    end do
                end if
                ! DerivConstrT1 stays the same throughout the iterations
            end do

        call halt_timer(CalcDerivConstr_Time)

    end subroutine CalcDerivConstr

    subroutine FindandUsetheForce(TotForce, TotDiffCoeffs, CoeffT2)

! This takes the current lambdas with the derivatives of the constraints and calculates a force
! for each cm, with an orthonormalisation correction.
! This is then used to rotate the coefficients by a defined timestep.

        integer :: a, m, Symm, w, SymMin, TempMaxOccVirt
        real(dp) :: TotForce, TotDiffCoeffs, CoeffT2(NoOrbs, NoOrbs)

        call set_timer(findandusetheforce_time, 30)

        if (tSeparateOccVirt) then
            TempMaxOccVirt = 2
        else
            TempMaxOccVirt = 1
        end if

        do w = 1, TempMaxOccVirt
            ! The force will be zero on those coefficients not being mixed,
            ! but still want to run over all, so that the diagonal 1 values
            ! are maintained.
            if (w == 1) then
                SymMin = 1
                MinMZ = 1
                if (tSeparateOccVirt) then
                    MaxMZ = NoOcc
                else
                    MaxMZ = NoOrbs
                end if
            else
                SymMin = 9
                MinMZ = NoOcc+1
                MaxMZ = NoOrbs
            end if

            do m = MinMZ, MaxMZ
                if (tStoreSpinOrbs) then
                    SymM = int(G1(SymLabelList2_rot(m))%sym%S)
                else
                    SymM = int(G1(SymLabelList2_rot(m)*2)%sym%S)
                end if
                do a = SymLabelCounts2_rot(1, SymM+SymMin), &
                        (SymLabelCounts2_rot(1, SymM+SymMin) + &
                            SymLabelCounts2_rot(2, SymM+SymMin)-1)
!               
                ! FIND THE FORCE 
                    ! find the corrected force. (in the case where the uncorrected force 
                    !is required, correction is set to 0.
                    ! DerivCoeff(m,a) is the derivative of the relevant potential energy w.r.t 
                    !cm without any constraints (no lambda terms).
                    ! ForceCorrect is then the latest force on coefficients.  This is 
                    !iteratively being corrected so that
                    ! it will finally move the coefficients so that they remain orthonormal.
                
                ! use THE FORCE
                    ForceCorrect(a,m) = DerivCoeff(a,m)-Correction(a,m)
                    CoeffT2(a,m) = CoeffT1(a,m)-(TimeStep*ForceCorrect(a,m))
                    ! Using the force to calculate the coefficients at time T2
                    ! (hopefully more orthonomal than those calculated in the
                    ! previous iteration).
                    
                    ! Calculate parameters for printing
                    TotForce = TotForce+ABS(ForceCorrect(a,m))
                    TotDiffCoeffs = TotDiffCoeffs+ABS(CoeffT2(a,m)-CoeffT1(a,m))
                end do
            end do
        end do

        TotForce = TotForce/(real(NoOrbs**2, dp))

        call halt_timer(findandusetheforce_time)

    end subroutine FindandUsetheForce

    subroutine CalcConstraints(CurrCoeff, Constraint, TotConstraints)  

        ! This calculates the value of each orthonomalisation constraint, using the shifted coefficients.
        ! Each of these should tend to 0 when the coefficients become orthonomal.

        integer :: l, i, j
        real(dp) :: CurrCoeff(NoOrbs, NoOrbs), TotConstraints, Constraint(TotNoConstraints) 

            TotConstraints = 0.0_dp
            do l = 1, TotNoConstraints
                i = lab(1,l)
                j = lab(2,l)
                if (i == j) then
                    Constraint(l) = Dot_Product(CurrCoeff(:,i), CurrCoeff(:,j))-1.0_dp
                else
                    Constraint(l) = Dot_Product(CurrCoeff(:,i), CurrCoeff(:,j))
                    ! Each of these components should tend towards 0 when the coefficients become orthonormal.
                end if
                TotConstraints = TotConstraints+ABS(Constraint(l))
            end do
     
    end subroutine CalcConstraints

    subroutine FullShake()

        ! This method calculates the lambdas by solving the full matrix
        ! equation.

        integer :: l, n, m, info, ipiv(TotNoConstraints)
        character(len=*), parameter :: this_routine = 'FullShake'

        call set_timer(FullShake_Time, 30)

! FULL MATRIX INVERSION METHOD

! Calculate matrix from the derivatives of the constraints w.r.t the the coefficients at t1 and t2. I.e. the initial 
! coefficients and those that have been moved by the corrected force.

            DerivConstrT1T2(:,:) = 0.0_dp
            do l = 1, TotNoConstraints
                do n = 1, TotNoConstraints
                    do m = 1, NoOrbs
                        ! Product of constraint i, j at time t1, mult by constraint l, n.
                        ! Add these over all m for a specific constraints to get matrix elements
                        DerivConstrT1T2(n,l) = DerivConstrT1T2(n,l)+(Dot_Product(DerivConstrT2(:, m, l), DerivConstrT1(:, m, n)))
                    end do
                end do
            end do       ! have filled up whole matrix

! Invert the matrix to calculate the lambda values.
! LU decomposition.
            call dgetrf(TotNoConstraints, TotNoConstraints, DerivConstrT1T2, TotNoConstraints, ipiv, info)
            if (info /= 0) then
                write(6,*) 'info ', info
                call Stop_All(this_routine,"The LU decomposition of matrix inversion failed...")
            end if

            do n = 1, TotNoConstraints
                ShakeLambdaNew(n) = Constraint(n)/(TimeStep*(-1))
            end do
            ! These are actually still the constraint values, but now Lambda(n)
            ! can go into dgetrs as the constraints (B in AX = B), and come out
            ! as the computed lambdas (X).

            call dgetrs('N', TotNoConstraints, 1, DerivConstrT1T2, TotNoConstraints, ipiv, ShakeLambdaNew, TotNoConstraints, info)
            if (info /= 0) call Stop_All(this_routine,"Error in dgetrs, solving for the lambdas...")

        call halt_timer(FullShake_Time)

    end subroutine FullShake

    subroutine ShakeApproximation()

        ! This is an approximation in which only the diagonal elements are
        ! considered in the matrix of the derivative of the constraints
        ! DerivConstrT1T2.

        integer :: m, l

        ! Use 'shake' algorithm in which the iterative scheme is applied to
        ! each constraint in succession.
        write(6,*) 'DerivConstrT1T2Diag calculated from the shake approx'
        
        DerivConstrT1T2Diag(:) = 0.0_dp
        do l = 1, TotNoConstraints 
            do m = 1, NoOrbs
                DerivConstrT1T2Diag(l) = DerivConstrT1T2Diag(l)+Dot_Product(DerivConstrT2(:, m, l), DerivConstrT1(:, m, l))
            end do
            ShakeLambdaNew(l) = Constraint(l)/((-1)*TimeStep*DerivConstrT1T2Diag(l))
            write(6,*) DerivConstrT1T2Diag(l)
        end do

    end subroutine ShakeApproximation

    subroutine UpdateLambdas()

        ! Use damping to update the lambdas, rather than completely replacing
        ! them with the new values.

        integer :: l

        do l = 1, TotNoConstraints
            ShakeLambda(l) = ShakeLambdaNew(l)
        end do

    end subroutine UpdateLambdas

    subroutine TestShakeConvergence(ConvergeCount, TotCorConstraints, ShakeIteration, tShakeNotConverged)  

    ! This calculates the value of each orthonomalisation constraint using the corrected coefficients.
    ! Each of these should tend to 0 when the coefficients become orthonomal.
    ! CovergeCount counts the number of constraints that individually have values below the specified
    ! convergence criteria.  If this  =  0, the shake is converged, else keep iterating.

        integer :: l, i, j, m, a, ConvergeCount
        real(dp) :: TotCorConstraints
        integer :: ShakeIteration
        logical :: tShakeNotConverged

        TotCorConstraints = 0.0_dp
        ConvergeCount = 0
        ConstraintCor(:) = 0.0_dp
        do l = 1, TotNoConstraints
            i = lab(1,l)
            j = lab(2,l)
            if (i == j) then
                ConstraintCor(l) = Dot_Product(CoeffCorT2(:,i), CoeffCorT2(:,j))-1.0_dp
            else
                ! Each of these components should tend towards 0 when the
                ! coefficients become orthonormal.
                ConstraintCor(l) = Dot_Product(CoeffCorT2(:,i), CoeffCorT2(:,j))
            end if
            
            TotCorConstraints = TotCorConstraints+ABS(ConstraintCor(l))
            ! Sum of all Contraint components - indication of overall
            ! orthonormality.
    
            if (ABS(ConstraintCor(l)) > ShakeConverged) ConvergeCount = ConvergeCount+1
            ! Count the number of constraints which are still well above 0.
            
        end do
        
        if (tShakeIter) then
            if (ShakeIteration == ShakeIterMax) then
                do m = 1, NoOrbs
                    do a = 1, NoOrbs
                        CoeffT1(a,m) = CoeffCorT2(a,m)
                    end do
                end do
                tShakeNotConverged = .false.
            end if
        elseif (ConvergeCount == 0) then
           tShakeNotConverged = .false.

            ! If convergence is reached, make the new coefficients coeff, to
            ! start the rotation iteration again.

            do m = 1, NoOrbs
                do a = 1, NoOrbs
                    CoeffT1(a,m) = CoeffCorT2(a,m)
                end do
            end do
        end if

    end subroutine TestShakeConvergence

    subroutine FinalizeNewOrbs()

! At the end of the orbital rotation, have a set of coefficients CoeffT1 which transform 
! the HF orbitals into a set of linear
! combinations ui which minimise |<ij|kl>|^2.  This is the final subroutine after 
! all iterations (but before the memory deallocation)
! that calculates the final 4 index integrals to be used in the NECI calculation.

        use sym_mod, only: GenSymStatePairs

        integer :: i, a, j
        real(dp) :: TotGSConstraints, GSConstraint(TotNoConstraints), CoeffTemp(SpatOrbs, SpatOrbs)
        
        ! First need to do a final explicit orthonormalisation. The orbitals
        ! are very close to being orthonormal, but not exactly. Need to make
        ! sure they are exact orthonormal using Gram Schmit.
        if (tStoreSpinOrbs .and. (.not.tMaxHLGap)) then
            CoeffTemp(:,:) = 0.0_dp
            do i = 1, SpatOrbs
                do j = 1, SpatOrbs
                    CoeffTemp(i,j) = CoeffT1(2*i, 2*j)
                end do
            end do

            call GRAMSCHMIDT(CoeffTemp, SpatOrbs)

            CoeffT1(:,:) = 0.0_dp
            do i = 1, SpatOrbs
                do j = 1, SpatOrbs
                    CoeffT1(2*i, 2*j) = CoeffTemp(i,j)
                    CoeffT1((2*i)-1,(2*j)-1) = CoeffTemp(i,j)
                end do
            end do
        elseif (.not.tMaxHLGap) then
            call GRAMSCHMIDT(CoeffT1, NoOrbs)
        end if
        
! Put routine in here that takes this rotation matrix, CoeffT1, and forms raises it to the power of a small number, alpha.
! Changeing this number allows us to see the change in plateau level with various rotations.

! Write out some final results of interest, like values of the constraints, values of new coefficients.
    
        write(6,*) 'The final transformation coefficients after gram schmidt orthonormalisation'
        do i = 1, NoOrbs
            do a = 1, NoOrbs
                write(6,'(F10.4)', advance='no') CoeffT1(a,i)
            end do
            write(6,*) ''
        end do

        call WriteTransformMat()
        
        call CalcConstraints(CoeffT1, GSConstraint, TotGSConstraints)  

        write(6,*) 'Final Potential Energy before orthogonalisation', PotEnergy

        call Transform2ElInts()

! Use these final coefficients to find the FourIndInts(i,j,k,l).
! These are now the <ij|kl> integrals we now want to use instead of the HF UMat.
! New potential energy is calculated in this routine using the orthogonalised coefficients.
! Compare to that before this, to make sure the orthogonalisation hasn't shifted them back to a non-minimal place.

        write(6,*) 'Final Potential Energy after orthogonalisation', PotEnergy

! Calculate the fock matrix, and print it out to see how much the off diagonal terms contribute.
! Also print out the sum of the diagonal elements to compare to the original value.
        call CalcFOCKMatrix()

        call RefillUMATandTMAT2D()        
! UMat is the 4 index integral matrix (2 electron), whereas TMAT2D is the 2 index integral (1 el) matrix
   
! This is the keyword that tells the NECI calculation that the orbitals are not HF.  It means that contributions to
! the energy from walkers on singly occupied determinants are included in the values printed.
! Making it true here allows us to go directly from a Rotation into a spawn if required.
        tRotatedOrbs = .true.

        call GENSymStatePairs(SpatOrbs,.false.)
        
    end subroutine FinalizeNewOrbs

    subroutine WriteSingHisttofile()

        integer :: i, j, k, BinNo, a, b, iunit
        real(dp) :: MaxFII, MinFII, BinIter, BinVal, SingExcit(NoOrbs, NoOrbs)

        ! <ik|jk> terms where all i, j and k are virtual
        ! Coulomb.
        ROHistSCijkVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(NoOcc+1, NoOcc+1, NoOcc+2, NoOcc+1)
        MinFII = FourIndInts(NoOcc+1, NoOcc+1, NoOcc+2, NoOcc+1)
        do i = NoOcc+1, NoOrbs
            do k = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,j,k) > MaxFII) MaxFII = FourIndInts(i,k,j,k)
                    if (FourIndInts(i,k,j,k) < MinFII) MinFII = FourIndInts(i,k,j,k)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSCijkVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = NoOcc+1, NoOrbs
            do k = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,j,k) /= 0.0_dp) then
                        BinNo = CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCijkVir(2, BinNo) = ROHistSCijkVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        ! Exchange.
        ROHistSEijkVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(NoOcc+1, NoOcc+1, NoOcc+1, NoOcc+2)
        MinFII = FourIndInts(NoOcc+1, NoOcc+1, NoOcc+1, NoOcc+2)
        do i = NoOcc+1, NoOrbs
            do k = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,k,j) > MaxFII) MaxFII = FourIndInts(i,k,k,j)
                    if (FourIndInts(i,k,k,j) < MinFII) MinFII = FourIndInts(i,k,k,j)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSEijkVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = NoOcc+1, NoOrbs
            do k = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,k,j) /= 0.0_dp) then
                        BinNo = CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEijkVir(2, BinNo) = ROHistSEijkVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        ! Antisymmetric.
        ROHistSASijkVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(NoOcc+1, NoOcc+1, NoOcc+2, NoOcc+1)-FourIndInts(NoOcc+1, NoOcc+1, NoOcc+1, NoOcc+2)
        MinFII = FourIndInts(NoOcc+1, NoOcc+1, NoOcc+2, NoOcc+1)-FourIndInts(NoOcc+1, NoOcc+1, NoOcc+1, NoOcc+2)
        do i = NoOcc+1, NoOrbs
            do k = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) > MaxFII) MaxFII = FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) < MinFII) MinFII = FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSASijkVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = NoOcc+1, NoOrbs
            do k = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) /= 0.0_dp) then
                        BinNo = CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASijkVir(2, BinNo) = ROHistSASijkVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        if (Iteration == 0) then
            iunit  =  get_free_unit()
            open(iunit, file='HistHFSingijkVir', status='unknown')
            do j = 1, 4002
                if ((ROHistSCijkVir(2,j) /= 0).or.(ROHistSEijkVir(2,j) /= 0).or.(ROHistSASijkVir(2,j) /= 0)) then
                    write(iunit,'(6F20.10)') ROHistSCijkVir(1,j), ROHistSCijkVir(2,j), ROHistSEijkVir(1,j), ROHistSEijkVir(2,j),&
                                                        &ROHistSASijkVir(1,j), ROHistSASijkVir(2,j)
                end if
            end do
            close(iunit)
        end if
        if ((Iteration > 1) .and. (.not.tNotConverged)) then
            iunit  =  get_free_unit()
            open(iunit, file='HistRotSingijkVir', status='unknown')
            do j = 1, 4002
                if ((ROHistSCijkVir(2,j) /= 0).or.(ROHistSEijkVir(2,j) /= 0).or.(ROHistSASijkVir(2,j) /= 0)) then
                    write(iunit,'(6F20.10)') ROHistSCijkVir(1,j), ROHistSCijkVir(2,j), ROHistSEijkVir(1,j), ROHistSEijkVir(2,j),&
                                                        &ROHistSASijkVir(1,j), ROHistSASijkVir(2,j)
                end if
            end do
            close(iunit)
        end if

        ! <ik|jk> where k is occupied, and i and j are both virtual.
        ! Coulomb.
        ROHistSCkOcijVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(NoOcc+1, 1, NoOcc+2, 1)
        MinFII = FourIndInts(NoOcc+1, 1, NoOcc+2, 1)
        do i = NoOcc+1, NoOrbs
            do k = 1, NoOcc
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,j,k) > MaxFII) MaxFII = FourIndInts(i,k,j,k)
                    if (FourIndInts(i,k,j,k) < MinFII) MinFII = FourIndInts(i,k,j,k)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSCkOcijVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = NoOcc+1, NoOrbs
            do k = 1, NoOcc
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,j,k) /= 0.0_dp) then
                        BinNo = CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCkOcijVir(2, BinNo) = ROHistSCkOcijVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        ! Exchange.
        ROHistSEkOcijVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(NoOcc+1, 1, 1, NoOcc+2)
        MinFII = FourIndInts(NoOcc+1, 1, 1, NoOcc+2)
        do i = NoOcc+1, NoOrbs
            do k = 1, NoOcc
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,k,j) > MaxFII) MaxFII = FourIndInts(i,k,k,j)
                    if (FourIndInts(i,k,k,j) < MinFII) MinFII = FourIndInts(i,k,k,j)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSEkOcijVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = NoOcc+1, NoOrbs
            do k = 1, NoOcc
                do j = i+1, NoOrbs
                    if (FourIndInts(i,k,k,j) /= 0.0_dp) then
                        BinNo = CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEkOcijVir(2, BinNo) = ROHistSEkOcijVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        !antisymmetric 
        ROHistSASkOcijVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(NoOcc+1, 1, NoOcc+2, 1)-FourIndInts(NoOcc+1, 1, 1, NoOcc+2)
        MinFII = FourIndInts(NoOcc+1, 1, NoOcc+2, 1)-FourIndInts(NoOcc+1, 1, 1, NoOcc+2)
        do i = NoOcc+1, NoOrbs
            do k = 1, NoOcc
                do j = i+1, NoOrbs
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) > MaxFII) MaxFII = FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) < MinFII) MinFII = FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSASkOcijVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = NoOcc+1, NoOrbs
            do k = 1, NoOcc
                do j = i+1, NoOrbs
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) /= 0.0_dp) then
                        BinNo = CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASkOcijVir(2, BinNo) = ROHistSASkOcijVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        if (Iteration == 0) then
            iunit  =  get_free_unit()
            open(iunit, file='HistHFSingkOcijVir', status='unknown')
            do j = 1, 4002
                if ((ROHistSCkOcijVir(2,j) /= 0).or.(ROHistSEkOcijVir(2,j) /= 0).or.(ROHistSASkOcijVir(2,j) /= 0)) then
                  write(iunit,'(6F20.10)') ROHistSCkOcijVir(1,j), ROHistSCkOcijVir(2,j), ROHistSEkOcijVir(1,j), &
                                                         ROHistSEkOcijVir(2,j), ROHistSASkOcijVir(1,j), ROHistSASkOcijVir(2,j)
                end if
            end do
            close(iunit)
        end if
        if ((Iteration > 1) .and. (.not.tNotConverged)) then
            iunit  =  get_free_unit()
            open(iunit, file='HistRotSingkOcijVir', status='unknown')
            do j = 1, 4002
                if ((ROHistSCkOcijVir(2,j) /= 0).or.(ROHistSEkOcijVir(2,j) /= 0).or.(ROHistSASkOcijVir(2,j) /= 0)) then
                  write(iunit,'(6F20.10)') ROHistSCkOcijVir(1,j), ROHistSCkOcijVir(2,j), ROHistSEkOcijVir(1,j), &
                                                         ROHistSEkOcijVir(2,j), ROHistSASkOcijVir(1,j), ROHistSASkOcijVir(2,j)
                end if
            end do
            close(iunit)
        end if

        ! <ik|jk> where i and k are both occupied, and j virtual.
        ! Coulomb.
        ROHistSCikOcjVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(1, 1, NoOcc+1, 1)
        MinFII = FourIndInts(1, 1, NoOcc+1, 1)
        do i = 1, NoOcc
            do k = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if (FourIndInts(i,k,j,k) > MaxFII) MaxFII = FourIndInts(i,k,j,k)
                    if (FourIndInts(i,k,j,k) < MinFII) MinFII = FourIndInts(i,k,j,k)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSCikOcjVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = 1, NoOcc
            do k = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if (FourIndInts(i,k,j,k) /= 0.0_dp) then
                        BinNo = CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCikOcjVir(2, BinNo) = ROHistSCikOcjVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        ! Exchange. 
        ROHistSEikOcjVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(1, 1, 1, NoOcc+1)
        MinFII = FourIndInts(1, 1, 1, NoOcc+1)
        do i = 1, NoOcc
            do k = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if (FourIndInts(i,k,k,j) > MaxFII) MaxFII = FourIndInts(i,k,k,j)
                    if (FourIndInts(i,k,k,j) < MinFII) MinFII = FourIndInts(i,k,k,j)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSEikOcjVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = 1, NoOcc
            do k = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if (FourIndInts(i,k,k,j) /= 0.0_dp) then
                        BinNo = CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEikOcjVir(2, BinNo) = ROHistSEikOcjVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        ! Antisymmetrised.
        ROHistSASikOcjVir(:,:) = 0.0_dp
        MaxFII = FourIndInts(1, 1, NoOcc+1, 1)-FourIndInts(1, 1, 1, NoOcc+1)
        MinFII = FourIndInts(1, 1, NoOcc+1, 1)-FourIndInts(1, 1, 1, NoOcc+1)
        do i = 1, NoOcc
            do k = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) > MaxFII) MaxFII = FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) < MinFII) MinFII = FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                end do
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSASikOcjVir(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do i = 1, NoOcc
            do k = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if ((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)) /= 0.0_dp) then
                        BinNo = CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASikOcjVir(2, BinNo) = ROHistSASikOcjVir(2, BinNo)+1.0         
                    end if
                end do
            end do
        end do

        if (Iteration == 0) then
            iunit  =  get_free_unit()
            open(iunit, file='HistHFSingikOcjVir', status='unknown')
            do j = 1, 4002
                if ((ROHistSCikOcjVir(2,j) /= 0).or.(ROHistSEikOcjVir(2,j) /= 0).or.(ROHistSASikOcjVir(2,j) /= 0)) then 
                  write(iunit,'(6F20.10)') ROHistSCikOcjVir(1,j), ROHistSCikOcjVir(2,j), ROHistSEikOcjVir(1,j), &
                                                         ROHistSEikOcjVir(2,j), ROHistSASikOcjVir(1,j), ROHistSASikOcjVir(2,j)
                end if
            end do
            close(iunit)
        end if
        if ((Iteration > 1) .and. (.not.tNotConverged)) then
            iunit  =  get_free_unit()
            open(iunit, file='HistRotSingikOcjVir', status='unknown')
            do j = 1, 4002
                if ((ROHistSCikOcjVir(2,j) /= 0).or.(ROHistSEikOcjVir(2,j) /= 0).or.(ROHistSASikOcjVir(2,j) /= 0)) then 
                  write(iunit,'(6F20.10)') ROHistSCikOcjVir(1,j), ROHistSCikOcjVir(2,j), ROHistSEikOcjVir(1,j), &
                                                         ROHistSEikOcjVir(2,j), ROHistSASikOcjVir(1,j), ROHistSASikOcjVir(2,j)
                end if
            end do
            close(iunit)
        end if

        ! Single excitations connected to the HF determinant.
        ROHistSing(:,:) = 0.0_dp
        MaxFII = 0.0_dp
        MinFII = 0.0_dp
        do j = NoOcc+1, NoOrbs
            do i = 1, NoOcc
                SingExcit(i,j) = 0.0_dp
                if (i == j) cycle
                a = SymLabelList2_rot(i)
                b = SymLabelList2_rot(j)
                do k = 1, NoOcc+1
                    SingExcit(i,j) = SingExcit(i,j)+real(TMAT2D(2*a, 2*b), dp)+((2*FourIndInts(i,k,j,k))-FourIndInts(i,k,k,j))
                end do
                if (SingExcit(i,j) > MaxFII) MaxFII = SingExcit(i,j)
                if (SingExcit(i,j) < MinFII) MinFII = SingExcit(i,j)
            end do
        end do
        BinIter = ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII = MaxFII+BinIter
        MinFII = MinFII-BinIter
        BinVal = MinFII
        do i = 1, 4002
            ROHistSing(1,i) = BinVal
            BinVal = BinVal+BinIter
        end do
        do j = NoOcc+1, NoOrbs
            do i = 1, NoOcc
                if (i == j) cycle
                if (SingExcit(i,j) /= 0.0_dp) then
                    BinNo = CEILING((SingExcit(i,j)-MinFII)*4002/(MaxFII-MinFII))
                    ROHistSing(2, BinNo) = ROHistSing(2, BinNo)+1.0         
                end if
            end do
        end do

        if (Iteration == 0) then
            iunit  =  get_free_unit()
            open(iunit, file='HistHFSingExcHF', status='unknown')
            do j = 1, 4002
                if (ROHistSing(2,j) /= 0) then
                    do i = 1, 2
                        write(iunit,'(F20.10)', advance='no') ROHistSing(i,j)
                    end do
                    write(iunit,*) ''
                end if
            end do
            close(iunit)
        end if
        if ((Iteration > 1) .and. (.not.tNotConverged)) then
            iunit  =  get_free_unit()
            open(iunit, file='HistRotSingExcHF', status='unknown')
            do j = 1, 4002
                if (ROHistSing(2,j) /= 0) then
                    do i = 1, 2
                        write(iunit,'(F20.10)', advance='no') ROHistSing(i,j)
                    end do
                    write(iunit,*) ''
                end if
            end do
            close(iunit)
        end if

    end subroutine WriteSingHisttofile 

    subroutine WriteDoubHisttofile()

        integer :: i, j, k, l, BinNo, iunit
        real(dp) :: MaxFII, MinFII, BinIter, OnePartOrbEnValue, BinVal

        ! Histogramming all coulomb terms <ij|ij> where i<j, and i and j are
        ! both virtual. In reality we are looking at i = <j, but the ERhistograms
        ! will show the i = j terms.
        if (tROHistVirtCoulomb) then
        
            ROHistDCijOcklVir(:,:) = 0.0_dp
            MinFII = FourIndInts(1, 2, NoOcc+1, NoOcc+2)
            MaxFII = FourIndInts(1, 2, NoOcc+1, NoOcc+2)
            do i = 1, NoOcc
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do l = NoOcc+1, NoOrbs
                            if (FourIndInts(i,j,k,l) < MinFII) MinFII = FourIndInts(i,j,k,l)
                            if (FourIndInts(i,j,k,l) > MaxFII) MaxFII = FourIndInts(i,j,k,l)
                        end do
                    end do
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistDCijOcklVir(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = 1, NoOcc
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do l = NoOcc+1, NoOrbs
                            if (FourIndInts(i,j,k,l) /= 0) then
                                BinNo = CEILING((FourIndInts(i,j,k,l)-MinFII)*4002/(MaxFII-MinFII))
                                ROHistDCijOcklVir(2, BinNo) = ROHistDCijOcklVir(2, BinNo)+1.0         
                            end if
                        end do
                    end do
                end do
            end do

            ! Antisymmetric.
            ROHistASijOcklVir(:,:) = 0.0_dp
            MinFII = FourIndInts(1, 2, NoOcc+1, NoOcc+2)-FourIndInts(1, 2, NoOcc+2, NoOcc+1)
            MaxFII = FourIndInts(1, 2, NoOcc+1, NoOcc+2)-FourIndInts(1, 2, NoOcc+2, NoOcc+1)
            do i = 1, NoOcc
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do l = NoOcc+1, NoOrbs
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) < MinFII) MinFII = (FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) > MaxFII) MaxFII = (FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                        end do
                    end do
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistASijOcklVir(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = 1, NoOcc
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do l = NoOcc+1, NoOrbs
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) /= 0) then
                                BinNo = CEILING(((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))-MinFII)*4002/(MaxFII-MinFII))
                                ROHistASijOcklVir(2, BinNo) = ROHistASijOcklVir(2, BinNo)+1.0         
                            end if
                        end do
                    end do
                end do
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFDoubijOcklVir', status='unknown')
                do j = 1, 4002
                    if ((ROHistDCijOcklVir(2,j) /= 0).or.(ROHistASijOcklVir(2,j) /= 0)) then
                        write(iunit,'(4F20.10)') ROHistDCijOcklVir(1,j), ROHistDCijOcklVir(2,j), &
                            ROHistASijOcklVir(1,j), ROHistASijOcklVir(2,j)
                    end if
                end do
                close(iunit)
            end if
            if ((.not.tNotConverged) .and. (Iteration > 1)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotDoubijOcklVir', status='unknown')
                do j = 1, 4002
                    if ((ROHistDCijOcklVir(2,j) /= 0).or.(ROHistASijOcklVir(2,j) /= 0)) then
                        write(iunit,'(4F20.10)') ROHistDCijOcklVir(1,j), ROHistDCijOcklVir(2,j), &
                            ROHistASijOcklVir(1,j), ROHistASijOcklVir(2,j)
                    end if
                end do
                close(iunit)
            end if

            ROHistDCijklVir(:,:) = 0.0_dp
            MinFII = FourIndInts(NoOrbs-1, NoOrbs, NoOrbs-1, NoOrbs)
            MaxFII = FourIndInts(NoOrbs-1, NoOrbs, NoOrbs-1, NoOrbs)
            do i = NoOcc+1, NoOrbs
                do j = NoOcc+1, NoOrbs
                    do k = i+1, NoOrbs
                        do l = j+1, NoOrbs
                            if (FourIndInts(i,j,k,l) < MinFII) MinFII = FourIndInts(i,j,k,l)
                            if (FourIndInts(i,j,k,l) > MaxFII) MaxFII = FourIndInts(i,j,k,l)
                        end do
                    end do
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistDCijklVir(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = NoOcc+1, NoOrbs
                do j = NoOcc+1, NoOrbs
                    do k = i+1, NoOrbs
                        do l = j+1, NoOrbs
                            if (FourIndInts(i,j,k,l) /= 0) then
                                BinNo = CEILING((FourIndInts(i,j,k,l)-MinFII)*4002/(MaxFII-MinFII))
                                ROHistDCijklVir(2, BinNo) = ROHistDCijklVir(2, BinNo)+1.0         
                            end if
                        end do
                    end do
                end do
            end do

            ! Antisymmetric.
            ROHistASijklVir(:,:) = 0.0_dp
            MinFII = FourIndInts(NoOrbs-3, NoOrbs-2, NoOrbs-1, NoOrbs)-FourIndInts(NoOrbs-3, NoOrbs-2, NoOrbs, NoOrbs-1)
            MaxFII = FourIndInts(NoOrbs-3, NoOrbs-2, NoOrbs-1, NoOrbs)-FourIndInts(NoOrbs-3, NoOrbs-2, NoOrbs, NoOrbs-1)
            do i = NoOcc+1, NoOrbs
                do j = NoOcc+1, NoOrbs
                    do k = i+1, NoOrbs
                        do l = j+1, NoOrbs
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) < MinFII) MinFII = (FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) > MaxFII) MaxFII = (FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                        end do
                    end do
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistASijklVir(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = NoOcc+1, NoOrbs
                do j = NoOcc+1, NoOrbs
                    do k = i+1, NoOrbs
                        do l = j+1, NoOrbs
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) /= 0) then
                                BinNo = CEILING(((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))-MinFII)*4002/(MaxFII-MinFII))
                                ROHistASijklVir(2, BinNo) = ROHistASijklVir(2, BinNo)+1.0         
                            end if
                        end do
                    end do
                end do
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFDoubijklVirt', status='unknown')
                do j = 1, 4002
                    if ((ROHistDCijklVir(2,j) /= 0).or.(ROHistASijklVir(2,j) /= 0)) then
                        write(iunit,'(4F20.10)') ROHistDCijklVir(1,j), ROHistDCijklVir(2,j), &
                            ROHistASijklVir(1,j), ROHistASijklVir(2,j)
                    end if
                end do
                close(iunit)
            end if
            if ((.not.tNotConverged) .and. (Iteration > 1)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotDoubijklVirt', status='unknown')
                do j = 1, 4002
                    if ((ROHistDCijklVir(2,j) /= 0).or.(ROHistASijklVir(2,j) /= 0)) then
                        write(iunit,'(4F20.10)') ROHistDCijklVir(1,j), ROHistDCijklVir(2,j), ROHistASijklVir(1,j), &
                            ROHistASijklVir(2,j)
                    end if
                end do
                close(iunit)
            end if
        end if

        ! Histogramming all one particle orbital energies (occupied and
        ! virtual) even though we are not changing occupied.  Would like to
        ! see HOMO-LUMO gap etc.
        if (tROHistOneElInts) then

            ROHistHijVirt(:,:) = 0.0_dp
            MinFII = TMAT2DRot(NoOcc+1, NoOcc+2)
            MaxFII = TMAT2DRot(NoOcc+1, NoOcc+2)
            do i = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if (TMAT2DRot(i,j) < MinFII) MinFII = TMAT2DRot(i,j)
                    if (TMAT2DRot(i,j) > MaxFII) MaxFII = TMAT2DRot(i,j)
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistHijVirt(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = NoOcc+1, NoOrbs
                do j = i+1, NoOrbs
                    if (TMAT2DRot(i,j) /= 0) then
                        BinNo = CEILING((TMAT2DRot(i,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistHijVirt(2, BinNo) = ROHistHijVirt(2, BinNo)+1.0         
                    end if
                end do
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFHijVirt', status='unknown')
                do j = 1, 4002
                    if (ROHistHijVirt(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistHijVirt(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
            if ((.not.tNotConverged) .and. (Iteration > 1)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotHijVirt', status='unknown')
                do j = 1, 4002
                    if (ROHistHijVirt(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistHijVirt(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
 
            ROHistHijOccVirt(:,:) = 0.0_dp
            MinFII = TMAT2DRot(1, NoOcc+1)
            MaxFII = TMAT2DRot(1, NoOcc+1)
            do i = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if (TMAT2DRot(i,j) < MinFII) MinFII = TMAT2DRot(i,j)
                    if (TMAT2DRot(i,j) > MaxFII) MaxFII = TMAT2DRot(i,j)
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistHijOccVirt(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = 1, NoOcc
                do j = NoOcc+1, NoOrbs
                    if (TMAT2DRot(i,j) /= 0) then
                        BinNo = CEILING((TMAT2DRot(i,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistHijOccVirt(2, BinNo) = ROHistHijOccVirt(2, BinNo)+1.0
                    end if
                end do
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFHijOccVirt', status='unknown')
                do j = 1, 4002
                    if (ROHistHijOccVirt(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistHijOccVirt(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
            if ((.not.tNotConverged) .and. (Iteration > 1)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotHijOccVirt', status='unknown')
                do j = 1, 4002
                    if (ROHistHijOccVirt(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistHijOccVirt(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
 
            ROHistHii(:,:) = 0.0_dp
            MinFII = TMAT2DRot(1,1)
            MaxFII = TMAT2DRot(1,1)
            do i = 1, NoOrbs
                if (TMAT2DRot(i,i) < MinFII) MinFII = TMAT2DRot(i,i)
                if (TMAT2DRot(i,i) > MaxFII) MaxFII = TMAT2DRot(i,i)
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistHii(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = 1, NoOrbs
                BinNo = CEILING((TMAT2DRot(i,i)-MinFII)*4002/(MaxFII-MinFII))
                ROHistHii(2, BinNo) = ROHistHii(2, BinNo)+1.0         
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFHii', status='unknown')
                do j = 1, 4002
                    if (ROHistHii(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistHii(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
            if ((.not.tNotConverged) .and. (Iteration > 1)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotHii', status='unknown')
                do j = 1, 4002
                    if (ROHistHii(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistHii(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
        end if
   
        if (tROHistOnePartOrbEn) then
            ROHistOnePartOrbEn(:,:) = 0.0_dp
            MaxFII = 0.0_dp
            MinFII = 0.0_dp
            do i = 1, NoOrbs
                OnePartOrbEnValue = 0.0_dp
                OnePartOrbEnValue = OnePartOrbEnValue+TMAT2DRot(i,i)
                do j = 1, NoOcc
                    OnePartOrbEnValue = OnePartOrbEnValue+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                end do
                if (OnePartOrbEnValue > MaxFII) MaxFII = OnePartOrbEnValue
                if (OnePartOrbEnValue < MinFII) MinFII = OnePartOrbEnValue
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistOnePartOrbEn(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = 1, NoOrbs
                OnePartOrbEnValue = 0.0_dp
                OnePartOrbEnValue = OnePartOrbEnValue+TMAT2DRot(i,i)
                do j = 1, NoOcc
                    OnePartOrbEnValue = OnePartOrbEnValue+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                end do
                BinNo = CEILING((OnePartOrbEnValue-MinFII)*4002/(MaxFII-MinFII))
                ROHistOnePartOrbEn(2, BinNo) = ROHistOnePartOrbEn(2, BinNo)+1.0         
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFOnePartOrbEn', status='unknown')
                do j = 1, 4002
                    if (ROHistOnePartOrbEn(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistOnePartOrbEn(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
            if ((Iteration > 1) .and. (.not.tNotConverged)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotOnePartOrbEn', status='unknown')
                do j = 1, 4002
                    if (ROHistOnePartOrbEn(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistOnePartOrbEn(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
        end if
  
        if (tROHistDoubExc) then
            ROHistDoubExc(:,:) = 0.0_dp
            MaxFII = 0.0_dp
            MinFII = 0.0_dp
            do l = NoOcc+1, NoOrbs
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do i = 1, NoOcc
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) > MaxFII) then
                                MaxFII = (FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            end if
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) < MinFII) then
                                MinFII = (FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            end if
                        end do
                    end do
                end do
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistDoubExc(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do l = NoOcc+1, NoOrbs
                do j = 1, NoOcc
                    do k = NoOcc+1, NoOrbs
                        do i = 1, NoOcc
                            if ((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)) /= 0) then
                                BinNo = CEILING(((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))-MinFII)*4002/(MaxFII-MinFII))
                                ROHistDoubExc(2, BinNo) = ROHistDoubExc(2, BinNo)+1.0         
                            end if
                        end do
                    end do
                end do
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHFDoubExc', status='unknown')
                do j = 1, 4002
                    if (ROHistDoubExc(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistDoubExc(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
            if ((Iteration > 1) .and. (.not.tNotConverged)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRotDoubExc', status='unknown')
                do j = 1, 4002
                    if (ROHistDoubExc(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistDoubExc(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if
        end if

        if (tROHistER) then
            ROHistER(:,:) = 0.0_dp
            MaxFII = 0.0_dp
            MinFII = 0.0_dp
            do i = 1, NoOrbs
                if (FourIndInts(i,i,i,i) > MaxFII) MaxFII = FourIndInts(i,i,i,i)
                if (FourIndInts(i,i,i,i) < MinFII) MinFII = FourIndInts(i,i,i,i)  
            end do
            BinIter = ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII = MaxFII+BinIter
            MinFII = MinFII-BinIter
            BinVal = MinFII
            do i = 1, 4002
                ROHistER(1,i) = BinVal
                BinVal = BinVal+BinIter
            end do
            do i = 1, NoOrbs
                BinNo = CEILING((FourIndInts(i,i,i,i)-MinFII)*4002/(MaxFII-MinFII))
                ROHistER(2, BinNo) = ROHistER(2, BinNo)+1.0         
            end do

            if (Iteration == 0) then 
                iunit  =  get_free_unit()
                open(iunit, file='HistHF-ER', status='unknown')
                do j = 1, 4002
                    if (ROHistER(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistER(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if

            if ((Iteration > 1) .and. (.not.tNotConverged)) then
                iunit  =  get_free_unit()
                open(iunit, file='HistRot-ER', status='unknown')
                do j = 1, 4002
                    if (ROHistER(2,j) /= 0) then
                        do i = 1, 2
                            write(iunit,'(F20.10)', advance='no') ROHistER(i,j)
                        end do
                        write(iunit,*) ''
                    end if
                end do
                close(iunit)
            end if

        end if

    end subroutine WriteDoubHisttofile 

    subroutine PrintIntegrals()

        integer :: i, j, k, l, io1, io2
        real(dp) :: DiagOneElPot, ERPot, ijVirtOneElPot, ijVirtCoulPot, ijVirtExchPot
        real(dp) :: singCoulijVirt, singExchijVirt, singCoulconHF, singExchconHF, ijklPot, ijklantisymPot
        real(dp) :: ijOccVirtOneElPot, ijOccVirtCoulPot, ijOccVirtExchPot

        io1 = 0
        io2 = 0
        if (tInitIntValues) then
            io1 = get_free_unit()
            open(io1, file='DiagIntegrals', status='unknown')
            write(io1,'(A10, 6A18)') "Iteration","<i|h|i> ivirt","<ii|ii> ivirt","<ij|ij> iOccjVirt","<ij|ji> iOccjVirt", &
                "<ij|ij> ijVirt","<ij|ji> ijVirt"

            io2 = get_free_unit()
            open(io2, file='SingExcIntegrals', status='unknown')
            write(io2,'(A10, 6A18)') "Iteration","<i|h|j> iOccjVirt","<i|h|j> ijVirt","<ik|jk> HFcon","<ik|kj> HFcon", &
                "<ik|jk> ijVirt","<ik|kj> ijVirt"

            DiagOneElPotInit = 0.0_dp
            ERPotInit = 0.0_dp
            ijVirtOneElPotInit = 0.0_dp
            ijVirtCoulPotInit = 0.0_dp
            ijVirtExchPotInit = 0.0_dp
            singCoulconHFInit = 0.0_dp
            singExchconHFInit = 0.0_dp
            singCoulijVirtInit = 0.0_dp
            singExchijVirtInit = 0.0_dp
            ijklPotInit = 0.0_dp
            ijklantisymPotInit = 0.0_dp
            ijOccVirtOneElPotInit = 0.0_dp
            ijOccVirtCoulPotInit = 0.0_dp
            ijOccVirtExchPotInit = 0.0_dp
            NoInts01 = 0
            NoInts02 = 0
            NoInts03 = 0
            NoInts04 = 0
            NoInts05 = 0
            NoInts06 = 0
            do i = 1, NoOrbs
                if (i > NoOcc) then
                    DiagOneElPotInit = DiagOneElPotInit+TMAT2DRot(i,i)
                    ERPotInit = ERPotInit+FourIndInts(i,i,i,i)
                    NoInts01 = NoInts01+1
                    do j = NoOcc+1, NoOrbs
                       ! The i, j terms with i and j both virtual.
                       if (j > i) then
                           ijVirtOneElPotInit = ijVirtOneElPotInit+TMAT2DRot(i,j)
                           ijVirtCoulPotInit = ijVirtCoulPotInit+FourIndInts(i,j,i,j)
                           ijVirtExchPotInit = ijVirtExchPotInit+FourIndInts(i,j,j,i)
                           NoInts02 = NoInts02+1
                       end if
                       do k = 1, NoOrbs
                           if (k > (NoOcc+1)) then
                               do l = NoOcc+1, NoOrbs
                                   if (l == j) cycle
                                   ijklPotInit = ijklPotInit+FourIndInts(i,j,k,l)
                                   ijklantisymPotInit = ijklantisymPotInit+FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)
                                   NoInts04 = NoInts04+1
                               end do
                            else
                                if (i == j) cycle
                                singCoulijVirtInit = singCoulijVirtInit+FourIndInts(i,k,j,k)
                                singExchijVirtInit = singExchijVirtInit+FourIndInts(i,k,k,j)
                                NoInts03 = NoInts03+1
                            end if
                       end do
                   end do
               else
                   do j = NoOcc+1, NoOrbs
                       do k = 1, NoOcc
                           singCoulconHFInit = singCoulconHFInit+FourIndInts(i,k,j,k)
                           singExchconHFInit = singExchconHFInit+FourIndInts(i,k,k,j)
                           NoInts06 = NoInts06+1
                       end do
                       ijOccVirtOneElPotInit = ijOccVirtOneElPotInit+TMAT2DRot(i,j)
                       ijOccVirtCoulPotInit = ijOccVirtCoulPotInit+FourIndInts(i,j,i,j)
                       ijOccVirtExchPotInit = ijOccVirtExchPotInit+FourIndInts(i,j,j,i)
                       NoInts05 = NoInts05+1
                   end do
               end if
            end do
        end if
        
        DiagOneElPot = 0.0_dp
        ERPot = 0.0_dp
        ijVirtOneElPot = 0.0_dp
        ijVirtCoulPot = 0.0_dp
        ijVirtExchPot = 0.0_dp
        singCoulconHF = 0.0_dp
        singExchconHF = 0.0_dp
        singCoulijVirt = 0.0_dp
        singExchijVirt = 0.0_dp
        ijklPot = 0.0_dp
        ijklantisymPot = 0.0_dp
        ijOccVirtOneElPot = 0.0_dp
        ijOccVirtCoulPot = 0.0_dp
        ijOccVirtExchPot = 0.0_dp
        do i = 1, NoOrbs
            if (i > NoOcc) then
                DiagOneElPot = DiagOneElPot+TMAT2DRot(i,i)
                ERPot = ERPot+FourIndInts(i,i,i,i)
                do j = NoOcc+1, NoOrbs
                   ! The i, j terms with i and j both virtual.
                   if (j > i) then
                       ijVirtOneElPot = ijVirtOneElPot+TMAT2DRot(i,j)
                       ijVirtCoulPot = ijVirtCoulPot+FourIndInts(i,j,i,j)
                       ijVirtExchPot = ijVirtExchPot+FourIndInts(i,j,j,i)
                   end if
                   do k = 1, NoOrbs
                       if (k > (NoOcc+1)) then
                           do l = NoOcc+1, NoOrbs
                               if (l == j) cycle
                               ijklPot = ijklPot+FourIndInts(i,j,k,l)
                               ijklantisymPot = ijklantisymPot+FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)
                           end do
                       else
                           if (i == j) cycle
                           singCoulijVirt = singCoulijVirt+FourIndInts(i,k,j,k)
                           singExchijVirt = singExchijVirt+FourIndInts(i,k,k,j)
                       end if
                   end do
               end do
           else
               do j = NoOcc+1, NoOrbs
                   do k = 1, NoOcc
                       singCoulconHF = singCoulconHF+FourIndInts(i,k,j,k)
                       singExchconHF = singExchconHF+FourIndInts(i,k,k,j)
                   end do
                   ijOccVirtOneElPot = ijOccVirtOneElPot+TMAT2DRot(i,j)
                   ijOccVirtCoulPot = ijOccVirtCoulPot+FourIndInts(i,j,i,j)
                   ijOccVirtExchPot = ijOccVirtExchPot+FourIndInts(i,j,j,i)
               end do
           end if
        end do

        DiagOneElPot = (DiagOneElPot-DiagOneElPotInit)/NoInts01
        ERPot = (ERPot-ERPotInit)/NoInts01
        ijVirtOneElPot = (ijVirtOneElPot-ijVirtOneElPotInit)/NoInts02
        ijVirtCoulPot = (ijVirtCoulPot-ijVirtCoulPotInit)/NoInts02
        ijVirtExchPot = (ijVirtExchPot-ijVirtExchPotInit)/NoInts02
        singCoulijVirt = (singCoulijVirt-singCoulijVirtInit)/NoInts03
        singExchijVirt = (singExchijVirt-singExchijVirtInit)/NoInts03
        singCoulconHF = (singCoulconHF-singCoulconHFInit)/NoInts06
        singExchconHF = (singExchconHF-singExchconHFInit)/NoInts06
        ijklPot = (ijklPot-ijklPotInit)/NoInts04
        ijklantisymPot = (ijklantisymPot-ijklantisymPotInit)/NoInts04
        ijOccVirtOneElPot = (ijOccVirtOneElPot-ijOccVirtOneElPotInit)/NoInts05
        ijOccVirtCoulPot = (ijOccVirtCoulPot-ijOccVirtCoulPotInit)/NoInts05
        ijOccVirtExchPot = (ijOccVirtExchPot-ijOccVirtExchPot)/NoInts05


        write(io1,'(I10, 6F18.10)') Iteration, DiagOneElPot, ERPot, ijOccVirtCoulPot, ijOccVirtExchPot, ijVirtCoulPot, &
            ijVirtExchPot
        write(io2,'(I10, 6F18.10)') Iteration, ijOccVirtOneElPot, ijVirtOneElPot, singCoulconHF, singExchconHF, singCoulijVirt, &
            singExchijVirt

        if ((.not.tNotConverged) .and. (.not.tInitIntValues)) then
            close(io1)
            close(io2)
        end if

    end subroutine PrintIntegrals

    subroutine CalcFOCKMatrix()

        use SystemData, only: nBasis
        use LoggingData, only: tRDMonfly

        integer :: i, j, k, l, a, b, ierr
        real(dp) :: FOCKDiagSumHF, FOCKDiagSumNew
        character(len=*), parameter :: this_routine = 'CalcFOCKMatrix'
        !NEED TO FIX THIS!

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

! When transforming the orbitals into approximate natural orbitals, we want to save memory, so don't bother
! calculating the whole matrix, just the diagonal elements that we actually need.

        if (tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) then
            allocate(ArrDiagNew(NoOrbs), stat = ierr)
            call LogMemAlloc('ArrDiagNew', NoOrbs, 8, this_routine, ArrDiagNewTag, ierr)
            ArrDiagNew(:) = 0.0_dp                     
        else
            allocate(ArrNew(NoOrbs, NoOrbs), stat = ierr)
            call LogMemAlloc('ArrNew', NoOrbs**2, 8, this_routine, ArrNewTag, ierr)
            ArrNew(:,:) = 0.0_dp                     
        end if

        ! First calculate the sum of the diagonal elements, ARR.
        ! Check if this is already being done.
        FOCKDiagSumHF = 0.0_dp
        do a = 1, nBasis        
            FOCKDiagSumHF = FOCKDiagSumHF+Arr(a,2)
        end do

        write(6,*) 'Sum of the fock matrix diagonal elements in the HF basis set = ', FOCKDiagSumHF

        FOCKDiagSumNew = 0.0_dp
        do j = 1, NoRotOrbs
            l = SymLabelList3_rot(j)
            if (tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) then
                do a = 1, NoOrbs
                    b = SymLabelList2_rot(a)
                    if (tStoreSpinOrbs.or.tTurnStoreSpinOff) then
                        ArrDiagNew(l) = ArrDiagNew(l)+(CoeffT1(a,j)*ARR(b,2)*CoeffT1(a,j))
                    else
                        ArrDiagNew(l) = ArrDiagNew(l)+(CoeffT1(a,j)*ARR(2*b, 2)*CoeffT1(a,j))
                    end if
                end do
                if (tStoreSpinOrbs.or.tTurnStoreSpinOff) then
                    FOCKDiagSumNew = FOCKDiagSumNew+(ArrDiagNew(l))
                else
                    FOCKDiagSumNew = FOCKDiagSumNew+(ArrDiagNew(l)*2)
                end if
            else
                do i = 1, NoRotOrbs
                    k = SymLabelList2_rot(i)
                    ArrNew(k,l) = 0.0_dp
                    do a = 1, NoOrbs
                        b = SymLabelList2_rot(a)
                        if (tStoreSpinOrbs.or.tTurnStoreSpinOff) then
                            ArrNew(k,l) = ArrNew(k,l)+(CoeffT1(a,i)*Arr(b,2)*CoeffT1(a,j))
                        else
                            ArrNew(k,l) = ArrNew(k,l)+(CoeffT1(a,i)*Arr(2*b, 2)*CoeffT1(a,j))
                        end if
                    end do
                end do
                if (tStoreSpinOrbs.or.tTurnStoreSpinOff) then
                    FOCKDiagSumNew = FOCKDiagSumNew+(ArrNew(l,l))
                else
                    FOCKDiagSumNew = FOCKDiagSumNew+(ArrNew(l,l)*2)
                end if
                ! Only running through spat orbitals, count each twice to
                ! compare to above.
            end if
        end do
        ! If we are truncation the virtual space, only the unfrozen entries
        ! will be transformed.
        
        write(6,*) 'Sum of the fock matrix diagonal elements in the transformed basis set = ', FOCKDiagSumNew

! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2) (ordered in terms of orbital number).
! ARR(:,2) needs to be ordered in terms of symmetry and then energy (like SymLabelList), so currently this ordering will not be 
! correct when reading in qchem intDUMPS as the orbital number ordering is by energy.

        ! If we are only writing out 1 ROFCIDUMP or we are not truncating at
        ! all - can refill ARR etc.
        if (NoDumpTruncs <= 1) then

            if (tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) then
                if (tStoreSpinOrbs.or.tTurnStoreSpinOff) then
                    do j = 1, NoOrbs
                        ARR(j,2) = ArrDiagNew(j)
                        ARR(j,1) = ArrDiagNew(BRR(j))
                    end do
                else
                    do j = 1, NoOrbs
                        ARR(2*j, 2) = ArrDiagNew(j)
                        ARR(2*j-1, 2) = ArrDiagNew(j)
                        ARR(2*j, 1) = ArrDiagNew(BRR(2*j)/2)
                        ARR(2*j-1, 1) = ArrDiagNew(BRR(2*j)/2)
                    end do
                end if
            else
                if (tStoreSpinOrbs.or.tTurnStoreSpinOff) then
                    do j = 1, NoRotOrbs
                        ARR(j,2) = ArrNew(j,j)
                        ARR(j,1) = ArrNew(BRR(j), BRR(j))
                    end do
                else
                    do j = 1, NoRotOrbs
                        ARR(2*j, 2) = ArrNew(j,j)
                        ARR(2*j-1, 2) = ArrNew(j,j)
                        ARR(2*j, 1) = ArrNew(BRR(2*j)/2, BRR(2*j)/2)
                        ARR(2*j-1, 1) = ArrNew(BRR(2*j)/2, BRR(2*j)/2)
                    end do
                end if
            end if

        end if

        if ((tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) .and. (NoDumpTruncs <= 1)) then
            deallocate(ArrDiagNew)
            call LogMemDealloc(this_routine, ArrDiagNewTag)
        elseif (NoDumpTruncs <= 1) then
            deallocate(ArrNew)
            call LogMemDealloc(this_routine, ArrNewTag)
        end if

        write(6,*) 'end of calcfockmatrix'
        call neci_flush(6)

    end subroutine CalcFOCKMatrix

    subroutine RefillUMATandTMAT2D()

        integer :: l, k, j, i, a, b, g, d, c, nBasis2, ierr
        integer(TagIntType) :: TMAT2DPartTag
        real(dp) :: NewTMAT
        real(dp), allocatable :: TMAT2DPart(:,:)
#ifdef __CMPLX
        call stop_all('RefillUMATandTMAT2D', 'Rotating orbitals not implemented for complex orbitals.')
#endif

        if (tStoreSpinOrbs) then
            allocate(TMAT2DPart((nBasis-NoFrozenVirt), nBasis), stat = ierr)
            call LogMemAlloc('TMAT2DPart',(nBasis-NoFrozenVirt)*nBasis, 8,'RefillUMAT', TMAT2DPartTag, ierr)
            if (NoDumpTruncs > 1) then
                allocate(TMAT2DNew((nBasis-NoFrozenVirt),(nBasis-NoFrozenVirt)), stat = ierr)
                call LogMemAlloc('TMAT2DNew',(nBasis-NoFrozenVirt)**2, 8,'RefillUMAT', TMAT2DNewTag, ierr)
                TMAT2DNew(:,:) = 0.0_dp
            end if
        else
            allocate(TMAT2DPart((nBasis-(NoFrozenVirt*2)), nBasis), stat = ierr)
            call LogMemAlloc('TMAT2DPart',(nBasis-(NoFrozenVirt*2))*nBasis, 8,'RefillUMAT', TMAT2DPartTag, ierr)
            if (NoDumpTruncs > 1) then
                allocate(TMAT2DNew((nBasis-NoFrozenVirt),(nBasis-NoFrozenVirt)), stat = ierr)
                call LogMemAlloc('TMAT2DNew',(nBasis-NoFrozenVirt)**2, 8,'RefillUMAT', TMAT2DNewTag, ierr)
                TMAT2DNew(:,:) = 0.0_dp
            end if
        end if
        TMAT2DPart(:,:) = 0.0_dp

        RefillUMAT_Time%timer_name = 'RefillUMATandTMAT'
        call set_timer(RefillUMAT_Time, 30)

        do i  =  1, nBasis
            write(6,*) SymLabelList2_rot(i), SymLabelList3_rot(i)
        end do

        ! Make the UMAT elements the four index integrals. These are calculated
        ! by transforming the HF orbitals using the coefficients that have been
        ! found.
        if (NoDumpTruncs <= 1) then
            do l = 1,(NoOrbs-(NoFrozenVirt))
                if (tTurnStoreSpinOff) then
                    d = CEILING(real(SymLabelList3_rot(l), dp)/2.0_dp)
                else
                    d = SymLabelList3_rot(l)
                end if
                do k = 1,(NoOrbs-(NoFrozenVirt))

                    if (tTurnStoreSpinOff) then
                        g = CEILING(real(SymLabelList3_rot(k), dp)/2.0_dp)
                    else
                        g = SymLabelList3_rot(k)
                    end if
     
                    do j = 1,(NoOrbs-(NoFrozenVirt))

                        if (tTurnStoreSpinOff) then
                            b = CEILING(real(SymLabelList3_rot(j), dp)/2.0_dp)
                        else
                            b = SymLabelList3_rot(j)
                        end if
                        do i = 1,(NoOrbs-(NoFrozenVirt))

                            if (tTurnStoreSpinOff) then
                                a = CEILING(real(SymLabelList3_rot(i), dp)/2.0_dp)
                            else
                                a = SymLabelList3_rot(i)
                            end if

                            if (tUseMP2VarDenMat.or.tFindCINatOrbs.or.tReadInCoeff) then
                                UMAT(UMatInd(a, b, g, d, 0, 0)) = (FourIndInts(i,k,j,l))
                            else
                                UMAT(UMatInd(a, b, g, d, 0, 0)) = (FourIndInts(i,j,k,l))
                            end if
                        end do
                    end do
                end do
            end do
        end if

        do a = 1, nBasis
            do k = 1, NoRotOrbs
                i = SymLabelList3_rot(k)
                NewTMAT = 0.0_dp
                do b = 1, NoOrbs
                    d = SymLabelList2_rot(b)
                    if (tStoreSpinOrbs) then
                        NewTMAT = NewTMAT+(CoeffT1(b,k)*real(TMAT2D(d,a), dp))
                    else
                        NewTMAT = NewTMAT+(CoeffT1(b,k)*real(TMAT2D(2*d, a), dp))
                    end if
                end do
                if (tStoreSpinOrbs) then
                    TMAT2DPart(i,a) = NewTMAT
                else
                    TMAT2DPart(2*i, a) = NewTMAT
                    TMAT2DPart(2*i-1, a) = NewTMAT
                end if
            end do
        end do

        if (tStoreSpinOrbs) then
            nBasis2 = nBasis-NoFrozenVirt
        else
            nBasis2 = nBasis-(NoFrozenVirt*2)
        end if
        do k = 1, nBasis2
            do l = 1, NoRotOrbs
                j = SymLabelList3_rot(l)
                NewTMAT = 0.0_dp
                do a = 1, NoOrbs
                    c = SymLabelList2_rot(a)
                    if (tStoreSpinOrbs) then
                        NewTMAT = NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,c))
                    else
                        NewTMAT = NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k, 2*c))
                    end if
                end do
                if (tStoreSpinOrbs) then
                    if (NoDumpTruncs > 1) then
                        TMAT2DNew(k,j) = NewTMAT
                    else
                        TMAT2D(k,j) = (NewTMAT)
                    end if
                else
                    if (NoDumpTruncs > 1) then
                        TMAT2DNew(k, 2*j) = NewTMAT
                        TMAT2DNew(k, 2*j-1) = NewTMAT
                    else
                        TMAT2D(k, 2*j) = (NewTMAT)
                        TMAT2D(k, 2*j-1) = (NewTMAT)
                    end if
                end if
            end do
        end do

        deallocate(TMAT2DPart)
        call LogMemDeAlloc('RefillUMAT', TMAT2DPartTag)


        if (tROHistSingExc) call WriteSingHisttofile()

        call set_timer(RefillUMAT_Time, 30)

        if (tTurnStoreSpinOff) then
            tStoreSpinOrbs  =  .false.
            NoOrbs  =  nBasis / 2
        end if

        write(6,'(A, I5, A)') ' Printing the new ROFCIDUMP file for a truncation of ', NoFrozenVirt,' orbitals.'
        if (tROFciDump .and. (NoDumpTruncs > 1)) then
            call PrintRepeatROFCIDUMP()
        elseif (tROFciDUmp) then
            call PrintROFCIDUMP()
        end if

    end subroutine RefillUMATandTMAT2D

    subroutine PrintROFCIDUMP()

        ! This prints out a new FCIDUMP file in the same format as the old one.

        integer :: i, j, k, l, iunit
        character(len=5) :: Label
        character(len=20) :: LabelFull

        PrintROFCIDUMP_Time%timer_name = 'PrintROFCIDUMP'
        call set_timer(PrintROFCIDUMP_Time, 30)

        Label = ''
        LabelFull = ''
        write(Label,'(I5)') NoFrozenVirt
        LabelFull = 'ROFCIDUMP-'//adjustl(Label)

        iunit  =  get_free_unit()
        open(iunit, file=LabelFull, status='unknown')
        
        write(iunit,'(2A6, I3, A7, I3, A5, I2, A)') '&FCI ','NORB = ',(NoOrbs-(NoFrozenVirt)),', NELEC = ', NEl,', MS2 = ', LMS,','
        write(iunit,'(A9)', advance='no') 'ORBSYM = '
        do i = 1,(NoOrbs-(NoFrozenVirt))
            if ((tUseMP2VarDenMat.or.tFindCINatOrbs) .and. (.not.lNoSymmetry) .and. tTruncRODump) then
                write(iunit,'(I1, A1)', advance='no') (SymOrbs_rot(i)+1),','
            else
                if (tStoreSpinOrbs) then
                    write(iunit,'(I1, A1)', advance='no') (int(G1(i)%sym%S)+1),','
                else
                    write(iunit,'(I1, A1)', advance='no') (int(G1(i*2)%sym%S)+1),','
                end if
            end if
        end do
        write(iunit,*) ''
        if (tStoreSpinOrbs) then
            write(iunit,'(A7, I1, A11)') 'ISYM = ', 1,' UHF = .TRUE.'
        else
            write(iunit,'(A7, I1)') 'ISYM = ', 1
        end if
        write(iunit,'(A5)') '&end'
       
        do i = 1,(NoOrbs-(NoFrozenVirt))
            do k = 1, i
                do j = 1,(NoOrbs-(NoFrozenVirt))
                    ! Potential to put symmetry in here, have currently taken
                    ! it out, because when we're only printing non-zero values,
                    ! it is kind of unnecessary - although it may be used to
                    ! speed things up.
                    do l = 1, j
                        if ((ABS(real(UMat(UMatInd(i, j, k, l, 0, 0)), dp))) /= 0.0_dp) &
                                        &write(iunit,'(F21.12, 4I3)') real(UMat(UMatInd(i, j, k, l, 0, 0)), dp), i, k, j, l 
                    end do
                end do
           end do
        end do

        ! TMAT2D stored as spin orbitals.
        do k = 1,(NoOrbs-(NoFrozenVirt))
            ! Symmetry?
            do i = k,(NoOrbs-(NoFrozenVirt))
                if (tStoreSpinOrbs) then
                    if ((real(TMAT2D(i,k), dp)) /= 0.0_dp) write(iunit,'(F21.12, 4I3)') real(TMAT2D(i,k), dp), i, k, 0, 0
                else
                    if ((real(TMAT2D(2*i, 2*k), dp)) /= 0.0_dp) write(iunit,'(F21.12, 4I3)') real(TMAT2D(2*i, 2*k), dp), i, k, 0, 0
                end if
            end do
        end do

! ARR has the energies of the orbitals (eigenvalues).  ARR(:,2) has ordering we want.
! ARR is stored as spin orbitals.

        do k = 1,(NoOrbs-(NoFrozenVirt))
            if (tStoreSpinOrbs) then
                write(iunit,'(F21.12, 4I3)') Arr(k,2), k, 0, 0, 0
            else
                write(iunit,'(F21.12, 4I3)') Arr(2*k, 2), k, 0, 0, 0
            end if
        end do

        write(iunit,'(F21.12, 4I3)') ECore, 0, 0, 0, 0
        
        call neci_flush(iunit)

        close(iunit)

        call halt_timer(PrintROFCIDUMP_Time)

    end subroutine PrintROFCIDUMP

    subroutine PrintRepeatROFCIDUMP()

        ! This prints out a new FCIDUMP file in the same format as the old one.

        integer :: i, j, k, l, ierr, a, b, g, d, iunit
        character(len=5) :: Label
        character(len=20) :: LabelFull
        character(len=*), parameter :: this_routine = 'PrintRepeatROFCIDUMP'

        PrintROFCIDUMP_Time%timer_name = 'PrintROFCIDUMP'
        call set_timer(PrintROFCIDUMP_Time, 30)

        Label = ''
        LabelFull = ''
        write(Label,'(I5)') NoFrozenVirt
        LabelFull = 'ROFCIDUMP-'//adjustl(Label)

        iunit = get_free_unit()
        open(iunit, file=LabelFull, status='unknown')
        
        write(iunit,'(2A6, I3, A7, I3, A5, I2, A)') '&FCI ','NORB = ',(NoOrbs-(NoFrozenVirt)),', NELEC = ', NEl,', MS2 = ', LMS,','
        write(iunit,'(A9)', advance='no') 'ORBSYM = '
        do i = 1,(NoOrbs-(NoFrozenVirt))
            if ((tUseMP2VarDenMat.or.tFindCINatOrbs) .and. (.not.lNoSymmetry) .and. tTruncRODump) then
                write(iunit,'(I1, A1)', advance='no') (SymOrbs_rot(i)+1),','
            else
                if (tStoreSpinOrbs) then
                    write(iunit,'(I1, A1)', advance='no') (int(G1(i)%sym%S)+1),','
                else
                    write(iunit,'(I1, A1)', advance='no') (int(G1(i*2)%sym%S)+1),','
                end if
            end if
        end do
        write(iunit,*) ''
        if (tStoreSpinOrbs) then
            write(iunit,'(A7, I1, A11)') 'ISYM = ', 1,' UHF = .TRUE.'
        else
            write(iunit,'(A7, I1)') 'ISYM = ', 1
        end if
        write(iunit,'(A5)') '&end'

        allocate(SymLabelList3_rotInv(NoOrbs), stat = ierr)
        call LogMemAlloc('SymLabelList3_rotInv', NoOrbs, 4, this_routine, SymLabelList3_rotInvTag, ierr)
        SymLabelList3_rotInv(:) = 0                     

        do i = 1, NoOrbs
            SymLabelList3_rotInv(SymLabelList3_rot(i)) = i
        end do
       
        do i = 1,(NoOrbs-(NoFrozenVirt))
            a = SymLabelList3_rotInv(i)
            do k = 1, i
                g = SymLabelList3_rotInv(k)
                do j = 1,(NoOrbs-(NoFrozenVirt))
                    b = SymLabelList3_rotInv(j)
                    ! Potential to put symmetry in here, have currently taken
                    ! it out, because when we're only printing non-zero values,
                    ! it is kind of unnecessary - although it may be used to
                    ! speed things up.
                    do l = 1, j
                        d = SymLabelList3_rotInv(l)
                        if ((ABS(FourIndInts(a,g,b,d))) /= 0.0_dp) &
                                        &write(iunit,'(F21.12, 4I3)') FourIndInts(a,g,b,d), i, k, j, l 
 
                    end do
                end do
           end do
        end do

        deallocate(SymLabelList3_rotInv)
        call LogMemDeAlloc(this_routine, SymLabelList3_rotInvTag)

        ! TMAT2D stored as spin orbitals.
        do k = 1,(NoOrbs-(NoFrozenVirt))
            ! Symmetry?
            do i = k,(NoOrbs-(NoFrozenVirt))
                if (tStoreSpinOrbs) then
                    if (TMAT2DNew(i,k) /= 0.0_dp) write(iunit,'(F21.12, 4I3)') TMAT2DNew(i,k), i, k, 0, 0
                else
                    if (TMAT2DNew(2*i, 2*k) /= 0.0_dp) write(iunit,'(F21.12, 4I3)') TMAT2DNew(2*i, 2*k), i, k, 0, 0
                end if
            end do
        end do

        ! ARR has the energies of the orbitals (eigenvalues). ARR(:,2) has
        ! ordering we want. ARR is stored as spin orbitals.

        if (tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs) then
            if (tStoreSpinOrbs) then
                do k = 1,(NoOrbs-(NoFrozenVirt))
                    write(iunit,'(F21.12, 4I3)') ArrDiagNew(k), k, 0, 0, 0
                end do
            else

                do k = 1,(NoOrbs-(NoFrozenVirt))

                    write(iunit,'(F21.12, 4I3)') ArrDiagNew(k), k, 0, 0, 0
                end do
            end if
        else
            if (tStoreSpinOrbs) then
                do k = 1,(NoOrbs-(NoFrozenVirt))
                    write(iunit,'(F21.12, 4I3)') ArrNew(k,k), k, 0, 0, 0
                end do
            else
                do k = 1,(NoOrbs-(NoFrozenVirt))
                    write(iunit,'(F21.12, 4I3)') ArrNew(k,k), k, 0, 0, 0
                end do
            end if
        end if

        write(iunit,'(F21.12, 4I3)') ECore, 0, 0, 0, 0
        
        call neci_flush(iunit)

        close(iunit)

        call halt_timer(PrintROFCIDUMP_Time)

    end subroutine PrintRepeatROFCIDUMP

    subroutine DeallocateMem()

        character(len=*), parameter :: this_routine = 'DeallocateMem'

        deallocate(Lab)
        call LogMemDealloc(this_routine, LabTag)
        deallocate(CoeffT1)
        call LogMemDealloc(this_routine, CoeffT1Tag)
        deallocate(CoeffCorT2)
        call LogMemDealloc(this_routine, CoeffCorT2Tag)
        deallocate(CoeffUncorT2)
        call LogMemDealloc(this_routine, CoeffUncorT2Tag)
        if (tLagrange) then
            deallocate(Lambdas)
            call LogMemDealloc(this_routine, LambdasTag)
            deallocate(DerivLambda)
            call LogMemDealloc(this_routine, DerivLambdaTag)
        end if 
        deallocate(DerivCoeff)
        call LogMemDealloc(this_routine, DerivCoeffTag)

        deallocate(DiagTMAT2Dfull)
        call LogMemDealloc(this_routine, DiagTMAT2DfullTag)
 
        deallocate(TwoIndInts01)
        call LogMemDealloc(this_routine, TwoIndInts01Tag)
        deallocate(ThreeIndInts02)
        call LogMemDealloc(this_routine, ThreeIndInts02Tag)
        deallocate(FourIndInts)
        call LogMemDealloc(this_routine, FourIndIntsTag)
        deallocate(FourIndInts02)
        call LogMemDealloc(this_routine, FourIndInts02Tag)

        if (tERLocalization .and. (.not.tStoreSpinOrbs)) then
            deallocate(TwoIndIntsER)
            call LogMemDeAlloc(this_routine, TwoIndIntsERTag)
            deallocate(ThreeIndInts01ER)
            call LogMemDeAlloc(this_routine, ThreeIndInts01ERTag)
            deallocate(ThreeIndInts02ER)
            call LogMemDeAlloc(this_routine, ThreeIndInts02ERTag)
            deallocate(FourIndIntsER)
            call LogMemDeAlloc(this_routine, FourIndIntsERTag)
        else
            deallocate(TMAT2DTemp)
            call LogMemDealloc(this_routine, TMAT2DTempTag)
            deallocate(TMAT2DPartRot01)
            call LogMemDealloc(this_routine, TMAT2DPartRot01Tag)
            deallocate(TMAT2DPartRot02)
            call LogMemDealloc(this_routine, TMAT2DPartRot02Tag)
            deallocate(TMAT2DRot)
            call LogMemDealloc(this_routine, TMAT2DRotTag)
     
            deallocate(TwoIndInts02)
            call LogMemDealloc(this_routine, TwoIndInts02Tag)
            deallocate(ThreeIndInts01)
            call LogMemDealloc(this_routine, ThreeIndInts01Tag)
            deallocate(ThreeIndInts03)
            call LogMemDealloc(this_routine, ThreeIndInts03Tag)
            deallocate(ThreeIndInts04)
            call LogMemDealloc(this_routine, ThreeIndInts04Tag)
            deallocate(UMATTemp02)
            call LogMemDealloc(this_routine, UMATTemp02Tag)
        end if 

        deallocate(UMATTemp01)
        call LogMemDealloc(this_routine, UMATTemp01Tag)
        deallocate(SymLabelList2_rot)
        call LogMemDealloc(this_routine, SymLabelList2_rotTag)
        deallocate(SymLabelCounts2_rot)
        call LogMemDealloc(this_routine, SymLabelCounts2_rotTag)
        deallocate(SymLabelListInv_rot)
        call LogMemDealloc(this_routine, SymLabelListInv_rotTag)

        if (tShake) then
            deallocate(ShakeLambda)
            call LogMemDealloc(this_routine, ShakeLambdaTag)
            deallocate(ShakeLambdaNew)
            call LogMemDealloc(this_routine, ShakeLambdaNewTag)
            deallocate(Constraint)
            call LogMemDealloc(this_routine, ConstraintTag)
            deallocate(ConstraintCor)
            call LogMemDealloc(this_routine, ConstraintCorTag)
            deallocate(DerivConstrT1)
            call LogMemDealloc(this_routine, DerivConstrT1Tag)
            deallocate(DerivConstrT2)
            call LogMemDealloc(this_routine, DerivConstrT2Tag)
            deallocate(ForceCorrect)
            call LogMemDealloc(this_routine, ForceCorrectTag)
            deallocate(Correction)
            call LogMemDealloc(this_routine, CorrectionTag)
            if (tShakeApprox) then
                deallocate(DerivConstrT1T2Diag)
                call LogMemDealloc(this_routine, DerivConstrT1T2DiagTag)
            else
                deallocate(DerivConstrT1T2)
                call LogMemDealloc(this_routine, DerivConstrT1T2Tag)
            end if
        end if

    end subroutine DeallocateMem
 
end module RotateOrbsMod
