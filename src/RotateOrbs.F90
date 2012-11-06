MODULE RotateOrbsMod

    USE Global_utilities
    USE Parallel_neci 
    USE IntegralsData , only : UMAT,nFrozen,ChemPot
    USE UMatCache , only : UMatInd
    use constants, only: dp, PI
    USE SystemData , only : ConvergedForce,TimeStep,tLagrange,tShake,tShakeApprox,ShakeConverged
    use SystemData, only : tROIteration,ROIterMax,tShakeIter,ShakeIterMax,OrbEnMaxAlpha
    USE SystemData , only : G1,ARR,NEl,nBasis,LMS,ECore,tSeparateOccVirt,Brr,nBasisMax,OrbOrder
    use SystemData , only : lNoSymmetry,tRotatedOrbs,tERLocalization,tRotateOccOnly
    USE SystemData, only : tOffDiagMin,DiagWeight,OffDiagWeight,tRotateVirtOnly,tOffDiagSqrdMax
    use SystemData, only : tOffDiagSqrdMin,tOffDiagMax,tDoubExcMin,tOneElIntMax,tOnePartOrbEnMax
    USE SystemData, only : tShakeDelay,ShakeStart,tVirtCoulombMax,tVirtExchangeMin,MaxMinFac
    use SystemData, only : tMaxHLGap,tHijSqrdMin,OneElWeight,DiagMaxMinFac,OneElMaxMinFac
    USE SystemData, only : tDiagonalizehij,tHFSingDoubExcMax,tSpinOrbs,tReadInCoeff,tUseMP2VarDenMat
    use SystemData, only : tStoreSpinOrbs,tROHF,tFindCINatOrbs,tUseHFOrbs,tUEG
    USE Logging , only : tROHistogramAll,tROFciDump,tROHistER,tROHistOffDiag,tROHistDoubExc, tPrintRODump
    use Logging, only : tROHistSingExc,tROHistOnePartOrbEn,tROHistOneElInts,tROHistVirtCoulomb
    USE Logging , only : tPrintInts,tTruncRODump,NoTruncOrbs,NoDumpTruncs,tTruncDumpbyVal,TruncEvalues,tWriteTransMat
    USE OneEInts , only : TMAT2D
    USE SymData , only : TwoCycleSymGens,SymLabelList,SymLabelCounts
    USE Timing_neci , only : end_timing,print_timing_report
    USE Soft_exit, only : test_SOFTEXIT
    USE RotateOrbsData 
    use sort_mod
    use util_mod, only: get_free_unit
    IMPLICIT NONE
    INTEGER , ALLOCATABLE :: Lab(:,:),LabVirtOrbs(:),LabOccOrbs(:),SymLabelList3_rotInv(:)
    real(dp) , ALLOCATABLE :: CoeffCorT2(:,:),CoeffUncorT2(:,:)
    real(dp) , ALLOCATABLE :: Lambdas(:,:),ArrNew(:,:),ArrDiagNew(:),TMAT2DTemp(:,:),TMAT2DRot(:,:),TMAT2DPartRot01(:,:)
    real(dp) , ALLOCATABLE :: TMAT2DPartRot02(:,:)
    real(dp) , ALLOCATABLE :: DerivCoeff(:,:),UMATTemp01(:,:,:,:),UMATTemp02(:,:,:,:)
    real(dp) , ALLOCATABLE :: DerivLambda(:,:),ForceCorrect(:,:),Correction(:,:),ShakeLambdaNew(:),ConstraintCor(:)
    real(dp) , ALLOCATABLE :: Constraint(:),ShakeLambda(:),DerivConstrT1(:,:,:),DerivConstrT2(:,:,:),DerivConstrT1T2(:,:)
    real(dp) , ALLOCATABLE :: DerivConstrT1T2Diag(:),FourIndInts(:,:,:,:)
    real(dp) , ALLOCATABLE :: TwoIndInts01(:,:,:,:),TwoIndInts02(:,:,:,:),ThreeIndInts01(:,:,:,:),FourIndInts02(:,:,:,:)
    real(dp) , ALLOCATABLE :: ThreeIndInts02(:,:,:,:),ThreeIndInts03(:,:,:,:),ThreeIndInts04(:,:,:,:)  
    real(dp) , ALLOCATABLE :: DiagTMAT2Dfull(:),TMAT2DNew(:,:) 
    real(dp) , ALLOCATABLE :: TwoIndIntsER(:,:,:),ThreeIndInts01ER(:,:),ThreeIndInts02ER(:,:),FourIndIntsER(:)
    INTEGER(TagIntType) :: TwoIndIntsERTag,ThreeIndInts01ERTag,ThreeIndInts02ERTag,FourIndIntsERTag
    INTEGER(TagIntType) :: TwoIndInts01Tag,TwoIndInts02Tag,ThreeIndInts01Tag,ThreeIndInts02Tag,ThreeIndInts03Tag,ThreeIndInts04Tag
    INTEGER(TagIntType) :: FourIndInts02Tag
    INTEGER(TagIntType) :: TMAT2DTempTag,TMAT2DRotTag,TMAT2DPartRot01Tag,TMAT2DPartRot02Tag
    INTEGER(TagIntType) :: LabTag,ForceCorrectTag,CorrectionTag,FourIndIntsTag,ArrDiagNewTag,ArrNewTag,UMATTemp01Tag,UMATTemp02Tag
    INTEGER :: ShakeIterInput,NoOcc,LowBound02,HighBound02,Iteration,TotNoConstraints
    INTEGER(TagIntType) :: CoeffCorT2Tag,CoeffUncorT2Tag,LambdasTag,DerivCoeffTag,DerivLambdaTag
    INTEGER(TagIntType) :: ShakeLambdaNewTag
    INTEGER(TagIntType) :: ShakeLambdaTag,ConstraintTag,ConstraintCorTag,DerivConstrT1Tag,DerivConstrT2Tag,DerivConstrT1T2Tag
    INTEGER(TagIntType) :: DerivConstrT1T2DiagTag
    INTEGER(TagIntType) :: LabVirtOrbsTag,LabOccOrbsTag
    INTEGER :: MinOccVirt,MaxOccVirt,MinMZ,MaxMZ,error,LowBound,HighBound
    INTEGER :: NoInts01,NoInts02,NoInts03,NoInts04,NoInts05,NoInts06
    INTEGER(TagIntType) :: DiagTMAT2DfullTag,TMAT2DNewTag,SymLabelList3_rotInvTag
    LOGICAL :: tNotConverged,tInitIntValues
    real(dp) :: OrthoNorm,ERPotEnergy,HijSqrdPotEnergy,OffDiagPotEnergy,CoulPotEnergy,PotEnergy,Force,TwoEInts,DistCs
    real(dp) :: OrthoForce,DistLs,LambdaMag,PEInts,PEOrtho
    real(dp) :: ForceInts,TotCorrectedForce
    real(dp) :: ijOccVirtPotEnergy,EpsilonMin,MaxTerm
    real(dp) :: DiagOneElPotInit,ERPotInit,ijVirtOneElPotInit,ijVirtCoulPotInit,ijVirtExchPotInit
    real(dp) :: singCoulijVirtInit,singExchijVirtInit,singCoulconHFInit,singExchconHFInit,ijklPotInit,ijklantisymPotInit
    real(dp) :: ijOccVirtOneElPotInit,ijOccVirtCoulPotInit,ijOccVirtExchPotInit
    real(dp) :: OrthoFac=1.0_dp,ROHistSing(2,4002),ROHistOffDiag(2,4002),ROHistDoubExc(2,4002),ROHistER(2,4002)
    real(dp) :: ROHistHijVirt(2,4002),ROHistHijOccVirt(2,4002),ROHistHii(2,4002)
    real(dp) :: ROHistOnePartOrbEn(2,4002),ROHistDCijOcklVir(2,4002),ROHistDEijOcklVir(2,4002),ROHistDCijklVir(2,4002)
    real(dp) :: ROHistDEijklVir(2,4002)
    real(dp) :: ROHistSCikOcjVir(2,4002),ROHistSEikOcjVir(2,4002),ROHistSCkOcijVir(2,4002),ROHistSEkOcijVir(2,4002)
    real(dp) :: ROHistSCijkVir(2,4002),ROHistSEijkVir(2,4002)
    real(dp) :: ROHistSASikOcjVir(2,4002),ROHistSASkOcijVir(2,4002),ROHistSASijkVir(2,4002),ROHistASijklVir(2,4002)
    real(dp) :: ROHistASijOcklVir(2,4002)
    TYPE(timer), save :: Rotation_Time,FullShake_Time,Shake_Time,Findtheforce_Time,Transform2ElInts_Time
    type(timer), save :: findandusetheforce_time,CalcDerivConstr_Time,TestOrthoConver_Time
    TYPE(timer), save :: RefillUMAT_Time,PrintROFCIDUMP_Time
! In this routine, alpha (a), beta (b), gamma (g) and delta (d) refer to the unrotated (HF) orbitals where possible such that < a b | g d > is an unrotated four index integral.   
! For the rotated orbitals, the letter i,j,k and l are generally used, i.e. < i j | k l > refers to a transformed four index integral.
! Differentiation of the potential energy (to find the force) is done with respect to coefficient c(z,m) (or c(a,m)), where zeta (z) or a refers to the HF index, and m to the rotated.
    

    contains

    SUBROUTINE RotateOrbs()


        IF(iProcIndex.eq.Root) THEN

! If we are reading in our own transformation matrix (coeffT1) don't need a lot of the initialisation stuff.
            IF(tReadInCoeff.or.tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs) THEN

                tNotConverged=.false.
                CALL FindNatOrbitals()

            ELSE
! Need to actually find the coefficient matrix and then use it.

                tNotConverged=.true.
                CALL InitLocalOrbs()        ! Set defaults, allocate arrays, write out headings for OUTPUT, set integarals to HF values.

                IF(tDiagonalizehij) THEN

                    CALL Diagonalizehij()
                    tNotConverged=.false.
                    Iteration=2

                ELSEIF(tMaxHLGap) THEN
                    CALL EquateDiagFock()
                    tNotConverged=.false.

                ELSE

                    tNotConverged=.true.

                    CALL WriteStats()           ! write out the original stats before any rotation.
                   
                    CALL set_timer(Rotation_Time,30)

                    do while(tNotConverged)     ! rotate the orbitals until the sum of the four index integral falls below a chose convergence value.

                        Iteration=Iteration+1
                        
                        CALL FindNewOrbs()      ! bulk of the calculation.
                                                ! do the actual transformations, moving the coefficients by a timestep according to the calculated force. 

                        CALL WriteStats()       ! write out the stats for this iteration.

                    enddo           

                    CALL halt_timer(Rotation_Time)
                    
                    WRITE(6,*) "Convergence criterion met. Finalizing new orbitals..."

                ENDIF


! Make symmetry, orbitals, one/two-electron integrals consistent with rest of NECI
                CALL FinalizeNewOrbs()


!        CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
                CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)

                CALL DeallocateMem()

            ENDIF            

            CALL neci_flush(6)
            CALL neci_flush(transform_unit)
        ENDIF


    END SUBROUTINE RotateOrbs


    SUBROUTINE FindNatOrbitals()
! This routine simply takes a transformation matrix and rotates the integrals to produce a new FCIDUMP file.
! In one case the transformation matrix is read in from a file TRANSFORMMAT.
! In the other, the transformation matrix is calculated from the MP2 variational density matrix.

! MP2VDM = D2_ab = sum_ijc [ t_ij^ac ( 2 t_ij^bc - t_ji^bc ) ]
! Where :  t_ij^ac = - < ab | ij > / ( E_a - E_i + E_b - Ej )
! Ref : J. Chem. Phys. 131, 034113 (2009) - note: in Eqn 1, the cb indices are the wrong way round (should be bc).
        USE NatOrbsMod , only : SetUpNatOrbLabels,FindNatOrbs,FillCoeffT1,DeallocateNatOrbs,PrintOccTable
        INTEGER :: i,a,ierr,MinReadIn,MaxReadIn, iunit
        CHARACTER(len=*) , PARAMETER :: this_routine='FindNatOrbitals'


        IF(tUseMP2VarDenMat) WRITE(6,*) '*** Transforming the HF orbitals into the MP2 approximate natural orbitals. ***'
        IF(tFindCINatOrbs) THEN
            WRITE(6,*) '*** Transforming the HF orbitals into approximate natural orbitals'
            WRITE(6,*) 'based on the one-electron density matrix found from the wavefunction calculated above. ***'
        ENDIF

        IF(tSpinOrbs) THEN
            IF(.not.tStoreSpinOrbs) THEN
                WRITE(6,*) "We want to use spin orbitals - turning on tStoreSpinOrbs."
                tStoreSpinOrbs=.true.
            ENDIF
        ENDIF

        IF(tROHF.and.tStoreSpinOrbs) CALL Stop_All(this_routine,"Cannot compress open shell systems into spatial " &
            & //"orbitals when rotating, turn off ROHF.")

        IF(tTruncRODump.and.(.not.tTruncDumpbyVal)) THEN 
            NoFrozenVirt=NoTruncOrbs(1)
        ELSEIF(tTruncRODump) THEN
            ! If the 'number of frozen orbitals' is given as a cutoff - take NoFrozenVirt to be 0 for all the allocation purposes - will set this later when
            ! we have the eigenvalues and know how many orbitals lie below it.
            NoFrozenVirt=0
            TruncEval=TruncEvalues(1)
        ELSE
            NoFrozenVirt=0
        ENDIF

        SpatOrbs=nBasis/2
        IF(tStoreSpinOrbs) THEN
            NoOrbs=nBasis
            NoOcc=NEl
            MinReadIn=1
            MaxReadIn=nBasis
            IF(tRotateVirtOnly) MinReadIn=NEl+1
            IF(tRotateOccOnly) MaxReadIn=NEl
            ! If tStoreSpinOrbs ARR(:,2) is not filled, but we want to use it later, so just fill it here.            
            do i=1,NoOrbs
                ARR(BRR(i),2)=ARR(i,1)
            enddo
            ALLOCATE(SymLabelCounts2_rot(2,32),stat=ierr)
            CALL LogMemAlloc('SymLabelCounts2_rot',2*32,4,this_routine,SymLabelCounts2_rotTag,ierr)
            SymLabelCounts2_rot(:,:)=0
            ! first 8 refer to the occupied, and the second to the virtual beta spin.
            ! third and fourth to the occupied and virtual alpha spin.
 
        ELSE
            NoOrbs=SpatOrbs
            NoOcc=NEl/2
            MinReadIn=1
            MaxReadIn=SpatOrbs
            IF(tRotateVirtOnly) MinReadIn=(NEl/2)+1
            IF(tRotateOccOnly) MaxReadIn=NEl/2
            ALLOCATE(SymLabelCounts2_rot(2,16),stat=ierr)
            CALL LogMemAlloc('SymLabelCounts2_rot',2*16,4,this_routine,SymLabelCounts2_rotTag,ierr)
            SymLabelCounts2_rot(:,:)=0
            ! first 8 refer to the occupied, and the second to the virtual.

        ENDIF
        NoRotOrbs=NoOrbs

        CALL ApproxMemReq()

!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
!        CALL neci_flush(6)
!        CALL Stop_All('','')


        ALLOCATE(SymLabelList2_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList2_rot',NoOrbs,4,this_routine,SymLabelList2_rotTag,ierr)
        SymLabelList2_rot(:)=0                     
        ALLOCATE(SymLabelList3_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3_rot',NoOrbs,4,this_routine,SymLabelList3_rotTag,ierr)
        SymLabelList3_rot(:)=0                     
 
        ALLOCATE(SymLabelListInv_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelListInv_rot',NoOrbs,4,this_routine,SymLabelListInv_rotTag,ierr)
        SymLabelListInv_rot(:)=0                     


        IF(tReadInCoeff.or.tUseHFOrbs) THEN
! No symmetry, so no reordering of the orbitals - symlabellist just goes from 1-NoOrbs.        
! When we are just reading in the coefficients and transforming, it does not matter about the ordering of the orbitals.
            do i=1,NoOrbs
                SymLabelList2_rot(i)=i
                SymLabelListInv_rot(i)=i
            enddo
!        ELSEIF(tUseMP2VarDenMat) THEN
! When we are calculating the MP2VDM ourselves, we need to be able to run over occupied and virtual separately, so the 
! orbitals are always separated.
! Also, we want to have the option of maintaining symmetry, so SymLabelList2_rot and SymLabelCounts2_rot are both constructed so 
! that the orbitals are labelled by symmetry within the occupied and virtual.
!            tSeparateOccVirt=.true.
!            CALL InitSymmArrays()

        ELSEIF(tFindCINatOrbs.or.tUseMP2VarDenMat) THEN

            CALL SetupNatOrbLabels() 

        ENDIF

!        OPEN(42,file="TRANSFORMMAT",status="old")
!        do i=MinReadIn,MaxReadIn
!            j=SymLabelList2_rot(i)
!            do a=MinReadIn,MaxReadIn
!                b=SymLabelList2_rot(a)
!                READ(42,*) CoeffT1(b,j)
!                READ(42,*) CoeffT1
!            enddo
!        enddo
!        CLOSE(42)
         
! Need to read to convert the UMAT matrix from UMATInd to the appropriate indexing for Transform2ElInts.        
! This just contains all the untransformed orbitals.
!        ALLOCATE(UMATTemp01(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!        CALL LogMemAlloc('UMATTemp01',NoOrbs**4,8,this_routine,UMATTemp01Tag,ierr)

!        CALL CopyAcrossUMAT()
        

! Yet another labelling system, SymLabelList3_rot is created here.
! This indicates the label of the transformed orbital.
! In the case where we are truncating the space, the transformed orbitals are ordered according to the size of the eigenvalues of the MP2VDM
! matrix when it is diagonalised.  We wish to keep them in this order when transforming the integrals etc, so that when we truncate the last 
! NoFrozenVirt orbitals, we are removing those with the smallest MP2VDM eigenvalues (occupation numbers).
! In the case where no truncation is made however, SymLabelList3_rot is the same as SymLabelList2_rot, so that the indexes remain the same as previously.
! This allows for the option of going straight into a spawning calc from the rotation, which is not possible when a truncation is performed 
! because of the messed up indices.
        IF(tTruncRODump) THEN
            IF(MOD(NoFrozenVirt,2).ne.0) CALL Stop_All(this_routine,"Must freeze virtual spin orbitals in pairs of 2.")
            IF(tStoreSpinOrbs) THEN
                NoRotOrbs=NoOrbs-NoFrozenVirt
            ELSE
                NoFrozenVirt=NoFrozenVirt/2
                NoRotOrbs=NoOrbs-NoFrozenVirt
            ENDIF            
            do i=1,NoOrbs
                SymLabelList3_rot(i)=i
            enddo
        ELSE
            do i=1,NoOrbs
                SymLabelList3_rot(i)=SymLabelList2_rot(i)
            enddo
        ENDIF


! The last two indices of these are the transformed and possibly truncated orbitals.        

!        ALLOCATE(TwoIndInts01(NoOrbs,NoOrbs,NoRotOrbs,NoRotOrbs),stat=ierr)
!        CALL LogMemAlloc('TwoIndInts01',(NoOrbs**4)*(NoRotOrbs**2),8,this_routine,TwoIndInts01Tag,ierr)

!        ALLOCATE(FourIndInts(NoRotOrbs,NoRotOrbs,NoOrbs,NoOrbs),stat=ierr)
!        CALL LogMemAlloc('FourIndInts',(NoOrbs**2)*(NoRotOrbs**2),8,this_routine,FourIndIntsTag,ierr)

        ALLOCATE(CoeffT1(NoOrbs,NoRotOrbs),stat=ierr)
        CALL LogMemAlloc(this_routine,NoRotOrbs*NoOrbs,8,this_routine,CoeffT1Tag,ierr)
        CoeffT1(:,:)=0.0_dp
        IF(tSeparateOccVirt) THEN
            do i=1,NoRotOrbs
                CoeffT1(i,i)=1.0_dp
            enddo
        ENDIF


        IF(tUEG) THEN

            CALL FindNatOrbs()

            CALL FillCoeffT1()


        ELSE
            IF(tReadInCoeff) THEN

                WRITE(6,'(A)') " Reading in the transformation matrix from TRANSFORMMAT, and using this to rotate the HF orbitals."

!                OPEN(72,FILE='TRANSFORMMAT',status='old')
!                READ(72,*) CoeffT1
!                CLOSE(72)

                
                iunit = get_free_unit()
                OPEN(iunit,FILE='TRANSFORMMAT',status='old')
                do i=1,NoOrbs
                    do a=1,NoOrbs
                        READ(iunit,*) CoeffT1(a,i)
                    enddo
                enddo
                CLOSE(iunit)
          
!                OPEN(78,FILE='TRANSFORMMATORIG',status='unknown')
!                do i=1,NoOrbs
!                    do a=1,NoOrbs
!                        WRITE(78,*) i,a,CoeffT1(i,a)
!                    enddo
!                enddo
!                stop

            ELSEIF(tFindCINatOrbs.or.tUseMP2VarDenMat.or.tUseHFOrbs) THEN

                
                IF(.not.tUseHFOrbs) CALL FindNatOrbs()
                
                ! Fill the coefficient matrix with the eigenvectors of the OneRDM.
                ! Find out the ordering ...need to read in according to SymLabelList2_rot, so that the transformation is all o.k.
                ! If both the HF and transformed indices are done this way, the symmetries and everything should be fine.
!                do i=1,NoRotOrbs
!                    i2=SymLabelList2_rot(i)
!                    do a=1,NoOrbs
!                        a2=SymLabelList2_rot(a)
!                        Coeff(a,i)=OneRDM(a2,i2)
!                    enddo
!                enddo
            

!           ELSEIF(tUseMP2VarDenMat) THEN
! This bit generates the MP2 variational density matrix, and uses this as the transformation matrix (CoeffT1).        
    
!                WRITE(6,*) "Calculating the MP2 vartiational density matrix, and using this to rotate the HF orbitals."
!                CALL CalcMP2VarDenMat()
            
                
                IF(tUseHFOrbs) THEN
                    CALL PrintOccTable()
                ELSE
                    CALL FillCoeffT1()
                ENDIF

            ENDIF

            IF(tPrintRODump) THEN
                ALLOCATE(FourIndInts(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
                CALL LogMemAlloc('FourIndInts',(NoOrbs**4),8,this_routine,FourIndIntsTag,ierr)

! Then, transform2ElInts
                WRITE(6,*) 'Transforming the four index integrals'
                CALL Transform2ElIntsMemSave()

                WRITE(6,*) 'Re-calculating the fock matrix'
                CALL CalcFOCKMatrix()

                WRITE(6,*) 'Refilling the UMAT and TMAT2D'
! The ROFCIDUMP is also printed out in here.        
                CALL RefillUMATandTMAT2D()        

                CALL neci_flush(6)


                IF((tFindCINatOrbs.or.tUseMP2VarDenMat).and.(NoDumpTruncs.gt.1)) CALL ReTruncROFciDump()

                IF((.not.tUseHFOrbs).and.(.not.tReadInCoeff)) CALL DeallocateNatOrbs()
            ENDIF

            IF(tWriteTransMat) CALL WriteTransformMat()

     
! If a truncation is being made, the new basis will not be in the correct energetic ordering - this does not matter, as we
! never go straight into a spawning and they will be reordered when the ROFCIDUMP file is read in again. 
            CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)

            DEALLOCATE(CoeffT1)
            CALL LogMemDeAlloc(this_routine,CoeffT1Tag)
            DEALLOCATE(SymLabelList2_rot)
            CALL LogMemDeAlloc(this_routine,SymLabelList2_rotTag)
            DEALLOCATE(SymLabelListInv_rot)
            CALL LogMemDeAlloc(this_routine,SymLabelListInv_rotTag)
            IF(tPrintRODump) THEN
                DEALLOCATE(FourIndInts)
                CALL LogMemDeAlloc(this_routine,FourIndIntsTag)
            ENDIF
        ENDIF

!        DEALLOCATE(UMATTemp01)
!        CALL LogMemDeAlloc(this_routine,UMATTemp01Tag)
!        DEALLOCATE(TwoIndInts01)
!        CALL LogMemDeAlloc(this_routine,TwoIndInts01Tag)

    END SUBROUTINE FindNatOrbitals 


    SUBROUTINE ReTruncROFciDump()
        USE NatOrbsMod , only : FillCoeffT1
        INTEGER :: i,j,ierr
        CHARACTER(len=*) , PARAMETER :: this_routine='ReTruncROFciDump'


        do i=2,NoDumpTruncs

            DEALLOCATE(ArrDiagNew)
            CALL LogMemDeAlloc(this_routine,ArrDiagNewTag)
            DEALLOCATE(CoeffT1)
            CALL LogMemDeAlloc(this_routine,CoeffT1Tag)
            DEALLOCATE(FourIndInts)
            CALL LogMemDeAlloc(this_routine,FourIndIntsTag)
            DEALLOCATE(SymOrbs_rot)
            CALL LogMemDeAlloc(this_routine,SymOrbs_rotTag)
            DEALLOCATE(TMAT2DNew)
            CALL LogMemDeAlloc(this_routine,TMAT2DNewTag)
            DEALLOCATE(EvaluesTrunc)
            CALL LogMemDeAlloc(this_routine,EvaluesTruncTag)


            IF(tTruncDumpbyVal) THEN
                NoFrozenVirt=0
                TruncEval=TruncEvalues(i)
            ELSE
                IF(tStoreSpinOrbs) THEN
                    NoFrozenVirt=NoTruncOrbs(i)
                ELSE
                    NoFrozenVirt=NoTruncOrbs(i)/2
                ENDIF            
            ENDIF
            NoRotOrbs=NoOrbs-NoFrozenVirt
 
            IF(MOD(NoFrozenVirt,2).ne.0) CALL Stop_All(this_routine,"Must freeze virtual spin orbitals in pairs of 2.")


            ALLOCATE(CoeffT1(NoOrbs,NoRotOrbs),stat=ierr)
            CALL LogMemAlloc(this_routine,NoRotOrbs*NoOrbs,8,this_routine,CoeffT1Tag,ierr)
            CoeffT1(:,:)=0.0_dp
            IF(tSeparateOccVirt) THEN
                do j=1,NoRotOrbs
                    CoeffT1(i,i)=1.0_dp
                enddo
            ENDIF

            CALL FillCoeffT1()


            ALLOCATE(FourIndInts(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('FourIndInts',(NoOrbs**4),8,this_routine,FourIndIntsTag,ierr)

! Then, transform2ElInts
            WRITE(6,*) 'Transforming the four index integrals.'
            CALL Transform2ElIntsMemSave()

            WRITE(6,*) 'Re-calculating the fock matrix.'
            CALL CalcFOCKMatrix()

            WRITE(6,*) 'Refilling the UMAT and TMAT2D.'
! The ROFCIDUMP is also printed out in here.        
            CALL RefillUMATandTMAT2D()        

            CALL neci_flush(6)

        enddo


    ENDSUBROUTINE ReTruncROFciDump



    SUBROUTINE ApproxMemReq()
! This routine makes a quick sum of the memory that will be require to transform the integrals from the HF to the new basis.

! Main arrays required are:
        MemAllocRot=0
        
! Symmetry/Labelling:
!   - SymLabelLists(NoOrbs) x 3 
!   - SymLabelCounts(32/16 - Spin/Spat)
        MemAllocRot=MemAllocRot+(3*NoOrbs*4)
        MemAllocRot=MemAllocRot+(32*4)
        
! Finding transformation matrices
!   - NatOrbsMat(NoOrbs,NoOrbs)
!   - Evalues(NoOrbs) x 2
        MemAllocRot=MemAllocRot+((NoOrbs**2)*8)
        MemAllocRot=MemAllocRot+(2*NoOrbs*8)


! Transformation of integrals
!   - CoeffT1(NoOrbs,NoRotOrbs) 
!   - FourIndInts(NoRotOrbs,NoRotOrbs,NoOrbs,NoOrbs)
!   - Temp4indints(NoRotOrbs,NoOrbs)
        IF(tPrintRODump) THEN
            MemAllocRot=MemAllocRot+(NoOrbs*NoRotOrbs*8*2)
            MemAllocRot=MemAllocRot+((NoRotOrbs**2)*(NoOrbs**2)*8)

! Transform fock
!   - ArrNew(NoOrbs) - reduce this?
            MemAllocRot=MemAllocRot+(NoOrbs*8)

! RefillTMAT2D
!   - TMAT2D(nBasis,nBasis) 
            MemAllocRot=MemAllocRot+((nBasis**2)*8)
        ENDIF

        WRITE(6,'(A72,F20.10,A15)') "Rough estimate of the memory required for the orbital transformation = " &
            ,REAL(MemAllocRot,dp)/1048576.0_dp," Mb/Processor"


    END SUBROUTINE ApproxMemReq



    SUBROUTINE WriteTransformMat()
        INTEGER :: w,x,i,a,b,iunit

! This file is printed to be used to produce cube files from QChem.
! Line 1 is the coefficients of HF spatial orbitals 1 2 3 ... which form transformed orbital 1 etc.

        iunit = get_free_unit()
        OPEN(iunit,FILE='MOTRANSFORM',FORM='UNFORMATTED',access='direct', recl=8)
! Need to put this back into the original order. 

        x = 0
        IF(tStoreSpinOrbs) THEN
            do i=1,NoOrbs-1,2
!                j=SymLabelListInv_rot(i)
                ! SymLabelList2_rot(i) gives the orbital label (from Dalton or QChem) corresponding to our
                ! label i.
                ! SymLabelListInv_rot(j) therefore gives the label used in CoeffT1 corresponding to the
                ! Qchem/Dalton label j.
                    
                do a=1,NoOrbs-1,2
                    b=SymLabelListInv_rot(a)
!                    WRITE(iunit,rec=x) CoeffT1(b,j)
                    WRITE(iunit,rec=x) CoeffT1(b,i)
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                enddo
            enddo
            do i=2,NoOrbs,2
                do a=2,NoOrbs,2
                    b=SymLabelListInv_rot(a)
!                    WRITE(iunit,rec=x) CoeffT1(b,j)
                    WRITE(iunit,rec=x) CoeffT1(b,i)
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                enddo
            enddo
        ELSE
            w=1
            x=1   !keep a counter of record number
            do while (w.le.2)
                do i=1,SpatOrbs
!                    j=SymLabelListInv_rot(i)
                    ! SymLabelList2_rot(i) gives the orbital label (from Dalton or QChem) corresponding to our
                    ! label i.
                    ! SymLabelListInv_rot(j) therefore gives the label used in CoeffT1 corresponding to the
                    ! Qchem/Dalton label j.
                    do a=1,SpatOrbs
                        b=SymLabelListInv_rot(a)
!                        WRITE(iunit,rec=x) CoeffT1(b,j)
                        WRITE(iunit,rec=x) CoeffT1(b,i)
                        x=x+1
                        ! a/b are the original (HF) orbitals, and i/j the transformed
                    enddo
                enddo
                w=w+1
                ! print the whole matrix twice, once for alpha spin, once for beta.
            enddo
        ENDIF
        CLOSE(iunit)
 
        OPEN(iunit,FILE='MOTRANSFORM02')
! Need to put this back into the original order. 
        w=1
        x=1   !keep a counter of record number
        do while (w.le.2)
            do i=1,SpatOrbs
!                j=SymLabelListInv_rot(i)
                ! SymLabelList2_rot(i) gives the orbital label (from Dalton or QChem) corresponding to our
                ! label i.
                ! SymLabelListInv_rot(j) therefore gives the label used in CoeffT1 corresponding to the
                ! Qchem/Dalton label j.
                do a=1,SpatOrbs
                    b=SymLabelListInv_rot(a)
!                    WRITE(iunit,'(F20.10)',advance='no') CoeffT1(b,j)
                    WRITE(iunit,'(F20.10)',advance='no') CoeffT1(b,i)
                    x=x+1
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                enddo
                WRITE(iunit,*) ''
            enddo
            w=w+1
            ! print the whole matrix twice, once for alpha spin, once for beta.
        enddo
        CLOSE(iunit)

        OPEN(iunit,FILE='TRANSFORMMAT',status='unknown')
        do i=1,NoOrbs
!            j=SymLabelListInv_rot(i)
            do a=1,NoOrbs
                b=SymLabelListInv_rot(a)
!                WRITE(iunit,*) CoeffT1(b,j)
                WRITE(iunit,*) CoeffT1(b,i)
            enddo
        enddo
        CALL neci_flush(iunit)
        CLOSE(iunit)
      


    END SUBROUTINE WriteTransformMat



   
    SUBROUTINE InitLocalOrbs()
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        INTEGER :: ierr

! Writing to output which PE is being maximised/minimised.        
        WRITE(6,*) '*****'
        IF(tERLocalization) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on Edmiston-Reudenberg localisation,"
            WRITE(6,*) "i.e. maximisation of the <ii|ii> integrals..."
            WRITE(6,*) "*****"
        ENDIF
        IF(tVirtCoulombMax) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on maximisation of the sum of the"
            WRITE(6,*) "<ij|ij> integrals, where i and j are both virtuals..."
            WRITE(6,*) "*****"
        ENDIF
        IF(tOffDiagSqrdMin) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on mimimisation "
            WRITE(6,*) "of <ij|kl>^2 integrals..."
            WRITE(6,*) "*****"
        ENDIF
        IF(tOffDiagMin) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on mimimisation "
            WRITE(6,*) "of <ij|kl> integrals..."
            WRITE(6,*) "*****"
        ENDIF
        IF(tDoubExcMin) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on mimimisation "
            WRITE(6,*) "of the double excitation hamiltonian elements."
            WRITE(6,*) "*****"
        ENDIF
        IF(tOnePartOrbEnMax) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on maximisation "
            WRITE(6,*) "of the virtual one particle orbital energies."
            WRITE(6,*) "*****"
        ELSEIF(tMaxHLGap) THEN
!This will transform all the orbitals within a particlar group to have the same diagonal fock matrix element.
            WRITE(6,*) "Transforming orbitals based on equating their diagonal fock matrix elements."
            WRITE(6,*) "*****"
        ENDIF

! Writing out which orthonormalisation method is being used...       
        IF(tLagrange) THEN
            IF(tShake) THEN
                CALL neci_flush(6)
                CALL Stop_All(this_routine,"ERROR. Both LAGRANGE and SHAKE keywords present in the input. &
                & These two orthonormalisation methods clash.")
            ENDIF
            WRITE(6,*) "Using a Lagrange multiplier to attempt to rotate orbitals in a way to maintain orthonormality"
        ELSEIF (tShake) THEN
            WRITE(6,*) "Using the shake algorithm to iteratively find lambdas which maintain "
            WRITE(6,*) "orthonormalisation with rotation"
        ELSE
            WRITE(6,*) "Explicity reorthonormalizing orbitals after each rotation."
        ENDIF
        
! Check for a few possible errors.
        IF(.not.TwoCycleSymGens) THEN
            CALL neci_flush(6)
            CALL Stop_All(this_routine,"ERROR. TwoCycleSymGens is false.  Symmetry is not abelian.") 
        ENDIF
        IF((tRotateOccOnly.or.tRotateVirtOnly).and.(.not.tSeparateOccVirt)) THEN
            tSeparateOccVirt=.true.
            WRITE(6,*) "NOTE. Cannot rotate only occupied or virtual without first separating them."
            WRITE(6,*) "SEPARATEOCCVIRT keyword is being turned on."
        ENDIF        
        IF((tOffDiagSqrdMax.and.tOffDiagSqrdMin).or.(tOffDiagMax.and.tOffDiagMin)) THEN
            CALL neci_flush(6)
            CALL Stop_All(this_routine,"ERROR. Cannot both maximise and minimise off diagonal elements simultaneously")
        ENDIF
        IF(tOnePartOrbEnMax.and.(.not.tSeparateOccVirt)) THEN
            CALL neci_flush(6)
            CALL Stop_All(this_routine, &
            "ERROR. Cannot currently maximise the one particle orbital energies without separating occupied and virtual.") 
        ENDIF
        WRITE(6,*) "*****"

!Zero values.
        OrthoNorm=0.0_dp
        ERPotEnergy=0.0_dp
        PotEnergy=0.0_dp
        Force=0.0_dp
        TwoEInts=0.0_dp
        PEInts=0.0_dp
        PEOrtho=0.0_dp
        ForceInts=0.0_dp
        DistCs=0.0_dp
        DistLs=0.0_dp
        LambdaMag=0.0_dp
        SpatOrbs=nBasis/2
        IF(tStoreSpinOrbs) THEN
            NoOrbs=nBasis
            NoOcc=NEl
        ELSE
            NoOrbs=SpatOrbs
            NoOcc=NEl/2
        ENDIF
        NoRotOrbs=NoOrbs
        Iteration=0
        OrthoForce=0.0_dp
        ShakeIterInput=ShakeIterMax
        TotNoConstraints=(NoOrbs*(NoOrbs+1))/2

!When maximising the one particle orbital energies, choose the zero value (Epsilon min).
        IF(tRotateVirtOnly.and.tOnePartOrbEnMax) THEN
            EpsilonMin=ARR(NEl+1,1)
            WRITE(6,*) 'Taking EpsilonMin to be the LUMO of the HF orbitals...'
            WRITE(6,*) 'EpsilonMin = ',EpsilonMin
        ELSEIF(tOnePartOrbEnMax) THEN
            EpsilonMin=ChemPot
            WRITE(6,*) 'Taking EpsilonMin to be the chemical potential (midway between HF HOMO and LUMO)...'
            WRITE(6,*) 'therefore EpsilonMin = ',EpsilonMin
        ENDIF

!Set timed routine names
        Rotation_Time%timer_name='RotateTime'
        Shake_Time%timer_name='ShakeTime'
        FullShake_Time%timer_name='FullShakeTime'
        Findtheforce_Time%timer_name='FindtheForceTime'
        Transform2ElInts_Time%timer_name='Transform2ElIntsTime'
        findandusetheforce_time%timer_name='Findandusetheforce'
        CalcDerivConstr_Time%timer_name='CalcDerivConstr'
        TestOrthoConver_Time%timer_name='TestOrthoConver'
       

!Allocate memory

!        ALLOCATE(CoeffT1Temp(NoOrbs,NoOrbs),stat=ierr)
!        CALL LogMemAlloc('CoeffT1Temp',NoOrbs**2,8,this_routine,CoeffT1TempTag,ierr)
        ALLOCATE(CoeffT1(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffT1',NoOrbs**2,8,this_routine,CoeffT1Tag,ierr)
        ALLOCATE(CoeffCorT2(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffCorT2',NoOrbs**2,8,this_routine,CoeffCorT2Tag,ierr)
        ALLOCATE(CoeffUncorT2(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffUncT2',NoOrbs**2,8,this_routine,CoeffUncorT2Tag,ierr)
        CoeffUncorT2(:,:)=0.0_dp
         
!        ALLOCATE(DerivCoeffTemp(NoOrbs,NoOrbs),stat=ierr)
!        CALL LogMemAlloc('DerivCoeffTemp',NoOrbs**2,8,this_routine,DerivCoeffTempTag,ierr)
        ALLOCATE(DerivCoeff(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('DerivCoeff',NoOrbs**2,8,this_routine,DerivCoeffTag,ierr)
  
        ALLOCATE(DiagTMAT2Dfull(NoOrbs-(NoOcc)),stat=ierr)
        CALL LogMemAlloc('DiagTMAT2Dfull',(NoOrbs-(NoOcc)),8,this_routine,DiagTMAT2DfullTag,ierr)
        ALLOCATE(UMATTemp01(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('UMATTemp01',NoOrbs**4,8,this_routine,UMATTemp01Tag,ierr)
        ALLOCATE(TwoIndInts01(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('TwoIndInts01',NoOrbs**4,8,this_routine,TwoIndInts01Tag,ierr)
!        ALLOCATE(ThreeIndInts02Temp(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!        CALL LogMemAlloc('ThreeIndInts02Temp',NoOrbs**4,8,this_routine,ThreeIndInts02TempTag,ierr)
        ALLOCATE(ThreeIndInts02(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts02',NoOrbs**4,8,this_routine,ThreeIndInts02Tag,ierr)
!        ALLOCATE(FourIndInts02Temp(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!        CALL LogMemAlloc('FourIndInts02Temp',NoOrbs**4,8,this_routine,FourIndInts02TempTag,ierr)
        ALLOCATE(FourIndInts(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts',NoOrbs**4,8,this_routine,FourIndIntsTag,ierr)
        ALLOCATE(FourIndInts02(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts02',NoOrbs**4,8,this_routine,FourIndInts02Tag,ierr)

        ! Partially transformed temporary arrays.
        IF(tERLocalization.and.(.not.tStoreSpinOrbs)) THEN
            ALLOCATE(TwoIndIntsER(NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('TwoIndIntsER',NoOrbs**3,8,this_routine,TwoIndIntsERTag,ierr)
            ALLOCATE(ThreeIndInts01ER(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ThreeIndInts01ER',NoOrbs**2,8,this_routine,ThreeIndInts01ERTag,ierr)
            ALLOCATE(ThreeIndInts02ER(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ThreeIndInts02ER',NoOrbs**2,8,this_routine,ThreeIndInts02ERTag,ierr)
            ALLOCATE(FourIndIntsER(NoOrbs),stat=ierr)
            CALL LogMemAlloc('FourIndIntsER',NoOrbs,8,this_routine,FourIndIntsERTag,ierr)
        ELSE
            ALLOCATE(TMAT2DTemp(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('TMAT2DTemp',NoOrbs**2,8,this_routine,TMAT2DTempTag,ierr)
            ALLOCATE(TMAT2DPartRot01(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('TMAT2DPartRot01',NoOrbs**2,8,this_routine,TMAT2DPartRot01Tag,ierr)
            ALLOCATE(TMAT2DPartRot02(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('TMAT2DPartRot02',NoOrbs**2,8,this_routine,TMAT2DPartRot02Tag,ierr)
            ALLOCATE(TMAT2DRot(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('TMAT2DRot',NoOrbs**2,8,this_routine,TMAT2DRotTag,ierr)
            ALLOCATE(UMATTemp02(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('UMATTemp02',NoOrbs**4,8,this_routine,UMATTemp02Tag,ierr)
 
!            ALLOCATE(TwoIndInts02Temp(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!            CALL LogMemAlloc('TwoIndInts02Temp',NoOrbs**4,8,this_routine,TwoIndInts02TempTag,ierr)
!            ALLOCATE(ThreeIndInts01Temp(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!            CALL LogMemAlloc('ThreeIndInts01Temp',NoOrbs**4,8,this_routine,ThreeIndInts01TempTag,ierr)
!            ALLOCATE(ThreeIndInts03Temp(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!            CALL LogMemAlloc('ThreeIndInts03Temp',NoOrbs**4,8,this_routine,ThreeIndInts03TempTag,ierr)
!            ALLOCATE(ThreeIndInts04Temp(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
!            CALL LogMemAlloc('ThreeIndInts04Temp',NoOrbs**4,8,this_routine,ThreeIndInts04TempTag,ierr)
        
            ! Partially transformed combined arrays.
            ALLOCATE(TwoIndInts02(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('TwoIndInts02',NoOrbs**4,8,this_routine,TwoIndInts02Tag,ierr)
            ALLOCATE(ThreeIndInts01(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ThreeIndInts01',NoOrbs**4,8,this_routine,ThreeIndInts01Tag,ierr)
            ALLOCATE(ThreeIndInts03(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ThreeIndInts03',NoOrbs**4,8,this_routine,ThreeIndInts03Tag,ierr)
            ALLOCATE(ThreeIndInts04(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ThreeIndInts04',NoOrbs**4,8,this_routine,ThreeIndInts04Tag,ierr)
        ENDIF

        ! Allocate according to orthonormalisation method being used.
        IF(tLagrange) THEN
            ALLOCATE(Lambdas(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('Lambdas',NoOrbs**2,8,this_routine,LambdasTag,ierr)
            ALLOCATE(DerivLambda(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('DerivLambda',NoOrbs**2,8,this_routine,DerivLambdaTag,ierr)
            Lambdas(:,:)=0.0_dp
            DerivLambda(:,:)=0.0_dp
        ENDIF

        IF(tShake) THEN
            ALLOCATE(ShakeLambda(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('ShakeLambda',TotNoConstraints,8,this_routine,ShakeLambdaTag,ierr)
            ShakeLambda(:)=0.0_dp                     
            ALLOCATE(ShakeLambdaNew(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('ShakeLambdaNew',TotNoConstraints,8,this_routine,ShakeLambdaNewTag,ierr)
            ShakeLambdaNew(:)=0.0_dp                     
            ALLOCATE(Constraint(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('Constraint',TotNoConstraints,8,this_routine,ConstraintTag,ierr)
            ALLOCATE(ConstraintCor(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('ConstraintCor',TotNoConstraints,8,this_routine,ConstraintCorTag,ierr)
            ALLOCATE(DerivConstrT1(NoOrbs,NoOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT1',NoOrbs*TotNoConstraints*NoOrbs,8,this_routine,DerivConstrT1Tag,ierr)
            DerivConstrT1(:,:,:)=0.0_dp
            ALLOCATE(DerivConstrT2(NoOrbs,NoOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT2',NoOrbs*TotNoConstraints*NoOrbs,8,this_routine,DerivConstrT2Tag,ierr)
            ALLOCATE(ForceCorrect(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ForceCorrect',NoOrbs**2,8,this_routine,ForceCorrectTag,ierr)
            ALLOCATE(Correction(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('Correction',NoOrbs**2,8,this_routine,CorrectionTag,ierr)
            IF(tShakeApprox) THEN
                ALLOCATE(DerivConstrT1T2Diag(TotNoConstraints),stat=ierr)
                CALL LogMemAlloc('DerivConstrT1T2Diag',TotNoConstraints,8,this_routine,DerivConstrT1T2DiagTag,ierr)
                DerivConstrT1T2Diag(:)=0.0_dp
            ELSE
                ALLOCATE(DerivConstrT1T2(TotNoConstraints,TotNoConstraints),stat=ierr)
                CALL LogMemAlloc('DerivConstrT1T2',TotNoConstraints**2,8,this_routine,DerivConstrT1T2Tag,ierr)
            ENDIF
        ENDIF               

        ! Indexing arrays.
        ALLOCATE(SymLabelList2_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList2_rot',NoOrbs,4,this_routine,SymLabelList2_rotTag,ierr)
        SymLabelList2_rot(:)=0                     
        ALLOCATE(SymLabelList3_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3_rot',NoOrbs,4,this_routine,SymLabelList3_rotTag,ierr)
        SymLabelList3_rot(:)=0                     
 
        ALLOCATE(SymLabelListInv_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelListInv_rot',NoOrbs,4,this_routine,SymLabelListInv_rotTag,ierr)
        SymLabelListInv_rot(:)=0                     
   
        ALLOCATE(Lab(2,TotNoConstraints),stat=ierr)
        CALL LogMemAlloc('Lab',2*TotNoConstraints,4,this_routine,LabTag,ierr)
        Lab(:,:)=0                     


! Do any initial calculations, and set up starting values for arrays used in rotation.        
        CALL InitRotCalc()


! Write out the headings for the results file.        
        transform_unit = get_free_unit()
        OPEN(transform_unit,FILE='Transform',STATUS='unknown')
        IF(tLagrange) THEN
            WRITE(transform_unit,"(A12,11A18)") "# Iteration","2.PotEnergy","3.PEInts","4.PEOrtho","5.Force","6.ForceInts", &
            "7.OrthoForce","8.Sum<ij|kl>^2",&
                        &"9.OrthoNormCondition","10.DistMovedbyCs","11.DistMovedByLs","12.LambdaMag"
            WRITE(6,"(A12,11A19)") "Iteration","2.PotEnergy","3.PEInts","4.PEOrtho","5.Force","6.ForceInts","7.OrthoForce", &
                    "8.Sum<ij|kl>^2",&
                    &"9.OrthoNormCondition","10.DistMovedbyCs","11.DistMovedbyLs","12.LambdaMag"
        ELSEIF(tERLocalization.and.tHijSqrdMin) THEN
            WRITE(transform_unit,"(A12,7A24)") "# Iteration","2.ERPotEnergy","3.HijSqrdPotEnergy","4.PotEnergy","5.Force", &
                "6.Totalcorrforce","7.OrthoNormCondition","8.DistMovedbyCs"
            WRITE(6,"(A12,7A24)") "# Iteration","2.ERPotEnergy","3.HijSqrdPotEnergy","4.PotEnergy","5.Force",              &
                "6.Totalcorrforce","7.OrthoNormCondition","8.DistMovedbyCs"
        ELSEIF(tERLocalization) THEN
            WRITE(transform_unit,"(A12,5A24)") "# Iteration","2.Sum_i<ii|ii>","3.Force","4.TotCorrForce",   &
                "5.OrthoNormCondition","6.DistMovedbyCs"
            WRITE(6,"(A12,5A24)") "Iteration","2.Sum_i<ii|ii>","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        ELSE
            WRITE(transform_unit,"(A12,5A24)") "# Iteration","2.PotEnergy","3.Force","4.Totalcorrforce",    &
                "5.OrthoNormCondition","6.DistMovedbyCs"
            WRITE(6,"(A12,5A24)") "Iteration","2.PotEnergy","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        ENDIF


    END SUBROUTINE InitLocalOrbs


    SUBROUTINE InitRotCalc()
! Sets up the initial arrays to be used in the orbital rotation.    
        CHARACTER(len=*) , PARAMETER :: this_routine='InitRotCalc'
        real(dp) :: RAN2
        INTEGER :: i,j,Const,iseed=-8,MinRot,MaxRot


        CALL InitSymmArrays()
! Creates an indexing system for each of the cases with symmetry on/off, and mixing all orbitals or separating
! the occuppied from virtual.
! The arrays used in this routine are labelled with a 2 (SymLabelList2_rot and SymLabelCount2), so as to not
! mess up the spawing/FCI calcs.
        do i=1,NoOrbs
            SymLabelList3_rot(i)=SymLabelList2_rot(i)
        enddo

! Set up constraint labels.  Constraint l is the dot product of i.j.
        Const=0
        do i=1,NoOrbs
            do j=i,NoOrbs
                Const=Const+1
                Lab(1,Const)=i
                Lab(2,Const)=j
            enddo
        enddo

! Just a check that the number of constraints labeled is the same as that calculated above.
        WRITE(6,*) 'Total number of constraints = ',TotNoConstraints
        IF(Const.ne.TotNoConstraints) THEN
            CALL Stop_all(this_routine,'ERROR in the number of constraints calculated.  lmax does not equal TotNoConstraints')
        ENDIF
 
! Zero/initialise the arrays
! In the case where symmetry is kept, the starting transformation matrix is just the identity.  Starting with a symmetric system
! means the symmetry is never broken.
! When we are breaking the symmetry, the starting transformation matrix is completely random, (and then orthonormalised).
! The ordering of the orbitals in CoeffT1 follow the ordering in SymLabelList2_rot.

        CoeffT1(:,:)=0.0_dp
        IF(tRotateOccOnly) THEN
            MinRot=1
            MaxRot=NoOcc
        ELSEIF(tRotateVirtOnly) THEN
            MinRot=NoOcc+1
            MaxRot=NoOrbs
        ELSE
            MinRot=1
            MaxRot=NoOrbs
        ENDIF
        do i=1,NoOrbs
            CoeffT1(i,i)=1.0_dp
        enddo
        ! If the symmetry is kept on, start with the symmetric identity matrix of coefficients, and it will be maintained.

! When only the occupied or virtual orbitals are rotated, the non rotated orbitals are still included in the transformation matrix,
! but start as the identity, and remain that way throughout.
        IF(lNoSymmetry) THEN
            do i=MinRot,MaxRot
                do j=MinRot,MaxRot
                    CoeffT1(j,i)=RAN2(iseed)*(1E-02)
                enddo
            enddo
! This bit is for when symmetry is being kept, but we want to start with random coefficients, rather than the HF orbitals.
! It is often used to check the same final potential energy is reached in both cases.
!        ELSE 
!            do w=MinOccVirt,MaxOccVirt
!                IF(w.eq.1) THEN
!                    SymMin=1
!                    MinMZ=1
!                    IF(tSeparateOccVirt) THEN
!                        MaxMZ=NoOcc
!                    ELSE
!                        MaxMZ=NoOrbs
!                    ENDIF
!                ELSE
!                    SymMin=9
!                    MinMZ=NoOcc+1
!                    MaxMZ=NoOrbs
!                ENDIF
!           
!                do m=MinMZ,MaxMZ
!                    SymM=INT(G1(SymLabelList2_rot(m)*2)%sym%S)
!                    do z=SymLabelCounts2_rot(1,SymM+SymMin),(SymLabelCounts2_rot(1,SymM+SymMin)+SymLabelCounts2_rot(2,SymM+SymMin)-1)
!                        CoeffT1(z,m)=RAN2(iseed)*(1E-01)
!                    enddo
!                enddo
!            enddo
!            WRITE(6,*) 'Starting from a randomised coefficient matrix (keeping symmetry)'
        ENDIF

! Ensures transformation matrix elements between the occupied and virtual orbitals are 0 (should be the case anyway though).
        IF(tSeparateOccVirt) CALL ZeroOccVirtElements(CoeffT1)

! Orthonormalise starting matrix.        
        CALL GRAMSCHMIDT(CoeffT1,NoOrbs)


!        WRITE(6,*) 'coefft1'
!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                WRITE(6,'(F15.10)',advance='no') CoeffT1(i,j)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        stop


! A UMATTemp is created (from the UMAT read in from the FCIDUMP) using the rotate orbs indexing system.
! i.e. in UMatTemp(1,1,1,1), the orbital involved is the first in SymLabelList2_rot.
! Doing this now, rather than using UMatInd in each transform2elint routine proved a lot faster.
        DerivCoeff(:,:)=0.0_dp
        UMATTemp01(:,:,:,:)=0.0_dp
        IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat) &
            .and.(.not.tFindCINatOrbs).and.(.not.tUseHFOrbs))&
            &.or.(tERLocalization.and.tStoreSpinOrbs)) UMATTemp02(:,:,:,:)=0.0_dp

        CALL CopyAcrossUMAT()

        CALL TestOrthonormality()
    

!        WRITE(6,*) 'i,j,TMAT2D'
!        do i=1,NoOrbs*2
!            do j=1,NoOrbs*2
!                WRITE(6,'(F20.10)',advance='no') REAL(TMAT2D(i,j),8)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'i,j,TMAT2DTemp'
!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                WRITE(6,'(F20.10)',advance='no') TMAT2DTemp(i,j)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        stop


! Find the lower and upper bounds for each processor when running over a=1,NoOrbs and then b=1,a.
! This is so that each processor is doing roughly the same amount of work.
!        do i=1,NoOrbs
!            SumNoOrbs=SumNoOrbs+i
!        enddo
!        NoOrbsPerProc=SumNoOrbs/nProcessors
!        j=1
!        k=0
!        do i=1,iProcIndex+1
!            NoOrbsRem=NoOrbsPerProc
!            do while (NoOrbsRem.gt.0)
!                k=k+1
!                NoOrbsRem=NoOrbsRem-k
!            enddo
!            IF(iProcIndex.eq.(i-1)) THEN
!                LowBound02=j
!                HighBound02=k-1
!            ENDIF
!            j=k
!        enddo
!        IF(iProcIndex.eq.(nProcessors-1)) HighBound02=NoOrbs



! With UMAT with the correct indexing and the starting coefficient, find the partially transformed 
! four index integrals (and hence the initial potential energy), and then the initial force.

        IF(tERLocalization.and.(.not.tStoreSpinOrbs)) THEN
            CALL Transform2ElIntsERlocal()
        ELSE
            CALL Transform2ElInts()
        ENDIF
        CALL FindTheForce()

        IF(tPrintInts) THEN
            tInitIntValues=.true.
            CALL PrintIntegrals()
            tInitIntValues=.false.
            ! This sets the initial values for the integral sums being printed.
            ! Values printed are then relative to these initial sums, per integral.
        ENDIF


    END SUBROUTINE InitRotCalc


    SUBROUTINE CopyAcrossUMAT()
        INTEGER :: a,b,g,d,i,j,k,l
        real(dp) :: s,t

        IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
        &.or.(tERLocalization.and.tStoreSpinOrbs)) TMAT2DTemp(:,:)=0.0_dp


! These loops can be sped up with spatial symmetry and pairwise permutation symmetry if needed.
        do a=1,NoOrbs
            i=SymLabelList2_rot(a)                  ! The spin orbital we are looking for.
            do g=1,a
                j=SymLabelList2_rot(g)
                IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
                &.or.(tERLocalization.and.tStoreSpinOrbs)) THEN
                    IF(tStoreSpinOrbs) THEN
                        s=REAL(TMAT2D(i,j),dp)
                        TMAT2DTemp(a,g)=s
                        TMAT2DTemp(g,a)=s
                    ELSE
                        s=REAL(TMAT2D(2*i,2*j),dp)
                        TMAT2DTemp(a,g)=s
                        TMAT2DTemp(g,a)=s
                    ENDIF
                ENDIF


                do b=1,NoOrbs
                    k=SymLabelList2_rot(b)
                    do d=1,b
                        l=SymLabelList2_rot(d)
                        t=REAL(UMAT(UMatInd(i,k,j,l,0,0)),dp)
                        UMATTemp01(a,g,b,d)=t                   !a,g,d,b chosen to make 'transform2elint' steps more efficient
                        UMATTemp01(g,a,b,d)=t
                        UMATTemp01(a,g,d,b)=t
                        UMATTemp01(g,a,d,b)=t
                        IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
                        &.or.(tERLocalization.and.tStoreSpinOrbs)) THEN
                            UMATTemp02(d,b,a,g)=t                   !d,b,a,g order also chosen to speed up the transformation.
                            UMATTemp02(d,b,g,a)=t
                            UMATTemp02(b,d,a,g)=t
                            UMATTemp02(b,d,g,a)=t
                        ENDIF
                    enddo
                enddo
            enddo
        enddo

!        do a=1,NoOrbs
!            do g=1,a
!                do b=1,NoOrbs
!                    do d=1,b
!                        WRITE(6,'(4I3,F20.10,4I3,F20.10)') a,b,g,d,UMATTemp01(a,g,b,d),a,b,d,g,UMATTemp01(a,d,b,g)
!                    enddo
!                enddo
!            enddo
!        enddo
!        stop
!        do a=1,NoOrbs
!            do g=1,a
!                WRITE(6,'(2I3,F20.10)') a,g,TMAT2DTemp(a,g)
!            enddo
!        enddo
!        stop


    END SUBROUTINE CopyAcrossUMAT        


    SUBROUTINE WriteStats()

        IF(tLagrange) THEN
            WRITE(6,"(I12,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts, &
                OrthoNorm,DistCs,DistLs,LambdaMag
            WRITE(transform_unit,"(I12,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce, &
                TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
        ELSEIF(tERLocalization.and.tHijSqrdMin) THEN
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(6,"(I12,7F24.10)") Iteration,ERPotEnergy,HijSqrdPotEnergy,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
                WRITE(transform_unit,"(I12,7F24.10)") Iteration,ERPotEnergy,HijSqrdPotEnergy,PotEnergy,Force, &
                    TotCorrectedForce,OrthoNorm,DistCs
            ENDIF
        ELSE
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(6,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
                WRITE(transform_unit,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
            ENDIF
        ENDIF
        CALL neci_flush(6)
        CALL neci_flush(transform_unit)

! after writing out stats, test for SOFTEXIT.
        if (test_SOFTEXIT()) then
            WRITE(6,*) 'SOFTEXIT detected, finalizing new orbitals.'
            tNotConverged=.false.
        endif


    END SUBROUTINE WriteStats
         



    SUBROUTINE InitSymmArrays()
! This routine creates indexing arrays for the cases with symmetry on/off, and either mixing all orbitals or 
! separating the occupied and virtuals.
! The arrays used specific to the orbital rotation are named with a 2. 

! The arrays produced are as follows...
! SymLabelList2_rot(NoOrbs) contains the spatial orbitals, ordered in groups of increasing symmetry label.
! - when the orbitals are being separated, the first NoOcc of SymLabelList2_rot are the occupied, and the rest are virtual.
! - essentially this array relates the orbital labelling used in the orbital rotation (1,2,3 according to the order 
! - in SymLabelList2_rot) to the labels used in arrays being fed in/out of this routine (UMAT etc).

! SymLabelCounts2_rot(1:Sym) is the index in SymLabelList where the symmetry block S starts
! SymLabelCounts2_rot(2:Sym) is the number of orbitals in symmetry block S.
! E.g. if symmetry S starts at index 2 and has 3 orbitals.
! SymLabelList2_rot(2)->SymLabelList2_rot(4) will give the indexes of these orbitals.
        use sym_mod, only: GenSymStatePairs
        INTEGER :: j,i,ierr
        CHARACTER(len=*) , PARAMETER :: this_routine='InitSymmArrays'

        IF(.not.tSeparateOccVirt) THEN
            SymLabelCounts(:,:)=0
            SymLabelList(:)=0
            IF(tStoreSpinOrbs) CALL Stop_All(this_routine,"There may be a problem with GENSymStatePairs when using spin orbitals.")
            CALL GENSymStatePairs(SpatOrbs,.false.)
        ENDIF
! Sets up the SymLabelList and SymLabelCounts arrays used in the spawing etc. (When the rotate
! orbs routine is called, this has not been done yet).
! If the symmetry is on, and all orbitals are being mixed, this will end up being the same as SymLabelList2_rot.


        IF(tSeparateOccVirt) THEN
            MinOccVirt=1
            MaxOccVirt=2
            IF(tRotateOccOnly) THEN
                MaxOccVirt=1
            ELSEIF(tRotateVirtOnly) THEN
                MinOccVirt=2
            ENDIF
            CALL InitOrbitalSeparation()
            ! rewrite all the symmetry lists to account for the separation and have simple option if
            ! symmetry is off.
        ELSE
            MinOccVirt=1
            MaxOccVirt=1
            ALLOCATE(SymLabelCounts2_rot(2,8),stat=ierr)
            CALL LogMemAlloc('SymLabelCounts2_rot',2*8,4,this_routine,SymLabelCounts2_rotTag,ierr)
            SymLabelCounts2_rot(:,:)=0                     
            do i=1,SpatOrbs   
                IF(tStoreSpinOrbs) THEN
                    SymLabelList2_rot(2*i)=2*SymLabelList(i)
                    SymLabelList2_rot(2*i-1)=(2*SymLabelList(i))-1
                ELSE
                    SymLabelList2_rot(i)=SymLabelList(i)
                ENDIF
            enddo
            IF(lNoSymmetry) THEN
                SymLabelCounts2_rot(1,1)=1
                SymLabelCounts2_rot(2,1)=NoOrbs
            ELSE
                do j=1,8
                    IF(tStoreSpinOrbs) THEN
                        SymLabelCounts2_rot(1,j)=(2*SymLabelCounts(1,j))-1
                        SymLabelCounts2_rot(2,j)=2*SymLabelCounts(2,j)
                    ELSE
                        do i=1,2
                            SymLabelCounts2_rot(i,j)=SymLabelCounts(i,j)
                        enddo
                    ENDIF
                enddo
            ENDIF
        ENDIF

        do i=1,NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i))=i
        enddo
        

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,16
!            WRITE(6,*) SymLabelCounts2_rot(1,i),SymLabelCounts2_rot(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
!        do i=1,NoOrbs
!            IF(tStoreSpinOrbs) THEN
!                WRITE(6,*) i,SymLabelList2_rot(i),INT(G1(SymLabelList2_rot(i))%sym%S)
!            ELSE
!                WRITE(6,*) i,SymLabelList2_rot(i),INT(G1(SymLabelList2_rot(i)*2)%sym%S)
!            ENDIF
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,NoOrbs
!            WRITE(6,*) SymLabelList2_rot(i),SymLabelListInv_rot(i)
!        enddo
!        CALL neci_flush(6)
!        CALL Stop_All('InitSymmArrays','Checking orbital labels.')


    ENDSUBROUTINE InitSymmArrays



    SUBROUTINE EquateDiagFock()
        INTEGER :: irr,NumInSym,Orbi,Orbj,w,i,j,k,ConjInd,OrbjConj
        real(dp) :: Angle,AngleConj,Check,Norm

        CoeffT1(:,:)=0.0_dp
!        MaxOccVirt=1
!        WRITE(6,*) MaxOccVirt,"***"

        do w=MinOccVirt,MaxOccVirt
!Do virtual and occupied orbitals seperately

            do irr=1,8
!Loop over irreps

                NumInSym=SymLabelCounts2_rot(2,(w-1)*8+irr)
!                WRITE(6,*) "NumInSym= ",NumInSym,irr-1

                do j=1,NumInSym
!Loop over the j-orthogonal vectors to create in this symmetry block

                    Orbj=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+j)

!See if this vector has already been done.
                    Check=0.0_dp
                    do i=1,NoOrbs
                        Check=Check+CoeffT1(i,Orbj)
                    enddo
                    IF(Check.ne.0.0_dp) THEN
!This vector is a conjugate pair of another vector and has already been worked out...
                        CYCLE
                    ENDIF

!Find out if we this vector will be complex. It will be real if j=N or j=N/2
                    IF(j.eq.NumInSym) THEN
!The vector will be the normalized 1,1,1 vector.

                        do i=1,NumInSym
                            Orbi=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)
                            CoeffT1(Orbi,Orbj)=1/SQRT(REAL(NumInSym,dp))
                        enddo

                    ELSEIF((mod(NumInSym,2).eq.0).and.(j.eq.(NumInSym/2))) THEN

                        do i=1,NumInSym
                            Orbi=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)
                            IF(mod(i,2).eq.1) THEN
                                CoeffT1(Orbi,Orbj)=-1/SQRT(REAL(NumInSym,dp))
                            ELSE
                                CoeffT1(Orbi,Orbj)=1/SQRT(REAL(NumInSym,dp))
                            ENDIF
                        enddo

                    ELSE
!Vector is complex - find its conjugate vector - do these at the same time.
                        ConjInd=NumInSym-j
                        OrbjConj=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+ConjInd)
                        
                        do i=1,NumInSym
                            
                            Orbi=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)

                            Angle=REAL(i*j*2,dp)*PI/REAL(NumInSym,dp)
                            AngleConj=REAL(i*ConjInd*2,dp)*PI/REAL(NumInSym,dp)

                            CoeffT1(Orbi,Orbj)=(1/SQRT(REAL(2*NumInSym,dp)))*(COS(Angle)+COS(AngleConj))
                            CoeffT1(Orbi,OrbjConj)=(1/SQRT(REAL(2*NumInSym,dp)))*(SIN(Angle)-SIN(AngleConj))

                        enddo

                    ENDIF

                enddo


!                    do i=1,NumInSym
!
!                        Orbi=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+i)
!                        WRITE(6,*) "Sym= ",irr-1, Orbj, Orbi
!!Coefficients are going to be C_jk = exp^(i j*k 2Pi/N), i.e roots of unity
!
!!                        IF(CoeffT1(i,j).ne.0.0_dp) CYCLE
!
!                        Prod=i*j
!                        IF(mod(Prod,NumInSym).eq.0) THEN
!! i*j = N or 2N, 3N, ...
!                            CoeffT1(Orbi,Orbj)=1.0_dp/SQRT(REAL(NumInSym,8))
!                            CoeffT1(Orbj,Orbi)=1.0_dp/SQRT(REAL(NumInSym,8))
!
!                        ELSEIF((mod(Prod*2,NumInSym).eq.0).and.(mod((2*Prod)/NumInSym,2).eq.1)) THEN
!! i*j = N/2 or 3N/2, 5N/2, ...
!                            CoeffT1(Orbi,Orbj)=-1.0_dp/SQRT(REAL(NumInSym,8))
!                            CoeffT1(Orbj,Orbi)=-1.0_dp/SQRT(REAL(NumInSym,8))
!
!                        ELSE
!!Here, the values will be complex. Therefore we need to take symmetric and antisymmetric combinations.
!! Symmetric is: C_jk        = 1/SQRT(2)  [ Phi_jk + Phi_(N-j)k ] = 1/SQRT(2)*2COS(2 j*k PI/N)
!! AntiSymm is:  C_(N-j)k    = 1/SQRT(-2) [ Phi_jk - Phi_(N-j)k ] = 1/SQRT(2)*2SIN(2 j*k PI/N)
!
!                            IF(i.gt.(NumInSym-i)) THEN
!
!!                                ConjInd=mod((NumInSym-mod(i+j-1,NumInSym))-(j-1),NumInSym)
!                                ConjInd=NumInSym-i
!                                ConjOrb=SymLabelList2_rot(SymLabelCounts2_rot(1,(w-1)*8+irr)-1+ConjInd)
!
!                                Angle=REAL(Prod*2,8)*3.141592654/REAL(NumInSym,8)
!                                CoeffT1(Orbi,Orbj)=1.0_dp/SQRT(REAL(2*NumInSym,8))*2.0_dp*COS(Angle)
!                                CoeffT1(ConjOrb,Orbj)=1.0_dp/SQRT(REAL(2*NumInSym,8))*2.0_dp*SIN(Angle)
!                                WRITE(6,*) "Ind = ",i," J= ",j, "ConjInd = ",ConjInd, " N = ",NumInSym, " Orbj = ",Orbj
!                                WRITE(6,*) "Angle = ",Angle, " CoeffT1(Orbi,Orbj) = ", CoeffT1(Orbi,Orbj), "CoeffT1(ConjOrb,Orbj) = ", CoeffT1(ConjOrb,Orbj)
!
!                            ENDIF
!
!                        ENDIF
!
!                    enddo
!                enddo
            enddo
        enddo

        do j=1,NoOrbs
            Norm=0.0_dp
            do i=1,NoOrbs
                Norm=Norm+(CoeffT1(i,j)**2)
            enddo
            IF(Norm.eq.0.0_dp) THEN
                CoeffT1(j,j)=1.0_dp
            ENDIF
        enddo

        do j=1,NoOrbs
            do i=1,NoOrbs
                WRITE(6,"(G13.5)",advance='no') CoeffT1(j,i)
            enddo
            WRITE(6,*) ""
        enddo

!Check normalization
        do j=1,NoOrbs
            Norm=0.0_dp
            do i=1,NoOrbs
                Norm=Norm+(CoeffT1(i,j)**2)
            enddo
            IF(abs(Norm-1.0_dp).gt.1.0e-7_dp_dp) THEN
                CALL Stop_All("EquateDiagFock","Rotation Coefficients not normalized")
            ENDIF
        enddo


!Check orthogonality
        do j=1,NoOrbs
            do i=1,NoOrbs
                IF(i.eq.j) CYCLE
                Norm=0.0_dp
                do k=1,NoOrbs
                    Norm=Norm+(CoeffT1(k,j)*CoeffT1(k,i))
                enddo
                IF(abs(Norm).gt.1.0e-7_dp_dp) THEN
                    WRITE(6,*) "COLUMNS: ",j,i
                    CALL Stop_All("EquateDiagFock","RotationCoefficients not orthogonal")
                ENDIF
            enddo
        enddo

    END SUBROUTINE EquateDiagFock


    SUBROUTINE InitOrbitalSeparation()
! This subroutine is called if the SEPARATEOCCVIRT keyword is present in the input, it sets up SymLabelList2_rot so that the first 
! NoOcc orbitals are the HF occupied, and the rest the virtual.  Within this separation, orbitals are ordered in symmetry 
! groups. 
! This means that two iterations of the rotate orbs routine will be performed, the first treats the occupied orbitals and the second
! the virtual.
        INTEGER :: i,j,ierr,SymCurr,Symi
        INTEGER(TagIntType) :: SymVirtOrbsTag,SymOccOrbsTag
        integer :: lo, hi
        INTEGER , ALLOCATABLE :: SymVirtOrbs(:),SymOccOrbs(:)
        CHARACTER(len=*) , PARAMETER :: this_routine='InitOrbitalSeparation'


        ALLOCATE(SymLabelCounts2_rot(2,16),stat=ierr)
        CALL LogMemAlloc('SymLabelCounts2_rot',2*16,4,this_routine,SymLabelCounts2_rotTag,ierr)
        SymLabelCounts2_rot(:,:)=0
        ! first 8 refer to the occupied, and the second to the virtual.

        ALLOCATE(LabVirtOrbs(NoOrbs-NoOcc),stat=ierr)
        CALL LogMemAlloc('LabVirtOrbs',(NoOrbs-NoOcc),4,this_routine,LabVirtOrbsTag,ierr)
        LabVirtOrbs(:)=0
        ALLOCATE(LabOccOrbs(NoOcc),stat=ierr)
        CALL LogMemAlloc('LabOccOrbs',(NoOcc),4,this_routine,LabOccOrbsTag,ierr)
        LabOccOrbs(:)=0
        ALLOCATE(SymVirtOrbs(NoOrbs-NoOcc),stat=ierr)
        CALL LogMemAlloc('SymVirtOrbs',(NoOrbs-NoOcc),4,this_routine,SymVirtOrbsTag,ierr)
        SymVirtOrbs(:)=0
        ALLOCATE(SymOccOrbs(NoOcc),stat=ierr)
        CALL LogMemAlloc('SymOccOrbs',(NoOcc),4,this_routine,SymOccOrbsTag,ierr)
        SymOccOrbs(:)=0


! First fill SymLabelList2_rot.

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index with the second lowest energy.

!        do i=1,nBasis
!            WRITE(6,*) BRR(i)
!        enddo

! this picks out the NoOcc lowest energy orbitals from BRR as these will be the occupied.
! these are then ordered according to symmetry, and the same done to the virtual.
        do i=1,NoOcc
            IF(tStoreSpinOrbs) THEN
                LabOccOrbs(i)=BRR(i)
                SymOccOrbs(i)=INT(G1(LabOccOrbs(i))%sym%S)
            ELSE
                LabOccOrbs(i)=(BRR(2*i))/2
                SymOccOrbs(i)=INT(G1(LabOccOrbs(i)*2)%sym%S)
            ENDIF
        enddo
        
        call sort (SymOccOrbs, LabOccOrbs)
        ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in terms of symmetry). 

        do i=1,NoOrbs-NoOcc
            IF(tStoreSpinOrbs) THEN
                LabVirtOrbs(i)=BRR(i+NEl)
                SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i))%sym%S)
            ELSE
                LabVirtOrbs(i)=(BRR((2*i)+NEl))/2
                SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i)*2)%sym%S)
            ENDIF
        enddo
        
        call sort (SymVirtOrbs, LabVirtOrbs)

! SymLabelList2_rot is then filled with the symmetry ordered occupied then virtual arrays.        
        do i=1,NoOcc
            SymLabelList2_rot(i)=LabOccOrbs(i)
        enddo
        j=0
        do i=NoOcc+1,NoOrbs
            j=j+1
            SymLabelList2_rot(i)=LabVirtOrbs(j)
        enddo

!        WRITE(6,*) 'symlabellist'
!        do i=1,NoOrbs
!            WRITE(6,'(2I4)') SymLabelList2_rot(i),INT(G1(SymLabelList2_rot(i)*2)%sym%S)
!        enddo
!        stop

!************
! Second fill SymLabelCounts2_rot.
! - the first 8 places of SymLabelCounts2_rot(1,:) and SymLabelCounts2_rot(2,:) refer to the occupied orbitals 
! - and the second 8 to the virtuals.

        IF(lNoSymmetry) THEN
            ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
            SymLabelCounts2_rot(1,1)=1
            SymLabelCounts2_rot(1,9)=NoOcc+1
            SymLabelCounts2_rot(2,1)=NoOcc
            SymLabelCounts2_rot(2,9)=NoOrbs-NoOcc
        ELSE 
            ! otherwise we run through the occupied orbitals, counting the number with each symmetry
            ! and noting where in SymLabelList2_rot each symmetry block starts.
            SymCurr=0
            SymLabelCounts2_rot(1,1)=1
            do i=1,NoOcc
                IF(tStoreSpinOrbs) THEN
                    Symi=INT(G1(SymLabelList2_rot(i))%sym%S)
                ELSE
                    Symi=INT(G1(SymLabelList2_rot(i)*2)%sym%S)
                ENDIF
                SymLabelCounts2_rot(2,(Symi+1))=SymLabelCounts2_rot(2,(Symi+1))+1
                IF(Symi.gt.SymCurr) THEN
                    SymLabelCounts2_rot(1,(Symi+1))=i
                    SymCurr=Symi
                ENDIF
            enddo
            ! the same is then done for the virtuals.
            SymCurr=0
            SymLabelCounts2_rot(1,9)=NoOcc+1
            do i=NoOcc+1,NoOrbs
                IF(tStoreSpinOrbs) THEN
                    Symi=INT(G1(SymLabelList2_rot(i))%sym%S)
                ELSE
                    Symi=INT(G1(SymLabelList2_rot(i)*2)%sym%S)
                ENDIF
                SymLabelCounts2_rot(2,(Symi+9))=SymLabelCounts2_rot(2,(Symi+9))+1
                IF(Symi.gt.SymCurr) THEN
                    SymLabelCounts2_rot(1,(Symi+9))=i
                    SymCurr=Symi
                ENDIF
            enddo
        ENDIF

        ! Go through each symmetry group, making sure the orbital pairs are ordered lowest to highest.
        do i=1,16
            IF(SymLabelCounts2_rot(2,i).ne.0) THEN
                lo = SymLabelCounts2_rot(1,i)
                hi = lo + SymLabelCounts2_rot(2,i) - 1
                call sort (SymLabelList2_rot(lo:hi))
            ENDIF
        enddo


! Deallocate the arrays just used in this routine.
        DEALLOCATE(LabOccOrbs)
        CALL LogMemDealloc(this_routine,LabOccOrbsTag)
        DEALLOCATE(LabVirtOrbs)
        CALL LogMemDealloc(this_routine,LabVirtOrbsTag)
        DEALLOCATE(SymOccOrbs)
        CALL LogMemDealloc(this_routine,SymOccOrbsTag)
        DEALLOCATE(SymVirtOrbs)
        CALL LogMemDealloc(this_routine,SymVirtOrbsTag)


    ENDSUBROUTINE InitOrbitalSeparation




    SUBROUTINE Diagonalizehij()
! This routine takes the original <i|h|j> matrix and diagonalises it.  The resulting coefficients from this process 
! are then the rotation coefficients to be applied to the four index integrals etc.
! This eliminates the <i|h|j> elements from the single excitations, and leaves only coulomb and exchange terms.
! In order to maintain the same HF energy, only the virtual elements are diagonalised, within symmetry blocks.
        INTEGER :: i,j,Sym,ierr,NoSymBlock,WorkSize,WorkCheck,SymStartInd
        INTEGER(TagIntType) WorkTag,DiagTMAT2DBlockTag,TMAT2DSymBlockTag
        REAL , ALLOCATABLE :: TMAT2DSymBlock(:,:),DiagTMAT2DBlock(:),Work(:)
        CHARACTER(len=*) , PARAMETER :: this_routine='Diagonalizehij'
 

        WRITE(6,*) 'The original coefficient matrix'
        do i=1,NoOrbs
            do j=1,NoOrbs
                WRITE(6,'(F20.10)',advance='no') CoeffT1(j,i)
            enddo
            WRITE(6,*) ''
        enddo

        WRITE(6,*) 'The original TMAT2D matrix'
        do i=1,NoOrbs
            do j=1,NoOrbs
                WRITE(6,'(F20.10)',advance='no') TMAT2DTemp(j,i)
            enddo
            WRITE(6,*) ''
        enddo
        TMAT2DRot(:,:)=0.0_dp
        DiagTMAT2Dfull(:)=0.0_dp

! The final <i|h|j> matrix will be TMAT2DRot, just copy accross the occupied elements, as these will not be changed.
!        do j=1,NoOcc
!            do i=1,NoOcc
!                TMAT2DRot(i,j)=TMAT2DTemp(i,j)
!            enddo
!        enddo
! Note, if decide to use the transform2elints to fill this, don't need this step.

! Now need to pick out symmetry blocks, from the virtual orbitals and diagonalize them.

! Take first symmetry, (0) and find the number of virtual orbitals with this symmetry.  If this is greater than 1, 
! take the block, diagonlize it, and put it into TMAT2DRot.
        
        Sym=0
        WorkSize=-1
        do while (Sym.le.7)

            NoSymBlock=SymLabelCounts2_rot(2,Sym+9)

            SymStartInd=SymLabelCounts2_rot(1,Sym+9)-1
            ! This is one less than the index that the symmetry starts, so that when we run through i=1,..., we can
            ! start at SymStartInd+i.

            IF(NoSymBlock.gt.1) THEN
                ALLOCATE(TMAT2DSymBlock(NoSymBlock,NoSymBlock),stat=ierr)
                CALL LogMemAlloc('TMAT2DSymBlock',NoSymBlock**2,8,this_routine,TMAT2DSymBlockTag,ierr)
                ALLOCATE(DiagTMAT2DBlock(NoSymBlock),stat=ierr)
                CALL LogMemAlloc('DiagTMAT2DBlock',NoSymBlock,8,this_routine,DiagTMAT2DBlockTag,ierr)

                WorkCheck=3*NoSymBlock+1
                WorkSize=WorkCheck
                ALLOCATE(Work(WorkSize),stat=ierr)
                CALL LogMemAlloc('Work',WorkSize,8,this_routine,WorkTag,ierr)

                do j=1,NoSymBlock
                    do i=1,NoSymBlock
                        TMAT2DSymBlock(i,j)=TMAT2DTemp(SymStartInd+i,SymStartInd+j)
                    enddo
                enddo

                WRITE(6,*) '*****'
                WRITE(6,*) 'Symmetry ',Sym,' has ',NoSymBlock,' orbitals .'
                WRITE(6,*) 'The TMAT2D for this symmetry block is '
                do i=1,NoSymBlock
                    do j=1,NoSymBlock
                        WRITE(6,'(F20.10)',advance='no') TMAT2DSymBlock(j,i)
                    enddo
                    WRITE(6,*) ''
                enddo

                CALL DSYEV('V','U',NoSymBlock,TMAT2DSymBlock,NoSymBlock,DiagTMAT2Dblock,Work,WorkSize,ierr)
                ! TMAT2DSymBlock goes in as the original TMAT2DSymBlock, comes out as the eigenvectors (Coefficients).
                ! TMAT2DBlock comes out as the eigenvalues in ascending order.
                IF(ierr.ne.0) THEN
                    WRITE(6,*) 'Problem with symmetry, ',Sym,' of TMAT2D'
                    CALL neci_flush(6)
                    CALL Stop_All(this_routine,"Diagonalization of TMAT2DSymBlock failed...")
                ENDIF

                WRITE(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are ,'
                do i=1,NoSymBlock
                    WRITE(6,'(F20.10)',advance='no') DiagTMAT2Dblock(i)
                enddo
                WRITE(6,*) ''
                WRITE(6,*) 'These go from orbital ,',SymStartInd+1,' to ',SymStartInd+NoSymBlock
               
                do i=1,NoSymBlock
                    DiagTMAT2Dfull(SymStartInd+i-NoOcc)=DiagTMAT2DBlock(i)
                enddo

!                do i=1,NoSymBlock
!                    TMAT2DRot(SymStartInd+i,SymStartInd+i)=DiagTMAT2DBlock(i)
!                enddo
                ! CAREFUL if eigenvalues are put in ascending order, this may not be correct, with the labelling system.
                ! may be better to just take coefficients and transform TMAT2DRot in transform2elints.
                ! a check that comes out as diagonal is a check of this routine anyway.

                WRITE(6,*) 'The eigenvectors (coefficients) for symmtry block ',Sym
                do i=1,NoSymBlock
                    do j=1,NoSymBlock
                        WRITE(6,'(F20.10)',advance='no') TMAT2DSymBlock(j,i)
                    enddo
                    WRITE(6,*) ''
                enddo

             
                do j=1,NoSymBlock
                    do i=1,NoSymBlock
                        CoeffT1(SymStartInd+i,SymStartInd+j)=TMAT2DSymBlock(i,j)
                    enddo
                enddo
                ! Directly fill the coefficient matrix with the eigenvectors from the diagonalization.

                DEALLOCATE(Work)
                CALL LogMemDealloc(this_routine,WorkTag)

                DEALLOCATE(DiagTMAT2DBlock)
                CALL LogMemDealloc(this_routine,DiagTMAT2DBlockTag)

                DEALLOCATE(TMAT2DSymBlock)
                CALL LogMemDealloc(this_routine,TMAT2DSymBlockTag)
            ELSEIF(NoSymBlock.eq.1) THEN
                DiagTMAT2Dfull(SymStartInd+1-NoOcc)=TMAT2DTemp(SymStartInd+1,SymStartInd+1)
                WRITE(6,*) '*****'
                WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
                WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to DiagTMAT2Dfull'
            ENDIF

            Sym=Sym+1
        enddo
 
        WRITE(6,*) '*****'
        WRITE(6,*) 'The final coefficient matrix'
        do i=1,NoOrbs
            do j=1,NoOrbs
                WRITE(6,'(F20.10)',advance='no') CoeffT1(j,i)
            enddo
            WRITE(6,*) ''
        enddo

        WRITE(6,*) '*****'
        WRITE(6,*) 'The diagonal elements of TMAT2D'
        do i=1,(NoOrbs-NoOcc)
            WRITE(6,*) DiagTMAT2Dfull(i)
        enddo



    ENDSUBROUTINE Diagonalizehij




    SUBROUTINE ZeroOccVirtElements(Coeff)
! This routine sets all the elements of the coefficient matrix that connect occupied and virtual orbitals to 0.
! This ensures that only occupied mix with occupied and virtual mix with virtual.
        real(dp) :: Coeff(NoOrbs,NoOrbs)
        INTEGER :: i,j

        do i=1,NoOcc
            do j=NoOcc+1,NoOrbs
                Coeff(i,j)=0.0_dp
                Coeff(j,i)=0.0_dp
            enddo
        enddo


    ENDSUBROUTINE ZeroOccVirtElements




    SUBROUTINE FindNewOrbs()
           
        IF(tERLocalization.and.(.not.tStoreSpinOrbs)) THEN
            CALL Transform2ElIntsERlocal()
        ELSE
            CALL Transform2ElInts()     ! Find the partially (and completely) transformed 4 index integrals to be used in further calcs.
        ENDIF


!Find derivatives of the c and lambda matrices and print the sum of off-diagonal matrix elements.
        CALL FindTheForce()
        ! This finds the unconstrained force (unless the lagrange keyword is present).
      
!Update coefficents by moving them in direction of force. Print sum of squared changes in coefficients. 
        IF(tShake) THEN
            CALL ShakeConstraints()
            ! Find the force that moves the coefficients while keeping them orthonormal, and use it 
            ! to get these new coefficients.
        ELSE
            CALL UseTheForce()
            ! This can be either completely unconstrained, or have the lagrange constraints imposed.
        ENDIF
!The coefficients coefft1(a,m) are now those that have been shifted by the time step.

!Test these for orthonomaility and then convergence.
!If they do not meet the convergence criteria, they go back into the previous step to produce another set of coefficients.

        call set_timer(testorthoconver_time,30)        

        CALL TestOrthonormality()
!Force should go to zero as we end in minimum - test for this

        CALL TestForConvergence()

        call halt_timer(testorthoconver_time)


    END SUBROUTINE FindNewOrbs


    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.
    SUBROUTINE Transform2ElInts()
        INTEGER :: i,j,k,l,a,b,g,d
        real(dp) :: t,Temp4indints(NoRotOrbs,NoOrbs)
        real(dp) :: Temp4indints02(NoRotOrbs,NoRotOrbs)  

        
        CALL set_timer(Transform2ElInts_time,30)

!Zero arrays from previous transform

        TwoIndInts01(:,:,:,:)=0.0_dp
        FourIndInts(:,:,:,:)=0.0_dp

        IF(tNotConverged) THEN
!            TwoIndInts02Temp(:,:,:,:)=0.0_dp
            TwoIndInts02(:,:,:,:)=0.0_dp
!            ThreeIndInts01Temp(:,:,:,:)=0.0_dp
            ThreeIndInts01(:,:,:,:)=0.0_dp
!            ThreeIndInts02Temp(:,:,:,:)=0.0_dp
            ThreeIndInts02(:,:,:,:)=0.0_dp
!            ThreeIndInts03Temp(:,:,:,:)=0.0_dp
            ThreeIndInts03(:,:,:,:)=0.0_dp
!            ThreeIndInts04Temp(:,:,:,:)=0.0_dp
            ThreeIndInts04(:,:,:,:)=0.0_dp
!            FourIndInts02Temp(:,:,:,:)=0.0_dp
            FourIndInts02(:,:,:,:)=0.0_dp
        ENDIF

! ************
!Transform the 1 electron, 2 index integrals (<i|h|j>).
        IF(tNotConverged) THEN
            TMAT2DRot(:,:)=0.0_dp
            TMAT2DPartRot01(:,:)=0.0_dp
            TMAT2DPartRot02(:,:)=0.0_dp

!            WRITE(6,*) 'coefft1'
!            do i=1,NoOrbs
!                do j=1,NoOrbs
!                    WRITE(6,'(2F20.10)',advance='no') CoeffT1(i,j)
!                enddo
!                WRITE(6,*) ''
!            enddo

!            WRITE(6,*) 'tmat2dtemp'
!            do i=1,NoOrbs
!                do j=1,NoOrbs
!                    WRITE(6,'(2F20.10)',advance='no') TMAT2DTemp(i,j)
!                enddo
!                WRITE(6,*) ''
!            enddo
!            stop


            CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,    &
                TMAT2DTemp(:,:),NoOrbs,0.0,TMAT2DPartRot01(:,:),NoOrbs)
            ! get TMAT2DPartRot(i,a) out of this.

            CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,    &
                TMAT2DTemp(:,:),NoOrbs,0.0,TMAT2DPartRot02(:,:),NoOrbs)
            ! get TMAT2DPartRot(a,j) out of this.
     
            CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,    &
                TMAT2DPartRot01(:,:),NoOrbs,0.0,TMAT2DRot(:,:),NoOrbs)
            ! get TMAT2DRot(i,j) out of this.

        ENDIF


!        WRITE(6,*) 'TMAT2DRot in the virtuals'
!        do j=NoOcc+1,NoOrbs
!            do i=NoOcc+1,NoOrbs 
!                WRITE(6,'(F20.10)',advance='no') TMAT2DRot(i,j)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        stop


! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i,j,k,l,0,0)

!        do b=LowBound02,HighBound02
!            do d=1,b
        do b=1,NoOrbs
            do d=1,b
                Temp4indints(:,:)=0.0_dp
                CALL DGEMM('T','N',NoRotOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMatTemp01(:,:,d,b),NoOrbs, &
                    0.0,Temp4indints(:,:),NoRotOrbs)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, we need the transpose of this.

                Temp4indints02(:,:)=0.0_dp
                CALL DGEMM('T','T',NoRotOrbs,NoRotOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(:,:),NoRotOrbs, &
                    0.0,Temp4indints02(:,:),NoRotOrbs)
                ! Get Temp4indits02(i,k)

                do i=1,NoRotOrbs
                    do k=1,i
                        TwoIndInts01(d,b,k,i)=Temp4indints02(k,i)
                        TwoIndInts01(b,d,k,i)=Temp4indints02(k,i)
                        TwoIndInts01(d,b,i,k)=Temp4indints02(k,i)
                        TwoIndInts01(b,d,i,k)=Temp4indints02(k,i)

                    enddo
                enddo
            enddo
        enddo
        

! These calculations are unnecessary when this routine is calculated to finalize the new orbs.
        IF(tNotConverged) THEN
!            do g=LowBound02,HighBound02
            do g=1,NoOrbs                
                do a=1,g
                    Temp4indints(:,:)=0.0_dp
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMatTemp02(:,:,a,g),NoOrbs, &
                        0.0,Temp4indints(:,:),NoOrbs)
                    ! Temp4indints(l,b) comes out of here, so need to use transpose of this to transform the b elements.

                    Temp4indints02(:,:)=0.0_dp
                    CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(:,:),NoOrbs, &
                        0.0,Temp4indints02(:,:),NoOrbs)
                    ! Temp4indints02(l,j) comes out of here

                    do l=1,NoOrbs
                        do j=1,l
                            TwoIndInts02(g,a,j,l)=Temp4indints02(j,l)
                            TwoIndInts02(a,g,j,l)=Temp4indints02(j,l)
                            TwoIndInts02(g,a,l,j)=Temp4indints02(j,l)
                            TwoIndInts02(a,g,l,j)=Temp4indints02(j,l)
                        enddo
                    enddo
                enddo
            enddo
!            CALL MPIDSumArr(TwoIndInts02Temp(:,:,:,:),NoOrbs**4,TwoIndInts02(:,:,:,:))
        ENDIF


! Calculating the 3 transformed, 4 index integrals. 01=a untransformed,02=b,03=g,04=d


!        do i=LowBound02,HighBound02
        do i=1,NoRotOrbs
            do k=1,i
                Temp4indints(:,:)=0.0_dp
                CALL DGEMM('T','N',NoRotOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts01(:,:,k,i),NoOrbs, &
                    0.0,Temp4indints(:,:),NoRotOrbs)

                IF(tNotConverged) THEN
                    do b=1,NoOrbs
                        do l=1,NoOrbs
                            ThreeIndInts02(i,k,l,b)=Temp4indints(l,b)
                            ThreeIndInts02(k,i,l,b)=Temp4indints(l,b)
                        enddo
                    enddo
                    Temp4indints02(:,:)=0.0_dp
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts01(:,:,k,i),NoOrbs, &
                        0.0,Temp4indints02(:,:),NoRotOrbs)
                    do d=1,NoOrbs
                        do j=1,NoOrbs
                            ThreeIndInts04(k,i,j,d)=Temp4indints02(j,d)
                            ThreeIndInts04(i,k,j,d)=Temp4indints02(j,d)
                        enddo
                    enddo
                ENDIF
                Temp4indints02(:,:)=0.0_dp
                CALL DGEMM('T','T',NoRotOrbs,NoRotOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(:,:),NoRotOrbs, &
                    0.0,Temp4indints02(:,:),NoRotOrbs)
                do l=1,NoRotOrbs
                    do j=1,l
                        FourIndInts(i,j,k,l)=Temp4indints02(j,l)
                        FourIndInts(i,l,k,j)=Temp4indints02(j,l)
                        FourIndInts(k,j,i,l)=Temp4indints02(j,l)
                        FourIndInts(k,l,i,j)=Temp4indints02(j,l)

                        IF(tNotConverged) THEN
                            FourIndInts02(j,k,l,i)=Temp4indints02(j,l)
                            FourIndInts02(j,i,l,k)=Temp4indints02(j,l)
                            FourIndInts02(l,k,j,i)=Temp4indints02(j,l)
                            FourIndInts02(l,i,j,k)=Temp4indints02(j,l)
                        ENDIF
                    enddo
                enddo
            enddo
        enddo

        IF(tNotConverged) THEN
!            do l=LowBound02,HighBound02
            do l=1,NoOrbs
                do j=1,l
                    Temp4indints(:,:)=0.0_dp
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts02(:,:,j,l),NoOrbs, &
                        0.0,Temp4indints(:,:),NoOrbs)
                    do a=1,NoOrbs
                        do k=1,NoOrbs
                            ThreeIndInts01(k,j,l,a)=Temp4indints(k,a)
                            ThreeIndInts01(k,l,j,a)=Temp4indints(k,a)
                        enddo
                    enddo
                    Temp4indints(:,:)=0.0_dp
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts02(:,:,j,l),NoOrbs, &
                        0.0,Temp4indints(:,:),NoOrbs)
                    do g=1,NoOrbs
                        do i=1,NoOrbs
                            ThreeIndInts03(i,l,j,g)=Temp4indints(i,g)
                            ThreeIndInts03(i,j,l,g)=Temp4indints(i,g)
                        enddo
                    enddo
                enddo
            enddo
!            CALL MPIDSumArr(ThreeIndInts01Temp(:,:,:,:),NoOrbs**4,ThreeIndInts01(:,:,:,:))
!            CALL MPIDSumArr(ThreeIndInts02Temp(:,:,:,:),NoOrbs**4,ThreeIndInts02(:,:,:,:))
!            CALL MPIDSumArr(ThreeIndInts03Temp(:,:,:,:),NoOrbs**4,ThreeIndInts03(:,:,:,:))
!            CALL MPIDSumArr(ThreeIndInts04Temp(:,:,:,:),NoOrbs**4,ThreeIndInts04(:,:,:,:))
!            CALL MPIDSumArr(FourIndInts02Temp(:,:,:,:),NoOrbs**4,FourIndInts02(:,:,:,:))
        ENDIF



! ***************************
! Calc the potential energies for this iteration (with these transformed integrals).        

! This can be sped up by merging the calculations of the potentials with the transformations, but while 
! we are playing around with different potentials, it is simpler to keep these separate.
    
        IF((.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs).and.(.not.tUseHFOrbs)) THEN
            
            PotEnergy=0.0_dp
            TwoEInts=0.0_dp
            PEInts=0.0_dp
            CALL CalcPotentials()

            IF(tPrintInts) CALL PrintIntegrals()
            IF((Iteration.eq.0).or.((.not.tNotConverged).and.(Iteration.gt.1))) CALL WriteDoubHisttofile()
            IF(tROHistSingExc.and.(Iteration.eq.0)) CALL WriteSingHisttofile()


! If doing Lagrange orthormalisations, find the change of the potential energy due to the orthonormality 
! of the orbitals...
            IF(tLagrange) THEN
                PEOrtho=0.0_dp
                do i=1,NoOrbs
                    do j=1,NoOrbs
                        t=0.0_dp
                        do a=1,NoOrbs
                            t=CoeffT1(a,i)*CoeffT1(a,j)
                        enddo
                        IF(i.eq.j) t=t-1.0_dp
                        PEOrtho=PEOrtho-Lambdas(i,j)*t
                        PotEnergy=PotEnergy-Lambdas(i,j)*t
                    enddo
                enddo
            ENDIF
        ENDIF

        CALL halt_timer(Transform2ElInts_Time)


    END SUBROUTINE Transform2ElInts


    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.
    SUBROUTINE Transform2ElIntsMemSave()
        INTEGER :: i,j,k,l,a,b,g,d,ierr,a2,b2,g2,d2
        INTEGER(TagIntType) Temp4indintsTag
        real(dp) , ALLOCATABLE :: Temp4indints(:,:)
#ifdef __CMPLX
        call stop_all('Transform2ElIntsMemSave', 'Rotating orbitals not implemented for complex orbitals.')
#endif
        
        Transform2ElInts_Time%timer_name='Transform2ElIntsTime'
        CALL set_timer(Transform2ElInts_time,30)

!Zero arrays from previous transform

 
        ALLOCATE(Temp4indints(NoRotOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('Temp4indints',NoRotOrbs*NoOrbs,8,'Transform2ElIntsMemSave',Temp4indintsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All('Transform2ElIntsMemSave','Problem allocating memory to Temp4indints.')
 
        FourIndInts(:,:,:,:)=0.0_dp

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i,j,k,l,0,0)

        do b=1,NoOrbs
            IF(tTurnStoreSpinOff) THEN
                b2=CEILING(REAL(SymLabelList2_rot(b))/2.0_dp)
            ELSE
                b2=SymLabelList2_rot(b)
            ENDIF
            do d=1,b
                IF(tTurnStoreSpinOff) THEN
                    d2=CEILING(REAL(SymLabelList2_rot(d))/2.0_dp)
                ELSE
                    d2=SymLabelList2_rot(d)
                ENDIF
                do a=1,NoOrbs
                    IF(tTurnStoreSpinOff) THEN
                        a2=CEILING(REAL(SymLabelList2_rot(a))/2.0_dp)
                    ELSE
                        a2=SymLabelList2_rot(a)
                    ENDIF
                    do g=1,a
                        IF(tTurnStoreSpinOff) THEN
                            g2=CEILING(REAL(SymLabelList2_rot(g))/2.0_dp)
                        ELSE
                            g2=SymLabelList2_rot(g)
                        ENDIF
                        FourIndInts(a,g,b,d)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),dp)
                        FourIndInts(g,a,b,d)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),dp)
                        FourIndInts(a,g,d,b)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),dp)
                        FourIndInts(g,a,d,b)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),dp)
                    enddo
                enddo
                Temp4indints(:,:)=0.0_dp
                CALL DGEMM('T','N',NoRotOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,FourIndInts(1:NoOrbs,1:NoOrbs,b,d), &
                    NoOrbs,0.0,Temp4indints(1:NoRotOrbs,1:NoOrbs),NoRotOrbs)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, we need the transpose of this.

                CALL DGEMM('T','T',NoRotOrbs,NoRotOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(1:NoRotOrbs,1:NoOrbs), &
                    NoRotOrbs,0.0,FourIndInts(1:NoRotOrbs,1:NoRotOrbs,b,d),NoRotOrbs)
                ! Get Temp4indits02(i,k)

                do i=1,NoRotOrbs
                    do k=1,i
                        FourIndInts(i,k,d,b)=FourIndInts(i,k,b,d)
                        FourIndInts(k,i,d,b)=FourIndInts(i,k,b,d)
                        FourIndInts(i,k,b,d)=FourIndInts(i,k,b,d)
                        FourIndInts(k,i,b,d)=FourIndInts(i,k,b,d)
                    enddo
                enddo
            enddo
        enddo
        

! Calculating the 3 transformed, 4 index integrals. 01=a untransformed,02=b,03=g,04=d
        do i=1,NoRotOrbs
            do k=1,i

                Temp4indints(:,:)=0.0_dp
                CALL DGEMM('T','N',NoRotOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,FourIndInts(i,k,1:NoOrbs,1:NoOrbs), &
                    NoOrbs,0.0,Temp4indints(1:NoRotOrbs,1:NoOrbs),NoRotOrbs)

                CALL DGEMM('T','T',NoRotOrbs,NoRotOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(1:NoRotOrbs,1:NoOrbs), &
                    NoRotOrbs,0.0,FourIndInts(i,k,1:NoRotOrbs,1:NoRotOrbs),NoRotOrbs)
                do l=1,NoRotOrbs
                    do j=1,l
                        FourIndInts(k,i,j,l)=FourIndInts(i,k,j,l)
                        FourIndInts(k,i,l,j)=FourIndInts(i,k,j,l)
                        FourIndInts(i,k,j,l)=FourIndInts(i,k,j,l)
                        FourIndInts(i,k,l,j)=FourIndInts(i,k,j,l)
                    enddo
                enddo
            enddo
        enddo
 
        DEALLOCATE(Temp4indints)
        CALL LogMemDeAlloc('Transform2ElIntsMemSave',Temp4indintsTag)
 
        CALL halt_timer(Transform2ElInts_Time)


    END SUBROUTINE Transform2ElIntsMemSave



   
! This is a transformation of the four index integrals for the ERlocalisation, in this only the <ii|ii> integrals are needed 
! therefore the process may be much simpler.
    SUBROUTINE Transform2ElIntsERlocal()
        INTEGER :: i,j,a,b,g,d,m
        real(dp) :: t,Temp4indints(NoOrbs,NoOrbs)
        real(dp) :: Temp4indints02(NoOrbs)  
 

        CALL set_timer(Transform2ElInts_time,30)

! Zero arrays from previous transform

        TwoIndIntsER(:,:,:)=0.0_dp

        ThreeIndInts01ER(:,:)=0.0_dp
        ThreeIndInts02ER(:,:)=0.0_dp

        FourIndIntsER(:)=0.0_dp

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i,j,k,l,0,0)


!        LowBound=iProcIndex*(NoOrbs/nProcessors)+1
!        HighBound=(iProcIndex+1)*(NoOrbs/nProcessors)
!        IF(iProcIndex.eq.(nProcessors-1)) HighBound=NoOrbs


! UMATTemp01(a,g,b,d)

!        do a=LowBound02,HighBound02
        do d=1,NoOrbs
            do b=1,d
                Temp4indints(:,:)=0.0_dp
                Temp4indints02(:)=0.0_dp
                CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMATTemp01(:,:,b,d), &
                    NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
                ! a -> m. Temp4indints(m,g) comes out of here.
                ! Want to transform g to m as well.

                do m=1,NoOrbs
                    do g=1,NoOrbs
                        Temp4indints02(m)=Temp4indints02(m)+(Temp4indints(m,g)*CoeffT1(g,m))
                    enddo
                enddo
                ! Now have Temp4indints(m,m) for each b and d.

                do m=1,NoOrbs
                    TwoIndIntsER(b,d,m)=Temp4indints02(m)
                    TwoIndIntsER(d,b,m)=Temp4indints02(m)
                enddo
            enddo
        enddo
        
! Now want to transform g to get one of the 3-transformed 4-index integrals <a m | m m>.
! These can be stored in 2-D arrays, as they can be specified by only m and z.

        do m=1,NoOrbs
            do b=1,NoOrbs
                do d=1,NoOrbs
                    ThreeIndInts01ER(b,m)=ThreeIndInts01ER(b,m)+(TwoIndIntsER(b,d,m)*CoeffT1(d,m))
                enddo
            enddo
        enddo
        ! ThreeIndInts01ER(z,m) is where z is alpha (a).


        TwoIndIntsER(:,:,:)=0.0_dp
!        do a=LowBound02,HighBound02
        do d=1,NoOrbs
            do b=1,NoOrbs
                Temp4indints(:,:)=0.0_dp
                Temp4indints02(:)=0.0_dp
                CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMATTemp01(:,:,b,d), &
                    NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
                ! a -> m. Temp4indints(m,g) comes out of here.
                ! Want to transform g to m as well.

                do m=1,NoOrbs
                    do g=1,NoOrbs
                        Temp4indints02(m)=Temp4indints02(m)+(Temp4indints(m,g)*CoeffT1(g,m))
                    enddo
                enddo
                ! Now have Temp4indints(m,m) for each a and g.

                do m=1,NoOrbs
                    TwoIndIntsER(b,d,m)=Temp4indints02(m)
                    TwoIndIntsER(d,b,m)=Temp4indints02(m)
                enddo
            enddo
        enddo

! Now want to transform g to get one of the 3-transformed 4-index integrals <a m | m m>.
! These can be stored in 2-D arrays, as they can be specified by only m and z.

        do m=1,NoOrbs
            do b=1,NoOrbs
                do d=1,NoOrbs
                    ThreeIndInts02ER(b,m)=ThreeIndInts02ER(b,m)+(TwoIndIntsER(b,d,m)*CoeffT1(d,m))
                enddo
            enddo
        enddo
        ! ThreeIndInts02ER(z,m) is where z is beta (b).

! Find the <ii|ii> integrals, to calculate the potential energy.

        do m=1,NoOrbs
            do a=1,NoOrbs
                FourIndIntsER(m)=FourIndIntsER(m)+(ThreeIndInts01ER(a,m)*CoeffT1(a,m))
            enddo
        enddo

! ***************************
! Calc the potential energies for this iteration (with these transformed integrals).        

! This can be sped up by merging the calculations of the potentials with the transformations, but while 
! we are playing around with different potentials, it is simpler to keep these separate.
        
        PotEnergy=0.0_dp
        TwoEInts=0.0_dp
        PEInts=0.0_dp
        CALL CalcPotentials()

        IF(tPrintInts) CALL PrintIntegrals()
        IF((Iteration.eq.0).or.((.not.tNotConverged).and.(Iteration.gt.1))) CALL WriteDoubHisttofile()
        IF(tROHistSingExc.and.(Iteration.eq.0)) CALL WriteSingHisttofile()


! If doing Lagrange orthormalisations, find the change of the potential energy due to the orthonormality 
! of the orbitals...
        IF(tLagrange) THEN
            PEOrtho=0.0_dp
            do i=1,NoOrbs
                do j=1,NoOrbs
                    t=0.0_dp
                    do a=1,NoOrbs
                        t=CoeffT1(a,i)*CoeffT1(a,j)
                    enddo
                    IF(i.eq.j) t=t-1.0_dp
                    PEOrtho=PEOrtho-Lambdas(i,j)*t
                    PotEnergy=PotEnergy-Lambdas(i,j)*t
                enddo
            enddo
        ENDIF

        CALL halt_timer(Transform2ElInts_Time)


    END SUBROUTINE Transform2ElIntsERlocal


    SUBROUTINE CalcPotentials()
    ! only temporarily like this, can tidy it up majorly
        INTEGER :: i,j,k,l,Starti,Finishi
        real(dp) :: MaxTerm

        l = 0
        IF(tERLocalization.and.(.not.tStoreSpinOrbs)) THEN
            ERPotEnergy=0.0_dp
            IF(tRotateVirtOnly) THEN 
                Starti=NoOcc+1
                Finishi=NoOrbs
            ELSEIF(tRotateOccOnly) THEN
                Starti=1
                Finishi=NoOcc
            ELSE
                Starti=1
                Finishi=NoOrbs
            ENDIF
            CoulPotEnergy=0.0_dp
            OffDiagPotEnergy=0.0_dp
!            do i=1,NoOrbs
!                IF((i.ge.Starti).and.(i.le.Finishi)) THEN
            do i=Starti,Finishi
                    ERPotEnergy=ERPotEnergy+FourIndIntsER(i)
                    IF(FourIndIntsER(i).lt.0) THEN
                        CALL neci_flush(6)
                        CALL Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    ENDIF
!                    WRITE(6,*) FourIndIntsER(i)
                    PotEnergy=PotEnergy+FourIndIntsER(i)
                    TwoEInts=TwoEInts+FourIndIntsER(i)
                    PEInts=PEInts+FourIndIntsER(i)
!                ENDIF
!                do k=i+1,NoOrbs
!                    IF((i.ge.Starti).and.(i.le.Finishi)) CoulPotEnergy=CoulPotEnergy+FourIndInts(i,k,i,k)
!                    do j=1,NoOrbs
!                        do l=j+1,NoOrbs
!                            OffDiagPotEnergy=OffDiagPotEnergy+FourIndInts(i,j,k,l)
!                        enddo
!                    enddo
!                enddo
            enddo
        ELSEIF(tERLocalization) THEN
            ERPotEnergy=0.0_dp
            PotEnergy=0.0_dp
            IF(tRotateVirtOnly) THEN 
                Starti=NoOcc+1
                Finishi=NoOrbs
            ELSEIF(tRotateOccOnly) THEN
                Starti=1
                Finishi=NoOcc
            ELSE
                Starti=1
                Finishi=NoOrbs
            ENDIF
            do i=Starti,Finishi
                IF(tStoreSpinOrbs) THEN
                    IF(MOD(i,2).eq.0) THEN
                        j=i-1
                    ELSE
                        j=i+1
                    ENDIF
                    ERPotEnergy=ERPotEnergy+FourIndInts(i,j,i,j)
                    IF((FourIndInts(i,j,i,j).lt.0).or.(FourIndInts(j,i,j,i).lt.0)) THEN
                        CALL neci_flush(6)
                        CALL Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    ENDIF
                    PotEnergy=PotEnergy+FourIndInts(i,j,i,j)
                    TwoEInts=TwoEInts+FourIndInts(i,j,i,j)
                    PEInts=PEInts+FourIndInts(i,j,i,j)
                ELSE
                    ERPotEnergy=ERPotEnergy+FourIndInts(i,i,i,i)
                    IF((FourIndInts(i,i,i,i).lt.0)) THEN
                        CALL neci_flush(6)
                        CALL Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    ENDIF
                    PotEnergy=PotEnergy+FourIndInts(i,i,i,i)
                    TwoEInts=TwoEInts+FourIndInts(i,i,i,i)
                    PEInts=PEInts+FourIndInts(i,i,i,i)
                ENDIF
            enddo
        ENDIF

        IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax.or.tOffDiagMin.or.tOffdiagMax) THEN
            do l=1,NoOrbs
                do j=1,l-1
!                     do k=1,NoOrbs
                    do k=1,j-1
                        do i=1,k-1
                            IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) THEN
                                IF(((i.ne.j).and.(j.ne.l)).and.((i.ne.k).or.(j.ne.l))) THEN
                                    PotEnergy=PotEnergy+(FourIndInts(i,j,k,l)**2)
                                    TwoEInts=TwoEInts+(FourIndInts(i,j,k,l)**2)
                                    PEInts=PEInts+(FourIndInts(i,j,k,l)**2)
                                ENDIF
                            ENDIF
                            IF(tOffDiagMin.or.tOffDiagMax) THEN
                                IF(.not.((k.eq.i).or.(j.eq.l))) THEN
                                    PotEnergy=PotEnergy+FourIndInts(i,j,k,l)
                                    TwoEInts=TwoEInts+FourIndInts(i,j,k,l)
                                    PEInts=PEInts+FourIndInts(i,j,k,l)
                                ENDIF
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo
        ENDIF
      
        IF(tDoubExcMin) THEN
            do i=1,NoOrbs
                do j=1,NoOrbs
                    do k=1,i-1
                        IF((k.eq.l).and.(k.eq.i)) CYCLE
                        do l=1,j-1
                            IF((j.eq.k).and.(j.eq.l)) CYCLE
                            IF((j.eq.k).and.(j.eq.i)) CYCLE
                            IF((j.eq.l).and.(j.eq.i)) CYCLE
                            PotEnergy=PotEnergy+(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            TwoEInts=TwoEInts+(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            PEInts=PEInts+(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                        enddo
                    enddo
                enddo
            enddo
        ENDIF
 
        IF(tOnePartOrbEnMax.or.tOneElIntMax) THEN
            do i=NoOcc+1,NoOrbs
                MaxTerm=0.0_dp
                MaxTerm=TMAT2DRot(i,i)
                IF(tOnePartOrbEnMax) THEN
                    do j=1,NoOcc
                        MaxTerm=MaxTerm+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                    enddo
                    MaxTerm=MaxTerm-EpsilonMin
                    MaxTerm=MaxTerm**OrbEnMaxAlpha
                ENDIF
                PotEnergy=PotEnergy+MaxTerm
            enddo
        ENDIF

        IF(tHijSqrdMin) THEN
            HijSqrdPotEnergy=0.0_dp
            do i=NoOcc+1,NoOrbs
                do j=NoOcc+1,NoOrbs
                    IF(j.gt.i) THEN
                        PotEnergy=PotEnergy+(TMAT2DRot(i,j)**2)
                        HijSqrdPotEnergy=HijSqrdPotEnergy+(TMAT2DRot(i,j)**2)
                    ENDIF
                enddo
            enddo
        ENDIF

        IF(tVirtCoulombMax) THEN
            ERPotEnergy=0.0_dp
            ijOccVirtPotEnergy=0.0_dp
            do i=1,NoOrbs
                IF(i.le.NoOcc) THEN
                    do j=NoOcc+1,NoOrbs
                        ijOccVirtPotEnergy=ijOccVirtPotEnergy+FourIndInts(i,j,i,j)
                    enddo
                ENDIF
                IF(i.gt.NoOcc) THEN
                    ERPotEnergy=ERPotEnergy+FourIndInts(i,i,i,i)
                    do j=NoOcc+1,NoOrbs
                        IF(j.le.i) CYCLE
                        PotEnergy=PotEnergy+FourIndInts(i,j,i,j)
                        TwoEInts=TwoEInts+FourIndInts(i,j,i,j)
                    enddo
                ENDIF
            enddo
        ENDIF

        IF(tHFSingDoubExcMax) THEN
            do i=1,NoOcc
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do l=NoOcc+1,NoOrbs
                            PotEnergy=PotEnergy+(FourIndInts(i,j,k,l)**2)
                        enddo
                        
                        !Sing excitations <ij|ik> where i and j are occ, k virt.
                        PotEnergy=PotEnergy+(FourIndInts(i,j,i,k)**2)
                    enddo
                enddo
            enddo
        ENDIF


    ENDSUBROUTINE CalcPotentials




    SUBROUTINE FindTheForce()
        INTEGER :: m,z,i,j,k,l,a,Symm,w,x,y,SymMin
        real(dp) :: OffDiagForcemz,DiagForcemz,OneElForcemz,LambdaTerm1,LambdaTerm2
        real(dp) :: NonDerivTerm,DerivPot
        LOGICAL :: leqm,jeqm,keqm
      
! Running over m and z, covers all matrix elements of the force matrix (derivative 
! of equation we are minimising, with respect to each translation coefficient) filling 
! them in as it goes.
        CALL set_timer(FindtheForce_time,30)
        
        DerivCoeff(:,:)=0.0_dp
        Force=0.0_dp
!        ForceTemp=0.0_dp
        ForceInts=0.0_dp
        OrthoForce=0.0_dp
        OffDiagForceMZ = 0
        SymMin = 0
        i = 0

        ! If the orbitals are being separated, do this whole loop twice, once for occupied and once for virtual
        ! i.e w = 1,2. Otherwise do them all at once.
        do w=MinOccVirt,MaxOccVirt
            IF(w.eq.1) THEN
                SymMin=1
                MinMZ=1
                IF(tSeparateOccVirt) THEN
                    MaxMZ=NoOcc
                ELSE
                    MaxMZ=NoOrbs
                ENDIF
            ELSE
                SymMin=9
                MinMZ=NoOcc+1
                MaxMZ=NoOrbs
            ENDIF
! If we are localising the occupied and virtual orbitals separately, the above block ensures that we loop over
! first the occupied then the virtual.  If we are not separating the orbitals we just run over all orbitals.

!            LowBound=iProcIndex*((MaxMZ-MinMZ)/nProcessors)+MinMZ
!            HighBound=(iProcIndex+1)*((MaxMZ-MinMZ)/nProcessors)+MinMZ-1
!            IF(iProcIndex.eq.(nProcessors-1)) HighBound=MaxMZ

!            do m=LowBound,HighBound
            do m=MinMZ,MaxMZ
                IF(tStoreSpinOrbs) THEN
                    SymM=INT(G1(SymLabelList2_rot(m))%sym%S)
                ELSE
                    SymM=INT(G1(SymLabelList2_rot(m)*2)%sym%S)
                ENDIF
                do z=SymLabelCounts2_rot(1,SymM+SymMin), &
                        (SymLabelCounts2_rot(1,SymM+SymMin) + &
                            SymLabelCounts2_rot(2,SymM+SymMin)-1)
 
                    ! Find the force on a coefficient c(m,z). 
                    OffDiagForcemz=0.0_dp
                    ! OffDiagForce is from any of the OffDiagMin/Max (Sqrd or not), or the double/single excitation
                    ! max/min, as only one of these terms may be used at once.

                    DiagForcemz=0.0_dp
                    ! DiagForce includes ER localisation, and the coulomb terms <ij|ij>.

                    OneElForcemz=0.0_dp
                    ! OneElForce includes that from the one electron integrals <i|h|j> and the one particle orbital 
                    ! energies.

                    ! DIAG TERMS
                    ! Maximise <ii|ii>, self interaction terms. 
                    IF(tERLocalization.and.(.not.tStoreSpinOrbs)) THEN
!                        DiagForcemz=DiagForcemz+ThreeIndInts01(m,m,m,z)+ThreeIndInts02(m,m,m,z)+ThreeIndInts03(m,m,m,z)+ThreeIndInts04(m,m,m,z)
                        DiagForcemz=DiagForcemz+(2*ThreeIndInts01ER(z,m))+(2*ThreeIndInts02ER(z,m))
                        ! Derivative of <ii|ii> only non-zero when i=m.
                        ! each of the four terms then correspond to zeta = a, b, g, then d in the unrotated basis.
                    ELSEIF(tERLocalization) THEN
!                    IF(tERLocalization) THEN
                        ! Looking at <ij|ij> terms where j=i+1 (i.e. i is alpha of spin orbital and j is beta - or vice versa).
                        IF(tStoreSpinOrbs) THEN
                            IF(MOD(m,2).eq.0) THEN      ! m = j
                                i=m-1
                                DiagForcemz=DiagForcemz+ThreeIndInts01(m,i,i,z)+ThreeIndInts01(z,i,i,m)+ &
                                    ThreeIndInts01(i,m,z,i)+ThreeIndInts01(i,z,m,i)
                            ELSE
                                j=m+1
                                DiagForcemz=DiagForcemz+ThreeIndInts01(m,j,j,z)+ThreeIndInts01(z,j,j,m)+ &
                                    ThreeIndInts01(j,m,z,j)+ThreeIndInts01(j,z,m,j)
                            ENDIF
                        ELSE
                            DiagForcemz=DiagForcemz+ThreeIndInts01(m,m,m,z)+ThreeIndInts02(m,m,m,z)+ &
                                ThreeIndInts03(m,m,m,z)+ThreeIndInts04(m,m,m,z)
                        ENDIF
                        ! First term when m=i and z=a, second when m=i and z=g.
                    ENDIF

 
                    ! Maximise <ij|ij>, coulomb terms, where i<j, i occ or virt, j virt only.
                    IF(tVirtCoulombMax) THEN
                        do i=1,NoOrbs
                            IF(i.eq.m) THEN
                                do j=NoOcc+1,NoOrbs
                                    IF(j.le.i) CYCLE        ! i<j.
                                    DiagForcemz=DiagForcemz+ThreeIndInts01(m,j,j,z)+ThreeIndInts03(m,j,j,z)
                                    ! First term for when m=i and z=a, second when m=i and z=g.
                                enddo
                            ENDIF
                            IF((m.gt.NoOcc).and.(m.gt.i)) DiagForcemz=DiagForcemz+ThreeIndInts02(i,i,m,z)+ThreeIndInts04(i,i,m,z)
                            ! This only contributes when j=m (no point in running over all j.
                            ! First term when m=j and z=b, second when m=j and z=d.
                        enddo
                    ENDIF

                    ! ONE ELECTRON TERMS
                    ! Minimise |<i|h|j>|^2 where either one or bot of i and j are virtual, but i<j.
                    IF(tHijSqrdMin) THEN
                        do j=NoOcc+1,NoOrbs
                            IF(m.ne.j) OneElForcemz=OneElForcemz+(2*TMAT2DRot(m,j)*TMAT2DPartRot02(z,j))
                            ! m=i and z=a.
                        enddo
                        do i=NoOcc+1,NoOrbs
                            IF(m.ne.i) OneElForcemz=OneElForcemz+(2*TMAT2DRot(i,m)*TMAT2DPartRot01(i,z))
                            ! m=j and z=b
                        enddo
                    ENDIF

                    ! OnePartOrbEnMax ; Maximisie sum_i [E_i - E_min]^Alpha
                    ! where E_i = <i|h|i> + sum_j <ij||ij> and E_min is either E_LUMO (rotating virtual only) or the chemical 
                    ! potential (midway between LUMO and HOMO, when rotating all), Alpha specified in input.
                    ! The derivative of the one part orb energies is then Alpha * NonDerivTerm * DerivPot^(Alpha-1)  

                    ! OneElIntMax ; Maximise <i|h|i>
                    IF(tOnePartOrbEnMax.or.tOneElIntMax) THEN
                        do i=NoOcc+1,NoOrbs
                            DerivPot=0.0_dp
                            DerivPot=DerivPot+TMAT2DPartRot02(z,m)+TMAT2DPartRot01(m,z)
                            ! First term when m=i and z=a, second when m=i and z=b.
                            ! This is all that is needed for OneElIntMax

                            IF(tOnePartOrbEnMax) THEN
                                NonDerivTerm=0.0_dp
                                IF(OrbEnMaxAlpha.ne.1.0_dp) THEN 
                                    ! The non-derived term in the chain rule, <i|h|i> + sum_j <ij||ij> - E_min.
                                    NonDerivTerm=NonDerivTerm+TMAT2DRot(i,i)-EpsilonMin
                                    do j=1,NoOcc
                                        NonDerivTerm=NonDerivTerm+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                                    enddo
                                    NonDerivTerm=OrbEnMaxAlpha*(NonDerivTerm**(OrbEnMaxAlpha-1))
                                ELSE
                                    ! If Alpha = 1, the NonDerivTerm will be raised to the power of 0, thus always 1.
                                    NonDerivTerm=1.0
                                ENDIF
                                IF(i.eq.m) THEN
                                    do j=1,NoOcc
                                        DerivPot=DerivPot+(2*ThreeIndInts01(m,j,j,z))-ThreeIndInts01(j,j,m,z)+ &
                                            (2*ThreeIndInts03(m,j,j,z))-ThreeIndInts03(m,j,z,j)
                                        ! First part is for when m=i and z=a, the second is for when m=i and z=g
                                    enddo
                                ENDIF
                                ! When m=j, for a particular i.
                                ! m and z run only over virtual, and j is over occupied. m will never = j.
                            ELSE
                                NonDerivTerm=1.0_dp
                            ENDIF

                            OneElForcemz=OneElForcemz+(NonDerivTerm*DerivPot)
                        enddo
                    ENDIF

                    ! OFFDIAGTERMS
                    ! Maximises the square of the single and double excitation integrals connected to the HF.
                    ! I.e maximises <ij|kl> where i,j are occupied and k,l are virtual (doubles), except k may be occuppied if
                    ! equal to i (<ij|il> singles).
                    ! Currently this is only used for rotating virtual only, so m can only equal k or l.
                    IF(tHFSingDoubExcMax) THEN
                        do i=1,NoOcc
                            do j=1,NoOcc
                                do k=NoOcc+1,NoOrbs
                                    IF(k.eq.m) THEN
                                        do l=NoOcc+1,NoOrbs
                                            OffDiagForcemz=OffDiagForcemz+(2*FourIndInts(i,j,m,l)*ThreeIndInts03(i,j,l,z))
                                            ! m=k and z=g.
                                        enddo
                                    ENDIF
                                    
                                    OffDiagForcemz=OffDiagForcemz+(2*FourIndInts(i,j,k,m)*ThreeIndInts04(i,k,j,z))
                                    ! m=l and z=d. 
                                    
                                    !Sing excitations <ij|il> where i and j are occ, l virt.
                                    OffDiagForcemz=OffDiagForcemz+(2*FourIndInts(i,j,i,m)*ThreeIndInts04(i,i,j,z))
                                    ! m=l
                                enddo
                            enddo
                        enddo
                    ENDIF

                    ! OffDiag Sqrd/notSqrd Min/Max treats the elements <ij|kl>
                    ! i<k and j<l.
                    IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax.or.tOffDiagMin.or.tOffDiagMax.or.tDoubExcMin) THEN
                        do l=1,NoOrbs
                            IF(l.eq.m) THEN
                                leqm=.true.
                            ELSE
                                leqm=.false.
                            ENDIF
                            do j=1,l-1                        
                                IF(j.eq.m) THEN
                                    jeqm=.true.
                                ELSE
                                    jeqm=.false.
                                ENDIF
                                do k=1,j-1                 
                                    IF(k.eq.l) CYCLE
                                    IF(k.eq.m) THEN
                                        keqm=.true.
                                    ELSE
                                        keqm=.false.
                                    ENDIF
!                                    Symi=IEOR(INT(G1(SymLabelList2_rot(k)*2)%sym%S),IEOR(INT(G1(SymLabelList2_rot(j)*2)%sym%S),INT(G1(SymLabelList2_rot(l)*2)%sym%S)))
                                    ! only i with symmetry equal to j x k x l will have integrals with overall
                                    ! symmetry A1 and therefore be non-zero.


                                    ! Running across i, ThreeIndInts01 only contributes if i.eq.m (which will happen once for each m)
                                    IF((m.le.k-1).and.(m.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
!                                    IF((m.le.k-1)) THEN
                                        IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+2* &
                                            (FourIndInts02(j,k,l,m)*ThreeIndInts01(k,j,l,z))
                                        IF(tOffDiagMin.or.tOffDiagMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts01(k,j,l,z)
                                        IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+ThreeIndInts01(k,j,l,z)- &
                                            ThreeIndInts01(l,j,k,z)
                                    ENDIF

                                    IF(jeqm) THEN
                                        do i=1,k-1
                                            IF((i.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
!                                        do i=SymLabelCounts2_rot(1,Symi+SymMin),(SymLabelCounts2_rot(1,Symi+SymMin)+SymLabelCounts2_rot(2,Symi+SymMin)-1)
                                                IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+2* &
                                                    (FourIndInts(i,j,k,l)*ThreeIndInts02(i,k,l,z))
                                                IF(tOffDiagMin.or.tOffDiagMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts02(i,k,l,z)
                                                IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+(ThreeIndInts02(i,k,l,z))- &
                                                    ThreeIndInts02(i,l,k,z)
                                            ENDIF
                                        enddo
                                    ENDIF

                                    IF(keqm) THEN
!                                        do i=SymLabelCounts2_rot(1,Symi+SymMin),(SymLabelCounts2_rot(1,Symi+SymMin)+SymLabelCounts2_rot(2,Symi+SymMin)-1)
                                        do i=1,k-1
                                            IF((i.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
                                                IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+2* &
                                                    (FourIndInts(i,j,k,l)*ThreeIndInts03(i,j,l,z))
                                                IF(tOffDiagMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+ &
                                                    ThreeIndInts03(i,j,l,z)
                                                IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+(ThreeIndInts03(i,j,l,z))- &
                                                    ThreeIndInts03(i,j,z,l)
                                            ENDIF
                                        enddo
                                    ENDIF

                                    IF(leqm) THEN
!                                        do i=SymLabelCounts2_rot(1,Symi+SymMin),(SymLabelCounts2_rot(1,Symi+SymMin)+SymLabelCounts2_rot(2,Symi+SymMin)-1)
                                        do i=1,k-1
                                            IF((i.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
                                                IF(tOffDiagSqrdMin.or.tOffDiagSqrdMin) OffDiagForcemz=OffDiagForcemz+2* &
                                                    (FourIndInts(i,j,k,l)*ThreeIndInts04(i,k,j,z))
                                                IF(tOffDiagMin.or.tOffDiagMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts04(i,k,j,z)
                                                IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+(ThreeIndInts04(i,k,j,z))- &
                                                    ThreeIndInts04(i,z,j,k)
                                            ENDIF
                                        enddo
                                    ENDIF
                                enddo
                            enddo
                        enddo
                    ENDIF
                    ! DerivCoeffTemp(z,m) then combines all the different forces on coefficient(m,z).
!                    IF(tStoreSpinOrbs) THEN
!                        WRITE(6,*) CEILING(m/2.0),CEILING(z/2.0),DiagForcemz
!                    ELSE
!                        WRITE(6,*) m,z,DiagForcemz
!                    ENDIF
                    DerivCoeff(z,m)=(MaxMinFac*OffDiagWeight*OffDiagForcemz)+(DiagMaxMinFac*DiagWeight*DiagForcemz)+ &
                        (OneElMaxMinFac*OneElWeight*OneElForcemz)
                    Force=Force+ABS(DerivCoeff(z,m))
                enddo
            enddo
        enddo
!        CALL MPIDSumArr(DerivCoeffTemp(:,:),NoOrbs**2,DerivCoeff(:,:))
!        CALL MPIDSum(ForceTemp,1,Force)

        Force=Force/REAL(NoOrbs**2,dp)


! Calculate the derivatives of orthogonalisation condition.
! Have taken this out of the m and z loop to make the shake faster, but can put it back in if start using it a lot.
        IF(tLagrange) THEN
            do x=MinMZ,MaxMZ
                m=SymLabelList2_rot(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM=INT(G1(m*2)%sym%S)
                do y=SymLabelCounts2_rot(1,SymM+SymMin), &
                        (SymLabelCounts2_rot(1,SymM+SymMin) + &
                            SymLabelCounts2_rot(2,SymM+SymMin)-1)
                    z=SymLabelList2_rot(y)
      
                    LambdaTerm1=0.0_dp
                    LambdaTerm2=0.0_dp
                    
                    do j=1,NoOrbs
                        LambdaTerm1=LambdaTerm1+(Lambdas(m,j)*CoeffT1(z,j))
                        LambdaTerm2=LambdaTerm2+(Lambdas(j,m)*CoeffT1(z,j))
                    enddo

! DerivCoeff is 'the force'.  I.e. the derivative of |<ij|kl>|^2 with 
! respect to each transformation coefficient.  It is the values of this matrix that will tend to 0 as
! we minimise the sum of the |<ij|kl>|^2 values.
! With the Lagrange keyword this includes orthonormality conditions, otherwise it is simply the unconstrained force.
                    DerivCoeff(z,m)=(2*OffDiagForcemz)-LambdaTerm1-LambdaTerm2
                    OrthoForce=OrthoForce-LambdaTerm1-LambdaTerm2
                enddo
            enddo
 
!If doing a lagrange calc we also need to find the force on the lambdas to ensure orthonormality...
            OrthoForce=OrthoForce/REAL(NoOrbs**2,dp)
            DerivLambda(:,:)=0.0_dp
            do i=1,NoOrbs
                do j=1,i
                    do a=1,NoOrbs
                        DerivLambda(i,j)=DerivLambda(i,j)+CoeffT1(a,i)*CoeffT1(a,j)
                    enddo
                    DerivLambda(j,i)=DerivLambda(i,j)
                enddo
            enddo
            do i=1,NoOrbs
                DerivLambda(i,i)=DerivLambda(i,i)-1.0_dp
            enddo
        ENDIF


        CALL halt_timer(FindtheForce_Time)


    END SUBROUTINE FindTheForce


    
    SUBROUTINE UseTheForce()
! This routine takes the old translation coefficients and Lambdas and moves them by a timestep in the direction 
! of the calculated force.
        INTEGER :: m,w,z,i,j,Symm,SymMin
        real(dp) :: NewCoeff,NewLambda

        DistCs=0.0_dp 
    
        do w=MinOccVirt,MaxOccVirt
            IF(w.eq.1) THEN
                SymMin=1
                MinMZ=1
                IF(tSeparateOccVirt) THEN
                    MaxMZ=NoOcc
                ELSE
                    MaxMZ=NoOrbs
                ENDIF
            ELSE
                SymMin=9
                MinMZ=NoOcc+1
                MaxMZ=NoOrbs
            ENDIF

!            LowBound=iProcIndex*((MaxMZ-MinMZ)/nProcessors)+MinMZ
!            HighBound=(iProcIndex+1)*((MaxMZ-MinMZ)/nProcessors)+MinMZ-1
!            IF(iProcIndex.eq.(nProcessors-1)) HighBound=MaxMZ
     
!            do m=LowBound,HighBound
            do m=MinMZ,MaxMZ
                IF(tStoreSpinOrbs) THEN
                    SymM=INT(G1(SymLabelList2_rot(m))%sym%S)
                ELSE
                    SymM=INT(G1(SymLabelList2_rot(m)*2)%sym%S)
                ENDIF

! Symmetry requirement that z must be from the same irrep as m
                do z=SymLabelCounts2_rot(1,SymM+SymMin), &
                        (SymLabelCounts2_rot(1,SymM+SymMin) + &
                            SymLabelCounts2_rot(2,SymM+SymMin)-1)
!               
                    ! Only coeffs with sym of m and z the same have non-zero coeffs.    
                    NewCoeff=0.0_dp
                    NewCoeff=CoeffT1(z,m)-(TimeStep*DerivCoeff(z,m))
                    DistCs=DistCs+abs(TimeStep*DerivCoeff(z,m))
                    CoeffT1(z,m)=NewCoeff
                enddo
            enddo
        enddo
!        CALL MPIDSumArr(CoeffT1Temp(:,:),NoOrbs**2,CoeffT1(:,:))
!        CALL MPIDSum(DistCsTemp,1,DistCs)


        DistCs=DistCs/(REAL(NoOrbs**2,dp))


        IF(tLagrange) THEN

            DistLs=0.0_dp
            LambdaMag=0.0_dp
            do i=1,NoOrbs
                do j=1,NoOrbs
                    NewLambda=0.0_dp
                    NewLambda=Lambdas(i,j)-(TimeStep*DerivLambda(i,j))  ! Timestep must be specified in the input file.
                    DistLs=DistLs+abs(TimeStep*DerivLambda(i,j))
                    Lambdas(i,j)=NewLambda
                    LambdaMag=LambdaMag+abs(NewLambda)
                enddo
            enddo
            DistLs=DistLs/(REAL(NoOrbs**2,dp))
            LambdaMag=LambdaMag/(REAL(NoOrbs**2,dp))

!        ELSE
!            CALL OrthoNormx(NoOrbs,NoOrbs,CoeffT1) !Explicitly orthonormalize the coefficient vectors.
        ENDIF


    ENDSUBROUTINE UseTheForce


   
    SUBROUTINE TestOrthonormality()
        INTEGER :: i,j
        real(dp) :: OrthoNormDP

        OrthoNorm=0.0_dp
        do i=1,NoOrbs
            do j=1,i
                OrthoNormDP=0.0_dp
                OrthoNormDP=Dot_Product(CoeffT1(:,i),CoeffT1(:,j))
                OrthoNorm=OrthoNorm+ABS(OrthoNormDP)
            enddo
        enddo
        OrthoNorm=OrthoNorm-real(NoOrbs,dp)
        OrthoNorm=(OrthoNorm*2.0_dp)/REAL((NoOrbs*(NoOrbs+1.0_dp)),dp)

    END SUBROUTINE TestOrthonormality



    SUBROUTINE TestForConvergence()
!This just tests the convergence on the grounds that the force is smaller that the input parameter: ConvergedForce

!     IF(Iteration.eq.500000) tNotConverged=.false.

        IF(tLagrange) THEN
            IF((abs(Force).lt.ConvergedForce).and.(abs(OrthoForce).lt.ConvergedForce)) THEN
                tNotConverged=.false.
            ENDIF
        ELSEIF(tROIteration) THEN
            IF(Iteration.eq.ROIterMax) THEN
                tNotConverged=.false.
            ENDIF
        ELSEIF(abs(TotCorrectedForce).lt.ConvergedForce) THEN
            tNotConverged=.false.
        ENDIF
! IF an ROIteration value is specified, use this to specify the end of the orbital rotation, otherwise use the 
! conversion limit (ConvergedForce).

    END SUBROUTINE TestForConvergence




    SUBROUTINE ShakeConstraints()
! DerivCoeff(k,a) is the unconstrained force on the original coefficients (CoeffT1(a,k)). 
        INTEGER :: w,l,a,m,ShakeIteration,ConvergeCount,SymM,SymMin
        real(dp) :: TotCorConstraints,TotConstraints,TotLambdas
        real(dp) :: TotUncorForce,TotDiffUncorCoeffs,TotDiffCorCoeffs
        LOGICAL :: tShakeNotConverged
        integer, save :: shake_io


!        WRITE(6,*) "Beginning shakeconstraints calculation"
        IF(Iteration.eq.1) THEN
            shake_io = get_free_unit()
            OPEN(shake_io,FILE='SHAKEstats',STATUS='unknown')
            WRITE(shake_io,'(A20,4A35,A20)') 'Shake Iteration','Sum Lambdas','Total of corrected forces', &
                & 'Sum unconstrained constraints',&
                                        &'Sum corrected constraints','Converge count'
        ENDIF
        IF(Mod(Iteration,10).eq.0) WRITE(shake_io,*) 'Orbital rotation iteration = ',Iteration


        ShakeIteration=0
        tShakeNotConverged=.true.

! Before we start iterating, take the current coefficients and find the derivative of the constraints with respect to them.

        CALL CalcDerivConstr(CoeffT1,DerivConstrT1)

!        WRITE(6,*) 'DerivContsrT1'
!            do l=1,TotNoConstraints
!                i=lab(1,l)
!                j=lab(2,l)
!                WRITE(6,*) i,j
!                do m=1,NoOrbs
!                    do a=1,NoOrbs
!                        WRITE(6,*) DerivConstrT1(a,m,l)
!                    enddo
!                enddo
!            enddo
!            stop

! Then find the coefficients at time t2, when moved by the completely unconstrained force and the values of the each 
! constraint at these positions.


        Correction(:,:)=0.0_dp
        CALL FindandUsetheForce(TotUncorForce,TotDiffUncorCoeffs,CoeffUncorT2)


        CALL CalcConstraints(CoeffUncorT2,Constraint,TotConstraints)

 
! Write stats from the beginning of the iteration to output.            
!        IF(Mod(Iteration,10).eq.0.or.Iteration.eq.1) THEN
!            CALL WriteShakeOUTstats01()
!        ENDIF

        CALL set_timer(Shake_Time,30)

        IF(tShakeDelay) THEN
            IF(Iteration.lt.ShakeStart) THEN
                ShakeIterMax=1
            ELSE
                ShakeIterMax=ShakeIterInput
            ENDIF
        ENDIF


! Actually starting the calculation.
        do while (tShakeNotConverged)

            ShakeIteration=ShakeIteration+1
            
            ForceCorrect(:,:)=0.0_dp          ! Zeroing terms that are re-calculated each iteration.
            CoeffCorT2(:,:)=0.0_dp
            ConstraintCor(:)=0.0_dp
            DerivConstrT2(:,:,:)=0.0_dp
            TotLambdas=0.0_dp
            TotCorrectedForce=0.0_dp
            TotDiffCorCoeffs=0.0_dp

            IF(ShakeIteration.ne.1) THEN
                CALL UpdateLambdas()
            ENDIF
            
            ShakeLambdaNew(:)=0.0_dp

            
            ! For a particular set of coefficients cm:
            ! Force(corrected)=Force(uncorrected)-Lambdas.DerivConstrT1
            ! Use these derivatives, and the current lambdas to find the trial corrected force.
            ! Then use this to get the (trial) shifted coefficients.

            ! Use the lambdas of this iteration to calculate the correction to the force due to the constraints.
            Correction(:,:)=0.0_dp
                

            do w=MinOccVirt,MaxOccVirt
                IF(w.eq.1) THEN
                    SymMin=1
                    MinMZ=1
                    IF(tSeparateOccVirt) THEN
                        MaxMZ=NoOcc
                    ELSE
                        MaxMZ=NoOrbs
                    ENDIF
                ELSE
                    SymMin=9
                    MinMZ=NoOcc+1
                    MaxMZ=NoOrbs
                ENDIF

                do m=MinMZ,MaxMZ
                    IF(tStoreSpinOrbs) THEN
                        SymM=INT(G1(SymLabelList2_rot(m))%sym%S)
                    ELSE
                        SymM=INT(G1(SymLabelList2_rot(m)*2)%sym%S)
                    ENDIF
                    do a=SymLabelCounts2_rot(1,SymM+SymMin), &
                            (SymLabelCounts2_rot(1,SymM+SymMin) + &
                                SymLabelCounts2_rot(2,SymM+SymMin)-1)
                        do l=1,TotNoConstraints
                            Correction(a,m)=Correction(a,m)+(ShakeLambda(l)*DerivConstrT1(a,m,l)) 
                        enddo
                    enddo
                enddo
            enddo

!            IF(tSeparateOccVirt) CALL ZeroOccVirtElements(Correction)


            CALL FindandUsetheForce(TotCorrectedForce,TotDiffCorCoeffs,CoeffCorT2)

! Use these new shifted coefficients to calculate the derivative of the constraints 
! (at time t2).
           
            CALL CalcDerivConstr(CoeffCorT2,DerivConstrT2) 
            
            
! Test for convergence, if convergence is reached, make the new coefficients the original ones to start the whole process again.
! Then exit out of this do loop and hence the subroutine.
            CALL TestShakeConvergence(ConvergeCount,TotCorConstraints,ShakeIteration,tShakeNotConverged)

! If the convergence criteria is met, exit out of this subroutine, a rotation has been made which keeps the coefficients 
! orthogonal.

! Write out stats of interest to output:
!            IF(Mod(Iteration,10).eq.0.or.Iteration.eq.1) THEN
!                CALL WriteShakeOUTstats02(ShakeIteration,TotLambdas,ConvergeCount) ! add correction to this.
!            ENDIF

! and to SHAKEstats file:
            CALL neci_flush(6)
            CALL neci_flush(shake_io)
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(shake_io,'(I20,4F35.20,I20)') ShakeIteration,TotLambdas,TotCorrectedForce,TotConstraints, &
                    TotCorConstraints,ConvergeCount 
            ENDIF

! If the convergence criteria is not met, use either the full matrix inversion method to find a new set of lambdas, or the shake algorithm 
! (in which case SHAKEAPPROX is required in the system block of the input).

            IF(tShakeApprox.and.tShakeNotConverged) THEN
!                WRITE(6,*) 'Using shake approximation to find new lambdas'
                CALL ShakeApproximation()
            ELSEIF(tShakeNotConverged) THEN
!                WRITE(6,*) 'Using the full diagonalisation shake method to find new lambdas'
                CALL FullShake()
            ELSE
                DistCs=TotDiffCorCoeffs
            ENDIF
   
    enddo

    CALL halt_timer(Shake_Time)

    ENDSUBROUTINE ShakeConstraints



    SUBROUTINE CalcDerivConstr(CurrCoeff,DerivConstr)
! This calculates the derivative of each of the orthonormalisation constraints, l, with respect
! to each set of coefficients cm.

        INTEGER :: l,i,j,a
        real(dp) :: CurrCoeff(NoOrbs,NoOrbs)
        real(dp) :: DerivConstr(NoOrbs,NoOrbs,TotNoConstraints)

        call set_timer(CalcDerivConstr_Time,30)


            DerivConstr(:,:,:)=0.0_dp
            do l=1,TotNoConstraints
                i=lab(1,l)
                j=lab(2,l)
                IF(i.eq.j) THEN
                    do a=1,NoOrbs
                        DerivConstr(a,i,l)=CurrCoeff(a,i)*2
                    enddo
                ELSE
                    do a=1,NoOrbs
                        DerivConstr(a,j,l)=CurrCoeff(a,i) 
                    enddo
                    do a=1,NoOrbs
                        DerivConstr(a,i,l)=CurrCoeff(a,j)
                    enddo
                ENDIF
                ! DerivConstrT1 stays the same throughout the iterations
            enddo

!            WRITE(6,*) 'DerivConstr'
!            do l=1,TotNoConstraints
!                i=Lab(1,l)
!                j=Lab(2,l)
!                WRITE(6,*) i,j
!                do m=1,NoOrbs
!                    do a=1,NoOrbs
!                        WRITE(6,'(F20.10)',advance='no') DerivConstr(a,m,l)
!                    enddo
!                    WRITE(6,*) ''
!                enddo
!            enddo

        call halt_timer(CalcDerivConstr_Time)

    ENDSUBROUTINE CalcDerivConstr



    SUBROUTINE FindandUsetheForce(TotForce,TotDiffCoeffs,CoeffT2)
! This takes the current lambdas with the derivatives of the constraints and calculates a force
! for each cm, with an orthonormalisation correction.
! This is then used to rotate the coefficients by a defined timestep.
        INTEGER :: a,m,Symm,w,SymMin,TempMaxOccVirt
        real(dp) :: TotForce,TotDiffCoeffs,CoeffT2(NoOrbs,NoOrbs)

!        WRITE(6,*) 'DerivCoeff'
!        do m=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(4F20.10)',advance='no') DerivCoeff(a,m)
!                IF(tStoreSpinOrbs) THEN
!                    WRITE(6,*) a,m,DerivCoeff(2*a,2*m)
!                ELSE
!                    WRITE(6,*) a,m,DerivCoeff(a,m)
!                ENDIF
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'Correction'
!        do m=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(4F20.10)',advance='no') Correction(a,m)
!                IF(tStoreSpinOrbs) THEN
!                    WRITE(6,*) a,m,Correction(2*a,2*m)
 !               ELSE
!                    WRITE(6,*) a,m,Correction(a,m)
!                ENDIF
!            enddo
!            WRITE(6,*) ''
!        enddo


        CALL set_timer(findandusetheforce_time,30)

        IF(tSeparateOccVirt) THEN
            TempMaxOccVirt=2
        ELSE
            TempMaxOccVirt=1
        ENDIF

!        do w=MinOccVirt,MaxOccVirt
        do w=1,TempMaxOccVirt
! the force will be zero on those coefficients not being mixed, but still want to run over all, so that the diagonal 1 values are maintained.
            IF(w.eq.1) THEN
                SymMin=1
                MinMZ=1
                IF(tSeparateOccVirt) THEN
                    MaxMZ=NoOcc
                ELSE
                    MaxMZ=NoOrbs
                ENDIF
            ELSE
                SymMin=9
                MinMZ=NoOcc+1
                MaxMZ=NoOrbs
            ENDIF

            do m=MinMZ,MaxMZ
                IF(tStoreSpinOrbs) THEN
                    SymM=INT(G1(SymLabelList2_rot(m))%sym%S)
                ELSE
                    SymM=INT(G1(SymLabelList2_rot(m)*2)%sym%S)
                ENDIF
                do a=SymLabelCounts2_rot(1,SymM+SymMin), &
                        (SymLabelCounts2_rot(1,SymM+SymMin) + &
                            SymLabelCounts2_rot(2,SymM+SymMin)-1)
!               
                ! FIND THE FORCE 
                    ! find the corrected force. (in the case where the uncorrected force is required, correction is set to 0.
                    ! DerivCoeff(m,a) is the derivative of the relevant potential energy w.r.t cm without any constraints (no lambda terms).
                    ! ForceCorrect is then the latest force on coefficients.  This is iteratively being corrected so that
                    ! it will finally move the coefficients so that they remain orthonormal.
                
                ! USE THE FORCE
                    ForceCorrect(a,m)=DerivCoeff(a,m)-Correction(a,m)
                    CoeffT2(a,m)=CoeffT1(a,m)-(TimeStep*ForceCorrect(a,m))
                    ! Using the force to calculate the coefficients at time T2 (hopefully more orthonomal than those calculated in
                    ! the previous iteration).
                    
                    ! Calculate parameters for printing
                    TotForce=TotForce+ABS(ForceCorrect(a,m))
                    TotDiffCoeffs=TotDiffCoeffs+ABS(CoeffT2(a,m)-CoeffT1(a,m))
                enddo
            enddo
        enddo


        TotForce=TotForce/(REAL(NoOrbs**2,dp))


!        WRITE(6,*) 'ForceCorrect'
!        do m=1,NoOrbs
!            do a=1,NoOrbs
!                WRITE(6,'(4F20.10)',advance='no') ForceCorrect(a,m)
!            enddo
!            WRITE(6,*) ''
!        enddo


        call halt_timer(findandusetheforce_time)

    ENDSUBROUTINE FindandUsetheForce



    SUBROUTINE CalcConstraints(CurrCoeff,Constraint,TotConstraints)  
! This calculates the value of each orthonomalisation constraint, using the shifted coefficients.
! Each of these should tend to 0 when the coefficients become orthonomal.
        INTEGER :: l,i,j
        real(dp) :: CurrCoeff(NoOrbs,NoOrbs),TotConstraints,Constraint(TotNoConstraints) 


            TotConstraints=0.0_dp
            do l=1,TotNoConstraints
                i=lab(1,l)
                j=lab(2,l)
                IF(i.eq.j) THEN
                    Constraint(l)=Dot_Product(CurrCoeff(:,i),CurrCoeff(:,j))-1.0_dp
                ELSE
                    Constraint(l)=Dot_Product(CurrCoeff(:,i),CurrCoeff(:,j))
                    ! Each of these components should tend towards 0 when the coefficients become orthonormal.
                ENDIF
                TotConstraints=TotConstraints+ABS(Constraint(l))
            enddo
     
!            do i=1,NoOrbs
!                do j=1,NoOrbs
!                    WRITE(6,'(I3)') i,j
!                    WRITE(6,'(F20.10)') CurrCoeff(i,j)
!                enddo
!                WRITE(6,*) ''
!            enddo
  
!            do l=1,TotNoConstraints
!                WRITE(6,*) Constraint(l)
!            enddo


    ENDSUBROUTINE CalcConstraints


    SUBROUTINE FullShake()
! This method calculates the lambdas by solving the full matrix equation.
        INTEGER :: l,n,m,info,ipiv(TotNoConstraints)
        CHARACTER(len=*), PARAMETER :: this_routine='FullShake'


        CALL set_timer(FullShake_Time,30)

! FULL MATRIX INVERSION METHOD

! Calculate matrix from the derivatives of the constraints w.r.t the the coefficients at t1 and t2. I.e. the initial 
! coefficients and those that have been moved by the corrected force.

            DerivConstrT1T2(:,:)=0.0_dp
            do l=1,TotNoConstraints
                do n=1,TotNoConstraints
                    do m=1,NoOrbs
                        ! Product of constraint i,j at time t1, mult by constraint l,n.
                        ! Add these over all m for a specific constraints to get matrix elements
                        DerivConstrT1T2(n,l)=DerivConstrT1T2(n,l)+(Dot_Product(DerivConstrT2(:,m,l),DerivConstrT1(:,m,n)))
                    enddo
                enddo
            enddo       ! have filled up whole matrix

!            do l=1,TotNoConstraints
!                do n=1,TotNoConstraints
!                    WRITE(6,'(F20.10)',advance='no') DerivConstrT1T2(n,l)
!                enddo
!                WRITE(6,*) ''
!            enddo

! Invert the matrix to calculate the lambda values.
! LU decomposition.
            call dgetrf(TotNoConstraints,TotNoConstraints,DerivConstrT1T2,TotNoConstraints,ipiv,info)
            if(info.ne.0) THEN
                WRITE(6,*) 'info ',info
                CALL Stop_All(this_routine,"The LU decomposition of matrix inversion failed...")
            endif

!            do n=1,TotNoConstraints
!                WRITE(6,*) Constraint(n)
!            enddo

            do n=1,TotNoConstraints
                ShakeLambdaNew(n)=Constraint(n)/(TimeStep*(-1))
            enddo
            ! These are actually still the constraint values, but now Lambda(n) can go into dgetrs as the constraints (B in AX=B), 
            ! and come out as the computed lambdas (X).

            call dgetrs('N',TotNoConstraints,1,DerivConstrT1T2,TotNoConstraints,ipiv,ShakeLambdaNew,TotNoConstraints,info)
            if(info.ne.0) CALL Stop_All(this_routine,"Error in dgetrs, solving for the lambdas...")


!            WRITE(6,*) 'Lambdas successfully calculated, beginning next shake iteration'
!            do n=1,TotNoConstraints
!                WRITE(6,*) ShakeLambdaNew(n)
!            enddo


        CALL halt_timer(FullShake_Time)


    ENDSUBROUTINE FullShake
 


    SUBROUTINE ShakeApproximation()
! This is an approximation in which only the diagonal elements are considered in the 
! matrix of the derivative of the constraints DerivConstrT1T2.
        INTEGER :: m,l


! Use 'shake' algorithm in which the iterative scheme is applied to each constraint in succession.
            WRITE(6,*) 'DerivConstrT1T2Diag calculated from the shake approx'
            
            DerivConstrT1T2Diag(:)=0.0_dp
            do l=1,TotNoConstraints 
                do m=1,NoOrbs
                    DerivConstrT1T2Diag(l)=DerivConstrT1T2Diag(l)+Dot_Product(DerivConstrT2(:,m,l),DerivConstrT1(:,m,l))
                enddo
                ShakeLambdaNew(l)=Constraint(l)/((-1)*TimeStep*DerivConstrT1T2Diag(l))
                WRITE(6,*) DerivConstrT1T2Diag(l)
            enddo


    ENDSUBROUTINE ShakeApproximation


    SUBROUTINE UpdateLambdas()
! Use damping to update the lambdas, rather than completely replacing them with the new values.
        INTEGER :: l
        

        do l=1,TotNoConstraints
            ShakeLambda(l)=ShakeLambdaNew(l)
        enddo

! DAMPING
!        do l=1,TotNoConstraints
!            ShakeLambda(l)=(0.9*ShakeLambda(l))+(0.1*ShakeLambdaNew(l))
!        enddo
! If decide to use this, make the 0.9 value a damping parameter in the input.


    ENDSUBROUTINE UpdateLambdas



    SUBROUTINE TestShakeConvergence(ConvergeCount,TotCorConstraints,ShakeIteration,tShakeNotConverged)  
! This calculates the value of each orthonomalisation constraint using the corrected coefficients.
! Each of these should tend to 0 when the coefficients become orthonomal.
! CovergeCount counts the number of constraints that individually have values below the specified
! convergence criteria.  If this = 0, the shake is converged, else keep iterating.
        INTEGER :: l,i,j,m,a,ConvergeCount
        real(dp) :: TotCorConstraints
        INTEGER :: ShakeIteration
        LOGICAL :: tShakeNotConverged


        TotCorConstraints=0.0_dp
        ConvergeCount=0
        ConstraintCor(:)=0.0_dp
        do l=1,TotNoConstraints
            i=lab(1,l)
            j=lab(2,l)
            IF(i.eq.j) THEN
                ConstraintCor(l)=Dot_Product(CoeffCorT2(:,i),CoeffCorT2(:,j))-1.0_dp
            ELSE
                ConstraintCor(l)=Dot_Product(CoeffCorT2(:,i),CoeffCorT2(:,j))
                ! Each of these components should tend towards 0 when the coefficients become orthonormal.
            ENDIF
            
            TotCorConstraints=TotCorConstraints+ABS(ConstraintCor(l))
            ! Sum of all Contraint components - indication of overall orthonormality.
    
            IF(ABS(ConstraintCor(l)).gt.ShakeConverged) ConvergeCount=ConvergeCount+1
            ! Count the number of constraints which are still well above 0.
            
        enddo
        
        IF(tShakeIter) THEN
            IF(ShakeIteration.eq.ShakeIterMax) THEN
                do m=1,NoOrbs
                    do a=1,NoOrbs
                        CoeffT1(a,m)=CoeffCorT2(a,m)
                    enddo
                enddo
!                WRITE(6,*) 'stopped at iteration, ',ShakeIteration
                tShakeNotConverged=.false.
            ENDIF
        ELSEIF(ConvergeCount.eq.0) THEN
           tShakeNotConverged=.false.
!           WRITE(6,*) 'Convergence reached in the shake algorithm'
!           WRITE(6,*) 'All constraints have values less than ',ShakeConverged

! If convergence is reached, make the new coefficients coeff, to start the rotation iteration again.

            do m=1,NoOrbs
                do a=1,NoOrbs
                    CoeffT1(a,m)=CoeffCorT2(a,m)
                enddo
            enddo
        ENDIF


    ENDSUBROUTINE TestShakeConvergence

  

    SUBROUTINE FinalizeNewOrbs()
! At the end of the orbital rotation, have a set of coefficients CoeffT1 which transform the HF orbitals into a set of linear
! combinations ui which minimise |<ij|kl>|^2.  This is the final subroutine after all iterations (but before the memory deallocation)
! that calculates the final 4 index integrals to be used in the NECI calculation.
        use sym_mod, only: GenSymStatePairs
        INTEGER :: i,a,j
        real(dp) :: TotGSConstraints,GSConstraint(TotNoConstraints),CoeffTemp(SpatOrbs,SpatOrbs)
        
!        WRITE(6,*) 'The final transformation coefficients before gram schmidt orthonormalisation'
!        do i=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(F20.10)',advance='no') CoeffT1(a,i)
!                IF(tStoreSpinOrbs) THEN
!                    WRITE(6,*) a,i,CoeffT1(2*a,2*i)
!                ELSE
!                    WRITE(6,*) a,i,CoeffT1(a,i)
!                ENDIF
!            enddo
!            WRITE(6,*) ''
!        enddo
!        stop

!        WRITE(6,*) 'The final values of the constraints with corrected coefficients'
!        do i=1,TotNoConstraints
!            WRITE(6,*) ConstraintCor(i)
!        enddo

! First need to do a final explicit orthonormalisation.  The orbitals are very close to being orthonormal, but not exactly.
! Need to make sure they are exact orthonormal using Gram Schmit.
        IF(tStoreSpinOrbs.and.(.not.tMaxHLGap)) THEN
            CoeffTemp(:,:)=0.0_dp
            do i=1,SpatOrbs
                do j=1,SpatOrbs
                    CoeffTemp(i,j)=CoeffT1(2*i,2*j)
                enddo
            enddo

            CALL GRAMSCHMIDT(CoeffTemp,SpatOrbs)

            CoeffT1(:,:)=0.0_dp
            do i=1,SpatOrbs
                do j=1,SpatOrbs
                    CoeffT1(2*i,2*j)=CoeffTemp(i,j)
                    CoeffT1((2*i)-1,(2*j)-1)=CoeffTemp(i,j)
                enddo
            enddo
        ELSEIF(.not.tMaxHLGap) THEN
            CALL GRAMSCHMIDT(CoeffT1,NoOrbs)
        ENDIF
        

! Put routine in here that takes this rotation matrix, CoeffT1, and forms raises it to the power of a small number, alpha.
! Changeing this number allows us to see the change in plateau level with various rotations.


! Write out some final results of interest, like values of the constraints, values of new coefficients.
    
        WRITE(6,*) 'The final transformation coefficients after gram schmidt orthonormalisation'
        do i=1,NoOrbs
            do a=1,NoOrbs
                WRITE(6,'(F10.4)',advance='no') CoeffT1(a,i)
            enddo
            WRITE(6,*) ''
        enddo

        CALL WriteTransformMat()
        
        CALL CalcConstraints(CoeffT1,GSConstraint,TotGSConstraints)  


!        WRITE(6,*) 'The values of the constraints after gram schmidt orthonormalisation'
        
!        do i=1,TotNoConstraints
!            WRITE(6,*) GSConstraint(i)
!        enddo

        
        WRITE(6,*) 'Final Potential Energy before orthogonalisation',PotEnergy

        CALL Transform2ElInts()
            

! Use these final coefficients to find the FourIndInts(i,j,k,l).
! These are now the <ij|kl> integrals we now want to use instead of the HF UMat.
! New potential energy is calculated in this routine using the orthogonalised coefficients.
! Compare to that before this, to make sure the orthogonalisation hasn't shifted them back to a non-minimal place.

        WRITE(6,*) 'Final Potential Energy after orthogonalisation',PotEnergy

! Calculate the fock matrix, and print it out to see how much the off diagonal terms contribute.
! Also print out the sum of the diagonal elements to compare to the original value.
        CALL CalcFOCKMatrix()

        CALL RefillUMATandTMAT2D()        
! UMat is the 4 index integral matrix (2 electron), whereas TMAT2D is the 2 index integral (1 el) matrix
   
! This is the keyword that tells the NECI calculation that the orbitals are not HF.  It means that contributions to
! the energy from walkers on singly occupied determinants are included in the values printed.
! Making it true here allows us to go directly from a Rotation into a spawn if required.
        tRotatedOrbs=.true.

        CALL GENSymStatePairs(SpatOrbs,.false.)
        

    ENDSUBROUTINE FinalizeNewOrbs


    SUBROUTINE WriteSingHisttofile()
        INTEGER :: i,j,k,BinNo,a,b,iunit
        real(dp) :: MaxFII,MinFII,BinIter,BinVal,SingExcit(NoOrbs,NoOrbs)


!<ik|jk> terms where all i,j and k are virtual
        !Coulomb
        ROHistSCijkVir(:,:)=0.0_dp
        MaxFII=FourIndInts(NoOcc+1,NoOcc+1,NoOcc+2,NoOcc+1)
        MinFII=FourIndInts(NoOcc+1,NoOcc+1,NoOcc+2,NoOcc+1)
        do i=NoOcc+1,NoOrbs
            do k=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,j,k).gt.MaxFII) MaxFII=FourIndInts(i,k,j,k)
                    IF(FourIndInts(i,k,j,k).lt.MinFII) MinFII=FourIndInts(i,k,j,k)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSCijkVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=NoOcc+1,NoOrbs
            do k=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,j,k).ne.0.0_dp) THEN
                        BinNo=CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCijkVir(2,BinNo)=ROHistSCijkVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Exchange
        ROHistSEijkVir(:,:)=0.0_dp
        MaxFII=FourIndInts(NoOcc+1,NoOcc+1,NoOcc+1,NoOcc+2)
        MinFII=FourIndInts(NoOcc+1,NoOcc+1,NoOcc+1,NoOcc+2)
        do i=NoOcc+1,NoOrbs
            do k=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,k,j).gt.MaxFII) MaxFII=FourIndInts(i,k,k,j)
                    IF(FourIndInts(i,k,k,j).lt.MinFII) MinFII=FourIndInts(i,k,k,j)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSEijkVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=NoOcc+1,NoOrbs
            do k=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,k,j).ne.0.0_dp) THEN
                        BinNo=CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEijkVir(2,BinNo)=ROHistSEijkVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !antisymmetric
        ROHistSASijkVir(:,:)=0.0_dp
        MaxFII=FourIndInts(NoOcc+1,NoOcc+1,NoOcc+2,NoOcc+1)-FourIndInts(NoOcc+1,NoOcc+1,NoOcc+1,NoOcc+2)
        MinFII=FourIndInts(NoOcc+1,NoOcc+1,NoOcc+2,NoOcc+1)-FourIndInts(NoOcc+1,NoOcc+1,NoOcc+1,NoOcc+2)
        do i=NoOcc+1,NoOrbs
            do k=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).gt.MaxFII) MaxFII=FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).lt.MinFII) MinFII=FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSASijkVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=NoOcc+1,NoOrbs
            do k=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).ne.0.0_dp) THEN
                        BinNo=CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASijkVir(2,BinNo)=ROHistSASijkVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo


        IF(Iteration.eq.0) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistHFSingijkVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCijkVir(2,j).ne.0).or.(ROHistSEijkVir(2,j).ne.0).or.(ROHistSASijkVir(2,j).ne.0)) THEN
                    WRITE(iunit,'(6F20.10)') ROHistSCijkVir(1,j),ROHistSCijkVir(2,j),ROHistSEijkVir(1,j),ROHistSEijkVir(2,j),&
                                                        &ROHistSASijkVir(1,j),ROHistSASijkVir(2,j)
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistRotSingijkVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCijkVir(2,j).ne.0).or.(ROHistSEijkVir(2,j).ne.0).or.(ROHistSASijkVir(2,j).ne.0)) THEN
                    WRITE(iunit,'(6F20.10)') ROHistSCijkVir(1,j),ROHistSCijkVir(2,j),ROHistSEijkVir(1,j),ROHistSEijkVir(2,j),&
                                                        &ROHistSASijkVir(1,j),ROHistSASijkVir(2,j)
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF


!<ik|jk> where k is occupied, and i and j are both virtual
        !Coulomb
        ROHistSCkOcijVir(:,:)=0.0_dp
        MaxFII=FourIndInts(NoOcc+1,1,NoOcc+2,1)
        MinFII=FourIndInts(NoOcc+1,1,NoOcc+2,1)
        do i=NoOcc+1,NoOrbs
            do k=1,NoOcc
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,j,k).gt.MaxFII) MaxFII=FourIndInts(i,k,j,k)
                    IF(FourIndInts(i,k,j,k).lt.MinFII) MinFII=FourIndInts(i,k,j,k)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSCkOcijVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=NoOcc+1,NoOrbs
            do k=1,NoOcc
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,j,k).ne.0.0_dp) THEN
                        BinNo=CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCkOcijVir(2,BinNo)=ROHistSCkOcijVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Exchange
        ROHistSEkOcijVir(:,:)=0.0_dp
        MaxFII=FourIndInts(NoOcc+1,1,1,NoOcc+2)
        MinFII=FourIndInts(NoOcc+1,1,1,NoOcc+2)
        do i=NoOcc+1,NoOrbs
            do k=1,NoOcc
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,k,j).gt.MaxFII) MaxFII=FourIndInts(i,k,k,j)
                    IF(FourIndInts(i,k,k,j).lt.MinFII) MinFII=FourIndInts(i,k,k,j)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSEkOcijVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=NoOcc+1,NoOrbs
            do k=1,NoOcc
                do j=i+1,NoOrbs
                    IF(FourIndInts(i,k,k,j).ne.0.0_dp) THEN
                        BinNo=CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEkOcijVir(2,BinNo)=ROHistSEkOcijVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !antisymmetric 
        ROHistSASkOcijVir(:,:)=0.0_dp
        MaxFII=FourIndInts(NoOcc+1,1,NoOcc+2,1)-FourIndInts(NoOcc+1,1,1,NoOcc+2)
        MinFII=FourIndInts(NoOcc+1,1,NoOcc+2,1)-FourIndInts(NoOcc+1,1,1,NoOcc+2)
        do i=NoOcc+1,NoOrbs
            do k=1,NoOcc
                do j=i+1,NoOrbs
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).gt.MaxFII) MaxFII=FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).lt.MinFII) MinFII=FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSASkOcijVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=NoOcc+1,NoOrbs
            do k=1,NoOcc
                do j=i+1,NoOrbs
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).ne.0.0_dp) THEN
                        BinNo=CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASkOcijVir(2,BinNo)=ROHistSASkOcijVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo


        IF(Iteration.eq.0) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistHFSingkOcijVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCkOcijVir(2,j).ne.0).or.(ROHistSEkOcijVir(2,j).ne.0).or.(ROHistSASkOcijVir(2,j).ne.0)) THEN
                  WRITE(iunit,'(6F20.10)') ROHistSCkOcijVir(1,j),ROHistSCkOcijVir(2,j),ROHistSEkOcijVir(1,j),ROHistSEkOcijVir(2,j),&
                                                        &ROHistSASkOcijVir(1,j),ROHistSASkOcijVir(2,j)
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistRotSingkOcijVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCkOcijVir(2,j).ne.0).or.(ROHistSEkOcijVir(2,j).ne.0).or.(ROHistSASkOcijVir(2,j).ne.0)) THEN
                  WRITE(iunit,'(6F20.10)') ROHistSCkOcijVir(1,j),ROHistSCkOcijVir(2,j),ROHistSEkOcijVir(1,j),ROHistSEkOcijVir(2,j),&
                                                        &ROHistSASkOcijVir(1,j),ROHistSASkOcijVir(2,j)
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF


! <ik|jk> where i and k are both occupied, and j virtual.
        ! Coulomb
        ROHistSCikOcjVir(:,:)=0.0_dp
        MaxFII=FourIndInts(1,1,NoOcc+1,1)
        MinFII=FourIndInts(1,1,NoOcc+1,1)
        do i=1,NoOcc
            do k=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(FourIndInts(i,k,j,k).gt.MaxFII) MaxFII=FourIndInts(i,k,j,k)
                    IF(FourIndInts(i,k,j,k).lt.MinFII) MinFII=FourIndInts(i,k,j,k)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSCikOcjVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=1,NoOcc
            do k=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(FourIndInts(i,k,j,k).ne.0.0_dp) THEN
                        BinNo=CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCikOcjVir(2,BinNo)=ROHistSCikOcjVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Exchange 
        ROHistSEikOcjVir(:,:)=0.0_dp
        MaxFII=FourIndInts(1,1,1,NoOcc+1)
        MinFII=FourIndInts(1,1,1,NoOcc+1)
        do i=1,NoOcc
            do k=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(FourIndInts(i,k,k,j).gt.MaxFII) MaxFII=FourIndInts(i,k,k,j)
                    IF(FourIndInts(i,k,k,j).lt.MinFII) MinFII=FourIndInts(i,k,k,j)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSEikOcjVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=1,NoOcc
            do k=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(FourIndInts(i,k,k,j).ne.0.0_dp) THEN
                        BinNo=CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEikOcjVir(2,BinNo)=ROHistSEikOcjVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Antisymmetrised
        ROHistSASikOcjVir(:,:)=0.0_dp
        MaxFII=FourIndInts(1,1,NoOcc+1,1)-FourIndInts(1,1,1,NoOcc+1)
        MinFII=FourIndInts(1,1,NoOcc+1,1)-FourIndInts(1,1,1,NoOcc+1)
        do i=1,NoOcc
            do k=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).gt.MaxFII) MaxFII=FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).lt.MinFII) MinFII=FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)
                enddo
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSASikOcjVir(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do i=1,NoOcc
            do k=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).ne.0.0_dp) THEN
                        BinNo=CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASikOcjVir(2,BinNo)=ROHistSASikOcjVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        IF(Iteration.eq.0) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistHFSingikOcjVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCikOcjVir(2,j).ne.0).or.(ROHistSEikOcjVir(2,j).ne.0).or.(ROHistSASikOcjVir(2,j).ne.0)) THEN 
                  WRITE(iunit,'(6F20.10)') ROHistSCikOcjVir(1,j),ROHistSCikOcjVir(2,j),ROHistSEikOcjVir(1,j),ROHistSEikOcjVir(2,j),&
                                                        &ROHistSASikOcjVir(1,j),ROHistSASikOcjVir(2,j)
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistRotSingikOcjVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCikOcjVir(2,j).ne.0).or.(ROHistSEikOcjVir(2,j).ne.0).or.(ROHistSASikOcjVir(2,j).ne.0)) THEN 
                  WRITE(iunit,'(6F20.10)') ROHistSCikOcjVir(1,j),ROHistSCikOcjVir(2,j),ROHistSEikOcjVir(1,j),ROHistSEikOcjVir(2,j),&
                                                        &ROHistSASikOcjVir(1,j),ROHistSASikOcjVir(2,j)
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF

!Single excitations connected to the HF determinant.
        ROHistSing(:,:)=0.0_dp
        MaxFII=0.0_dp
        MinFII=0.0_dp
        do j=NoOcc+1,NoOrbs
            do i=1,NoOcc
                SingExcit(i,j)=0.0_dp
                IF(i.eq.j) CYCLE
                a=SymLabelList2_rot(i)
                b=SymLabelList2_rot(j)
                do k=1,NoOcc+1
!                    IF(k.eq.j) CYCLE
!                    IF(k.eq.i) CYCLE
                    SingExcit(i,j)=SingExcit(i,j)+REAL(TMAT2D(2*a,2*b),dp)+((2*FourIndInts(i,k,j,k))-FourIndInts(i,k,k,j))
                enddo
                IF(SingExcit(i,j).gt.MaxFII) MaxFII=SingExcit(i,j)
                IF(SingExcit(i,j).lt.MinFII) MinFII=SingExcit(i,j)
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.0_dp
        MaxFII=MaxFII+BinIter
        MinFII=MinFII-BinIter
        BinVal=MinFII
        do i=1,4002
            ROHistSing(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do j=NoOcc+1,NoOrbs
            do i=1,NoOcc
                IF(i.eq.j) CYCLE
                IF(SingExcit(i,j).ne.0.0_dp) THEN
                    BinNo=CEILING((SingExcit(i,j)-MinFII)*4002/(MaxFII-MinFII))
                    ROHistSing(2,BinNo)=ROHistSing(2,BinNo)+1.0         
                ENDIF
            enddo
        enddo

        IF(Iteration.eq.0) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistHFSingExcHF',STATUS='unknown')
            do j=1,4002
                IF(ROHistSing(2,j).ne.0) THEN
                    do i=1,2
                        WRITE(iunit,'(F20.10)',advance='no') ROHistSing(i,j)
                    enddo
                    WRITE(iunit,*) ''
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE='HistRotSingExcHF',STATUS='unknown')
            do j=1,4002
                IF(ROHistSing(2,j).ne.0) THEN
                    do i=1,2
                        WRITE(iunit,'(F20.10)',advance='no') ROHistSing(i,j)
                    enddo
                    WRITE(iunit,*) ''
                ENDIF
            enddo
            CLOSE(iunit)
        ENDIF



    ENDSUBROUTINE WriteSingHisttofile 




    SUBROUTINE WriteDoubHisttofile()
        INTEGER :: i,j,k,l,BinNo, iunit
        real(dp) :: MaxFII,MinFII,BinIter,OnePartOrbEnValue,BinVal


!        OPEN(34,FILE='FourIndInts',STATUS='unknown')
!        WRITE(34,'(A19,A20,A19,A20)') 'i,j,k,l','','i,j,l,k',''
!        do l=1,NoOrbs
!            do k=1,l
!                do j=1,NoOrbs
!                    do i=1,NoOrbs
!                        WRITE(34,'(I10,A1,I2,A1,I2,A1,I2,F20.10,I10,A1,I2,A1,I2,A1,I2,F20.10)') i,',',j,',',k,',',l,FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k),&
!                        &i,',',j,',',l,',',k,FourIndInts(i,j,l,k)-FourIndInts(i,j,k,l)
!                    enddo
!                enddo
!            enddo
!        enddo
!        CLOSE(34)

! Histogramming all coulomb terms <ij|ij> where i<j, and i and j are both virtual.
! In reality we are looking at i=<j, but the ERhistograms will show the i=j terms.
        IF(tROHistVirtCoulomb) THEN
        
            ROHistDCijOcklVir(:,:)=0.0_dp
            MinFII=FourIndInts(1,2,NoOcc+1,NoOcc+2)
            MaxFII=FourIndInts(1,2,NoOcc+1,NoOcc+2)
            do i=1,NoOcc
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do l=NoOcc+1,NoOrbs
                            IF(FourIndInts(i,j,k,l).lt.MinFII) MinFII=FourIndInts(i,j,k,l)
                            IF(FourIndInts(i,j,k,l).gt.MaxFII) MaxFII=FourIndInts(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistDCijOcklVir(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOcc
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do l=NoOcc+1,NoOrbs
                            IF(FourIndInts(i,j,k,l).ne.0) THEN
                                BinNo=CEILING((FourIndInts(i,j,k,l)-MinFII)*4002/(MaxFII-MinFII))
                                ROHistDCijOcklVir(2,BinNo)=ROHistDCijOcklVir(2,BinNo)+1.0         
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo

            !antisymmetric
            ROHistASijOcklVir(:,:)=0.0_dp
            MinFII=FourIndInts(1,2,NoOcc+1,NoOcc+2)-FourIndInts(1,2,NoOcc+2,NoOcc+1)
            MaxFII=FourIndInts(1,2,NoOcc+1,NoOcc+2)-FourIndInts(1,2,NoOcc+2,NoOcc+1)
            do i=1,NoOcc
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do l=NoOcc+1,NoOrbs
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).lt.MinFII) MinFII=(FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).gt.MaxFII) MaxFII=(FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistASijOcklVir(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOcc
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do l=NoOcc+1,NoOrbs
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).ne.0) THEN
                                BinNo=CEILING(((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))-MinFII)*4002/(MaxFII-MinFII))
                                ROHistASijOcklVir(2,BinNo)=ROHistASijOcklVir(2,BinNo)+1.0         
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFDoubijOcklVir',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijOcklVir(2,j).ne.0).or.(ROHistASijOcklVir(2,j).ne.0)) THEN
                        WRITE(iunit,'(4F20.10)') ROHistDCijOcklVir(1,j),ROHistDCijOcklVir(2,j), &
                            ROHistASijOcklVir(1,j),ROHistASijOcklVir(2,j)
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotDoubijOcklVir',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijOcklVir(2,j).ne.0).or.(ROHistASijOcklVir(2,j).ne.0)) THEN
                        WRITE(iunit,'(4F20.10)') ROHistDCijOcklVir(1,j),ROHistDCijOcklVir(2,j), &
                            ROHistASijOcklVir(1,j),ROHistASijOcklVir(2,j)
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF



            ROHistDCijklVir(:,:)=0.0_dp
            MinFII=FourIndInts(NoOrbs-1,NoOrbs,NoOrbs-1,NoOrbs)
            MaxFII=FourIndInts(NoOrbs-1,NoOrbs,NoOrbs-1,NoOrbs)
            do i=NoOcc+1,NoOrbs
                do j=NoOcc+1,NoOrbs
                    do k=i+1,NoOrbs
                        do l=j+1,NoOrbs
                            IF(FourIndInts(i,j,k,l).lt.MinFII) MinFII=FourIndInts(i,j,k,l)
                            IF(FourIndInts(i,j,k,l).gt.MaxFII) MaxFII=FourIndInts(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistDCijklVir(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=NoOcc+1,NoOrbs
                do j=NoOcc+1,NoOrbs
                    do k=i+1,NoOrbs
                        do l=j+1,NoOrbs
                            IF(FourIndInts(i,j,k,l).ne.0) THEN
                                BinNo=CEILING((FourIndInts(i,j,k,l)-MinFII)*4002/(MaxFII-MinFII))
                                ROHistDCijklVir(2,BinNo)=ROHistDCijklVir(2,BinNo)+1.0         
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo

            !antisymmetric
            ROHistASijklVir(:,:)=0.0_dp
            MinFII=FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs-1,NoOrbs)-FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs,NoOrbs-1)
            MaxFII=FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs-1,NoOrbs)-FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs,NoOrbs-1)
            do i=NoOcc+1,NoOrbs
                do j=NoOcc+1,NoOrbs
                    do k=i+1,NoOrbs
                        do l=j+1,NoOrbs
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).lt.MinFII) MinFII=(FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).gt.MaxFII) MaxFII=(FourIndInts(i,j,k,l)- &
                                FourIndInts(i,j,l,k))
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistASijklVir(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=NoOcc+1,NoOrbs
                do j=NoOcc+1,NoOrbs
                    do k=i+1,NoOrbs
                        do l=j+1,NoOrbs
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).ne.0) THEN
                                BinNo=CEILING(((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))-MinFII)*4002/(MaxFII-MinFII))
                                ROHistASijklVir(2,BinNo)=ROHistASijklVir(2,BinNo)+1.0         
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFDoubijklVirt',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijklVir(2,j).ne.0).or.(ROHistASijklVir(2,j).ne.0)) THEN
                        WRITE(iunit,'(4F20.10)') ROHistDCijklVir(1,j),ROHistDCijklVir(2,j), &
                            ROHistASijklVir(1,j),ROHistASijklVir(2,j)
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotDoubijklVirt',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijklVir(2,j).ne.0).or.(ROHistASijklVir(2,j).ne.0)) THEN
                        WRITE(iunit,'(4F20.10)') ROHistDCijklVir(1,j),ROHistDCijklVir(2,j),ROHistASijklVir(1,j), &
                            ROHistASijklVir(2,j)
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
        ENDIF
   


! Histogramming all one particle orbital energies (occupied and virtual) even though we are not changing occupied.  Would like to see HOMO-LUMO gap etc.
        IF(tROHistOneElInts) THEN

            ROHistHijVirt(:,:)=0.0_dp
            MinFII=TMAT2DRot(NoOcc+1,NoOcc+2)
            MaxFII=TMAT2DRot(NoOcc+1,NoOcc+2)
            do i=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(TMAT2DRot(i,j).lt.MinFII) MinFII=TMAT2DRot(i,j)
                    IF(TMAT2DRot(i,j).gt.MaxFII) MaxFII=TMAT2DRot(i,j)
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistHijVirt(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(TMAT2DRot(i,j).ne.0) THEN
                        BinNo=CEILING((TMAT2DRot(i,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistHijVirt(2,BinNo)=ROHistHijVirt(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFHijVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistHijVirt(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotHijVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistHijVirt(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
 


            ROHistHijOccVirt(:,:)=0.0_dp
            MinFII=TMAT2DRot(1,NoOcc+1)
            MaxFII=TMAT2DRot(1,NoOcc+1)
            do i=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(TMAT2DRot(i,j).lt.MinFII) MinFII=TMAT2DRot(i,j)
                    IF(TMAT2DRot(i,j).gt.MaxFII) MaxFII=TMAT2DRot(i,j)
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistHijOccVirt(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(TMAT2DRot(i,j).ne.0) THEN
                        BinNo=CEILING((TMAT2DRot(i,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistHijOccVirt(2,BinNo)=ROHistHijOccVirt(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFHijOccVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijOccVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistHijOccVirt(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotHijOccVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijOccVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistHijOccVirt(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
 


            ROHistHii(:,:)=0.0_dp
            MinFII=TMAT2DRot(1,1)
            MaxFII=TMAT2DRot(1,1)
            do i=1,NoOrbs
                IF(TMAT2DRot(i,i).lt.MinFII) MinFII=TMAT2DRot(i,i)
                IF(TMAT2DRot(i,i).gt.MaxFII) MaxFII=TMAT2DRot(i,i)
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistHii(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOrbs
                BinNo=CEILING((TMAT2DRot(i,i)-MinFII)*4002/(MaxFII-MinFII))
                ROHistHii(2,BinNo)=ROHistHii(2,BinNo)+1.0         
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFHii',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHii(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistHii(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotHii',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHii(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistHii(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
        ENDIF
   

        IF(tROHistOnePartOrbEn) THEN
            ROHistOnePartOrbEn(:,:)=0.0_dp
            MaxFII=0.0_dp
            MinFII=0.0_dp
            do i=1,NoOrbs
                OnePartOrbEnValue=0.0_dp
                OnePartOrbEnValue=OnePartOrbEnValue+TMAT2DRot(i,i)
                do j=1,NoOcc
                    OnePartOrbEnValue=OnePartOrbEnValue+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                enddo
                IF(OnePartOrbEnValue.gt.MaxFII) MaxFII=OnePartOrbEnValue
                IF(OnePartOrbEnValue.lt.MinFII) MinFII=OnePartOrbEnValue
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistOnePartOrbEn(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOrbs
                OnePartOrbEnValue=0.0_dp
                OnePartOrbEnValue=OnePartOrbEnValue+TMAT2DRot(i,i)
                do j=1,NoOcc
                    OnePartOrbEnValue=OnePartOrbEnValue+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                enddo
                BinNo=CEILING((OnePartOrbEnValue-MinFII)*4002/(MaxFII-MinFII))
                ROHistOnePartOrbEn(2,BinNo)=ROHistOnePartOrbEn(2,BinNo)+1.0         
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFOnePartOrbEn',STATUS='unknown')
                do j=1,4002
                    IF(ROHistOnePartOrbEn(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistOnePartOrbEn(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotOnePartOrbEn',STATUS='unknown')
                do j=1,4002
                    IF(ROHistOnePartOrbEn(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistOnePartOrbEn(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
        ENDIF
  

        IF(tROHistDoubExc) THEN
            ROHistDoubExc(:,:)=0.0_dp
            MaxFII=0.0_dp
            MinFII=0.0_dp
            do l=NoOcc+1,NoOrbs
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do i=1,NoOcc
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).gt.MaxFII) THEN
                                MaxFII=(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            ENDIF
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).lt.MinFII) THEN
                                MinFII=(FourIndInts(i,j,k,l))-FourIndInts(i,j,l,k)
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistDoubExc(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do l=NoOcc+1,NoOrbs
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do i=1,NoOcc
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).ne.0) THEN
                                BinNo=CEILING(((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))-MinFII)*4002/(MaxFII-MinFII))
                                ROHistDoubExc(2,BinNo)=ROHistDoubExc(2,BinNo)+1.0         
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHFDoubExc',STATUS='unknown')
                do j=1,4002
                    IF(ROHistDoubExc(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistDoubExc(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
            IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRotDoubExc',STATUS='unknown')
                do j=1,4002
                    IF(ROHistDoubExc(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistDoubExc(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF
        ENDIF


        IF(tROHistER) THEN
            ROHistER(:,:)=0.0_dp
            MaxFII=0.0_dp
            MinFII=0.0_dp
            do i=1,NoOrbs
                IF(FourIndInts(i,i,i,i).gt.MaxFII) MaxFII=FourIndInts(i,i,i,i)
                IF(FourIndInts(i,i,i,i).lt.MinFII) MinFII=FourIndInts(i,i,i,i)  
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.0_dp
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistER(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOrbs
                BinNo=CEILING((FourIndInts(i,i,i,i)-MinFII)*4002/(MaxFII-MinFII))
                ROHistER(2,BinNo)=ROHistER(2,BinNo)+1.0         
            enddo

            IF(Iteration.eq.0) THEN 
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistHF-ER',STATUS='unknown')
                do j=1,4002
                    IF(ROHistER(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistER(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF

            IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
                iunit = get_free_unit()
                OPEN(iunit,FILE='HistRot-ER',STATUS='unknown')
                do j=1,4002
                    IF(ROHistER(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(iunit,'(F20.10)',advance='no') ROHistER(i,j)
                        enddo
                        WRITE(iunit,*) ''
                    ENDIF
                enddo
                CLOSE(iunit)
            ENDIF

        ENDIF


    ENDSUBROUTINE WriteDoubHisttofile 




    SUBROUTINE PrintIntegrals()
        INTEGER :: i,j,k,l, io1, io2
        real(dp) :: DiagOneElPot,ERPot,ijVirtOneElPot,ijVirtCoulPot,ijVirtExchPot
        real(dp) :: singCoulijVirt,singExchijVirt,singCoulconHF,singExchconHF,ijklPot,ijklantisymPot
        real(dp) :: ijOccVirtOneElPot,ijOccVirtCoulPot,ijOccVirtExchPot

        io1 = 0
        io2 = 0
        IF(tInitIntValues) THEN
            io1 = get_free_unit()
            OPEN(io1,FILE='DiagIntegrals',STATUS='unknown')
            WRITE(io1,'(A10,6A18)') "Iteration","<i|h|i> ivirt","<ii|ii> ivirt","<ij|ij> iOccjVirt","<ij|ji> iOccjVirt", &
                "<ij|ij> ijVirt","<ij|ji> ijVirt"

            io2 = get_free_unit()
            OPEN(io2,FILE='SingExcIntegrals',STATUS='unknown')
            WRITE(io2,'(A10,6A18)') "Iteration","<i|h|j> iOccjVirt","<i|h|j> ijVirt","<ik|jk> HFcon","<ik|kj> HFcon", &
                "<ik|jk> ijVirt","<ik|kj> ijVirt"


            DiagOneElPotInit=0.0_dp
            ERPotInit=0.0_dp
            ijVirtOneElPotInit=0.0_dp
            ijVirtCoulPotInit=0.0_dp
            ijVirtExchPotInit=0.0_dp
            singCoulconHFInit=0.0_dp
            singExchconHFInit=0.0_dp
            singCoulijVirtInit=0.0_dp
            singExchijVirtInit=0.0_dp
            ijklPotInit=0.0_dp
            ijklantisymPotInit=0.0_dp
            ijOccVirtOneElPotInit=0.0_dp
            ijOccVirtCoulPotInit=0.0_dp
            ijOccVirtExchPotInit=0.0_dp
            NoInts01=0
            NoInts02=0
            NoInts03=0
            NoInts04=0
            NoInts05=0
            NoInts06=0
            do i=1,NoOrbs
                IF(i.gt.NoOcc) THEN
                    DiagOneElPotInit=DiagOneElPotInit+TMAT2DRot(i,i)
                    ERPotInit=ERPotInit+FourIndInts(i,i,i,i)
                    NoInts01=NoInts01+1
                    do j=NoOcc+1,NoOrbs
                       ! The i,j terms with i and j both virtual.
                       IF(j.gt.i) THEN
                           ijVirtOneElPotInit=ijVirtOneElPotInit+TMAT2DRot(i,j)
                           ijVirtCoulPotInit=ijVirtCoulPotInit+FourIndInts(i,j,i,j)
                           ijVirtExchPotInit=ijVirtExchPotInit+FourIndInts(i,j,j,i)
                           NoInts02=NoInts02+1
                       ENDIF
                       do k=1,NoOrbs
                           IF(k.gt.(NoOcc+1)) THEN
                               do l=NoOcc+1,NoOrbs
                                   IF(l.eq.j) CYCLE
                                   ijklPotInit=ijklPotInit+FourIndInts(i,j,k,l)
                                   ijklantisymPotInit=ijklantisymPotInit+FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)
                                   NoInts04=NoInts04+1
                               enddo
                            ELSE
                                IF(i.eq.j) CYCLE
                                singCoulijVirtInit=singCoulijVirtInit+FourIndInts(i,k,j,k)
                                singExchijVirtInit=singExchijVirtInit+FourIndInts(i,k,k,j)
                                NoInts03=NoInts03+1
                            ENDIF
                       enddo
                   enddo
               ELSE
                   do j=NoOcc+1,NoOrbs
                       do k=1,NoOcc
                           singCoulconHFInit=singCoulconHFInit+FourIndInts(i,k,j,k)
                           singExchconHFInit=singExchconHFInit+FourIndInts(i,k,k,j)
                           NoInts06=NoInts06+1
                       enddo
                       ijOccVirtOneElPotInit=ijOccVirtOneElPotInit+TMAT2DRot(i,j)
                       ijOccVirtCoulPotInit=ijOccVirtCoulPotInit+FourIndInts(i,j,i,j)
                       ijOccVirtExchPotInit=ijOccVirtExchPotInit+FourIndInts(i,j,j,i)
                       NoInts05=NoInts05+1
                   enddo
               ENDIF
            enddo
        ENDIF
        

        DiagOneElPot=0.0_dp
        ERPot=0.0_dp
        ijVirtOneElPot=0.0_dp
        ijVirtCoulPot=0.0_dp
        ijVirtExchPot=0.0_dp
        singCoulconHF=0.0_dp
        singExchconHF=0.0_dp
        singCoulijVirt=0.0_dp
        singExchijVirt=0.0_dp
        ijklPot=0.0_dp
        ijklantisymPot=0.0_dp
        ijOccVirtOneElPot=0.0_dp
        ijOccVirtCoulPot=0.0_dp
        ijOccVirtExchPot=0.0_dp
        do i=1,NoOrbs
            IF(i.gt.NoOcc) THEN
                DiagOneElPot=DiagOneElPot+TMAT2DRot(i,i)
                ERPot=ERPot+FourIndInts(i,i,i,i)
                do j=NoOcc+1,NoOrbs
                   ! The i,j terms with i and j both virtual.
                   IF(j.gt.i) THEN
                       ijVirtOneElPot=ijVirtOneElPot+TMAT2DRot(i,j)
                       ijVirtCoulPot=ijVirtCoulPot+FourIndInts(i,j,i,j)
                       ijVirtExchPot=ijVirtExchPot+FourIndInts(i,j,j,i)
                   ENDIF
                   do k=1,NoOrbs
                       IF(k.gt.(NoOcc+1)) THEN
                           do l=NoOcc+1,NoOrbs
                               IF(l.eq.j) CYCLE
                               ijklPot=ijklPot+FourIndInts(i,j,k,l)
                               ijklantisymPot=ijklantisymPot+FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)
                           enddo
                       ELSE
                           IF(i.eq.j) CYCLE
                           singCoulijVirt=singCoulijVirt+FourIndInts(i,k,j,k)
                           singExchijVirt=singExchijVirt+FourIndInts(i,k,k,j)
                       ENDIF
                   enddo
               enddo
           ELSE
               do j=NoOcc+1,NoOrbs
                   do k=1,NoOcc
                       singCoulconHF=singCoulconHF+FourIndInts(i,k,j,k)
                       singExchconHF=singExchconHF+FourIndInts(i,k,k,j)
                   enddo
                   ijOccVirtOneElPot=ijOccVirtOneElPot+TMAT2DRot(i,j)
                   ijOccVirtCoulPot=ijOccVirtCoulPot+FourIndInts(i,j,i,j)
                   ijOccVirtExchPot=ijOccVirtExchPot+FourIndInts(i,j,j,i)
               enddo
           ENDIF
        enddo

!        WRITE(6,*) 'ijOccVirtExchPot',ijOccVirtExchPot,ijOccVirtExchPotInit

        DiagOneElPot=(DiagOneElPot-DiagOneElPotInit)/NoInts01
        ERPot=(ERPot-ERPotInit)/NoInts01
        ijVirtOneElPot=(ijVirtOneElPot-ijVirtOneElPotInit)/NoInts02
        ijVirtCoulPot=(ijVirtCoulPot-ijVirtCoulPotInit)/NoInts02
        ijVirtExchPot=(ijVirtExchPot-ijVirtExchPotInit)/NoInts02
        singCoulijVirt=(singCoulijVirt-singCoulijVirtInit)/NoInts03
        singExchijVirt=(singExchijVirt-singExchijVirtInit)/NoInts03
        singCoulconHF=(singCoulconHF-singCoulconHFInit)/NoInts06
        singExchconHF=(singExchconHF-singExchconHFInit)/NoInts06
        ijklPot=(ijklPot-ijklPotInit)/NoInts04
        ijklantisymPot=(ijklantisymPot-ijklantisymPotInit)/NoInts04
        ijOccVirtOneElPot=(ijOccVirtOneElPot-ijOccVirtOneElPotInit)/NoInts05
        ijOccVirtCoulPot=(ijOccVirtCoulPot-ijOccVirtCoulPotInit)/NoInts05
        ijOccVirtExchPot=(ijOccVirtExchPot-ijOccVirtExchPot)/NoInts05


        WRITE(io1,'(I10,6F18.10)') Iteration,DiagOneElPot,ERPot,ijOccVirtCoulPot,ijOccVirtExchPot,ijVirtCoulPot,ijVirtExchPot
        WRITE(io2,'(I10,6F18.10)') Iteration,ijOccVirtOneElPot,ijVirtOneElPot,singCoulconHF,singExchconHF,singCoulijVirt, &
            singExchijVirt

        IF((.not.tNotConverged).and.(.not.tInitIntValues)) THEN
            CLOSE(io1)
            CLOSE(io2)
        ENDIF


    ENDSUBROUTINE PrintIntegrals


    SUBROUTINE CalcFOCKMatrix()
        USE SystemData , only : nBasis
        USE Logging , only : tRDMonfly
        INTEGER :: i,j,k,l,a,b,ierr
        real(dp) :: FOCKDiagSumHF,FOCKDiagSumNew
        CHARACTER(len=*) , PARAMETER :: this_routine='CalcFOCKMatrix'
        !NEED TO FIX THIS!

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

! When transforming the orbitals into approximate natural orbitals, we want to save memory, so don't bother
! calculating the whole matrix, just the diagonal elements that we actually need.

    
        IF(tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) THEN
            ALLOCATE(ArrDiagNew(NoOrbs),stat=ierr)
            CALL LogMemAlloc('ArrDiagNew',NoOrbs,8,this_routine,ArrDiagNewTag,ierr)
            ArrDiagNew(:)=0.0_dp                     
        ELSE
            ALLOCATE(ArrNew(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ArrNew',NoOrbs**2,8,this_routine,ArrNewTag,ierr)
            ArrNew(:,:)=0.0_dp                     
        ENDIF

!        WRITE(6,*) 'The diagonal fock elements in the HF basis set'
!        do a=1,nBasis
!            WRITE(6,'(F20.10)',advance='no') Arr(a,2)
!        enddo


! First calculate the sum of the diagonal elements, ARR.
! Check if this is already being done.
        FOCKDiagSumHF=0.0_dp
        do a=1,nBasis        
            FOCKDiagSumHF=FOCKDiagSumHF+Arr(a,2)
        enddo

        WRITE(6,*) 'Sum of the fock matrix diagonal elements in the HF basis set = ',FOCKDiagSumHF

!        WRITE(6,*) 'Coeffs'
!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                WRITE(6,'(F20.10)',advance='no') CoeffT1(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

! Then calculate the fock matrix in the transformed basis, and the sum of the new diagonal elements.
! Our Arr in spin orbitals.
!        do j=1,NoOrbs
!            ArrNew(j,j)=Arr(2*j,2)
!        enddo

        FOCKDiagSumNew=0.0_dp
        do j=1,NoRotOrbs
            l=SymLabelList3_rot(j)
            IF(tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) THEN
                do a=1,NoOrbs
                    b=SymLabelList2_rot(a)
                    IF(tStoreSpinOrbs.or.tTurnStoreSpinOff) THEN
                        ArrDiagNew(l)=ArrDiagNew(l)+(CoeffT1(a,j)*ARR(b,2)*CoeffT1(a,j))
                    ELSE
                        ArrDiagNew(l)=ArrDiagNew(l)+(CoeffT1(a,j)*ARR(2*b,2)*CoeffT1(a,j))
                    ENDIF
                enddo
                IF(tStoreSpinOrbs.or.tTurnStoreSpinOff) THEN
                    FOCKDiagSumNew=FOCKDiagSumNew+(ArrDiagNew(l))
                ELSE
                    FOCKDiagSumNew=FOCKDiagSumNew+(ArrDiagNew(l)*2)
                ENDIF
            ELSE
                do i=1,NoRotOrbs
                    k=SymLabelList2_rot(i)
                    ArrNew(k,l)=0.0_dp
                    do a=1,NoOrbs
                        b=SymLabelList2_rot(a)
                        IF(tStoreSpinOrbs.or.tTurnStoreSpinOff) THEN
                            ArrNew(k,l)=ArrNew(k,l)+(CoeffT1(a,i)*Arr(b,2)*CoeffT1(a,j))
                        ELSE
                            ArrNew(k,l)=ArrNew(k,l)+(CoeffT1(a,i)*Arr(2*b,2)*CoeffT1(a,j))
                        ENDIF
                    enddo
                enddo
                IF(tStoreSpinOrbs.or.tTurnStoreSpinOff) THEN
                    FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(l,l))
                ELSE
                    FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(l,l)*2)
                ENDIF
                !only running through spat orbitals, count each twice to compare to above.
            ENDIF
        enddo
        ! If we are truncation the virtual space, only the unfrozen entries will be transformed.
        

        WRITE(6,*) 'Sum of the fock matrix diagonal elements in the transformed basis set = ',FOCKDiagSumNew

!        WRITE(6,*) 'The fock matrix for the transformed orbitals'
!        do j=1,NoOrbs
!            do i=1,NoOrbs
!                WRITE(6,'(F20.10)',advance='no') ArrNew(i,j)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'BRR then ARR before being changed',nBasis
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2),ArrDiagNew(i)
!        enddo
       
!        WRITE(6,*) 'to here',NoDumpTruncs,NoOrbs 
!        CALL neci_flush(6)

! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2) (ordered in terms of orbital number).
! ARR(:,2) needs to be ordered in terms of symmetry and then energy (like SymLabelList), so currently this ordering will not be 
! correct when reading in qchem INTDUMPS as the orbital number ordering is by energy.

        IF(NoDumpTruncs.le.1) THEN
! If we are only writing out 1 ROFCIDUMP or we are not truncating at all - can refill ARR etc.            

            IF(tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly) THEN
                IF(tStoreSpinOrbs.or.tTurnStoreSpinOff) THEN
                    do j=1,NoOrbs
                        ARR(j,2)=ArrDiagNew(j)
                        ARR(j,1)=ArrDiagNew(BRR(j))
                    enddo
                ELSE
                    do j=1,NoOrbs
                        ARR(2*j,2)=ArrDiagNew(j)
                        ARR(2*j-1,2)=ArrDiagNew(j)
                        ARR(2*j,1)=ArrDiagNew(BRR(2*j)/2)
                        ARR(2*j-1,1)=ArrDiagNew(BRR(2*j)/2)
                    enddo
                ENDIF
            ELSE
                IF(tStoreSpinOrbs.or.tTurnStoreSpinOff) THEN
                    do j=1,NoRotOrbs
                        ARR(j,2)=ArrNew(j,j)
                        ARR(j,1)=ArrNew(BRR(j),BRR(j))
                    enddo
                ELSE
                    do j=1,NoRotOrbs
                        ARR(2*j,2)=ArrNew(j,j)
                        ARR(2*j-1,2)=ArrNew(j,j)
                        ARR(2*j,1)=ArrNew(BRR(2*j)/2,BRR(2*j)/2)
                        ARR(2*j-1,1)=ArrNew(BRR(2*j)/2,BRR(2*j)/2)
                    enddo
                ENDIF
            ENDIF

        ENDIF

!        WRITE(6,*) 'BRR then ARR after being changed'
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
!        CALL neci_flush(6)
!        stop       

        IF((tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs.or.tRDMonfly).and.(NoDumpTruncs.le.1)) THEN
            DEALLOCATE(ArrDiagNew)
            CALL LogMemDealloc(this_routine,ArrDiagNewTag)
        ELSEIF(NoDumpTruncs.le.1) THEN
            DEALLOCATE(ArrNew)
            CALL LogMemDealloc(this_routine,ArrNewTag)
        ENDIF


        WRITE(6,*) 'end of calcfockmatrix'
        call neci_flush(6)

    ENDSUBROUTINE CalcFOCKMatrix


    SUBROUTINE RefillUMATandTMAT2D()
        INTEGER :: l,k,j,i,a,b,g,d,c,nBasis2,ierr
        INTEGER(TagIntType) :: TMAT2DPartTag
        real(dp) :: NewTMAT
        real(dp) , ALLOCATABLE :: TMAT2DPart(:,:)
#ifdef __CMPLX
        call stop_all('RefillUMATandTMAT2D', 'Rotating orbitals not implemented for complex orbitals.')
#endif

        IF(tStoreSpinOrbs) THEN
            ALLOCATE(TMAT2DPart((nBasis-NoFrozenVirt),nBasis),stat=ierr)
            CALL LogMemAlloc('TMAT2DPart',(nBasis-NoFrozenVirt)*nBasis,8,'RefillUMAT',TMAT2DPartTag,ierr)
            IF(NoDumpTruncs.gt.1) THEN
                ALLOCATE(TMAT2DNew((nBasis-NoFrozenVirt),(nBasis-NoFrozenVirt)),stat=ierr)
                CALL LogMemAlloc('TMAT2DNew',(nBasis-NoFrozenVirt)**2,8,'RefillUMAT',TMAT2DNewTag,ierr)
                TMAT2DNew(:,:)=0.0_dp
            ENDIF
        ELSE
            ALLOCATE(TMAT2DPart((nBasis-(NoFrozenVirt*2)),nBasis),stat=ierr)
            CALL LogMemAlloc('TMAT2DPart',(nBasis-(NoFrozenVirt*2))*nBasis,8,'RefillUMAT',TMAT2DPartTag,ierr)
            IF(NoDumpTruncs.gt.1) THEN
                ALLOCATE(TMAT2DNew((nBasis-NoFrozenVirt),(nBasis-NoFrozenVirt)),stat=ierr)
                CALL LogMemAlloc('TMAT2DNew',(nBasis-NoFrozenVirt)**2,8,'RefillUMAT',TMAT2DNewTag,ierr)
                TMAT2DNew(:,:)=0.0_dp
            ENDIF
        ENDIF
        TMAT2DPart(:,:)=0.0_dp

        RefillUMAT_Time%timer_name='RefillUMATandTMAT'
        CALL set_timer(RefillUMAT_Time,30)

        do i = 1,nBasis
            WRITE(6,*) SymLabelList2_rot(i),SymLabelList3_rot(i)
        enddo


! Make the UMAT elements the four index integrals.  These are calculated by transforming the HF orbitals using the
! coefficients that have been found
        IF(NoDumpTruncs.le.1) THEN
            do l=1,(NoOrbs-(NoFrozenVirt))
                IF(tTurnStoreSpinOff) THEN
                    d=CEILING(REAL(SymLabelList3_rot(l))/2.0_dp)
                ELSE
                    d=SymLabelList3_rot(l)
                ENDIF
                do k=1,(NoOrbs-(NoFrozenVirt))

                    IF(tTurnStoreSpinOff) THEN
                        g=CEILING(REAL(SymLabelList3_rot(k))/2.0_dp)
                    ELSE
                        g=SymLabelList3_rot(k)
                    ENDIF
     
                    do j=1,(NoOrbs-(NoFrozenVirt))

                        IF(tTurnStoreSpinOff) THEN
                            b=CEILING(REAL(SymLabelList3_rot(j))/2.0_dp)
                        ELSE
                            b=SymLabelList3_rot(j)
                        ENDIF
                        do i=1,(NoOrbs-(NoFrozenVirt))

                            IF(tTurnStoreSpinOff) THEN
                                a=CEILING(REAL(SymLabelList3_rot(i))/2.0_dp)
                            ELSE
                                a=SymLabelList3_rot(i)
                            ENDIF

                            IF(tUseMP2VarDenMat.or.tFindCINatOrbs.or.tReadInCoeff) THEN
                                UMAT(UMatInd(a,b,g,d,0,0))=(FourIndInts(i,k,j,l))
                            ELSE
                                UMAT(UMatInd(a,b,g,d,0,0))=(FourIndInts(i,j,k,l))
                            ENDIF
                        enddo
                    enddo
                enddo
            enddo
        ENDIF

! Also calculate the 2 index integrals, and make these the elements of the TMAT2D matrix.
! TMAT2D is in spin orbitals.

!        WRITE(6,*) 'TMAT2D before transformation' 
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l),8)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'SpatOrbs',SpatOrbs
!        WRITE(6,*) 'NoRotOrbs',NoRotOrbs
!        CALL neci_flush(6)

!        do l=1,NoRotOrbs
!            j=SymLabelList3_rot(l)
!            do k=1,NoRotOrbs
!                i=SymLabelList3_rot(k)
!                NewTMAT02=0.0_dp
!
!                do b=1,NoOrbs
!                    d=SymLabelList2_rot(b)
!                    NewTMAT=0.0_dp
!                    ! for a particular beta, find the transformed integral <i|h|b>
!                    do a=1,NoOrbs
!                        c=SymLabelList2_rot(a)
!                        IF(tStoreSpinOrbs) THEN
!                            NewTMAT=NewTMAT+(CoeffT1(a,k)*REAL(TMAT2D(c,d),8))
!                        ELSE
!                            NewTMAT=NewTMAT+(CoeffT1(a,k)*REAL(TMAT2D(2*c,2*d),8))
!                        ENDIF
!                    enddo
                    ! NewTMAT is then <i|h|b> for a particular i and b.

                    ! then transform the beta part as well.
!                    NewTMAT02=NewTMAT02+(CoeffT1(b,l)*NewTMAT)
                    ! NewTMAT02 become <i|h|j> for a particular i and j.
!                enddo
!                IF(tStoreSpinOrbs) THEN
!                    TMAT2D(i,j)=(NewTMAT02)
!                ELSE
!                    TMAT2D(2*i,2*j)=(NewTMAT02)
!                    TMAT2D(2*i-1,2*j-1)=(NewTMAT02)
!                ENDIF
!            enddo
!        enddo



        do a=1,nBasis
            do k=1,NoRotOrbs
                i=SymLabelList3_rot(k)
                NewTMAT=0.0_dp
                do b=1,NoOrbs
                    d=SymLabelList2_rot(b)
                    IF(tStoreSpinOrbs) THEN
                        NewTMAT=NewTMAT+(CoeffT1(b,k)*REAL(TMAT2D(d,a),dp))
                    ELSE
                        NewTMAT=NewTMAT+(CoeffT1(b,k)*REAL(TMAT2D(2*d,a),dp))
                    ENDIF
                enddo
                IF(tStoreSpinOrbs) THEN
                    TMAT2DPart(i,a)=NewTMAT
                ELSE
                    TMAT2DPart(2*i,a)=NewTMAT
                    TMAT2DPart(2*i-1,a)=NewTMAT
                ENDIF
            enddo
        enddo



        IF(tStoreSpinOrbs) THEN
            nBasis2=nBasis-NoFrozenVirt
        ELSE
            nBasis2=nBasis-(NoFrozenVirt*2)
        ENDIF
        do k=1,nBasis2
            do l=1,NoRotOrbs
                j=SymLabelList3_rot(l)
                NewTMAT=0.0_dp
                do a=1,NoOrbs
                    c=SymLabelList2_rot(a)
                    IF(tStoreSpinOrbs) THEN
                        NewTMAT=NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,c))
                    ELSE
                        NewTMAT=NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,2*c))
                    ENDIF
                enddo
                IF(tStoreSpinOrbs) THEN
                    IF(NoDumpTruncs.gt.1) THEN
                        TMAT2DNew(k,j)=NewTMAT
                    ELSE
                        TMAT2D(k,j)=(NewTMAT)
                    ENDIF
                ELSE
                    IF(NoDumpTruncs.gt.1) THEN
                        TMAT2DNew(k,2*j)=NewTMAT
                        TMAT2DNew(k,2*j-1)=NewTMAT
                    ELSE
                        TMAT2D(k,2*j)=(NewTMAT)
                        TMAT2D(k,2*j-1)=(NewTMAT)
                    ENDIF
                ENDIF
            enddo
        enddo

    
!        WRITE(6,*) 'TMAT2D after transformation'
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l),8)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        CALL neci_flush(6)
!        stop

        DEALLOCATE(TMAT2DPart)
        CALL LogMemDeAlloc('RefillUMAT',TMAT2DPartTag)


        IF(tROHistSingExc) CALL WriteSingHisttofile()

        CALL set_timer(RefillUMAT_Time,30)

        IF(tTurnStoreSpinOff) THEN
            tStoreSpinOrbs = .false.
            NoOrbs = nBasis / 2
        ENDIF

        WRITE(6,'(A,I5,A)') ' Printing the new ROFCIDUMP file for a truncation of ',NoFrozenVirt,' orbitals.'
        IF(tROFciDump.and.(NoDumpTruncs.gt.1)) THEN
            CALL PrintRepeatROFCIDUMP()
        ELSEIF(tROFciDUmp) THEN
            CALL PrintROFCIDUMP()
        ENDIF



    ENDSUBROUTINE RefillUMATandTMAT2D


    

    SUBROUTINE PrintROFCIDUMP()
!This prints out a new FCIDUMP file in the same format as the old one.
        INTEGER :: i,j,k,l,iunit
        CHARACTER(len=5) :: Label
        CHARACTER(len=20) :: LabelFull


        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
        CALL set_timer(PrintROFCIDUMP_Time,30)

        Label=''
        LabelFull=''
        WRITE(Label,'(I5)') NoFrozenVirt
        LabelFull='ROFCIDUMP-'//adjustl(Label)

        iunit = get_free_unit()
        OPEN(iunit,FILE=LabelFull,STATUS='unknown')
        
        WRITE(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',(NoOrbs-(NoFrozenVirt)),',NELEC=',NEl,',MS2=',LMS,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,(NoOrbs-(NoFrozenVirt))
            IF((tUseMP2VarDenMat.or.tFindCINatOrbs).and.(.not.lNoSymmetry).and.tTruncRODump) THEN
                WRITE(iunit,'(I1,A1)',advance='no') (SymOrbs_rot(i)+1),','
            ELSE
                IF(tStoreSpinOrbs) THEN
                    WRITE(iunit,'(I1,A1)',advance='no') (INT(G1(i)%sym%S)+1),','
                ELSE
                    WRITE(iunit,'(I1,A1)',advance='no') (INT(G1(i*2)%sym%S)+1),','
                ENDIF
            ENDIF
        enddo
        WRITE(iunit,*) ''
        IF(tStoreSpinOrbs) THEN
            WRITE(iunit,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        ELSE
            WRITE(iunit,'(A7,I1)') 'ISYM=',1
        ENDIF
        WRITE(iunit,'(A5)') '&END'
       
        do i=1,(NoOrbs-(NoFrozenVirt))
            do k=1,i
                do j=1,(NoOrbs-(NoFrozenVirt))
!                    Sym=IEOR(INT(G1(j*2)%sym%S),IEOR(INT(G1(k*2)%sym%S),INT(G1(i*2)%sym%S)))
                    ! Potential to put symmetry in here, have currently taken it out, because when we're only printing non-zero values,
                    ! it is kind of unnecessary - although it may be used to speed things up.
                    do l=1,j
!                        Syml=INT(G1(l*2)%sym%S)
!                        IF((Syml.eq.Sym).and.((REAL(UMat(UMatInd(i,j,k,l,0,0)),8)).ne.0.0_dp)) &
                        IF((ABS(REAL(UMat(UMatInd(i,j,k,l,0,0)),dp))).ne.0.0_dp) &
                                        &WRITE(iunit,'(F21.12,4I3)') REAL(UMat(UMatInd(i,j,k,l,0,0)),dp),i,k,j,l 
                    enddo
                enddo
           enddo
        enddo

       
!        do i=1,SpatOrbs
!            do j=1,SpatOrbs
!                do l=j,SpatOrbs
!                    Sym=IEOR(INT(G1(l*2)%sym%S),IEOR(INT(G1(j*2)%sym%S),INT(G1(i*2)%sym%S)))
!                    do l=SymLabelCounts2_rot(1,Sym+SymMin),(SymLabelCounts2_rot(1,Sym+SymMin)+SymLabelCounts2_rot(2,Sym+SymMin)-1)
                   ! Potential to put symmetry in here.
!                    do k=i,SpatOrbs
!                        Symk=INT(G1(k*2)%sym%S)
!                        IF(Symk.eq.Sym) WRITE(iunit,'(F21.12,4I3)') REAL(UMat(UMatInd(i,j,k,l,0,0)),8),i,k,j,l 
!                    enddo
!                enddo
!            enddo
!        enddo


! TMAT2D stored as spin orbitals
        do k=1,(NoOrbs-(NoFrozenVirt))
            ! Symmetry?
            do i=k,(NoOrbs-(NoFrozenVirt))
                IF(tStoreSpinOrbs) THEN
                    IF((REAL(TMAT2D(i,k),dp)).ne.0.0_dp) WRITE(iunit,'(F21.12,4I3)') REAL(TMAT2D(i,k),dp),i,k,0,0
                ELSE
                    IF((REAL(TMAT2D(2*i,2*k),dp)).ne.0.0_dp) WRITE(iunit,'(F21.12,4I3)') REAL(TMAT2D(2*i,2*k),dp),i,k,0,0
                ENDIF
            enddo
        enddo

! ARR has the energies of the orbitals (eigenvalues).  ARR(:,2) has ordering we want.
! ARR is stored as spin orbitals.

        do k=1,(NoOrbs-(NoFrozenVirt))
            IF(tStoreSpinOrbs) THEN
                WRITE(iunit,'(F21.12,4I3)') Arr(k,2),k,0,0,0
            ELSE
                WRITE(iunit,'(F21.12,4I3)') Arr(2*k,2),k,0,0,0
            ENDIF
        enddo

        WRITE(iunit,'(F21.12,4I3)') ECore,0,0,0,0
        
        CALL neci_flush(iunit)

        CLOSE(iunit)

        CALL halt_timer(PrintROFCIDUMP_Time)


    ENDSUBROUTINE PrintROFCIDUMP


    SUBROUTINE PrintRepeatROFCIDUMP()
!This prints out a new FCIDUMP file in the same format as the old one.
        INTEGER :: i,j,k,l,ierr,a,b,g,d, iunit
        CHARACTER(len=5) :: Label
        CHARACTER(len=20) :: LabelFull
        CHARACTER(len=*) , PARAMETER :: this_routine='PrintRepeatROFCIDUMP'


        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
        CALL set_timer(PrintROFCIDUMP_Time,30)

        Label=''
        LabelFull=''
        WRITE(Label,'(I5)') NoFrozenVirt
        LabelFull='ROFCIDUMP-'//adjustl(Label)

        iunit = get_free_unit()
        OPEN(iunit,FILE=LabelFull,STATUS='unknown')
        
        WRITE(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',(NoOrbs-(NoFrozenVirt)),',NELEC=',NEl,',MS2=',LMS,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,(NoOrbs-(NoFrozenVirt))
            IF((tUseMP2VarDenMat.or.tFindCINatOrbs).and.(.not.lNoSymmetry).and.tTruncRODump) THEN
                WRITE(iunit,'(I1,A1)',advance='no') (SymOrbs_rot(i)+1),','
            ELSE
                IF(tStoreSpinOrbs) THEN
                    WRITE(iunit,'(I1,A1)',advance='no') (INT(G1(i)%sym%S)+1),','
                ELSE
                    WRITE(iunit,'(I1,A1)',advance='no') (INT(G1(i*2)%sym%S)+1),','
                ENDIF
            ENDIF
        enddo
        WRITE(iunit,*) ''
        IF(tStoreSpinOrbs) THEN
            WRITE(iunit,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        ELSE
            WRITE(iunit,'(A7,I1)') 'ISYM=',1
        ENDIF
        WRITE(iunit,'(A5)') '&END'

 
        ALLOCATE(SymLabelList3_rotInv(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3_rotInv',NoOrbs,4,this_routine,SymLabelList3_rotInvTag,ierr)
        SymLabelList3_rotInv(:)=0                     

        do i=1,NoOrbs
            SymLabelList3_rotInv(SymLabelList3_rot(i))=i
        enddo
       
        do i=1,(NoOrbs-(NoFrozenVirt))
            a=SymLabelList3_rotInv(i)
            do k=1,i
                g=SymLabelList3_rotInv(k)
                do j=1,(NoOrbs-(NoFrozenVirt))
                    b=SymLabelList3_rotInv(j)
!                    Sym=IEOR(INT(G1(j*2)%sym%S),IEOR(INT(G1(k*2)%sym%S),INT(G1(i*2)%sym%S)))
                    ! Potential to put symmetry in here, have currently taken it out, because when we're only printing non-zero values,
                    ! it is kind of unnecessary - although it may be used to speed things up.
                    do l=1,j
                        d=SymLabelList3_rotInv(l)
!                        Syml=INT(G1(l*2)%sym%S)
!                        IF((Syml.eq.Sym).and.((REAL(UMat(UMatInd(i,j,k,l,0,0)),8)).ne.0.0_dp)) &
!                        WRITE(6,*) i,a,k,g,j,b,l,d,FourIndInts(a,b,g,d)
                        IF((ABS(FourIndInts(a,g,b,d))).ne.0.0_dp) &
                                        &WRITE(iunit,'(F21.12,4I3)') FourIndInts(a,g,b,d),i,k,j,l 
 
                    enddo
                enddo
           enddo
        enddo

        DEALLOCATE(SymLabelList3_rotInv)
        CALL LogMemDeAlloc(this_routine,SymLabelList3_rotInvTag)


! TMAT2D stored as spin orbitals
        do k=1,(NoOrbs-(NoFrozenVirt))
            ! Symmetry?
            do i=k,(NoOrbs-(NoFrozenVirt))
                IF(tStoreSpinOrbs) THEN
                    IF(TMAT2DNew(i,k).ne.0.0_dp) WRITE(iunit,'(F21.12,4I3)') TMAT2DNew(i,k),i,k,0,0
                ELSE
                    IF(TMAT2DNew(2*i,2*k).ne.0.0_dp) WRITE(iunit,'(F21.12,4I3)') TMAT2DNew(2*i,2*k),i,k,0,0
                ENDIF
            enddo
        enddo

! ARR has the energies of the orbitals (eigenvalues).  ARR(:,2) has ordering we want.
! ARR is stored as spin orbitals.

        IF(tUseMP2VarDenMat.or.tFindCINatOrbs.or.tUseHFOrbs) THEN
            IF(tStoreSpinOrbs) THEN
                do k=1,(NoOrbs-(NoFrozenVirt))
                    WRITE(iunit,'(F21.12,4I3)') ArrDiagNew(k),k,0,0,0
                enddo
            ELSE

                do k=1,(NoOrbs-(NoFrozenVirt))

                    WRITE(iunit,'(F21.12,4I3)') ArrDiagNew(k),k,0,0,0
                enddo
            ENDIF
        ELSE
            IF(tStoreSpinOrbs) THEN
                do k=1,(NoOrbs-(NoFrozenVirt))
                    WRITE(iunit,'(F21.12,4I3)') ArrNew(k,k),k,0,0,0
                enddo
            ELSE
                do k=1,(NoOrbs-(NoFrozenVirt))
                    WRITE(iunit,'(F21.12,4I3)') ArrNew(k,k),k,0,0,0
                enddo
            ENDIF
        ENDIF


        WRITE(iunit,'(F21.12,4I3)') ECore,0,0,0,0
        
        CALL neci_flush(iunit)

        CLOSE(iunit)

        CALL halt_timer(PrintROFCIDUMP_Time)


    ENDSUBROUTINE PrintRepeatROFCIDUMP



    SUBROUTINE DeallocateMem()
        CHARACTER(len=*) , PARAMETER :: this_routine='DeallocateMem'


        DEALLOCATE(Lab)
        CALL LogMemDealloc(this_routine,LabTag)
!        DEALLOCATE(CoeffT1Temp)
!        CALL LogMemDealloc(this_routine,CoeffT1TempTag)
        DEALLOCATE(CoeffT1)
        CALL LogMemDealloc(this_routine,CoeffT1Tag)
        DEALLOCATE(CoeffCorT2)
        CALL LogMemDealloc(this_routine,CoeffCorT2Tag)
        DEALLOCATE(CoeffUncorT2)
        CALL LogMemDealloc(this_routine,CoeffUncorT2Tag)
        IF(tLagrange) THEN
            DEALLOCATE(Lambdas)
            CALL LogMemDealloc(this_routine,LambdasTag)
            DEALLOCATE(DerivLambda)
            CALL LogMemDealloc(this_routine,DerivLambdaTag)
        ENDIF 
        DEALLOCATE(DerivCoeff)
        CALL LogMemDealloc(this_routine,DerivCoeffTag)

        DEALLOCATE(DiagTMAT2Dfull)
        CALL LogMemDealloc(this_routine,DiagTMAT2DfullTag)
!        DEALLOCATE(ThreeIndInts02Temp)
!        CALL LogMemDealloc(this_routine,ThreeIndInts02TempTag)
!        DEALLOCATE(FourIndInts02Temp)
!        CALL LogMemDealloc(this_routine,FourIndInts02TempTag)
 
        DEALLOCATE(TwoIndInts01)
        CALL LogMemDealloc(this_routine,TwoIndInts01Tag)
        DEALLOCATE(ThreeIndInts02)
        CALL LogMemDealloc(this_routine,ThreeIndInts02Tag)
        DEALLOCATE(FourIndInts)
        CALL LogMemDealloc(this_routine,FourIndIntsTag)
        DEALLOCATE(FourIndInts02)
        CALL LogMemDealloc(this_routine,FourIndInts02Tag)

        IF(tERLocalization.and.(.not.tStoreSpinOrbs)) THEN
            DEALLOCATE(TwoIndIntsER)
            CALL LogMemDeAlloc(this_routine,TwoIndIntsERTag)
            DEALLOCATE(ThreeIndInts01ER)
            CALL LogMemDeAlloc(this_routine,ThreeIndInts01ERTag)
            DEALLOCATE(ThreeIndInts02ER)
            CALL LogMemDeAlloc(this_routine,ThreeIndInts02ERTag)
            DEALLOCATE(FourIndIntsER)
            CALL LogMemDeAlloc(this_routine,FourIndIntsERTag)

        ELSE
            DEALLOCATE(TMAT2DTemp)
            CALL LogMemDealloc(this_routine,TMAT2DTempTag)
            DEALLOCATE(TMAT2DPartRot01)
            CALL LogMemDealloc(this_routine,TMAT2DPartRot01Tag)
            DEALLOCATE(TMAT2DPartRot02)
            CALL LogMemDealloc(this_routine,TMAT2DPartRot02Tag)
            DEALLOCATE(TMAT2DRot)
            CALL LogMemDealloc(this_routine,TMAT2DRotTag)
     
!            DEALLOCATE(TwoIndInts02Temp)
!            CALL LogMemDealloc(this_routine,TwoIndInts02TempTag)
!            DEALLOCATE(ThreeIndInts01Temp)
!            CALL LogMemDealloc(this_routine,ThreeIndInts01TempTag)
!            DEALLOCATE(ThreeIndInts03Temp)
!            CALL LogMemDealloc(this_routine,ThreeIndInts03TempTag)
!            DEALLOCATE(ThreeIndInts04Temp)
!            CALL LogMemDealloc(this_routine,ThreeIndInts04TempTag)
            DEALLOCATE(TwoIndInts02)
            CALL LogMemDealloc(this_routine,TwoIndInts02Tag)
            DEALLOCATE(ThreeIndInts01)
            CALL LogMemDealloc(this_routine,ThreeIndInts01Tag)
            DEALLOCATE(ThreeIndInts03)
            CALL LogMemDealloc(this_routine,ThreeIndInts03Tag)
            DEALLOCATE(ThreeIndInts04)
            CALL LogMemDealloc(this_routine,ThreeIndInts04Tag)
            DEALLOCATE(UMATTemp02)
            CALL LogMemDealloc(this_routine,UMATTemp02Tag)
        ENDIF 

        DEALLOCATE(UMATTemp01)
        CALL LogMemDealloc(this_routine,UMATTemp01Tag)
        DEALLOCATE(SymLabelList2_rot)
        CALL LogMemDealloc(this_routine,SymLabelList2_rotTag)
        DEALLOCATE(SymLabelCounts2_rot)
        CALL LogMemDealloc(this_routine,SymLabelCounts2_rotTag)
        DEALLOCATE(SymLabelListInv_rot)
        CALL LogMemDealloc(this_routine,SymLabelListInv_rotTag)

        IF(tShake) THEN
            DEALLOCATE(ShakeLambda)
            CALL LogMemDealloc(this_routine,ShakeLambdaTag)
            DEALLOCATE(ShakeLambdaNew)
            CALL LogMemDealloc(this_routine,ShakeLambdaNewTag)
            DEALLOCATE(Constraint)
            CALL LogMemDealloc(this_routine,ConstraintTag)
            DEALLOCATE(ConstraintCor)
            CALL LogMemDealloc(this_routine,ConstraintCorTag)
            DEALLOCATE(DerivConstrT1)
            CALL LogMemDealloc(this_routine,DerivConstrT1Tag)
            DEALLOCATE(DerivConstrT2)
            CALL LogMemDealloc(this_routine,DerivConstrT2Tag)
            DEALLOCATE(ForceCorrect)
            CALL LogMemDealloc(this_routine,ForceCorrectTag)
            DEALLOCATE(Correction)
            CALL LogMemDealloc(this_routine,CorrectionTag)
            IF(tShakeApprox) THEN
                DEALLOCATE(DerivConstrT1T2Diag)
                CALL LogMemDealloc(this_routine,DerivConstrT1T2DiagTag)
            ELSE
                DEALLOCATE(DerivConstrT1T2)
                CALL LogMemDealloc(this_routine,DerivConstrT1T2Tag)
            ENDIF
        ENDIF
       

    END SUBROUTINE DeallocateMem
 

    SUBROUTINE WriteShakeOUTstats01()
! Debbugging option
        INTEGER :: i,j,l,m,a

        WRITE(6,*) 'Original coefficients'
        do m=1,NoOrbs
            do a=1,NoOrbs
                WRITE(6,'(4F20.10)',advance='no') CoeffT1(a,m)
            enddo
            WRITE(6,*) ''
        enddo
        
        WRITE(6,*) 'Uncorrected Force'
        do m=1,NoOrbs
            do a=1,NoOrbs
                WRITE(6,'(4F20.10)',advance='no') DerivCoeff(a,m)
            enddo
            WRITE(6,*) ''
        enddo


        WRITE(6,*) 'Coefficients at t2, having been shifted by the uncorrected force'
        do m=1,NoOrbs
            do a=1,NoOrbs
                WRITE(6,'(4F20.10)',advance='no') CoeffUncorT2(a,m)
            enddo
            WRITE(6,*) ''
        enddo

        WRITE(6,*) 'The value of each constraint with these new unconstrained coefficients'
        do l=1,TotNoConstraints
            i=Lab(1,l)
            j=Lab(2,l)
            WRITE(6,'(F20.10)',advance='no') Constraint(l)
            WRITE(6,*) l,i,j 
        enddo
        CALL neci_flush(6) 

    ENDSUBROUTINE WriteShakeOUTstats01



    SUBROUTINE WriteShakeOUTstats02(ShakeIteration,TotLambdas,ConvergeCount)
! Debugging option    
        INTEGER :: m,a,l,ConvergeCount,ShakeIteration
        real(dp) :: TotLambdas 
    
            WRITE(6,*) 'Iteration number ,', ShakeIteration

            WRITE(6,*) 'Lambdas used for this iteration'
            do l=1,TotNoConstraints
                WRITE(6,'(10F20.10)') ShakeLambda(l)
                TotLambdas=TotLambdas+ShakeLambda(l)
            enddo

            WRITE(6,*) 'Corrected Force'
            do m=1,NoOrbs
                do a=1,NoOrbs
                    WRITE(6,'(4F20.10)',advance='no') ForceCorrect(a,m)
                enddo
                WRITE(6,*) ''
            enddo
    
            WRITE(6,*) 'Coefficients having been shifted by the corrected force (at t2)'
            do m=1,NoOrbs
                do a=1,NoOrbs
                    WRITE(6,'(4F20.10)',advance='no') CoeffCorT2(a,m)
                enddo
                WRITE(6,*) ''
            enddo
            WRITE(6,*) 'total corrected force = ',TotCorrectedForce 

            WRITE(6,*) 'The value of each constraint with the recent constrained coefficients'
            do l=1,TotNoConstraints
                WRITE(6,'(10F20.10)',advance='no') ConstraintCor(l)
            enddo

            WRITE(6,*) 'The number of constraints with values less than the convergence limit'
            WRITE(6,*) ConvergeCount


!            WRITE(6,*) 'Derivative of constraints w.r.t original coefficients (t1)'
!            do m=1,NoOrbs
!                WRITE(6,*) 'm equals ', m 
!                do l=1,TotNoConstraints
!                    WRITE(6,*) 'i, j equals ',Lab(1,l),Lab(2,l)
!                    do a=1,NoOrbs
!                        WRITE(6,'(4F20.10)',advance='no') DerivConstrT1(a,m,l)
!                    enddo
!                enddo
!            enddo

!            WRITE(6,*) 'Derivative of constraints w.r.t shifted coefficients (t2)'
!            do m=1,NoOrbs
!                WRITE(6,*) 'm equals ', m  
!                do l=1,TotNoConstraints
!                    WRITE(6,*) 'i,j equals ',Lab(1,l),Lab(2,l)
!                    do a=1,NoOrbs
!                        WRITE(6,'(4F20.10)',advance='no') DerivConstrT2(a,m,l)
!                    enddo
!                enddo
!            enddo


    ENDSUBROUTINE WriteShakeOUTstats02
   

END MODULE RotateOrbsMod
