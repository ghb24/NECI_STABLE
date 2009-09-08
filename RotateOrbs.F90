MODULE RotateOrbsMod

    USE Global_utilities
    USE Parallel 
    USE IntegralsData , only : UMAT,nFrozen,ChemPot
    USE UMatCache , only : UMatInd
    USE HElem , only : HElement
    USE SystemData , only : ConvergedForce,TimeStep,tLagrange,tShake,tShakeApprox,ShakeConverged,tROIteration,ROIterMax,tShakeIter,ShakeIterMax,OrbEnMaxAlpha
    USE SystemData , only : G1,ARR,NEl,nBasis,LMS,ECore,tSeparateOccVirt,Brr,nBasisMax,OrbOrder,lNoSymmetry,tRotatedOrbs,tERLocalization,tRotateOccOnly
    USE SystemData, only : tOffDiagMin,DiagWeight,OffDiagWeight,tRotateVirtOnly,tOffDiagSqrdMax,tOffDiagSqrdMin,tOffDiagMax,tDoubExcMin,tOneElIntMax,tOnePartOrbEnMax
    USE SystemData, only : tShakeDelay,ShakeStart,tVirtCoulombMax,tVirtExchangeMin,MaxMinFac,tMaxHLGap,tHijSqrdMin,OneElWeight,DiagMaxMinFac,OneElMaxMinFac
    USE SystemData, only : tDiagonalizehij,tHFSingDoubExcMax,tSpinOrbs,tReadInCoeff,tUseMP2VarDenMat,tStoreSpinOrbs,tROHF,tFindCINatOrbs
    USE Logging , only : tROHistogramAll,tROFciDump,tROHistER,tROHistOffDiag,tROHistDoubExc,tROHistSingExc,tROHistOnePartOrbEn,tROHistOneElInts,tROHistVirtCoulomb
    USE Logging , only : tPrintInts,tTruncRODump,NoFrozenVirt
    USE OneEInts , only : TMAT2D
    USE SymData , only : TwoCycleSymGens,SymLabelList,SymLabelCounts
    USE Timing , only : end_timing,print_timing_report
    USE Soft_exit, only : test_SOFTEXIT
    IMPLICIT NONE
    INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
    INTEGER , ALLOCATABLE :: Lab(:,:),LabVirtOrbs(:),LabOccOrbs(:),SymLabelCounts2(:,:),SymLabelList2(:),SymLabelListInv(:),SymOrbs(:),SymLabelList3(:)
    REAL*8 , ALLOCATABLE :: CoeffT1(:,:),CoeffCorT2(:,:),CoeffUncorT2(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:),ArrNew(:,:),TMAT2DTemp(:,:),TMAT2DRot(:,:),TMAT2DPartRot01(:,:),TMAT2DPartRot02(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:),UMATTemp01(:,:,:,:),UMATTemp02(:,:,:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:),ForceCorrect(:,:),Correction(:,:),ShakeLambdaNew(:),ConstraintCor(:)
    REAL*8 , ALLOCATABLE :: Constraint(:),ShakeLambda(:),DerivConstrT1(:,:,:),DerivConstrT2(:,:,:),DerivConstrT1T2(:,:),DerivConstrT1T2Diag(:),FourIndInts(:,:,:,:)
    REAL*8 , ALLOCATABLE :: TwoIndInts01(:,:,:,:),TwoIndInts02(:,:,:,:),ThreeIndInts01(:,:,:,:),FourIndInts02(:,:,:,:)
    REAL*8 , ALLOCATABLE :: ThreeIndInts02(:,:,:,:),ThreeIndInts03(:,:,:,:),ThreeIndInts04(:,:,:,:)  
    REAL*8 , ALLOCATABLE :: DiagTMAT2Dfull(:)   
    REAL*8 , ALLOCATABLE :: TwoIndIntsER(:,:,:),ThreeIndInts01ER(:,:),ThreeIndInts02ER(:,:),FourIndIntsER(:)
    INTEGER :: TwoIndIntsERTag,ThreeIndInts01ERTag,ThreeIndInts02ERTag,FourIndIntsERTag,SymLabelList3Tag,SymOrbsTag
    INTEGER :: TwoIndInts01Tag,TwoIndInts02Tag,ThreeIndInts01Tag,ThreeIndInts02Tag,ThreeIndInts03Tag,ThreeIndInts04Tag,FourIndInts02Tag
    INTEGER :: LowBound02,HighBound02,TMAT2DTempTag,TMAT2DRotTag,TMAT2DPartRot01Tag,TMAT2DPartRot02Tag
    INTEGER :: LabTag,ForceCorrectTag,CorrectionTag,SpatOrbs,FourIndIntsTag,ArrNewTag,UMATTemp01Tag,UMATTemp02Tag,ShakeIterInput,NoOrbs,NoOcc
    INTEGER :: CoeffT1Tag,CoeffCorT2Tag,CoeffUncorT2Tag,LambdasTag,DerivCoeffTag,DerivLambdaTag,Iteration,TotNoConstraints,ShakeLambdaNewTag
    INTEGER :: ShakeLambdaTag,ConstraintTag,ConstraintCorTag,DerivConstrT1Tag,DerivConstrT2Tag,DerivConstrT1T2Tag,DerivConstrT1T2DiagTag
    INTEGER :: LabVirtOrbsTag,LabOccOrbsTag,MinOccVirt,MaxOccVirt,MinMZ,MaxMZ,SymLabelCounts2Tag,SymLabelList2Tag,SymLabelListInvTag,error,LowBound,HighBound
    INTEGER :: NoInts01,NoInts02,NoInts03,NoInts04,NoInts05,NoInts06,DiagTMAT2DfullTag,NoRotOrbs
    LOGICAL :: tNotConverged,tInitIntValues
    REAL*8 :: OrthoNorm,ERPotEnergy,HijSqrdPotEnergy,OffDiagPotEnergy,CoulPotEnergy,PotEnergy,Force,TwoEInts,DistCs,OrthoForce,DistLs,LambdaMag,PEInts,PEOrtho
    REAL*8 :: ForceInts,TotCorrectedForce
    REAL*8 :: ijOccVirtPotEnergy,EpsilonMin,MaxTerm
    REAL*8 :: DiagOneElPotInit,ERPotInit,ijVirtOneElPotInit,ijVirtCoulPotInit,ijVirtExchPotInit
    REAL*8 :: singCoulijVirtInit,singExchijVirtInit,singCoulconHFInit,singExchconHFInit,ijklPotInit,ijklantisymPotInit,ijOccVirtOneElPotInit,ijOccVirtCoulPotInit,ijOccVirtExchPotInit
    REAL*8 :: OrthoFac=1.D0,ROHistSing(2,4002),ROHistOffDiag(2,4002),ROHistDoubExc(2,4002),ROHistER(2,4002),ROHistHijVirt(2,4002),ROHistHijOccVirt(2,4002),ROHistHii(2,4002)
    REAL*8 :: ROHistOnePartOrbEn(2,4002),ROHistDCijOcklVir(2,4002),ROHistDEijOcklVir(2,4002),ROHistDCijklVir(2,4002),ROHistDEijklVir(2,4002)
    REAL*8 :: ROHistSCikOcjVir(2,4002),ROHistSEikOcjVir(2,4002),ROHistSCkOcijVir(2,4002),ROHistSEkOcijVir(2,4002),ROHistSCijkVir(2,4002),ROHistSEijkVir(2,4002)
    REAL*8 :: ROHistSASikOcjVir(2,4002),ROHistSASkOcijVir(2,4002),ROHistSASijkVir(2,4002),ROHistASijklVir(2,4002),ROHistASijOcklVir(2,4002)
    TYPE(timer), save :: Rotation_Time,FullShake_Time,Shake_Time,Findtheforce_Time,Transform2ElInts_Time,findandusetheforce_time,CalcDerivConstr_Time,TestOrthoConver_Time
! In this routine, alpha (a), beta (b), gamma (g) and delta (d) refer to the unrotated (HF) orbitals where possible such that < a b | g d > is an unrotated four index integral.   
! For the rotated orbitals, the letter i,j,k and l are generally used, i.e. < i j | k l > refers to a transformed four index integral.
! Differentiation of the potential energy (to find the force) is done with respect to coefficient c(z,m) (or c(a,m)), where zeta (z) or a refers to the HF index, and m to the rotated.
    

    contains

    SUBROUTINE RotateOrbs()
        INTEGER :: i,j   


        IF(iProcIndex.eq.Root) THEN

! If we are reading in our own transformation matrix (coeffT1) don't need a lot of the initialisation stuff.
            IF(tReadInCoeff.or.tUseMP2VarDenMat.or.tFindCINatOrbs) THEN

                tNotConverged=.false.
                CALL UseMP2VarDenMat()

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

            CALL FLUSH(6)
            CALL FLUSH(12)
        ENDIF


    END SUBROUTINE RotateOrbs


    SUBROUTINE UseMP2VarDenMat()
! This routine simply takes a transformation matrix and rotates the integrals to produce a new FCIDUMP file.
! In one case the transformation matrix is read in from a file TRANSFORMMAT.
! In the other, the transformation matrix is calculated from the MP2 variational density matrix.

! MP2VDM = D2_ab = sum_ijc [ t_ij^ac ( 2 t_ij^bc - t_ji^bc ) ]
! Where :  t_ij^ac = - < ab | ij > / ( E_a - E_i + E_b - Ej )
! Ref : J. Chem. Phys. 131, 034113 (2009) - note: in Eqn 1, the cb indices are the wrong way round (should be bc).
        INTEGER :: i,j,k,l,a,ierr,MinReadIn,MaxReadIn
        CHARACTER(len=*) , PARAMETER :: this_routine='UseMP1VarDenMat'


        IF(tStoreSpinOrbs) THEN 
            IF(.not.tSpinOrbs) THEN
                WRITE(6,*) "FCIDUMP and therefore UMAT are in spin orbitals - turning on SPINORBS keyword."
                tSpinOrbs=.true.
            ENDIF
            IF(tROHF) CALL Stop_All(this_routine,"Cannot compress open shell systems into spatial orbitals when rotating, turn off ROHF.")
        ENDIF

!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
!        CALL FLUSH(6)

        SpatOrbs=nBasis/2
        IF(tSpinOrbs) THEN
            NoOrbs=nBasis
            NoOcc=NEl
            MinReadIn=1
            MaxReadIn=nBasis
            IF(tRotateVirtOnly) MinReadIn=NEl+1
            IF(tRotateOccOnly) MaxReadIn=NEl
            do i=1,NoOrbs
                ARR(BRR(i),2)=ARR(i,1)
            enddo
        ELSE
            NoOrbs=SpatOrbs
            NoOcc=NEl/2
            MinReadIn=1
            MaxReadIn=SpatOrbs
            IF(tRotateVirtOnly) MinReadIn=(NEl/2)+1
            IF(tRotateOccOnly) MaxReadIn=NEl/2
        ENDIF
        NoRotOrbs=NoOrbs

! Only virtual orbitals need to be transformed.

        ALLOCATE(SymLabelList2(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList2',NoOrbs,4,this_routine,SymLabelList2Tag,ierr)
        SymLabelList2(:)=0                     
        ALLOCATE(SymLabelList3(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3',NoOrbs,4,this_routine,SymLabelList3Tag,ierr)
        SymLabelList3(:)=0                     
 
        ALLOCATE(SymLabelListInv(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelListInv',NoOrbs,4,this_routine,SymLabelListInvTag,ierr)
        SymLabelListInv(:)=0                     


        IF(tReadInCoeff) THEN
! No symmetry, so no reordering of the orbitals - symlabellist just goes from 1-NoOrbs.        
! When we are just reading in the coefficients and transforming, it does not matter about the ordering of the orbitals.
            do i=1,NoOrbs
                SymLabelList2(i)=i
                SymLabelListInv(i)=i
            enddo
        ELSEIF(tUseMP2VarDenMat) THEN
! When we are calculating the MP2VDM ourselves, we need to be able to run over occupied and virtual separately, so the 
! orbitals are always separated.
! Also, we want to have the option of maintaining symmetry, so SymLabelList2 and SymLabelCounts2 are both constructed so 
! that the orbitals are labelled by symmetry within the occupied and virtual.
            tSeparateOccVirt=.true.
            CALL InitSymmArrays()

        ELSEIF(tFindCINatOrbs) THEN

            ! Doesn't like being allocated outside this module.
            IF(tSpinOrbs) THEN
                ALLOCATE(SymLabelCounts2(2,32),stat=ierr)
                CALL LogMemAlloc('SymLabelCounts2',2*32,4,this_routine,SymLabelCounts2Tag,ierr)
                SymLabelCounts2(:,:)=0
                ! first 8 refer to the occupied, and the second to the virtual beta spin.
                ! third and fourth to the occupied and virtual alpha spin.

            ELSE
                ALLOCATE(SymLabelCounts2(2,16),stat=ierr)
                CALL LogMemAlloc('SymLabelCounts2',2*16,4,this_routine,SymLabelCounts2Tag,ierr)
                SymLabelCounts2(:,:)=0
                ! first 8 refer to the occupied, and the second to the virtual.
            ENDIF

!            tSeparateOccVirt=.true.
            CALL SetupNatOrbLabels() 

        ENDIF

!        WRITE(6,*) 'SymLabelList'
!        do i=1,NoOrbs
!            WRITE(6,*) SymLabelList2(i)
!        enddo
!        CALL FLUSH(6)
!        stop

!            OPEN(42,file="TRANSFORMMAT",status="old")
!            do i=MinReadIn,MaxReadIn
!                j=SymLabelList2(i)
!                do a=MinReadIn,MaxReadIn
!                    b=SymLabelList2(a)
!                    READ(42,*) CoeffT1(b,j)
!                    READ(42,*) CoeffT1
!                enddo
!            enddo
!            CLOSE(42)
         
! Need to read to convert the UMAT matrix from UMATInd to the appropriate indexing for Transform2ElInts.        
! This just contains all the untransformed orbitals.
        ALLOCATE(UMATTemp01(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('UMATTemp01',NoOrbs**4,8,this_routine,UMATTemp01Tag,ierr)

        CALL CopyAcrossUMAT()
        

! Yet another labelling system, SymLabelList3 is created here.
! This indicates the label of the transformed orbital.
! In the case where we are truncating the space, the transformed orbitals are ordered according to the size of the eigenvalues of the MP2VDM
! matrix when it is diagonalised.  We wish to keep them in this order when transforming the integrals etc, so that when we truncate the last 
! NoFrozenVirt orbitals, we are removing those with the smallest MP2VDM eigenvalues (occupation numbers).
! In the case where no truncation is made however, SymLabelList3 is the same as SymLabelList2, so that the indexes remain the same as previously.
! This allows for the option of going straight into a spawning calc from the rotation, which is not possible when a truncation is performed 
! because of the messed up indices.
        IF(tTruncRODump) THEN
            IF(MOD(NoFrozenVirt,2).ne.0) CALL Stop_All(this_routine,"Must freeze virtual spin orbitals in pairs of 2.")
            IF(tSpinOrbs) THEN
                NoRotOrbs=NoOrbs-NoFrozenVirt
            ELSE
                NoFrozenVirt=NoFrozenVirt/2
                NoRotOrbs=NoOrbs-NoFrozenVirt
            ENDIF            
            do i=1,NoOrbs
                SymLabelList3(i)=i
            enddo
        ELSE
            do i=1,NoOrbs
                SymLabelList3(i)=SymLabelList2(i)
            enddo
        ENDIF


        ALLOCATE(CoeffT1(NoOrbs,NoRotOrbs),stat=ierr)
        CALL LogMemAlloc(this_routine,NoRotOrbs*NoOrbs,8,this_routine,CoeffT1Tag,ierr)
        CoeffT1(:,:)=0.D0
        IF(tSeparateOccVirt) THEN
            do i=1,NoRotOrbs
                CoeffT1(i,i)=1.D0
            enddo
        ENDIF


! The last two indices of these are the transformed and possibly truncated orbitals.        

        ALLOCATE(TwoIndInts01(NoOrbs,NoOrbs,NoRotOrbs,NoRotOrbs),stat=ierr)
        CALL LogMemAlloc('TwoIndInts01',(NoOrbs**4)*(NoRotOrbs**2),8,this_routine,TwoIndInts01Tag,ierr)

        ALLOCATE(FourIndInts(NoRotOrbs,NoRotOrbs,NoRotOrbs,NoRotOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts',NoRotOrbs**4,8,this_routine,FourIndIntsTag,ierr)


        IF(tReadInCoeff) THEN

            WRITE(6,*) "Reading in the transformation matrix from TRANSFORMMAT, and using this to rotate the HF orbitals."

            OPEN(72,FILE='TRANSFORMMAT',status='old')
            READ(72,*) CoeffT1
            CLOSE(72)

!            OPEN(78,FILE='TRANSFORMMATORIG',status='unknown')
!            do i=1,NoOrbs
!                do a=1,NoOrbs
!                    WRITE(78,*) i,a,CoeffT1(i,a)
!                enddo
!            enddo
!            stop

        ELSEIF(tFindCINatOrbs) THEN
            
            CALL FindNatOrbs()
            
            ! Fill the coefficient matrix with the eigenvectors of the OneRDM.
            ! Find out the ordering ...need to read in according to symlabellist2, so that the transformation is all o.k.
            ! If both the HF and transformed indices are done this way, the symmetries and everything should be fine.
!            do i=1,NoRotOrbs
!                i2=SymLabelList2(i)
!                do a=1,NoOrbs
!                    a2=SymLabelList2(a)
!                    Coeff(a,i)=OneRDM(a2,i2)
!                enddo
!            enddo
            

        ELSEIF(tUseMP2VarDenMat) THEN
! This bit generates the MP2 variational density matrix, and uses this as the transformation matrix (CoeffT1).        
    
            WRITE(6,*) "Calculating the MP2 vartiational density matrix, and using this to rotate the HF orbitals."
            CALL CalcMP2VarDenMat()

        ENDIF

! Then, transform2ElInts
        WRITE(6,*) 'Transforming the four index integrals'
        CALL Transform2ElInts()

        WRITE(6,*) 'Re-calculating the fock matrix'
        CALL CalcFOCKMatrix()

        WRITE(6,*) 'Refilling the UMAT and TMAT2D'
! The ROFCIDUMP is also printed out in here.        
        CALL RefillUMATandTMAT2D()        
 
        CALL FLUSH(6)
 
! If a truncation is being made, the new basis will not be in the correct energetic ordering - this does not matter, as we
! never go straight into a spawning and they will be reordered when the ROFCIDUMP file is read in again. 
        CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)

        DEALLOCATE(CoeffT1)
        CALL LogMemDeAlloc(this_routine,CoeffT1Tag)
        DEALLOCATE(SymLabelList2)
        CALL LogMemDeAlloc(this_routine,SymLabelList2Tag)
        DEALLOCATE(SymLabelListInv)
        CALL LogMemDeAlloc(this_routine,SymLabelListInvTag)
        DEALLOCATE(UMATTemp01)
        CALL LogMemDeAlloc(this_routine,UMATTemp01Tag)
        DEALLOCATE(TwoIndInts01)
        CALL LogMemDeAlloc(this_routine,TwoIndInts01Tag)
        DEALLOCATE(FourIndInts)
        CALL LogMemDeAlloc(this_routine,FourIndIntsTag)


    END SUBROUTINE UseMP2VarDenMat 

    

    SUBROUTINE CalcMP2VarDenMat()
        INTEGER :: w,k,l,a,b,c,i,j,a2,b2,c2,i2,j2,ierr,LWORK2,INFO,MP2VDMTag,WORK2Tag,NoSpins,Skip,NoSpinArray,SymLabelListSpinTag
        INTEGER :: MP2VDMSymTag,MP2VDMTempTag,EvaluesSymTag,NoSymBlock,SymStartInd,Sym,EvaluesSpinTag,MP2VDMSpinTag,SymOrbsSpinTag,FinalSpin
        INTEGER :: OrigNoSymBlock
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        REAL*8 :: Evalues(NoOrbs),MP2VDMSum,SumTrace,SumDiagTrace
        REAL*8 , ALLOCATABLE :: MP2VDM(:,:),MP2VDMSym(:,:),MP2VDMTemp(:,:),EvaluesSym(:),WORK2(:),MP2VDMSpin(:,:),EvaluesSpin(:)
        INTEGER , ALLOCATABLE :: SymOrbsSpin(:),SymLabelListSpin(:)

        ALLOCATE(MP2VDM(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('MP2VDM',NoOrbs**2,8,this_routine,MP2VDMTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Error assigning memory to MP2VDM.')
        MP2VDM(:,:)=0.D0

        Evalues(:)=0.D0

! Calculating the MP2VDM (D2_ab) matrix whose eigenvectors become the transformation matrix.        
! The eigenvalues are the occupation numbers of the new orbitals.  These should decrease exponentially so that when we remove the 
! orbitals with small occupation numbers we should have little affect on the energy.
        do a2=NoOcc+1,NoOrbs
            a=SymLabelList2(a2)
            do b2=NoOcc+1,NoOrbs
                b=SymLabelList2(b2)
                MP2VDMSum=0.D0

                do c2=NoOcc+1,NoOrbs
                    c=SymLabelList2(c2)
                    do i2=1,NoOcc
                        i=SymLabelList2(i2)
                        do j2=1,NoOcc
                            j=SymLabelList2(j2)
                            IF(tSpinOrbs) THEN
                                IF((ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)).eq.0.D0) THEN
                                    IF((REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)).ne.0.D0) THEN
                                        WRITE(6,*) i,j,a,c,REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)
                                        CALL Stop_All(this_routine,"Dividing a non-zero by zero.")
                                    ENDIF
                                ENDIF
!                                CALL Stop_All(this_routine,"have not set UseMP2VarDenMat up for spin orbitals yet.")
                                MP2VDMSum=MP2VDMSum+&
                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (2.D0*(REAL(UMAT(UMatInd(b,c,i,j,0,0))%v,8))) )/&
                                            &( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) * (ARR(i,2)+ARR(j,2)-ARR(b,2)-ARR(c,2)) ) )
                                MP2VDMSum=MP2VDMSum-&                                            
                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (REAL(UMAT(UMatInd(c,b,i,j,0,0))%v,8)) )/&
                                            &( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) * (ARR(i,2)+ARR(j,2)-ARR(c,2)-ARR(b,2)) ) )
 
                            ELSE
!                                MP2VDMSum=MP2VDMSum+&
!                                            &(UMATTemp01(i,a,j,c)/(ARR(2*i,1)-ARR(2*a,1)+ARR(2*j,1)-ARR(2*c,1)))*&
!                                            &(2.D0*UMATTemp01(i,b,j,c)/(ARR(2*i,1)-ARR(2*b,1)+ARR(2*j,1)-ARR(2*c,1)))&
!                                            &-(UMATTemp01(i,a,j,c)/(ARR(2*i,1)-ARR(2*a,1)+ARR(2*j,1)-ARR(2*c,1)))*&
!                                            &(UMATTemp01(j,b,i,c)/(ARR(2*i,1)-ARR(2*c,1)+ARR(2*j,1)-ARR(2*b,1)))

                                MP2VDMSum=MP2VDMSum+&
                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (2.D0*(REAL(UMAT(UMatInd(b,c,i,j,0,0))%v,8))) )/&
                                            &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*b,2)-ARR(2*c,2)) ) )
                                MP2VDMSum=MP2VDMSum-&                                            
                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (REAL(UMAT(UMatInd(c,b,i,j,0,0))%v,8)) )/&
                                            &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*c,2)-ARR(2*b,2)) ) )
                            ENDIF
                        enddo
                    enddo
                enddo
                MP2VDM(a2,b2)=MP2VDMSum
                MP2VDM(b2,a2)=MP2VDMSum

            enddo
        enddo
        WRITE(6,*) 'Finished filling MP2VDM'

        SumTrace=0.D0
        do i=1,NoOrbs
            SumTrace=SumTrace+MP2VDM(i,i)
            do j=1,NoOrbs
                WRITE(6,*) MP2VDM(j,i)
            enddo
        enddo
        CALL FLUSH(6)

        
!        OPEN(83,file='D2MATRIX-sym',status='unknown')
!        do a=1,NoOrbs
!            do b=1,NoOrbs
!                WRITE(83,'(4I3,F20.10)') a,INT(G1(SymLabelList2(a)*2)%sym%S,4),b,INT(G1(SymLabelList2(b)*2)%sym%S,4),MP2VDM(a,b)
!            enddo
!        enddo
!        CLOSE(83)

        WRITE(6,*) 'Calculating eigenvectors and eigenvalues of MP2VDM'
        CALL FLUSH(6)
        IF(lNoSymmetry) THEN
            LWORK2=3*(NoOrbs-NoOcc)+1

            ALLOCATE(WORK2(LWORK2),stat=ierr)
            CALL LogMemAlloc('WORK2',LWORK2,8,this_routine,WORK2Tag,ierr)

            CALL DSYEV( 'V', 'L', (NoOrbs-NoOcc), MP2VDM((NoOcc+1):NoOrbs,(NoOcc+1):NoOrbs), (NoOrbs-NoOcc), Evalues((NoOcc+1):NoOrbs), WORK2, LWORK2, INFO )
            ! MP2VDM goes in as the D2 matrix, comes out as the eigenvectors.
            ! Evalues come out in ascending order.
            ! The transformed orbitals are now all jumbled up as the eigenvectors are ordered according to the eigenvalues.  When symmetry is turned off, this
            ! doesn't matter.  However we like to see the eigenvalues from largest to smallest so that we can see the occupation numbers decrease, so they are 
            ! reordered to be in descending order below.

            DEALLOCATE(WORK2)
            CALL LogMemDealloc(this_routine,WORK2Tag)

            WRITE(6,*) 'Matrix diagonalised'
            CALL FLUSH(6)

            ! Eigenvalues come out in descending order, with the corresponding eigenvectors.
            CALL SortEvecbyEval((NoOrbs-NoOcc),Evalues((NoOcc+1):NoOrbs),(NoOrbs-NoOcc),MP2VDM((NoOcc+1):NoOrbs,(NoOcc+1):NoOrbs))

        ELSE
! If we want to maintain the symmetry, we cannot have all the orbitals jumbled up when the diagonaliser reorders the eigenvectors.
! Must instead feed each symmetry block in separately.
! This means that although the transformed orbitals are jumbled within the symmetry blocks, the symmetry labels are all that are relevant and this is unaffected.
     
            ALLOCATE(MP2VDMTemp(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('MP2VDMTemp',NoOrbs**2,8,this_routine,MP2VDMTempTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Error allocating memory to MP2VDMTemp.')
            MP2VDMTemp(:,:)=0.D0

            ! If we are using spin orbitals, need to feed in the alpha and beta spins separately.
            ! Otherwise these jumble up and the final ordering is uncorrect. 
            ! There should be no non-zero values between these, but can put a check in for this.
            IF(tSpinOrbs) THEN
                ! Odd indices are alpha, even beta.
                NoSpins=2
                Skip=2
            ELSE
                NoSpins=1
                Skip=1
            ENDIF
          
            Sym=0
            LWORK2=-1
            do while (Sym.le.7)

                NoSymBlock=SymLabelCounts2(2,Sym+9)

                SymStartInd=SymLabelCounts2(1,Sym+9)-1
                ! This is one less than the index that the symmetry starts, so that when we run through i=1,..., we can
                ! start at SymStartInd+i.

                OrigNoSymBlock=NoSymBlock

                ! Run through twice (once for each spin) for spin orbitals, and just once otherwise.
                do w=1,NoSpins
                    IF(tSpinOrbs) THEN
                        IF(w.eq.1) THEN
                            NoSymBlock=CEILING(REAL(OrigNoSymBlock)/2.D0)
                            FinalSpin=OrigNoSymBlock
                        ELSEIF(w.eq.2) THEN
                            NoSymBlock=FLOOR(REAL(OrigNoSymBlock)/2.D0)
                            SymStartInd=SymStartInd+1
                            FinalSpin=OrigNoSymBlock-1
                        ENDIF
                        ! If OrigNoSymBlock is even, NoSymBlock will be the same for both w.
                    ELSE
                        FinalSpin=NoSymBlock
                    ENDIF

                    IF(NoSymBlock.gt.1) THEN

                        ALLOCATE(MP2VDMSym(NoSymBlock,NoSymBlock),stat=ierr)
                        CALL LogMemAlloc('MP2VDMSym',NoSymBlock**2,8,this_routine,MP2VDMSymTag,ierr)
                        IF(ierr.ne.0) CALL Stop_All(this_routine,'Error allocating memory to MP2VDMSym.')
                        MP2VDMSym(:,:)=0.D0

                        ALLOCATE(EvaluesSym(NoSymBlock),stat=ierr)
                        CALL LogMemAlloc('EvaluesSym',NoSymBlock,8,this_routine,EvaluesSymTag,ierr)
                        IF(ierr.ne.0) CALL Stop_All(this_routine,'Error allocating memory to EvaluesSym.')
                        EvaluesSym(:)=0

                        LWORK2=3*NoSymBlock+1
                        ALLOCATE(WORK2(LWORK2),stat=ierr)
                        CALL LogMemAlloc('WORK2',LWORK2,8,this_routine,WORK2Tag,ierr)
                        IF(ierr.ne.0) CALL Stop_All(this_routine,'Error allocating memory to WORK2.')

                        l=0
                        do j=1,FinalSpin,Skip
                            l=l+1
                            k=0
                            do i=1,FinalSpin,Skip
                                k=k+1
                                MP2VDMSym(k,l)=MP2VDM(SymStartInd+i,SymStartInd+j)
                            enddo
                        enddo


                        WRITE(6,*) '*****'
                        WRITE(6,*) 'Symmetry ',Sym, 'with w ',w,' has ',NoSymBlock,' orbitals .'
                        WRITE(6,*) 'The MP2VDM for this symmetry block is '
                        do i=1,NoSymBlock
                            do j=1,NoSymBlock
                                WRITE(6,'(F20.10)',advance='no') MP2VDMSym(j,i)
                            enddo
                            WRITE(6,*) ''
                        enddo

                        CALL DSYEV('V','L',NoSymBlock,MP2VDMSym,NoSymBlock,EvaluesSym,WORK2,LWORK2,ierr)
                        ! MP2VDMSym goes in as the original MP2VDMSym, comes out as the eigenvectors (Coefficients).
                        ! EvaluesSym comes out as the eigenvalues in ascending order.

                        WRITE(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are ,'
                        do i=1,NoSymBlock
                            WRITE(6,'(F20.10)',advance='no') EvaluesSym(i)
                        enddo
                        WRITE(6,*) ''
                        WRITE(6,*) 'These go from orbital ,',SymStartInd+1,' to ',SymStartInd+NoSymBlock
                       
                        k=0
                        do i=1,FinalSpin,Skip
                            k=k+1
                            Evalues(SymStartInd+i)=EvaluesSym(k)
                        enddo

                        ! CAREFUL if eigenvalues are put in ascending order, this may not be correct, with the labelling system.
                        ! may be better to just take coefficients and transform TMAT2DRot in transform2elints.
                        ! a check that comes out as diagonal is a check of this routine anyway.

                        WRITE(6,*) 'The eigenvectors (coefficients) for symmtry block ',Sym
                        do i=1,NoSymBlock
                            do j=1,NoSymBlock
                                WRITE(6,'(F20.10)',advance='no') MP2VDMSym(j,i)
                            enddo
                            WRITE(6,*) ''
                        enddo
                     
                        l=0
                        do j=1,FinalSpin,Skip
                            l=l+1
                            k=0
                            do i=1,FinalSpin,Skip
                                k=k+1
                                MP2VDMTemp(SymStartInd+i,SymStartInd+j)=MP2VDMSym(k,l)
                            enddo
                        enddo
                        ! Directly fill the coefficient matrix with the eigenvectors from the diagonalization.

                        DEALLOCATE(WORK2)
                        CALL LogMemDealloc(this_routine,WORK2Tag)

                        DEALLOCATE(MP2VDMSym)
                        CALL LogMemDealloc(this_routine,MP2VDMSymTag)

                        DEALLOCATE(EvaluesSym)
                        CALL LogMemDealloc(this_routine,EvaluesSymTag)

                    ELSEIF(NoSymBlock.eq.1) THEN
                        ! The eigenvalue is the lone value, while the eigenvector is 1.

                        Evalues(SymStartInd+1)=MP2VDM(SymStartInd+1,SymStartInd+1)
                        MP2VDMTemp(SymStartInd+1,SymStartInd+1)=1.D0
                        WRITE(6,*) '*****'
                        WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
                        WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to MP2VDM'
                    ENDIF

                enddo

                Sym=Sym+1
            enddo

            WRITE(6,*) 'Matrix diagonalised'
            CALL FLUSH(6)


! Here, if symmetry is kept, we are going to have to reorder the eigenvectors according to the size of the eigenvalues, while taking
! the orbital labels (and therefore symmetries) with them. This will be put back into MP2VDM from MP2VDMTemp.

! Want to reorder the eigenvalues from largest to smallest, taking the eigenvectors with them and the symmetry as well.  
    
            IF(tTruncRODump) THEN
                ! If we are truncating, the orbitals stay in this order, so we want to take their symmetries with them.

                ALLOCATE(SymOrbs(NoOrbs),stat=ierr)
                CALL LogMemAlloc('SymOrbs',NoOrbs,4,this_routine,SymOrbsTag,ierr)
                SymOrbs(:)=0

                do i=1,NoOrbs
                    SymOrbs(i)=INT(G1(SymLabelList2(i)*2)%sym%S,4)
                enddo

                IF(tSpinOrbs) THEN
                    MP2VDM(:,:)=0.D0
                    do w=1,2
                        IF(w.eq.1) NoSpinArray=CEILING(REAL(NoOrbs-NoOcc)/2.D0)
                        IF(w.eq.2) NoSpinArray=FLOOR(REAL(NoOrbs-NoOcc)/2.D0)
         
                        ALLOCATE(EvaluesSpin(NoSpinArray),stat=ierr)
                        CALL LogMemAlloc('EvaluesSpin',NoSpinArray,8,this_routine,EvaluesSpinTag,ierr)
                        EvaluesSpin(:)=0.D0
 
                        ALLOCATE(MP2VDMSpin(NoSpinArray,NoSpinArray),stat=ierr)
                        CALL LogMemAlloc('MP2VDMSpin',NoSpinArray**2,8,this_routine,MP2VDMSpinTag,ierr)
                        MP2VDMSpin(:,:)=0.D0
 
                        ALLOCATE(SymOrbsSpin(NoSpinArray),stat=ierr)
                        CALL LogMemAlloc('SymOrbsSpin',NoSpinArray,4,this_routine,SymOrbsSpinTag,ierr)
                        SymOrbsSpin(:)=0

                        IF(MOD((NoOcc+w),2).ne.0) THEN
                            FinalSpin=NoOrbs-1
                        ELSE
                            FinalSpin=NoOrbs
                        ENDIF

                        k=0
                        do i=NoOcc+w,FinalSpin,2
                            k=k+1
                            l=0
                            do j=NoOcc+w,FinalSpin,2
                                l=l+1
                                MP2VDMSpin(l,k)=MP2VDMTemp(j,i)
                            enddo
                            EvaluesSpin(k)=Evalues(i)
                            SymOrbsSpin(k)=SymOrbs(i)
                        enddo
      
                        CALL SortEvecbyEvalPlus1(NoSpinArray,EvaluesSpin(1:NoSpinArray),NoSpinArray,MP2VDMSpin(1:NoSpinArray,&
                                                &1:NoSpinArray),SymOrbsSpin(1:NoSpinArray))

                        k=0
                        do i=NoOcc+w,FinalSpin,2
                            k=k+1
                            l=0
                            do j=NoOcc+w,FinalSpin,2
                                l=l+1
                                MP2VDM(j,i)=MP2VDMSpin(l,k)
                            enddo
                            Evalues(i)=EvaluesSpin(k)
                            SymOrbs(i)=SymOrbsSpin(k)
                        enddo
         
                        DEALLOCATE(EvaluesSpin)
                        CALL LogMemDeAlloc(this_routine,EvaluesSpinTag)
                        DEALLOCATE(MP2VDMSpin)
                        CALL LogMemDeAlloc(this_routine,MP2VDMSpinTag)
                        DEALLOCATE(SymOrbsSpin)
                        CALL LogMemDeAlloc(this_routine,SymOrbsSpinTag)

                    enddo
                ELSE 
                    CALL SortEvecbyEvalPlus1((NoOrbs-NoOcc),Evalues((NoOcc+1):NoOrbs),(NoOrbs-NoOcc),MP2VDMTemp((NoOcc+1):NoOrbs,&
                                            &(NoOcc+1):NoOrbs),SymOrbs((NoOcc+1):NoOrbs))
                    MP2VDM(:,:)=0.D0                                            
                    MP2VDM(:,:)=MP2VDMTemp(:,:)                                            
                ENDIF

            ELSE
                ! If we are not truncating, they orbitals get put back into their original order, so the symmetry information is still 
                ! correct, no need for the SymOrbs array.
                ! Instead, just take the labels of SymLabelList3 with them.

                IF(tSpinOrbs) THEN
                ! Need to reorder the alpha and beta separately, and put them back into the total array so that they alternate.                    
                    MP2VDM(:,:)=0.D0
                    do w=1,2
                        
                        IF(w.eq.1) NoSpinArray=CEILING(REAL(NoOrbs-NoOcc)/2.D0)
                        IF(w.eq.2) NoSpinArray=FLOOR(REAL(NoOrbs-NoOcc)/2.D0)
         
!                        WRITE(6,*) 'w',w
!                        WRITE(6,*) 'NoSpinArray',NoSpinArray

                        ALLOCATE(EvaluesSpin(NoSpinArray),stat=ierr)
                        CALL LogMemAlloc('EvaluesSpin',NoSpinArray,8,this_routine,EvaluesSpinTag,ierr)
                        EvaluesSpin(:)=0.D0
                        ALLOCATE(MP2VDMSpin(NoSpinArray,NoSpinArray),stat=ierr)
                        CALL LogMemAlloc('MP2VDMSpin',NoSpinArray**2,8,this_routine,MP2VDMSpinTag,ierr)
                        MP2VDMSpin(:,:)=0.D0
                        ALLOCATE(SymLabelListSpin(NoSpinArray),stat=ierr)
                        CALL LogMemAlloc('SymLabelListSpin',NoSpinArray,4,this_routine,SymLabelListSpinTag,ierr)
                        SymLabelListSpin(:)=0

                        ! IF we start on an odd number, end on an odd number.
                        IF(MOD((NoOcc+w),2).ne.0) THEN
                            FinalSpin=NoOrbs-1
                        ELSE
                            FinalSpin=NoOrbs
                        ENDIF
                        k=0
                        do i=NoOcc+w,FinalSpin,2
                            k=k+1
                            l=0
                            do j=NoOcc+w,FinalSpin,2
                                l=l+1
                                MP2VDMSpin(l,k)=MP2VDMTemp(j,i)
                            enddo
                            EvaluesSpin(k)=Evalues(i)
                            SymLabelListSpin(k)=SymLabelList3(i)
                        enddo
                        
                        CALL SortEvecbyEvalPlus1(NoSpinArray,EvaluesSpin(1:NoSpinArray),NoSpinArray,MP2VDMSpin(1:NoSpinArray,&
                                            &1:NoSpinArray),SymLabelListSpin(1:NoSpinArray))

                        k=0
                        do i=NoOcc+w,FinalSpin,2
                            k=k+1
                            l=0
                            do j=NoOcc+w,FinalSpin,2
                                l=l+1
                                MP2VDM(j,i)=MP2VDMSpin(l,k)
                            enddo
                            Evalues(i)=EvaluesSpin(k)
                            SymLabelList3(i)=SymLabelListSpin(k)
                        enddo

                        DEALLOCATE(EvaluesSpin)
                        CALL LogMemDeAlloc(this_routine,EvaluesSpinTag)
                        DEALLOCATE(MP2VDMSpin)
                        CALL LogMemDeAlloc(this_routine,MP2VDMSpinTag)
                        DEALLOCATE(SymLabelListSpin)
                        CALL LogMemDeAlloc(this_routine,SymLabelListSpinTag)

                    enddo
 
                ELSE
                    CALL SortEvecbyEvalPlus1((NoOrbs-NoOcc),Evalues((NoOcc+1):NoOrbs),(NoOrbs-NoOcc),MP2VDMTemp((NoOcc+1):NoOrbs,&
                                            &(NoOcc+1):NoOrbs),SymLabelList3((NoOcc+1):NoOrbs))
                    MP2VDM(:,:)=0.D0                                            
                    MP2VDM(:,:)=MP2VDMTemp(:,:)                                            
                ENDIF 
            ENDIF

! If we are not truncating, there is no need to order the orbitals in terms of the eigenvalue, so they may just stay as they are.             
! Accept that we want to look at the eigenvalues in descending order, so just order these to write out, and copy over the MP2VDMTemp
! to MP2VDM.


        ENDIF

        SumDiagTrace=0.D0
        do i=1,NoOrbs
            SumDiagTrace=SumDiagTrace+Evalues(i)
        enddo
        IF((ABS(SumDiagTrace-SumTrace)).gt.1E-10) THEN
            WRITE(6,*) 'The sum of the diagonal MP2VDM values : ', SumTrace
            WRITE(6,*) 'The sum of the eigenvalues of the MP2VDM matrix : ',SumDiagTrace
            CALL Stop_All('CalcMP2VarDenMat','The sum of the trace of the MP2VDM before diagonalisation does not equal that after.')
        ENDIF

!        do i=1,NoOrbs
!            WRITE(6,*) Evalues(i)
!        enddo
!        do i=1,NoOrbs
!            WRITE(6,*) MP2VDM(:,i)
!        enddo
!        CALL FLUSH(6)
!        stop

        do i=1,NoOcc
            Evalues(i)=1.D0
            MP2VDM(i,i)=1.D0
        enddo

        OPEN(73,FILE='D2EVS',status='unknown')
        WRITE(73,*) NoOrbs
        do i=1,NoOcc
            WRITE(73,*) NoOrbs-i+1,Evalues(i)
        enddo
        do i=NoOcc+1,NoOrbs
            WRITE(73,*) NoOrbs-i+1,Evalues(i)
        enddo
        CLOSE(73)

        CALL HistEvalues(Evalues,MP2VDM)

        WRITE(6,*) 'Eigen-values: '
        do i=1,NoOrbs
            WRITE(6,*) Evalues(i)
        enddo

        WRITE(6,*) 'MP2VDM matrix'
        do i=1,NoOrbs
            WRITE(6,*) MP2VDM(:,i)
        enddo

!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                MP2VDMInv(i,j)=MP2VDM(j,i)
!                CoeffT1(i,j)=MP2VDMInv(i,j)
!            enddo
!        enddo
!        CoeffT1(:,:)=MP2VDMInv(:,:)
        CoeffT1(:,:)=0.D0
        do i=1,NoRotOrbs
            CoeffT1(:,i)=MP2VDM(:,i)
        enddo



        OPEN(74,FILE='TRANSFORMMAT',status='unknown')
        do i=1,NoOrbs
            do j=1,NoRotOrbs
                WRITE(74,*) i,j,CoeffT1(i,j)
            enddo
        enddo
        CLOSE(74)
        
        DEALLOCATE(MP2VDM)
        CALL LogMemDeAlloc(this_routine,MP2VDMTag)

        IF(.not.lNoSymmetry) THEN
            DEALLOCATE(MP2VDMTemp)
            CALL LogMemDeAlloc(this_routine,MP2VDMTempTag)
        ENDIF



    END SUBROUTINE CalcMP2VarDenMat


    SUBROUTINE HistEvalues(Evalues,MP2VDM)
        INTEGER :: i,k,NoEvalues,a,b,x
        REAL*8 :: EvaluesCount((NoOrbs-NoOcc),2),Evalues(1:NoOrbs)
        REAL*8 :: OrbEnergies(1:NoOrbs),EvalueEnergies(1:NoOrbs)
        REAL*8 :: MP2VDM(1:NoOrbs,1:NoOrbs),SumEvalues

        EvaluesCount(:,:)=0.D0

        k=1
        EvaluesCount(k,1)=Evalues(NoOcc+1)
        EvaluesCount(k,2)=1.D0
        do i=NoOcc+2,NoOrbs
            IF((ABS(Evalues(i)-Evalues(i-1))).ge.(1E-10)) THEN
                k=k+1
                EvaluesCount(k,1)=Evalues(i)
                EvaluesCount(k,2)=1.D0
            ELSE
                EvaluesCount(k,2)=EvaluesCount(k,2)+1.D0
            ENDIF
        enddo
        NoEvalues=k

        
        OPEN(73,FILE='D2EVS-plot',status='unknown')
        do i=1,NoEvalues
            WRITE(73,*) EvaluesCount(i,1),Evaluescount(i,2)
        enddo
        CLOSE(73)

        
        OPEN(73,FILE='D2EVS-plot-rat',status='unknown')
        k=0
        do i=NoOcc+1,NoOrbs
            k=k+1
            WRITE(73,*) REAL(k)/REAL(NoOrbs-NoOcc),Evalues(i)
        enddo
        CLOSE(73)


        OrbEnergies(:)=0.D0
        EvalueEnergies(:)=0.D0
        do i=1,NoOrbs
            EvalueEnergies(i)=Evalues(i)
! We are only interested in the diagonal elements.            
            do a=1,NoOrbs
                b=SymLabelList2(a)
                IF(tStoreSpinOrbs) THEN
                    OrbEnergies(i)=OrbEnergies(i)+(MP2VDM(a,i)*ARR(b,2)*MP2VDM(a,i))
                ELSE
                    OrbEnergies(i)=OrbEnergies(i)+(MP2VDM(a,i)*ARR(2*b,2)*MP2VDM(a,i))
                ENDIF
            enddo
        enddo
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

        WRITE(6,*) 'OrbEnergies'
        do i=1,NoOrbs
            WRITE(6,*) OrbEnergies(i)
        enddo
        WRITE(6,*) 'EvalueEnergies'
        SumEvalues=0.D0
        do i=1,NoOrbs
            WRITE(6,*) EvalueEnergies(i)
            IF(tSpinOrbs) THEN
                SumEvalues=SumEvalues+EvalueEnergies(i)
            ELSE
                SumEvalues=SumEvalues+(2*EvalueEnergies(i))
            ENDIF
        enddo

        CALL Sort2Real(NoOrbs,OrbEnergies(1:NoOrbs),EvalueEnergies(1:NoOrbs))

        OPEN(73,FILE='EVALUES-ENERGY',status='unknown')
        do i=1,NoOrbs
            WRITE(73,*) OrbEnergies(NoOrbs-i+1),EvalueEnergies(NoOrbs-i+1)
        enddo
        WRITE(73,*) 'Sum of occupation numbers (eigenvalues) = ',SumEvalues
        WRITE(73,*) 'Number of electrons = ',NEl
        CALL FLUSH(73)
        CLOSE(73)
        CALL FLUSH(6)


        OPEN(73,FILE='OccupationTable',status='unknown')
        x=1
        do while (x.le.NoOrbs)
            WRITE(73,'(A16,A5)',advance='no') 'HF Orb En    ','Sym'
            do i=x,x+9
                IF(i.gt.NoOrbs) THEN
                    WRITE(73,*) ''
                    EXIT
                ENDIF
                WRITE(73,'(ES16.6)',advance='no') Evalues(i)
            enddo
            WRITE(73,*) ''

            do a=1,NoOrbs
                b=SymLabelListInv(a)
                IF(tStoreSpinOrbs) THEN
                    WRITE(73,'(F16.10,I5)',advance='no') ARR(a,1),INT(G1(a)%sym%S,4)
                ELSE
                    WRITE(73,'(F16.10,I5)',advance='no') ARR(2*a,1),INT(G1(2*a)%sym%S,4)
                ENDIF
                do i=x,x+9
                    IF(i.gt.NoOrbs) THEN
                        WRITE(73,*) ''
                        EXIT
                    ENDIF
                    WRITE(73,'(F16.10)',advance='no') MP2VDM(b,i)
                enddo
                WRITE(73,*) ''
            enddo
            WRITE(73,*) ''
            x=x+10
        enddo
        CALL FLUSH(73)
        CLOSE(73)
        CALL FLUSH(6)


    END SUBROUTINE HistEvalues


    SUBROUTINE InitLocalOrbs()
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        INTEGER :: ierr

! Writing to output which PE is being maximised/minimised.        
        WRITE(6,*) '*****'
        IF(tERLocalization) THEN
!            IF(tSpinOrbs) CALL Stop_All(this_routine,"Cannot do an ER localisation with spin orbitals.")
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
                CALL FLUSH(6)
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
            CALL FLUSH(6)
            CALL Stop_All(this_routine,"ERROR. TwoCycleSymGens is false.  Symmetry is not abelian.") 
        ENDIF
        IF((tRotateOccOnly.or.tRotateVirtOnly).and.(.not.tSeparateOccVirt)) THEN
            tSeparateOccVirt=.true.
            WRITE(6,*) "NOTE. Cannot rotate only occupied or virtual without first separating them."
            WRITE(6,*) "SEPARATEOCCVIRT keyword is being turned on."
        ENDIF        
        IF((tOffDiagSqrdMax.and.tOffDiagSqrdMin).or.(tOffDiagMax.and.tOffDiagMin)) THEN
            CALL FLUSH(6)
            CALL Stop_All(this_routine,"ERROR. Cannot both maximise and minimise off diagonal elements simultaneously")
        ENDIF
        IF(tOnePartOrbEnMax.and.(.not.tSeparateOccVirt)) THEN
            CALL FLUSH(6)
            CALL Stop_All(this_routine,"ERROR. Cannot currently maximise the one particle orbital energies without separating occupied and virtual.") 
        ENDIF
        WRITE(6,*) "*****"

!Zero values.
        OrthoNorm=0.D0
        ERPotEnergy=0.D0
        PotEnergy=0.D0
        Force=0.D0
        TwoEInts=0.D0
        PEInts=0.D0
        PEOrtho=0.D0
        ForceInts=0.D0
        DistCs=0.D0
        DistLs=0.D0
        LambdaMag=0.D0
        SpatOrbs=nBasis/2
        IF(tSpinOrbs) THEN
            NoOrbs=nBasis
            NoOcc=NEl
        ELSE
            NoOrbs=SpatOrbs
            NoOcc=NEl/2
        ENDIF
        NoRotOrbs=NoOrbs
        Iteration=0
        OrthoForce=0.D0
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
        CoeffUncorT2(:,:)=0.D0
         
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
        IF(tERLocalization.and.(.not.tSpinOrbs)) THEN
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
            Lambdas(:,:)=0.D0
            DerivLambda(:,:)=0.D0
        ENDIF

        IF(tShake) THEN
            ALLOCATE(ShakeLambda(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('ShakeLambda',TotNoConstraints,8,this_routine,ShakeLambdaTag,ierr)
            ShakeLambda(:)=0.D0                     
            ALLOCATE(ShakeLambdaNew(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('ShakeLambdaNew',TotNoConstraints,8,this_routine,ShakeLambdaNewTag,ierr)
            ShakeLambdaNew(:)=0.D0                     
            ALLOCATE(Constraint(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('Constraint',TotNoConstraints,8,this_routine,ConstraintTag,ierr)
            ALLOCATE(ConstraintCor(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('ConstraintCor',TotNoConstraints,8,this_routine,ConstraintCorTag,ierr)
            ALLOCATE(DerivConstrT1(NoOrbs,NoOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT1',NoOrbs*TotNoConstraints*NoOrbs,8,this_routine,DerivConstrT1Tag,ierr)
            DerivConstrT1(:,:,:)=0.D0
            ALLOCATE(DerivConstrT2(NoOrbs,NoOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT2',NoOrbs*TotNoConstraints*NoOrbs,8,this_routine,DerivConstrT2Tag,ierr)
            ALLOCATE(ForceCorrect(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('ForceCorrect',NoOrbs**2,8,this_routine,ForceCorrectTag,ierr)
            ALLOCATE(Correction(NoOrbs,NoOrbs),stat=ierr)
            CALL LogMemAlloc('Correction',NoOrbs**2,8,this_routine,CorrectionTag,ierr)
            IF(tShakeApprox) THEN
                ALLOCATE(DerivConstrT1T2Diag(TotNoConstraints),stat=ierr)
                CALL LogMemAlloc('DerivConstrT1T2Diag',TotNoConstraints,8,this_routine,DerivConstrT1T2DiagTag,ierr)
                DerivConstrT1T2Diag(:)=0.D0
            ELSE
                ALLOCATE(DerivConstrT1T2(TotNoConstraints,TotNoConstraints),stat=ierr)
                CALL LogMemAlloc('DerivConstrT1T2',TotNoConstraints**2,8,this_routine,DerivConstrT1T2Tag,ierr)
            ENDIF
        ENDIF               

        ! Indexing arrays.
        ALLOCATE(SymLabelList2(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList2',NoOrbs,4,this_routine,SymLabelList2Tag,ierr)
        SymLabelList2(:)=0                     
        ALLOCATE(SymLabelList3(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3',NoOrbs,4,this_routine,SymLabelList3Tag,ierr)
        SymLabelList3(:)=0                     
 
        ALLOCATE(SymLabelListInv(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelListInv',NoOrbs,4,this_routine,SymLabelListInvTag,ierr)
        SymLabelListInv(:)=0                     
   
        ALLOCATE(Lab(2,TotNoConstraints),stat=ierr)
        CALL LogMemAlloc('Lab',2*TotNoConstraints,4,this_routine,LabTag,ierr)
        Lab(:,:)=0                     


! Do any initial calculations, and set up starting values for arrays used in rotation.        
        CALL InitRotCalc()


! Write out the headings for the results file.        
        OPEN(12,FILE='Transform',STATUS='unknown')
        IF(tLagrange) THEN
            WRITE(12,"(A12,11A18)") "# Iteration","2.PotEnergy","3.PEInts","4.PEOrtho","5.Force","6.ForceInts","7.OrthoForce","8.Sum<ij|kl>^2",&
                                    &"9.OrthoNormCondition","10.DistMovedbyCs","11.DistMovedByLs","12.LambdaMag"
            WRITE(6,"(A12,11A19)") "Iteration","2.PotEnergy","3.PEInts","4.PEOrtho","5.Force","6.ForceInts","7.OrthoForce","8.Sum<ij|kl>^2",&
                                    &"9.OrthoNormCondition","10.DistMovedbyCs","11.DistMovedbyLs","12.LambdaMag"
        ELSEIF(tERLocalization.and.tHijSqrdMin) THEN
            WRITE(12,"(A12,7A24)") "# Iteration","2.ERPotEnergy","3.HijSqrdPotEnergy","4.PotEnergy","5.Force","6.Totalcorrforce","7.OrthoNormCondition","8.DistMovedbyCs"
            WRITE(6,"(A12,7A24)") "# Iteration","2.ERPotEnergy","3.HijSqrdPotEnergy","4.PotEnergy","5.Force","6.Totalcorrforce","7.OrthoNormCondition","8.DistMovedbyCs"
        ELSEIF(tERLocalization) THEN
            WRITE(12,"(A12,5A24)") "# Iteration","2.Sum_i<ii|ii>","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
            WRITE(6,"(A12,5A24)") "Iteration","2.Sum_i<ii|ii>","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        ELSE
            WRITE(12,"(A12,5A24)") "# Iteration","2.PotEnergy","3.Force","4.Totalcorrforce","5.OrthoNormCondition","6.DistMovedbyCs"
            WRITE(6,"(A12,5A24)") "Iteration","2.PotEnergy","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        ENDIF


    END SUBROUTINE InitLocalOrbs


    SUBROUTINE InitRotCalc()
! Sets up the initial arrays to be used in the orbital rotation.    
        CHARACTER(len=*) , PARAMETER :: this_routine='InitRotCalc'
        REAL*8 :: t,RAN2,s,EpsilonHOMO,EpsilonLUMO,EpsilonMidVirt
        INTEGER :: x,y,i,j,k,l,ierr,Const,iseed=-8,a,b,g,d,MinRot,MaxRot,MinA
        INTEGER :: SumNoOrbs,NoOrbsPerProc,NoOrbsRem,m,w,z,SymM,SymMin


        CALL InitSymmArrays()
! Creates an indexing system for each of the cases with symmetry on/off, and mixing all orbitals or separating
! the occuppied from virtual.
! The arrays used in this routine are labelled with a 2 (SymLabelList2 and SymLabelCount2), so as to not
! mess up the spawing/FCI calcs.
        do i=1,NoOrbs
            SymLabelList3(i)=SymLabelList2(i)
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
! The ordering of the orbitals in CoeffT1 follow the ordering in SymLabelList2.

        CoeffT1(:,:)=0.D0
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
            CoeffT1(i,i)=1.D0
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
!                    SymM=INT(G1(SymLabelList2(m)*2)%sym%S,4)
!                    do z=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
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
! i.e. in UMatTemp(1,1,1,1), the orbital involved is the first in SymLabelList2.
! Doing this now, rather than using UMatInd in each transform2elint routine proved a lot faster.
        DerivCoeff(:,:)=0.D0
        UMATTemp01(:,:,:,:)=0.D0
        IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
        &.or.(tERLocalization.and.tSpinOrbs)) UMATTemp02(:,:,:,:)=0.D0

        CALL CopyAcrossUMAT()

        CALL TestOrthonormality()
    

!        WRITE(6,*) 'i,j,TMAT2D'
!        do i=1,NoOrbs*2
!            do j=1,NoOrbs*2
!                WRITE(6,'(F20.10)',advance='no') REAL(TMAT2D(i,j)%v,8)
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

        IF(tERLocalization.and.(.not.tSpinOrbs)) THEN
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
        INTEGER :: a,b,g,d,i,j,k,l,i2,j2,k2,l2,Spina,Spinb,Sping,Spind
        REAL*8 :: s,t


        IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
        &.or.(tERLocalization.and.tSpinOrbs)) TMAT2DTemp(:,:)=0.D0


! These loops can be sped up with spatial symmetry and pairwise permutation symmetry if needed.
        do a=1,NoOrbs
            i=SymLabelList2(a)                  ! The spin orbital we are looking for.
!            IF(tSpinOrbs) THEN
!                i2=CEILING(i/2.0)     ! The spatial orbital we are looking for.
!                IF(MOD(i,2).eq.0) THEN
!                    Spina=2                         ! Alpha orbital has spin down.
!                ELSE
!                    Spina=1                         ! Alpha orbital has spin up.
!                ENDIF
!            ENDIF
            do g=1,a
                j=SymLabelList2(g)
!                IF(tSpinOrbs) THEN
!                    j2=CEILING(j/2.0)
!                    IF(MOD(j,2).eq.0) THEN
!                        Sping=2
!                    ELSE
!                        Sping=1
!                    ENDIF
!                ENDIF
                IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
                &.or.(tERLocalization.and.tSpinOrbs)) THEN
                    IF(tSpinOrbs) THEN
                        s=REAL(TMAT2D(i,j)%v,8)
                        TMAT2DTemp(a,g)=s
                        TMAT2DTemp(g,a)=s
                    ELSE
                        s=REAL(TMAT2D(2*i,2*j)%v,8)
                        TMAT2DTemp(a,g)=s
                        TMAT2DTemp(g,a)=s
                    ENDIF
                ENDIF


                do b=1,NoOrbs
                    k=SymLabelList2(b)
!                    IF(tSpinOrbs) THEN
!                        IF(b.eq.a) CYCLE    ! alpha and beta spin orbitals the same - <ab|gd> = 0
!                        k2=CEILING(k/2.0)
!                        IF(MOD(k,2).eq.0) THEN
!                            Spinb=2
!                        ELSE
!                            Spinb=1
!                        ENDIF
!                    ENDIF
                    do d=1,b
                        l=SymLabelList2(d)
!                        IF(tSpinOrbs) THEN
!                            IF(d.eq.g) CYCLE        ! gamma and delta spin orbitals the same.
!                            l2=CEILING(l/2.0)
!                            IF(MOD(l,2).eq.0) THEN
!                                Spind=2
!                            ELSE
!                                Spind=1
!                            ENDIF

!                            t=REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0))%v,8)
!                            IF((Spina.eq.Sping).and.(Spinb.eq.Spind)) THEN
!                                UMATTemp01(a,g,b,d)=t
!                                UMATTemp01(g,a,b,d)=t
!                                UMATTemp01(a,g,d,b)=t
!                                UMATTemp01(g,a,d,b)=t
!                                IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat))&
!                                &.or.(tERLocalization.and.tSpinOrbs)) THEN
!                                    UMATTemp02(b,d,a,g)=t
!                                    UMATTemp02(b,d,g,a)=t
!                                    UMATTemp02(d,b,a,g)=t
!                                    UMATTemp02(d,b,g,a)=t
!                                ENDIF
!                            ENDIF
!                        ELSE
                            t=REAL(UMAT(UMatInd(i,k,j,l,0,0))%v,8)
                            UMATTemp01(a,g,b,d)=t                   !a,g,d,b chosen to make 'transform2elint' steps more efficient
                            UMATTemp01(g,a,b,d)=t
                            UMATTemp01(a,g,d,b)=t
                            UMATTemp01(g,a,d,b)=t
                            IF(((.not.tERLocalization).and.(.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs))&
                            &.or.(tERLocalization.and.tSpinOrbs)) THEN
                                UMATTemp02(d,b,a,g)=t                   !d,b,a,g order also chosen to speed up the transformation.
                                UMATTemp02(d,b,g,a)=t
                                UMATTemp02(b,d,a,g)=t
                                UMATTemp02(b,d,g,a)=t
                            ENDIF
!                        ENDIF
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
            WRITE(6,"(I12,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
            WRITE(12,"(I12,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
        ELSEIF(tERLocalization.and.tHijSqrdMin) THEN
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(6,"(I12,7F24.10)") Iteration,ERPotEnergy,HijSqrdPotEnergy,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
                WRITE(12,"(I12,7F24.10)") Iteration,ERPotEnergy,HijSqrdPotEnergy,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
            ENDIF
        ELSE
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(6,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
                WRITE(12,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
            ENDIF
        ENDIF
        CALL FLUSH(6)
        CALL FLUSH(12)

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
! SymLabelList2(NoOrbs) contains the spatial orbitals, ordered in groups of increasing symmetry label.
! - when the orbitals are being separated, the first NoOcc of SymLabelList2 are the occupied, and the rest are virtual.
! - essentially this array relates the orbital labelling used in the orbital rotation (1,2,3 according to the order 
! - in SymLabelList2) to the labels used in arrays being fed in/out of this routine (UMAT etc).

! SymLabelCounts2(1:Sym) is the index in SymLabelList where the symmetry block S starts
! SymLabelCounts2(2:Sym) is the number of orbitals in symmetry block S.
! E.g. if symmetry S starts at index 2 and has 3 orbitals.
! SymLabelList2(2)->SymLabelList2(4) will give the indexes of these orbitals.
        INTEGER :: j,i,ierr,SymSum
        CHARACTER(len=*) , PARAMETER :: this_routine='InitSymmArrays'

        SymLabelCounts(:,:)=0
        SymLabelList(:)=0

        IF(.not.tSeparateOccVirt) THEN
            IF(tSpinOrbs) CALL Stop_All(this_routine,"There may be a problem with GENSymStatePairs when using spin orbitals.")
            CALL GENSymStatePairs(SpatOrbs,.false.)
        ENDIF
! Sets up the SymLabelList and SymLabelCounts arrays used in the spawing etc. (When the rotate
! orbs routine is called, this has not been done yet).
! If the symmetry is on, and all orbitals are being mixed, this will end up being the same as SymLabelList2.


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
            ALLOCATE(SymLabelCounts2(2,8),stat=ierr)
            CALL LogMemAlloc('SymLabelCounts2',2*8,4,this_routine,SymLabelCounts2Tag,ierr)
            SymLabelCounts2(:,:)=0                     
            do i=1,SpatOrbs   
                IF(tSpinOrbs) THEN
                    SymLabelList2(2*i)=2*SymLabelList(i)
                    SymLabelList2(2*i-1)=(2*SymLabelList(i))-1
                ELSE
                    SymLabelList2(i)=SymLabelList(i)
                ENDIF
            enddo
            IF(lNoSymmetry) THEN
                SymLabelCounts2(1,1)=1
                SymLabelCounts2(2,1)=NoOrbs
            ELSE
                do j=1,8
                    IF(tSpinOrbs) THEN
                        SymLabelCounts2(1,j)=(2*SymLabelCounts(1,j))-1
                        SymLabelCounts2(2,j)=2*SymLabelCounts(2,j)
                    ELSE
                        do i=1,2
                            SymLabelCounts2(i,j)=SymLabelCounts(i,j)
                        enddo
                    ENDIF
                enddo
            ENDIF
        ENDIF

        do i=1,NoOrbs
            SymLabelListInv(SymLabelList2(i))=i
        enddo
        

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,16
!            WRITE(6,*) SymLabelCounts2(1,i),SymLabelCounts2(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
!        do i=1,NoOrbs
!            IF(tSpinOrbs) THEN
!                WRITE(6,*) i,SymLabelList2(i),INT(G1(SymLabelList2(i))%sym%S,4)
!            ELSE
!                WRITE(6,*) i,SymLabelList2(i),INT(G1(SymLabelList2(i)*2)%sym%S,4)
!            ENDIF
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,NoOrbs
!            WRITE(6,*) SymLabelList2(i),SymLabelListInv(i)
!        enddo
!        CALL FLUSH(6)
!        stop


    ENDSUBROUTINE InitSymmArrays



    SUBROUTINE EquateDiagFock()
        INTEGER :: irr,NumInSym,Orbi,Orbj,w,i,j,Prod,ConjOrb,k,ConjInd,OrbjConj
        REAL*8 :: Angle,AngleConj,Check,Norm
        REAL*8 , PARAMETER :: PI=3.14159265358979323846264338327950288419716939937510D0

        CoeffT1(:,:)=0.D0
!        MaxOccVirt=1
!        WRITE(6,*) MaxOccVirt,"***"

        do w=MinOccVirt,MaxOccVirt
!Do virtual and occupied orbitals seperately

            do irr=1,8
!Loop over irreps

                NumInSym=SymLabelCounts2(2,(w-1)*8+irr)
!                WRITE(6,*) "NumInSym= ",NumInSym,irr-1

                do j=1,NumInSym
!Loop over the j-orthogonal vectors to create in this symmetry block

                    Orbj=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+j)

!See if this vector has already been done.
                    Check=0.D0
                    do i=1,NoOrbs
                        Check=Check+CoeffT1(i,Orbj)
                    enddo
                    IF(Check.ne.0.D0) THEN
!This vector is a conjugate pair of another vector and has already been worked out...
                        CYCLE
                    ENDIF

!Find out if we this vector will be complex. It will be real if j=N or j=N/2
                    IF(j.eq.NumInSym) THEN
!The vector will be the normalized 1,1,1 vector.

                        do i=1,NumInSym
                            Orbi=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+i)
                            CoeffT1(Orbi,Orbj)=1/SQRT(REAL(NumInSym,8))
                        enddo

                    ELSEIF((mod(NumInSym,2).eq.0).and.(j.eq.(NumInSym/2))) THEN

                        do i=1,NumInSym
                            Orbi=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+i)
                            IF(mod(i,2).eq.1) THEN
                                CoeffT1(Orbi,Orbj)=-1/SQRT(REAL(NumInSym,8))
                            ELSE
                                CoeffT1(Orbi,Orbj)=1/SQRT(REAL(NumInSym,8))
                            ENDIF
                        enddo

                    ELSE
!Vector is complex - find its conjugate vector - do these at the same time.
                        ConjInd=NumInSym-j
                        OrbjConj=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+ConjInd)
                        
                        do i=1,NumInSym
                            
                            Orbi=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+i)

                            Angle=REAL(i*j*2,8)*PI/REAL(NumInSym,8)
                            AngleConj=REAL(i*ConjInd*2,8)*PI/REAL(NumInSym,8)

                            CoeffT1(Orbi,Orbj)=(1/SQRT(REAL(2*NumInSym,8)))*(COS(Angle)+COS(AngleConj))
                            CoeffT1(Orbi,OrbjConj)=(1/SQRT(REAL(2*NumInSym,8)))*(SIN(Angle)-SIN(AngleConj))

                        enddo

                    ENDIF

                enddo


!                    do i=1,NumInSym
!
!                        Orbi=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+i)
!                        WRITE(6,*) "Sym= ",irr-1, Orbj, Orbi
!!Coefficients are going to be C_jk = exp^(i j*k 2Pi/N), i.e roots of unity
!
!!                        IF(CoeffT1(i,j).ne.0.D0) CYCLE
!
!                        Prod=i*j
!                        IF(mod(Prod,NumInSym).eq.0) THEN
!! i*j = N or 2N, 3N, ...
!                            CoeffT1(Orbi,Orbj)=1.D0/SQRT(REAL(NumInSym,8))
!                            CoeffT1(Orbj,Orbi)=1.D0/SQRT(REAL(NumInSym,8))
!
!                        ELSEIF((mod(Prod*2,NumInSym).eq.0).and.(mod((2*Prod)/NumInSym,2).eq.1)) THEN
!! i*j = N/2 or 3N/2, 5N/2, ...
!                            CoeffT1(Orbi,Orbj)=-1.D0/SQRT(REAL(NumInSym,8))
!                            CoeffT1(Orbj,Orbi)=-1.D0/SQRT(REAL(NumInSym,8))
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
!                                ConjOrb=SymLabelList2(SymLabelCounts2(1,(w-1)*8+irr)-1+ConjInd)
!
!                                Angle=REAL(Prod*2,8)*3.141592654/REAL(NumInSym,8)
!                                CoeffT1(Orbi,Orbj)=1.D0/SQRT(REAL(2*NumInSym,8))*2.D0*COS(Angle)
!                                CoeffT1(ConjOrb,Orbj)=1.D0/SQRT(REAL(2*NumInSym,8))*2.D0*SIN(Angle)
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
            Norm=0.D0
            do i=1,NoOrbs
                Norm=Norm+(CoeffT1(i,j)**2)
            enddo
            IF(Norm.eq.0.D0) THEN
                CoeffT1(j,j)=1.D0
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
            Norm=0.D0
            do i=1,NoOrbs
                Norm=Norm+(CoeffT1(i,j)**2)
            enddo
            IF(abs(Norm-1.D0).gt.1.D-7) THEN
                CALL Stop_All("EquateDiagFock","Rotation Coefficients not normalized")
            ENDIF
        enddo


!Check orthogonality
        do j=1,NoOrbs
            do i=1,NoOrbs
                IF(i.eq.j) CYCLE
                Norm=0.D0
                do k=1,NoOrbs
                    Norm=Norm+(CoeffT1(k,j)*CoeffT1(k,i))
                enddo
                IF(abs(Norm).gt.1.D-7) THEN
                    WRITE(6,*) "COLUMNS: ",j,i
                    CALL Stop_All("EquateDiagFock","RotationCoefficients not orthogonal")
                ENDIF
            enddo
        enddo

    END SUBROUTINE EquateDiagFock


    SUBROUTINE InitOrbitalSeparation()
! This subroutine is called if the SEPARATEOCCVIRT keyword is present in the input, it sets up SymLabelList2 so that the first 
! NoOcc orbitals are the HF occupied, and the rest the virtual.  Within this separation, orbitals are ordered in symmetry 
! groups. 
! This means that two iterations of the rotate orbs routine will be performed, the first treats the occupied orbitals and the second
! the virtual.
        INTEGER :: i,j,ierr,t,SymCurr,Symi,SymVirtOrbsTag,SymOccOrbsTag,Temp
        INTEGER , ALLOCATABLE :: SymVirtOrbs(:),SymOccOrbs(:)
        CHARACTER(len=*) , PARAMETER :: this_routine='InitOrbitalSeparation'


        ALLOCATE(SymLabelCounts2(2,16),stat=ierr)
        CALL LogMemAlloc('SymLabelCounts2',2*16,4,this_routine,SymLabelCounts2Tag,ierr)
        SymLabelCounts2(:,:)=0
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


! First fill SymLabelList2.

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index with the second lowest energy.

!        do i=1,nBasis
!            WRITE(6,*) BRR(i)
!        enddo

! this picks out the NoOcc lowest energy orbitals from BRR as these will be the occupied.
! these are then ordered according to symmetry, and the same done to the virtual.
        do i=1,NoOcc
            IF(tSpinOrbs) THEN
                LabOccOrbs(i)=BRR(i)
                SymOccOrbs(i)=INT(G1(LabOccOrbs(i))%sym%S,4)
            ELSE
                LabOccOrbs(i)=(BRR(2*i))/2
                SymOccOrbs(i)=INT(G1(LabOccOrbs(i)*2)%sym%S,4)
            ENDIF
        enddo
        
        CALL NECI_SORT2I(NoOcc,SymOccOrbs,LabOccOrbs)
        ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in terms of symmetry). 

        do i=1,NoOrbs-NoOcc
            IF(tSpinOrbs) THEN
                LabVirtOrbs(i)=BRR(i+NEl)
                SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i))%sym%S,4)
            ELSE
                LabVirtOrbs(i)=(BRR((2*i)+NEl))/2
                SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i)*2)%sym%S,4)
            ENDIF
        enddo
        
        CALL NECI_SORT2I((NoOrbs-NoOcc),SymVirtOrbs,LabVirtOrbs)

! SymLabelList2 is then filled with the symmetry ordered occupied then virtual arrays.        
        do i=1,NoOcc
            SymLabelList2(i)=LabOccOrbs(i)
        enddo
        j=0
        do i=NoOcc+1,NoOrbs
            j=j+1
            SymLabelList2(i)=LabVirtOrbs(j)
        enddo

!        WRITE(6,*) 'symlabellist'
!        do i=1,NoOrbs
!            WRITE(6,'(2I4)') SymLabelList2(i),INT(G1(SymLabelList2(i)*2)%sym%S,4)
!        enddo
!        stop

!************
! Second fill SymLabelCounts2.
! - the first 8 places of SymLabelCounts2(1,:) and SymLabelCounts2(2,:) refer to the occupied orbitals 
! - and the second 8 to the virtuals.

        IF(lNoSymmetry) THEN
            ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
            SymLabelCounts2(1,1)=1
            SymLabelCounts2(1,9)=NoOcc+1
            SymLabelCounts2(2,1)=NoOcc
            SymLabelCounts2(2,9)=NoOrbs-NoOcc
        ELSE 
            ! otherwise we run through the occupied orbitals, counting the number with each symmetry
            ! and noting where in SymLabelList2 each symmetry block starts.
            SymCurr=0
            SymLabelCounts2(1,1)=1
            do i=1,NoOcc
                IF(tSpinOrbs) THEN
                    Symi=INT(G1(SymLabelList2(i))%sym%S,4)
                ELSE
                    Symi=INT(G1(SymLabelList2(i)*2)%sym%S,4)
                ENDIF
                SymLabelCounts2(2,(Symi+1))=SymLabelCounts2(2,(Symi+1))+1
                IF(Symi.gt.SymCurr) THEN
                    SymLabelCounts2(1,(Symi+1))=i
                    SymCurr=Symi
                ENDIF
            enddo
            ! the same is then done for the virtuals.
            SymCurr=0
            SymLabelCounts2(1,9)=NoOcc+1
            do i=NoOcc+1,NoOrbs
                IF(tSpinOrbs) THEN
                    Symi=INT(G1(SymLabelList2(i))%sym%S,4)
                ELSE
                    Symi=INT(G1(SymLabelList2(i)*2)%sym%S,4)
                ENDIF
                SymLabelCounts2(2,(Symi+9))=SymLabelCounts2(2,(Symi+9))+1
                IF(Symi.gt.SymCurr) THEN
                    SymLabelCounts2(1,(Symi+9))=i
                    SymCurr=Symi
                ENDIF
            enddo
        ENDIF

        ! Go through each symmetry group, making sure the orbital pairs are ordered lowest to highest.
        do i=1,16
            IF(SymLabelCounts2(2,i).ne.0) THEN
                CALL NECI_SORTI(SymLabelCounts2(2,i),SymLabelList2(SymLabelCounts2(1,i):(SymLabelCounts(1,i)+SymLabelCounts2(2,i)-1)))
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
        INTEGER :: i,j,Sym,ierr,NoSymBlock,TMAT2DSymBlockTag,WorkSize,WorkCheck,WorkTag,DiagTMAT2DBlockTag,SymStartInd
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
        TMAT2DRot(:,:)=0.D0
        DiagTMAT2Dfull(:)=0.D0

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

            NoSymBlock=SymLabelCounts2(2,Sym+9)

            SymStartInd=SymLabelCounts2(1,Sym+9)-1
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
                    CALL FLUSH(6)
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
        REAL*8 :: Coeff(NoOrbs,NoOrbs)
        INTEGER :: i,j

        do i=1,NoOcc
            do j=NoOcc+1,NoOrbs
                Coeff(i,j)=0.D0
                Coeff(j,i)=0.D0
            enddo
        enddo


    ENDSUBROUTINE ZeroOccVirtElements




    SUBROUTINE FindNewOrbs()
           
        IF(tERLocalization.and.(.not.tSpinOrbs)) THEN
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
        REAL*8 :: t,Temp4indints(NoRotOrbs,NoOrbs)
        REAL*8 :: Temp4indints02(NoRotOrbs,NoRotOrbs)  

        
        CALL set_timer(Transform2ElInts_time,30)

!Zero arrays from previous transform

        TwoIndInts01(:,:,:,:)=0.D0
        FourIndInts(:,:,:,:)=0.D0

        IF(tNotConverged) THEN
!            TwoIndInts02Temp(:,:,:,:)=0.D0
            TwoIndInts02(:,:,:,:)=0.D0
!            ThreeIndInts01Temp(:,:,:,:)=0.D0
            ThreeIndInts01(:,:,:,:)=0.D0
!            ThreeIndInts02Temp(:,:,:,:)=0.D0
            ThreeIndInts02(:,:,:,:)=0.D0
!            ThreeIndInts03Temp(:,:,:,:)=0.D0
            ThreeIndInts03(:,:,:,:)=0.D0
!            ThreeIndInts04Temp(:,:,:,:)=0.D0
            ThreeIndInts04(:,:,:,:)=0.D0
!            FourIndInts02Temp(:,:,:,:)=0.D0
            FourIndInts02(:,:,:,:)=0.D0
        ENDIF

! ************
!Transform the 1 electron, 2 index integrals (<i|h|j>).
        IF(tNotConverged) THEN
            TMAT2DRot(:,:)=0.D0
            TMAT2DPartRot01(:,:)=0.D0
            TMAT2DPartRot02(:,:)=0.D0

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


            CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TMAT2DTemp(:,:),NoOrbs,0.0,TMAT2DPartRot01(:,:),NoOrbs)
            ! get TMAT2DPartRot(i,a) out of this.

            CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TMAT2DTemp(:,:),NoOrbs,0.0,TMAT2DPartRot02(:,:),NoOrbs)
            ! get TMAT2DPartRot(a,j) out of this.
     
            CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TMAT2DPartRot01(:,:),NoOrbs,0.0,TMAT2DRot(:,:),NoOrbs)
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
                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',NoRotOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMatTemp01(:,:,d,b),NoOrbs,0.0,Temp4indints(:,:),NoRotOrbs)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, we need the transpose of this.

                Temp4indints02(:,:)=0.D0
                CALL DGEMM('T','T',NoRotOrbs,NoRotOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(:,:),NoRotOrbs,0.0,Temp4indints02(:,:),NoRotOrbs)
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
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMatTemp02(:,:,a,g),NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
                    ! Temp4indints(l,b) comes out of here, so need to use transpose of this to transform the b elements.

                    Temp4indints02(:,:)=0.D0
                    CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(:,:),NoOrbs,0.0,Temp4indints02(:,:),NoOrbs)
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
                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',NoRotOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts01(:,:,k,i),NoOrbs,0.0,Temp4indints(:,:),NoRotOrbs)

                IF(tNotConverged) THEN
                    do b=1,NoOrbs
                        do l=1,NoOrbs
                            ThreeIndInts02(i,k,l,b)=Temp4indints(l,b)
                            ThreeIndInts02(k,i,l,b)=Temp4indints(l,b)
                        enddo
                    enddo
                    Temp4indints02(:,:)=0.D0
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts01(:,:,k,i),NoOrbs,0.0,Temp4indints02(:,:),NoRotOrbs)
                    do d=1,NoOrbs
                        do j=1,NoOrbs
                            ThreeIndInts04(k,i,j,d)=Temp4indints02(j,d)
                            ThreeIndInts04(i,k,j,d)=Temp4indints02(j,d)
                        enddo
                    enddo
                ENDIF
                Temp4indints02(:,:)=0.D0
                CALL DGEMM('T','T',NoRotOrbs,NoRotOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,Temp4indints(:,:),NoRotOrbs,0.0,Temp4indints02(:,:),NoRotOrbs)
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
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts02(:,:,j,l),NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
                    do a=1,NoOrbs
                        do k=1,NoOrbs
                            ThreeIndInts01(k,j,l,a)=Temp4indints(k,a)
                            ThreeIndInts01(k,l,j,a)=Temp4indints(k,a)
                        enddo
                    enddo
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,TwoIndInts02(:,:,j,l),NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
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
    
        IF((.not.tReadInCoeff).and.(.not.tUseMP2VarDenMat).and.(.not.tFindCINatOrbs)) THEN
            
            PotEnergy=0.D0
            TwoEInts=0.D0
            PEInts=0.D0
            CALL CalcPotentials()

            IF(tPrintInts) CALL PrintIntegrals()
            IF((Iteration.eq.0).or.((.not.tNotConverged).and.(Iteration.gt.1))) CALL WriteDoubHisttofile()
            IF(tROHistSingExc.and.(Iteration.eq.0)) CALL WriteSingHisttofile()


! If doing Lagrange orthormalisations, find the change of the potential energy due to the orthonormality 
! of the orbitals...
            IF(tLagrange) THEN
                PEOrtho=0.D0
                do i=1,NoOrbs
                    do j=1,NoOrbs
                        t=0.D0
                        do a=1,NoOrbs
                            t=CoeffT1(a,i)*CoeffT1(a,j)
                        enddo
                        IF(i.eq.j) t=t-1.D0
                        PEOrtho=PEOrtho-Lambdas(i,j)*t
                        PotEnergy=PotEnergy-Lambdas(i,j)*t
                    enddo
                enddo
            ENDIF
        ENDIF

        CALL halt_timer(Transform2ElInts_Time)


    END SUBROUTINE Transform2ElInts

   
! This is a transformation of the four index integrals for the ERlocalisation, in this only the <ii|ii> integrals are needed 
! therefore the process may be much simpler.
    SUBROUTINE Transform2ElIntsERlocal()
        INTEGER :: i,j,k,l,a,b,g,d,m
        REAL*8 :: t,Temp4indints(NoOrbs,NoOrbs)
        REAL*8 :: Temp4indints02(NoOrbs)  
 

        CALL set_timer(Transform2ElInts_time,30)

! Zero arrays from previous transform

        TwoIndIntsER(:,:,:)=0.D0

        ThreeIndInts01ER(:,:)=0.D0
        ThreeIndInts02ER(:,:)=0.D0

        FourIndIntsER(:)=0.D0

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
                Temp4indints(:,:)=0.D0
                Temp4indints02(:)=0.D0
                CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMATTemp01(:,:,b,d),NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
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


        TwoIndIntsER(:,:,:)=0.D0
!        do a=LowBound02,HighBound02
        do d=1,NoOrbs
            do b=1,NoOrbs
                Temp4indints(:,:)=0.D0
                Temp4indints02(:)=0.D0
                CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,CoeffT1(:,:),NoOrbs,UMATTemp01(:,:,b,d),NoOrbs,0.0,Temp4indints(:,:),NoOrbs)
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
        
        PotEnergy=0.D0
        TwoEInts=0.D0
        PEInts=0.D0
        CALL CalcPotentials()

        IF(tPrintInts) CALL PrintIntegrals()
        IF((Iteration.eq.0).or.((.not.tNotConverged).and.(Iteration.gt.1))) CALL WriteDoubHisttofile()
        IF(tROHistSingExc.and.(Iteration.eq.0)) CALL WriteSingHisttofile()


! If doing Lagrange orthormalisations, find the change of the potential energy due to the orthonormality 
! of the orbitals...
        IF(tLagrange) THEN
            PEOrtho=0.D0
            do i=1,NoOrbs
                do j=1,NoOrbs
                    t=0.D0
                    do a=1,NoOrbs
                        t=CoeffT1(a,i)*CoeffT1(a,j)
                    enddo
                    IF(i.eq.j) t=t-1.D0
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
        REAL*8 :: MaxTerm


        IF(tERLocalization.and.(.not.tSpinOrbs)) THEN
            ERPotEnergy=0.D0
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
            CoulPotEnergy=0.D0
            OffDiagPotEnergy=0.D0
!            do i=1,NoOrbs
!                IF((i.ge.Starti).and.(i.le.Finishi)) THEN
            do i=Starti,Finishi
                    ERPotEnergy=ERPotEnergy+FourIndIntsER(i)
                    IF(FourIndIntsER(i).lt.0) THEN
                        CALL FLUSH(6)
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
            ERPotEnergy=0.D0
            PotEnergy=0.D0
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
                IF(tSpinOrbs) THEN
                    IF(MOD(i,2).eq.0) THEN
                        j=i-1
                    ELSE
                        j=i+1
                    ENDIF
                    ERPotEnergy=ERPotEnergy+FourIndInts(i,j,i,j)
                    IF((FourIndInts(i,j,i,j).lt.0).or.(FourIndInts(j,i,j,i).lt.0)) THEN
                        CALL FLUSH(6)
                        CALL Stop_All('CalcPotentials','A <ii|ii> value is less than 0.')
                    ENDIF
                    PotEnergy=PotEnergy+FourIndInts(i,j,i,j)
                    TwoEInts=TwoEInts+FourIndInts(i,j,i,j)
                    PEInts=PEInts+FourIndInts(i,j,i,j)
                ELSE
                    ERPotEnergy=ERPotEnergy+FourIndInts(i,i,i,i)
                    IF((FourIndInts(i,i,i,i).lt.0)) THEN
                        CALL FLUSH(6)
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
                MaxTerm=0.D0
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
            HijSqrdPotEnergy=0.D0
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
            ERPotEnergy=0.D0
            ijOccVirtPotEnergy=0.D0
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
        INTEGER :: m,z,i,j,k,l,a,b,g,d,Symm,Symb,Symb2,Symi,Symd,w,x,y,SymMin,n,o,p,q
        REAL*8 :: OffDiagForcemz,DiagForcemz,OneElForcemz,t1,t2,t3,t4,LambdaTerm1,LambdaTerm2,t1Temp,ForceTemp
        REAL*8 :: NonDerivTerm,DerivPot
        CHARACTER(len=*) , PARAMETER :: this_routine='FindtheForce'
        LOGICAL :: leqm,jeqm,keqm,ieqm
      
! Running over m and z, covers all matrix elements of the force matrix (derivative 
! of equation we are minimising, with respect to each translation coefficient) filling 
! them in as it goes.
        CALL set_timer(FindtheForce_time,30)
        
        DerivCoeff(:,:)=0.D0
        Force=0.D0
!        ForceTemp=0.D0
        ForceInts=0.D0
        OrthoForce=0.D0

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
                IF(tSpinOrbs) THEN
                    SymM=INT(G1(SymLabelList2(m))%sym%S,4)
                ELSE
                    SymM=INT(G1(SymLabelList2(m)*2)%sym%S,4)
                ENDIF
                do z=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
 
                    ! Find the force on a coefficient c(m,z). 
                    OffDiagForcemz=0.D0
                    ! OffDiagForce is from any of the OffDiagMin/Max (Sqrd or not), or the double/single excitation
                    ! max/min, as only one of these terms may be used at once.

                    DiagForcemz=0.D0
                    ! DiagForce includes ER localisation, and the coulomb terms <ij|ij>.

                    OneElForcemz=0.D0
                    ! OneElForce includes that from the one electron integrals <i|h|j> and the one particle orbital 
                    ! energies.

                    ! DIAG TERMS
                    ! Maximise <ii|ii>, self interaction terms. 
                    IF(tERLocalization.and.(.not.tSpinOrbs)) THEN
!                        DiagForcemz=DiagForcemz+ThreeIndInts01(m,m,m,z)+ThreeIndInts02(m,m,m,z)+ThreeIndInts03(m,m,m,z)+ThreeIndInts04(m,m,m,z)
                        DiagForcemz=DiagForcemz+(2*ThreeIndInts01ER(z,m))+(2*ThreeIndInts02ER(z,m))
                        ! Derivative of <ii|ii> only non-zero when i=m.
                        ! each of the four terms then correspond to zeta = a, b, g, then d in the unrotated basis.
                    ELSEIF(tERLocalization) THEN
!                    IF(tERLocalization) THEN
                        ! Looking at <ij|ij> terms where j=i+1 (i.e. i is alpha of spin orbital and j is beta - or vice versa).
                        IF(tSpinOrbs) THEN
                            IF(MOD(m,2).eq.0) THEN      ! m = j
                                i=m-1
                                DiagForcemz=DiagForcemz+ThreeIndInts01(m,i,i,z)+ThreeIndInts01(z,i,i,m)+ThreeIndInts01(i,m,z,i)+ThreeIndInts01(i,z,m,i)
                            ELSE
                                j=m+1
                                DiagForcemz=DiagForcemz+ThreeIndInts01(m,j,j,z)+ThreeIndInts01(z,j,j,m)+ThreeIndInts01(j,m,z,j)+ThreeIndInts01(j,z,m,j)
                            ENDIF
                        ELSE
                            DiagForcemz=DiagForcemz+ThreeIndInts01(m,m,m,z)+ThreeIndInts02(m,m,m,z)+ThreeIndInts03(m,m,m,z)+ThreeIndInts04(m,m,m,z)
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
                            DerivPot=0.D0
                            DerivPot=DerivPot+TMAT2DPartRot02(z,m)+TMAT2DPartRot01(m,z)
                            ! First term when m=i and z=a, second when m=i and z=b.
                            ! This is all that is needed for OneElIntMax

                            IF(tOnePartOrbEnMax) THEN
                                NonDerivTerm=0.D0
                                IF(OrbEnMaxAlpha.ne.1.D0) THEN 
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
                                        DerivPot=DerivPot+(2*ThreeIndInts01(m,j,j,z))-ThreeIndInts01(j,j,m,z)+(2*ThreeIndInts03(m,j,j,z))-ThreeIndInts03(m,j,z,j)
                                        ! First part is for when m=i and z=a, the second is for when m=i and z=g
                                    enddo
                                ENDIF
                                ! When m=j, for a particular i.
                                ! m and z run only over virtual, and j is over occupied. m will never = j.
                            ELSE
                                NonDerivTerm=1.D0
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
!                                    Symi=IEOR(INT(G1(SymLabelList2(k)*2)%sym%S,4),IEOR(INT(G1(SymLabelList2(j)*2)%sym%S,4),INT(G1(SymLabelList2(l)*2)%sym%S,4)))
                                    ! only i with symmetry equal to j x k x l will have integrals with overall
                                    ! symmetry A1 and therefore be non-zero.


                                    ! Running across i, ThreeIndInts01 only contributes if i.eq.m (which will happen once for each m)
                                    IF((m.le.k-1).and.(m.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
!                                    IF((m.le.k-1)) THEN
                                        IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+2*(FourIndInts02(j,k,l,m)*ThreeIndInts01(k,j,l,z))
                                        IF(tOffDiagMin.or.tOffDiagMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts01(k,j,l,z)
                                        IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+ThreeIndInts01(k,j,l,z)-ThreeIndInts01(l,j,k,z)
                                    ENDIF

                                    IF(jeqm) THEN
                                        do i=1,k-1
                                            IF((i.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
!                                        do i=SymLabelCounts2(1,Symi+SymMin),(SymLabelCounts2(1,Symi+SymMin)+SymLabelCounts2(2,Symi+SymMin)-1)
                                                IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+2*(FourIndInts(i,j,k,l)*ThreeIndInts02(i,k,l,z))
                                                IF(tOffDiagMin.or.tOffDiagMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts02(i,k,l,z)
                                                IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+(ThreeIndInts02(i,k,l,z))-ThreeIndInts02(i,l,k,z)
                                            ENDIF
                                        enddo
                                    ENDIF

                                    IF(keqm) THEN
!                                        do i=SymLabelCounts2(1,Symi+SymMin),(SymLabelCounts2(1,Symi+SymMin)+SymLabelCounts2(2,Symi+SymMin)-1)
                                        do i=1,k-1
                                            IF((i.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
                                                IF(tOffDiagSqrdMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+2*(FourIndInts(i,j,k,l)*ThreeIndInts03(i,j,l,z))
                                                IF(tOffDiagMin.or.tOffDiagSqrdMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts03(i,j,l,z)
                                                IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+(ThreeIndInts03(i,j,l,z))-ThreeIndInts03(i,j,z,l)
                                            ENDIF
                                        enddo
                                    ENDIF

                                    IF(leqm) THEN
!                                        do i=SymLabelCounts2(1,Symi+SymMin),(SymLabelCounts2(1,Symi+SymMin)+SymLabelCounts2(2,Symi+SymMin)-1)
                                        do i=1,k-1
                                            IF((i.ne.j).and.((i.ne.k).or.(j.ne.l))) THEN
                                                IF(tOffDiagSqrdMin.or.tOffDiagSqrdMin) OffDiagForcemz=OffDiagForcemz+2*(FourIndInts(i,j,k,l)*ThreeIndInts04(i,k,j,z))
                                                IF(tOffDiagMin.or.tOffDiagMax) OffDiagForcemz=OffDiagForcemz+ThreeIndInts04(i,k,j,z)
                                                IF(tDoubExcMin) OffDiagForcemz=OffDiagForcemz+(ThreeIndInts04(i,k,j,z))-ThreeIndInts04(i,z,j,k)
                                            ENDIF
                                        enddo
                                    ENDIF
                                enddo
                            enddo
                        enddo
                    ENDIF
                    ! DerivCoeffTemp(z,m) then combines all the different forces on coefficient(m,z).
!                    IF(tSpinOrbs) THEN
!                        WRITE(6,*) CEILING(m/2.0),CEILING(z/2.0),DiagForcemz
!                    ELSE
!                        WRITE(6,*) m,z,DiagForcemz
!                    ENDIF
                    DerivCoeff(z,m)=(MaxMinFac*OffDiagWeight*OffDiagForcemz)+(DiagMaxMinFac*DiagWeight*DiagForcemz)+(OneElMaxMinFac*OneElWeight*OneElForcemz)
                    Force=Force+ABS(DerivCoeff(z,m))
                enddo
            enddo
        enddo
!        CALL MPIDSumArr(DerivCoeffTemp(:,:),NoOrbs**2,DerivCoeff(:,:))
!        CALL MPIDSum(ForceTemp,1,Force)

        Force=Force/REAL(NoOrbs**2,8)


! Calculate the derivatives of orthogonalisation condition.
! Have taken this out of the m and z loop to make the shake faster, but can put it back in if start using it a lot.
        IF(tLagrange) THEN
            do x=MinMZ,MaxMZ
                m=SymLabelList2(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM=INT(G1(m*2)%sym%S,4)
                do y=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
                    z=SymLabelList2(y)
      
                    LambdaTerm1=0.D0
                    LambdaTerm2=0.D0
                    
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
            OrthoForce=OrthoForce/REAL(NoOrbs**2,8)
            DerivLambda(:,:)=0.D0
            do i=1,NoOrbs
                do j=1,i
                    do a=1,NoOrbs
                        DerivLambda(i,j)=DerivLambda(i,j)+CoeffT1(a,i)*CoeffT1(a,j)
                    enddo
                    DerivLambda(j,i)=DerivLambda(i,j)
                enddo
            enddo
            do i=1,NoOrbs
                DerivLambda(i,i)=DerivLambda(i,i)-1.D0
            enddo
        ENDIF


        CALL halt_timer(FindtheForce_Time)


    END SUBROUTINE FindTheForce


    
    SUBROUTINE UseTheForce()
! This routine takes the old translation coefficients and Lambdas and moves them by a timestep in the direction 
! of the calculated force.
        INTEGER :: m,w,x,y,z,i,j,Symm,SymMin
        REAL*8 :: NewCoeff,NewLambda

        DistCs=0.D0 
    
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
                IF(tSpinOrbs) THEN
                    SymM=INT(G1(SymLabelList2(m))%sym%S,4)
                ELSE
                    SymM=INT(G1(SymLabelList2(m)*2)%sym%S,4)
                ENDIF

! Symmetry requirement that z must be from the same irrep as m
                do z=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
!               
                    ! Only coeffs with sym of m and z the same have non-zero coeffs.    
                    NewCoeff=0.D0
                    NewCoeff=CoeffT1(z,m)-(TimeStep*DerivCoeff(z,m))
                    DistCs=DistCs+abs(TimeStep*DerivCoeff(z,m))
                    CoeffT1(z,m)=NewCoeff
                enddo
            enddo
        enddo
!        CALL MPIDSumArr(CoeffT1Temp(:,:),NoOrbs**2,CoeffT1(:,:))
!        CALL MPIDSum(DistCsTemp,1,DistCs)


        DistCs=DistCs/(REAL(NoOrbs**2,8))


        IF(tLagrange) THEN

            DistLs=0.D0
            LambdaMag=0.D0
            do i=1,NoOrbs
                do j=1,NoOrbs
                    NewLambda=0.D0
                    NewLambda=Lambdas(i,j)-(TimeStep*DerivLambda(i,j))  ! Timestep must be specified in the input file.
                    DistLs=DistLs+abs(TimeStep*DerivLambda(i,j))
                    Lambdas(i,j)=NewLambda
                    LambdaMag=LambdaMag+abs(NewLambda)
                enddo
            enddo
            DistLs=DistLs/(REAL(NoOrbs**2,8))
            LambdaMag=LambdaMag/(REAL(NoOrbs**2,8))

!        ELSE
!            CALL OrthoNormx(NoOrbs,NoOrbs,CoeffT1) !Explicitly orthonormalize the coefficient vectors.
        ENDIF


    ENDSUBROUTINE UseTheForce


   
    SUBROUTINE TestOrthonormality()
        INTEGER :: i,j,a
        REAL*8 :: OrthoNormDP

        OrthoNorm=0.D0
        do i=1,NoOrbs
            do j=1,i
                OrthoNormDP=0.D0
                OrthoNormDP=Dot_Product(CoeffT1(:,i),CoeffT1(:,j))
                OrthoNorm=OrthoNorm+ABS(OrthoNormDP)
            enddo
        enddo
        OrthoNorm=OrthoNorm-real(NoOrbs,8)
        OrthoNorm=(OrthoNorm*2.D0)/REAL((NoOrbs*(NoOrbs+1.D0)),8)

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
        INTEGER :: w,x,y,i,j,l,a,m,ShakeIteration,ConvergeCount,ierr,Const,SymM,SymMin
        REAL*8 :: TotCorConstraints,TotConstraints,TotLambdas
        REAL*8 :: TotUncorForce,TotDiffUncorCoeffs,TotDiffCorCoeffs
        LOGICAL :: tShakeNotConverged
        CHARACTER(len=*), PARAMETER :: this_routine='ShakeConstraints'


!        WRITE(6,*) "Beginning shakeconstraints calculation"
        IF(Iteration.eq.1) THEN
            OPEN(8,FILE='SHAKEstats',STATUS='unknown')
            WRITE(8,'(A20,4A35,A20)') 'Shake Iteration','Sum Lambdas','Total of corrected forces','Sum unconstrained constraints',&
                                        &'Sum corrected constraints','Converge count'
        ENDIF
        IF(Mod(Iteration,10).eq.0) WRITE(8,*) 'Orbital rotation iteration = ',Iteration


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


        Correction(:,:)=0.D0
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
            
            ForceCorrect(:,:)=0.D0          ! Zeroing terms that are re-calculated each iteration.
            CoeffCorT2(:,:)=0.D0
            ConstraintCor(:)=0.D0
            DerivConstrT2(:,:,:)=0.D0
            TotLambdas=0.D0
            TotCorrectedForce=0.D0
            TotDiffCorCoeffs=0.D0

            IF(ShakeIteration.ne.1) THEN
                CALL UpdateLambdas()
            ENDIF
            
            ShakeLambdaNew(:)=0.D0

            
            ! For a particular set of coefficients cm:
            ! Force(corrected)=Force(uncorrected)-Lambdas.DerivConstrT1
            ! Use these derivatives, and the current lambdas to find the trial corrected force.
            ! Then use this to get the (trial) shifted coefficients.

            ! Use the lambdas of this iteration to calculate the correction to the force due to the constraints.
            Correction(:,:)=0.D0
                

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
                    IF(tSpinOrbs) THEN
                        SymM=INT(G1(SymLabelList2(m))%sym%S,4)
                    ELSE
                        SymM=INT(G1(SymLabelList2(m)*2)%sym%S,4)
                    ENDIF
                    do a=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
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
            CALL FLUSH(6)
            CALL FLUSH(8)
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(8,'(I20,4F35.20,I20)') ShakeIteration,TotLambdas,TotCorrectedForce,TotConstraints,TotCorConstraints,ConvergeCount 
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

        INTEGER :: l,m,i,j,a,Symm
        REAL*8 :: CurrCoeff(NoOrbs,NoOrbs)
        REAL*8 :: DerivConstr(NoOrbs,NoOrbs,TotNoConstraints)

        call set_timer(CalcDerivConstr_Time,30)


            DerivConstr(:,:,:)=0.D0
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
        INTEGER :: a,l,m,Symm,x,y,w,SymMin,TempMaxOccVirt
        REAL*8 :: TotForce,TotDiffCoeffs,CoeffT2(NoOrbs,NoOrbs)

!        WRITE(6,*) 'DerivCoeff'
!        do m=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(4F20.10)',advance='no') DerivCoeff(a,m)
!                IF(tSpinOrbs) THEN
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
!                IF(tSpinOrbs) THEN
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
                IF(tSpinOrbs) THEN
                    SymM=INT(G1(SymLabelList2(m))%sym%S,4)
                ELSE
                    SymM=INT(G1(SymLabelList2(m)*2)%sym%S,4)
                ENDIF
                do a=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
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


        TotForce=TotForce/(REAL(NoOrbs**2,8))


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
        INTEGER :: l,i,j,m
        REAL*8 :: CurrCoeff(NoOrbs,NoOrbs),TotConstraints,Constraint(TotNoConstraints) 


            TotConstraints=0.D0
            do l=1,TotNoConstraints
                i=lab(1,l)
                j=lab(2,l)
                IF(i.eq.j) THEN
                    Constraint(l)=Dot_Product(CurrCoeff(:,i),CurrCoeff(:,j))-1.D0
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
        INTEGER :: l,n,m,info,ipiv(TotNoConstraints),work,ierr
        CHARACTER(len=*), PARAMETER :: this_routine='FullShake'


        CALL set_timer(FullShake_Time,30)

! FULL MATRIX INVERSION METHOD

! Calculate matrix from the derivatives of the constraints w.r.t the the coefficients at t1 and t2. I.e. the initial 
! coefficients and those that have been moved by the corrected force.

            DerivConstrT1T2(:,:)=0.D0
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
        INTEGER :: m,l,ierr
        CHARACTER(len=*), PARAMETER :: this_routine='ShakeApproximation'


! Use 'shake' algorithm in which the iterative scheme is applied to each constraint in succession.
            WRITE(6,*) 'DerivConstrT1T2Diag calculated from the shake approx'
            
            DerivConstrT1T2Diag(:)=0.D0
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
        REAL*8 :: TotCorConstraints
        INTEGER :: ShakeIteration
        LOGICAL :: tShakeNotConverged


        TotCorConstraints=0.D0
        ConvergeCount=0
        ConstraintCor(:)=0.D0
        do l=1,TotNoConstraints
            i=lab(1,l)
            j=lab(2,l)
            IF(i.eq.j) THEN
                ConstraintCor(l)=Dot_Product(CoeffCorT2(:,i),CoeffCorT2(:,j))-1.D0
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
        INTEGER :: w,x,y,i,a,j,b
        REAL*8 :: TotGSConstraints,GSConstraint(TotNoConstraints),CoeffTemp(SpatOrbs,SpatOrbs)
        REAL*8 , ALLOCATABLE :: Ident(:,:)
        
!        WRITE(6,*) 'The final transformation coefficients before gram schmidt orthonormalisation'
!        do i=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(F20.10)',advance='no') CoeffT1(a,i)
!                IF(tSpinOrbs) THEN
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
        IF(tSpinOrbs.and.(.not.tMaxHLGap)) THEN
            CoeffTemp(:,:)=0.D0
            do i=1,SpatOrbs
                do j=1,SpatOrbs
                    CoeffTemp(i,j)=CoeffT1(2*i,2*j)
                enddo
            enddo

            CALL GRAMSCHMIDT(CoeffTemp,SpatOrbs)

            CoeffT1(:,:)=0.D0
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


! This file is printed to be used to produce cube files from QChem.
! Line 1 is the coefficients of HF spatial orbitals 1 2 3 ... which form transformed orbital 1 etc.

        OPEN(66,FILE='MOTRANSFORM',FORM='UNFORMATTED',access='direct', recl=8)
! Need to put this back into the original order. 
        w=1
        x=1   !keep a counter of record number
        do while (w.le.2)
            do i=1,SpatOrbs
                j=SymLabelListInv(i)
                ! SymLabelList2(i) gives the orbital label (from Dalton or QChem) corresponding to our
                ! label i.
                ! SymLabelListInv(j) therefore gives the label used in CoeffT1 corresponding to the
                ! Qchem/Dalton label j.
                do a=1,SpatOrbs
                    b=SymLabelListInv(a)
                    WRITE(66,rec=x) CoeffT1(b,j)
                    x=x+1
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                enddo
            enddo
            w=w+1
            ! print the whole matrix twice, once for alpha spin, once for beta.
        enddo
        CLOSE(66)
 
        OPEN(67,FILE='MOTRANSFORM02')
! Need to put this back into the original order. 
        w=1
        x=1   !keep a counter of record number
        do while (w.le.2)
            do i=1,SpatOrbs
                j=SymLabelListInv(i)
                ! SymLabelList2(i) gives the orbital label (from Dalton or QChem) corresponding to our
                ! label i.
                ! SymLabelListInv(j) therefore gives the label used in CoeffT1 corresponding to the
                ! Qchem/Dalton label j.
                do a=1,SpatOrbs
                    b=SymLabelListInv(a)
                    WRITE(67,'(F20.10)',advance='no') CoeffT1(b,j)
                    x=x+1
                    ! a/b are the original (HF) orbitals, and i/j the transformed
                enddo
                WRITE(67,*) ''
            enddo
            w=w+1
            ! print the whole matrix twice, once for alpha spin, once for beta.
        enddo
        CLOSE(67)

!        OPEN(71,FILE='TRANSFORMMAT',status='unknown')
!        do i=1,NoOrbs
!            j=SymLabelListInv(i)
!            do a=1,NoOrbs
!                b=SymLabelListInv(a)
!                WRITE(71,*) CoeffT1(b,j)
!            enddo
!        enddo
!        CALL FLUSH(71)
!        CLOSE(71)
      
        
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
        INTEGER :: i,j,k,BinNo,a,b
        REAL*8 :: MaxFII,MinFII,BinIter,BinVal,SingExcit(NoOrbs,NoOrbs)


!<ik|jk> terms where all i,j and k are virtual
        !Coulomb
        ROHistSCijkVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF(FourIndInts(i,k,j,k).ne.0.D0) THEN
                        BinNo=CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCijkVir(2,BinNo)=ROHistSCijkVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Exchange
        ROHistSEijkVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF(FourIndInts(i,k,k,j).ne.0.D0) THEN
                        BinNo=CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEijkVir(2,BinNo)=ROHistSEijkVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !antisymmetric
        ROHistSASijkVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).ne.0.D0) THEN
                        BinNo=CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASijkVir(2,BinNo)=ROHistSASijkVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo


        IF(Iteration.eq.0) THEN
            OPEN(74,FILE='HistHFSingijkVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCijkVir(2,j).ne.0).or.(ROHistSEijkVir(2,j).ne.0).or.(ROHistSASijkVir(2,j).ne.0)) THEN
                    WRITE(74,'(6F20.10)') ROHistSCijkVir(1,j),ROHistSCijkVir(2,j),ROHistSEijkVir(1,j),ROHistSEijkVir(2,j),&
                                                        &ROHistSASijkVir(1,j),ROHistSASijkVir(2,j)
                ENDIF
            enddo
            CLOSE(74)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            OPEN(75,FILE='HistRotSingijkVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCijkVir(2,j).ne.0).or.(ROHistSEijkVir(2,j).ne.0).or.(ROHistSASijkVir(2,j).ne.0)) THEN
                    WRITE(75,'(6F20.10)') ROHistSCijkVir(1,j),ROHistSCijkVir(2,j),ROHistSEijkVir(1,j),ROHistSEijkVir(2,j),&
                                                        &ROHistSASijkVir(1,j),ROHistSASijkVir(2,j)
                ENDIF
            enddo
            CLOSE(75)
        ENDIF


!<ik|jk> where k is occupied, and i and j are both virtual
        !Coulomb
        ROHistSCkOcijVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF(FourIndInts(i,k,j,k).ne.0.D0) THEN
                        BinNo=CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCkOcijVir(2,BinNo)=ROHistSCkOcijVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Exchange
        ROHistSEkOcijVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF(FourIndInts(i,k,k,j).ne.0.D0) THEN
                        BinNo=CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEkOcijVir(2,BinNo)=ROHistSEkOcijVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !antisymmetric 
        ROHistSASkOcijVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).ne.0.D0) THEN
                        BinNo=CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASkOcijVir(2,BinNo)=ROHistSASkOcijVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo


        IF(Iteration.eq.0) THEN
            OPEN(70,FILE='HistHFSingkOcijVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCkOcijVir(2,j).ne.0).or.(ROHistSEkOcijVir(2,j).ne.0).or.(ROHistSASkOcijVir(2,j).ne.0)) THEN
                    WRITE(70,'(6F20.10)') ROHistSCkOcijVir(1,j),ROHistSCkOcijVir(2,j),ROHistSEkOcijVir(1,j),ROHistSEkOcijVir(2,j),&
                                                        &ROHistSASkOcijVir(1,j),ROHistSASkOcijVir(2,j)
                ENDIF
            enddo
            CLOSE(70)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            OPEN(68,FILE='HistRotSingkOcijVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCkOcijVir(2,j).ne.0).or.(ROHistSEkOcijVir(2,j).ne.0).or.(ROHistSASkOcijVir(2,j).ne.0)) THEN
                    WRITE(68,'(6F20.10)') ROHistSCkOcijVir(1,j),ROHistSCkOcijVir(2,j),ROHistSEkOcijVir(1,j),ROHistSEkOcijVir(2,j),&
                                                        &ROHistSASkOcijVir(1,j),ROHistSASkOcijVir(2,j)
                ENDIF
            enddo
            CLOSE(68)
        ENDIF


! <ik|jk> where i and k are both occupied, and j virtual.
        ! Coulomb
        ROHistSCikOcjVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF(FourIndInts(i,k,j,k).ne.0.D0) THEN
                        BinNo=CEILING((FourIndInts(i,k,j,k)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSCikOcjVir(2,BinNo)=ROHistSCikOcjVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Exchange 
        ROHistSEikOcjVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF(FourIndInts(i,k,k,j).ne.0.D0) THEN
                        BinNo=CEILING((FourIndInts(i,k,k,j)-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSEikOcjVir(2,BinNo)=ROHistSEikOcjVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        !Antisymmetrised
        ROHistSASikOcjVir(:,:)=0.D0
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
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                    IF((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j)).ne.0.D0) THEN
                        BinNo=CEILING(((FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))-MinFII)*4002/(MaxFII-MinFII))
                        ROHistSASikOcjVir(2,BinNo)=ROHistSASikOcjVir(2,BinNo)+1.0         
                    ENDIF
                enddo
            enddo
        enddo

        IF(Iteration.eq.0) THEN
            OPEN(61,FILE='HistHFSingikOcjVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCikOcjVir(2,j).ne.0).or.(ROHistSEikOcjVir(2,j).ne.0).or.(ROHistSASikOcjVir(2,j).ne.0)) THEN 
                    WRITE(61,'(6F20.10)') ROHistSCikOcjVir(1,j),ROHistSCikOcjVir(2,j),ROHistSEikOcjVir(1,j),ROHistSEikOcjVir(2,j),&
                                                        &ROHistSASikOcjVir(1,j),ROHistSASikOcjVir(2,j)
                ENDIF
            enddo
            CLOSE(61)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            OPEN(62,FILE='HistRotSingikOcjVir',STATUS='unknown')
            do j=1,4002
                IF((ROHistSCikOcjVir(2,j).ne.0).or.(ROHistSEikOcjVir(2,j).ne.0).or.(ROHistSASikOcjVir(2,j).ne.0)) THEN 
                    WRITE(62,'(6F20.10)') ROHistSCikOcjVir(1,j),ROHistSCikOcjVir(2,j),ROHistSEikOcjVir(1,j),ROHistSEikOcjVir(2,j),&
                                                        &ROHistSASikOcjVir(1,j),ROHistSASikOcjVir(2,j)
                ENDIF
            enddo
            CLOSE(62)
        ENDIF

!Single excitations connected to the HF determinant.
        ROHistSing(:,:)=0.D0
        MaxFII=0.D0
        MinFII=0.D0
        do j=NoOcc+1,NoOrbs
            do i=1,NoOcc
                SingExcit(i,j)=0.D0
                IF(i.eq.j) CYCLE
                a=SymLabelList2(i)
                b=SymLabelList2(j)
                do k=1,NoOcc+1
!                    IF(k.eq.j) CYCLE
!                    IF(k.eq.i) CYCLE
                    SingExcit(i,j)=SingExcit(i,j)+REAL(TMAT2D(2*a,2*b)%v,8)+((2*FourIndInts(i,k,j,k))-FourIndInts(i,k,k,j))
                enddo
                IF(SingExcit(i,j).gt.MaxFII) MaxFII=SingExcit(i,j)
                IF(SingExcit(i,j).lt.MinFII) MinFII=SingExcit(i,j)
            enddo
        enddo
        BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                IF(SingExcit(i,j).ne.0.D0) THEN
                    BinNo=CEILING((SingExcit(i,j)-MinFII)*4002/(MaxFII-MinFII))
                    ROHistSing(2,BinNo)=ROHistSing(2,BinNo)+1.0         
                ENDIF
            enddo
        enddo

        IF(Iteration.eq.0) THEN
            OPEN(56,FILE='HistHFSingExcHF',STATUS='unknown')
            do j=1,4002
                IF(ROHistSing(2,j).ne.0) THEN
                    do i=1,2
                        WRITE(56,'(F20.10)',advance='no') ROHistSing(i,j)
                    enddo
                    WRITE(56,*) ''
                ENDIF
            enddo
            CLOSE(56)
        ENDIF
        IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
            OPEN(60,FILE='HistRotSingExcHF',STATUS='unknown')
            do j=1,4002
                IF(ROHistSing(2,j).ne.0) THEN
                    do i=1,2
                        WRITE(60,'(F20.10)',advance='no') ROHistSing(i,j)
                    enddo
                    WRITE(60,*) ''
                ENDIF
            enddo
            CLOSE(60)
        ENDIF



    ENDSUBROUTINE WriteSingHisttofile 




    SUBROUTINE WriteDoubHisttofile()
        INTEGER :: i,j,k,l,BinNo
        REAL*8 :: MaxFII,MinFII,BinIter,OnePartOrbEnValue,MinHFTMAT,BinVal


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
        
            ROHistDCijOcklVir(:,:)=0.D0
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
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
            ROHistASijOcklVir(:,:)=0.D0
            MinFII=FourIndInts(1,2,NoOcc+1,NoOcc+2)-FourIndInts(1,2,NoOcc+2,NoOcc+1)
            MaxFII=FourIndInts(1,2,NoOcc+1,NoOcc+2)-FourIndInts(1,2,NoOcc+2,NoOcc+1)
            do i=1,NoOcc
                do j=1,NoOcc
                    do k=NoOcc+1,NoOrbs
                        do l=NoOcc+1,NoOrbs
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).lt.MinFII) MinFII=(FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).gt.MaxFII) MaxFII=(FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(54,FILE='HistHFDoubijOcklVir',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijOcklVir(2,j).ne.0).or.(ROHistASijOcklVir(2,j).ne.0)) THEN
                        WRITE(54,'(4F20.10)') ROHistDCijOcklVir(1,j),ROHistDCijOcklVir(2,j),ROHistASijOcklVir(1,j),ROHistASijOcklVir(2,j)
                    ENDIF
                enddo
                CLOSE(54)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                OPEN(55,FILE='HistRotDoubijOcklVir',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijOcklVir(2,j).ne.0).or.(ROHistASijOcklVir(2,j).ne.0)) THEN
                        WRITE(55,'(4F20.10)') ROHistDCijOcklVir(1,j),ROHistDCijOcklVir(2,j),ROHistASijOcklVir(1,j),ROHistASijOcklVir(2,j)
                    ENDIF
                enddo
                CLOSE(55)
            ENDIF



            ROHistDCijklVir(:,:)=0.D0
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
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
            ROHistASijklVir(:,:)=0.D0
            MinFII=FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs-1,NoOrbs)-FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs,NoOrbs-1)
            MaxFII=FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs-1,NoOrbs)-FourIndInts(NoOrbs-3,NoOrbs-2,NoOrbs,NoOrbs-1)
            do i=NoOcc+1,NoOrbs
                do j=NoOcc+1,NoOrbs
                    do k=i+1,NoOrbs
                        do l=j+1,NoOrbs
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).lt.MinFII) MinFII=(FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))
                            IF((FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k)).gt.MaxFII) MaxFII=(FourIndInts(i,j,k,l)-FourIndInts(i,j,l,k))
                        enddo
                    enddo
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(51,FILE='HistHFDoubijklVirt',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijklVir(2,j).ne.0).or.(ROHistASijklVir(2,j).ne.0)) THEN
                        WRITE(51,'(4F20.10)') ROHistDCijklVir(1,j),ROHistDCijklVir(2,j),ROHistASijklVir(1,j),ROHistASijklVir(2,j)
                    ENDIF
                enddo
                CLOSE(51)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                OPEN(52,FILE='HistRotDoubijklVirt',STATUS='unknown')
                do j=1,4002
                    IF((ROHistDCijklVir(2,j).ne.0).or.(ROHistASijklVir(2,j).ne.0)) THEN
                        WRITE(52,'(4F20.10)') ROHistDCijklVir(1,j),ROHistDCijklVir(2,j),ROHistASijklVir(1,j),ROHistASijklVir(2,j)
                    ENDIF
                enddo
                CLOSE(52)
            ENDIF
        ENDIF
   


! Histogramming all one particle orbital energies (occupied and virtual) even though we are not changing occupied.  Would like to see HOMO-LUMO gap etc.
        IF(tROHistOneElInts) THEN

            ROHistHijVirt(:,:)=0.D0
            MinFII=TMAT2DRot(NoOcc+1,NoOcc+2)
            MaxFII=TMAT2DRot(NoOcc+1,NoOcc+2)
            do i=NoOcc+1,NoOrbs
                do j=i+1,NoOrbs
                    IF(TMAT2DRot(i,j).lt.MinFII) MinFII=TMAT2DRot(i,j)
                    IF(TMAT2DRot(i,j).gt.MaxFII) MaxFII=TMAT2DRot(i,j)
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(27,FILE='HistHFHijVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(27,'(F20.10)',advance='no') ROHistHijVirt(i,j)
                        enddo
                        WRITE(27,*) ''
                    ENDIF
                enddo
                CLOSE(27)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                OPEN(28,FILE='HistRotHijVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(28,'(F20.10)',advance='no') ROHistHijVirt(i,j)
                        enddo
                        WRITE(28,*) ''
                    ENDIF
                enddo
                CLOSE(28)
            ENDIF
 


            ROHistHijOccVirt(:,:)=0.D0
            MinFII=TMAT2DRot(1,NoOcc+1)
            MaxFII=TMAT2DRot(1,NoOcc+1)
            do i=1,NoOcc
                do j=NoOcc+1,NoOrbs
                    IF(TMAT2DRot(i,j).lt.MinFII) MinFII=TMAT2DRot(i,j)
                    IF(TMAT2DRot(i,j).gt.MaxFII) MaxFII=TMAT2DRot(i,j)
                enddo
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(25,FILE='HistHFHijOccVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijOccVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(25,'(F20.10)',advance='no') ROHistHijOccVirt(i,j)
                        enddo
                        WRITE(25,*) ''
                    ENDIF
                enddo
                CLOSE(25)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                OPEN(26,FILE='HistRotHijOccVirt',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHijOccVirt(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(26,'(F20.10)',advance='no') ROHistHijOccVirt(i,j)
                        enddo
                        WRITE(26,*) ''
                    ENDIF
                enddo
                CLOSE(26)
            ENDIF
 


            ROHistHii(:,:)=0.D0
            MinFII=TMAT2DRot(1,1)
            MaxFII=TMAT2DRot(1,1)
            do i=1,NoOrbs
                IF(TMAT2DRot(i,i).lt.MinFII) MinFII=TMAT2DRot(i,i)
                IF(TMAT2DRot(i,i).gt.MaxFII) MaxFII=TMAT2DRot(i,i)
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(21,FILE='HistHFHii',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHii(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(21,'(F20.10)',advance='no') ROHistHii(i,j)
                        enddo
                        WRITE(21,*) ''
                    ENDIF
                enddo
                CLOSE(21)
            ENDIF
            IF((.not.tNotConverged).and.(Iteration.gt.1)) THEN
                OPEN(20,FILE='HistRotHii',STATUS='unknown')
                do j=1,4002
                    IF(ROHistHii(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(20,'(F20.10)',advance='no') ROHistHii(i,j)
                        enddo
                        WRITE(20,*) ''
                    ENDIF
                enddo
                CLOSE(20)
            ENDIF
        ENDIF
   

        IF(tROHistOnePartOrbEn) THEN
            ROHistOnePartOrbEn(:,:)=0.D0
            MaxFII=0.D0
            MinFII=0.D0
            do i=1,NoOrbs
                OnePartOrbEnValue=0.D0
                OnePartOrbEnValue=OnePartOrbEnValue+TMAT2DRot(i,i)
                do j=1,NoOcc
                    OnePartOrbEnValue=OnePartOrbEnValue+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                enddo
                IF(OnePartOrbEnValue.gt.MaxFII) MaxFII=OnePartOrbEnValue
                IF(OnePartOrbEnValue.lt.MinFII) MinFII=OnePartOrbEnValue
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
            MaxFII=MaxFII+BinIter
            MinFII=MinFII-BinIter
            BinVal=MinFII
            do i=1,4002
                ROHistOnePartOrbEn(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,NoOrbs
                OnePartOrbEnValue=0.D0
                OnePartOrbEnValue=OnePartOrbEnValue+TMAT2DRot(i,i)
                do j=1,NoOcc
                    OnePartOrbEnValue=OnePartOrbEnValue+(2*FourIndInts(i,j,i,j))-FourIndInts(i,j,j,i)
                enddo
                BinNo=CEILING((OnePartOrbEnValue-MinFII)*4002/(MaxFII-MinFII))
                ROHistOnePartOrbEn(2,BinNo)=ROHistOnePartOrbEn(2,BinNo)+1.0         
            enddo

            IF(Iteration.eq.0) THEN 
                OPEN(23,FILE='HistHFOnePartOrbEn',STATUS='unknown')
                do j=1,4002
                    IF(ROHistOnePartOrbEn(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(23,'(F20.10)',advance='no') ROHistOnePartOrbEn(i,j)
                        enddo
                        WRITE(23,*) ''
                    ENDIF
                enddo
                CLOSE(23)
            ENDIF
            IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
                OPEN(22,FILE='HistRotOnePartOrbEn',STATUS='unknown')
                do j=1,4002
                    IF(ROHistOnePartOrbEn(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(22,'(F20.10)',advance='no') ROHistOnePartOrbEn(i,j)
                        enddo
                        WRITE(22,*) ''
                    ENDIF
                enddo
                CLOSE(22)
            ENDIF
        ENDIF
  

        IF(tROHistDoubExc) THEN
            ROHistDoubExc(:,:)=0.D0
            MaxFII=0.D0
            MinFII=0.D0
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
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(36,FILE='HistHFDoubExc',STATUS='unknown')
                do j=1,4002
                    IF(ROHistDoubExc(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(36,'(F20.10)',advance='no') ROHistDoubExc(i,j)
                        enddo
                        WRITE(36,*) ''
                    ENDIF
                enddo
                CLOSE(36)
            ENDIF
            IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
                OPEN(24,FILE='HistRotDoubExc',STATUS='unknown')
                do j=1,4002
                    IF(ROHistDoubExc(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(24,'(F20.10)',advance='no') ROHistDoubExc(i,j)
                        enddo
                        WRITE(24,*) ''
                    ENDIF
                enddo
                CLOSE(24)
            ENDIF
        ENDIF


        IF(tROHistER) THEN
            ROHistER(:,:)=0.D0
            MaxFII=0.D0
            MinFII=0.D0
            do i=1,NoOrbs
                IF(FourIndInts(i,i,i,i).gt.MaxFII) MaxFII=FourIndInts(i,i,i,i)
                IF(FourIndInts(i,i,i,i).lt.MinFII) MinFII=FourIndInts(i,i,i,i)  
            enddo
            BinIter=ABS(MaxFII-MinFII)/4000.D0
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
                OPEN(32,FILE='HistHF-ER',STATUS='unknown')
                do j=1,4002
                    IF(ROHistER(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(32,'(F20.10)',advance='no') ROHistER(i,j)
                        enddo
                        WRITE(32,*) ''
                    ENDIF
                enddo
                CLOSE(32)
            ENDIF

            IF((Iteration.gt.1).and.(.not.tNotConverged)) THEN
                OPEN(34,FILE='HistRot-ER',STATUS='unknown')
                do j=1,4002
                    IF(ROHistER(2,j).ne.0) THEN
                        do i=1,2
                            WRITE(34,'(F20.10)',advance='no') ROHistER(i,j)
                        enddo
                        WRITE(34,*) ''
                    ENDIF
                enddo
                CLOSE(34)
            ENDIF

        ENDIF


    ENDSUBROUTINE WriteDoubHisttofile 




    SUBROUTINE PrintIntegrals()
        INTEGER :: i,j,k,l
        REAL*8 :: DiagOneElPot,ERPot,ijVirtOneElPot,ijVirtCoulPot,ijVirtExchPot
        REAL*8 :: singCoulijVirt,singExchijVirt,singCoulconHF,singExchconHF,ijklPot,ijklantisymPot,ijOccVirtOneElPot,ijOccVirtCoulPot,ijOccVirtExchPot

        IF(tInitIntValues) THEN
            OPEN(17,FILE='DiagIntegrals',STATUS='unknown')
            WRITE(17,'(A10,6A18)') "Iteration","<i|h|i> ivirt","<ii|ii> ivirt","<ij|ij> iOccjVirt","<ij|ji> iOccjVirt","<ij|ij> ijVirt","<ij|ji> ijVirt"

            OPEN(18,FILE='SingExcIntegrals',STATUS='unknown')
            WRITE(18,'(A10,6A18)') "Iteration","<i|h|j> iOccjVirt","<i|h|j> ijVirt","<ik|jk> HFcon","<ik|kj> HFcon","<ik|jk> ijVirt","<ik|kj> ijVirt"


            DiagOneElPotInit=0.D0
            ERPotInit=0.D0
            ijVirtOneElPotInit=0.D0
            ijVirtCoulPotInit=0.D0
            ijVirtExchPotInit=0.D0
            singCoulconHFInit=0.D0
            singExchconHFInit=0.D0
            singCoulijVirtInit=0.D0
            singExchijVirtInit=0.D0
            ijklPotInit=0.D0
            ijklantisymPotInit=0.D0
            ijOccVirtOneElPotInit=0.D0
            ijOccVirtCoulPotInit=0.D0
            ijOccVirtExchPotInit=0.D0
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
        

        DiagOneElPot=0.D0
        ERPot=0.D0
        ijVirtOneElPot=0.D0
        ijVirtCoulPot=0.D0
        ijVirtExchPot=0.D0
        singCoulconHF=0.D0
        singExchconHF=0.D0
        singCoulijVirt=0.D0
        singExchijVirt=0.D0
        ijklPot=0.D0
        ijklantisymPot=0.D0
        ijOccVirtOneElPot=0.D0
        ijOccVirtCoulPot=0.D0
        ijOccVirtExchPot=0.D0
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


        WRITE(17,'(I10,6F18.10)') Iteration,DiagOneElPot,ERPot,ijOccVirtCoulPot,ijOccVirtExchPot,ijVirtCoulPot,ijVirtExchPot
        WRITE(18,'(I10,6F18.10)') Iteration,ijOccVirtOneElPot,ijVirtOneElPot,singCoulconHF,singExchconHF,singCoulijVirt,singExchijVirt

        IF((.not.tNotConverged).and.(.not.tInitIntValues)) THEN
            CLOSE(17)
            CLOSE(18)
        ENDIF


    ENDSUBROUTINE PrintIntegrals





    SUBROUTINE CalcFOCKMatrix()
        USE SystemData , only : nBasis
        INTEGER :: i,j,k,l,a,b,ierr
        REAL*8 :: FOCKDiagSumHF,FOCKDiagSumNew,ArrTemp(nBasis)
        CHARACTER(len=*) , PARAMETER :: this_routine='CalcFOCKMatrix'
        !NEED TO FIX THIS!

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

    
        ALLOCATE(ArrNew(nBasis,nBasis),stat=ierr)
        CALL LogMemAlloc('ArrNew',nBasis**2,8,this_routine,ArrNewTag,ierr)
        ArrNew(:,:)=0.D0                     

!        WRITE(6,*) 'The diagonal fock elements in the HF basis set'
!        do a=1,nBasis
!            WRITE(6,'(F20.10)',advance='no') Arr(a,2)
!        enddo


! First calculate the sum of the diagonal elements, ARR.
! Check if this is already being done.
        FOCKDiagSumHF=0.D0
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

        FOCKDiagSumNew=0.D0
        do j=1,NoRotOrbs
            l=SymLabelList3(j)
            do i=1,NoRotOrbs
                k=SymLabelList3(i)
                ArrNew(k,l)=0.D0
                do a=1,NoOrbs
                    b=SymLabelList2(a)
                    IF(tSpinOrbs) THEN
                        ArrNew(k,l)=ArrNew(k,l)+(CoeffT1(a,i)*Arr(b,2)*CoeffT1(a,j))
                    ELSE
                        ArrNew(k,l)=ArrNew(k,l)+(CoeffT1(a,i)*Arr(2*b,2)*CoeffT1(a,j))
                    ENDIF
                enddo
            enddo
            IF(tSpinOrbs) THEN
                FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(j,j))
            ELSE
                FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(j,j)*2)
            ENDIF
            !only running through spat orbitals, count each twice to compare to above.
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

!        WRITE(6,*) 'BRR then ARR before being changed'
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
       

! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2) (ordered in terms of orbital number).
! ARR(:,2) needs to be ordered in terms of symmetry and then energy (like SymLabelList), so currently this ordering will not be 
! correct when reading in qchem INTDUMPS as the orbital number ordering is by energy.

        IF(tSpinOrbs) THEN
            do j=1,NoRotOrbs
                ARR(j,2)=ArrNew(j,j)
                ARRTemp(j)=ArrNew(BRR(j),BRR(j))
            enddo
        ELSE
            do j=1,NoRotOrbs
                ARR(2*j,2)=ArrNew(j,j)
                ARR(2*j-1,2)=ArrNew(j,j)
                ARRTemp(2*j)=ArrNew(BRR(2*j)/2,BRR(2*j)/2)
                ARRTemp(2*j-1)=ArrNew(BRR(2*j)/2,BRR(2*j)/2)
            enddo
        ENDIF
! fill both ARR(:,1) and ARR(:,2) with the diagonal values, then sort ARR(:,1) according to the values (energies), 
! taking the orbital labels (BRR) with them.

!        WRITE(6,*) 'ARRtemp'
!        do i=1,nBasis
!            WRITE(6,*) BRR(i),ARRTemp(i)
!        enddo

!        CALL NECI_SORT2(nBasis,ARRTemp,BRR)

        IF(tSpinOrbs) THEN
            do j=1,(nBasis-NoFrozenVirt)
                ARR(j,1)=ArrTemp(j)
            enddo
        ELSE
            do j=1,(nBasis-(NoFrozenVirt*2))
                ARR(j,1)=ArrTemp(j)
            enddo
        ENDIF


!        WRITE(6,*) 'BRR then ARR after being changed'
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
!        CALL FLUSH(6)
!        stop       

        DEALLOCATE(ArrNew)
        CALL LogMemDealloc(this_routine,ArrNewTag)


    ENDSUBROUTINE CalcFOCKMatrix


    SUBROUTINE RefillUMATandTMAT2D()
        INTEGER :: l,k,j,i,a,b,g,d,c,BinNo,l2,k2,j2,i2,nBasis2
        REAL*8 :: NewTMAT,TMAT2DPart((nBasis-(NoFrozenVirt*2)),nBasis)


! Make the UMAT elements the four index integrals.  These are calculated by transforming the HF orbitals using the
! coefficients that have been found
        do l=1,(NoOrbs-(NoFrozenVirt))
!            IF(tSpinOrbs) THEN
!                l2=l*2-1
!                d=CEILING(SymLabelList3(2*l)/2.0)
!            ELSE
                d=SymLabelList3(l)
!            ENDIF
            do k=1,(NoOrbs-(NoFrozenVirt))
!                IF(tSpinOrbs) THEN
!                    k2=k*2
!                    g=CEILING(SymLabelList3(2*k)/2.0)
!                ELSE
                    g=SymLabelList3(k)
!                ENDIF
 
                do j=1,(NoOrbs-(NoFrozenVirt))
!                    IF(tSpinOrbs) THEN
!                        j2=j*2-1
!                        b=CEILING(SymLabelList3(2*j)/2.0)
!                    ELSE
                        b=SymLabelList3(j)
!                    ENDIF
                    do i=1,(NoOrbs-(NoFrozenVirt))
!                        IF(tSpinOrbs) THEN
!                            i2=i*2
!                            a=CEILING(SymLabelList3(2*i)/2.0)
!                            UMAT(UMatInd(a,b,g,d,0,0))=HElement(FourIndInts(i2,j2,k2,l2))
!                        ELSE
                            a=SymLabelList3(i)
                            UMAT(UMatInd(a,b,g,d,0,0))=HElement(FourIndInts(i,j,k,l))
!                        ENDIF
                    enddo
                enddo
            enddo
        enddo

! Also calculate the 2 index integrals, and make these the elements of the TMAT2D matrix.
! TMAT2D is in spin orbitals.

!        WRITE(6,*) 'TMAT2D before transformation' 
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l)%v,8)
!            enddo
!            WRITE(6,*) ''
!        enddo


        do a=1,nBasis
            do k=1,NoRotOrbs
                i=SymLabelList3(k)
                NewTMAT=0.D0
                do b=1,NoOrbs
                    d=SymLabelList2(b)
                    IF(tSpinOrbs) THEN
                        NewTMAT=NewTMAT+(CoeffT1(b,k)*REAL(TMAT2D(d,a)%v,8))
                    ELSE
                        NewTMAT=NewTMAT+(CoeffT1(b,k)*REAL(TMAT2D(2*d,a)%v,8))
                    ENDIF
                enddo
                IF(tSpinOrbs) THEN
                    TMAT2DPart(i,a)=NewTMAT
                ELSE
                    TMAT2DPart(2*i,a)=NewTMAT
                    TMAT2DPart(2*i-1,a)=NewTMAT
                ENDIF
            enddo
        enddo
        IF(tSpinOrbs) THEN
            nBasis2=nBasis-NoFrozenVirt
        ELSE
            nBasis2=nBasis-(NoFrozenVirt*2)
        ENDIF
        do k=1,nBasis2
            do l=1,NoRotOrbs
                j=SymLabelList3(l)
                NewTMAT=0.D0
                do a=1,NoOrbs
                    c=SymLabelList2(a)
                    IF(tSpinOrbs) THEN
                        NewTMAT=NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,c))
                    ELSE
                        NewTMAT=NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,2*c))
                    ENDIF
                enddo
                IF(tSpinOrbs) THEN
                    TMAT2D(k,j)=HElement(NewTMAT)
                ELSE
                    TMAT2D(k,2*j)=HElement(NewTMAT)
                    TMAT2D(k,2*j-1)=HElement(NewTMAT)
                ENDIF
            enddo
        enddo

    
!        WRITE(6,*) 'TMAT2D after transformation'
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l)%v,8)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        CALL FLUSH(6)
!        stop


        IF(tROHistSingExc) CALL WriteSingHisttofile()

        WRITE(6,*) 'Printing the new ROFCIDUMP file.'
        IF(tROFciDump) CALL PrintROFCIDUMP()


    ENDSUBROUTINE RefillUMATandTMAT2D


    

    SUBROUTINE PrintROFCIDUMP()
!This prints out a new FCIDUMP file in the same format as the old one.
        INTEGER :: i,j,k,l,Sym,Syml,Symk


        OPEN(48,FILE='ROFCIDUMP',STATUS='unknown')
        
        WRITE(48,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',(NoOrbs-(NoFrozenVirt)),',NELEC=',NEl,',MS2=',LMS,','
        WRITE(48,'(A9)',advance='no') 'ORBSYM='
        do i=1,(NoOrbs-(NoFrozenVirt))
            IF((tUseMP2VarDenMat.or.tFindCINatOrbs).and.(.not.lNoSymmetry).and.tTruncRODump) THEN
                WRITE(48,'(I1,A1)',advance='no') (SymOrbs(i)+1),','
            ELSE
                IF(tSpinOrbs) THEN
                    WRITE(48,'(I1,A1)',advance='no') (INT(G1(i)%sym%S,4)+1),','
                ELSE
                    WRITE(48,'(I1,A1)',advance='no') (INT(G1(i*2)%sym%S,4)+1),','
                ENDIF
            ENDIF
        enddo
        WRITE(48,*) ''
        IF(tStoreSpinOrbs) THEN
            WRITE(48,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        ELSE
            WRITE(48,'(A7,I1)') 'ISYM=',1
        ENDIF
        WRITE(48,'(A5)') '&END'
       
        do i=1,(NoOrbs-(NoFrozenVirt))
            do k=1,i
                do j=1,(NoOrbs-(NoFrozenVirt))
!                    Sym=IEOR(INT(G1(j*2)%sym%S,4),IEOR(INT(G1(k*2)%sym%S,4),INT(G1(i*2)%sym%S,4)))
                    ! Potential to put symmetry in here, have currently taken it out, because when we're only printing non-zero values,
                    ! it is kind of unnecessary - although it may be used to speed things up.
                    do l=1,j
!                        Syml=INT(G1(l*2)%sym%S,4)
!                        IF((Syml.eq.Sym).and.((REAL(UMat(UMatInd(i,j,k,l,0,0))%v,8)).ne.0.D0)) &
                        IF((ABS(REAL(UMat(UMatInd(i,j,k,l,0,0))%v,8))).ne.0.D0) &
                                        &WRITE(48,'(F21.12,4I3)') REAL(UMat(UMatInd(i,j,k,l,0,0))%v,8),i,k,j,l 
                    enddo
                enddo
           enddo
        enddo

       
!        do i=1,SpatOrbs
!            do j=1,SpatOrbs
!                do l=j,SpatOrbs
!                    Sym=IEOR(INT(G1(l*2)%sym%S,4),IEOR(INT(G1(j*2)%sym%S,4),INT(G1(i*2)%sym%S,4)))
!                    do l=SymLabelCounts2(1,Sym+SymMin),(SymLabelCounts2(1,Sym+SymMin)+SymLabelCounts2(2,Sym+SymMin)-1)
                   ! Potential to put symmetry in here.
!                    do k=i,SpatOrbs
!                        Symk=INT(G1(k*2)%sym%S,4)
!                        IF(Symk.eq.Sym) WRITE(48,'(F21.12,4I3)') REAL(UMat(UMatInd(i,j,k,l,0,0))%v,8),i,k,j,l 
!                    enddo
!                enddo
!            enddo
!        enddo


! TMAT2D stored as spin orbitals
        do k=1,(NoOrbs-(NoFrozenVirt))
            ! Symmetry?
            do i=k,(NoOrbs-(NoFrozenVirt))
                IF(tSpinOrbs) THEN
                    IF((REAL(TMAT2D(i,k)%v,8)).ne.0.D0) WRITE(48,'(F21.12,4I3)') REAL(TMAT2D(i,k)%v,8),i,k,0,0
                ELSE
                    IF((REAL(TMAT2D(2*i,2*k)%v,8)).ne.0.D0) WRITE(48,'(F21.12,4I3)') REAL(TMAT2D(2*i,2*k)%v,8),i,k,0,0
                ENDIF
            enddo
        enddo

! ARR has the energies of the orbitals (eigenvalues).  ARR(:,2) has ordering we want.
! ARR is stored as spin orbitals.

        do k=1,(NoOrbs-(NoFrozenVirt))
            IF(tSpinOrbs) THEN
                WRITE(48,'(F21.12,4I3)') Arr(k,2),k,0,0,0
            ELSE
                WRITE(48,'(F21.12,4I3)') Arr(2*k,2),k,0,0,0
            ENDIF
        enddo

        WRITE(48,'(F21.12,4I3)') ECore,0,0,0,0
        
        CLOSE(48)


    ENDSUBROUTINE PrintROFCIDUMP



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

        IF(tERLocalization.and.(.not.tSpinOrbs)) THEN
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
        DEALLOCATE(SymLabelList2)
        CALL LogMemDealloc(this_routine,SymLabelList2Tag)
        DEALLOCATE(SymLabelCounts2)
        CALL LogMemDealloc(this_routine,SymLabelCounts2Tag)
        DEALLOCATE(SymLabelListInv)
        CALL LogMemDealloc(this_routine,SymLabelListInvTag)

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
        INTEGER :: i,j,x,y,l,m,a

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
        CALL FLUSH(6) 

    ENDSUBROUTINE WriteShakeOUTstats01



    SUBROUTINE WriteShakeOUTstats02(ShakeIteration,TotLambdas,ConvergeCount)
! Debugging option    
        INTEGER :: m,a,l,ConvergeCount,ShakeIteration
        REAL*8 :: TotLambdas 
    
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
