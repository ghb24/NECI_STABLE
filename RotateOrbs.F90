!TO DO
!Spatial symmetry
!Finalize orbs to fully test
!Permutational symmetry where possible
!Parallelize

MODULE RotateOrbsMod

    USE Global_utilities
    USE IntegralsData , only : UMAT,nFrozen
    USE UMatCache , only : UMatInd
    USE HElem , only : HElement
    USE SystemData , only : ConvergedForce,TimeStep,tLagrange,tShake,tShakeApprox,ShakeConverged,tROIteration,ROIterMax,tShakeIter,ShakeIterMax
    USE SystemData , only : G1,ARR,NEl,nBasis,LMS,ECore,tSeparateOccVirt,Brr,nBasisMax,OrbOrder,lNoSymmetry,tRotatedOrbs,tERLocalization
    USE Logging , only : tROHistogram,tROFciDump,tERHist
    USE OneEInts , only : TMAT2D
    USE SymData , only : TwoCycleSymGens,SymLabelList,SymLabelCounts
    USE Timing , only : end_timing,print_timing_report
    USE Soft_exit, only : test_SOFTEXIT
    IMPLICIT NONE
    INTEGER , ALLOCATABLE :: Lab(:,:),LabVirtOrbs(:),LabOccOrbs(:),SymLabelCounts2(:,:)
    REAL*8 , ALLOCATABLE :: CoeffT1(:,:),CoeffCorT2(:,:),CoeffUncorT2(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:),ArrNew(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:),UMATTemp01(:,:,:,:),UMATTemp02(:,:,:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:),ForceCorrect(:,:),Correction(:,:),ShakeLambdaNew(:),ConstraintCor(:)
    REAL*8 , ALLOCATABLE :: Constraint(:),ShakeLambda(:),DerivConstrT1(:,:,:),DerivConstrT2(:,:,:),DerivConstrT1T2(:,:),DerivConstrT1T2Diag(:),FourIndInts(:,:,:,:)
!    REAL*8 , ALLOCATABLE :: PartTransfInts01(:,:,:,:),PartTransfInts02(:,:,:,:),PartTransfInts03(:,:,:,:),PartTransfInts04(:,:,:,:),PartTransfInts05(:,:,:,:)
    REAL*8 , ALLOCATABLE :: OneIndInts01(:,:,:,:),OneIndInts02(:,:,:,:),TwoIndInts01(:,:,:,:),TwoIndInts02(:,:,:,:),ThreeIndInts01(:,:,:,:),FourIndInts02(:,:,:,:)
    REAL*8 , ALLOCATABLE :: ThreeIndInts02(:,:,:,:),ThreeIndInts03(:,:,:,:),ThreeIndInts04(:,:,:,:),ThreeIndInts(:,:,:,:)   
    !These are the arrays to store the
!partially transformed integrals for each iteration. The FourIndInts are the full <ij|kl> integrals for the current iteration.
!OneIndInts = <i b|g d> ; TwoIndInts = <i b|k d> ; ThreeIndInts = <i j|k d>
!    INTEGER :: PartTransfInts01Tag,PartTransfInts02Tag,PartTransfInts03Tag,PartTransfInts04Tag,PartTransfInts05Tag
    INTEGER :: OneIndInts01Tag,OneIndInts02Tag,TwoIndInts01Tag,TwoIndInts02Tag,ThreeIndInts01Tag,ThreeIndInts02Tag,ThreeIndInts03Tag,ThreeIndInts04Tag,FourIndInts02Tag
    INTEGER :: LabTag,ForceCorrectTag,CorrectionTag,SpatOrbs,FourIndIntsTag,ArrNewTag,UMATTemp01Tag,UMATTemp02Tag,ThreeIndIntsTag
    INTEGER :: CoeffT1Tag,CoeffCorT2Tag,CoeffUncorT2Tag,LambdasTag,DerivCoeffTag,DerivLambdaTag,Iteration,TotNoConstraints,ShakeLambdaNewTag
    INTEGER :: ShakeLambdaTag,ConstraintTag,ConstraintCorTag,DerivConstrT1Tag,DerivConstrT2Tag,DerivConstrT1T2Tag,DerivConstrT1T2DiagTag
    INTEGER :: LabVirtOrbsTag,LabOccOrbsTag,OccVirt,MinMZ,MaxMZ,SymLabelCounts2Tag
    LOGICAL :: tNotConverged
    REAL*8 :: OrthoNorm,PotEnergy,Force,TwoEInts,DistCs,OrthoForce,DistLs,LambdaMag,PEInts,PEOrtho,ForceInts,TotCorrectedForce
    REAL*8 :: OrthoFac=1.D0,ROHistSing(2,2000),ROHistDoub(2,2000),ROHistER(2,2000)

    TYPE(timer), save :: Rotation_Time,FullShake_Time,Shake_Time,Findtheforce_Time,Transform2ElInts_Time,findandusetheforce_time,CalcDerivConstr_Time,TestOrthoConver_Time


    contains

    SUBROUTINE RotateOrbs()


!Sym of spatial orbital i is INT(G1(2*i)%sym%S,4) - should go from 0 -> 7

        CALL InitLocalOrbs()        ! Set defaults, allocate arrays, write out headings for OUTPUT, set integarals to HF values.

        CALL WriteStats()           ! write out the original stats before any rotation.

        CALL set_timer(Rotation_Time,30)


!        IF(tROFciDump) THEN
!           CALL CalcFOCKMatrix()
!           CALL RefillUMATandTMAT2D()        
!           CALL PrintROFCIDUMP()
!        ENDIF
!        stop

!        CALL InitOrbitalSeparation()        
        
        tNotConverged=.true.

        do while(tNotConverged)     ! rotate the orbitals until the sum of the four index integral falls below a chose convergence value.

            Iteration=Iteration+1
            
            CALL FindNewOrbs()      ! bulk of the calculation.
                                    ! do the actual transformations, moving the coefficients by a timestep according to the calculated force. 

            CALL WriteStats()       ! write out the stats for this iteration.

        enddo           

        CALL halt_timer(Rotation_Time)

        WRITE(6,*) "Convergence criterion met. Finalizing new orbitals..."

!Make symmetry, orbitals, one/two-electron integrals consistent with rest of neci
        CALL FinalizeNewOrbs()

!        CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
        CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)

        CALL DeallocateMem()

        CALL FLUSH(6)
        CALL FLUSH(12)
!        CALL end_timing()
!        CALL print_timing_report()
!        CALL Stop_All("RotateOrbs","This code is still in the testing phase")

    END SUBROUTINE RotateOrbs



    SUBROUTINE FindNewOrbs()
           

        CALL Transform2ElInts()     ! Find the partially (and completely) transformed 4 index integrals to be used in further calcs.
       

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

    
    
    SUBROUTINE TestOrthonormality()
        INTEGER :: i,j,a
        REAL*8 :: OrthoNormDP

        OrthoNorm=0.D0
        do i=1,SpatOrbs
            do j=1,i
                OrthoNormDP=0.D0
                OrthoNormDP=Dot_Product(CoeffT1(:,i),CoeffT1(:,j))
                OrthoNorm=OrthoNorm+ABS(OrthoNormDP)
            enddo
        enddo
        OrthoNorm=OrthoNorm-real(SpatOrbs,8)
        OrthoNorm=(OrthoNorm*2.D0)/REAL((SpatOrbs*(SpatOrbs+1.D0)),8)

    END SUBROUTINE TestOrthonormality




    SUBROUTINE FindTheForce
        INTEGER :: m,z,i,j,k,l,a,b,g,d,Symm,Symb,Symb2,Symi,Symd,w,x,y,SymMin
        REAL*8 :: Deriv4indint,Deriv4indintsqrd,t1,t2,t3,t4,LambdaTerm1,LambdaTerm2,t1Temp
        CHARACTER(len=*) , PARAMETER :: this_routine='FindtheForce'
        LOGICAL :: leqm,jeqm,keqm,ieqm
      
! Running over m and z, covers all matrix elements of the force matrix (derivative 
! of equation we are minimising, with respect to each translation coefficient) filling 
! them in as it goes.
        CALL set_timer(FindtheForce_time,30)
        

        DerivCoeff(:,:)=0.D0
        Force=0.D0
        ForceInts=0.D0
        OrthoForce=0.D0


        do w=1,OccVirt
            IF(w.eq.1) THEN
                SymMin=1
                MinMZ=1
                IF(tSeparateOccVirt) THEN
                    MaxMZ=NEl/2
                ELSE
                    MaxMZ=SpatOrbs
                ENDIF
            ELSE
                SymMin=9
                MinMZ=(NEl/2)+1
                MaxMZ=SpatOrbs
            ENDIF
! If we are localising the occupied and virtual orbitals separately, the above block ensures that we loop over
! first the occupied then the virtual.  If we are not separating the orbitals we just run over all orbitals.

            do x=MinMZ,MaxMZ
                m=SymLabelList(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM=INT(G1(m*2)%sym%S,4)
                do y=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
                    z=SymLabelList(y)
 
!                    WRITE(6,*) 'x,m,y,z,SymM,SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)'
!                    WRITE(6,*) x,m,y,z,SymM,SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
!                    stop

                    Deriv4indintsqrd=0.D0

                    IF(tERLocalization) THEN
                        ! This is actually just the derivative of <ii|ii> (not squared) but I am just keeping
                        ! the naming consistent.
                        Deriv4indintsqrd=Deriv4indintsqrd+ThreeIndInts01(m,m,m,z)+ThreeIndInts02(m,m,m,z) &
                                                         +ThreeIndInts03(m,m,m,z)+ThreeIndInts04(m,m,m,z)

                    ELSE
! This runs over the full <ij|kl> integrals from the previous iteration. 
! In the future we can take advantage of the permutational symmetry of the matrix elements
                        do l=1,SpatOrbs
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
                                do k=1,SpatOrbs                 
                                    IF(k.eq.m) THEN
                                        keqm=.true.
                                    ELSE
                                        keqm=.false.
                                    ENDIF
!                                    Symi=IEOR(INT(G1(k*2)%sym%S,4),IEOR(INT(G1(j*2)%sym%S,4),INT(G1(l*2)%sym%S,4)))

                                    ! only i with symmetry equal to j x k x l will have integrals with overall
                                    ! symmetry A1 and therefore be non-zero.

                                    ! running across i, ThreeIndInts01 only contributes if i.eq.m (which will happen once for each m)
                                    IF(m.lt.k) Deriv4indintsqrd=Deriv4indintsqrd+2*(FourIndInts02(j,k,l,m)*ThreeIndInts01(k,j,l,z))

                                    IF(jeqm) THEN
                                        do i=1,k-1
                                            Deriv4indintsqrd=Deriv4indintsqrd+2*(FourIndInts(i,j,k,l)*ThreeIndInts02(i,k,l,z))
                                        enddo
                                    ENDIF

                                    IF(keqm) THEN
                                        do i=1,k-1
                                            Deriv4indintsqrd=Deriv4indintsqrd+2*(FourIndInts(i,j,k,l)*ThreeIndInts03(i,j,l,z))
                                        enddo
                                    ENDIF

                                    IF(leqm) THEN
                                        do i=1,k-1
                                            Deriv4indintsqrd=Deriv4indintsqrd+2*(FourIndInts(i,j,k,l)*ThreeIndInts04(i,k,j,z))
                                        enddo
                                    ENDIF

! Already calculated the four index integrals and some partial integrals
! in Transform2ElInts routine of the previous iteration.

! Deriv4indintsqrd is the derivative of the overall expression for the sum of the squares of the <ij|kl> matrix.
! This accumulates as the loop sums over i,j,k and l.
                                       
                                enddo
                            enddo
                        enddo
                    ENDIF
                    DerivCoeff(z,m)=Deriv4indintsqrd  
                    Force=Force+DerivCoeff(z,m)
                enddo
            enddo
        enddo

        Force=Force/REAL(SpatOrbs**2,8)


!        WRITE(6,*) 'found the force'
!        do m=1,SpatOrbs
!            do z=1,SpatOrbs
!                WRITE(6,'(F15.10)',advance='no') DerivCoeff(z,m)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        stop

! Calculate the derivatives of orthogonalisation condition.
! Have taken this out of the m and z loop to make the shake faster, but can put it back in if start using it a lot.
        IF(tLagrange) THEN
            do x=MinMZ,MaxMZ
                m=SymLabelList(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM=INT(G1(m*2)%sym%S,4)
                do y=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
                    z=SymLabelList(y)
      
                    LambdaTerm1=0.D0
                    LambdaTerm2=0.D0
                    
                    do j=1,SpatOrbs
                        LambdaTerm1=LambdaTerm1+(Lambdas(m,j)*CoeffT1(z,j))
                        LambdaTerm2=LambdaTerm2+(Lambdas(j,m)*CoeffT1(z,j))
                    enddo

! DerivCoeff is 'the force'.  I.e. the derivative of |<ij|kl>|^2 with 
! respect to each transformation coefficient.  It is the values of this matrix that will tend to 0 as
! we minimise the sum of the |<ij|kl>|^2 values.
! With the Lagrange keyword this includes orthonormality conditions, otherwise it is simply the unconstrained force.
                    DerivCoeff(z,m)=(2*Deriv4indintsqrd)-LambdaTerm1-LambdaTerm2
                    OrthoForce=OrthoForce-LambdaTerm1-LambdaTerm2
                enddo
            enddo
 
!If doing a lagrange calc we also need to find the force on the lambdas to ensure orthonormality...
            OrthoForce=OrthoForce/REAL(SpatOrbs**2,8)
            DerivLambda(:,:)=0.D0
            do i=1,SpatOrbs
                do j=1,i
                    do a=1,SpatOrbs
                        DerivLambda(i,j)=DerivLambda(i,j)+CoeffT1(a,i)*CoeffT1(a,j)
                    enddo
                    DerivLambda(j,i)=DerivLambda(i,j)
                enddo
            enddo
            do i=1,SpatOrbs
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


        do w=1,OccVirt
            IF(w.eq.1) THEN
                SymMin=1
                MinMZ=1
                IF(tSeparateOccVirt) THEN
                    MaxMZ=NEl/2
                ELSE
                    MaxMZ=SpatOrbs
                ENDIF
            ELSE
                SymMin=9
                MinMZ=(NEl/2)+1
                MaxMZ=SpatOrbs
            ENDIF

            do x=MinMZ,MaxMZ
                m=SymLabelList(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM=INT(G1(m*2)%sym%S,4)
                do y=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
                    z=SymLabelList(y)
                    ! Only coeffs with sym of m and z the same have non-zero coeffs.    
                    NewCoeff=0.D0
                    IF(tERLocalization) THEN
                        NewCoeff=CoeffT1(z,m)+(TimeStep*DerivCoeff(z,m))
                        ! Edmiston Reudenberg localisation MAXimises <ii|ii>, whereas in the other case we are minimising.
                    ELSE
                        NewCoeff=CoeffT1(z,m)-(TimeStep*DerivCoeff(z,m))
                    ENDIF
                    DistCs=DistCs+abs(TimeStep*DerivCoeff(z,m))
                    CoeffT1(z,m)=NewCoeff
                enddo
            enddo
        enddo

        DistCs=DistCs/(REAL(SpatOrbs**2,8))

        
!        IF(tSeparateOccVirt) CALL ZeroOccVirtElements(CoeffT1)


        IF(tLagrange) THEN

            DistLs=0.D0
            LambdaMag=0.D0
            do i=1,SpatOrbs
                do j=1,SpatOrbs
                    NewLambda=0.D0
                    NewLambda=Lambdas(i,j)-(TimeStep*DerivLambda(i,j))  ! Timestep must be specified in the input file.
                    DistLs=DistLs+abs(TimeStep*DerivLambda(i,j))
                    Lambdas(i,j)=NewLambda
                    LambdaMag=LambdaMag+abs(NewLambda)
                enddo
            enddo
            DistLs=DistLs/(REAL(SpatOrbs**2,8))
            LambdaMag=LambdaMag/(REAL(SpatOrbs**2,8))

!        ELSE
!            CALL OrthoNormx(SpatOrbs,SpatOrbs,CoeffT1) !Explicitly orthonormalize the coefficient vectors.
        ENDIF

    ENDSUBROUTINE UseTheForce


    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.
    SUBROUTINE Transform2ElInts()
        INTEGER :: i,j,k,l,a,b,g,d
        REAL*8 :: t,FourIndIntMax,Temp4indints(SpatOrbs,SpatOrbs)
        
        
        CALL set_timer(Transform2ElInts_time,30)


!Zero arrays from previous transform
        OneIndInts01(:,:,:,:)=0.D0
        OneIndInts02(:,:,:,:)=0.D0

        TwoIndInts01(:,:,:,:)=0.D0
        TwoIndInts02(:,:,:,:)=0.D0

        ThreeIndInts(:,:,:,:)=0.D0
        ThreeIndInts01(:,:,:,:)=0.D0
        ThreeIndInts02(:,:,:,:)=0.D0
        ThreeIndInts03(:,:,:,:)=0.D0
        ThreeIndInts04(:,:,:,:)=0.D0


        FourIndInts(:,:,:,:)=0.D0


!<alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i,j,k,l,0,0)

!        do i=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(F20.10)',advance='no') CoeffT1(a,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

! Calculating the one-transformed, four index integrals.
        do b=1,SpatOrbs
            do d=1,b
                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,UMatTemp01(:,:,d,b),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                do g=1,SpatOrbs
                    do i=1,SpatOrbs
                        OneIndInts01(g,i,b,d)=Temp4indints(i,g)
                        OneIndInts01(g,i,d,b)=Temp4indints(i,g)
                    enddo
                enddo
            enddo
        enddo

! Only the fourindints are required once this routine is called to finalize the new orbs, these calculations are unnecessary. 
        IF(tNotConverged) THEN
            do g=1,SpatOrbs
                do a=1,g
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,UMatTemp02(:,:,a,g),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                    do b=1,SpatOrbs
                        do l=1,SpatOrbs
                            OneIndInts02(b,l,g,a)=Temp4indints(l,b)
                            OneIndInts02(b,l,a,g)=Temp4indints(l,b)
                        enddo
                    enddo
                enddo
            enddo
        ENDIF




! Calculating the two-transformed, four index integrals.        
        
        do b=1,SpatOrbs
            do d=1,b
                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,OneIndInts01(:,:,d,b),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                do i=1,SpatOrbs
                    do k=1,i
                        TwoIndInts01(d,b,k,i)=Temp4indints(k,i)
                        TwoIndInts01(b,d,k,i)=Temp4indints(k,i)
                        TwoIndInts01(d,b,i,k)=Temp4indints(k,i)
                        TwoIndInts01(b,d,i,k)=Temp4indints(k,i)
                    enddo
                enddo
            enddo
        enddo


        IF(tNotConverged) THEN
            do g=1,SpatOrbs
                do a=1,g
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,OneIndInts02(:,:,a,g),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                    do l=1,SpatOrbs
                        do j=1,l
                            TwoIndInts02(g,a,j,l)=Temp4indints(j,l)
                            TwoIndInts02(a,g,j,l)=Temp4indints(j,l)
                            TwoIndInts02(g,a,l,j)=Temp4indints(j,l)
                            TwoIndInts02(a,g,l,j)=Temp4indints(j,l)
                        enddo
                    enddo
                enddo
            enddo
        ENDIF



! Calculating the 3 transformed, 4 index integrals. 01=a untransformed,02=b,03=g,04=d

        do i=1,SpatOrbs
            do k=1,i
                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,TwoIndInts01(:,:,k,i),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                do b=1,SpatOrbs
                    do l=1,SpatOrbs
                        ThreeIndInts(b,l,k,i)=Temp4indints(l,b)
                        ThreeIndInts(b,l,i,k)=Temp4indints(l,b)
                        ThreeIndInts02(i,k,l,b)=Temp4indints(l,b)
                        ThreeIndInts02(k,i,l,b)=Temp4indints(l,b)
                    enddo
                enddo
                IF(tNotConverged) THEN
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,TwoIndInts01(:,:,k,i),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                    do d=1,SpatOrbs
                        do j=1,SpatOrbs
                            ThreeIndInts04(k,i,j,d)=Temp4indints(j,d)
                            ThreeIndInts04(i,k,j,d)=Temp4indints(j,d)
                        enddo
                    enddo
                ENDIF
            enddo
        enddo


        IF(tNotConverged) THEN
            do l=1,SpatOrbs
                do j=1,l
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,TwoIndInts02(:,:,j,l),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                    do a=1,SpatOrbs
                        do k=1,SpatOrbs
                            ThreeIndInts01(k,j,l,a)=Temp4indints(k,a)
                            ThreeIndInts01(k,l,j,a)=Temp4indints(k,a)
                        enddo
                    enddo
                    Temp4indints(:,:)=0.D0
                    CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,TwoIndInts02(:,:,j,l),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                    do g=1,SpatOrbs
                        do i=1,SpatOrbs
                            ThreeIndInts03(i,l,j,g)=Temp4indints(i,g)
                            ThreeIndInts03(i,j,l,g)=Temp4indints(i,g)
                        enddo
                    enddo
                enddo
            enddo
        ENDIF

                    
!        do j=1,3,2
!            do k=1,SpatOrbs,2
!                do i=1,k
!                    do d=1,SpatOrbs
!                        WRITE(6,'(4I5,2F20.10)') d,i,k,j,ThreeIndInts04(d,i,k,j)
!                    enddo
!                enddo
!            enddo
!        enddo

!        CALL FLUSH(6)
!        CALL end_timing()
!        CALL print_timing_report()
!        stop


 


! We now have all the 3-transformed, 4 index integrals needed to calculate the force.


        PotEnergy=0.D0
        TwoEInts=0.D0
        PEInts=0.D0
        FourIndIntMax=0.D0
        do i=1,SpatOrbs
            do k=1,i
                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',SpatOrbs,SpatOrbs,SpatOrbs,1.0,CoeffT1(:,:),SpatOrbs,ThreeIndInts(:,:,k,i),SpatOrbs,0.0,Temp4indints(:,:),SpatOrbs)
                do l=1,SpatOrbs
                    do j=1,l
                        FourIndInts(i,j,k,l)=Temp4indints(j,l)
                        FourIndInts(i,l,k,j)=Temp4indints(j,l)
                        FourIndInts(k,j,i,l)=Temp4indints(j,l)
                        FourIndInts(k,l,i,j)=Temp4indints(j,l)
                        FourIndInts02(j,k,l,i)=Temp4indints(j,l)
                        FourIndInts02(j,i,l,k)=Temp4indints(j,l)
                        FourIndInts02(l,k,j,i)=Temp4indints(j,l)
                        FourIndInts02(l,i,j,k)=Temp4indints(j,l)
                        IF(.not.tERLocalization) THEN
                            IF(.not.((k.eq.i).or.(j.eq.l))) THEN
                                PotEnergy=PotEnergy+(Temp4indints(j,l)**2)
                                TwoEInts=TwoEInts+(Temp4indints(j,l)**2)
                                PEInts=PEInts+(Temp4indints(j,l)**2)
                            ENDIF
                        ENDIF
                   enddo
                enddo
            enddo
            IF(tERLocalization) THEN
                PotEnergy=PotEnergy+FourIndInts(i,i,i,i)
                TwoEInts=TwoEInts+FourIndInts(i,i,i,i)
                PEInts=PEInts+FourIndInts(i,i,i,i)
            ENDIF
        enddo

        IF((tROHistogram.or.tERHist).and.((Iteration.eq.1).or.(.not.tNotConverged))) CALL WriteDoubHisttofile()
        IF(tROHistogram.and.(Iteration.eq.1)) CALL WriteSingHisttofile()

! Now find the change of the potential energy due to the orthonormality of the orbitals...
        IF(tLagrange) THEN
            PEOrtho=0.D0
            do i=1,SpatOrbs
                do j=1,SpatOrbs
                    t=0.D0
                    do a=1,SpatOrbs
                        t=CoeffT1(a,i)*CoeffT1(a,j)
                    enddo
                    IF(i.eq.j) t=t-1.D0
                    PEOrtho=PEOrtho-Lambdas(i,j)*t
                    PotEnergy=PotEnergy-Lambdas(i,j)*t
                enddo
            enddo
        ENDIF

        CALL halt_timer(Transform2ElInts_Time)


    END SUBROUTINE Transform2ElInts



    SUBROUTINE WriteDoubHisttofile()
        INTEGER :: i,j,k,l,BinNo
        REAL*8 :: MaxFII,MinFII,BinIter,BinVal


        ROHistDoub(:,:)=0.D0
        MaxFII=0.D0
        MinFII=0.D0
        do l=1,SpatOrbs
            do j=1,l-1
                do k=1,SpatOrbs
                    do i=1,k-1
                        IF(FourIndInts(i,j,k,l).gt.MaxFII) MaxFII=FourIndInts(i,j,k,l)
                        IF(FourIndInts(i,j,k,l).lt.MinFII) MinFII=FourIndInts(i,j,k,l)  
                    enddo
                enddo
            enddo
        enddo
        BinIter=(MaxFII-MinFII)/2000
        BinVal=MinFII
        do i=1,2000
            ROHistDoub(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo
        do l=1,SpatOrbs
            do j=1,l-1
                do k=1,SpatOrbs
                    do i=1,k-1
                        BinNo=CEILING((FourIndInts(i,j,k,l)-MinFII)*2000/(MaxFII-MinFII))
                        ROHistDoub(2,BinNo)=ROHistDoub(2,BinNo)+1.0         
                    enddo
                enddo
            enddo
        enddo

        IF(Iteration.eq.1) THEN 
            OPEN(36,FILE='ROHistHFDoub',STATUS='unknown')
            do j=1,2000
                do i=1,2
                    WRITE(36,'(F20.10)',advance='no') ROHistDoub(i,j)
                enddo
                WRITE(36,*) ''
            enddo
            CLOSE(36)
        ENDIF
        

        IF(.not.tNotConverged) THEN
            OPEN(24,FILE='ROHistDoubRotd',STATUS='unknown')
            do j=1,2000
                do i=1,2
                    WRITE(24,'(F20.10)',advance='no') ROHistDoub(i,j)
                enddo
                WRITE(24,*) ''
            enddo
            CLOSE(24)
        ENDIF


        IF(tERHist) THEN
            ROHistER(:,:)=0.D0
            MaxFII=0.D0
            MinFII=0.D0
            do i=1,SpatOrbs
                IF(FourIndInts(i,i,i,i).gt.MaxFII) MaxFII=FourIndInts(i,i,i,i)
                IF(FourIndInts(i,i,i,i).lt.MinFII) MinFII=FourIndInts(i,i,i,i)  
            enddo
            BinIter=(MaxFII-MinFII)/2000
            BinVal=MinFII
            do i=1,2000
                ROHistER(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
            do i=1,SpatOrbs
                BinNo=CEILING((FourIndInts(i,i,i,i)-MinFII)*2000/(MaxFII-MinFII))
                ROHistER(2,BinNo)=ROHistER(2,BinNo)+1.0         
            enddo

            IF(Iteration.eq.1) THEN 
                OPEN(32,FILE='ROHistHF-ER',STATUS='unknown')
                do j=1,2000
                    do i=1,2
                        WRITE(32,'(F20.10)',advance='no') ROHistER(i,j)
                    enddo
                    WRITE(32,*) ''
                enddo
                CLOSE(32)
            ENDIF
            

            IF(.not.tNotConverged) THEN
                OPEN(34,FILE='ROHistRotd-ER',STATUS='unknown')
                do j=1,2000
                    do i=1,2
                        WRITE(34,'(F20.10)',advance='no') ROHistER(i,j)
                    enddo
                    WRITE(34,*) ''
                enddo
                CLOSE(34)
            ENDIF

        ENDIF


    ENDSUBROUTINE WriteDoubHisttofile 



    SUBROUTINE WriteSingHisttofile()
        INTEGER :: i,j,k,BinNo
        REAL*8 :: MaxFII,MinFII,BinIter,BinVal,SingExcit(SpatOrbs,SpatOrbs)


        ROHistSing(:,:)=0.D0
        MaxFII=0.D0
        MinFII=0.D0
        do j=1,SpatOrbs
            do i=1,SpatOrbs
                SingExcit(i,j)=0.D0
                IF(i.eq.j) CYCLE
                do k=1,SpatOrbs
                    IF(k.eq.j) CYCLE
                    SingExcit(i,j)=SingExcit(i,j)+REAL(TMAT2D(i,j)%v,8)+(FourIndInts(i,k,j,k)-FourIndInts(i,k,k,j))
                enddo
                IF(SingExcit(i,j).gt.MaxFII) MaxFII=SingExcit(i,j)
                IF(SingExcit(i,j).lt.MinFII) MinFII=SingExcit(i,j)
            enddo
        enddo

        BinIter=(MaxFII-MinFII)/2000
        BinVal=MinFII
        do i=1,2000
            ROHistSing(1,i)=BinVal
            BinVal=BinVal+BinIter
        enddo

        do j=1,SpatOrbs
            do i=1,SpatOrbs
                IF(i.eq.j) CYCLE
                BinNo=CEILING((SingExcit(i,j)-MinFII)*2000/(MaxFII-MinFII))
                ROHistSing(2,BinNo)=ROHistSing(2,BinNo)+1.0         
            enddo
        enddo


        IF(Iteration.eq.1) THEN
            OPEN(56,FILE='ROHistHFSing',STATUS='unknown')
            do j=1,2000
                do i=1,2
                    WRITE(56,'(F20.10)',advance='no') ROHistSing(i,j)
                enddo
                WRITE(56,*) ''
            enddo
            CLOSE(56)
        ENDIF
        
        IF(.not.tNotConverged) THEN
            OPEN(60,FILE='ROHistRotSing',STATUS='unknown')
            do j=1,2000
                do i=1,2
                    WRITE(60,'(F20.10)',advance='no') ROHistSing(i,j)
                enddo
                WRITE(60,*) ''
            enddo
            CLOSE(60)
        ENDIF



    ENDSUBROUTINE WriteSingHisttofile 




!This just tests the convergence on the grounds that the force is smaller that the input parameter: ConvergedForce
    SUBROUTINE TestForConvergence()

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



    SUBROUTINE InitLocalOrbs()
        USE SystemData , only : nBasis
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        REAL*8 :: t,RAN2
        INTEGER :: x,y,i,j,k,l,ierr,Const,iseed=-7,a,b,g,d

        IF(tERLocalization) THEN
            WRITE(6,*) "Calculating new molecular orbitals based on Edmiston-Reudenberg localisation,"
            WRITE(6,*) "i.e. maximisation of the <ii|ii> integrals..."
            WRITE(6,*) "*****"
        ELSE
            WRITE(6,*) "Calculating new molecular orbitals based on mimimisation "
            WRITE(6,*) "of <ij|kl>^2 integrals..."
            WRITE(6,*) "*****"
        ENDIF

       
! Check for possible errors.
        IF(.not.TwoCycleSymGens) THEN
            CALL Stop_All(this_routine,"ERROR. TwoCycleSymGens is false.  Symmetry is not abelian.") 
        ENDIF
!        IF((.not.lNoSymmetry).and.tSeparateOccVirt) THEN
!            CALL Stop_All(this_routine,"ERROR. This routine is only set up to maintain symmetry if all orbitals are to be mixed.") 
!        ENDIF

        IF(tLagrange) THEN
            IF(tShake) THEN
                CALL FLUSH(6)
                CALL FLUSH(12)
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
        WRITE(6,*) "*****"


        OrthoNorm=0.D0
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
        Iteration=0
        OrthoForce=0.D0


        
!Init symm arrays in here... for each of the sepoccorbs and mixall situations. 

        CALL GENSymStatePairs(SpatOrbs,.false.)
        ! if the symmetry is on, and all orbitals are being mixed, this is all that is required.

        IF(tSeparateOccVirt) THEN
            OccVirt=2
            CALL InitOrbitalSeparation()
            ! rewrite all the symmetry lists to account for the separation and have simple option if
            ! symmetry is off.
        ELSE
            OccVirt=1
            ALLOCATE(SymLabelCounts2(2,8),stat=ierr)
            CALL LogMemAlloc('SymLabelCounts2',2*8,4,this_routine,SymLabelCounts2Tag,ierr)
            SymLabelCounts2(:,:)=0                     
            IF(lNoSymmetry) THEN
                SymLabelCounts2(1,1)=1
                SymLabelCounts2(2,1)=SpatOrbs
            ELSE
                do j=1,8
                    do i=1,2
                        SymLabelCounts2(i,j)=SymLabelCounts(i,j)
                    enddo
                enddo
            ENDIF
        ENDIF

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,16
!            WRITE(6,*) SymLabelCounts2(1,i),SymLabelCounts2(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
!        do i=1,SpatOrbs
!            WRITE(6,*) SymLabelList(i),INT(G1(SymLabelList(i)*2)%sym%S,4)
!        enddo

!Set timed routine names
        Rotation_Time%timer_name='RotateTime'
        Shake_Time%timer_name='ShakeTime'
        FullShake_Time%timer_name='FullShakeTime'
        Findtheforce_Time%timer_name='FindtheForceTime'
        Transform2ElInts_Time%timer_name='Transform2ElIntsTime'
        findandusetheforce_time%timer_name='Findandusetheforce'
        CalcDerivConstr_Time%timer_name='CalcDerivConstr'
        TestOrthoConver_Time%timer_name='TestOrthoConver'
       
! Set up constraint labels.
        TotNoConstraints=(SpatOrbs*(SpatOrbs+1))/2

        ALLOCATE(Lab(2,TotNoConstraints),stat=ierr)
        CALL LogMemAlloc('Lab',2*TotNoConstraints,4,this_routine,LabTag,ierr)
        Lab(:,:)=0                     

        Const=0
        do i=1,SpatOrbs
            do j=i,SpatOrbs
                Const=Const+1
                Lab(1,Const)=i
                Lab(2,Const)=j
!                Lab2(i,j)=Const
            enddo
        enddo

!        WRITE(6,*) 'constraint labels'
!        do l=1,TotNoConstraints
!            i=Lab(1,l)
!            j=Lab(2,l)
!            WRITE(6,*) l,i,j
!        enddo
       

        WRITE(6,*) 'Total number of constraints = ',TotNoConstraints
!        WRITE(6,*) 'Const = ',Const
        IF(Const.ne.TotNoConstraints) THEN
            CALL Stop_all(this_routine,'ERROR in the number of constraints calculated.  lmax does not equal TotNoConstraints')
        ENDIF
            

!Allocate memory

        ALLOCATE(CoeffT1(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffT1',SpatOrbs**2,8,this_routine,CoeffT1Tag,ierr)
        ALLOCATE(CoeffCorT2(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffCorT2',SpatOrbs**2,8,this_routine,CoeffCorT2Tag,ierr)
        ALLOCATE(CoeffUncorT2(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffUncT2',SpatOrbs**2,8,this_routine,CoeffUncorT2Tag,ierr)
        CoeffUncorT2(:,:)=0.D0
        ALLOCATE(DerivCoeff(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('DerivCoeff',SpatOrbs**2,8,this_routine,DerivCoeffTag,ierr)

        ALLOCATE(OneIndInts01(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('OneIndInts01',SpatOrbs**4,8,this_routine,OneIndInts01Tag,ierr)
        ALLOCATE(OneIndInts02(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('OneIndInts02',SpatOrbs**4,8,this_routine,OneIndInts02Tag,ierr)

        ALLOCATE(TwoIndInts01(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('TwoIndInts01',SpatOrbs**4,8,this_routine,TwoIndInts01Tag,ierr)
        ALLOCATE(TwoIndInts02(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('TwoIndInts02',SpatOrbs**4,8,this_routine,TwoIndInts02Tag,ierr)

        ALLOCATE(ThreeIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts',SpatOrbs**4,8,this_routine,ThreeIndIntsTag,ierr)
        ALLOCATE(ThreeIndInts01(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts01',SpatOrbs**4,8,this_routine,ThreeIndInts01Tag,ierr)
        ALLOCATE(ThreeIndInts02(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts02',SpatOrbs**4,8,this_routine,ThreeIndInts02Tag,ierr)
        ALLOCATE(ThreeIndInts03(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts03',SpatOrbs**4,8,this_routine,ThreeIndInts03Tag,ierr)
        ALLOCATE(ThreeIndInts04(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts04',SpatOrbs**4,8,this_routine,ThreeIndInts04Tag,ierr)

        ALLOCATE(FourIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts',SpatOrbs**4,8,this_routine,FourIndIntsTag,ierr)
        ALLOCATE(FourIndInts02(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts02',SpatOrbs**4,8,this_routine,FourIndInts02Tag,ierr)
 

        ALLOCATE(UMATTemp01(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('UMATTemp01',SpatOrbs**4,8,this_routine,UMATTemp01Tag,ierr)
        ALLOCATE(UMATTemp02(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('UMATTemp02',SpatOrbs**4,8,this_routine,UMATTemp02Tag,ierr)

        IF(tLagrange) THEN
            ALLOCATE(Lambdas(SpatOrbs,SpatOrbs),stat=ierr)
            CALL LogMemAlloc('Lambdas',SpatOrbs**2,8,this_routine,LambdasTag,ierr)
            ALLOCATE(DerivLambda(SpatOrbs,SpatOrbs),stat=ierr)
            CALL LogMemAlloc('DerivLambda',SpatOrbs**2,8,this_routine,DerivLambdaTag,ierr)
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
            ALLOCATE(DerivConstrT1(SpatOrbs,SpatOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT1',SpatOrbs*TotNoConstraints*SpatOrbs,8,this_routine,DerivConstrT1Tag,ierr)
            DerivConstrT1(:,:,:)=0.D0
            ALLOCATE(DerivConstrT2(SpatOrbs,SpatOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT2',SpatOrbs*TotNoConstraints*SpatOrbs,8,this_routine,DerivConstrT2Tag,ierr)
            ALLOCATE(ForceCorrect(SpatOrbs,SpatOrbs),stat=ierr)
            CALL LogMemAlloc('ForceCorrect',SpatOrbs**2,8,this_routine,ForceCorrectTag,ierr)
            ALLOCATE(Correction(SpatOrbs,SpatOrbs),stat=ierr)
            CALL LogMemAlloc('Correction',SpatOrbs**2,8,this_routine,CorrectionTag,ierr)
            IF(tShakeApprox) THEN
                ALLOCATE(DerivConstrT1T2Diag(TotNoConstraints),stat=ierr)
                CALL LogMemAlloc('DerivConstrT1T2Diag',TotNoConstraints,8,this_routine,DerivConstrT1T2DiagTag,ierr)
                DerivConstrT1T2Diag(:)=0.D0
            ELSE
                ALLOCATE(DerivConstrT1T2(TotNoConstraints,TotNoConstraints),stat=ierr)
                CALL LogMemAlloc('DerivConstrT1T2',TotNoConstraints**2,8,this_routine,DerivConstrT1T2Tag,ierr)
            ENDIF
        ENDIF               

!Zero/initialise the arrays

        CoeffT1(:,:)=0.D0
        IF(lNoSymmetry) THEN
            do i=1,SpatOrbs
                do j=1,SpatOrbs
                    CoeffT1(j,i)=RAN2(iseed)*(1E-06)
                enddo
            enddo
        ENDIF
        do i=1,SpatOrbs
            CoeffT1(i,i)=1.D0
        enddo
        ! If the symmetry is kept on, start with the symmetric identity matrix of coefficients, and it will be maintained.
        
        IF(tSeparateOccVirt) CALL ZeroOccVirtElements(CoeffT1)

        CALL GRAMSCHMIDT(CoeffT1,SpatOrbs)


!        WRITE(6,*) 'coefft1'
!        do i=1,SpatOrbs
!            do j=1,SpatOrbs
!                WRITE(6,'(F15.10)',advance='no') CoeffT1(i,j)
!            enddo
!            WRITE(6,*) ''
!        enddo

        DerivCoeff(:,:)=0.D0
        UMATTemp01(:,:,:,:)=0.D0
        UMATTemp02(:,:,:,:)=0.D0

!These loops can be sped up with spatial symmetry and pairwise permutation symmetry if needed.
!Fill all the original partially transformed integral arrays with the original two-electron integrals.
        do a=1,SpatOrbs
            do g=1,a
                do b=1,SpatOrbs
                    do d=1,b
                        t=REAL(UMAT(UMatInd(a,b,g,d,0,0))%v,8)
                        IF(.not.((g.eq.a).or.(b.eq.d))) THEN
                            PotEnergy=PotEnergy+(t**2)          !Potential energy starts as this since the orbitals are orthonormal by construction.
                            TwoEInts=TwoEInts+(t**2)
                        ENDIF
                        UMATTemp01(a,g,b,d)=t                   !a,g,d,b chosen to make 'transform2elint' steps more efficient
                        UMATTemp01(g,a,b,d)=t
                        UMATTemp01(a,g,d,b)=t
                        UMATTemp01(g,a,d,b)=t

                        UMATTemp02(d,b,a,g)=t                   !d,b,a,g order also chosen to speed up the transformation.
                        UMATTemp02(d,b,g,a)=t
                        UMATTemp02(b,d,a,g)=t
                        UMATTemp02(b,d,g,a)=t
                    enddo
                enddo
            enddo
        enddo
        PEInts=PotEnergy
        CALL TestOrthonormality()


        OPEN(12,FILE='Transform',STATUS='unknown')
        IF(tLagrange) THEN
!We want to write out: Iteration, Potential energy, Force, Sum<ij|kl>^2, OrthonormalityCondition
            WRITE(12,"(A)") "# Iteration   2.PotEnergy   3.PEInts   4.PEOrtho    5.Force   6.ForceInts   7.OrthoForce    8.Sum<ij|kl>^2   9.OrthoNormCondition   10.DistMovedbyCs   11.DistMovedByLs   12.LambdaMag"
            WRITE(6,"(A)") "Iteration   2.PotEnergy   3.PEInts   4.PEOrtho   5.Force   6.ForceInts   7.OrthoForce    8.Sum<ij|kl>^2   9.OrthoNormCondition   10.DistMovedbyCs   11.DistMovedbyLs   12.LambdaMag"
        ELSE
            WRITE(12,"(A12,5A24)") "# Iteration","2.PotEnergy","3.Force","4.Totalcorrforce","5.OrthoNormCondition","6.DistMovedbyCs"
            WRITE(6,"(A12,5A24)") "Iteration","2.PotEnergy","3.Force","4.TotCorrForce","5.OrthoNormCondition","6.DistMovedbyCs"
        ENDIF

    END SUBROUTINE InitLocalOrbs

    SUBROUTINE DeallocateMem()
        CHARACTER(len=*) , PARAMETER :: this_routine='DeallocateMem'


        DEALLOCATE(Lab)
        CALL LogMemDealloc(this_routine,LabTag)
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
        
        DEALLOCATE(OneIndInts01)
        CALL LogMemDealloc(this_routine,OneIndInts01Tag)
        DEALLOCATE(OneIndInts02)
        CALL LogMemDealloc(this_routine,OneIndInts02Tag)
        DEALLOCATE(TwoIndInts01)
        CALL LogMemDealloc(this_routine,TwoIndInts01Tag)
        DEALLOCATE(TwoIndInts02)
        CALL LogMemDealloc(this_routine,TwoIndInts02Tag)
        DEALLOCATE(ThreeIndInts)
        CALL LogMemDealloc(this_routine,ThreeIndIntsTag)
        DEALLOCATE(ThreeIndInts01)
        CALL LogMemDealloc(this_routine,ThreeIndInts01Tag)
        DEALLOCATE(ThreeIndInts02)
        CALL LogMemDealloc(this_routine,ThreeIndInts02Tag)
        DEALLOCATE(ThreeIndInts03)
        CALL LogMemDealloc(this_routine,ThreeIndInts03Tag)
        DEALLOCATE(ThreeIndInts04)
        CALL LogMemDealloc(this_routine,ThreeIndInts04Tag)
        DEALLOCATE(FourIndInts)
        CALL LogMemDealloc(this_routine,FourIndIntsTag)
        DEALLOCATE(FourIndInts02)
        CALL LogMemDealloc(this_routine,FourIndInts02Tag)
        DEALLOCATE(UMATTemp01)
        CALL LogMemDealloc(this_routine,UMATTemp01Tag)
        DEALLOCATE(UMATTemp02)
        CALL LogMemDealloc(this_routine,UMATTemp02Tag)


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
    


    SUBROUTINE WriteStats()

        IF(tLagrange) THEN
            WRITE(6,"(I7,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
            WRITE(12,"(I7,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
        ELSE
            IF(Mod(Iteration,10).eq.0) THEN
                WRITE(6,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
                WRITE(12,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
            ENDIF
        ENDIF
        CALL FLUSH(6)
        CALL FLUSH(12)

! after writing out stats, test for softexit.
        if (test_SOFTEXIT()) then
            WRITE(6,*) 'SOFTEXIT detected, finalizing new orbitals.'
            tNotConverged=.false.
        endif


    END SUBROUTINE WriteStats


    SUBROUTINE InitSymmArrays()
        INTEGER :: i,ierr,SymSum
        CHARACTER(len=*) , PARAMETER :: this_routine='InitSymmArrays'


        ! SymLabelCounts2(1:Sym) is the index in SymLabelList where the symmetry block S starts
        ! SymLabelCounts2(2:Sym) is the number of orbitals in symmetry block S.
        ! E.g. Symmetry S may start at index 2 and have 3 orbitals.
        ! SymLabelList(2)->SymLabelList(4) will give the indexes of these orbitals.

        WRITE(6,*) 'Sym Label Counts'
        do i=1,8
            WRITE(6,*) SymLabelCounts2(1,i),SymLabelCounts2(2,i)
        enddo

        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
        do i=1,SpatOrbs
            WRITE(6,*) SymLabelList(i),INT(G1(SymLabelList(i)*2)%sym%S,4)
        enddo
        stop

    ENDSUBROUTINE InitSymmArrays


    SUBROUTINE InitOrbitalSeparation()
! This subroutine is called if the SEPARATEOCCVIRT keyword is present in the input.
! This means that two iterations of the rotate orbs routine will be performed, one in which the occupied orbtials are localised,
! and one in which the virtual orbitals are occupied.
! First need to reorder orbitals so that the first NEl are the occupied HF (ordered by sym then energy).
        INTEGER :: i,j,ierr,t,SymCurr,Symi
        CHARACTER(len=*) , PARAMETER :: this_routine='InitOrbitalSeparation'

 
! Kind of pick out the NEl lowest energy, reorder them, then reorder the virtual ones.

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index with the second lowest energy.

!        do i=1,nBasis
!            WRITE(6,*) BRR(i)
!        enddo

!        DEALLOCATE(SymLabelCounts)

        ALLOCATE(SymLabelCounts2(2,16),stat=ierr)
        CALL LogMemAlloc('SymLabelCounts2',2*16,4,this_routine,SymLabelCounts2Tag,ierr)
        SymLabelCounts2(:,:)=0
        ! first 8 refer to the occupied, and the second to the virtual.

        SymLabelList(:)=0
        ALLOCATE(LabVirtOrbs((nBasis-NEl)/2),stat=ierr)
        CALL LogMemAlloc('LabVirtOrbs',((nBasis-NEl)/2),4,this_routine,LabVirtOrbsTag,ierr)
        LabVirtOrbs(:)=0
        ALLOCATE(LabOccOrbs(NEl/2),stat=ierr)
        CALL LogMemAlloc('LabOccOrbs',(NEl/2),4,this_routine,LabOccOrbsTag,ierr)
        LabOccOrbs(:)=0

! this next part picks out the NEl/2 lowest energy orbitals from BRR as these will be the occupied.
! these occcupied orbitals are then ordered according to symmetry, and the same done to the virual.
        do i=1,NEl/2
            LabOccOrbs(i)=(BRR(2*i))/2
        enddo
        CALL NECI_SORTI(NEl/2,LabOccOrbs)
        ! Sorts LabOrbs numerically and therefore into symmetry

        do i=1,(nBasis-NEl)/2
            LabVirtOrbs(i)=(BRR((2*i)+NEl))/2
        enddo
        CALL NECI_SORTI((nBasis-NEl)/2,LabVirtOrbs)

        do i=1,NEl/2
            SymLabelList(i)=LabOccOrbs(i)
        enddo
        j=0
        do i=(NEl/2)+1,SpatOrbs
            j=j+1
            SymLabelList(i)=LabVirtOrbs(j)
        enddo
        ! a new symlabellist is now created with the occupied then virtual orbitals separately ordered
        ! according to symmetry.


        IF(lNoSymmetry) THEN
            SymLabelCounts2(1,1)=1
            SymLabelCounts2(1,9)=(NEl/2)+1
            SymLabelCounts2(2,1)=NEl/2
            SymLabelCounts2(2,9)=(nBasis-NEl)/2
            ! basically filling these to say that all spatial orbitals have 0 symmetry.
            ! SymLabelList will already have all 0's in it.
        ELSE 
            SymCurr=0
            SymLabelCounts2(1,1)=1
            do i=1,NEl/2
                Symi=INT(G1(SymLabelList(i)*2)%sym%S,4)
                SymLabelCounts2(2,(Symi+1))=SymLabelCounts2(2,(Symi+1))+1
                IF(Symi.gt.SymCurr) THEN
                    SymLabelCounts2(1,(Symi+1))=i
                    SymCurr=Symi
                ENDIF
            enddo
            SymCurr=0
            SymLabelCounts2(1,9)=(NEl/2)+1
            do i=(NEl/2)+1,SpatOrbs
                Symi=INT(G1(SymLabelList(i)*2)%sym%S,4)
                SymLabelCounts2(2,(Symi+9))=SymLabelCounts2(2,(Symi+9))+1
                IF(Symi.gt.SymCurr) THEN
                    SymLabelCounts2(1,(Symi+9))=i
                    SymCurr=Symi
                ENDIF
            enddo
        ENDIF

        DEALLOCATE(LabOccOrbs)
        CALL LogMemDealloc(this_routine,LabOccOrbsTag)
        DEALLOCATE(LabVirtOrbs)
        CALL LogMemDealloc(this_routine,LabVirtOrbsTag)

!        WRITE(6,*) 'labocc'
!        do i=1,NEl/2
!            WRITE(6,*) LabOccOrbs(i)
!        enddo
!        WRITE(6,*) 'labvirt'
!        do i=1,(nBasis-NEl)/2
!            WRITE(6,*) LabVirtOrbs(i)
!        enddo
!        WRITE(6,*) 'symlabellist'
!        do i=1,SpatOrbs
!            WRITE(6,*) SymLabelList(i)
!        enddo


! After this routine, the array SymLabelList contains the labels of the occupied or virtual orbitals
! in order of symmetry then energy.


    ENDSUBROUTINE InitOrbitalSeparation




    SUBROUTINE ZeroOccVirtElements(Coeff)
! This routine sets all the elements of the coefficient matrix that connect occupied and virtual orbitals to 0.
! This ensures that only occupied mix with occupied and virtual mix with virtual.
! It is not the most efficient method but is o.k for now.
        REAL*8 :: Coeff(SpatOrbs,SpatOrbs)
        INTEGER :: i,j


        do i=1,NEl/2
            do j=NEl/2+1,SpatOrbs
                Coeff(SymLabelList(i),SymLabelList(j))=0.D0
                Coeff(SymLabelList(j),SymLabelList(i))=0.D0
            enddo
        enddo

!        WRITE(6,*) 'labels virt occ'
!        do i=1,NEl/2
!            WRITE(6,*) i,LabOccOrbs(i),AOSymm(LabOccOrbs(i))
!        enddo
!        do j=1,((nBasis-NEl)/2)
!            WRITE(6,*) j,LabVirtOrbs(j),AOSymm(LabVirtOrbs(j))
!        enddo

!        stop

! LabOccOrbs are the orbital indices of the occupied orbitals.



    ENDSUBROUTINE ZeroOccVirtElements




! DerivCoeff(k,a) is the unconstrained force on the original coefficients (CoeffT1(a,k)). 
    SUBROUTINE ShakeConstraints()
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
!                do m=1,SpatOrbs
!                    do a=1,SpatOrbs
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
!        IF(Mod(Iteration,10).eq.0.or.Iteration.eq.1) CALL WriteShakeOUTstats01()

        CALL set_timer(Shake_Time,30)

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
                

            do w=1,OccVirt
                IF(w.eq.1) THEN
                    SymMin=1
                    MinMZ=1
                    IF(tSeparateOccVirt) THEN
                        MaxMZ=NEl/2
                    ELSE
                        MaxMZ=SpatOrbs
                    ENDIF
                ELSE
                    SymMin=9
                    MinMZ=(NEl/2)+1
                    MaxMZ=SpatOrbs
                ENDIF

                do x=MinMZ,MaxMZ
                    m=SymLabelList(x)
                    SymM=INT(G1(m*2)%sym%S,4)
                    do y=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
                        a=SymLabelList(y)
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
!            IF(Mod(Iteration,10).eq.0.or.Iteration.eq.1) CALL WriteShakeOUTstats02(ShakeIteration,TotLambdas,ConvergeCount) ! add correction to this.


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



    SUBROUTINE WriteShakeOUTstats01()
        INTEGER :: i,j,x,y,l,m,a
            

            WRITE(6,*) 'Original coefficients'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    WRITE(6,'(4F20.10)',advance='no') CoeffT1(a,m)
                enddo
                WRITE(6,*) ''
            enddo
            
            WRITE(6,*) 'Uncorrected Force'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    WRITE(6,'(4F20.10)',advance='no') DerivCoeff(a,m)
                enddo
                WRITE(6,*) ''
            enddo


            WRITE(6,*) 'Coefficients at t2, having been shifted by the uncorrected force'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
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
        INTEGER :: m,a,l,ConvergeCount,ShakeIteration
        REAL*8 :: TotLambdas 
    
            WRITE(6,*) 'Iteration number ,', ShakeIteration

            WRITE(6,*) 'Lambdas used for this iteration'
            do l=1,TotNoConstraints
                WRITE(6,'(10F20.10)') ShakeLambda(l)
                TotLambdas=TotLambdas+ShakeLambda(l)
            enddo

            WRITE(6,*) 'Corrected Force'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    WRITE(6,'(4F20.10)',advance='no') ForceCorrect(a,m)
                enddo
                WRITE(6,*) ''
            enddo
    
            WRITE(6,*) 'Coefficients having been shifted by the corrected force (at t2)'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
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
!            do m=1,SpatOrbs
!                WRITE(6,*) 'm equals ', m 
!                do l=1,TotNoConstraints
!                    WRITE(6,*) 'i, j equals ',Lab(1,l),Lab(2,l)
!                    do a=1,SpatOrbs
!                        WRITE(6,'(4F20.10)',advance='no') DerivConstrT1(a,m,l)
!                    enddo
!                enddo
!            enddo

!            WRITE(6,*) 'Derivative of constraints w.r.t shifted coefficients (t2)'
!            do m=1,SpatOrbs
!                WRITE(6,*) 'm equals ', m  
!                do l=1,TotNoConstraints
!                    WRITE(6,*) 'i,j equals ',Lab(1,l),Lab(2,l)
!                    do a=1,SpatOrbs
!                        WRITE(6,'(4F20.10)',advance='no') DerivConstrT2(a,m,l)
!                    enddo
!                enddo
!            enddo


    ENDSUBROUTINE WriteShakeOUTstats02



    SUBROUTINE CalcDerivConstr(CurrCoeff,DerivConstr)
    ! This calculates the derivative of each of the orthonormalisation constraints, l, with respect
    ! to each set of coefficients cm.

        INTEGER :: l,m,i,j,a,Symm
        REAL*8 :: CurrCoeff(SpatOrbs,SpatOrbs)
        REAL*8 :: DerivConstr(SpatOrbs,SpatOrbs,TotNoConstraints)

        call set_timer(CalcDerivConstr_Time,30)


            DerivConstr(:,:,:)=0.D0
            do l=1,TotNoConstraints
                i=lab(1,l)
                j=lab(2,l)
                IF(i.eq.j) THEN
                    do a=1,SpatOrbs
                        DerivConstr(a,i,l)=CurrCoeff(a,i)*2
                    enddo
                ELSE
                    do a=1,SpatOrbs
                        DerivConstr(a,j,l)=CurrCoeff(a,i) 
                    enddo
                    do a=1,SpatOrbs
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
!                do m=1,SpatOrbs
!                    do a=1,SpatOrbs
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
        INTEGER :: a,l,m,Symm,x,y,w,SymMin
        REAL*8 :: TotForce,TotDiffCoeffs,CoeffT2(SpatOrbs,SpatOrbs)

!        WRITE(6,*) 'DerivCoeff'
!        do m=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(4F20.10)',advance='no') DerivCoeff(a,m)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'Correction'
!        do m=1,SpatOrbs
!            do a=1,SpatOrbs
!                WRITE(6,'(4F20.10)',advance='no') Correction(a,m)
!            enddo
!            WRITE(6,*) ''
!        enddo


        CALL set_timer(findandusetheforce_time,30)


        do w=1,OccVirt
            IF(w.eq.1) THEN
                SymMin=1
                MinMZ=1
                IF(tSeparateOccVirt) THEN
                    MaxMZ=NEl/2
                ELSE
                    MaxMZ=SpatOrbs
                ENDIF
            ELSE
                SymMin=9
                MinMZ=(NEl/2)+1
                MaxMZ=SpatOrbs
            ENDIF

            do x=MinMZ,MaxMZ
                m=SymLabelList(x)
! Symmetry requirement that z must be from the same irrep as m
                SymM=INT(G1(m*2)%sym%S,4)
                do y=SymLabelCounts2(1,SymM+SymMin),(SymLabelCounts2(1,SymM+SymMin)+SymLabelCounts2(2,SymM+SymMin)-1)
                    a=SymLabelList(y)
                ! FIND THE FORCE 
                    IF(tERLocalization) THEN
                        ForceCorrect(a,m)=DerivCoeff(a,m)+Correction(a,m)
                    ! find the corrected force. (in the case where the uncorrected force is required, correction is set to 0.
                    ! DerivCoeff(m,a) is the derivative of |<ij|kl>|^2 w.r.t cm without any constraints (no lambda terms).
                    ! ForceCorrect is then the latest force on coefficients.  This is iteratively being corrected so that
                    ! it will finally move the coefficients so that they remain orthonormal.
                
                ! USE THE FORCE
                        CoeffT2(a,m)=CoeffT1(a,m)+(TimeStep*ForceCorrect(a,m))
                    ELSE
                        ForceCorrect(a,m)=DerivCoeff(a,m)-Correction(a,m)
                        CoeffT2(a,m)=CoeffT1(a,m)-(TimeStep*ForceCorrect(a,m))
                    ENDIF
                    ! Using the force to calculate the coefficients at time T2 (hopefully more orthonomal than those calculated in
                    ! the previous iteration).
                    
                    ! Calculate parameters for printing
                    TotForce=TotForce+ABS(ForceCorrect(a,m))
                    TotDiffCoeffs=TotDiffCoeffs+ABS(CoeffT2(a,m)-CoeffT1(a,m))
                enddo
            enddo
        enddo

!            IF(tSeparateOccVirt) THEN
!                CALL ZeroOccVirtElements(CoeffT2)
!                do m=1,SpatOrbs
!                    Symm=AOSymm(m) 
!                    do a=(SymInd(Symm+1)),(SymInd(Symm+1)+MBas(Symm+1)-1)
!                        TotForce=TotForce+ForceCorrect(a,m)
!                        TotDiffCoeffs=TotDiffCoeffs+ABS(CoeffT2(a,m)-CoeffT1(a,m))
!                    enddo
!                enddo
!            ENDIF

!        WRITE(6,*) 'ForceCorrect'
!        do m=1,SpatOrbs
!            do a=1,SpatOrbs
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
        REAL*8 :: CurrCoeff(SpatOrbs,SpatOrbs),TotConstraints,Constraint(TotNoConstraints) 


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
     
!            do i=1,SpatOrbs
!                do j=1,SpatOrbs
!                    WRITE(6,'(I3)') i,j
!                    WRITE(6,'(F20.10)') CurrCoeff(i,j)
!                enddo
!                WRITE(6,*) ''
!            enddo
  
!            do l=1,TotNoConstraints
!                WRITE(6,*) Constraint(l)
!            enddo


    ENDSUBROUTINE CalcConstraints


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
                ! Sum of all Contraint componenets - indication of overall orthonormality.
        
                IF(ABS(ConstraintCor(l)).gt.ShakeConverged) ConvergeCount=ConvergeCount+1
                ! Count the number of constraints which are still well above 0.
                
            enddo
            
            IF(tShakeIter) THEN
                IF(ShakeIteration.eq.ShakeIterMax) THEN
                    do m=1,SpatOrbs
                        do a=1,SpatOrbs
                            CoeffT1(a,m)=CoeffCorT2(a,m)
                        enddo
                    enddo
!                    WRITE(6,*) 'stopped at iteration, ',ShakeIteration
                    tShakeNotConverged=.false.
                ENDIF
            ELSEIF(ConvergeCount.eq.0) THEN
               tShakeNotConverged=.false.
!                WRITE(6,*) 'Convergence reached in the shake algorithm'
!                WRITE(6,*) 'All constraints have values less than ',ShakeConverged

! If convergence is reached, make the new coefficients coeff, to start the rotation iteration again.

                do m=1,SpatOrbs
                    do a=1,SpatOrbs
                        CoeffT1(a,m)=CoeffCorT2(a,m)
                    enddo
                enddo
            ENDIF


    ENDSUBROUTINE TestShakeConvergence



    SUBROUTINE ShakeApproximation()
    ! This is an approximation in which only the diagonal elements are considered in the 
    ! matrix of the derivative of the constraints DerivConstrT1T2.
        INTEGER :: m,l,ierr
        CHARACTER(len=*), PARAMETER :: this_routine='ShakeApproximation'


! Use 'shake' algorithm in which the iterative scheme is applied to each constraint in succession.
            WRITE(6,*) 'DerivConstrT1T2Diag calculated from the shake approx'
            
            DerivConstrT1T2Diag(:)=0.D0
            do l=1,TotNoConstraints 
                do m=1,SpatOrbs
                    DerivConstrT1T2Diag(l)=DerivConstrT1T2Diag(l)+Dot_Product(DerivConstrT2(:,m,l),DerivConstrT1(:,m,l))
                enddo
                ShakeLambdaNew(l)=Constraint(l)/((-1)*TimeStep*DerivConstrT1T2Diag(l))
                WRITE(6,*) DerivConstrT1T2Diag(l)
            enddo


    ENDSUBROUTINE ShakeApproximation




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
                    do m=1,SpatOrbs
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
   

    SUBROUTINE FinalizeNewOrbs()
! At the end of the orbital rotation, have a set of coefficients CoeffT1 which transform the HF orbitals into a set of linear
! combinations ui which minimise |<ij|kl>|^2.  This is the final subroutine after all iterations (but before the memory deallocation)
! that calculates the final 4 index integrals to be used in the NECI calculation.
        INTEGER :: x,y,i,a
        REAL*8 :: TotGSConstraints,GSConstraint(TotNoConstraints)
        
        WRITE(6,*) 'The final transformation coefficients before gram schmidt orthonormalisation'
        do i=1,SpatOrbs
            do a=1,SpatOrbs
                WRITE(6,'(F20.10)',advance='no') CoeffT1(a,i)
            enddo
            WRITE(6,*) ''
        enddo

        WRITE(6,*) 'The final values of the constraints with corrected coefficients'
        do i=1,TotNoConstraints
            WRITE(6,*) ConstraintCor(i)
        enddo

! First need to do a final explicit orthonormalisation.  The orbitals are very close to being orthonormal, but not exactly.
! Need to make sure they are exact orthonormal using Gram Schmit.
        CALL GRAMSCHMIDT(CoeffT1,SpatOrbs)

! Write out some final results of interest, like values of the constraints, values of new coefficients.

        WRITE(6,*) 'The final transformation coefficients after gram schmidt orthonormalisation'
        do i=1,SpatOrbs
            do a=1,SpatOrbs
                WRITE(6,'(F20.10)',advance='no') CoeffT1(a,i)
            enddo
            WRITE(6,*) ''
        enddo

! This file is printed to be used to produce cube files from QChem.
! Line 1 is the coefficients of HF spatial orbitals 1 2 3 ... which form transformed orbital 1 etc.

        OPEN(66,FILE='MOTRANSFORM',FORM='UNFORMATTED')

        do a=1,SpatOrbs
            do i=1,SpatOrbs
                WRITE(66) CoeffT1(a,i),CoeffT1(a,i)
            enddo
!            WRITE(66) ''
            do i=1,SpatOrbs
               WRITE(66) CoeffT1(a,i),CoeffT1(a,i)
            enddo
!            WRITE(66) ''
        enddo

        CLOSE(66)
      
 
!        OPEN(67,FILE='MOTRANSFORM02',status='unknown')

!        do a=1,SpatOrbs
!            do i=1,SpatOrbs
!                WRITE(67,'(2F15.10)',advance='no') CoeffT1(a,i),CoeffT1(a,i)
!            enddo
!            WRITE(67,*) ''
!            do i=1,SpatOrbs
!               WRITE(67,'(2F15.10)',advance='no') CoeffT1(a,i),CoeffT1(a,i)
!            enddo
!            WRITE(67,*) ''
!        enddo

!        CLOSE(67)
        


        CALL CalcConstraints(CoeffT1,GSConstraint,TotGSConstraints)  


        WRITE(6,*) 'The values of the constraints after gram schmidt orthonormalisation'
        
        do i=1,TotNoConstraints
            WRITE(6,*) GSConstraint(i)
        enddo

        
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



    SUBROUTINE RefillUMATandTMAT2D()
        INTEGER :: l,k,j,i,a,b,BinNo
        REAL*8 :: NewTMAT,TMAT2DPart(nBasis,nBasis)


! Make the UMAT elements the four index integrals.  These are calculated by transforming the HF orbitals using the
! coefficients that have been found
        do l=1,SpatOrbs
            do k=1,SpatOrbs
                do j=1,SpatOrbs
                    do i=1,SpatOrbs
                        UMAT(UMatInd(i,j,k,l,0,0))=HElement(FourIndInts(i,j,k,l))
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
            do k=1,SpatOrbs
                NewTMAT=0.D0
                do b=1,SpatOrbs
                    NewTMAT=NewTMAT+(CoeffT1(b,k)*REAL(TMAT2D(2*b,a)%v,8))
                enddo
                TMAT2DPart(2*k,a)=NewTMAT
                TMAT2DPart(2*k-1,a)=NewTMAT
            enddo
        enddo
        do k=1,nBasis
            do l=1,SpatOrbs
                NewTMAT=0.D0
                do a=1,SpatOrbs
                    NewTMAT=NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,2*a))
                enddo
                TMAT2D(k,2*l)=HElement(NewTMAT)
                TMAT2D(k,2*l-1)=HElement(NewTMAT)
            enddo
        enddo

    
!        WRITE(6,*) 'TMAT2D after transformation'
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l)%v,8)
!            enddo
!            WRITE(6,*) ''
!        enddo


        IF(tROHistogram) CALL WriteSingHisttofile()

        IF(tROFciDump) CALL PrintROFCIDUMP()


    ENDSUBROUTINE RefillUMATandTMAT2D


    

    SUBROUTINE PrintROFCIDUMP()
!This prints out a new FCIDUMP file in the same format as the old one.
        INTEGER :: i,j,k,l,Sym


        OPEN(48,FILE='ROFCIDUMP',STATUS='unknown')
        

        WRITE(48,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB= ',SpatOrbs,',NELEC=',NEl,',MS2=',LMS,','
        WRITE(48,'(A9)',advance='no') 'ORBSYM='
        i=1
        do while (i.le.SpatOrbs)
            WRITE(48,'(A2)',advance='no') '1,'
            i=i+1
        enddo
        WRITE(48,*) ''
        WRITE(48,'(A7,I1)') 'ISYM=',1
        WRITE(48,'(A5)') '&END'
       
        do i=1,SpatOrbs
            do k=1,i
                do j=1,SpatOrbs
                    ! Potential to put symmetry in here.
                    do l=1,j
                        WRITE(48,'(F21.12,4I3)') REAL(UMat(UMatInd(i,j,k,l,0,0))%v,8),i,k,j,l 
                    enddo
                enddo
            enddo
        enddo

! TMAT2D stored as spin orbitals
        do k=1,SpatOrbs
            ! Symmetry?
            do i=k,SpatOrbs
                WRITE(48,'(F21.12,4I3)') REAL(TMAT2D(2*i,2*k)%v,8),i,k,0,0
            enddo
        enddo

! ARR has the energies of the orbitals (eigenvalues).  ARR(:,2) has ordering we want.
! ARR is stored as spin orbitals.

        do k=1,SpatOrbs
            WRITE(48,'(F21.12,4I3)') Arr(2*k,2),k,0,0,0
        enddo

        WRITE(48,'(F21.12,4I3)') ECore,0,0,0,0
        
        CLOSE(48)


    ENDSUBROUTINE PrintROFCIDUMP




    SUBROUTINE CalcFOCKMatrix()
        USE SystemData , only : nBasis
        INTEGER :: i,j,a,ierr
        REAL*8 :: FOCKDiagSumHF,FOCKDiagSumNew,ArrTemp(nBasis)
        CHARACTER(len=*) , PARAMETER :: this_routine='CalcFOCKMatrix'

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

    
        ALLOCATE(ArrNew(nBasis,nBasis),stat=ierr)
        CALL LogMemAlloc('ArrNew',nBasis**2,8,this_routine,ArrNewTag,ierr)
        ArrNew(:,:)=0.D0                     

        WRITE(6,*) 'The diagonal fock elements in the HF basis set'
        do a=1,nBasis
            WRITE(6,'(F20.10)',advance='no') Arr(a,2)
        enddo


! First calculate the sum of the diagonal elements, ARR.
! Check if this is already being done.
        FOCKDiagSumHF=0.D0
        do a=1,nBasis        
            FOCKDiagSumHF=FOCKDiagSumHF+Arr(a,2)
        enddo

        WRITE(6,*) 'Sum of the diagonal elements of the fock matrix in the HF basis set = ',FOCKDiagSumHF

!        WRITE(6,*) 'Coeffs'
!        do i=1,SpatOrbs
!            do j=1,SpatOrbs
!                WRITE(6,'(F20.10)',advance='no') CoeffT1(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

! Then calculate the fock matrix in the transformed basis, and the sum of the new diagonal elements.
! Our Arr in spin orbitals.
!        do j=1,SpatOrbs
!            ArrNew(j,j)=Arr(2*j,2)
!        enddo

        FOCKDiagSumNew=0.D0
        do j=1,SpatOrbs
            do i=1,SpatOrbs
                ArrNew(i,j)=0.D0
                do a=1,SpatOrbs
                    ArrNew(i,j)=ArrNew(i,j)+(CoeffT1(a,i)*Arr(2*a,2)*CoeffT1(a,j))
                enddo
!                ArrNew(2*i-1,2*j-1)=ArrNewTemp
!                ArrNew(2*i,2*j)=ArrNewTemp
            enddo
            FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(j,j)*2)
            !only running through spat orbitals, count each twice to compare to above.
        enddo
        
!        do j=1,SpatOrbs
!            FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(j,j)*2)
!        enddo

        WRITE(6,*) 'Sum of the diagonal elements of the fock matrix in the transformed basis set = ',FOCKDiagSumNew

        WRITE(6,*) 'The fock matrix for the transformed orbitals'
        do j=1,SpatOrbs
            do i=1,SpatOrbs
                WRITE(6,'(F20.10)',advance='no') ArrNew(i,j)
            enddo
            WRITE(6,*) ''
        enddo

!        WRITE(6,*) 'BRR then ARR before being changed'
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
       

! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2) (ordered in terms of orbital number).


        do j=1,SpatOrbs
            ARR(2*j,2)=ArrNew(j,j)
            ARR(2*j-1,2)=ArrNew(j,j)
            ARRTemp(2*j)=ArrNew(BRR(2*j)/2,BRR(2*j)/2)
            ARRTemp(2*j-1)=ArrNew(BRR(2*j)/2,BRR(2*j)/2)
        enddo
! fill both ARR(:,1) and ARR(:,2) with the diagonal values, then sort ARR(:,1) according to the values (energies), 
! taking the orbital labels (BRR) with them.

!        WRITE(6,*) 'ARRtemp'
!        do i=1,nBasis
!            WRITE(6,*) BRR(i),ARRTemp(i)
!        enddo

        CALL NECI_SORT2(nBasis,ARRTemp,BRR)

        do j=1,nBasis
            ARR(j,1)=ArrTemp(j)
        enddo

!        WRITE(6,*) 'BRR then ARR after being changed'
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1)
!        enddo
       
        CALL FLUSH(6)

        DEALLOCATE(ArrNew)
        CALL LogMemDealloc(this_routine,ArrNewTag)
       

    ENDSUBROUTINE CalcFOCKMatrix



END MODULE RotateOrbsMod
