!TO DO
!Spatial symmetry
!Finalize orbs to fully test
!Permutational symmetry where possible
!Parallelize

MODULE RotateOrbsMod

    USE Global_utilities
    USE IntegralsData , only : UMAT
    USE UMatCache , only : UMatInd
    USE HElem , only : HElement
    USE SystemData , only : ConvergedForce,TimeStep,tLagrange,tShake,tShakeApprox,ShakeConverged,tROIteration,ROIterMax,tShakeIter,ShakeIterMax
    USE SystemData , only : G1,ARR
    USE OneEInts , only : TMAT2D
    USE SymData , only : TwoCycleSymGens
    USE Timing , only : end_timing,print_timing_report
    IMPLICIT NONE
    INTEGER , ALLOCATABLE :: Lab(:,:),AOSymm(:)
    REAL*8 , ALLOCATABLE :: CoeffT1(:,:),CoeffCorT2(:,:),CoeffUncorT2(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:),ForceCorrect(:,:),Correction(:,:),ShakeLambdaNew(:),ConstraintCor(:)
    REAL*8 , ALLOCATABLE :: Constraint(:),ShakeLambda(:),DerivConstrT1(:,:,:),DerivConstrT2(:,:,:),DerivConstrT1T2(:,:),DerivConstrT1T2Diag(:)
    REAL*8 , ALLOCATABLE :: OneIndInts(:,:,:,:), TwoIndInts(:,:,:,:), ThreeIndInts(:,:,:,:), FourIndInts(:,:,:,:)   !These are the arrays to store the
!partially transformed integrals for each iteration. The FourIndInts are the full <ij|kl> integrals for the current iteration.
!OneIndInts = <i b|g d> ; TwoIndInts = <i b|k d> ; ThreeIndInts = <i j|k d>
    INTEGER :: SpatOrbs,OneIndIntsTag,TwoIndIntsTag,ThreeIndIntsTag,FourIndIntsTag,LabTag,ForceCorrectTag,CorrectionTag,SymInd(8),AOSymmTag,MBas(8)
    INTEGER :: CoeffT1Tag,CoeffCorT2Tag,CoeffUncorT2Tag,LambdasTag,DerivCoeffTag,DerivLambdaTag,Iteration,TotNoConstraints,ShakeLambdaNewTag
    INTEGER :: ShakeLambdaTag,ConstraintTag,ConstraintCorTag,DerivConstrT1Tag,DerivConstrT2Tag,DerivConstrT1T2Tag,DerivConstrT1T2DiagTag
    LOGICAL :: tNotConverged
    REAL*8 :: OrthoNorm,PotEnergy,Force,TwoEInts,DistCs,OrthoForce,DistLs,LambdaMag,PEInts,PEOrtho,ForceInts,TotCorrectedForce
    REAL*8 :: OrthoFac=1.D0

    TYPE(timer), save :: Rotation_Time,Shake_Time,Findtheforce_Time

    contains

    SUBROUTINE RotateOrbs()


!Sym of spatial orbital i is INT(G1(2*i)%sym%S,4) - should go from 0 -> 7

        CALL InitLocalOrbs()        ! Set defaults, allocate arrays, write out headings for OUTPUT, set integarals to HF values.

        CALL WriteStats()           ! write out the original stats before any rotation.

        CALL set_timer(Rotation_Time,30)

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

!option to create new FCIDUMP?

        CALL DeallocateMem()

        CLOSE(12)
        
        CALL FLUSH(6)
        CALL FLUSH(12)
        CALL end_timing()
        CALL print_timing_report()

        CALL Stop_All("RotateOrbs","This code is still in the testing phase")

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

        CALL TestOrthonormality()
!Force should go to zero as we end in minimum - test for this
        CALL TestForConvergence()

    END SUBROUTINE FindNewOrbs

    
    
    SUBROUTINE TestOrthonormality()
        INTEGER :: i,j,a

        OrthoNorm=0.D0
        do i=1,SpatOrbs
            do j=1,i
                do a=1,SpatOrbs
                    OrthoNorm=OrthoNorm+CoeffT1(a,i)*CoeffT1(a,j)
                enddo
            enddo
        enddo
        OrthoNorm=OrthoNorm-real(SpatOrbs,8)
        OrthoNorm=(OrthoNorm*2.D0)/REAL((SpatOrbs*(SpatOrbs+1.D0)),8)

    END SUBROUTINE TestOrthonormality




    SUBROUTINE FindTheForce
        INTEGER :: m,z,i,j,k,l,a,b,g,d,Symm,Symb,Symb2,Symi,Symd
        REAL*8 :: Deriv4indint,Deriv4indintsqrd,t1,t2,t3,t4,LambdaTerm1,LambdaTerm2
        CHARACTER(len=*) , PARAMETER :: this_routine='FindtheForce'
      
! Running over m and z, covers all matrix elements of the force matrix (derivative 
! of equation we are minimising, with respect to each translation coefficient) filling 
! them in as it goes.
        CALL set_timer(FindtheForce_time,30)

        DerivCoeff(:,:)=0.D0
        Force=0.D0
        Deriv4indintsqrd=0.D0
        ForceInts=0.D0
        OrthoForce=0.D0
       

        
        do m=1,SpatOrbs
! To include: symmetry requirement that z must be from the same irrep as m
            Symm=AOSymm(m)                                             
            ! Gives the symmetry of the orbital m (0 -> 7)
            do z=(SymInd(Symm+1)),(SymInd(Symm+1)+MBas(Symm+1)-1)          
            ! Runs only over those z with the same symm. 
            
                Deriv4indintsqrd=0.D0
                    
! This runs over the full <ij|kl> integrals from the previous iteration. 
! In the future we can take advantage of the permutational symmetry of the matrix elements
                do l=1,SpatOrbs
                    do j=1,l-1                          ! j < l
                        do k=1,SpatOrbs                 
                            Symi=IEOR(AOSymm(k),IEOR(AOSymm(j),AOSymm(l)))
                            ! only i with symmetry equal to j x k x l will have integrals with overall
                            ! symmetry A1 and therefore be non-zero.

                            do i=(SymInd(Symi+1)),(SymInd(Symi+1)+MBas(Symi+1)-1)
                                IF(i.lt.k) THEN         ! i < k                            

!                                IF((((i-1)*SpatOrbs)+k).ge.(((j-1)*SpatOrbs)+l)) THEN
                                   

! Already calculated the four index integrals and some partial integrals
! in Transform2ElInts routine of the previous iteration.

! Deriv4indint is the derivative of <ij|kl> with respect to a coefficient
! c(m,z), at the current m,z,i,j,k,l of the loop.
! This involves 4 terms in which the derivative is taken w.r.t the i,j,k or l position in <ij|kl>. 
                                t1=0.D0
                                t2=0.D0
                                t3=0.D0
                                t4=0.D0

! Alternative method, but actually appears to take longer...                                
!                                IF (m.eq.l) t4=ThreeIndInts(i,j,k,z)
!                                do d=1,SpatOrbs
!                                    IF(m.eq.j) t2=t2+(CoeffT1(d,l)*TwoIndInts(i,z,k,d))
!                                    do b=1,SpatOrbs
!                                        IF(m.eq.k) t3=t3+(CoeffT1(b,j)*CoeffT1(d,l)*OneIndInts(i,b,z,d))
!                                        do g=1,SpatOrbs
!                                            IF(m.eq.i) t1=t1+(Coefft1(b,j)*Coefft1(g,k)*Coefft1(d,l)*real(Umat(Umatind(z,b,g,d,0,0))%v,8))
!                                        enddo
!                                    enddo
!                                enddo


                                IF (m.eq.i) THEN
                                    do d=1,SpatOrbs
                                        do g=1,SpatOrbs
                                            Symb=IEOR(AOSymm(g),IEOR(AOSymm(d),AOSymm(z)))
                                            do b=(SymInd(Symb+1)),(SymInd(Symb+1)+MBas(Symb+1)-1)
                                                t1=t1+(coefft1(b,j)*coefft1(g,k)*coefft1(d,l)*real(umat(umatind(z,b,g,d,0,0))%v,8))
                                            enddo
                                        enddo
                                    enddo
                                ENDIF
                                IF (m.eq.j) THEN
                                    Symd=IEOR(AOSymm(i),IEOR(AOSymm(k),AOSymm(z)))
                                    do d=(SymInd(Symd+1)),(SymInd(Symd+1)+MBas(Symd+1)-1)
                                        t2=t2+(CoeffT1(d,l)*TwoIndInts(i,z,k,d))
                                    enddo
                                ENDIF
                                IF (m.eq.k) THEN
                                    do d=1,SpatOrbs
                                        Symb2=IEOR(AOSymm(d),IEOR(AOSymm(i),AOSymm(z)))
                                        do b=(SymInd(Symb2+1)),(SymInd(Symb2+1)+MBas(Symb2+1)-1)
                                            t3=t3+CoeffT1(b,j)*CoeffT1(d,l)*OneIndInts(i,b,z,d)
                                        enddo
                                    enddo
                                ENDIF
                                IF (m.eq.l) THEN
                                    t4=ThreeIndInts(i,j,k,z)
                                ENDIF

                                Deriv4indint=0.D0
                                Deriv4indint=(t1 + t2 + t3 + t4)
                                ! At a particular i, j, k and l; t1, t2, t3 and t4 are only non-zero if m is equal to 
                                ! i, j, k and/or l respectively.

! Deriv4indintsqrd is the derivative of the overall expression for the sum of the squares of the <ij|kl> matrix.
! This accumulates as the loop sums over i,j,k and l.
                                Deriv4indintsqrd=Deriv4indintsqrd+(FourIndInts(i,j,k,l)*Deriv4indint) 
                                ENDIF  
!                            ENDIF
                                
                            enddo
                        enddo
                    enddo
                enddo
                DerivCoeff(z,m)=(2*Deriv4indintsqrd)  
                Force=Force+DerivCoeff(z,m)
                ForceInts=ForceInts+2*Deriv4indintsqrd
            enddo
        enddo

 
        Force=Force/REAL(SpatOrbs**2,8)
        ForceInts=ForceInts/REAL(SpatOrbs**2,8)


! Calculate the derivatives of orthogonalisation condition.
! Have taken this out of the m and z loop to make the shake faster, but can put it back in if start using it a lot.
        IF(tLagrange) THEN
            do m=1,SpatOrbs
                Symm=AOSymm(m)                                             
                ! Gives the symmetry of the orbital m (0 -> 7)
                do z=(SymInd(Symm+1)),(SymInd(Symm+1)+MBas(Symm+1)-1)          
                ! Runs only over those z with the same symm. 
      
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
        INTEGER :: m,z,i,j,Symm
        REAL*8 :: NewCoeff,NewLambda

        DistCs=0.D0 
        do m=1,SpatOrbs
            Symm=AOSymm(m)                                             
            ! Gives the symmetry of the orbital m (0 -> 7)
            do z=(SymInd(Symm+1)),(SymInd(Symm+1)+MBas(Symm+1)-1)          
            ! Only coeffs with sym of m and z the same have non-zero coeffs.    
                NewCoeff=0.D0
                NewCoeff=CoeffT1(z,m)-(TimeStep*DerivCoeff(z,m))
                DistCs=DistCs+abs(TimeStep*DerivCoeff(z,m))
                CoeffT1(z,m)=NewCoeff
            enddo
        enddo
        DistCs=DistCs/(REAL(SpatOrbs**2,8))

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
        INTEGER :: i,j,k,l,a,b,g,d,Symd,Symi
        REAL*8 :: t

!Zero arrays from previous transform
        OneIndInts(:,:,:,:)=0.D0
        TwoIndInts(:,:,:,:)=0.D0
        ThreeIndInts(:,:,:,:)=0.D0
        FourIndInts(:,:,:,:)=0.D0
!<alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i,j,k,l,0,0)

!loop over i's in each symmetry irrep
        do i=1,SpatOrbs
!loop over beta,gamma,delta
            do b=1,SpatOrbs
                do g=1,SpatOrbs
                    Symd=IEOR(AOSymm(g),IEOR(AOSymm(b),AOSymm(i)))
                    do d=(SymInd(Symd+1)),(SymInd(Symd+1)+MBas(Symd+1)-1)
!The <i b | g d> must be symmetry allowed...
                        t=0.D0
                        Symi=AOSymm(i)
                        do a=(SymInd(Symi+1)),(SymInd(Symi+1)+MBas(Symi+1)-1)     
                        ! only CoeffT1(a,i) will only be non-zero if i and a same symmetry.
                            t=t+CoeffT1(a,i)*REAL(UMAT(UMatInd(a,b,g,d,0,0))%v,8)
!this is a sum over alpha to store <i beta | gamma delta> - OneIndInts
                        enddo
                        OneIndInts(i,b,g,d)=t
                        OneIndInts(i,d,g,b)=t
                    enddo
                enddo
            enddo
        enddo

        do k=1,SpatOrbs !loop over k orbitals in each symmetry irrep
            do i=1,k    !i =< k
                do b=1,SpatOrbs
                    Symd=IEOR(AOSymm(b),IEOR(AOSymm(i),AOSymm(k)))
                    do d=(SymInd(Symd+1)),(SymInd(Symd+1)+MBas(Symd+1)-1)
                    !Test that <i b | k d> is symmetry allowed
                        !Can test that (i,k) < (b,d) here
                        t=0.D0
!this is a sum over gamma to store <i beta | j delta> - TwoIndInts
                        do g=1,SpatOrbs
                            t=t+CoeffT1(g,k)*OneIndInts(i,b,g,d)
                        enddo
                        TwoIndInts(i,b,k,d)=t
                        TwoIndInts(k,b,i,d)=t
                        TwoIndInts(i,d,k,b)=t
                        TwoIndInts(k,d,i,b)=t
                    enddo
                enddo
            enddo
        enddo

        do k=1,SpatOrbs !loop over j orbitals in each symmetry irrep
            do j=1,SpatOrbs
                do i=1,k                                !i =< k
                    Symd=IEOR(AOSymm(i),IEOR(AOSymm(j),AOSymm(k)))
                    do d=(SymInd(Symd+1)),(SymInd(Symd+1)+MBas(Symd+1)-1)
                        t=0.D0
!this is a sum over beta to store <i j | k delta> - ThreeIndInts
                        do b=1,SpatOrbs
                            t=t+CoeffT1(b,j)*TwoIndInts(i,b,k,d)
                        enddo
                        ThreeIndInts(i,j,k,d)=t
                        ThreeIndInts(k,j,i,d)=t
                    enddo
                enddo
            enddo
        enddo

        PotEnergy=0.D0
        TwoEInts=0.D0
        PEInts=0.D0
        do l=1,SpatOrbs
            do j=1,l                                    ! j =< l
                do k=1,SpatOrbs
                    Symi=IEOR(AOSymm(k),IEOR(AOSymm(j),AOSymm(l)))
                    do i=(SymInd(Symi+1)),(SymInd(Symi+1)+MBas(Symi+1)-1)    ! i =< k
                        IF(i.le.k) THEN
                            !Test that (i,k) =< (j,l) here
                            !Test that integral is symmetry allowed
                            t=0.D0
!this is a sum over delta to store <i j | k l> - FourIndInts
                            do d=1,SpatOrbs
                                t=t+CoeffT1(d,l)*ThreeIndInts(i,j,k,d)
                            enddo
                        ENDIF 
                        FourIndInts(i,j,k,l)=t
                        FourIndInts(k,j,i,l)=t
                        FourIndInts(i,l,k,j)=t
                        FourIndInts(k,l,i,j)=t
                        IF(.not.((k.eq.i).or.(j.eq.l))) THEN
                            PotEnergy=PotEnergy+(t**2)
                            TwoEInts=TwoEInts+(t**2)
                            PEInts=PEInts+(t**2)
                        ENDIF
                    enddo
                enddo
            enddo
        enddo

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

    END SUBROUTINE Transform2ElInts


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
        ELSEIF(abs(Force).lt.ConvergedForce) THEN
            tNotConverged=.false.
        ENDIF
! IF an ROIteration value is specified, use this to specify the end of the orbital rotation, otherwise use the 
! conversion limit (ConvergedForce).

    END SUBROUTINE TestForConvergence



    SUBROUTINE FinalizeNewOrbs()
! At the end of the orbital rotation, have a set of coefficients CoeffT1 which transform the HF orbitals into a set of linear
! combinations ui which minimise |<ij|kl>|^2.  This is the final subroutine after all iterations (but before the memory deallocation)
! that calculates the final 4 index integrals to be used in the NECI calculation.
        INTEGER :: i,a


! First need to do a final explicit orthonormalisation.  The orbitals are very close to being orthonormal, but not exactly.
! Need to make sure they are exact orthonormal using Gram Schmit maybe...

! Write out some final results of interest, like values of the constraints, values of new coefficients.
        WRITE(6,*) 'The final transformation coefficients'
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


        CALL Transform2ElInts()
! Use these final coefficients to find the FourIndInts(i,j,k,l).
! These are now the <ij|kl> integrals we now want to use instead of the HF UMat.

!        CALL RefillUMATandTMAT2D()        
! UMat is the 4 index integral matrix (2 electron), whereas TMAT2D is the 2 index integral (1 el) matrix

! Calculate the fock matrix, and print it out to see how much the off diagonal terms contribute.
! Also print out the sum of the diagonal elements to compare to the original value.
        CALL CalcFOCKMatrix()
        stop

    ENDSUBROUTINE FinalizeNewOrbs



    SUBROUTINE RefillUMATandTMAT2D()
        INTEGER :: l,k,j,i,a,b
        REAL :: NewTMAT2D


        ! Make the UMAT elements the four index integrals.  These are calculated by transforming the HF orbitals using the
        ! coefficients that have been found
        do l=1,SpatOrbs
            do k=1,SpatOrbs
                do j=1,SpatOrbs
                    do i=1,SpatOrbs
                        UMAT(UMatInd(i,j,k,l,0,0))=HElement(FourIndInts(i,j,k,l))
                    enddo
                enddo
                ! Also calculate the 2 index integrals, and make these the elements of the TMAT2D matrix.
                NewTMAT2D=0.D0
                do a=1,SpatOrbs
                    do b=1,SpatOrbs
                        NewTMAT2D=NewTMAT2D+(CoeffT1(b,k)*CoeffT1(a,l)*REAL(TMAT2D(b,a)%v,8))
                    enddo
                enddo
                TMAT2D(k,l)=HElement(NewTMAT2D)
            enddo
        enddo


    ENDSUBROUTINE RefillUMATandTMAT2D


    SUBROUTINE CalcFOCKMatrix()
        USE SystemData , only : nBasis
        INTEGER :: i,j,a,ArrNewTag,ierr
        REAL*8 , ALLOCATABLE :: ArrNew(:,:)
        REAL*8 :: FOCKDiagSumHF,FOCKDiagSumNew
        CHARACTER(len=*) , PARAMETER :: this_routine='CalcFOCKMatrix'

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

    
        ALLOCATE(ArrNew(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ArrNew',SpatOrbs**2,8,this_routine,ArrNewTag,ierr)
        ArrNew(:,:)=0.D0                     


! First calculate the sum of the diagonal elements of ARR.
! Check if this is already being done.
        FOCKDiagSumHF=0.D0
        do a=1,nBasis        
            FOCKDiagSumHF=FOCKDiagSumHF+Arr(a,2)
        enddo

! Then calculate the fock matrix in the transformed basis, and the sum of the new diagonal elements.
! Our ArrNew is in spatial orbitals, count each twice to compare to above.
        FOCKDiagSumNew=0.D0
        do j=1,SpatOrbs
            do i=1,SpatOrbs
                do a=1,SpatOrbs
                    ArrNew(i,j)=ArrNew(i,j)+(CoeffT1(a,i)*Arr(a,2)*CoeffT1(a,j))
                enddo
            enddo
            FOCKDiagSumNew=FOCKDiagSumNew+(ArrNew(j,j)*2)
        enddo

        WRITE(6,*) 'Sum of the diagonal elements of the fock matrix in the HF basis set = ',FOCKDiagSumHF
        WRITE(6,*) 'Sum of the diagonal elements of the fock matrix in the transformed basis set = ',FOCKDiagSumNew

        WRITE(6,*) 'The fock matrix for the transformed orbitals'
        do j=1,SpatOrbs
            do i=1,SpatOrbs
                WRITE(6,'(F20.10)',advance='no') ArrNew(i,j)
            enddo
            WRITE(6,*) ''
        enddo

        DEALLOCATE(ArrNew)
        CALL LogMemDealloc(this_routine,ArrNewTag)
       

    ENDSUBROUTINE CalcFOCKMatrix



    SUBROUTINE InitLocalOrbs()
        USE SystemData , only : nBasis
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        REAL*8 :: t
        INTEGER :: i,j,k,l,ierr,Const

        WRITE(6,*) "Calculating new molecular orbitals based on mimimisation of <ij|kl>^2 integrals..."
       
! Check for possible errors.
        IF(.not.TwoCycleSymGens) THEN
            CALL Stop_All(this_routine,"ERROR. TwoCycleSymGens is false.  Symmetry is not abelian.") 
        ENDIF

        IF(tLagrange.and.tShake) THEN
            CALL FLUSH(6)
            CALL FLUSH(12)
            CALL Stop_All(this_routine,"ERROR. Both LAGRANGE and SHAKE keywords present in the input. &
            & These two orthonormalisation methods clash.")
        ENDIF

        IF(tLagrange) THEN
            WRITE(6,*) "Using a Lagrange multiplier to attempt to rotate orbitals in a way to maintain orthonormality"
        ELSEIF (tShake) THEN
            WRITE(6,*) "Using the shake algorithm to iteratively find lambdas which maintain orthonormalisation with rotation"
        ELSE
            WRITE(6,*) "Explicity reorthonormalizing orbitals after each rotation."
        ENDIF


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
        TotNoConstraints=(SpatOrbs*(SpatOrbs+1))/2

        
!Set timed routine names
        Rotation_Time%timer_name='RotateTime'
        Shake_Time%timer_name='ShakeTime'
        Findtheforce_Time%timer_name='FindtheForceTime'


 
! Set up constraint labels.
            
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
       
        WRITE(6,*) 'TotNoConstraints = ',TotNoConstraints
        WRITE(6,*) 'Const = ',Const
        IF(Const.ne.TotNoConstraints) THEN
            CALL Stop_all(this_routine,'ERROR in the number of constraints calculated.  lmax does not equal TotNoConstraints')
        ENDIF
            
! Set up symmetry arrays.

        CALL InitSymmArrays()


!Allocate memory

        ALLOCATE(CoeffT1(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffT1',SpatOrbs**2,8,this_routine,CoeffT1Tag,ierr)
        ALLOCATE(CoeffCorT2(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffCorT2',SpatOrbs**2,8,this_routine,CoeffCorT2Tag,ierr)
        ALLOCATE(CoeffUncorT2(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffUncT2',SpatOrbs**2,8,this_routine,CoeffUncorT2Tag,ierr)
        ALLOCATE(DerivCoeff(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('DerivCoeff',SpatOrbs**2,8,this_routine,DerivCoeffTag,ierr)
        ALLOCATE(OneIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('OneIndInts',SpatOrbs**4,8,this_routine,OneIndIntsTag,ierr)
        ALLOCATE(TwoIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('TwoIndInts',SpatOrbs**4,8,this_routine,TwoIndIntsTag,ierr)
        ALLOCATE(ThreeIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts',SpatOrbs**4,8,this_routine,ThreeIndIntsTag,ierr)
        ALLOCATE(FourIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts',SpatOrbs**4,8,this_routine,FourIndIntsTag,ierr)
        
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
        do i=1,SpatOrbs
            CoeffT1(i,i)=1.D0
        enddo
        DerivCoeff(:,:)=0.D0
!These loops can be speeded up with spatial symmetry and pairwise permutation symmetry if needed.
!Fill all the original partially transformed integral arrays with the original two-electron integrals.
        do i=1,SpatOrbs
            do k=1,i
                do j=1,SpatOrbs
                    do l=1,j
                        t=REAL(UMAT(UMatInd(i,j,k,l,0,0))%v,8)
                        IF(.not.((k.eq.i).or.(j.eq.l))) THEN
                            PotEnergy=PotEnergy+(t**2)       !Potential energy starts as this since the orbitals are orthonormal by construction.
                            TwoEInts=TwoEInts+(t**2)
                        ENDIF
                        FourIndInts(i,j,k,l)=t
                        FourIndInts(k,j,i,l)=t
                        FourIndInts(i,l,k,j)=t
                        FourIndInts(k,l,i,j)=t
                        ThreeIndInts(i,j,k,l)=t
                        ThreeIndInts(k,j,i,l)=t
                        ThreeIndInts(i,l,k,j)=t
                        ThreeIndInts(k,l,i,j)=t
                        TwoIndInts(i,j,k,l)=t
                        TwoIndInts(k,j,i,l)=t
                        TwoIndInts(i,l,k,j)=t
                        TwoIndInts(k,l,i,j)=t
                        OneIndInts(i,j,k,l)=t
                        OneIndInts(k,j,i,l)=t
                        OneIndInts(i,l,k,j)=t
                        OneIndInts(k,l,i,j)=t
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
        DEALLOCATE(AOSymm)
        CALL LogMemDealloc(this_routine,AOSymmTag)
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
        DEALLOCATE(OneIndInts)
        CALL LogMemDealloc(this_routine,OneIndIntsTag)
        DEALLOCATE(TwoIndInts)
        CALL LogMemDealloc(this_routine,TwoIndIntsTag)
        DEALLOCATE(ThreeIndInts)
        CALL LogMemDealloc(this_routine,ThreeIndIntsTag)
        DEALLOCATE(FourIndInts)
        CALL LogMemDealloc(this_routine,FourIndIntsTag)

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
            IF(Mod(Iteration,1000).eq.0) THEN
                WRITE(6,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
                WRITE(12,"(I12,5F24.10)") Iteration,PotEnergy,Force,TotCorrectedForce,OrthoNorm,DistCs
            ENDIF
        ENDIF
        CALL FLUSH(6)
        CALL FLUSH(12)

    END SUBROUTINE WriteStats


    SUBROUTINE InitSymmArrays()
        INTEGER :: i,ierr,SymSum
        CHARACTER(len=*) , PARAMETER :: this_routine='InitSymmArrays'
         
        
        ALLOCATE(AOSymm(SpatOrbs),stat=ierr)
        CALL LogMemAlloc('AOSymm',SpatOrbs,4,this_routine,AOSymmTag,ierr)
            
! AOSymm holds the symmetry of each of the spatial orbitals.    
        do i=1,SpatOrbs
            AOSymm(i)=INT(G1(2*i)%sym%S,4)
            ! The symmetry (0 -> 7) of each of the spatial orbitals.
            MBas(AOSymm(i)+1)=MBas(AOSymm(i)+1)+1
            ! Position 1 of SymInd gives the number of orbitals in symmetry 0, position 2 the number in symm 1 etc.
        enddo
        SymInd(1)=1
        SymSum=0
        do i=2,8 
            SymSum=SymSum+MBas(i-1)
            SymInd(i)=SymSum+1
        enddo

!        WRITE(6,*) 'AOSymm'
!        do i=1,SpatOrbs
!            WRITE(6,*) AOSymm(i)
!        enddo
!        WRITE(6,*) 'SymInd'
!        do i=1,8
!            WRITE(6,*) SymInd(i)
!        enddo
!        WRITE(6,*) 'MBas'
!        do i=1,8
!            WRITE(6,*) MBas(i)
!        enddo


    ENDSUBROUTINE InitSymmArrays



! DerivCoeff(k,a) is the unconstrained force on the original coefficients (CoeffT1(a,k)). 
    SUBROUTINE ShakeConstraints()
        INTEGER :: i,j,l,a,m,ShakeIteration,ConvergeCount,ierr,Const
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
        IF(Mod(Iteration,1000).eq.0) WRITE(8,*) 'Orbital rotation iteration = ',Iteration

        ShakeIteration=0
        tShakeNotConverged=.true.


! Before we start iterating, take the current coefficients and find the derivative of the constraints with respect to them.
        
        CALL CalcDerivConstr(CoeffT1,DerivConstrT1)



! Then find the coefficients at time t2, when moved by the completely unconstrained force and the values of the each 
! constraint at these positions.

        Correction(:,:)=0.D0
        CALL FindandUsetheForce(TotUncorForce,TotDiffUncorCoeffs,CoeffUncorT2)

        CALL CalcConstraints(CoeffUncorT2,Constraint,TotConstraints)


! Write stats from the beginning of the iteration to output.            
!        IF(Mod(Iteration,1000).eq.0.or.Iteration.eq.1) CALL WriteShakeOUTstats01()

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
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    do l=1,TotNoConstraints
                        Correction(a,m)=Correction(a,m)+(ShakeLambda(l)*DerivConstrT1(a,m,l)) 
                    enddo
                enddo
            enddo

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
!            IF(Mod(Iteration,1000).eq.0.or.Iteration.eq.1) CALL WriteShakeOUTstats02(ShakeIteration,TotLambdas,ConvergeCount) ! add correction to this.


! and to SHAKEstats file:
            CALL FLUSH(6)
            CALL FLUSH(8)
            IF(Mod(Iteration,1000).eq.0) THEN
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
        INTEGER :: l,m,a
            

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
                WRITE(6,'(10F20.10)',advance='no') Constraint(l)
            enddo
 

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
 
        REAL*8 :: CurrCoeff(SpatOrbs,SpatOrbs)
        REAL*8 :: DerivConstr(SpatOrbs,SpatOrbs,TotNoConstraints)
        INTEGER :: l,m,i,j,a,Symm
       
        DerivConstr(:,:,:)=0.D0
        do l=1,TotNoConstraints
            i=lab(1,l)
            j=lab(2,l)
            do m=1,SpatOrbs 
!                Symm=AOSymm(m) 
!                do a=(SymInd(Symm+1)),(SymInd(Symm+1)+MBas(Symm+1)-1)
                IF (m.eq.i.and.m.eq.j) THEN
                    do a=1,SpatOrbs
                        DerivConstr(a,m,l)=CurrCoeff(a,j)*2
                    enddo
                ELSEIF (m.eq.j) THEN
                    do a=1,SpatOrbs
                        DerivConstr(a,m,l)=CurrCoeff(a,i) 
                    enddo
                ELSEIF (m.eq.i) THEN
                    do a=1,SpatOrbs
                        DerivConstr(a,m,l)=CurrCoeff(a,j)
                    enddo
                ENDIF
                ! DerivConstrT1 stays the same throughout the iterations
            enddo
        enddo
                  
    ENDSUBROUTINE CalcDerivConstr



    SUBROUTINE FindandUsetheForce(TotForce,TotDiffCoeffs,CoeffT2)
    ! This takes the current lambdas with the derivatives of the constraints and calculates a force
    ! for each cm, with an orthonormalisation correction.
    ! This is then used to rotate the coefficients by a defined timestep.
        INTEGER :: a,l,m,Symm
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

                
        do m=1,SpatOrbs
            ! FIND THE FORCE 
            Symm=AOSymm(m) 
            do a=(SymInd(Symm+1)),(SymInd(Symm+1)+MBas(Symm+1)-1)

!            do a=1,SpatOrbs

                ForceCorrect(a,m)=DerivCoeff(a,m)-Correction(a,m)
                ! find the corrected force. (in the case where the uncorrected force is required, correction is set to 0.
                ! DerivCoeff(m,a) is the derivative of |<ij|kl>|^2 w.r.t cm without any constraints (no lambda terms).
                ! ForceCorrect is then the latest force on coefficients.  This is iteratively being corrected so that
                ! it will finally move the coefficients so that they remain orthonormal.
            
            ! USE THE FORCE
                CoeffT2(a,m)=CoeffT1(a,m)-(TimeStep*ForceCorrect(a,m))
                ! Using the force to calculate the coefficients at time T2 (hopefully more orthonomal than those calculated in
                ! the previous iteration).
                
                
                ! Calculate parameters for printing
                TotForce=TotForce+ForceCorrect(a,m)
                TotDiffCoeffs=TotDiffCoeffs+ABS(CoeffT2(a,m)-CoeffT1(a,m))

            enddo
        enddo
            

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


!            WRITE(6,*) 'DerivConstrT1T2 '
!            do l=1,TotNoConstraints
!                do n=1,TotNoConstraints
!                    WRITE(6,'(10F20.10)',advance='no') DerivConstrT1T2(n,l)
!                enddo
!                write(6,*) ''
!            enddo
           

! Invert the matrix to calculate the lambda values.
! LU decomposition.
            call dgetrf(TotNoConstraints,TotNoConstraints,DerivConstrT1T2,TotNoConstraints,ipiv,info)

            if(info.ne.0) THEN
                WRITE(6,*) 'info ',info
                CALL Stop_All(this_routine,"The LU decomposition of matrix inversion failed...")
            endif

            do n=1,TotNoConstraints
                ShakeLambdaNew(n)=Constraint(n)/(TimeStep*(-1))
            enddo
            ! These are actually still the constraint values, but now Lambda(n) can go into dgetrs as the constraints (B in AX=B), 
            ! and come out as the computed lambdas (X).

            call dgetrs('N',TotNoConstraints,1,DerivConstrT1T2,TotNoConstraints,ipiv,ShakeLambdaNew,TotNoConstraints,info)
            if(info.ne.0) CALL Stop_All(this_routine,"Error in dgetrs, solving for the lambdas...")

!            WRITE(6,*) 'Lambdas successfully calculated, beginning next shake iteration'
          
!    enddo

    ENDSUBROUTINE FullShake
   


END MODULE RotateOrbsMod
