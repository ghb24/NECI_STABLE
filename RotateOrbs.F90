!TO DO
!Spatial symmetry
!Finalize orbs to fully test
!Permutational symmetry where possible
!Parallelize

MODULE RotateOrbsMod

    USE Global_utilities
    USE IntegralsData , only : UMAT
    USE UMatCache , only : UMatInd
    USE SystemData , only : ConvergedForce,TimeStep,tLagrange,tShake,tShakeApprox,ShakeConverged
    IMPLICIT NONE
    REAL*8 , ALLOCATABLE :: Coeff(:,:),CoeffT2(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:),Lab(:,:),ForceCorrect(:,:),Correction(:)
    REAL*8 , ALLOCATABLE :: Constraint(:),Lambda(:),DerivConstrT1(:,:,:),DerivConstrT2(:,:,:),DerivConstrT1T2(:,:),DerivConstrT1T2Diag(:)
    REAL*8 , ALLOCATABLE :: OneIndInts(:,:,:,:), TwoIndInts(:,:,:,:), ThreeIndInts(:,:,:,:), FourIndInts(:,:,:,:)   !These are the arrays to store the
!partially transformed integrals for each iteration. The FourIndInts are the full <ij|kl> integrals for the current iteration.
!OneIndInts = <i b|g d> ; TwoIndInts = <i b|k d> ; ThreeIndInts = <i j|k d>
    INTEGER :: SpatOrbs,OneIndIntsTag,TwoIndIntsTag,ThreeIndIntsTag,FourIndIntsTag,LabTag,ForceCorrectTag,CorrectionTag
    INTEGER :: CoeffTag,CoeffT2Tag,LambdasTag,DerivCoeffTag,DerivLambdaTag,Iteration,TotNoConstraints
    INTEGER :: LambdaTag,ConstraintTag,DerivConstrT1Tag,DerivConstrT2Tag,DerivConstrT1T2Tag,DerivConstrT1T2DiagTag
    LOGICAL :: tNotConverged
    REAL*8 :: OrthoNorm,PotEnergy,Force,TwoEInts,DistCs,OrthoForce,DistLs,LambdaMag,PEInts,PEOrtho,ForceInts
    REAL*8 :: OrthoFac=1.D0

    contains

    SUBROUTINE RotateOrbs()

        CALL InitLocalOrbs()        ! Set defaults, allocate arrays, write out headings for OUTPUT, set integarals to HF values.

        CALL WriteStats()           ! write out the original stats before any rotation.

        tNotConverged=.true.
        do while(tNotConverged)     ! rotate the orbitals until the sum of the four index integral falls below a chose convergence value.

            Iteration=Iteration+1

            CALL FindNewOrbs()      ! bulk of the calculation.
                                    ! do the actual transformations, moving the coefficients by a timestep according to the calculated force. 

            CALL WriteStats()       ! write out the stats for this iteration.

        enddo           

        WRITE(6,*) "Convergence criterion met. Finalizing new orbitals..."

!Make symmetry, orbitals, one/two-electron integrals consistent with rest of neci
!        CALL FinalizeNewOrbs()

!option to create new FCIDUMP?

        CALL DeallocateMem()

        CLOSE(12)

        CALL Stop_All("RotateOrbs","This code is still in the testing phase")

    END SUBROUTINE RotateOrbs

    SUBROUTINE FindNewOrbs()
        
        CALL Transform2ElInts()     ! Find the partially (and completely) transfored 4 index integrals to be used in further calcs.

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
!The coefficients coeff(a,m) are now those that have been shifted by the time step.

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
                    OrthoNorm=OrthoNorm+Coeff(a,i)*Coeff(a,j)
                enddo
            enddo
        enddo
        OrthoNorm=OrthoNorm-real(SpatOrbs,8)
        OrthoNorm=(OrthoNorm*2.D0)/REAL((SpatOrbs*(SpatOrbs+1.D0)),8)

    END SUBROUTINE TestOrthonormality

    SUBROUTINE FindTheForce
        INTEGER :: m,z,i,j,k,l,a,b,g,d
        REAL*8 :: Deriv4indint,Deriv4indintsqrd,t1,t2,t3,t4,LambdaTerm1,LambdaTerm2
        CHARACTER(len=*) , PARAMETER :: this_routine='FindtheForce'
      
! Running over m and z, covers all matrix elements of the force matrix (derivative 
! of equation we are minimising, with respect to each translation coefficient) filling 
! them in as it goes.

        DerivCoeff(:,:)=0.D0
        Force=0.D0
        Deriv4indintsqrd=0.D0
        ForceInts=0.D0
        OrthoForce=0.D0
        
        do m=1,SpatOrbs
! To include: symmetry requirement that z must be from the same irrep as m
            do z=1,SpatOrbs

                Deriv4indintsqrd=0.D0
                    
! This runs over the full <ij|kl> integrals from the previous iteration. 
! In the future we can take advantage of the permutational symmetry of the matrix elements
                do i=1,SpatOrbs
                    do k=1,i-1                       ! i < k
                        do j=1,SpatOrbs
                            do l=1,j-1               ! j < l

! To include: If statement to check the symmetry of i x j x k x l is A1.
!                                If((((i-1)*nBasis)+k).ge.(((j-1)*nBasis)+l)) THEN

! Already calculated the four index integrals and some partial integrals
! in Transform2ElInts routine of the previous iteration.

! Deriv4indint is the derivative of <ij|kl> with respect to a coefficient
! c(m,z), at the current m,z,i,j,k,l of the loop.
! This involves 4 terms in which the derivative is taken w.r.t the i,j,k or l position in <ij|kl>. 
                                t1=0.D0
                                t2=0.D0
                                t3=0.D0
                                t4=0.D0
                                
                                IF (m.eq.i) THEN
                                    do b=1,SpatOrbs
                                        do g=1,SpatOrbs
                                            do d=1,SpatOrbs
                                                t1=t1+(Coeff(b,j)*Coeff(g,k)*Coeff(d,l)*REAL(UMAT(UMatInd(z,b,g,d,0,0))%v,8))
                                            enddo
                                        enddo
                                    enddo
                                ENDIF

                                IF (m.eq.j) THEN
                                    do d=1,SpatOrbs     
                                        t2=t2+(Coeff(d,l)*TwoIndInts(i,z,k,d))
                                    enddo
                                ENDIF

                                IF (m.eq.k) THEN
                                    do b=1,SpatOrbs
                                        do d=1,SpatOrbs
                                            t3=t3+Coeff(b,j)*Coeff(d,l)*OneIndInts(i,b,z,d)
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
                            enddo
                        enddo
                    enddo
                enddo

! Calculate the derivatives of orthogonalisation condition.

                IF(tLagrange) THEN
                       
                    LambdaTerm1=0.D0
                    LambdaTerm2=0.D0
                    
                    do j=1,SpatOrbs
                        LambdaTerm1=LambdaTerm1+(Lambdas(m,j)*Coeff(z,j))
                        LambdaTerm2=LambdaTerm2+(Lambdas(j,m)*Coeff(z,j))
                    enddo

! DerivCoeff is 'the force'.  I.e. the derivative of |<ij|kl>|^2 with 
! respect to each transformation coefficient.  It is the values of this matrix that will tend to 0 as
! we minimise the sum of the |<ij|kl>|^2 values.
! With the Lagrange keyword this includes orthonormality conditions, otherwise it is simply the unconstrained force.
                    DerivCoeff(z,m)=(2*Deriv4indintsqrd)-LambdaTerm1-LambdaTerm2
                    OrthoForce=OrthoForce-LambdaTerm1-LambdaTerm2
                ELSE
                    DerivCoeff(z,m)=(2*Deriv4indintsqrd)  
                    Force=Force+DerivCoeff(z,m)
                    ForceInts=ForceInts+2*Deriv4indintsqrd
                ENDIF

            enddo
        enddo

        Force=Force/REAL(SpatOrbs**2,8)
        ForceInts=ForceInts/REAL(SpatOrbs**2,8)
        OrthoForce=OrthoForce/REAL(SpatOrbs**2,8)

!If doing a lagrange calc we also need to find the force on the lambdas to ensure orthonormality...
        IF(tLagrange) THEN
            DerivLambda(:,:)=0.D0
            do i=1,SpatOrbs
                do j=1,i
                    do a=1,SpatOrbs
                        DerivLambda(i,j)=DerivLambda(i,j)+Coeff(a,i)*Coeff(a,j)
                    enddo
                    DerivLambda(j,i)=DerivLambda(i,j)
                enddo
            enddo
            do i=1,SpatOrbs
                DerivLambda(i,i)=DerivLambda(i,i)-1.D0
            enddo
        ENDIF


    END SUBROUTINE FindTheForce

    
    SUBROUTINE UseTheForce()
! This routine takes the old translation coefficients and Lambdas and moves them by a timestep in the direction 
! of the calculated force.
        INTEGER :: m,z,i,j
        REAL*8 :: NewCoeff,NewLambda

        DistCs=0.D0 
        do m=1,SpatOrbs
            do z=1,SpatOrbs
                NewCoeff=0.D0
                NewCoeff=Coeff(z,m)-(TimeStep*DerivCoeff(z,m))
                DistCs=DistCs+abs(TimeStep*DerivCoeff(z,m))
                Coeff(z,m)=NewCoeff
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
!            CALL OrthoNormx(SpatOrbs,SpatOrbs,Coeff) !Explicitly orthonormalize the coefficient vectors.
        ENDIF

    ENDSUBROUTINE UseTheForce


    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.
    SUBROUTINE Transform2ElInts()
        INTEGER :: i,j,k,l,a,b,g,d
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
                    do d=1,b    !delta =< beta
!If statement here - the <i b | g d> must be symmetry allowed...
                        t=0.D0
                        do a=1,SpatOrbs     !This loop only wants to be over the atomic orbitals in the symmetry class of the current i loop
!sum over alpha to store <i beta | gamma delta> - OneIndInts
                            t=t+Coeff(a,i)*REAL(UMAT(UMatInd(a,b,g,d,0,0))%v,8)
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
                    do d=1,b    !d =< b
                        !Can test that (i,k) < (b,d) here
                        !Can test that <i b | k d> is symmetry allowed here
                        t=0.D0
!sum over gamma to store <i beta | j delta> - TwoIndInts
                        do g=1,SpatOrbs
                            t=t+Coeff(g,k)*OneIndInts(i,b,g,d)
                        enddo
                        TwoIndInts(i,b,k,d)=t
                        TwoIndInts(k,b,i,d)=t
                        TwoIndInts(i,d,k,b)=t
                        TwoIndInts(k,d,i,b)=t
                    enddo
                enddo
            enddo
        enddo

        do j=1,SpatOrbs !loop over j orbitals in each symmetry irrep
            do i=1,SpatOrbs
                do k=1,i    !k =< i
                    do d=1,SpatOrbs
                        !Test is symmetry allowed here
                        t=0.D0
!sum over beta to store <i j | k delta> - ThreeIndInts
                        do b=1,SpatOrbs
                            t=t+Coeff(b,j)*TwoIndInts(i,b,k,d)
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
        do l=1,SpatOrbs !loop over l orbitals in each symmetry irrep
            do j=1,l    ! j =< l
                do i=1,SpatOrbs
                    do k=1,i    ! k =< i
                        !Test that (i,k) =< (j,l) here
                        !Test that integral is symmetry allowed
                        t=0.D0
!sum over delta to store <i j | k l> - FourIndInts
                        do d=1,SpatOrbs
                            t=t+Coeff(d,l)*ThreeIndInts(i,j,k,d)
                        enddo
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
                        t=Coeff(a,i)*Coeff(a,j)
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

    IF(Iteration.eq.5) THEN
        tNotConverged=.false.
    ENDIF
!        IF(tLagrange) THEN
!            IF((abs(Force).lt.ConvergedForce).and.(abs(OrthoForce).lt.ConvergedForce)) THEN
!                tNotConverged=.false.
!            ENDIF
!        ELSE
!            IF(abs(Force).lt.ConvergedForce) THEN
!                tNotConverged=.false.
!            ENDIF
!        ENDIF

    END SUBROUTINE TestForConvergence


    SUBROUTINE InitLocalOrbs()
        USE SystemData , only : nBasis
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        REAL*8 :: t
        INTEGER :: i,j,k,l,ierr,Const

        WRITE(6,*) "Calculating new molecular orbitals based on mimimisation of <ij|kl>^2 integrals..."
        IF(tLagrange) THEN
            IF(tShake) CALL Stop_All(this_routine,"ERROR. Both LAGRANGE and SHAKE keywords present in the input. &
            & These two orthonormalisation methods clash.")
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
     
 
! Set up constraint labels.
            
        ALLOCATE(Lab(2,TotNoConstraints),stat=ierr)
        CALL LogMemAlloc('Lab',2*TotNoConstraints,4,this_routine,LabTag,ierr)
        Lab(:,:)=0.D0                     

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
            WRITE(6,*) 'ERROR in the number of constraints calculated.  lmax does not equal TotNoConstraints'
            STOP
        ENDIF
            

!Allocate memory

        ALLOCATE(Coeff(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('Coeff',SpatOrbs**2,8,this_routine,CoeffTag,ierr)
        ALLOCATE(CoeffT2(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('CoeffT2',SpatOrbs**2,8,this_routine,CoeffT2Tag,ierr)
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
            ALLOCATE(Lambda(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('Lambda',TotNoConstraints,8,this_routine,LambdaTag,ierr)
            Lambda(:)=0.D0                     
            ALLOCATE(Constraint(TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('Constraint',TotNoConstraints,8,this_routine,ConstraintTag,ierr)
            ALLOCATE(DerivConstrT1(SpatOrbs,SpatOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT1',SpatOrbs*TotNoConstraints*SpatOrbs,8,this_routine,DerivConstrT1Tag,ierr)
            DerivConstrT1(:,:,:)=0.D0
            ALLOCATE(DerivConstrT2(SpatOrbs,SpatOrbs,TotNoConstraints),stat=ierr)
            CALL LogMemAlloc('DerivConstrT2',SpatOrbs*TotNoConstraints*SpatOrbs,8,this_routine,DerivConstrT2Tag,ierr)
            ALLOCATE(ForceCorrect(SpatOrbs,SpatOrbs),stat=ierr)
            CALL LogMemAlloc('ForceCorrect',SpatOrbs**2,8,this_routine,ForceCorrectTag,ierr)
            ALLOCATE(Correction(SpatOrbs),stat=ierr)
            CALL LogMemAlloc('Correction',SpatOrbs,8,this_routine,CorrectionTag,ierr)
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
        Coeff(:,:)=0.D0
        do i=1,SpatOrbs
            Coeff(i,i)=1.D0
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
            WRITE(12,"(A)") "# Iteration   2.PotEnergy   3.Force      4.OrthoNormCondition   5.DistMovedbyCs"
            WRITE(6,"(A)") "Iteration   2.PotEnergy   3.Force   4.OrthoNormCondition   5.DistMovedbyCs"
        ENDIF

    END SUBROUTINE InitLocalOrbs

    SUBROUTINE DeallocateMem()
        CHARACTER(len=*) , PARAMETER :: this_routine='DeallocateMem'
            
        WRITE(6,*) "1999"
        CALL FLUSH(6)
        DEALLOCATE(Lab)
        CALL LogMemDealloc(this_routine,LabTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)

        DEALLOCATE(Coeff)
        CALL LogMemDealloc(this_routine,CoeffTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(CoeffT2)
        CALL LogMemDealloc(this_routine,CoeffT2Tag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(Lambdas)
        CALL LogMemDealloc(this_routine,LambdasTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(DerivCoeff)
        CALL LogMemDealloc(this_routine,DerivCoeffTag)
        WRITE(6,*) "15"
        CALL FLUSH(6)
        DEALLOCATE(DerivLambda)
        CALL LogMemDealloc(this_routine,DerivLambdaTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(OneIndInts)
        CALL LogMemDealloc(this_routine,OneIndIntsTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(TwoIndInts)
        CALL LogMemDealloc(this_routine,TwoIndIntsTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(ThreeIndInts)
        CALL LogMemDealloc(this_routine,ThreeIndIntsTag)
        WRITE(6,*) "1"
        CALL FLUSH(6)
        DEALLOCATE(FourIndInts)
        CALL LogMemDealloc(this_routine,FourIndIntsTag)
        WRITE(6,*) "199"
        CALL FLUSH(6)

        IF(tShake) THEN
            DEALLOCATE(Lambda)
            CALL LogMemDealloc(this_routine,LambdaTag)
        WRITE(6,*) "19"
        CALL FLUSH(6)
            DEALLOCATE(Constraint)
            CALL LogMemDealloc(this_routine,ConstraintTag)
        WRITE(6,*) "11111111"
        CALL FLUSH(6)
            DEALLOCATE(DerivConstrT1)
            CALL LogMemDealloc(this_routine,DerivConstrT1Tag)
        WRITE(6,*) "111111"
        CALL FLUSH(6)
            DEALLOCATE(DerivConstrT2)
            CALL LogMemDealloc(this_routine,DerivConstrT2Tag)
        WRITE(6,*) "11111"
        CALL FLUSH(6)
            DEALLOCATE(ForceCorrect)
            CALL LogMemDealloc(this_routine,ForceCorrectTag)
        WRITE(6,*) "1111"
        CALL FLUSH(6)
            DEALLOCATE(Correction)
            CALL LogMemDealloc(this_routine,CorrectionTag)
        WRITE(6,*) "111"
        CALL FLUSH(6)
            IF(tShakeApprox) THEN
                DEALLOCATE(DerivConstrT1T2)
                CALL LogMemDealloc(this_routine,DerivConstrT1T2Tag)
            ELSE
                DEALLOCATE(DerivConstrT1T2Diag)
                CALL LogMemDealloc(this_routine,DerivConstrT1T2DiagTag)
        WRITE(6,*) "11"
        CALL FLUSH(6)
            ENDIF
        ENDIF


    END SUBROUTINE DeallocateMem
        
    SUBROUTINE WriteStats()

        IF(tLagrange) THEN
            WRITE(6,"(I7,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
            WRITE(12,"(I7,11F18.10)") Iteration,PotEnergy,PEInts,PEOrtho,Force,ForceInts,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
        ELSE
            WRITE(6,"(I7,4F18.10)") Iteration,PotEnergy,Force,OrthoNorm,DistCs
            WRITE(12,"(I7,4F18.10)") Iteration,PotEnergy,Force,OrthoNorm,DistCs
        ENDIF
        CALL FLUSH(6)
        CALL FLUSH(12)

    END SUBROUTINE WriteStats




! DerivCoeff(k,a) is the unconstrained force on the original coefficients (Coeff(a,k)). 
    SUBROUTINE ShakeConstraints()
        INTEGER :: i,j,l,a,m,ShakeIteration,ConvergeCount,ierr,Const
        REAL*8 :: TotConstraints,TotCorrectedForce,TotDiffCoeffs,TotLambdas
        LOGICAL :: tShakeNotConverged
        CHARACTER(len=*), PARAMETER :: this_routine='ShakeConstraints'


        WRITE(6,*) "Beginning shakeconstraints calculation"
        IF(Iteration.eq.1) THEN
            OPEN(8,FILE='SHAKEstats',STATUS='unknown')
            WRITE(8,'(A20,4A35,A20)') 'Shake Iteration','Sum Lambdas','Total of corrected forces','Tot movement of coefficients',&
                                        &'Sum constraint contributions','Converge count'
        ENDIF
        WRITE(8,*) 'Orbital rotation iteration = ',Iteration

        ShakeIteration=0
        tShakeNotConverged=.true.
! Actually starting the calculation.
        do while (tShakeNotConverged)

            ShakeIteration=ShakeIteration+1
            
            ForceCorrect(:,:)=0.D0          ! Zeroing terms that are re-calculated each iteration.
            CoeffT2(:,:)=0.D0
            Constraint(:)=0.D0
            DerivConstrT2(:,:,:)=0.D0
            TotLambdas=0.D0

! Write stats from the beginning of the iteration to output.            
            CALL WriteShakeOUTstats01(ShakeIteration,TotLambdas)

            do m=1,SpatOrbs                     ! Run through all vectors cm, finding the corrected force each.
                
                ! Calculate the derivative of the constraints w.r.t the original coefficients (at t1).
                ! This only needs to be done once, because the original coefficients are the same in each iteration.
                IF(ShakeIteration.eq.1) THEN
                    do l=1,TotNoConstraints
                        CALL CalcDerivConstr(m,l,Coeff,DerivConstrT1)
                    enddo
                ENDIF

                ! For a particular set of coefficients cm:
                ! Force(corrected)=Force(uncorrected)-Lambdas.DerivConstrT1
                ! Use these derivatives, and the current lambdas to find the trial corrected force.
                ! Then use this to get the (trial) shifted coefficients.

                CALL FindandUsetheForce(m,TotCorrectedForce,TotDiffCoeffs)
                
            enddo


! Use these new shifted coefficients to calculate the values of each constraint, and the derivative of the constraints 
! (at time t2).
            
            CALL CalcConstraints(ConvergeCount,TotConstraints)


! Write out stats of interest to SHAKEstats file:
            CALL FLUSH(6)
            CALL FLUSH(8)
            WRITE(8,'(I20,4F35.20,I20)') ShakeIteration,TotLambdas,TotCorrectedForce,TotDiffCoeffs,TotConstraints,ConvergeCount 
! and to output:
            CALL WriteShakeOUTstats02(ConvergeCount)


            
! Test for convergence, if convergence is reached, make the new coefficients the original ones to start the whole process again.
! Then exit out of this do loop and hence the subroutine.
            CALL TestShakeConvergence(ShakeIteration,ConvergeCount,tShakeNotConverged)

! If the convergence criteria is met, exit out of this subroutine, a rotation has been made which keeps the coefficients 
! orthogonal.

! Otherwise, use either the full matrix inversion method to find a new set of lambdas, or the shake algorithm (in which case 
! SHAKEAPPROX is required in the system block of the input.

            IF(tShakeApprox) THEN
                WRITE(6,*) 'Using shake approximation to find new lambdas'
                CALL ShakeApproximation()
            ELSE
                WRITE(6,*) 'Using the full diagonalisation shake method to find new lambdas'
                CALL FullShake()
            ENDIF
    
    enddo


    ENDSUBROUTINE ShakeConstraints




    SUBROUTINE WriteShakeOUTstats01(ShakeIteration,TotLambdas)
        INTEGER :: ShakeIteration,l,m,a
        REAL*8 :: TotLambdas
            
    
            WRITE(6,*) 'Iteration number ,', ShakeIteration
            WRITE(6,*) 'Original coefficients'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    WRITE(6,'(4F20.10)',advance='no') Coeff(a,m)
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

            WRITE(6,*) 'Lambdas used for this iteration iteration'
            do l=1,TotNoConstraints
                WRITE(6,'(10F20.10)') Lambda(l)
                TotLambdas=TotLambdas+Lambda(l)
            enddo

    ENDSUBROUTINE WriteShakeOUTstats01


    SUBROUTINE WriteShakeOUTstats02(ConvergeCount)
        INTEGER :: m,a,l,ConvergeCount
    
            WRITE(6,*) 'Corrected Force'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    WRITE(6,'(4F20.10)',advance='no') ForceCorrect(a,m)
                enddo
                WRITE(6,*) ''
            enddo
    
            WRITE(6,*) 'Coefficients having been shifted by the force (at t2)'
            do m=1,SpatOrbs
                do a=1,SpatOrbs
                    WRITE(6,'(4F20.10)',advance='no') CoeffT2(a,m)
                enddo
                WRITE(6,*) ''
            enddo

            WRITE(6,*) 'The value of each constraint with these new coefficients'
            do l=1,TotNoConstraints
                WRITE(6,'(10F20.10)',advance='no') Constraint(l)
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



    SUBROUTINE CalcDerivConstr(m,l,CurrCoeff,DerivConstr)
    ! This calculates the derivative of each of the orthonormalisation constraints, l, with respect
    ! to each set of coefficients cm.
 
        REAL*8 :: CurrCoeff(SpatOrbs,SpatOrbs)
        REAL*8 :: DerivConstr(SpatOrbs,SpatOrbs,TotNoConstraints)
        INTEGER :: l,m,i,j,a
        
       
        i=lab(1,l)
        j=lab(2,l)
        do a=1,SpatOrbs
            IF (m.eq.i.and.m.eq.j) THEN
                DerivConstr(a,m,l)=CurrCoeff(a,j)*2
            ELSEIF (m.eq.j) THEN
                DerivConstr(a,m,l)=CurrCoeff(a,i) 
            ELSEIF (m.eq.i) THEN
                DerivConstr(a,m,l)=CurrCoeff(a,j)
            ELSE
                DerivConstr(a,m,l)=0.D0
            ENDIF
            ! DerivConstrT1 stays the same throughout the iterations
        enddo

    
    ENDSUBROUTINE CalcDerivConstr



    SUBROUTINE FindandUsetheForce(m,TotCorrectedForce,TotDiffCoeffs)
    ! This takes the current lambdas with the derivatives of the constraints and calculates a force
    ! for each cm, with an orthonormalisation correction.
    ! This is then used to rotate the coefficients by a defined timestep.
        INTEGER :: a,l,m
        REAL*8 :: TotCorrectedForce,TotDiffCoeffs

        IF(m.eq.1) THEN
            WRITE(6,*) 'm value then'
            WRITE(6,*) 'correction terms'
        ENDIF
        WRITE(6,*) m 

        ! FIND THE FORCE
        ! Use the lambdas of this iteration to calculate the correction to the force due to the constraints.
        Correction(:)=0.D0              
        do a=1,SpatOrbs
            do l=1,TotNoConstraints
                Correction(a)=Correction(a)+(Lambda(l)*DerivConstrT1(a,m,l)) 
                ! Correction is specific to this m, overwrite every time.
            enddo
       
            ForceCorrect(a,m)=DerivCoeff(a,m)-Correction(a)
            ! find the corrected force.
            ! DerivCoeff(m,a) is the derivative of |<ij|kl>|^2 w.r.t cm without any constraints (no lambda terms).
            ! ForceCorrect is then the latest force on coefficients.  This is iteratively being corrected so that
            ! it will finally move the coefficients so that they remain orthonormal.
        
        ! USE THE FORCE
            CoeffT2(a,m)=Coeff(a,m)-(TimeStep*ForceCorrect(a,m))
            ! Using the force to calculate the coefficients at time T2 (hopefully more orthonomal than those calculated in
            ! the previous iteration).
            
            
            ! Calculate parameters for printing
            TotCorrectedForce=TotCorrectedForce+ForceCorrect(a,m)
            TotDiffCoeffs=TotDiffCoeffs+ABS(CoeffT2(a,m)-Coeff(a,m))

            WRITE(6,'(4F20.10)') Correction(a)

        enddo
        

    ENDSUBROUTINE FindandUsetheForce



    SUBROUTINE CalcConstraints(ConvergeCount,TotConstraints)  
    ! This calculates the value of each orthonomalisation constraint, using the shifted coefficients.
    ! Each of these should tend to 0 when the coefficients become orthonomal.
        INTEGER :: l,i,j,m,ConvergeCount
        REAL*8 :: TotConstraints


            TotConstraints=0.D0
            ConvergeCount=0

            do l=1,TotNoConstraints
                i=lab(1,l)
                j=lab(2,l)
                IF(i.eq.j) THEN
                    Constraint(l)=Dot_Product(CoeffT2(:,i),CoeffT2(:,j))-1.D0
                ELSE
                    Constraint(l)=Dot_Product(CoeffT2(:,i),CoeffT2(:,j))
                    ! Each of these components should tend towards 0 when the coefficients become orthonormal.
                ENDIF


! Find the derivative of the constraints with respect to the new coefficients (the coefficients at time t2).
                do m=1,SpatOrbs
                    CALL CalcDerivConstr(m,l,CoeffT2,DerivConstrT2)
                enddo
                
                TotConstraints=TotConstraints+ABS(Constraint(l))
                ! Sum of all Contraint componenets - indication of overall orthonormality.
        
                IF(ABS(Constraint(l)).gt.ShakeConverged) ConvergeCount=ConvergeCount+1
                ! Count the number of constraints which are still well above 0.
                
            enddo

    ENDSUBROUTINE CalcConstraints




    SUBROUTINE TestShakeConvergence(ShakeIteration,ConvergeCount,tShakeNotConverged)
    ! ConvergeCount is the number of constraints that individually have values below
    ! the specified convergence criteria.
    ! If this = 0, the shake is converged, else keep iterating.
        INTEGER :: ConvergeCount,ShakeIteration
        LOGICAL :: tShakeNotConverged


!            IF(ConvergeCount.eq.0) THEN
!                tShakeNotConverged=.false.
!                WRITE(6,*) 'Convergence reached in the shake algorithm'
!                WRITE(6,*) 'All constraints have values less than ',ShakeConverged

! If convergence is reached, make the new coefficients coeff, to start the rotation iteration again.

!                do m=1,SpatOrbs
!                    do a=1,SpatOrbs
!                        Coeff(a,m)=CoeffT2(a,m)
!                    enddo
!                enddo
!            ENDIF
!            IF(.not.tShakeNotConverged) EXIT

! Temporary criteria to stop, to limit the number of iterations while debugging.
            IF(ShakeIteration.eq.10) THEN
                WRITE(6,*) 'stopped at iteration, ',ShakeIteration
                tShakeNotConverged=.false.
            endif


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
                Lambda(l)=Constraint(l)/((-1)*TimeStep*DerivConstrT1T2Diag(l))
                WRITE(6,*) DerivConstrT1T2Diag(l)
            enddo
            
            WRITE(6,*) 'New Lambdas calculated by approx method'
            do l=1,TotNoConstraints
                WRITE(6,*) Lambda(l)
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
                        DerivConstrT1T2(l,n)=DerivConstrT1T2(l,n)+(Dot_Product(DerivConstrT1(:,m,l),DerivConstrT2(:,m,n)))
                    enddo
                enddo
            enddo       ! have filled up whole matrix


            WRITE(6,*) 'DerivConstrT1T2 '
            do l=1,TotNoConstraints
                do n=1,TotNoConstraints
                    WRITE(6,'(10F20.10)',advance='no') DerivConstrT1T2(l,n)
                enddo
                write(6,*) ''
            enddo
            
! Invert the matrix to calculate the lambda values.
! LU decomposition.
            call dgetrf(TotNoConstraints,TotNoConstraints,DerivConstrT1T2,TotNoConstraints,ipiv,info)

            if(info.ne.0) THEN
                WRITE(6,*) 'info ',info
                CALL Stop_All(this_routine,"The LU decomposition of matrix inversion failed...")
            endif

            do l=1,TotNoConstraints
                Lambda(l)=Constraint(l)/(TimeStep*(-1))
            enddo
            ! These are actually still the constraint values, but now Lambda(l) can go into dgetrs as the constraints (B in AX=B), 
            ! and come out as the computed lambdas (X).

            call dgetrs('N',TotNoConstraints,1,DerivConstrT1T2,TotNoConstraints,ipiv,Lambda,TotNoConstraints,info)
            if(info.ne.0) CALL Stop_All(this_routine,"Error in dgetrs, solving for the lambdas...")

            WRITE(6,*) 'Lambdas successfully calculated, beginning next shake iteration'

   
!    enddo

    ENDSUBROUTINE FullShake
   


END MODULE RotateOrbsMod
