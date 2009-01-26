!TO DO
!Spatial symmetry
!Finalize orbs to fully test
!Permutational symmetry where possible
!Parallelize

MODULE LocalizeOrbsMod

    USE Global_utilities
    USE IntegralsData , only : UMAT
    USE UMatCache , only : UMatInd
    USE SystemData , only : ConvergedForce,TimeStep
    IMPLICIT NONE
    REAL*8 , ALLOCATABLE :: Coeff(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:)
    REAL*8 , ALLOCATABLE :: OneIndInts(:,:,:,:), TwoIndInts(:,:,:,:), ThreeIndInts(:,:,:,:), FourIndInts(:,:,:,:)   !These are the arrays to store the
!partially transformed integrals for each iteration. The FourIndInts are the full <ij|kl> integrals for the current iteration.
!OneIndInts = <i b|g d> ; TwoIndInts = <i b|k d> ; ThreeIndInts = <i j|k d>
    INTEGER :: SpatOrbs,OneIndIntsTag,TwoIndIntsTag,ThreeIndIntsTag,FourIndIntsTag
    INTEGER :: CoeffTag,LambdasTag,DerivCoeffTag,DerivLambdaTag,Iteration
    LOGICAL :: tNotConverged
    REAL*8 :: OrthoNorm,PotEnergy,Force,TwoEInts,DistCs,OrthoForce,DistLs,LambdaMag
    REAL*8 :: OrthoFac=1.D0

    contains

    SUBROUTINE LocalizeOrbs()

        CALL InitLocalOrbs()

        WRITE(6,"(I7,8F25.10)") Iteration,PotEnergy,Force,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
        WRITE(12,"(I7,8F25.10)") Iteration,PotEnergy,Force,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
        tNotConverged=.true.
        do while(tNotConverged)

            Iteration=Iteration+1

            CALL RotateOrbs() 

            WRITE(6,"(I7,8F25.10)") Iteration,PotEnergy,Force,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
            WRITE(12,"(I7,8F25.10)") Iteration,PotEnergy,Force,OrthoForce,TwoEInts,OrthoNorm,DistCs,DistLs,LambdaMag
            CALL FLUSH(6)
            CALL FLUSH(12)

        enddo

        WRITE(6,*) "Convergence criterion met. Finalizing new orbitals..."

!Make symmetry, orbitals, one/two-electron integrals consistent with rest of neci
!        CALL FinalizeNewOrbs()

!option to create new FCIDUMP?

        CALL DeallocateMem()

        CLOSE(12)

        CALL Stop_All("LocalizeOrbs","This code is still in the testing phase")

    END SUBROUTINE LocalizeOrbs

    SUBROUTINE RotateOrbs()

!Find derivatives of the c and lambda matrices and print the sum of off-diagonal matrix elements.
        CALL FindTheForce()

!Update coefficents by moving them in direction of force. Also update lambda matrix. Print sum of squared changes in coefficients. 
        CALL UseTheForce()

        CALL TestOrthonormality()

!Force should go to zero as we end in minimum - test for this
        CALL TestForConvergence()

        CALL Transform2ElInts()

    END SUBROUTINE RotateOrbs

!This just tests the convergence on the grounds that the force is smaller that the input parameter: ConvergedForce
    SUBROUTINE TestForConvergence()

        IF(abs(Force).lt.ConvergedForce) THEN
            tNotConverged=.false.
        ENDIF

    END SUBROUTINE TestForConvergence

    SUBROUTINE TestOrthonormality()
        INTEGER :: i,j,a

        OrthoNorm=0.D0
        do i=1,SpatOrbs
            do j=1,i
                do a=1,SpatOrbs
                    OrthoNorm=OrthoNorm+Coeff(i,a)*Coeff(a,j)
                enddo
            enddo
        enddo
        OrthoNorm=OrthoNorm-real(SpatOrbs,8)
        OrthoNorm=(OrthoNorm*2.D0)/REAL((SpatOrbs*(SpatOrbs+1.D0)),8)

    END SUBROUTINE TestOrthonormality

    SUBROUTINE FindTheForce
        INTEGER :: m,z,i,j,k,l,a,b,g,d
        REAL*8 :: Deriv4indint,Deriv4indintsqrd,t1,t2,t3,LambdaTerm1,LambdaTerm2

! Running over m and z, covers all matrix elements of the force matrix (derivative 
! of equation we are minimising, with respect to each translation coefficient) filling 
! them in as it goes.

    DerivCoeff(:,:)=0.D0
    Force=0.D0
    Deriv4indintsqrd=0.D0
    
    do m=1,SpatOrbs
! To include: symmetry requirement that z must be from the same irrep as m
        do z=1,SpatOrbs
                
! This runs over the full <ij|kl> integrals from the previous iteration. 
! In the future we can take advantage of the permutational symmetry of the matrix elements
! however this only takes into account 1/4 of the off-diagonal plus the diagonal elements
! and this must be accounted for when calculating Deriv4indintsqrd.
            do i=1,SpatOrbs
                do k=1,SpatOrbs
                    do j=1,SpatOrbs
                        do l=1,SpatOrbs

! To include: If statement to check the symmetry of i x j x k x l is A1.
!                                If((((i-1)*nBasis)+k).ge.(((j-1)*nBasis)+l)) THEN

! Already calculated the four index integrals and some partial integrals
! in Transform2ElInts routine of the previous iteration.

! Deriv4indint is the derivative of <ij|kl> with respect to a coefficient
! c(m,z), at the current m,z,i,j,k,l of the loop.
! This involves 4 terms in which the derivative is taken w.r.t the i,j,k or l position in <ij|kl>. 

                            t1=0.D0    
                            do b=1,SpatOrbs
                                do g=1,SpatOrbs
                                    do d=1,SpatOrbs
                                        t1=t1+(Coeff(j,b)*Coeff(k,g)*Coeff(l,d)*REAL(UMAT(UMatInd(z,b,g,d,0,0))%v,8))
                                    enddo
                                enddo
                            enddo
                            
                            t2=0.D0    
                            do d=1,SpatOrbs     
                                t2=t2+(Coeff(l,d)*TwoIndInts(i,z,k,d))
                            enddo

                            t3=0.D0
                            do b=1,SpatOrbs
                                do d=1,SpatOrbs
                                    t3=t3+Coeff(j,b)*Coeff(l,d)*OneIndInts(i,b,z,d)
                                enddo
                            enddo

                            Deriv4indint=(t1 + t2 + t3 + ThreeIndInts(i,j,k,z))

! Deriv4indintsqrd is the derivative of the overall expression for the sum of the squares of the <ij|kl> matrix.
! This accumulates as the loop sums over i,j,k and l.
!                            IF(Iteration.lt.1000) THEN
                                Deriv4indintsqrd=Deriv4indintsqrd+(FourIndInts(i,j,k,l)*Deriv4indint)   
!                            ELSE
!                                Deriv4indintsqrd=0.D0
!                            ENDIF

!                                endif
                        enddo
                    enddo
                enddo
            enddo

! Calculate the derivatives of orthogonalisation condition.

            LambdaTerm1=0.D0
            LambdaTerm2=0.D0
            do j=1,SpatOrbs
                LambdaTerm1=LambdaTerm1+(Lambdas(m,j)*Coeff(z,j))
            enddo

            do i=1,SpatOrbs
                LambdaTerm2=LambdaTerm2+(Lambdas(i,m)*Coeff(i,z))
            enddo

! DerivCoeff is 'the force'.  I.e. the derivative of |<ij|kl>|^2 (with orthogonalisation constraints) with 
! respect to each transformation coefficient.  It is the values of this matrix that will tend to 0 as
! we minimise the sum of the |<ij|kl>|^2 values.
!            IF(Iteration.lt.1000) THEN
                DerivCoeff(m,z)=(2*Deriv4indintsqrd)-OrthoFac*(LambdaTerm1-LambdaTerm2)
!            ELSE
!                DerivCoeff(m,z)=-LambdaTerm1-LambdaTerm2
!            ENDIF
            Force=Force+DerivCoeff(m,z)

        enddo
    enddo

    Force=Force/REAL(SpatOrbs**2,8)
    OrthoForce=(LambdaTerm1+LambdaTerm2)/REAL(SpatOrbs**2,8)

!We also need to find the force on the lambdas to ensure orthonormality...
    DerivLambda(:,:)=0.D0
    do i=1,SpatOrbs
        do j=1,i
            do a=1,SpatOrbs
                DerivLambda(i,j)=DerivLambda(i,j)+Coeff(i,a)*Coeff(a,j)
            enddo
            DerivLambda(j,i)=DerivLambda(i,j)
        enddo
    enddo
    do i=1,SpatOrbs
        DerivLambda(i,i)=DerivLambda(i,i)-1.D0
    enddo


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
            NewCoeff=Coeff(m,z)-(TimeStep*DerivCoeff(m,z))
            DistCs=DistCs+abs(TimeStep*DerivCoeff(m,z))
            Coeff(m,z)=NewCoeff
        enddo
    enddo
    DistCs=DistCs/(REAL(SpatOrbs**2,8))

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
                            t=t+Coeff(i,a)*REAL(UMAT(UMatInd(a,b,g,d,0,0))%v,8)
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
                            t=t+Coeff(k,g)*OneIndInts(i,b,g,d)
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
                            t=t+Coeff(j,b)*TwoIndInts(i,b,k,d)
                        enddo
                        ThreeIndInts(i,j,k,d)=t
                        ThreeIndInts(k,j,i,d)=t
                    enddo
                enddo
            enddo
        enddo

        PotEnergy=0.D0
        TwoEInts=0.D0
        do l=1,SpatOrbs !loop over l orbitals in each symmetry irrep
            do j=1,l    ! j =< l
                do i=1,SpatOrbs
                    do k=1,i    ! k =< i
                        !Test that (i,k) =< (j,l) here
                        !Test that integral is symmetry allowed
                        t=0.D0
!sum over delta to store <i j | k l> - FourIndInts
                        do d=1,SpatOrbs
                            t=t+Coeff(l,d)*ThreeIndInts(i,j,k,d)
                        enddo
                        FourIndInts(i,j,k,l)=t
                        FourIndInts(k,j,i,l)=t
                        FourIndInts(i,l,k,j)=t
                        FourIndInts(k,l,i,j)=t
                        PotEnergy=PotEnergy+(t**2)
                        TwoEInts=TwoEInts+(t**2)
                    enddo
                enddo
            enddo
        enddo

! Now find the change of the potential energy due to the orthonormality of the orbitals...
        do i=1,SpatOrbs
            do j=1,SpatOrbs
                t=0.D0
                do a=1,SpatOrbs
                    t=Coeff(i,a)*Coeff(a,j)
                enddo
                IF(i.eq.j) t=t-1.D0
                PotEnergy=PotEnergy-Lambdas(i,j)*t
            enddo
        enddo

    END SUBROUTINE Transform2ElInts


    SUBROUTINE InitLocalOrbs()
        USE SystemData , only : nBasis
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        REAL*8 :: t
        INTEGER :: i,j,k,l,ierr

        WRITE(6,*) "Calculating new molecular orbitals based on mimimisation of <ij|kl>^2 integrals..."

        OrthoNorm=0.D0
        PotEnergy=0.D0
        Force=0.D0
        TwoEInts=0.D0
        DistCs=0.D0
        DistLs=0.D0
        LambdaMag=0.D0
        SpatOrbs=nBasis/2
        Iteration=0
        OrthoForce=0.D0

!Allocate memory
        ALLOCATE(Coeff(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('Coeff',SpatOrbs**2,8,this_routine,CoeffTag,ierr)
        ALLOCATE(Lambdas(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('Lambdas',SpatOrbs**2,8,this_routine,LambdasTag,ierr)
        ALLOCATE(DerivCoeff(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('DerivCoeff',SpatOrbs**2,8,this_routine,DerivCoeffTag,ierr)
        ALLOCATE(DerivLambda(SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('DerivLambda',SpatOrbs**2,8,this_routine,DerivLambdaTag,ierr)
        ALLOCATE(OneIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('OneIndInts',SpatOrbs**4,8,this_routine,OneIndIntsTag,ierr)
        ALLOCATE(TwoIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('TwoIndInts',SpatOrbs**4,8,this_routine,TwoIndIntsTag,ierr)
        ALLOCATE(ThreeIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('ThreeIndInts',SpatOrbs**4,8,this_routine,ThreeIndIntsTag,ierr)
        ALLOCATE(FourIndInts(SpatOrbs,SpatOrbs,SpatOrbs,SpatOrbs),stat=ierr)
        CALL LogMemAlloc('FourIndInts',SpatOrbs**4,8,this_routine,FourIndIntsTag,ierr)

!Zero/initialise the arrays
        Coeff(:,:)=0.D0
        do i=1,SpatOrbs
            Coeff(i,i)=1.D0
        enddo
        Lambdas(:,:)=0.D0
        DerivCoeff(:,:)=0.D0
        DerivLambda(:,:)=0.D0
!These loops can be speeded up with spatial symmetry and pairwise permutation symmetry if needed.
!Fill all the original partially transformed integral arrays with the original two-electron integrals.
        do i=1,SpatOrbs
            do k=1,i
                do j=1,SpatOrbs
                    do l=1,j
                        t=REAL(UMAT(UMatInd(i,j,k,l,0,0))%v,8)
                        PotEnergy=PotEnergy+(t**2)       !Potential energy starts as this since the orbitals are orthonormal by construction.
                        TwoEInts=TwoEInts+(t**2)
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
        CALL TestOrthonormality()

        OPEN(12,FILE='Transform',STATUS='unknown')
!We want to write out: Iteration, Potential energy, Force, Sum<ij|kl>^2, OrthonormalityCondition
        WRITE(12,"(A)") "# Iteration   PotEnergy   Force   OrthoForce    Sum<ij|kl>^2   OrthoNormCondition   DistMovedbyCs   DistMovedByLs   LambdaMag"
        WRITE(6,"(A)") "Iteration   PotEnergy   Force   OrthoForce    Sum<ij|kl>^2   OrthoNormCondition   DistMovedbyCs   DistMovedbyLs   LambdaMag"

    END SUBROUTINE InitLocalOrbs

    SUBROUTINE DeallocateMem()
        CHARACTER(len=*) , PARAMETER :: this_routine='DeallocateMem'

        DEALLOCATE(Coeff)
        CALL LogMemDealloc(this_routine,CoeffTag)
        DEALLOCATE(Lambdas)
        CALL LogMemDealloc(this_routine,LambdasTag)
        DEALLOCATE(DerivCoeff)
        CALL LogMemDealloc(this_routine,DerivCoeffTag)
        DEALLOCATE(DerivLambda)
        CALL LogMemDealloc(this_routine,DerivLambdaTag)
        DEALLOCATE(OneIndInts)
        CALL LogMemDealloc(this_routine,OneIndIntsTag)
        DEALLOCATE(TwoIndInts)
        CALL LogMemDealloc(this_routine,TwoIndIntsTag)
        DEALLOCATE(ThreeIndInts)
        CALL LogMemDealloc(this_routine,ThreeIndIntsTag)
        DEALLOCATE(FourIndInts)
        CALL LogMemDealloc(this_routine,FourIndIntsTag)

    END SUBROUTINE DeallocateMem
        

END MODULE LocalizeOrbsMod
