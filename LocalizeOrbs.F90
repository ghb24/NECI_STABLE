MODULE LocalizeOrbsMod

    REAL*8 , ALLOCATABLE :: Coeff(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:)
    REAL*8 , ALLOCATABLE :: OneIndInts(:,:,:,:), TwoIndInts(:,:,:,:), ThreeIndInts(:,:,:,:), FourIndInts(:,:,:,:)   !These are the arrays to store the
!partially transformed integrals for each iteration. The FourIndInts are the full <ij|kl> integrals for the current iteration.
!OneIndInts = <i b|g d> ; TwoIndInts = <i b|k d> ; ThreeIndInts = <i j|k d>
    INTEGER :: SpatOrbs,OneIndIntsTag,TwoIndIntsTag,ThreeIndIntsTag,FourIndIntsTag
    INTEGER :: CoeffTag,LambdasTag,DerivCoeffTag,DerivLambdaTag
    LOGICAL :: tNotConverged

    contains

    SUBROUTINE LocalizeOrbs()

        CALL InitLocalOrbs()

        tNotConverged=.true.
        do while(tNotConverged)

            CALL RotateOrbs() 

        enddo

!Make symmetry, orbitals, one/two-electron integrals consistent with rest of neci
        CALL FinalizeNewOrbs()

!option to create new FCIDUMP?

        CALL DeallocateMem()

        CLOSE(12)

    END SUBROUTINE LocalizeOrbs

    SUBROUTINE RotateOrbs()

!Find derivatives of the c and lambda matrices and print the sum of off-diagonal matrix elements.
        CALL FindTheForce()

!Update coefficents by moving them in direction of force. Also update lambda matrix. Print sum of squared changes in coefficients. 
        CALL UseTheForce()

!Force should go to zero as we end in minimum - test for this
        CALL TestForConvergence()

        CALL Transform2ElInts()

    END SUBROUTINE RotateOrbs


    SUBROUTINE FindTheForce

!loop over m,z

    !loop over i,j,k,l

        !seperate loop over bgd/agd/abd...

    !loop over j

    !loop over i

    END SUBROUTINE FindTheForce

    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.
    SUBROUTINE Transform2EInts()
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
                            t=t+Coeff(i,a)*UMAT(UMatInd(a,b,g,d,0,0))
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
                    enddo
                enddo
            enddo
        enddo

    END SUBROUTINE Transform2EInts


    SUBROUTINE InitLocalOrbs()
        USE SystemData , only : nBasis
        CHARACTER(len=*) , PARAMETER :: this_routine='InitLocalOrbs'
        REAL*8 :: t
        INTEGER :: i,j,k,l,ierr

        SpatOrbs=nBasis/2

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
                        t=UMAT(UMatInd(i,j,k,l,0,0))
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

        OPEN(12,FILE='Transform.Dat',STATUS='unknown')

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
