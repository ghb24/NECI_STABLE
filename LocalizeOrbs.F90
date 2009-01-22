MODULE LocalizeOrbsMod

    REAL*8 , ALLOCATABLE :: Coeff(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:)

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

        CALL DeallocAllMem()

    END SUBROUTINE LocalizeOrbs

    SUBROUTINE RotateOrbs()

!Find derivatives of the c and lambda matrices and print the sum of off-diagonal matrix elements.
        CALL FindTheForce()

!Update coefficents by moving them in direction of force. Also update lambda matrix. Print sum of squared changes in coefficients. 
        CALL UseTheForce()

!Force should go to zero as we end in minimum - test for this
        CALL TestForConvergence()

    END SUBROUTINE RotateOrbs


    SUBROUTINE FindTheForce

!loop over m,z

    !loop over i,j,k,l

        !loop over a,b,g,d to find 'new' matrix element

        !seperate loop over bgd/agd/abd...

    !loop over j

    !loop over i

    END SUBROUTINE FindTheForce

        

END MODULE LocalizeOrbsMod
