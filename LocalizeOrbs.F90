MODULE LocalizeOrbsMod

    REAL*8 , ALLOCATABLE :: Coeff(:,:)
    REAL*8 , ALLOCATABLE :: Lambdas(:,:)
    REAL*8 , ALLOCATABLE :: DerivCoeff(:,:)
    REAL*8 , ALLOCATABLE :: DerivLambda(:,:)
    REAL*8 , ALLOCATABLE :: OneIndInts(:,:), TwoIndInts(:,:), ThreeIndInts(:,:), FourIndInts(:,:)   !These are the arrays to store the
                                                                                                    !partially transformed integrals for
                                                                                                    !each iteration. The FourIndInts are
                                                                                                    !the full <ij|kl> integrals for the
                                                                                                    !current iteration.

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

        CALL Transform2ElInts()

    END SUBROUTINE RotateOrbs


    SUBROUTINE FindTheForce

!loop over m,z

    !loop over i,j,k,l

        !seperate loop over bgd/agd/abd...

    !loop over j

    !loop over i

    END SUBROUTINE FindTheForce

    SUBROUTINE Transform2EInts()

!Zero arrays from previous transform

!loop over i
    !loop over beta,gamma,delta
        !sum over alpha to store <i beta | gamma delta> - OneIndInts

!loop over j
    !loop over i, gamma, delta
        !sum over beta to store <i j | gamma delta> - TwoIndInts

!loop over k
    !loop over i, j, delta
        !sum over gamma to store <i j | k delta> - ThreeIndInts

!loop over l
    !loop over i,j,k
        !sum over delta to store <i j | k l> - FourIndInts

    END SUBROUTINE Transform2EInts

        

END MODULE LocalizeOrbsMod
