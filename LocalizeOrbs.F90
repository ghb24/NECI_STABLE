MODULE LocalizeOrbsMod

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

        

END MODULE LocalizeOrbsMod
