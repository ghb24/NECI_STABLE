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

        

END MODULE LocalizeOrbsMod
