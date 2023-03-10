INTEGER FUNCTION ICMPDETS(D1, D2, NEL)
    IMPLICIT NONE
    INTEGER NEL, I
    INTEGER D1(1:NEL), D2(1:NEL)
    DO I = 1, NEL
        IF (D1(I) < D2(I)) THEN
            ICMPDETS = -1
            RETURN
        ELSE IF (D1(I) > D2(I)) THEN
            ICMPDETS = 1
            RETURN
        END IF
    END DO
    ICMPDETS = 0
    RETURN
END

INTEGER FUNCTION IGETEXCITLEVEL(NI, NJ, NEL)
    IMPLICIT NONE
    INTEGER I, J, NEL, NI(NEL), NJ(NEL), IC
    IC = 0
    DO I = 1, NEL
        DO J = 1, NEL
            IF (NI(I) == NJ(J)) THEN
                IC = IC + 1
                EXIT
            END IF
        END DO
    END DO
    IGETEXCITLEVEL = NEL - IC
    RETURN
END

!This routine is the same as IGetExcitLevel, but for two differences.
!First, this will only work for ordered lists of occupied orbitals
!Secondly, this will only find the number of excitation levels apart up
!to a maximum excitation level of MaxExcit. If the distance away is more than this
!it will simply return MaxExcit+1. i.e. if only want to know if it is a double
!excitation or less, then CALL IGetExcitLevel_2(NI,NJ,NEL,2). This will return
!0, 1, 2, or 3 if the distance is more than a double excitation.
INTEGER FUNCTION IGetExcitLevel_2(NI, NJ, NEL, MaxExcit)
    IMPLICIT NONE
    INTEGER :: NEL, NI(NEl), NJ(NEl), MaxExcit, I, IC, j, MaxIC
    LOGICAL :: CommonOrb

    IF (MaxExcit >= NEl) THEN
        MaxIC = NEl - 1
    ELSE
        MaxIC = MaxExcit
    END IF
    I = 1   !Index for electron in NI
    IC = 0  !Number of excitation levels between determinants
    DO WHILE ((IC <= MaxIC) .and. (I <= NEl))
        CommonOrb = .false.
        do j = 1, NEL
        IF (NJ(j) == NI(I)) THEN
            CommonOrb = .true.
            EXIT
        ELSEIF (NJ(j) > NI(I)) THEN
            EXIT
        END IF
        end do
        IF (.not. CommonOrb) THEN
            IC = IC + 1
        END IF
        I = I + 1
    END DO
    IGetExcitLevel_2 = IC
    RETURN
END

!.. Very unnecessarily complex, but would be faster for large numbers of
!.. electrons.

INTEGER FUNCTION IGETEXCITLEVEL_(NI, NJ, NEL)
    IMPLICIT NONE
    INTEGER I, J, NEL, NI(NEL), NJ(NEL), IC
    IC = 0
    I = 1
    J = 1
!.. We only count differences from I to J
    DO WHILE (I <= NEL .AND. J <= NEL)
    DO WHILE (NI(I) < NJ(J) .AND. I <= NEL)
        I = I + 1
        IC = IC + 1
    END DO
    DO WHILE (NI(I) > NJ(J) .AND. I <= NEL .AND. J <= NEL)
        J = J + 1
!               IC=IC+1
    END DO
    IF (NI(I) == NJ(J)) THEN
        I = I + 1
        J = J + 1
    END IF
    END DO
    IC = IC + (NEL + 1 - I)
    IGETEXCITLEVEL_ = IC
    RETURN
END

!.. GETHELEMENT
!.. Get matrix element of the hamiltonian
FUNCTION GETHELEMENT(II, IJ, HAMIL, LAB, NROW, NDET)
    use constants, only: dp
    IMPLICIT NONE
    real(dp) HAMIL(*), GETHELEMENT
    INTEGER LAB(*)
    INTEGER NDET, NROW(NDET), IJ, II, I, J, INDXROW, IMAX, K
!.. We only have half of H, so if J<I, return the symmetrical (J,I) element
!.. Or if we have the whole H, it's quicker to look closer to its beginning
    IF (IJ < II) THEN
        I = IJ
        J = II
    ELSE
        I = II
        J = IJ
    END IF
    GETHELEMENT = 0.0_dp
    INDXROW = 1
!.. Find the Ith row
    DO K = 1, I - 1
        INDXROW = INDXROW + NROW(K)
    END DO
    IMAX = INDXROW + NROW(I) - 1
    DO K = INDXROW, IMAX
        IF (LAB(K) > J) RETURN
        IF (LAB(K) == J) THEN
            GETHELEMENT = HAMIL(K)
            RETURN
        END IF
    END DO
    RETURN
END
