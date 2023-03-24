#include "macros.h"

module rhodiag_mod

    use HElem, only: HElement_t_size
    use constants, only: dp, int32
    use global_utilities, only: timer, set_timer, halt_timer
    use util_mod, only: near_zero
    use error_handling_neci, only: stop_all
    use blas_interface_mod, only: dsyev, zheev
    better_implicit_none
    private

contains

    FUNCTION RHODIAG_CP(RHOIJ, I_P, I_V)
        INTEGER I_P, I_V
        type(timer), save :: proc_timer
        real(dp) RHOIJ(0:I_V, 0:I_V)
        real(dp) RIJMAT(I_V, I_V), RHODIAG_CP
        real(dp) WLIST(I_V), WORK(3 * I_V)
        INTEGER(int32) :: INFO
        INTEGER I, J
        real(dp) SI
        real(dp) RII
        character(*), parameter :: this_routine = 'RHODIAG_CP'
        RII = RHOIJ(0, 0)
!.. Diagonalize
        proc_timer%timer_name = 'RHODIAG_CP'
        call set_timer(proc_timer)
        RIJMAT(1:I_V, 1:I_V) = 0.0_dp
        DO I = 1, I_V
        DO J = I, I_V
            RIJMAT(I, J) = RHOIJ(I - 1, J - 1) + 0.0_dp
            RIJMAT(I, J) = RIJMAT(I, J) + 0.0_dp
            SI = RHOIJ(I - 1, J - 1)
        END DO
        END DO
!         RIJMAT(1,1)=1.0_dp
!         RIJMAT(1,2)=1.0_dp
!         RIJMAT(2,1)=1.0_dp
!         RIJMAT(2,2)=1.0_dp
!         WRITE(stdout,*) ((RIJMAT(I,J),J=1,I_V),I=1,I_V)
        CALL DSYEV('V', 'U', I_V, RIJMAT(1, 1), I_V, WLIST(1), WORK(1), 3 * I_V, INFO)
!         WRITE(stdout,*) ((RIJMAT(I,J),J=1,I_V),I=1,I_V)
!         WRITE(stdout,*) (WLIST(I),I=1,I_V)
        IF (INFO /= 0) THEN
            WRITE(stdout, *) 'DYSEV error: ', INFO
            call stop_all(this_routine, "DSYEV error")
        END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
        SI = 0.0_dp
!.. divide through by RHOII^P
        DO I = 1, I_V
            SI = SI + RIJMAT(1, I) * RIJMAT(1, I) * ((WLIST(I) / RII)**I_P)
!           WRITE(stdout,"(I5,3G25.16)") I,WLIST(I+1),RIJMAT(I*NLIST+1),WI
        END DO
        RHODIAG_CP = SI
!         DO I=0,NLIST-1
!            DO J=0,NLIST-1
!               WRITE(34,"(2I5)",advance='no') I,J
!               WRITE(34,"(A)",advance='no') "  ("
!               DO K=1,NEL
!                  WRITE(34,"(I3,A)",advance='no') LSTE(K,J),","
!               ENDDO
!            WRITE(34,"(A,G25.16") ")",RIJMAT(I*NLIST+J+1)
!            ENDDO
!         ENDDO

        call halt_timer(proc_timer)
        RETURN
    END

!.. 29/6/06 Based on RHODIAG_CPP, instead of diagonalizing
!.. a matrix of RHOIJ elements, and raising the result ^P,
!.. we diagonalize the HIJ elements and work out e^-beta lambda.

!.. See 1/7/04
    RECURSIVE FUNCTION HDIAG_CPP(HIJ, I_P, I_V, IMISS, TSUB, BETA, DLWDB, HIJS) result(HDIAG_CPPRES)
        INTEGER I_P, I_V
        type(timer), save :: proc_timer
        HElement_t(dp) HIJ(I_V + 1, I_V + 1), RIJMAT(I_V, I_V)
        real(dp) WLIST(I_V), WORK(3 * I_V)
        HElement_t(dp) NWORK(4 * I_V)
        INTEGER(int32) :: INFO
        INTEGER I, J, IMISS, II, IJ
        HElement_t(dp) :: SI, SI2
!.. do we subtract out lower vertices here or later?
        LOGICAL TSUB
        real(dp) BETA
        HElement_t(dp) HIJS(I_V + 1), HIJS2(I_V), DLWT, T, U, HDIAG_CPPRES
        HElement_t(dp) :: DLWDB, R, RII, S, DD2
        character(*), parameter :: this_routine = 'HDIAG_CPP'
! Optimise the 1V case
        IF (I_V == 1) THEN
            DLWDB = HIJ(1, 1)
            HDIAG_CPPRES = 1.0_dp
            RETURN
        END IF
        R = HIJ(1, 1)
        S = -BETA
        R = R * S
        RII = EXP(R)
!.. Diagonalize
!         WRITE(stdout,*) "...",I_V,IMISS
        proc_timer%timer_name = 'HDIAG_CPP '
        call set_timer(proc_timer, 55)
        RIJMAT(1:I_V, 1:I_V) = (0.0_dp)
        DLWDB = 0.0_dp
        II = 0
        DO I = 1, I_V + 1
        IF (I /= IMISS) THEN
            IJ = II
            II = II + 1
            DO J = I, I_V + 1
            IF (J /= IMISS) THEN
                IJ = IJ + 1
                RIJMAT(II, IJ) = HIJ(I, J)
            END IF
            END DO
            HIJS2(II) = HIJS(I)
        END IF
        END DO
        SI = 0.0_dp
!.. Now subtract out the smaller submatrices first
!.. In order to count the subsets only once, we need to only
!.. remove up to IMISS
        IF (TSUB) THEN
        DO I = 2, IMISS - 1
!IMISS-1
            DD2 = 0.0_dp
            SI = SI - HDIAG_CPP(RIJMAT, I_P, I_V - 1, I, TSUB, BETA, DD2, HIJS2)
            DLWDB = DLWDB - DD2
        END DO
        END IF
        IF (HElement_t_size == 1) THEN
            CALL DSYEV('V', 'U', I_V, RIJMAT, I_V, WLIST, WORK, 3 * I_V, INFO)
            IF (INFO /= 0) THEN
                WRITE(stdout, *) 'DYSEV error: ', INFO
                call stop_all(this_routine, "DSYEV error")
            END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
!.. now calculate exp(-beta lambda) for each eigenvalue, with the
!.. appropriate projection onto the root
            SI2 = 0.0_dp
            DLWT = 0.0_dp
            DO I = 1, I_V
!            WRITE(stdout,*) WLIST(I),RIJMAT(1,I)
                R = HIJ(1, 1)
                R = EXP(-BETA * (WLIST(I) - R))
                S = RIJMAT(1, I) * RIJMAT(1, I)
                SI2 = SI2 + S * R
!/RII
                T = R
!/RII
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                DO J = 1, I_V
                    U = HIJS2(J) * RIJMAT(J, I) * RIJMAT(1, I)
                    DLWT = DLWT + U * T
!                 WRITE(stdout,*) I,J,HIJS(J),RIJMAT(J,I),RIJMAT(1,I),U,T,DLWT
                END DO
            END DO
        ELSE
!.. The complex case
            CALL ZHEEV('V', 'U', I_V, RIJMAT, I_V, WLIST, NWORK, 4 * I_V, WORK, INFO)
            IF (INFO /= 0) THEN
                WRITE(stdout, *) 'ZHEEV error: ', INFO
                call stop_all(this_routine, "ZHEEV error")
            END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
!.. now calculate exp(-beta lambda) for each eigenvalue, with the
!.. appropriate projection onto the root
            SI2 = 0.0_dp
            DLWT = 0.0_dp
            DO I = 1, I_V
!            WRITE(stdout,*) WLIST(I),RIJMAT(1,I)
                S = abs(RIJMAT(1, I))**2
                R = HIJ(1, 1)
                R = EXP(-BETA * (WLIST(I) - R))
                SI2 = SI2 + S * R
!%/RII
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                U = R
!/RII
                DO J = 1, I_V
#ifdef CMPLX_
                    T = HIJS2(J) * RIJMAT(J, I) * conjg(RIJMAT(1, I))
#else
                    T = HIJS2(J) * RIJMAT(J, I) * (RIJMAT(1, I))
#endif
                    DLWT = DLWT + T * U
!                 WRITE(stdout,*) I,J,HIJS(J),RIJMAT(J,I),RIJMAT(1,I),T,U,DLWT
                END DO
            END DO
        END IF
        S = DLWT
        DLWDB = DLWDB + S
        HDIAG_CPPRES = SI + SI2
        call halt_timer(proc_timer)
        RETURN
    END
!.. See 1/7/04
    RECURSIVE FUNCTION RHODIAG_CPP(RHOIJ, I_P, I_V, IMISS, TSUB, DBETA, DLWDB, HIJS, tLogWeight) RESULT(RHODIAG_CPPRES)
        INTEGER I_P, I_V
        type(timer), save :: proc_timer
        HElement_t(dp) RHOIJ(I_V + 1, I_V + 1), HIJS(I_V + 1), HIJS2(I_V)
        HElement_t(dp) RIJMAT(I_V, I_V), NWORK(4 * I_V), S, DLWT, S2
        real(dp) WORK(3 * I_V)
        HElement_t(dp) :: R
        real(dp) WLIST(I_V)
        INTEGER(int32) :: INFO
        INTEGER I, J, IMISS, II, IJ
        real(dp) SS2
        HElement_t(dp) :: SI, SS, SI2
!.. do we subtract out lower vertices here or later?
        LOGICAL TSUB
        HElement_t(dp) :: DLWDB, RII, DD2
        HElement_t(dp) RhoDiag_CPPRES
        real(dp) DBETA
        LOGICAL tLogWeight
        character(*), parameter :: this_routine = 'RHODIAG_CPP'
! Optimise the 1V case
        IF (I_V == 1) THEN
            DLWDB = HIJS(1)
            if (tLogWeight) then
                RHODIAG_CPPRES = 0.0_dp
            else
                RHODIAG_CPPRES = 1.0_dp
            END IF
            RETURN
        END IF
        RII = RHOIJ(1, 1)
!.. Diagonalize
!         WRITE(stdout,*) "...",I_V,IMISS
        proc_timer%timer_name = 'RHODIAG_C2'
        call set_timer(proc_timer, 55)
        RIJMAT(1:I_V, 1:I_V) = (0.0_dp)
        IF (.not. near_zero(DBETA)) DLWDB = 0.0_dp
        II = 0
        DO I = 1, I_V + 1
        IF (I /= IMISS) THEN
            IJ = II
            II = II + 1
            DO J = I, I_V + 1
            IF (J /= IMISS) THEN
                IJ = IJ + 1
                RIJMAT(II, IJ) = RHOIJ(I, J)
            END IF
            END DO
            HIJS2(II) = HIJS(I)
        END IF
        END DO
        SI = 0.0_dp
!.. Now subtract out the smaller submatrices first
!.. In order to count the subsets only once, we need to only
!.. remove up to IMISS
        IF (TSUB) THEN
        DO I = 2, IMISS - 1
            DD2 = 0.0_dp
            IF (tLogWeight) THEN
! we return e~' instead of w' E~'
! e~'[G]=e~[G]-sum_{g in G} e~'[g]
!  we ignore w' returned in SI
!  DD2 returns e~'[g]
!  where g is this graph missing out vertex I.
                SI = SI - RHODIAG_CPP(RIJMAT, I_P, I_V - 1, I, TSUB, DBETA, DD2, HIJS2, tLogWeight)
                IF (.not. near_zero(DBETA)) DLWDB = DLWDB - DD2
            ELSE
! we return w' E~'
! w'[G]E~'[G]=w[G]E~[G]-sum_{g in G} w'[g] E~'[g]
!SI is w' and DLWDB is E~'
!  The next line subtracts out the contribution from this graph, but missing out vertex I.
!   done recursively, this makes w'[G]E~'[G]
                SI = SI - RHODIAG_CPP(RIJMAT, I_P, I_V - 1, I, TSUB, DBETA, DD2, HIJS2, tLogWeight)
                IF (.not. near_zero(DBETA)) DLWDB = DLWDB - DD2
            END IF
!            WRITE(stdout,*) "SI=",SI
        END DO
        END IF
!         WRITE(stdout,*) I_V
        IF (HElement_t_size == 1) THEN
            CALL DSYEV('V', 'U', I_V, RIJMAT, I_V, WLIST, WORK, 3 * I_V, INFO)
            IF (INFO /= 0) THEN
                WRITE(stdout, *) 'DYSEV error: ', INFO
                call stop_all(this_routine, "DSYEV error")
            END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
!.. divide through by RHOII^P
            DLWT = 0.0_dp
            SI2 = 0.0_dp
            DO I = 1, I_V
!               WRITE(stdout,*) WLIST(I),RIJMAT(1,I)
                R = (RII)
                R = ((WLIST(I) / R)**I_P)
                S = R
                SS = S * RIJMAT(1, I) * RIJMAT(1, I)
                SI2 = SI2 + SS
                IF (.not. near_zero(DBETA)) THEN
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                    DO J = 1, I_V
                        S2 = S * RIJMAT(J, I) * RIJMAT(1, I)
                        DLWT = DLWT + S2 * HIJS2(J)
                    END DO
                END IF
            END DO
! DLWT is the value of w[G] E~[G]
! SI2 is w[G]
!DLWDB contains the subtracted out subgraphs
            if (tLogWeight) then
! we return e~' instead of w' E~'
! e~'[G]=e~[G]-sum_{g in G} e~'[g]
                S = DLWT
                S2 = SI2
                S = S / S2
                SS = S
                DLWDB = SS + DLWDB
            else
                SS = DLWT
                DLWDB = DLWDB + SS
            end if
        ELSE
!.. The complex case
            CALL ZHEEV('V', 'U', I_V, RIJMAT, I_V, WLIST, NWORK, 4 * I_V, WORK, INFO)
            IF (INFO /= 0) THEN
                WRITE(stdout, *) 'ZHEEV error: ', INFO
                call stop_all(this_routine, "ZHEEV error")
            END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
!.. divide through by RHOII^P
            SI2 = 0.0_dp
            DLWT = 0.0_dp
            DO I = 1, I_V
                R = (RII)
                SS = ((WLIST(I) / R)**I_P)
                SS2 = abs(RIJMAT(1, I))**2
                SI2 = SI2 + SS * SS2
                IF (.not. near_zero(DBETA)) THEN
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                    DO J = 1, I_V
                        S = SS
#ifdef CMPLX_
                        S = S * RIJMAT(J, I) * conjg(RIJMAT(1, I))
#else
                        S = S * RIJMAT(J, I) * (RIJMAT(1, I))
#endif
                        DLWT = DLWT + S * HIJS2(J)
                    END DO
                END IF
            END DO
! DLWT is the value of w[G] E~[G]
! SI2 is w[G]
!DLWDB contains the subtracted out subgraphs
            if (tLogWeight) then
! we return e~' instead of w' E~'
! e~'[G]=e~[G]-sum_{g in G} e~'[g]

                S = DLWT
                S2 = SI2
                S = S / S2
                SS = S
                DLWDB = SS + DLWDB
            else
                SS = DLWT
                DLWDB = DLWDB + SS
            end if
        END IF

        if (tLogWeight) then
!  We pass the log x[G] around
            SI2 = LOG(SI2)
            RHODIAG_CPPRES = SI + SI2
!  We return zero as the weight of the graph
        ELSE
            RHODIAG_CPPRES = SI + SI2
        END IF

        call halt_timer(proc_timer)
        RETURN
    END
!  As rhodiag_cpp, but only deal with this vertex level, and subtract out the two-vertex star contribtion from this vertex level.
    FUNCTION RHODIAG_CPPS2VS(RHOIJ, I_P, I_V, DBETA, DLWDB, HIJS)
        INTEGER I_P, I_V
        type(timer), save :: proc_timer
        HElement_t(dp) RHOIJ(I_V + 1, I_V + 1), HIJS(I_V + 1), HIJS2(I_V)
        HElement_t(dp) RIJMAT(I_V, I_V), NWORK(4 * I_V), S, DLWT, S2
        real(dp) WORK(3 * I_V)
        HElement_t(dp) :: R
        real(dp) WLIST(I_V)
        INTEGER(int32) :: INFO
        INTEGER I, J
        real(dp) SI, SS2
        HElement_t(dp) :: SS, SI2
!.. do we subtract out lower vertices here or later?
        HElement_t(dp) DLWDB, RII
        HElement_t(dp) RhoDiag_CPPS2VS
        real(dp) DBETA
        character(*), parameter :: this_routine = 'RHODIAG_CPPS2VS'
! Optimise the 1V case
        IF (I_V == 1) THEN
            RII = HIJS(1)
            DLWDB = DLWDB - RII
            RHODIAG_CPPS2VS = -1.0_dp
            RETURN
        END IF
        RII = RHOIJ(1, 1)
!.. Diagonalize
!         WRITE(stdout,*) "...",I_V,IMISS
        proc_timer%timer_name = 'RHODIAG_C2'
        call set_timer(proc_timer, 55)
        RIJMAT(1:I_V, 1:I_V) = (0.0_dp)
        DO J = 1, I_V
            RIJMAT(1, J) = RHOIJ(1, J)
            RIJMAT(J, J) = RHOIJ(J, J)
            HIJS2(J) = HIJS(J)
        END DO
        SI = 0.0_dp
!  No need to deal with small submatrices
        IF (HElement_t_size == 1) THEN
            CALL DSYEV('V', 'U', I_V, RIJMAT, I_V, WLIST, WORK, 3 * I_V, INFO)
            IF (INFO /= 0) THEN
                WRITE(stdout, *) 'DYSEV error: ', INFO
                call stop_all(this_routine, "DSYEV error")
            END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
!.. divide through by RHOII^P
            DLWT = 0.0_dp
            SI2 = 0.0_dp
            DO I = 1, I_V
!               WRITE(stdout,*) WLIST(I),RIJMAT(1,I)
                R = (RII)
                R = ((WLIST(I) / R)**I_P)
                S = R
                SS = S * RIJMAT(1, I) * RIJMAT(1, I)
                SI2 = SI2 + SS
                IF (.not. near_zero(DBETA)) THEN
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                    DO J = 1, I_V
                        S2 = S * RIJMAT(J, I) * RIJMAT(1, I)
                        DLWT = DLWT + S2 * HIJS2(J)
                    END DO
                END IF
            END DO
!            S=DLWDB
!            S=S+DLWT
            SS = DLWT
!We're doing a subtraction
            DLWDB = DLWDB - SS
!            WRITE(stdout,*) I_V,DLWDB
        ELSE
!.. The complex case
            CALL ZHEEV('V', 'U', I_V, RIJMAT, I_V, WLIST, NWORK, 4 * I_V, WORK, INFO)
            IF (INFO /= 0) THEN
                WRITE(stdout, *) 'ZHEEV error: ', INFO
                call stop_all(this_routine, "ZHEEV error")
            END IF
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues
!.. divide through by RHOII^P
            SI2 = 0.0_dp
            DLWT = 0.0_dp
            DO I = 1, I_V
                R = (RII)
                SS = ((WLIST(I) / R)**I_P)
                SS2 = abs(RIJMAT(1, I))**2
                SI2 = SI2 + SS * SS2
                IF (.not. near_zero(DBETA)) THEN
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
                    DO J = 1, I_V
                        S = SS
#ifdef CMPLX_
                        S = S * RIJMAT(J, I) * conjg(RIJMAT(1, I))
#else
                        S = S * RIJMAT(J, I) * (RIJMAT(1, I))
#endif
                        DLWT = DLWT + S * HIJS2(J)
                    END DO
                END IF
            END DO
            S = DLWDB
! In subtraction mode
            S = S - DLWT
            DLWDB = S
        END IF
        RHODIAG_CPPS2VS = SI - SI2
        call halt_timer(proc_timer)
        RETURN
    END

end module
