module frsblk_mod
    use matmul_mod, only: my_hpsi
    implicit none
    private
    public :: NECI_FRSBLKH, NECI_WRITE_MATRIX

contains

    SUBROUTINE NECI_FRSBLKH(M, ICMAX, N, H, LAB, V0, VS, NKRY, NKRY1, NBLOCK, NROW, &
                            LSCR, LISCR, A, W, V, AM, BM, T, WT, SCR, ISCR, INDEX, &
                            NCYCLE, B2LIMIT, PRINTOUT, TLargest, tDie)
        use constants, only: dp

        IMPLICIT NONE
    ! ==----------------------------------------------------------------==
    ! ==    M is NDET : the number of rows in Hamil                     ==
    ! ==    ICMAX : the MAX. number of columns in Hamil                 ==
    ! ==    N is NEVAL : the number of columns in the wavevector        ==
    ! ==    PRINTOUT is a logical testing whether to print results      ==
    ! ==    tDie indicates whether STOP should be called if it fails
    ! ==    to diagonalise
    ! ==----------------------------------------------------------------==

        LOGICAL :: PRINTOUT, TLargest, tDie2, tFail
        LOGICAL :: tDie
        integer :: LISCR, ISCR(LISCR), LSCR, NBLOCK, NKRY1, NKRY, N, ICMAX, M
        REAL(dp) :: V0(M, N), VS(M, N), A(N, N), W(N)
        REAL(dp) :: V(M * NBLOCK * NKRY1)
        REAL(dp) :: AM(NBLOCK * NBLOCK * NKRY1), BM(NBLOCK * NBLOCK * NKRY)
        REAL(dp) :: T(3 * NBLOCK * NKRY * NBLOCK * NKRY), WT(NBLOCK * NKRY)
        real(dp) :: SCR(LSCR), B2LIMIT
        integer :: NCYCLE, NDIAG
        integer :: INDEX(N)
    !..temp
        real(dp) :: H(M, ICMAX), ReturnNan
    !..temp
        integer :: LAB(M, ICMAX), NROW(M)
    ! ==----------------------------------------------------------------==
        tDie2 = tDie
        ReturnNan = -1.0_dp
        tFail = .false.

    !..Initialise V0
        CALL NECI_SETUP_MATRIX(M, N, V0, .FALSE.)
        CALL NECI_MGS(M, N, V0, M, A, N, tDie2, tFail)
        if (tFail) then
            W(1) = sqrt(ReturnNan)
            return
        end if

    !..Calculate top 80% of N
        NDIAG = N - INT(N * 0.2_dp)
    !..
        IF (PRINTOUT) THEN
            WRITE(stdout, 20000) M, N, NKRY, NBLOCK, NDIAG, LSCR, B2LIMIT
        END IF
    20000 FORMAT(2X, 'M:', I7 / 2X, 'N:', I7 / 2X, 'NKRY:', I7 / 2X, 'NBLOCK:', I7 / 2X, 'NDIAG:', I7 / 2X, 'LSCR:', I7 / 2X, 'B2LIMIT:', E10.2)
    ! ==----------------------------------------------------------------==
        CALL NECI_FRSBLK(M, N, NKRY, NBLOCK, V0, VS, A, V, AM, BM, T, W, WT, INDEX, &
                         SCR, LSCR, ISCR, LISCR, NDIAG, B2LIMIT, H, ICMAX, LAB, NROW, NCYCLE, &
                         PRINTOUT, TLargest, tDie2, tFail)
        if (tFail) then
            W(1) = sqrt(ReturnNan)
            return
        end if
    ! ==----------------------------------------------------------------==
    !..Uncomment to test to see if routine is working
    !..Exact eigenstates
    !      T1 = neci_etime()
    !      CALL DSYEV('V','U',M,H,M,WH,WORK2,3*M,INFO)
    !      WRITE(stdout,'(//14X,''Exact'',15X,''Lanczos'',10X,''Residual'')')
    !      DO I=1,N
    !        AUX=ABS(DDOT(M,V0(1,I),1,H(1,M-I+1),1))
    !        AUX=1.0_dp-AUX
    !        WRITE(stdout,'(6X,I3,2E19.11,2X,E10.3)') I,WH(M-I+1),W(I),AUX
    !      ENDDO
    !      T2 = neci_etime()
    !      T3=(T2-T1)
    !      WRITE(stdout,'(//5X,''TIME FOR EXACT DIAGONALISATION'',F10.2)')
    !     &       T3/1000.0_dp
    !      IF(PRINTOUT) THEN
    !          WRITE(stdout,'(//10X,''Neval'',15X,''Eigenvalue'')')
    !          DO I=1,N
    !              WRITE(stdout,'(10X,I3,15X,F19.11)') I,-1.0_dp*W(I)
    !          ENDDO
    !      ENDIF
    ! ==----------------------------------------------------------------==
    END
    ! ====================================================================
    SUBROUTINE NECI_FRSBLK(M, N, NKRY, NBLOCK, V0, VS, A, V, AM, BM, T, W, WT, &
                           INDEX, SCR, LSCR, ISCR, LISCR, NDIAG, B2LIMIT, &
                           H, ICMAX, LAB, NROW, NCYCLE, PRINTOUT, TLargest, tDie2, tFail)
        use constants, only: dp, sp
        use global_utilities
        use util_mod, only: neci_etime
        use error_handling_neci, only: stop_all
        IMPLICIT NONE
        LOGICAL :: PRINTOUT, TLargest, tDie2, tFail
        integer :: LISCR, ISCR(LISCR), I, J, NK, LL, NHPSI, N, NBLK, M, NBLEFF, NBL
        integer :: NKRY, NBLOCK, ICMAX, LSCR, NLEFT, NDIAG, NCURR, NCYCLE, ICYCLE
        real(dp) :: V0(M, N), VS(M, N), A(N, N), W(N)
        real(dp) :: V(M * NBLOCK * (NKRY + 1))
        real(dp) :: AM(NBLOCK * NBLOCK * (NKRY + 1)), BM(NBLOCK * NBLOCK * NKRY)
        real(dp) :: T(3 * NBLOCK * NKRY * NBLOCK * NKRY), WT(NBLOCK * NKRY)
        real(dp) :: SCR(LSCR), AUX, B2LIMIT, B2, DDOT, B2MIN, B2MAX
        integer :: NCONV
        integer :: INDEX(N), NROW(M), LAB(M, ICMAX), INFO
        type(timer), save :: proc_timer
        real(dp) t1, t2, t3, tarr(2)
    !..temp
        real(dp) ::  H(*)
        character(*), parameter :: this_routine = 'NECI_FRSBLK'

    !..temp
    ! ==-----------------------------------------------------------------==
        proc_timer%timer_name = '    FRSBLK'
        call set_timer(proc_timer)
        NK = NBLOCK * NKRY
    !..scratch space for divide-and-conquer banded matrix routine
    !      K=INT(LOG(DFLOAT(NK)))+1
    !      LSCR1=1+4*NK+2*NK*K+3*NK*NK
    !      LSCR2=2+5*NK
    !      LL=MAX(3*NK,LSCR1+LSCR2)
    !.. scratch space for banded matrix diagonaliser
        LL = 3 * NK
        IF (LSCR < LL) THEN
            WRITE(stdout, *) ' LL:', LL
    !        WRITE(stdout,*) 'LSCR1:',LSCR1
    !        WRITE(stdout,*) 'LSCR2:',LSCR2
            WRITE(stdout, *) 'LSCR:', LSCR
            call stop_all(this_routine, ' LSCR TOO SMALL ')
        END IF
        T1 = neci_etime(tarr)
    ! ====================================================================
        INDEX(1:N) = 0
        NHPSI = 0
        ICYCLE = 0
        NCONV = 0
        IF (PRINTOUT) THEN
            WRITE(stdout, 10000) ICYCLE, NCONV, B2MAX, B2MIN, NHPSI
        END IF
    10000 FORMAT(2X, 'ICYCLE:', I3, 1X, 'NCONV:', I3, 2X, 'B2MAX:', F10.5, 6X, 'B2MIN:', F10.5, 5X, 'NHPSI:', I3)
    !..VS=H.V0
    !..My matrix multiplication routine implimented 8/11/02 DCT
        CALL MY_HPSI(M, N, NROW, LAB, H, V0, VS, TLargest)
    !      CALL NECI_HSPI(M,N,H,V0,VS)
        NHPSI = NHPSI + N
    !..   Ovlap: V0^T VS.
        CALL NECI_OVLAP(M, N, A, V0, VS)
    !..   AY=YE
        CALL DSYEV('V', 'U', N, A, N, W, SCR, LSCR, INFO)
        CALL NECI_REORDER(N, N, W, A)
    !..Rotate: V0 -> V0.Y
        CALL NECI_ROTATE(M, N, V0, A, SCR, M * N)
    !..Rotate: VS -> VS.Y
        CALL NECI_ROTATE(M, N, VS, A, SCR, M * N)
    !=====================================================================
        NCONV = 0
        DO ICYCLE = 1, NCYCLE
    !        CALL NECI_WRITE_MATRIX('   W:   ',N,1,W)
    ! ====================================================================
    !..Residual: H(VY) - (VY)E and test for convergence
            B2MAX = 0.0_dp
            B2MIN = 1.D30
            NCURR = NCONV + 1
            DO J = NCURR, N
                CALL DCOPY(M, VS(1, J), 1, SCR, 1)
                CALL DAXPY(M, -W(J), V0(1, J), 1, SCR, 1)
                B2 = DDOT(M, SCR, 1, SCR, 1)
                IF (B2 < B2LIMIT) THEN
                    NCONV = NCONV + 1
                    INDEX(J) = NCONV
                END IF
                IF (B2 > B2MAX) B2MAX = B2
                IF (B2 < B2MIN) B2MIN = B2
            END DO
    !..
            IF (PRINTOUT) THEN
                WRITE(stdout, '(5X,I4,2X,I4,2X,2(E10.3,3X),F6.2)') ICYCLE, NCONV, B2MAX, B2MIN, real(NHPSI, dp) / real(N, dp)
            END IF
    !..Order states
            DO I = NCURR, N
                J = INDEX(I)
                IF (J /= 0 .AND. J /= I) THEN
                    CALL DSWAP(M, V0(1, I), 1, V0(1, J), 1)
                    CALL DSWAP(M, VS(1, I), 1, VS(1, J), 1)
                    INDEX(J) = J
                    INDEX(I) = 0
                    AUX = W(I)
                    W(I) = W(J)
                    W(J) = AUX
                END IF
            END DO
            NCURR = NCONV + 1
            NLEFT = N - NCONV
            IF (NCURR > NDIAG) GOTO 100
            NBL = MIN(NLEFT, NBLOCK)
            DO I = NCURR, N, NBL
                NBLEFF = MIN(NBL, N - I + 1)
    ! ==------------------------------------------------------------------==
    ! Preparation for refinement
    ! ==------------------------------------------------------------------==
                CALL NECI_PRPKRV(M, NBLEFF, NKRY, NCONV, V0(1, 1), V0(1, I), VS(1, I), &
                                 V, AM, BM, H, W(I), SCR, LSCR, NHPSI, LAB, NROW, TLargest, tDie2, tFail)
                if (tFail) return

    ! Refinement Loop
    ! ==------------------------------------------------------------------==
                NBLK = NBLEFF * NKRY
                CALL NECI_KRYREF(M, NBLEFF, NKRY, NCONV, V0, V, AM, BM, T, NBLK * NBLK, &
                                 WT, NBLK, H, SCR, LSCR, ISCR, LISCR, NHPSI, LAB, NROW, TLargest, tDie2, tFail)
                if (tFail) return
    ! ==------------------------------------------------------------------==
    !..V=[V_1 V_2.. V_L] Y'
                CALL DGEMM('N', 'N', M, NBLEFF, NBLK, 1.0_dp, V, M, T, NBLK, 0.0_dp, V0(1, I), M)
    !==-------------------------------------------------------------------==
    !==  End of refinement over states                                    ==
    !==-------------------------------------------------------------------==
            END DO
    !..GS orthogonalisation on refined states
            CALL NECI_MY_GSORTHO(M, V0, NCURR - 1, V0(1, NCURR), NLEFT, A, tDie2, tFail)
            if (tFail) return
    !..HPSI: H V0. Enter only refined states
            CALL MY_HPSI(M, NLEFT, NROW, LAB, H, V0(1, NCURR), VS(1, NCURR), TLargest)
    !        CALL NECI_HSPI(M,NLEFT,H,V0(1,NCURR),VS(1,NCURR))
            NHPSI = NHPSI + NLEFT
    !..Ovlap: V0^T VS
            CALL NECI_OVLAP(M, N, A, V0, VS)
    !..AY=YE
            CALL DSYEV('V', 'U', N, A, N, W, SCR, LSCR, INFO)
            CALL NECI_REORDER(N, N, W, A)
    !..Rotate: V -> VY
            CALL NECI_ROTATE(M, N, V0, A, SCR, M * N)
    !..Rotate: HV -> HVY
            CALL NECI_ROTATE(M, N, VS, A, SCR, M * N)
    ! ==============================================================
        END DO                     ! Loop over icycle
    !.. End of Lanczos diagonalisation
    100 CONTINUE
        IF (PRINTOUT) THEN
            WRITE(stdout, '(//''    NCONV:'',I5)') NCONV
        END IF
        T2 = neci_etime(tarr)
        T3 = (T2 - T1)
        IF (PRINTOUT) THEN
            WRITE(stdout, '(//5X,''TIME FOR LANCZOS DIAGONALISATION'',F10.2)') T3 / 1000.0_dp
        END IF
        call halt_timer(proc_timer)
    !     ================================================================
        RETURN
    END
    ! ======================================================================
    SUBROUTINE NECI_MGS(M, N, A, LDA, R, LDR, tDie2, tFail)
    !     ==--------------------------------------------------------------==
        use global_utilities
        use constants, only: dp, sp
    !      IMPLICIT real(dp) (A-H,O-Z)
        IMPLICIT NONE
        integer :: LDA, LDR, N, M
        integer(sp) :: INFO
        real(dp) :: A(LDA, *), R(LDR, *)
        LOGICAL :: tDie2, tfail
        type(timer), save :: proc_timer
    !     ==================================================================
        proc_timer%timer_name = 'NECI_MGS'
        call set_timer(proc_timer)
        CALL NECI_RGS(M, N, A, M, R, tDie2, tFail)
        if (tFail) return
        CALL DTRTRI('U', 'N', N, R, N, INFO)
        call halt_timer(proc_timer)
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_SETUP_MATRIX(M, N, A, TSYM)
        use constants, only: dp
        use dSFMT_interface, only: genrand_real2_dSFMT
        IMPLICIT NONE
    !..M = NDET
    !..N = NEVAL
        integer :: M, N, I, J
        real(dp) :: A(M, N)
        LOGICAL TSYM
        integer iseed
        IF (TSYM) THEN
        DO I = 1, N
        DO J = I, M
            A(J, I) = genrand_real2_dSFMT()
            A(I, J) = A(J, I)
        END DO
        END DO
        ELSE
        DO I = 1, N
        DO J = 1, M
            A(J, I) = genrand_real2_dSFMT()
        END DO
        END DO
        END IF
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_WRITE_MATRIX(CHAR, M, N, A)
    !      IMPLICIT real(dp) (A-H,O-Z)
        use constants, only: dp
        use error_handling_neci, only: stop_all
        IMPLICIT NONE
        CHARACTER(*) CHAR
        INTEGER :: N, M, I, J
        real(dp) :: A(M, N)
        WRITE(stdout, *) CHAR
        DO I = 1, M
            WRITE(stdout, 1000)(A(I, J), J=1, N)
        END DO
    1000 FORMAT(12E15.6)
        RETURN
    END
    ! ======================================================================
    SUBROUTINE NECI_PUTTMAT(CHAR, T, M, N, A, MA, NA, IBEG, JBEG)
        use constants, only: dp
        use error_handling_neci, only: stop_all
        IMPLICIT NONE
        integer :: N, M, MA, NA, IBEG, JBEG, I, J
        real(dp) :: T(M, N), A(MA, NA)
        CHARACTER(1) CHAR
        character(*), parameter :: this_routine = 'NECI_PUTTMAT'
        IF (MA * IBEG > M) call stop_all(this_routine, ' MA+IBEG.GT.M')
        IF (NA * JBEG > N) call stop_all(this_routine, ' MA+IBEG.GT.M')
        IF (CHAR == 'N') THEN
    !        DO J=1,NA
    !          DO I=1,MA
    !            T(I+MA*(IBEG-1),J+NA*(JBEG-1))=A(I,J)
    !          ENDDO
    !        ENDDO
            T((ibeg - 1) * MA + 1:ibeg * MA,(jbeg - 1) * MA + 1:jbeg * MA) = A
        ELSEIF (CHAR == 'T') THEN
    !        DO J=1,NA
    !          DO I=1,MA
    !            T(I+MA*(IBEG-1),J+NA*(JBEG-1))=A(J,I)
    !          ENDDO
    !        ENDDO
            T((ibeg - 1) * MA + 1:ibeg * MA,(jbeg - 1) * MA + 1:jbeg * MA) = transpose(A)
        ELSE
            call stop_all(this_routine, 'ILLEGAL CHAR')
        END IF
        RETURN
    END
    ! ======================================================================
    SUBROUTINE NECI_REORDER(M, N, W, A)
        use constants, only: dp
    !      IMPLICIT real(dp) (A-H,O-Z)
        IMPLICIT NONE
        integer :: N, M, J
        real(dp) :: W(N), A(M, N)
        real(dp) :: AUX
        DO J = 1, N / 2
            CALL DSWAP(M, A(1, J), 1, A(1, N - J + 1), 1)
            AUX = W(N - J + 1)
            W(N - J + 1) = W(J)
            W(J) = AUX
        END DO
        RETURN
    END
    ! ======================================================================
    SUBROUTINE NECI_KRYREF( &
        M, N, NKRY, NCONV, V0, V, AM, BM, T, LDT, WT, &
        LDWT, H, SCR, LSCR, ISCR, LISCR, NHPSI, LAB, NROW, TLargest, tDie2, tFail)
    !     ==================================================================
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp) (A-H,O-Z)
        IMPLICIT NONE
        integer :: N, M, LAB(M, *), NROW(M), LSCR, LDWT, LDT, NKRY
        integer :: LISCR, J, NCONV, NHPSI, I, ISCR(LISCR)
        real(dp) :: V0(M, *), V(M, N, NKRY + 1), AM(N, N, NKRY + 1), BM(N, N, NKRY)
        real(dp) :: WT(LDWT), H(*), SCR(LSCR)
        real(dp) :: T(LDT)
        LOGICAL :: TLargest, tDie2, tFail
        type(timer), save :: proc_timer
    !     ==================================================================
        proc_timer%timer_name = '    KRYREF               '
        call set_timer(proc_timer)
        DO I = 3, NKRY
    !..V_I = V_I - V_(I-1) A_(I-1)
            CALL NECI_RSDBLK('N', M, N, V(1, 1, I - 1), AM(1, 1, I - 1), V(1, 1, I))
    !..V_I B_I = V_I
            CALL NECI_MGS(M, N, V(1, 1, I), M, BM(1, 1, I), N, tDie2, tFail)
            if (tFail) return
    !..HPSI:  H V_I
            CALL MY_HPSI(M, N, NROW, LAB, H, V(1, 1, I), V(1, 1, I + 1), TLargest)
    !        CALL NECI_HSPI(M,N,H,V(1,1,I),V(1,1,I+1))
            NHPSI = NHPSI + N
    !.. project out converged states
            CALL NECI_PRJCNV(M, N, NCONV, V0, V(1, 1, I + 1), SCR)
    !..V_(I+1) = H V_I - V_(I-1)B_I
            CALL NECI_RSDBLK('T', M, N, V(1, 1, I - 1), BM(1, 1, I), V(1, 1, I + 1))
    !..Ovlap:  A_I=V_I^T V_(I+1)
            CALL NECI_OVLAP(M, N, AM(1, 1, I), V(1, 1, I), V(1, 1, I + 1))
    !..Setup T-matrix. First diagonal terms
    !        CALL AZZERO(T,N*NKRY*N*NKRY)
    !        CALL NECI_PUTTMAT('N',T,I*N,I*N,AM(1,1,1),N,N,1,1)
    !        DO J=2,I
    !          CALL NECI_PUTTMAT('N',T,I*N,I*N,AM(1,1,J),N,N,J,J)
    !          CALL NECI_PUTTMAT('T',T,I*N,I*N,BM(1,1,J),N,N,J-1,J)
    !          CALL NECI_PUTTMAT('N',T,I*N,I*N,BM(1,1,J),N,N,J,J-1)
    !        ENDDO
    !..TY=YE
    !        CALL DSYEV('V','U',I*N,T,I*N,WT,SCR,LSCR,INFO)
    !        CALL NECI_REORDER(I*N,WT,T)
    !.. Compute error V_3 Y'', where Y'' is the last N rows of Y
    !        CALL GETMAT('N',T,I*N,I*N,AM(1,1,I+1),N,N,I,1)
    !        CALL DGEMM('N','N',M,N,N,1.0_dp,V(1,1,I+1),M,AM(1,1,I+1),
    !     &       N,0.0_dp,SCR,M)
    !        DO J=1,N
    !          AUX=DNRM2(M,SCR(M*(J-1)+1),1)**2
    !          IF(AUX.LT.B2MIN) B2MIN=AUX
    !          IF(AUX.GT.B2MAX) B2MAX=AUX
    !        ENDDO
        END DO                     !  End loop over ikry
    ! ==================================================================
    !..Setup T-matrix. First diagonal terms

        I = NKRY
        T(1:(I * N)**2) = 0
        CALL NECI_PUTTMAT('N', T, I * N, I * N, AM(1, 1, 1), N, N, 1, 1)
        DO J = 2, I
            CALL NECI_PUTTMAT('N', T, I * N, I * N, AM(1, 1, J), N, N, J, J)
            CALL NECI_PUTTMAT('T', T, I * N, I * N, BM(1, 1, J), N, N, J - 1, J)
            CALL NECI_PUTTMAT('N', T, I * N, I * N, BM(1, 1, J), N, N, J, J - 1)
        END DO
        IF (N > 0) THEN
    !..Full matrix diag.
            CALL NECI_JACOBI(I * N, N, T, WT, SCR, LSCR, ISCR, 5 * I * N, ISCR(5 * I * N + 1))
        ELSE
    !..banded matrix diagonalisation. This seems to be slower than
    !..full diag.
            CALL NECI_BANDM(I * N, N, T, WT, SCR, LSCR, ISCR, 5 * I * N, ISCR(5 * I * N + 1))
        END IF
    !..
    ! ===================================================================
        call halt_timer(proc_timer)
        RETURN
    END
    ! ==================================================================
    SUBROUTINE NECI_PRPKRV(M, N, NKRY, NCONV, VCONV, V0, VS, V, AM, &
                           BM, H, W, SCR, LSCR, NHPSI, LAB, NROW, TLargest, tDie2, tFail)
    ! ==------------------------------------------------------------------==
    ! == Returns A_1,A_2,B_2 and V in a form suitable for krylov          ==
    ! == refinement                                                       ==
    ! ==------------------------------------------------------------------==
        use constants, only: dp
        use global_utilities
    !      IMPLICIT real(dp) (A-H,O-Z)
        IMPLICIT NONE
        integer :: M, N, LAB(M, *), NROW(M), LSCR, I, NHPSI, NKRY, NCONV
        real(dp) :: VCONV(M, *), V0(M, N), VS(M, N), V(M, N, NKRY + 1)
        real(dp) :: AM(N, N, NKRY + 1), BM(N, N, NKRY), W(N), H(*), SCR(LSCR)
        LOGICAL :: TLargest, tDie2, tFail
        type(timer), save :: proc_timer
    ! ==------------------------------------------------------------------==
        proc_timer%timer_name = 'NECI_PRPKRV'
        call set_timer(proc_timer)
    ! ==------------------------------------------------------------------==
        CALL DCOPY(M * N, V0, 1, V(1, 1, 1), 1)
        CALL DCOPY(M * N, VS, 1, V(1, 1, 2), 1)
    !..Compute residual
        DO I = 1, N
            CALL DAXPY(M, -W(I), V(1, I, 1), 1, V(1, I, 2), 1)
        END DO
    !..V2 B2=R
        BM = 0.0_dp
        CALL NECI_MGS(M, N, V(1, 1, 2), M, BM(1, 1, 2), N, tDie2, tFail)
        if (tFail) return
    !..   Setup A_1=diag(w)
        AM(:, :, 1) = 0.0_dp
        DO I = 1, N
            AM(I, I, 1) = W(I)
        END DO
    !..V_3=H.V2
        CALL MY_HPSI(M, N, NROW, LAB, H, V(1, 1, 2), V(1, 1, 3), TLargest)
    !      CALL NECI_HSPI(M,N,H,V(1,1,2),V(1,1,3))
        NHPSI = NHPSI + N
    !.. project out converged states
        CALL NECI_PRJCNV(M, N, NCONV, VCONV, V(1, 1, 3), SCR)
    !..V_3=V_3-V_1 B_2^T
        CALL NECI_RSDBLK('T', M, N, V(1, 1, 1), BM(1, 1, 2), V(1, 1, 3))
    !..Ovlap
        CALL NECI_OVLAP(M, N, AM(1, 1, 2), V(1, 1, 2), V(1, 1, 3))
    ! ==------------------------------------------------------------------==
        call halt_timer(proc_timer)
        RETURN
    END
    !     ==================================================================
    SUBROUTINE NECI_MY_GSORTHO(M, C0, N0, CP, NP, SMAT, tDie2, tFail)
    !     ==--------------------------------------------------------------==
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp) (A-H,O-Z)
        IMPLICIT NONE
    !     Arguments
        LOGICAL :: tDie2, tFail
        INTEGER N0, NP, M
        real(dp) C0(M, *), CP(M, *)
        real(dp) SMAT(NP, *)             !MAX(N0,NP)
    !     Variables
        type(timer), save :: proc_timer
    !     ==--------------------------------------------------------------==
        proc_timer%timer_name = 'MY_GSORTHO'
        call set_timer(proc_timer)
        IF (NP < 1) RETURN
        IF (N0 > 0) THEN
    !..SMAT=CP.C0
            CALL DGEMM('T', 'N', NP, N0, M, 1.0_dp, CP, M, C0, M, 0.0_dp, SMAT, NP)
    !..CP -> CP-C0*SMAT
            CALL DGEMM('N', 'T', M, NP, N0, -1.0_dp, C0, M, SMAT, NP, 1.0_dp, CP, M)
        END IF
        CALL NECI_MGS(M, NP, CP, M, SMAT, NP, tDie2, tFail)
        call halt_timer(proc_timer)
    !     ==--------------------------------------------------------------==
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_PUTTAB(N, KD, A, SCR, LSCR)
        use error_handling_neci, only: stop_all
        use constants, only: dp
        IMPLICIT NONE
        integer :: N, LSCR, KD, I, IBEG, J
        real(dp) :: A(N, N), SCR(LSCR)
        character(*), parameter :: this_routine = 'NECI_PUTTAB'
    !..
        IF (LSCR < N) THEN
            WRITE(stdout, *) ' LSCR:', LSCR
            WRITE(stdout, *) 'N:', N
            call stop_all(this_routine, 'LSCR LT N ')
        END IF
        SCR(1:N) = 0.0_dp
        DO J = 1, N
            CALL DCOPY(N, A(1, J), 1, SCR, 1)
            A(1:N, J) = 0.0_dp
            IBEG = MAX(J - KD, 1)
            DO I = IBEG, J
                A(1 + KD + I - J, J) = SCR(I)
            END DO
        END DO
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_HPSI(M, N, H, V0, VS)
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp)(A-H,O-Z)
        IMPLICIT NONE
        integer :: N, M
        real(dp) :: H(M, M), V0(M, N), VS(M, N)
        type(timer), save :: proc_timer
        proc_timer%timer_name = '      NECI_HSPI'
        call set_timer(proc_timer)
        CALL DGEMM('N', 'N', M, N, M, 1.0_dp, H, M, V0, M, 0.0_dp, VS, M)
        call halt_timer(proc_timer)
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_OVLAP(M, N, A, V0, VS)
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp)(A-H,O-Z)
        IMPLICIT NONE
        integer :: N, M
        real(dp) :: A(N, N), V0(M, N), VS(M, N)
        type(timer), save :: proc_timer
        proc_timer%timer_name = 'NECI_OVLAP'
        call set_timer(proc_timer)
        CALL DGEMM('T', 'N', N, N, M, 1.0_dp, V0, M, VS, M, 0.0_dp, A, N)
        call halt_timer(proc_timer)
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_RSDBLK(CHAR, M, N, V1, A, V2)
        use global_utilities
        use constants, only: dp
        use util_mod, only: stop_all
        IMPLICIT NONE
        CHARACTER(1) CHAR
        integer :: N, M
        real(dp) :: A(N, N), V1(M, N), V2(M, N)
        type(timer), save :: proc_timer
        character(*), parameter :: this_routine = 'NECI_RSDBLK'
    !..V_I = V_I - V_(I-1) A_(I-1)
        proc_timer%timer_name = '    NECI_RSDBLK'
        call set_timer(proc_timer)
        IF (CHAR == 'N') THEN
            CALL DGEMM('N', 'N', M, N, N, -1.0_dp, V1, M, A, N, 1.0_dp, V2, M)
        ELSEIF (CHAR == 'T') THEN
            CALL DGEMM('N', 'T', M, N, N, -1.0_dp, V1, M, A, N, 1.0_dp, V2, M)
        ELSE
            call stop_all(this_routine, ' CHAR ILLEGAL IN NECI_RSDBLK ')
        END IF
        call halt_timer(proc_timer)
        RETURN
    END
    ! =======================================================================
    SUBROUTINE NECI_ROTATE(M, N, V0, A, SCR, LSCR)
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp)(A-H,O-Z)
        IMPLICIT NONE
        integer :: N, M, LSCR
        REAL(dp) :: A(N, N), V0(M, N), SCR(LSCR)
        type(timer), save :: proc_timer
        proc_timer%timer_name = '    NECI_ROTATE'
        call set_timer(proc_timer)
        CALL DGEMM('N', 'N', M, N, N, 1.0_dp, V0, M, A, N, 0.0_dp, SCR, M)
        CALL DCOPY(M * N, SCR, 1, V0, 1)
        call halt_timer(proc_timer)
        RETURN
    END

    SUBROUTINE NECI_JACOBI(M, N, T, WT, SCR, LSCR, ISCR, LISCR, IFAIL)
    ! ==------------------------------------------------------------------==
    ! Returns the n largest eigenvalues of MxM symtric matrix T           ==
    ! ==------------------------------------------------------------------==
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp) (A-H,O-Z)
        IMPLICIT NONE
        INTEGER :: IFAIL(*), INFO, M, N, LISCR, LSCR, MEVAL, ISCR(LISCR)
        REAL(dp) :: T(M, M, *), WT(*), SCR(LSCR)
        type(timer), save :: proc_timer
        proc_timer%timer_name = '    JACOBI'
        call set_timer(proc_timer)
        IF (M == N) THEN
            CALL DSYEV('V', 'U', N, T, N, WT, SCR, LSCR, INFO)
        ELSE
            CALL DSYEVX('V', 'I', 'U', M, T, M, 1.0_dp, 1.0_dp, M - N + 1, M, 0.0_dp, MEVAL, WT, T(1, 1, 2), M, SCR, LSCR, ISCR, IFAIL, INFO)
            IF (MEVAL < N) THEN
                WRITE(stdout, *) ' WARNING| DSYEVX RETURNED MEVAL < N', MEVAL, N
            END IF
            CALL DCOPY(M * N, T(1, 1, 2), 1, T(1, 1, 1), 1)
        END IF
        CALL NECI_REORDER(M, N, WT, T)
        call halt_timer(proc_timer)
        RETURN
    END
    !     ==================================================================
    SUBROUTINE NECI_BANDM(IN, N, T, WT, SCR, LSCR, ISCR, LISCR, IFAIL)
    ! ==------------------------------------------------------------------==
    ! Returns the n largest eigenvalues of MxM symtric banded matrix T    ==
    ! ==------------------------------------------------------------------==
        use global_utilities
        use constants, only: dp
        use error_handling_neci, only: stop_all
        IMPLICIT NONE
        type(timer), save :: proc_timer
        INTEGER :: IFAIL(*), IN, LSCR, LISCR, N, MEVAL, INFO, ISCR(LISCR)
        REAL(dp) :: T(IN, IN, *), WT(IN), SCR(LSCR)
        character(*), parameter :: this_routine = 'NECI_BANDM'
    !     ==================================================================
        proc_timer%timer_name = '     BANDM'
        call set_timer(proc_timer)
    !..Put T is form suitable for banded matrix diagonalisation
        CALL NECI_PUTTAB(IN, N, T, SCR, LSCR)
        IF (LSCR < 7 * IN) THEN
            WRITE(stdout, *) ' 7*IN:', 7 * IN
            WRITE(stdout, *) 'LSCR:', LSCR
            call stop_all(this_routine, ' LSCR TOO SMALL IN BANDM')
        END IF
        IF (LISCR < 5 * IN) THEN
            WRITE(stdout, *) ' 5*IN:', 5 * IN
            WRITE(stdout, *) 'LISCR:', LISCR
            call stop_all(this_routine, ' LISCR TOO SMALL IN BANDM')
        END IF
        CALL DSBEVX('V', 'I', 'U', IN, N, T, IN, T(1, 1, 2), IN, 1.0_dp, 1.0_dp, &
                    IN - N + 1, IN, 1.0e-10_dp, MEVAL, WT, T(1, 1, 3), IN, SCR, ISCR, IFAIL, INFO)
        IF (MEVAL < N) THEN
            WRITE(stdout, *) ' WARNING| DSBEVX RETURNED MEVAL < N', MEVAL, N
        END IF
        CALL DCOPY(IN * IN, T(1, 1, 3), 1, T(1, 1, 1), 1)
        CALL NECI_REORDER(IN, N, WT, T)
        call halt_timer(proc_timer)
        RETURN
    END
    ! =====================================================================
    SUBROUTINE NECI_PRJCNV(M, NP, N0, V0, V, SMAT)
        use global_utilities
        use constants, only: dp
    !      IMPLICIT real(dp)(A-H,O-Z)
        IMPLICIT NONE
        INTEGER :: M, N, N0, NP
        REAL(dp) :: V(M, NP), V0(M, *), SMAT(N0, *)
        type(timer), save :: proc_timer
        IF (N0 == 0) RETURN
        proc_timer%timer_name = '    PRJCNV'
        call set_timer(proc_timer)
    !..SMAT=V0^T.V
        CALL DGEMM('T', 'N', N0, NP, M, 1.0_dp, V0, M, V, M, 0.0_dp, SMAT, N0)
    !deb      DO J=1,N0
    !deb        CALL DSCAL(NP,W(J),SMAT(J,1),N0)
    !deb      ENDDO
    !..   V -> V-V0*SMAT
        CALL DGEMM('N', 'N', M, NP, N0, -1.0_dp, V0, M, SMAT, N0, 1.0_dp, V, M)
        call halt_timer(proc_timer)
        RETURN
    END
    !     =================================================================
    SUBROUTINE NECI_RGS(M, N, CP, LDCP, SMAT, tDie2, tFail)
    !     ==--------------------------------------------------------------==
    !     ==  GRAM-SCHMIDT ORTHOGONALIZATION                              ==
    !     ==--------------------------------------------------------------==
        use global_utilities
        use constants, only: dp
        use error_handling_neci, only: warning_neci
        IMPLICIT NONE
    !     Arguments
        INTEGER N, M, LDCP
        real(dp) SMAT(N, N)
        real(dp) CP(LDCP, N)
        logical :: tDie2, tFail
    !     VARIABLES
        type(timer), save :: proc_timer
    !     ==--------------------------------------------------------------==
        IF (N <= 0) RETURN
        proc_timer%timer_name = 'NECI_RGS'
        call set_timer(proc_timer)
        SMAT = 0
        CALL DSYRK('U', 'T', N, M, 1.0_dp, CP, M, 0.0_dp, SMAT, N)
        CALL NECI_UINV('U', SMAT, N, N, tDie2, tFail)
        if (tFail) return
        CALL DTRMM('R', 'U', 'N', 'N', M, N, 1.0_dp, SMAT, N, CP, M)
        call halt_timer(proc_timer)
    !     ==--------------------------------------------------------------==
        RETURN
    END
    !     ==================================================================
    SUBROUTINE NECI_UINV(UPLO, SMAT, LDA, N, tDie2, tFail)
    !     ==--------------------------------------------------------------==
        use constants, only: dp, sp
        use error_handling_neci, only: stop_all, warning_neci
        IMPLICIT NONE
    !     ARGUMENTS
        CHARACTER(1) UPLO
        INTEGER LDA, N
        real(dp) SMAT(LDA, N)
    !     VARIABLES
        INTEGER(sp) INFO
        LOGICAL :: tDie2, tFail
    !     ==--------------------------------------------------------------==
        CALL DPOTRF(UPLO, N, SMAT, LDA, INFO)
        IF (INFO /= 0) THEN
            WRITE(stdout, *) "INFO is : ", INFO
        END IF
        IF (INFO /= 0) THEN
        if (tDie2) then
            CALL Stop_All('UINV', 'ILLEGAL RESULTS DGETRF')
        else
            CALL Warning_neci('UINV', 'ILLEGAL RESULTS DGETRF')
            tFail = .true.
            return
        end if
        end if
        CALL DTRTRI(UPLO, 'N', N, SMAT, LDA, INFO)
        IF (INFO /= 0) then
        if (tDie2) then
            CALL Stop_All('UINV', 'ILLEGAL RESULTS DTRTRI')
        else
            CALL Warning_neci('UINV', 'ILLEGAL RESULTS DTRTRI')
            tFail = .true.
        end if
        end if
    !     ==--------------------------------------------------------------==
        RETURN
    END
    !     ==================================================================

end module
