      SUBROUTINE NECI_FRSBLKH(M,ICMAX,N,H,LAB,V0,VS,NKRY,NKRY1,
     &    NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT,SCR,ISCR,INDEX,
     &    NCYCLE,B2LIMIT,PRINTOUT,TLargest,TSeeded,tDie)
!      IMPLICIT real(dp) (A-H,O-Z) 
      use constants, only: dp
      IMPLICIT NONE
C ==----------------------------------------------------------------==
C ==    M is NDET : the number of rows in Hamil                     ==
C ==    ICMAX : the MAX. number of columns in Hamil                 ==
C ==    N is NEVAL : the number of columns in the wavevector        ==
C ==    PRINTOUT is a logical testing whether to print results      ==
C ==    TSeeded if true, means that a seed for RAN2 has been set,   ==
C ==    and shouldn't be reset.
C ==    tDie indicates whether STOP should be called if it fails
C ==    to diagonalise
C ==----------------------------------------------------------------==
      
      LOGICAL :: PRINTOUT,TLargest,TSeeded,tDie2,tFail
      LOGICAL :: tDie
      integer :: LISCR,ISCR(LISCR),LSCR,NBLOCK,NKRY1,NKRY,N,ICMAX,M
      REAL(dp) :: V0(M,N),VS(M,N),A(N,N),W(N) 
      REAL(dp) :: V(M*NBLOCK*NKRY1)
      REAL(dp) :: AM(NBLOCK*NBLOCK*NKRY1),BM(NBLOCK*NBLOCK*NKRY)
      REAL(dp) :: T(3*NBLOCK*NKRY*NBLOCK*NKRY),WT(NBLOCK*NKRY)
      real(dp) :: SCR(LSCR),B2LIMIT
      integer :: NCYCLE,NDIAG
      integer :: INDEX(N)
C..temp
      real(dp) :: H(M,ICMAX),ReturnNan
C..temp
      integer :: LAB(M,ICMAX),NROW(M)
C ==----------------------------------------------------------------==
      tDie2 = tDie
      ReturnNan=-1.0_dp
      tFail=.false.

C..Initialise V0 
      CALL NECI_SETUP_MATRIX(M,N,V0,.FALSE.,TSeeded)
      CALL NECI_MGS(M,N,V0,M,A,N,tDie2,tFail)
      if(tFail) then
          W(1) = sqrt(ReturnNan)
          return
      endif

C..Calculate top 80% of N 
      NDIAG=N-INT(N*0.2) 
C..
      IF(PRINTOUT) THEN
          WRITE(6,20000) M,N,NKRY,NBLOCK,NDIAG,LSCR,B2LIMIT
      ENDIF
20000 FORMAT(2X,'M:',I7/2X,'N:',I7/2X,'NKRY:',I7/2X,'NBLOCK:',I7/2X,
     &    'NDIAG:',I7/2X,'LSCR:',I7/2X,'B2LIMIT:',E10.2)
C ==----------------------------------------------------------------==
      CALL NECI_FRSBLK(M,N,NKRY,NBLOCK,V0,VS,A,V,AM,BM,T,W,WT,INDEX,
     &     SCR,LSCR,ISCR,LISCR,NDIAG,B2LIMIT,H,ICMAX,LAB,NROW,NCYCLE,
     &     PRINTOUT,TLargest,tDie2,tFail)
      if(tFail) then
          W(1) = sqrt(ReturnNan)
          return
      endif
C ==----------------------------------------------------------------==
C..Uncomment to test to see if routine is working
C..Exact eigenstates
C      T1 = neci_etime()
C      CALL DSYEV('V','U',M,H,M,WH,WORK2,3*M,INFO)
C      WRITE(6,'(//14X,''Exact'',15X,''Lanczos'',10X,''Residual'')')
C      DO I=1,N
C        AUX=ABS(DDOT(M,V0(1,I),1,H(1,M-I+1),1))
C        AUX=1.0_dp-AUX
C        WRITE(6,'(6X,I3,2E19.11,2X,E10.3)') I,WH(M-I+1),W(I),AUX
C      ENDDO
C      T2 = neci_etime()
C      T3=(T2-T1)
C      WRITE(6,'(//5X,''TIME FOR EXACT DIAGONALISATION'',F10.2)')
C     &       T3/1000.0_dp
C      IF(PRINTOUT) THEN
C          WRITE(6,'(//10X,''Neval'',15X,''Eigenvalue'')')
C          DO I=1,N
C              WRITE(6,'(10X,I3,15X,F19.11)') I,-1.0_dp*W(I)
C          ENDDO
C      ENDIF
C ==----------------------------------------------------------------==
      END
C ====================================================================
      SUBROUTINE NECI_FRSBLK(M,N,NKRY,NBLOCK,V0,VS,A,V,AM,BM,T,W,WT,
     &     INDEX,SCR,LSCR,ISCR,LISCR,NDIAG,B2LIMIT,H,ICMAX,LAB,NROW,
     &     NCYCLE,PRINTOUT,TLargest,tDie2,tFail)
      use constants, only: dp,sp
      use global_utilities
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      LOGICAL :: PRINTOUT,TLargest,tDie2,tFail
      integer :: LISCR,ISCR(LISCR),I,J,NK,LL,NHPSI,N,NBLK,M,NBLEFF,NBL
      integer :: NKRY,NBLOCK,ICMAX,LSCR,NLEFT,NDIAG,NCURR,NCYCLE,ICYCLE
      real(dp) :: V0(M,N),VS(M,N),A(N,N),W(N) 
      real(dp) :: V(M*NBLOCK*(NKRY+1))
      real(dp) :: AM(NBLOCK*NBLOCK*(NKRY+1)),BM(NBLOCK*NBLOCK*NKRY)
      real(dp) :: T(3*NBLOCK*NKRY*NBLOCK*NKRY),WT(NBLOCK*NKRY)
      real(dp) :: SCR(LSCR),AUX,B2LIMIT,B2,DDOT,B2MIN,B2MAX
      integer :: NCONV
      integer :: INDEX(N),NROW(M),LAB(M,ICMAX),INFO 
      type(timer), save :: proc_timer
      real(sp) t1,t2,t3,tarr(2),neci_etime
C..temp
      real(dp) ::  H(*)
C..temp
C ==-----------------------------------------------------------------==
      proc_timer%timer_name='    FRSBLK'
      call set_timer(proc_timer)
      NK=NBLOCK*NKRY
C..scratch space for divide-and-conquer banded matrix routine
C      K=INT(LOG(DFLOAT(NK)))+1
C      LSCR1=1+4*NK+2*NK*K+3*NK*NK
C      LSCR2=2+5*NK
C      LL=MAX(3*NK,LSCR1+LSCR2)
C.. scratch space for banded matrix diagonaliser
      LL=3*NK      
      IF(LSCR.LT.LL) THEN 
        WRITE(6,*) ' LL:',LL
C        WRITE(6,*) 'LSCR1:',LSCR1
C        WRITE(6,*) 'LSCR2:',LSCR2
        WRITE(6,*) 'LSCR:',LSCR
        STOP ' LSCR TOO SMALL '  
      ENDIF
      T1 = neci_etime(tarr) 
C ====================================================================
      INDEX(1:N)=0
      NHPSI=0
      IF(PRINTOUT) THEN
          WRITE(6,10000) ICYCLE,NCONV,B2MAX,B2MIN,NHPSI
      ENDIF
10000 FORMAT(2X,'ICYCLE:',I3,1X,'NCONV:',I3,2X,'B2MAX:',F10.5,
     &  6X,'B2MIN:',F10.5,5X,'NHPSI:',I3)
C..VS=H.V0 
C..My matrix multiplication routine implimented 8/11/02 DCT
      CALL MY_HPSI(M,N,NROW,LAB,H,V0,VS,TLargest)
C      CALL NECI_HSPI(M,N,H,V0,VS)
      NHPSI=NHPSI+N 
C..   Ovlap: V0^T VS. 
      CALL NECI_OVLAP(M,N,A,V0,VS)
C..   AY=YE
      CALL DSYEV('V','U',N,A,N,W,SCR,LSCR,INFO)
      CALL NECI_REORDER(N,N,W,A)
C..Rotate: V0 -> V0.Y 
      CALL NECI_ROTATE(M,N,V0,A,SCR,M*N)
C..Rotate: VS -> VS.Y 
      CALL NECI_ROTATE(M,N,VS,A,SCR,M*N)
C=====================================================================
      NCONV=0 
      DO ICYCLE=1,NCYCLE
C        CALL NECI_WRITE_MATRIX('   W:   ',N,1,W)
C ====================================================================
C..Residual: H(VY) - (VY)E and test for convergence
        B2MAX=0.0_dp
        B2MIN=1.D30
        NCURR=NCONV+1 
        DO J=NCURR,N
          CALL DCOPY(M,VS(1,J),1,SCR,1)
          CALL DAXPY(M,-W(J),V0(1,J),1,SCR,1)
          B2=DDOT(M,SCR,1,SCR,1)
          IF(B2.LT.B2LIMIT) THEN 
            NCONV=NCONV+1
            INDEX(J)=NCONV 
          ENDIF
          IF(B2.GT.B2MAX) B2MAX=B2 
          IF(B2.LT.B2MIN) B2MIN=B2  
        ENDDO
C..
        IF(PRINTOUT) THEN
            WRITE(6,'(5X,I4,2X,I4,2X,2(E8.3,3X),F6.2)') 
     &       ICYCLE,NCONV,B2MAX,B2MIN,real(NHPSI,dp)/real(N,dp)
        ENDIF
C..Order states         
        DO I=NCURR,N
          J=INDEX(I)
          IF(J.NE.0.AND.J.NE.I) THEN 
            CALL DSWAP(M,V0(1,I),1,V0(1,J),1)
            CALL DSWAP(M,VS(1,I),1,VS(1,J),1)
            INDEX(J)=J
            INDEX(I)=0 
            AUX=W(I)
            W(I)=W(J)
            W(J)=AUX
          ENDIF
        ENDDO
        NCURR=NCONV+1 
        NLEFT=N-NCONV
        IF(NCURR.GT.NDIAG) GOTO 100 
        NBL=MIN(NLEFT,NBLOCK) 
        DO I=NCURR,N,NBL
          NBLEFF=MIN(NBL,N-I+1)
C ==------------------------------------------------------------------==
C Preparation for refinement
C ==------------------------------------------------------------------==
          CALL NECI_PRPKRV(M,NBLEFF,NKRY,NCONV,V0(1,1),V0(1,I),VS(1,I),
     &         V,AM,BM,H,W(I),SCR,LSCR,NHPSI,LAB,NROW,TLargest,tDie2,
     &          tFail)
          if(tFail) return

C Refinement Loop
C ==------------------------------------------------------------------==
          NBLK=NBLEFF*NKRY
          CALL NECI_KRYREF(M,NBLEFF,NKRY,NCONV,V0,V,AM,BM,T,NBLK*NBLK,
     &      WT,NBLK,H,SCR,LSCR,ISCR,LISCR,NHPSI,LAB,NROW,TLargest,tDie2,
     &      tFail)
          if(tFail) return
C ==------------------------------------------------------------------==
C..V=[V_1 V_2.. V_L] Y' 
          CALL DGEMM('N','N',M,NBLEFF,NBLK,1.0_dp,V,M,T,NBLK,0.0_dp,
     &         V0(1,I),M)
C==-------------------------------------------------------------------==
C==  End of refinement over states                                    ==
C==-------------------------------------------------------------------==
        ENDDO
C..GS orthogonalisation on refined states
        CALL NECI_MY_GSORTHO(M,V0,NCURR-1,V0(1,NCURR),NLEFT,A,tDie2,
     &          tFail)
        if(tFail) return
C..HPSI: H V0. Enter only refined states 
        CALL MY_HPSI(M,NLEFT,NROW,LAB,H,
     &        V0(1,NCURR),VS(1,NCURR),TLargest)
C        CALL NECI_HSPI(M,NLEFT,H,V0(1,NCURR),VS(1,NCURR))
        NHPSI=NHPSI+NLEFT 
C..Ovlap: V0^T VS  
        CALL NECI_OVLAP(M,N,A,V0,VS)
C..AY=YE
        CALL DSYEV('V','U',N,A,N,W,SCR,LSCR,INFO)
        CALL NECI_REORDER(N,N,W,A)
C..Rotate: V -> VY 
        CALL NECI_ROTATE(M,N,V0,A,SCR,M*N)
C..Rotate: HV -> HVY 
        CALL NECI_ROTATE(M,N,VS,A,SCR,M*N)
C ==============================================================
      ENDDO                     ! Loop over icycle 
C.. End of Lanczos diagonalisation  
 100  CONTINUE
      IF(PRINTOUT) THEN
          WRITE(6,'(//''    NCONV:'',I5)') NCONV
      ENDIF
      T2 = neci_etime(tarr) 
      T3=(T2-T1)
      IF(PRINTOUT) THEN
          WRITE(6,'(//5X,''TIME FOR LANCZOS DIAGONALISATION'',F10.2)')
     &     T3/1000.0_dp
      ENDIF
      call halt_timer(proc_timer)
C     ================================================================  
      RETURN
      END 
C ======================================================================
      SUBROUTINE NECI_MGS(M,N,A,LDA,R,LDR,tDie2,tFail)
C     ==--------------------------------------------------------------==
      use global_utilities
      use constants, only : dp,sp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: LDA,LDR,N,M
      integer(sp) :: INFO
      real(dp) :: A(LDA,*),R(LDR,*)
      LOGICAL :: tDie2,tfail
      type(timer), save :: proc_timer
C     ==================================================================
      proc_timer%timer_name='NECI_MGS'
      call set_timer(proc_timer)
      CALL NECI_RGS(M,N,A,M,R,tDie2,tFail)
      if(tFail) return
      CALL DTRTRI('U','N',N,R,N,INFO)
      call halt_timer(proc_timer)
      RETURN 
      END 
C =======================================================================
      SUBROUTINE NECI_SETUP_MATRIX(M,N,A,TSYM,TSeeded)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only : dp
      IMPLICIT NONE
C..M = NDET
C..N = NEVAL
      integer :: M,N,I,J
      real(dp) :: A(M,N),RAN2
      LOGICAL TSYM,TSeeded 
      integer iseed
      IF(TSeeded) THEN
          iseed=0
      ELSE
          iseed=1
      ENDIF
      IF(TSYM) THEN 
        DO I=1,N
          DO J=I,M 
            A(J,I)=RAN2(iseed)
            A(I,J)=A(J,I)
          ENDDO
        ENDDO
      ELSE       
        DO I=1,N
          DO J=1,M 
C            IF(J.EQ.NUN) THEN
C              IF(I.EQ.1) THEN
C                A(J,I)=1.0_dp
C              ELSE
C               A(J,I)=0.0_dp
C              ENDIF
C            ENDIF
            A(J,I)=RAN2(iseed)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
C =======================================================================
      SUBROUTINE NECI_WRITE_MATRIX(CHAR,M,N,A)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      CHARACTER(*) CHAR
      INTEGER :: N,M,I,J
      real(dp) :: A(M,N)
      WRITE(6,*) CHAR
      DO I=1,M
        WRITE(6,1000) (A(I,J),J=1,N)
      ENDDO
 1000 FORMAT(12E15.6)
      RETURN
      END
C ======================================================================
      SUBROUTINE NECI_PUTTMAT(CHAR,T,M,N,A,MA,NA,IBEG,JBEG)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: N,M,MA,NA,IBEG,JBEG,I,J
      real(dp) :: T(M,N),A(MA,NA)
      CHARACTER(1) CHAR
      IF(MA*IBEG.GT.M) STOP ' MA+IBEG.GT.M'
      IF(NA*JBEG.GT.N) STOP ' MA+IBEG.GT.M'
      IF(CHAR.EQ.'N') THEN 
        DO J=1,NA
          DO I=1,MA 
            T(I+MA*(IBEG-1),J+NA*(JBEG-1))=A(I,J)
          ENDDO
        ENDDO
      ELSEIF(CHAR.EQ.'T') THEN 
        DO J=1,NA
          DO I=1,MA 
            T(I+MA*(IBEG-1),J+NA*(JBEG-1))=A(J,I)
          ENDDO
        ENDDO
      ELSE
        STOP 'ILLEGAL CHAR'
      ENDIF
      RETURN
      END
C ======================================================================
      SUBROUTINE NECI_GETMAT(CHAR,T,M,N,A,MA,NA,IBEG,JBEG)
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: N,M,MA,NA,IBEG,JBEG,I,J
      real(dp) :: T(M,N),A(MA,NA)
      CHARACTER(1) :: CHAR
      IF(MA*IBEG.GT.M) STOP ' MA+IBEG.GT.M'
      IF(NA*JBEG.GT.N) STOP ' MA+IBEG.GT.M'
      IF(CHAR.EQ.'N') THEN 
        DO J=1,NA
          DO I=1,MA 
            A(I,J)=T(I+MA*(IBEG-1),J+NA*(JBEG-1))
          ENDDO
        ENDDO
      ELSEIF(CHAR.EQ.'T') THEN 
        DO J=1,NA
          DO I=1,MA 
            A(J,I)=T(I+MA*(IBEG-1),J+NA*(JBEG-1))
          ENDDO
        ENDDO
      ELSE
        STOP 'ILLEGAL CHAR'
      ENDIF
      RETURN
      END
C ======================================================================
      SUBROUTINE NECI_REORDER(M,N,W,A)
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: N,M,J
      real(dp) :: W(N),A(M,N)
      real(dp) :: AUX
      DO J=1,N/2
        CALL DSWAP(M,A(1,J),1,A(1,N-J+1),1)
        AUX=W(N-J+1)
        W(N-J+1)=W(J)
        W(J)=AUX
      ENDDO
      RETURN 
      END
C ======================================================================
      SUBROUTINE NECI_KRYREF(M,N,NKRY,NCONV,V0,V,AM,BM,T,LDT,WT,
     &       LDWT,H,SCR,LSCR,ISCR,LISCR,NHPSI,LAB,NROW,TLargest,tDie2
     &       ,tFail)
C     ==================================================================
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: N,M,LAB(M,*),NROW(M),LSCR,LDWT,LDT,NKRY
      integer :: LISCR,J,NCONV,NHPSI,I,ISCR(LISCR)
      real(dp) :: V0(M,*),V(M,N,NKRY+1),AM(N,N,NKRY+1),BM(N,N,NKRY)
      real(dp) :: T(LDT),WT(LDWT),H(*),SCR(LSCR)
      LOGICAL :: TLargest,tDie2,tFail
      type(timer), save :: proc_timer
C     ==================================================================
      proc_timer%timer_name='    KRYREF'
      call set_timer(proc_timer)
      DO I=3,NKRY
C..V_I = V_I - V_(I-1) A_(I-1) 
        CALL NECI_RSDBLK('N',M,N,V(1,1,I-1),AM(1,1,I-1),V(1,1,I))
C..V_I B_I = V_I
        CALL NECI_MGS(M,N,V(1,1,I),M,BM(1,1,I),N,tDie2,tFail)
        if(tFail) return
C..HPSI:  H V_I
        CALL MY_HPSI(M,N,NROW,LAB,H,V(1,1,I),V(1,1,I+1),TLargest)
C        CALL NECI_HSPI(M,N,H,V(1,1,I),V(1,1,I+1))
        NHPSI=NHPSI+N 
C.. project out converged states
        CALL NECI_PRJCNV(M,N,NCONV,V0,V(1,1,I+1),SCR)
C..V_(I+1) = H V_I - V_(I-1)B_I 
        CALL NECI_RSDBLK('T',M,N,V(1,1,I-1),BM(1,1,I),V(1,1,I+1))
C..Ovlap:  A_I=V_I^T V_(I+1)
        CALL NECI_OVLAP(M,N,AM(1,1,I),V(1,1,I),V(1,1,I+1))
C..Setup T-matrix. First diagonal terms
C        CALL AZZERO(T,N*NKRY*N*NKRY) 
C        CALL NECI_PUTTMAT('N',T,I*N,I*N,AM(1,1,1),N,N,1,1)
C        DO J=2,I
C          CALL NECI_PUTTMAT('N',T,I*N,I*N,AM(1,1,J),N,N,J,J)
C          CALL NECI_PUTTMAT('T',T,I*N,I*N,BM(1,1,J),N,N,J-1,J)
C          CALL NECI_PUTTMAT('N',T,I*N,I*N,BM(1,1,J),N,N,J,J-1)
C        ENDDO  
C..TY=YE
C        CALL DSYEV('V','U',I*N,T,I*N,WT,SCR,LSCR,INFO)
C        CALL NECI_REORDER(I*N,WT,T)
C.. Compute error V_3 Y'', where Y'' is the last N rows of Y 
C        CALL GETMAT('N',T,I*N,I*N,AM(1,1,I+1),N,N,I,1)
C        CALL DGEMM('N','N',M,N,N,1.0_dp,V(1,1,I+1),M,AM(1,1,I+1),
C     &       N,0.0_dp,SCR,M)
C        DO J=1,N
C          AUX=DNRM2(M,SCR(M*(J-1)+1),1)**2
C          IF(AUX.LT.B2MIN) B2MIN=AUX
C          IF(AUX.GT.B2MAX) B2MAX=AUX
C        ENDDO
      ENDDO                     !  End loop over ikry
C ==================================================================
C..Setup T-matrix. First diagonal terms      
      I=NKRY
      T(1:(I*N)**2)=0
      CALL NECI_PUTTMAT('N',T,I*N,I*N,AM(1,1,1),N,N,1,1)
      DO J=2,I
        CALL NECI_PUTTMAT('N',T,I*N,I*N,AM(1,1,J),N,N,J,J)
        CALL NECI_PUTTMAT('T',T,I*N,I*N,BM(1,1,J),N,N,J-1,J)
        CALL NECI_PUTTMAT('N',T,I*N,I*N,BM(1,1,J),N,N,J,J-1)
      ENDDO 
      IF(N.GT.0) THEN 
C..Full matrix diag. 
        CALL NECI_JACOBI(I*N,N,T,WT,SCR,LSCR,ISCR,5*I*N,ISCR(5*I*N+1))
      ELSE
C..banded matrix diagonalisation. This seems to be slower than 
C..full diag.  
        CALL NECI_BANDM(I*N,N,T,WT,SCR,LSCR,ISCR,5*I*N,ISCR(5*I*N+1))
      ENDIF
C..
C ===================================================================
      call halt_timer(proc_timer)
      RETURN
      END
C ==================================================================
      SUBROUTINE NECI_PRPKRV(M,N,NKRY,NCONV,VCONV,V0,VS,V,AM,BM,H,W,
     &     SCR,LSCR,NHPSI,LAB,NROW,TLargest,tDie2,tFail)
C ==------------------------------------------------------------------==
C == Returns A_1,A_2,B_2 and V in a form suitable for krylov          ==
C == refinement                                                       ==
C ==------------------------------------------------------------------==
      use constants, only: dp
      use global_utilities
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: M,N,LAB(M,*),NROW(M),LSCR,I,NHPSI,NKRY,NCONV
      real(dp) :: VCONV(M,*),V0(M,N),VS(M,N),V(M,N,NKRY+1)
      real(dp) :: AM(N,N,NKRY+1),BM(N,N,NKRY),W(N),H(*),SCR(LSCR)
      LOGICAL :: TLargest,tDie2,tFail
      type(timer), save :: proc_timer
C ==------------------------------------------------------------------==
      proc_timer%timer_name='NECI_PRPKRV'
      call set_timer(proc_timer)
C ==------------------------------------------------------------------==
      CALL DCOPY(M*N,V0,1,V(1,1,1),1)
      CALL DCOPY(M*N,VS,1,V(1,1,2),1)
C..Compute residual 
      DO I=1,N
        CALL DAXPY(M,-W(I),V(1,I,1),1,V(1,I,2),1)
      ENDDO
C..V2 B2=R  
      BM=0.0_dp
      CALL NECI_MGS(M,N,V(1,1,2),M,BM(1,1,2),N,tDie2,tFail)
      if(tFail) return
C..   Setup A_1=diag(w)
      AM(:,:,1)=0.0_dp
      DO I=1,N 
        AM(I,I,1)=W(I)
      ENDDO
C..V_3=H.V2
      CALL MY_HPSI(M,N,NROW,LAB,H,V(1,1,2),V(1,1,3),TLargest)
C      CALL NECI_HSPI(M,N,H,V(1,1,2),V(1,1,3))
      NHPSI=NHPSI+N 
C.. project out converged states
        CALL NECI_PRJCNV(M,N,NCONV,VCONV,V(1,1,3),SCR)
C..V_3=V_3-V_1 B_2^T
      CALL NECI_RSDBLK('T',M,N,V(1,1,1),BM(1,1,2),V(1,1,3))
C..Ovlap 
      CALL NECI_OVLAP(M,N,AM(1,1,2),V(1,1,2),V(1,1,3))
C ==------------------------------------------------------------------==
      call halt_timer(proc_timer)
      RETURN
      END
C     ==================================================================
      SUBROUTINE NECI_MY_GSORTHO(M,C0,N0,CP,NP,SMAT,tDie2,tFail)
C     ==--------------------------------------------------------------==
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
C     Arguments
      LOGICAL :: tDie2,tFail
      INTEGER    N0,NP,M
      real(dp)     C0(M,*),CP(M,*)
      real(dp)     SMAT(NP,*)             !MAX(N0,NP)
C     Variables
      type(timer), save :: proc_timer
C     ==--------------------------------------------------------------==
      proc_timer%timer_name='MY_GSORTHO'
      call set_timer(proc_timer)
      IF(NP.LT.1) RETURN
      IF(N0.GT.0) THEN
C..SMAT=CP.C0 
        CALL DGEMM('T','N',NP,N0,M,1.0_dp,CP,M,C0,M,0.0_dp,SMAT,NP)
C..CP -> CP-C0*SMAT
        CALL DGEMM('N','T',M,NP,N0,-1.0_dp,C0,M,SMAT,NP,
     *             1.0_dp,CP,M)
      ENDIF
      CALL NECI_MGS(M,NP,CP,M,SMAT,NP,tDie2,tFail)
      call halt_timer(proc_timer)
C     ==--------------------------------------------------------------==
      RETURN
      END
C =======================================================================
      SUBROUTINE NECI_PUTTAB(N,KD,A,SCR,LSCR)
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: N,LSCR,KD,I,IBEG,J
      real(dp) :: A(N,N),SCR(LSCR)
C..
      IF(LSCR.LT.N) THEN 
        WRITE(6,*) ' LSCR:',LSCR
        WRITE(6,*) 'N:',N
        STOP 'LSCR LT N '
      ENDIF
      SCR(1:N)=0.0_dp
      DO J=1,N
        CALL DCOPY(N,A(1,J),1,SCR,1)
        A(1:N,J)=0.0_dp
        IBEG=MAX(J-KD,1)
        DO I=IBEG,J
          A(1+KD+I-J,J)=SCR(I)
        ENDDO
      ENDDO
      RETURN
      END 
C =======================================================================
      SUBROUTINE NECI_HPSI(M,N,H,V0,VS)
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp)(A-H,O-Z)
      IMPLICIT NONE
      integer :: N,M
      real(dp) :: H(M,M),V0(M,N),VS(M,N)
      type(timer), save :: proc_timer
      proc_timer%timer_name='      NECI_HSPI'
      call set_timer(proc_timer)
      CALL DGEMM('N','N',M,N,M,1.0_dp,H,M,V0,M,0.0_dp,VS,M)
      call halt_timer(proc_timer)
      RETURN
      END 
C =======================================================================
      SUBROUTINE NECI_OVLAP(M,N,A,V0,VS)
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp)(A-H,O-Z)
      IMPLICIT NONE
      integer :: N,M
      real(dp) :: A(N,N),V0(M,N),VS(M,N)
      type(timer), save :: proc_timer
      proc_timer%timer_name='NECI_OVLAP'
      call set_timer(proc_timer)
      CALL DGEMM('T','N',N,N,M,1.0_dp,V0,M,VS,M,0.0_dp,A,N)
      call halt_timer(proc_timer)
      RETURN
      END 
C =======================================================================
      SUBROUTINE NECI_RSDBLK(CHAR,M,N,V1,A,V2)
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp)(A-H,O-Z)
      IMPLICIT NONE
      CHARACTER(1) CHAR 
      integer :: N,M
      real(dp) :: A(N,N),V1(M,N),V2(M,N)
      type(timer), save :: proc_timer
C..V_I = V_I - V_(I-1) A_(I-1)
      proc_timer%timer_name='    NECI_RSDBLK'
      call set_timer(proc_timer)
      IF(CHAR.EQ.'N') THEN 
        CALL DGEMM('N','N',M,N,N,-1.0_dp,V1,M,A,N,1.0_dp,V2,M)
      ELSEIF(CHAR.EQ.'T') THEN 
        CALL DGEMM('N','T',M,N,N,-1.0_dp,V1,M,A,N,1.0_dp,V2,M)
      ELSE
        STOP ' CHAR ILLEGAL IN NECI_RSDBLK ' 
      ENDIF
      call halt_timer(proc_timer)
      RETURN
      END
C =======================================================================
      SUBROUTINE NECI_ROTATE(M,N,V0,A,SCR,LSCR)
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp)(A-H,O-Z)
      IMPLICIT NONE
      integer :: N,M,LSCR
      REAL(dp) :: A(N,N),V0(M,N),SCR(LSCR)
      type(timer), save :: proc_timer
      proc_timer%timer_name='    NECI_ROTATE'
      call set_timer(proc_timer)
      CALL DGEMM('N','N',M,N,N,1.0_dp,V0,M,A,N,0.0_dp,SCR,M)
      CALL DCOPY(M*N,SCR,1,V0,1)
      call halt_timer(proc_timer)
      RETURN 
      END 


      SUBROUTINE NECI_JACOBI(M,N,T,WT,SCR,LSCR,ISCR,LISCR,IFAIL)
C ==------------------------------------------------------------------==
C Returns the n largest eigenvalues of MxM symtric matrix T           ==
C ==------------------------------------------------------------------==
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER :: IFAIL(*),INFO,M,N,LISCR,LSCR,MEVAL,ISCR(LISCR)
      REAL(dp) :: T(M,M,*),WT(*),SCR(LSCR)
      type(timer), save :: proc_timer
      proc_timer%timer_name='    JACOBI'
      call set_timer(proc_timer)
      IF(M.EQ.N) THEN 
        CALL DSYEV('V','U',N,T,N,WT,SCR,LSCR,INFO)
      ELSE
        CALL DSYEVX('V','I','U',M,T,M,1.0_dp,1.0_dp,M-N+1,M,0.0_dp,
     &       MEVAL,WT,T(1,1,2),M,SCR,LSCR,ISCR,IFAIL,INFO)
        IF(MEVAL.LT.N) THEN 
          WRITE(6,*) ' WARNING| DSYEVX RETURNED MEVAL < N',MEVAL,N
        ENDIF
        CALL DCOPY(M*N,T(1,1,2),1,T(1,1,1),1)
      ENDIF
      CALL NECI_REORDER(M,N,WT,T)
      call halt_timer(proc_timer)
      RETURN 
      END 
C     ==================================================================
      SUBROUTINE NECI_BANDM(IN,N,T,WT,SCR,LSCR,ISCR,LISCR,IFAIL)
C ==------------------------------------------------------------------==
C Returns the n largest eigenvalues of MxM symtric banded matrix T    ==
C ==------------------------------------------------------------------==
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      type(timer), save :: proc_timer
      INTEGER :: IFAIL(*),IN,LSCR,LISCR,N,MEVAL,INFO,ISCR(LISCR)
      REAL(dp) :: T(IN,IN,*),WT(IN),SCR(LSCR)
C     ==================================================================
      proc_timer%timer_name='     BANDM'
      call set_timer(proc_timer)
C..Put T is form suitable for banded matrix diagonalisation 
      CALL NECI_PUTTAB(IN,N,T,SCR,LSCR)
      IF(LSCR.LT.7*IN) THEN
        WRITE(6,*) ' 7*IN:',7*IN
        WRITE(6,*) 'LSCR:',LSCR
        STOP ' LSCR TOO SMALL IN BANDM'
      ENDIF
      IF(LISCR.LT.5*IN) THEN
        WRITE(6,*) ' 5*IN:',5*IN
        WRITE(6,*) 'LISCR:',LISCR
        STOP ' LISCR TOO SMALL IN BANDM'
      ENDIF
C        CALL DSBEV('V','U',IN,N,T,IN,WT,T(1,2),IN,SCR,INFO)
      CALL DSBEVX('V','I','U',IN,N,T,IN,T(1,1,2),IN,1.0_dp,1.0_dp,
     &     IN-N+1,IN,1.D-10,MEVAL,WT,T(1,1,3),IN,SCR,ISCR,IFAIL,INFO)
      IF(MEVAL.LT.N) THEN 
        WRITE(6,*) ' WARNING| DSBEVX RETURNED MEVAL < N',MEVAL,N
      ENDIF
      CALL DCOPY(IN*IN,T(1,1,3),1,T(1,1,1),1)
      CALL NECI_REORDER(IN,N,WT,T)
      call halt_timer(proc_timer)
      RETURN 
      END
C =====================================================================
      SUBROUTINE NECI_PRJCNV(M,NP,N0,V0,V,SMAT)
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp)(A-H,O-Z)
      IMPLICIT NONE
      INTEGER :: M,N,N0,NP
      REAL(dp) :: V(M,NP),V0(M,*),SMAT(N0,*)
      type(timer), save :: proc_timer
      IF(N0.EQ.0) RETURN
      proc_timer%timer_name='    PRJCNV'
      call set_timer(proc_timer)
C..SMAT=V0^T.V 
      CALL DGEMM('T','N',N0,NP,M,1.0_dp,V0,M,V,M,0.0_dp,SMAT,N0)
Cdeb      DO J=1,N0
Cdeb        CALL DSCAL(NP,W(J),SMAT(J,1),N0)
Cdeb      ENDDO
C..   V -> V-V0*SMAT
      CALL DGEMM('N','N',M,NP,N0,-1.0_dp,V0,M,SMAT,N0,1.0_dp,V,M)
      call halt_timer(proc_timer)
      RETURN 
      END 
C     =================================================================
      SUBROUTINE NECI_RGS(M,N,CP,LDCP,SMAT,tDie2,tFail)
C     ==--------------------------------------------------------------==
C     ==  GRAM-SCHMIDT ORTHOGONALIZATION                              ==
C     ==--------------------------------------------------------------==
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
C     Arguments
      INTEGER    N,M,LDCP
      real(dp)     SMAT(N,N)
      real(dp)     CP(LDCP,N)
      logical :: tDie2,tFail
C     VARIABLES
      type(timer), save :: proc_timer
C     ==--------------------------------------------------------------==
      IF(N.LE.0) RETURN
      proc_timer%timer_name='NECI_RGS'
      call set_timer(proc_timer)
      SMAT=0
      CALL DSYRK('U','T',N,M,1.0_dp,CP,M,0.0_dp,SMAT,N)
      CALL NECI_UINV('U',SMAT,N,N,tDie2,tFail)
      if(tFail) return
      CALL DTRMM('R','U','N','N',M,N,1.0_dp,SMAT,N,CP,M)
      call halt_timer(proc_timer)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE NECI_UINV(UPLO,SMAT,LDA,N,tDie2,tFail)
C     ==--------------------------------------------------------------==
      use constants, only: dp,sp
      IMPLICIT NONE
C     ARGUMENTS
      CHARACTER(1) UPLO
      INTEGER   LDA,N
      real(dp)    SMAT(LDA,N)
C     VARIABLES
      INTEGER(sp)   INFO
      LOGICAL :: tDie2,tFail
C     ==--------------------------------------------------------------==
      CALL DPOTRF(UPLO,N,SMAT,LDA,INFO)
      IF(INFO.NE.0) THEN
          WRITE(6,*) "INFO is : ", INFO
      ENDIF
      IF(INFO.NE.0) THEN
          if(tDie2) then
              CALL Stop_All('UINV','ILLEGAL RESULTS DGETRF')
          else
              CALL Warning_neci('UINV','ILLEGAL RESULTS DGETRF')
              tFail=.true.
              return
          endif
      endif
      CALL DTRTRI(UPLO,'N',N,SMAT,LDA,INFO)
      IF(INFO.NE.0) then
          if(tDie2) then
              CALL Stop_All('UINV','ILLEGAL RESULTS DTRTRI')
          else
              CALL Warning_neci('UINV','ILLEGAL RESULTS DTRTRI')
              tFail = .true.
          endif
      endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
