      SUBROUTINE FRSBLKH(M,ICMAX,N,H,LAB,V0,VS,NKRY,NKRY1,
     &	NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT,SCR,ISCR,INDEX,
     &	WH,WORK2,V2,NCYCLE,B2LIMIT)
      IMPLICIT REAL*8 (A-H,O-Z) 
C ==----------------------------------------------------------------==
C ==	M is NDET : the number of rows in Hamil			    ==
C ==    ICMAX : the MAX. number of columns in Hamil		    ==
C ==    N is NEVAL : the number of columns in the wavevector        ==
C ==----------------------------------------------------------------==
      DIMENSION V0(M,N),VS(M,N),A(N,N),W(N) 
      DIMENSION V(M*NBLOCK*NKRY1)
      DIMENSION AM(NBLOCK*NBLOCK*NKRY1),BM(NBLOCK*NBLOCK*NKRY)
      DIMENSION T(3*NBLOCK*NKRY*NBLOCK*NKRY),WT(NBLOCK*NKRY)
      DIMENSION SCR(LSCR),ISCR(LISCR)
      DIMENSION INDEX(N)
C..temp
      DIMENSION H(M,ICMAX),WH(M),WORK2(3*M)
      DIMENSION V2(M,N)  
C..temp
      DIMENSION LAB(M,ICMAX),NROW(M)
C ==----------------------------------------------------------------==
C..Initialise V0 
      CALL SETUP_MATRIX(M,N,V0,.FALSE.)
      CALL MGS(M,N,V0,M,A,N)
C..Calculate top 80% of N 
      NDIAG=N-INT(N*0.2) 
C..
      WRITE(6,20000) M,N,NKRY,NBLOCK,NDIAG,LSCR,B2LIMIT
20000 FORMAT(2X,'M:',I7/2X,'N:',I7/2X,'NKRY:',I7/2X,'NBLOCK:',I7/2X,
     &	'NDIAG:',I7/2X,'LSCR:',I7/2X,'B2LIMIT:',E10.2)
C ==----------------------------------------------------------------==
      CALL FRSBLK(M,N,NKRY,NBLOCK,V0,VS,A,V,AM,BM,T,W,WT,INDEX,
     &     SCR,LSCR,ISCR,LISCR,NDIAG,B2LIMIT,H,ICMAX,LAB,NROW,NCYCLE)
C ==----------------------------------------------------------------==
C..Uncomment to test to see if routine is working
C      CALL TIPRI
C..Exact eigenstates
C      T1 = TIMEC()
C      CALL DSYEV('V','U',M,H,M,WH,WORK2,3*M,INFO)
C      WRITE(6,'(//14X,''Exact'',15X,''Lanczos'',10X,''Residual'')')
C      DO I=1,N
C        AUX=ABS(DDOT(M,V0(1,I),1,H(1,M-I+1),1))
C        AUX=1.D0-AUX
C        WRITE(6,'(6X,I3,2E19.11,2X,E10.3)') I,WH(M-I+1),W(I),AUX
C      ENDDO
C      T2 = TIMEC()
C      T3=(T2-T1)
C      WRITE(6,'(//5X,''TIME FOR EXACT DIAGONALISATION'',F10.2)')
C     &	   T3/1000.D0
      WRITE(6,'(//10X,''Neval'',15X,''Eigenvalue'')')
      DO I=1,N
        WRITE(6,'(10X,I3,15X,F19.11)') I,-1.D0*W(I)
      ENDDO
C ==----------------------------------------------------------------==
      END
C ====================================================================
      SUBROUTINE FRSBLK(M,N,NKRY,NBLOCK,V0,VS,A,V,AM,BM,T,W,WT,INDEX,
     &     SCR,LSCR,ISCR,LISCR,NDIAG,B2LIMIT,H,ICMAX,LAB,NROW,NCYCLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V0(M,N),VS(M,N),A(N,N),W(N) 
      DIMENSION V(M*NBLOCK*(NKRY+1))
      DIMENSION AM(NBLOCK*NBLOCK*(NKRY+1)),BM(NBLOCK*NBLOCK*NKRY)
      DIMENSION T(3*NBLOCK*NKRY*NBLOCK*NKRY),WT(NBLOCK*NKRY)
      DIMENSION SCR(LSCR),ISCR(LISCR)
      DIMENSION INDEX(N),NROW(M),LAB(M,ICMAX) 
C..temp
      DIMENSION H(*)
C..temp
C ==-----------------------------------------------------------------==
      CALL TISET('    FRSBLK',ISUB)
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
      T1 = TIMEC() 
C ====================================================================
      CALL IAZZERO(INDEX,N)
      NHPSI=0
      WRITE(6,10000) ICYCLE,NCONV,B2MAX,B2MIN,NHPSI
10000 FORMAT(2X,'ICYCLE:',I3,1X,'NCONV:',I3,2X,'B2MAX:',F10.5,
     &		6X,'B2MIN:',F10.5,5X,'NHPSI:',I3)
C..VS=H.V0 
C..My matrix multiplication routine implimented 8/11/02 DCT
      CALL MY_HPSI(M,ICMAX,N,NROW,LAB,H,V0,VS)
C      CALL HPSI(M,N,H,V0,VS)
      NHPSI=NHPSI+N 
C..   Ovlap: V0^T VS. 
      CALL OVLAP(M,N,A,V0,VS)
C..   AY=YE
      CALL DSYEV('V','U',N,A,N,W,SCR,LSCR,INFO)
      CALL REORDER(N,N,W,A)
C..Rotate: V0 -> V0.Y 
      CALL ROTATE(M,N,V0,A,SCR,M*N)
C..Rotate: VS -> VS.Y 
      CALL ROTATE(M,N,VS,A,SCR,M*N)
C=====================================================================
      NCONV=0 
      DO ICYCLE=1,NCYCLE
C        CALL WRITE_MATRIX('   W:   ',N,1,W)
C ====================================================================
C..Residual: H(VY) - (VY)E and test for convergence
        B2MAX=0.D0
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
        WRITE(6,'(5X,I4,2X,I4,2X,2(E8.3,3X),F6.2)') 
     &       ICYCLE,NCONV,B2MAX,B2MIN,DFLOAT(NHPSI)/DFLOAT(N)
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
          CALL PRPKRV(M,NBLEFF,NKRY,NCONV,V0(1,1),V0(1,I),VS(1,I),V,
     &         AM,BM,H,W(I),SCR,LSCR,NHPSI,ICMAX,LAB,NROW)

C Refinement Loop
C ==------------------------------------------------------------------==
          NBLK=NBLEFF*NKRY
          CALL KRYREF(M,NBLEFF,NKRY,NCONV,V0,W,V,AM,BM,T,NBLK*NBLK,WT,
     &         NBLK,H,SCR,LSCR,ISCR,LISCR,NHPSI,ICMAX,LAB,NROW)
C ==------------------------------------------------------------------==
C..V=[V_1 V_2.. V_L] Y' 
          CALL DGEMM('N','N',M,NBLEFF,NBLK,1.D0,V,M,T,NBLK,0.D0,
     &         V0(1,I),M)
C==-------------------------------------------------------------------==
C==  End of refinement over states                                    ==
C==-------------------------------------------------------------------==
        ENDDO
C..GS orthogonalisation on refined states
        CALL MY_GSORTHO(M,V0,NCURR-1,V0(1,NCURR),NLEFT,A)
C..HPSI: H V0. Enter only refined states 
        CALL MY_HPSI(M,ICMAX,NLEFT,NROW,LAB,H,
     &		V0(1,NCURR),VS(1,NCURR))
C        CALL HPSI(M,NLEFT,H,V0(1,NCURR),VS(1,NCURR))
        NHPSI=NHPSI+NLEFT 
C..Ovlap: V0^T VS  
        CALL OVLAP(M,N,A,V0,VS)
C..AY=YE
        CALL DSYEV('V','U',N,A,N,W,SCR,LSCR,INFO)
        CALL REORDER(N,N,W,A)
C..Rotate: V -> VY 
        CALL ROTATE(M,N,V0,A,SCR,M*N)
C..Rotate: HV -> HVY 
        CALL ROTATE(M,N,VS,A,SCR,M*N)
C ==============================================================
      ENDDO                     ! Loop over icycle 
C.. End of Lanczos diagonalisation  
 100  CONTINUE 
      WRITE(6,'(//''    NCONV:'',I5)') NCONV
      T2 = TIMEC() 
      T3=(T2-T1)
      WRITE(6,'(//5X,''TIME FOR LANCZOS DIAGONALISATION'',F10.2)')
     &     T3/1000.D0
      CALL TIHALT('    FRSBLK',ISUB)
C     ================================================================  
      RETURN
      END 
C ======================================================================
      SUBROUTINE MGS(M,N,A,LDA,R,LDR)
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(LDA,*),R(LDR,*)
C     ==================================================================
      CALL TISET('       MGS:',ISUB)
      CALL RGS(M,N,A,M,R,N,I,1)
      CALL DTRTRI('U','N',N,R,N,INFO)
      CALL TIHALT('       MGS:',ISUB)
      RETURN 
      END 
C =======================================================================
      SUBROUTINE SETUP_MATRIX(M,N,A,TSYM)
      IMPLICIT REAL*8 (A-H,O-Z)
C..M = NDET
C..N = NEVAL
      DIMENSION A(M,N)
      LOGICAL TSYM 
      IF(TSYM) THEN 
        DO I=1,N
          DO J=I,M 
            A(J,I)=RAND()
            A(I,J)=A(J,I)
          ENDDO
        ENDDO
      ELSE       
        DO I=1,N
          DO J=1,M 
C            IF(J.EQ.NUN) THEN
C              IF(I.EQ.1) THEN
C                A(J,I)=1.D0
C              ELSE
C               A(J,I)=0.D0
C              ENDIF
C            ENDIF
            A(J,I)=RAND()
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
C =======================================================================
      SUBROUTINE WRITE_MATRIX(CHAR,M,N,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*10 CHAR
      DIMENSION A(M,N)
      WRITE(6,*) CHAR
      DO I=1,M
        WRITE(6,1000) (A(I,J),J=1,N)
      ENDDO
 1000 FORMAT(12E15.6)
      RETURN
      END
C ======================================================================
      SUBROUTINE PUTTMAT(CHAR,T,M,N,A,MA,NA,IBEG,JBEG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(M,N),A(MA,NA)
      CHARACTER*1 CHAR
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
      SUBROUTINE GETMAT(CHAR,T,M,N,A,MA,NA,IBEG,JBEG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(M,N),A(MA,NA)
      CHARACTER*1 CHAR
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
      SUBROUTINE REORDER(M,N,W,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(N),A(M,N)
      DO J=1,N/2
        CALL DSWAP(M,A(1,J),1,A(1,N-J+1),1)
        AUX=W(N-J+1)
        W(N-J+1)=W(J)
        W(J)=AUX
      ENDDO
      RETURN 
      END
C ======================================================================
      SUBROUTINE KRYREF(M,N,NKRY,NCONV,V0,W,V,AM,BM,T,LDT,WT,
     &	   LDWT,H,SCR,LSCR,ISCR,LISCR,NHPSI,ICMAX,LAB,NROW)
C     ==================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V0(M,*),V(M,N,NKRY+1),AM(N,N,NKRY+1),BM(N,N,NKRY)
      DIMENSION T(LDT,*),WT(LDWT),H(*),W(N),SCR(LSCR),ISCR(LISCR)
      DIMENSION LAB(M,*),NROW(M)
C     ==================================================================
      CALL TISET('    KRYREF',ISUB)
      DO I=3,NKRY
C..V_I = V_I - V_(I-1) A_(I-1) 
        CALL RSDBLK('N',M,N,V(1,1,I-1),AM(1,1,I-1),V(1,1,I))
C..V_I B_I = V_I
        CALL MGS(M,N,V(1,1,I),M,BM(1,1,I),N)
C..HPSI:  H V_I
        CALL MY_HPSI(M,ICMAX,N,NROW,LAB,H,V(1,1,I),V(1,1,I+1))
C        CALL HPSI(M,N,H,V(1,1,I),V(1,1,I+1))
        NHPSI=NHPSI+N 
C.. project out converged states
        CALL PRJCNV(M,N,NCONV,V0,W,V(1,1,I+1),SCR)
C..V_(I+1) = H V_I - V_(I-1)B_I 
        CALL RSDBLK('T',M,N,V(1,1,I-1),BM(1,1,I),V(1,1,I+1))
C..Ovlap:  A_I=V_I^T V_(I+1)
        CALL OVLAP(M,N,AM(1,1,I),V(1,1,I),V(1,1,I+1))
C..Setup T-matrix. First diagonal terms
C        CALL AZZERO(T,N*NKRY*N*NKRY) 
C        CALL PUTTMAT('N',T,I*N,I*N,AM(1,1,1),N,N,1,1)
C        DO J=2,I
C          CALL PUTTMAT('N',T,I*N,I*N,AM(1,1,J),N,N,J,J)
C          CALL PUTTMAT('T',T,I*N,I*N,BM(1,1,J),N,N,J-1,J)
C          CALL PUTTMAT('N',T,I*N,I*N,BM(1,1,J),N,N,J,J-1)
C        ENDDO  
C..TY=YE
C        CALL DSYEV('V','U',I*N,T,I*N,WT,SCR,LSCR,INFO)
C        CALL REORDER(I*N,WT,T)
C.. Compute error V_3 Y'', where Y'' is the last N rows of Y 
C        CALL GETMAT('N',T,I*N,I*N,AM(1,1,I+1),N,N,I,1)
C        CALL DGEMM('N','N',M,N,N,1.D0,V(1,1,I+1),M,AM(1,1,I+1),
C     &       N,0.D0,SCR,M)
C        DO J=1,N
C          AUX=DNRM2(M,SCR(M*(J-1)+1),1)**2
C          IF(AUX.LT.B2MIN) B2MIN=AUX
C          IF(AUX.GT.B2MAX) B2MAX=AUX
C        ENDDO
      ENDDO                     !  End loop over ikry
C ==================================================================
C..Setup T-matrix. First diagonal terms      
      I=NKRY
      CALL AZZERO(T,I*N*I*N) 
      CALL PUTTMAT('N',T,I*N,I*N,AM(1,1,1),N,N,1,1)
      DO J=2,I
        CALL PUTTMAT('N',T,I*N,I*N,AM(1,1,J),N,N,J,J)
        CALL PUTTMAT('T',T,I*N,I*N,BM(1,1,J),N,N,J-1,J)
        CALL PUTTMAT('N',T,I*N,I*N,BM(1,1,J),N,N,J,J-1)
      ENDDO 
      IF(N.GT.0) THEN 
C..Full matrix diag. 
        CALL JACOBI(I*N,N,T,WT,SCR,LSCR,ISCR,5*I*N,ISCR(5*I*N+1))
      ELSE
C..banded matrix diagonalisation. This seems to be slower than 
C..full diag.  
        CALL BANDM(I*N,N,T,WT,SCR,LSCR,ISCR,5*I*N,ISCR(5*I*N+1))
      ENDIF
C..
C ===================================================================
      CALL TIHALT('    KRYREF',ISUB)
      RETURN
      END
C ==================================================================
      SUBROUTINE PRPKRV(M,N,NKRY,NCONV,VCONV,V0,VS,V,AM,BM,H,W,
     &     SCR,LSCR,NHPSI,ICMAX,LAB,NROW)
C ==------------------------------------------------------------------==
C == Returns A_1,A_2,B_2 and V in a form suitable for krylov          ==
C == refinement                                                       ==
C ==------------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VCONV(M,*),V0(M,N),VS(M,N),V(M,N,NKRY+1),
     &     AM(N,N,NKRY+1),BM(N,N,NKRY),SCR(LSCR),W(N),H(*)
      DIMENSION LAB(M,*),NROW(M)
C ==------------------------------------------------------------------==
      CALL TISET('    PRPKRV',ISEUB)
C ==------------------------------------------------------------------==
      CALL DCOPY(M*N,V0,1,V(1,1,1),1)
      CALL DCOPY(M*N,VS,1,V(1,1,2),1)
C..Compute residual 
      DO I=1,N
        CALL DAXPY(M,-W(I),V(1,I,1),1,V(1,I,2),1)
      ENDDO
C..V2 B2=R  
      CALL AZZERO(BM,N*N*NKRY)
      CALL MGS(M,N,V(1,1,2),M,BM(1,1,2),N)
C..   Setup A_1=diag(w)
      CALL AZZERO(AM,N*N)
      DO I=1,N 
        AM(I,I,1)=W(I)
      ENDDO
C..V_3=H.V2
      CALL MY_HPSI(M,ICMAX,N,NROW,LAB,H,V(1,1,2),V(1,1,3))
C      CALL HPSI(M,N,H,V(1,1,2),V(1,1,3))
      NHPSI=NHPSI+N 
C.. project out converged states
        CALL PRJCNV(M,N,NCONV,VCONV,W,V(1,1,3),SCR)
C..V_3=V_3-V_1 B_2^T
      CALL RSDBLK('T',M,N,V(1,1,1),BM(1,1,2),V(1,1,3))
C..Ovlap 
      CALL OVLAP(M,N,AM(1,1,2),V(1,1,2),V(1,1,3))
C ==------------------------------------------------------------------==
      CALL TIHALT('    PRPKRV',ISEUB)
      RETURN
      END
C     ==================================================================
      SUBROUTINE MY_GSORTHO(M,C0,N0,CP,NP,SMAT)
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
C     Arguments
      INTEGER    N0,NP
      REAL*8     C0(M,*),CP(M,*)
      REAL*8     SMAT(NP,*)             !MAX(N0,NP)
C     Variables
      INTEGER    ISUB
C     ==--------------------------------------------------------------==
      CALL TISET('MY_GSORTHO',ISUB)
      IF(NP.LT.1) RETURN
      IF(N0.GT.0) THEN
C..SMAT=CP.C0 
        CALL DGEMM('T','N',NP,N0,M,1.D0,CP,M,C0,M,0.D0,SMAT,NP)
C..CP -> CP-C0*SMAT
        CALL DGEMM('N','T',M,NP,N0,-1.D0,C0,M,SMAT,NP,
     *             1.D0,CP,M)
      ENDIF
      CALL MGS(M,NP,CP,M,SMAT,NP)
      CALL TIHALT('MY_GSORTHO',ISUB)
C     ==--------------------------------------------------------------==
      RETURN
      END
C =======================================================================
      SUBROUTINE PUTTAB(N,KD,A,SCR,LSCR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),SCR(LSCR)
C..
      IF(LSCR.LT.N) THEN 
        WRITE(6,*) ' LSCR:',LSCR
        WRITE(6,*) 'N:',N
        STOP 'LSCR LT N '
      ENDIF
      CALL AZZERO(SCR,N)
      DO J=1,N
        CALL DCOPY(N,A(1,J),1,SCR,1)
        CALL AZZERO(A(1,J),N)
        IBEG=MAX(J-KD,1)
        DO I=IBEG,J
          A(1+KD+I-J,J)=SCR(I)
        ENDDO
      ENDDO
      RETURN
      END 
C =======================================================================
      SUBROUTINE HPSI(M,N,H,V0,VS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION H(M,M),V0(M,N),VS(M,N)
      CALL TISET('      HPSI',ISUB)
      CALL DGEMM('N','N',M,N,M,1.D0,H,M,V0,M,0.D0,VS,M)
      CALL TIHALT('      HPSI',ISUB)
      RETURN
      END 
C =======================================================================
      SUBROUTINE OVLAP(M,N,A,V0,VS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),V0(M,N),VS(M,N)
      CALL TISET('     OVLAP',ISUB)
      CALL DGEMM('T','N',N,N,M,1.D0,V0,M,VS,M,0.D0,A,N)
      CALL TIHALT('     OVLAP',ISUB)
      RETURN
      END 
C =======================================================================
      SUBROUTINE RSDBLK(CHAR,M,N,V1,A,V2)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 CHAR 
      DIMENSION A(N,N),V1(M,N),V2(M,N)
C..V_I = V_I - V_(I-1) A_(I-1)
      CALL TISET('    RSDBLK',ISUB)
      IF(CHAR.EQ.'N') THEN 
        CALL DGEMM('N','N',M,N,N,-1.D0,V1,M,A,N,1.D0,V2,M)
      ELSEIF(CHAR.EQ.'T') THEN 
        CALL DGEMM('N','T',M,N,N,-1.D0,V1,M,A,N,1.D0,V2,M)
      ELSE
        STOP ' CHAR ILLEGAL IN RSDBLK ' 
      ENDIF
      CALL TIHALT('    RSDBLK',ISUB)
      RETURN
      END
C =======================================================================
      SUBROUTINE ROTATE(M,N,V0,A,SCR,LSCR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),V0(M,N),SCR(LSCR)
      CALL TISET('    ROTATE',ISUB)
      CALL DGEMM('N','N',M,N,N,1.D0,V0,M,A,N,0.D0,SCR,M)
      CALL DCOPY(M*N,SCR,1,V0,1)
      CALL TIHALT('    ROTATE',ISUB)
      RETURN 
      END 
C     ==================================================================
      BLOCK DATA TIMER
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
      DATA NTIMER /0/
      DATA TNAM /200*'          '/
      DATA CTSUM,CTSTA,WTSUM,WTSTA,NCALLS
     *     /200*0.0D0,200*0.0D0,200*0.0D0,200*0.0D0,200*0/
      END
C     ==================================================================
      SUBROUTINE TISET(SNAME,NSTIM)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      INTEGER   NSTIM
      CHARACTER SNAME*10
C     Variables
      INTEGER   I
      REAL*8    TIMEC,TIMEF
C     ==--------------------------------------------------------------==
C     ==  START THE CPU AND WALL CLOCK TIME FOR ROUTINE SNAME         ==
C     ==--------------------------------------------------------------==
      DO I=1,NTIMER
        IF(SNAME.EQ.TNAM(I)) THEN
          NSTIM=I
          GOTO 100
        ENDIF
      ENDDO
      NTIMER=NTIMER+1
      IF(NTIMER.GT.NTIMX) THEN
        WRITE(*,*) ' PARAMETER NTIMX IN INCLUDE FILE time.inc TOO ',
     *               'SMALL'
        CALL STOPGM('TISET',' PARAMETER')
      ENDIF
      NSTIM=NTIMER
      TNAM(NSTIM)=SNAME
  100 CONTINUE
      CTSTA(NSTIM)=TIMEC()*0.001D0
      WTSTA(NSTIM)=TIMEF()*0.001D0
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TIHALT(SNAME,NSTIM)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      INTEGER   NSTIM
      CHARACTER SNAME*10
C     Variables
      INTEGER   I
      REAL*8    DCT,DWT,
     &          TIMEC,TIMEF
C     ==--------------------------------------------------------------==
C     ==  STOP THE CPU AND WALL CLOCK TIME FOR ROUTINE SNAME          ==
C     ==--------------------------------------------------------------==
      IF(NSTIM.EQ.0) THEN
        DO I=1,NTIMER
          IF(SNAME.EQ.TNAM(I)) THEN
            NSTIM=I
            GOTO 100
          ENDIF
        ENDDO
        WRITE(*,*) ' NO SUCH TIME HANDLE ',SNAME
        CALL STOPGM('TIHALT',' HANDLE   ')
  100   CONTINUE
      ELSE
        IF(SNAME.NE.TNAM(NSTIM)) THEN
          WRITE(*,*) ' INCONSISTENT HANDLE IN TIHALT ',SNAME,TNAM(NSTIM)
          CALL STOPGM('TIHALT',' HANDLE   ')
        ENDIF
      ENDIF
      DCT=TIMEC()*0.001D0-CTSTA(NSTIM)
      DWT=TIMEF()*0.001D0-WTSTA(NSTIM)
      CTSUM(NSTIM)=CTSUM(NSTIM)+DCT
      WTSUM(NSTIM)=WTSUM(NSTIM)+DWT
      NCALLS(NSTIM)=NCALLS(NSTIM)+1
      DO I=1,NTIMER
        CTSTA(I)=CTSTA(I)+DCT
        WTSTA(I)=WTSTA(I)+DWT
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TIPRI
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Variables
      INTEGER J,IMAX,I
      REAL*8  CPMAX,CPTHRS,TOTALCT,TOTALWT
C     ==--------------------------------------------------------------==
C     ==  PRINT TIMING INFO                                           ==
C     ==--------------------------------------------------------------==
      WRITE(*,'(/,/,1X,64(1H*))')
      WRITE(*,'(" *",62X,"*")')
      WRITE(*,'(" *",27X,A,27X,"*")') ' TIMING '
      WRITE(*,'(" *",62X,"*")')
      WRITE(*,'(1X,64(1H*))')
      WRITE(*,'(A,A)') ' SUBROUTINE',
     *      '            CALLS         CPU TIME        ELAPSED TIME'
      TOTALCT=0.D0
      TOTALWT=0.D0
      DO J=1,NTIMER
        IMAX=0
        CPMAX=-0.1D0
        DO I=1,NTIMER
          IF(CTSUM(I).GT.CPMAX) THEN
            IMAX=I
            CPMAX=CTSUM(I)
          ENDIF
        ENDDO
        IF(J.EQ.1) CPTHRS=CPMAX/1000.D0
Cdeb        IF(J.EQ.1) CPTHRS=CPMAX/100.D0
        IF(CPMAX.LT.CPTHRS) GOTO 100
        I=IMAX
        WRITE(*,'(1X,A10,9X,I8,5X,F12.2,8X,F12.2)') 
     *                    TNAM(I),NCALLS(I),CTSUM(I),WTSUM(I)
        TOTALCT=TOTALCT+CTSUM(I)
        TOTALWT=TOTALWT+WTSUM(I)
        CTSUM(I)=-1.0
      ENDDO
  100 CONTINUE
      WRITE(6,'(1X,64(1H-))')
C      WRITE(6,*) "TOTAL TIME",TOTALCT,TOTALWT
      WRITE(6,'(1X,"TOTAL TIME",T34,F12.2,8X,F12.2)')
     &                    TOTALCT,TOTALWT
C      WRITE(6,'(1X,64(1H*),/)')
C     ==--------------------------------------------------------------==
      WRITE(*,'(1X,64(1H*),/)')
C     ==--------------------------------------------------------------==
      RETURN
      END




C     =============AJWT=========
      REAL*8 FUNCTION  TIGET(SNAME,NSTIM)
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      INTEGER   NSTIM
      CHARACTER SNAME*10
C     Variables
      INTEGER   I
      REAL*8    TIMEC,TIMEF
C     ==--------------------------------------------------------------==
      DO I=1,NTIMER
        IF(SNAME.EQ.TNAM(I)) THEN
          NSTIM=I
          GOTO 100
        ENDIF
      ENDDO
      WRITE(*,*) ' TIMER',SNAME,' NOT FOUND'
      CALL STOPGM('TISET',' PARAMETER')
  100 CONTINUE
      TIGET=CTSUM(NSTIM)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ============/AJWT=====================



C     ==================================================================
C      SUBROUTINE STOPGM(A,B)
C     ==--------------------------------------------------------------==
C      IMPLICIT NONE
C     Arguments
C      CHARACTER*(*) A,B
C     Variables
C     ==---------------------------------------------------------==
C      WRITE(*,'(//,A,A,A,A)') ' PROGRAM STOPS IN SUBROUTINE ',A,'| ',B
C      STOP 999
C     ==---------------------------------------------------------==
C      END
C     =============================================================
      FUNCTION TIMEC()
C     ==---------------------------------------------------------==
C     ==  GENERALL CPU TIME ROUTINE                              ==
C     ==  CPU TIME IS RETURNED IN msec                           ==
C     ==---------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 TIMEC

      REAL*4 ETIME,TARRAY(2),D
      EXTERNAL ETIME
C     ==---------------------------------------------------------==
      D=ETIME(TARRAY)
      TIMEC = DBLE(TARRAY(1))*1000.D0

C     ==---------------------------------------------------------==
      RETURN
      END
C==================================================================
      FUNCTION TIMEF()
C     ==--------------------------------------------------------==
C     ==  GENERAL TIME ROUTINE FOR ELAPSED TIME			==
C     ==  ELAPSED TIME IS RETURNED IN msec			==
C     ==--------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 TIMEF

      REAL*4 TARRAY(2),ETIME
C==-------------------------------------------------------------==
      TIMEF=ETIME(TARRAY)*1000.D0

C==-------------------------------------------------------------==
      RETURN
      END
C=================================================================

      SUBROUTINE JACOBI(M,N,T,WT,SCR,LSCR,ISCR,LISCR,IFAIL)
C ==------------------------------------------------------------------==
C Returns the n largest eigenvalues of MxM symtric matrix T           ==
C ==------------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(M,M,*),WT(*),SCR(LSCR),ISCR(LISCR),IFAIL(*)
      CALL TISET('    JACOBI',ISUB)
      IF(M.EQ.N) THEN 
        CALL DSYEV('V','U',N,T,N,WT,SCR,LSCR,INFO)
      ELSE
        CALL DSYEVX('V','I','U',M,T,M,1.D0,1.D0,M-N+1,M,0.D0,
     &       MEVAL,WT,T(1,1,2),M,SCR,LSCR,ISCR,IFAIL,INFO)
        IF(MEVAL.LT.N) THEN 
          WRITE(6,*) ' WARNING| DSYEVX RETURNED MEVAL < N',MEVAL,N
        ENDIF
        CALL DCOPY(M*N,T(1,1,2),1,T(1,1,1),1)
      ENDIF
      CALL REORDER(M,N,WT,T)
      CALL TIHALT('    JACOBI',ISUB)
      RETURN 
      END 
C     ==================================================================
      SUBROUTINE BANDM(IN,N,T,WT,SCR,LSCR,ISCR,LISCR,IFAIL)
C ==------------------------------------------------------------------==
C Returns the n largest eigenvalues of MxM symtric banded matrix T    ==
C ==------------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(IN,IN,*),WT(IN),SCR(LSCR),ISCR(LISCR),IFAIL(*)
C     ==================================================================
      CALL TISET('     BANDM',ISUB)
C..Put T is form suitable for banded matrix diagonalisation 
      CALL PUTTAB(IN,N,T,SCR,LSCR)
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
      CALL DSBEVX('V','I','U',IN,N,T,IN,T(1,1,2),IN,1.D0,1.D0,
     &     IN-N+1,IN,1.D-10,MEVAL,WT,T(1,1,3),IN,SCR,ISCR,IFAIL,INFO)
      IF(MEVAL.LT.N) THEN 
        WRITE(6,*) ' WARNING| DSBEVX RETURNED MEVAL < N',MEVAL,N
      ENDIF
      CALL DCOPY(IN*IN,T(1,1,3),1,T(1,1,1),1)
      CALL REORDER(IN,N,WT,T)
      CALL TIHALT('     BANDM',ISUB)
      RETURN 
      END
C =====================================================================
      SUBROUTINE PRJCNV(M,NP,N0,V0,W,V,SMAT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(M,NP),V0(M,*),SMAT(N0,*),W(*) 
      IF(N0.EQ.0) RETURN
      CALL TISET('    PRJCNV',ISUB)
C..SMAT=V0^T.V 
      CALL DGEMM('T','N',N0,NP,M,1.D0,V0,M,V,M,0.D0,SMAT,N0)
Cdeb      DO J=1,N0
Cdeb        CALL DSCAL(NP,W(J),SMAT(J,1),N0)
Cdeb      ENDDO
C..   V -> V-V0*SMAT
      CALL DGEMM('N','N',M,NP,N0,-1.D0,V0,M,SMAT,N0,1.D0,V,M)
      CALL TIHALT('    PRJCNV',ISUB)
      RETURN 
      END 
C     =================================================================
      SUBROUTINE RGS(M,N,CP,LDCP,SMAT,ISCR,LSCR)
C     ==--------------------------------------------------------------==
C     ==  GRAM-SCHMIDT ORTHOGONALIZATION                              ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
C     Arguments
      INTEGER    N,ISCR,LSCR(ISCR)
      REAL*8     SMAT(N,N)
      REAL*8     CP(LDCP,N)
C     VARIABLES
      INTEGER    ISUB
C     ==--------------------------------------------------------------==
      IF(N.LE.0) RETURN
      CALL TISET('       RGS',ISUB)
      CALL AZZERO(SMAT,N*N)
      CALL DSYRK('U','T',N,M,1.D0,CP,M,0.D0,SMAT,N)
      CALL UINV('U',SMAT,N,N,LSCR,1)
      CALL DTRMM('R','U','N','N',M,N,1.D0,SMAT,N,CP,M)
      CALL TIHALT('       RGS',ISUB)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE UINV(UPLO,SMAT,LDA,N,ISCR,LSCR)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     ARGUMENTS
      CHARACTER UPLO*1
      INTEGER   LDA,N,LSCR,ISCR(LSCR)
      REAL*8    SMAT(LDA,N)
C     VARIABLES
      INTEGER   INFO
C     ==--------------------------------------------------------------==
      CALL DPOTRF(UPLO,N,SMAT,LDA,INFO)
      IF(INFO.NE.0) CALL STOPGM('UINV','ILLEGAL RESULTS DGETRF')
      CALL DTRTRI(UPLO,'N',N,SMAT,LDA,INFO)
      IF(INFO.NE.0) CALL STOPGM('UINV','ILLEGAL RESULTS DTRTRI')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
