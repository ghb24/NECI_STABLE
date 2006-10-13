C     ==================================================================
#ifndef __HP
      SUBROUTINE AZZERO(A,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      REAL*8 A(N)
C     Variables
#if defined(__VECTOR)
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,N
        A(I)=0.0D0
      ENDDO
#else
      INTEGER I,II
C     ==--------------------------------------------------------------==
      II=MOD(N,4)
      DO I=1,II
        A(I)=0.0D0
      ENDDO
      DO I=1+II,N,4
        A(I)  =0.0D0
        A(I+1)=0.0D0
        A(I+2)=0.0D0
        A(I+3)=0.0D0
      ENDDO
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
      SUBROUTINE IAZZERO(M,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      INTEGER M(N)
C     Variables
#if defined(__VECTOR)
      INTEGER I
#else
      INTEGER I,II
#endif
C     ==--------------------------------------------------------------==
#if defined(__VECTOR)
      DO I=1,N
        M(I)=0
      ENDDO
#else
      II=MOD(N,4)
      DO I=1,II
        M(I)=0
      ENDDO
      DO I=1+II,N,4
        M(I)  =0
        M(I+1)=0
        M(I+2)=0
        M(I+3)=0
      ENDDO
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#if ! defined(__NOINT8)
      SUBROUTINE I8AZZERO(M,N)
C     ==--------------------------------------------------------------==
C     == IAZZERO version for INTEGER*8                                ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      INTEGER*8 M(N)
C     Variables
#if defined(__VECTOR)
      INTEGER I
#else
      INTEGER I,II
#endif
C     ==--------------------------------------------------------------==
#if defined(__VECTOR)
      DO I=1,N
        M(I)=0
      ENDDO
#else
      II=MOD(N,4)
      DO I=1,II
        M(I)=0
      ENDDO
      DO I=1+II,N,4
        M(I)  =0
        M(I+1)=0
        M(I+2)=0
        M(I+3)=0
      ENDDO
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
      FUNCTION XRANF()
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 XRANF
C     Variables
      INTEGER M,KONST
      DATA    M/100001/,KONST/125/
      SAVE    M
C     ==--------------------------------------------------------------==
      M=M*KONST
      M=M-2796203*(M/2796203)
      XRANF=DFLOAT(M)/2796203.D0
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE MPRAND(SEED,N,A)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      REAL*8  A(N),SEED
C     Variables
      REAL*8  XRANF
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,N
        A(I)=XRANF()
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE UNITMX(A,N)
C     ==--------------------------------------------------------------==
C     Arguments
      INTEGER N
      REAL*8  A(N,N)
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      CALL AZZERO(A,N*N)
      DO I=1,N
         A(I,I)=1.D0
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SIMPSN(N,INTE,SUM)
C     ==--------------------------------------------------------------==
C     == COMPUTES ONE-DIMENSIONAL INTEGRALS BY THE SIMPSON METHOD     ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER   N
      REAL*8    INTE(N),SUM
C     Variables
      REAL*8    C1,C2,C3,C4
      INTEGER   I
      PARAMETER (C1=109./48.,C2=-5./48., C3=63./48.,C4=49./48.)
C     ==--------------------------------------------------------------==
      INTE(1)   = INTE(1)*C1
      INTE(2)   = INTE(2)*C2
      INTE(3)   = INTE(3)*C3
      INTE(4)   = INTE(4)*C4
      INTE(N-1) = INTE(N-1)*C1
      INTE(N-2) = INTE(N-2)*C2
      INTE(N-3) = INTE(N-3)*C3
      INTE(N-4) = INTE(N-4)*C4
      SUM=0.D0
      DO I=1,N-1
         SUM=SUM+INTE(I)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE ZCLEAN(A,N,NGW)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,NGW
      REAL*8  A(2,NGW,N)
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,N
         A(2,1,I)=0.D0
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#if !(defined(CRAY))
      SUBROUTINE ICOPY(N,A,IA,B,IB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,IA,IB
      INTEGER A(IA*N),B(IB*N)
C     Variables
      INTEGER I,IAX,IBX
C     ==--------------------------------------------------------------==
      DO I=1,N
        IAX=(I-1)*IA + 1
        IBX=(I-1)*IB + 1
        B(IBX) = A(IAX)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#ifdef POINTER8
      SUBROUTINE I8COPY(N,A,IA,B,IB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER   N,IA,IB
      INTEGER*8 A(IA*N),B(IB*N)
C     Variables
      INTEGER I,IAX,IBX
C     ==--------------------------------------------------------------==
      DO I=1,N
        IAX=(I-1)*IA + 1
        IBX=(I-1)*IB + 1
        B(IBX) = A(IAX)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
#endif
      SUBROUTINE SYMMA(A,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      REAL*8  A(N,N)
C     Variables
      INTEGER I,J
C     ==--------------------------------------------------------------==
      DO I=1,N
        DO J=I+1,N
          A(I,J)=0.5D0*(A(I,J)+A(J,I))
          A(J,I)=A(I,J)
        ENDDO
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE MATMOV(N,M,A,LDA,B,LDB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    N,M,LDA,LDB
      COMPLEX*16 A(LDA,*),B(LDB,*)
C     Variables
#if defined(T3D) || defined(T3E)
      COMPLEX*16 CONE,CZERO
      PARAMETER  (CONE=(1.0D0,0.0D0),CZERO=(0.0D0,0.0D0))
C     ==--------------------------------------------------------------==
      CALL CGESUM('N',N,M,CONE,A,LDA,CZERO,B,LDB)
#elif defined(__VECTOR)
      INTEGER    I,J
C     ==--------------------------------------------------------------==
CMIC$ DO ALL PRIVATE(I,J) SHARED(M,N,A,B)
      DO I=1,M
        DO J=1,N
          B(J,I)=A(J,I)
        ENDDO
      ENDDO
#elif defined(__HP)
      INTEGER    I,J
C     ==--------------------------------------------------------------==
      DO I=1,M
        DO J=1,N
          B(J,I)=A(J,I)
        ENDDO
      ENDDO
#else
      INTEGER    I,J,JJ
C     ==--------------------------------------------------------------==
      DO I=1,M
        JJ=MOD(N,4)
        DO J=1,JJ
          B(J,I)=A(J,I)
        ENDDO
        DO J=1+JJ,N,4
          B(J  ,I)  =  A(J  ,I)
          B(J+1,I)  =  A(J+1,I)
          B(J+2,I)  =  A(J+2,I)
          B(J+3,I)  =  A(J+3,I)
        ENDDO
      ENDDO
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#if !(defined(T3D) || defined(T3E))
      SUBROUTINE CGESUM(TRANSA,N,M,CA,A,LDA,CB,B,LDB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER  TRANSA*1
      INTEGER    N,M,LDA,LDB
      COMPLEX*16 A(LDA,*),B(LDB,*),CA,CB
C     Variables
      INTEGER    I,J
C     ==--------------------------------------------------------------==
      CALL TEST_SCR('CGESUM','LDB',LDB,N)
      IF(TRANSA.EQ.'N'.OR.TRANSA.EQ.'n') THEN
        CALL TEST_SCR('CGESUM','LDA[N]',LDA,N)
        DO I=1,M
          DO J=1,N
            B(J,I)=CB*B(J,I)+CA*A(J,I)
          ENDDO
        ENDDO
      ELSE
        CALL TEST_SCR('CGESUM','LDA[M]',LDA,M)
        DO I=1,M
          DO J=1,N
            B(J,I)=CB*B(J,I)+CA*A(I,J)
          ENDDO
        ENDDO
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
#endif
C     ==================================================================
      SUBROUTINE DATUM(DATX)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER DATX*(*)
C     Variables
#if defined(CRAY)
C     ==--------------------------------------------------------------==
      WRITE(DATX,'(A8,'' ON '',A8)') CLOCK(),DATE()
#elif defined(__NEC)
      CHARACTER*8 A1,A2
C     ==--------------------------------------------------------------==
      CALL DATIM(A1,A2)
      WRITE(DATX,'(A8,'' ON '',A8)') A2,A1
#elif defined(__IBM)
C     ==--------------------------------------------------------------==
      CALL FDATE_(DATX)
#else
C     ==--------------------------------------------------------------==
      WRITE(DATX,'(50X)')
      CALL FDATE(DATX)
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE NUMCPUS(NCPUS)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER NCPUS
C     Variables
      INTEGER IA,IE,MCPUS
      LOGICAL ERROR
      CHARACTER*80 STRING
      SAVE MCPUS
      DATA MCPUS /0/
C     ==--------------------------------------------------------------==
      IF(MCPUS.NE.0) THEN
        NCPUS=MCPUS
      ELSE
        CALL MY_GETENV('NCPUS',STRING)
        CALL XSTRING(STRING,IA,IE)
        IF(IA.GT.IE.OR.STRING(IA:IE).EQ.' ') THEN
          NCPUS=1
        ELSE
          CALL READSI(STRING,IA,IE,NCPUS,ERROR)
        ENDIF
        NCPUS=MAX(NCPUS,1)
        MCPUS=NCPUS
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE INVMAT(N,A,B,INFO)
C     ==--------------------------------------------------------------==
C     == Inverse Matrix A(N,N)                                        ==
C     ==--------------------------------------------------------------==
C     == INPUT:  N Dimension                                          ==
C     ==         A(N,N) Matrix                                        ==
C     == OUTPUT: A(N,N) inverse matrix of A                           ==
C     ==         B is a work array                                    ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,INFO
      REAL*8  A(N,N),B(N,N)
C     Variables
      INTEGER LB
C     ==--------------------------------------------------------------==
      IF(N.EQ.1) THEN
        A(1,1)=1.D0/A(1,1)
      ELSE
        LB=(N-1)*N
C       Compute an LU factorization
        CALL DGETRF(N,N,A,N,B(1,1),INFO)
        IF(INFO.EQ.0) THEN
C         Compute the inverse of a matrix using the LU factorization
          CALL DGETRI(N,A,N,B(1,1),B(1,2),LB,INFO)
        ENDIF
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE INVERSEMAT(N,A,LDA,B,INFO)
C     ==--------------------------------------------------------------==
C     == Inverse Matrix A(LDA,N)                                      ==
C     ==--------------------------------------------------------------==
C     == INPUT:  N   Dimension                                          ==
C     ==         LDA Leading dimension of A                           ==
C     ==         A(LDA,N) Matrix                                      ==
C     == OUTPUT: A(LDA,N) inverse matrix of A                         ==
C     ==         B is a work array                                    ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,LDA,INFO
      REAL*8  A(LDA,N),B(N,N)
C     Variables
      INTEGER LB
C     ==--------------------------------------------------------------==
      IF(N.EQ.1) THEN
        A(1,1)=1.D0/A(1,1)
      ELSE
        LB=(N-1)*N
C       Compute an LU factorization
        CALL DGETRF(N,N,A,LDA,B(1,1),INFO)
        IF(INFO.EQ.0) THEN
C         Compute the inverse of a matrix using the LU factorization
          CALL DGETRI(N,A,LDA,B(1,1),B(1,2),LB,INFO)
        ENDIF
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE DSPEVY(IOPT,AP,W,Z,LDZ,N,AUX,NAUX)
C     ==--------------------------------------------------------------==
C     == DIAGONALIZATION ROUTINE: FOLLOW THE ESSL CONVENTION          ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER IOPT,LDZ,N,NAUX
      REAL*8    AP(N*(N+1)/2),W(N),Z(LDZ,N),AUX(NAUX)
C     Variables
      INTEGER INFO
C     ==--------------------------------------------------------------==
      CALL TEST_SCR('DSPEVY','NAUX',NAUX,3*N)
      INFO=0
      IF(IOPT.EQ.0) THEN
        CALL DSPEV('N','L',N,AP,W,Z,LDZ,AUX,INFO)
      ELSEIF(IOPT.EQ.1) THEN
        CALL DSPEV('V','L',N,AP,W,Z,LDZ,AUX,INFO)
      ELSEIF(IOPT.EQ.20) THEN
        CALL DSPEV('N','U',N,AP,W,Z,LDZ,AUX,INFO)
      ELSEIF(IOPT.EQ.21) THEN
        CALL DSPEV('V','U',N,AP,W,Z,LDZ,AUX,INFO)
      ENDIF
      IF(INFO.NE.0) CALL STOPGM('DSPEVY','FAILED TO DIAGONALIZE')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#if !(defined(ESSL))
C     ==================================================================
      SUBROUTINE ZGETMO(A,LDA,M,N,B,LDB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    LDA,M,N,LDB
      COMPLEX*16 A(LDA,N),B(LDB,M)
C     Variables
      INTEGER    I,J
C     ==--------------------------------------------------------------==
      DO I=1,M
        DO J=1,N
          B(J,I)=A(I,J)
        ENDDO
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE ZGTHR(N,A,B,IND)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    N,IND(N)
      COMPLEX*16 A(*),B(N)
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,N
        B(I)=A(IND(I))
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE ZSCTR(N,A,IND,B)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    N,IND(N)
      COMPLEX*16 A(N),B(*)
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,N
        B(IND(I))=A(I)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE DURAND(SEED,N,A)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      REAL*8  A(N),SEED,XRANF
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,N
        A(I)=XRANF()
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE DAXPYI(NZ,ALPHA,X,INDX,Y)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER NZ,INDX(NZ)
      REAL*8  ALPHA,X(NZ),Y(*)
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,NZ
        Y(INDX(I))=Y(INDX(I))+ALPHA*X(I)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#if   !(defined(__SGI))
      FUNCTION IZAMAX(N,A,INC)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER IZAMAX
C     Arguments
      INTEGER    N,INC
      COMPLEX*16 A(*)
C     Variables
      INTEGER    I,II
      REAL*8     XMAX,Y
C     ==--------------------------------------------------------------==
      IZAMAX=0
      XMAX=-1.0
      DO I=1,N
        II=(I-1)*INC+1
        Y=ABS(A(II))
        IF(Y.GT.XMAX) THEN
          XMAX=Y
          IZAMAX=I
        ENDIF
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
#endif
#if !(defined(CRAY) || defined(__IBM))
      FUNCTION IDAMIN(N,X,INX)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER IDAMIN
C     Arguments
      INTEGER N,INX
      REAL*8 X(N)
C     Variables
      INTEGER I,IJ,II
      REAL*8 XMIN
C     ==--------------------------------------------------------------==
      XMIN=X(1)
      IJ=1
      DO I=2,N
        II=(I-1)*INX+1
        IF(XMIN.GE.ABS(X(II))) IJ=I
      ENDDO
      IDAMIN=IJ
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
#if defined(T3D) || defined(T3E)
      SUBROUTINE MY_GETARG(M,STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER M
      CHARACTER STRING*(*)
C     Variables
      INTEGER   L,IERROR
C     ==--------------------------------------------------------------==
      CALL PXFGETARG(M,STRING,L,IERROR)
      IF(IERROR.NE.0) STRING=' '
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE MY_GETENV(EVAR,STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER EVAR*(*),STRING*(*)
C     Variables
      INTEGER   L,IERROR
C     ==--------------------------------------------------------------==
      CALL PXFGETENV(EVAR,15,STRING,L,IERROR)
      IF(IERROR.NE.0) STRING=' '
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE MY_GETLOG(STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER STRING*(*)
C     Variables
      INTEGER L,IERROR
C     ==--------------------------------------------------------------==
      CALL PXFGETLOGIN(STRING,L,IERROR)
      IF(IERROR.NE.0) STRING=' '
C     ==--------------------------------------------------------------==
      RETURN
      END
C     =================================================================
#elif defined(__SR2201)
      SUBROUTINE MY_GETARG(M,STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER M
      CHARACTER*(*) STRING
C     ==--------------------------------------------------------------==
      CALL GETARG(M+1,STRING)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     =================================================================
      SUBROUTINE MY_GETENV(EVAR,STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER LEN1,LEN2
      CHARACTER*(*) EVAR,STRING
C     ==--------------------------------------------------------------==
      LEN1=LEN(EVAR)
      LEN2=80
      CALL GETENV(EVAR,LEN1,STRING,LEN2)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     =================================================================
      SUBROUTINE MY_GETLOG(STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER*(*) STRING
C     ==--------------------------------------------------------------==
      CALL GETLOG(STRING)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     =================================================================
#else
      SUBROUTINE MY_GETARG(M,STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER M
      CHARACTER*(*) STRING
C     ==--------------------------------------------------------------==
      CALL GETARG(M,STRING)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     =================================================================
      SUBROUTINE MY_GETENV(EVAR,STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER*(*) EVAR,STRING
C     ==--------------------------------------------------------------==
      CALL GETENV(EVAR,STRING)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     =================================================================
      SUBROUTINE MY_GETLOG(STRING)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      CHARACTER*(*) STRING
C     ==--------------------------------------------------------------==
#if defined(__IBM)
      CALL GETLOG_(STRING)
#elif defined(_vpp_)
      CALL GETENV('USER',STRING)
#else
      CALL GETLOG(STRING)
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
      SUBROUTINE ZCLEAN_K(A,N,NGW)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    N,NGW
      COMPLEX*16 A(2*NGW,N)
C     Variables
      INTEGER    I
C     ==--------------------------------------------------------------==
      DO I=1,N
        A(NGW+1,I)=DCMPLX(0.D0,0.D0)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SETKWF(NGW,NSTATE,A)
C     ==--------------------------------------------------------------==
C     == A USES INVERSION SYMMETRY.                                   ==
C     == WE TRANSFORM A IN ARRAY WITHOUT INVERSION SYMMETRY.          ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'geq0.inc'
C     Arguments
      INTEGER    NGW,NSTATE
      COMPLEX*16 A(NGW,2*NSTATE)
C     Variables
      INTEGER    I,II,J
C     ==--------------------------------------------------------------==
      DO I=NSTATE,2,-1
        CALL DCOPY(2*NGW,A(1,I),1,A(1,2*I-1),1)
      ENDDO
      DO I=1,NSTATE
        II=2*I
        CALL DCOPY(2*NGW,A(1,II-1),1,A(1,II),1)
        DO J=1,NGW
          A(J,II)=DCONJG(A(J,II))
        ENDDO
      ENDDO
      IF(GEQ0) CALL ZCLEAN_K(A,NSTATE,NGW)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE FSKIP(NF,NREC)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
C     Arguments
      INTEGER NF,NREC
C     Variables
      INTEGER I
C     ==--------------------------------------------------------------==
      DO I=1,NREC
        READ(NF,END=20,ERR=30)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
C     ==--------------------------------------------------------------==
 20   CONTINUE
      CALL STOPGM('FSKIP','END OF FILE')
 30   CONTINUE
      CALL STOPGM('FSKIP','READ ERROR')
      END
C     ==================================================================
      FUNCTION DGIVE(A,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 DGIVE
C     Arguments
      INTEGER N
      REAL*8  A(N)
C     ==--------------------------------------------------------------==
      DGIVE=A(N)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      FUNCTION ZGIVE(A,N)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      COMPLEX*16 ZGIVE
C     Arguments
      INTEGER    N
      COMPLEX*16 A(N)
C     ==--------------------------------------------------------------==
      ZGIVE=A(N)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
C     == USED TO DEBUG.                                               ==
C     ==================================================================
      SUBROUTINE ROUNDOFF(N,A)
C     ==--------------------------------------------------------------==
C     == ROUNDOFF ERROR NEAR AN INTEGER VALUE                         ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER   N
      REAL*8    A(N)
C     Variables
      REAL*8    EPS
      PARAMETER (EPS=1.D-15)
      INTEGER   I
C     ==--------------------------------------------------------------==
      DO I=1,N
        A(I) = A(I)- DMOD(A(I),EPS)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE DEBROT(A,N,NSTATE)
C     ==--------------------------------------------------------------==
C     == G=0 COMPONENT OF C0 IS REAL                                  ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    N,NSTATE
      COMPLEX*16 A(N,NSTATE)
C     Variables
      REAL*8     TOLL
      PARAMETER  (TOLL=1.D-13)
      INTEGER    I,II
      COMPLEX*16 ZS
C     ==--------------------------------------------------------------==
      DO I=1,NSTATE
        DO II=1,N
          IF (ABS(A(II,I)).GT.TOLL) GOTO 10
        ENDDO
        II=N
 10     CONTINUE
        IF(ABS(DIMAG(A(II,I))).GT.TOLL) THEN
          ZS=ABS(A(II,I))/A(II,I)
        ELSE
          ZS=DCMPLX(1.D0,0.D0)
        ENDIF
        CALL ZSCAL(N,ZS,A(1,I),1)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
C      SUBROUTINE STORE(C0,NUM,IOPT,IUNIT)
C     ==--------------------------------------------------------------==
C     == STORE IN A FILE AN ARRAY TO DEBUG                            ==
C     ==--------------------------------------------------------------==
C      IMPLICIT NONE
C     Arguments
C      INTEGER NUM,IOPT,IUNIT
C      REAL*8  C0(2,*)
C     Variables
C      INTEGER I,J
C      REAL*8  DNRM2,DOTP
C     ==--------------------------------------------------------------==
C      OPEN(UNIT=IUNIT,STATUS= 'UNKNOWN')
C      IF(IOPT.EQ.1) THEN
C        WRITE(IUNIT,*) 0, SQRT(DOTP(NUM,C0(1,1),C0(1,1))),0.D0
C      ELSEIF(IOPT.EQ.2) THEN
C        WRITE(IUNIT,*) 0, DNRM2(2*IOPT*NUM,C0(1,1),1),0.D0
C      ENDIF
C      IF(IOPT.LE.2) THEN
C        DO I=1,NUM*IOPT
C          WRITE(IUNIT,*) I, C0(1,I), C0(2,I)
C        ENDDO
C      ELSE
C        DO I=1,NUM/2
C          DO J=1,2
C            WRITE(IUNIT,*) 2*(I-1)+J, C0(J,I)
C          ENDDO
C        ENDDO
C        IF(2*(NUM/2).NE.NUM) WRITE(IUNIT,*) NUM, C0(1,NUM/2+1)
C      ENDIF
C      CLOSE(UNIT=IUNIT)
CC     ==--------------------------------------------------------------==
C      RETURN
C      END
C     ==================================================================
C     == USED TO DEBUG.                                               ==
C     ==================================================================
      SUBROUTINE DETERM(DMAT,N,ISCR,LSCR,DD)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER    N,ISCR,LSCR
      COMPLEX*16 DMAT(N,N),DD
C     Variables
      INTEGER    I,INFO
C     ==--------------------------------------------------------------==
      IF(LSCR.LT.N) CALL STOPGM('DETERM','SCRATCH TOO SMALL')
      CALL ZGETRF(N,N,DMAT,N,ISCR,INFO)
      IF(INFO.NE.0) CALL STOPGM('DETERM','ILLEGAL RESULTS')
      DD=1.D0
      DO I=1,N
        DD=DD*DMAT(I,I)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      FUNCTION NXXFUN(NSTATE)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'system.h'
C     Arguments
      INTEGER NXXFUN,NSTATE
C     Variables
      INTEGER N1,N2,I,NS
C     ==--------------------------------------------------------------==
      N1=0
      N2=0
      DO I=0,NPROC-1
        NS=NST12(I,2)-NST12(I,1)+1
        N1=MAX(N1,NS)
        N2=MAX(N2,SPARM(3,I))
      ENDDO
      NXXFUN=N1*N2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WRITEVEC(CHAR,VEC,N,IUNIT)
      CHARACTER*7 CHAR
      REAL*8 VEC(N)
      OPEN(IUNIT,FILE=CHAR,STATUS='UNKNOWN')
      DO I=1,N
        WRITE(IUNIT,'(I6,F19.9)') I,VEC(I)
      ENDDO
      WRITE(6,'(A,A,A,I3)') ' FILE ', CHAR,' WRITTEN OUT'
      CLOSE(IUNIT)
      RETURN 
      END 
      SUBROUTINE READVEC(CHAR,VEC,N,IUNIT)
      CHARACTER*7 CHAR
      REAL*8 VEC(N)
      OPEN(IUNIT,FILE=CHAR,STATUS='OLD')
      DO I=1,N
        READ(IUNIT,*) II,VEC(I)
      ENDDO
      WRITE(6,'(A,A,A,I3)') ' FILE ', CHAR,' READ IN'
      CLOSE(IUNIT)
      RETURN 
      END 

C.. A fudge of the unportrable routine DGEMA to add two matrices
      SUBROUTINE DGEMA(TRANSA,TRANSB,M,N,ALPHA,A,LDA,BETA,B,LDB,C,LDC)
         IMPLICIT NONE
         CHARACTER*1 TRANSA,TRANSB
         INTEGER*4 M,N,LDA,LDB,LDC
         REAL*8 ALPHA,A(M,N),B(M,N),BETA,C(M,N)
         INTEGER I,J
         IF(LDA.NE.M) STOP 'LDA NE M not supported'
         IF(LDB.NE.M) STOP 'LDB NE M not supported'
         IF(LDC.NE.M) STOP 'LDC NE M not supported'
         IF(TRANSA.NE.'N'.AND.TRANSA.NE.'n')
     &      STOP 'TRANSA not "N" not supported'
         IF(TRANSB.NE.'N'.AND.TRANSB.NE.'n')
     &      STOP 'TRANSB not "N" not supported'
         DO I=1,M
            DO J=1,N
               C(I,J)=ALPHA*A(I,J)+BETA*B(I,J)
            ENDDO
         ENDDO
         RETURN
      END


C.. RAN2 - generate a decent random number from Numerical Recipes
C.. 15th Aug '05.  It would appear this is not Numerical Recipes'
C.. RAN2, but an implementation of something different with the same
C.. interface.

C.. Initialize with a seed <0

C.. Modified this to use REAL*8s.
      REAL*8 FUNCTION ORAN2(IDUM)
      IMPLICIT NONE
      PARAMETER (M=714025,IA=1366,IC=150889)
      INTEGER IR(97),IFF,J,M,IY,IA,IC,IDUM
      REAL*8 RM
      DATA IFF /0/
      RM=1.4005112D-6
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      ORAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END

C.. Numerical Recipes RAN2

      REAL*8 FUNCTION OORAN2(idum)
        IMPLICIT NONE
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL*8 AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.D0/IM1,
     &          IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,
     &          IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,
     &          EPS=1.2D-7,RNMX=1.D0-EPS)
        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
!        INTEGER C
        DATA idum2/123456789/, iv/NTAB*0/,iy/0/
!        SAVE C
!        C=C+1
!        IF(MOD(C,10000).EQ.0) WRITE(6,*) "PRNG: ",C
        if (idum.le.0) then
           idum=max(-idum,1)
           idum2=idum
           do j=NTAB+8,1,-1
              k=idum/IQ1
              idum=IA1*(idum-k*IQ1)-k*IR1
              if (idum.lt.0) idum=idum+IM1
              if (j.le.NTAB) iv(j)=idum
           enddo
           iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1) iy=iy+IMM1
        ooran2=min(AM*iy,RNMX)
        return
      END

!  Following the same interface as RAN2, but using the RANLUX generator (lux level 3).
!  If ISEED.NE.0, it is used as a seed, and reset to 0.
!  Generates a 32-bit float (REAL*4), and stores this in a real*8
      REAL*8 FUNCTION RAN2(ISEED)
         IMPLICIT NONE
         INTEGER ISEED
         REAL*8 r
         IF(ISEED.NE.0) THEN
!  Init RANLUX
            IF(ISEED.LT.0) ISEED=-ISEED
            CALL RLUXGO(3,ISEED,0,0)
            ISEED=0
         ENDIF
         CALL RANLUX(r,1)
         RAN2=r
         RETURN
      END
