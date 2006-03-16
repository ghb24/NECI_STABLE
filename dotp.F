C     ==================================================================
      FUNCTION DOTP(N,A,B)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'geq0.inc'
C     Arguments
      INTEGER N
      REAL*8 DOTP,A(2*N),B(2*N)
C     Variables
#if defined (_vpp_)
      INTEGER II
      REAL*8  DDO
#else
      REAL*8 DDOT
#endif
C     ==--------------------------------------------------------------==
      IF(N.EQ.0) THEN
        DOTP=0.D0
        RETURN
      ENDIF
      IF(GEQ0) THEN
        DOTP=A(1)*B(1)
      ELSE
        DOTP=2.0D0*(A(1)*B(1)+A(2)*B(2))
      ENDIF
      IF(N.GT.1) THEN
#if defined (_vpp_)
        DDO=0.D0
        DO II=3,2*N
          DDO=DDO+A(II)*B(II)
        ENDDO
        DOTP=DOTP+DDO+DDO
#else
        DOTP=DOTP+2.0D0*DDOT(2*N-2,A(3),1,B(3),1)
#endif
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
