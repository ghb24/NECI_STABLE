      FUNCTION ERF(X)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This has been lifted from numerical recipes
      IF(X.LT.0.)THEN
        ERF=-GAMMP(0.5D0,X**2)
      ELSE 
        ERF=GAMMP(0.5D0,X**2)
      ENDIF
C..
      RETURN
      END
C ==-------------------------------------------------------------==
      FUNCTION GAMMP(A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
C..
      IF(X.LT.0.D0.OR.A.LE.0.D0) STOP
      IF(X.LT.A+1.D0)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.D0-GAMMCF
      ENDIF
C..
      RETURN
      END
C ==-------------------------------------------------------------==
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.D0)THEN
        IF(X.LT.0.D0)STOP
        GAMSER=0.D0
        RETURN
      ENDIF
      AP=A
      SUM=1.D0/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.D0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11      CONTINUE
        STOP 'A too large, ITMAX too small'
1       GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
C..
      RETURN
      END
C ==-------------------------------------------------------------==
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      IMPLICIT REAL*8 (A-H,O-Z)
C..
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.D0
      A0=1.D0
      A1=X
      B0=0.D0
      B1=1.D0
      FAC=1.D0
      DO 11 N=1,ITMAX
	AN=FLOAT(N)
	ANA=AN-A
	A0=(A1+A0*ANA)*FAC
	B0=(B1+B0*ANA)*FAC
	ANF=AN*FAC
	A1=X*A0+ANF*A1
	B1=X*B0+ANF*B1
	IF(A1.NE.0.D0)THEN
	  FAC=1./A1
	  G=B1*FAC
	  IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
	  GOLD=G
	ENDIF
11    CONTINUE
      STOP 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
C..
      RETURN
      END
