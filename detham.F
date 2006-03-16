      SUBROUTINE DETHAM(NDET,NEL,NMRKS,NBASISMAX,
     &	NHG,HAMIL,G1,LAB,NROW,TCOUNT,NMSH,FCK,ZIA,
     &	NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BRR)
      IMPLICIT NONE
C      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 HAMIL(*)
      INTEGER NDET,NEL,NMSH,NMAX
      INTEGER G1(4,*),LAB(*)
      INTEGER NMRKS(NEL,*)
      INTEGER NROW(NDET),GC
C..Cube Arrays
      COMPLEX*16 FCK(NMSH,NMSH,NMSH)
      COMPLEX*16 ZIA(-NMSH/2:NMSH/2,NMAX,NMAX)
      REAL*8 ALAT(3),ECORE
      INTEGER NBASISMAX(*),NHG,BRR(*)
      REAL*8 UMAT(*)
C..
      INTEGER ICMAX,ISUB,KI,IBEG,IBEGJ,KJ,IC,I,J,IMAX,IDAMAX
      LOGICAL TCOUNT,TMC
      REAL*8 GETHELEMENT2
      REAL*8 SUM
      INTEGER IGETEXCITLEVEL
C ==-------------------------------------------------------------------==
      CALL TISET('    DETHAM',ISUB)
C ==-------------------------------------------------------------------==
C..Global counter
      GC=0
      CALL IAZZERO(NROW,NDET)
C..   Now we need to match up any two determinants
      DO KI=1,NDET
         IF(KI.EQ.1) THEN
           IBEG=0
         ELSE
           IBEG=IBEG+NROW(KI-1)
         ENDIF
         IF(TMC) THEN 
           IBEGJ=1
	 ELSE
	   IBEGJ=KI
	 ENDIF
         DO KJ=IBEGJ,NDET
C..   IC is the number of pairs which are the *SAME* in both determinants
            IC=0
            IF(KI.EQ.KJ) THEN
              IC=NEL
            ELSE
              IC=NEL-IGETEXCITLEVEL(NMRKS(1,KI),NMRKS(1,KJ),NEL)
C..   If determinants differ by more than 2 spin orbitals IC < 2
              IF(IC.LT.(NEL-2)) GOTO 500
            ENDIF
C..
            SUM=GETHELEMENT2(NMRKS(1:NEL,KI),NMRKS(1:NEL,KJ),NEL,
     &         NBASISMAX,G1,NHG,BRR,NMSH,FCK,ZIA,NMAX,
     &         ALAT,UMAT,NEL-IC,ECORE)
C            CALL SLTCND(NEL,NBASISMAX,NHG,
C     &		NMRKS(1:NEL,KI),NMRKS(1:NEL,KJ),G1,IC,NMSH,
C     &		FCK,ZIA,NMAX,ALAT,UMAT,SUM)
            IF(ABS(SUM).LT.1.D-10) SUM=0.D0
            IF(ABS(SUM).GT.0.D0.OR.KI.EQ.KJ) THEN
               GC=GC+1
C..   Stores the number of non-zero elements in each row
               NROW(KI)=NROW(KI)+1
               IF(.NOT.TCOUNT) THEN
                 LAB(IBEG+NROW(KI))=KJ
                 HAMIL(IBEG+NROW(KI))=SUM
               ENDIF
            ENDIF
 500        CONTINUE
         ENDDO
      ENDDO
C..No. of columns
      IF(TCOUNT) THEN
        IMAX=IDAMAX(NDET,DFLOAT(NROW),1)
        ICMAX=NROW(IMAX)
        WRITE(6,*) ' MAXIMUM WIDTH OF HAMIL : ' , ICMAX
        WRITE(6,*) ' TOTAL NUMBER OF NON-ZERO ELEMENTS : ' , GC
      ENDIF
C      IF(.NOT.TCOUNT) THEN
C        WRITE(6,*)
C        WRITE(6,*) ' SECTION OF MAP OF HAMILTONIAN : '
C        WRITE(6,*)
C        WRITE(6,*) NDET,ICMAX
C        DO I=1,NDET
C          WRITE(6,*) (LAB(J,I),J=1,ICMAX)
C          WRITE(6,*) I,(HAMIL(J,I),J=1,ICMAX)
C        ENDDO
C      ENDIF
C ==-------------------------------------------------------------------==
      CALL TIHALT('    DETHAM',ISUB)
C ==-------------------------------------------------------------------==
      RETURN
      END
C ==---------------------------------------------------------------==
      SUBROUTINE CHK_MS(ND1,ND2,NBASISMAX,NEL,IFLG,G1)
      IMPLICIT REAL*8 (A-H,O-Z)
C..
      DIMENSION ND1(NEL),ND2(NEL)
      INTEGER G1(4,*)
C..
      MS1=0
      MS2=0
      DO I=1,NEL
        MS1=MS1+G1(4,ND1(I))
        MS2=MS2+G1(4,ND2(I))
      ENDDO
      MS1=MS1/2
      MS2=MS2/2
C..
      IF(MS1.NE.MS2) THEN
        IF(IFLG.EQ.1) STOP ' PROBLEM IN TWO DIFFER '
        IF(IFLG.EQ.2) STOP ' PROBLEM IN ONE DIFFER '
        IF(IFLG.EQ.3) STOP ' PROBLEM IN SAME '
      ENDIF 
C..
      RETURN
      END
