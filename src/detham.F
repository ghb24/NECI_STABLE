      SUBROUTINE DETHAM(NDET,NEL,NMRKS,NBASISMAX,
     &	NHG,HAMIL,G1,LAB,NROW,TCOUNT,NMSH,FCK,
     &	NMAX,ALAT,UMAT,ICMAX,GC,TMC,ECORE,BRR)
      use global_utilities
      USE HElem
      Use Determinants, only: get_helement
      use SystemData , only : BasisFN
      use CalcData , only : TStar
      IMPLICIT NONE
C      IMPLICIT REAL*8 (A-H,O-Z)
      TYPE(HElement) HAMIL(*)
      INTEGER NDET,NEL,NMSH,NMAX
      INTEGER LAB(*)
      TYPE(BasisFn) G1(*)
      INTEGER NMRKS(NEL,*)
      INTEGER NROW(NDET),GC
C..Cube Arrays
      COMPLEX*16 FCK(NMSH,NMSH,NMSH)
      REAL*8 ALAT(3),ECORE
      INTEGER nBasisMax(5,*),NHG,BRR(*)
      TYPE(HElement) UMat(*)
C..
      INTEGER ICMAX,KI,IBEG,IBEGJ,KJ,IC,I,J,IMAX,IDAMAX
      LOGICAL TCOUNT,TMC
      INTEGER IGETEXCITLEVEL_2,STEP,IMAXJ
      TYPE(HElement) totSUM
      type(timer), save :: proc_timer
!      LOGICAL TSTAR
C ==-------------------------------------------------------------------==
      proc_timer%timer_name='    DETHAM'
      call set_timer(proc_timer)
C ==-------------------------------------------------------------------==
C..Global counter
      GC=0
      NROW(1:NDET)=0
C..   Now we need to match up any two determinants
      DO KI=1,NDET
         IF(mod(KI,1000).eq.0) WRITE(6,*) KI
         IF(KI.EQ.1) THEN
           IBEG=0
         ELSE
           IBEG=IBEG+NROW(KI-1)
         ENDIF
         IF(TMC) THEN 
           IBEGJ=1
           STEP=1
           IMAXJ=NDET
	 ELSE
           IF(TSTAR) THEN
               IF(KI.EQ.1) THEN
                  IBEGJ=1
                  IMAXJ=NDET
                  STEP=1
               ELSE
                  IBEGJ=KI
                  IMAXJ=KI
                  STEP=1
               ENDIF
           ELSE
               IBEGJ=KI
               IMAXJ=NDET
               STEP=1
           ENDIF
	 ENDIF
         IF(STEP.EQ.0) STEP=1
         DO KJ=IBEGJ,IMAXJ,STEP
C..   IC is the number of pairs which are the *SAME* in both determinants
            IC=0
            IF(KI.EQ.KJ) THEN
              IC=NEL
            ELSE
C..              IC=NEL-IGETEXCITLEVEL(NMRKS(1,KI),NMRKS(1,KJ),NEL)
              IC=NEL-IGETEXCITLEVEL_2(NMRKS(1,KI),NMRKS(1,KJ),NEL,2)
C..   If determinants differ by more than 2 spin orbitals IC < 2
              IF(IC.LT.(NEL-2)) GOTO 500
            ENDIF
C..
!            CALL WRITEDET(6,NMRKS(1,KI),NEL,.FALSE.)
!            CALL WRITEDET(6,NMRKS(1,KJ),NEL,.TRUE.)
            totsum = get_helement(NMRKS(:,KI), NMRKS(:,KJ),nel-ic)
C            CALL SLTCND(NEL,NBASISMAX,NHG,
C     &		NMRKS(1:NEL,KI),NMRKS(1:NEL,KJ),G1,IC,NMSH,
C     &		FCK,NMAX,ALAT,UMAT,TOTSUM)
            IF(ABS(TOTSUM%v).LT.1.D-10) TOTSUM%v=0.D0
            IF(TOTSUM%v.NE.0.D0.OR.KI.EQ.KJ) THEN
               GC=GC+1
C..   Stores the number of non-zero elements in each row
               NROW(KI)=NROW(KI)+1
               IF(.NOT.TCOUNT) THEN
                 LAB(IBEG+NROW(KI))=KJ
                 HAMIL(IBEG+NROW(KI))=TOTSUM
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
      call halt_timer(proc_timer)
C ==-------------------------------------------------------------------==
      RETURN
      END
C ==---------------------------------------------------------------==
      SUBROUTINE CHK_MS(ND1,ND2,NBASISMAX,NEL,IFLG,G1)
      use SystemData, only: BasisFN
      IMPLICIT REAL*8 (A-H,O-Z)
C..
      DIMENSION ND1(NEL),ND2(NEL)
      TYPE(BASISFN) G1(*)
C..
      MS1=0
      MS2=0
      DO I=1,NEL
        MS1=MS1+G1(ND1(I))%MS
        MS2=MS2+G1(ND2(I))%MS
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