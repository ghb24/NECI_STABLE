      SUBROUTINE READ_PSI(BOX,BOA,COA,NDET,NEVAL,
     &		NBASISMAX,NEL,CK,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER nBasisMax(5,*)
C..
      DIMENSION CK(NDET,NEVAL)
      DIMENSION W(NEVAL)
C..
      OPEN(10,FILE='PSI_LONG',STATUS='OLD',ERR=10)
      READ(10,*) NMAXCK,NMAYCK,NMAZCK,NELCK,NEVALCK,NDETCK
      READ(10,*) BOXCK,BOACK,COACK
      IF(    NMAXCK.NE.NBASISMAX(1,2)
     &   .OR.NMAYCK.NE.NBASISMAX(2,2)
     &   .OR.NMAZCK.NE.NBASISMAX(3,2)) 
     &		       STOP ' !!! DIFFERENT SIZE BASIS !!! '
      IF(NELCK.NE.NEL) STOP ' !!! DIFFERENT NUMBER OF ELECTRONS !!! '
      IF(NEVALCK.NE.NEVAL) STOP ' !!! NEVAL DIFFERENT !!! '
      IF(NDETCK.NE.NDET) STOP ' !!! DIFFERENT NUMBER OF DETS !!! '
      IF(BOXCK.NE.BOX.OR.BOACK.NE.BOA.OR.COACK.NE.COA) 
     &		       STOP ' !!! DIFFERENT SIZE CUBES !!! ' 
C..
      DO I=1,NEVAL
        DO J=1,NDET
          READ(10,*) CK(J,I)
        ENDDO
      ENDDO
      DO I=1,NEVAL
         READ(10,*) W(I)
      ENDDO
      CLOSE(10)
C..
      WRITE(6,'(A)') ' HAVE READ IN PSI FROM FILE ./PSI_LONG '
C..
      RETURN
10    CONTINUE
C.. PSI_LONG doesn't exist.  what about PSI_COMP
      OPEN(10,FILE='PSI_COMP',STATUS='OLD')
      READ(10,*) NMAXCK,NMAYCK,NMAZCK,NELCK,NEVALCK,NDETCK
      READ(10,*) BOXCK,BOACK,COACK
      IF(    NMAXCK.NE.NBASISMAX(1,2)
     &   .OR.NMAYCK.NE.NBASISMAX(2,2)
     &   .OR.NMAZCK.NE.NBASISMAX(3,2)) 
     &		       STOP ' !!! DIFFERENT SIZE BASIS !!! '
      IF(NELCK.NE.NEL) STOP ' !!! DIFFERENT NUMBER OF ELECTRONS !!! '
      IF(NEVALCK.NE.NEVAL) STOP ' !!! NEVAL DIFFERENT !!! '
      IF(NDETCK.NE.NDET) STOP ' !!! DIFFERENT NUMBER OF DETS !!! '
      IF(BOXCK.NE.BOX.OR.BOACK.NE.BOA.OR.COACK.NE.COA) 
     &		       STOP ' !!! DIFFERENT SIZE CUBES !!! ' 
C..
      I=1
      DO WHILE(I.NE.0)
        READ(10,*) I,J,C
        IF(I.NE.0) CK(J,I)=C
      ENDDO
      DO I=1,NEVAL
         READ(10,*) W(I)
      ENDDO
      CLOSE(10)
      
      END
C ==----------------------------------------------------------------==   
      SUBROUTINE WRITE_PSI(BOX,BOA,COA,NDET,NEVAL,
     &	NBASISMAX,NEL,CK,W)
      USE HElem
      IMPLICIT NONE
      INTEGER nBasisMax(5,*)
C..
      TYPE(HElement) CK(NDET,NEVAL)
      REAL*8 W(NEVAL),BOA,BOX,COA
      INTEGER I,K,NEVAL,NDET,NEL,J
C..
      OPEN(10,FILE='PSI_LONG',STATUS='UNKNOWN')
      WRITE(10,*) NBASISMAX(1,2),NBASISMAX(2,2),NBASISMAX(3,2),
     &      NEL,NEVAL,NDET 
      WRITE(10,*) BOX,BOA,COA
      DO I=1,NEVAL
        DO J=1,NDET
          WRITE(10,*) CK(J,I),ABS(CK(J,I))
        ENDDO
      ENDDO
      DO I=1,NEVAL
         WRITE(10,*) W(I)
      ENDDO
      CLOSE(10)
C..
      RETURN
      END   

      SUBROUTINE WRITE_PSI_COMP(BOX,BOA,COA,NDET,NEVAL,
     &	NBASISMAX,NEL,CK,W)
      USE HElem
      IMPLICIT NONE
      INTEGER nBasisMax(5,*)
C..
      TYPE(HElement) CK(NDET,NEVAL)
      REAL*8 W(NEVAL),BOA,BOX,COA
      INTEGER I,K,NEVAL,NDET,NEL,J
C..
      OPEN(10,FILE='PSI_COMP',STATUS='UNKNOWN')
      WRITE(10,*) NBASISMAX(1,2),NBASISMAX(2,2),NBASISMAX(3,2),
     &      NEL,NEVAL,NDET 
      WRITE(10,*) BOX,BOA,COA
      DO I=1,NEVAL
        DO J=1,NDET
          IF(CK(J,I).AGT. 1.D-2) THEN
            WRITE(10,"(2I5)",advance='no') I,J
            WRITE(10,*) CK(J,I),ABS(CK(J,I))
         ENDIF
        ENDDO
      ENDDO
      WRITE(10,*) 0,0,0
      DO I=1,NEVAL
         WRITE(10,*) W(I)
      ENDDO
      CLOSE(10)
C..
      RETURN
      END   
