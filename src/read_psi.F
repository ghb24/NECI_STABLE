      SUBROUTINE READ_PSI(BOX,BOA,COA,NDET,NEVAL,
     &        NBASISMAX,NEL,CK,W)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER nBasisMax(5,*)
      integer :: NEVAL,NDET,I,J,NEL,NDETCK,NEVALCK,NELCK
      real(dp) :: CK(NDET,NEVAL)
      real(dp) :: W(NEVAL)
      integer :: NMAXCK,NMAYCK,NMAZCK
      real(dp) :: BOX,BOA,COA,C,BOXCK,BOACK,COACK
      character(*), parameter :: this_routine = 'read_psi'
C..
      OPEN(10,FILE='PSI_LONG',STATUS='OLD',ERR=10)
      READ(10,*) NMAXCK,NMAYCK,NMAZCK,NELCK,NEVALCK,NDETCK
      READ(10,*) BOXCK,BOACK,COACK
      IF(    NMAXCK.NE.NBASISMAX(1,2)
     &   .OR.NMAYCK.NE.NBASISMAX(2,2)
     &   .OR.NMAZCK.NE.NBASISMAX(3,2)) 
     &     call stop_all(this_routine, ' !!! DIFFERENT SIZE BASIS !!! ')
      IF(NELCK.NE.NEL) call stop_all(this_routine, 
     &                        ' !!! DIFFERENT NUMBER OF ELECTRONS !!! ')
      IF(NEVALCK.NE.NEVAL) call stop_all(this_routine,
     &                                      ' !!! NEVAL DIFFERENT !!! ')
      IF(NDETCK.NE.NDET) call stop_all(this_routine, 
     &                             ' !!! DIFFERENT NUMBER OF DETS !!! ')
      IF(BOXCK.NE.BOX.OR.BOACK.NE.BOA.OR.COACK.NE.COA) 
     &    call stop_all(this_routine, ' !!! DIFFERENT SIZE CUBES !!! ' )
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
     &     call stop_all(this_routine, ' !!! DIFFERENT SIZE BASIS !!! ')
      IF(NELCK.NE.NEL) call stop_all(this_routine,
     &                        ' !!! DIFFERENT NUMBER OF ELECTRONS !!! ')
      IF(NEVALCK.NE.NEVAL) call stop_all(this_routine,
     &                                      ' !!! NEVAL DIFFERENT !!! ')
      IF(NDETCK.NE.NDET) call stop_all(this_routine,
     &                             ' !!! DIFFERENT NUMBER OF DETS !!! ')
      IF(BOXCK.NE.BOX.OR.BOACK.NE.BOA.OR.COACK.NE.COA) 
     &    call stop_all(this_routine, ' !!! DIFFERENT SIZE CUBES !!! ' )
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
     &    NBASISMAX,NEL,CK,W)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER nBasisMax(5,*)
      INTEGER I,NEVAL,NDET,NEL,J
      HElement_t(dp) CK(NDET,NEVAL)
      real(dp) W(NEVAL),BOA,BOX,COA
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
     &    NBASISMAX,NEL,CK,W)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER I,NEVAL,NDET,NEL,J
      INTEGER nBasisMax(5,*)
      HElement_t(dp) CK(NDET,NEVAL)
      real(dp) W(NEVAL),BOA,BOX,COA
C..
      OPEN(10,FILE='PSI_COMP',STATUS='UNKNOWN')
      WRITE(10,*) NBASISMAX(1,2),NBASISMAX(2,2),NBASISMAX(3,2),
     &      NEL,NEVAL,NDET 
      WRITE(10,*) BOX,BOA,COA
      DO I=1,NEVAL
        DO J=1,NDET
          IF(abs(CK(J,I)).gt. 1.0e-2_dp) THEN
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
