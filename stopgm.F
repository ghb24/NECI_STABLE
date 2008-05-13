C     ==================================================================
      SUBROUTINE STOPGM(A,B)
C     ==--------------------------------------------------------------==
#ifdef PARALLEL
         Use Parallel, only : me => iProcIndex, StopAll =>MPIStopAll
#endif
      IMPLICIT NONE
C     Arguments
      CHARACTER*(*) A,B
C     Variables
#ifdef PARALLEL
      INTEGER       NC
#endif
C     ==--------------------------------------------------------------==
C     CALL END_SWAP
#ifdef PARALLEL
      WRITE(*,'(//,A,A,A,A,A,I4,A)')
     &     ' PROGRAM STOPS IN SUBROUTINE ',A,'| ',B,' [PROC=',ME,']'
      NC=999
      CALL STOPALL(NC)
#else
      WRITE(*,'(//,A,A,A,A)') ' PROGRAM STOPS IN SUBROUTINE ',A,'| ',B
      STOP 999
#endif
C     ==--------------------------------------------------------------==
      END
C     ==================================================================
      SUBROUTINE TEST_SCR(A,TAG,LSCR,LENGTH)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER*(*) A,TAG
      INTEGER       LSCR,LENGTH
C     ==--------------------------------------------------------------==
      IF(LSCR.LT.LENGTH) THEN
        WRITE(*,'("LSCR=",I9," ",A,"=",I9)') LSCR,TAG,LENGTH
        CALL STOPGM(A,'SCRATCH ARRAY TOO SMALL')
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      
