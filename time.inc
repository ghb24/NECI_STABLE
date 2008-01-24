C     ==================================================================
C     == INCLUDE FILE FOR TIMING                                      ==
C     ==================================================================
      INTEGER        NTIMX
      PARAMETER      (NTIMX=200)
C     ==--------------------------------------------------------------==
      CHARACTER*10   TNAM
      COMMON /TNAMENECI/ TNAM(NTIMX)
C     ==================================================================
      REAL*8         CTSUM(NTIMX),CTSTA(NTIMX),WTSUM(NTIMX),WTSTA(NTIMX)
      COMMON /TIMBUNECI/ CTSUM,CTSTA,WTSUM,WTSTA
C     ==================================================================
      INTEGER        NTIMER,NCALLS(NTIMX)
! now in LOGREAD module      INTEGER        iGlobalTimerLevel
      COMMON /TIMINNECI/ NTIMER,NCALLS!,iGlobalTimerLevel
C     ==================================================================
