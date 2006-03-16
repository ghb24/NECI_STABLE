#ifndef __HP
C     ==================================================================
      FUNCTION TIMEC()
C     ==--------------------------------------------------------------==
C     ==  GENERALL CPU TIME ROUTINE                                   ==
C     ==  CPU TIME IS RETURNED IN msec                                ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 TIMEC
#if   defined(YMP)
      REAL*8 SECOND
C     ==--------------------------------------------------------------==
      TIMEC = SECOND()*1000.D0
#elif defined(T3D)
      INTEGER IRTC
C     DCLOCK is the cpu clock.
      REAL*8 DCLOCK
      PARAMETER (DCLOCK=150.D0)
C     ==--------------------------------------------------------------==
      TIMEC = DBLE(IRTC()) * (0.001D0/DCLOCK)
#elif defined(T3E)
C     DCLOCK is the cpu clock.
      REAL TSECND
      EXTERNAL TSECND
C     ==--------------------------------------------------------------==
      TIMEC = TSECND() * 1.0E3
#elif defined(__NEC)
      REAL*8 D
C     ==--------------------------------------------------------------==
      CALL CLOCK(D)
      TIMEC=D*1000.D0
#elif defined(__alpha) || defined(__SUN)
      REAL*4 ETIME,TARRAY(2),D
      EXTERNAL ETIME
C     ==--------------------------------------------------------------==
      D=ETIME(TARRAY)
      TIMEC = DBLE(TARRAY(1))*1000.D0
#elif defined(__SGI)
      INTEGER MCLOCK
C     ==--------------------------------------------------------------==
      TIMEC = DBLE(MCLOCK())*10.D0
#elif defined(__SR2201)
      REAL*8 D
C     ==--------------------------------------------------------------==
      CALL XCLOCK(D,8)
      TIMEC = D*1000.D0
#elif defined(_vpp_)
      REAL*8 D
C     ==--------------------------------------------------------------==
      CALL CLOCK(D,1,2)
      TIMEC=D
#else 
      INTEGER MCLOCK
C     ==--------------------------------------------------------------==
      TIMEC = DBLE(MCLOCK())*10.D0
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#if ! defined(__IBM) && ! defined(CRAY)
      FUNCTION TIMEF()
C     ==--------------------------------------------------------------==
C     ==  GENERAL TIME ROUTINE FOR ELAPSED TIME                       ==
C     ==  ELAPSED TIME IS RETURNED IN msec                            ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      REAL*8 TIMEF
#if defined(__SGI) || defined(_SUN) || defined(__alpha) || defined(_vpp_)
      REAL*4 TARRAY(2),ETIME
C     ==--------------------------------------------------------------==
      TIMEF=ETIME(TARRAY)*1000.D0
#elif defined(__NEC)
      REAL*8 D
C     ==--------------------------------------------------------------==
      CALL ETIME(D)
      TIMEF=D*1000.D0
#elif defined(__SR2201)
      REAL*8 D
C     ==--------------------------------------------------------------==
      CALL XCLOCK(D,8)
      TIMEF = D*1000.D0
#elif defined(__PGI)
      INTEGER MCLOCK
      TIMEF = DBLE(MCLOCK())*10.D0
#else
      REAL*8 TIME
      EXTERNAL TIME
C     ==--------------------------------------------------------------==
      TIMEF=TIME()
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
#endif
