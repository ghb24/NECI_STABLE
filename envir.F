C     ==================================================================
      SUBROUTINE ENVIR
C     ==--------------------------------------------------------------==
C     ==  LOOK FOR THE ENVIRONMENT THIS JOB IS RUNNING                ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'envj.inc'
      CHARACTER FFORMAT*9
      INTEGER I1,I2
#if defined(CRAY) && !defined(CRAYX1)
      INTEGER ISHOST,ISDIR
      INTEGER HOSTNM,GETCWD,GETPID,GETUID
#elif defined(CRAYX1)
      INTEGER*8 ISHOST
      INTEGER*8 pathlen,IERROR
#elif defined(__SGI)
      INTEGER ISHOST,ISDIR
      INTEGER GETHOSTNAME
      INTEGER GETPID,GETUID
#elif defined(__SR2201) || defined(__SR8000)
      INTEGER ISHOST,ISDIR
      INTEGER GETPID,GETUID,GETCWD,HOSTNM
#elif defined(__alpha)
      INTEGER ISDIR,ISHOST
      INTEGER GETCWD,GETPID,GETUID,HOSTNM
#elif defined(__SUN) || defined(__Linux) || defined(__WINNT)
      INTEGER ISDIR
      INTEGER GETCWD,GETPID,GETUID
#elif defined(_vpp_) || defined(__NEC)
      INTEGER ISHOST,ISDIR
      INTEGER GETPID,GETUID,GETCWD,HOSTNM
#else
      INTEGER ISHOST,ISDIR
      INTEGER HOSTNM_,GETCWD_,GETPID_,GETUID_
#endif
#if defined(JOBLIMIT)
      INTEGER GETJOBLIMIT
#endif
C     ==--------------------------------------------------------------==
C     ==  LOOK FOR THE ENVIRONMENT THIS JOB IS RUNNING                ==
C     ==--------------------------------------------------------------==
      WRITE(HNAME,'(64(" "))')
#ifdef CRAYX1
      HNAME = ' '
      CALL PXFGETENV('HOST',0,HNAME,0,IERROR)
      CURDIR = ' '
      CALL PXFGETCWD(CURDIR,pathlen,IERROR)
      CALL PXFGETPID(MY_PID,IERROR)
      CALL PXFGETUID(MY_UID,IERROR)
      USER = ' '
      CALL MY_GETLOG(USER)
#elif defined(CRAY)
      HNAME = ' '
      ISHOST = HOSTNM(HNAME)
      CURDIR = ' '
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
      USER = ' '
      CALL MY_GETLOG(USER)
C..this causes many problems for many implementations
C     TJLIMIT=GETJOBLIMIT()
#elif defined(__SGI)
      ISHOST = GETHOSTNAME(HNAME,32)
      CALL GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
      CALL MY_GETLOG(USER)
#elif defined(__SR2201) || defined(__SR8000)
      ISHOST = HOSTNM(HNAME)
      CURDIR = ' '
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
      USER = ' '
      CALL MY_GETLOG(USER)
#elif defined(__alpha)
      ISHOST = HOSTNM(HNAME)
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
      CALL MY_GETLOG(USER)
#elif defined(__SUN) || defined(__Linux)
      CALL GETENV('HOSTNAME',HNAME)
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
      CALL MY_GETLOG(USER)
#elif defined(__WINNT)
      CALL GETENV('HOSTNAME',HNAME)
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = 0
      CALL MY_GETLOG(USER)
#elif defined(__PGI) 
      CALL GETENV('HOSTNAME',HNAME)
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
      CALL MY_GETLOG(USER)
 
#elif defined(_vpp_) || defined(__NEC)
      ISHOST = HOSTNM(HNAME)
      ISDIR = GETCWD(CURDIR)
      MY_PID = GETPID()
      MY_UID = GETUID()
#else
      ISHOST = HOSTNM_(HNAME)
      ISDIR = GETCWD_(CURDIR)
      MY_PID = GETPID_()
      MY_UID = GETUID_()
      CALL MY_GETLOG(USER)
#endif
#if defined(JOBLIMIT)
      TJLIMIT=GETJOBLIMIT()
      IF(TJLIMIT.GT.8640000D0) TJLIMIT=0.D0
#else
      TJLIMIT=0.D0
#endif
      CALL MY_GETENV('TMPDIR',TMPDIR)
      IF(TMPDIR(1:1).EQ.' ') TMPDIR=CURDIR

      CALL XSTRING(HNAME,I1,I2)
      WRITE(FFORMAT,'(A,I2,A)') '(A,T',MAX(20,65-(I2-I1)),',A)'
      WRITE(*,FFORMAT)   ' THIS JOB RUNS ON: ',HNAME(I1:I2)
      WRITE(*,'(A)')   ' THE CURRENT DIRECTORY IS: '
      CALL XSTRING(CURDIR,I1,I2)
      WRITE(FFORMAT,'(A,I2,A)') '(T',MAX(2,65-(I2-I1)),',A)'
      WRITE(*,FFORMAT)   CURDIR(I1:I2)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
