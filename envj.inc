C     ==================================================================
C     == ENVIRONMENT VARIABLES                                        ==
C     ==================================================================
C     == USER   User name                                             ==
C     == HNAME  Computer name                                         ==
C     == CURDIR Current Directory                                     ==
C     == TMPDIR Temporary Directory                                   ==
C     == MY_PID Process ID Number                                     ==
C     == MY_UID User ID Number                                        ==
C     ==--------------------------------------------------------------==
      CHARACTER    USER*10,HNAME*32,CURDIR*60,TMPDIR*80
      INTEGER      MY_PID,MY_UID
      COMMON/ENVJ/ MY_PID,MY_UID,USER,HNAME,CURDIR,TMPDIR
C     ==================================================================
C     == TJLIMIT Job or CPU limit time (seconds)                      ==
C     ==--------------------------------------------------------------==
      REAL*8       TJLIMIT
      COMMON/ENVR/ TJLIMIT
C     ==================================================================
