C     ==================================================================
C     == ENVIRONMENT                                                  ==
C     ==================================================================
C     == USER   User name                                             ==
C     == HNAME  Computer name                                         ==
C     == CURDIR Current Directory                                     ==
C     == TMPDIR Temporary Directory                                   ==
C     ==--------------------------------------------------------------==
C     == MY_PID Process ID Number                                     ==
C     == MY_UID User ID Number                                        ==
C     ==--------------------------------------------------------------==
C     == TJLIMIT Job or CPU limit time (seconds)                      ==
C     ==--------------------------------------------------------------==
      CHARACTER     USER*10,HNAME*64,CURDIR*255,TMPDIR*255
      INTEGER       MY_PID,MY_UID
      REAL*8        TJLIMIT
      COMMON/ENVJC/ USER,HNAME,CURDIR,TMPDIR
      COMMON/ENVJI/ MY_PID,MY_UID
      COMMON/ENVJR/ TJLIMIT
C     ==================================================================
