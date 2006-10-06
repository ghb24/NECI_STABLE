C     ==================================================================
      BLOCK DATA TIMER
C     ==--------------------------------------------------------------==
C     == INITIALIZATION ARRAY FOR SUBROUTINE TIMING                   ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
      DATA NTIMER /0/
      DATA TNAM /200*'          '/
      DATA CTSUM,CTSTA,WTSUM,WTSTA,NCALLS
     *     /200*0.0D0,200*0.0D0,200*0.0D0,200*0.0D0,200*0/
      END


C     =============AJWT=========
      REAL*8 FUNCTION  TIGET(NSTIM)
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      INTEGER   NSTIM
      CHARACTER SNAME*10
C     Variables
      INTEGER   I
      REAL*8    TIMEC,TIMEF
      TIGET=CTSUM(NSTIM)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ============/AJWT=====================    

C     ==================================================================
      SUBROUTINE TISET(SNAME,NSTIM)
C     ==--------------------------------------------------------------==
C     ==  START THE CPU AND WALL CLOCK TIME FOR ROUTINE SNAME         ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      INTEGER   NSTIM
      CHARACTER SNAME*10
C     Variables
      INTEGER   I
      REAL*8    TIMEC,TIMEF
C     ==--------------------------------------------------------------==
      DO I=1,NTIMER
        IF(SNAME.EQ.TNAM(I)) THEN
          NSTIM=I
          GOTO 100
        ENDIF
      ENDDO
      NTIMER=NTIMER+1
      IF(NTIMER.GT.NTIMX) THEN
        WRITE(*,*) ' PARAMETER NTIMX IN INCLUDE FILE time.inc TOO ',
     *               'SMALL'
        CALL STOPGM('TISET',' PARAMETER')
      ENDIF
      NSTIM=NTIMER
      TNAM(NSTIM)=SNAME
  100 CONTINUE
      CTSTA(NSTIM)=TIMEC()*0.001D0
      WTSTA(NSTIM)=TIMEF()*0.001D0
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TIHALT(SNAME,NSTIM)
C     ==--------------------------------------------------------------==
C     ==  STOP THE CPU AND WALL CLOCK TIME FOR ROUTINE SNAME          ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      INTEGER   NSTIM
      CHARACTER SNAME*10
C     Variables
      INTEGER   I
      REAL*8    DCT,DWT,
     &          TIMEC,TIMEF
C     ==--------------------------------------------------------------==
      IF(NSTIM.EQ.0) THEN
        DO I=1,NTIMER
          IF(SNAME.EQ.TNAM(I)) THEN
            NSTIM=I
            GOTO 100
          ENDIF
        ENDDO
        WRITE(*,*) ' NO SUCH TIME HANDLE ',SNAME
        CALL STOPGM('TIHALT',' HANDLE   ')
  100   CONTINUE
      ELSE
        IF(SNAME.NE.TNAM(NSTIM)) THEN
          WRITE(*,*) ' INCONSISTENT HANDLE IN TIHALT ',SNAME,TNAM(NSTIM)
          CALL STOPGM('TIHALT',' HANDLE   ')
        ENDIF
      ENDIF
      DCT=TIMEC()*0.001D0-CTSTA(NSTIM)
      DWT=TIMEF()*0.001D0-WTSTA(NSTIM)
      CTSUM(NSTIM)=CTSUM(NSTIM)+DCT
      WTSUM(NSTIM)=WTSUM(NSTIM)+DWT
      NCALLS(NSTIM)=NCALLS(NSTIM)+1
      DO I=1,NTIMER
        CTSTA(I)=CTSTA(I)+DCT
        WTSTA(I)=WTSTA(I)+DWT
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TIPRI
C     ==--------------------------------------------------------------==
C     ==  PRINT TIMING INFO                                           ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Variables
      INTEGER J,IMAX,I
      REAL*8  CPMAX,CPTHRS,TOTALCT,TOTALWT
C     ==--------------------------------------------------------------==
      WRITE(6,'(/,/,1X,64(1H*))')
      WRITE(6,'(" *",62X,"*")')
      WRITE(6,'(" *",27X,A,27X,"*")') ' TIMING '
      WRITE(6,'(" *",62X,"*")')
      WRITE(6,'(1X,64(1H*))')
      WRITE(6,'(A,A)') ' SUBROUTINE',
     *      '            CALLS         CPU TIME        ELAPSED TIME'
      TOTALCT=0.D0
      TOTALWT=0.D0
      DO J=1,NTIMER
        IMAX=0
        CPMAX=-0.1D0
        DO I=1,NTIMER
          IF(CTSUM(I).GT.CPMAX) THEN
            IMAX=I
            CPMAX=CTSUM(I)
          ENDIF
        ENDDO
        IF(J.EQ.1) CPTHRS=CPMAX/100.D0
        IF(CPMAX.LT.CPTHRS) GOTO 100
        I=IMAX
        WRITE(6,'(1X,A10,7X,I13,5X,F12.2,8X,F12.2)') 
     *                    TNAM(I),NCALLS(I),CTSUM(I),WTSUM(I)
        TOTALCT=TOTALCT+CTSUM(I)
        TOTALWT=TOTALWT+WTSUM(I)
        CTSUM(I)=-1.0
      ENDDO
  100 CONTINUE
      WRITE(6,'(1X,64(1H-))')
      WRITE(6,*) "TOTAL TIME",TOTALCT,TOTALWT
C      WRITE(6,'(1X,"TOTAL TIME",T34,F12.2,8X,F12.2)')
C     &                    TOTALCT,TOTALWT
C      WRITE(6,'(1X,64(1H*),/)')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TISTART(TIME1)
C     ==--------------------------------------------------------------==
C     == CALLED ONCE AT THE BEGINNING IN ORDER TO START THE TIMER     ==
C     == AND DO SOME INITIALISATION [CRAY]                            ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'
C     Arguments
      REAL*8 TIME1
C     Variables
      REAL*8  TIMEC
      INTEGER IFIRST
      DATA    IFIRST /0/
C     ==--------------------------------------------------------------==
      IF(IFIRST.EQ.0) THEN
        TIME1=TIMEC()*1.D-3
        IFIRST=1
        TJSTART=TIME1
        TJSLOOP=TJSTART
      ELSE
        CALL STOPGM('TISTART','CALLED ONCE')
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TILIMIT(TSTOP)
C     ==--------------------------------------------------------------==
C     == .TRUE. IF THE JOB TIME IS TOO SHORT FOR ANOTHER LOOP         ==
C     == A LOOP IS DEFINED BY 2 CALLS TO TILIMIT (TESTEX)             ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'        !TJSTART
      INCLUDE 'envj.inc'        !TJLIMIT
C     Arguments
      LOGICAL TSTOP
C     Variables
      REAL*8  TJELOOP,TREMAIN,
     &        TIMEC
C     ==--------------------------------------------------------------==
      IF(TIMESTOP) THEN
        TSTOP=.TRUE.
        RETURN
      ENDIF
      TSTOP=.FALSE.
      TIMESTOP=.FALSE.
      IF(TJLIMIT.EQ.0) RETURN
C     End of time for the last loop
      TJELOOP=TIMEC()*1.D-3
C     Total time for the last loop
      TIMELOOP=TJELOOP-TJSLOOP
C     Remaining job time 
      TREMAIN=TJLIMIT-(TJELOOP-TJSTART)
C     We add 10%.
      IF(TREMAIN.LT.(1.01D0*TIMELOOP)) THEN
        TSTOP=.TRUE.
        TIMESTOP=.TRUE.
      ENDIF
      TJSLOOP=TJELOOP
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE TILIMEX
C     ==--------------------------------------------------------------==
C     == EXIT IF THE JOB TIME LIMIT IS REACHED                        ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'time.inc'        !TJSTART
      INCLUDE 'envj.inc'        !TJLIMIT
C     Variables
      REAL*8    TREMAIN,TIME1,
     &          TIMEC
      CHARACTER DATX*26
C     ==--------------------------------------------------------------==
      IF(.NOT.TIMESTOP) RETURN
C     Remaining job time
      TIME1=TIMEC()*1.D-3
      TREMAIN=TJLIMIT-(TIME1-TJSTART)
      CALL DATUM(DATX)
      WRITE(6,'(1X,64(1H*))') 
      WRITE(6,'(1X,"*",62X,"*")')
      WRITE(6,'(1X,"*",12X,A,T65,"*")') 
     *     'JOB LIMIT TIME EXCEEDED FOR A NEW LOOP'
      WRITE(6,'(1X,"*",62X,"*")')
      WRITE(6,'(1X,64(1H*))') 
      WRITE(6,'(1X,"*",A,T43,F13.3," SECONDS *")')
     &     ' THE JOB TIME LIMIT IS:',TJLIMIT
      WRITE(6,'(1X,"*",A,T43,F13.3," SECONDS *")')
     &     ' THE TIME OF THE LAST LOOP IS:',TIMELOOP
      WRITE(6,'(1X,"*",A,T43,F13.3," SECONDS *")')
     &     ' THE REMAINING JOB TIME IS:',TREMAIN
      WRITE(6,'(A,A,A,A)')
     *  ' *     ',' THE COMMAND WAS ISSUED AT ',DATX(1:24),'      *'
      WRITE(6,'(1X,"*",62X,"*")')
      WRITE(6,'(1X,64(1H*))')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
