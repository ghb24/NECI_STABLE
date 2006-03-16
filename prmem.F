C     ==================================================================
      SUBROUTINE PRMEM(SUBR)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
      INCLUDE 'envj.inc'
C     Arguments
      CHARACTER SUBR*10
C     Variables
      CHARACTER STRING*140
#if defined(CRAY)
      INTEGER ISS,
     &        IHPSTAT
C     ==--------------------------------------------------------------==
      ISS=IHPSTAT(1)/1000
      WRITE(6,'(A4,A10,A24,I16,A11)') 
     *    ' ***',SUBR,'|   CURRENT HEAP LENGTH ',ISS,' KWords ***'
#elif defined(__alpha)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)')
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS   '
      WRITE(STRING,'(A,I10,A)') 'ps -o rssize=RSS -p ',MY_PID,
     *  ' | sed -e "/RSS/d" | sed -e "s/$/Bytes    ***/" '
      CALL SYSTEM(STRING)
#elif defined(__SUN)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)')
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS '
      WRITE(STRING,'(A,I10,A)') 'ps -o rss -p ',MY_PID,
     *  ' | awk ''NR == 2 { printf("%7d KBYTES ***\\n", $1) }'' '
      CALL SYSTEM(STRING)
#elif defined(__NEC)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)')
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS '
      CALL FLUSH(6)
      WRITE(STRING,'(A,I10,A,A)') 'ps -lp ',MY_PID,
     *  ' | grep -v PID |  cut -c49-55 | ',
     *  ' awk '' { print $1/25 " MBytes ***"}'''
      CALL SYSTEM(STRING)
      WRITE(6,*)
#elif defined(_vpp_)
#if defined(_vpp_old)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)')
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS '
      WRITE(STRING,'(A,I10,A)') 'ps -lc -p ',MY_PID,
     *' | sed -e "/ SZ/d" | cut -c42-48 - | sed -e "s/$/(4kB) pages/"'
      CALL SYSTEM(STRING)
#else
      INTEGER VRTSIZE,          !OUT process virtual size
     &        STKSIZE,          !OUT current stack size
     &        BRKSIZE           !OUT current heap  size
C     ==--------------------------------------------------------------==
      CALL GETMEMSIZE(VRTSIZE,STKSIZE,BRKSIZE)
      VRTSIZE=VRTSIZE*32
      WRITE(6,'(A,A,A,I7,A)')
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS ',VRTSIZE,
     &    ' kBYTES ***'
#endif
#elif defined(__SGI)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)') 
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS '
      CALL FLUSH(6)
      WRITE(STRING,'(A,I10,A,A,$)') 'ps -l -p ',MY_PID,
     *  ' | awk ''NR == 2 { print $10 }'' ',
     *  ' | awk -F: ''{ printf("%6d KBYTES  ***", $2*16) }'' '
C     *  ' | awk -F: ''{ print $1*16 " KBYTES   ***" }'' '
      CALL SYSTEM(STRING)
      WRITE(6,*)
#elif defined(__HP)
      real*8 rdata,vdata,rstack,vstack
      call memme(vdata,vstack)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,F7.2,A,A,F7.2,A)')
     *    ' ***',SUBR,'| PROGSIZE DATA=',vdata,' MByte',
     *    '  /  STACK=',vstack,' MByte'
#elif defined(__SR2201)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)')
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS '
      WRITE(STRING,'(A,I10,A,A)') 'ps -o rss -p ',MY_PID,
     *  ' | sed -e "/RSS/d" | sed -e "s/K.*$/ KBytes ***/" ',
     *  ' | sed -e "s/M.*$/ MBytes ***/" '
      CALL SYSTEM(STRING)
#elif defined(__IBM)
C     ==--------------------------------------------------------------==
      WRITE(6,'(A,A,A,$)') 
     *    ' ***',SUBR,'| THE NEW SIZE OF THE PROGRAM IS '
      WRITE(STRING,'(A,I12,A)') 'ps -Fvsz -p ',MY_PID,
     *  ' | sed -e "/VSZ/d" | sed -e "s/$/ KBytes ***/" '
      CALL SYSTEM(STRING)
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
