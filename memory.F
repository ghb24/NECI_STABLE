C     ==================================================================
      SUBROUTINE N_MEMORY(IP_POINT,LENGTH,TAG)
C     ==--------------------------------------------------------------==
      use MemoryManager, only: LogMemAlloc,LeaveMemoryManager
      use MemoryManager, only: LookupPointer,li
      IMPLICIT NONE
      INCLUDE 'memc.inc'
      INCLUDE 'irat.inc'
C     Arguments
      REAL*8    POINTDUMMY
      POINTER   (IP_POINT,POINTDUMMY)
      INTEGER   LENGTH
      CHARACTER TAG*(*),TAG2*(10)
C     Variables
#if defined(POINTER8)
      INTEGER*8 I_AD1,I_AD2
      INTEGER*8 MALLOC,LOC
#elif defined(__NEC)
      INTEGER   I_AD1,I_AD2
      INTEGER   MALLOCW,LOC
#else
      INTEGER   I_AD1,I_AD2
      INTEGER   MALLOC,LOC
#endif
      REAL*8    XM(0:LENGTH+1)
      POINTER   (IP_XM,XM)
#ifdef CRAY
      INTEGER   IERR1,IERR2
#endif
      CHARACTER LINE*100
      INTEGER   I,IERROR,LEN,IA,IE,iTag,
     &          I_FIRST
      DATA      I_FIRST /0/
      SAVE      I_FIRST


C      WRITE(14,*) TAG,LENGTH


C     ==--------------------------------------------------------------==
C     ==  ALLOCATION OF LENGTH * REAL*8 WORDS OF MEMORY               ==
C     ==--------------------------------------------------------------==
      IF(I_FIRST.EQ.0) THEN
        I_BLOCK=0
        I_FIRST=1
        LB_MIN=0
        IB_MIN=1
        TAG2="       "
        DO I=1,MBM
          I_LENB(I)=0
          I_TAG(I)=TAG2  
        ENDDO
        NUMPEAK=0
        MEMPEAK=0
        MEMTOT=0
      ENDIF
      IERROR=0
      LEN=LENGTH+2
      IF(LEN.LT.0) THEN
        WRITE(*,*) ' MEMORY| MEMORY REQUIRED:',LEN,' WORDS'
        CALL Stop_All('MEMORY','NEGATIVE OR TOO BIG VALUE')
      ENDIF
#ifdef CRAY
C     Integer=integer*8
      CALL HPALLOC(IP_XM,LEN,IERR1,0)
      IF(IERR1.NE.0) THEN
        CALL HPSHRINK
        CALL HPALLOC(IP_XM,LEN,IERR2,0)
        IF(IERR2.NE.0) IERROR=1
      ENDIF
#elif defined(__NEC)
C     Integer=integer*8
      LEN=MAX(LEN,400)
      IP_XM=MALLOCW(LEN)
#else
C     Check if 8 * LEN is in the range of integer.
      IF(IRAT.EQ.2) THEN
C       Integer=integer*4 : 2^31=2147483648, 2^28=268435456
        IF(LEN.GT.2**28) THEN
          WRITE(*,*) ' MEMORY| MEMORY REQUIRED:',LEN,' WORDS'
          CALL Stop_All('MEMORY','TOO BIG VALUE')
        ENDIF
      ENDIF
      LEN = 8 * LEN
      LEN=MAX(LEN,400)
#if defined(__IBM) || defined(__SR2201)
      IP_XM=MALLOC(%VAL(LEN))
#else
      IP_XM=MALLOC(LEN)
#endif
CC AJWT 5/6/06 from CPMD 3.9.11
C detect malloc failure. comparing a pointer with a possibly signed
C integer
C may be a problem. on some machines the 'heap' grows down so successful
C allocations may result in a 'negative' value.
C we are being careful here. on the following machines a failed malloc
C _is_
C indicated by a NULL pointer, for the rest we fail on negative as well.
C AK 2005/04/16.
CC /AJWT
#if defined(__HP) || defined(__Linux) || defined(__alpha) || defined(__OSX) || defined(__WINNT)
      IF(IP_XM.EQ.0) THEN
        IERROR=1
      ENDIF
#else
      IF(IP_XM.LE.0) THEN
        IERROR=1
      ENDIF
#endif
#endif
      IF(IERROR.NE.0) THEN
        WRITE(*,'(/,1X,64("*"))')
!Used to be processor ME
        WRITE(*,'(1X,A9,I6,A14,I12,A28,A10)') 'PROCESSOR',1,
     *             ' ALLOCATION OF',LENGTH,' WORDS OF MEMORY FAILED
     *             FOR ',TAG
        WRITE(*,'(1X,64("*"),/)')
        CALL LeaveMemoryManager
        CALL N_XSTRING(TAG,IA,IE)
        WRITE(LINE,'(A,A,A)') 'ALLOCATION FAILED (',TAG(IA:IE),')'
        IE=IE+19+1
        CALL Stop_All('MEMORY',LINE(1:IE))
      ENDIF
      IP_POINT=LOC(XM(1))
C..Statistics and memory checks
      XM(0)=ZFIRST
      XM(LENGTH+1)=ZLAST
      I_BLOCK=I_BLOCK+1
      IF(I_BLOCK.GT.MEM_LEN) I_BLOCK=MEM_LEN
      IF(NUMPEAK.LT.I_BLOCK) NUMPEAK=I_BLOCK
      I_AD1(I_BLOCK)=LOC(XM(0))
      I_AD2(I_BLOCK)=LOC(XM(LENGTH+1))
      I_LEN(I_BLOCK)=LENGTH
      I_TAG(I_BLOCK)=TAG
      ! Store in the new memory book-keeping system.
      call LogMemAlloc(Tag,Length,8,'unknown',iTag)
      if (iTag.ne.-1) LookupPointer(iTag)=int(I_AD1(I_BLOCK),kind=li)
      MEMTOT=MEMTOT+LENGTH
      IF(MEMPEAK.LT.MEMTOT) MEMPEAK=MEMTOT
      IF(LENGTH.GT.LB_MIN) THEN
        I_LENB(IB_MIN)=LENGTH
        IB_TAG(IB_MIN)=TAG
        LB_MIN=LENGTH
        DO I=MBM,1,-1
          IF(LB_MIN.GE.I_LENB(I)) THEN
            LB_MIN=I_LENB(I)
            IB_MIN=I
          ENDIF
        ENDDO
      ENDIF
      IF(MEMCHK) CALL MEMORY_CHECK
C      WRITE(6,"(A10,Z)") TAG,IP_POINT
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

!  Used to be in readsr.F
C     ==================================================================
      SUBROUTINE N_XSTRING(STRING,IA,IE)
C     ==--------------------------------------------------------------==
C     == Give index of the first and last non-blank character         ==
C     == blank characters are " " and CHAR(9) (tab)                   ==
C     == more blank characters are char(13) (cr) and CHAR(10) (newl)  ==
C     == STRING (INPUT)                                               ==
C     == IA,IE  (OUTPUT)                                              ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER STRING*(*)
      INTEGER IA,IE
C     Variables
      INTEGER SLEN,I
C     ==--------------------------------------------------------------==
      SLEN=LEN(STRING)
      IA=1
      DO I=1,SLEN
         IF(STRING(I:I).NE.' '.AND.
     &       STRING(I:I).NE.CHAR(9).AND.      !tab character
     &       STRING(I:I).NE.CHAR(10).AND.     !newline character
     &       STRING(I:I).NE.CHAR(13) ) THEN   !carriage return character
           IA=I
           GOTO 10
         ENDIF
      ENDDO
 10   CONTINUE
      DO I=IA,SLEN
         IF(STRING(I:I).EQ.' '.OR.
     &       STRING(I:I).EQ.CHAR(0).OR. !\0 character
     &       STRING(I:I).EQ.CHAR(9)) THEN !tab character
           IE=I-1
           GOTO 20
         ENDIF
      ENDDO
      IE=SLEN
 20   CONTINUE
C     ==--------------------------------------------------------------==
      RETURN
      END
