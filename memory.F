C     ==================================================================
      SUBROUTINE MEMORY(IP_POINT,LENGTH,TAG)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
!      INCLUDE 'system.h'
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
      INTEGER*8 MALLOC,GET_ADDR
#elif defined(__NEC)
      INTEGER   I_AD1,I_AD2
      INTEGER   MALLOCW,GET_ADDR
#else
      INTEGER   I_AD1,I_AD2
      INTEGER   MALLOC,GET_ADDR
#endif
      REAL*8    XM(0:LENGTH+1)
      POINTER   (IP_XM,XM)
#ifdef CRAY
      INTEGER   IERR1,IERR2
#endif
      CHARACTER LINE*100
      INTEGER   I,IERROR,LEN,IA,IE,
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
        CALL STOPGM('MEMORY','NEGATIVE OR TOO BIG VALUE')
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
          CALL STOPGM('MEMORY','TOO BIG VALUE')
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
        CALL PRMEM('    MEMORY')
        CALL PMAXMEM
        CALL XSTRING(TAG,IA,IE)
        WRITE(LINE,'(A,A,A)') 'ALLOCATION FAILED (',TAG(IA:IE),')'
        IE=IE+19+1
        CALL STOPGM('MEMORY',LINE(1:IE))
      ENDIF
      IP_POINT=GET_ADDR(XM(1))
C..Statistics and memory checks
      XM(0)=ZFIRST
      XM(LENGTH+1)=ZLAST
      I_BLOCK=I_BLOCK+1
      IF(I_BLOCK.GT.MEM_LEN) I_BLOCK=MEM_LEN
      IF(NUMPEAK.LT.I_BLOCK) NUMPEAK=I_BLOCK
      I_AD1(I_BLOCK)=GET_ADDR(XM(0))
      I_AD2(I_BLOCK)=GET_ADDR(XM(LENGTH+1))
      I_LEN(I_BLOCK)=LENGTH
      I_TAG(I_BLOCK)=TAG
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


!! AJWT 26/3/07 A routine to integrate Fortran 90 memory and array management with MEMORY
!! NB This does not allocate memory.
!! this should be called after statements like
!!
!!  TYPE(HElement), Pointer :: UMat(:,:)
!!  INTEGER ierr
!!  ALLOCATE(UMat(10,20),STAT=ierr
!!  CALL MemAlloc(ierr,UMat,10*20*HElementSize,"UMat")
!!

C     ==================================================================
      SUBROUTINE MemAlloc(IERROR,XM,LENGTH,TAG)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
!      INCLUDE 'system.h'
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
      INTEGER*8 MALLOC,GET_ADDR
#elif defined(__NEC)
      INTEGER   I_AD1,I_AD2
      INTEGER   MALLOCW,GET_ADDR
#else
      INTEGER   I_AD1,I_AD2
      INTEGER   MALLOC,GET_ADDR
#endif
      REAL*8    XM(LENGTH)
#ifdef CRAY
      INTEGER   IERR1,IERR2
#endif
      CHARACTER LINE*100
      INTEGER   I,IERROR,LEN,IA,IE,
     &          I_FIRST
      DATA      I_FIRST /0/
      SAVE      I_FIRST


C      WRITE(14,*) TAG,LENGTH


      IF(IERROR.NE.0) THEN
        WRITE(*,'(/,1X,64("*"))')
!Used to be processor ME
        WRITE(*,'(1X,A9,I6,A14,I12,A28,A10)') 'PROCESSOR',1,
     *             ' ALLOCATION OF',LENGTH,' WORDS OF MEMORY FAILED
     *             FOR ',TAG
        WRITE(*,'(1X,64("*"),/)')
        CALL PRMEM('    MEMORY')
        CALL PMAXMEM
        CALL XSTRING(TAG,IA,IE)
        WRITE(LINE,'(A,A,A)') 'ALLOCATION FAILED (',TAG(IA:IE),')'
        IE=IE+19+1
        CALL STOPGM('MEMORY',LINE(1:IE))
      ENDIF
      IP_POINT=GET_ADDR(XM(1))
C..Statistics and memory checks
      I_BLOCK=I_BLOCK+1
      IF(I_BLOCK.GT.MEM_LEN) I_BLOCK=MEM_LEN
      IF(NUMPEAK.LT.I_BLOCK) NUMPEAK=I_BLOCK
      I_AD1(I_BLOCK)=GET_ADDR(XM(1))
!GET_ADDR(XM(1))
      I_AD2(I_BLOCK)=0
!GET_ADDR(XM(LENGTH+1))
      I_LEN(I_BLOCK)=LENGTH
      I_TAG(I_BLOCK)=TAG
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
