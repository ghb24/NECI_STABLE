C     ==================================================================
C     == MEM_LEN  MAXIMUM NUMBER OF ALLOCATED ARRAYS                  ==
C     == MBM      NUMBER OF BIGGEST ALLOCATIONS (TO DISPLAY)          ==
C     == ZFIRST   VALUE OF XM(0) TO CHECK MEMORY                      ==
C     == ZLAST    VALUE OF XM(N+1)                                    ==
C     ==--------------------------------------------------------------==
      INTEGER   MEM_LEN,MBM
      PARAMETER (MEM_LEN=65536)
      PARAMETER (MBM=10)
      REAL*8    ZFIRST,ZLAST
C     PARAMETER (ZFIRST=Z'EFFFFFFFFFFFFFFF')
C     PARAMETER (ZLAST=Z'0FFFFFFFFFFFFFFF')
      PARAMETER (ZFIRST=1234567890.1234567E0)
      PARAMETER (ZLAST=9876543210.9876543E0)
C     ==--------------------------------------------------------------==
C     == I_LEN   LENGTH          FOR EACH ALLOCATED ARRAY             ==
C     == I_AD1   INITIAL ADDRESS FOR EACH ALLOCATED ARRAY             ==
C     == I_AD2   FINAL   ADDRESS FOR EACH ALLOCATED ARRAY             ==
C     == I_BLOCK NUMBER OF ALLOCATIONS                                ==
C     == I_LENB  LENGTH FOR THE MBM BIGGEST ALLOCATIONS               ==
C     == IB_MIN  INDEX OF THE SMALLEST ALLOCATION IN THE MBM BIGGEST  ==
C     == LB_MIN  LENGTH OF THE SMALLEST ALLOCATION IN THE MBM BIGGEST ==
C     == NUMPEAK PEAK NUMBER OF ALLOCATIONS                           ==
C     == MEMTOT  TOTAL ALLOCATED MEMORY                               ==
C     == MEMPEAK PEAK ALLOCATED MEMORY   = MAX(MEMTOT)                ==
C     ==--------------------------------------------------------------==
      INTEGER      I_LEN,I_LENB,I_BLOCK,IB_MIN,LB_MIN,
     &             NUMPEAK,MEMTOT,MEMPEAK
      DIMENSION    I_AD1(MEM_LEN),I_AD2(MEM_LEN),I_LEN(MEM_LEN),
     &             I_LENB(MBM)
      COMMON/MEMCNECI2/I_AD1,I_AD2,I_LEN,I_LENB,
     &             I_BLOCK,IB_MIN,LB_MIN,NUMPEAK,MEMTOT,MEMPEAK
C     ==--------------------------------------------------------------==
C     == I_TAG   TAG FOR EACH ALLOCATION                              ==
C     == IB_TAG  TAG FOR THE MBM BIGGEST ALLOCATIONS                  ==
C     ==--------------------------------------------------------------==
      CHARACTER*10 I_TAG,IB_TAG
      COMMON/MEMCNECI3/I_TAG(MEM_LEN),IB_TAG(MBM)
C     ==--------------------------------------------------------------==
C     == MEMCHK .TRUE. MEMORY_CHECK DISPLAY MESSAGES                  ==
C     ==--------------------------------------------------------------==
      LOGICAL MEMCHK
      COMMON /MEMCNECI1/ MEMCHK
C     ==================================================================
