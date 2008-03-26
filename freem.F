C     ==================================================================
      SUBROUTINE FREEM(IP_POINT)
C     ==--------------------------------------------------------------==
      use MemoryManager, only: LogMemDealloc
      use MemoryManager, only: LookupPointer
      IMPLICIT NONE
      REAL*8  POINTDUMMY
      POINTER (IP_POINT,POINTDUMMY)
#ifdef POINTER8
      INTEGER*8 I_AD1,I_AD2
#else
      INTEGER  I_AD1,I_AD2
#endif
!      INCLUDE 'system.h'
      INCLUDE 'memc.inc'
      LOGICAL ERROR
      REAL*8  XM(0:*)
      POINTER (IP_XM,XM)
#ifdef CRAY
      INTEGER IERR
#endif
      INTEGER I,II,LENGTH,iTag
C     ==--------------------------------------------------------------==
C     ==  FREE MEMORY ASSOCIATED WITH POINTER IP_POINT                ==
C     ==--------------------------------------------------------------==
#if defined(CRAY) || defined(__NEC)
      IP_XM=IP_POINT-1
#else
      IP_XM=IP_POINT-8
#endif
      ! Log using the new memory book-keeping system.  Need to find out
      ! the  position in the log first.
      ! Assume it wasn't logged correctly (e.g. due to a full log), 
      ! and correct this assumption if was logged in LookupPointer.
      iTag=-1
      do i=1,size(LookupPointer)
          if (LookupPointer(i).eq.IP_XM) then
              iTag=i
              exit
          end if
      end do
      call LogMemDealloc('unknown',iTag)
C..Statistics and memory checks
      DO I=1,I_BLOCK
        IF(IP_XM.EQ.I_AD1(I)) THEN
          II=I
          GOTO 100
        ENDIF
      ENDDO
      WRITE(*,*) ' FREEM | Could not find pointer ',IP_XM
      CALL PRINT_MEMORY
      CALL STOPGM('FREEM',' NO SUCH POINTER')
  100 CONTINUE
      LENGTH=I_LEN(II)
      ERROR=ZFIRST.NE.XM(0).OR.ZLAST.NE.XM(LENGTH+1)
      IF(ERROR) CALL MEMORY_CHECK
      IF(MEMCHK) CALL MEMORY_CHECK
      DO I=II,I_BLOCK-1
        I_AD1(I)=I_AD1(I+1)
        I_AD2(I)=I_AD2(I+1)
        I_LEN(I)=I_LEN(I+1)
        I_TAG(I)=I_TAG(I+1)
      ENDDO
      I_BLOCK=I_BLOCK-1
      MEMTOT=MEMTOT-LENGTH
c..free memory
#ifdef CRAY
      CALL HPDEALLC(IP_XM,IERR,0)
      IF(IERR.NE.0) THEN
        WRITE(*,*) ' FREEM | Error in deallocation of memory '
        CALL STOPGM('FREEM',' ')
      ENDIF
#elif __IBM
      CALL FREE(%VAL(IP_XM))
#elif __SR2201
      CALL FREE(%VAL(IP_XM))
#elif __NEC
      CALL FREEW(IP_XM)
#else
      CALL FREE(IP_XM)
#endif
C     ==--------------------------------------------------------------==
      RETURN
      END

C     ==================================================================
      SUBROUTINE MEMORY_CHECK
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
#ifdef POINTER8
      INTEGER*8 I_AD1,I_AD2
#else
      INTEGER   I_AD1,I_AD2
#endif
!      INCLUDE 'system.h'
      INCLUDE 'memc.inc'
      LOGICAL ERROR
      REAL*8  X1(*),X2(*)
      POINTER (IP_X1,X1),(IP_X2,X2)
      INTEGER I
C     ==--------------------------------------------------------------==
C     ==  CHECK BORDERS OF ALLOCATED BLOCKS OF MEMORY                 ==
C     ==--------------------------------------------------------------==
      ERROR=.FALSE.
      DO I=1,I_BLOCK
        IP_X1=I_AD1(I)
        IP_X2=I_AD2(I)
        IF(IP_X2.NE.0) THEN
            IF(ZFIRST.NE.X1(1)) ERROR=.TRUE.
            IF(ZLAST.NE.X2(1)) ERROR=.TRUE.
        ENDIF
      ENDDO
      IF(.NOT.ERROR) THEN
        IF(MEMCHK)
     *  WRITE(6,'(A)') ' MEMORY_CHECK| ALL ARRAY BOUNDARIES ARE OK '
!        IF(PARENT.AND.MEMCHK)
!     *  WRITE(6,'(A)') ' MEMORY_CHECK| ALL ARRAY BOUNDARIES ARE OK '
      ELSE
        CALL PRINT_MEMORY
        CALL STOPGM('MEMORY_CHECK',' CORRUPTED MEMORY')
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END

C     ==================================================================
      SUBROUTINE PRINT_MEMORY
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
#ifdef POINTER8
      INTEGER*8 I_AD1,I_AD2
#else
      INTEGER   I_AD1,I_AD2
#endif
!      INCLUDE 'system.h'
      INCLUDE 'memc.inc'
#ifdef POINTER8
      INTEGER*8 I_SORT(MEM_LEN)
#else
      INTEGER   I_SORT(MEM_LEN)
#endif      
      INTEGER   INDEX(MEM_LEN), I,II,IC1,IC2
      REAL*8    X1(*),X2(*)
      POINTER   (IP_X1,X1),(IP_X2,X2)
      CHARACTER TTT(2)*10
#ifdef PARALLEL
      INTEGER   IFORM
      CHARACTER FORMFREEM*100
#endif
      DATA TTT /'        OK',' CORRUPTED'/
C     ==--------------------------------------------------------------==
#ifdef POINTER8
      CALL I8COPY(I_BLOCK,I_AD1,1,I_SORT,1)
      CALL KB07AI8(I_SORT,I_BLOCK,INDEX)
#else
      CALL ICOPY(I_BLOCK,I_AD1,1,I_SORT,1)
      CALL KB07AI(I_SORT,I_BLOCK,INDEX)
#endif
#ifdef PARALLEL
      IFORM=INT(LOG10(REAL(NPROC)))+1
      WRITE(FORMFREEM,'(A,I1,A,I1,A,I2,A)')
     &     '(//,1X,"==",I',IFORM,'":",I',IFORM,',',63-2*IFORM,'("="))'
      WRITE(6,FORMFREEM) ME,NPROC
#else
      WRITE(6,'(//,1x,66("="))')
#endif
      WRITE(6,'(A,A)') ' NAME            FIRST ADDRESS    ',
     *     '       LAST ADDRESS        LENGTH  '
      DO I=1,I_BLOCK
        II=INDEX(I)
        IP_X1=I_AD1(II)
        IP_X2=I_AD2(II)
        IF(IP_X2.NE.0) THEN
           IF(ZFIRST.NE.X1(1)) THEN
             IC1=2
           ELSE
             IC1=1
           ENDIF
           IF(ZLAST.NE.X2(1)) THEN
             IC2=2
           ELSE
             IC2=1
           ENDIF
        ELSE
         IC1=1
         IC2=1
        ENDIF
      
        WRITE(6,'(1X,A10,I14,A10,I14,A10,I8)') 
     *      I_TAG(II),I_AD1(II),TTT(IC1),I_AD2(II),TTT(IC2),I_LEN(II)
      ENDDO
      WRITE(6,'(1x,66("-"))')
      WRITE(6,'(1X,A12,I5,A1,6X,A12,4X,I10,A2,F7.1,A7)')
     &     '[PEAK NUMBER',NUMPEAK,']',
     &     'PEAK MEMORY',MEMPEAK,' =',MEMPEAK*8/1.D6,' MBytes'
      WRITE(6,'(1X,A12,I5,A1,6X,A12,4X,I10,A2,F7.1,A7)')
     &     '[ALL. NUMBER',I_BLOCK,']',
     &     'TOTAL MEMORY',MEMTOT,' =',MEMTOT*8/1.D6,' MBytes'
      WRITE(6,'(1x,66("="))')
C     ==--------------------------------------------------------------==
      RETURN
      END

C     ==================================================================
      SUBROUTINE PMAXMEM
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
#ifdef POINTER8
      INTEGER*8 I_AD1,I_AD2
#else
      INTEGER   I_AD1,I_AD2
#endif
!      INCLUDE 'system.h'
      INCLUDE 'memc.inc'
      INTEGER INDEX(MBM),I_SORT(MBM), I,I1,I2
      INTEGER IMAX
C     ==--------------------------------------------------------------==
C      IF(I_BLOCK.LT.MBM) THEN
C         IMAX=I_BLOCK
C      ELSE
         IMAX=MBM
C      ENDIF
      CALL ICOPY(IMAX,I_LENB,1,I_SORT,1)
      CALL KB07AI(I_SORT,IMAX,INDEX)
      WRITE(6,'(/,1x,64("="))')
      WRITE(6,'(20X,A)') '  BIG MEMORY ALLOCATIONS '
      DO I=IMAX,1,-2
        I1=INDEX(I)
        I2=INDEX(I-1)
        WRITE(6,'(1X,A10,4X,I10,16X,A10,4X,I10)') IB_TAG(I1),I_LENB(I1),
     *          IB_TAG(I2),I_LENB(I2)
      ENDDO
      WRITE(6,'(1x,64("-"))')
      WRITE(6,'(1X,A12,I5,A1,6X,A11,3X,I10,A2,F7.1,A7)')
     &     '[PEAK NUMBER',NUMPEAK,']',
     &     'PEAK MEMORY',MEMPEAK,' =',MEMPEAK*8/1.D6,' MBytes'
      WRITE(6,'(1x,64("="))')
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
!! AJWT 26/3/07 A routine to integrate Fortran 90 memory and array management with FREEM
!! NB This does not deallocate memory.
!! this should be called in statements like
!!
!!  TYPE(HElement), Pointer :: UMat(:,:)
!!  CALL MemDealloc(UMat)
!!  DEALLOCATE(UMat)
!!
!!  the memory should have been ALLOCATEd and MemAlloc called previously.
C     ==================================================================
      SUBROUTINE MemDealloc(XM)
C     ==--------------------------------------------------------------==
      use MemoryManager, only: LogMemDealloc
      use MemoryManager, only: LookupPointer
      IMPLICIT NONE
      REAL*8  POINTDUMMY
      POINTER (IP_POINT,POINTDUMMY)
#ifdef POINTER8
      INTEGER*8 GET_ADDR
      INTEGER*8 I_AD1,I_AD2
      INTEGER*8 IP_XM
#else
      INTEGER GET_ADDR
      INTEGER  I_AD1,I_AD2
      INTEGER*8 IP_XM
#endif
!      INCLUDE 'system.h'
      INCLUDE 'memc.inc'
      LOGICAL ERROR
      REAL*8  XM(0:*)
#ifdef CRAY
      INTEGER IERR
#endif
      INTEGER I,II,LENGTH,iTag
C     ==--------------------------------------------------------------==
C     ==  FREE MEMORY ASSOCIATED WITH POINTER IP_POINT                ==
C     ==--------------------------------------------------------------==
      IP_XM=GET_ADDR(XM)
      ! Log using the new memory book-keeping system.  Need to find out
      ! the  position in the log first.
      do i=1,size(LookupPointer)
          if (LookupPointer(i).eq.IP_XM) then
              iTag=i
              exit
          end if
      end do
C..Statistics and memory checks
      DO I=1,I_BLOCK
        IF(IP_XM.EQ.I_AD1(I)) THEN
          II=I
          GOTO 100
        ENDIF
      ENDDO
      WRITE(*,*) ' FREEM | Could not find pointer ',IP_XM
      CALL PRINT_MEMORY
      CALL STOPGM('FREEM',' NO SUCH POINTER')
  100 CONTINUE
      LENGTH=I_LEN(II)
      ERROR=.FALSE.
      IF(ERROR) CALL MEMORY_CHECK
      IF(MEMCHK) CALL MEMORY_CHECK
      DO I=II,I_BLOCK-1
        I_AD1(I)=I_AD1(I+1)
        I_AD2(I)=I_AD2(I+1)
        I_LEN(I)=I_LEN(I+1)
        I_TAG(I)=I_TAG(I+1)
      ENDDO
      I_BLOCK=I_BLOCK-1
      MEMTOT=MEMTOT-LENGTH
      RETURN
      END
