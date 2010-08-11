C     ==================================================================
      SUBROUTINE N_FREEM(IP_POINT)
C     ==--------------------------------------------------------------==
      use MemoryManager, only: LogMemDealloc
      use MemoryManager, only: LookupPointer
      use legacy_memory_data
      IMPLICIT NONE
      REAL*8  POINTDUMMY
      POINTER (IP_POINT,POINTDUMMY)
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
      LookupPointer(iTag)=0
      call LogMemDealloc('unknown',iTag)
C..Statistics and memory checks
      DO I=1,I_BLOCK
        IF(IP_XM.EQ.I_AD1(I)) THEN
          II=I
          GOTO 100
        ENDIF
      ENDDO
      write (6,*) ' N_FREEM | Could not find pointer ',IP_XM
      CALL N_PRINT_MEMORY
      CALL Stop_All('N_FREEM',' NO SUCH POINTER')
  100 CONTINUE
      LENGTH=I_LEN(II)
      ERROR=ZFIRST.NE.XM(0).OR.ZLAST.NE.XM(LENGTH+1)
      IF(ERROR) CALL N_MEMORY_CHECK
      IF(MEMCHK) CALL N_MEMORY_CHECK
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
        write (6,*) ' N_FREEM | Error in deallocation of memory '
        CALL Stop_All('N_FREEM',' ')
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
      SUBROUTINE N_MEMORY_CHECK
C     ==--------------------------------------------------------------==
      use legacy_memory_data
      IMPLICIT NONE
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
     *  WRITE(6,'(A)') ' N_MEMORY_CHECK| ALL ARRAY BOUNDARIES ARE OK '
!        IF(PARENT.AND.MEMCHK)
!     *  WRITE(6,'(A)') ' N_MEMORY_CHECK| ALL ARRAY BOUNDARIES ARE OK '
      ELSE
        CALL N_PRINT_MEMORY
        CALL Stop_All('N_MEMORY_CHECK',' CORRUPTED MEMORY')
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END

C     ==================================================================
      SUBROUTINE N_PRINT_MEMORY
C     ==--------------------------------------------------------------==
      use legacy_memory_data
      use util_mod, only: NECI_ICOPY
#ifdef PARALLEL
         use Parallel, only: me => iProcIndex, nProc => nProcessors
#endif
      IMPLICIT NONE
#ifdef POINTER8
      INTEGER*8 I_SORT(MEM_LEN)
#else
      INTEGER   I_SORT(MEM_LEN)
#endif      
      INTEGER   I,II,IC1,IC2
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
      CALL NECI_I8COPY(I_BLOCK,I_AD1,1,I_SORT,1)
#else
      CALL NECI_ICOPY(I_BLOCK,I_AD1,1,I_SORT,1)
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
        II=I
!NDEX(I)
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

      SUBROUTINE NECI_I8COPY(N,A,IA,B,IB)

         IMPLICIT NONE
!        Arguments
         INTEGER   N,IA,IB
         INTEGER*8 A(IA*N),B(IB*N)
!        Variables
         INTEGER I,IAX,IBX

         DO I=1,N
           IAX=(I-1)*IA + 1
           IBX=(I-1)*IB + 1
           B(IBX) = A(IAX)
         ENDDO

         RETURN
      END SUBROUTINE NECI_I8COPY
