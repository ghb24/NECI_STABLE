      SUBROUTINE NECI_ICOPY(N,A,IA,B,IB)
         != Copy elements from integer array A to B.
         != Simple version of BLAS routine ICOPY.
         != In:
         !=    N: number of elements in A.
         !=    A: vector to be copied.
         !=    IA: increment between elements to be copied in A.  
         !=        IA=1 for continuous data blocks.
         !=    IB: increment between elements to be copied to in B.  
         !=        IB=1 for continuous data blocks.
         != Out:
         !=    B: result vector.
         IMPLICIT NONE
!        Arguments
         INTEGER N,IA,IB
         INTEGER A(IA*N),B(IB*N)
!        Variables
         INTEGER I,IAX,IBX
     
         DO I=1,N
           IAX=(I-1)*IA + 1
           IBX=(I-1)*IB + 1
           B(IBX) = A(IAX)
         ENDDO
     
         RETURN
      END SUBROUTINE NECI_ICOPY

!   A fudge of the unportrable routine DGEMA to add two matrices
      SUBROUTINE DGEMA(TRANSA,TRANSB,M,N,ALPHA,A,LDA,BETA,B,LDB,C,LDC)
         IMPLICIT NONE
         CHARACTER*1 TRANSA,TRANSB
         INTEGER*4 M,N,LDA,LDB,LDC
         REAL*8 ALPHA,A(M,N),B(M,N),BETA,C(M,N)
         INTEGER I,J
         IF(LDA.NE.M) STOP 'LDA NE M not supported'
         IF(LDB.NE.M) STOP 'LDB NE M not supported'
         IF(LDC.NE.M) STOP 'LDC NE M not supported'
         IF(TRANSA.NE.'N'.AND.TRANSA.NE.'n')
     &      STOP 'TRANSA not "N" not supported'
         IF(TRANSB.NE.'N'.AND.TRANSB.NE.'n')
     &      STOP 'TRANSB not "N" not supported'
         DO I=1,M
            DO J=1,N
               C(I,J)=ALPHA*A(I,J)+BETA*B(I,J)
            ENDDO
         ENDDO
         RETURN
      END


!  Following the same interface as Numerical Recipes' RAN2, but using the RANLUX generator (lux level 3).
!  If ISEED.NE.0, it is used as a seed, and reset to 0.
!  Generates a 32-bit float (REAL*4), and stores this in a real*8
      REAL*8 FUNCTION RAN2(ISEED)
         IMPLICIT NONE
         INTEGER ISEED
!         REAL*8 r
         IF(ISEED.NE.0) THEN
!  Init RANLUX
            IF(ISEED.LT.0) ISEED=-ISEED
            CALL RLUXGO(3,ISEED,0,0)
            ISEED=0
         ENDIF
         CALL RANLUX(RAN2,1)
!         RAN2=r
         RETURN
      END FUNCTION RAN2

      integer function record_length(bytes)
         ! Some compilers use record lengths in units of bytes.
         ! Some compilers use record lengths in units of words.
         ! This is an utter *pain* for reading unformatted files,
         ! where you must specify the record length.
         !
         ! In:
         !    bytes: number of bytes in record type of interest (should
         !    be a multiple of 4).
         !
         ! Returns:
         !    record_length: size of record in units of the compiler's
         !    choice.
         integer bytes
         inquire(iolength=record_length) bytes
         record_length = (bytes/4)*record_length 
      end function record_length
