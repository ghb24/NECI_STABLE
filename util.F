#if !(defined(CRAY))
      SUBROUTINE ICOPY(N,A,IA,B,IB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,IA,IB
      INTEGER A(IA*N),B(IB*N)
C     Variables
      INTEGER I,IAX,IBX
C     ==--------------------------------------------------------------==
      DO I=1,N
        IAX=(I-1)*IA + 1
        IBX=(I-1)*IB + 1
        B(IBX) = A(IAX)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#ifdef POINTER8
      SUBROUTINE I8COPY(N,A,IA,B,IB)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER   N,IA,IB
      INTEGER*8 A(IA*N),B(IB*N)
C     Variables
      INTEGER I,IAX,IBX
C     ==--------------------------------------------------------------==
      DO I=1,N
        IAX=(I-1)*IA + 1
        IBX=(I-1)*IB + 1
        B(IBX) = A(IAX)
      ENDDO
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
#endif
#endif

C.. A fudge of the unportrable routine DGEMA to add two matrices
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


C.. RAN2 - generate a decent random number from Numerical Recipes
C.. 15th Aug '05.  It would appear this is not Numerical Recipes'
C.. RAN2, but an implementation of something different with the same
C.. interface.

C.. Initialize with a seed <0

C.. Modified this to use REAL*8s.
      REAL*8 FUNCTION ORAN2(IDUM)
      IMPLICIT NONE
      INTEGER IR(97),IFF,J,M,IY,IA,IC,IDUM
      PARAMETER (M=714025,IA=1366,IC=150889)
      REAL*8 RM
      DATA IFF /0/
      RM=1.4005112D-6
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1) STOP
      IY=IR(J)
      ORAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END

C.. Numerical Recipes RAN2

      REAL*8 FUNCTION OORAN2(idum)
        IMPLICIT NONE
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL*8 AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.D0/IM1,
     &          IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,
     &          IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,
     &          EPS=1.2D-7,RNMX=1.D0-EPS)
        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
!        INTEGER C
        DATA idum2/123456789/, iv/NTAB*0/,iy/0/
!        SAVE C
!        C=C+1
!        IF(MOD(C,10000).EQ.0) WRITE(6,*) "PRNG: ",C
        if (idum.le.0) then
           idum=max(-idum,1)
           idum2=idum
           do j=NTAB+8,1,-1
              k=idum/IQ1
              idum=IA1*(idum-k*IQ1)-k*IR1
              if (idum.lt.0) idum=idum+IM1
              if (j.le.NTAB) iv(j)=idum
           enddo
           iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1) iy=iy+IMM1
        ooran2=min(AM*iy,RNMX)
        return
      END

!  Following the same interface as RAN2, but using the RANLUX generator (lux level 3).
!  If ISEED.NE.0, it is used as a seed, and reset to 0.
!  Generates a 32-bit float (REAL*4), and stores this in a real*8
      REAL*8 FUNCTION RAN2(ISEED)
         IMPLICIT NONE
         INTEGER ISEED
         REAL*8 r
         IF(ISEED.NE.0) THEN
!  Init RANLUX
            IF(ISEED.LT.0) ISEED=-ISEED
            CALL RLUXGO(3,ISEED,0,0)
            ISEED=0
         ENDIF
         CALL RANLUX(r,1)
         RAN2=r
         RETURN
      END


      REAL*8 FUNCTION FACTRL(N)
          != Return the factorial of N, N!.
          != This is not done in the most efficient way possible (i.e. use
          != with care if N is large or if called many times!).
          != The reader is referred to
          != http://www.luschny.de/math/factorial/FastFactorialFunctions.htm
          != if a more efficient procedure is desired.
          implicit none
          integer :: i,N
          FACTRL=1.d0
          do i=2,N
              FACTRL=FACTRL*i
          end do
      end function FACTRL
