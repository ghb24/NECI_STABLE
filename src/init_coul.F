C==========================================================
C.. Routines needed for the initialisation of the Fourier 
C.. Method of evaluating the Coulomb integrals
C  ========================================================              
      SUBROUTINE INITFOU(NMESHX,CK,NMAX,
     &     A,TALPHA,ALPHA,OMEGA,ZIA)
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp,int64,sp
      use util_mod, only: neci_etime
      IMPLICIT NONE
      integer :: NMESHX,NMAX
      complex(dp) CK(NMESHX,NMESHX,NMESHX)
      complex(dp) ZIA(-NMESHX/2:NMESHX/2,NMAX,NMAX) 
      real(dp) :: ALPHA,OMEGA
      real(dp) A(3),t3
      real(dp) t(2),t1,t2
      LOGICAL TALPHA
      INTEGER , SAVE :: IFIRST=0
C..
      IF(IFIRST.EQ.1) RETURN 
      IFIRST=1
C..      
      T1 = neci_etime(t)
      CALL GEN_CK_FFT(NMESHX,CK,A,TALPHA,ALPHA,OMEGA)
      CALL GEN_ZIA(NMESHX,NMAX,ZIA)
      T2 = neci_etime(t)
      T3=(T2-T1)
      WRITE(6,*) 'V0=',CK(NMESHX/2,NMESHX/2,NMESHX/2)
      WRITE(6,*) ' TIME FOR INITIALISATION:',T3/1000.
C..      
      RETURN 
      END
C =========================================================
      SUBROUTINE GEN_CK_FFT(N,DIST,A,TALPHA,ALPHA,OMEGA)
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp,int64,pi

      IMPLICIT NONE
#if !defined(__alpha) || !defined(__SGI) 
      !INTEGER FFTW_R2HC
      !PARAMETER (FFTW_R2HC=0)
      !INTEGER FFTW_HC2R
      !PARAMETER (FFTW_HC2R=1)
      !INTEGER FFTW_DHT
      !PARAMETER (FFTW_DHT=2)
      !INTEGER FFTW_REDFT00
      !PARAMETER (FFTW_REDFT00=3)
      !INTEGER FFTW_REDFT01
      !PARAMETER (FFTW_REDFT01=4)
      !INTEGER FFTW_REDFT10
      !PARAMETER (FFTW_REDFT10=5)
      !INTEGER FFTW_REDFT11
      !PARAMETER (FFTW_REDFT11=6)
      !INTEGER FFTW_RODFT00
      !PARAMETER (FFTW_RODFT00=7)
      !INTEGER FFTW_RODFT01
      !PARAMETER (FFTW_RODFT01=8)
      !INTEGER FFTW_RODFT10
      !PARAMETER (FFTW_RODFT10=9)
      !INTEGER FFTW_RODFT11
      !PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      !INTEGER FFTW_BACKWARD
      !PARAMETER (FFTW_BACKWARD=+1)
      !INTEGER FFTW_MEASURE
      !PARAMETER (FFTW_MEASURE=0)
      !INTEGER FFTW_DESTROY_INPUT
      !PARAMETER (FFTW_DESTROY_INPUT=1)
      !INTEGER FFTW_UNALIGNED
      !PARAMETER (FFTW_UNALIGNED=2)
      !INTEGER FFTW_CONSERVE_MEMORY
      !PARAMETER (FFTW_CONSERVE_MEMORY=4)
      !INTEGER FFTW_EXHAUSTIVE
      !PARAMETER (FFTW_EXHAUSTIVE=8)
      !INTEGER FFTW_PRESERVE_INPUT
      !PARAMETER (FFTW_PRESERVE_INPUT=16)
      !INTEGER FFTW_PATIENT
      !PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      !INTEGER FFTW_ESTIMATE_PATIENT
      !PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      !INTEGER FFTW_BELIEVE_PCOST
      !PARAMETER (FFTW_BELIEVE_PCOST=256)
      !INTEGER FFTW_DFT_R2HC_ICKY
      !PARAMETER (FFTW_DFT_R2HC_ICKY=512)
      !INTEGER FFTW_NONTHREADED_ICKY
      !PARAMETER (FFTW_NONTHREADED_ICKY=1024)
      !INTEGER FFTW_NO_BUFFERING
      !PARAMETER (FFTW_NO_BUFFERING=2048)
      !INTEGER FFTW_NO_INDIRECT_OP
      !PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      !INTEGER FFTW_ALLOW_LARGE_GENERIC
      !PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      !INTEGER FFTW_NO_RANK_SPLITS
      !PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      !INTEGER FFTW_NO_VRANK_SPLITS
      !PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      !INTEGER FFTW_NO_VRECURSE
      !PARAMETER (FFTW_NO_VRECURSE=65536)
      !INTEGER FFTW_NO_SIMD
      !PARAMETER (FFTW_NO_SIMD=131072)
#endif
      integer :: N,I,J,K,N1,N2,N3
      complex(dp) DIST(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1)
      integer(int64) PLAN
      real(dp) A(3),ALPHA2,ALPHA,HSTEPX,HSTEPY,HSTEPZ,DERF
      real(dp) :: X,Y,Z,AUX,SUM,G2,GX,GY,GZ,OMEGA,OMEGAP,R
      LOGICAL TALPHA
      character(len=*), parameter :: t_r="GEN_CK_FFT"
C..
!      PI=ACOS(-1.0_dp)
      ALPHA2=ALPHA*ALPHA
C..
      HSTEPX=2.0_dp*A(1)/N
      HSTEPY=2.0_dp*A(2)/N
      HSTEPZ=2.0_dp*A(3)/N
C..
#if defined(NAGF95) || defined(__GFORTRAN__) || defined(BLUEGENE_HACKS)
      call stop_all(t_r,"No ERF in NAG/GFortran?")
#endif
      DO I=-N/2,N/2-1
        DO J=-N/2,N/2-1
          DO K=-N/2,N/2-1
            X=real(I,dp)*HSTEPX
            Y=real(J,dp)*HSTEPY
            Z=real(K,dp)*HSTEPZ
            AUX=X*X+Y*Y+Z*Z
            IF(AUX.NE.0.0_dp) THEN
              IF(TALPHA) THEN 
                R=SQRT(AUX) 
#if !defined(NAGF95) && !defined(__GFORTRAN__) && !defined(BLUEGENE_HACKS)
                AUX=1.0_dp/R*(-1)**(I+J+K)*DERF(R/ALPHA)
#endif                
              ELSE
                AUX=1.0_dp/SQRT(AUX)*(-1)**(I+J+K)
              ENDIF
            ELSE
              IF(TALPHA) THEN 
                AUX=2.0_dp/SQRT(PI)/ALPHA
              ELSE
                AUX=real(N/2,dp)*(-1)**(I+J+K)
              ENDIF
            ENDIF
            DIST(I,J,K)=CMPLX(AUX,0.0_dp,dp)
c            write(16,'(3i3,2e19.8)') i,j,k,dist(i,j,k)
          ENDDO
        ENDDO
      ENDDO
c      write(6,*) ' inverse transform:'
#ifdef __SGI
C..FFT parameters
      N1=N
      N2=N
      N3=N
      LA1=N1
      LA2=N2 
C..Initialise FFT
      CALL ZFFT3DI( N1, N2, N3, COEFF)
C..FORWARD TRANSFORM 
      CALL ZFFT3D(-1, N1, N2, N3, DIST, LA1, LA2, COEFF )
#elif __alpha 
      N1=N
      N2=N
      N3=N
      LA1=N1
      LA2=N2 
      CALL ZFFT_3D('C','C','F',DIST,DIST,N1,N2,N3,LA1,LA2,1,1,1)
#else
      N1=N
      N2=N
      N3=N
#ifndef DISABLE_FFTW
      CALL DFFTW_PLAN_DFT_3D(PLAN,N1,N2,N3,DIST,
     &  DIST,
     &  FFTW_FORWARD,FFTW_ESTIMATE)
      CALL DFFTW_EXECUTE(PLAN)
      CALL DFFTW_DESTROY_PLAN(PLAN)
#else
      call stop_all("gen_ck_fft","FFTW disabled")
#endif
#endif
C..Shift origin and normalise 
      DO I=-N/2,N/2-1
        DO J=-N/2,N/2-1
          DO K=-N/2,N/2-1
            DIST(I,J,K)=DIST(I,J,K)*(-1)**(i+j+k)/real(n1*n2*n3,dp)
          enddo
        enddo
      enddo      
C..the short-range correction for alpha 
      IF(TALPHA) THEN 
        OMEGAP=8.0_dp*OMEGA  
        DO I=-N/2,N/2-1
          DO J=-N/2,N/2-1
            DO K=-N/2,N/2-1
              GX=PI*I/A(1)
              GY=PI*J/A(2)
              GZ=PI*K/A(3)
              G2=GX*GX+GY*GY+GZ*GZ
              IF(G2.NE.0.0_dp) THEN 
!                AUX=EXP(-G2*ALPHA2/4.0_dp)
                SUM=4.0_dp*PI*(1.0_dp-AUX)/OMEGAP/G2
              ELSE
                SUM=PI*ALPHA2/OMEGAP 
              ENDIF
              DIST(I,J,K)=DIST(I,J,K)+CMPLX(SUM,0.0_dp,dp)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C..Inverse transform 
c      CALL ZFFT3D( 1, N1, N2, N3, DIST, LA1, LA2, COEFF )
C.test fft 
c      X= .2
c      Y= 0.9
c      Z= -0.2
c      UIMAG=DCMPLX(0.0_dp,1.0_dp)
c      CSUM=DCMPLX(0.0_dp,0.0_dp)
c      DO I=-N/2,N/2-1
c        DO J=-N/2,N/2-1
c          DO K=-N/2,N/2-1
c            GX=PI*I
c            GY=PI*J
c            GZ=PI*K
c            CDOT=DCMPLX(0.0_dp,GX*X+GY*Y+GZ*Z)
c            CAUX=ZEXP(CDOT)*DIST(I,J,K)
cc*(-1)**(i+j+k)
c            CSUM=CSUM+CAUX
c            write(6,'(3i3,2e19.8)') i,j,k,dist(i,j,k)/(n1*n2*n3)
c          enddo
c        enddo
c      enddo
c      csum=csum/dfloat(n1*n2*n3)
c      sum=sqrt(real(csum)**2+dimag(csum)**2)
c      write(6,*) csum,1.0_dp/sqrt(x*x+y*y+z*z) 
c      write(6,*) sum
c      st--op
      RETURN
      END
C ============================================================
      SUBROUTINE GEN_ZIA(KMAX,NMAX,ZIA)
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: K,N,M,KMAX,NMAX
      complex(dp) ZIA(-KMAX/2:KMAX/2,NMAX,NMAX),ZIO
      DO K=-KMAX/2,KMAX/2
        DO N=1,NMAX
          DO M=1,NMAX
            ZIA(K,N,M)=0.25_dp*(ZIO(N-M+K)+ZIO(-N+M+K)
     &           -ZIO(N+M+K)-ZIO(-N-M+K))
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      FUNCTION ZIO(K)
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only : dp,pi
      IMPLICIT NONE
      complex(dp) ZIO
!      INTEGER , DATA , SAVE :: IFIRST/0/ 
!      real(dp) :: SAVE PI
      INTEGER :: K
!      IF(IFIRST.EQ.0) THEN 
!        PI=ACOS(-1.0_dp)
!        IFIRST=1
!      ENDIF
      IF(K.EQ.0) THEN 
        ZIO=(1.0_dp,0.0_dp)
      ELSE
        ZIO=CMPLX(0.0_dp,-((-1.0_dp)**K-1.0_dp)/(REAL(K,dp)*PI),dp)
      ENDIF
      RETURN
      END
