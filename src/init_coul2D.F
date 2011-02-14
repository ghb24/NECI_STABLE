C==========================================================
C.. Routines needed for the initialisation of the Fourier 
C.. Method of evaluating the Coulomb integrals
C  ========================================================              
      SUBROUTINE INITFOU2D(NMESHX,CK,NMAX,
     &     A,TALPHA,ALPHA,OMEGA,ZIA)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CK(NMESHX,NMESHX,NMESHX)
      COMPLEX*16 ZIA(-NMESHX/2:NMESHX/2,NMAX,NMAX) 
      REAL*8 A(3)
      real*4 t1,t2,t3,t(2),etime
      LOGICAL TALPHA
      DATA IFIRST/0/
      SAVE IFIRST 
C..
      IF(IFIRST.EQ.1) RETURN 
      IFIRST=1
C..      
      T1 = etime(t)
      CALL GEN_CK_FFT2D(NMESHX,CK,A,TALPHA,ALPHA,OMEGA)
      CALL GEN_ZIA(NMESHX,NMAX,ZIA)
      T2 = etime(t)
      T3=(T2-T1)
      WRITE(6,*) 'V0=',CK(NMESHX/2+1,NMESHX/2+1,NMESHX/2+1)
      WRITE(6,*) ' TIME FOR INITIALISATION:',T3/1000.
C..      
      RETURN 
      END
C =========================================================
      SUBROUTINE GEN_CK_FFT2D(N,DIST,A,TALPHA,ALPHA,OMEGA)
      IMPLICIT REAL*8(A-H,O-Z)
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
      COMPLEX*16 DIST(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1)
      INTEGER*8 PLAN
      REAL*8 A(3) 
      LOGICAL TALPHA
      COMPLEX*16 D(-N/2:N/2-1,-N/2:N/2-1)
C..
      PI=ACOS(-1.D0)
      ALPHA2=ALPHA*ALPHA
C..
      HSTEPX=2.D0*A(1)/N
      HSTEPY=2.D0*A(2)/N
      HSTEPZ=0
    
C..
      DO I=-N/2,N/2-1
        DO J=-N/2,N/2-1
c          DO K=-N/2,N/2-1
            K=0
            X=DFLOAT(I)*HSTEPX
            Y=DFLOAT(J)*HSTEPY
            Z=DFLOAT(K)*HSTEPZ
            AUX=X*X+Y*Y+Z*Z
            IF(AUX.NE.0.D0) THEN
              IF(TALPHA) THEN 
                R=SQRT(AUX) 
                AUX=1.D0/R*(-1)**(I+J)*DERF(R/ALPHA)
              ELSE
                AUX=1.D0/SQRT(AUX)*(-1)**(I+J)
              ENDIF
            ELSE
              IF(TALPHA) THEN 
                AUX=(2.D0/SQRT(PI)/ALPHA)*(-1)**(I+J)
              ELSE
                AUX=DFLOAT(N/2)*(-1)**(I+J)
              ENDIF
            ENDIF
            D(I,J)=DCMPLX(AUX,0.D0)
c            write(16,'(3i3,2e19.8)') i,j,k,d(i,j)

c          ENDDO
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
C      CALL DFFTW_PLAN_DFT_3D(PLAN,N1,N2,N3,DIST,
C     &  DIST,
C     &  FFTW_FORWARD,FFTW_ESTIMATE)
C      CALL DFFTW_EXECUTE(PLAN)
C      CALL DFFTW_DESTROY_PLAN(PLAN)
c      CALL DFFTW_PLAN_DFT_2D(PLAN,N1,N2,
c     &  DIST(-N/2:N/2-1,-N/2:N/2-1,0),
c     &  DIST(-N/2:N/2-1,-N/2:N/2-1,0),
c     &  FFTW_FORWARD,FFTW_ESTIMATE)
      CALL DFFTW_PLAN_DFT_2D(PLAN,N1,N2,
     &  D,
     &  DIST(-N/2,-N/2,0),
     &  FFTW_FORWARD,FFTW_ESTIMATE)
c      call dfftw_print_plan(plan)
      CALL DFFTW_EXECUTE(PLAN)
      CALL DFFTW_DESTROY_PLAN(PLAN)
#endif
#endif
C..Shift origin and normalise 
      DO I=-N/2,N/2-1
        DO J=-N/2,N/2-1
C          DO K=-N/2,N/2-1
            K=0
            DIST(I,J,0)=DIST(I,J,0)*(-1)**(i+j+k)/dfloat(n1*n2)
C/dfloat(n1*n2)
C          enddo
        enddo
      enddo      
C..the short-range correction for alpha 
      IF(TALPHA) THEN 

       OMEGAP=4.D0*OMEGA  
        DO I=-N/2,N/2-1
          DO J=-N/2,N/2-1
C           DO K=-N/2,N/2-1
            K=0
              GX=PI*I/A(1)
              GY=PI*J/A(2)
              GZ=0
C..PI*K/A(3)
              G=SQRT(GX*GX+GY*GY+GZ*GZ)
              IF(G.NE.0.D0) THEN 
c sum=2 Pi * erf(a G/2)/G/omegap
                SUM=2.D0*PI*DERF(ALPHA*G/2)/OMEGAP/G
C                AUX=EXP(-G2*ALPHA2/4.D0)
C                SUM=4.D0*PI*(1.D0-AUX)/OMEGAP/G2
              ELSE
                 SUM=ALPHA/(SQRT(PI)*OMEGAP)
C                SUM=PI*ALPHA2/OMEGAP 
              ENDIF
              DIST(I,J,K)=DIST(I,J,K)+DCMPLX(SUM,0.D0)
C            ENDDO
          ENDDO
        ENDDO
      ENDIF
C..Inverse transform 
c      CALL ZFFT3D( 1, N1, N2, N3, DIST, LA1, LA2, COEFF )
C.test fft
C       X=  0.2
C       Y=  0.9
C       Z=  0.0
C       UIMAG=DCMPLX(0.D0,1.D0)
C       CSUM=DCMPLX(0.D0,0.D0)
C       DO I=-N/2,N/2-1
C         DO J=-N/2,N/2-1
C             k=0
C             GX=PI*I
C             GY=PI*J
C             GZ=PI*K
C             CDOT=DCMPLX(0.D0,GX*X+GY*Y+GZ*Z)
C             CAUX= EXP(CDOT)*DIST(I,J,K)
C             CSUM=CSUM+CAUX
C               write(6,'(3i3,2e19.8)') i,j,k,dist(i,j,k)
C         enddo
C       enddo
C       sum=sqrt(real(csum)**2+dimag(csum)**2)
C       write(6,*) X,csum,1.D0/sqrt(x*x*a(1)**2+y*y*a(2)**2)
C       stop
      RETURN
      END
