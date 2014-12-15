C==========================================================
C.. Routines needed for the initialisation of the Fourier 
C.. Method of evaluating the Coulomb integrals
C  ========================================================              
      SUBROUTINE INITFOU(NMESHX,CK,COEFF,NMAX,
     &     A,TALPHA,ALPHA,OMEGA,ZIA)
         STOP 'FFTW library not installed for CPMD NECI'
      RETURN 
      END
C =========================================================
      SUBROUTINE GEN_CK_FFT(N,DIST,COEFF,A,TALPHA,ALPHA,OMEGA)
         STOP 'FFTW library not installed for CPMD NECI'
      RETURN
      END
C ============================================================
      SUBROUTINE GEN_ZIA(KMAX,NMAX,ZIA)
         STOP 'FFTW library not installed for CPMD NECI'
      RETURN
      END
      complex FUNCTION ZIO(K)
         STOP 'FFTW library not installed for CPMD NECI'
         zio = 0
      RETURN
      END
C==========================================================
C.. Routines needed for the initialisation of the Fourier 
C.. Method of evaluating the Coulomb integrals
C  ========================================================              
      SUBROUTINE INITFOU2D(NMESHX,CK,COEFF,NMAX,
     &     A,TALPHA,ALPHA,OMEGA)
         STOP 'FFTW library not installed for CPMD NECI'
      RETURN 
      END
C =========================================================
      SUBROUTINE GEN_CK_FFT2D(N,DIST,COEFF,A,TALPHA,ALPHA,OMEGA)
         STOP 'FFTW library not installed for CPMD NECI'
      RETURN
      END
