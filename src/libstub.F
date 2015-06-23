! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
C==========================================================
C.. Routines needed for the initialisation of the Fourier 
C.. Method of evaluating the Coulomb integrals
C  ========================================================              
      SUBROUTINE INITFOU(NMESHX,CK,COEFF,NMAX,
     &     A,TALPHA,ALPHA,OMEGA,ZIA)
         call stop_all('initfou',
     &                'FFTW library not installed for CPMD NECI')
      RETURN 
      END
C =========================================================
      SUBROUTINE GEN_CK_FFT(N,DIST,COEFF,A,TALPHA,ALPHA,OMEGA)
         call stop_all('gen_ck_fft',
     &              'FFTW library not installed for CPMD NECI')
      RETURN
      END
C ============================================================
      SUBROUTINE GEN_ZIA(KMAX,NMAX,ZIA)
         call stop_all('gen_zia',
     &              'FFTW library not installed for CPMD NECI')
      RETURN
      END
      complex FUNCTION ZIO(K)
         call stop_all('zio',
     &         'FFTW library not installed for CPMD NECI')
         zio = 0
      RETURN
      END
C==========================================================
C.. Routines needed for the initialisation of the Fourier 
C.. Method of evaluating the Coulomb integrals
C  ========================================================              
      SUBROUTINE INITFOU2D(NMESHX,CK,COEFF,NMAX,
     &     A,TALPHA,ALPHA,OMEGA)
         call stop_all('initfou2d',
     &            'FFTW library not installed for CPMD NECI')
      RETURN 
      END
C =========================================================
      SUBROUTINE GEN_CK_FFT2D(N,DIST,COEFF,A,TALPHA,ALPHA,OMEGA)
         call stop_all('gen_ck_fft2d',
     &               'FFTW library not installed for CPMD NECI')
      RETURN
      END
