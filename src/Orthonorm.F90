! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
SUBROUTINE OrthoNormx(n,m,a)
  use constants, only: dp,sp
   implicit none
   INTEGER :: i, j, k,n,m,lda,lwork
   real(dp) :: work(n),tau(n),a(m,n)
   INTEGER(sp) info
   real(dp) , ALLOCATABLE :: aTa(:,:)

!Input the number of vectors, n, the dimensionality of the space, m, and the matrix of vectors, a(m,n). 
!a is returned as n orthonormal vectors.

!    DGEQRF implicitly computes the QR factorization of an M by N 
!    matrix A:

!      A(MxN) = Q(MxK) * R(KxN)

!    where K = min ( M, N ).  For our purposes, it should always
!    be the case that N < M, so that K = N.

!    DORGQR explicitly forms the Q matrix.
!First, compute the QR factorization.
   k = n
   lda = m
   lwork = n

   call dgeqrf ( m, n, a, lda, tau, work, lwork, info )
   if ( info /= 0 ) then
        write (6, '(a)' ) ' '
        write (6, '(a)' ) 'Q_FACTOR - Warning:'
        write (6, '(a,i8)' ) '  DGEQRF returned nonzero INFO = ', info
        stop
    end if

!  Construct Q explicitly.
    
    call dorgqr ( m, n, k, a, lda, tau, work, lwork, info )

    if ( info /= 0 ) then
        write (6, '(a)' ) ' '
        write (6, '(a)' ) 'Q_FACTOR - Warning:'
        write (6, '(a,i8)' ) '  DORGQR returned nonzero INFO = ', info
        stop
    end if


!  Perform orthonormality test.

    allocate ( aTa(n,n) )

    aTa = matmul ( transpose ( a ), a )

!   aTa should be the identity...

    do i=1,n
        do j=1,n
            WRITE(6,"(F15.7)",advance='no') aTa(j,i)
        enddo
        WRITE(6,*) ""
    enddo


    deallocate ( aTa )

END SUBROUTINE OrthoNormx

!Lowdin Orthoganalize
!for any non-singular R, let S=R RT
!P = S^(-1/2) R is orthogonal.
!MAT is NxN and is returned as an orthogal matrix
!R1 and R2 are NxN workspaces
      SUBROUTINE LOWDIN_ORTH(MAT,N,R1,R2,WORK)
         use constants, only: dp,sp
         use HElem
         IMPLICIT NONE
         INTEGER N
         HElement_t MAT(N,N),R1(N,N),R2(N,N)
         real(dp) L(N),LL,RWORK(3*N)
         HElement_t WORK(3*N)
         INTEGER I,J
         integer(sp) info
!R=MAT
!S= R1=1.0_dp * R * RT + 0.0_dp*R1
         IF(HElement_t_size.EQ.1) THEN
            CALL DGEMM('N','T',N,N,N,1.0_dp,MAT,N,MAT,N,0.0_dp,R1,N)
         ELSE
            CALL ZGEMM('N','C',N,N,N,(1.0_dp,0.0_dp),MAT,N,MAT,N,(0.0_dp,0.0_dp),R1,N)
         ENDIF
!         CALL Write_HEMatrix("R",N,N,R1)
! Diagonalize S=R1 into eigenvectors U=R1 and eigenvalues L
         IF(HElement_t_size.EQ.1) THEN
            CALL DSYEV('V','U',N,R1,N,L,WORK,N*3,INFO)
         ELSE
            CALL ZHEEV('V','U',N,R1,N,L,WORK,N*3,RWORK,INFO)
         ENDIF
! eigenvector 1 is given by R1(I,1)
         IF(INFO.NE.0) THEN
            WRITE(6,*) "INFO=",INFO," on diag in LOWDIN_ORTH. Stopping"
            STOP 'Error in LOWDIN_ORTH.'
         ENDIF
! Calculate P = S^(-1/2) R = U L^(-1/2) UT R.  U=R1
! First let R2=U R. U=R1.  R=MAT
! R2=1.0_dp * R1 * MAT + 0.0_dp*R1
         IF(HElement_t_size.EQ.1) THEN
            CALL DGEMM('T','N',N,N,N,1.0_dp,R1,N,MAT,N,0.0_dp,R2,N)
         ELSE
            CALL ZGEMM('C','N',N,N,N,(1.0_dp,0.0_dp),R1,N,MAT,N,(0.0_dp,0.0_dp),R2,N)
         ENDIF
! Now let R2=L^(-1/2) (U R) = L^(-1/2) R2
! row I is multiplied by (L(I))^(-1/2)
         DO I=1,N
            LL=L(I)**(-0.5_dp)
            DO J=1,N
               R2(I,J)=R2(I,J)*(LL)
            ENDDO
         ENDDO
! Now let MAT = P = U (L^(-1/2) UT R) = U R2.  U=R1
! MAT=1.0_dp * U * R2 + 0.0_dp*MAT
         IF(HElement_t_size.EQ.1) THEN
            CALL DGEMM('N','N',N,N,N,1.0_dp,R1,N,R2,N,0.0_dp,MAT,N)
         ELSE
            CALL ZGEMM('N','N',N,N,N,(1.0_dp,0.0_dp),R1,N,R2,N,(0.0_dp,0.0_dp),MAT,N)
         ENDIF
! and we should be done, with an orthoganal matrix in MAT
      END SUBROUTINE LOWDIN_ORTH                      
 
      SUBROUTINE GRAMSCHMIDT(MAT,LEN)
! MAT(IELEMENT,IVECTOR)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER LEN
         HElement_t MAT(LEN,LEN),DOT
         real(dp) NORM,SNORM
         INTEGER I,J,K
         DO I=1,LEN
! First dot with all lower vectors, and remove their components
            DO J=1,I-1
               DOT=0.0_dp
               NORM=0.0_dp
               DO K=1,LEN
#ifdef __CMPLX
                  DOT=DOT+conjg(MAT(K,J))*MAT(K,I)
#else
                  DOT=DOT+(MAT(K,J))*MAT(K,I)
#endif
               ENDDO
               DO K=1,LEN
                  MAT(K,I)=MAT(K,I)-MAT(K,J)*DOT
               ENDDO 
            ENDDO
            NORM=0.0_dp
            DO K=1,LEN
               NORM=NORM+abs(MAT(K,I))**2
            ENDDO        
            SNORM=SQRT(NORM)    
!            WRITE(6,*) NORM
            DO K=1,LEN
!               WRITE(6,*) MAT(K,I),MAT(K,I)/SNORM
               MAT(K,I)=MAT(K,I)/(SNORM)
            ENDDO
         ENDDO
         RETURN
      END SUBROUTINE GRAMSCHMIDT


