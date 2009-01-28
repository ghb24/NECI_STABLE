SUBROUTINE OrthoNormx(n,m,a)

   REAL*8 :: work(n),tau(n),a(m,n)
   INTEGER :: k,n,m,lda,lwork,info
   REAL*8 , ALLOCATABLE :: aTa(:,:)

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
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Q_FACTOR - Warning:'
        write ( *, '(a,i8)' ) '  DGEQRF returned nonzero INFO = ', info
        stop
    end if

!  Construct Q explicitly.
    
    call dorgqr ( m, n, k, a, lda, tau, work, lwork, info )

    if ( info /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Q_FACTOR - Warning:'
        write ( *, '(a,i8)' ) '  DORGQR returned nonzero INFO = ', info
        stop
    end if


!  Perform orthonormality test.
!
!    allocate ( aTa(n,n) )
!
!    aTa = matmul ( transpose ( a ), a )
!
!!   aTa should be the identity...
!
!    do i=1,n
!        do j=1,n
!            WRITE(6,"(F15.7)",advance='no') aTa(j,i)
!        enddo
!        WRITE(6,*) ""
!    enddo
!
!
!    deallocate ( aTa )

END SUBROUTINE OrthoNormx

