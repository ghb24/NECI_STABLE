! Lowdin Orthoganalize
! for any non-singular R, let S=R RT
! P = S^(-1/2) R is orthogonal.
! MAT is NxN and is returned as an orthogal matrix
! R1 and R2 are NxN workspaces
      SUBROUTINE LOWDIN_ORTH(MAT,N,R1,R2,WORK)
         use constants, only: dp,sp
         use HElem
         IMPLICIT NONE
         INTEGER N
         HElement_t(dp) MAT(N,N),R1(N,N),R2(N,N)
         real(dp) L(N),LL,RWORK(3*N)
         HElement_t(dp) WORK(3*N)
         INTEGER I,J
         integer(sp) info
         character(*), parameter :: this_routine = 'LOWDIN_ORTH'
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
            call stop_all(this_routine, 'Error in LOWDIN_ORTH.')
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
 
      SUBROUTINE GRAMSCHMIDT_NECI(MAT,LEN)
! MAT(IELEMENT,IVECTOR)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER LEN
         HElement_t(dp) MAT(LEN,LEN),DOT
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
      END SUBROUTINE GRAMSCHMIDT_NECI


