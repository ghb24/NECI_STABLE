      SUBROUTINE MY_HPSI(NDET,ICMAX,NEVAL,NROW,LAB,HAMIL,CK,CKN,
     &                      TLargest)
      use global_utilities
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HAMIL(*),CKN(NDET,NEVAL)
      INTEGER LAB(*),NROW(NDET)
      DIMENSION CK(NDET,NEVAL)
      LOGICAL :: TLargest
      type(timer), save :: proc_timer
C ==-----------------------------------------------------------------==
      proc_timer%timer_name=' MY_HPSI  '
      call set_timer(proc_timer)
C ==-----------------------------------------------------------------==
C..Run over rows
      CKN=0.d0
      DO I=1,NDET
C..Run over columns
         IF(I.EQ.1) THEN 
            IBEG=0
         ELSE
            IBEG=IBEG+NROW(I-1)
         ENDIF
         DO K=1,NEVAL
C..
            J=LAB(IBEG+1)
            CKN(I,K)=CKN(I,K)+HAMIL(IBEG+1)*CK(J,K)
            DO L=2,NROW(I)
               J=LAB(IBEG+L)
               AUX=HAMIL(IBEG+L)
               CKN(I,K)=CKN(I,K)+AUX*CK(J,K)
               CKN(J,K)=CKN(J,K)+AUX*CK(I,K)
            ENDDO
C..
         ENDDO
      ENDDO
C..Need to do the following mult. by -1 so that eigenvalues are calc. in
C..ascending order
      IF(.NOT.TLargest) THEN
          CALL DSCAL(NEVAL*NDET,-1.D0,CKN,1)
      ENDIF
c      DO I=1,NEVAL
c        DO J=1,NDET
c          CKN(J,I)=-1.D0*CKN(J,I)
c        ENDDO
c      ENDDO
C ==-----------------------------------------------------------------==
      call halt_timer(proc_timer)
C ==-----------------------------------------------------------------==
      RETURN
      END
C ==-----------------------------------------------------------------==
      SUBROUTINE TST_MATMUL(NDET,NEVAL,HTEMP,CK,CKN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 HTEMP(NDET,NDET),CK(NDET,NEVAL),CKN(NDET,NEVAL)
C..
      CALL DGEMM('N','N',NDET,NEVAL,NDET,1.D0,HTEMP,NDET,
     &    CK,NDET,1.D0,CKN,NDET)
C      DO I=1,NDET
C         DO K=1,NEVAL
C            WRITE(6,*) CKN(I,K)
C         ENDDO
C      ENDDO 
C      DO I=1,NDET
C        DO J=1,NDET
C          IF(HTEMP(I,J).NE.0.D0) THEN
C            WRITE(6,*) HTEMP(I,J),I,J
C          ENDIF
C        ENDDO
C      ENDDO
      RETURN
      END
