! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
      SUBROUTINE MY_HPSI(NDET,NEVAL,NROW,LAB,HAMIL,CK,CKN,
     &                      TLargest)
      use global_utilities
      use constants, only : dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      integer :: NEVAL,NDET,I,J,K,L,IBEG
      real(dp) :: HAMIL(*),CKN(NDET,NEVAL)
      INTEGER LAB(*),NROW(NDET)
      real(dp) :: CK(NDET,NEVAL),AUX
      LOGICAL :: TLargest
      type(timer), save :: proc_timer
C ==-----------------------------------------------------------------==
      proc_timer%timer_name=' MY_HPSI  '
      call set_timer(proc_timer)
C ==-----------------------------------------------------------------==
C..Run over rows
      CKN=0.0_dp
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
          CALL DSCAL(NEVAL*NDET,-1.0_dp,CKN,1)
      ENDIF
c      DO I=1,NEVAL
c        DO J=1,NDET
c          CKN(J,I)=-1.0_dp*CKN(J,I)
c        ENDDO
c      ENDDO
C ==-----------------------------------------------------------------==
      call halt_timer(proc_timer)
C ==-----------------------------------------------------------------==
      RETURN
      END
C ==-----------------------------------------------------------------==
      SUBROUTINE TST_MATMUL(NDET,NEVAL,HTEMP,CK,CKN)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only : dp
      IMPLICIT NONE
      integer :: NDET,NEVAL
      real(dp) HTEMP(NDET,NDET),CK(NDET,NEVAL),CKN(NDET,NEVAL)
C..
      CALL DGEMM('N','N',NDET,NEVAL,NDET,1.0_dp,HTEMP,NDET,
     &    CK,NDET,1.0_dp,CKN,NDET)
C      DO I=1,NDET
C         DO K=1,NEVAL
C            WRITE(6,*) CKN(I,K)
C         ENDDO
C      ENDDO 
C      DO I=1,NDET
C        DO J=1,NDET
C          IF(HTEMP(I,J).NE.0.0_dp) THEN
C            WRITE(6,*) HTEMP(I,J),I,J
C          ENDIF
C        ENDDO
C      ENDDO
      RETURN
      END
