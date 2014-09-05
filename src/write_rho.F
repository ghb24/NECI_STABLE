! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
      SUBROUTINE WRITE_RHO(NC,FILENM,RHO,NMESHX,NMESHY,
     &    NMESHZ,ALAT,TMOVIE,TGOPENMOL,RS)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: NMESHX,NMESHY,NMESHZ,I,J,K,NC
      real(dp) :: RHO(NMESHX,NMESHY,NMESHZ),ALAT(3),SUM
      real(dp) :: rs,rs3
      LOGICAL TMOVIE,TGOPENMOL
      CHARACTER FILENM*10,FILENM2*16,FILENM3*14
      OPEN(10,FILE=FILENM,STATUS='UNKNOWN')
      WRITE(10,*) NMESHX,NMESHY,NMESHZ
      IF(TMOVIE) THEN
        FILENM2=FILENM(1:NC)//'.movie'
        OPEN(11,FILE=FILENM2,STATUS='UNKNOWN')
C..Compute mean value of function for plotting function
        SUM =0
        DO K=1,NMESHZ
          DO J=1,NMESHY
          DO I=1,NMESHX
            SUM=SUM+RHO(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      SUM=SUM/(NMESHX*NMESHY*NMESHZ)
C..     
      WRITE(11,*) 1,0
      WRITE(11,*) 2000,0,0,0
      WRITE(11,*) 1
      WRITE(11,*) 3
      WRITE(11,*)  2*SUM ,0
      WRITE(11,*)  SUM ,0
      WRITE(11,*)  SUM/2 ,0
      WRITE(11,*)  1,NMESHX,NMESHY,NMESHZ
      WRITE(11,*) 0
      WRITE(11,*) '0.   0.  0.'
      WRITE(11,*) ALAT(1),0.0_dp,ALAT(2),0.0000_dp,0.0000_dp,ALAT(3)
      WRITE(11,*) '8    12'
      WRITE(11,*) 0,   0,   0
      WRITE(11,*) ALAT(1),   0,   0
      WRITE(11,*) 0,   ALAT(2),   0
      WRITE(11,*) 0,   0,   ALAT(3)
      WRITE(11,*) ALAT(1),   ALAT(2),    0
      WRITE(11,*) ALAT(1),   0,   ALAT(3)
      WRITE(11,*) 0,   ALAT(2),ALAT(3)
      WRITE(11,*) ALAT(1),ALAT(2),ALAT(3)
      WRITE(11,*) '1   2'
      WRITE(11,*) '1   3'
      WRITE(11,*) '1   4 '
      WRITE(11,*) '2   5 '
      WRITE(11,*) '2   6'
      WRITE(11,*) '3   5'
      WRITE(11,*) '3   7'
      WRITE(11,*) '5   8 '
      WRITE(11,*) '4   6'
      WRITE(11,*) '4   7'
      WRITE(11,*) '7   8 '
      WRITE(11,*) '6   8 '
      WRITE(11,*) '0'
      ENDIF
      IF(TGOPENMOL) THEN
        FILENM3=FILENM(1:NC)//'.txt'
c        FILENM3=FILENM(1:NC)//'.plt'
        OPEN(12,FILE=FILENM3,STATUS='UNKNOWN')
c        OPEN(12,FILE=FILENM3,STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(12,'(I1,1X,I3)') 3,2
c        WRITE(12,*) 3,2
        WRITE(12,'(3I3)') NMESHX,NMESHY,NMESHZ
c        WRITE(12,*) NMESHX,NMESHY,NMESHZ
        WRITE(12,'(6F13.6)') -ALAT(3)/2.0_dp,ALAT(3)/2.0_dp,
     &    -alat(2)/2.0_dp,ALAT(2)/2.0_dp,-alat(1)/2.0_dp,ALAT(1)/2.0_dp
c        WRITE(12,*)-ALAT(1)/2.,ALAT(1)/2.0_dp,-alat(2)/2.0_dp,
c     &   ALAT(2)/2.0_dp,-alat(3)/2.0_dp,ALAT(3)/2.0_dp
      ENDIF
      DO K=1,NMESHZ
        DO J=1,NMESHY
        DO I=1,NMESHX
          WRITE(10,*) I,J,K,RHO(I,J,K)
        ENDDO
      ENDDO
      ENDDO
      CLOSE(10)
      WRITE(6,*) ' WRITTEN FILE : ' , FILENM
      IF(TMOVIE) THEN
        DO I=1,NMESHX
          DO J=1,NMESHY
            DO K=1,NMESHZ
            WRITE(11,*) RHO(I,J,K)
          ENDDO
        ENDDO
        ENDDO
        CLOSE(11)
        WRITE(6,*) ' WRITTEN FILE : ' , FILENM2
      ENDIF
      IF(TGOPENMOL) THEN
        RS3=RS*RS*RS
          DO K=1,NMESHZ
            DO J=1,NMESHY
              WRITE(12,'(6E13.6)') (RHO(I,J,K)*RS3,I=1,NMESHX)
          ENDDO
        ENDDO
c        WRITE(12,'(6E13.6)') RHO
c        WRITE(12) RHO
        CLOSE(12)
        WRITE(6,*) ' WRITTEN FILE : ' , FILENM3
      ENDIF
C ==--------------------------------------------------------------==
      RETURN
      END
