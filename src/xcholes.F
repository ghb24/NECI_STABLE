      SUBROUTINE XCHOLES(EVEC,PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     &    NHG,RHO,XCHOLE,SPAC,RS,ALAT,OMEGA,NMRKS,NDET,NEL,NEVAL)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      INTEGER :: NEL,NMSH,NMAX,NEVAL,NDET,NHG
      real(dp) :: EVEC(NDET,NEVAL),PSIR(-NMSH:NMSH)
      real(dp) :: RHO(NMSH,NMSH,NMSH),XCHOLE(NMSH,NMSH,NMSH)
      real(dp) :: SITAB(NMSH,NMAX),ALAT(3),RS,OMEGA,SPAC
      integer :: NMRKS(NEL,*)
      INTEGER G1(5,NHG),LEN
      CHARACTER(50) :: FILENM
      CHARACTER(50) :: CHART
      integer :: NC,II,IY,IZ,IX
C..
      WRITE(6,*) 
     & '       ==--------------------------------------------------== '
      WRITE(6,*) '          COMPUTING EXCHANGE-CORRELATION HOLES '
      WRITE(6,*) 
     & '       ==--------------------------------------------------== '
      CHART='XCHD'
      NC=4
C..
      II=1
C..100
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,1,0,0,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,1,0,0,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..010
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,0,1,0,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,0,1,0,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..001
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,0,0,1,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,0,0,1,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..011
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,0,1,1,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,0,1,1,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..101
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,1,0,1,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,1,0,1,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..110
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,1,1,0,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,1,1,0,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..111             
      CALL GEN_XCHOLE(EVEC(1,II),PSIR,IX,IY,IZ,G1,SITAB,NMAX,NMSH,
     & NHG,1,1,1,RHO,.FALSE.,XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,
     & NEL)
      CALL GET_FILENM(NC,CHART,II,1,1,1,FILENM,LEN)
      CALL WRITE_LINE(LEN,FILENM,PSIR,-NMSH,NMSH,IX,IY,IZ,SPAC,RS)
C..
      RETURN
      END
C ==---------------------------------------------------------------------==
      SUBROUTINE GET_FILENM(L1,FILENM1,II,I,J,K,FILENM,LEN)
      IMPLICIT NONE
      CHARACTER FILENM1*50,FILEAD1*1,FILEAD2*2,FILEAD3*3,
     &       FILE3*3,FILENM*50,FILENM2*8
      integer :: LEN,I,J,K,II,L1
      WRITE(FILE3,'(3I1)') I,J,K
      FILENM2=FILENM1(1:L1)//FILE3//'.'
      IF(II.LT.9) THEN
        WRITE(FILEAD1,'(I1)') II
        FILENM=FILENM2//FILEAD1
        LEN=L1+4+1
      ELSEIF(II.LT.99) THEN
        WRITE(FILEAD2,'(I2)') II
        FILENM=FILENM2//FILEAD2
        LEN=L1+4+2
      ELSEIF(II.LT.999) THEN
        WRITE(FILEAD3,'(I3)') II
        FILENM=FILENM2//FILEAD3
      LEN=L1+4+3
      ENDIF
      RETURN
      END
C ==---------------------------------------------------------------------==
      SUBROUTINE WRITE_LINE(NC,FILENM,FUNC,IBEG,IEND,IX,IY,IZ,
     &       SPAC,RS)
!      IMPLICIT real(dp)(A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: IBEG,IEND,IX,IY,IZ,NC,I
      real(dp) :: FUNC(IBEG:IEND),SPAC,RS,RS3
      CHARACTER FILENM*50
C..
      RS3=RS*RS*RS
      OPEN(10,FILE=FILENM(1:NC),STATUS='UNKNOWN')
      IF(IX.LT.0) GOTO 100
      WRITE(10,'(A,3I3)') '#Observation point:',IX,IY,IZ
 100  CONTINUE
      DO I=IBEG,IEND
        IF(FUNC(I).NE.0.0_dp.OR.IBEG.EQ.1.OR.IEND.EQ.1) THEN
          WRITE(10,*) I*SPAC,FUNC(I)*RS3
        ENDIF
      ENDDO
      CLOSE(10)
      WRITE(6,*)  
     &    ' WROTE FILE : ' , FILENM(1:NC)
      RETURN
      END
C ==---------------------------------------------------------------------==
      SUBROUTINE PLANARAV(RHO,NMESHX,DLINE,I,J,K,SPAC,ALAT)
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: NMESHX,I,J,K,IX,IY,IZ
      real(dp) :: RHO(NMESHX,NMESHX,NMESHX)
      real(dp) :: DLINE(NMESHX),ALAT(3)
      real(dp) :: spac,SUM
C..
      DLINE=0.0_dp
      IF(I.EQ.0.AND.J.EQ.0.AND.K.EQ.1) THEN
         DO IZ=1,NMESHX
            SUM=0.0_dp
            SPAC=ALAT(3)/NMESHX
            DO IY=1,NMESHX
               DO IX=1,NMESHX
                  SUM=SUM+RHO(IX,IY,IZ)
               ENDDO
            ENDDO
            DLINE(IZ)=SUM/real(NMESHX*NMESHX,dp)
         ENDDO
      ELSEIF(I.EQ.1.AND.J.EQ.0.AND.K.EQ.0) THEN
         SPAC=ALAT(1)/NMESHX
         DO IX=1,NMESHX
            SUM=0.0_dp
            DO IY=1,NMESHX
               DO IZ=1,NMESHX
                  SUM=SUM+RHO(IX,IY,IZ)
               ENDDO
            ENDDO
            DLINE(IX)=SUM/real(NMESHX*NMESHX,dp)
         ENDDO
      ELSEIF(I.EQ.0.AND.J.EQ.1.AND.K.EQ.0) THEN
         SPAC=ALAT(1)/NMESHX
         DO IY=1,NMESHX
            SUM=0.0_dp
            DO IZ=1,NMESHX
               DO IX=1,NMESHX
                  SUM=SUM+RHO(IX,IY,IZ)
               ENDDO
            ENDDO
            DLINE(IY)=SUM/real(NMESHX*NMESHX,dp)
         ENDDO
      ELSE
         WRITE(6,*) ' AVERAGING DIRECTION NOT PROGRAMMED!!'
      ENDIF
      RETURN
      END
