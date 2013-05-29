! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
      SUBROUTINE RLSXC1_SC(NEL,NMSH,NMAX,SITAB,IOBS,JOBS,KOBS,
     &  NDET1,NDET2,G1,NHG,AUXC)
      use SystemData, only: BasisFN
!      IMPLICIT real(dp) (A-H,O-Z)
      use constants, only: dp
      IMPLICIT NONE
      integer :: NEL,NHG,NMSH,NMAX
      real(dp) :: AUXC(NMSH,NMSH,NMSH)
      integer :: NDET1(NEL),NDET2(NEL)
C..Cube array
      real(dp) :: SITAB(NMSH,NMAX) 
      TYPE(BasisFN) G1(NHG)
      INTEGER TMP1(NEL),TMP2(NEL)
      INTEGER IDS1(NEL-1),IDS2(NEL-1)
      EXTERNAL CHK,U1
      real(dp) :: U2ND2R2,U1ND2R2,U1ND1R2,U2ND1R2
      real(dp) :: U2ND2R1,U1ND2R1,U1ND1R1,U2ND1R1,U1KIR1,U2KIR2,U1
      integer :: IQ,I,J,K,ND1,II,KI,KOBS,JOBS,IOBS,IDS,IX,IY,IZ,ND2
C..
      TMP1(1:NEL)=0
      TMP2(1:NEL)=0
      AUXC=0.0_dp
C..Compare NDET1 to NDET2
      IQ=0
      DO I=1,NEL
        K=0
        DO J=1,NEL
          IF(NDET1(I).EQ.NDET2(J)) CONTINUE
          IF(NDET1(I).NE.NDET2(J)) K=K+1
        ENDDO
        IF(K.EQ.NEL) THEN
          IQ=IQ+1
          ND1=NDET1(I)
          IF(IQ.GE.1) GOTO 100
        ENDIF
      ENDDO
100   CONTINUE
C..Now we compare NDET2 to NDET1
      IQ=0
      DO I=1,NEL
        K=0
        DO J=1,NEL
          IF(NDET2(I).EQ.NDET1(J)) CONTINUE
          IF(NDET2(I).NE.NDET1(J)) K=K+1
        ENDDO
        IF(K.EQ.NEL) THEN
          IQ=IQ+1
          ND2=NDET2(I)
          IF(IQ.GE.1) GOTO 200
        ENDIF
      ENDDO
200   CONTINUE
C..Here we make sure that the last elements of TMP12 are ND1 and ND2
      TMP1(NEL)=ND1
      TMP2(NEL)=ND2
C..Add in rest of NDET1 and NDET2 to TMP12
      II=0
      DO I=1,NEL
        IF(NDET1(I).NE.ND1) THEN
          II=II+1
          TMP1(II)=NDET1(I)
        ENDIF
      ENDDO
C..
      II=0
      DO I=1,NEL
        IF(NDET2(I).NE.ND2) THEN
          II=II+1
          TMP2(II)=NDET2(I)
        ENDIF
      ENDDO
C..Kronecker deltas
      IQ=0
      DO I=1,NEL
        IF(NDET1(I).NE.ND1) THEN
          IQ=IQ+1
          CALL CHK(G1(NDET1(I))%Ms,G1(ND1)%Ms,IDS1(IQ))
          CALL CHK(G1(NDET2(I))%Ms,G1(ND2)%Ms,IDS2(IQ))
          IF(IQ.GE.(NEL-1)) GOTO 500
        ENDIF
      ENDDO
500   CONTINUE
C..
      CALL CHK(G1(ND1)%Ms,G1(ND2)%Ms,IDS)
      DO KI=1,NEL-1
       
        DO IZ=1,NMSH
          DO IY=1,NMSH 
            DO IX=1,NMSH
C..
              U1ND1R1=U1(ND1,NMAX,NMSH,NHG,G1,SITAB,IOBS,JOBS,KOBS)
              U1ND2R1=U1(ND2,NMAX,NMSH,NHG,G1,SITAB,IOBS,JOBS,KOBS)
              U2KIR2=U1(TMP1(KI),NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
              U1KIR1=U1(TMP1(KI),NMAX,NMSH,NHG,G1,SITAB,
     &                 IOBS,JOBS,KOBS)    
              U2ND1R2=U1(ND1,NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)
              U2ND2R2=U1(ND2,NMAX,NMSH,NHG,G1,SITAB,IX,IY,IZ)

              AUXC(IX,IY,IZ)=AUXC(IX,IY,IZ)+
     &          U1ND1R1*U1ND2R1*U2KIR2*U2KIR2*real(IDS,dp)  
     &         +U1KIR1*U1KIR1*U2ND1R2*U2ND2R2*real(IDS,dp)   
     &         -IDS1(KI)*IDS2(KI)*U1ND1R1*U1KIR1*U2KIR2*U2ND2R2   
     &         -IDS1(KI)*IDS2(KI)*U1KIR1*U1ND2R1*U2ND1R2*U2KIR2    
            ENDDO
          ENDDO
        ENDDO
C..
      ENDDO  
      RETURN
      END 
