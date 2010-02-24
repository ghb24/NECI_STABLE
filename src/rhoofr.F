      SUBROUTINE NECI_RHOOFR(NHG,CK,G1,RHO,NMSH,SITAB,NMAX,NMRKS,
     &                       NEL,NDET,NEVAL,RS,ALAT,OMEGA)
      use global_utilities
      IMPLICIT REAL*8 (A-H,O-Z)
C..
      real*8, allocatable, save :: CFS(:) ! NDET
      integer, save :: tagCFS=0
C..
      INTEGER G1(*)
      DIMENSION CK(NDET,NEVAL),NMRKS(NEL,*)
      DIMENSION RHO(NMSH,NMSH,NMSH),AUX(NMSH,NMSH,NMSH)
      DIMENSION SITAB(NMSH,NMAX),ALAT(3)
      DATA EPS/1.D-8/
      DATA IFIRST/0/
      SAVE IFIRST
C..
      IF(IFIRST.EQ.0) THEN
        IFIRST=1
        SITAB=0.d0
        CALL GEN_SITAB(NMAX,NMSH,SITAB)
      ENDIF
C..
      allocate(CFS(NDET))
      call LogMemAlloc('CFS',NDET,8,'NECI_RHOOFR',tagCFS)
      CFS=0.d0
      CALL DCOPY(NDET,CK(1,1),1,CFS,1)
C..
      AUX=0.d0
      RHO=0.d0
C..
      DO KI=1,NDET
        IF(ABS(CFS(KI)).LT.EPS) GOTO 300
        DO KJ=KI,NDET
          IF(ABS(CFS(KJ)).LT.EPS) GOTO 500
C..IC is the number of pairs which are the *SAME* in both determinants
          ICC=2
          IF(KI.EQ.KJ) ICC=1
C..
          IC=0
          DO I=1,NEL
            DO J=1,NEL
              IF(NMRKS(I,KJ).EQ.NMRKS(J,KI)) IC=IC+1
            ENDDO
          ENDDO
C..If determinants differ by more than 1 spin orbital IC < (NEL-1) 
          IF(IC.LT.(NEL-1)) GOTO 500
          CALL RLS(NEL,NMSH,NHG,NDET,NMRKS(1:NEL,KI),
     &             NMRKS(1:NEL,KJ),G1,SITAB,NMAX,IC,AUX)
          CF=CFS(KI)*CFS(KJ)*DFLOAT(ICC)
          DO K=1,NMSH-1
            DO J=1,NMSH-1
              DO I=1,NMSH-1
                RHO(I,J,K)=RHO(I,J,K)+CF*AUX(I,J,K)
              ENDDO
            ENDDO
          ENDDO
500     CONTINUE
        ENDDO
300   CONTINUE
      ENDDO
      CST=1.D0/OMEGA
      CALL DSCAL(NMSH*NMSH*NMSH,CST,RHO,1)
C..Calculate integral of RHO 
      SUM=0.D0
      DO K=1,NMSH-1
        DO J=1,NMSH-1
          DO I=1,NMSH-1
            SUM=SUM+RHO(I,J,K) 
          ENDDO
        ENDDO
      ENDDO
      SUM=SUM*OMEGA/DFLOAT(NMSH*NMSH*NMSH)
      WRITE(6,'(A,F12.6)') ' INTEGRATED DENSITY : ' , SUM
      CALL WRITE_RHO(9,'DENSITY.1',RHO,NMSH,NMSH,NMSH,ALAT,
     &    .FALSE.,.TRUE.,RS)
      RETURN
      END
C ==---------------------------------------------------------==
      SUBROUTINE RLS(NEL,NMSH,NHG,NDET,NDET1,NDET2,G1,
     &     SITAB,NMAX,IC,AUX)
      use SystemData, only: BasisFN
      IMPLICIT REAL*8 (A-H,O-Z)
C..
      DIMENSION NDET1(NEL),NDET2(NEL),AUX(NMSH,NMSH,NMSH)
      DIMENSION SITAB(NMSH,NMAX),AUXX(NEL),AUXY(NEL),AUXZ(NEL)
      TYPE(BasisFN) G1(NHG)
      INTEGER N1E(NHG),N2E(NHG)
C..
      N1E(1:NHG)=0
      N2E(1:NHG)=0
      CALL LINEUP(NEL,NDET1,NDET2,NHG,N1E,N2E,IFLAG,IPTOT)
      ND1=0
      ND2=0
C..
      AUXX=0.d0
      AUXY=0.d0
      AUXZ=0.d0
C..
      IF(IC.EQ.(NEL-1)) THEN
        AUX=0.d0
C..The determinants differ by ONE spin orbital
        ND1=0
        ND2=0
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
100     CONTINUE
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
200     CONTINUE
C..Check spin indices
        CALL CKRN(NHG,G1,ND1,ND2,ID)
        IF(ID.EQ.1) THEN
          DO K=1,NMSH-1
            AUXZ1=SITAB(K,G1(ND1)%k(3))
            AUXZ2=SITAB(K,G1(ND2)%k(3))
            DO J=1,NMSH-1
              AUXY1=SITAB(J,G1(ND1)%k(2))
              AUXY2=SITAB(J,G1(ND2)%k(2))
              DO I=1,NMSH-1
                AUXX1=SITAB(I,G1(ND1)%k(1))
                AUXX2=SITAB(I,G1(ND2)%k(1))
C..
                AUX(I,J,K)=AUXZ1*AUXZ2*AUXY1*AUXY2*
     &                        AUXX1*AUXX2
C..
              ENDDO
            ENDDO
          ENDDO
        ELSE
          AUX=0.d0
          RETURN
        ENDIF
      ELSEIF(IC.EQ.NEL) THEN
        AUX=0.d0
C..Differ by none
        DO K=1,NMSH-1
          DO IK=1,NEL
            AUXZ(IK)=SITAB(K,G1(NDET1(IK))%k(3))
          ENDDO
          DO J=1,NMSH-1
            DO IJ=1,NEL
              AUXY(IJ)=SITAB(J,G1(NDET1(IJ))%k(2)) 
            ENDDO
            DO I=1,NMSH-1
              DO II=1,NEL
                AUXX(II)=SITAB(I,G1(NDET1(II))%k(1))
              ENDDO
C..
              DO IJK=1,NEL
                AUX(I,J,K)=AUX(I,J,K)+
     &            AUXZ(IJK)*AUXZ(IJK)*
     &            AUXY(IJK)*AUXY(IJK)*
     &            AUXX(IJK)*AUXX(IJK)
              ENDDO
C..
            ENDDO
          ENDDO
        ENDDO
      ELSE
        AUX=0.d0
        RETURN
      ENDIF
C..
      CALL DSCAL(NMSH*NMSH*NMSH,DFLOAT(IPTOT),AUX,1)
      RETURN
      END
C ==---------------------------------------------------------==
      SUBROUTINE GEN_SITAB(NMAX,NMSH,SITAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SITAB(NMSH,NMAX)
C..
      PI=ACOS(-1.D0)
      DO I=1,NMAX
        DO J=1,NMSH
          X=DFLOAT(J)/DFLOAT(NMSH)
          SITAB(J,I)=SIN(PI*I*X)*SQRT(2.D0)
        ENDDO
      ENDDO
      RETURN
      END
C ==---------------------------------------------------------==
      SUBROUTINE CKRN(NHG,G1,ND1,ND2,ID)
      use SystemData, only: BasisFN
      IMPLICIT REAL*8 (A-H,O-Z)
      TYPE(BasisFN) G1(NHG)
C..
      ID=0
      IF(G1(ND1)%Ms.EQ.G1(ND2)%Ms) ID=1
      RETURN
      END
