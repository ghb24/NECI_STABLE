C.. This is based upon the GEN_COUL for electrons in a cubic box
C.. but with a simpler coulomb integral
      SUBROUTINE GEN_COUL_HUBNPBC(NEL,NBASISMAX,
     &   NHG,G1,NMSH,NMAX,FCK,ZIA,UMAT,ISS,THUB,UHUB,OMEGA)
      USE UMatCache , only : UMatInd
      IMPLICIT NONE
      INCLUDE 'basis.inc'
      TYPE(BASISFN) G1(NHG)
      INTEGER NBASISMAX(5,3),ISS,NEL
      REAL*8 UMAT(NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
      REAL*8 UHUB,OMEGA
C..Cube arrays
      COMPLEX*16 FCK(-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1)
      COMPLEX*16 ZIA(-NMSH/2:NMSH/2,NMAX,NMAX)
      LOGICAL THUB,CHKMOMEQ,LCHKSYM
      INTEGER NHG,NMSH,NMAX
      INTEGER KA(4),KB(4),M
      INTEGER A,B,C,LL(3),KK(3)
      INTEGER ISUB,II,ISPINSKIP
      REAL*8 SUMS(3),SUM
      INTEGER I,J,K,L,ID1,ID2,ID3,ID4
      INCLUDE 'cons.inc'
C..This routine generates *ALL* possible combinations of Coulomb integrals
C..stored in the form (u1 u2 | U | u1' u2') = UMAT(n1 n2 n3 n4)
C..
C..This first call calculates the inner integral 
C..The call to SCOUL calculates the outer integral 
C ==--------------------------------------------------------------------==
      CALL TISET('GENCOULUEG',ISUB)
C ==--------------------------------------------------------------------==
      OPEN(10,FILE='UMAT',STATUS='UNKNOWN')
      II=0
      ISPINSKIP=NBASISMAX(2,3)
      LL(1)=NBASISMAX(1,2)-NBASISMAX(1,1)+1
      LL(2)=NBASISMAX(2,2)-NBASISMAX(2,1)+1
      LL(3)=NBASISMAX(3,2)-NBASISMAX(3,1)+1
      DO I=1,NHG,ISPINSKIP
        DO J=I,NHG,ISPINSKIP
          DO K=1,I,ISPINSKIP
            DO L=1,J,ISPINSKIP
              SUM=1.D0
C..Original call
C              CALL SLATCOULFOU(G1(1,I),G1(1,J),
C     &		G1(1,K),G1(1,L),NMSH,FCK,ZIA,NMAX,SUM)
              CALL GTID(NBASISMAX,I,ID1)
              CALL GTID(NBASISMAX,J,ID2)
              CALL GTID(NBASISMAX,K,ID3)
              CALL GTID(NBASISMAX,L,ID4)
              SUM=UHUB
              DO A=1,3
                SUMS(A)=0.D0
                DO B=1,LL(A)
                 SUMS(A)=SUMS(A)
     &            +SIN(B*PI*(G1(I)%K(A)-NBASISMAX(A,1)+1)/(LL(A)+1))
     &            *SIN(B*PI*(G1(J)%K(A)-NBASISMAX(A,1)+1)/(LL(A)+1))
     &            *SIN(B*PI*(G1(K)%K(A)-NBASISMAX(A,1)+1)/(LL(A)+1))
     &            *SIN(B*PI*(G1(L)%K(A)-NBASISMAX(A,1)+1)/(LL(A)+1))
C     &                          *SIN(B*PI*G1(A,J)/(LL(A)+1))
C     &                          *SIN(B*PI*G1(A,K)/(LL(A)+1))
C     &                          *SIN(B*PI*G1(A,L)/(LL(A)+1))
                ENDDO
                SUM=SUM*SUMS(A)*(2.D0/(LL(A)+1))**2
              ENDDO

C.. UMAT is stored as just spatial orbitals (not spinorbitals)
C..   we store <IJ|KL>  (
! Get the index of physical order UMAT element <IJ|KL>.  Indices are internally reordered such that I>=K, J>=L,(I,K)>=(J,L)
              UMAT(UMatInd(ID1,ID2,ID3,ID4))=SUM
              IF(ABS(SUM).GT.1.D-10)
     &            WRITE(10,'(4I7,F19.9)') ID1,ID2,ID3,ID4,SUM
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(10)
      IF(THUB) THEN
C.. V0 is subtracted from the diagonal elements of the hamiltonian for
C.. the UEG, but we don't want this for the hubbard model
         FCK(0,0,0)=(0.D0,0.D0)
      ENDIF
      WRITE(6,*) ' !!! FINISHED CALCULATING ALL 2E INTEGRALS !!! '
C ==--------------------------------------------------------------------==
      CALL TIHALT('GENCOULUEG',ISUB)
C ==--------------------------------------------------------------------==
      RETURN
      END
C ==------------------------------------------------------------------==
C.. This is based upon the GEN_COUL for electrons in a cubic box
C.. but with a simpler coulomb integral
      SUBROUTINE GEN_COUL_UEG(NEL,NBASISMAX,
     &   NHG,G1,NMSH,NMAX,FCK,ZIA,UMAT,ISS,THUB,UHUB,OMEGA,ALAT)
      USE UMatCache , only : UMatInd
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'basis.inc'
      TYPE(BASISFN) G1(NHG)
      INTEGER NBASISMAX(5,5),ISS 
      DIMENSION UMAT(NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
C..Cube arrays
      COMPLEX*16 FCK(-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1)
      COMPLEX*16 ZIA(-NMSH/2:NMSH/2,NMAX,NMAX)
      LOGICAL THUB,CHKMOMEQ,LCHKSYM
      INTEGER KA(4),KB(4),M
      INTEGER A,B,C,D,TX,TY
      REAL*8 LX,LY,P,Q,R,T2
      REAL*8 ALAT(3)
      COMPLEX*16 S,CI
      INCLUDE 'cons.inc'
C..This routine generates *ALL* possible combinations of Coulomb integrals
C..stored in the form (u1 u2 | U | u1' u2') = UMAT(n1 n2 n3 n4)
C..
C..This first call calculates the inner integral 
C..The call to SCOUL calculates the outer integral 
C ==--------------------------------------------------------------------==
      CALL TISET('GENCOULUEG',ISUB)
C ==--------------------------------------------------------------------==
      OPEN(10,FILE='UMAT',STATUS='UNKNOWN')
      II=0
      ISPINSKIP=NBASISMAX(2,3)
      DO I=1,NHG,ISPINSKIP
        DO J=I,NHG,ISPINSKIP
          DO K=1,I,ISPINSKIP
            DO L=1,J,ISPINSKIP
              SUM=0.D0
C..Original call
C              CALL SLATCOULFOU(G1(1,I),G1(1,J),
C     &		G1(1,K),G1(1,L),NMSH,FCK,ZIA,NMAX,SUM)
              CALL GTID(NBASISMAX,I,ID1)
              CALL GTID(NBASISMAX,J,ID2)
              CALL GTID(NBASISMAX,K,ID3)
              CALL GTID(NBASISMAX,L,ID4)
              IF(THUB) THEN
                 CALL SETUPSYM(NBASISMAX,KA)
                 CALL SETUPSYM(NBASISMAX,KB)
                 CALL ADDELECSYM(K,G1,NBASISMAX,KA)
                 CALL ADDELECSYM(L,G1,NBASISMAX,KA)
                 CALL ADDELECSYM(I,G1,NBASISMAX,KB)
                 CALL ADDELECSYM(J,G1,NBASISMAX,KB)
                 CALL ROUNDSYM(KA,NBASISMAX)
                 CALL ROUNDSYM(KB,NBASISMAX)
C                 IF(I.EQ.1.AND.J.EQ.1.AND.K.LE.L)
C     &          WRITE(40,"(6I3,$)") (G1(II,K),II=1,2),(G1(II,L),II=1,2),
C     &                  (KA(II),II=1,2)
                 IF(LCHKSYM(KA,KB)) THEN
C                 IF(CHKMOMEQ(KA,KB,NBASISMAX,3)) THEN
C..Omega for the Hubbard model is the number of sites in the lattice
C,(G1(II,K),II=1,2),
C     &                  (G1(II,L),II=1,2),(KB(II),II=1,2)
C                    WRITE(40,"(A10,$)") "Y"
                    SUM=UHUB/OMEGA
                 ELSE
                    SUM=0.D0
                 ENDIF
C                 WRITE(40,*)
             IF(.FALSE.) THEN
                 S=0.D0
                 CI=(0.D0,1.D0)         
C.. Check the Coulomb matrix elements by summing manually
                 C=G1(K)%K(1)+G1(L)%K(1)-G1(I)%K(1)-G1(J)%K(1)
                 D=G1(K)%K(2)+G1(L)%K(2)-G1(I)%K(2)-G1(J)%K(2)
                 TX=NBASISMAX(1,4)
                 TY=NBASISMAX(2,4)
                 T2=TX*TX+TY*TY
                 LX=NBASISMAX(1,5)
                 LY=NBASISMAX(2,5)
                 DO A=1,NHG,2
                  Q=(G1(A)%K(1)*TX+G1(A)%K(2)*TY)/T2
                  P=(G1(A)%K(1)*-TY+G1(A)%K(2)*TX)/T2
                  S=S+EXP(CI*(P*C/LX+Q*D/LY)
     &            *2*PI)*UHUB/(OMEGA**2)
                 ENDDO
                 WRITE(41,"(2I3,$)") C,D
                 IF(ABS(REAL(S)-SUM).GT.1.D-7) THEN
                    WRITE(41,*) SUM,S
                 ELSE
                    WRITE(41,*) SUM,SUM
                 ENDIF
C                 SUM=REAL(S)
         ENDIF
              ELSE
C.. we splice in the simple version of the fourier transform

C.. New, but not needed
C                 CALL SETUPSYM(NBASISMAX,KA)
C                 CALL SETUPSYM(NBASISMAX,KB)
C                 CALL ADDELECSYM(K,G1,NBASISMAX,KA)
C                 CALL ADDELECSYM(L,G1,NBASISMAX,KA)
C                 CALL ADDELECSYM(I,G1,NBASISMAX,KB)
C                 CALL ADDELECSYM(J,G1,NBASISMAX,KB)
C                 CALL ROUNDSYM(KA,NBASISMAX)
C                 CALL ROUNDSYM(KB,NBASISMAX)
C                 IF(LCHKSYM(KA,KB)) THEN
                   A=G1(I)%K(1)-G1(K)%K(1)
                   B=G1(I)%K(2)-G1(K)%K(2)
                   C=G1(I)%K(3)-G1(K)%K(3)
                   IF(    (G1(L)%K(1)-G1(J)%K(1)).EQ.A
     &               .AND.(G1(L)%K(2)-G1(J)%K(2)).EQ.B
     &               .AND.(G1(L)%K(3)-G1(J)%K(3)).EQ.C
     &               .AND.((A.NE.0).OR.(B.NE.0).OR.(C.NE.0))) THEN
C                     WRITE(6,*) "(",I,J,"|",K,L,")",A,B,C
                     
CC  AJWT  <IJ|r_12^-1|KL> = v_(G_I-G_K) delta_((G_I-G_K)-(G_L-G_J)
CC  v_G = 4 Pi/ G**2.  G=2 Pi/L(nx,ny,nx) etc.
                     SUM=((A/ALAT(1))**2+(B/ALAT(2))**2)
                     IF(ALAT(3).NE.0.D0) SUM=SUM+(C/ALAT(3))**2
                     SUM=1/(PI*SUM*OMEGA)
                   ELSE
                     SUM=0.D0
                   ENDIF

CC Old UEG - doesn't seem to work
C                 CALL SLATCOULFOUCOU(G1(1,I),G1(1,J),
C     &            G1(1,K),G1(1,L),NMSH,FCK,ZIA,NMAX,SUM)
C                 ENDIF
               ENDIF
C.. UMAT is stored as just spatial orbitals (not spinorbitals)
C..   we store <IJ|KL>  (
! Get the index of physical order UMAT element <IJ|KL>.  Indices are internally reordered such that I>=K, J>=L,(I,K)>=(J,L)
              UMAT(UMatInd(ID1,ID2,ID3,ID4))=SUM
              IF(ABS(SUM).GT.1.D-10)
     &            WRITE(10,'(4I7,F19.9)') ID1,ID2,ID3,ID4,SUM
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(10)
      IF(THUB) THEN
C.. V0 is subtracted from the diagonal elements of the hamiltonian for
C.. the UEG, but we don't want this for the hubbard model
         FCK(0,0,0)=(0.D0,0.D0)
      ENDIF
      WRITE(6,*) ' !!! FINISHED CALCULATING ALL 2E INTEGRALS !!! '
C ==--------------------------------------------------------------------==
      CALL TIHALT('GENCOULUEG',ISUB)
C ==--------------------------------------------------------------------==
      RETURN
      END
C ==------------------------------------------------------------------==
C ========================================================
      SUBROUTINE SLATCOULFOUCOU(G1,G2,G1P,G2P,N,
     &     CK,ZIA,NMAX,OUT)
C=========================================================
C.. Returns the Couloumb integral between the Slater
C.. determinants of plane wave basis, using the Fourier
C.. method.
C=========================================================
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER P,Q
      INTEGER G1(4),G2(4),G1P(4),G2P(4)
      COMPLEX*16 ZIA(-N/2:N/2,NMAX,NMAX)
      COMPLEX*16 CK(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1)
      INTEGER GD1,GD2,I,GD(3)
      LOGICAL T,T2
C.. check if the spins are correct.
      IF(G1(4).NE.G1P(4).OR.G2(4).NE.G2P(4)) THEN
         OUT=0.D0
         RETURN
      ENDIF
      T=.TRUE.
      T2=.TRUE.
      TOT=0
      DO I=1,3
         GD1=G1(I)-G1P(I)
         GD2=G2(I)-G2P(I)
C.. The 2* is to account for the fact that the UEG has basis fns which have
C.. k= 2 pi n/L rather than pi n/L
         GD(I)=2*GD1
         IF(GD1.NE.GD2) T=.FALSE.
         IF(GD1.NE.0) T2=.FALSE.
      ENDDO
      IF(T2) THEN
C.. GD1=0,0,0, which has to be removed as this cancels the background
C..   term
         T=.FALSE.
C         WRITE(6,*) CK(GD(1),GD(2),GD(3)),GD(1),GD(2),GD(3)
      ENDIF
C      IF(T) THEN
         OUT=REAL(CK(2*GD(1),2*GD(2),2*GD(3)))
c         OUT=4*3.1415926535D0/(GD(1)*GD(1)+GD(2)*GD(2)+GD(3)*GD(3))
C      ELSE
C         OUT=0.D0
C      ENDIF
      RETURN
      END
C
      LOGICAL FUNCTION CHKMOMEQ(K1,K2,NBASISMAX,KDIM)
C.. NB the third column of NBASISMAX tells us whether it is tilted
         INTEGER KDIM,K1(KDIM),K2(KDIM),NBASISMAX(5,3)
         INTEGER J,LDIM,KK1,KK2
         LOGICAL LCMP
         LCMP=.TRUE.
         DO J=1,3
            IF(NBASISMAX(1,3).EQ.1) THEN
C.. we have a tilted lattice
               IF(J.EQ.1) THEN
                  KK1=K1(1)+K1(2)
                  KK2=K2(1)+K2(2)
                  LDIM=NBASISMAX(1,2)*2
               ELSEIF(J.EQ.2) THEN
                  KK1=K1(1)-K1(2)
                  KK2=K2(1)-K2(2)
                  LDIM=NBASISMAX(1,2)*2
               ELSE
                  KK1=K1(J)
                  KK2=K2(J)
                  LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
               ENDIF

            ELSE
               KK1=K1(J)
               KK2=K2(J)
               LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
            ENDIF
            KK1=MOD(KK1,LDIM)
            IF(KK1.LT.NBASISMAX(J,1)) KK1=KK1+LDIM
            IF(KK1.GT.NBASISMAX(J,2)) KK1=KK1-LDIM
            KK2=MOD(KK2,LDIM)
            IF(KK2.LT.NBASISMAX(J,1)) KK2=KK2+LDIM
            IF(KK2.GT.NBASISMAX(J,2)) KK2=KK2-LDIM
            IF(KK1.NE.KK2) LCMP=.FALSE.
         ENDDO
         IF(KDIM.EQ.4.AND.K1(4).NE.K2(4)) LCMP=.FALSE.
         CHKMOMEQ=LCMP
         RETURN
      END
 
      SUBROUTINE GetHubUMatEl(I,J,K,L,UMAT,NBASISMAX,G1,SUM)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BASISFN) G1(*),KA,KB
         INTEGER NBASISMAX(5,5)
C..Cube arrays
         LOGICAL CHKMOMEQ,LCHKSYM
         INTEGER I,J,K,L
         REAL*8 SUM,UMAT
                 CALL SETUPSYM(NBASISMAX,KA)
                 CALL SETUPSYM(NBASISMAX,KB)
                 CALL ADDELECSYM(K*2,G1,NBASISMAX,KA)
                 CALL ADDELECSYM(L*2,G1,NBASISMAX,KA)
                 CALL ADDELECSYM(I*2,G1,NBASISMAX,KB)
                 CALL ADDELECSYM(J*2,G1,NBASISMAX,KB)
                 CALL ROUNDSYM(KA,NBASISMAX)
                 CALL ROUNDSYM(KB,NBASISMAX)
                 IF(LCHKSYM(KA,KB)) THEN
                    SUM=UMAT
                 ELSE
                    SUM=0.D0
                 ENDIF
      RETURN
      END
