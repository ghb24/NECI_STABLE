! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
C.. GETHELEMENT2T
C.. Get matrix element of the hamiltonian
C.. IC is the number of basis fns that differ in NI and NJ (or -1 if not
Cknown)
C.. ECORE is the uniform background energy.  It has been renamed from
C.. ETRIAL, but ETRIAL not longer has a function
      FUNCTION GETHELEMENT2T(NI,NJ,NEL,NBASISMAX,
     &   NBASIS,ECORE,IC2,BTRANS,NSBASIS,NSPINS)
         use constants, only: dp
         IMPLICIT NONE      
         INTEGER NBASIS,NSBASIS,NSPINS
         INTEGER NEL,NI(NEL),NJ(NEL),IC,nBasisMax(5,*),IC2
         real(dp) TOTSUM,ECORE
         INTEGER IGETEXCITLEVEL
         real(dp) BTRANS(*),GETHELEMENT2T
         IC=IC2
         GETHELEMENT2T=0.0_dp
         IF(IC.LT.0) IC=IGETEXCITLEVEL(NI,NJ,NEL)
C.. if we differ by more than 2 spin orbital, then the hamiltonian
Celement is 0         
         IF(IC.GT.2) RETURN
C.. SLTCND has IC is # electrons the same in 2 dets
         CALL SLTCNDT(NEL,NBASISMAX,
     &         NBASIS,NI,NJ,NEL-IC,BTRANS,NSBASIS,
     &         NSPINS,TOTSUM)
         GETHELEMENT2T=TOTSUM
         IF(IC.EQ.0) GETHELEMENT2T=GETHELEMENT2T+ECORE
         RETURN
      END

C.. The Slater Condon Rules for a transformed basis
      SUBROUTINE SLTCNDT(NEL,NBASISMAX,
     &   NHG,NDET1,NDET2,IC,
     &   BTRANS,NSBASIS,NSPINS,TOTSUM)
         use constants, only: dp
      IMPLICIT NONE
C..
C.. 
      INTEGER IC,NEL,NHG
      INTEGER NDET1(NEL),NDET2(NEL)
      INTEGER N1E(NHG),N2E(NHG)
      INTEGER nBasisMax(5,*),IFLAG,IPTOT
      real(dp) COULCOUPLE,BTRANS(*),TOTSUM
      INTEGER ISPINSKIP,NSBASIS,NSPINS
      ISPINSKIP=NBASISMAX(2,3)
C..
C..Routine to read in UMAT
C..Need to get the determinants in the form of max. coincidence
      N1E(1:NHG)=0
      N2E(1:NHG)=0
C..Also returns the sign change IPTOT
      CALL LINEUP(NEL,NDET1,NDET2,NHG,N1E,N2E,IFLAG,IPTOT)
      TOTSUM=0.0_dp
C..We now check to see which Slater-Condon rule to apply
      COULCOUPLE=1.0_dp
      IF(IC.EQ.(NEL-2)) THEN
C..Differ by two 
        CALL SCR2T(NEL,NDET1,NDET2,
     &         COULCOUPLE,BTRANS,
     &   NSBASIS,NSPINS,TOTSUM)
      ELSEIF(IC.EQ.(NEL-1)) THEN
C..Differ by one 
        CALL SCR1T(NEL,NDET1,NDET2,COULCOUPLE,BTRANS,NSBASIS,NSPINS,
     &             TOTSUM)
      ELSEIF(IC.EQ.NEL) THEN
C..Differ by none 
        CALL SCR0T(NEL,NDET1,COULCOUPLE,BTRANS,NSBASIS,
     &             NSPINS,TOTSUM)
C.. If we're in the UEG, we need to add in the part of the background
C.. contribution which doesn't cancel
C         IF(NBASISMAX(1,3).EQ.-1) THEN
C            WRITE(6,*) FCK(0,0,0)
C            TOTSUM=TOTSUM-NEL*REAL(FCK(0,0,0))/2

C         ENDIF
      ELSE
        TOTSUM=0.0_dp
        RETURN
      ENDIF
      TOTSUM=TOTSUM*real(IPTOT,dp)
      RETURN
      END

C.. Slater-Condon Rules in a transformed basis.   based on SCR0 and SCR1
C.. TMAT has sometimes been hidden in ZIA
      SUBROUTINE SCR0T(NEL,NDET1,
     &   COULCOUPLE,BTRANS,
     &   NSBASIS,NSPINS,TOTSUM)
         use constants, only: dp
         IMPLICIT NONE
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         INTEGER NEL,NSBASIS,NSPINS
         INTEGER NDET1(NEL)
         real(dp) TRUMAT,TRTMAT
         real(dp) COULCOUPLE,TOTSUM
         INTEGER I,J
         real(dp) SUM1,SUM2,PI

C.. Basis transformation matrix.  
         real(dp) BTRANS(NSBASIS,NSBASIS,NSPINS)
         PI=ACOS(-1.0_dp)
C..This is the Slater-Condon rule when the two determinants are the
C*SAME*
         SUM1=0.0_dp
         SUM2=0.0_dp
C..We must first do the one electron integrals
         DO I=1,NEL
C.. we extract the KEs from TMAT
            CALL GETTRTMATEL(NDET1(I),NDET1(I),BTRANS,NSBASIS,
     &               NSPINS,TRTMAT)
            SUM1=SUM1+TRTMAT
         ENDDO
C..Now for the two electron part
         DO I=1,NEL-1
            DO J=I+1,NEL
               CALL GETTRUMATEL(NDET1(I),NDET1(J),NDET1(I),NDET1(J),
     &            BTRANS,NSBASIS,NSPINS,TRUMAT)
C.. the spin checking occurs in GETTRUMATEL
               SUM2=SUM2+TRUMAT
C.. and also the exchange
               CALL GETTRUMATEL(NDET1(I),NDET1(J),NDET1(J),NDET1(I),
     &            BTRANS,NSBASIS,NSPINS,TRUMAT)
               SUM2=SUM2-TRUMAT
CUMAT(ID(I),ID(J),ID(J),ID(I))
           ENDDO
         ENDDO
C..
         TOTSUM=SUM1+SUM2*COULCOUPLE
         RETURN
      END



      SUBROUTINE SCR1T(NEL,NDET1,NDET2,
     &   COULCOUPLE,BTRANS,
     &   NSBASIS,NSPINS,TOTSUM)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL
         INTEGER NSBASIS,NSPINS
         INTEGER NDET1(NEL),NDET2(NEL)
         INTEGER ID(NEL-1)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         real(dp) TRTMAT,TRUMAT
         real(dp) COULCOUPLE,TOTSUM,SUM2
         INTEGER ND1,ND2,IQ,I,J,K
         real(dp) BTRANS(NSBASIS,NSBASIS,NSPINS)
C.. NB in HF we need to include TMAT element

C..In this routine the determinants differ by only 1 spin orbital
C..The integers ND1 and ND2 store the spin orbitals that are different 
C..in each determinant
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
100      CONTINUE
C.. ND1 is the basis fn which differs in NDET1 to NDET2
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
200      CONTINUE 
C.. ND2 is the different basis fn in NDET2
C..The two electron integrals are to be done 
         IQ=0
         DO I=1,NEL
           IF(NDET1(I).NE.ND1) THEN
             IQ=IQ+1
C          ID(IQ) = GTID(NDET1(I))
             ID(IQ)=NDET1(I)
             IF(IQ.GE.(NEL-1)) GOTO 500
           ENDIF 
         ENDDO 
500      CONTINUE
         IF(ND1.EQ.ND2) STOP ' !!! PROBLEM IN SCR1T !!! '
C         ID1 = GTID(ND1)
C         ID2 = GTID(ND2)
C..
         SUM2=0.0_dp
         DO I=1,NEL-1
C.. non-transformed version
C        SUM2=SUM2+(UMAT(ID1,ID(I),ID2,ID(I))*DFLOAT(IDS)-
C     &      DFLOAT(IDS1(I)*IDS2(I))*UMAT(ID1,ID(I),ID(I),ID2))

            CALL GETTRUMATEL(ND1,ID(I),ND2,ID(I),BTRANS,NSBASIS,
     &                       NSPINS,TRUMAT)
C.. the spin checking occurs in GETTRUMATEL
            SUM2=SUM2+TRUMAT
C.. and also the exchange
            CALL GETTRUMATEL(ND1,ID(I),ID(I),ND2,
     &         BTRANS,NSBASIS,NSPINS,TRUMAT)

             SUM2=SUM2-TRUMAT
C            SUM2=SUM2-UMAT(ID(I),ID(J),ID(J),ID(I))
         ENDDO
C..  We then add in the non-diagonal part of the kinetic energy -
C..  <ph_a|T|ph_a'> where a and a' are the only basis fns that differ in
C..  NDET1 and NDET2
         CALL GETTRTMATEL(ND1,ND2,BTRANS,NSBASIS,
     &            NSPINS,TRTMAT)
         TOTSUM=SUM2*COULCOUPLE+TRTMAT
         RETURN
      END

      SUBROUTINE SCR2T(NEL,NDET1,NDET2,
     &     COULCOUPLE,BTRANS,
     &   NSBASIS,NSPINS,TOTSUM)
         use constants, only: dp
C..This routine is the EASIEST of the Slater-Condon rules.
C..The determinants differ by two spin orbitals.
C..The arrays ND1 and ND2 contain the integers from determinant 1 and 2
C..that are NOT common to both.  
         IMPLICIT NONE
         INTEGER NEL
         INTEGER NSBASIS,NSPINS
         INTEGER NDET1(NEL)
         INTEGER NDET2(NEL)
         INTEGER ND1(2),ND2(2)
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         real(dp) COULCOUPLE,TOTSUM,SUM2,TRUMAT
         real(dp) BTRANS(NSBASIS,NSBASIS,NSPINS)
         INTEGER I,J,K,IQ
         TOTSUM=0.0_dp
C..Here we zero ND1 and ND2
         DO I=1,2
            ND1(I)=0
            ND2(I)=0
         ENDDO
C..IQ ensures that only two differences are recorded
C..This first routine compares NDET1 to NDET2
         IQ=0
         DO I=1,NEL
            K=0
            DO J=1,NEL
               IF(NDET1(I).EQ.NDET2(J)) CONTINUE
               IF(NDET1(I).NE.NDET2(J)) K=K+1
            ENDDO
            IF(K.EQ.NEL) THEN
               IQ=IQ+1
               ND1(IQ)=NDET1(I)
               IF(IQ.GE.2) GOTO 100
            ENDIF
         ENDDO
 100     CONTINUE
C..This routine compares NDET2 to NDET1
         IQ=0
         DO I=1,NEL
            K=0
            DO J=1,NEL
               IF(NDET2(I).EQ.NDET1(J)) CONTINUE
               IF(NDET2(I).NE.NDET1(J)) K=K+1
            ENDDO
            IF(K.EQ.NEL) THEN
               IQ=IQ+1
               ND2(IQ)=NDET2(I)
               IF(IQ.GE.2) GOTO 200
            ENDIF
         ENDDO
 200     CONTINUE
C      ID2 = GTID(ND1(1))
C      ID1 = GTID(ND1(2))
C..
C      ID4 = GTID(ND2(1))
C      ID3 = GTID(ND2(2))
C..Check the spin combinations of N', N'-1, N and N-1 
C      IDS(1:4)=0
C      CALL GENKRON(G1,NBASISMAX,ND1,ND2,IDS)
C..N'>N'-1 and N>N-1 *See Levine for rule 3*
C      TOTSUM=UMAT(ID1,ID2,ID3,ID4)*DFLOAT(IDS(1)*IDS(2))-
C     &      UMAT(ID1,ID2,ID4,ID3)*DFLOAT(IDS(3)*IDS(4))
C      IF(TOTSUM.NE.0.0_dp) CALL CHK_MS(NDET1,NDET2,NBASISMAX,
C     &   NEL,3,G)
      
         CALL GETTRUMATEL(ND1(2),ND1(1),ND2(2),ND2(1),
     &         BTRANS,NSBASIS,NSPINS,TRUMAT)
C.. the spin checking occurs in GETTRUMATEL
         SUM2=TOTSUM+TRUMAT
C.. and also the exchange
         CALL GETTRUMATEL(ND1(2),ND1(1),ND2(1),ND2(2),
     &         BTRANS,NSBASIS,NSPINS,TRUMAT)
         TOTSUM=TOTSUM-TRUMAT
C         SUM2=SUM2-UMAT(ID(I),ID(J),ID(J),ID(I))
         TOTSUM=TOTSUM*COULCOUPLE
         RETURN
      END 

C.. Calculate <phi_a phi_b | U | phi_c phi_d> from <u_i u_j|U|u_k u_l>
      SUBROUTINE GETTRUMATEL(A,B,C,D,BTRANS,NSBASIS,NSPINS,TOTSUM)
         use constants, only: dp
         use procedure_pointers, only: get_umat_el
         USE UMatCache , Only : GTID
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER A,B,C,D
         INTEGER I,J,K,L
C.. was (NHG/ISS,NHG/ISS,NHG/ISS,NHG/ISS)
         INTEGER NSPINS,NSBASIS
         real(dp) BTRANS(NSBASIS,NSBASIS,NSPINS)
         HElement_t TOTSUM,SUM1,SUM2,SUM3
         INTEGER ASPN,BSPN,CSPN,DSPN,AE,BE,CE,DE
         INTEGER IDI,IDJ,IDK,IDL
         ASPN=MOD(A-1,NSPINS)+1
         BSPN=MOD(B-1,NSPINS)+1
         CSPN=MOD(C-1,NSPINS)+1
         DSPN=MOD(D-1,NSPINS)+1
C.. AE is the index if new basis fn A within the ASPN part of the BTRANS
C..  matrix
         AE=(A-1)/NSPINS+1
         BE=(B-1)/NSPINS+1
         CE=(C-1)/NSPINS+1
         DE=(D-1)/NSPINS+1
         IF(ASPN.NE.CSPN.AND.BSPN.NE.DSPN) THEN
            TOTSUM=0.0_dp
            RETURN
         ENDIF
         TOTSUM=0.0_dp
         DO I=1,NSBASIS
C.. We have to generate the normal basis number corresponding to I
C.. and then from it get the ID within UMAT
            IDI = GTID((I-1)*NSPINS+1+ASPN-1)
            SUM1=0.0_dp
            DO J=1,NSBASIS
               IDJ = GTID((J-1)*NSPINS+1+BSPN-1)
               SUM2=0.0_dp
               DO K=1,NSBASIS
                  IDK = GTID((K-1)*NSPINS+1+CSPN-1)
                  SUM3=0.0_dp
                  DO L=1,NSBASIS
                     IDL = GTID((L-1)*NSPINS+1+DSPN-1)
                     SUM3=SUM3+(BTRANS(DE,L,DSPN))*
     &       get_umat_el(IDI,IDJ,IDK,IDL)
C     &                  UMAT(IDI,IDJ,IDK,IDL)
                  ENDDO
                  SUM2=SUM2+SUM3*(BTRANS(CE,K,CSPN))
               ENDDO
               SUM1=SUM1+SUM2*(BTRANS(BE,J,BSPN))
            ENDDO
            TOTSUM=TOTSUM+SUM1*(BTRANS(AE,I,ASPN))
         ENDDO
         RETURN
      END

C.. Calculate <phi_a | T | phi_c > from <u_i |T|u_k>
      SUBROUTINE GETTRTMATEL(A,C,BTRANS,NSBASIS,
     &            NSPINS,TOTSUM)
         use constants, only: dp
         USE OneEInts, Only : GetTMatEl
         IMPLICIT NONE
         INTEGER A,C
         INTEGER I,K
         INTEGER NSPINS,NSBASIS
         real(dp) BTRANS(NSBASIS,NSBASIS,NSPINS)
         real(dp) TOTSUM,SUM2
         INTEGER ASPN,CSPN,AE,CE
         INTEGER IDI,IDK
         ASPN=MOD(A-1,NSPINS)+1
         CSPN=MOD(C-1,NSPINS)+1
C.. AE is the index if new basis fn A within the ASPN part of the BTRANS
C..  matrix
         AE=(A-1)/NSPINS+1
         CE=(C-1)/NSPINS+1
         IF(ASPN.NE.CSPN) THEN
            TOTSUM=0.0_dp
            RETURN
         ENDIF
         TOTSUM=0.0_dp
         DO I=1,NSBASIS
C.. We have to generate the normal basis number corresponding to I
            IDI=(I-1)*NSPINS+1+ASPN-1
            SUM2=0.0_dp
            DO K=1,NSBASIS
               IDK=(K-1)*NSPINS+1+CSPN-1
               SUM2=SUM2+BTRANS(CE,K,CSPN)*
     &            (GetTMATEl(IDI,IDK))
            ENDDO
            TOTSUM=TOTSUM+SUM2*BTRANS(AE,I,ASPN)
         ENDDO
         RETURN
      END
