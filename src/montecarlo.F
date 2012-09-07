      INTEGER FUNCTION IFINDBASISFN(NQNS,G1,NBASIS)
         use SystemData, only: BasisFN
!         IMPLICIT real(dp) (A-H,O-Z)
         IMPLICIT NONE
         INTEGER NBASIS,I,J
         TYPE(BasisFN) G1(NBASIS),NQNS
         LOGICAL L
         DO I=1,NBASIS
            L=.TRUE.
            DO J=1,3
               IF(G1(I)%k(J).NE.NQNS%k(J)) L=.FALSE.
            ENDDO
            IF(G1(I)%Sym%s.NE.NQNS%Sym%s) L=.FALSE.
            IF(G1(I)%Ms.NE.NQNS%Ms) L=.FALSE.
            IF(L) THEN
               IFINDBASISFN=I
               RETURN
            ENDIF
         ENDDO
         IFINDBASISFN=0
         RETURN
      END

      SUBROUTINE GENNEWDET2(NI,NJ,NEL,NBASIS,NBASISMAX,G1,ISEED,NCH,BRR)
         use SystemData, only: BasisFN
         use sort_mod
         use util_mod, only: NECI_ICOPY
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NBASIS,nBasisMax(5,*),NI(NEL),NJ(NEL),NQN(5)
         TYPE(BasisFN) G1(NBASIS)
         INTEGER NCH,ND,BRR(NBASIS),ISEED
         LOGICAL L
         INTEGER IEL,IQN,ISGN,I,IT,IFINDBASISFN
         real(dp) RAN2
         IF(NBASISMAX(3,3).EQ.1) THEN
C.. Generic spatial basis
         call stop_all('GENNEWDET2',
     &         'G1 variable needs to be made into a BasisFN.')
!            CALL GENNEWDET3(NI,NJ,NEL,NBASIS,G1,ISEED,NCH,BRR)
            RETURN
         ENDIF
C.. Copy NI to NJ
         CALL NECI_ICOPY(NEL,NI,1,NJ,1)
         L=.TRUE.
         ND=NCH
C.. we need to keep on doing this until we make a valid choice
         DO WHILE(L)
            IEL=INT(RAN2(ISEED)*NEL)+1
            !CALL NECI_ICOPY(4,G1(NJ(IEL)),1,NQN,1)
            call stop_all('GENNEWDET2',
     &                  'NQN variable needs to be made into a BasisFN.')
            IQN=3+(NBASISMAX(4,2)-NBASISMAX(4,1))/2
C.. Sometimes no spin freedom
            IQN=INT(RAN2(ISEED)*IQN)+1
            ISGN=INT(RAN2(ISEED)*(NBASISMAX(IQN,2)-NBASISMAX(IQN,1)+1))
     &         +NBASISMAX(IQN,1)
            IF(IQN.EQ.4) THEN
C.. change the spin
               NQN(4)=-NQN(4)
            ELSE
               NQN(IQN)=NQN(IQN)+ISGN
            ENDIF
            CALL HUB_PBC(NQN,NBASISMAX)
!            IT=IFINDBASISFN(NQN,G1,NBASIS)
            IT=0
            IF(IT.EQ.0) THEN 
               L=.TRUE.
            ELSE
               L=.FALSE.
            ENDIF
C.. make sure the new basis fn isn't already in the det.  if so, go round again
            DO I=1,NEL
               IF(NJ(I).EQ.IT) L=.TRUE.
            ENDDO
            IF(.NOT.L) THEN
               ND=ND-1
               NJ(IEL)=IT
               call sort (nJ)
C.. if we've more electrons to change, we go round again
               IF(ND.GT.0) L=.TRUE.
            ENDIF
         ENDDO
         RETURN
      END

      SUBROUTINE HUB_PBC(NQN, NBASISMAX)
         use constants, only: dp,sizeof_int
         IMPLICIT NONE
         INTEGER NQN(5), nBasisMax(5,*)
         INTEGER MPN,MMN,IQN,NMAXX,MX,MY,I,J,KX,KY,NMINX
         real(dp) :: XX,YY,AX,AY
         LOGICAL TALLOW
         IF(NBASISMAX(1,3).EQ.0) THEN
C.. PBC non-tilted
            DO IQN=1,3
               IF(NQN(IQN).GT.NBASISMAX(IQN,2))
     &            NQN(IQN)=NBASISMAX(IQN,1)
               IF(NQN(IQN).LT.NBASISMAX(IQN,1))
     &            NQN(IQN)=NBASISMAX(IQN,2)
            ENDDO
         ELSEIF(NBASISMAX(1,3).EQ.1) THEN
C.. PBC tilted
            TALLOW=.FALSE.
            DO WHILE(.NOT.TALLOW)
               TALLOW=.TRUE.
               NMAXX=NBASISMAX(1,5)
               AX=NBASISMAX(1,4)
               AY=NBASISMAX(2,4)
C.. (XX,YY) is the position of the bottom right corner of the unit cell
               XX=((AX+AY)/2.D0)*real(NMAXX,dp)
               YY=((AY-AX)/2.D0)*real(NMAXX,dp)
               MX=int(XX*AX+YY*AY,sizeof_int)
               MY=int(XX*AY-YY*AX,sizeof_int)
               I=NQN(1)
               J=NQN(2)
               KX=int(real(I,dp)*AX+real(J,dp)*AY,sizeof_int)
               KY=int(real(I,dp)*AY-real(J,dp)*AX,sizeof_int)
               IF(KX.GT.MX) THEN
                  TALLOW=.FALSE.
                  NQN(1)=NQN(1)-int(AX*real(NMAXX,dp),sizeof_int)
                  NQN(2)=NQN(2)-int(AY*real(NMAXX,dp),sizeof_int)
               ENDIF
               IF(KY.GT.MY) THEN
                  TALLOW=.FALSE.
                  NQN(1)=NQN(1)-int(AY*real(NMAXX,dp),sizeof_int)
                  NQN(2)=NQN(2)+int(AX*real(NMAXX,dp),sizeof_int)
               ENDIF
               IF(KX.LE.-MX) THEN
                  TALLOW=.FALSE.
                  NQN(1)=NQN(1)+int(AX*real(NMAXX,dp),sizeof_int)
                  NQN(2)=NQN(2)+int(AY*real(NMAXX,dp),sizeof_int)
               ENDIF
               IF(KY.LE.-MY) THEN
                  TALLOW=.FALSE.
                  NQN(1)=NQN(1)+int(AY*real(NMAXX,dp),sizeof_int)
                  NQN(2)=NQN(2)-int(AX*real(NMAXX,dp),sizeof_int)
               ENDIF
            ENDDO
         ELSEIF(NBASISMAX(1,3).EQ.3) THEN
C.. non-PBC tilted
            NMAXX=NBASISMAX(1,2)
            NMINX=-NMAXX
            MPN=NQN(1)+NQN(2)
            MMN=NQN(1)-NQN(2)
            NQN(1)=(MPN+MMN)/2
            NQN(2)=(MPN-MMN)/2
         ELSEIF(NBASISMAX(1,3).EQ.4) THEN            
C.. non-PBC non-tilted - leave things as they are
         ENDIF
         RETURN
      END

C.. gen a new det from this one in a generic spatial basis.
C.. We do this by exciting an electron into a spinorbital of close energy.
C.. This does not take into account any symmetries.
      SUBROUTINE GENNEWDET3(NI,NJ,NEL,NBASIS,G1,ISEED,NCH,BRR)
        use sort_mod
        use constants, only: dp
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL,NBASIS,NI(NEL),NJ(NEL),G1(*)
         INTEGER NCH,ND,BRR(NBASIS),ISEED
         LOGICAL L
         real(dp) RN
         INTEGER IGETENEIGHBOUR,IGETSNEIGHBOUR
         INTEGER IEL,ISGN,IT,I,IDIST
         real(dp) RAN2
C.. Copy NI to NJ
         CALL NECI_ICOPY(NEL,NI,1,NJ,1)
         L=.TRUE.
         ND=NCH
C.. we need to keep on doing this until we make a valid choice
C.. We keep spin fixed here.
        call stop_all('GENNEWDET3',
     &         'G1 variable needs to be made into a BasisFN.')
         IT=0
         DO WHILE(L)
            L=.FALSE.
            IEL=INT(RAN2(ISEED)*NEL)+1
            ISGN=INT(RAN2(ISEED)+0.5)*2-1
            IF(RAN2(ISEED).LT.0.25) THEN
C.. 1/4 prob of changing the spin
!               IT=IGETSNEIGHBOUR(NI(IEL),NBASIS,BRR,ISGN,G1)
            ELSE
C.. change the energy level.  Have an exponential distribution of
C.. distance
               RN=0.D0
               DO WHILE(RN.EQ.0.D0)
                  RN=RAN2(ISEED)
               ENDDO
               IDIST=0
               DO WHILE(RN.LE.1)
                  RN=RN*2
                  IDIST=IDIST+1
               ENDDO
!               IT=IGETENEIGHBOUR(NI(IEL),NBASIS,BRR,ISGN,G1,IDIST)
            ENDIF
C.. make sure the new basis fn isn't already in the det.  if so, go round again
            DO I=1,NEL
               IF(NJ(I).EQ.IT) L=.TRUE.
            ENDDO
            IF(.NOT.L) THEN
               ND=ND-1
               NJ(IEL)=IT
               call sort(nJ(1:nel))
C.. if we've more electrons to change, we go round again
               IF(ND.GT.0) L=.TRUE.
            ENDIF
         ENDDO
         RETURN
      END

      SUBROUTINE CLEARCACHE(LSTDET,NEL,NFIRST)
         INTEGER LSTDET(NEL+2,*)
C.. Clear the first node, and end the chain
         LSTDET(NEL+2,NFIRST)=0
         LSTDET(1,NFIRST)=0
         RETURN
      END

      INTEGER FUNCTION IGETSNEIGHBOUR(IEL,NBASIS,BRR,ISGN,G1)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER IEL,NBASIS,BRR(NBASIS),ISGN
         TYPE(BasisFN) G1(NBASIS)
         INTEGER I,J
         J=0
         DO I=1,NBASIS
            IF(BRR(I).EQ.IEL) J=I
         ENDDO
C.. make sure we have an orbital of the correct spin
C.. See if there's an orbital with the same energy
         DO WHILE(G1(BRR(J))%Ms.EQ.G1(IEL)%Ms.OR.IEL.EQ.BRR(J))
            J=J+ISGN
            IF(J.LT.1) J=NBASIS
            IF(J.GT.NBASIS) J=1
         ENDDO 
         IGETSNEIGHBOUR=BRR(J)  
         RETURN
      END


      INTEGER FUNCTION IGETENEIGHBOUR(IEL,NBASIS,BRR,ISGN,G1,IDIST)
         use SystemData, only: BasisFN
         IMPLICIT NONE
         INTEGER IEL,NBASIS,BRR(NBASIS),ISGN
         TYPE(BasisFN) G1(NBASIS)
         INTEGER I,J,IDIST,IDIST2
         J=0
         DO I=1,NBASIS
            IF(BRR(I).EQ.IEL) J=I
         ENDDO
         IDIST2=IDIST
C.. make sure we have an orbital of the correct spin
         DO WHILE(IDIST2.NE.0)
            J=J+ISGN
            IF(J.LT.1) J=NBASIS
            IF(J.GT.NBASIS) J=1
            IF(G1(BRR(J))%Ms.EQ.G1(IEL)%Ms) IDIST2=IDIST2-1
         ENDDO 
         IGETENEIGHBOUR=BRR(J)  
         RETURN
      END

