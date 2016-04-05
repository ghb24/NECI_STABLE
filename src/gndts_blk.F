! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
      SUBROUTINE GNDTS_BLK(NEL,NBASIS,BRR,NBASISMAX,NMRKS,TCOUNT,
     &         NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS,
     &         TPARITY,SymRestrict,IFDET,TGENFDET,NDETTOT,BLOCKSYM
     &         )
         use SystemData, only: BasisFN,BasisFNSize
         use sym_mod
         use util_mod, only: NECI_ICOPY
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),nBasisMax(5,*),NDET
         INTEGER NMRKS(NEL,NDET)
         INTEGER II,NLMAX
         INTEGER NBLOCKS,NBLOCKSTARTS(NBLOCKS+1),OII
         TYPE(BASISFN) G1(NBASIS),ISYM,KJ,BLOCKSYM(NBLOCKS),IMAX(2)
         TYPE(BasisFN) SymRestrict
         LOGICAL TCOUNT
         INTEGER I,NI(NEL),NJ(NEL)
         LOGICAL TSPN,TPARITY,TDONE
         INTEGER LMS
         LOGICAL TGENFDET
         INTEGER IFDET,NDETTOT,IDEG
         real(dp) DETSC,TDETSC
         real(dp) CALCT
         DETSC=1D200
         II=0
         IF(TCOUNT) THEN
            NLMAX=0
            NBLOCKS=0
         ELSE
            NLMAX=NDET
         ENDIF
         I=0
         OII=0
         NDETTOT=0
C.. set the comparison det to an invalid one, so all dets are counted
         NI(1)=0
         IF(TCOUNT) OPEN(14,FILE="BLOCKS",STATUS="UNKNOWN")
         CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &      TPARITY,SymRestrict,.TRUE.,TDONE,IMAX,ISYM)
         DO WHILE(.NOT.TDONE)
            CALL SETUPSYM(KJ)
            CALL  GENSYMDETS_R(NI,ISYM,NEL,G1,BRR,NBASIS,NMRKS,
     &            II,NLMAX,NJ,KJ,1,1,NBASISMAX)
            IF(II.NE.OII) THEN
C.. we've found an occupied block 
               I=I+1
               CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
               NDETTOT=NDETTOT+(II-OII)*IDEG
               IF(TCOUNT) THEN
                  WRITE(14,"(I5)",advance='no') I
                  CALL WRITEALLSYM(14,ISYM,.FALSE.)
                  WRITE(14,"(2I10)")
     &               II-OII,
     &               (II-OII)*IDEG
               ELSE
                  NBLOCKSTARTS(I)=OII+1
                  BLOCKSYM(I) = ISYM
                  IF(TGENFDET) THEN
                     TDETSC=CALCT(NMRKS(1:NEL,OII+1),NEL)
                     IF(TDETSC.LT.DETSC) THEN
                        IFDET=OII+1
                        DETSC=TDETSC
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            OII=II
            CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &         TPARITY,SymRestrict,.FALSE.,TDONE,IMAX,ISYM)
         ENDDO
         NBLOCKS=I
         IF(.NOT.TCOUNT) NBLOCKSTARTS(I+1)=II+1
         IF(TCOUNT) CLOSE(14)
         IF(.NOT.TGENFDET) IFDET=1
         RETURN
      END         

 

      RECURSIVE SUBROUTINE GENSYMDETSSD_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                     NLIST,NLMAX,NJ,KJ,NELEC,NBF,NBASISMAX)
         use SystemData, only: BasisFN,BasisFNSize
         use sort_mod
         use sym_mod
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, NLMAX,NLIST
         TYPE(BASISFN) G1(NBASIS),KI,KJ,KJ2
         INTEGER LSTE(NEL,NLMAX),NJ(NEL),NELEC,NBF
         INTEGER I,BRR(NBASIS),NN(NEL),nBasisMax(5,*)
         INTEGER IGETEXCITLEVEL,ICE
         DO I=NBF,NBASIS
            NJ(NELEC)=BRR(I)
            KJ2=KJ
C.. Check if we've filled all the electrons
            IF(NELEC.EQ.NEL) THEN
               CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
               CALL ROUNDSYM(KJ2,NBASISMAX)
               IF(LCHKSYM(KJ2,KI)) THEN
                CALL NECI_ICOPY(NEL,NJ,1,NN,1)
                call sort (nN)
                ICE=IGETEXCITLEVEL(NI,NN,NEL)
                IF(ICE.GT.0.AND.ICE.LE.2) THEN
C.. we've found a det with the right sym.
                  NLIST=NLIST+1
                  IF(NLIST.LE.NLMAX) THEN
C.. if there's space, we save it
                     CALL NECI_ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
                  ENDIF
                ENDIF
               ENDIF
            ELSE
C.. otherwise we need to add more electrons:
               CALL GENSYMDETSSD_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                  NLIST,NLMAX,NJ,KJ2,NELEC+1,I+1,NBASISMAX)
            ENDIF
         ENDDO
         RETURN
      END

C.. Get the determinants of the same symmetry as KI, but which are only
C.. singles and doubles of NI
      SUBROUTINE GENSYMDETSSD(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &         NLIST,NBASISMAX)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry,NullBasisFn
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST)
         TYPE(BASISFN) G1(NBASIS),KI,KJ
         INTEGER NJ(NEL),nBasisMax(5,*)
         INTEGER NLMAX
         KJ=NullBasisFn
         NJ(1:NEL)=0
         NLMAX=NLIST
         NLIST=0
         CALL GENSYMDETSSD_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX)
         RETURN
      END
 

C.. Generate determinants of the same symmetry as NI (including NI).
C.. As they are generated, they are placed in LSTE, which is of length
C.. NLIST, up until LSTE is full.  NLIST at the end contains the total
C.. number of determinants available (even if LSTE doesn't hold them all).
      SUBROUTINE GENSYMDETS(NI,NEL,G1,BRR,NBASIS,LSTE,NLIST,NBASISMAX)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry,NullBasisFn
         use sym_mod, only: getsym
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST),nBasisMax(5,*)
         TYPE(BASISFN) G1(NBASIS),KI,KJ
         INTEGER NJ(NEL)
         INTEGER NLMAX
         KJ=NullBasisFn
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KI)
         NLMAX=NLIST
         NLIST=0
         CALL GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX)
         RETURN
      END

      RECURSIVE SUBROUTINE GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                     NLIST,NLMAX,NJ,KJ,NELEC,NBF,NBASISMAX)
         use SystemData, only: BasisFN,BasisFNSize
         use sort_mod
         use sym_mod
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, NLMAX,NLIST
         INTEGER LSTE(NEL,NLMAX),NJ(NEL),NELEC,NBF
         INTEGER I,J,BRR(NBASIS),NN(NEL),nBasisMax(5,*)
         LOGICAL LCMP
         TYPE(BASISFN) G1(NBASIS),KI,KJ,KJ2
         DO I=NBF,NBASIS
            NJ(NELEC)=BRR(I)
            KJ2=KJ
C.. Check if we've filled all the electrons
            IF(NELEC.EQ.NEL) THEN
               CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
               IF(LCHKSYM(KJ2,KI)) THEN
C                  CALL ROUNDSYM(KJ3,NBASISMAX)
C.. we've found a det with the right sym.
                  NLIST=NLIST+1
                  CALL NECI_ICOPY(NEL,NJ,1,NN,1)
                  call sort (nN)
C.. Just check to see if it's our original det
                  LCMP=.TRUE.
                  DO J=1,NEL
                     IF(NN(J).NE.NI(J)) LCMP=.FALSE.
                  ENDDO
                  IF(LCMP) THEN
C.. Roll back NLIST
                     NLIST=NLIST-1
                  ELSEIF(NLIST.LE.NLMAX) THEN
C.. if there's space, we save it
                     CALL NECI_ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
                  ENDIF
               ENDIF
            ELSE
C.. otherwise we need to add more electrons:
               CALL GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                  NLIST,NLMAX,NJ,KJ2,NELEC+1,I+1,NBASISMAX)
            ENDIF
         ENDDO
         RETURN
      END

      
C.. Generate determinants with a given symmetry, given by KI, as GENSYMDETS
      SUBROUTINE GENSYMDETSS(KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,NBASISMAX)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry,NullBasisFn
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST)
         TYPE(BASISFN) G1(NBASIS),KI,KJ
         INTEGER NJ(NEL),nBasisMax(5,*)
         INTEGER NLMAX
         KJ=NullBasisFn
         NI(1:NEL)=0
         NLMAX=NLIST
         NLIST=0
         CALL GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX)
         RETURN
      END
