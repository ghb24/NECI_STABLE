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


      RECURSIVE SUBROUTINE GENSYMDETSSDN_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                NLIST,NLMAX,NJ,KJ,NELEC,NBF,NBASISMAX,BETA,I_P,
     &       NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,PP,RHOEPS)
         use SystemData, only: BasisFN,BasisFNSize
         use sym_mod
         use sort_mod
         use util_mod, only: NECI_ICOPY
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, NLMAX,NLIST
         INTEGER LSTE(NEL,NLMAX),NJ(NEL),NELEC,NBF
         INTEGER I,BRR(NBASIS),NN(NEL),nBasisMax(5,*)
         INTEGER IGETEXCITLEVEL,ICE
         real(dp) RN,RP,RH,PP,RH2,RHOEPS
         real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE
         INTEGER I_P,NMSH,NMAX,NTAY
         TYPE(BASISFN) G1(NBASIS),KI,KJ,KJ2
         DO I=NBF,NBASIS
C.. We continue if we're not searching for a specific det, or if we've
C.. yet to find a specific det.
            IF(PP.EQ.0.0_dp.OR.NLIST.LT.NLMAX) THEN
               NJ(NELEC)=BRR(I)
               KJ2=KJ
C.. Check if we've filled all the electrons
               IF(NELEC.EQ.NEL) THEN
                  CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
                  IF(LCHKSYM(KJ2,KI)) THEN
                   CALL NECI_ICOPY(NEL,NJ,1,NN,1)
                   call sort (nN)
                   ICE=IGETEXCITLEVEL(NI,NN,NEL)
C.. Check the connectivity
                   IF(ICE.GT.0.AND.ICE.LE.2) THEN
                     NLIST=NLIST+1
                     IF(RP.NE.0.0_dp) THEN
                        CALL CALCRHO2(NI,NN,BETA,I_P,NEL,
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &             RH2,NTAY,ICE,ECORE)
                        IF(RP.GE.0.0_dp) THEN
                           RN=RN+RH2**RP
                        ELSEIF(ABS(RH2).GT.RHOEPS) THEN
                           CALL CALCRHO2(NN,NN,BETA,I_P,NEL,
     &             G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,
     &             RH,NTAY,0,ECORE)
                           RN=RN+RH**ABS(RP)
                        ENDIF
                     ENDIF
                     IF(NLIST.LE.NLMAX) THEN
C.. if there's space, we save it
                        IF(PP.NE.0.0_dp) THEN
                           IF(RN.GE.PP) THEN
                              CALL NECI_ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
                           ELSE
                              NLIST=NLIST-1
                           ENDIF
                        ELSE
                           CALL NECI_ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
                        ENDIF
                     ENDIF
                   ENDIF
                  ENDIF
               ELSE
C.. otherwise we need to add more electrons:
                  CALL GENSYMDETSSDN_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &           NLIST,NLMAX,NJ,KJ2,NELEC+1,I+1,NBASISMAX,BETA,I_P,
     &         NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,PP,RHOEPS)
               ENDIF
            ENDIF
         ENDDO
         RETURN
      END

C.. Get the determinants of the same symmetry as KI, but which are only
C.. singles and doubles of NI.  If RP!=0, these are weighted:
C.. if ABS(RHO_IJ)<RHOEPS the weight is zero.
C.. Otherwise if (RP<0) the weight is RHO_JJ**ABS(RP).
C.. If RP>0, the weight is RHO_IJ**ABS(RP)
C.. RN is returned as the total of the weights of all the dets (if
C.. RP!=0.0_dp).  
C.. 
C.. If PP>0.0_dp, then as RN is summed, if a det's contribution to RN
C.. exceeds PP, then the procedure ends.
C..
C.. NLIST on input is set to the number of dets which can fit into 
C.. LSTE
C.. On exit, it contains the number of dets generated.
C.. If there are more dets generated than the input NLIST, then only
C.. the first NLIST are stored.

      SUBROUTINE GENSYMDETSSDN(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &         NLIST,NBASISMAX,
     &         BETA,I_P,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,
     &         PP,RHOEPS)
         use SystemData, only: BasisFN,BasisFNSize,Symmetry,NullBasisFn
         use constants, only: dp
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST)
         INTEGER NJ(NEL),nBasisMax(5,*)
         INTEGER NLMAX
         real(dp) RN,RP,PP,RHOEPS
         real(dp) BETA,FCK(*),ALAT(*),UMAT(*),ECORE
         INTEGER I_P,NMSH,NMAX,NTAY
         TYPE(BASISFN) G1(NBASIS),KI,KJ
         KJ=NullBasisFn
         NJ(1:NEL)=0
         NLMAX=NLIST
         NLIST=0
         RN=0.0_dp
         CALL GENSYMDETSSDN_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX,
     &         BETA,I_P,NMSH,FCK,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,
     &         PP,RHOEPS)
         RETURN
      END
C.. Get the next symmetric determinant in sequence.
C.. TDONE will be set to TRUE if there are no more dets.
C.. KI is the sym, NJ contains the working
      SUBROUTINE GENNEXTSYMDET(KI,NEL,G1,BRR,NBASIS,
     &            TDONE,NJ,NN,INDJ,NELEC,NBASISMAX,nActiveBasis)
         use SystemData, only: BasisFN
         use sort_mod
         use sym_mod
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL, NBASIS
         TYPE(BasisFN) G1(nBasis),kI,kJ2
         INTEGER NJ(NEL),NELEC,INDJ(NEL),NLAST
         INTEGER I,BRR(NBASIS),NN(NEL),nBasisMax(5,*)
         LOGICAL LCONT,TDONE
         INTEGER nActiveBasis(2)
         TDONE=.FALSE.
         LCONT=.TRUE.
         DO WHILE (LCONT)
C.. Remove the current electron in this slot
C.. Only do this if there was one there before.
C            IF((NELEC.GT.1.AND.INDJ(NELEC).GT.(INDJ(NELEC-1)))
C     &               .OR.(NELEC.EQ.1.AND.INDJ(1).GT.0)) THEN
C.. we don't even have to set it to zero, as it'll be overwritten
C            ENDIF
            NLAST=NActiveBasis(2)-(NEL-NELEC)
C.. if we've past the last allowable electron in our position
            IF(INDJ(NELEC).GE.NLAST) THEN
C.. Move back a slot
               NELEC=NELEC-1
C.. If we move back past the beginning
               IF(NELEC.LT.nActiveBasis(1)) THEN
C.. we're done
                  TDONE=.TRUE.
                  LCONT=.FALSE.
               ENDIF
            ELSE
C..Increment this slot
               I=INDJ(NELEC)+1
               INDJ(NELEC)=I
               IF(I.LE.nActiveBasis(2)) NJ(NELEC)=BRR(I)
C.. Check if we've filled all the electrons
               IF(NELEC.EQ.NEL) THEN
                  CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
                  IF(I.LE.nActiveBasis(2).AND.LCHKSYM(KI,KJ2)) THEN
                     CALL NECI_ICOPY(NEL,NJ,1,NN,1)
                     call sort (nN)
                     LCONT=.FALSE.
                  ENDIF
               ELSE
C.. Otherwise we need to go to the next slot
                  NELEC=NELEC+1
                  INDJ(NELEC)=I
               ENDIF
            ENDIF
         ENDDO
         RETURN
      END

C.. WORK(2*NEL+18)
      SUBROUTINE GENNEXTDET(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &      TPARITY,SymRestrict,ISYM,NI,TSETUP,TDONE,Work,nActiveBasis)
         use SystemData, only: BasisFN,BasisFNSize
         use gnd_work_type
         use sym_mod
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),nBasisMax(5,*)
         INTEGER LMS,NI(NEL)
C         INTEGER WORK(5*NEL+18)
         LOGICAL TSPN,TPARITY,TSETUP,TDONE
         TYPE(GNDWork) Work
         TYPE(BasisFN) G1(nBasis),ISym,SymRestrict
         INTEGER nActiveBasis(2)
C.. WORK(7*NEL+18)
C..  1-4                NSWORK
C..  5-14               IMAX
C..  15-19              ISYM
C..  20                 NELEC
C..  21-NEL+20          NIWORK
C..  NEL+21 - 2*NEL+20  NINDJ
         CALL GENNEXTDET_(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &      TPARITY,SymRestrict,NI,TSETUP,TDONE,Work,nActiveBasis)
CWORK(5),
C     &      WORK(15),WORK(20),WORK(21),WORK(NEL+21))
         ISYM = Work%ISym
         RETURN
      END

      SUBROUTINE GENNEXTDET_(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &      TPARITY,SymRestrict,NI,TSETUP,TDONE,Work,nActiveBasis)
C,IMAX,ISYM,NELEC,
C     &      NIWORK,NINDJ)
         use SystemData, only: BasisFN,BasisFNSize
         use gnd_work_type
         use sym_mod
         IMPLICIT NONE
         TYPE(GNDWork) Work 
         INTEGER NEL,NBASIS,BRR(NBASIS),nBasisMax(5,*)
         TYPE (BasisFN) G1(NBASIS),SymRestrict
         INTEGER LMS,NI(NEL)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE
         INTEGER nActiveBasis(2)
C..  1-4                NSWORK
C..  5-12               IMAX
C..  13-17              ISYM
C..  18                 NELEC
C..  19-NEL+18          NIWORK
C..  NEL+19 - 2*NEL+18  NINDJ
         INTEGER NINDJ(NEL)
         IF(TSETUP) THEN
            CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &         TPARITY,SymRestrict,TSETUP,TDONE,Work%IMax,Work%ISym)
C,IMAX,ISYM)
            Work%nElec=1
            Work%niWork(nEl+1)=0
            RETURN
         ENDIF
C.. Go to the next det.
         TMORE=.TRUE.
         DO WHILE(TMORE)
            CALL GENNEXTSYMDET(Work%ISYM,NEL,G1,BRR,NBASIS,
     &           TDONE,Work%niWork,NI,Work%niWork(nEl+1)
     &            ,Work%nElec,NBASISMAX,nActiveBasis)
            IF(TDONE) THEN
               CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &            TPARITY,SymRestrict,TSETUP,TMORE,Work%IMax,Work%iSym)
               TMORE=.NOT.TMORE
               Work%nElec=1
               NINDJ(1)=0
            ELSE
               TMORE=.FALSE.
            ENDIF
         ENDDO
      END
