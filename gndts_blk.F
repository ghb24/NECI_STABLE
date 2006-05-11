      SUBROUTINE GNDTS_BLK(NEL,NBASIS,BRR,NBASISMAX,NMRKS,TCOUNT,
     &         NDET,G1,II,NBLOCKSTARTS,NBLOCKS,TSPN,LMS,
     &         TPARITY,IPARITY,IFDET,TMAT,TGENFDET,NDETTOT,BLOCKSYM
     &         )
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),NBASISMAX(5,2),NDET
         INTEGER NMRKS(NEL,NDET)
         INTEGER II,G1(5,NBASIS),NLMAX
         LOGICAL TCOUNT
         INTEGER ISYM(5),I,J,NFAC,NI(NEL),KJ(5),NJ(NEL)
         INTEGER KI1,KI2,KI3,KI4
         INTEGER NBLOCKS,NBLOCKSTARTS(NBLOCKS+1),OII
         LOGICAL TSPN,TPARITY,TDONE
         INTEGER LMS,IPARITY(3),LMIN,LMAX,BLOCKSYM(5,NBLOCKS)
         LOGICAL TGENFDET
         INTEGER IFDET,IC,IMAX(5,2),NDETTOT,IDEG
         REAL*8 DETSC,TDETSC,TMAT(*)
         REAL*8 CALCT
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
     &      TPARITY,IPARITY,.TRUE.,TDONE,IMAX,ISYM)
         DO WHILE(.NOT.TDONE)
            CALL SETUPSYM(NBASISMAX,KJ)
            CALL  GENSYMDETS_R(NI,ISYM,NEL,G1,BRR,NBASIS,NMRKS,
     &            II,NLMAX,NJ,KJ,1,1,NBASISMAX)
            IF(II.NE.OII) THEN
C.. we've found an occupied block 
               I=I+1
               CALL GETSYMDEGEN(ISYM,NBASISMAX,IDEG)
               NDETTOT=NDETTOT+(II-OII)*IDEG
               IF(TCOUNT) THEN
                  WRITE(14,"(5I5,Z5,2I10)") I,(ISYM(J),J=1,5),II-OII,
     &               (II-OII)*IDEG
               ELSE
                  NBLOCKSTARTS(I)=OII+1
                  CALL ICOPY(4,ISYM,1,BLOCKSYM(1,I),1)
                  IF(TGENFDET) THEN
                     TDETSC=CALCT(NMRKS(1:NEL,OII+1),NEL,
     &                 G1,TMAT,NBASIS)
                     IF(TDETSC.LT.DETSC) THEN
                        IFDET=OII+1
                        DETSC=TDETSC
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            OII=II
            CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &         TPARITY,IPARITY,.FALSE.,TDONE,IMAX,ISYM)
         ENDDO
         NBLOCKS=I
         IF(.NOT.TCOUNT) NBLOCKSTARTS(I+1)=II+1
         IF(TCOUNT) CLOSE(14)
         IF(.NOT.TGENFDET) IFDET=1
         RETURN
      END         

 

      RECURSIVE SUBROUTINE GENSYMDETSSD_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                     NLIST,NLMAX,NJ,KJ,NELEC,NBF,NBASISMAX)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, G1(5,NBASIS),NLMAX,NLIST
         INTEGER LSTE(NEL,NLMAX),KI(5),KJ(5),NJ(NEL),NELEC,NBF
         INTEGER I,J,KJ2(5),BRR(NBASIS),NN(NEL),NBASISMAX(5,2)
         LOGICAL LCMP,LCHKSYM,TTILT
         INTEGER IGETEXCITLEVEL,ICE
         DO I=NBF,NBASIS
            NJ(NELEC)=BRR(I)
            DO J=1,5
               KJ2(J)=KJ(J)
            ENDDO
C.. Check if we've filled all the electrons
            IF(NELEC.EQ.NEL) THEN
               CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
               CALL ROUNDSYM(KJ2,NBASISMAX)
               IF(LCHKSYM(KJ2,KI)) THEN
                CALL ICOPY(NEL,NJ,1,NN,1)
                CALL SORTI(NEL,NN)
                ICE=IGETEXCITLEVEL(NI,NN,NEL)
                IF(ICE.GT.0.AND.ICE.LE.2) THEN
C.. we've found a det with the right sym.
                  NLIST=NLIST+1
                  IF(NLIST.LE.NLMAX) THEN
C.. if there's space, we save it
                     CALL ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
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
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,G1(5,NBASIS),BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST)
         INTEGER KI(5),KJ(5),NJ(NEL),NBASISMAX(5,2)
         INTEGER I,J,NLMAX
         DO J=1,5
            KJ(J)=0
         ENDDO
         DO J=1,NEL
            NJ(J)=0
         ENDDO
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
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,G1(5,NBASIS),BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST),NBASISMAX(5,2)
         INTEGER KI(5),KJ(5),NJ(NEL)
         INTEGER I,J,NLMAX
         DO J=1,5
            KJ(J)=0
         ENDDO
         CALL GETSYM(NI,NEL,G1,NBASISMAX,KI)
         NLMAX=NLIST
         NLIST=0
         CALL GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX)
         RETURN
      END

      RECURSIVE SUBROUTINE GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                     NLIST,NLMAX,NJ,KJ,NELEC,NBF,NBASISMAX)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, G1(5,NBASIS),NLMAX,NLIST
         INTEGER LSTE(NEL,NLMAX),KI(5),KJ(5),NJ(NEL),NELEC,NBF
         INTEGER I,J,KJ2(5),BRR(NBASIS),NN(NEL),NBASISMAX(5,2)
         LOGICAL LCMP,LCHKSYM,TTILT
C         INTEGER KJ3(4)
         DO I=NBF,NBASIS
            NJ(NELEC)=BRR(I)
            DO J=1,5
               KJ2(J)=KJ(J)
C               KJ3(J)=KJ(J)
            ENDDO
C.. Check if we've filled all the electrons
            IF(NELEC.EQ.NEL) THEN
               CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
C               CALL WRITEDET(16,NJ,NEL,.FALSE.)
C                WRITE(16,*) KJ2(1),KJ2(2),KJ2(3)
C               WRITE(6,*) KJ2
C               WRITE(6,*) KI
C               WRITE(6,*)
               IF(LCHKSYM(KJ2,KI)) THEN
C                  CALL ROUNDSYM(KJ3,NBASISMAX)
C.. we've found a det with the right sym.
                  NLIST=NLIST+1
                  CALL ICOPY(NEL,NJ,1,NN,1)
                  CALL SORTI(NEL,NN)
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
                     CALL ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
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
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,G1(5,NBASIS),BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST)
         INTEGER KI(5),KJ(5),NJ(NEL),NBASISMAX(5,2)
         INTEGER I,J,NLMAX
         DO J=1,5
            KJ(J)=0
         ENDDO
         DO J=1,NEL
            NI(J)=0
         ENDDO
         NLMAX=NLIST
         NLIST=0
         CALL GENSYMDETS_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX)
         RETURN
      END


      RECURSIVE SUBROUTINE GENSYMDETSSDN_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &                NLIST,NLMAX,NJ,KJ,NELEC,NBF,NBASISMAX,BETA,I_P,
     &       NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,PP,RHOEPS)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, G1(5,NBASIS),NLMAX,NLIST
         INTEGER LSTE(NEL,NLMAX),KI(5),KJ(5),NJ(NEL),NELEC,NBF
         INTEGER I,J,KJ2(5),BRR(NBASIS),NN(NEL),NBASISMAX(5,2)
         LOGICAL LCMP,LCHKSYM,TTILT
         INTEGER IGETEXCITLEVEL,ICE
         REAL*8 RN,RP,RH,PP,RH2,RHOEPS
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         INTEGER I_P,NMSH,NMAX,NTAY
         DO I=NBF,NBASIS
C.. We continue if we're not searching for a specific det, or if we've
C.. yet to find a specific det.
            IF(PP.EQ.0.D0.OR.NLIST.LT.NLMAX) THEN
               NJ(NELEC)=BRR(I)
               DO J=1,5
                  KJ2(J)=KJ(J)
               ENDDO
C.. Check if we've filled all the electrons
               IF(NELEC.EQ.NEL) THEN
                  CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
                  IF(LCHKSYM(KJ2,KI)) THEN
                   CALL ICOPY(NEL,NJ,1,NN,1)
                   CALL SORTI(NEL,NN)
                   ICE=IGETEXCITLEVEL(NI,NN,NEL)
C.. Check the connectivity
                   IF(ICE.GT.0.AND.ICE.LE.2) THEN
                     NLIST=NLIST+1
                     IF(RP.NE.0.D0) THEN
                        CALL CALCRHO2(NI,NN,BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH2,NTAY,ICE,ECORE)
                        IF(RP.GE.0.D0) THEN
                           RN=RN+RH2**RP
                        ELSEIF(ABS(RH2).GT.RHOEPS) THEN
                           CALL CALCRHO2(NN,NN,BETA,I_P,NEL,
     &             NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,
     &             RH,NTAY,0,ECORE)
                           RN=RN+RH**ABS(RP)
                        ENDIF
                     ENDIF
                     IF(NLIST.LE.NLMAX) THEN
C.. if there's space, we save it
                        IF(PP.NE.0.D0) THEN
                           IF(RN.GE.PP) THEN
                              CALL ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
                           ELSE
                              NLIST=NLIST-1
                           ENDIF
                        ELSE
                           CALL ICOPY(NEL,NN,1,LSTE(1,NLIST),1)
                        ENDIF
                     ENDIF
                   ENDIF
                  ENDIF
               ELSE
C.. otherwise we need to add more electrons:
                  CALL GENSYMDETSSDN_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &           NLIST,NLMAX,NJ,KJ2,NELEC+1,I+1,NBASISMAX,BETA,I_P,
     &         NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,PP,RHOEPS)
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
C.. RP!=0.D0).  
C.. 
C.. If PP>0.D0, then as RN is summed, if a det's contribution to RN
C.. exceeds PP, then the procedure ends.
C..
C.. NLIST on input is set to the number of dets which can fit into 
C.. LSTE
C.. On exit, it contains the number of dets generated.
C.. If there are more dets generated than the input NLIST, then only
C.. the first NLIST are stored.

      SUBROUTINE GENSYMDETSSDN(NI,KI,NEL,G1,BRR,NBASIS,LSTE,
     &         NLIST,NBASISMAX,
     &         BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,
     &         PP,RHOEPS)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,G1(5,NBASIS),BRR(NBASIS)
         INTEGER NLIST,LSTE(NEL,NLIST)
         INTEGER KI(5),KJ(5),NJ(NEL),NBASISMAX(5,2)
         INTEGER I,J,NLMAX
         REAL*8 RN,RP,PP,RHOEPS
         REAL*8 BETA,FCK(*),ZIA(*),ALAT(*),UMAT(*),ECORE
         INTEGER I_P,NMSH,NMAX,NTAY
         DO J=1,5
            KJ(J)=0
         ENDDO
         DO J=1,NEL
            NJ(J)=0
         ENDDO
         NLMAX=NLIST
         NLIST=0
         RN=0.D0
         CALL GENSYMDETSSDN_R(NI,KI,NEL,G1,BRR,NBASIS,LSTE,NLIST,
     &                        NLMAX,NJ,KJ,1,1,NBASISMAX,
     &         BETA,I_P,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,NTAY,ECORE,RN,RP,
     &         PP,RHOEPS)
         RETURN
      END

C.  Irrep symmetries are specified in SYM(5).
C.. if SYM(5)=0, we assume it's totally symmetric
C.. Other irreps contributing to the symmetry have bits set in 
C.. SYM.
C.. e.g. if irreps are a1,a2,b1,b2
      LOGICAL FUNCTION LCHKSYM(ISYM,JSYM)
         IMPLICIT NONE
         INTEGER I,ISYM(5),JSYM(5),IS,JS
         INTEGER SYMPROD
         LOGICAL LSYMSYM
         LCHKSYM=.TRUE.
         DO I=1,4
            IF(ISYM(I).NE.JSYM(I)) LCHKSYM=.FALSE.
         ENDDO
C.. if the symmetry product of I and J doesn't contain the totally
C.. symmetric irrep, we set sym to .FALSE.
         LCHKSYM=LCHKSYM.AND.LSYMSYM(SYMPROD(ISYM(5),JSYM(5)))
      RETURN
      END
      
      LOGICAL FUNCTION LCHKSYMD(NI,NJ,NEL,G1,NBASISMAX)
         IMPLICIT NONE
         INTEGER ISYM(5),JSYM(5)
         INTEGER NEL,NI(NEL),NJ(NEL),G1(5,*),NBASISMAX(5,3)
         LOGICAL LCHKSYM
         CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
         CALL GETSYM(NJ,NEL,G1,NBASISMAX,JSYM)
         LCHKSYMD=LCHKSYM(ISYM,JSYM)
         RETURN
      END
C.. NBASISMAX descriptor (1,3)
C
C HUBBARD:
C  BITS
C. 0 Tilted
C. 1 non-pbc
C. 2 real-space
C.. which effects to values
C.. MOM SPACE
C 0 Non-Tilted Lattice - pbc
C 1 Tilted Lattice - pbc
C 2 Non-Tilted lattice - no pbc
C 3 Tilted Lattice - no pbc
C.. four following are REAL
C 4 Non-Tilted Lattice - pbc
C 5 Tilted Lattice - pbc
C 6 Non-Tilted lattice - no pbc
C 7 Tilted Lattice - no pbc
C
C (3,3)
C -2 Particle in a box
C -1 UEG
C 0 Hubbard
C 1 Generic spatial


C.. This only works for momentum variables - 1-4
      SUBROUTINE ADDELECSYM(IEL,G1,NBASISMAX,ISYM)
         IMPLICIT NONE
         INTEGER IELEC,G1(5,*),NBASISMAX(5,3)
         INTEGER I,IEL,SSYM
         INTEGER ISYM(5)
         CALL GETUNCSFELEC(IEL,IELEC,SSYM)
        IF(NBASISMAX(1,3).LT.4) THEN
C.. Momentum space
            DO I=1,4
               ISYM(I)=ISYM(I)+G1(I,IELEC)
            ENDDO
C.. Symmetry space
         ELSEIF(NBASISMAX(3,3).EQ.0.AND.NBASISMAX(1,3).GE.4) THEN
C.. We have no symmetries, so do nothing. (we're in real space)
C.. except Ms
            ISYM(4)=ISYM(4)+G1(4,IELEC)
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
C.. deal with momentum
            DO I=1,4
               ISYM(I)=ISYM(I)+G1(I,IELEC)
            ENDDO
         ENDIF
C.. SSYM keeps track of the total S change on adding this electron
C.. (it is +/-CSF_NSBASIS)
         ISYM(4)=ISYM(4)+SSYM
         RETURN
      END 
      
      SUBROUTINE ROUNDSYM(ISYM,NBASISMAX)
         IMPLICIT NONE
         INTEGER ISYM(5),NBASISMAX(5,3)
         INTEGER I
         IF(NBASISMAX(3,3).EQ.-2) THEN
C.. particle in a box
C.. parity symmetries
            DO I=1,3
               ISYM(I)=MOD(ISYM(I),2)
            ENDDO
         ELSEIF(NBASISMAX(3,3).EQ.-1) THEN
C.. UEG (can't remember the symmetries of that
C.. probably momentum  conservation)
         ELSEIF(NBASISMAX(3,3).EQ.0) THEN
C.. Hubbard model
            IF(NBASISMAX(1,3).LT.2) THEN
C.. momentum conservation - various PBC issues
               CALL MOMPBCSYM(ISYM,NBASISMAX)
C            ELSEIF(NBASISMAX(1,3).EQ.2) THEN
C.. non-pbc mom space has parity symmetry
C               DO I=1,3
C                  ISYM(I)=MOD(ISYM(I),2)
C               ENDDO
            ELSEIF(NBASISMAX(1,3).GE.2) THEN
C.. we're in real space so no sym
               DO I=1,3
                  ISYM(I)=0
               ENDDO
            ENDIF
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
C.. Generic spatial symmetries
C..         We need do nothing.
C.. However, there is still momentum conservation - various PBC issues
               CALL MOMPBCSYM(ISYM,NBASISMAX)
         ENDIF
         RETURN 
      END

C.. NBASISMAX descriptor (1,3)
C
C HUBBARD:
C  BITS
C. 0 Tilted
C. 1 non-pbc
C. 2 real-space
C.. which effects to values
C.. MOM SPACE
C 0 Non-Tilted Lattice - pbc
C 1 Tilted Lattice - pbc
C 2 Non-Tilted lattice - no pbc
C 3 Tilted Lattice - no pbc
C.. four following are REAL
C 4 Non-Tilted Lattice - pbc
C 5 Tilted Lattice - pbc
C 6 Non-Tilted lattice - no pbc
C 7 Tilted Lattice - no pbc
C
C (3,3)
C -2 Particle in a box
C -1 UEG
C 0 Hubbard
C 1 Generic spatial


      SUBROUTINE MOMPBCSYM(K1,NBASISMAX)
C.. NB the third column of NBASISMAX tells us whether it is tilted
         IMPLICIT NONE
         INTEGER K1(3),NBASISMAX(5,5)
         INTEGER J,LDIM,AX,AY,LENX,LENY,KK2,T1,T2
         REAL*8 R1,R2,NORM
         AX=NBASISMAX(1,4)
         AY=NBASISMAX(2,4)
         LENX=NBASISMAX(1,5)
         LENY=NBASISMAX(2,5)
         IF(NBASISMAX(1,3).EQ.0.OR.NBASISMAX(1,3).EQ.0) THEN
C.. A non-tilted lattice with PBC
            DO J=1,3
C..  non-tilted
               KK2=K1(J)
               LDIM=NBASISMAX(J,2)-NBASISMAX(J,1)+1
               KK2=MOD(KK2,LDIM)
               IF(KK2.LT.NBASISMAX(J,1)) KK2=KK2+LDIM
               IF(KK2.GT.NBASISMAX(J,2)) KK2=KK2-LDIM
               K1(J)=KK2 
            ENDDO
         ELSEIF(NBASISMAX(1,3).EQ.1) THEN
C.. we have a tilted lattice with PBC
C.. we want the a1,a2 components of k
            NORM=AX*AX+AY*AY
            R1=(AX*K1(1)+AY*K1(2))/NORM
            R2=(AX*K1(2)-AY*K1(1))/NORM
            R1=R1/LENX+0.5D0
            R2=R2/LENY+0.5D0
C.. T1= highest integer less than R1
            T1=INT(ABS(R1))
            IF(R1.LT.0.D0) THEN
               T1=-T1
               IF(T1.NE.R1) T1=T1-1
            ENDIF
            T2=INT(ABS(R2))
            IF(R2.LT.0.D0) THEN
               T2=-T2
               IF(T2.NE.R2) T2=T2-1
            ENDIF
            IF(R1.EQ.T1) T1=T1-1
            IF(R1.GT.1.D0.OR.R1.LE.0.D0) R1=R1-T1
            IF(R2.GE.1.D0.OR.R2.LT.0.D0) R2=R2-T2
            R1=(R1-0.5D0)*LENX
            R2=(R2-0.5D0)*LENY
            K1(1)=NINT(R1*AX-R2*AY)
            K1(2)=NINT(R1*AY+R2*AX)
         ENDIF
         RETURN
      END

C.. Get the next symmetric determinant in sequence.
C.. TDONE will be set to TRUE if there are no more dets.
C.. KI is the sym, NJ contains the working
      SUBROUTINE GENNEXTSYMDET(KI,NEL,G1,BRR,NBASIS,
     &            TDONE,NJ,NN,INDJ,NELEC,NBASISMAX)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL), NBASIS, G1(5,NBASIS),KJ2(5)
         INTEGER KI(5),NJ(NEL),NELEC,INDJ(NEL),NLAST
         INTEGER I,J,BRR(NBASIS),NN(NEL),NBASISMAX(5,2)
         LOGICAL LCHKSYM,LCONT,TDONE
         TDONE=.FALSE.
         LCONT=.TRUE.
         DO WHILE (LCONT)
C.. Remove the current electron in this slot
C.. Only do this if there was one there before.
C            IF((NELEC.GT.1.AND.INDJ(NELEC).GT.(INDJ(NELEC-1)))
C     &               .OR.(NELEC.EQ.1.AND.INDJ(1).GT.0)) THEN
C.. we don't even have to set it to zero, as it'll be overwritten
C            ENDIF
            NLAST=NBASIS-(NEL-NELEC)
C.. if we've past the last allowable electron in our position
            IF(INDJ(NELEC).GE.NLAST) THEN
C.. Move back a slot
               NELEC=NELEC-1
C.. If we move back past the beginning
               IF(NELEC.EQ.0) THEN
C.. we're done
                  TDONE=.TRUE.
                  LCONT=.FALSE.
               ENDIF
            ELSE
C..Increment this slot
               I=INDJ(NELEC)+1
               INDJ(NELEC)=I
               IF(I.LE.NBASIS) NJ(NELEC)=BRR(I)
C.. Check if we've filled all the electrons
               IF(NELEC.EQ.NEL) THEN
                  CALL GETSYM(NJ,NEL,G1,NBASISMAX,KJ2)
                  IF(I.LE.NBASIS.AND.LCHKSYM(KI,KJ2)) THEN
                     CALL ICOPY(NEL,NJ,1,NN,1)
                     CALL SORTI(NEL,NN)
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
     &      TPARITY,IPARITY,ISYM,NI,TSETUP,TDONE,WORK)
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),NBASISMAX(5,3),G1(5,NBASIS)
         INTEGER LMS,IPARITY(3),ISYM(5),NIWORK(NEL),NI(NEL)
         INTEGER WORK(5*NEL+18)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE
C.. WORK(7*NEL+18)
C..  1-4                NSWORK
C..  5-14               IMAX
C..  15-19              ISYM
C..  20                 NELEC
C..  21-NEL+20          NIWORK
C..  NEL+21 - 2*NEL+20  NINDJ
         CALL GENNEXTDET_(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &      TPARITY,IPARITY,NI,TSETUP,TDONE,WORK(5),
     &      WORK(15),WORK(20),WORK(21),WORK(NEL+21))
         CALL ICOPY(5,WORK(15),1,ISYM,1)
         RETURN
      END
      SUBROUTINE GENNEXTDET_(NEL,NBASIS,BRR,NBASISMAX,G1,TSPN,LMS,
     &      TPARITY,IPARITY,NI,TSETUP,TDONE,IMAX,ISYM,NELEC,
     &      NIWORK,NINDJ)
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),NBASISMAX(5,3),G1(5,NBASIS)
         INTEGER LMS,IPARITY(3),ISYM(5),NIWORK(NEL),NI(NEL)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE,KALLOWED,TMORE2
         INTEGER ILEV,IMAX(5,2)
C..  1-4                NSWORK
C..  5-12               IMAX
C..  13-17              ISYM
C..  18                 NELEC
C..  19-NEL+18          NIWORK
C..  NEL+19 - 2*NEL+18  NINDJ
         INTEGER NELEC,NINDJ(NEL)
         IF(TSETUP) THEN
            CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &         TPARITY,IPARITY,TSETUP,TDONE,IMAX,ISYM)
            NELEC=1
            NINDJ(1)=0
            RETURN
         ENDIF
C.. Go to the next det.
         TMORE=.TRUE.
         DO WHILE(TMORE)
            CALL GENNEXTSYMDET(ISYM,NEL,G1,BRR,NBASIS,
     &           TDONE,NIWORK,NI,NINDJ,NELEC,NBASISMAX)
            IF(TDONE) THEN
               CALL GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &            TPARITY,IPARITY,TSETUP,TMORE,IMAX,ISYM)
               TMORE=.NOT.TMORE
               NELEC=1
               NINDJ(1)=0
            ELSE
               TMORE=.FALSE.
            ENDIF
         ENDDO
      END
      SUBROUTINE GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &      TPARITY,IPARITY,TSETUP,TDONE,IMAX,ISYM)
         IMPLICIT NONE
         INTEGER NEL,NBASISMAX(5,3)
         INTEGER LMS,IPARITY(5),ISYM(5)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE,KALLOWED,TMORE2
         INTEGER ILEV,IMAX(5,2)
         INCLUDE 'sym.inc'
         IF(TSETUP) THEN
            DO ILEV=1,3
               IF(TPARITY) THEN
                  IMAX(ILEV,1)=IPARITY(ILEV)
                  IMAX(ILEV,2)=IPARITY(ILEV)
               ELSE
                  IMAX(ILEV,1)=NBASISMAX(ILEV,1)
                  IMAX(ILEV,2)=NBASISMAX(ILEV,2)
C                  IF(NBASISMAX(1,3).EQ.2) THEN
C.. hubbard non-pbc mom space
C                     IMAX(ILEV,1)=IMAX(ILEV,1)*NEL
C                     IMAX(ILEV,2)=IMAX(ILEV,2)*NEL
C                  ENDIF
               ENDIF
            ENDDO
            IF(TSPN) THEN
               IMAX(4,1)=LMS
               IMAX(4,2)=LMS
            ELSE
               IMAX(4,1)=NBASISMAX(4,1)*NEL
               IMAX(4,2)=NBASISMAX(4,2)*NEL
            ENDIF
C.. If we're specifying a sym (TPARITY) in IPARITY(5), and
C.. we have a system with all 1D reducible orbs, then we put
C.. that into IMAX
            IF(NBASISMAX(5,2).NE.0) THEN
               IF(TPARITY) THEN
                  IMAX(5,1)=2**IPARITY(5)
                  IMAX(5,2)=IMAX(5,1)
               ELSE
                  IMAX(5,1)=1
                  IMAX(5,2)=2**NBASISMAX(5,2)
               ENDIF
            ELSE
C.. we've got a sym system with polydimensional irreps, which leads to
C.. dets with combinations of irreps, so we cannot put sym into blocks
               IMAX(5,1)=0
               IMAX(5,2)=0
            ENDIF
            TDONE=.FALSE.
            CALL DOSYMLIMDEGEN(IMAX,NBASISMAX)
            DO ILEV=1,5
               ISYM(ILEV)=IMAX(ILEV,1)
            ENDDO
         ENDIF
         IF(TSETUP.AND.KALLOWED(ISYM,NBASISMAX)) RETURN
C.. Go to the next sym.
         TMORE2=.TRUE.
         TMORE=.TRUE.
         ILEV=5
         DO WHILE(TMORE2)
            DO WHILE (ILEV.GT.0)
               IF(ILEV.EQ.5.AND.ISYM(ILEV).NE.0) THEN
C.. symmetry specifiers are incremented by multiplying*2 (unless there are no syms counted)
                  ISYM(ILEV)=ISYM(ILEV)*2
               ELSE
                  ISYM(ILEV)=ISYM(ILEV)+1
               ENDIF
               IF(ISYM(ILEV).GT.IMAX(ILEV,2)) THEN
                  ISYM(ILEV)=IMAX(ILEV,1)
                  ILEV=ILEV-1
                  IF(ILEV.EQ.0) THEN
                     TMORE2=.FALSE.
C.. If we've run out of syms, we give up
                     TMORE=.FALSE.
                  ENDIF
               ELSEIF(ILEV.LT.4) THEN
C.. We've just incremented one of the higher columns, now go down to the
C.. lower ones.
                  ILEV=ILEV+1
                  ISYM(ILEV)=IMAX(ILEV,1)-1
                  
               ELSEIF(KALLOWED(ISYM,NBASISMAX)) THEN
                  TMORE2=.FALSE.
                  ILEV=0
               ENDIF
            ENDDO
         ENDDO
         TDONE=.NOT.TMORE
      END
      SUBROUTINE DOSYMLIMDEGEN(IMAX,NBASISMAX)
         IMPLICIT NONE
         INTEGER IMAX(5,2),NBASISMAX(5,5),I
         IF(NBASISMAX(3,3).EQ.0) THEN
            DO I=1,3
               IF(IMAX(I,2).NE.IMAX(I,1)) IMAX(I,1)=0
            ENDDO
         ENDIF
C.. always a spin symmetry
         IF(IMAX(4,1).NE.IMAX(4,2)) IMAX(4,1)=0
      END
      SUBROUTINE GETSYMDEGEN(ISYM,NBASISMAX,IDEGEN)
         IMPLICIT NONE
         INTEGER ISYM(5),NBASISMAX(5,5),IDEGEN,I,ISYM2(5),J
         LOGICAL KALLOWED,TDO
         IDEGEN=0
         IF(NBASISMAX(3,3).EQ.0) THEN
C.. Hubbard
            DO I=0,7
               TDO=.TRUE.
               DO J=1,3
                  IF(BTEST(I,J-1).EQ.0) THEN
                     ISYM2(J)=ISYM(J)
                  ELSE
                     ISYM2(J)=-ISYM(J)
                     IF(ISYM(J).EQ.0) TDO=.FALSE.
                  ENDIF
               ENDDO
               IF(TDO.AND.KALLOWED(ISYM2,NBASISMAX)) IDEGEN=IDEGEN+1
            ENDDO
         ELSE
            IDEGEN=1
         ENDIF
C.. Spin
         IF(ISYM(4).NE.0) IDEGEN=IDEGEN*2
      END

C.. Initialize symmetry to take into account the core electrons
      SUBROUTINE SETUPSYM(NBASISMAX,ISYM)
         IMPLICIT NONE
         INTEGER NBASISMAX(5,6),ISYM(5)
         CALL ICOPY(5,NBASISMAX(1,6),1,ISYM,1)
         RETURN
      END


