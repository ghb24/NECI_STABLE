      LOGICAL FUNCTION AREDETSEXCITS(D1,D2,NEL,NBASISMAX,ICE)
         IMPLICIT NONE
         INTEGER ICE,NEL,D1(NEL),D2(NEL),NBASISMAX(*)
         LOGICAL ISUHFDET
         INTEGER IGETEXCITLEVEL
         ICE=IGETEXCITLEVEL(D1,D2,NEL)
         IF(ICE.LE.2) THEN
C.. The new node is connected to this node
            IF(ISUHFDET(D1,NEL,NBASISMAX).OR.
     &         ISUHFDET(D2,NEL,NBASISMAX)) THEN
               IF(ICE.NE.1) THEN
                  AREDETSEXCITS=.TRUE.
               ELSE
                  AREDETSEXCITS=.FALSE.
               ENDIF
            ELSE
               AREDETSEXCITS=.TRUE.
            ENDIF
         ELSE
            AREDETSEXCITS=.FALSE.
         ENDIF
         RETURN
      END

      LOGICAL FUNCTION ISVALIDDET(DET,NEL)
         IMPLICIT NONE
         INTEGER NEL,DET(NEL)
         INTEGER I
         IF(DET(1).LT.1) THEN
            ISVALIDDET=.FALSE.
            RETURN
         ENDIF
         DO I=2,NEL
            IF(DET(I-1).GE.DET(I)) THEN
               ISVALIDDET=.FALSE.
               RETURN
            ENDIF
         ENDDO
         ISVALIDDET=.TRUE.
         RETURN
      END

      SUBROUTINE EXCIT_DUMP(UNT,NI,NEL,G1,NBASIS,NBASISMAX,
     &         NMEM,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER NMEM(*),UNT
         INTEGER IFRZ(0:NBASIS)
C.. now call a routine to reset the current type's lists
         CALL GSETE_DUMP(UNT,NI,NEL,G1,NBASIS,NBASISMAX,
     &        NMEM(1),NMEM(2),NMEM(3),NMEM(4),NMEM(5),
     &         NMEM(7),NMEM(11),
     &         NMEM(11+(NEL+4)*4),
     &         NMEM(11+(NEL+NBASIS+8)*4),IFRZ)
         RETURN
      END
      SUBROUTINE GSETE_DUMP(UNT,NI,NEL,G1,NBASIS,NBASISMAX,
     &         NTYPES,NTYPE,TFULL,TFREE,TTYPE,SPNS,
     &         NFULL,NFREE,
     &         TYPES,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,NTYPES),NTYPES,NTYPE
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         INTEGER SPNS(2,2),IC
         LOGICAL TFULL,TFREE,TTYPE
         INTEGER I,J
         INTEGER IFRZ(0:NBASIS),UNT
         
         IF(NTYPE.EQ.0) RETURN
C..Debugging junk
         CALL WRITEDET(UNT,NI,NEL,.FALSE.)
         WRITE(UNT,"(5Z5,I)") (TYPES(I,NTYPE),I=1,5)
        WRITE(UNT,"(9Z5,I)") (TYPES(I,NTYPE),I=1,5),SPNS(1,1),SPNS(2,1),
     &      SPNS(1,2),SPNS(2,2),NTYPE
         DO J=1,2
         WRITE(UNT,"(A,$)") "FULL"
         DO I=-3,NFULL(J,1,0)
            WRITE(UNT,"(I3,$)") NFULL(J,1,I)
         ENDDO
         WRITE(UNT,"(A,$)") ':'
         DO I=-3,NFULL(J,2,0)
            WRITE(UNT,"(I3,$)") NFULL(J,2,I)
         ENDDO
         WRITE(UNT,*)
         WRITE(UNT,"(A,$)") "FREE"
         DO I=-3,NFREE(J,1,0)
            WRITE(UNT,"(I3,$)") NFREE(J,1,I)
         ENDDO
         WRITE(UNT,"(A,$)") ':'
         DO I=-3,NFREE(J,2,0)
            WRITE(UNT,"(I3,$)") NFREE(J,2,I)
         ENDDO
         WRITE(UNT,*)
         ENDDO
         RETURN
      END
      SUBROUTINE RESETEXIT(NI,NEL,G1,NBASIS,NBASISMAX,
     &         NMEM,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6),NJ(NEL)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER NMEM(*)
         INTEGER IFRZ(0:NBASIS)
C.. now call a routine to reset the current type's lists
         CALL GSETE_RSTTYP(NI,NEL,G1,NBASIS,NBASISMAX,
     &        NMEM(1),NMEM(2),NMEM(3),NMEM(4),NMEM(5),
     &         NMEM(7),NMEM(11),
     &         NMEM(11+(NEL+4)*4),
     &         NMEM(11+(NEL+NBASIS+8)*4),IFRZ)
         RETURN
      END
      SUBROUTINE GSETE_RSTTYP(NI,NEL,G1,NBASIS,NBASISMAX,
     &         NTYPES,NTYPE,TFULL,TFREE,TTYPE,SPNS,
     &         NFULL,NFREE,
     &         TYPES,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,NTYPES),NTYPES,NTYPE
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         INTEGER SPNS(2,2),IC
         LOGICAL TFULL,TFREE,TTYPE
         INTEGER IEL,ISPN,I,J
         LOGICAL TDONE
         INTEGER IFRZ(0:NBASIS)
C.. NXXXX(a,b,-3) is  is the spin (either 1 or 2) for the first electron to be 
C.. excited (so NFREE(NFREE(1,-3),-2) is the 1st electron to be excited)
C.. [???(8/4/05)]

C.. NXXXX(a,b,-2) corresponds to the first basis fn to be excited, and -1 to the
C.. second.

C.. NXXXX(a,b,0) is the number of basisfns available in this group
C.. NXXXX(a,b,1:...) are the basisfns in this group.  for FREE, it is
C..                  the actual basis fn, and FULL it is the index
C..                  within NI.

C.. a corresponds to the sym to be excited - 1 for the 1st and 2 for 2nd
C.. b is the spin of the appropriate symmetry - beta=1, alpha=2

C.. Go through all lists for this type and remove the any frozen
C.. orbitals.  We do not add orbitals if they have become unfrozen
         
         DO IEL=1,2
            DO ISPN=1,2
C.. full
               I=1
               DO WHILE(I.LE.NFULL(IEL,ISPN,0))
                  IF(IFRZ(NI(NFULL(IEL,ISPN,I))).NE.0) THEN
C.. if we've been frozen we remove the orbital
                     NFULL(IEL,ISPN,0)=NFULL(IEL,ISPN,0)-1
                     DO J=I,NFULL(IEL,ISPN,0)
                        NFULL(IEL,ISPN,J)=NFULL(IEL,ISPN,J+1)
                     ENDDO
                     IF(NFULL(IEL,ISPN,-2).GT.I)
     &                  NFULL(IEL,ISPN,-2)=NFULL(IEL,ISPN,-2)-1
                     IF(NFULL(IEL,ISPN,-1).GT.I)
     &                  NFULL(IEL,ISPN,-1)=NFULL(IEL,ISPN,-1)-1
                  ENDIF
                  I=I+1
               ENDDO
C.. If we've now run out of excitations flag this.
               IF(  NFULL(IEL,ISPN,-2).GT.NFULL(IEL,ISPN,0)
     &          .OR.NFULL(IEL,ISPN,-1).GT.NFULL(IEL,ISPN,0)) 
     &            TFULL=.FALSE.
C.. Now the empty ones
               IF(.NOT.TFULL) TTYPE=.FALSE.
               I=1
               DO WHILE(I.LE.NFREE(IEL,ISPN,0))
                  IF(IFRZ(NFREE(IEL,ISPN,I)).NE.0) THEN
C.. if we've been frozen we remove the orbital
                     NFREE(IEL,ISPN,0)=NFREE(IEL,ISPN,0)-1
                     DO J=I,NFREE(IEL,ISPN,0)
                        NFREE(IEL,ISPN,J)=NFREE(IEL,ISPN,J+1)
                     ENDDO
                     IF(NFREE(IEL,ISPN,-2).GT.I)
     &                  NFREE(IEL,ISPN,-2)=NFREE(IEL,ISPN,-2)-1
                     IF(NFREE(IEL,ISPN,-1).GT.I)
     &                  NFREE(IEL,ISPN,-1)=NFREE(IEL,ISPN,-1)-1
                  ENDIF
                  I=I+1
               ENDDO
               IF(  NFREE(IEL,ISPN,-2).GT.NFREE(IEL,ISPN,0)
     &          .OR.NFREE(IEL,ISPN,-1).GT.NFREE(IEL,ISPN,0)) 
     &            TFREE=.FALSE.
            ENDDO
         ENDDO
         RETURN
      END
      
C.. Get the orbitals which are excited in going from I to J
C.. IEX(1,*) are in I, and IEX(2,*) are in J
      SUBROUTINE GETEXCITATION(NI,NJ,NEL,EX)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),EX(2,2)
         INTEGER I,J
         INTEGER IC1,IC2
         CALL IAZZERO(EX,2*2)
         IC1=0
         IC2=0
         I=1
         J=1
C.. We only count differences from I to J
         DO WHILE(I.LE.NEL.AND.J.LE.NEL)
            DO WHILE(NI(I).LT.NJ(J).AND.I.LE.NEL)
               IC1=IC1+1
               IF(IC1.LE.2) EX(1,IC1)=NI(I)
               I=I+1
            ENDDO
            DO WHILE(NI(I).GT.NJ(J).AND.I.LE.NEL.AND.J.LE.NEL)
               IC2=IC2+1
               IF(IC2.LE.2) EX(2,IC2)=NJ(J)
               J=J+1
            ENDDO
            IF(I.LE.NEL.AND.J.LE.NEL.AND.NI(I).EQ.NJ(J)) THEN
               I=I+1
               J=J+1
            ENDIF
         ENDDO
         DO WHILE(I.LE.NEL)
            IC1=IC1+1
            IF(IC1.LE.2) EX(1,IC1)=NI(I)
            I=I+1
         ENDDO
         DO WHILE(J.LE.NEL)
            IC2=IC2+1
            IF(IC2.LE.2) EX(2,IC2)=NJ(J)
            J=J+1
         ENDDO
         RETURN
      END



C.. Generate excitations of up to 2 electrons within the symmetry
C.. (including Ms) group of NI (but not using momentum sym)
C..
C.. If TSYM is set, then each excitation is checked against the sym
C.. ISYM, and discarded if it has the wrong symmetry
      SUBROUTINE GENRANDSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,ISEED,IC,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6),NJ(NEL)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER ISYM(5)
         LOGICAL TSYM
         LOGICAL TSETUP
         INTEGER NMEM(*)
         INTEGER ISEED
         INTEGER IFRZ(0:NBASIS)
         STOP 'Excitation generators not functional'
         IF(TSETUP) THEN
C.. we need to be called with TSETUP=.T. first, to work out how much
C.. memory needs to be allocated for temp storage.  We return this in 
C.. NMEM
            NTYPES=0
            CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &         IFRZ)
C.. NTYPES
            NMEM(1)=1
C.. NTYPE
            NMEM(1)=NMEM(1)+1
C.. TFULL
            NMEM(1)=NMEM(1)+1
C.. TFREE
            NMEM(1)=NMEM(1)+1
C.. TTYPE
            NMEM(1)=NMEM(1)+1
C.. ICOUNT
            NMEM(1)=NMEM(1)+1
C.. SPNS(2,2)
            NMEM(1)=NMEM(1)+4
C.. NFULL and NFREE
            NMEM(1)=NMEM(1)+(NEL+4+NBASIS+4)*4
C.. the list of excitation types (each type is 5 long)
            NMEM(1)=NMEM(1)+5*NTYPES
            RETURN
         ENDIF
         IF(NMEM(1).EQ.0) THEN
C.. this is the first time we've been called with allocated memory.
C.. NTYPES=0, so we need to find the excitation types again
C.. The first time counts them again (inelegant, but not too
C.. time-consuming), and the second time saves them
            CALL IAZZERO(NMEM(1),10)
C            WRITE(6,*) "SETUP"
C            WRITE(6,*) NMEM(6)
C            CALL WRITEDET(6,NI,NEL,.TRUE.)
            CALL COUNTSYMEXCIT2(NI,NEL,G1,NBASIS,NBASISMAX,NMEM(6),
     &      NMEM(11+(NEL+NBASIS+8)*4),NMEM(1),IFRZ)
         ENDIF
         I=NMEM(6)
C.. Decide upon which type of excitation to choose
         CALL GSETE_GETRANDEXCITTYPE(I,NMEM(2),
     &       NMEM(11+(NEL+NBASIS+8)*4),NMEM(1),ISEED)
C.. now call a routine to generate a random excitation
         CALL GSETE_GENRAND(NI,NJ,IC,NEL,G1,NBASIS,NBASISMAX,
     &        TSYM,ISYM,NMEM(1),NMEM(2),NMEM(3),NMEM(4),NMEM(5),
     &         NMEM(7),NMEM(11),
     &         NMEM(11+(NEL+4)*4),
     &         NMEM(11+(NEL+NBASIS+8)*4),ISEED,IFRZ)
         IC=I
         RETURN
      END
C.. Return a count of the number of up to double excitations of NI
C.. available within the same symmetry and Sz.  This does not enumerate them all.
      SUBROUTINE COUNTSYMEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,NEX,IFRZ)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,G1(5,NBASIS),NBASIS,NBASISMAX(5,6)
         INTEGER NEX

         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER IFRZ(0:NBASIS)
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         STOP 'Excitation generators not functional'
         NTYPES=0
         CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         CALL MEMORY(IP_TYPES,5*NTYPES,'TYPES')
         CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         IC=0
         DO I=1,NTYPES
           CALL GETSYMEXTYPECOUNT(NI,NEL,G1,NBASIS,NBASISMAX,TYPES,I,
     &         NFULL,NFREE,IFRZ)
           IC=IC+TYPES(5,I)
         ENDDO
         NEX=IC
         CALL FREEM(IP_TYPES)
         RETURN
      END

C.. Generate excitations of up to 2 electrons within the symmetry
C.. (including Ms) group of NI (but not using momentum sym)
C..
C.. If TSYM is set, then each excitation is checked against the sym
C.. ISYM, and discarded if it has the wrong symmetry
      SUBROUTINE GENSYMEXCITIT(NI,NEL,G1,NBASIS,NBASISMAX,TSYM,ISYM,
     &         TSETUP,NMEM,NJ,IC,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6),NJ(NEL)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER ISYM(5)
         LOGICAL TSYM
         LOGICAL TSETUP
         INTEGER NMEM(*)
         INTEGER IFRZ(0:NBASIS)
         STOP 'Excitation generators not functional'
         IF(TSETUP) THEN
C.. we need to be called with TSETUP=.T. first, to work out how much
C.. memory needs to be allocated for temp storage.  We return this in 
C.. NMEM
            NTYPES=0
            CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &         IFRZ)
C.. NTYPES
            NMEM(1)=1
C.. NTYPE
            NMEM(1)=NMEM(1)+1
C.. TFULL
            NMEM(1)=NMEM(1)+1
C.. TFREE
            NMEM(1)=NMEM(1)+1
C.. TTYPE
            NMEM(1)=NMEM(1)+1
C.. ICOUNT
            NMEM(1)=NMEM(1)+1
C.. SPNS(2,2)
            NMEM(1)=NMEM(1)+4
C.. NFULL and NFREE
            NMEM(1)=NMEM(1)+(NEL+4+NBASIS+4)*4
C.. the list of excitation types (each type is 5 long)
            NMEM(1)=NMEM(1)+5*NTYPES
            RETURN
         ENDIF
         IF(NMEM(1).EQ.0) THEN
C.. this is the first time we've been called with allocated memory.
C.. NTYPES=0, so we need to find the excitation types again
C.. The first time counts them again (inelegant, but not too
C.. time-consuming), and the second time saves them
           CALL GETSYMEXTYPES(NBASISMAX,NMEM,NMEM,NI,NEL,G1,NBASIS,IFRZ)
           CALL GETSYMEXTYPES(NBASISMAX,NMEM(11+(NEL+NBASIS+8)*4),NMEM
     &      ,NI,NEL,G1,NBASIS,IFRZ)
           CALL IAZZERO(NMEM(2),9)
         ENDIF
C.. now call a routine to enumerate the excitations
         CALL GSETE_GET(NI,NJ,IC,NEL,G1,NBASIS,NBASISMAX,
     &        TSYM,ISYM,NMEM(1),NMEM(2),NMEM(3),NMEM(4),NMEM(5),
     &         NMEM(7),NMEM(11),
     &         NMEM(11+(NEL+4)*4),
     &         NMEM(11+(NEL+NBASIS+8)*4),IFRZ)
         NMEM(6)=NMEM(6)+1
         RETURN
      END
C.. Recurse through all possible allowed excitations, and put them (in a
C.. non-sorted order, which makes this unsuitable for MCPATHSR3) into
C.. list LIST, and their excitation level into list ICLIST.
      SUBROUTINE GENSYMEXCIT(NI,NEL,G1,NBASIS,NBASISMAX,TSYM,ISYM,
     &            LIST,ICLIST,NLIST,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,LIST(NEL,NLIST),ICLIST(NLIST)
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER ISYM(5)
         LOGICAL TSYM
         INTEGER IFRZ(0:NBASIS)
         STOP 'Excitation generators not functional'
         NTYPES=0
         CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         CALL MEMORY(IP_TYPES,5*NTYPES,'TYPES')
         CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         IC=0
C         CALL WRITEDET(60,NI,NEL,.TRUE.)
         DO I=1,NTYPES
            CALL GETSYMEXTYPEEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSYM,ISYM,
     &         TYPES(1,I),LIST(1,IC+1),ICLIST(IC+1),NLIST-IC,IFRZ)
            IC=IC+TYPES(5,I)
            IF(IC.GT.NLIST) STOP 'LIST FOR EXCITATIONS NOT LONG ENOUGH'
         ENDDO
         NLIST=IC
         CALL FREEM(IP_TYPES)
         RETURN
      END

      SUBROUTINE GSETE_GETRANDEXCITTYPE(IC,NTYPE,TYPES,NTYPES,ISEED)
         IMPLICIT NONE
         INTEGER IC,NTYPE,TYPES(5,NTYPES),NTYPES,ISEED
         INTEGER I
         REAL*8 RAN2
C.. Pick a random excit type such that all excitations have equal prob
         I=RAN2(ISEED)*IC
         NTYPE=0
C         WRITE(6,"(2I,$)"),I,IC
         DO WHILE(I.GE.0)
            NTYPE=NTYPE+1
            I=I-TYPES(5,NTYPE)
         ENDDO
C         WRITE(6,*) NTYPE
      END

C.. Get a count of the number of excitations in this excitor.  
C.. If it hasn't yet been generated, it is generated, otherwise
C.. it is just returned
      SUBROUTINE GETSYMEXCITCOUNT(NI,NEL,G1,NBASIS,NBASISMAX,
     &         NMEM,IC,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,*),NTYPES
         POINTER (IP_TYPES, TYPES)
         INTEGER I,IC
         INTEGER ISYM(5)
         LOGICAL TSYM
         LOGICAL TSETUP
         INTEGER NMEM(*)
         INTEGER IFRZ(0:NBASIS)
         IF(NMEM(1).EQ.0) THEN
C.. this is the first time we've been called with allocated memory.
C.. NTYPES=0, so we need to find the excitation types again
C.. The first time counts them again (inelegant, but not too
C.. time-consuming), and the second time saves them
            CALL IAZZERO(NMEM(1),10)
C            WRITE(6,*) "SETUP"
C            WRITE(6,*) NMEM(6)
C            CALL WRITEDET(6,NI,NEL,.TRUE.)
            CALL COUNTSYMEXCIT2(NI,NEL,G1,NBASIS,NBASISMAX,NMEM(6),
     &      NMEM(11+(NEL+NBASIS+8)*4),NMEM(1),IFRZ)
         ENDIF
         IC=NMEM(6)
         RETURN
      END

C.. COunt the total number of symmetric excitations - it stores
C.. each individual type's count in the TYPES(5,*)
      SUBROUTINE COUNTSYMEXCIT2(NI,NEL,G1,NBASIS,NBASISMAX,NEX,
     &      TYPES,NTYPES,IFRZ)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,G1(5,NBASIS),NBASIS,NBASISMAX(5,6)
         INTEGER NEX

         INTEGER TYPES(5,*),NTYPES
         INTEGER I,IC
         INTEGER IFRZ(0:NBASIS)
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         NTYPES=0
         CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         CALL GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         IC=0
         DO I=1,NTYPES
           CALL GETSYMEXTYPECOUNT(NI,NEL,G1,NBASIS,NBASISMAX,TYPES,I,
     &         NFULL,NFREE,IFRZ)
           IC=IC+TYPES(5,I)
C            WRITE(6,*) I,IC
         ENDDO
         NEX=IC
         RETURN
      END
C.. Generate a random excitation of TYPE TYPES(1:5,NTYPE)
      SUBROUTINE GSETE_GENRAND(NI,NJ,IC,NEL,G1,NBASIS,NBASISMAX,
     &         TSYM,ISYM,NTYPES,NTYPE,TFULL,TFREE,TTYPE,SPNS,
     &         NFULL,NFREE,
     &         TYPES,ISEED,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,NTYPES),NTYPES,NTYPE
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         INTEGER SPNS(2,2),IC
         LOGICAL TFULL,TFREE,TTYPE
         INTEGER ISYM(5),I,J
         LOGICAL TSYM
         LOGICAL TDONE
         INTEGER ISEED
         INTEGER IFRZ(0:NBASIS)
         TDONE=.FALSE.
         NJ(1)=0
C.. Initialize this type
         CALL GSETE_INIT(NI,NEL,G1,NBASIS,NFULL,NFREE,
     &               TYPES,NTYPE,NBASISMAX,.FALSE.,IFRZ)
C.. Find a random full set to excite
         CALL GSETE_RANDFLFR(TYPES(1,NTYPE),SPNS(1,1),NFULL,
     &               TFULL,ISEED,IC,.TRUE.)
C.. Find a random empty set to excite
         CALL GSETE_RANDFLFR(TYPES(1,NTYPE),SPNS(1,2),NFREE,
     &                  TFREE,ISEED,IC,.FALSE.)
C.. Get the next excitation
         CALL GSETE_NEXTDET(SPNS,NFREE,NFULL,TSYM,ISYM,
     &                  NI,NJ,NEL,IC,G1,NBASISMAX)

C.. Debugging
C         CALL GSETE_DUMP(56,NI,NEL,G1,NBASIS,NBASISMAX,
C     &         NTYPES,NTYPE,TFULL,TFREE,TTYPE,SPNS,
C     &         NFULL,NFREE,
C     &         TYPES,IFRZ)
         RETURN
      END
C.. Get in NJ the next excitation from the state variables input
C.. If we've run out of excitations in this type, we move to the next
C.. If we've run out of types, we set NJ(1) to 0.
C.. Best not to call directly, but let  GENSYMEXCITIT() do the work
      SUBROUTINE GSETE_GET(NI,NJ,IC,NEL,G1,NBASIS,NBASISMAX,
     &         TSYM,ISYM,NTYPES,NTYPE,TFULL,TFREE,TTYPE,SPNS,
     &         NFULL,NFREE,
     &         TYPES,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYPES(5,NTYPES),NTYPES,NTYPE
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         INTEGER SPNS(2,2),IC
         LOGICAL TFULL,TFREE,TTYPE
         INTEGER ISYM(5),I,J
         LOGICAL TSYM
         LOGICAL TDONE
         INTEGER IFRZ(0:NBASIS)
         TDONE=.FALSE.
         NJ(1)=0
         DO WHILE(.NOT.TDONE)
            IF(.NOT.TTYPE) THEN
C.. We have no more full orbitals to excite, so we need another sym
C.. type
               TFULL=.FALSE.
               IF(NTYPE.LT.NTYPES) THEN
                  NTYPE=NTYPE+1
                  CALL GSETE_INIT(NI,NEL,G1,NBASIS,NFULL,NFREE,
     &               TYPES,NTYPE,NBASISMAX,.TRUE.,IFRZ)
                  TTYPE=.TRUE.
   
               ELSE
                  TDONE=.TRUE.
               ENDIF
            ENDIF
            IF(TTYPE) THEN
C.. Debugging printout
C         WRITE(6,"(9Z5,I)") (TYPES(I,NTYPE),I=1,5),SPNS(1,1),SPNS(2,1),
C     &      SPNS(1,2),SPNS(2,2),NTYPE
C         DO J=1,2
C         WRITE(6,"(A,$)") "FULL"
C         DO I=-3,NFULL(J,1,0)
C            WRITE(6,"(I3,$)") NFULL(J,1,I)
C         ENDDO
C         WRITE(6,"(A,$)") ':'
C         DO I=-3,NFULL(J,2,0)
C            WRITE(6,"(I3,$)") NFULL(J,2,I)
C         ENDDO
C         WRITE(6,*)
C         WRITE(6,"(A,$)") "FREE"
C         DO I=-3,NFREE(J,1,0)
C            WRITE(6,"(I3,$)") NFREE(J,1,I)
C         ENDDO
C         WRITE(6,"(A,$)") ':'
C         DO I=-3,NFREE(J,2,0)
C            WRITE(6,"(I3,$)") NFREE(J,2,I)
C         ENDDO
C         WRITE(6,*)
C         ENDDO
               IF(.NOT.TFULL) THEN
C.. we need to setup a new full set
                  CALL GSETE_SETUPFLFR(TYPES(1,NTYPE),SPNS(1,1),NFULL,
     &               TFULL,.TRUE.)
                  TFREE=.FALSE.
               ENDIF
               IF(TFULL) THEN
                  IF(.NOT.TFREE) THEN
C.. we need to setup the free (virtual) orbitals.
                    CALL GSETE_SETUPFLFR(TYPES(1,NTYPE),SPNS(1,2),NFREE,
     &                  TFREE,.FALSE.)
                  ENDIF
                  IF(TFREE) THEN
C.. Get the next excitation
                     CALL GSETE_NEXTDET(SPNS,NFREE,NFULL,TSYM,ISYM,
     &                  NI,NJ,NEL,IC,G1,NBASISMAX)
C.. Now move to the next free
                    CALL INCCOUNTS(NFREE,SPNS(1,2),TFREE,TYPES(1,NTYPE))
C                     TFREE=.FALSE.
                  ELSE
C.. otherwise mark that we haven't got an excitation
                     NJ(1)=0
                  ENDIF
                  IF(.NOT.TFREE) THEN
C.. We've run out of excitations for this free set, so go round again 
C.. with a new full set.
                    CALL INCCOUNTS(NFULL,SPNS(1,1),TFULL,TYPES(1,NTYPE))
                     IF(.NOT.TFULL) TTYPE=.FALSE.
                  ENDIF
                  IF(NJ(1).NE.0) THEN
C.. otherwise we've found a det
                     TDONE=.TRUE.
                  ENDIF
               ELSE
C.. there're no full available so we move to the next type
                  TTYPE=.FALSE.
               ENDIF
            ENDIF
         ENDDO
C.. If we're here and NJ(1)=0, then we've run out of dets
         RETURN
      END

C.. Initialize state variables for a new symmetry type.
C.. If the symmetry specifications haven't changed (i.e. we're just
C.. considering a different set of spins - workd out by looking at
C.. TYPES(NTYPE-1), and if TOPT is true, we don't need to do anything,
C.. so just return.      
      SUBROUTINE GSETE_INIT(NI,NEL,G1,NBASIS,NFULL,NFREE,TYPES,NTYPE,
     &      NBASISMAX,TOPT,IFRZ)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,G1(5,*),NBASIS,TYPES(5,*),NTYPE,TYP
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         INTEGER NSYM(2),NSYME(2),IEL,I,NSPN,NBASISMAX(5,6),NS,NSE
         LOGICAL TNS,TOPT
         INTEGER ISYM
         INTEGER IFRZ(0:NBASIS)
C.. -2 corresponds to the first electrom to be excited, and -1 to the
C.. second. also for NFREE
C.. NFREE(1,-3) is the spin (either 1 or 2) for the first electron to be 
C.. excited (so NFREE(NFREE(1,-3),-2) is the 1st electron to be excited)
C.. This will also be NSPN1
C.. NFREE(2,-2) will be the second electron - NSPN2
C.. We store the positions of the correct symmetry electrons within NI
C.. in NFULL and NFREE
         IF(NTYPE.GT.1) THEN
            IF(TYPES(2,NTYPE).EQ.TYPES(2,NTYPE-1).AND.TOPT) THEN
C.. The symmetry hasn't changed, so just return
               RETURN
            ENDIF
         ENDIF
         TYP=TYPES(2,NTYPE)

         CALL IAZZERO(NFULL,4*(NEL+3))
         CALL IAZZERO(NFREE,4*(NBASIS+3))
C.. The coding for TYP(2) is:
C.. SSSSRRRR QQQQPPPP for pq->rs (sym1,sym2->syme1,syme2)
C.. 
         NSYM(1)=IAND(TYP,15)-1
         NSYM(2)=IAND(TYP,15*16)/16-1
         NSYME(1)=IAND(TYP,15*256)/256-1
         NSYME(2)=IAND(TYP,15*4096)/4096-1

C.. No symmetry, so all basis fns are allowed 
         IF(TYP.EQ.0) TNS=.TRUE.
C.. 
         IEL=1
C.. Go through the basis looking for electrons of the right sym which
C.. are either filled, or virtuals which can be.
         DO I=1,NBASIS
            DO WHILE(NI(IEL).LT.I.AND.IEL.LT.NEL)
               IEL=IEL+1
            ENDDO
            IF(IAND(IFRZ(0),8).EQ.0.OR.IFRZ(I).EQ.0) THEN
             NSPN=(G1(4,I)+3)/2
             IF(NI(IEL).EQ.I) THEN
C.. a filled electron
               DO ISYM=1,2
                  IF(G1(5,I).EQ.NSYM(ISYM).OR.TNS) THEN
                     NFULL(ISYM,NSPN,0)=NFULL(ISYM,NSPN,0)+1
                     NFULL(ISYM,NSPN,NFULL(ISYM,NSPN,0))=IEL
C                     WRITE(6,*) "S",ISYM,NSPN,NFULL(ISYM,NSPN,0)
                  ENDIF
               ENDDO
             ELSE
C.. a virtual orbital
               DO ISYM=1,2
                  IF(G1(5,I).EQ.NSYME(ISYM).OR.TNS) THEN
                     NFREE(ISYM,NSPN,0)=NFREE(ISYM,NSPN,0)+1
                     NFREE(ISYM,NSPN,NFREE(ISYM,NSPN,0))=I
C                     WRITE(6,*) "S",ISYM,NSPN,NFREE(ISYM,NSPN,0)
                  ENDIF
               ENDDO
             ENDIF
            ENDIF
         ENDDO
         RETURN
      END
C.. Setup the full or free basis fn list pointers.
C.. See GSETE_NEXTDET for the structure
C.. TFULLORFREE is true for full and false for free
C.. NLISTS is either NFULL or NFREE.
C.. TCONT is either TFULL or TFREE.
      SUBROUTINE GSETE_SETUPFLFR(TYP,SPNS,NLISTS,TCONT,TFULLORFREE)
         IMPLICIT NONE
         INTEGER TYP(5),SPNS(2)
         INTEGER NLISTS(2,2,-3:0)
         LOGICAL TCONT
         INTEGER I,J,K
         LOGICAL TFULLORFREE
C.. Free spins are coded as the top four bits.
C.. Full spins are the bottom four
         IF(TFULLORFREE) THEN
            I=1
         ELSE
            I=16
         ENDIF
         SPNS(1)=IAND(TYP(4),3*I)/I
         SPNS(2)=IAND(TYP(4),12*I)/(4*I)
         TCONT=.TRUE.
C.. Zero everything
         DO I=1,2
            DO J=1,2
               DO K=-3,-1
                  NLISTS(I,J,K)=0
               ENDDO
            ENDDO
         ENDDO
         IF(SPNS(1).EQ.SPNS(2).AND.TYP(3).EQ.1) THEN
C.. both excitations from the same spin and sym
            NLISTS(1,SPNS(1),-2)=1
            NLISTS(2,SPNS(2),-1)=2
C.. Flag if there are no possible combinations
            IF(NLISTS(2,SPNS(2),0).LT.2) TCONT=.FALSE.
         ELSEIF(SPNS(2).EQ.0) THEN
            NLISTS(1,SPNS(1),-2)=1
C.. Flag if there are no possible combinations
            IF(NLISTS(1,SPNS(1),0).LT.1) TCONT=.FALSE.
         ELSE
C.. One of one and one of the other
            NLISTS(1,SPNS(1),-2)=1
            NLISTS(2,SPNS(2),-1)=1
C.. Flag if there are no possible combinations
            IF(NLISTS(1,SPNS(1),0).LT.1) TCONT=.FALSE.
            IF(NLISTS(2,SPNS(2),0).LT.1) TCONT=.FALSE.
         ENDIF
         RETURN
      END
C.. Setup the full or free basis fn list pointers to a random value.
C.. See GSETE_NEXTDET for the structure
C.. TFULLORFREE is true for full and false for free
C.. NLISTS is either NFULL or NFREE.
C.. TCONT is either TFULL or TFREE.
      SUBROUTINE GSETE_RANDFLFR(TYP,SPNS,NLISTS,TCONT,ISEED,IC,
     &      TFULLORFREE)
         IMPLICIT NONE
         INTEGER TYP(5),SPNS(2)
         INTEGER NLISTS(2,2,-3:0)
         LOGICAL TCONT,TFULLORFREE
         INTEGER I,J,K,ISEED,IC
         REAL*8 RAN2
C.. Free spins are coded as the top four bits.
C.. Full spins are the bottom four
         IF(TFULLORFREE) THEN
            I=1
         ELSE
            I=16
         ENDIF
         SPNS(1)=IAND(TYP(4),3*I)/I
         SPNS(2)=IAND(TYP(4),12*I)/(4*I)
         TCONT=.TRUE.
C.. We don't need to flag if there are no possible combinations
C.. As we shouldn't be here if there are not.
         DO I=1,2
            DO J=1,2
               DO K=-3,-1
                  NLISTS(I,J,K)=0
               ENDDO
            ENDDO
         ENDDO
         IF(SPNS(1).EQ.SPNS(2).AND.TYP(3).EQ.1) THEN
C.. both excitations from the same spin and sym
C.. so we need a more convoluted algorithm
            J=NLISTS(1,SPNS(1),0)
            IC=J*(J-1)/2
            I=RAN2(ISEED)*J*(J-1)/2+1
            J=J-1
            K=1
            DO WHILE(I.GT.J)
               K=K+1
               I=I-J
               J=J-1
            ENDDO
            I=I+K
C.. K,I is the new pair
            NLISTS(1,SPNS(1),-2)=K
            NLISTS(2,SPNS(2),-1)=I
C            WRITE(6,*) "IC1"
         ELSEIF(SPNS(2).EQ.0) THEN
            NLISTS(1,SPNS(1),-2)=RAN2(ISEED)*NLISTS(1,SPNS(1),0)+1
            IC=NLISTS(1,SPNS(1),0)
C            WRITE(6,*) "IC2"
         ELSE
C.. One of one and one of the other
            NLISTS(1,SPNS(1),-2)=RAN2(ISEED)*NLISTS(1,SPNS(1),0)+1
            NLISTS(2,SPNS(2),-1)=RAN2(ISEED)*NLISTS(2,SPNS(2),0)+1
            IC=NLISTS(1,SPNS(1),0)*NLISTS(2,SPNS(2),0)
C            WRITE(6,*) "IC3"
         ENDIF
         RETURN
      END
C.. Given list pointers, get the next excited determinant, or set
C.. NJ(1)=0 if we generate one whose symmetry's wrong (if TSYM is set)
C..
C.. NXXXX(a,b,-3) is  is the spin (either 1 or 2) for the first electron to be 
C.. excited (so NFREE(NFREE(1,-3),-2) is the 1st electron to be excited)

C.. NXXXX(a,b,-2) corresponds to the first basis fn to be excited, and -1 to the
C.. second.

C.. NXXXX(a,b,0) is the number of basisfns available in this group
C.. NXXXX(a,b,1:...) are the basisfns in this group.  for FREE, it is
C..                  the actual basis fn, and FULL it is the index
C..                  within NI.

C.. a corresponds to the sym to be excited - 1 for the 1st and 2 for 2nd
C.. b is the spin of the appropriate symmetry - beta=1, alpha=2

C.. This will also be NSPN1
C.. NFREE(2,-2) will be the second electron - NSPN2
C.. We store the positions of the correct symmetry electrons within NI
C.. in NFULL and NFREE
      SUBROUTINE GSETE_NEXTDET(SPNS,NFREE,NFULL,TSYM,ISYM,
     &            NI,NJ,NEL,IC,G1,NBASISMAX)
         IMPLICIT NONE
         INTEGER SPNS(2,2)
         INTEGER NFULL(2,2,-3:0)
         INTEGER NFREE(2,2,-3:0)
         LOGICAL TSYM
         INTEGER ISYM(5),ISYM2(5)
         INTEGER NJ(NEL),NI(NEL),NEL
         LOGICAL LCHKSYM
         INTEGER IC,G1(*),NBASISMAX(*),I
         LOGICAL ISVALIDDET
         CALL ICOPY(NEL,NI,1,NJ,1)
C         WRITE(56,"(I4,A,I4,$)") 
C     &                  NJ(NFULL(1,SPNS(1,1),NFULL(1,SPNS(1,1),-2)))
C     &            ,'->',NFREE(1,SPNS(1,2),NFREE(1,SPNS(1,2),-2))
         NJ(NFULL(1,SPNS(1,1),NFULL(1,SPNS(1,1),-2)))
     &                  =NFREE(1,SPNS(1,2),NFREE(1,SPNS(1,2),-2))
         IF(SPNS(2).NE.0) THEN
            NJ(NFULL(2,SPNS(2,1),NFULL(2,SPNS(2,1),-1)))
     &                  =NFREE(2,SPNS(2,2),NFREE(2,SPNS(2,2),-1))
C         WRITE(56,"(I4,A,I4,$)") 
C     &                  NJ(NFULL(2,SPNS(2,1),NFULL(2,SPNS(2,1),-1)))
C     &            ,'->',NFREE(2,SPNS(2,2),NFREE(2,SPNS(2,2),-1))
            IC=2
         ELSE
            IC=1
         ENDIF
         CALL SORTI(NEL,NJ)
C         IF(NJ(1).EQ.0) WRITE(6,*) "ERROR"
         IF(.NOT.ISVALIDDET(NJ,NEL)) THEN
            WRITE(6,*) "Invalid det in GSETE_NEXTDET"
            CALL WRITEDET(6,NJ,NEL,.TRUE.)
            STOP "Invalid Determinant in GSETE_NEXTDET"
         ENDIF
         IF(TSYM) THEN
C.. if the symmetry's wrong, we remove the excitation
            CALL GETSYM(NJ,NEL,G1,NBASISMAX,ISYM2)
C.. by setting the first electron to be 0
C         CALL WRITEDET(6,NI,NEL,.FALSE.)
C         CALL WRITEDET(6,NJ,NEL,.TRUE.)
C            WRITE(6,"(A,8I3)") "SYM:", (ISYM(I),I=1,4),
C     &         (ISYM2(I),I=1,4)
            IF(.NOT.LCHKSYM(ISYM,ISYM2)) NJ(1)=0
         ENDIF
         RETURN
      END

C.. For a given Excitation symmetry type, generate all excitatiosn
      SUBROUTINE GETSYMEXTYPEEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSYM,ISYM,
     &         TYP,LIST,ICLIST,NLIST,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS,LIST(NEL,NLIST),ICLIST(NLIST)
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER TYP(5),NSPN,NJ(NEL)
         INTEGER NFULL(2,-3:NEL),NFREE(2,-3:NBASIS),NSYM,NSYME
         INTEGER SPN1,SPN2,IC,ICARRY
         LOGICAL TFULL,TFREE
         INTEGER IEL,I
         INTEGER ISYM(5),ISYM2(5)
         LOGICAL TSYM
         LOGICAL LCHKSYM
         INTEGER IFRZ(0:NBASIS)
C.. -2 corresponds to the first electrom to be excited, and -1 to the
C.. second. also for NFREE
C.. NFREE(1,-3) is the spin (either 1 or 2) for the first electron to be 
C.. excited (so NFREE(NFREE(1,-3),-2) is the 1st electron to be excited)
C.. This will also be NSPN1
C.. NFREE(2,-2) will be the second electron - NSPN2
C.. We store the positions of the correct symmetry electrons within NI
C.. in NFULL and NFREE





C         CALL GSETE_INIT(NI,NEL,G1,NBASIS,NFULL,NFREE,TYP,IFRZ)
         CALL GSETE_INIT(NI,NEL,G1,NBASIS,NFULL,NFREE,TYP,1
     &               ,NBASISMAX,.FALSE.,IFRZ)
C.. Setup the full lists
         CALL GSETE_SETUPFLFR(TYP,SPN1,SPN2,NFULL,TFULL,.TRUE.)

C.. Now work through both lists of full and free, generating all
C.. excitations
         IC=0
         DO WHILE(TFULL)
C.. Setup the free lists
            CALL GSETE_SETUPFLFR(TYP,SPN1,SPN2,NFREE,TFREE,.FALSE.)
            DO WHILE(TFREE)
               IC=IC+1
               IF(IC.LE.NLIST) THEN
                  CALL GSETE_NEXTDET(SPN1,SPN2,NFREE,NFULL,TSYM,ISYM,
     &               NI,LIST(1,IC),NEL,ICLIST(IC),G1,NBASISMAX)
C.. if the symmetry's wrong, we remove the excitation
                  IF(LIST(1,IC).EQ.0) IC=IC-1
               ENDIF
C.. Now move to the next free
               CALL INCCOUNTS(NFREE,SPN1,SPN2,TFREE)
            ENDDO
            CALL INCCOUNTS(NFULL,SPN1,SPN2,TFULL)
         ENDDO
         TYP(5)=IC
         RETURN
      END

C.. Increment the pointer for either full (in which case NFREE=NFULL) or
C.. FREE (in which case NFREE=NFREE) orbitals, setting TFREE(or TFULL)
C.. to false if there are no more combinations of these orbitals
C.. allowed.
      SUBROUTINE INCCOUNTS(NFREE,SPNS,TFREE,TYP)
         IMPLICIT NONE
         INTEGER SPNS(2),NFREE(2,2,-3:0),TYP(5)
         LOGICAL TFREE
         INTEGER ICARRY
C.. NFREE could be NFREE or NFULL
         IF(SPNS(1).EQ.SPNS(2).AND.TYP(3).EQ.1) THEN
C.. spins and sym same
C.. so we have to ensure the two pointers are different:
C.. (1,2),(1,3),...,(1,N),(2,3),(2,4),...,(2,N),...,(N-1,N)
            ICARRY=0
            NFREE(2,SPNS(2),-1)=NFREE(2,SPNS(2),-1)+1
            IF(NFREE(2,SPNS(2),-1).GT.NFREE(2,SPNS(2),0)) THEN
               NFREE(1,SPNS(1),-2)=NFREE(1,SPNS(1),-2)+1
               NFREE(2,SPNS(2),-1)=NFREE(1,SPNS(1),-2)+1
            ENDIF
            IF(NFREE(1,SPNS(1),-2).GT.NFREE(1,SPNS(1),0)-1) THEN
               TFREE=.FALSE.
            ENDIF
         ELSE
C.. Spins and sym different, so we can just go through
C..(1,1),(1,2),...,(1,N),(2,1),(2,2),...,(2,N),...(N,N)
            IF(SPNS(2).NE.0) THEN
               ICARRY=0
               NFREE(2,SPNS(2),-1)=NFREE(2,SPNS(2),-1)+1
               IF(NFREE(2,SPNS(2),-1).GT.NFREE(2,SPNS(2),0)) THEN
                  NFREE(2,SPNS(2),-1)=1
                  ICARRY=1
               ENDIF
            ELSE
               ICARRY=1
            ENDIF
            NFREE(1,SPNS(1),-2)=NFREE(1,SPNS(1),-2)+ICARRY
            IF(NFREE(1,SPNS(1),-2).GT.NFREE(1,SPNS(1),0)) THEN
               TFREE=.FALSE.
            ENDIF
         ENDIF
      END
C.. Enumerate all possible classes of excitations of a det, and store
C.. them in TYPES (if NTYPES>0 on entry), or just count them if
C.. NTYPES=0.
C
C.. each type is a 5-integer set
C.. TYPE(1,*)=0
C.. TYPE(2,*)=symmetries of ab->cd
C..             coded in binary as              DDDDCCCCBBBBAAAA
C..           or 0 if no symmetry at all
C.. TYPE(3,*)=0 normally and 1 is both excitations are from the same sym
C..                                   (and therefore to the same sym)
C.. TYPE(4,*)=spin excitation type
C..                     we code SPN1,SPN2-> SPNE1,SPNE2 (I,J->K,L)
C..                     as LL KK JJ II in binary (beta=1, alpha=2)
C.. TYPE(5,*)=number of excitations of this type (calculated later)
      SUBROUTINE GETSYMEXTYPES(NBASISMAX,TYPES,NTYPES,NI,NEL,G1,NBASIS,
     &      IFRZ)
         IMPLICIT NONE
         INTEGER NBASISMAX(5,6),NBASIS
         INTEGER NTYPES,TYPES(5,NTYPES) 

         INTEGER NTY2
         INTEGER I,J,K,L,ITYPE
         INTEGER TYP(5)
         INTEGER JMIN,JMAX
         INTEGER IEXC
         INTEGER SYMES(2,2,0:8),NSPN,IEL
         INTEGER NI(NEL),NEL,G1(5,*)
         LOGICAL TALL
         LOGICAL ISUHFDET,NOSING
         INTEGER IFRZ(0:NBASIS)
         NTY2=0
C.. Loop over all four possible symmetry labels
C.. (in principle, but given we only take into account spatial and spin
C.. symmetries here, we only loop over indices 2 and 4
         NOSING=ISUHFDET(NI,NEL,NBASISMAX).OR.(IAND(IFRZ(0),4).NE.0)
         TYP(1)=0
         TYP(3)=0
         TYP(5)=0
         IF(NBASISMAX(3,3).EQ.1) THEN
C.. Generic spatial
            JMAX=8
            JMIN=1
C.. J corresponds to the symmetry of orbitals being excited.
C.. For spatial symmetries, either a single a1 (J=0) excitation is
C.. allowed, or doubles of any other symmetry
         ELSE
            JMAX=0
            JMIN=0
         ENDIF
C.. Get all the syms
         CALL IAZZERO(SYMES,36)
         IEL=1
         DO I=1,NBASIS  
            DO WHILE(NI(IEL).LT.I.AND.IEL.LT.NEL)
               IEL=IEL+1
            ENDDO
            IF(NBASISMAX(3,3).EQ.1) THEN
               J=G1(5,I)+1
            ELSE
               J=0
            ENDIF
            IF(NI(IEL).EQ.I) THEN
               L=1
            ELSE
               L=2
            ENDIF
            K=(G1(4,I)+3)/2
            IF(IAND(IFRZ(0),8).EQ.0.OR.IFRZ(I).EQ.0)
     &            SYMES(L,K,J)=SYMES(L,K,J)+1
         ENDDO       
C.. J is the sym of elec 1, and K of elec 2
         DO J=JMIN,JMAX
            DO K=JMIN,JMAX
C.. we need to find all possible symmetry pair combinations
C.. ab->cd
C.. Singles are easy (a->a) ITYPE=0
C.. Doubles for non-sym are classified as ITYPE=-1
C.. Doubles for sym are ITYPE>0
               DO ITYPE=-1,JMAX
                 TALL=.TRUE.
                 IF(ITYPE.LT.1) THEN
C.. ITYPE=-1 doesn't apply to spatial symmetries
                  IF(JMAX.NE.0.AND.ITYPE.EQ.-1) TALL=.FALSE.
C.. Disallow if the occupied and virtual orbital have different sym
                  IF(ITYPE.EQ.0.AND.(J.NE.K.OR.NOSING)) TALL=.FALSE.
                  TYP(2)=256*K+J
                 ELSE
C.. spatial symmetries.
C.. we consider ab->cd.  Here a=J, b=K, c=ITYPE, and d is worked out
C.. DDDDCCCCBBBBAAAA
                  L=IEOR(IEOR(J-1,K-1),ITYPE-1)+1
                  TYP(2)=4096*L+256*ITYPE+16*K+J
                  IF(L.LT.ITYPE.OR.K.LT.J) TALL=.FALSE.
                 ENDIF

                 IF(TALL) THEN
C.. For each possible symmetry subgrouping (e.g. a symmetry element, or
C.. alpha/beta.  We do not consider momentum conservation here, as it is
C.. rather complicated to do.
                  TYP(3)=0
                  IF(ITYPE.EQ.0) THEN
C.. Single excits a0->a0

C.. Beta -> Beta 00010001
                     TYP(4)=17
                     CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
C.. Alpha -> Alpha        00100010
                     TYP(4)=34
                     CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
                  ELSE
C.. ab->cd
C.. if both excitations from the same place, they end up in the same
C.. place, and we need to tell the excitation routines
                     IF(K.EQ.J) TYP(3)=1
C.. we code SPN1,SPN2-> SPNE1,SPNE2 (I,J->K,L)
C.. as LL KK JJ II in binary (beta=1, alpha=2)

C.. beta, beta->beta,beta 01010101
                     TYP(4)='55'X
                     CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
C.. alpha, alpha->alpha alpha 10101010
                     TYP(4)='AA'X
                     CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
C.. beta,alpha->  beta,alpha 10011001
                     TYP(4)='99'X
                     CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
                     IF(TYP(3).NE.1) THEN
C.. the two sets of excitors are different, so we need to consider 
C.. beta, alpha-> alpha, beta  01101001
                        TYP(4)='69'X
                        CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
C.. the two sets of excitors are different, so we need to consider 
C.. alpha,beta-> alpha,beta 01100110
                        TYP(4)='66'X
                        CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
C.. the two sets of excitors are different, so we need to consider 
C.. alpha beta -> beta, alpha 10010110 
                        TYP(4)='96'X
                        CALL GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
                     ENDIF
                  ENDIF
                 ENDIF
               ENDDO
            ENDDO
         ENDDO
         NTYPES=NTY2
         RETURN
      END
C.. If we have space in TYPES, and we calculate that there is a non-zero
C.. number of excitations here, we add the type to TYPES
      SUBROUTINE GSET_ADDTYP(TYP,TYPES,NTY2,NTYPES,SYMES)
         IMPLICIT NONE
         INTEGER TYP(5),TYPES(5,NTYPES),NTY2,NTYPES,SYMES(2,2,0:8)
         INTEGER SPN1,SPN2,NSYM(2),NSYME(2),N,SPNE1,SPNE2
         NSYM(1)=IAND(TYP(2),15)
         NSYM(2)=IAND(TYP(2),15*16)/16
         NSYME(1)=IAND(TYP(2),15*256)/256
         NSYME(2)=IAND(TYP(2),15*4096)/4096
C         WRITE(6,*) "[",NSYM(1),NSYM(2),NSYME(1),NSYME(2),"]"
         N=1
         SPN1=IAND(TYP(4),3)
         SPN2=IAND(TYP(4),12)/4
         SPNE1=IAND(TYP(4),3*16)/16
         SPNE2=IAND(TYP(4),3*64)/64
         IF(NSYM(1).EQ.NSYM(2).AND.SPN1.EQ.SPN2) N=2
         IF(SYMES(1,SPN1,NSYM(1)).GE.N.AND.
     &      SYMES(2,SPNE1,NSYME(1)).GE.N.AND.
     &      (SPN2.EQ.0.OR.
     &      (SYMES(1,SPN2,NSYM(2)).GE.N.AND.
     &       SYMES(2,SPNE2,NSYME(2)).GE.N)))    THEN
            NTY2=NTY2+1
C            WRITE(6,*) "YES"
            IF(NTY2.LE.NTYPES) 
     &         CALL ICOPY(5,TYP,1,TYPES(1,NTY2),1)
         ENDIF
      END

C.. Work out the number of excitations in TYPES(1:5,NTYPE) and store in
C.. TYPES(5,NTYPE).
C.. This does not take into account TSYM exclusions.
      SUBROUTINE GETSYMEXTYPECOUNT(NI,NEL,G1,NBASIS,NBASISMAX,TYPES,
     &   NTYPE,NFULL,NFREE,IFRZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),NBASIS
         INTEGER NLIST,NBASISMAX(5,6)
         INTEGER G1(5,NBASIS)
         INTEGER NFULL(2,2,-3:NEL),NFREE(2,2,-3:NBASIS)
         INTEGER SPNS(2,2),IC
         LOGICAL TFULL,TFREE,TTYPE
         INTEGER ISYM(5),I,J
         LOGICAL TSYM
         LOGICAL TDONE
         INTEGER ISEED,IC1,IC2,TYPES(5,*),NTYPE
         INTEGER IFRZ(0:NBASIS)
         TDONE=.FALSE.
         NJ(1)=0
C.. We have no more full orbitals to excite, so we need another sym
C.. type
         CALL GSETE_INIT(NI,NEL,G1,NBASIS,NFULL,NFREE,
     &               TYPES,NTYPE,NBASISMAX,.TRUE.,IFRZ)
C         CALL GSETE_DUMP(6,NI,NEL,G1,NBASIS,NBASISMAX,
C     &         0,NTYPE,TFULL,TFREE,TTYPE,SPNS,
C     &         NFULL,NFREE,
C     &         TYPES,IFRZ)
C.. we need to setup a new full set
         CALL GSETE_RANDFLFR(TYPES(1,NTYPE),SPNS(1,1),NFULL,
     &               TFULL,ISEED,IC1,.TRUE.)
C.. we need to setup the excited orbitals.
         CALL GSETE_RANDFLFR(TYPES(1,NTYPE),SPNS(1,2),NFREE,
     &                  TFREE,ISEED,IC2,.FALSE.)
C.. Get the next excitation
C      WRITE(6,*) IC1,IC2
         TYPES(5,NTYPE)=IC1*IC2
         RETURN
      END

C.. Generate excitations of order at most NORDER from NI (excluding NI
Citself) in order
C.. ICLIST is the list of orders of the excitations
C.. NMIN is the minimum level of excitation
C.. IF NORDER=0, we SET NORDER=NEL

C.. THis does not use the excitation generators, but enumerates all possible
C.. determinants within double excitations, which will end up fairly inefficient.
      SUBROUTINE GENEXCIT(NI,NORDER,NBASIS,NEL,LIST,ICLIST,NLIST,NMIN,
     &            G1,TSYM,NBASISMAX,TCOUNT)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NORDER,NBASIS,LIST(NEL,*),ICLIST(*)
         INTEGER NEXCIT(NEL),NLIST,NBASISMAX(5,7)
         INTEGER ISYM(5),NO,IFDET
         LOGICAL TSYM,TCOUNT
         INTEGER G1(5,NBASIS),NMIN

         INTEGER STORE(6),ICOUNT,ILEVEL
         INCLUDE 'sym.inc'


         NO=NORDER
         IF(NO.EQ.0) NO=NEL
         IF(TSYM) CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
         IF(NBASISMAX(4,7).EQ.1) THEN
            NLIST=1
C.. If we're actually using CSFs
            CALL GNCSFS(NEL,NBASIS,0,NBASISMAX,
     &         LIST,TCOUNT,G1,TSYM,ISYM(4),TSYM,ISYM,NLIST,IFDET,
     &         .TRUE.,NI,ICLIST,.TRUE.,NO)
         ELSE
            IF(TNEWEXCITATIONS.AND.NORDER.LE.2) THEN
               IF(NORDER.EQ.1) ILEVEL=1
               IF(NORDER.EQ.2) THEN
                  IF(NMIN.EQ.1) THEN
                     ILEVEL=3
                  ELSE
                     ILEVEL=2
                  ENDIF
               ENDIF
C.. use a new excitation routine - THIS DOES NOT PRODUCE AN ORDERED LIST
               CALL SYMSETUPEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &            TCOUNT,ICOUNT,ILEVEL)
               IF(.NOT.TCOUNT) THEN
                  IF(ICOUNT.GT.NLIST) THEN
                     WRITE(6,*) "Exctiation list size",NLIST
                     WRITE(6,*) "Requires size",ICOUNT
                     STOP "EXCITATIONS LIST NOT LARGE ENOUGH"
                  ENDIF
                  CALL SYMGENEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &               LIST,ICLIST,ICOUNT)
               ENDIF
               NLIST=ICOUNT
            ELSE
               NLIST=1
               CALL GENEXCIT_R(NI,NO,NEL,1,NBASIS,LIST,ICLIST,
     &            NLIST,NEXCIT,NO,1,NMIN,G1,ISYM,TSYM,NBASISMAX,TCOUNT)
               NLIST=NLIST-1
            ENDIF
         ENDIF
         RETURN
      END

C.. AT 29/1/04 - Looks like this won't scale very well to lots of
Celectrons
C..      as there will be very few which are within NORDER excitations
Cof
C..      our original det. (Or for that matter a large basis)
C.. Recursively go through each electron ( we're on NELEC), filling
CLIST.
C.. with up to NORDER excitations.
C.. NEXCIT is the determinant being constructed.
C.. NLEFT is the number of possible excitations left
C.. NSTARTFN is the basis fn to start with for this electron
      RECURSIVE SUBROUTINE GENEXCIT_R(NI,NORDER,NEL,NELEC,NBASIS,LIST,
     &         ICLIST,NLISTPOS,NEXCIT,NLEFT,NSTARTFN,NMIN,G1,ISYM,TSYM,
     &         NBASISMAX,TCOUNT)
         IMPLICIT NONE
         INTEGER NEL,NELEC,NORDER,NLISTPOS,NLEFT,NBASIS,NSTARTFN,NNLEFT
         INTEGER NI(NEL),NEXCIT(NEL),LIST(NEL,*),ICLIST(*)
         INTEGER I,J,I1,I2,NMIN,NMAXEX,ISYM2(5),K
         LOGICAL LISINOLD,LSYM,TSYM,LCHKSYM,TCOUNT
         INTEGER G1(5,NBASIS),ISYM(5),NBASISMAX(5,3)
C.. NMAXEX is the maximum number of excitations we're allowed to  have
Cleft
         NMAXEX=NORDER-NMIN
C..  I is the new basis fn for electron NELEC
C.. We see if it's in the original determinant
C.. we need to leave at least NEL-NELEC basis fns for the following
Celectrons
         DO I=NSTARTFN,NBASIS-(NEL-NELEC)
            LISINOLD=.FALSE.
            DO J=1,NEL
               IF(NI(J).EQ.I) LISINOLD=.TRUE.
            ENDDO
            NEXCIT(NELEC)=I
            NNLEFT=NLEFT
            IF(.NOT.LISINOLD) NNLEFT=NNLEFT-1
C.. NNLEFT is the number of excitations left
C.. check to see if we're actually exciting an electron, and if we've
Cany
C.. spare excitations left
C.. if we've allocated the last electron, we need to see if the det
C.. is allowed, and if so, store it
            IF(NELEC.EQ.NEL) THEN
C.. if we've excited at least one elec
               IF(NNLEFT.LE.NMAXEX.AND.NNLEFT.GE.0) THEN
C.. we check whether we're allowed this excitation owing to sym
C                  CALL WRITEDET(6,NEXCIT,NEL,.TRUE.)
                  IF(TSYM) THEN
                     CALL GETSYM(NEXCIT,NEL,G1,NBASISMAX,ISYM2)
                     LSYM=LCHKSYM(ISYM,ISYM2)
                  ELSE
                     LSYM=.TRUE.
                  ENDIF
C                  WRITE(6,*) ISYM,ISYM2,LSYM
                  IF(LSYM) THEN
                     IF(.NOT.TCOUNT) THEN
                        CALL ICOPY(NEL,NEXCIT,1,LIST(1:NEL,NLISTPOS),1)
                        ICLIST(NLISTPOS)=NORDER-NNLEFT
                     ENDIF
                     NLISTPOS=NLISTPOS+1
                  ENDIF
               ENDIF
            ELSEIF(NNLEFT.GE.0) THEN
C.. we start the next electron on the basis fn after this one
               CALL GENEXCIT_R(NI,NORDER,NEL,NELEC+1,NBASIS,LIST,ICLIST,
     &            NLISTPOS,NEXCIT,NNLEFT,I+1,NMIN,G1,ISYM,TSYM,
     &            NBASISMAX,TCOUNT)
            ENDIF
         ENDDO
         RETURN
      END


