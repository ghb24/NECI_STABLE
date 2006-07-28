C.. Symmetries can be unset(=0), or have bits set for the different irreps included
C.. bit 0 corresponds to totally symmetric.
C.. To multiply symmetries, for each bit set in each of the two symmetries, we look up the 
C.. product in the symmetry table, and OR that with the product.
      INTEGER FUNCTION SYMPROD(ISYM1,ISYM2)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER ISYM1, ISYM2
         INTEGER IS1,IS2
         INTEGER I,J
         IF(ISYM1.EQ.0.OR.ISYM2.EQ.0) THEN
            SYMPROD=0
            RETURN
         ENDIF
         IF(IP_SYMTABLE.EQ.0) STOP 'SYMMETRY TABLE NOT ALLOCATED'
         IS1=ISYM1
         I=1
         SYMPROD=0
         DO WHILE(IS1.NE.0)
            IF(IAND(IS1,1)) THEN
               IS2=ISYM2
               J=1
               DO WHILE(IS2.NE.0)
                  IF(IAND(IS2,1)) THEN
                     SYMPROD=IOR(SYMPROD,SYMTABLE(I,J))
                  ENDIF
                  IS2=RSHIFT(IS2,1)
                  J=J+1
               ENDDO
            ENDIF
            IS1=RSHIFT(IS1,1)
            I=I+1
         ENDDO
         RETURN
      END

      SUBROUTINE WRITESYMTABLE(IUNIT)
         IMPLICIT NONE
         INTEGER IUNIT,I,J
         INCLUDE 'sym.inc'
         DO I=0,NSYM-1
            DO J=0,NSYM-1
               WRITE(IUNIT,"(Z5,$)")SYMTABLE(I+1,J+1)
            ENDDO
            WRITE(IUNIT,*)
         ENDDO
      END

      LOGICAL FUNCTION LSYMSYM(ISYM)
         INTEGER ISYM
         LSYMSYM=(ISYM.EQ.0.OR.BTEST(ISYM,0))
         RETURN
      END


C.. Generate a symmetry table with molpro symmetries.
C.. irreps are simply made out of up to three generators.
C.. MOLPRO classifies irreps with bits 0-3 corresponding to those generators.
C.. symmetry products are merely exclusive ors of the molpro irrep numbers
C.. We set each of the MOLPRO irreps to a bit in our symmetry specifier.
C.. A1 corresponds to bit 0 (i.e. irrep 1)
      SUBROUTINE GENMOLPSYMTABLE(NSYMMAX,G1,NBASIS,ARR,BRR)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER NSYMMAX
         INTEGER I,J
         INTEGER G1(5,NBASIS),NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS)
         INTEGER SYMT2(NSYMMAX,NSYMMAX)
         POINTER (IP_SYMT2,SYMT2)
         

C.. Setup the symmetry product table
         NSYM=NSYMMAX

         CALL MEMORY(IP_SYMTABLE,NSYM*NSYM,'SYMTABLE')
         IP_SYMT2=IP_SYMTABLE
         CALL IAZZERO(SYMTABLE,NSYM*NSYM)
         DO I=0,NSYM-1
            DO J=0,NSYM-1
               SYMT2(I+1,J+1)=IBSET(0,IEOR(I,J))
            ENDDO
         ENDDO

C.. Now generate a list of sym labels.
         CALL MEMORY(IP_SYMLABELS,NBASIS*2,"SYMLABELS")
         NSYMLABELS=NSYM
         DO I=1,NBASIS,2
C.. place the sym label of each state in SYMLABELS(1,ISTATE).  For molp sym, this is 
C.. the log_2 of the symmetry bit string
            SYMLABELS(1,(I+1)/2)=INT(DLOG(G1(5,I)+0.D0)/DLOG(2.D0)+1.4)
         ENDDO
C.. list the symmetry string of each sym label
         DO I=1,NSYM
            SYMLABELS(2,I)=2**(I-1)
         ENDDO
      END

C.. Freeze the SYM LABELS
C.. NHG is the old number of orbs
C.. NBASIS is the new number of orbs
C.. GG(I) is the new index of the (old) orb I

      SUBROUTINE FREEZESYMLABELS(NHG,NBASIS,GG)
         IMPLICIT NONE
         INTEGER NHG,NBASIS,GG(NHG)
         INTEGER I
         INTEGER NSL(NBASIS)
         INCLUDE 'sym.inc'
C.. SYMLABELS is used to classify all states which transform with the same symmetry
C.. for the excitation generation routines
C.. Each state's symmetry falls into a class SYMLABELS(1,ISTATE).
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SYMLABELS(2,ISYMLABEL)
C.. The characters of this class are stored in SYMLABELCHARS(1:NROT, SYMLABELS(1,ISTATE))
C.. The total number of symmetry labels is NSYMLABELS
         DO I=1,NHG,2
            IF(GG(I).NE.0) THEN
               NSL((GG(I)+1)/2)=SYMLABELS(1,(I+1)/2)
            ENDIF
         ENDDO
         DO I=1,NBASIS/2
            SYMLABELS(1,I)=NSL(I)
C            WRITE(6,*) "SL",I,SYMLABELS(1,I)
         ENDDO
      END

      SUBROUTINE GENMOLPSYMREPS(NSYMMAX,G1,NBASIS,ARR,BRR)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER NSYMMAX
         INTEGER I,J
         INTEGER G1(5,NBASIS),NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS)
         INTEGER SYMT2(NSYMMAX,NSYMMAX)
         POINTER (IP_SYMT2,SYMT2)
         
C.. now work out which reps are degenerate and label them
         CALL MEMORY(IP_SYMREPS,2*NBASIS,'SYMREPS')
         J=0
         DO I=1,NBASIS
            IF(I.GT.1.AND.ARR(I).EQ.ARR(I-1)
     &         .AND.G1(5,BRR(I)).EQ.G1(5,BRR(I-1))) THEN
C.. we have the same degenerate rep as the previous entry
               SYMREPS(2,J)=SYMREPS(2,J)+1
            ELSE
C.. we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ENDIF
            SYMREPS(1,I)=J
         ENDDO
         CALL GENSYMLABELPAIRS(NBASIS/2)
      END

C.. delete a symmetry table if one existed.
      SUBROUTINE ENDSYM()
         INCLUDE 'sym.inc'
         IF(IP_SYMTABLE.NE.0) CALL FREEM(IP_SYMTABLE)
         IP_SYMTABLE=0
         IF(IP_SYMREPS.NE.0) CALL FREEM(IP_SYMREPS)
         IP_SYMREPS=0
         IF(IP_SYMLABELS.NE.0) CALL FREEM(IP_SYMLABELS)
         IP_SYMLABELS=0
         IF(IP_SYMLABELCHARS.NE.0) CALL FREEM(IP_SYMLABELCHARS)
         IP_SYMLABELCHARS=0
         IF(IP_IRREPCHARS.NE.0) CALL FREEM(IP_IRREPCHARS)
         IP_IRREPCHARS=0
         IF(IP_SYMLABELPAIRS.NE.0) CALL FREEM(IP_SYMLABELPAIRS)
         IP_SYMLABELPAIRS=0
         IF(IP_SYMLABELPAIRLIST.NE.0) CALL FREEM(IP_SYMLABELPAIRLIST)
         IP_SYMLABELPAIRLIST=0
         IF(IP_SYMLABELPAIRPRODS.NE.0) CALL FREEM(IP_SYMLABELPAIRPRODS)
         IP_SYMLABELPAIRPRODS=0
         IF(IP_SYMLABELLIST.NE.0) CALL FREEM(IP_SYMLABELLIST)
         IP_SYMLABELLIST=0
         IF(IP_SYMLABELCOUNTS.NE.0) CALL FREEM(IP_SYMLABELCOUNTS)
         IP_SYMLABELLIST=0
      END

C.. Precompute a list of the symmetry product of all pairs of symmetry labels
      SUBROUTINE GENSYMLABELPAIRS(NSTATES)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,IDECOMP,TOT,NPRODS,NSTATES
         INTEGER TEMPLIST(NSTATES)
         CALL MEMORY(IP_SYMLABELLIST,NSTATES,'SYMLABELLIST')
         CALL MEMORY(IP_SYMLABELCOUNTS,2*NSYMLABELS,'SYMLABELCOUNTS')
         DO I=1,NSTATES
            SYMLABELLIST(I)=I
            TEMPLIST(I)=SYMLABELS(1,I)
         ENDDO
C.. order according to sym label, so SYMLABELLIST gets a list of states grouped under SYMLABEL
         CALL SORT2I(NSTATES,TEMPLIST,SYMLABELLIST)
         CALL IAZZERO(SYMLABELCOUNTS,2*NSYMLABELS)
         SYMLABELCOUNTS(1,1)=1
         SYMLABELCOUNTS(2,1)=1
         DO I=2,NSTATES
            IF(TEMPLIST(I).NE.TEMPLIST(I-1)) THEN
C.. add a new sym label
               SYMLABELCOUNTS(2,TEMPLIST(I))=1
               SYMLABELCOUNTS(1,TEMPLIST(I))=I
C.. sort the symlabellist
               CALL SORTI(SYMLABELCOUNTS(2,TEMPLIST(I-1)),
     &            SYMLABELLIST(SYMLABELCOUNTS(1,TEMPLIST(I-1))))
            ELSE
               SYMLABELCOUNTS(2,TEMPLIST(I))=
     &            SYMLABELCOUNTS(2,TEMPLIST(I))+1
            ENDIF
         ENDDO
         CALL SORTI(SYMLABELCOUNTS(2,TEMPLIST(I-1)),
     &         SYMLABELLIST(SYMLABELCOUNTS(1,TEMPLIST(I-1))))
C         DO I=1,NSYMLABELS
C            WRITE(6,*) "NSL",I,SYMLABELCOUNTS(1,I),SYMLABELCOUNTS(2,I),
C     &         SYMLABELS(2,I)
C         ENDDO
C.. First allocate the list where we will store the numbers and positions of the sym pairs in the sym pair table
C.. This can get very big if there are very many different reps
         CALL MEMORY(IP_SYMLABELPAIRLIST,2*(2**NSYM),"SYMPAIRLIST")
         CALL IAZZERO(SYMLABELPAIRLIST,2*(2**NSYM))
C.. Now enumerate all pairs, and classify their product, but don't store them.
         NSYMPAIRPRODS=0
         CALL GENALLSYMLABELPAIRS(NSTATES,.FALSE.)
         TOT=0
         DO I=0,2**NSYM-1
            SYMLABELPAIRLIST(1,I)=TOT
            TOT=TOT+SYMLABELPAIRLIST(2,I)
            SYMLABELPAIRLIST(2,I)=0
         ENDDO
C.. I=SYMLABELPAIRLIST(1,SYMPROD) is the index of the set of symlabels with sym prod SYMPROD in
C.. SYMLABELPAIRS
C.. SYMLABELPAIRLIST(2,SYMPROD) is 0 (but will be the total number of pairs with that SYMPROD)
C.. [ SYMLABELPAIRS(1,I) , SYMLABELPAIRS(2,I) ] is the pair of states whose prod is the symlabel.
C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available.

         CALL MEMORY(IP_SYMLABELPAIRS,2*TOT,'SYMLABPS')
         CALL MEMORY(IP_SYMPAIRPRODS,NSYMPAIRPRODS,'SYMPAIRPRODS')
         WRITE(6,*) TOT," SYM LABEL PAIRS"
         WRITE(6,*) NSYMPAIRPRODS," DISTINCT ORBITAL PAIR PRODUCTS"
         NSYMPAIRPRODS=0
         CALL GENALLSYMLABELPAIRS(NSTATES,.TRUE.)
C         DO I=1,NSYMPAIRPRODS
C            WRITE(6,*) "NSPP",I,SYMLABELPAIRLIST(1,SYMPAIRPRODS(I)),
C     &         SYMLABELPAIRLIST(2,SYMPAIRPRODS(I))
C         ENDDO
      END

      SUBROUTINE GENALLSYMLABELPAIRS(NSTATES,TSTORE)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         LOGICAL TSTORE
         INTEGER I,J,PROD
         INTEGER SYMPROD
         INTEGER NSTATES
         DO I=1,NSTATES
CSYMLABELS
            DO J=I,NSTATES
CSYMLABELS
C               WRITE(6,*) I,J,SYMLABELS(1,I),SYMLABELS(1,J)
               PROD=SYMPROD(SYMLABELS(2,SYMLABELS(1,I)),
     &               SYMLABELS(2,SYMLABELS(1,J)))
               IF(TSTORE) THEN
C.. put the pair into the list of pairs.
                  SYMLABELPAIRS(1,SYMLABELPAIRLIST(1,PROD)+
     &               SYMLABELPAIRLIST(2,PROD))=I
                  SYMLABELPAIRS(2,SYMLABELPAIRLIST(1,PROD)+
     &               SYMLABELPAIRLIST(2,PROD))=J
C                  WRITE(6,"(2I4,3Z4,I4)") I,J,
C     &               SYMLABELS(2,SYMLABELS(1,I)),
C     &               SYMLABELS(2,SYMLABELS(1,J)),PROD,
C     &               SYMLABELPAIRLIST(2,PROD)
               ENDIF
C.. If the counter in the pairlist is zero, then we've found a new sym label pair product.
               IF(SYMLABELPAIRLIST(2,PROD).EQ.0) THEN
                  NSYMPAIRPRODS=NSYMPAIRPRODS+1
                  IF(TSTORE) SYMPAIRPRODS(NSYMPAIRPRODS)=PROD
C                  WRITE(6,*) "SPP",NSYMPAIRPRODS,PROD
               ENDIF
C.. incrememnt the counter in the pairlist
               SYMLABELPAIRLIST(2,PROD)=SYMLABELPAIRLIST(2,PROD)+1
            ENDDO
         ENDDO
      END


C..   SYMREPS(2,*) has two sets of data:
C.. SYMREPS(1,IBASISFN) contains the numnber of the representation
C.. of which IBASISFN is a part.
C.. SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
C.. The new method does the following:
C.. Identify all the completely filled reps.
C..   use ADDELECSYM to add together the momenta of these.
C..   These together are totally symmetric
C.. Identify all part-filled degenerate non-reduced representations
C..   Use SYMPROD and ADDELECSYM to generate the resultant symmetry of these
      SUBROUTINE GETSYM(NI2,NEL,G1,NBASISMAX,ISYM)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),G1(5,*),NBASISMAX(5,3),ISYM(5)
         INTEGER I,J,NI2(NEL)
         INTEGER NREPS(NEL),NELECS(NEL),SSYM
         INTEGER SYMPROD
         INCLUDE 'sym.inc'
         I=1
         CALL IAZZERO(NREPS,NEL)
         CALL SETUPSYM(NBASISMAX,ISYM)
         DO I=1,NEL
            CALL GETUNCSFELEC(NI2(I),NI(I),SSYM)
         ENDDO
         DO I=1,NEL
C.. Count all electrons in each rep
C.. NREPS(J) is the rep, and NELECS(J) is the number of electrons in that rep
            J=1
            DO WHILE(NREPS(J).NE.0.AND.J.LT.NEL)
               IF(NREPS(J).EQ.SYMREPS(1,NI(I))) THEN
C.. We've found the slot for the rep.  increment it and leave.
                  NELECS(J)=NELECS(J)+1
                  J=NEL
               ENDIF
               J=J+1
            ENDDO
            IF(J.LE.NEL) THEN
C.. need to put the new rep in a new space
               NREPS(J)=SYMREPS(1,NI(I))
               NELECS(J)=1
            ENDIF
         END DO
C.. now go through and see which are closed and which open
         DO I=1,NEL
            J=1
            DO WHILE(NREPS(J).NE.SYMREPS(1,NI(I)))
               J=J+1
            ENDDO
C.. electron NI(I) is in rep NREPS(J)
            IF(NELECS(J).NE.SYMREPS(2,NREPS(J))) THEN
C.. we don't have a closed shell
C.. add the sym product
               ISYM(5)=SYMPROD(ISYM(5),G1(5,NI(I)))
            ENDIF
C.. add the momentum
            CALL ADDELECSYM(NI(I),G1,NBASISMAX,ISYM)
         ENDDO
C.. round the momentum
         CALL ROUNDSYM(ISYM,NBASISMAX)
         RETURN
      END

C.. Given a set of characters of states, generate all relevant irreps which span the set of characters.
      SUBROUTINE GENIRREPS(IMPROPER_OP,NROTOP)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J,K
         LOGICAL LDO,LDO2
         INTEGER IDECOMP
         INTEGER NEXTSYMLAB
         LOGICAL GETIRREPDECOMP
         COMPLEX*16 REPCHARS(NROT,NSYMLABELS)
         INTEGER NREPS,NROTOP
         REAL*8 NORM
         LOGICAL INV,IMPROPER_OP(NROTOP)
         NREPS=0
C.. Initialize the table with the totally symmetric rep.
         INV=.FALSE.
         DO I=1,NROT
            IRREPCHARS(I,1)=1
            IF(IMPROPER_OP(MOD(I-1,NROTOP)+1)) INV=.TRUE.
         ENDDO
         NSYM=1
         IF(INV) THEN
            WRITE(6,*) "Inversion centre detected"
            NSYM=NSYM+1
C.. There's an inversion centre, so we can immediately create an A1u irrep
            DO I=1,NROT
               IF(IMPROPER_OP(MOD(I-1,NROTOP)+1)) THEN
                  IRREPCHARS(I,NSYM)=-1
               ELSE
                  IRREPCHARS(I,NSYM)=1
               ENDIF
            ENDDO
C            CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         ENDIF
         LDO=.TRUE.
         NEXTSYMLAB=1
         LDO2=.TRUE.
         DO WHILE(LDO.OR.LDO2)
C            CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
C            WRITE(6,*) NREPS," non-reducible"
C            CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
C.. First see if all the products of chars are decomposable
            LDO=.FALSE.
         lp1:DO I=1,NSYM
               DO J=I,NSYM
                  NREPS=NREPS+1
                  IF(NREPS.GT.NSYMLABELS) STOP 'TOO MANY REPS'
                  DO K=1,NROT
                     REPCHARS(K,NREPS)=DCONJG(IRREPCHARS(K,I))
     &                                    *IRREPCHARS(K,J)
                  ENDDO
C                  WRITE(6,*) NREPS,"PROD",I,J
C                  CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDPRD")
                  IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                     IF(ABS(NORM-NROT).LE.1.D-2) THEN
C.. if it's an irrep
                        NSYM=NSYM+1
                        IF(NSYM.GT.32) STOP "MORE than 32 irreps"
                        DO K=1,NROT
                           IRREPCHARS(K,NSYM)=REPCHARS(K,NREPS)
                        ENDDO
C                        CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
                        NREPS=NREPS-1
                        LDO=.TRUE.
                        EXIT lp1
                     ELSE
C                        WRITE(6,*) "IDECOMP:", IDECOMP,NORM,"SYMS:",NSYM
C                      CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"REMAIN")
C.. It's not an irrep, but we cannot reduce it.  Store only if we think we've got all the irreps.
                        IF(LDO2) NREPS=NREPS-1
                        NREPS=NREPS-1
                     ENDIF
                  ELSE
C                     WRITE(6,*) "IDECOMP:", IDECOMP
                     NREPS=NREPS-1
                  ENDIF
               END DO
            END DO lp1
C            WRITE(6,*) LDO,NEXTSYMLAB,NSYMLABELS
            IF(LDO) CYCLE
C.. Check to see if the next symlabel's char is decomposable
        lp2: DO WHILE (NEXTSYMLAB.LE.NSYMLABELS)
               NREPS=NREPS+1
               IF(NREPS.GT.NSYMLABELS) STOP 'TOO MANY REPS'
               DO I=1,NROT
                  REPCHARS(I,NREPS)=SYMLABELCHARS(I,NEXTSYMLAB)
               ENDDO
C               CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDST ")
               IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                  IF(ABS(NORM-NROT).LE.1.D-2) THEN
C.. if it's an irrep
                     NSYM=NSYM+1
                     IF(NSYM.GT.32) STOP "MORE than 32 irreps"
                     DO I=1,NROT
                        IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
                     ENDDO
C                     CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
                     NREPS=NREPS-1
                     LDO=.TRUE.
                     EXIT lp2
                  ELSE
C                     WRITE(6,*) "IDECOMP:", IDECOMP,NORM,"SYMS:",NSYM
C                     CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"REMAIN")
C.. It's not an irrep, but we cannot reduce it.  Store only if we think we've got all the irreps.
                     IF(LDO2) NREPS=NREPS-1
                  ENDIF
               ELSE
C                  WRITE(6,*) "IDECOMP:", IDECOMP
                  NREPS=NREPS-1
               ENDIF
               NEXTSYMLAB=NEXTSYMLAB+1
               IF(.NOT.LDO) THEN
C.. We've not manage to add any more irreps, so we have achieved self-consistency.  Do one more pass to check, saving all C.. non-reducible reps
                  LDO=.TRUE.
                  LDO2=.FALSE.
                  NREPS=0
               ENDIF
            END DO lp2
         ENDDO
C.. 
         WRITE(6,*) "IRREP TABLE"
         CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         IF(NREPS.GT.0) THEN
            WRITE(6,*) NREPS," non-reducible"
               CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
C            IF(NREPS.GT.1) THEN
               STOP "More than 1 non-reducible reps found."
C            ENDIF
C.. we can cope with a single reducible rep.
C            NSYM=NSYM+1
C            DO I=1,NROT
C               IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
C            ENDDO
         ENDIF
C..   Classify each of the symlabels with its decomposition into irreps
         DO I=1,NSYMLABELS
            CALL DECOMPOSEREP(SYMLABELCHARS(1,I),IDECOMP)
            SYMLABELS(2,I)=IDECOMP
         ENDDO
      END


C.. Display irrep table      
      SUBROUTINE WRITEIRREPTAB(IUNIT,CHARS,NROT,NSYM)
         IMPLICIT NONE
         INTEGER IUNIT,NROT,NSYM
         COMPLEX*16 CHARS(NROT,NSYM)
         CHARACTER*6 STR
         INTEGER I,J
         LOGICAL LCOMP,LREAL
         LCOMP=.FALSE.
         LREAL=.FALSE.
         DO I=1,NSYM
            DO J=1,NROT
               IF(ABS(REAL(CHARS(J,I))).GT.1.D-2.AND.
     &            ABS(DIMAG(CHARS(J,I))).GT.1.D-2) LCOMP=.TRUE.
              IF(ABS(REAL(CHARS(J,I))-NINT(REAL(CHARS(J,I))).GT.1.D-2))
     &            LREAL=.TRUE.
            IF(ABS(DIMAG(CHARS(J,I))-NINT(DIMAG(CHARS(J,I))).GT.1.D-2))
     &            LREAL=.TRUE.
            ENDDO
         ENDDO
         DO I=1,NSYM
            WRITE(STR,"(A,I3)"), "SYM", I
            CALL WRITECHARSF(IUNIT,CHARS(1,I),NROT,STR,LCOMP,LREAL)
         ENDDO
         WRITE(IUNIT,*)
      END 
C.. Display a line of characters
      SUBROUTINE WRITECHARS(IUNIT,CHARS,NROT,STR)
         IMPLICIT NONE
         INTEGER IUNIT,NROT
         COMPLEX*16 CHARS(NROT)
         INTEGER I,J
         CHARACTER*6 STR
         LOGICAL LCOMP,LREAL
C.. First do a check for the format
            LCOMP=.FALSE.
            LREAL=.FALSE.
            DO J=1,NROT
               IF(ABS(REAL(CHARS(J))).GT.1.D-2.AND.
     &            ABS(DIMAG(CHARS(J))).GT.1.D-2) LCOMP=.TRUE.
               IF(ABS(REAL(CHARS(J))-NINT(REAL(CHARS(J))).GT.1.D-2))
     &            LREAL=.TRUE.
              IF(ABS(DIMAG(CHARS(J))-NINT(DIMAG(CHARS(J))).GT.1.D-2))
     &            LREAL=.TRUE.
            ENDDO
            CALL WRITECHARSF(IUNIT,CHARS,NROT,STR,LCOMP,LREAL)
      END
      SUBROUTINE WRITECHARSF(IUNIT,CHARS,NROT,STR,LCOMP,LREAL)
         IMPLICIT NONE
         INTEGER IUNIT,NROT
         COMPLEX*16 CHARS(NROT)
         INTEGER I,J
         CHARACTER*6 STR
         LOGICAL LCOMP,LREAL
            WRITE(IUNIT,"(A6,A,$)"), STR,":   "
            DO J=1,NROT
               IF(LCOMP) THEN
                  IF(LREAL) THEN
                     WRITE(IUNIT,"(A,2G16.9,A,$)") "(",
     &                  NINT(REAL(CHARS(J))*1000)/1000.D0,
     &                  NINT(DIMAG(CHARS(J))*1000)/1000.D0
     &                  ,")"
                  ELSE
                     WRITE(IUNIT,"(A,2F6.3,A,$)") "(",CHARS(J),")"
                  ENDIF
               ELSE
                  IF(DIMAG(CHARS(J)).GT.1.D-2) THEN
C.. write in terms of I.
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G14.9,$)") " "
                     ELSE                        
                        IF(ABS(DIMAG(CHARS(J))-1.D0).LT.1.D-2) THEN
                           WRITE(IUNIT,"(A,$)") " -I "
                        ELSEIF(ABS(DIMAG(CHARS(J))+1.D0).LT.1.D-2) THEN
                           WRITE(IUNIT,"(A,$)") "  I "
                        ELSE 
                         WRITE(IUNIT,"(I2,A,$)") NINT(DIMAG(CHARS(J))),
     &                     "I "
                        ENDIF
                     ENDIF
                  ELSE
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G14.9,$)") REAL(CHARS(J)), "    "
                     ELSE
                        WRITE(IUNIT,"(I3,$)") NINT(REAL(CHARS(J)))
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            WRITE(IUNIT,*)
      END

C.. Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
C.. CHARS at the end contains the remainder after the decomposition.
      SUBROUTINE DECOMPOSEREP(CHARSIN,IDECOMP)
         IMPLICIT NONE
         INTEGER IDECOMP
         COMPLEX*16 CHARS(NROT),CHARSIN(NROT),TOT
         REAL*8 CNORM
         INTEGER I,J
         REAL*8 NORM,DIFF
         INCLUDE 'sym.inc'
         IDECOMP=0
         CALL DCOPY(NROT*2,CHARSIN,1,CHARS,1)
C         WRITE(6,*) "Decompose Rep"
C         CALL WRITECHARS(6,CHARS,NROT,"REP   ")
C,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+CHARS(J)*CHARS(J)
         ENDDO
         DO I=1,NSYM
            TOT=0
            DO J=1,NROT
               TOT=TOT+DCONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
            IF(TOT.NE.0) THEN
C.. Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+DCONJG(IRREPCHARS(J,I))*IRREPCHARS(J,I)
               ENDDO
C               WRITE(6,*) "IRREP ",I,(TOT+0.D0)/NORM
               DIFF=ABS(TOT-NINT(ABS(TOT/NORM))*NORM)
               IF(DIFF.GT.1.D-2) THEN
                  WRITE(6,*) 'Symmetry decomposition not complete'
                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IRREP ")
                  CALL WRITECHARS(6,CHARS,NROT,"CHARS ")
                  WRITE(6,*) "Dot product: ",(TOT+0.D0)/NORM,TOT,NORM
                  STOP 'Incomplete symmetry decomposition'
C.. The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
C.. Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
               ELSEIF(ABS(TOT).GT.1.D-2) THEN
C.. We've found an (ir)rep which is wholly in CHARS
                  IDECOMP=IBSET(IDECOMP,I-1)
                  CNORM=0
C                  WRITE(6,*) I,DIFF,TOT,TOT/NORM
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+DCONJG(CHARS(J))*CHARS(J)
                  ENDDO
C                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"DIRREP")
C                  CALL WRITECHARS(6,CHARS,NROT,"DCHARS")
               ENDIF
            ENDIF
         ENDDO
      END
   
 
C.. Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
C.. CHARS at the end contains the remainder after the decomposition.
C.. Return .FALSE. if the decomposition is complete and CHARS contains only 0.
C.. This is used internally in the symmetry routine and destroys CHARS.  For general decomposition,
C,, use DECOMPOSEREP
      LOGICAL FUNCTION GETIRREPDECOMP(CHARS,IRREPCHARS,NIRREPS,NROT,
     &         IDECOMP,CNORM)
         IMPLICIT NONE
         INTEGER NIRREPS, NROT,IDECOMP
         COMPLEX*16 IRREPCHARS(NROT,NIRREPS),CHARS(NROT)
         REAL*8 CNORM, NORM,DIFF
         COMPLEX*16 TOT
         INTEGER I,J
         IDECOMP=0
C,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+DCONJG(CHARS(J))*CHARS(J)
         ENDDO
         DO I=1,NIRREPS
            TOT=0
            DO J=1,NROT
               TOT=TOT+DCONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
            IF(ABS(TOT).GE.1.D-2) THEN
C.. Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+DCONJG(IRREPCHARS(J,I))*IRREPCHARS(J,I)
               ENDDO
C               WRITE(6,*) "IRREP ",I,(TOT+0.D0)/NORM
C                CALL WRITECHARS(6,CHARS,NROT,"REP   ")
C                CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IRREP ")
               DIFF=ABS(TOT-NINT(ABS(TOT/NORM))*NORM)
               IF(DIFF.GE.1.D-2
     &            .AND.CNORM.EQ.NROT) THEN
C.. The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
C.. Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
C                  DO J=1,NROT
C                    IRREPCHARS(J,I)=IRREPCHARS(J,I)-CHARS(J)*TOT/CNORM
C                  ENDDO
C                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"NOW   ")
               ELSEIF(DIFF.LT.1.D-2) THEN
C.. We've found an (ir)rep which is wholly in CHARS
                  IDECOMP=IBSET(IDECOMP,I-1)
                  CNORM=0
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+DCONJG(CHARS(J))*CHARS(J)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         GETIRREPDECOMP=.FALSE.
         DO J=1,NROT
            IF(ABS(CHARS(J)).GT.1.D-2) GETIRREPDECOMP=.TRUE.
         ENDDO
      END


      SUBROUTINE GENSYMTABLE
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J,K
         LOGICAL GETIRREPDECOMP
         COMPLEX*16 CHARS(NROT)
         INTEGER IDECOMP
         REAL*8 CNORM
         CALL MEMORY(IP_SYMTABLE,NSYM*NSYM,'SYMTABLE')
         DO I=1,NSYM
            DO J=I,NSYM
               DO K=1,NROT
                  CHARS(K)=DCONJG(IRREPCHARS(K,I))*IRREPCHARS(K,J)
               ENDDO
               IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,
     &            IDECOMP,CNORM)) THEN
                  WRITE(6,*) "Multiplication of SYMS ",I,J,
     &               " not reducible,"
                  CALL WRITECHARS(6,CHARS,NROT,"REMAIN")
                  STOP "Symmetry table element not reducible"
               ENDIF
               SYMTABLE(I,J)=IDECOMP
               SYMTABLE(J,I)=IDECOMP
C               WRITE(6,"(2I3,B12)") I,J,IDECOMP
            ENDDO
         ENDDO
      END

      SUBROUTINE GENSYMREPS(G1,NBASIS,ARR,BRR,DEGENTOL)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J
         INTEGER G1(5,NBASIS),NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS,2)
         REAL*8 DEGENTOL

C.. now work out which reps are degenerate and label them
         CALL MEMORY(IP_SYMREPS,2*NBASIS,'SYMREPS')
         J=0
         DO I=1,NBASIS
C            WRITE(6,*) "SR2",I
            IF(I.GT.1.AND.ABS(ARR(I,2)-ARR(I-1,2)).LT.DEGENTOL
     &         .AND.G1(5,I).EQ.G1(5,I-1)) THEN
C.. we have the same degenerate rep as the previous entry
               SYMREPS(2,J)=SYMREPS(2,J)+1
            ELSE
C.. we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ENDIF
            SYMREPS(1,I)=J
         ENDDO
C         DO I=1,NBASIS
C            WRITE(6,*) "SR1",SYMREPS(1,I),SYMREPS(2,I)
C         ENDDO   
      END

