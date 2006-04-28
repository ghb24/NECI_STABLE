C.. Symmetries can be unset(=0), or have bits set for the different irreps included
C.. bit 0 corresponds to totally symmetric.
C.. To multiply symmetries, for each bit set in each of the two symmetries, we look up the 
C.. product in the symmetry table, and set that bit.
      INTEGER FUNCTION SYMPROD(ISYM1,ISYM2)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER ISYM1, ISYM2
         INTEGER IS1,IS2
         INTEGER I,J
         IF(IS1.EQ.0.OR.IS2.EQ.0) THEN
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
                     SYMPROD=IBSET(SYMPROD,SYMTABLE(I,J))
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

      LOGICAL FUNCTION LSYMSYM(ISYM)
         INTEGER ISYM
         LSYMSYM=(ISYM.EQ.0.OR.ISYM.EQ.1)
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
         

C.. Setup the symmetry product table
         NSYM=NSYMMAX

         CALL MEMORY(IP_SYMTABLE,NSYM*NSYM,'SYMTABLE')
         DO I=0,NSYM-1
            DO J=0,NSYM-1
               SYMTABLE(I+1,J+1)=IBSET(0,IEOR(I,J))
               SYMTABLE(J+1,I+1)=IBSET(0,IEOR(I,J))
               WRITE(6,"(I4,$)")SYMTABLE(I+1,J+1)
            ENDDO
            WRITE(6,*)
         ENDDO

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

      SUBROUTINE FAKEREADCHARTABLE(NHG)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
C.. First set up the labels and chars
         INTEGER CHARS(16,7)
         DATA ((CHARS(J,I),J=1,16),I=1,7)
     &               /1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &               2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2,
     &               2, 0, 0,-2, 0, 0, 0, 0,-2, 0, 0, 2, 0, 0, 0, 0,
     &               1,-1,-1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1, 1,-1,-1,
     &               4, 0, 0,-4, 0, 0, 0, 0,-4, 0, 0, 4, 0, 0, 0, 0,
     &               1,-1,-1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,
     &               2,-2,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2/
         INTEGER LABS(19)
         DATA LABS /1,2,2,1,3,3,4,2,2,1,5,5,5,5,
     &              3,3,6,7,7/
         INTEGER NHG,I,J
         INTEGER ARRAY(*)
         POINTER (IP_ARRAY,ARRAY)
         NSYMLABELS=7
         NROT=16
         NHG=19
         CALL MEMORY(IP_SYMLABELS,NHG,"SYMLABELS")
         CALL MEMORY(IP_SYMLABELCHARS,NROT*NSYMLABELS,"SYMLABELCH")
         IP_ARRAY=IP_SYMLABELCHARS
         CALL MEMORY(IP_IRREPCHARS,NROT*NSYMLABELS,"IRREPCH")
         DO I=1,7
            DO J=1,16
               ARRAY((I-1)*16+J)=CHARS(J,I)
            ENDDO
         ENDDO
         WRITE(6,*) "SYMMETRY CLASSES"
         CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
      END

C.. Given a set of characters of states, generate all relevant irreps which span the set of characters.
      SUBROUTINE GENIRREPS
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J,K
         LOGICAL LDO,LDO2
         INTEGER IDECOMP
         INTEGER NEXTSYMLAB
         LOGICAL GETIRREPDECOMP
         INTEGER REPCHARS(NROT,NSYMLABELS),NREPS,NORM
         NREPS=0
C.. Initialize the table with the totally symmetric rep.
         DO I=1,NROT
            IRREPCHARS(I,1)=1
         ENDDO
         NSYM=1
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
                     REPCHARS(K,NREPS)=IRREPCHARS(K,I)*IRREPCHARS(K,J)
                  ENDDO
C                  CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDPRD")
                  IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                     IF(NORM.EQ.NROT) THEN
C.. if it's an irrep
                        NSYM=NSYM+1
                        DO K=1,NROT
                           IRREPCHARS(K,NSYM)=REPCHARS(K,NREPS)
                        ENDDO
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
                  IF(NORM.EQ.NROT) THEN
C.. if it's an irrep
                     NSYM=NSYM+1
                     DO I=1,NROT
                        IRREPCHARS(I,NSYM)=REPCHARS(I,NREPS)
                     ENDDO
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
               ENDIF
            END DO lp2
         ENDDO
C.. 
         WRITE(6,*) "IRREP TABLE"
         CALL WRITEIRREPTAB(6,IRREPCHARS,NROT,NSYM)
         IF(NREPS.GT.0) THEN
            WRITE(6,*) NREPS," non-reducible"
            CALL WRITEIRREPTAB(6,REPCHARS,NROT,NREPS)
            STOP "Non-reducible reps found."
         ENDIF
      END


C.. Display irrep table      
      SUBROUTINE WRITEIRREPTAB(IUNIT,IRREPCHARS,NROT,NSYM)
         IMPLICIT NONE
         INTEGER IUNIT,NROT,NSYM,IRREPCHARS(NROT,NSYM)
         INTEGER I,J
         DO I=1,NSYM
            WRITE(IUNIT,"(A,I3,A,$)"), "SYM ", I,":   "
            DO J=1,NROT
               WRITE(IUNIT,"(I3,$)") IRREPCHARS(J,I)
            ENDDO
            WRITE(IUNIT,*)
         ENDDO
         WRITE(IUNIT,*)
      END 
C.. Display a line of characters
      SUBROUTINE WRITECHARS(IUNIT,CHARS,NROT,STR)
         IMPLICIT NONE
         INTEGER IUNIT,NROT,CHARS(NROT)
         INTEGER I,J
         CHARACTER*6 STR
            WRITE(IUNIT,"(A6,A,$)"), STR,":   "
            DO J=1,NROT
               WRITE(IUNIT,"(I3,$)") CHARS(J)
            ENDDO
            WRITE(IUNIT,*)
      END 
C.. Decompose rep CHARS into irreps in IRREPCHARS.  Bit 0 in IDECOMP corresponds to the first irrep etc.
C.. CHARS at the end contains the remainder after the decomposition.
C.. Return .FALSE. if the decomposition is complete and CHARS contains only 0.
      LOGICAL FUNCTION GETIRREPDECOMP(CHARS,IRREPCHARS,NIRREPS,NROT,
     &         IDECOMP,CNORM)
         IMPLICIT NONE
         INTEGER NIRREPS, NROT,IDECOMP,IRREPCHARS(NROT,NIRREPS)
         INTEGER CHARS(NROT),CNORM
         INTEGER I,J,TOT,NORM
         IDECOMP=0
C,. First check norm of this state
         CNORM=0
         DO J=1,NROT
            CNORM=CNORM+CHARS(J)*CHARS(J)
         ENDDO
         DO I=1,NIRREPS
            TOT=0
            DO J=1,NROT
               TOT=TOT+IRREPCHARS(J,I)*CHARS(J)
            ENDDO
            IF(TOT.NE.0) THEN
C.. Calculate the normalization of the state I which matches (if it's an irrep, this will be 1)
               NORM=0
               DO J=1,NROT
                  NORM=NORM+IRREPCHARS(J,I)*IRREPCHARS(J,I)
               ENDDO
C               WRITE(6,*) "IRREP ",I,(TOT+0.D0)/NORM
               IF(MOD(TOT,NORM).NE.0.AND.CNORM.EQ.NROT) THEN
C.. The given representation CHARS has fewer irreps in it than the one in IRREPCHARS, and is an irrep
C.. Hurrah!  Remove it from the one in IRREPCHARS, and keep on going)
                  DO J=1,NROT
                    IRREPCHARS(J,I)=IRREPCHARS(J,I)-CHARS(J)*TOT/CNORM
                  ENDDO
C                  CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"NOW   ")
               ELSEIF(MOD(TOT,NORM).EQ.0) THEN
C.. We've found an (ir)rep which is wholly in CHARS
                  IDECOMP=IBSET(IDECOMP,I-1)
                  CNORM=0
                  DO J=1,NROT
                     CHARS(J)=CHARS(J)-(IRREPCHARS(J,I)*TOT)/NORM
                     CNORM=CNORM+CHARS(J)*CHARS(J)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         GETIRREPDECOMP=.FALSE.
         DO J=1,NROT
            IF(CHARS(J).NE.0) GETIRREPDECOMP=.TRUE.
         ENDDO
      END


      SUBROUTINE GENSYMTABLE
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J,K
         LOGICAL GETIRREPDECOMP
         INTEGER CHARS(NROT),IDECOMP,CNORM
         CALL MEMORY(IP_SYMTABLE,NSYM*NSYM,'SYMTABLE')
         DO I=1,NSYM
            DO J=I,NSYM
               DO K=1,NROT
                  CHARS(K)=IRREPCHARS(K,I)*IRREPCHARS(K,J)
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

