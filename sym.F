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


