C.. Symmetries can be unset(=0), or have bits set for the different irreps included
C.. bit 0 corresponds to totally symmetric.
      
! JSS: if Abelian symmetry (currently just k-point symmetry), then
! symmetry=0 corresponds to the symmetric representation.  It is not
! possible to unset symmetries for k-point jobs.
      
C.. To multiply symmetries, for each bit set in each of the two symmetries, we look up the 
C.. product in the symmetry table, and OR that with the product.

! JSS: In Abelian (k-point) symmetry, a representation can be described by "quantum" numbers.  
! The multiplication of two irreps is equal to the sum of such vectors
! describing the irreps, subject to a modulo of the periodic conditions
! due to the size of the symmetry cell (which corresponds to the k-point
! mesh).
         FUNCTION SYMPROD(ISYM1,ISYM2)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         TYPE(Symmetry) ISYM1, ISYM2
         TYPE(Symmetry) SYMPROD
         TYPE(Symmetry) IS1,IS2
         INTEGER I,J,Abel1(3),Abel2(3)
         integer*8  ComposeAbelianSym
         if (TAbelian) then
             call DecomposeAbelianSym(ISym1%s,Abel1)
             call DecomposeAbelianSym(ISym2%s,Abel2)
             do i=1,3
               Abel1(i)=modulo(Abel1(i)+Abel2(i),NProp(i))
             end do
             SymProd%s=ComposeAbelianSym(Abel1)
         else
             IF(ISYM1%s.EQ.0.OR.ISYM2%s.EQ.0) THEN
                SYMPROD%s=0
                RETURN
             ENDIF
             IF(IP_SYMTABLE.EQ.0) STOP 'SYMMETRY TABLE NOT ALLOCATED'
             IS1=ISYM1
             I=1
             SYMPROD%s=0
             DO WHILE(IS1%s.NE.0)
                IF(IAND(IS1%s,1)) THEN
                   IS2=ISYM2
                   J=1
                   DO WHILE(IS2%s.NE.0)
                      IF(IAND(IS2%s,1)) THEN
                         SYMPROD%s=IOR(SYMPROD%s,SYMTABLE(I,J)%s)
                      ENDIF
                      IS2%s=RSHIFT(IS2%s,1)
                      J=J+1
                   ENDDO
                ENDIF
                IS1%s=RSHIFT(IS1%s,1)
                I=I+1
             ENDDO
         end if
         RETURN
      END

      FUNCTION SymConj(s2)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         TYPE(Symmetry) s,SymConj,s2
         INTEGER i,AbelConj(3)
         integer*8 ComposeAbelianSym
         if (TAbelian) then
             ! K-point symmetry has k_-i=k_i.  We store k-vectors from 0
             ! to N rather than -N/2 to N.  Hence k_-i=mod(-k_i+N,N) for
             ! each component of the vector.
             call DecomposeAbelianSym(s2%s,AbelConj)
             do i=1,3
               AbelConj(i)=modulo(-AbelConj(i)+NProp(i),NProp(i))
             end do
             SymConj%s=ComposeAbelianSym(AbelConj)
         else
             SymConj%s=0
             s=s2
             i=1
             DO WHILE(s%s.NE.0)
                IF(IAND(s%s,1)) SymConj%s=IOR(SymConj%s,
     &         LSHIFT(1,SymConjTab(I)-1))
                s%s=RSHIFT(s%s,1)
                i=i+1
             ENDDO
         end if
      END

      SUBROUTINE WRITESYMTABLE(IUNIT)
         IMPLICIT NONE
         INTEGER IUNIT,I,J
         INCLUDE 'sym.inc'
         
         DO I=0,NSYM-1
            DO J=0,NSYM-1
               CALL WRITESYM(IUNIT,SYMTABLE(I+1,J+1),.FALSE.)
            ENDDO
            WRITE(IUNIT,*)
         ENDDO
         IF(NSYM.EQ.0) THEN
            WRITE(6,*) "No Symmetry table found."
         ENDIF 
      END

      LOGICAL FUNCTION LSYMSYM(SYM)
         include 'sym.inc'
         Type(Symmetry) SYM
         if (TAbelian) then
             LSymSym=Sym%s.eq.0
         else
             LSYMSYM=(SYM%s.EQ.0.OR.BTEST(SYM%s,0))
         end if
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
         TYPE(BasisFN) G1(*)
         INTEGER NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS)
         TYPE(Symmetry) SYMT2(NSYMMAX,NSYMMAX)
         POINTER (IP_SYMT2,SYMT2)
        
         WRITE(6,"(A,I3,A)") "Generating abelian symmetry table with",  &
     &      INT(DLOG(NSYMMAX+0.D0)/DLOG(2.D0)+.4), " generators" 

C.. Setup the symmetry product table
         NSYM=NSYMMAX
         CALL MEMORY(IP_SYMTABLE,NSYM*NSYM*SymmetrySize,'SYMTABLE')
         IP_SYMT2=IP_SYMTABLE
         CALL MEMORY(IP_SymConjTab,NSYM,'SymConjTab')
         CALL IAZZERO(SYMTABLE,NSYM*NSYM*SymmetrySize)
         DO I=0,NSYM-1
            SymConjTab(I+1)=I+1
            DO J=0,NSYM-1
               SYMT2(I+1,J+1)%s=IBSET(0,IEOR(I,J))
            ENDDO
         ENDDO
C.. Now generate a list of sym labels.
         CALL MEMORY(IP_SYMLABELS,NSYM*SymmetrySize,"SYMLABELS")
         CALL MEMORY(IP_SymClasses,NBASIS,"SymClasses")
         NSYMLABELS=NSYM
         DO I=1,NBASIS,2
C.. place the sym label of each state in SymClasses(ISTATE).  For molp sym, this is 
C.. the log_2 of the symmetry bit string
            IF(G1(I)%Sym%s.EQ.0) THEN
C.. we don't have symmetry, so fake it.
               SymClasses((I+1)/2)=1
            ELSE
          SymClasses((I+1)/2)=INT(DLOG(G1(I)%Sym%s+0.D0)/DLOG(2.D0)+1.4)
            ENDIF
         ENDDO
C.. list the symmetry string of each sym label
         DO I=1,NSYM
            SYMLABELS(I)%s=2**(I-1)
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
C.. Each state's symmetry falls into a class SymClasses(ISTATE).
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SymLabels(ISYMLABEL)
C.. The characters of this class are stored in SYMLABELCHARS(1:NROT, SymClasses(ISTATE))
C.. The total number of symmetry labels is NSYMLABELS
         DO I=1,NHG,2
            IF(GG(I).NE.0) THEN
               NSL((GG(I)+1)/2)=SymClasses((I+1)/2)
            ENDIF
         ENDDO
         DO I=1,NBASIS/2
            SymClasses(I)=NSL(I)
C            WRITE(6,*) "SL",I,SymClasses(I)
         ENDDO
      END

      SUBROUTINE GENMOLPSYMREPS(NSYMMAX,G1,NBASIS,ARR,BRR)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER NSYMMAX
         INTEGER I,J
         TYPE(BasisFN) G1(NBASIS)
         INTEGER NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS)
         INTEGER SYMT2(NSYMMAX,NSYMMAX)
         POINTER (IP_SYMT2,SYMT2)
         
C.. now work out which reps are degenerate and label them
         CALL MEMORY(IP_SYMREPS,2*NBASIS,'SYMREPS')
         J=0
         DO I=1,NBASIS
            IF(I.GT.1.AND.(ABS(ARR(I)-ARR(I-1)).LT.1.D-5
     &         .AND.G1(BRR(I))%Sym%s.EQ.G1(BRR(I-1))%Sym%s)) THEN
C.. we have the same degenerate rep as the previous entry
               SYMREPS(2,J)=SYMREPS(2,J)+1
            ELSE
C.. we have a new rep
               J=J+1
               SYMREPS(2,J)=1
            ENDIF
            SYMREPS(1,BRR(I))=J
         ENDDO
      END

C.. delete a symmetry table if one existed.
      SUBROUTINE ENDSYM()
         INCLUDE 'sym.inc'
         IF(IP_SYMTABLE.NE.0) CALL FREEM(IP_SYMTABLE)
         IP_SYMTABLE=0
         IF(IP_SymConjTab.NE.0) CALL FREEM(IP_SymConjTab)
         IP_SymConjTab=0
         IF(IP_SYMREPS.NE.0) CALL FREEM(IP_SYMREPS)
         IP_SYMREPS=0
         IF(IP_SymClasses.NE.0) CALL FREEM(IP_SymClasses)
         IP_SymClasses=0
         IF(IP_SYMLABELS.NE.0) CALL FREEM(IP_SYMLABELS)
         IP_SYMLABELS=0
         IF(IP_SYMLABELCHARS.NE.0) CALL FREEM(IP_SYMLABELCHARS)
         IP_SYMLABELCHARS=0
         IF(IP_IRREPCHARS.NE.0) CALL FREEM(IP_IRREPCHARS)
         IP_IRREPCHARS=0
         IF(IP_SymStatePairs.NE.0) CALL FREEM(IP_SymStatePairs)
         IP_SymStatePairs=0
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
      SUBROUTINE GENSymStatePairs(NSTATES)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,IDECOMP,TOT,NPRODS,NSTATES
         INTEGER TEMPLIST(NSTATES)
         CALL MEMORY(IP_SYMLABELLIST,NSTATES,'SYMLABELLIST')
         CALL MEMORY(IP_SYMLABELCOUNTS,2*NSYMLABELS,'SYMLABELCOUNTS')
C.. First deal with listing single states
         DO I=1,NSTATES
            SYMLABELLIST(I)=I
            TEMPLIST(I)=SymClasses(I)
         ENDDO
C.. order according to sym label, so SYMLABELLIST gets a list of states grouped under SYMLABEL
         CALL SORT2I(NSTATES,TEMPLIST,SYMLABELLIST)
         CALL IAZZERO(SYMLABELCOUNTS,2*NSYMLABELS)
         SYMLABELCOUNTS(1,TEMPLIST(1))=1
         SYMLABELCOUNTS(2,TEMPLIST(1))=1
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
!         DO I=1,NSYMLABELS
!            WRITE(6,*) "NSL",I,SYMLABELCOUNTS(1,I),SYMLABELCOUNTS(2,I),
!     &         SymLabels(I)
!         ENDDO

C.. Now deal with pairs of states


         CALL MEMORY(IP_SymPairProds,
     &         SymPairProdSize*nSymLabels*nSymLabels,
     &      "SymPairProds")
         CALL IAZZERO(SymPairProds,
     &         SymPairProdSize*nSymLabels*nSymLabels)
C.. Now enumerate all pairs, and classify their product, but don't store them.
         nSymPairProds=0
         CALL GENALLSymStatePairs(NSTATES,.FALSE.)

C.. Now sort the SymPairProds into order
         CALL SORTSYMARR(nSymPairProds,SymPairProds,SymPairProdSize)
         TOT=0
         DO I=1,nSymPairProds
            SymPairProds(I)%nIndex=TOT
            TOT=TOT+SymPairProds(I)%nPairs
            SymPairProds(I)%nPairs=0
         ENDDO
         WRITE(6,*) TOT," STATE PAIRS"
         WRITE(6,*) NSYMPAIRPRODS," DISTINCT ORBITAL PAIR PRODUCT SYMS"
         CALL MEMORY(IP_SymStatePairs,
     &         TOT*2,
     &      "SymStatePairs")
         CALL IAZZERO(SymStatePairs,TOT*2)
         CALL GENALLSymStatePairs(NSTATES,.TRUE.)
      END

      SUBROUTINE GENALLSymStatePairs(NSTATES,TSTORE)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         LOGICAL TSTORE
         INTEGER I,J
         TYPE(Symmetry) PROD
         TYPE(Symmetry) SYMPROD
         INTEGER NSTATES,iProd
         DO I=1,NSTATES
CSYMLABELS
            DO J=I,NSTATES
CSYMLABELS
C               WRITE(6,*) I,J,SymClasses(I),SymClasses(J)
               PROD=SYMPROD(SymLabels(SymClasses(I)),
     &               SymLabels(SymClasses(J)))
               CALL FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
               IF(TSTORE) THEN
C.. put the pair into the list of pairs.
                  SymStatePairs(1,SymPairProds(iProd)%nIndex
     &                            +SymPairProds(iProd)%nPairs) =I
                  SymStatePairs(2,SymPairProds(iProd)%nIndex
     &                            +SymPairProds(iProd)%nPairs) =J
               ENDIF
               IF(iProd.EQ.nSymPairProds+1) THEN
                  nSymPairProds=nSymPairProds+1
                  SymPairProds(iProd)%Sym=Prod
                  SymPairProds(iProd)%nIndex=0
                  SymPairProds(iProd)%nPairs=0
               ENDIF
C.. incrememnt the counter in the pairlist
               SymPairProds(iProd)%nPairs=SymPairProds(iProd)%nPairs+1
            ENDDO
         ENDDO
      END

      SUBROUTINE FindSymProd(Prod,SymPairProds,nSymPairProds,iProd)
         INCLUDE 'symtypes.inc'
         TYPE(SymPairProd) SymPairProds(nSymPairProds)
         INTEGER nSymPairProds,iProd
         TYPE(Symmetry) Prod
         LOGICAL SYMEQ
         DO iProd=1,nSymPairProds
            IF(SYMEQ(SymPairProds(iProd)%Sym,Prod)) EXIT
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
         INCLUDE 'sym.inc'
         INTEGER NEL,NI(NEL),NBASISMAX(5,3)
         TYPE(BasisFn) G1(*),ISym
         INTEGER I,J,NI2(NEL)
         INTEGER NREPS(NEL),NELECS(NEL),SSYM
         TYPE(Symmetry) SYMPROD
         LOGICAL ISCSF,ISC
         I=1
         CALL IAZZERO(NREPS,NEL)
         CALL SETUPSYM(NBASISMAX,ISYM)
         ISC=ISCSF(NI2,NEL)
         IF(ISC) THEN
            DO I=1,NEL
               CALL GETUNCSFELEC(NI2(I),NI(I),SSYM)
            ENDDO
         ELSE
            CALL ICOPY(NEL,NI2,1,NI,1)
         ENDIF
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
               ISYM%Sym=SYMPROD(ISYM%Sym,G1(NI(I))%Sym)
            ENDIF
C.. add the momentum
            CALL ADDELECSYM(NI(I),G1,NBASISMAX,ISYM)
         ENDDO
C.. round the momentum
         CALL ROUNDSYM(ISYM,NBASISMAX)
         IF(ISC) CALL CSFGETSPIN(NI2,NEL,ISYM%Ms) 
         RETURN
      END

C.. Given a set of characters of states, generate all relevant irreps which span the set of characters.
      SUBROUTINE GENIRREPS(TKP,IMPROPER_OP,NROTOP)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J,K
         LOGICAL LDO,LDO2
         TYPE(Symmetry) iDecomp
         INTEGER NEXTSYMLAB
         LOGICAL GETIRREPDECOMP
         COMPLEX*16 REPCHARS(NROT,NSYMLABELS*10)
         INTEGER NREPS,NROTOP
         REAL*8 NORM
         LOGICAL TKP,INV,IMPROPER_OP(NROTOP)
         NREPS=0
C.. Initialize the table with the totally symmetric rep.
         INV=.FALSE.
         DO I=1,NROT
            IRREPCHARS(I,1)=1
            IF(IMPROPER_OP(MOD(I-1,NROTOP)+1).and..not.TKP) INV=.TRUE.
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
            NREPS=0
         lp1:DO I=1,NSYM
               DO J=I,NSYM
                  NREPS=NREPS+1
                  IF(NREPS.GT.NSYMLABELS*10) STOP 'TOO MANY REPS'
                  DO K=1,NROT
                     REPCHARS(K,NREPS)=DCONJG(IRREPCHARS(K,I))
     &                                    *IRREPCHARS(K,J)
                  ENDDO
                  
!                  WRITE(6,*) NREPS,"PROD",I,J
                  CALL MEMORY_CHECK
!                  CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDPRD")
                  IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM,TAbelian)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                     IF(ABS(NORM-NROT).LE.1.D-2) THEN
C.. if it's an irrep
                        NSYM=NSYM+1
                        IF(NSYM.GT.64) STOP "MORE than 64 irreps"
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
!                        WRITE(6,*) "NR",NREPS,LDO2
                        IF(LDO2) NREPS=NREPS-1
!                        NREPS=NREPS-1
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
               IF(NREPS.GT.NSYMLABELS*10) STOP 'TOO MANY REPS'
               DO I=1,NROT
                  REPCHARS(I,NREPS)=SYMLABELCHARS(I,NEXTSYMLAB)
               ENDDO
C               CALL WRITECHARS(6,REPCHARS(1,NREPS),NROT,"ADDST ")
               IF(GETIRREPDECOMP(REPCHARS(1,NREPS),IRREPCHARS,
     &               NSYM,NROT,IDECOMP,NORM,TAbelian)) THEN
C.. CHARWORK now contains the remainder, which will be a new irrep (or combination or irreps), which we need to add
                  IF(ABS(NORM-NROT).LE.1.D-2) THEN
C.. if it's an irrep
                     NSYM=NSYM+1
                     IF(NSYM.GT.64) STOP "MORE than 64 irreps"
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
            SymLabels(I)=IDECOMP
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
                  IF(ABS(DIMAG(CHARS(J))).GT.1.D-2) THEN
C.. write in terms of I.
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G14.9,A,$)") CHARS(J)," "
                     ELSE                        
                        IF(ABS(DIMAG(CHARS(J))+1.D0).LT.1.D-2) THEN
                           WRITE(IUNIT,"(A,$)") " -I "
                        ELSEIF(ABS(DIMAG(CHARS(J))-1.D0).LT.1.D-2) THEN
                           WRITE(IUNIT,"(A,$)") "  I "
                        ELSE 
                         WRITE(IUNIT,"(I2,A,$)") NINT(DIMAG(CHARS(J))),
     &                     "I "
                        ENDIF
                     ENDIF
                  ELSE
                     IF(LREAL) THEN
                        WRITE(IUNIT,"(G21.9,$)") REAL(CHARS(J)), "    "
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
         INCLUDE 'sym.inc'
         TYPE(Symmetry) IDECOMP
         COMPLEX*16 CHARS(NROT),CHARSIN(NROT),TOT
         REAL*8 CNORM
         INTEGER I,J
         REAL*8 NORM,DIFF
         if (TAbelian) then
             ! We shouldn't be here!  Using symmetry "quantum" numbers
             ! rather than irreps.
             stop "Should not be decomposing irreps with Abelian "//
     &             "symmetry option on."
         end if
         IDECOMP%s=0
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
!            CALL WRITECHARS(6,IRREPCHARS(1,I),NROT,"IR")
!            CALL WRITECHARS(6,CHARS(1),NROT,"CH")
            DO J=1,NROT
               TOT=TOT+DCONJG(IRREPCHARS(J,I))*CHARS(J)
            ENDDO
!            WRITE(6,*) I,TOT
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
                  IDECOMP%s=IBSET(IDECOMP%s,I-1)
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
     &         IDECOMP,CNORM,TAbelian)
         IMPLICIT NONE
         INTEGER NIRREPS, NROT
         INCLUDE 'basis.inc'
         TYPE(Symmetry) IDECOMP
         COMPLEX*16 IRREPCHARS(NROT,NIRREPS),CHARS(NROT)
         REAL*8 CNORM, NORM,DIFF
         COMPLEX*16 TOT
         INTEGER I,J
         logical TAbelian
         if (TAbelian) then
             ! We shouldn't be here!  Using symmetry "quantum" numbers
             ! rather than irreps.
             stop "Should not be decomposing irreps with Abelian "//
     &             "symmetry option on."
         end if
         IDECOMP%s=0
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
                  IDECOMP%s=IBSET(IDECOMP%s,I-1)
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
         TYPE(Symmetry) IDECOMP
         REAL*8 CNORM
         CALL MEMORY(IP_SYMTABLE,NSYM*NSYM*SymmetrySize,'SYMTABLE')
         CALL MEMORY(IP_SymConjTab,nSym,'SymConjTab')
         DO I=1,NSYM
            DO K=1,NROT
               CHARS(K)=DCONJG(IRREPCHARS(K,I))
            ENDDO
            IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,
     &         IDECOMP,CNORM,TAbelian)) THEN
               WRITE(6,*) "Conjugate of SYM ",I,
     &            " not reducible,"
               CALL WRITECHARS(6,CHARS,NROT,"REMAIN")
               STOP "Symmetry table element not conjugable"
            ENDIF
            K=0
            DO WHILE(.NOT.BTEST(IDECOMP%s,0))
               K=K+1
               IDECOMP%s=RSHIFT(IDECOMP%s,1)
            ENDDO
            IF(IDECOMP%s.NE.1) THEN
               WRITE(6,*) "Conjugate of SYM ",I,
     &            " not a single SYM,"
               STOP
            ENDIF
            SymConjTab(I)=K+1
            DO J=I,NSYM
               DO K=1,NROT
                  CHARS(K)=IRREPCHARS(K,I)*IRREPCHARS(K,J)
               ENDDO
               IF(GETIRREPDECOMP(CHARS,IRREPCHARS,NSYM,NROT,
     &            IDECOMP,CNORM,TAbelian)) THEN
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
         WRITE(6,*) "Symmetry, Symmetry Conjugate"
         DO I=1,NSYM
            WRITE(6,*) I,SymConjTab(I)
         ENDDO
      END

      SUBROUTINE GENSYMREPS(G1,NBASIS,ARR,BRR,DEGENTOL)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER I,J
         TYPE(BasisFN) G1(nBasis)
         INTEGER NBASIS,BRR(NBASIS)
         REAL*8 ARR(NBASIS,2)
         REAL*8 DEGENTOL
         LOGICAL SYMEQ

C.. now work out which reps are degenerate and label them
         CALL MEMORY(IP_SYMREPS,2*NBASIS,'SYMREPS')
         J=0
         DO I=1,NBASIS
C            WRITE(6,*) "SR2",I
            IF(I.GT.1.AND.ABS(ARR(I,2)-ARR(I-1,2)).LT.DEGENTOL
     &         .AND.SYMEQ(G1(I)%Sym,G1(I-1))) THEN
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

C.  Irrep symmetries are specified in SYM(5).
C.. if SYM(5)=0, we assume it's totally symmetric
C.. Other irreps contributing to the symmetry have bits set in 
C.. SYM.
C.. e.g. if irreps are a1,a2,b1,b2
      LOGICAL FUNCTION LCHKSYM(ISYM,JSYM)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BASISFN) ISYM,JSYM
         INTEGER I,IS,JS
         TYPE(SYMMETRY) SYMPROD,SymConj
         LOGICAL LSYMSYM
         LCHKSYM=.TRUE.
         DO I=1,4
            IF(ISYM%K(I).NE.JSYM%K(I)) LCHKSYM=.FALSE.
         ENDDO
C.. if the symmetry product of I and J doesn't contain the totally
C.. symmetric irrep, we set sym to .FALSE.
        LCHKSYM=LCHKSYM.AND.LSYMSYM(SYMPROD(SymConj(ISYM%SYM),JSYM%SYM))
      RETURN
      END
      
      LOGICAL FUNCTION LCHKSYMD(NI,NJ,NEL,G1,NBASISMAX)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BASISFN) ISYM,JSYM,G1(*)
         INTEGER NEL,NI(NEL),NJ(NEL),NBASISMAX(5,3)
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
         INCLUDE 'basis.inc'
         TYPE(BASISFN) ISYM,G1(*)
         INTEGER IELEC,NBASISMAX(5,3)
         INTEGER I,IEL,SSYM
         CALL GETUNCSFELEC(IEL,IELEC,SSYM)
        IF(NBASISMAX(1,3).LT.4) THEN
C.. Momentum space
            DO I=1,3
               ISYM%K(I)=ISYM%K(I)+G1(IELEC)%K(I)
            ENDDO
C.. Symmetry space
         ELSEIF(NBASISMAX(3,3).EQ.0.AND.NBASISMAX(1,3).GE.4) THEN
C.. We have no symmetries, so do nothing. (we're in real space)
C.. except Ms
         ELSEIF(NBASISMAX(3,3).EQ.1) THEN
C.. deal with momentum
            DO I=1,3
               ISYM%k(I)=ISYM%k(I)+G1(IELEC)%k(I)
            ENDDO
         ENDIF
         ISYM%MS=ISYM%MS+G1(IELEC)%MS
C.. SSYM keeps track of the total S change on adding this electron
C.. (it is +/-CSF_NSBASIS)
         I=ISYM%MS+0
         ISYM%Ms=I+SSYM
         RETURN
      END 
      
      SUBROUTINE ROUNDSYM(ISYM,NBASISMAX)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BasisFN) ISYM
         INTEGER NBASISMAX(5,3)
         INTEGER I
         IF(NBASISMAX(3,3).EQ.-2) THEN
C.. particle in a box
C.. parity symmetries
            DO I=1,3
               ISYM%k(I)=MOD(ISYM%k(I),2)
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
                  ISYM%k(I)=0
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

      LOGICAL FUNCTION SYMLT(A,B)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(Symmetry) A,B
         IF(A%s.GE.0) THEN
            IF(B%s.GE.0) THEN
               SYMLT=A%s.LT.B%s
            ELSE
               SYMLT=.TRUE.
            ENDIF
         ELSE
            IF(B%s.GE.0) THEN
               SYMLT=.FALSE.
            ELSE
               SYMLT=A%s.LT.B%s
            ENDIF
         ENDIF    
         RETURN
      END
      LOGICAL FUNCTION SYMNE(A,B)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(Symmetry) A,B
         SYMNE=A%s.NE.B%s
         RETURN
      END
      LOGICAL FUNCTION SYMEQ(A,B)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(Symmetry) A,B
         SYMEQ=A%s.EQ.B%s
         RETURN
!Need to cope with 'unsigned integers'
      END
      LOGICAL FUNCTION SYMGT(A,B)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(Symmetry) A,B
         IF(A%s.GE.0) THEN
            IF(B%s.GE.0) THEN
               SYMGT=A%s.GT.B%s
            ELSE
               SYMGT=.FALSE.
            ENDIF
         ELSE
            IF(B%s.GE.0) THEN
               SYMGT=.TRUE.
            ELSE
               SYMGT=A%s.GT.B%s
            ENDIF
         ENDIF    
         RETURN
      END
C.. A binary search to find VAL in TAB.  TAB is sorted

      SUBROUTINE BINARYSEARCHSYM(VAL,TAB,ROWLEN,LEN,LOC)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(Symmetry) VAL
         INTEGER LOC,LEN,ROWLEN
         INTEGER TAB(ROWLEN,LEN)
         INTEGER I,J,IFIRST,N,ILAST
         LOGICAL SYMNE,SYMLT,SYMGT,SYMEQ
         I=1
         J=LEN
         IFIRST=I
         ILAST=J
         DO WHILE(J-I.GE.1)
            N=(I+J)/2
C            WRITE(6,"(3I4,$)") I,J,N
C            CALL WRITESYM(6,TAB(1,I),.FALSE.)
C            CALL WRITESYM(6,TAB(1,J),.FALSE.)
C            CALL WRITESYM(6,TAB(1,N),.TRUE.)
            IF(SYMLT(TAB(1,N),VAL).AND.I.NE.N) THEN
               IF(SYMNE(TAB(1,N),TAB(1,IFIRST))) IFIRST=N
C.. reset the lower limit
               I=N
            ELSEIF(SYMGT(TAB(1,N),VAL)) THEN
               IF(SYMNE(TAB(1,N),TAB(ILAST))) ILAST=N
C.. reset the upper limit
               J=N
            ELSEIF(SYMEQ(TAB(1,N),VAL)) THEN
C.. bingo, we've got it!
               LOC=N
               RETURN
            ELSE
C.. we've reached a situation where I and J's entries have the same value, and it's
C.. not the one we want.  Leave the loop.
               I=J
            ENDIF
         ENDDO
         IF(SYMEQ(TAB(1,I),VAL)) THEN
            LOC=I
         ELSEIF(SYMEQ(TAB(1,J),VAL)) THEN
            LOC=J
         ELSE
C.. Failure
            LOC=0
         ENDIF
      END
      SUBROUTINE GENNEXTSYM(NEL,NBASISMAX,TSPN,LMS,
     &      TPARITY,IPARITY,TSETUP,TDONE,IMAX,ISYM)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER NEL,NBASISMAX(5,3)
         INTEGER LMS
         TYPE(BasisFN) IPARITY,ISYM,IMax(2)
         LOGICAL TSPN,TPARITY,TSETUP,TMORE,TDONE,KALLOWED,TMORE2
         INTEGER ILEV
         integer*8 TotSymRep
         IF(TSETUP) THEN
            DO ILEV=1,3
               IF(TPARITY) THEN
                  IMAX(1)%k(iLev)=IPARITY%k(ILEV)
                  IMAX(2)%k(ILEV)=IPARITY%k(ILEV)
               ELSE
                  IMAX(1)%k(iLev)=NBASISMAX(ILEV,1)
                  IMAX(2)%k(iLev)=NBASISMAX(ILEV,2)
C                  IF(NBASISMAX(1,3).EQ.2) THEN
C.. hubbard non-pbc mom space
C                     IMAX(ILEV,1)=IMAX(ILEV,1)*NEL
C                     IMAX(ILEV,2)=IMAX(ILEV,2)*NEL
C                  ENDIF
               ENDIF
            ENDDO
            IF(TSPN) THEN
               IMAX(1)%Ms=LMS
               IMAX(2)%Ms=LMS
            ELSE
               IMAX(1)%Ms=NBASISMAX(4,1)*NEL
               IMAX(2)%Ms=NBASISMAX(4,2)*NEL
            ENDIF
C.. If we're specifying a sym (TPARITY) in IPARITY(5), and
C.. we have a system with all 1D reducible orbs, then we put
C.. that into IMAX
            IF(NBASISMAX(5,2).NE.0) THEN
               IF(TPARITY) THEN
                  IMAX(1)%Sym%s=IPARITY%Sym%s
                  IMAX(2)%Sym%s=IMAX(1)%Sym%s
               ELSE
                  IMAX(1)%Sym%s=TotSymRep()
                  IMAX(2)%Sym%s=2**NBASISMAX(5,2)
               ENDIF
            ELSE
C.. we've got a sym system with polydimensional irreps, which leads to
C.. dets with combinations of irreps, so we cannot put sym into blocks
                
!  JSS: if only 1D symmetries, then a determinant can only interact with
!  other determinants of the same symmetry.  This applies to Abelian
!  groups.  If there are multi-dimensional irreps, then this is no
!  longer the case, so we set the symmetries to be 0 (i.e. ignore
!  symmetry when generating determinants which interact).  This is not
!  equivalent to setting %s=0 if the Abelian case (which corresponds to
!  the totally symmetric irrep).
               IMAX(1)%Sym%s=0
               IMAX(2)%Sym%s=0
            ENDIF
            TDONE=.FALSE.
            CALL DOSYMLIMDEGEN(IMAX,NBASISMAX)
            ISym=IMax(1)
         ENDIF
         IF(TSETUP.AND.KALLOWED(ISYM,NBASISMAX)) RETURN
C.. Go to the next sym.
         TMORE2=.TRUE.
         TMORE=.TRUE.
         ILEV=5
         DO WHILE(TMORE2)
            DO WHILE (ILEV.GT.0)
               IF(ILEV.EQ.5.AND.ISYM%k(ILEV).NE.0) THEN
C.. symmetry specifiers are incremented by multiplying*2 (unless there are no syms counted)
                  ISYM%k(ILEV)=ISYM%k(ILEV)*2
               ELSE
                  ISYM%k(ILEV)=ISYM%k(ILEV)+1
               ENDIF
               IF(ISYM%k(ILEV).GT.IMAX(2)%k(ILEV)) THEN
                  ISYM%k(ILEV)=IMAX(1)%k(ILEV)
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
                  ISYM%k(ILEV)=IMAX(1)%k(ILEV)-1
                  
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
         INCLUDE 'basis.inc'
         TYPE(BasisFN) IMax(2)
         INTEGER NBASISMAX(5,5),I
         IF(NBASISMAX(3,3).EQ.0) THEN
            DO I=1,3
               IF(IMax(2)%k(I).NE.IMAX(1)%k(I)) IMAX(1)%k(I)=0
            ENDDO
         ENDIF
C.. always a spin symmetry
         IF(IMAX(1)%Ms.NE.IMAX(2)%Ms) IMAX(1)%Ms=0
      END
      SUBROUTINE GETSYMDEGEN(ISYM,NBASISMAX,IDEGEN)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BasisFN) ISym,ISym2
         INTEGER NBASISMAX(5,5),IDEGEN,I,J
         LOGICAL KALLOWED,TDO
         IDEGEN=0
         IF(NBASISMAX(3,3).EQ.0) THEN
C.. Hubbard
            DO I=0,7
               TDO=.TRUE.
               DO J=1,3
                  IF(BTEST(I,J-1).EQ.0) THEN
                     ISYM2%k(J)=ISYM%k(J)
                  ELSE
                     ISYM2%k(J)=-ISYM%k(J)
                     IF(ISYM%k(J).EQ.0) TDO=.FALSE.
                  ENDIF
               ENDDO
               IF(TDO.AND.KALLOWED(ISYM2,NBASISMAX)) IDEGEN=IDEGEN+1
            ENDDO
         ELSE
            IDEGEN=1
         ENDIF
C.. Spin
         IF(ISYM%k(4).NE.0) IDEGEN=IDEGEN*2
      END

C.. Initialize symmetry to take into account the core electrons
      SUBROUTINE SETUPSYM(NBASISMAX,ISYM)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER NBASISMAX(5,6)
         TYPE(BasisFN) ISym
         ISym=FrozenSym
         RETURN
      END

      SUBROUTINE WRITEALLSYM(IUNIT,SYM,LTERM)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         INTEGER IUNIT
         TYPE(BASISFN) SYM
         LOGICAL LTERM
         INTEGER J
         WRITE(IUNIT,"(4I5,$)") SYM%K(1),SYM%K(2),SYM%K(3),SYM%MS
         CALL WRITESYM(IUNIT,SYM%SYM,LTERM)
      END
      SUBROUTINE WRITESYM(IUNIT,SYM,LTERM)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INTEGER IUNIT
         TYPE(SYMMETRY) SYM
         LOGICAL LTERM
         INTEGER J
         IF(TAbelian) THEN
            WRITE(IUNIT,"(Z12,$)") SYM
         ELSEIF(NSYM.LE.16) THEN
            WRITE(IUNIT,"(Z5,$)") SYM
         ELSEIF(NSYM.LE.24) THEN
            WRITE(IUNIT,"(Z7,$)") SYM
         ELSEIF(NSYM.LE.32) THEN
            WRITE(IUNIT,"(Z9,$)") SYM
         ELSEIF(NSYM.LE.40) THEN
            WRITE(IUNIT,"(Z11,$)") SYM
         ELSEIF(NSYM.LE.48) THEN
            WRITE(IUNIT,"(Z13,$)") SYM
         ELSEIF(NSYM.LE.56) THEN
            WRITE(IUNIT,"(Z15,$)") SYM
         ELSE
            WRITE(IUNIT,"(Z17,$)") SYM
         ENDIF
         IF(LTERM) WRITE(IUNIT,*)
      END
      SUBROUTINE SetupFreezeAllSym(Sym)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         TYPE(BasisFN) Sym
C.. Set to be totally symmetric
         FrozenSym=Sym
      END
      SUBROUTINE SetupFreezeSym(Sym)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         TYPE(Symmetry) Sym
         CALL IAZZERO(FrozenSym,BasisFNSize)
C.. Set to be totally symmetric
         FrozenSym%Sym=Sym
      END
 
!Deal with K-point symmetries, using translational symmetry operations.
      ! JSS: use Abelian symmetry formulation (allows us to go beyond 64
      ! symmetry operations, and hence deal with larger k-point meshes).
      SUBROUTINE GenKPtIrreps(cKptChars,nTranslat,nKps,KpntInd,nStates)
      IMPLICIT NONE
      COMPLEX*16 cKptChars(nTranslat,nKps)
      INTEGER nTranslat,nKps,KpntInd(nStates)
      INTEGER I,J,nStates
      integer*8  ComposeAbelianSym
      INCLUDE 'sym.inc'
      type(Symmetry) SymConj,SymProd
      nSymLabels=nKps
      nRot=nTranslat
      CALL MEMORY(IP_SYMLABELCHARS,nSymLabels*nRot*2,'SYMLABELCH')
      CALL AZZERO(SYMLABELCHARS,nSymLabels*nRot*2)
      CALL MEMORY(IP_SYMLABELS,nSymLabels,'SYMLABELS')
      CALL MEMORY(IP_SymClasses,nStates,'SymClasses')
!     DO I=1,nKps
!        DO J=1,nRot
!           SymLabelChars(J,I)=cKptChars(J,I)
!        ENDDO
!     ENDDO
      DO I=1,nStates
        SymClasses(I)=KpntInd(I)
        SymLabels(KPntInd(I))%s=ComposeAbelianSym(KpntSym(:,KPntInd(I)))
      END DO
      write (6,*) 
      write(6,'(a11," |",a9," |",a10)')' K-vector','  Label','Conjugate'
      write (6,'(34("-"))')
      do i=1,nSymLabels
        write (6,'("(",3i3,")"," | ",z8," | ",z8)') 
     &                   KpntSym(:,I),SymLabels(I),SymConj(SymLabels(I))
      end do
      write (6,'(/,a)') 'Symmetry Multiplication Table'
      do i=1,nSymLabels
        do j=1,nSymLabels
          write (6,'(z12,$)') SymProd(SymLabels(I),SymLabels(J))
        end do
        write (6,*) 
      end do
      write (6,'(/)') 
!        WRITE(6,*) "SYMMETRY CLASSES"
!        CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
!.. Allocate memory gor irreps.
!.. Assume there will be no more than 64 irreps
!        CALL MEMORY(IP_IRREPCHARS,NROT*64*2,"IRREPCH")
      END

      subroutine  DecomposeAbelianSym(ISym,AbelSym)
      ! Store the symmetry index as integer*8.  For Abelian symmetry
      ! we need to have 3 numbers stored in this index.  We store
      ! according to isym=\sum_i AbelSym(i)*32768**(i-1).
      ! This allows point groups with more than 64 irreps to be used in
      ! the point group is Abelian (as all translational/k-point
      ! symmetries are).
      ! Decompose the symmetry label back into the appropriate "quantum"
      ! numbers.
      ! Store the symmetry index as integer*8.  For Abelian symmetry
      ! we need to have 3 numbers stored in this index.  We store
      ! according to isym=1+\sum_i AbelSym(i)*32768**(i-1).
      ! Decompose the symmetry label back into the appropriate
      ! numbers...
      implicit none
      integer*8 Isym
      integer  AbelSym(3)
      AbelSym(3)=RShift(Isym,30)
      AbelSym(2)=Iand(RShift(ISym,15),32767)
      AbelSym(1)=Iand(Isym,32767)
      return
      end subroutine DecomposeAbelianSym

      integer*8 function ComposeAbelianSym(AbelSym)
          implicit none
          integer  AbelSym(3)
          integer*8 TempVar
          TempVar=AbelSym(3)
          ComposeAbelianSym=LShift(Tempvar,30)+LShift(AbelSym(2),15)
     &                      +AbelSym(1)
      end function ComposeAbelianSym


      integer*8 function TotSymRep
          ! Our definition of the totally symmetric representation
          ! changes according to whether we're using Abelian/k-point
          ! symmetry or the standard symmetry.  It's just a matter of
          ! convenience, rather than some deep theoretical insight!
          include 'sym.inc'
          if (TAbelian) then
              TotSymRep=0
          else
              TotSymRep=1
          end if
      end function TotSymRep
