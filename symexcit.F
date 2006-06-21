
C.. Setup for symmetry routine below.
      SUBROUTINE SYMSETUPEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,TCOUNT,
     &      ICOUNT)
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL)
         INTEGER G1(5,NBASIS),NBASISMAX(*),NBASIS
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER STORE(6),ICOUNT
         LOGICAL TCOUNT
         CALL IAZZERO(STORE,6)
         CALL MEMORY(IP_DSTORE,NEL+NBASIS/32+1,'DSTORE')
         STORE(1)=IP_DSTORE
         CALL SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,TCOUNT,
     &      ICOUNT,DSTORE(1),
     &      DSTORE(NEL+1))
C.. If we're just counting, we don't need to keep DSTORE
         IF(TCOUNT) CALL FREEM(IP_DSTORE)
      END


C.. A test subroutine to trial a new symmetry excitation generation algorithm
C.. Algorithm: First list all the different symmetry classes in NI.
C.. For each pair of symmetry classes, determine its symmetry product.
C.. Classify and store each possible pair of orbitals under its symmetry product.
C.. With each of the symmetry products, calculate []' such that []x[]'
C.. contains the totally symmetric rep.
C.. We do this by checking whether any of the []' in the 
C.. global symprods table multiply by [] to give the sym rep.

      SUBROUTINE SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &   TCOUNT,ICOUNT,CLASSES,ILUT)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL)
         INTEGER G1(5,NBASIS),NBASISMAX(*),NBASIS
         LOGICAL TCOUNT
         INTEGER STORE(*)

         INTEGER CLASSES(NEL),CLASSCOUNT(2,NEL)
         INTEGER SYMPRODS(0:NEL*NEL),NPR,PR
         INTEGER SYMPRODCOUNT(3,0:NEL*NEL)
         INTEGER NPAIRS,ICC1,ICC2,ICC3,ICC,ICC4
         INTEGER I,J,NCL,K,L
         LOGICAL TNEW
         INTEGER ORBPAIRS(2,*)
         POINTER (IP_ORBPAIRS,ORBPAIRS)
         INTEGER SYMPRODIND(2,3,0:*)
         POINTER (IP_SYMPRODIND,SYMPRODIND)
         INTEGER NALLOWPPS(3,*),SPP
         POINTER(IP_NALLOWPPS,NALLOWPPS)
         LOGICAL L1B,L1A,L2B,L2A
         INTEGER NEXCITS,NSING,NDOUB,NEXCITTYPES

         INTEGER EXCITTYPES(5,*)
         POINTER (IP_EXCITTYPES,EXCITTYPES)

         INTEGER ICOUNT
         INTEGER SYMPROD
         LOGICAL LSYMSYM
         INTEGER ISPN
         INTEGER ILUT(0:NBASIS/32)
         INTEGER ISUB
         
         CALL TISET('SYMSUEXCIT',ISUB)


C.. List the symmetry classes in NI.
C.. For each orb, just see if its class is in the list, and if not add it
C.. Each state's symmetry falls into a class SYMLABELS(1,ISTATE).
         NCL=0
         DO I=1,NEL
            TNEW=.TRUE.
            DO J=1,NCL
               IF(SYMLABELS(1,(NI(I)+1)/2).EQ.CLASSES(J)) THEN
                  ISPN=(G1(4,NI(I))+3)/2
                  CLASSCOUNT(ISPN,J)=CLASSCOUNT(ISPN,J)+1
                  TNEW=.FALSE.
                  EXIT
               ENDIF
            ENDDO
            IF(TNEW) THEN
C.. add the new class
C.. ISPN=1 is beta, ISPN=2 is alpha
               ISPN=(G1(4,NI(I))+3)/2
               NCL=NCL+1
               CLASSES(NCL)=SYMLABELS(1,(NI(I)+1)/2)
               CLASSCOUNT(ISPN,NCL)=1
               CLASSCOUNT(3-ISPN,NCL)=0
C               WRITE(6,"(A,2I4,Z4,2I4)") "CLASS ",NCL,CLASSES(NCL),
C     &            SYMLABELS(2,CLASSES(NCL)),I
            ENDIF
         ENDDO
C         WRITE(6,*) NCL," Symmetry Classes"
C.. For each pair of symmetry classes, determine its symmetry product.
C.. We list them and then sort them
         CALL IAZZERO(SYMPRODCOUNT,(NEL*NEL+1)*3)
         NPR=0
         NPAIRS=0
         DO I=1,NCL
            DO J=I,NCL
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SYMLABELS(2,ISYMLABEL)
               PR=SYMPROD( SYMLABELS(2,CLASSES(I)),
     &                     SYMLABELS(2,CLASSES(J)))
               TNEW=.TRUE.
               DO K=1,NPR
                  IF(SYMPRODS(K).EQ.PR) THEN
                     TNEW=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF(TNEW) THEN
                  NPR=NPR+1
                  K=NPR
                  SYMPRODS(NPR)=PR
               ENDIF
C.. 1 corresponds to beta beta, 2 to alpha/beta and beta/alpha, and 3 to alpha, alpha
               IF(I.EQ.J) THEN
                  ICC1=CLASSCOUNT(1,I)*(CLASSCOUNT(1,I)-1)/2
                  ICC2=CLASSCOUNT(1,I)*CLASSCOUNT(2,I)
                  ICC3=CLASSCOUNT(2,I)*(CLASSCOUNT(2,I)-1)/2
               ELSE
                  ICC1=CLASSCOUNT(1,I)*CLASSCOUNT(1,J)
                  ICC2=CLASSCOUNT(1,I)*CLASSCOUNT(2,J)
     &                +CLASSCOUNT(2,I)*CLASSCOUNT(1,J)
                  ICC3=CLASSCOUNT(2,I)*CLASSCOUNT(2,J)
               ENDIF
               SYMPRODCOUNT(1,K)=SYMPRODCOUNT(1,K)+ICC1
               SYMPRODCOUNT(2,K)=SYMPRODCOUNT(2,K)+ICC2
               SYMPRODCOUNT(3,K)=SYMPRODCOUNT(3,K)+ICC3
     
               NPAIRS=NPAIRS+ICC1+ICC2+ICC3
C                  WRITE(6,"(A,I6,Z4,5I6)") "SYMPROD ",K,SYMPRODS(K),I,J,
C     &               ICC1,ICC2,ICC3
            ENDDO
         ENDDO
         SYMPRODS(0)=0
         SYMPRODCOUNT(1,0)=0
         SYMPRODCOUNT(2,0)=0
         SYMPRODCOUNT(3,0)=0
C         DO I=1,NPR
C            WRITE(6,"(Z4,3I4)") SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
C         ENDDO
         CALL SORTIIARR(NPR+1,SYMPRODS,SYMPRODCOUNT,3)
C         WRITE(6,*) NPR," Symmetry Products"
C         WRITE(6,*) NPAIRS," Orbital Pairs"
C.. Create a cumulative sum.
         DO I=1,NPR
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)+SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)+SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)+SYMPRODCOUNT(3,I)
C            WRITE(6,*) SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
         ENDDO


C.. Allocate enough memory to store the index
         IF(STORE(5).EQ.0) THEN
            CALL MEMORY(IP_SYMPRODIND,2*3*(NPR+1),'SYMPRODIND')   
         ELSE
C.. we've acually got space to store it somewhere
            IP_SYMPRODIND=STORE(5)
         ENDIF
         CALL IAZZERO(SYMPRODIND,2*3*(NPR+1))

C.. Now shift this such that SYMPRODCOUNT(ISPN,I) is the index of
C.. the first orbital pair of SYMPROD(I) with spin ISPN in ORBPAIRS
C.. Store in SYMPRODIND(1,ISPN,I) too
         DO I=NPR,1,-1
            SYMPRODCOUNT(3,I)=SYMPRODCOUNT(2,I)
            SYMPRODCOUNT(2,I)=SYMPRODCOUNT(1,I)
            SYMPRODCOUNT(1,I)=SYMPRODCOUNT(3,I-1)
            SYMPRODIND(1,3,I)=SYMPRODCOUNT(3,I)
            SYMPRODIND(1,2,I)=SYMPRODCOUNT(2,I)
            SYMPRODIND(1,1,I)=SYMPRODCOUNT(1,I)
         ENDDO
C         DO I=1,NPR
C            WRITE(6,"(Z4,3I4)") SYMPRODS(I),SYMPRODCOUNT(1,I),
C     &         SYMPRODCOUNT(2,I),SYMPRODCOUNT(3,I)
C         ENDDO
C.. We allocate enough memory to store all the pairs.
C.. Each pair consists of (ORB1,ORB2) where ORB1<ORB2
         IF(STORE(4).EQ.0) THEN
            CALL MEMORY(IP_ORBPAIRS,2*NPAIRS/IRAT+1,'ORBPAIRS')
         ELSE
            IP_ORBPAIRS=STORE(4)
         ENDIF

C.. Classify and store each possible pair of orbitals under its symmetry product.
C.. With each symmetry product (denoted []), calculate []'' which is the 
C.. remainder after having removed the orbitals.


C.. Now classify each pair of orbs under its sym prod
C.. SYMPRODS(I) specifies a sym prod.  Orbitals for this are stored in
C.. ORBPAIRS starting at index SYMPRODIND(1,I-1)+1.
C.. The first free space is at SYMPRODCOUNT(I-1)+1.
         CALL IAZZERO(ORBPAIRS,2*NPAIRS)
         L=0
         DO I=1,NEL
            
            DO J=I+1,NEL
               L=L+1
               PR=SYMPROD(G1(5,NI(I)),G1(5,NI(J)))
C.. Find the product in the sorted symprod list
               CALL BINARYSEARCH(PR,SYMPRODS(1),1,NPR+1,K,ICC1,ICC2)
               ISPN=(G1(4,NI(I))+G1(4,NI(J)))/2+2

C.. Increment the offset for the next one
               SYMPRODIND(2,ISPN,K)=SYMPRODIND(2,ISPN,K)+1
               ICC=SYMPRODIND(1,ISPN,K)+SYMPRODIND(2,ISPN,K)
C..                 The index          +  offset of the element

               ORBPAIRS(1,ICC)=NI(I)
               ORBPAIRS(2,ICC)=NI(J)
               
C               WRITE(6,"(A,2I4,2Z4,2I4)") 
C     &            "SP",ICC,ISPN,PR,SYMPRODS(K),NI(I),NI(J)
C               WRITE(6,*) L,NI(I),NI(J),PR,K
            ENDDO
         ENDDO

C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements

C.. Now classify each pair of orbs under its sym prod
C.. SYMPRODS(I) specifies a sym prod.  Orbitals for this are stored in
C.. ORBPAIRS starting at index SYMPRODIND(1,I-1)+1 and have number SYMPRODIND(2,I-1)

C         DO I=1,L
C            WRITE(6,*) "OP",I,ORBPAIRS(1,I),ORBPAIRS(2,I)
C         ENDDO

C.. With each of the symmetry products, calculate []' such that []x[]'
C.. contains the totally symmetric rep.
C.. We do this by checking whether any of the []' in the 
C.. global symprods table multiply by [] to give the sym rep.

C.. I=SYMLABELPAIRLIST(1,SYMPROD) is the index of the set of symlabels with sym prod SYMPROD in
C.. SYMLABELPAIRS
C.. SYMLABELPAIRLIST(2,SYMPROD) is he total number of pairs [of states] with that SYMPROD
C.. [ SYMLABELPAIRS(1,I) , SYMLABELPAIRS(2,I) ] is the pair of states whose prod is the symlabel.
C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available.

C         WRITE(6,*) NPR, NSYMPAIRPRODS
C.. Now go through the list of all pairs, finding out how many pairs are to be excluded as they contain some of the
C.. orbitals in this determinant.
C.. We first create a quick lookup table to enable us to quickly check whether a given orbital is in this determinant
C.. in order 1, not order NEL.
C         WRITE(6,*) NBASIS
         CALL IAZZERO(ILUT,NBASIS/32+1)
         DO I=1,NEL
            ILUT(NI(I)/32)=IBSET(ILUT(NI(I)/32),MOD(NI(I),32)-1)
         ENDDO
C         DO I=0,NBASIS/32
C            WRITE(6,"(A,Z10)") "LUT: ",ILUT(I)
C         ENDDO
C.. Now look through the list of our pairs.  For each pair sym of the complete list which has a 
C.. symmetric product with any of our pair syms, we work out how many allowed pairs there are in
C.. the complete list, and store that value in NALLOWPPS
         IF(STORE(3).EQ.0) THEN
            CALL MEMORY(IP_NALLOWPPS,3*NSYMPAIRPRODS/IRAT+1,'NALLOWPPS')
         ELSE
            IP_NALLOWPPS=STORE(3)
         ENDIF
         NDOUB=0
         NEXCITTYPES=0
         DO I=1,NPR
            DO J=1,NSYMPAIRPRODS
C.. If it's the first time around, we flag NALLOWPPS as invalid
               IF(I.EQ.1) THEN
                  NALLOWPPS(1,J)=-1
                  NALLOWPPS(2,J)=-1
                  NALLOWPPS(3,J)=-1
               ENDIF
C.. For each sympairprod, check if [] x []' contains the sym rep
C               WRITE(6,*) SYMPRODS(I),SYMPAIRPRODS(J)
               SPP=SYMPAIRPRODS(J)
               IF(LSYMSYM(SYMPROD(SYMPRODS(I), SPP))) THEN
C.. Check if we've worked out the allowed number of excitation pairs, otherwise work it out
                  IF(NALLOWPPS(1,J).EQ.-1) THEN
C.. Zero the counters
                     NALLOWPPS(1,J)=0
                     NALLOWPPS(2,J)=0
                     NALLOWPPS(3,J)=0

                     DO K=SYMLABELPAIRLIST(1,SPP),
     &                 SYMLABELPAIRLIST(1,SPP)+SYMLABELPAIRLIST(2,SPP)-1

C.. WRONG if there are definitely no more disallowed pairs in the list, stop the loop
C                        IF(SYMLABELPAIRS(1,K).GT.NI(NEL)) EXIT

C.. Now check according to ISPN
C.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
C.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
                        ICC1=SYMLABELPAIRS(1,K)*2-1
                        ICC2=ICC1+1
                        ICC3=SYMLABELPAIRS(2,K)*2-1
                        ICC4=ICC3+1
                        L1B=.NOT.BTEST(ILUT(ICC1/32),
     &                        MOD(ICC1,32)-1)
                        L1A=.NOT.BTEST(ILUT(ICC2/32),
     &                        MOD(ICC2,32)-1)
                        L2B=.NOT.BTEST(ILUT(ICC3/32),
     &                        MOD(ICC3,32)-1)
                        L2A=.NOT.BTEST(ILUT(ICC4/32),
     &                        MOD(ICC4,32)-1)
C.. L1B is set if the beta of the first virtual is not in NI, i.e. is allowed
C.. The virtuals we are exciting to are the same state.  We're only allowed an AB excitation
                        IF(ICC1.EQ.ICC3) THEN
                           NALLOWPPS(1,J)=NALLOWPPS(1,J)+1
                           NALLOWPPS(3,J)=NALLOWPPS(3,J)+1
                         IF(L1B.AND.L2A) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                        ELSE
                         IF(L1B.AND.L2B) NALLOWPPS(1,J)=NALLOWPPS(1,J)+1
                         IF(L1B.AND.L2A) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                         IF(L1A.AND.L2B) NALLOWPPS(2,J)=NALLOWPPS(2,J)+1
                         IF(L1A.AND.L2A) NALLOWPPS(3,J)=NALLOWPPS(3,J)+1
                        ENDIF
C                        WRITE(6,*) SPP,SYMLABELPAIRS(1,K)*2-1,
C     &                     SYMLABELPAIRS(2,K)*2-1
C     &                     ,L1B,L1A,L2B,L2A
                     ENDDO
                  ENDIF
                  DO ISPN=1,3
                     K=NALLOWPPS(ISPN,J)
                     L=SYMPRODIND(2,ISPN,I)
                     ICC=L*K
C                     WRITE(6,"(2I5,2Z6,4I5)") I,J,SYMPRODS(I),
C     &                  SYMPAIRPRODS(J),ISPN,L,K,ICC
                     NDOUB=NDOUB+ICC
                     IF(ICC.GT.0) NEXCITTYPES=NEXCITTYPES+1
                  ENDDO
               ENDIF
            ENDDO
         ENDDO       

C         WRITE(6,*) "Number of double excitations: ",NDOUB
C.. Now calculate the number of single excitations

C.. SYMLABELLIST holds a list of states grouped under symmlabel
C.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
C.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I

         NSING=0
         DO I=1,NCL
C.. For each class, see what sym prods it can interact with
            DO J=1,NSYMLABELS
               IF(LSYMSYM(SYMPROD(SYMLABELS(2,J),
     &                            SYMLABELS(2,CLASSES(I))))) THEN
                  DO ISPN=1,2
                     ICC=CLASSCOUNT(ISPN,I)*
     &                  (SYMLABELCOUNTS(2,J)-CLASSCOUNT(ISPN,I))
                     NSING=NSING+ICC
                     IF(ICC.NE.0) NEXCITTYPES=NEXCITTYPES+1
C                     WRITE(6,"(2Z4,4I)") SYMLABELS(2,CLASSES(I)),
C     &                  SYMLABELS(2,J),ISPN,ICC,
C     &                  CLASSCOUNT(ISPN,I),
C     &                  (SYMLABELCOUNTS(2,J)-CLASSCOUNT(ISPN,I))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO         

C         WRITE(6,*) "Number of single excitations: ",NSING
         NEXCITS=NSING+NDOUB
         ICOUNT=NEXCITS
C         WRITE(6,*) "Total number of excitations: ",NEXCITS
C         WRITE(6,*) "Total number of excitation types: ",NEXCITTYPES

         IF(TCOUNT) THEN
C.. If we're just counting, we're done, so we get rid of some pointers.
            IP_EXCITTYPES=0
C.. However, we do save the length of the memory required.
C.. EXCITTYPES
            STORE(2)=(NEXCITTYPES*5)
C.. NALLOWPPS
            STORE(3)=3*NSYMPAIRPRODS
C.. ORBPAIRS
            STORE(4)=2*NPAIRS
C.. SYMPRODIND
            STORE(5)=2*3*(NPR+1)
C.. indicate that these are lengths
            STORE(6)=0
            
            CALL FREEM(IP_NALLOWPPS)
            CALL FREEM(IP_ORBPAIRS)
            CALL FREEM(IP_SYMPRODIND)
         ELSE
C.. Now allocate memory to store all the excitation types
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SYMLABELS(2,J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
         IF(STORE(2).EQ.0) THEN
          CALL MEMORY(IP_EXCITTYPES,(NEXCITTYPES*5)/IRAT+1,'EXCITTYPES')
         ELSE
            IP_EXCITTYPES=STORE(2)
         ENDIF

C.. Now store the singles:
         NEXCITTYPES=0
         DO I=1,NCL
C.. For each class, see what sym prods it can interact with
            DO J=1,NSYMLABELS
               IF(LSYMSYM(SYMPROD(SYMLABELS(2,J),
     &                            SYMLABELS(2,CLASSES(I))))) THEN
                  DO ISPN=1,2
                     ICC=CLASSCOUNT(ISPN,I)*
     &                  (SYMLABELCOUNTS(2,J)-CLASSCOUNT(ISPN,I))
                     IF(ICC.NE.0) THEN
                        NEXCITTYPES=NEXCITTYPES+1
                        EXCITTYPES(1,NEXCITTYPES)=1
                        EXCITTYPES(2,NEXCITTYPES)=ISPN
                        EXCITTYPES(3,NEXCITTYPES)=I
                        EXCITTYPES(4,NEXCITTYPES)=J
                        EXCITTYPES(5,NEXCITTYPES)=ICC
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO         
C.. And the doubles
         DO I=1,NPR
            DO J=1,NSYMPAIRPRODS
               SPP=SYMPAIRPRODS(J)
               IF(LSYMSYM(SYMPROD(SYMPRODS(I), SPP))) THEN
                  DO ISPN=1,3
                     K=NALLOWPPS(ISPN,J)
                     L=SYMPRODIND(2,ISPN,I)
                     ICC=L*K
                     IF(ICC.NE.0) THEN
                        NEXCITTYPES=NEXCITTYPES+1
                        EXCITTYPES(1,NEXCITTYPES)=2
                        EXCITTYPES(2,NEXCITTYPES)=ISPN
                        EXCITTYPES(3,NEXCITTYPES)=I
                        EXCITTYPES(4,NEXCITTYPES)=J
                        EXCITTYPES(5,NEXCITTYPES)=ICC
C                     WRITE(6,"(2I5,2Z6,4I5)") I,J,SYMPRODS(I),
C     &                  SYMPAIRPRODS(J),ISPN,L,K,ICC
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
C.. Store all the pointers we need
         STORE(2)=IP_EXCITTYPES
         STORE(3)=IP_NALLOWPPS
         STORE(4)=IP_ORBPAIRS
         STORE(5)=IP_SYMPRODIND
         STORE(6)=NEXCITTYPES
         ENDIF
C..      ENDIF(.NOT.TCOUNT)
         CALL TIHALT('SYMSUEXCIT',ISUB)
      END

      SUBROUTINE SYMGENEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,LSTE,
     &   ICLIST,ICOUNT)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,G1(5,NBASIS),NBASIS,NBASISMAX(*)
         INTEGER STORE(6)
         INTEGER ORBPAIRS(2,*)
         POINTER (IP_ORBPAIRS,ORBPAIRS)
         INTEGER SYMPRODIND(2,3,0:*)
         POINTER (IP_SYMPRODIND,SYMPRODIND)
         INTEGER NALLOWPPS(3,*),SPP
         POINTER(IP_NALLOWPPS,NALLOWPPS)

         INTEGER EXCITTYPES(5,*)
         POINTER (IP_EXCITTYPES,EXCITTYPES)
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER ISUB,I
         INTEGER LSTE(NEL,*),ICLIST(*),ICOUNT
         CALL TISET('SYMGNEXCIT',ISUB)

         IP_DSTORE=STORE(1)
         IP_EXCITTYPES=STORE(2)
         IP_NALLOWPPS=STORE(3)
         IP_ORBPAIRS=STORE(4)
         IP_SYMPRODIND=STORE(5)
         
         CALL SYMGENALLEXCITS(NI,NEL,EXCITTYPES,STORE(6),DSTORE(1),
     &            SYMPRODIND,DSTORE(NEL+1),ORBPAIRS,G1,NBASISMAX,
     &            LSTE,ICLIST,ICOUNT)
   
         CALL FREEM(IP_EXCITTYPES)
         CALL FREEM(IP_NALLOWPPS)
         CALL FREEM(IP_ORBPAIRS)
         CALL FREEM(IP_SYMPRODIND)
         CALL FREEM(IP_DSTORE)

         CALL TIHALT('SYMGNEXCIT',ISUB)
 
      END
      SUBROUTINE SYMGENALLEXCITS(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
     &               SYMPRODIND,ILUT,ORBPAIRS,G1,NBASISMAX,
     &               LSTE,ICLIST,ICOUNT)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES),NEXCITTYPES
         
         INTEGER I,J,K,L
         INTEGER NJ(NEL),NK(NEL)
         INTEGER CLASSES(*)
         INTEGER IEXCIT
         INTEGER IFROM,ITO,ISPN
         INTEGER IFROMSL,ITOSL
         INTEGER SPP
         LOGICAL L1,L2
         INTEGER ICC1,ICC2,ICC3,ICC4
         LOGICAL L1B,L1A,L2B,L2A

         INTEGER ORBPAIRS(2,*)
         INTEGER ILUT(0:*)
         INTEGER SYMPRODIND(2,3,0:*)
         INTEGER ISYM(5)
         INTEGER G1(*),NBASISMAX(*)
         INTEGER LSTE(NEL,*),ICLIST(*),ICOUNT
         INCLUDE 'sym.inc'

         ICOUNT=0

C.. Now work out all the excitations
         DO IEXCIT=1,NEXCITTYPES
            IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
C.. a single excitation
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SYMLABELS(2,J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
               ISPN=EXCITTYPES(2,IEXCIT)-2
               IFROM=CLASSES(EXCITTYPES(3,IEXCIT))
               ITO=EXCITTYPES(4,IEXCIT)
               J=0
               I=1
C               WRITE(6,*) ISPN,IFROM,ITO,EXCITTYPES(5,IEXCIT)
C.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
C.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I
               DO WHILE (I.LE.NEL)
                  L1=J.LT.SYMLABELCOUNTS(2,IFROM)
                  DO WHILE(L1)
                     IFROMSL=(SYMLABELLIST(SYMLABELCOUNTS(1,IFROM)+J)*2
     &                  +ISPN)
C                     WRITE(6,*) I,NI(I),J,IFROM,IFROMSL
                     IF(IFROMSL.LT.NI(I)) J=J+1
                     L1=IFROMSL.LT.NI(I)
     &                  .AND.J.LT.SYMLABELCOUNTS(2,IFROM)
                  ENDDO
C                  WRITE(6,*) I,J,SYMLABELCOUNTS(2,IFROM)
                  IF(J.LT.SYMLABELCOUNTS(2,IFROM).AND.IFROMSL.EQ.NI(I))
     &               THEN
C.. We've found an orb in NI with the correct sym.  Now go through the list of possible excitations of it
                     K=0
                     L=0
                     DO WHILE(K.LT.SYMLABELCOUNTS(2,ITO))
                        L2=.TRUE.
                        DO WHILE (L2)
                           L=L+1
                           ITOSL=2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)
     &                              +ISPN
                           L2=L.LE.NEL.AND.NI(L).LT.ITOSL
                        ENDDO
                        IF(L.GT.NEL.OR.NI(L).NE.ITOSL) THEN
C.. We've found a virtual into which we can excite our occupied orb.
                           CALL ICOPY(NEL,NI,1,NK,1)
                           NK(I)=ITOSL
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=1
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM(5) 
                        ENDIF
                        K=K+1
                     ENDDO
                  ENDIF
                  I=I+1
               ENDDO
            ELSE
C.. a double excitation
C.. We store each excitation type as:
C.. 1   TYPE (single=1, double=2)
C.. 2   SPIN (for single, 1=beta, 2=alpha.  For double, 1=beta/beta; 2=alpha/beta; 3=alpha/alpha;)
C.. 3   FROM (for single, I in CLASSES(I); for double, I in SYMPRODS(I) )
C.. 4   TO   (for single, J in SYMLABELS(2,J); for double, J in SYMPAIRPRODS(J) )
C.. 5  COUNT (Total number of excitations in this category)
               ISPN=EXCITTYPES(2,IEXCIT)
               IFROM=EXCITTYPES(3,IEXCIT)
               ITO=EXCITTYPES(4,IEXCIT)
C               WRITE(6,*) "EXC",IEXCIT,ISPN,IFROM,ITO,SYMPRODIND(2,ISPN,
C     &            IFROM)
C.. Go through the list of pairs with a given symprod.
C.. SYMPRODIND(1,ISPN,I)+1 contains the index of the first element of spin ISPN of sym
C.. SYMPRODS(I) in ORBPAIRS
C.. SYMPRODIND(2,ISPN,I) contains the number of such elements
               DO I=1,SYMPRODIND(2,ISPN,IFROM)
C.. Now go through the list of virtual pairs, excluding those with orbitals in NI 
                  SPP=SYMPAIRPRODS(ITO)
C                  WRITE(6,*) "SPP:", SPP
C.. I=SYMLABELPAIRLIST(1,SYMPROD) is the index of the set of symlabels with sym prod SYMPROD in
C.. SYMLABELPAIRS
C.. SYMLABELPAIRLIST(2,SYMPROD) is he total number of pairs [of states] with that SYMPROD
C.. [ SYMLABELPAIRS(1,I) , SYMLABELPAIRS(2,I) ] is the pair of states whose prod is the symlabel.
C.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available.
                  DO K=SYMLABELPAIRLIST(1,SPP),
     &              SYMLABELPAIRLIST(1,SPP)+SYMLABELPAIRLIST(2,SPP)-1
C.. Now check according to ISPN
C.. ICC1 is the beta orbital corresponding to the first state, and ICC2 the alpha
C.. ICC3 is the beta orbital corresponding to the  state, and ICC4 the alpha
                     ICC1=SYMLABELPAIRS(1,K)*2-1
                     ICC2=ICC1+1
                     ICC3=SYMLABELPAIRS(2,K)*2-1
                     ICC4=ICC3+1
C                     WRITE(6,*) ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
C     &                      ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)
C                     WRITE(6,*) ISPN,ICC1,ICC3
                     L1B=BTEST(ILUT(ICC1/32),
     &                     MOD(ICC1,32)-1)
                     L1A=BTEST(ILUT(ICC2/32),
     &                     MOD(ICC2,32)-1)
                     L2B=BTEST(ILUT(ICC3/32),
     &                     MOD(ICC3,32)-1)
                     L2A=BTEST(ILUT(ICC4/32),
     &                     MOD(ICC4,32)-1)
C                     WRITE(6,*) L1B,L1A,L2B,L2A
C.. L1B is set if the beta of the first virtual is in NI, i.e. is disallowed
                     IF(ISPN.EQ.1) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                        IF(ICC1.NE.ICC3.AND..NOT.(L1B.OR.L2B)) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC1
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC3
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM(5) 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ELSEIF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
                        IF(.NOT.(L1B.OR.L2A)) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC1
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC4
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM(5) 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
                        IF(.NOT.(L1A.OR.L2B).AND.ICC1.NE.ICC3) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC2
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC3
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM(5) 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ELSEIF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                        IF(ICC1.NE.ICC3.AND..NOT.(L1A.OR.L2A)) THEN
                           CALL ICOPY(NEL,NI,1,NK,1)
                           DO J=1,NEL
                              IF(NI(J).EQ.ORBPAIRS(1,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC2
                              ENDIF
                              IF(NI(J).EQ.ORBPAIRS(2,
     &                           SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                                 NK(J)=ICC4
                              ENDIF
                           ENDDO
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                           CALL SORTI(NEL,NK)
                           ICOUNT=ICOUNT+1
                           ICLIST(ICOUNT)=2
                           CALL ICOPY(NEL,NK,1,LSTE(1,ICOUNT),1)
C                           CALL WRITEDET(26,NK,NEL,.FALSE.)
C                           CALL GETSYM(NK,NEL,G1,NBASISMAX,ISYM)
C                           WRITE(26,'(Z6)') ISYM(5) 
C                           CALL WRITEDET(6,NK,NEL,.TRUE.)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO 
      END

C.. IF(TSETUP) Generate an iterator which allows up to double excitations to be generated
C.. one at a time (in an unordered fashion) from a given det.  THis needs to be called twice,
C.. first with STORE(1)=0, and then again.  Finally it can be called with TSETUP=.FALSE. 
C.. to actually generate the excitations
      SUBROUTINE GENSYMEXCITIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &         TSETUP,NMEM,NJ,IC,IFRZ,STORE)
         IMPLICIT NONE
         INCLUDE 'irat.inc'
         INTEGER NEL,NI(NEL)
         INTEGER G1(5,NBASIS),NBASISMAX(*),NBASIS
         POINTER (IP_DSTORE,DSTORE)
         INTEGER DSTORE(*)
         INTEGER STORE(6),ICOUNT
         INTEGER NMEM(*),NJ(NEL),IC,IFRZ,I
         LOGICAL TSETUP
         INTEGER GET_ADDR
         IF(TSETUP) THEN
C.. This is the first time we've been called for setup.
            IF(STORE(1).EQ.0) THEN
               CALL IAZZERO(STORE,6)
               CALL MEMORY(IP_DSTORE,NEL+NBASIS/32+1,'DSTORE')
               STORE(1)=IP_DSTORE
C.. Just count.
               CALL SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &            .TRUE.,ICOUNT,DSTORE(1), DSTORE(NEL+1))
               CALL FREEM(IP_DSTORE)
               STORE(1)=NEL+NBASIS/32+1
               NMEM(1)=6+STORE(1)+STORE(2)+STORE(3)+STORE(4)+STORE(5)
            ELSE
C.. The second setup.  Now NMEM is allocated, we store all the info
C.. NMEM is as follows:
C..   1           -  5           STORE(1:5)
C..   6           -  6           STORE(6) = NEXCITTYPES
C.. Data for the iterators
C..   7  IEXCIT
C..   8  ISPN
C..   9  IFROM
C..   10 ITO
C..   11 I
C..   12 J
C..   13 K
C..   14 L
C..   15 ICC(1:4)
C..   19 LS(1:2,1:2)  ([1 or 2], [A or B])
C..   23...
C..  (STORE(1)=23)-  STORE(2)-1  DSTORE
C..   STORE(2)    -  STORE(3)-1  EXCITTYPES
C..   STORE(3)    -  STORE(4)-1  NALLOWPPS
C..   STORE(4)    -  STORE(5)-1  ORBPAIRS
C..   STORE(5)    -  ...         SYMPRODIND

C..   DSTORE(1)   -  DSTORE(NEL) CLASSES
C..   DSTORE(NEL+1) - ,,,        ILUT

               CALL IAZZERO(NMEM,22)
               NMEM(11)=-1
               NMEM(7)=0
               ICOUNT=23
C.. Put the addresses in store
               DO I=1,5
                  NMEM(I)=GET_ADDR(NMEM(ICOUNT))
                  ICOUNT=ICOUNT+STORE(I)
               ENDDO
               IP_DSTORE=NMEM(1)
               CALL SYMSETUPEXCITS2(NI,NEL,G1,NBASIS,NBASISMAX,NMEM(1),
     &            .FALSE.,ICOUNT,DSTORE(1), DSTORE(NEL+1))
C.. and now instead the indices within NMEM
               ICOUNT=23
               DO I=1,5
                  NMEM(I)=ICOUNT
                  ICOUNT=ICOUNT+STORE(I)
               ENDDO
C.. Second setup finally complete.
            ENDIF
         ELSE
C.. Actually generate a det
            IP_DSTORE=GET_ADDR(NMEM(NMEM(1)))
            CALL SYMGENEXCITIT(NI,NEL,NMEM(NMEM(2)),NMEM(6),DSTORE(1),
     &         NMEM(NMEM(5)),DSTORE(NEL+1),NMEM(NMEM(4)),
     &         NMEM(7),NMEM(8),NMEM(9),NMEM(10),NMEM(11),NMEM(12),
     &         NMEM(13),NMEM(14),NMEM(15),NMEM(19),NJ,IC)
         ENDIF
      END
      SUBROUTINE RESETEXIT2(NI,NEL,G1,NBASIS,NBASISMAX,
     &            NMEM,IFRZ)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,G1(*),NBASIS,NBASISMAX(*)
         INTEGER NMEM(*),IFRZ(*)
         CALL IAZZERO(NMEM,22)
         NMEM(11)=-1
         NMEM(7)=0
      END
C                  CALL RESETEXIT(IPATH(1,LOCTAB(3,IVLEVEL)),NEL,G1,
C     &               NBASIS,NBASISMAX,CURGEN,IFRZ2)

      SUBROUTINE SYMGENEXCITIT(NI,NEL,EXCITTYPES,NEXCITTYPES,CLASSES,
     &               SYMPRODIND,ILUT,ORBPAIRS,IEXCIT,ISPN,IFROM,ITO,
     &               I,J,K,L,ICC,LS,
     &               NK,IC)
      IMPLICIT NONE
      INTEGER NEL,NI(NEL),EXCITTYPES(5,NEXCITTYPES),NEXCITTYPES
      
      INTEGER I,J,K,L
      INTEGER NK(NEL)
      INTEGER CLASSES(*)
      INTEGER IEXCIT
      INTEGER IFROM,ITO,ISPN
      INTEGER IFROMSL,ITOSL
      INTEGER SPP
      LOGICAL L1,L2
      INTEGER ICC(4)
C.. 1,1= 1B, 1,2=1A; 2,1=2B, 2,2=2A.
      LOGICAL LS(2,2)

      INTEGER ORBPAIRS(2,*)
      INTEGER ILUT(0:*)
      INTEGER SYMPRODIND(2,3,0:*)
      INTEGER IC
      INCLUDE 'sym.inc'

      DO WHILE(1)
C.. see if we need a new EXCIT
         IF(I.LT.0) THEN
C.. move to the next excitation
            IEXCIT=IEXCIT+1
            IF(IEXCIT.GT.NEXCITTYPES) THEN
C.. We're done
               NK(1)=0
               RETURN
            ENDIF
            IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
C.. a single
               ISPN=EXCITTYPES(2,IEXCIT)-2
               IFROM=CLASSES(EXCITTYPES(3,IEXCIT))
               ITO=EXCITTYPES(4,IEXCIT)
               L=0
               K=SYMLABELCOUNTS(2,ITO)
               J=0
               I=0
            ELSE
C.. a double
               ISPN=EXCITTYPES(2,IEXCIT)
               IFROM=EXCITTYPES(3,IEXCIT)
               ITO=EXCITTYPES(4,IEXCIT)
               I=0
               L=5
               K=-2
            ENDIF
         ENDIF
         IF(EXCITTYPES(1,IEXCIT).EQ.1) THEN
C.. singles
C.. We always need a new K.  K and L run in parallel
            K=K+1
C.. We've stored IFROMSL in ICC1
            IFROMSL=ICC(1)
C.. For each K, we check to see if it's valid
            IF(K.GE.SYMLABELCOUNTS(2,ITO)) THEN
C.. no more possible K, so we get a new I
               K=-1
               L=0
C.. I and J run in parallel, and for each J, we check whether it's in 
C.. the det.  If it is, then we carry on.
               L2=I.LT.NEL
               L1=.TRUE.
               DO WHILE(L2)
                  I=I+1
                  L1=J.LT.SYMLABELCOUNTS(2,IFROM)
                  DO WHILE(L1)
                     IFROMSL=(SYMLABELLIST(SYMLABELCOUNTS(1,IFROM)+J)*2
     &                     +ISPN)
                     ICC(1)=IFROMSL
C                     WRITE(6,*) I,NI(I),J,IFROM,IFROMSL
                     IF(IFROMSL.LT.NI(I)) J=J+1
                     L1=IFROMSL.LT.NI(I)
     &                     .AND.J.LT.SYMLABELCOUNTS(2,IFROM)
                  ENDDO
C.. see if we need a new I.
                  L2=IFROMSL.NE.NI(I).AND.I.LT.NEL
               ENDDO
C.. If we've gone too far, signal a new excit
               IF(L1.OR.(I.EQ.NEL.AND.IFROMSL.NE.NI(I))) I=-1
               CYCLE
            ENDIF
C.. it's a valid K, but is it in the det already?
C.. ITOSL is the orb it corresponds to
            ITOSL=2*SYMLABELLIST(SYMLABELCOUNTS(1,ITO)+K)
     &                              +ISPN
C.. Check if it's in the det
            L2=.TRUE.
            DO WHILE (L2)
               L=L+1
               L2=L.LE.NEL.AND.NI(L).LT.ITOSL
            ENDDO
            IF(L.LE.NEL.AND.NI(L).EQ.ITOSL) THEN
C.. We've found an L in the det which is the same as ITOSL. 
C.. we go round again, getting another K
               CYCLE
            ENDIF
C.. hoorah!  We've got an I in the det, and a K not in the det.  Create an excitation
            CALL ICOPY(NEL,NI,1,NK,1)
            NK(I)=ITOSL
            CALL SORTI(NEL,NK)
            IC=1
C.. quit the do loop
            EXIT
         ELSE
C.. doubles
            SPP=SYMPAIRPRODS(ITO)
C.. See if we need a new K.  L is the spin label, and goes from 1..4.
            IF(L.GT.4) THEN
               L=1
               K=K+1
               IF(K.LT.0.OR.
     &            K.GE.SYMLABELPAIRLIST(1,SPP)+SYMLABELPAIRLIST(2,SPP))
     &            THEN
C.. K is invalid.  Get a new I
                  K=1
                  I=I+1
                  IF(I.GT.SYMPRODIND(2,ISPN,IFROM)) THEN
C.. I is now invalid
                     I=-1
                     CYCLE
                  ELSE
C.. reset K
                     K=SYMLABELPAIRLIST(1,SPP)
                  ENDIF
               ENDIF
C.. We've got a new K, so we need to reset some variables
               ICC(1)=SYMLABELPAIRS(1,K)*2-1
               ICC(2)=ICC(1)+1
               ICC(3)=SYMLABELPAIRS(2,K)*2-1
               ICC(4)=ICC(3)+1
C                     WRITE(6,*) ORBPAIRS(1,SYMPRODIND(1,ISPN,IFROM)+I),
C     &                      ORBPAIRS(2,SYMPRODIND(1,ISPN,IFROM)+I)
C                     WRITE(6,*) ISPN,ICC1,ICC3
               LS(1,1)=BTEST(ILUT(ICC(1)/32),
     &                     MOD(ICC(1),32)-1)
               LS(1,2)=BTEST(ILUT(ICC(2)/32),
     &                     MOD(ICC(2),32)-1)
               LS(2,1)=BTEST(ILUT(ICC(3)/32),
     &                     MOD(ICC(3),32)-1)
               LS(2,2)=BTEST(ILUT(ICC(4)/32),
     &                     MOD(ICC(4),32)-1)
            ENDIF
C.. Now check for an excitation
            IF(L.EQ.1) THEN
               L=2
               IF(ISPN.EQ.1.AND.ICC(1).NE.ICC(3).AND.
     &            .NOT.(LS(1,1).OR.LS(1,2))) THEN
                  CALL ICOPY(NEL,NI,1,NK,1)
                  DO J=1,NEL
                     IF(NI(J).EQ.ORBPAIRS(1,
     &                  SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                        NK(J)=ICC(1)
                     ENDIF
                     IF(NI(J).EQ.ORBPAIRS(2,
     &                  SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                        NK(J)=ICC(3)
                     ENDIF
                  ENDDO
                  CALL SORTI(NEL,NK)
                  IC=2
                  EXIT
               ENDIF
            ENDIF
            IF(L.EQ.2) THEN
               L=3
               IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, then allow
                  IF(.NOT.(LS(1,1).OR.LS(2,2))) THEN
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           NK(J)=ICC(1)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           NK(J)=ICC(4)
                        ENDIF
                     ENDDO
                     CALL SORTI(NEL,NK)
                     IC=2
                     EXIT
                  ENDIF
               ENDIF
            ENDIF
            IF(L.EQ.3) THEN
               L=4
               IF(ISPN.EQ.2) THEN
C.. If neither virtuals are in NI, and they're not the same(which would give
C.. us the same excitation as previously), then allow
                 IF(.NOT.(LS(1,2).OR.LS(2,1)).AND.ICC(1).NE.ICC(3)) THEN
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           NK(J)=ICC(2)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           NK(J)=ICC(3)
                        ENDIF
                     ENDDO
                     CALL SORTI(NEL,NK)
                     IC=2
                     EXIT
                  ENDIF
               ENDIF
            ENDIF
            IF(L.EQ.4) THEN
               L=5
               IF(ISPN.EQ.3) THEN
C.. If both virtuals aren't the samem and neither are in NI, then allow
                 IF(ICC(1).NE.ICC(3).AND..NOT.(LS(1,2).OR.LS(2,2))) THEN
                     CALL ICOPY(NEL,NI,1,NK,1)
                     DO J=1,NEL
                        IF(NI(J).EQ.ORBPAIRS(1,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           NK(J)=ICC(2)
                        ENDIF
                        IF(NI(J).EQ.ORBPAIRS(2,
     &                     SYMPRODIND(1,ISPN,IFROM)+I)) THEN
                           NK(J)=ICC(4)
                        ENDIF
                     ENDDO
                     CALL SORTI(NEL,NK)
                     IC=2
                     EXIT
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO 
      END

