         INCLUDE 'symtypes.inc'
         TYPE(BasisFN) FrozenSym
!.. The number of symmetries (irreps), and a product table of irreps (the result being an irrep bitmask)
         INTEGER NSYM
!   For translational symmetry groups, we need to know the periodic condition of propogation, as the 
!   multiplication of two irreps is equivalent to the addition of their propogators subject to the modulus 
!   of the identity operation (ie supercell dimension), and the k-vectors in the dimensions of the symmetry 
!   supercell.  The number of bits each property takes up.
         integer Nprop(3),KPntSym(3,*),PropBitLen
         INTEGER SymConjTab(nSym)
         POINTER (IP_SymConjTab,SymConjTab)
         TYPE(Symmetry) SYMTABLE(NSYM,NSYM)
         POINTER (IP_SYMTABLE,SYMTABLE),(IP_KPntSym,KPntSym)
!.. SYMREPS is used to group together degenerate sets of orbitals of the same sym
!.. (e.g. the six orbitals which might make up a T2g set), and is used for working 
!.. out the symmetry of a determinant in GETSYM
!.. It uses that fact that even for non-abelian groups a completely filled degenerate symmetry set is totally symmetric.
!..  Thus each member of a set of states which when completely filled gives a totally symmetric det should be labelled with the same symrep
!..   SYMREPS(2,*) has two sets of data:
!.. SYMREPS(1,IBASISFN) contains the numnber of the representation
!.. of which IBASISFN is a part.
!.. SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
         INTEGER SYMREPS(2,*)
         POINTER (IP_SYMREPS,SYMREPS)
!.. SymClasses is used to classify all states which transform with the same symmetry
!.. for the excitation generation routines
!.. Each state's symmetry falls into a class ISYMLABEL=SymClasses(ISTATE).
!.. The symmetry bit string, decomposing the sym label into its component irreps is in 
!.. SYMLABELS(ISYMLABEL)
!.. The characters of this class are stored in SYMLABELCHARS(1:NROT, SymClasses(ISTATE))
!.. The total number of symmetry labels is NSYMLABELS
         Type(Symmetry) SYMLABELS(*)
         POINTER (IP_SYMLABELS,SYMLABELS)
         INTEGER StateSymMap(*),StateSymMap2(*)
         POINTER (IP_StateSymMap,StateSymMap)
         Pointer (IP_StateSymMap2,StateSymMap2)
         INTEGER SymClasses(*)
         POINTER (IP_SymClasses,SymClasses)
         INTEGER SymClasses2(*)
         POINTER (IP_SymClasses2,SymClasses2)
         COMPLEX*16 SYMLABELCHARS(NROT,NSYMLABELS)
         INTEGER NSYMLABELS
         POINTER (IP_SYMLABELCHARS,SYMLABELCHARS)
!.. SYMLABELLIST holds a list of states grouped under symmlabel
         INTEGER SYMLABELLIST(*)
!.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
!.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I
!.. SYMLABELCOUNTSCUM(I) is the cumulative number of states with symlabel I
!.. SYMLABELINTSCUM(I) is the cumulative number of one-electron integrals with symlabel I

         INTEGER SYMLABELCOUNTS(2,*)
         POINTER (IP_SYMLABELLIST,SYMLABELLIST)
         POINTER (IP_SYMLABELCOUNTS,SYMLABELCOUNTS)
         INTEGER SYMLABELCOUNTSCUM(*)
         POINTER (IP_SYMLABELCOUNTSCUM,SYMLABELCOUNTSCUM)
         INTEGER SYMLABELINTSCUM(*)
         POINTER (IP_SYMLABELINTSCUM,SYMLABELINTSCUM)
!.. These ones are for when freezing orbitals
         INTEGER SYMLABELCOUNTSCUM2(*)
         POINTER (IP_SYMLABELCOUNTSCUM2,SYMLABELCOUNTSCUM2)
         INTEGER SYMLABELINTSCUM2(*)
         POINTER (IP_SYMLABELINTSCUM2,SYMLABELINTSCUM2)
         !.. NROT is the number of symmetry operations
         INTEGER NROT
!.. All symmetries are decomposable into component irreps.
!.. The characters corresponding to each irrep are in IRREPCHARS
         COMPLEX*16 IRREPCHARS(NROT,NSYM)
         POINTER (IP_IRREPCHARS,IRREPCHARS)

!.. SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs available, the number of pairs of
!.. states (listed in SymStatePairs), and the index of the start of this list
!.. For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
!.. [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose prod is of that symmetry.
         INTEGER SymStatePairs(2,0:*)
         TYPE(SymPairProd) SymPairProds(*)
         INTEGER nSymPairProds


         POINTER (IP_SymStatePairs,SymStatePairs)
         POINTER (IP_SYMPAIRPRODS,SYMPAIRPRODS)
         LOGICAL TAbelian  ! TAbelian for Abelian point groups (specifically k-point symmetry).
         COMMON/SYM/ IP_SYMTABLE,NSYM,Nprop,IP_SYMREPS,NROT,NSYMLABELS, &
     &               IP_SYMLABELS,IP_SYMLABELCHARS,IP_IRREPCHARS,       &
     &               IP_SymStatePairs,IP_KPntSym,                       &
     &               IP_SYMPAIRPRODS,                                   &
     &               IP_SYMLABELLIST,IP_SYMLABELCOUNTS,                 &
     &               IP_SymConjTab,IP_SYMLABELCOUNTSCUM,IP_SymClasses,  &
     &               IP_SYMLABELINTSCUM,IP_SYMLABELCOUNTSCUM2,          &
     &               IP_StateSymMap,IP_StateSymMap2,                    &
     &               IP_SYMLABELINTSCUM2,IP_SYMCLASSES2,                &
     &               FrozenSym  ,NSYMPAIRPRODS,TAbelian,PropBitLen
