C.. The number of symmetries (irreps), and a product table of irreps (the result being an irrep bitmask)
         INTEGER NSYM,SYMTABLE(NSYM,NSYM)
         POINTER (IP_SYMTABLE,SYMTABLE)
C.. SYMREPS is used to group together degenerate sets of orbitals of the same sym
C.. (e.g. the six orbitals which might make up a T2g set), and is used for working 
C.. out the symmetry of a determinant in GETSYM
C..   SYMREPS(2,*) has two sets of data:
C.. SYMREPS(1,IBASISFN) contains the numnber of the representation
C.. of which IBASISFN is a part.
C.. SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
         INTEGER SYMREPS(2,*)
         POINTER (IP_SYMREPS,SYMREPS)
C.. SYMLABELS is used to classify all states which transform with the same symmetry
C.. for the excitation generation routines
C.. Each state's symmetry falls into a class SYMLABELS(1,ISTATE).
C.. The symmetry bit string, decomposing the sym label into its component irreps is in 
C.. SYMLABELS(2,ISYMLABEL)
C.. The characters of this class are stored in SYMLABELCHARS(1:NROT, SYMLABELS(1,ISTATE))
C.. The total number of symmetry labels is NSYMLABELS
         INTEGER SYMLABELS(2,*)
         COMPLEX*16 SYMLABELCHARS(NROT,NSYMLABELS)
         INTEGER NSYMLABELS
         POINTER (IP_SYMLABELS,SYMLABELS)
         POINTER (IP_SYMLABELCHARS,SYMLABELCHARS)
C.. SYMLABELLIST holds a list of states grouped under symmlabel
         INTEGER SYMLABELLIST(*)
C.. SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state of symlabel I
C.. SYMLABELCOUNTS(2,I) is the number of states with symlabel I
         INTEGER SYMLABELCOUNTS(2,*)
         POINTER (IP_SYMLABELLIST,SYMLABELLIST)
         POINTER (IP_SYMLABELCOUNTS,SYMLABELCOUNTS)
C.. NROT is the number of symmetry operations
         INTEGER NROT
C.. All symmetries are decomposable into component irreps.
C.. The characters corresponding to each irrep are in IRREPCHARS
         COMPLEX*16 IRREPCHARS(NROT,NSYM)
         POINTER (IP_IRREPCHARS,IRREPCHARS)
         INTEGER SYMLABELPAIRLIST(2,0:*),SYMLABELPAIRS(2,0:*)
         INTEGER SYMPAIRPRODS(*),NSYMPAIRPRODS
         POINTER (IP_SYMLABELPAIRLIST,SYMLABELPAIRLIST)
         POINTER (IP_SYMLABELPAIRS,SYMLABELPAIRS)
         POINTER (IP_SYMPAIRPRODS,SYMPAIRPRODS)
         LOGICAL TNEWEXCITATIONS
         COMMON/SYM/ IP_SYMTABLE,NSYM,IP_SYMREPS,NROT,NSYMLABELS,
     &               IP_SYMLABELS,IP_SYMLABELCHARS,IP_IRREPCHARS,
     &               IP_SYMLABELPAIRLIST,IP_SYMLABELPAIRS,
     &               IP_SYMPAIRPRODS,NSYMPAIRPRODS,
     &               IP_SYMLABELLIST,IP_SYMLABELCOUNTS,
     &               TNEWEXCITATIONS
