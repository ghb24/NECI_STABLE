C.. The number of symmetries (irreps), and a product table of irreps (the result being an irrep bitmask)
         INTEGER NSYM,SYMTABLE(NSYM,NSYM)
         POINTER (IP_SYMTABLE,SYMTABLE)
C..   SYMREPS(2,*) has two sets of data:
C.. SYMREPS(1,IBASISFN) contains the numnber of the representation
C.. of which IBASISFN is a part.
C.. SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
         INTEGER SYMREPS(2,*)
         POINTER (IP_SYMREPS,SYMREPS)
C.. Each basis fn's symmetry falls into a class SYMLABELS(IBASISFN).
C.. The characters of this class are stored in SMYLABELCHARS(1:NROT, SYMLABELS(IBASISFN))
C.. The total numebr of symmetry labels is NSYMLABELS
         INTEGER SYMLABELS(*)
         INTEGER SYMLABELCHARS(NROT,NSYMLABELS)
         INTEGER NSYMLABELS
         POINTER (IP_SYMLABELS,SYMLABELS)
         POINTER (IP_SYMLABELCHARS,SYMLABELCHARS)
C.. NROT is the number of symmetry operations
         INTEGER NROT
C.. All symmetries are decomposable into component irreps.
C.. The characters corresponding to each irrep are in IRREPCHARS
         INTEGER IRREPCHARS(NROT,NSYM)
         POINTER (IP_IRREPCHARS,IRREPCHARS)
         COMMON/SYM/ IP_SYMTABLE,NSYM,IP_SYMREPS,NROT,NSYMLABELS,
     &               IP_SYMLABELS,IP_SYMLABELCHARS,IP_IRREPCHARS
