! 1. comment.
! 2. pointers-> allocatable
! 3. use statements.
! 4. memory {,de}allocation
module SymData

    use SystemData, only: BasisFn,Symmetry,SymmetrySize

    implicit none

    save

    !Hold information about pairs of orbitals according to the symmetry of their product.
    TYPE SymPairProd 
        TYPE(Symmetry) Sym  !The symmetry of a set of pairs
        INTEGER nPairs      !The number of pairs where the direct product of the syms of the orbs in the pair is Sym
        INTEGER nIndex      !The index of the first of these pairs in the complete list of pairs (SymStatePairs)
        INTEGER nPairsStateSS !For abelian, SymPairProds actually holds info on pairs of sym classes not individual
                              !states.  This is the number of pairs of states which could be generated.            
                              !  This is if the states have the same spin 
        INTEGER nPairsStateOS !For abelian, SymPairProds actually holds info on pairs of sym classes not individual 
                              !states.  This is the number of pairs of states which could be generated.          
                              !  This is if the states have opposite spin                                     
    ENDTYPE
    INTEGER, PARAMETER :: SymPairProdSize=SymmetrySize+2

    !Used for  SymSetupExcits* to hold internal data.  Each symclass corresponds to a single symmetry of an orbital which is present in the determinant being excited from
    TYPE SymClass 
        TYPE(Symmetry) SymRem    ! The symmetry remaining from the determinant we excite from when one orbital from this symmetry is removed
        INTEGER        SymLab    ! The symmetry label (i.e. from SymClasses) corresponding to the symmetry of the orbital being removed.
        INTEGER        spacer    ! The spacer is there to make sure we have a structure which is a multiple of 8-bytes for 64-bit machines.
    ENDTYPE
    INTEGER, PARAMETER :: SymClassSize=2+SymmetrySize

    !  The resultant symmetry of any frozen orbitals
    TYPE(BasisFN) :: FrozenSym

    ! The number of symmetries (irreps), and a product table of irreps (the
    ! result being an irrep bitmask)
    INTEGER :: NSYM

    ! For translational symmetry groups:
    ! We need to know the periodic condition of propogation, as the
    ! multiplication of two irreps is equivalent to the addition of their
    ! propogators subject to the modulus of the identity operation (ie
    ! supercell dimension)
    INTEGER :: Nprop(3)
    ! and the k-vectors in the dimensions of the symmetry supercell  
    INTEGER, ALLOCATABLE :: KPntSym(:,:) ! size=3,nKP
    ! and the number of bits each property takes up.
    INTEGER :: PropBitLen

    ! The symmetry conjugate of each symmetry (i.e. its index in the list of irreps)
    INTEGER, ALLOCATABLE :: SymConjTab(:) ! length=nSym

    TYPE(Symmetry), ALLOCATABLE :: SYMTABLE(:,:) ! size=NSYM,NSYM

    ! SYMREPS is used to group together degenerate sets of orbitals of the same
    ! sym (e.g. the six orbitals which might make up a T2g set), and is used
    ! for working out the symmetry of a determinant in GETSYM It uses that fact
    ! that even for non-abelian groups a completely filled degenerate symmetry
    ! set is totally symmetric.  Thus each member of a set of states which when
    ! completely filled gives a totally symmetric det should be labelled with
    ! the same symrep
    ! SYMREPS(2,:) has two sets of data:
    !     SYMREPS(1,IBASISFN) contains the numnber of the representation
    !                         of which IBASISFN is a part.
    !     SYMPREPS(2,IREP) contains the degeneracy of the rep IREP
    INTEGER, ALLOCATABLE :: SYMREPS(:,:) ! size=2,

    ! The total number of symmetry labels is NSYMLABELS
    INTEGER :: NSYMLABELS
    ! The symmetry bit string, decomposing the sym label into its component
    ! irreps is in SYMLABELS(ISYMLABEL).
    Type(Symmetry), ALLOCATABLE :: SYMLABELS(:)
    INTEGER, ALLOCATABLE ::  StateSymMap(:),StateSymMap2(:)
    ! SymClasses is used to classify all states which transform with the same
    ! symmetry for the excitation generation routines.
    ! Each state's symmetry falls into a class ISYMLABEL=SymClasses(ISTATE).
    INTEGER, POINTER ::  SymClasses(:)
    INTEGER, POINTER ::  SymClasses2(:)
    ! The characters of this class are stored in 
    ! SYMLABELCHARS(1:NROT, SymClasses(ISTATE)).
    COMPLEX*16, ALLOCATABLE ::  SYMLABELCHARS(:,:) ! size=NROT,NSYMLABELS

    !.. SYMLABELLIST holds a list of states grouped under symmlabel
    INTEGER, ALLOCATABLE ::  SYMLABELLIST(:)

    ! SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state
    ! of symlabel I
    ! SYMLABELCOUNTS(2,I) is the number of states with symlabel I
    ! SYMLABELCOUNTSCUM(I) is the cumulative number of states with symlabel I
    ! SYMLABELINTSCUM(I) is the cumulative number of one-electron integrals with symlabel I
    INTEGER, ALLOCATABLE ::  SYMLABELCOUNTS(:,:) ! size=2,
    INTEGER, POINTER ::  SYMLABELCOUNTSCUM(:)
    INTEGER, POINTER ::  SYMLABELINTSCUM(:)

    ! These ones are for when freezing orbitals
    INTEGER, POINTER ::  SYMLABELCOUNTSCUM2(:)
    INTEGER, POINTER ::  SYMLABELINTSCUM2(:)

    ! NROT is the number of symmetry operations
    INTEGER :: NROT
    ! All symmetries are decomposable into component irreps.
    ! The characters corresponding to each irrep are in IRREPCHARS
    COMPLEX*16, ALLOCATABLE ::  IRREPCHARS(:,:) ! size=NROT,NSYM

    ! SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs
    ! available, the number of pairs of states (listed in SymStatePairs), and
    ! the index of the start of this list.
    TYPE(SymPairProd), ALLOCATABLE :: SymPairProds(:)
    INTEGER :: nSymPairProds
    ! For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
    ! [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose
    ! prod is of that symmetry.
    INTEGER, ALLOCATABLE ::  SymStatePairs(:,:) ! shape=2,0:*


    LOGICAL TAbelian  ! TAbelian for Abelian point groups (specifically k-point
                      ! symmetry).

    LOGICAL tAbelianFastExcitGen
    ! tAbelianFastExcitGen is a temporary flag. 
    !  It is used to indicate that if an Abelian symmetry group is present
    !   the excitation generators should use optimized routines
    !   to take this into account.  Not all excitation generator functions
    !   currently work with this.  USE WITH CARE


    ! Memory logging tags.
    INTEGER :: tagKPntSym 
    INTEGER :: tagSymConjTab
    INTEGER :: tagSYMTABLE
    INTEGER :: tagSYMREPS 
    INTEGER :: tagSYMLABELS
    INTEGER :: tagStateSymMap
    INTEGER :: tagStateSymMap2
    INTEGER :: tagSymClasses
    INTEGER :: tagSymClasses2
    INTEGER :: tagSYMLABELCHARS
    INTEGER :: tagSYMLABELLIST
    INTEGER :: tagSYMLABELCOUNTS
    INTEGER :: tagSYMLABELCOUNTSCUM
    INTEGER :: tagSYMLABELINTSCUM
    INTEGER :: tagSYMLABELCOUNTSCUM2
    INTEGER :: tagSYMLABELINTSCUM2
    INTEGER :: tagIRREPCHARS
    INTEGER :: tagSymStatePairs
    INTEGER :: tagSymPairProds

end module SymData
