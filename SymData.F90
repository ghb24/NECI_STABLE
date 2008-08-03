! 1. comment.
! 2. pointers-> allocatable
! 3. use statements.
! 4. memory {,de}allocation
module SymData

    use System, only: BasisFn,Symmetry,SymmetrySize

    implicit none

    TYPE SymPairProd 
        TYPE(Symmetry) Sym
        INTEGER nPairs
        INTEGER nIndex
    ENDTYPE
    ! The spacer is there to make sure we have a structure which is a multiple
    ! of 8-bytes for 64-bit machines.
    INTEGER, PARAMETER  :: SymPairProdSize=SymmetrySize+2

    TYPE SymClass 
        TYPE(Symmetry) SymRem
        INTEGER        SymLab
        INTEGER        spacer
    ENDTYPE
    ! The spacer is there to make sure we have a structure which is a multiple
    ! of 8-bytes for 64-bit machines.
    INTEGER, PARAMETER :: SymClassSize=2+SymmetrySize

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

    INTEGER :: SymConjTab(:) ! length=nSym

    TYPE(Symmetry) SYMTABLE(NSYM,NSYM)

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
    INTEGER SYMREPS(:,:) ! size=2,

    ! SymClasses is used to classify all states which transform with the same
    ! symmetry for the excitation generation routines.
    ! Each state's symmetry falls into a class ISYMLABEL=SymClasses(ISTATE).
    ! The symmetry bit string, decomposing the sym label into its component
    ! irreps is in SYMLABELS(ISYMLABEL).
    ! The characters of this class are stored in 
    ! SYMLABELCHARS(1:NROT, SymClasses(ISTATE)).
    ! The total number of symmetry labels is NSYMLABELS
    INTEGER NSYMLABELS
    Type(Symmetry) SYMLABELS(:)
    INTEGER StateSymMap(:),StateSymMap2(:)
    INTEGER SymClasses(:)
    INTEGER SymClasses2(:)
    COMPLEX*16 SYMLABELCHARS(NROT,NSYMLABELS)

    !.. SYMLABELLIST holds a list of states grouped under symmlabel
    INTEGER SYMLABELLIST(:)

    ! SYMLABELCOUNTS(1,I) is the index within SYMLABELLIST of the first state
    ! of symlabel I
    ! SYMLABELCOUNTS(2,I) is the number of states with symlabel I
    ! SYMLABELCOUNTSCUM(I) is the cumulative number of states with symlabel I
    ! SYMLABELINTSCUM(I) is the cumulative number of one-electron integrals with symlabel I
    INTEGER SYMLABELCOUNTS(:,:) ! size=2,
    INTEGER SYMLABELCOUNTSCUM(:)
    INTEGER SYMLABELINTSCUM(:)

    ! These ones are for when freezing orbitals
    INTEGER SYMLABELCOUNTSCUM2(:)
    INTEGER SYMLABELINTSCUM2(:)

    ! NROT is the number of symmetry operations
    INTEGER NROT
    ! All symmetries are decomposable into component irreps.
    ! The characters corresponding to each irrep are in IRREPCHARS
    COMPLEX*16 IRREPCHARS(:,:) ! size=NROT,NSYM

    ! SYMPAIRPRODS(1:NSYMPAIRPRODS) contains the list of all SYMPRODs
    ! available, the number of pairs of states (listed in SymStatePairs), and
    ! the index of the start of this list.
    ! For a given (unique) SymPairProds(J)%Sym, I=SymPairProds(J)%Index.
    ! [ SymStatePairs(1,I) , SymStatePairs(2,I) ] is the pair of states whose
    ! prod is of that symmetry.
    INTEGER SymStatePairs(:,:) ! shape=2,0:*
    TYPE(SymPairProd) SymPairProds(:)
    INTEGER nSymPairProds


    LOGICAL TAbelian  ! TAbelian for Abelian point groups (specifically k-point
                      ! symmetry).

end module SymData
