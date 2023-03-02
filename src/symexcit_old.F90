!! NB SYMSETUPEXCITS uses a store which has integer(int64) to cope with 64-bit machines, and passes this directly to SYMSETUPEXCITS2
!!  Iterators which use GENSYMEXCITIT2 use an INTEGER STORE, and create an integer(int64) STORE2 internally, and pass this to SYMSETUPEXCITS2

!! SYMSETUPEXCITS takes an integer(int64) STORE(6) which stores the pointers used by SYMSETUPEXCITS2
!!   This STORE can then be passed to SYMGENEXCITS.
!.. Setup for symmetry routine below.
SUBROUTINE SYMSETUPEXCITS(NI, NEL, NBASIS, STORE, TCOUNT, ICOUNT, ILEVEL, iMinElec1, iMaxElec1)
    use SystemData, only: Symmetry, SymmetrySize, g1, BasisFN
    use SymData, only: SymClassSize
    IMPLICIT NONE
    INTEGER NEL, NI(NEL), NBASIS
    integer, POINTER :: DSTORE(:)
    INTEGER STORE(6)
    INTEGER ICOUNT
    LOGICAL TCOUNT
    INTEGER ILEVEL
    INTEGER iMinElec1, iMaxElec1

    external :: SYMSETUPEXCITS3

    STORE(1:6) = 0
    allocate(DSTORE(SymClassSize * NEL + (nBasis / 32) + 1 + SymmetrySize * (NEL * NEL + 1)))
    STORE(1) = 1 !Meaning the first element of DSTORE
! STORE(1) -->
    CALL SYMSETUPEXCITS3(NI, NEL, G1, NBASIS, STORE, STORE(1), STORE(1), STORE(1), &
                         STORE(1), TCOUNT, ICOUNT, DSTORE(1), DSTORE(SymClassSize * NEL + 1), &
                         DSTORE(SymClassSize * NEL + 1 + (nBasis / 32) + 1), &
                         ILEVEL, iMinElec1, iMaxElec1)
!.. If we're just counting, we don't need to keep DSTORE
    IF (TCOUNT) deallocate(dstore)
END

!.. IF(TSETUP) Generate an iterator which allows up to double excitations to be generated
!.. one at a time (in an unordered fashion) from a given det.  THis needs to be called twice,
!.. first with STORE(1)=0, and then again.  Finally it can be called with TSETUP=.FALSE.
!.. to actually generate the excitations
!! First run.  TSETUP=.TRUE.  STORE(1)=0.  NMEM can just be an INTEGER
!!  requires INTEGER STORE(6) to hold initialization data.  NMEM is returned containing the memory required to hold the excitation generator (call this nLength)
!! Second run,  TSETUP=.TRUE.  STORE(:) contains data from first run.  NMEM now is an array of nLenght INTEGERs.  The excitation generator data is stored in this.
!! Third run,   TSETUP=.FALSE.  will generate an excitation of NI, and put it in NEL.

SUBROUTINE GENSYMEXCITIT2(NI, NEL, G1, NBASIS, TSETUP, NMEM, NJ, IC, STORE, ILEVEL)
    use SystemData, only: Symmetry, BasisFN
    IMPLICIT NONE
    INTEGER NEL, NI(NEL), NBASIS
    INTEGER G1(nBasis)
!  STORE contains lengths of various components of the excitation generator
    INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
    INTEGER, target :: NMEM(*)
    INTEGER NJ(NEL), IC
    LOGICAL TSETUP
    INTEGER ILEVEL
    external :: GenSymExcitIt2Par
    CALL GenSymExcitIt2Par(NI, NEL, G1, NBASIS, TSETUP, NMEM, NJ, IC, STORE, ILEVEL, 1, NEL)
END
