#include "macros.h"

!! NB SYMSETUPEXCITS uses a store which has integer(int64) to cope with 64-bit machines, and passes this directly to SYMSETUPEXCITS2
!!  Iterators which use GENSYMEXCITIT2 use an INTEGER STORE, and create an integer(int64) STORE2 internally, and pass this to SYMSETUPEXCITS2

!! SYMSETUPEXCITS takes an integer(int64) STORE(6) which stores the pointers used by SYMSETUPEXCITS2
!!   This STORE can then be passed to SYMGENEXCITS.
!.. Setup for symmetry routine below.
SUBROUTINE SYMSETUPEXCITS(NI, NEL, NBASIS, STORE, TCOUNT, ICOUNT, ILEVEL, iMinElec1, iMaxElec1)
    use SystemData, only: Symmetry, BasisFN
    use error_handling_neci, only: stop_all
    IMPLICIT NONE
    INTEGER NEL, NI(NEL), NBASIS
    INTEGER STORE(6)
    INTEGER ICOUNT
    LOGICAL TCOUNT
    INTEGER ILEVEL
    INTEGER iMinElec1, iMaxElec1
    routine_name("SYMSETUPEXCITS")
    unused_var(NI)
    unused_var(NEL)
    unused_var(NBASIS)
    unused_var(STORE)
    unused_var(TCOUNT)
    unused_var(ICOUNT)
    unused_var(ILEVEL)
    unused_var(iMinElec1)
    unused_var(iMaxElec1)
    call stop_all(this_routine, 'See if this routine is actually used.')
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
    use error_handling_neci, only: stop_all
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
