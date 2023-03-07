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
end subroutine
