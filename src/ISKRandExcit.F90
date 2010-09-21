#include "macros.h"

MODULE ISKRandExcit 
!ISK (Inversion-Symmetry K-point) wavefunctions are a linear combinations of two determinants, 
!where all orbitals are exchanged with 
!the equivalent orbital where the k-points of all occupied orbitals have been switched.
!In simple notation we will consider an excitation where determinants i and j are in the original ISK,
!and determinants a and b in the excited ISK.
!As with HPHF wavefunctions, a large simplification occurs when it is realised that P(i->a) = P(j->b),
!therefore all excitation are connected to the determinantal excitations of just one of the constituent determinants.
!This means that only one constituent determinant will be considered in the space.

    use SystemData, only: nel, tCSF, Alat, G1, nbasis, nbasismax, nmsh, arr
    use IntegralsData, only: UMat, fck, nMax
    use SymData, only: nSymLabels
    use dSFMT_interface, only : genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: gen_rand_excit, ConstructClassCounts, &
                                    CalcNonUniPGen, ScratchSize 
    use DetBitOps, only: DetBitLT, DetBitEQ, FindExcitBitDet, &
                         FindBitExcitLevel,MaskAlpha,MaskBeta
    use FciMCData, only: pDoubles
    use constants, only: dp,n_int
    use HElem
    use sltcnd_mod, only: sltcnd_excit
    use bit_reps, only: NIfD, NIfDBO, NIfTot
    use sort_mod
    IMPLICIT NONE

    contains

