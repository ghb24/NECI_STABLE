#include "macros.h"
module gndts_blk_mod
    use SystemData, only: BasisFN, Symmetry, NullBasisFn, BasisFNSize
    use constants, only: dp
    use sort_mod, only: sort
    use sym_mod, only: GENNEXTSYM, SETUPSYM, GETSYMDEGEN, &
        GETSYM, LCHKSYM, WRITEALLSYM, getsym, ROUNDSYM, LCHKSYM
    use util_mod, only: NECI_ICOPY
    better_implicit_none
    private
    public :: gndts_blk, gensymdetss

contains

    SUBROUTINE GNDTS_BLK(NEL, NBASIS, BRR, NBASISMAX, NMRKS, TCOUNT, &
                         NDET, G1, II, NBLOCKSTARTS, NBLOCKS, TSPN, LMS, &
                         TPARITY, SymRestrict, IFDET, TGENFDET, NDETTOT, BLOCKSYM)
        INTEGER NEL, NBASIS, BRR(NBASIS), nBasisMax(5, *), NDET
        INTEGER NMRKS(NEL, NDET)
        INTEGER II, NLMAX
        INTEGER NBLOCKS, NBLOCKSTARTS(NBLOCKS + 1), OII
        TYPE(BASISFN) G1(NBASIS), ISYM, KJ, BLOCKSYM(NBLOCKS), IMAX(2)
        TYPE(BasisFN) SymRestrict
        LOGICAL TCOUNT
        INTEGER I, NI(NEL), NJ(NEL)
        LOGICAL TSPN, TPARITY, TDONE
        INTEGER LMS
        LOGICAL TGENFDET
        INTEGER IFDET, NDETTOT, IDEG
        real(dp) DETSC, TDETSC
        real(dp) CALCT
        DETSC = 1D200
        II = 0
        IF (TCOUNT) THEN
            NLMAX = 0
            NBLOCKS = 0
        ELSE
            NLMAX = NDET
        END IF
        I = 0
        OII = 0
        NDETTOT = 0
!.. set the comparison det to an invalid one, so all dets are counted
        NI(1) = 0
        IF (TCOUNT) OPEN(14, FILE="BLOCKS", STATUS="UNKNOWN")
        CALL GENNEXTSYM(NEL, NBASISMAX, TSPN, LMS, TPARITY, SymRestrict, .TRUE., TDONE, IMAX, ISYM)
        DO WHILE (.NOT. TDONE)
            CALL SETUPSYM(KJ)
            CALL GENSYMDETS_R(NI, ISYM, NEL, G1, BRR, NBASIS, NMRKS, II, NLMAX, NJ, KJ, 1, 1, NBASISMAX)
            IF (II /= OII) THEN
!.. we've found an occupied block
                I = I + 1
                CALL GETSYMDEGEN(ISYM, NBASISMAX, IDEG)
                NDETTOT = NDETTOT + (II - OII) * IDEG
                IF (TCOUNT) THEN
                    WRITE(14, "(I5)", advance='no') I
                    CALL WRITEALLSYM(14, ISYM, .FALSE.)
                    WRITE(14, "(2I10)") II - OII,(II - OII) * IDEG
                ELSE
                    NBLOCKSTARTS(I) = OII + 1
                    BLOCKSYM(I) = ISYM
                    IF (TGENFDET) THEN
                        TDETSC = CALCT(NMRKS(1:NEL, OII + 1), NEL)
                        IF (TDETSC < DETSC) THEN
                            IFDET = OII + 1
                            DETSC = TDETSC
                        END IF
                    END IF
                END IF
            END IF
            OII = II
            CALL GENNEXTSYM(NEL, NBASISMAX, TSPN, LMS, TPARITY, SymRestrict, .FALSE., TDONE, IMAX, ISYM)
        END DO
        NBLOCKS = I
        IF (.NOT. TCOUNT) NBLOCKSTARTS(I + 1) = II + 1
        IF (TCOUNT) CLOSE(14)
        IF (.NOT. TGENFDET) IFDET = 1
        RETURN
    END

    RECURSIVE SUBROUTINE GENSYMDETSSD_R(NI, KI, NEL, G1, BRR, NBASIS, LSTE, NLIST, NLMAX, NJ, KJ, NELEC, NBF, NBASISMAX)
        INTEGER NEL, NI(NEL), NBASIS, NLMAX, NLIST
        TYPE(BASISFN) G1(NBASIS), KI, KJ, KJ2
        INTEGER LSTE(NEL, NLMAX), NJ(NEL), NELEC, NBF
        INTEGER I, BRR(NBASIS), NN(NEL), nBasisMax(5, *)
        INTEGER IGETEXCITLEVEL, ICE
        DO I = NBF, NBASIS
            NJ(NELEC) = BRR(I)
            KJ2 = KJ
!.. Check if we've filled all the electrons
            IF (NELEC == NEL) THEN
                CALL GETSYM(NJ, NEL, G1, NBASISMAX, KJ2)
                CALL ROUNDSYM(KJ2, NBASISMAX)
                IF (LCHKSYM(KJ2, KI)) THEN
                    CALL NECI_ICOPY(NEL, NJ, 1, NN, 1)
                    call sort(nN)
                    ICE = IGETEXCITLEVEL(NI, NN, NEL)
                    IF (ICE > 0 .AND. ICE <= 2) THEN
!.. we've found a det with the right sym.
                        NLIST = NLIST + 1
                        IF (NLIST <= NLMAX) THEN
!.. if there's space, we save it
                            CALL NECI_ICOPY(NEL, NN, 1, LSTE(1, NLIST), 1)
                        END IF
                    END IF
                END IF
            ELSE
!.. otherwise we need to add more electrons:
                CALL GENSYMDETSSD_R(NI, KI, NEL, G1, BRR, NBASIS, LSTE, NLIST, NLMAX, NJ, KJ2, NELEC + 1, I + 1, NBASISMAX)
            END IF
        END DO
        RETURN
    END

    RECURSIVE SUBROUTINE GENSYMDETS_R(NI, KI, NEL, G1, BRR, NBASIS, LSTE, NLIST, NLMAX, NJ, KJ, NELEC, NBF, NBASISMAX)
        INTEGER NEL, NI(NEL), NBASIS, NLMAX, NLIST
        INTEGER LSTE(NEL, NLMAX), NJ(NEL), NELEC, NBF
        INTEGER I, J, BRR(NBASIS), NN(NEL), nBasisMax(5, *)
        LOGICAL LCMP
        TYPE(BASISFN) G1(NBASIS), KI, KJ, KJ2
        DO I = NBF, NBASIS
            NJ(NELEC) = BRR(I)
            KJ2 = KJ
!.. Check if we've filled all the electrons
            IF (NELEC == NEL) THEN
                CALL GETSYM(NJ, NEL, G1, NBASISMAX, KJ2)
                IF (LCHKSYM(KJ2, KI)) THEN
!                  CALL ROUNDSYM(KJ3,NBASISMAX)
!.. we've found a det with the right sym.
                    NLIST = NLIST + 1
                    CALL NECI_ICOPY(NEL, NJ, 1, NN, 1)
                    call sort(nN)
!.. Just check to see if it's our original det
                    LCMP = .TRUE.
                    DO J = 1, NEL
                        IF (NN(J) /= NI(J)) LCMP = .FALSE.
                    END DO
                    IF (LCMP) THEN
!.. Roll back NLIST
                        NLIST = NLIST - 1
                    ELSEIF (NLIST <= NLMAX) THEN
!.. if there's space, we save it
                        CALL NECI_ICOPY(NEL, NN, 1, LSTE(1, NLIST), 1)
                    END IF
                END IF
            ELSE
!.. otherwise we need to add more electrons:
                CALL GENSYMDETS_R(NI, KI, NEL, G1, BRR, NBASIS, LSTE, NLIST, NLMAX, NJ, KJ2, NELEC + 1, I + 1, NBASISMAX)
            END IF
        END DO
        RETURN
    END

!.. Generate determinants with a given symmetry, given by KI, as GENSYMDETS
    SUBROUTINE GENSYMDETSS(KI, NEL, G1, BRR, NBASIS, LSTE, NLIST, NBASISMAX)
        INTEGER NEL, NI(NEL), NBASIS, BRR(NBASIS)
        INTEGER NLIST, LSTE(NEL, NLIST)
        TYPE(BASISFN) G1(NBASIS), KI, KJ
        INTEGER NJ(NEL), nBasisMax(5, *)
        INTEGER NLMAX
        KJ = NullBasisFn
        NI(1:NEL) = 0
        NLMAX = NLIST
        NLIST = 0
        CALL GENSYMDETS_R(NI, KI, NEL, G1, BRR, NBASIS, LSTE, NLIST, NLMAX, NJ, KJ, 1, 1, NBASISMAX)
        RETURN
    END

end module gndts_blk_mod
