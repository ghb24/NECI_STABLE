#include "macros.h"
module hfbasis_mod
    use Determinants, only: get_helement, orderbasis
    use DeterminantData, only: write_det, write_det_len
    use HElem, only: HElement_t_size
    use IntegralsData, only: DMatEpsilon
    use LoggingData, only: HFLogLevel
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use OneEInts, only: GetTMATEl, TMAT2D, TMAT2D2
    use SystemData, only: BasisFN, tMolpro, IPARITY, nbasisMax, symmetry, &
        symmax, brr, g1, nbasis, lms, nel, tHub, tGUGA, t_k_space_hubbard, &
        get_basisfn
    use UMatCache, only: UMatInd, GTID
    use constants, only: dp, int32, stdout
    use dSFMT_interface, only: genrand_real2_dSFMT
    use error_handling_neci, only: neci_flush, stop_all
    use global_utilities, only: timer, set_timer, halt_timer
    use lattice_mod, only: lat
    use procedure_pointers, only: get_umat_el
    use sort_mod, only: sort
    use sym_mod, only: GENMOLPSYMTABLE, GENMOLPSYMREPS, TotSymRep
    use sym_mod, only: checkMomentumInvalidity, MomPbcSym
    use sym_mod, only: getsym, writesym, symprod
    use basic_float_math, only: near_zero
    use util_mod, only: int_fmt

    use scrtransf_mod, only: gethelement2t

    use blas_interface_mod, only: dcopy, dgemm, dsyev, zheev
    use Orthonorm_mod, only: GRAMSCHMIDT_NECI, LOWDIN_ORTH

    use scrtransf_mod, only: GETTRTMATEL, GETTRUMATEL

    use frsblk_mod, only: neci_write_matrix
    better_implicit_none

    private
    public :: CALCHFTMAT, setuphfbasis, readhftmat, &
        readhfumat, orderbasishf, calchfumat, readhfbasis, iFindBasisFn

#ifndef CMPLX_
    public :: calchfbasis
#endif

contains

!.. Ti=Sum_a,b(c_ia* c_ib <u_a(1)|h(1)|u_b(1)>)
!..   =Sum_a(|c_ia|^2 <u_a|(p^2)/2|u_a>)
    SUBROUTINE CALCHFTMAT(NBASIS, HFBASIS, NORBUSED)
        INTEGER :: NBASIS
        real(dp) HFBASIS(NBASIS, NBASIS)
        INTEGER I, J, A, B, NORBUSED
        HElement_t(dp) SUM1
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
        OPEN(10, FILE='TMAT2', STATUS='UNKNOWN')
        DO I = 1, NORBUSED
        DO J = 1, NORBUSED
            SUM1 = 0.0_dp
            DO A = 1, NBASIS
            DO B = 1, NBASIS
                SUM1 = SUM1 + HFBASIS(I, A) * HFBASIS(J, B) * (GetTMATEl(A, B))
            END DO
            END DO
            TMAT2D2(I, J) = SUM1
            IF (ABS(SUM1) > 1.0e-10_dp) WRITE(10, *) I, J, TMAT2D2(I, J)
        END DO
        END DO
        CLOSE(10)
        RETURN
    END
!.. Ti=Sum_a,b(c_ia* c_ib <u_a(1)|h(1)|u_b(1)>)
!..   =Sum_a(|c_ia|^2 <u_a|(p^2)/2|u_a>)
    SUBROUTINE READHFTMAT(NBASIS)
        INTEGER NBASIS
        real(dp) R
        INTEGER I, J
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
        OPEN(10, FILE='TMAT2', STATUS='OLD')
        I = 0
        DO WHILE (.NOT.(I == NBASIS .AND. J == NBASIS))
            READ(10, *, END=11) I, J, R
            TMAT2D2(I, J) = R
            TMAT2D2(J, I) = R
        END DO
11      CLOSE(10)
        RETURN
    END
    SUBROUTINE READHFUMAT(UMAT2, NBASIS)
        INTEGER NBASIS
        real(dp) UMAT2(*)
        real(dp) SUM1
        INTEGER A, B, C, D
        INTEGER NHG
        NHG = NBASIS
        WRITE(6, *) 'READING HF UMAT'
! ==--------------------------------------------------------------------=
        OPEN(10, FILE='UMAT2', STATUS='OLD')
        DO WHILE (.TRUE.)

            READ(10, *, END=12) A, B, C, D, SUM1
            UMAT2(UMatInd(A, B, C, D)) = SUM1
!..Symmetries not needed by UMatInd
!              UMAT2(C,D,A,B)=SUM1
!              UMAT2(B,A,D,C)=SUM1
!              UMAT2(D,C,B,A)=SUM1
        END DO
12      CLOSE(10)
    END

    SUBROUTINE SETUPHFBASIS(NBASISMAX, G1, NBASIS, HFE, ARR, BRR)
        INTEGER NBASIS, nBasisMax(5, *)
        TYPE(BasisFN) G1(nBasis)
        real(dp) ARR(NBASIS, 2)
        real(dp) HFE(NBASIS)
        INTEGER BRR(NBASIS), ORBORDER(8, 2)
        INTEGER I
#ifdef CMPLX_
        call stop_all('SETUPHFBASIS', 'HF not implemented for complex orbitals.')
#endif

!.. We now need to modify G1.  Pretend that all basis functions have
!.. a different X quantum number, and all have the same spin.
        DO I = 1, NBASIS
            G1(I)%k(1) = 0
            G1(I)%Ms = -2 * (MOD(I, 2)) + 1
            G1(I)%Sym = TotSymRep()
            ARR(I, 1) = HFE(I)
            ARR(I, 2) = HFE(I)
            BRR(I) = I
        END DO
!.. Now modify NBASISMAX
        NBASISMAX(4, 2) = 1
        NBASISMAX(4, 1) = -1
        NBASISMAX(3, 2) = 0
        NBASISMAX(3, 1) = 0
        NBASISMAX(2, 1) = 0
        NBASISMAX(2, 2) = 0
        NBASISMAX(1, 1) = 0
!         NBASISMAX(1,2)=NBASIS-1
        NBASISMAX(5, 2) = 0

!.. Generic spatial symmetry
        NBASISMAX(3, 3) = 1
        CALL GENMOLPSYMTABLE(1, G1, NBASIS)
        CALL ORDERBASIS(NBASIS, ARR, BRR, ORBORDER, NBASISMAX, G1)
        CALL GENMOLPSYMREPS()

        RETURN
    END
    SUBROUTINE CALCHFUMAT(UMAT, UMAT2, NBASIS, HFBASIS, ISS, NORBUSED)
        INTEGER NBASIS, ISS
        real(dp) UMAT(*)
        real(dp) UMAT2(*)
        real(dp) UMATT((((NBASIS / ISS) * (NBASIS / ISS - 1)) / 2)**2)
        real(dp) HFBASIS(NBASIS, NBASIS), SUM1
        INTEGER I, J, K, L, A, B, C, D
        type(timer), save :: proc_timer
        INTEGER ID1, ID2, ID3, ID4, NHG, NORBUSED
        LOGICAL LSPN
        character(*), parameter :: this_routine = 'CALCHFUMAT'
        NHG = NBASIS
        proc_timer%timer_name = 'CALCHFUMAT'
        call set_timer(proc_timer)
        WRITE(6, *) 'CALCULATING HF UMAT'
        call stop_all(this_routine, "HF UMAT calc broken through U/TMAT reindexing. Please fix")
#ifdef CMPLX_
        call stop_all('CALCHFUMAT', 'HF not implemented for complex orbitals.')
#endif
! ==--------------------------------------------------------------------==
        OPEN(10, FILE='UMAT2', STATUS='UNKNOWN')
!.. A, B, C, D denote basis fns in the HF basis.  These basis fns
!.. alternate in spin i.e. fn 1 has alpha, fn 2 beta, fn 3 alpha etc.
!.. We need to take into account the spin when writing out the new U
!.. matrix.  For <A(1) B(2) |U| C(1) D(2)> to be non-zero, A and C must
!.. have the same spin, and B and D must have the same spin.
!.. Thus the sum1 A+C must be even, as must B+D
        WRITE(6, *) "Index 1..."
        DO A = 1, NORBUSED
        DO J = 1, NHG / ISS
        DO K = 1, NHG / ISS
        DO L = 1, NHG / ISS
            SUM1 = 0.0_dp
            DO I = 1, NHG
                ID1 = GTID(I)
                SUM1 = SUM1 + HFBASIS(A, I) * UMAT(UMatInd(ID1, J, K, L))
            END DO
            UMATT(UMatInd(J, K, L, A)) = SUM1
!             IF(ABS(SUM1).GT.1.0e-9_dp) WRITE(6,*) J,K,L,A,SUM1
        END DO
        END DO
        END DO
        END DO
        WRITE(6, *) "Index 2..."
        DO A = 1, NORBUSED
        DO B = 1, NORBUSED
        DO K = 1, NHG / ISS
        DO L = 1, NHG / ISS
            SUM1 = 0.0_dp
            DO J = 1, NHG
                ID2 = GTID(J)
                SUM1 = SUM1 + HFBASIS(B, J) * UMATT(UMatInd(ID2, K, L, A))
            END DO
            UMAT2(UMatInd(K, L, A, B)) = SUM1
!             IF(ABS(SUM1).GT.1.0e-9_dp) WRITE(6,*) K,L,A,B,SUM1
        END DO
        END DO
        END DO
        END DO
        WRITE(6, *) "Index 3..."
        DO A = 1, NORBUSED
        DO B = 1, NORBUSED
        DO C = 1, NORBUSED
        DO L = 1, NHG / ISS
            SUM1 = 0.0_dp
            DO K = 1, NHG
                ID3 = GTID(K)
                SUM1 = SUM1 + HFBASIS(C, K) * UMAT2(UMatInd(ID3, L, A, B))
            END DO
            UMATT(UMatInd(L, A, B, C)) = SUM1
!             IF(ABS(SUM1).GT.1.0e-9_dp) WRITE(6,*) L,A,B,C,SUM1
        END DO
        END DO
        END DO
        END DO
        UMAT2(1:nBasis**4) = 0.0_dp
        WRITE(6, *) "Index 4..."
        DO A = 1, NORBUSED
        DO B = 1, NORBUSED
        DO C = 1, NORBUSED
        DO D = 1, NORBUSED
            LSPN = (MOD(A + C, 2) == 0) .AND. (MOD(B + D, 2) == 0)
            IF (LSPN .AND. near_zero(UMAT2(UMatInd(A, B, C, D)))) THEN
                SUM1 = 0.0_dp
                DO L = 1, NHG
                    ID4 = GTID(L)
                    SUM1 = SUM1 + HFBASIS(D, L) * UMATT(UMatInd(ID4, A, B, C))
                END DO
                UMAT2(UMatInd(A, B, C, D)) = SUM1
!..Symmetries not needed by UMatInd
!              UMAT2(C,D,A,B)=SUM1
!              UMAT2(B,A,D,C)=SUM1
!              UMAT2(D,C,B,A)=SUM1
                IF (ABS(SUM1) > 1.0e-10_dp) WRITE(10, '(4I7,F19.9)') A, B, C, D, SUM1
            END IF
        END DO
        END DO
        END DO
        END DO
        CLOSE(10)
        WRITE(6, *) ' !!! FINISHED CALCULATING HF UMAT !!! '
        call halt_timer(proc_timer)
        RETURN
    END
    SUBROUTINE READHFBASIS(HFBASIS, HFE, G1, NBASIS)
        INTEGER NBASIS, NQNS(5), NN
        TYPE(BasisFN) G1(nBasis)
        real(dp) HFBASIS(NBASIS, NBASIS), HFE(NBASIS)
        INTEGER I, L, J, NB, NE, IG
        real(dp) VAL
        character(*), parameter :: this_routine = 'READHFBASIS'
        WRITE(6, *) "Loading HF BASIS"
        OPEN(10, FILE='HFBASIS', STATUS='OLD')
        READ(10, *)
        READ(10, *) NB, NE
!.. NE is NEVAL, and NB is NBASIS/2
!.. NBASIS is the number of orbitals, so *2 to get # spinorbitals
!         IF(NE.NE.NEL) call stop_all(this_routine, 'NEL in HFBASIS <> NEL')

        IF (NE /= NB) call stop_all(this_routine, 'NEVAL <> NBASIS in HFBASIS not supported')
        IF (NB * 2 /= NBASIS) call stop_all(this_routine, 'NBASIS in HFBASIS <> NHG')
        DO I = 1, NB
        DO L = -1, 1, 2
            READ(10, *)
            READ(10, *) HFE(I * 2 + (L - 1) / 2)
            NQNS(4) = L
            DO J = 1, NB
                READ(10, *) NN, NQNS(1), NQNS(2), NQNS(3), VAL
                IG = IFINDBASISFN(get_basisfn(NQNS), G1, NBASIS)
                HFBASIS(I * 2 + (L - 1) / 2, J * 2 + (L - 1) / 2) = VAL
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
            END DO
        END DO
        END DO
        CLOSE(10)
        RETURN
    END


#ifndef CMPLX_
    SUBROUTINE CALCHFBASIS(NBASIS, NBASISMAX, G1, BRR, ECORE, &
                           UMAT, HFE, HFBASIS, NHFIT, NEL, MS, HFMIX, EDELTA, CDELTA, TRHF, &
                           IHFMETHOD, TREADHF, FRAND, HFDET, ILOGGING)
        INTEGER NSPINS, NSBASIS
        INTEGER NBASIS, nBasisMax(5, *)
        TYPE(BasisFN) G1(nBasis)
        real(dp) UMAT(*)
        real(dp) ECORE
        INTEGER BRR(NBASIS)
        HElement_t(dp) HFBASIS(NBASIS, NBASIS)
        real(dp) HFE(NBASIS)
        HElement_t(dp), allocatable :: FMAT(:), OFMAT(:)
        HElement_t(dp), allocatable :: DMAT(:), ODMAT(:)
        HElement_t(dp), allocatable :: WORK(:)
        real(dp), allocatable :: HFES(:)
        HElement_t(dp), allocatable :: R1(:), R2(:)
        integer(TagIntType), save :: tagR1 = 0, tagR2 = 0, tagHFES = 0
        integer(TagIntType), save :: tagFMAT = 0, tagOFMat = 0
        integer(TagIntType), save :: tagWork = 0, tagDMAT = 0, tagODMat = 0
        INTEGER NHFIT, NEL
        real(dp) HFMIX, EDELTA, CDELTA
        INTEGER MS, IHFMETHOD
        LOGICAL TRHF, TREADHF
        real(dp) FRAND
        INTEGER HFDET(NEL)
        INTEGER ILOGGING
        character(*), parameter :: this_routine = 'CALCHFBASIS'
        NSPINS = 1 + (NBASISMAX(4, 2) - NBASISMAX(4, 1)) / 2
        NSBASIS = NBASIS / NSPINS
        allocate(FMAT(NSBASIS * NSBASIS * NSPINS))
        call LogMemAlloc('FMAT', NSBASIS * NSBASIS * NSPINS, HElement_t_size * 8, this_routine, tagFMAT)
        allocate(OFMAT(NSBASIS * NSBASIS * NSPINS))
        call LogMemAlloc('OFMAT', NSBASIS * NSBASIS * NSPINS, HElement_t_size * 8, this_routine, tagOFMAT)
        allocate(DMAT(NSBASIS * NSBASIS * NSPINS))
        call LogMemAlloc('DMAT', NSBASIS * NSBASIS * NSPINS, HElement_t_size * 8, this_routine, tagDMAT)
        allocate(ODMAT(NSBASIS * NSBASIS * NSPINS))
        call LogMemAlloc('ODMAT', NSBASIS * NSBASIS * NSPINS, HElement_t_size * 8, this_routine, tagODMAT)
        allocate(WORK(NSBASIS * 3))
        call LogMemAlloc('WORK', NSBASIS * 3, HElement_t_size * 8, this_routine, tagWORK)
        allocate(R1(NSBASIS * NSBASIS))
        call LogMemAlloc('R1', NSBASIS * NSBASIS, HElement_t_size * 8, this_routine, tagR1)
        allocate(R2(NSBASIS * NSBASIS))
        call LogMemAlloc('R2', NSBASIS * NSBASIS, HElement_t_size * 8, this_routine, tagR2)
        allocate(HFES(NSBASIS * NSPINS))
        call LogMemAlloc('HFES', NSBASIS * NSPINS, 8, this_routine, tagHFES)
!.. Generate initial HFBASIS vectors as the energy ordered single
!.. particle basis fns, separated into up and down blocks
        FMAT = (0.0_dp)
        DMAT = (0.0_dp)
        IF (IHFMETHOD == 0 .OR. IHFMETHOD == -1) THEN
            CALL UHFSCF(NBASIS, G1, BRR, ECORE, HFE, HFBASIS, NHFIT, NEL, MS, FMAT, &
                        DMAT, ODMAT, WORK, NSPINS, NSBASIS, HFES, OFMAT, &
                        HFMIX, EDELTA, CDELTA, TRHF, R1, R2, &
                        IHFMETHOD, TREADHF, FRAND, HFDET, ILOGGING)
        ELSE
            CALL UHFGRADDESC(NBASIS, NBASISMAX, G1, BRR, ECORE, &
                             UMAT, HFE, HFBASIS, NHFIT, NEL, MS, NSPINS, NSBASIS, HFES, &
                             HFMIX, FMAT, OFMAT, DMAT, ODMAT, EDELTA, CDELTA, R1, R2, WORK, TRHF, &
                             IHFMETHOD, TREADHF, FRAND, HFDET, ILOGGING)
        END IF

        deallocate(FMAT, OFMAT, DMAT, ODMAT, WORK, R1, R2, HFES)
        call LogMemDealloc(this_routine, tagFMAT)
        call LogMemDealloc(this_routine, tagOFMAT)
        call LogMemDealloc(this_routine, tagDMAT)
        call LogMemDealloc(this_routine, tagODMAT)
        call LogMemDealloc(this_routine, tagWORK)
        call LogMemDealloc(this_routine, tagR1)
        call LogMemDealloc(this_routine, tagR2)
        call LogMemDealloc(this_routine, tagHFES)

    END subroutine CALCHFBASIS
#endif


!.. Unrestricted HF SCF code
#ifndef CMPLX_
    SUBROUTINE UHFSCF(NBASIS, G1, BRR, ECORE, HFE, HFBASIS, NHFIT, NEL, MS, &
                      FMAT, DMAT, ODMAT, WORK, NSPINS, NSBASIS, HFES, OFMAT, HFMIX, &
                      EDELTA, CDELTA, TRHF, R1, R2, IHFMETHOD, TREADHF, &
                      FRAND, HFDET, ILOGGING)
        INTEGER NSPINS, NSBASIS
        INTEGER NBASIS
        TYPE(BasisFn) G1(*)
        real(dp) ECORE
        INTEGER BRR(NBASIS)
        HElement_t(dp) HFBASIS(NBASIS, NBASIS)
        real(dp) HFES(NSBASIS, NSPINS), HFE(NBASIS)
        HElement_t(dp) FMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) OFMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) DMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) ODMAT(NSBASIS, NSBASIS, NSPINS), WORK(NBASIS * 3)
        real(dp) HFMIX
        HElement_t(dp) R1(NSBASIS, NSBASIS), R2(NSBASIS, NSBASIS)
        INTEGER NHFIT, NEL
        INTEGER I, J, K, L, ISPN, NELS(NSPINS)
        real(dp) ELAST, EDELTA, ECUR, CDELTA
        INTEGER IHFIT
        INTEGER MS, IHFMETHOD, IRHFB
        real(dp) F
        real(dp) RMSD, FRAND
        LOGICAL BR
        LOGICAL TRHF, TREADHF
        INTEGER HFDET(*)
        INTEGER ILOGGING
        F = HFMIX
!         EDELTA=1.0e-8_dp
        ELAST = 1.D20
        WRITE(6, *) "Performing Hartree-Fock SCF diagonalisation..."
        IF (IHFMETHOD == -1) THEN
            WRITE(6, *) "Method -1:Rotational Mixing "
        ELSEIF (IHFMETHOD == 0) THEN
            WRITE(6, *) "Method 0:Linear Mixing"
        END IF
        WRITE(6, *) "HF Mixing", HFMIX
        WRITE(6, *) "E Thresh:", EDELTA
        WRITE(6, *) "RMSD Thresh:", CDELTA
        IF (NSPINS == 2) THEN
            NELS(2) = (MS + NEL) / 2
            NELS(1) = NEL - NELS(2)
            WRITE(6, *) " Beta, Alpha: ", NELS(1), NELS(2)
        ELSE
            NELS(1) = NEL
        END IF
        IF (TREADHF) THEN
            CALL READHFFMAT(NBASIS, FMAT, HFES, G1, NSPINS, NSBASIS, .FALSE.)
        ELSE
            CALL GENHFGUESS(FMAT, NSPINS, NSBASIS, BRR, G1, .FALSE., MS, FRAND, NELS, HFDET)
        END IF
        WRITE(6, *) "Iteration   Energy     MSD"
        BR = .TRUE.
        IHFIT = 1
        IRHFB = 0
        IF (TRHF .AND. NSPINS > 1) THEN
        IF (NELS(2) > NELS(1)) THEN
            IRHFB = 2
        ELSE
            IRHFB = 1
        END IF
        END IF
        DO WHILE (BR)
        IF (IRHFB > 0) THEN
            CALL DCOPY(NSBASIS * NSBASIS * HElement_t_size, FMAT(1, 1, IRHFB), 1, FMAT(1, 1, 3 - IRHFB), 1)
            CALL DCOPY(NSBASIS * NSBASIS * HElement_t_size, DMAT(1, 1, IRHFB), 1, DMAT(1, 1, 3 - IRHFB), 1)
        END IF
        CALL DCOPY(NSBASIS * NSBASIS * NSPINS * HElement_t_size, DMAT, 1, ODMAT, 1)
        CALL DCOPY(NSBASIS * NSBASIS * NSPINS * HElement_t_size, FMAT, 1, OFMAT, 1)
!.. Construct the Density Matrix
        CALL GENDMAT(NSPINS, NSBASIS, NELS, FMAT, DMAT, .FALSE.)
!.. See how much our density matrix has changed from last time
        RMSD = 0.0_dp
        DO ISPN = 1, NSPINS
        DO I = 1, NSBASIS
        DO J = 1, NSBASIS
            IF (.NOT. TRHF .OR. TRHF .AND. ISPN == IRHFB) RMSD = RMSD + abs(DMAT(I, J, ISPN) - ODMAT(I, J, ISPN))**2
        END DO
        END DO
        END DO
        RMSD = SQRT(RMSD / (NSBASIS * NSBASIS * NSPINS))
!.. replace our HF orbitals in FMAT with the Fock matrix
!.. FMAT just stores the results and the HF orbitals (currently) in it
!.. are not used in calculating the F matrix
        CALL GENFMAT(FMAT, DMAT, NSBASIS, NSPINS)
        CALL DIAGFMAT(NSPINS, NSBASIS, NELS, FMAT, DMAT, HFES, WORK, ECORE, ECUR)
        IF (IRHFB > 0) THEN
            CALL DCOPY(NSBASIS * NSBASIS, FMAT(1, 1, IRHFB), 1, FMAT(1, 1, 3 - IRHFB), 1)
            CALL DCOPY(NSBASIS, HFES(1, IRHFB), 1, HFES(1, 3 - IRHFB), 1)
        END IF
!.. FMAT now contains HF orbitals again, and ECUR the Fock Energy
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector
        WRITE(6, *) IHFIT, ECUR, RMSD
!.. Now add back in some of our original F matrix
        IF (IHFMETHOD == -1) THEN
            CALL HFROTMIX(FMAT, OFMAT, NSPINS, NSBASIS, F, R1, R2, WORK)
        ELSE
            CALL HFLINMIX(FMAT, OFMAT, NSPINS, NSBASIS, F, R1, R2, WORK)
        END IF
        IHFIT = IHFIT + 1
        IF (IHFIT > NHFIT) THEN
            WRITE(6, *) "*** WARNING Hartree-Fock did not converge ***"
            BR = .FALSE.
        END IF
        IF (ABS(ECUR - ELAST) <= EDELTA .AND. RMSD <= CDELTA .AND. IHFIT > 5) THEN
            WRITE(6, *) "*** Hartree-Fock converged in ", IHFIT, " iterations."
            WRITE(6, *) "*** HF ENERGY=", ECUR
            BR = .FALSE.
        END IF
        ELAST = ECUR

        END DO
        IF (BTEST(ILOGGING, 11)) then
            CALL WRITEHFPSIALL(NBASIS, FMAT, HFES, G1, NSPINS, NSBASIS, .FALSE.)
        end if
!.. We write out HFMAT
        HFBASIS = (0.0_dp)
        DO I = 1, NSBASIS
        DO ISPN = 1, NSPINS
            K = (I - 1) * NSPINS + ISPN
            DO J = 1, NSBASIS
                L = (J - 1) * NSPINS + ISPN
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
                HFBASIS(K, L) = FMAT(J, I, ISPN)
            END DO
            HFE(K) = HFES(I, ISPN)
        END DO
        END DO
    END subroutine UHFSCF
#endif

#ifndef CMPLX_
    SUBROUTINE GENDMAT(NSPINS, NSBASIS, NELS, FMAT, DMAT, LTRANS)
        INTEGER NSPINS, NSBASIS, NELS(NSPINS)
        HElement_t(dp) FMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) DMAT(NSBASIS, NSBASIS, NSPINS)
        LOGICAL LTRANS

        INTEGER ISPN, I, J, K
        HElement_t(dp) TOT
!.. Construct the Density Matrix
        DO ISPN = 1, NSPINS
        DO I = 1, NSBASIS
        DO J = I, NSBASIS
            TOT = 0.0_dp
!.. Sum over occupied orbitals for our spin
            DO K = 1, NELS(ISPN)
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector
#ifdef CMPLX_
                IF (LTRANS) THEN
                    TOT = TOT + CONJG(FMAT(K, I, ISPN)) * FMAT(K, J, ISPN)
                ELSE
                    TOT = TOT + FMAT(I, K, ISPN) * CONJG(FMAT(J, K, ISPN))
                END IF
#else
                IF (LTRANS) THEN
                    TOT = TOT + (FMAT(K, I, ISPN)) * FMAT(K, J, ISPN)
                ELSE
                    TOT = TOT + FMAT(I, K, ISPN) * (FMAT(J, K, ISPN))
                END IF
#endif
            END DO
            DMAT(I, J, ISPN) = TOT
            DMAT(J, I, ISPN) = TOT
        END DO
        END DO
        IF (HFLogLevel > 0) CALL WRITE_HEMATRIX("DMAT", NSBASIS, NSBASIS, DMAT(1, 1, ISPN))
        END DO
        RETURN
    END subroutine
#endif

#ifndef CMPLX_
    SUBROUTINE GENFMAT(FMAT, DMAT, NSBASIS, NSPINS)
        INTEGER NSBASIS, NSPINS
        HElement_t(dp) FMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) DMAT(NSBASIS, NSBASIS, NSPINS)
!         real(dp) TMAT(NSBASIS*NSPINS,NSBASIS*NSPINS)
        HElement_t(dp) TOT
        INTEGER I, J, K, L, ISPN, ISPN2, ID1, ID2, ID3, ID4
        real(dp) RHFMUL
!         IF(TRHF) THEN
!            RHFMUL=2.0_dp
!         ELSE
        RHFMUL = 1.0_dp
!         ENDIF
!.. Construct the Fock Matrix
        FMAT = (0.0_dp)
!. <ui|F|uj>=tij+Sum_kl(Dkl <ui uk|U|uj ul>-s<ui uk|U|ul uj>) where s=1 if i,k same spin
        DO ISPN = 1, NSPINS
        DO I = 1, NSBASIS
        DO J = I, NSBASIS
            TOT = GetTMATEl((I - 1) * NSPINS + ISPN,(J - 1) * NSPINS + ISPN)
!                  WRITE(6,*) I,J,TOT
!.. Now sum in the alpha and beta block of u matrix
            ID1 = GTID((I - 1) * NSPINS + ISPN)
            ID2 = GTID((J - 1) * NSPINS + ISPN)
            DO ISPN2 = 1, NSPINS
            DO K = 1, NSBASIS
            DO L = 1, NSBASIS
                ID3 = GTID((K - 1) * NSPINS + ISPN2)
                ID4 = GTID((L - 1) * NSPINS + ISPN2)
                IF (abs(DMAT(K, L, ISPN2)) > DmatEpsilon) then
                    TOT = TOT + (RHFMUL) * DMAT(K, L, ISPN2) * get_umat_el(ID1, ID3, ID2, ID4)
                end if
                IF (ISPN2 == ISPN .AND. (abs(DMAT(K, L, ISPN2)) > DmatEpsilon)) then
                    TOT = TOT - DMAT(K, L, ISPN2) * get_umat_el(ID1, ID3, ID4, ID2)
                end if
            END DO
            END DO
            END DO
            FMAT(I, J, ISPN) = TOT
#ifdef CMPLX_
            FMAT(J, I, ISPN) = CONJG(TOT)
#endif
        END DO
        END DO
        IF (HFLogLevel > 0) CALL WRITE_HEMATRIX("FMAT", NSBASIS, NSBASIS, FMAT(1, 1, ISPN))
        END DO
    END subroutine
#endif

#ifndef CMPLX_
    SUBROUTINE DIAGFMAT(NSPINS, NSBASIS, NELS, FMAT, DMAT, HFES, WORK, ECORE, ECUR)
        INTEGER NSPINS, NSBASIS, NELS(NSPINS)
        real(dp) HFES(NSBASIS, NSPINS)
        HElement_t(dp) FMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) DMAT(NSBASIS, NSBASIS, NSPINS)
!         real(dp) TMAT(NSBASIS*NSPINS,NSBASIS*NSPINS)
        real(dp) ECUR, ECORE
        HElement_t(dp) WORK(3 * NSBASIS), RWORK(3 * NSBASIS)

        INTEGER(int32) :: INFO
        HElement_t(dp) TOT, TOT2
        INTEGER ISPN, I, J, K
        character(*), parameter :: this_routine = 'DIAGFMAT'

!.. First calculate the HF energy double counting contrib
        TOT = 0.0_dp
        DO ISPN = 1, NSPINS
        DO J = 1, NSBASIS
        DO K = 1, NSBASIS
            TOT = TOT + DMAT(J, K, ISPN) * (FMAT(J, K, ISPN) - GetTMATEl((J - 1) * NSPINS + ISPN,(K - 1) * NSPINS + ISPN))
        END DO
        END DO
        END DO

!.. Now diagonalize the Fock matrix
        DO ISPN = 1, NSPINS
        IF (HElement_t_size == 1) THEN
            CALL DSYEV('V', 'U', NSBASIS, FMAT(1, 1, ISPN), NSBASIS, HFES(1, ISPN), WORK, 3 * NSBASIS, INFO)
        ELSE
            CALL ZHEEV('V', 'U', NSBASIS, FMAT(1, 1, ISPN), NSBASIS, HFES(1, ISPN), WORK, 3 * NSBASIS, RWORK, INFO)
        END IF
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector
        IF (INFO /= 0) THEN
            WRITE(6, *) 'DYSEV error: ', INFO
            call stop_all(this_routine, "DYSEV error")
        END IF
        IF (HFLogLevel > 0) CALL WRITE_HEMATRIX("EIGVEC", NSBASIS, NSBASIS, FMAT(1, 1, ISPN))
        IF (HFLogLevel > 0) CALL NECI_WRITE_MATRIX("EIGVAL", 1, NSBASIS, HFES(1, ISPN))
        END DO
!.. now calculate the sum of the occupied Fock orbitals
        TOT2 = 0.0_dp
        DO ISPN = 1, NSPINS
        DO I = 1, NELS(ISPN)
            TOT2 = TOT2 + (HFES(I, ISPN))
        END DO
        END DO
!.. subtract out the double counting term, and add in the core energy
        ECUR = (TOT2) - (TOT) / 2.0_dp + ECORE
    END subroutine diagfmat
#endif



#ifndef CMPLX_
    SUBROUTINE WRITEHFPSIALL(NBASIS, FMAT, HFES, G1, NSPINS, NSBASIS, TRANSP)
        INTEGER NBASIS, NSBASIS, NSPINS, I, J, K
        TYPE(BasisFN) G1(nBasis)
        real(dp) FMAT(NSBASIS, NSBASIS, NSPINS), HFES(NSBASIS, NSPINS)
        LOGICAL TRANSP
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector unless transposed
        OPEN(10, FILE='HFBASIS', STATUS='UNKNOWN')
        WRITE(10, *) 'NBASIS,NEVAL'
        WRITE(10, *) NSBASIS, NSBASIS
        DO K = 1, NSBASIS
        IF (NSPINS > 1) THEN
            WRITE(10, *) ' BETA ELECTRON NO: ', K
            WRITE(10, *) HFES(K, 1)
            DO I = 1, NSBASIS
            IF (TRANSP) THEN
                WRITE(10, '(I7,1X,3I7,F19.9)') I,(G1((I - 1) * NSPINS + 1)%k(J), J=1, 3), FMAT(K, I, 1)
            ELSE
                WRITE(10, '(I7,1X,3I7,F19.9)') I,(G1((I - 1) * NSPINS + 1)%k(J), J=1, 3), FMAT(I, K, 1)
            END IF
            END DO
        END IF
        WRITE(10, *) 'ALPHA ELECTRON NO: ', K
        WRITE(10, *) HFES(K, 2)
        DO I = 1, NSBASIS
        IF (TRANSP) THEN
            WRITE(10, '(I7,1X,3I7,F19.9)') I,(G1((I - 1) * NSPINS + 2)%k(J), J=1, 3), FMAT(K, I, 2)
        ELSE
            WRITE(10, '(I7,1X,3I7,F19.9)') I,(G1((I - 1) * NSPINS + 2)%k(J), J=1, 3), FMAT(I, K, 2)
        END IF
        END DO
        END DO
        CLOSE(10)
    END subroutine
#endif

#ifndef CMPLX_
!.. Generate initial density matrix, as well as a guess at the HF DET
    SUBROUTINE GENHFGUESS(FMAT, NSPINS, NSBASIS, BRR, G1, TRANS, LMS, FRAND, NELS, HFDET)
        INTEGER NSPINS, NSBASIS
        HElement_t(dp) FMAT(NSBASIS, NSBASIS, NSPINS)
        real(dp) PI, R
        INTEGER ISPN, I, J
        INTEGER BRR(NSBASIS * NSPINS), NELR, IREAL, IS
        TYPE(BasisFN) G1(*)
        INTEGER LMS
!.. Working space
        real(dp) FRAND
        LOGICAL TRANS
        INTEGER NELS(NSPINS)
        INTEGER HFDET(*), NEL
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector unless transposed
        PI = 3.141592653589793_dp
        R = genrand_real2_dSFMT()
        FMAT = (0.0_dp)
        WRITE(6, *) "Generating HF Guess..."
        NEL = 0
        DO IS = 1, NSPINS
        IF (LMS < 0) THEN
            ISPN = IS
        ELSE
            ISPN = NSPINS + 1 - IS
        END IF
        IREAL = 1
        WRITE(6, "(A,I2,A)", advance='no') "Spin ", IS, ":"

        DO I = 1, NSBASIS
        DO WHILE (G1(BRR(IREAL))%Ms /= (-3 + 2 * ISPN))
            IREAL = IREAL + 1
        END DO
        NELR = (BRR(IREAL) - 1) / NSPINS + 1
        IF (TRANS) THEN
            FMAT(I, NELR, ISPN) = (1.0_dp)
        ELSE
            FMAT(NELR, I, ISPN) = (1.0_dp)
        END IF
        IF (I <= NELS(IS)) THEN
            WRITE(6, "(I4,A)", advance='no') BRR(IREAL), ","
            NEL = NEL + 1
            HFDET(NEL) = BRR(IREAL)
        END IF
        IREAL = IREAL + 1
        END DO
        WRITE(6, *)
        DO I = 1, NSBASIS
        DO J = 1, NSBASIS
            FMAT(I, J, ISPN) = FMAT(I, J, ISPN) + (FRAND * genrand_real2_dSFMT())
        END DO
        END DO
        CALL GRAMSCHMIDT_NECI(FMAT(1, 1, ISPN), NSBASIS)
!            CALL LOWDIN_ORTH(FMAT(1,1,ISPN),NSBASIS,R1,R2,WORK)
        END DO
        call sort(HFDet(1:nel))
        call write_det(6, HFDET, .true.)
    END subroutine
#endif


#ifndef CMPLX_
    SUBROUTINE HFLINMIX(FMAT, OFMAT, NSPINS, N, FMIX, R1, R2, WORK)
        INTEGER N, NSPINS
        real(dp) FMIX
        HElement_t(dp) FMAT(N, N, NSPINS), OFMAT(N, N, NSPINS)
        HElement_t(dp) R1(*), R2(*), WORK(*)
        INTEGER I, J, ISPN
        DO ISPN = 1, NSPINS
        DO I = 1, N
        DO J = 1, N
            FMAT(I, J, ISPN) = (FMIX) * FMAT(I, J, ISPN) + (1.0_dp - FMIX) * OFMAT(I, J, ISPN)
        END DO
        END DO
        IF (N < 3) THEN
!.. For N<2, Lowdin will return in FMAT exactly what we started with in OFMAT, which is
!.. rather pointless in mixing, so we use Gram-Schmidt to mix things up a
!.. bit.
            CALL GRAMSCHMIDT_NECI(FMAT(1, 1, ISPN), N)
        ELSE
            CALL LOWDIN_ORTH(FMAT(1, 1, ISPN), N, R1, R2, WORK)
        END IF
        END DO
    END subroutine
#endif

#ifndef CMPLX_
    SUBROUTINE HFROTMIX(FMAT, OFMAT, NSPINS, N, FMIX, R1, R2, WORK)
        INTEGER NSPINS, N
        real(dp) FMIX
        HElement_t(dp) :: FMAT(N, N, NSPINS)
        real(dp) OFMAT(N, N, NSPINS)
        HElement_t(dp) R1(N, N), R2(N, N), WORK(3 * N)
        INTEGER I, J, ISPN
!.. OFMAT is the old HF orbitals
!.. FMAT is the new HF orbitals
!.. as HF orbitals make an orthoganal set, then FMAT will merely be a
!.. rotation of OFMAT around some axis
!.. If we want to only include part of the new HF orbitals, we should
!.. only rotate OFMAT by part of that amount.

!.. F' = OFMAT.  F=FMAT.  R is the rotation
!.. F=R F'.  We want F'' (the new FMAT we're to generate) to be
!.. F'' = R^a F'  (where a=FMIX is between 0 and 1)
!.. R^a = (I +(R-I))^a = I+a(R-I)+a(a-1)/2 (R-I)^2 + ... (Taylor)

!.. R = F F'T as F and F' are orthogonal matrices
        DO ISPN = 1, NSPINS

!.. Work out R = R1=1.0_dp * F * F'T + 0.0_dp*R1
            CALL DGEMM('N', 'T', N, N, N, 1.0_dp, FMAT(1, 1, ISPN), N, OFMAT(1, 1, ISPN), N, 0.0_dp, R1, N)
!.. now let P=R1=R-I
            DO I = 1, N
                R1(I, I) = R1(I, I) - 1.0_dp
            END DO
!.. R^a = I+aP(I+(a-1)/2 P(I+ (a-2)/3 P (I+...) ) ) F'
!.. Let FMAT be the accumulator
            CALL DCOPY(N * N, OFMAT(1, 1, ISPN), 1, FMAT(1, 1, ISPN), 1)
!.. Go up to 2nd order taylor
            DO J = 1, 0, -1
!.. Set R2=I
                R2 = 0.0_dp
                DO I = 1, N
                    R2(I, I) = 1.0_dp
                END DO
!.. Work out R2=(a-J)/(J+1.0_dp) * R1 * FMAT + 1.0_dp*R2.  a=FMIX
                CALL DGEMM('N', 'N', N, N, N,(FMIX - J) / (J + 1.0_dp), R1, N, FMAT(1, 1, ISPN), N, 1.0_dp, R2, N)
                CALL DCOPY(N * N, R2, 1, FMAT(1, 1, ISPN), 1)
            END DO
!.. Now reorthoganalise, as this is just an approximation
            CALL LOWDIN_ORTH(FMAT(1, 1, ISPN), N, R1, R2, WORK)
        END DO
    END subroutine
#endif

!.. Unrestricted Hartree Fock with a gradient descent method.
#ifndef CMPLX_
    SUBROUTINE UHFGRADDESC(NBASIS, NBASISMAX, G1, BRR, ECORE, UMAT, HFE, HFBASIS, &
                           NHFIT, NEL, MS, NSPINS, NSBASIS, HFES, HFMIX, CMAT, OCMAT, DEDCIJ, &
                           DMAT, EDELTA, CDELTA, R1, R2, WORK, TRHF, IHFMETHOD, &
                           TREADHF, FRAND, HFDET, ILOGGING)
        TYPE(BasisFn) G1(*)
        INTEGER NSPINS, NSBASIS
        INTEGER NBASIS, nBasisMax(5, *)
        real(dp) UMAT(*)
        real(dp) ECORE
        INTEGER BRR(NBASIS)
        real(dp) HFBASIS(NBASIS, NBASIS), HFE(NBASIS)
        real(dp) HFES(NSBASIS, NSPINS)
        HElement_t(dp) CMAT(NSBASIS, NSBASIS, NSPINS)
        real(dp) OCMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) DEDCIJ(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) DMAT(NSBASIS, NSBASIS, NSPINS)
        HElement_t(dp) WORK(NBASIS * 3)
        INTEGER INORDER(100, 2), ILOGGING
        real(dp) EORDER(100, 2)
!,HFMIX
        HElement_t(dp) R1(NSBASIS, NSBASIS), R2(NSBASIS, NSBASIS)
        INTEGER NHFIT, NEL, IHFMETHOD, NELEX2
!,NSTART(NEL)
        INTEGER I, J, K, L, ISPN, NELS(NSPINS)
        real(dp) ELAST, HFMIX, ECUR, FRAND
!         INTEGER INDS(NBASIS)
!         real(dp) TOT,ELAST,EDELTA,ECUR,TOT2,CDELTA
!         INTEGER ID1,ID2,ID3,ID4,INFO,
        INTEGER IHFIT
        INTEGER MS
!         real(dp) F
!         real(dp) RMSD
        LOGICAL BR
        LOGICAL TRHF, TREADHF
        real(dp) TOT, MIX, TOT2
        INTEGER NDET1(0:NEL + 1), NSPN(NSPINS), NELEX
        INTEGER IRHFB, JSPN
        real(dp) RMSD, EDELTA, CDELTA, EN, ECUR2
        INTEGER HFDET(*)

        character(*), parameter :: this_routine = 'UHFGRADDESC'
        irhfb = 0
        IF (NSBASIS > 99) call stop_all(this_routine, 'ERROR - hardcoded NSBASIS limit of 100')
        ELAST = 1.D20
        WRITE(6, *) "Performing Hartree-Fock Gradient Descent..."
        IF (IHFMETHOD == 1) THEN
            WRITE(6, *) "Method 1:Singles replacement "
        ELSEIF (IHFMETHOD == 2) THEN
            WRITE(6, *) "Method 2:Explicit differential"
        END IF
        IF (NSPINS == 2) THEN
            NELS(2) = (MS + NEL) / 2
            NELS(1) = NEL - NELS(2)
            WRITE(6, *) " Beta, Alpha: ", NELS(1), NELS(2)
        ELSE
            NELS(1) = NEL
        END IF
!.. Cij - i corresponds to rows and new basis functions, phi_i
!..       j corresponds to columns and old basis functions, u_j
!.. phi_i=sum_j=1,M cij u_j

!         VMAT=0.0_dp
        IF (TREADHF) THEN
            CALL READHFFMAT(NBASIS, CMAT, HFES, G1, NSPINS, NSBASIS, .TRUE.)
        ELSE
            CALL GENHFGUESS(CMAT, NSPINS, NSBASIS, BRR, G1, .TRUE., MS, FRAND, NELS, HFDET)
        END IF
!         DO ISPN=1,NSPINS
!         DO I=1,NSBASIS
!            WRITE(6,*) (FMAT(I,J,ISPN),J=1,NSBASIS)
!         ENDDO
!         ENDDO

!.. Initialize our HF det in NDET1
        DO I = 1, NSPINS
            NSPN(I) = 0
        END DO
        I = 1
        NDET1(0) = 0
        NDET1(NEL + 1) = NBASIS + 1
        DO WHILE (I <= NEL)
        DO ISPN = 1, NSPINS
        IF (NSPN(ISPN) < NELS(ISPN)) THEN
            NDET1(I) = NSPN(ISPN) * NSPINS + ISPN
            NSPN(ISPN) = NSPN(ISPN) + 1
            I = I + 1
        END IF
        END DO
        END DO
        call write_det(6, NDET1(1), .true.)
        BR = .TRUE.
        WRITE(6, *) "Iteration   Energy     MSD   Fock Energy"
        IHFIT = 0
        IF (TRHF .AND. NSPINS > 1) THEN
        IF (NELS(2) > NELS(1)) THEN
            IRHFB = 2
        ELSE
            IRHFB = 1
        END IF
        END IF
        DO WHILE (BR)
        IF (IRHFB > 0) THEN
            CALL DCOPY(NSBASIS * NSBASIS, CMAT(1, 1, IRHFB), 1, CMAT(1, 1, 3 - IRHFB), 1)
        END IF
        CALL DCOPY(NSBASIS * NSBASIS * NSPINS, CMAT, 1, OCMAT, 1)
        IHFIT = IHFIT + 1
!.. First Calculate dE/dcij
!.. dE/dcij is automatically 0 if i>N as the HF det only depends on
!.. phi_1 to phi_N.  All values of j must be iterated as each phi_i is
!.. dependent on all u_j
        ECUR = GETHELEMENT2T(NDET1(1), NDET1(1), NEL, NBASISMAX, NBASIS, ECORE, 0, CMAT, NSBASIS, NSPINS)
!.. Calculate the Gradient
        IF (IHFMETHOD == 1) THEN
            CALL CALCDEDCIJ(CMAT, DEDCIJ, NDET1, NSPINS, NSBASIS, ECORE, NBASIS, NBASISMAX, NELS, NEL, ECUR)
        ELSEIF (IHFMETHOD == 2) THEN
            CALL CALCDEDCIJ2(CMAT, DEDCIJ, NDET1, NSPINS, NSBASIS, UMAT, NELS, NEL)
        END IF
!.. DEDCIJ now comtains all elements of dE/dcij
!.. To move down the slope, we subtract a small amount of this from cij,
!.. and re-orthogonalise
!.. HFMIX is -ve
!              WRITE(6,*) ((CMAT(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
!              WRITE(6,*)
!              WRITE(6,*) ((DEDCIJ(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
!.. modify the velocity.
!              R1 = VMAT(:,:,ISPN) + MIX*DEDCIJ(:,:,ISPN)
!              CALL DCOPY(NSBASIS*NSBASIS,R1,1,VMAT(1,1,ISPN),1)
!              R1 = CMAT(:,:,ISPN) + 0.1_dp*VMAT(:,:,ISPN)

!.. Remove the projection of the "force" already in the direction of
!.. the coefficients
        DO ISPN = 1, NSPINS
            R1 = 0.0_dp
            DO I = 1, NELS(ISPN)
                TOT = 0.0_dp
                DO J = 1, NSBASIS
                    TOT = TOT + DEDCIJ(I, J, ISPN) * CMAT(I, J, ISPN)
                END DO
                TOT2 = 0.0_dp
                DO J = 1, NSBASIS
                    R1(I, J) = DEDCIJ(I, J, ISPN) - TOT * CMAT(I, J, ISPN)
                    TOT2 = TOT2 + R1(I, J)**2
                END DO
                TOT2 = SQRT(TOT2)
                DO J = 1, NSBASIS
!                  R1(I,J)=R1(I,J)/TOT2
                END DO
            END DO
            MIX = -HFMIX
!/ABS(ECUR)
!              IF(ABS(ECUR).LT.1.0e-4_dp) MIX=HFMIX
            IF (IHFIT > NHFIT) BR = .FALSE.
            R2 = CMAT(:, :, ISPN) + R1
            CALL DCOPY(NSBASIS * NSBASIS, R2, 1, CMAT(1, 1, ISPN), 1)
!              WRITE(6,*)
!              WRITE(6,*) ((CMAT(I,J,ISPN),I=1,NSBASIS),J=1,NSBASIS)
!              WRITE(6,*)
!              WRITE(6,*)
            CALL LOWDIN_ORTH(CMAT(1, 1, ISPN), NSBASIS, R1, R2, WORK)
        END DO
        RMSD = 0.0_dp
        DO ISPN = 1, NSPINS
            DO I = 1, NSBASIS
                DO J = 1, NSBASIS
                    IF (.NOT. TRHF .OR. (TRHF .AND. ISPN == IRHFB)) then
                        RMSD = RMSD + (CMAT(I, J, ISPN) - OCMAT(I, J, ISPN))**2
                    end if
                END DO
            END DO
        END DO
        RMSD = SQRT(RMSD / (NSBASIS * NSBASIS * NSPINS))
        IF (IHFIT > NHFIT) THEN
            WRITE(6, *) "** WARNING Hartree-Fock did not converge **"
            BR = .FALSE.
        END IF
        IF (ABS(ECUR - ELAST) < EDELTA .AND. RMSD < CDELTA .AND. IHFIT > 5) THEN
            WRITE(6, *) "*** Hartree-Fock converged in ", IHFIT, " iterations."
            WRITE(6, *) "*** HF ENERGY=", ECUR
            BR = .FALSE.
        END IF
        ELAST = ECUR

!.. Construct the Density Matrix
        CALL GENDMAT(NSPINS, NSBASIS, NELS, CMAT, DMAT, .TRUE.)
!.. Use the Density Matrix to generate the Fock matrix (in DEDCIJ)
        CALL GENFMAT(DEDCIJ, DMAT, NSBASIS, NSPINS)
        CALL DIAGFMAT(NSPINS, NSBASIS, NELS, DEDCIJ, DMAT, HFES, WORK, ECORE, ECUR2)

        WRITE(6, "(I6)", advance='no') IHFIT
        WRITE(6, *) ECUR, RMSD, ECUR2
!.. DEDCIJ now contains HF orbitals, and ECUR the Fock Energy
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector

!.. calculate the orbital energies every time
        DO ISPN = 1, NSPINS
        DO I = 1, NSBASIS
            NELEX = (I - 1) * NSPINS + 1 + ISPN - 1
            EN = 0.0_dp
#ifdef CMPLX_
            call stop_all(this_routine, "not implemented for complex")
#else
            CALL GETTRTMATEL(NELEX, NELEX, CMAT, NSBASIS, NSPINS, EN)
#endif
            HFES(I, ISPN) = EN
            DO JSPN = 1, NSPINS
            DO J = 1, NELS(JSPN)
                NELEX2 = (J - 1) * NSPINS + 1 + JSPN - 1
                IF (NELEX /= NELEX2) THEN
!.. we're not allowed to count the current electron again
#ifdef CMPLX_
                    call stop_all(this_routine, "not implemented for complex")
#else
                    CALL GETTRUMATEL(NELEX, NELEX2, NELEX, NELEX2, CMAT, NSBASIS, NSPINS, EN)
#endif

!                        HFES(I,ISPN)=HFES(I,ISPN)+EN
                    EN = 0.0_dp
                    IF (ISPN == JSPN) then
#ifdef CMPLX_
                        call stop_all(this_routine, "not implemented for complex")
#else
                        CALL GETTRUMATEL(NELEX, NELEX2, NELEX2, NELEX, CMAT, NSBASIS, NSPINS, EN)
#endif
                    endif
!                        HFES(I,ISPN)=HFES(I,ISPN)-EN
                END IF
            END DO
            END DO
        END DO
        DO I = 1, NSBASIS
            INORDER(I, ISPN) = I
        END DO
        CALL DCOPY(NSBASIS, HFES(1, ISPN), 1, EORDER(1, ISPN), 1)
        call sort(eorder(1:nsBasis, iSpn), inOrder(1:nsBasis, iSpn))
        END DO
        IF (BR) THEN
!.. Now re-order the orbitals
!.. although we don't actually use this info
            DO ISPN = 1, NSPINS
                R1 = 0.0_dp
                DO I = 1, NSBASIS
!.. If R1(2,1) is occupied, then postmultiplying C by R1 moves
!.. column 2 in C to column 1.
!.. INORDER(1,ISPN) is the old index of the lowest energy orb
                    R1(I, INORDER(I, ISPN)) = 1
                END DO
!.. Work out R2=1.0_dp CMAT*R1+ 0.0_dp*R2
                CALL DGEMM('N', 'N', NSBASIS, NSBASIS, NSBASIS, 1.0_dp, R1, NSBASIS, CMAT(1, 1, ISPN), NSBASIS, 0.0_dp, R2, NSBASIS)
!            CALL DCOPY(NSBASIS*NSBASIS,R2,1,CMAT(1,1,ISPN),1)
                CALL DGEMM('N', 'N', NSBASIS, NSBASIS, NSBASIS, 1.0_dp, R1, NSBASIS, OCMAT(1, 1, ISPN), NSBASIS, 0.0_dp, R2, NSBASIS)
!            CALL DCOPY(NSBASIS*NSBASIS,R2,1,OCMAT(1,1,ISPN),1)
            END DO

            DO I = 1, NSBASIS
            DO ISPN = 1, NSPINS
!                  WRITE(6,*) I,ISPN*2-3,HFES(I,ISPN)
!                  WRITE(6,*) I,ISPN*2-3,HFES(INORDER(I,ISPN),ISPN)
!,
!     &               INORDER(I,ISPN)
            END DO
            END DO
            DO ISPN = 1, NSPINS
            DO I = 1, NSBASIS
!                  WRITE(6,*) I,ISPN*2-3,EORDER(I,ISPN),INORDER(I,ISPN)
!HFES(INORDER(I,ISPN),ISPN),
!     &               INORDER(I,ISPN)
            END DO
            END DO
        END IF
        IF (MOD(IHFIT, 10) == 0 .AND. BTEST(ILOGGING, 11)) then
            CALL WRITEHFPSIALL(NBASIS, CMAT, HFES, G1, NSPINS, NSBASIS, .TRUE.)
        end if
        END DO
        IF (BTEST(ILOGGING, 11)) then
            CALL WRITEHFPSIALL(NBASIS, CMAT, HFES, G1, NSPINS, NSBASIS, .TRUE.)
        end if
!.. We write out HFMAT
        HFBASIS = 0.0_dp
        DO I = 1, NSBASIS
        DO ISPN = 1, NSPINS
            K = (I - 1) * NSPINS + ISPN
            DO J = 1, NSBASIS
                L = (J - 1) * NSPINS + ISPN
!.. eigenvector N is in FMAT(i,N,ISPN), where i is the component of the
!.. vector
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
                HFBASIS(K, L) = CMAT(I, J, ISPN)
            END DO
            HFE(K) = HFES(I, ISPN)
        END DO
        END DO
        RETURN
    END subroutine UHFGRADDESC
#endif

#ifndef CMPLX_
    SUBROUTINE CALCDEDCIJ(CMAT, DEDCIJ, NDET1, NSPINS, NSBASIS, ECORE, NBASIS, NBASISMAX, NELS, NEL, ECUR)
        INTEGER NSBASIS, NSPINS, NBASIS, NEL
        real(dp) CMAT(NSBASIS, NSBASIS, NSPINS)
        real(dp) DEDCIJ(NSBASIS, NSBASIS, NSPINS)
        INTEGER NDET1(0:NEL + 1), NDET2(NEL), NELS(NSPINS)
        real(dp) ECORE
        INTEGER NBASISMAX(*)

        INTEGER I, J, K, L, M, ISPN
        INTEGER NELEX, NELNEW, IPOSO, IPOSN, ISGNCH
        real(dp) ECUR, SCRRES, TOT
        DEDCIJ = 0.0_dp
        iposo = 0
        iposn = 0
        DO ISPN = 1, NSPINS
        DO I = 1, NELS(ISPN)
!.. First find out which new basis function index this is
            NELEX = (I - 1) * NSPINS + 1 + ISPN - 1
            DO J = 1, NSBASIS
!.. first part is d_ji <Psi|H|Psi>.  d=cT
                TOT = ECUR * CMAT(I, J, ISPN)
!.. now add in factor for all single replacements of NELEX
!.. We only need to look at the same spin, as we multiply by dki, which
!.. is zero if i and k have different spins
                DO K = NELS(ISPN) + 1, NSBASIS
                    NELNEW = (K - 1) * NSPINS + 1 + ISPN - 1
                    M = 1
                    DO L = 1, NEL
                    IF (NDET1(L) == NELEX) THEN
                        IPOSO = L
                    ELSEIF (NDET1(L - 1) < NELNEW .AND. NDET1(L) > NELNEW) THEN
!.. We slot in the new det here
                        NDET2(M) = NELNEW
                        IF (M < NEL) NDET2(M + 1) = NDET1(L)
                        IPOSN = M
                        M = M + 2
                    ELSE
                        NDET2(M) = NDET1(L)
                        M = M + 1
                    END IF
                    END DO
!.. If we haven't yet put it the new electron, we put it at the end
                    IF (M == NEL) NDET2(NEL) = NELNEW

!.. calculate how many positions move from Old to new electron, and work
!.. out appropriate sign change
                    ISGNCH = (-1)**(IPOSN - IPOSO + NEL)
!.. At this point, NDET2 contains a det which is Psi with  phi_i
!.. replaced by phi_k and rerdered.
!.. Now calculate <NDET2 | H | PSI>
                    SCRRES = GETHELEMENT2T(NDET1(1), NDET2, NEL, NBASISMAX, NBASIS, ECORE, 1, CMAT, NSBASIS, NSPINS)
!.. We multiply by the amount of phi_k in u_j, as well as the sign
                    TOT = TOT + ISGNCH * SCRRES * CMAT(K, J, ISPN)
                END DO
!.. TOT now contains dE/dcij, so we store this
                DEDCIJ(I, J, ISPN) = 2 * TOT
            END DO
        END DO
        END DO
    END subroutine
#endif

#ifndef CMPLX_
    SUBROUTINE CALCDEDCIJ2(CMAT, DEDCIJ, NDET1, NSPINS, NSBASIS, UMAT, NELS, NEL)
        INTEGER NSBASIS, NSPINS, NEL
        real(dp) CMAT(NSBASIS, NSBASIS, NSPINS)
        real(dp) DEDCIJ(NSBASIS, NSBASIS, NSPINS)
        INTEGER NDET1(0:NEL + 1), NELS(NSPINS)
        real(dp) UMAT(*)
        INTEGER I, J, K, A, B, C, F, JSPN, KSPN, JJ
        INTEGER IDA, IDB, IDC, IDF
        real(dp) TOT, TOT1, TOT2, TOT1B
!.. We calculate dE/dcij as
!.. dE/dc_kf = 2(Sum_a c_ka <a|h|f>
!.. +Sum_j Sum_abc c_ja c_kb c_jc (<af|U|cb>-<af|U|bc>+<ab|U|cf>-<ab|U|fc>))
        DEDCIJ = 0.0_dp
        DO KSPN = 1, NSPINS
!.. K is an HF orbital
            DO K = 1, NELS(KSPN)
!.. F is the basis orbital
                DO F = 1, NSBASIS
!.. TMAT ID
                    IDF = (F - 1) * NSPINS + 1 + KSPN - 1
                    TOT = 0.0_dp
!.. deal with the one-electron integrals first
                    DO A = 1, NSBASIS
                        IDA = (A - 1) * NSPINS + 1 + KSPN - 1
                        TOT = TOT + CMAT(K, A, KSPN) * (GetTMATEl(IDA, IDF))
                    END DO
!.. UMAT ID
                    IDF = GTID((F - 1) * NSPINS + 1 + KSPN - 1)
                    DO I = 1, NEL
                        JJ = NDET1(I)
                        JSPN = MOD(JJ - 1, NSPINS) + 1
                        J = (JJ - 1) / NSPINS + 1
!.. Here A and C correspond to J and B corresponds to K
                        DO A = 1, NSBASIS
                            IDA = GTID((A - 1) * NSPINS + 1 + JSPN - 1)
                            TOT1 = 0.0_dp
                            DO B = 1, NSBASIS
                                IDB = GTID((B - 1) * NSPINS + 1 + KSPN - 1)
                                TOT1B = 0.0_dp
                                DO C = 1, NSBASIS
                                    TOT2 = 0.0_dp
                                    IDC = GTID((C - 1) * NSPINS + 1 + JSPN - 1)
                                    TOT2 = TOT2 + UMAT(UMatInd(IDA, IDF, IDC, IDB))
                                    TOT2 = TOT2 + UMAT(UMatInd(IDA, IDB, IDC, IDF))
                                    IF (KSPN == JSPN) THEN
                                        TOT2 = TOT2 - UMAT(UmatInd(IDA, IDF, IDB, IDC))
                                        TOT2 = TOT2 - UMAT(UMatInd(IDA, IDB, IDF, IDC))
                                    END IF
                                    TOT1B = TOT1B + TOT2 * CMAT(J, C, JSPN)
                                END DO
                                TOT1 = TOT1 + TOT1B * CMAT(K, B, KSPN)
                            END DO
                            TOT = TOT + TOT1 * CMAT(J, A, JSPN)
                        END DO
                    END DO
                    DEDCIJ(K, F, KSPN) = 2 * TOT
                END DO
            END DO
        END DO
        RETURN
    END
#endif

#ifndef CMPLX_
    SUBROUTINE READHFFMAT(NBASIS, FMAT, HFES, G1, NSPINS, NSBASIS, TRANSP)
        TYPE(BasisFn) G1(*)
        INTEGER NBASIS, NQNS(5), NN, NSPINS, NSBASIS
        real(dp) FMAT(NSBASIS, NSBASIS, NSPINS), HFES(NSBASIS, NSPINS)
        INTEGER I, L, J, NB, NE, IG
        real(dp) VAL
        LOGICAL TRANSP
        character(*), parameter :: this_routine = 'READHFFMAT'
        WRITE(6, *) "Loading HF BASIS"
        OPEN(10, FILE='HFBASIS', STATUS='OLD')
        READ(10, *)
        READ(10, *) NB, NE
!.. NE is NEVAL, and NB is NBASIS/2
!.. NBASIS is the number of orbitals, so *2 to get # spinorbitals
!         IF(NE.NE.NEL) call stop_all(this_routine, 'NEL in HFBASIS <> NEL')

        IF (NE /= NB) call stop_all(this_routine, 'NEVAL <> NBASIS in HFBASIS not supported')
        IF (NB * 2 /= NBASIS) call stop_all(this_routine, 'NBASIS in HFBASIS <> NHG')
        DO I = 1, NB
        DO L = -1, 1, 2
            READ(10, *)
            READ(10, *) HFES(I,(L + 3) / 2)
            NQNS(4) = L
            DO J = 1, NB
                READ(10, *) NN, NQNS(1), NQNS(2), NQNS(3), VAL
                IG = IFINDBASISFN(get_basisfn(NQNS), G1, NBASIS)
                IF (TRANSP) THEN
                    FMAT(I, J,(L + 3) / 2) = VAL
                ELSE
                    FMAT(J, I,(L + 3) / 2) = VAL
                END IF
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
            END DO
        END DO
        END DO
        CLOSE(10)
    END subroutine
#endif


    SUBROUTINE ORDERBASISHF(ARR, BRR, HFE, HFBASIS, NBASIS, FDET, NEL)
        INTEGER NBASIS, BRR(NBASIS)
        real(dp) ARR(NBASIS, 2), HFE(NBASIS), HFBASIS(NBASIS, NBASIS)
!.. HFBASIS(HFBASISFN,PRIMBASISFN) has PRIMBASISFN varying slowest
        INTEGER I, J, IBF, ITOT, NEL, FDET(NEL), ICUR
        real(dp) MX, OEN
        character(*), parameter :: this_routine = 'ORDERBASISHF'
        ICUR = 1
        DO I = 1, NBASIS
        if (ICUR < NEL) then
            if (FDET(ICUR + 1) == I) ICUR = ICUR + 1
        end if
        MX = 0.0_dp
        IBF = 0
        DO J = 1, NBASIS
        IF (ABS(HFBASIS(I, J)) > MX) THEN
            MX = ABS(HFBASIS(I, J))
            IBF = J
        END IF
        END DO
        IF (I == FDET(ICUR) .AND. MX < 0.95_dp) THEN
!.. The largest element isn't big enough, so we abort
            WRITE(6, *) "Largest coeff of HF basis fn ", I, " is ", MX
            WRITE(6, *) "Aborting ORDERBASISHF"
            CALL neci_flush(stdout)
            call stop_all(this_routine, "ORDERBASISHF failed - HF Basis not converged")
        END IF
        ARR(I, 1) = HFE(I)
        ARR(IBF, 2) = HFE(I)
        BRR(I) = IBF
        END DO
!.. We need to now go through each set of degenerate orbitals, and make
!.. the correct ones are paired together in BRR otherwise bad things
!.. happen in FREEZEBASIS
!.. We do this by ensuring that within a degenerate set, the BRR are in
!.. ascending order
        OEN = ARR(1, 1)
        J = 1
!         G1(3,BRR(1))=J
        ITOT = 1
        DO I = 2, NBASIS
        IF (ABS(ARR(I, 1) - OEN) > 1.0e-4_dp) THEN
!.. We don't have degenerate orbitals
!.. First deal with the last set of degenerate orbitals
!.. We sort them into order of BRR
            call sort(brr(i - itot:i - 1), arr(i - itot:i - 1, 1))
!               CALL SORT2_(ITOT,BRR(I-ITOT),ARR(I-ITOT,1))
!.. now setup the new degenerate set.
            J = J + 1
            ITOT = 1
        ELSE
            ITOT = ITOT + 1
        END IF
        OEN = ARR(I, 1)
!.. If we've got a generic spatial sym or hf we mark degeneracies
!               G1(3,BRR(I))=J
        END DO
        RETURN
    END

#ifndef CMPLX_
    SUBROUTINE Write_HEMatrix(CHAR, M, N, A)
        CHARACTER(*) CHAR
        Integer I, M, N, J
        HElement_t(dp) A(M, N)
        WRITE(6, *) CHAR
        DO I = 1, M
            IF (HElement_t_size == 1) THEN
                WRITE(6, "(12E15.6)")(A(I, J), J=1, N)
            ELSE
                WRITE(6, "(6('(',E15.6,',',E15.6,')'))")(A(I, J), J=1, N)
            END IF
        END DO
! 1000 FORMAT(12E15.6)
    END subroutine
#endif


      INTEGER FUNCTION IFINDBASISFN(NQNS,G1,NBASIS)
         INTEGER NBASIS,I,J
         TYPE(BasisFN) G1(NBASIS),NQNS
         LOGICAL L
         DO I=1,NBASIS
            L=.TRUE.
            DO J=1,3
               IF(G1(I)%k(J).NE.NQNS%k(J)) L=.FALSE.
            ENDDO
            IF(G1(I)%Sym%s.NE.NQNS%Sym%s) L=.FALSE.
            IF(G1(I)%Ms.NE.NQNS%Ms) L=.FALSE.
            IF(L) THEN
               IFINDBASISFN=I
               RETURN
            ENDIF
         ENDDO
         IFINDBASISFN=0
         RETURN
      END FUNCTION IFINDBASISFN


end module
