#include "macros.h"

module hubbard_mod
    USE OneEInts, only: TMat2D, TMATSYM, SetupTMAT
    USE Parallel_neci, only: iProcIndex
    use SymData, only: SymReps, tagSymReps, &
        nSym, SymConjTab, SymClasses, SymLabels, nSymLabels, tAbelian, &
        SymTable, tagSymConjTab, tagSymClasses, tagSymLabels, tagSymTable
    use SystemData, only: BasisFN, BasisFNSize, BasisFNSizeB, &
        t_open_bc_x, t_open_bc_y, Symmetry, SymmetrySize, SymmetrySizeB, &
        t_k_space_hubbard, twisted_bc, nn_bhub
    use UMatCache, only: UMatInd
    use constants, only: Pi, Pi2, THIRD, dp, int64, sizeof_int
    use global_utilities, only: timer, set_timer, halt_timer, &
        LogMemAlloc, LogMemDealloc
    use hfbasis_mod, only: iFindBasisFn
    use lattice_mod, only: lat
    use sym_mod, only: RoundSym, AddElecSym, mompbcsym
    use util_mod, only: stop_all
    better_implicit_none
    private
    public :: genhubmomirrepssymtable, genhubsymreps, &
        hubkin, hubkinn, setbasislim_hubtilt, setbasislim_hub, calctmathub, &
        calcumathubreal, write_kspace_umat

contains

    SUBROUTINE CALCUMATHUBREAL(NBASIS, UHUB, UMAT)
        INTEGER NBASIS
        HElement_t(dp) UMAT(*)
        real(dp) UHUB
        INTEGER I, IC
        OPEN(10, FILE='UMAT', STATUS='UNKNOWN')
        IC = 1
!Only do this for each spatial orbital.
        DO I = 1, NBASIS / 2
            UMAT(UMatInd(I, I, I, I)) = h_cast(UHUB)
            WRITE(10, *) I, I, I, I, UHUB, UMatInd(I, I, I, I)
        END DO
        CLOSE(10)
        RETURN
    END

    SUBROUTINE HUBKIN(I, J, K, NBASISMAX, BHUB, TTILT, TOTSUM, TREAL)
! returns the non-interacting energy of state with
!..quatum numbers (i,j,k) for Hubbard model of
!..lengths (LX,LY,LZ) PBC.
!..Returned in sum
        INTEGER nBasisMax(5, *), AX, AY
        INTEGER LX, LY, LZ, K, J, I
        LOGICAL TTILT, TREAL
        real(dp) :: TOTSUM, BHUB
        real(dp) :: tb_x, tb_y, nn
        IF (TREAL) THEN
            TOTSUM = 0.0_dp
        ELSE
            tb_x = twisted_bc(1)
            tb_y = twisted_bc(2)

            LX = NBASISMAX(1, 2) - NBASISMAX(1, 1) + 1
            LY = NBASISMAX(2, 2) - NBASISMAX(2, 1) + 1
            LZ = NBASISMAX(3, 2) - NBASISMAX(3, 1) + 1
            AX = NBASISMAX(1, 4)
            AY = NBASISMAX(2, 4)
            IF (TTILT) THEN
!CCC.. NBASISMAX goes from -NMAXX+1 to MAXX so LX=2MAXX
                LX = NBASISMAX(1, 5)
                LX = (LX * (AX * AX + AY * AY))
                totsum = cos(2 * pi * ((I + tb_x) * ax + (j + tb_y) * ay) / lx)
                totsum = totsum + cos(2 * pi * ((i + tb_x) * ay - (j + tb_y) * ax) / lx)

                nn = cos(2 * pi * (i + tb_x) * AX / lx) + cos((J + tb_y) * ay / lx)
            ELSE
                IF (LX > 1) then
                    TOTSUM = COS(2 * PI * (I + tb_x) / LX)
                    nn = cos(2 * pi * (I + tb_x + J * tb_y) / LX)
                end if
                IF (LY > 1) then
                    TOTSUM = TOTSUM + COS(2 * PI * (J + tb_y) / LY)
                    nn = cos(2 * pi * (I + tb_x - (J + tb_y)) / LY)
                end if
                IF (LZ > 1) TOTSUM = TOTSUM + COS(2 * PI * K / (LZ))
            END IF
            TOTSUM = 2.0_dp * (TOTSUM * BHUB + nn * nn_bhub)
        END IF
        RETURN
    END

    SUBROUTINE HUBKINN(I, J, K, NBASISMAX, BHUB, TTILT, TOTSUM, TREAL)
! returns the non-interacting energy of state with
!..quatum numbers (i,j,k) for Hubbard model of
!..lengths (LX,LY,LZ).  NON-PBC
!..Returned in sum
        INTEGER nBasisMax(5, *), AX, AY
        LOGICAL TTILT, TREAL
        INTEGER II, JJ, KK, K, J, I, LZ, LX, LY
        real(dp) :: BHUB, TOTSUM
        if (tReal) then
            TOTSUM = 0.0_dp
            return
        else

            LX = NBASISMAX(1, 2) - NBASISMAX(1, 1) + 1
            LY = NBASISMAX(2, 2) - NBASISMAX(2, 1) + 1
            LZ = NBASISMAX(3, 2) - NBASISMAX(3, 1) + 1
            AX = NBASISMAX(1, 4)
            AY = NBASISMAX(2, 4)

            IF (TTILT) THEN
!CCC.. NBASISMAX goes from -NMAXX+1 to MAXX so LX=2MAXX
                LX = NBASISMAX(1, 5)
                LX = (LX * (AX * AX + AY * AY))
                TOTSUM = COS(PI * (I * AX - J * AY) / (LX + 1)) + COS(PI * (I * AY + J * AX) / (LX + 1))
            ELSE
                II = I - NBASISMAX(1, 1) + 1
                JJ = J - NBASISMAX(2, 1) + 1
                KK = K - NBASISMAX(3, 1) + 1
                IF (LX > 1) TOTSUM = COS(PI * II / (LX + 1))
                IF (LY > 1) TOTSUM = TOTSUM + COS(PI * JJ / (LY + 1))
                IF (LZ > 1) TOTSUM = TOTSUM + COS(PI * KK / (LZ + 1))
!         IF(LY.GT.1) TOTSUM=TOTSUM+COS(PI*J/(LY+1))
!         IF(LZ.GT.1) TOTSUM=TOTSUM+COS(PI*K/(LZ+1))
            END IF
            TOTSUM = TOTSUM * 2.0_dp * BHUB
            RETURN
        end if
    END

    SUBROUTINE CALCTMATHUB(NBASIS, NBASISMAX, BHUB, TTILT, G1, TREAL, TPBC)
        INTEGER NBASIS, nBasisMax(5, *)
        TYPE(BasisFN) G1(nBasis)
        real(dp) BHUB
        INTEGER iSize
        INTEGER I, J
        INTEGER DX, DY, DZ, LX, LY, LZ
        LOGICAL TTILT, TREAL, TPBC
        real(dp) TOTSUM
        integer temp_k(3)
        IF (iProcIndex == 0) OPEN(10, FILE='TMAT', STATUS='UNKNOWN')
        CALL SetupTMAT(NBASIS, 2, iSize)
        TOTSUM = 0.0_dp
        LX = NBASISMAX(1, 2) - NBASISMAX(1, 1)
        LY = NBASISMAX(2, 2) - NBASISMAX(2, 1)
        LZ = NBASISMAX(3, 2) - NBASISMAX(3, 1)
! here i have to make some changes in the real-space tilted
! case to construct the correct TMAT
        IF (LY == 0) LY = -1
        IF (LZ == 0) LZ = -1
        IF (TREAL) THEN
        DO I = 1, NBASIS
        DO J = 1, NBASIS
            DX = ABS(G1(I)%k(1) - G1(J)%k(1))
            DY = ABS(G1(I)%k(2) - G1(J)%k(2))
            DZ = ABS(G1(I)%k(3) - G1(J)%k(3))
            temp_k = [DX, DY, DZ]
            call mompbcsym(temp_k, NBASISMAX)
            temp_k = abs(temp_k)
            TOTSUM = BHUB
            IF (TPBC) THEN
!This bit is only for if the hubbard lattice only has one site in a certain dimension
                IF (DX == LX .and. .not. t_open_bc_x) THEN
                    DX = 1
                    IF (LX == 1) TOTSUM = TOTSUM + BHUB
                END IF
                IF (DY == LY .and. .not. t_open_bc_y) THEN
                    DY = 1
                    IF (LY == 1) TOTSUM = TOTSUM + BHUB
                END IF
                IF (DZ == LZ) THEN
                    DZ = 1
                    IF (LZ == 1) TOTSUM = TOTSUM + BHUB
                END IF
            END IF
! make this more general to also use open BC and a
! tilted real-space lattice
            if (.not. t_open_bc_x .and. .not. t_open_bc_y) then
            IF (DX + DY + DZ == 1 .or. sum(temp_k) == 1) then
            if (G1(I)%Ms == G1(J)%MS) THEN
!This is for if the site can interact with a periodic image.
                TMAT2D(I, J) = TOTSUM
                IF (iProcIndex == 0) WRITE(10, *) I, J, TOTSUM
            end if
            END IF
            else
            if (DX + DY + DZ == 1 .and. G1(I)%Ms == G1(J)%MS) THEN
                TMAT2D(I, J) = TOTSUM
                IF (iProcIndex == 0) WRITE(10, *) I, J, TOTSUM
            end if
            end if
        END DO
        END DO
        ELSE
        IF (TPBC) THEN
        DO I = 1, NBASIS
            CALL HUBKIN(G1(I)%k(1), G1(I)%k(2), G1(I)%k(3), NBASISMAX, BHUB, TTILT, TOTSUM, TREAL)
            TMAT2D(I, I) = TOTSUM
            IF (iProcIndex == 0) WRITE(10, *) I, I, TMAT2D(I, I)
        END DO
        ELSE
        DO I = 1, NBASIS
            CALL HUBKINN(G1(I)%k(1), G1(I)%k(2), G1(I)%k(3), NBASISMAX, BHUB, TTILT, TOTSUM, TREAL)
            TMAT2D(I, I) = TOTSUM
            IF (iProcIndex == 0) WRITE(10, *) I, I, TMAT2D(I, I)
        END DO
        END IF
        END IF
        IF (iProcIndex == 0) CLOSE(10)
        RETURN
    END
!.. NBASISMAX descriptor (1,3)
!
! HUBBARD:
! 0 Non-Tilted Lattice - pbc
! 1 Tilted Lattice - pbc
! 2 Non-Tilted lattice - no pbc
! 3 Tilted Lattice - no pbc
!.. four following are REAL
! 4 Non-Tilted Lattice - pbc
! 5 Tilted Lattice - pbc
! 6 Non-Tilted lattice - no pbc
! 7 Tilted Lattice - no pbc
!
    SUBROUTINE SETBASISLIM_HUB(NBASISMAX, NMAXX, NMAXY, NMAXZ, LEN, TPBC, TREAL)
        INTEGER nBasisMax(5, *), NMAXX, NMAXY, NMAXZ, LEN
        LOGICAL TPBC, TREAL
        IF (TPBC) THEN
            NBASISMAX(1, 3) = 0
        ELSE
!.. Non-tilted, not pbc
            NBASISMAX(1, 3) = 2
        END IF
        IF (TREAL) NBASISMAX(1, 3) = NBASISMAX(1, 3) + 4
!         IF(.NOT.TPBC.AND..NOT.TREAL) THEN
!.. non-pbc has Huckel MOs starting k from 1
!            NBASISMAX(1,1)=1
!            NBASISMAX(2,1)=1
!            NBASISMAX(3,1)=1
!            NBASISMAX(1,2)=NMAXX
!            NBASISMAX(2,2)=NMAXY
!            NBASISMAX(3,2)=NMAXZ
!         ELSE
        IF (MOD(NMAXX, 2) == 0) THEN
            NBASISMAX(1, 2) = NMAXX / 2
            NBASISMAX(1, 1) = -NMAXX / 2 + 1
        ELSE
            NBASISMAX(1, 2) = NMAXX / 2
            NBASISMAX(1, 1) = -NMAXX / 2
        END IF
        IF (MOD(NMAXY, 2) == 0) THEN
            NBASISMAX(2, 2) = NMAXY / 2
            NBASISMAX(2, 1) = -NMAXY / 2 + 1
        ELSE
            NBASISMAX(2, 2) = NMAXY / 2
            NBASISMAX(2, 1) = -NMAXY / 2
        END IF
        IF (MOD(NMAXZ, 2) == 0) THEN
            NBASISMAX(3, 2) = NMAXZ / 2
            NBASISMAX(3, 1) = -NMAXZ / 2 + 1
        ELSE
            NBASISMAX(3, 2) = NMAXZ / 2
            NBASISMAX(3, 1) = -NMAXZ / 2
        END IF
!         ENDIF
        NBASISMAX(1, 4) = 0
        NBASISMAX(2, 4) = 1
        NBASISMAX(1, 5) = NMAXX
        NBASISMAX(2, 5) = NMAXY
        NBASISMAX(3, 5) = NMAXZ
        LEN = NMAXX * NMAXY * NMAXZ * ((NBASISMAX(4, 2) - NBASISMAX(4, 1)) / 2 + 1)
    END

    SUBROUTINE SETBASISLIM_HUBTILT(NBASISMAX, NMAXX, NMAXY, NMAXZ, LEN, TPBC, ITILTX, ITILTY)
        INTEGER nBasisMax(5, *), NMAXX, NMAXY, NMAXZ, LEN
        LOGICAL TPBC
        INTEGER ITILTX, ITILTY
        character(*), parameter :: this_routine = 'SETBASISLIM_HUBTILT'
        IF (TPBC) THEN
!.. Indicate tilted
            NBASISMAX(1, 3) = 1
        ELSE
!.. Not periodic boundaries
            NBASISMAX(1, 3) = 3
        END IF
        NBASISMAX(1, 2) = int(real(NMAXX, dp) * (real(ITILTX + ITILTY, dp) / 2.0_dp), sizeof_int)
        NBASISMAX(1, 1) = -NBASISMAX(1, 2)
!+1-MOD(NMAXX,2)
        IF (NMAXY /= NMAXX) call stop_all(this_routine, 'CANNOT HANDLE NON-SQUARE TILTED HUBBARD')
        NBASISMAX(2, 2) = int(real(NMAXX, dp) * (real(ITILTY + ITILTX, dp) / 2.0_dp), sizeof_int)
        NBASISMAX(2, 1) = -NBASISMAX(2, 2)
!+1-MOD(NMAXX,2)
        IF (NMAXZ > 1) call stop_all(this_routine, 'CANNOT HANDLE TILTED 3D HUBBARD')
        NBASISMAX(3, 2) = 0
        NBASISMAX(3, 1) = 0
        NBASISMAX(1, 4) = ITILTX
        NBASISMAX(2, 4) = ITILTY
        NBASISMAX(1, 5) = NMAXX
        NBASISMAX(2, 5) = NMAXY
        NBASISMAX(3, 5) = NMAXZ
!Len is number of basis functions
        LEN = NMAXX * NMAXY * (ITILTX * ITILTX + ITILTY * ITILTY) * ((NBASISMAX(4, 2) - NBASISMAX(4, 1)) / 2 + 1)
!         LEN=NMAXX*NMAXY*2
!     &         *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
    END

!Generate the Sym table for a Hubbard Lattice
    Subroutine GenHubMomIrrepsSymTable(G1, nBasis, nBasisMax)
        INTEGER nBasis, nBasisMax(*)
        TYPE(BasisFN) G1(nBasis)
        INTEGER J, nSyms, K
        integer(int64) :: I
        TYPE(BasisFN) NQNS, S
        character(*), parameter:: this_routine = 'GenHubMomIrrepsSymTable'
        nSym = nBasis / 2
        WRITE(6, "(A,I3,A)") "Generating abelian symmetry table with", nBasis / 2, " generators for Hubbard momentum"

!.. Now generate a list of sym labels.
        write(6, *) 'SIZES', nSymLabels, nBasis, allocated(symlabels), associated(symclasses), allocated(symconjtab), allocated(symtable)
        if (allocated(SymLabels)) then
            write(6, '(a/a)') 'Warning: symmetry info already allocated.', 'Deallocating and reallocating.'
            deallocate(SymLabels)
            call LogMemDealloc(this_routine, tagSymLabels)
        end if
        allocate(SymLabels(nSym))
        call LogMemAlloc('SymLabels', nSym, SymmetrySize, this_routine, tagSymLabels)
        if (associated(SymClasses)) then
            deallocate(SymClasses)
            call LogMemDealloc(this_routine, tagSymClasses)
        end if
        allocate(SymClasses(nBasis))
        call LogMemAlloc('SymClasses', nBasis, 4, this_routine, tagSymClasses)
        NSYMLABELS = NSYM
        nSyms = 1
        tAbelian = .false.
        DO I = 1, NBASIS / 2
!.. place the sym label of each state in SymClasses(ISTATE).
            IF (ALL(G1(I * 2)%k == 0)) THEN
                SymClasses(I) = 1
                NQNS = G1(I * 2)
            ELSE
                nSyms = nSyms + 1
                SymClasses(I) = nSyms
            END IF
!.. list the symmetry string of each sym label
!.. try the new symmetry encoding here
            if (t_k_space_hubbard) then
                SymLabels(I)%s = i
            else
                SymLabels(I)%s = 2_int64**(I - 1_int64)
            end if
        END DO
!.. Setup the symmetry product table
        if (allocated(SymTable)) then
            deallocate(SymTable)
            call LogMemDealloc(this_routine, tagSymTable)
        end if
        allocate(SymTable(nSym, nSym))
        call LogMemAlloc('SymTable', nSym**2, SymmetrySize, this_routine, tagSymTable)
        if (allocated(SymConjTab)) then
            deallocate(SymConjTab)
            call LogMemDealloc(this_routine, tagSymConjTab)
        end if
        allocate(SymConjTab(nSym))
        call LogMemAlloc('SymConjTable', nSym, 4, this_routine, tagSymConjTab)
        SYMTABLE(1:NSYM, 1:NSYM) = Symmetry(0)

        DO I = 1, nBasis / 2
            NQNS%k = -G1(I * 2)%k
            CALL RoundSym(NQNS, nBasisMax)
            if (t_k_space_hubbard) then
! i could actually include the symmetry rounding in the
! get_orb_from_k_vec functionality!
! i was doing this incorrectly.. i want the spin orbital
! j
                J = lat%get_orb_from_k_vec(NQNS%k, 2)
            else
                J = iFindBasisFn(NQNS, G1, nBasis)
            end if
            IF (J == 0) THEN
                WRITE(6, *) "Cannot Find symmetry conjugate to basis fn ", I * 2
                call stop_all(this_routine, "Cannot find symmetry conjugate.")
            END IF
            SymConjTab(SymClasses(I)) = SymClasses(J / 2)
            DO J = 1, nBasis / 2
                S = G1(I * 2)
                CALL AddElecSym(J * 2, G1, nBasisMax, S)
                CALL RoundSym(S, nBasisMax)
                NQNs%k = S%k
                if (t_k_space_hubbard) then
                    k = lat%get_orb_from_k_vec(S%k, 2)
                else
                    K = iFindBasisFn(NQNS, G1, nBasis)
                end if

                IF (K == 0) THEN
                    WRITE(6, *) "Cannot find symmetry product of basis fns ", I * 2, J * 2
                    call stop_all(this_routine, "Cannot find symmetry product.")
                END IF
                SymTable(SymClasses(I), SymClasses(J)) = SymLabels(SymClasses(K / 2))
            END DO
        END DO
        DO I = 1, nBasis / 2
            G1(I * 2 - 1)%Sym = SymLabels(SymClasses(I))
            G1(I * 2)%Sym = SymLabels(SymClasses(I))
        END DO
        WRITE(6, *) "Symmetry, Symmetry Conjugate"
        DO I = 1, NSYM
            WRITE(6, *) I, SymConjTab(I)
        END DO
    End

!Hubbard Sym Reps are different from normal ones.  In hubbard, we have split degenerate sets into 1D subcomponents,
!  but a complete degenerate set will count as a sym rep (i.e. if it is filled, its sym can be discounted),
!  not just one of the subcomponents
    SUBROUTINE GENHUBSYMREPS(NBASIS, ARR, BRR)
        INTEGER I, J
        INTEGER NBASIS, BRR(NBASIS)
        real(dp) ARR(NBASIS)
        character(*), parameter :: this_routine = 'GenHubSymReps'

!.. now work out which reps are degenerate and label them
        if (allocated(SymReps)) deallocate(symreps)
        allocate(SymReps(2, nBasis))
        call LogMemAlloc('SymReps', 2 * nBasis, 4, this_routine, tagSymReps)
! this does not make sense..
        SYMREPS = 0
!.. we have a new rep
        J = 1
        SYMREPS(1, BRR(1)) = 1
        SYMREPS(2, 1) = 1
        DO I = 2, NBASIS
        IF (ABS(ARR(I) - ARR(I - 1)) < 1.0e-5_dp) THEN
!.. we have the same degenerate rep as the previous entry
            SYMREPS(2, J) = SYMREPS(2, J) + 1
        ELSE
!.. we have a new rep
            J = J + 1
            SYMREPS(2, J) = 1
        END IF
        SYMREPS(1, BRR(I)) = J
        END DO

    END


    subroutine write_kspace_umat()
! subroutine to output the umat also in the case of the
! momentum space hubbard model, to be easier able to compare it
! with the DMRG calculations for the GUGA matrix elements!
        use SystemData, only: G1, nSpatOrbs, nBasisMax, uHub, omega
        use sym_mod, only: mompbcsym
        use UMatCache, only: gtid

        integer :: i, j, k, l, k_in(3), k_out(3)
        open(10, file="UMAT", status="unknown")
! do it really naively and just loop over all the indices and
! check if the momentum consercation is fullfilled!
! actually it is more efficient to loop over the spatial orbitals
! only! and check the momentum for the spin-orbs!
        do i = 1, nSpatOrbs
        do j = 1, nSpatOrbs
        do k = 1, nSpatOrbs
        do l = 1, nSpatOrbs
! have to figure out which orbitals to compare in
! the physicist notation!
! and convert to fake spin-orbtitals
            k_in = G1(2 * i)%k + G1(2 * k)%k
            k_out = G1(2 * j)%k + G1(2 * l)%k
! apply periodic BC
            call mompbcsym(k_in, nBasisMax)
            call mompbcsym(k_out, nBasisMax)
            if (all(k_in == k_out)) then
                write(10, '(4i7, f19.9)') i, j, k, l, uHub / Omega
            end if
        end do
        end do
        end do
        end do
        close(10)

    end subroutine write_kspace_umat

end module
