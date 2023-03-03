#include "macros.h"
MODULE HFCalc
    use constants, only: dp, int64, MPIArg
    use util_mod, only: stop_all, neci_flush

    use hfbasis_mod, only: setuphfbasis, calchfbasis, readhftmat, readhfumat, calchfbasis, orderbasishf, calchfumat, calchftmat, &
        readhfbasis


    implicit none
    save
contains
    subroutine HFDoCalc()
        Use global_utilities
        use SystemData, only: BasisFN
        use IntegralsData, only: tHFBasis, tHFCalc, iHFMethod, tReadHF, nHFIt, HFMix, HFCDelta, HFEDelta
        use IntegralsData, only: HFRand, tRHF, ntFrozen, tReadTUMat
        use SystemData, only: tCPMD, tHFOrder, nBasisMax, G1, Arr, Brr, ECore, nEl, nBasis, iSpinSkip, LMS
        use SystemData, only: tHub, lmsbasis
        Use LoggingData, only: iLogging
        Use Determinants, only: FDet, nUHFDet, write_det
        use IntegralsData, only: UMat, tagUMat, umat_win
        Use UMatCache, only: GetUMatSize
        Use OneEInts, only: TMat2D, TMat2D2, SetupTMat2, DestroyTMat
        use sort_mod
        use shared_memory_mpi
        use HElem, only: helement_t_size, helement_t_sizeb
        use MemoryManager, only: TagIntType
        character(25), parameter :: this_routine = 'HFDoCalc'
        HElement_t(dp), ALLOCATABLE :: HFBASIS(:), HFE(:)
        HElement_t(dp), pointer :: UMat2(:)
        INTEGER(MPIArg):: umat2_win
        integer i
        integer nOrbUsed
        integer TMatInt
        integer(int64) :: UMatInt
        integer(TagIntType), save :: tagUMat2 = 0, tagHFE = 0, tagHFBasis = 0

!C.. If we are using an HF basis instead of our primitive basis, we need
!C.. to load in the coeffs of the HF eigenfunctions in terms of the
!C.. primitive basis.
!C.. We load the coeffs from a file HFBASIS
        IF (THFBASIS .OR. THFCALC .OR. (THFORDER .AND. .NOT. TCPMD)) THEN
            ! NOTE: while HFBasis and HFE are declared to be  arrays,
            ! many of the following routines assume them to be real arrays.
            ! These need to be changed for use with complex code.
            allocate(HFBasis(nBasis * nBasis))
            call LogMemAlloc('HFBASIS', nBasis * nBasis, HElement_t_size * 8, this_routine, tagHFBasis)
            HFBASIS = (0.0_dp)
!C.. Allocate an array to store the HF Energies
            allocate(HFE(nBasis))
            call LogMemAlloc('HFE', nBasis, HElement_t_size * 8, this_routine, tagHFE)
            HFE = (0.0_dp)
            IF (THFORDER .AND. .NOT. THFBASIS) THEN
!C.. If we're not using HF, but just calculating the HF order
!C.. We generate the HF energies (this has no mixing or randomisation, so should jsut
!C.. re-order the orbitals and give us some energy)
!C.. HF basis is NOT using the LMS value set in the input
                CALL CALCHFBASIS(NBASIS, NBASISMAX, G1, BRR, ECORE, UMAT, HFE, HFBASIS, 1, NEL, LMSBASIS, 1.0_dp, HFEDELTA, &
                                 HFCDELTA, .TRUE., 0, TREADHF, 0.0_dp, FDET, ILOGGING)
                CALL ORDERBASISHF(ARR, BRR, HFE, HFBASIS, NBASIS, FDET, NEL)
            else if (THFCALC) THEN
                CALL CALCHFBASIS(NBASIS, NBASISMAX, G1, BRR, ECORE, UMAT, HFE, HFBASIS, NHFIT, NEL, LMS, HFMIX, HFEDELTA, &
                                 HFCDELTA, TRHF, IHFMETHOD, TREADHF, HFRAND, FDET, ILOGGING)
                CALL SETUPHFBASIS(NBASISMAX, G1, NBASIS, HFE, ARR, BRR)
            else if (THFBASIS) THEN
                CALL READHFBASIS(HFBASIS, HFE, G1, NBASIS)
                CALL SETUPHFBASIS(NBASISMAX, G1, NBASIS, HFE, ARR, BRR)
            end if
            write(stdout, *) "FINAL HF BASIS"
            CALL WRITEBASIS(6, G1, nBasis, ARR, BRR)

            write(stdout, "(A)", advance='no') " Fermi det (D0):"
            call write_det(6, FDET, .true.)
            CALL neci_flush(6)
!C.. If in Hubbard, we generate site-spin occupations
            IF (THUB) THEN
!  Don't think this works
!               CALL GENSITESPINOCC(NBASIS,NBASIS/ISPINSKIP,ISPINSKIP,NBASISMAX,G1,NEL,LMS,BRR,HFBASIS)
            end if
!C.. We now generate a new U matrix corresponding to the HF basis fns
!C.. This requires a new matrix UMAT2 to be set up, as our HF basis is
!C.. unrestricted
            !THIS ROUTINE NO LONGER WORKS WITH NEW TMAT/UMAT MODULARISATION
            IF (THFBASIS) THEN
                write(stdout, *) "Allocating TMAT2"
                CALL SetupTMAT2(nBasis, 2, TMATINT)
                NORBUSED = NBASIS - NTFROZEN
                IF (TREADTUMAT) THEN
                    CALL READHFTMAT(NBASIS)
                ELSE
                    CALL CALCHFTMAT(NBASIS, HFBASIS, NORBUSED)
                end if
                CALL DestroyTMAT(.false.)
                TMAT2D => TMAT2D2
                NULLIFY (TMAT2D2)
!C.. Allocate the new matrix
                CALL GetUMatSize(nBasis, UMATINT)
                call shared_allocate_mpi(umat2_win, umat2, (/UMatInt/))
                !allocate(UMat2(UMatInt), stat=ierr)
                LogAlloc(ierr, 'UMAT2', int(UMatInt), HElement_t_SizeB, tagUMat2)
                UMAT2 = 0.0_dp
!C.. We need to pass the TMAT to CALCHFUMAT as TMAT is no longer diagona
!C.. This also modified G1, ARR, BRR
                IF (TREADTUMAT) THEN
                    CALL READHFUMAT(UMAT2, NBASIS)
                ELSE
                    CALL CALCHFUMAT(UMAT, UMAT2, NBASIS, HFBASIS, ISPINSKIP, NORBUSED)
                end if
!C.. Now we can remove the old UMATRIX, and set the pointer UMAT to point
!C.. to UMAT2
                call shared_sync_mpi(umat2_win)
                LogDealloc(tagUMat)
                call shared_deallocate_mpi(umat_win, umat)

                !Deallocate(UMat)
                umat_win = umat2_win
                UMat => UMat2
                nullify (UMat2)
                tagUMat = tagUMat2
                tagUMat2 = 0
!C.. spinskip values
                ISPINSKIP = 1
                NBASISMAX(2, 3) = 1
!C.. Indicate that we're in UHF and that <D0|H|D1>=0 for D1 being a
!C.. single excitation
                IF (LMS == 0) THEN
                    NBASISMAX(4, 5) = 1
                    DO I = 1, NEL
                        NUHFDET(I) = BRR(I)
                    end do
                    call sort(nUHFDet(1:nel))
                ELSE
                    NBASISMAX(4, 5) = 2
                end if
            end if
!C.. Now deallocate the HF arrays
            deallocate(HFE, HFBasis)
            call LogMemDealloc(this_routine, tagHFE)
            call LogMemDealloc(this_routine, tagHFBasis)
            CALL neci_flush(6)
        end if
    End Subroutine HFDoCalc
End Module HFCalc
