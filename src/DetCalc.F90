#include "macros.h"
MODULE DetCalc
    use constants, only: dp,n_int
    use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB, tStoreSpinOrbs, &
         t_non_hermitian_2_body
    use sort_mod

    use bit_reps, only: writebitdet

    use DetCalcData

    use MemoryManager, only: TagIntType

    use gndts_mod, only: gndts

    use UMatCache, only: UMat2D, tUMat2D, tDeferred_UMat2D, SetupUMat2d_dense

    use procedure_pointers, only: get_umat_el

    use excit_mod, only: isvaliddet, genexcit

    use gndts_blk_mod, only: gndts_blk

    IMPLICIT NONE
    save

!From input
    INTEGER DETINV !The index in the list of dets of a det to investigate
    INTEGER IOBS, JOBS, KOBS

    LOGICAL TCALCHMAT, TENERGY, TREAD, TBLOCK
    LOGICAL tFindDets           !Set if we are to enumerate all determinants within given constraints
    LOGICAL tCompressDets       !Set if once we've found the dets we compress to bit format

    TYPE(BasisFN), pointer :: BLOCKSYM(:)  !The Symmetry of each block.  nBlocks elements
    INTEGER(TagIntType) :: tagBlockSym
    INTEGER, ALLOCATABLE :: NBLOCKSTARTS(:) !Index of the first det of different symmetry blocks in the complete list of dets
    INTEGER(TagIntType) :: tagNBLOCKSTARTS = 0
    INTEGER NBLOCKS                        !Number of Symmetry blocks
    INTEGER iFDet                       ! The index of the Fermi det in the list of dets.
    HElement_t(dp), pointer :: CKN(:, :) !  (nDet,nEval)  Temporary storage for the Lanczos routine
    INTEGER(TagIntType) :: tagCKN = 0

    real(dp), ALLOCATABLE :: ExpandedHamil(:, :)    ! (NDet,NDet) This is the hamiltonian in expanded form,
    !so that it can be histogrammed against.

    INTEGER iExcitLevel                 ! The excitation level at which determinants are cut off.

CONTAINS
    Subroutine DetCalcInit

        Use global_utilities
        Use Determinants, only: FDet, specdet, tSpecDet, tDefineDet, &
                                DefDet, write_det
        Use IntegralsData, only: NFROZEN
        use SystemData, only: lms, lms2, nBasis, nBasisMax, nEl, SymRestrict
        use SystemData, only: Alat, arr, brr, boa, box, coa, ecore, g1, Beta
        use SystemData, only: tParity, tSpn, Symmetry, STot, NullBasisFn, tUHF, tMolpro
        use sym_mod
        use LoggingData, only: tLogDets
        use HElem
        use util_mod, only: get_free_unit, NECI_ICOPY
        Type(BasisFn) ISym

        integer i, j, ii, iunit
        integer ierr, norb
        integer nDetTot

        character(25), parameter :: this_routine = 'DetCalcInit'

        IF (.NOT. TCALCHMAT) THEN
            write(stdout, *) "Not storing the H matrix."
            IF (TENERGY .AND. .NOT. TBLOCK) THEN
                write(stdout, *) "Cannot calculate energies without blocking the Hamiltonian."
                TENERGY = .FALSE.
            end if
            IF (TENERGY .AND. NBLK /= 0) THEN
!C.. We're doing a Lanczos without calculating the H mat
                write(stdout, *) "Cannot perform Lanczos without Hamiltonian"
                TENERGY = .FALSE.
            end if
        end if

        ! If we want to have UMat2D, and it isn't yet filled in, generate it
        ! here. All of the integrals setup/freezing etc is done...
        if (tDeferred_Umat2d .and. .not. tUMat2D) then

            ASSERT(.not. btest(nbasis, 0))

            ! And fill in the array
            call SetupUMat2d_dense(nBasis)
        end if

!Copied Specdet information from Calc.F, so if inspect is present, but no determinant/csf specified, it will still run.
        if (TSPECDET) then
            if (.not. associated(specdet)) then
                !specdet not allocated. Allocate it and copy fdet
                allocate(specdet(nel))
                write(stdout, *) "TSPECDET set, but not allocated.  using FDET"
                CALL NECI_ICOPY(NEL, FDET, 1, SPECDET, 1)
            else if (.not. ISVALIDDET(SPECDET, NEL)) then
                write(stdout, *) "TSPECDET set, but invalid.  using FDET"
                CALL NECI_ICOPY(NEL, FDET, 1, SPECDET, 1)
            end if
        end if

!C      IF(TCALCHMAT.OR.NPATHS.NE.0.OR.DETINV.GT.0.OR.TBLOCK) THEN
        iExcitLevel = ICILEVEL
        IF (tFindDets) THEN
!C..Need to determine the determinants
            IF (iExcitLevel /= 0) THEN
                write(stdout, *) "Performing truncated CI at level ", iExcitLevel
                IF (TSPECDET) THEN
                    write(stdout, *) "Using SPECDET:"
                    call write_det(6, SPECDET, .true.)!
                    CALL NECI_ICOPY(NEL, SPECDET, 1, FDET, 1)
                ELSE
                    write(stdout, *) "Using Fermi DET:"
                    call write_det(6, FDET, .true.)
                end if
!C.. if we're doing a truncated CI expansion
                CALL GENEXCIT(FDET, iExcitLevel, NBASIS, NEL, reshape([0], [1, 1]), &
                              reshape([0], [1, 1]), NDET, 1, G1, .TRUE., NBASISMAX, .TRUE.)
                write(stdout, *) "NDET out of GENEXCIT ", NDET
!C.. We need to add in the FDET
                NDET = NDET + 1
                II = NDET
                NBLOCKS = 1
            else if (TBLOCK) THEN
                write(stdout, *) "Determining determinants and blocks."
                IF (TPARITY) THEN
                    write(stdout, *) "Using symmetry restriction:"
                    CALL WRITEALLSYM(6, SymRestrict, .TRUE.)
                end if
                IF (TSPN) THEN
                    write(stdout, *) "Using spin restriction:", LMS
                end if
                if (tUHF .and. tMolpro) then
                    !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                    CALL GNDTS_BLK(NEL, nBasis, BRR, NBASISMAX, NMRKS, .TRUE.,             &
         &                NDET, G1, II, NBLOCKSTARTS, NBLOCKS, TSPN, -LMS2, TPARITY,        &
         &               SymRestrict, IFDET,.NOT. TREAD, NDETTOT, BLOCKSYM)
                else
                    CALL GNDTS_BLK(NEL, nBasis, BRR, NBASISMAX, NMRKS, .TRUE.,             &
         &                NDET, G1, II, NBLOCKSTARTS, NBLOCKS, TSPN, LMS2, TPARITY,        &
         &               SymRestrict, IFDET,.NOT. TREAD, NDETTOT, BLOCKSYM)
                end if
                write(stdout, *) "NBLOCKS:", NBLOCKS
                write(stdout, *) "Determining determinants."
                IF (TPARITY) THEN
                    write(stdout, *) "Using symmetry restriction:"
                    CALL WRITEALLSYM(6, SymRestrict, .TRUE.)
                end if
                IF (TSPN) THEN
                    write(stdout, *) "Using spin restriction:", LMS
                end if
                if (tUHF .and. tMolpro) then
                    !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                    CALL GNDTS(NEL, nBasis, BRR, NBASISMAX, NMRKS, .TRUE., G1, TSPN, -LMS, TPARITY, SymRestrict, II, IFDET)
                else
                    CALL GNDTS(NEL, nBasis, BRR, NBASISMAX, NMRKS, .TRUE., G1, TSPN, LMS, TPARITY, SymRestrict, II, IFDET)
                end if
                NBLOCKS = 1
                NDET = II
            end if
!C..
            IF (II == 0) THEN
                write(stdout, *) "No determinants found.  Cannot continue"
                call stop_all(this_routine, "No determinants found.  Cannot continue")
            end if
!C.. NEL now only includes active electrons
            write(stdout, *) "Number of determinants found to be: ", II
            write(stdout, *) "Allocating initial memory for calculation of energy..."
            CALL neci_flush(6)
            allocate(NMrks(nEl, II), stat=ierr)
            LogAlloc(ierr, 'NMRKS', NEL * II, 4, tagNMRKS)
            NMRKS(1:NEL, 1:II) = 0
            allocate(NBLOCKSTARTS(NBLOCKS + 1), stat=ierr)
            call LogMemAlloc('NBLOCKSTARTS', NBLOCKS + 1, 4, this_routine, tagNBLOCKSTARTS, ierr)
            NBLOCKSTARTS(1:NBLOCKS + 1) = 0
            allocate(BlockSym(NBLOCKS + 1), stat=ierr)
            LogAlloc(ierr, 'BLOCKSYM', NBLOCKS + 1, BasisFNSizeB, tagBlockSym)

            BLOCKSYM(1:NBLOCKS) = NullBasisFn
!C..

            NDET = II
            IF (iExcitLevel /= 0) THEN
!C.. Use HAMIL to temporarily hold a list of excitation levels
                CALL NECI_ICOPY(NEL, FDET, 1, NMRKS, 1)
                allocate(Hamil(II), stat=ierr)
                LogAlloc(ierr, 'HAMIL', II, HElement_t_sizeB, tagHamil)
                NDET = 0
                CALL GENEXCIT(&
                        FDET, iExcitLevel, NBASIS, NEL, &
                        NMRKS(:1, :2), int(HAMIL), NDET, 1, G1, .TRUE., NBASISMAX, .FALSE.)
                Deallocate(Hamil)
                LogDealloc(tagHamil)
                NDET = NDET + 1
                NBLOCKSTARTS(1) = 1
                NBLOCKSTARTS(2) = II + 1
                IFDET = 1
            else if (TBLOCK) THEN
                if (tUHF .and. tMolpro) then
                    !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                    CALL GNDTS_BLK(NEL, nBasis, BRR, NBASISMAX, NMRKS, .FALSE., NDET, G1, II, NBLOCKSTARTS, NBLOCKS, TSPN, -LMS2, TPARITY, &
         &               SymRestrict, IFDET,.NOT. TREAD, NDETTOT, BLOCKSYM)
                else
                    CALL GNDTS_BLK(NEL, nBasis, BRR, NBASISMAX, NMRKS, .FALSE., NDET, G1, II, NBLOCKSTARTS, NBLOCKS, TSPN, LMS2, TPARITY, &
         &               SymRestrict, IFDET,.NOT. TREAD, NDETTOT, BLOCKSYM)
                end if
            ELSE
                if (tUHF .and. tMolpro) then
                    !When breaking spin symmetry in molpro, it is important to occupy alpha orbs preferentially
                    CALL GNDTS(NEL, nBasis, BRR, NBASISMAX, NMRKS, .FALSE., G1, TSPN, -LMS, TPARITY, SymRestrict, II, IFDET)
                else
                    CALL GNDTS(NEL, nBasis, BRR, NBASISMAX, NMRKS, .FALSE., G1, TSPN, LMS, TPARITY, SymRestrict, II, IFDET)
                end if
                NBLOCKSTARTS(1) = 1
                NBLOCKSTARTS(2) = II + 1
            end if
            if (tLogDets) THEN
                iunit = get_free_unit()
                open(iunit, FILE='DETS', STATUS='UNKNOWN')
                DO I = 1, NDET
                    call write_det(iunit, NMRKS(:, I), .false.)
                    CALL GETSYM(NMRKS(:, I), NEL, G1, NBASISMAX, ISYM)
                    CALL WRITESYM(iunit, ISym%Sym, .TRUE.)
                end do
                close(iunit)
            end if

!Update: 14.03.2018, K.Ghanem
!FDET has already been assigned in DetPreFreezeInit.
!No idea why it is overwirtten here.
!Therefore, I comment out the following code and hope for the best.
!Instead, I look for FDET in the list of determinants NMRKS and assign the index to IFDET.

!C.. Now generate the fermi determiant
!C.. Work out the fermi det
            IFDET = 0
            DO I = 1, NDET
                IF (ALL(NMRKS(:, I) == FDET)) THEN
                    IFDET = I
                    Exit
                END IF
            END DO
            IF (IFDET == 0) call stop_all("DetCalcInit", "Fermi determinant is not found in NMRKS!")

            write(stdout, *) ' NUMBER OF SYMMETRY UNIQUE DETS ', NDET

!C         write(stdout,*) ' TOTAL NUMBER OF DETS.' , NDETTOT
            IF (NEVAL == 0) THEN
                write(stdout, *) 'NEVAL=0.  Setting NEVAL=NDET'
                NEVAL = NDET
            end if
            IF (NEVAL > NDET) THEN
                write(stdout, *) 'NEVAL>NDET.  Setting NEVAL=NDET'
                NEVAL = NDET
            end if

            IF (ABS(DETINV) > NDET) THEN
                write(stdout, *) 'DETINV=', DETINV, '>NDET=', NDET
                write(stdout, *) 'Setting DETINV to 0'
                DETINV = 0
            end if
            CALL neci_flush(6)

!C ==----------------------------------------------------------------==
!C..Set up memory for c's, nrow and the label
            IF (TCALCHMAT) THEN
                write(stdout, *) "CK Size", NDET * NEVAL * HElement_t_size
                allocate(CkN(nDet, nEval), stat=ierr)
                LogAlloc(ierr, 'CKN', nDet * nEval, HElement_t_sizeB, tagCKN)
                CKN = (0.0_dp)
!C..
                allocate(Ck(nDet, nEval), stat=ierr)
                LogAlloc(ierr, 'CK', nDet * nEval, HElement_t_sizeB, tagCK)
                CK = (0.0_dp)
!C..
                allocate(W(nEval), stat=ierr)
                LogAlloc(ierr, 'W', nEval, 8, tagW)
                W = 0.0_dp
            end if
!C..
            IF (TREAD) THEN
                CALL READ_PSI(BOX, BOA, COA, NDET, NEVAL, NBASISMAX, NEL, CK, W)
            end if
        end if

!      TMC=TCALCHMAT.AND.(.NOT.TENERGY)

    End Subroutine DetCalcInit

    Subroutine DoDetCalc
        Use global_utilities
        use util_mod, only: get_free_unit
        use Determinants, only: get_helement, FDet, DefDet, tDefineDet
        use SystemData, only: Alat, arr, brr, boa, box, coa, ecore, g1, Beta
        use SystemData, only: t_new_real_space_hubbard
        use SystemData, only: nBasis, nBasisMax, nEl, nMsh, LzTot, TSPN, LMS
        use IntegralsData, only: FCK, NMAX, UMat
        Use LoggingData, only: iLogging, tLogDets, tCalcVariationalEnergy
        use Parallel_neci, only: iProcIndex
        use DetBitops, only: DetBitEQ, EncodeBitDet, FindBitExcitLevel
        use bit_rep_data, only: nifd, NIfTot, NIfD
        use bit_reps, only: decode_bit_det
        use sym_mod
        use HElem
        use MemoryManager, only: TagIntType
        use hist_data, only: tHistSpawn
        use CalcData, only: tFCIDavidson
        !use davidson_neci, only: DavidsonCalcType, DestroyDavidsonCalc
        !use davidson_neci, only: perform_davidson
        !use sparse_arrays, only: calculate_sparse_hamiltonian, deallocate_sparse_ham, &
        !sparse_ham, hamil_diag, HDiagTag, SparseHamilTags
        !use hamiltonian_linalg, only: sparse_hamil_type
        !use constants, only: size_n_int

        use davidson_neci, only: DavidsonCalcType, DestroyDavidsonCalc
        use davidson_neci, only: davidson_direct_ci_init, davidson_direct_ci_end, perform_davidson
        use hamiltonian_linalg, only: direct_ci_type, tCalcHFIndex
        use FCIMCData, only: davidson_ras, davidson_classes
        use ras, only: generate_entire_ras_space
        use real_space_hubbard, only: init_real_space_hubbard

        real(dp), ALLOCATABLE :: TKE(:), A(:, :), V(:), AM(:), BM(:), T(:), WT(:), SCR(:), WH(:), WORK2(:), V2(:, :), FCIGS(:)
        HElement_t(dp), ALLOCATABLE :: WORK(:)
        INTEGER, ALLOCATABLE :: INDEX(:), ISCR(:), Temp(:)
        integer(TagIntType) :: TKETag = 0, ATag = 0, VTag = 0, AMTag = 0, BMTag = 0, TTag = 0
        INTEGER(TagIntType) :: WTTag = 0, SCRTag = 0, ISCRTag = 0, INDEXTag = 0, WHTag = 0, Work2Tag = 0, V2Tag = 0, WorkTag = 0
        integer :: ierr, Lz
        character(25), parameter :: this_routine = 'DoDetCalc'
        real(dp) EXEN, GSEN
        Type(BasisFn) ISym, IHFSYM
        INTEGER GC, I, ICMAX
        INTEGER IN, IND, INDZ
        INTEGER NBLOCK!,OpenOrbs,OpenOrbsSym,Ex(2,NEl)
        INTEGER nKry1
        INTEGER(KIND=n_int) :: ilut(0:NIfTot), ilut_temp(0:NIfTot)
        INTEGER J, JR, iGetExcitLevel_2, ExcitLevel, iunit
        INTEGER LSCR, LISCR, MaxIndex
        LOGICAL tMC!,TestClosedShellDet,Found,tSign
        real(dp) :: GetHElement, calct, calcmcen, calcdlwdb, norm, temp_hel
        external :: GetHElement
        integer:: ic, TempnI(NEl), MomSymDet(NEl), ICSym, ICConnect, PairedUnit, SelfInvUnit
        integer(n_int) :: iLutMomSym(0:NIfTot)
        logical :: tSuccess
        type(DavidsonCalcType) :: davidsonCalc
        integer(n_int), allocatable, dimension(:, :) :: davidson_ilut
        integer, allocatable, dimension(:) :: davidson_parities
        integer :: nI(nel)
        integer :: davidson_size
        !character (len=*), parameter :: t_r = "DoDetCalc"
        !integer(TagIntType) :: IlutTag

        IF (tEnergy) THEN
            write(stdout, '(1X,A,E19.3)') ' B2LIMIT : ', B2L
            write(stdout, *) ' NBLK : ', NBLK
            write(stdout, *) ' NKRY : ', NKRY
            write(stdout, *) ' NEVAL : ', NEVAL

            write(stdout, *) ' NCYCLE : ', NCYCLE
            write(stdout, *) ' TENERGY : ', TENERGY
            write(stdout, *) ' IOBS : ', IOBS
            write(stdout, *) ' JOBS : ', JOBS
            write(stdout, *) ' KOBS : ', KOBS
            write(stdout, *) ' NMSH : ', NMSH
            IF (IOBS > NMSH .OR. IOBS <= 0 .OR. JOBS > NMSH .OR. JOBS <= 0 .OR. KOBS > NMSH .OR. KOBS <= 0) THEN
                call stop_all(this_routine, ' !!! REFERENCE PARTICLE NOT IN BOX !!! ')
            end if
        end if

!C.. now back to the storing H
        IF (TCALCHMAT) THEN
            if (t_new_real_space_hubbard) then
                call init_real_space_hubbard()
            end if
            write(stdout, *) "Calculating H matrix"
!C..We need to measure HAMIL and LAB first
            allocate(NROW(NDET), stat=ierr)
            CALL LogMemAlloc('NROW', NDET, 4, this_routine, NROWTag, ierr)
            NROW(1:NDET) = 0
            ICMAX = 1
!Falsify tMC
            TMC = .FALSE.
            allocate(HAMIL(0), LAB(0))
            CALL DETHAM(NDET, NEL, NMRKS, HAMIL, LAB, NROW, .TRUE., ICMAX, GC, TMC)
            deallocate(HAMIL, LAB)
            write(stdout, *) ' FINISHED COUNTING '
            write(stdout, *) "Allocating memory for hamiltonian: ", GC * 2
            CALL neci_flush(6)
!C..Now we know size, allocate memory to HAMIL and LAB
            LENHAMIL = GC
            allocate(Hamil(LenHamil), stat=ierr)
            LogAlloc(ierr, 'HAMIL', LenHamil, HElement_t_sizeB, tagHamil)
            HAMIL = (0.0_dp)
!C..
            allocate(LAB(LENHAMIL), stat=ierr)
            CALL LogMemAlloc('LAB', LenHamil, 4, this_routine, LabTag, ierr)

            LAB(1:LENHAMIL) = 0
!C..Now we store HAMIL and LAB
            CALL DETHAM(NDET, NEL, NMRKS, HAMIL, LAB, NROW, .FALSE., ICMAX, GC, TMC)

            IF (BTEST(ILOGGING, 7)) THEN
!C.. we write out H now
                iunit = get_free_unit()
                open(iunit, FILE='HAMIL', STATUS='UNKNOWN')
                J = 0
                JR = 0
!C            HMAX=-dflmax()
!C            HMIN=dflmax()
                DO I = 1, LENHAMIL
                    DO WHILE (I > J)
                        JR = JR + 1
                        J = J + NROW(JR)
                    end do
                    write(iunit, "(2I12)", advance='no') JR, LAB(I)
                    IF (HElement_t_size == 1) THEN
                        write(iunit, *) HAMIL(I)
                    ELSE
                        write(iunit, *) HAMIL(I), ABS(HAMIL(I))
                    end if
!C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
!C               write(14,"(A)",advance='no'),"|"
!C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
!C              write(14,"(F27.20)") HAMIL(I)
!C               CALL WRITEDET(14,NMRKS(1,LAB(I)),NEL,.FALSE.)
!C               write(14,"(A)",advance='no'),"|"
!C               CALL WRITEDET(14,NMRKS(1,JR),NEL,.FALSE.)
!C               write(14,"(F27.20)") HAMIL(I)

!C               IF(HAMIL(I).GT.HMAX) HMAX=HAMIL(I)
!C               IF(HAMIL(I).LT.HMIN) HMIN=HAMIL(I)
                end do
                close(iunit)
            end if
            temp_hel = real(GETHELEMENT(IFDET, IFDET, HAMIL, LAB, NROW, NDET), dp)
            write(stdout, *) '<D0|H|D0>=', temp_hel
            write(stdout, *) '<D0|T|D0>=', CALCT(NMRKS(1, IFDET), NEL)
            CALL neci_flush(6)
!CC         CALL HAMHIST(HMIN,HMAX,LENHAMIL,NHISTBOXES)
        end if
!C.. We've now finished calculating H if we were going to.
!C.. IF ENERGY CALC (for which we need to have calced H)
        IF (TENERGY) THEN
            IF (NBLK /= 0) THEN
!C..Things needed for Friesner-Pollard diagonalisation
                IF (TMC) call stop_all(this_routine, 'TMC and TENERGY set - Stopping')
                IF (HElement_t_size /= 1) call stop_all(this_routine, 'Cannot do Lanczos on Complex orbitals.')
                NKRY1 = NKRY + 1
                NBLOCK = MIN(NEVAL, NBLK)
                LSCR = MAX(NDET * NEVAL, 8 * NBLOCK * NKRY)
                LISCR = 6 * NBLOCK * NKRY
!C..
!            write(stdout,'(/,/,8X,64(1H*))')
                write(stdout, '(7X," *",62X,"*")')
                write(stdout, '(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
                write(stdout, '(7X," *",62X,"*")')
!            write(stdout,'(7X,1X,64(1H*))')
!C..Set up memory for FRSBLKH

                allocate(A(NEVAL, NEVAL), stat=ierr)
                CALL LogMemAlloc('A', NEVAL**2, 8, this_routine, ATag, ierr)
                A = 0.0_dp
!C..
!C,, W is now allocated with CK
!C..
                allocate(V(NDET * NBLOCK * NKRY1), stat=ierr)
                CALL LogMemAlloc('V', NDET * NBLOCK * NKRY1, 8, this_routine, VTag, ierr)
                V = 0.0_dp
!C..
                allocate(AM(NBLOCK * NBLOCK * NKRY1), stat=ierr)
                CALL LogMemAlloc('AM', NBLOCK * NBLOCK * NKRY1, 8, this_routine, AMTag, ierr)
                AM = 0.0_dp
!C..
                allocate(BM(NBLOCK * NBLOCK * NKRY), stat=ierr)
                CALL LogMemAlloc('BM', NBLOCK * NBLOCK * NKRY, 8, this_routine, BMTag, ierr)
                BM = 0.0_dp
!C..
                allocate(T(3 * NBLOCK * NKRY * NBLOCK * NKRY), stat=ierr)
                CALL LogMemAlloc('T', 3 * NBLOCK * NKRY * NBLOCK * NKRY, 8, this_routine, TTag, ierr)
                T = 0.0_dp
!C..
                allocate(WT(NBLOCK * NKRY), stat=ierr)
                CALL LogMemAlloc('WT', NBLOCK * NKRY, 8, this_routine, WTTag, ierr)
                WT = 0.0_dp
!C..
                allocate(SCR(LScr), stat=ierr)
                CALL LogMemAlloc('SCR', LScr, 8, this_routine, SCRTag, ierr)
                SCR = 0.0_dp
                allocate(ISCR(LIScr), stat=ierr)
                CALL LogMemAlloc('IScr', LIScr, 4, this_routine, IScrTag, ierr)
                ISCR(1:LISCR) = 0
                allocate(INDEX(NEVAL), stat=ierr)
                CALL LogMemAlloc('INDEX', NEVAL, 4, this_routine, INDEXTag, ierr)
                INDEX(1:NEVAL) = 0
!C..
                allocate(WH(NDET), stat=ierr)
                CALL LogMemAlloc('WH', NDET, 8, this_routine, WHTag, ierr)
                WH = 0.0_dp
                allocate(WORK2(3 * NDET), stat=ierr)
                CALL LogMemAlloc('WORK2', 3 * NDET, 8, this_routine, WORK2Tag, ierr)
                WORK2 = 0.0_dp
                allocate(V2(NDET, NEVAL), stat=ierr)
                CALL LogMemAlloc('V2', NDET * NEVAL, 8, this_routine, V2Tag, ierr)
                V2 = 0.0_dp
!C..Lanczos iterative diagonalising routine
                if (t_non_hermitian_2_body) then
                    call stop_all(this_routine, &
                                  "NECI_FRSBLKH not adapted for non-hermitian Hamiltonians")
                end if
                CALL NECI_FRSBLKH(NDET, ICMAX, NEVAL, HAMIL, LAB, CK, CKN, NKRY, NKRY1, NBLOCK, NROW, LSCR, LISCR, A, W, V, AM, BM, T, WT, &
         &  SCR, ISCR, INDEX, NCYCLE, B2L, .true., .false., .true.)

!Multiply all eigenvalues by -1.
                CALL DSCAL(NEVAL, -1.0_dp, W, 1)
            ELSE
!C.. We splice in a non-Lanczos diagonalisin routine if NBLOCK=0
                IF (NEVAL /= NDET) THEN
                    write(stdout, *) "NEVAL.NE.NDET.", NEVAL, NDET, " Cannot exactly diagonalize."
                    call stop_all(this_routine, "Cannot exactly diagonalise")
                end if
                write(stdout, *) "NBLK=0.  Doing exact diagonalization."
                IF (TCALCHMAT) THEN
                    allocate(WORK(4 * NDET), stat=ierr)
                    CALL LogMemAlloc('WORK', 4 * NDET, 8 * HElement_t_size, this_routine, WorkTag, ierr)
                    allocate(WORK2(3 * NDET), stat=ierr)
                    CALL LogMemAlloc('WORK2', 3 * NDET, 8, this_routine, WORK2Tag, ierr)
                    if (t_non_hermitian_2_body) then
                        call stop_all(this_routine, &
                                      "HDIAG_nec is not setup for non-hermitian Hamiltonians")
                    end if
                    CALL HDIAG_neci(NDET, HAMIL, LAB, NROW, CK, W, WORK2, WORK, NBLOCKSTARTS, NBLOCKS)
                end if
            end if
!C..
!  Since we no longer use HAMIL or LAB, we deallocate
            if (.not. tCalcVariationalEnergy) then
                LogDealloc(tagHamil)
                Deallocate(Hamil)
                DEallocate(LAB)
                CALL LogMemDealloc(this_routine, LabTag)
            end if
            allocate(TKE(NEVAL), stat=ierr)
            CALL LogMemAlloc('TKE', NEVAL, 8, this_routine, TKETag, ierr)

            EXEN = CALCMCEN(NEVAL, W, BETA)
            write(stdout, "(A,F19.9)") "EXACT E(BETA)=", EXEN
            GSEN = CALCDLWDB(IFDET, NDET, NEVAL, CK, W, BETA)
            write(stdout, "(A,F19.9)") "EXACT DLWDB(D0)=", GSEN
            write(stdout, "(A,F19.9)") "GROUND E=", W(1)
!C.. END ENERGY CALC
!      end if
        else if (tFCIDavidson) THEN
            if (.not. TSPN .or. LMS /= 0) then
                call stop_all("DoDetCalc", "FCI-Davidson only works for closed shell systems.")
            end if
            davidsonCalc = davidson_direct_ci_init()
            davidson_size = davidsonCalc%super%space_size
            allocate(davidson_ilut(0:NIfTot, davidson_size))
            allocate(davidson_parities(davidson_size))
            call generate_entire_ras_space(davidson_ras, davidson_classes, davidson_size, davidson_ilut, davidson_parities)
            !Find HF index
            !Set this flag, otherwise hfindex will be overwritten
            tCalcHFIndex = .false.
            davidsonCalc%super%hfindex = 0
            CALL EncodeBitDet(FDet, iLut(0:nifd))
            do i = 1, davidson_size
                if (DetBitEq(davidson_ilut(:, i), ilut)) then
                    davidsonCalc%super%hfindex = i
                    exit
                end if
            end do
            IF (davidsonCalc%super%hfindex == 0) call stop_all("DoDetCalc", "Fermi determinant is not found in RAS space!")
            if (t_non_hermitian_2_body) then
                call stop_all(this_routine, &
                              "perform_davidson not adapted for non-hermitian Hamiltonians!")
            end if
            call perform_davidson(davidsonCalc, direct_ci_type, .true.)
        end if

        call neci_flush(6)
!C.. If we're calculating rhos (for which we have to have calced H
!No longer used
!      IF(TRHOIJ) THEN
!         IF((.NOT.TENERGY).AND.(.NOT.TREADRHO)) THEN
!            write(stdout,*) "Calculating approx RHOs"
!            write(stdout,*) "Using Trotter decomposition? ",TTROT
!            write(stdout,*) "Order of Taylor: ",ABS(NTAY)
!            CALL CALCAPPROXRHOIJ(BETA,I_P,HAMIL,LAB,NROW,NDET,RHOMIN,RHOMAX,NRHOS,RHOEPS,TTROT,NTAY)
!         end if
!      end if
!C..Free HAMIL AND LAB memory if we no longer need them
!      IF(TCALCHMAT.AND..NOT.(TMONTE.AND.TMC)) THEN
!      end if

!C.. IF we want to compress the found determinants for use later...
        IF (tFindDets) THEN
            IF (tCompressDets) THEN
!Need to find symmetry of the reference determinant, so that we can only look for determinants of the correct symmetry.
                CALL GETSYM(FDET, NEL, G1, NBASISMAX, IHFSYM)
                IF (.not. associated(NMRKS)) THEN
                    write(stdout, *) "NMRKS not allocated"
                    CALL neci_flush(6)
                    CALL Stop_All("DoDetCalc", "NMRKS not allocated so cannot compress dets.")
                end if
!First, we want to count the number of determinants of the correct symmetry...
                Det = 0
                norm = 0.0_dp
                if (tFCIDavidson) then
                    do i = 1, davidson_size
                        CALL decode_bit_det(nI, davidson_ilut(:, i))
                        CALL GETSYM(nI, NEL, G1, NBASISMAX, ISYM)
                        !CALL GetLz(nI,NEL,Lz)
                        IF (ISym%Sym%S == IHFSYM%Sym%S) THEN
                            Det = Det + 1
                            norm = norm + (davidsonCalc%davidson_eigenvector(i))**2
                        end if
                    end do
                else
                    do i = 1, NDET
                        CALL GETSYM(NMRKS(:, i), NEL, G1, NBASISMAX, ISYM)
                        IF (ISym%Sym%S == IHFSYM%Sym%S) THEN
                            Det = Det + 1
                            IF (tEnergy) THEN
                                norm = norm + (REAL(CK(i, 1), dp))**2
                            end if
                        end if
                    end do
                end if
                write(stdout, "(I25,A,I4,A)") Det, " determinants of symmetry ", IHFSym%Sym%S, " found."
                write(stdout, *) "Normalization of eigenvector 1 is: ", norm
                CALL neci_flush(6)

                allocate(FCIDets(0:NIfTot, Det), stat=ierr)
                IF (ierr /= 0) CALL Stop_All("DetCalc", "Cannot allocate memory to hold vector")
                allocate(FCIDetIndex(0:NEl + 1), stat=ierr) !+1 so we can store the end of the array too
                allocate(Temp(Det), stat=ierr)
                IF (ierr /= 0) CALL Stop_All("DetCalc", "Cannot allocate memory to hold vector")
                IF (tEnergy .or. tFCIDavidson) THEN
                    allocate(FCIGS(Det), stat=ierr)
                    IF (ierr /= 0) CALL Stop_All("DetCalc", "Cannot allocate memory to hold vector")
                end if
                if (tCalcVariationalEnergy) then
                    !This allows us to resort to get back to the hamiltonian ordering
                    allocate(ReIndex(Det), stat=ierr)
                    if (ierr /= 0) CALL Stop_All("DetCalc", "Cannot allocate memory to hold vector")
                    ReIndex(:) = 0
                end if

                Det = 0
                FCIDetIndex(:) = 0
                if (tFCIDavidson) then
                    do i = 1, davidson_size
                        CALL decode_bit_det(nI, davidson_ilut(:, i))
                        CALL GETSYM(nI, NEL, G1, NBASISMAX, ISYM)
                        !CALL GetLz(nI,NEL,Lz)
                        !IF((ISym%Sym%S.eq.IHFSYM%Sym%S).and.(Lz.eq.LzTot)) THEN
                        IF (ISym%Sym%S == IHFSYM%Sym%S) THEN
                            Det = Det + 1
                            ExcitLevel = FindBitExcitLevel(davidson_ilut(:, i), iLut) !iGetExcitLevel_2(FDet,nI,NEl,NEl)
                            FCIDetIndex(ExcitLevel + 1) = FCIDetIndex(ExcitLevel + 1) + 1
                            Temp(Det) = ExcitLevel    !Temp will now temporarily hold the excitation level of the determinant.
                            FCIDets(:, Det) = davidson_ilut(:, i)
                            FCIGS(Det) = davidson_parities(i) * davidsonCalc%davidson_eigenvector(i) / norm
                        end if
                        if (tCalcVariationalEnergy) ReIndex(i) = i
                    end do
                else
                    do i = 1, NDet
                        CALL GETSYM(NMRKS(:, i), NEL, G1, NBASISMAX, ISYM)
                        IF (ISym%Sym%S == IHFSYM%Sym%S) THEN
                            Det = Det + 1
                            ExcitLevel = iGetExcitLevel_2(FDet, NMRKS(:, i), NEl, NEl)
                            ! FCIDetIndex is off by one, for later cumulative indexing
                            FCIDetIndex(ExcitLevel + 1) = FCIDetIndex(ExcitLevel + 1) + 1
                            Temp(Det) = ExcitLevel    !Temp will now temporarily hold the excitation level of the determinant.
                            CALL EncodeBitDet(NMRKS(:, i), FCIDets(0:NIfTot, Det))
                            IF (tEnergy) THEN
                                FCIGS(Det) = REAL(CK(i, 1), dp) / norm
                            end if
                        end if
                        if (tCalcVariationalEnergy) ReIndex(i) = i
                    end do
                end if

                IF (iExcitLevel <= 0) THEN
                    MaxIndex = NEl
                ELSE
                    MaxIndex = MIN(iExcitLevel, NEl)
                end if
!NB FCIDetIndex is off by 1 for later cumulation
                do i = 1, MaxIndex
                    write(stdout, *) "Number at excitation level: ", i, " is: ", FCIDetIndex(i + 1)
                end do

                ! We now want to sort the determinants according to the
                ! excitation level (stored in Temp)
                IF (.not. tEnergy .and. .not. tFCIDavidson) THEN
                    call sort(temp(1:Det), FCIDets(:, 1:Det))
                ELSE
                    if (tCalcVariationalEnergy) then
                        call sort(temp(1:Det), FCIDets(:, 1:Det), FCIGS(1:Det), ReIndex(1:Det))
                    else
                        call sort(temp(1:Det), FCIDets(:, 1:Det), FCIGS(1:Det))
                    end if
                end if

!Test that HF determinant is the first determinant
                CALL EncodeBitDet(FDet, iLut(0:nifd))
                do i = 0, nifd
                    IF (iLut(i) /= FCIDets(i, 1)) THEN
                        CALL Stop_All("DetCalc", "Problem with ordering the determinants by excitation level")
                    end if
                end do

!Change it so that FCIDetIndex indexes the start of the block of that excitation level.
                FCIDetIndex(0) = 1
                do i = 1, NEl + 1
                    FCIDetIndex(i) = FCIDetIndex(i - 1) + FCIDetIndex(i)
                end do
                FCIDetIndex(nEl + 1) = Det + 1
                IF (FCIDetIndex(nEl + 1) /= Det + 1) THEN
                    write(stdout, *) "FCIDETIndex:", FCIDetIndex(:)
                    CALL Stop_All("DetCalc", "Error in the indexing of determinant excitation level")
                end if

!We now need to sort within the excitation level by the "number" of the determinant
                do i = 1, MaxIndex
                    IF (.not. tEnergy .and. .not. tFCIDavidson) THEN
                        call sort(FCIDets(:, FCIDetIndex(i):FCIDetIndex(i + 1) - 1), &
                                  temp(FCIDetIndex(i):FCIDetIndex(i + 1) - 1))
                    ELSE
                        if (tCalcVariationalEnergy) then
                            call sort(FCIDets(:, FCIDetIndex(i):FCIDetIndex(i + 1) - 1), &
                                      temp(FCIDetIndex(i):FCIDetIndex(i + 1) - 1), &
                                      FCIGS(FCIDetIndex(i):FCIDetIndex(i + 1) - 1), &
                                      ReIndex(FCIDetIndex(i):FCIDetIndex(i + 1) - 1))
                        else
                            call sort(FCIDets(:, FCIDetIndex(i):FCIDetIndex(i + 1) - 1), &
                                      temp(FCIDetIndex(i):FCIDetIndex(i + 1) - 1), &
                                      FCIGS(FCIDetIndex(i):FCIDetIndex(i + 1) - 1))
                        end if
                    end if
                end do

                IF (tEnergy .or. tFCIDavidson) THEN
                    IF (tLogDETS .and. iProcIndex == 0) THEN
                        iunit = get_free_unit()
                        open(iunit, FILE='SymDETS', STATUS='UNKNOWN')

                        do i = 1, Det
                            write(iunit, "(2I17)", advance='no') i, temp(i)
                            do j = 0, nifd
                                write(iunit, "(I17)", advance='no') FCIDets(j, i)
                            end do
                            write(iunit, "(A,G25.16,A)", advance='no') " ", FCIGS(i), "  "
                            Call WriteBitDet(iunit, FCIDets(:, i), .true.)
                        end do
                        close(iunit)
                    end if
                    DEallocate(FCIGS)
                ELSE
                    IF (tLogDETS .and. iProcIndex == 0) THEN
                        iunit = get_free_unit()
                        open(iunit, FILE='SymDETS', STATUS='UNKNOWN')
                        write(iunit, *) "FCIDETIndex: ", FCIDetIndex(:)
                        write(iunit, *) "***"
                        do i = 1, Det
                            write(iunit, "(2I13)", advance='no') i, temp(i)
                            do j = 0, nifd
                                write(iunit, "(I13)", advance='no') FCIDets(j, i)
                            end do
                            write(iunit, "(A)", advance='no') " "
                            Call WriteBitDet(iunit, FCIDets(:, i), .true.)
                        end do
                        close(iunit)
                    end if
                end if !tEnergy - for dumping compressed ordered GS wavefunction
                DEallocate(Temp)
            end if !tCompressDets
        end if !tFindDets
!C..
        IF (TEnergy) THEN
            CALL CFF_CHCK(NDET, NEVAL, NMRKS, NEL, G1, CK, TKE)
            IF (BTEST(ILOGGING, 7)) CALL WRITE_PSI(BOX, BOA, COA, NDET, NEVAL, NBASISMAX, NEL, CK, W)
            IF (BTEST(ILOGGING, 8)) CALL WRITE_PSI_COMP(BOX, BOA, COA, NDET, NEVAL, NBASISMAX, NEL, CK, W)
            write(stdout, *) '       ====================================================== '
            write(stdout, '(A5,5X,A15,1X,A18,1x,A20)') 'STATE', 'KINETIC ENERGY', 'COULOMB ENERGY', 'TOTAL ENERGY'
            iunit = get_free_unit()
            open(iunit, FILE='ENERGIES', STATUS='UNKNOWN')
            DO IN = 1, NEVAL
                write(stdout, '(I5,2X,3(F19.11,2x))') IN, TKE(IN), W(IN) - TKE(IN), W(IN)
                !            write(iunit,"(I7)",advance='no') IN
                !            CALL WRITEDET(iunit,NMRKS(1,IN),NEL,.FALSE.)
                write(iunit, "(F19.11)") W(IN)
            end do
            close(iunit)
            write(stdout, *) '       ====================================================== '
        else if (tFCIDavidson) THEN
            open(iunit, FILE='ENERGIES', STATUS='UNKNOWN')
            write(iunit, "(F19.11)") davidsonCalc%davidson_eigenvalue
            close(iunit)
!C., END energy calc
        end if

        IF (tFCIDavidson) THEN
            call davidson_direct_ci_end(davidsonCalc)
            call DestroyDavidsonCalc(davidsonCalc)
            DEallocate(davidson_ilut)
            DEallocate(davidson_parities)
        end if
!C.. Jump to here if just read Psi in
        CONTINUE

        !Now deallocate NMRKS if tFindDets and not tEnergy
        if (tFindDets .and. tCompressDets .and. (.not. tEnergy)) then
            DEallocate(NMRKS)
            CALL LogMemDealloc(this_routine, tagNMRKS)
        end if

    End Subroutine DoDetCalc

END MODULE DetCalc

!  Given an exact calculation of eigen-vectors and -values, calculate the expectation value of E(Beta)
FUNCTION CALCMCEN(NEVAL, W, BETA)
    use constants, only: dp
    IMPLICIT NONE
    INTEGER NEVAL, IK
    real(dp) W(NEVAL), BETA, DNORM, EN, CALCMCEN
    EN = 0.0_dp
    DNORM = 0.0_dp
    DO IK = 1, NEVAL
        EN = EN + (W(IK)) * EXP(-(W(IK) - W(1)) * BETA)
        DNORM = DNORM + EXP(-(W(IK) - W(1)) * BETA)
    end do
    CALCMCEN = EN / DNORM
    RETURN
END

!  Given an exact calculation of eigen-vectors and -values, calculate the expectation value of E~(Beta)_I for det I
FUNCTION CALCDLWDB(I, NDET, NEVAL, CK, W, BETA)
    use constants, only: dp
    IMPLICIT NONE
    INTEGER NDET, NEVAL, IK, I
    HElement_t(dp) CK(NDET, NEVAL)
    real(dp) W(NEVAL), BETA, DNORM, EN, CALCDLWDB
    EN = 0.0_dp
    DNORM = 0.0_dp
    DO IK = 1, NEVAL
        EN = EN + abs(CK(I, IK))**2 * (W(IK)) * EXP(-(W(IK) - W(1)) * BETA)
        DNORM = DNORM + abs(CK(I, IK))**2 * EXP(-(W(IK) - W(1)) * BETA)
    end do
    CALCDLWDB = EN / DNORM
    RETURN
END

SUBROUTINE CFF_CHCK(NDET, NEVAL, NM, NEL, G1, CG, TKE)
    use constants, only: dp
    use util_mod, only: get_free_unit
    USE OneEInts, only: GetTMATEl
    use SystemData, only: BasisFN
    use HElem
    IMPLICIT NONE
    INTEGER NEL, NM(NEL, *), NDET, NEVAL, iunit
    HElement_t(dp) CG(NDET, NEVAL)
    real(dp) TKE(NEVAL)
    TYPE(BASISFN) G1(*)
    real(dp) PI, S, SUM1
    real(dp) AUX
    INTEGER I, J, IN, IEL, L
!..Calculate the expectation value of the kinetic energy
!..<Psi|T|Psi>
    PI = ACOS(-1.0_dp)
    DO IN = 1, NEVAL
        TKE(IN) = 0.0_dp
        DO I = 1, NDET
            SUM1 = 0.0_dp
            DO J = 1, NEL
                AUX = GetTMATEl(NM(J, I), NM(J, I))
!((ALAT(1)**2)*((G1(1,NM(J,I))**2)/(ALAT(1)**2)+
!     &     (G1(2,NM(J,I))**2)/(ALAT(2)**2)+
!     &     (G1(3,NM(J,I))**2)/(ALAT(3)**2)))
                SUM1 = SUM1 + (AUX)
            end do
!..Cube multiplier
!          CST=PI*PI/(2.0_dp*ALAT(1)*ALAT(1))
!.. Deal with the UEG
!          IF(NBASISMAX(1,1).LE.0) CST=CST*4.0_dp
!          SUM1=CST*SUM1
            TKE(IN) = TKE(IN) + SUM1 * abs(CG(I, IN))**2
        end do
    end do
! ==--------------------------------------------------------------==
    IF (.FALSE.) THEN
!      IF(BTEST(ILOGGING,7)) THEN
        iunit = get_free_unit()
        open(iunit, FILE='PSI', STATUS='UNKNOWN')
        DO J = 1, NEVAL
            IF (J == 1) THEN
                write(iunit, *) ' GROUND STATE COEFFICIENTS  '
            ELSE
                write(iunit, *) ' COEFFICIENTS FOR EXCITED STATE NUMBER : ', J
            end if
            S = 0.0_dp
            DO I = 1, NDET
                IF (abs(CG(I, J)) > 1.0e-15_dp) THEN
                    DO IEL = 1, NEL
                        write(iunit, "(I3,I3,2I3,2X)", advance='no') (G1(NM(1, IEL))%K(L), L=1, 3)
                    end do
                    IF (HElement_t_size == 1) THEN
                        write(iunit, "(F19.9,1X,I7)") CG(I, J), I
                    ELSE
                        write(iunit, "(F19.9,1X,I7)") CG(I, J), I
                    end if
                end if
                S = S + abs(CG(I, J))**2
            end do
            write(iunit, '(A,F19.5)') ' SQUARE OF COEFFICIENTS : ', S
            write(iunit, *)
        end do
        close(iunit)
    end if
    RETURN
END
