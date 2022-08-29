module CAS_distribution_init
    use SystemData, only: tHPHFInts, tHPHF, lms, lztot, t_non_hermitian, tspn

    use CalcData, only: DiagSft, InitialPart, InitWalkers, OccCasorbs, RealCoeffExcitThresh,&
                        tAllRealCoeff, tReadPops, tRealCoeffByExcitLevel, &
                        tRestartHighPop, tStartSinglePart, tTruncInitiator, &
                        VirtCASorbs

    use Determinants, only: get_helement

    use hphf_integrals, only: hphf_diag_helement

    use DeterminantData, only: write_det, write_det_len

    use DetCalcData, only: NKRY, NBLK, B2L, det, ncycle

    use bit_rep_data, only: NIfTot, nIfD

    use bit_reps, only: encode_det, clear_all_flags, encode_sign

    use hash, only: FindWalkerHash

    use load_balance_calcnodes, only: DetermineDetNode

    use DetBitOps, only: FindBitExcitLevel, TestClosedShellDet, &
                         IsAllowedHPHF, DetBitEq, &
                         EncodeBitDet

    use initial_trial_states, only: calc_trial_states_lanczos, &
                                    set_trial_populations, set_trial_states, calc_trial_states_direct
    use global_det_data, only: global_determinant_data, set_det_diagH, &
                               set_det_offdiagH, clean_global_det_data, &
                               init_global_det_data, set_spawn_rate, &
                               store_decoding
    use semi_stoch_gen, only: init_semi_stochastic, end_semistoch, &
                              enumerate_sing_doub_kpnt
    use semi_stoch_procs, only: return_mp1_amp_and_mp2_energy
    use initiator_space_procs, only: init_initiator_space
    use kp_fciqmc_data_mod, only: tExcitedStateKP
    use sym_general_mod, only: ClassCountInd
    use trial_wf_gen, only: init_trial_wf, end_trial_wf
    use load_balance, only: clean_load_balance, init_load_balance
    use matel_getter, only: get_diagonal_matel, get_off_diagonal_matel
    use ueg_excit_gens, only: gen_ueg_excit
    use gndts_mod, only: gndts
    use excit_gen_5, only: gen_excit_4ind_weighted2
    use tc_three_body_excitgen, only: gen_excit_mol_tc, setup_mol_tc_excitgen
    use pcpp_excitgen, only: gen_rand_excit_pcpp, init_pcpp_excitgen, finalize_pcpp_excitgen

    use tau_search, only: init_tau_search, max_death_cpt

    use fcimc_helper, only: CalcParentFlag, update_run_reference

    use cont_time_rates, only: spawn_rate_full, oversample_factors, &
                               secondary_gen_store, ostag

    use Parallel_neci, only: iProcIndex, nNodes, mpisumall

    use util_mod, only: operator(.isclose.)

    use FciMCData, only: ll_node, HFSym, ProjEDet, tSinglePartPhase, NoatHF, TotParts, &
        iLutRef, CurrentDets, OldAllHFCyc, AllTotParts, iter_data_fciqmc, AllNoAbortedOld, &
        AllNoatHf, AllTotPartsOld, AllTotWalkers, AllTotWalkersOld, Hii, &
        nWalkerHashes, OldAllNoatHF, TotPartsOld, TotWalkers, TotWalkersOld, &
        HashIndex, OldAllAvWalkersCyc

    use dSFMT_interface, only: genrand_real2_dSFMT

    use sym_mod

    use constants

    implicit none
    private
    public :: InitFCIMC_CAS

contains

    subroutine InitFCIMC_CAS()

        ! Routine to initialise the particle distribution according to a CAS diagonalisation.
        ! This hopefully will help with close-lying excited states of the same sym.

        type(BasisFN) :: CASSym
        integer :: i, ierr, nEval, NKRY1, NBLOCK, LSCR, LISCR, DetIndex
        integer :: iNode, nBlocks, nBlockStarts(2), DetHash
        integer :: CASSpinBasisSize, nCASDet, ICMax, GC, LenHamil, iInit
        integer :: nHPHFCAS, iCasDet, ExcitLevel
        real(dp) :: NoWalkers
        integer, allocatable :: CASBrr(:), CASDet(:), CASFullDets(:, :), nRow(:), Lab(:), ISCR(:), INDEX(:)
        integer, pointer :: CASDetList(:, :) => null()
        integer(n_int) :: iLutnJ(0:NIfTot)
        logical :: tMC, tHPHF_temp, tHPHFInts_temp
        HElement_t(dp) :: HDiagTemp, HOffDiagTemp
        HElement_t(dp), allocatable :: Hamil(:)
        real(dp), allocatable :: CK(:, :), W(:), CKN(:, :), A_Arr(:, :), V(:), BM(:), T(:), WT(:)
        real(dp), allocatable :: SCR(:), WH(:), Work2(:), V2(:, :), AM(:)
        real(dp), allocatable :: Work(:)
        integer(TagIntType) :: ATag = 0, VTag = 0, BMTag = 0, TTag = 0, WTTag = 0, SCRTag = 0, WHTag = 0, Work2Tag = 0, V2Tag = 0
        integer(TagIntType) :: ISCRTag = 0, IndexTag = 0, AMTag = 0
        integer(TagIntType) :: WorkTag = 0
        real(dp) :: CASRefEnergy, TotWeight, PartFac, amp, rat, r, GetHElement
        external :: GetHElement
        real(dp), dimension(lenof_sign) :: temp_sign
        real(dp) :: energytmp(nel), max_wt
        integer  :: tmp_det(nel), det_max, run
        type(ll_node), pointer :: TempNode
        character(len=*), parameter :: this_routine = 'InitFCIMC_CAS'
#ifdef CMPLX_
        call stop_all(this_routine, "StartCAS currently does not work with complex walkers")
#endif
        if (tReadPops) call stop_all(this_routine, "StartCAS cannot work with with ReadPops")
        if (tStartSinglePart) call stop_all(this_routine, "StartCAS cannot work with StartSinglePart")
        if (tRestartHighPop) call stop_all(this_routine, "StartCAS cannot with with dynamically restarting calculations")

        write(stdout, *) "Initialising walkers proportional to a CAS diagonalisation..."
        write(stdout, '(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (" &
            , (OccCASOrbs + VirtCASOrbs) / 2, ",", OccCASOrbs, ")"
        do I = NEl - OccCASorbs + 1, NEl
            write(stdout, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(stdout, G1(BRR(I))%SYM, .FALSE.)
            write(stdout, '(I4)', advance='no') G1(BRR(I))%Ml
            write(stdout, '(2F19.9)') ARR(I, 1), ARR(BRR(I), 2)
        end do
        write(stdout, '(A)') " ================================================================================================="
        do I = NEl + 1, NEl + VirtCASOrbs
            write(stdout, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(stdout, G1(BRR(I))%SYM, .FALSE.)
            write(stdout, '(I4)', advance='no') G1(BRR(I))%Ml
            write(stdout, '(2F19.9)') ARR(I, 1), ARR(BRR(I), 2)
        end do

        CASSpinBasisSize = OccCASorbs + VirtCASorbs
        allocate(CASBrr(1:CASSpinBasisSize))
        allocate(CASDet(1:OccCasOrbs))
        do i = 1, CASSpinBasisSize
            !Run through the cas space, and create an array which will map these orbtials to the
            !orbitals they actually represent.
            CASBrr(i) = BRR(i + (NEl - OccCasorbs))
        end do

        !Calculate symmetry of CAS determinants, and check that this will be the same as the reference determinant
        !for the rest of the FCIMC calculations.
        CASDet = CasBRR(1:OccCasOrbs)
        call sort(CasDet)

        write(stdout, *) "CAS Det is: "
        call write_det_len(stdout, CASDet, OccCASOrbs, .true.)
        call GetSym(CASDet, OccCASOrbs, G1, nBasisMax, CASSym)
        write(stdout, *) "Spatial symmetry of CAS determinants: ", CASSym%Sym%S
        write(stdout, *) "Ms of CAS determinants: ", CASSym%Ms
        if (tFixLz) then
            write(stdout, *) "Ml of CAS determinants: ", CASSym%Ml
        end if
        call neci_flush(stdout)

        if (CASSym%Ml /= LzTot) call stop_all(this_routine, "Ml of CAS ref det does not match Ml of full reference det")
        if (CASSym%Ms /= 0) call stop_all(this_routine, "CAS diagonalisation can only work with closed shell CAS spaces initially")
        if (CASSym%Sym%S /= HFSym%Sym%S) then
            call stop_all(this_routine, "Sym of CAS ref det does not match Sym of full reference det")
        end if

        !First, we need to generate all the excitations.
        call gndts(OccCASorbs, CASSpinBasisSize, CASBrr, nBasisMax, CASDetList, .true., G1, tSpn, LMS, .true., CASSym, nCASDet, iCASDet)

        if (nCASDet == 0) call stop_all(this_routine, "No CAS determinants found.")
        write(stdout, *) "Number of symmetry allowed CAS determinants found to be: ", nCASDet
        allocate(CASDetList(OccCASorbs, nCASDet), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating CASDetList")
        CASDetList(:, :) = 0

        !Now fill up CASDetList...
        call gndts(OccCASorbs, CASSpinBasisSize, CASBrr, nBasisMax, &
                   CASDetList, .false., G1, tSpn, LMS, .true., CASSym, &
                   nCASDet, iCASDet)

        !We have a complication here. If we calculate the hamiltonian from these CAS determinants, then we are not
        !including the mean-field generated from the other occupied orbitals. We need to either 'freeze' the occupied
        !orbitals and modify the 1 & two electron integrals, or add the other electrons back into the list. We do the latter.
        allocate(CASFullDets(NEl, nCASDet), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating CASFullDets")
        CASFullDets(:, :) = 0

        ! Get the first part of a determinant with the lowest energy, rather
        ! than lowest index number orbitals
        energytmp = ARR(ProjEDet(:, 1), 2)
        tmp_det = ProjEDet(:, 1)
        call sort(energytmp, tmp_det)

        ! Construct the determinants resulting from the CAS expansion.
        do i = 1, nCASDet
            CASFullDets(1:nel - OccCASorbs, i) = tmp_det(1:nel - OccCASOrbs)
            CASFullDets(nel - OccCASorbs + 1:nel, i) = CASDetList(1:OccCASorbs, i)
            call sort(CASFullDets(:, i))
        end do
        deallocate(CASDetList)

        write(stdout, *) "First CAS determinant in list is: "
        call write_det(stdout, CASFullDets(:, 1), .true.)

        if (nCASDet > 1300) then
            !Do lanczos
            nEval = 4
        else
            nEval = nCASDet
        end if
        write(stdout, "(A,I4,A)") "Calculating lowest ", nEval, " eigenstates of CAS Hamiltonian..."
        allocate(Ck(nCASDet, nEval), stat=ierr)
        Ck = 0.0_dp
        allocate(W(nEval), stat=ierr)    !Eigenvalues
        W = 0.0_dp
        if (ierr /= 0) call stop_all(this_routine, "Error allocating")

        write(stdout, *) "Calculating hamiltonian..."
        allocate(nRow(nCASDet), stat=ierr)
        nRow = 0
        ICMax = 1
        tMC = .false.

        !HACK ALERT!! Need to fill up array in space of determinants, not HPHF functions.
        !Turn off tHPHFInts and tHPHF and turn back on after the hamiltonian constructed.
        tHPHF_temp = tHPHF
        tHPHFInts_temp = tHPHFInts
        tHPHF = .false.
        tHPHFInts = .false.

        ! do not pass an unallocated array, so dummy-allocate
        allocate(Hamil(0))
        allocate(Lab(0))
        CALL Detham(nCASDet, NEl, CASFullDets, Hamil, Lab, nRow, .true., ICMax, GC, tMC)
        deallocate(Lab)
        deallocate(Hamil)
        LenHamil = GC
        write(stdout, *) "Allocating memory for hamiltonian: ", LenHamil * 2
        allocate(Hamil(LenHamil), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating Hamil")
        Hamil = 0.0_dp
        allocate(Lab(LenHamil), stat=ierr)
        if (ierr /= 0) call stop_all(this_routine, "Error allocating Lab")
        Lab = 0
        call Detham(nCASDet, NEl, CASFullDets, Hamil, Lab, nRow, .false., ICMax, GC, tMC)

        ! Assuming that this routine does not work with complex
        CASRefEnergy = GETHELEMENT(1, 1, HAMIL, LAB, NROW, NCASDET)
        write(stdout, *) "Energy of first CAS det is: ", CASRefEnergy

        ! Turn back on HPHFs if needed.
        tHPHF = tHPHF_temp
        tHPHFInts = tHPHFInts_temp

!        if(abs(CASRefEnergy-Hii).gt.1.0e-7_dp) then
!            call stop_all(this_routine,"CAS reference energy does not match reference energy of full space")
!        end if

        if (nCASDet > 1300) then
            !Lanczos
            NKRY1 = NKRY + 1
            NBLOCK = MIN(NEVAL, NBLK)
            LSCR = MAX(nCASDet * NEVAL, 8 * NBLOCK * NKRY)
            LISCR = 6 * NBLOCK * NKRY
            allocate(A_Arr(NEVAL, NEVAL), stat=ierr)
            CALL LogMemAlloc('A_Arr', NEVAL**2, 8, this_routine, ATag, ierr)
            A_Arr = 0.0_dp
            allocate(V(nCASDet * NBLOCK * NKRY1), stat=ierr)
            CALL LogMemAlloc('V', nCASDet * NBLOCK * NKRY1, 8, this_routine, VTag, ierr)
            V = 0.0_dp
            allocate(AM(NBLOCK * NBLOCK * NKRY1), stat=ierr)
            CALL LogMemAlloc('AM', NBLOCK * NBLOCK * NKRY1, 8, this_routine, AMTag, ierr)
            AM = 0.0_dp
            allocate(BM(NBLOCK * NBLOCK * NKRY), stat=ierr)
            CALL LogMemAlloc('BM', NBLOCK * NBLOCK * NKRY, 8, this_routine, BMTag, ierr)
            BM = 0.0_dp
            allocate(T(3 * NBLOCK * NKRY * NBLOCK * NKRY), stat=ierr)
            CALL LogMemAlloc('T', 3 * NBLOCK * NKRY * NBLOCK * NKRY, 8, this_routine, TTag, ierr)
            T = 0.0_dp
            allocate(WT(NBLOCK * NKRY), stat=ierr)
            CALL LogMemAlloc('WT', NBLOCK * NKRY, 8, this_routine, WTTag, ierr)
            WT = 0.0_dp
            allocate(SCR(LScr), stat=ierr)
            CALL LogMemAlloc('SCR', LScr, 8, this_routine, SCRTag, ierr)
            SCR = 0.0_dp
            allocate(ISCR(LIScr), stat=ierr)
            CALL LogMemAlloc('IScr', LIScr, 4, this_routine, IScrTag, ierr)
            ISCR(1:LISCR) = 0
            allocate(INDEX(NEVAL), stat=ierr)
            CALL LogMemAlloc('INDEX', NEVAL, 4, this_routine, INDEXTag, ierr)
            INDEX(1:NEVAL) = 0
            allocate(WH(nCASDet), stat=ierr)
            CALL LogMemAlloc('WH', nCASDet, 8, this_routine, WHTag, ierr)
            WH = 0.0_dp
            allocate(WORK2(3 * nCASDet), stat=ierr)
            CALL LogMemAlloc('WORK2', 3 * nCASDet, 8, this_routine, WORK2Tag, ierr)
            WORK2 = 0.0_dp
            allocate(V2(nCASDet, NEVAL), stat=ierr)
            CALL LogMemAlloc('V2', nCASDet * NEVAL, 8, this_routine, V2Tag, ierr)
            V2 = 0.0_dp
            allocate(CkN(nCASDet, nEval), stat=ierr)
            CkN = 0.0_dp
            !C..Lanczos iterative diagonalising routine
            if (t_non_hermitian) then
                call stop_all(this_routine, &
                              "NECI_FRSBLKH not adapted for non-hermitian Hamiltonians!")
            end if
         CALL NECI_FRSBLKH(nCASDet, ICMAX, NEVAL, HAMIL, LAB, CK, CKN, NKRY, NKRY1, NBLOCK, NROW, LSCR, LISCR, A_Arr, W, V, AM, BM, T, WT, &
         &  SCR, ISCR, INDEX, NCYCLE, B2L, .false., .false., .true.)
            !Multiply all eigenvalues by -1.
            CALL DSCAL(NEVAL, -1.0_dp, W, 1)
            if (CK(1, 1) < 0.0_dp) then
                do i = 1, nCASDet
                    CK(i, 1) = -CK(i, 1)
                end do
            end if

            deallocate(CKN, A_Arr, V, BM, T, WT, SCR, WH, V2, iscr, index, AM)
            call logmemdealloc(this_routine, ATag)
            call logmemdealloc(this_routine, VTag)
            call logmemdealloc(this_routine, BMTag)
            call logmemdealloc(this_routine, TTag)
            call logmemdealloc(this_routine, WTTag)
            call logmemdealloc(this_routine, SCRTag)
            call logmemdealloc(this_routine, WHTag)
            call logmemdealloc(this_routine, V2Tag)
            call logmemdealloc(this_routine, iscrTag)
            call logmemdealloc(this_routine, indexTag)
            call logmemdealloc(this_routine, AMTag)
        else
            !complete diagonalisation
            allocate(Work(4 * nCASDet), stat=ierr)
            call LogMemAlloc('Work', 4 * nCASDet, 8, this_routine, WorkTag, ierr)
            allocate(Work2(3 * nCASDet), stat=ierr)
            call logMemAlloc('Work2', 3 * nCASDet, 8, this_routine, Work2Tag, ierr)
            nBlockStarts(1) = 1
            nBlockStarts(2) = nCASDet + 1
            nBlocks = 1
            if (t_non_hermitian) then
                call stop_all(this_routine, &
                              "HDIAG_neci is not set up for non-hermitian Hamiltonians!")
            end if
            call HDIAG_neci(nCASDet, Hamil, Lab, nRow, CK, W, Work2, Work(1), nBlockStarts, nBlocks)
            deallocate(Work)
            call LogMemDealloc(this_routine, WorkTag)
        end if
        !Deallocate all the lanczos arrays now.
        deallocate(nrow, lab, work2)
        call logmemdealloc(this_routine, Work2Tag)

        write(stdout, *) "Diagonalisation complete. Lowest energy CAS eigenvalues/corr E are: "
        do i = 1, min(NEval, 10)
            write(stdout, *) i, W(i), W(i) - CASRefEnergy
        end do

        TotWeight = 0.0_dp
        nHPHFCAS = 0
        max_wt = 0
        det_max = 0
        do i = 1, nCASDet
            if (tHPHF) then
                !Only allow valid HPHF functions
                call EncodeBitDet(CASFullDets(:, i), iLutnJ)
                if (IsAllowedHPHF(iLutnJ)) then
                    nHPHFCAS = nHPHFCAS + 1
                    if (.not. TestClosedShellDet(iLutnJ)) then
                        !Open shell. Weight is sqrt(2) of det weight.
                        TotWeight = TotWeight + (abs(CK(i, 1)) * sqrt(2.0_dp))
                        !Return this new weight to the CK array, so that we do not need to do this a second time.
                        CK(i, 1) = CK(i, 1) * sqrt(2.0_dp)
                    else
                        !Closed Shell
                        TotWeight = TotWeight + abs(CK(i, 1))
                    end if
                end if
            else
                TotWeight = TotWeight + abs(CK(i, 1))
            end if

            ! Find the maximum weighted determinant.
            if (abs(ck(i, 1)) > max_wt) then
                max_wt = abs(ck(i, 1))
                det_max = i
            end if
        end do

        ! Output details
        write(stdout, *) "Total weight of lowest eigenfunction: ", TotWeight
        write(stdout, *) "Maximum weighted det: ", det_max, max_wt

        !If the reference det is not the maximum weighted det, suggest that
        !we change it!
        if (.not. all(CASFullDets(:, det_max) == ProjEDet(:, 1))) then
            write(stdout, *) 'The specified reference determinant is not the &
                       &maximum weighted determinant in the CAS expansion'
            write(stdout, *) 'Use following det as reference:'
            call write_det(6, CASFullDets(:, det_max), .true.)
            call warning_neci(this_routine, "Poor reference chosen")
        end if

        if (tHPHF) write(stdout, *) "Converting into HPHF space. Total HPHF CAS functions: ", nHPHFCAS

        if ((InitialPart.isclose.1._dp) .or. (InitialPart >= (InitWalkers * nNodes) - 50)) then
            !Here, all the walkers will be assigned to the CAS wavefunction.
            !InitialPart = 1 by default
            write(stdout, "(A)") "All walkers specified in input will be distributed according to the CAS wavefunction."
            write(stdout, "(A)") "Shift will be allowed to vary from the beginning"
            write(stdout, "(A,F20.9)") "Setting initial shift to equal CAS correlation energy", W(1) - CASRefEnergy
            DiagSft = W(1) - CASRefEnergy
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac = (real(InitWalkers, dp) * real(nNodes, dp)) / TotWeight
        else
            !Here, not all walkers allowed will be initialised to the CAS wavefunction.
            write(stdout, "(A,I15,A)") "Initialising ", int(InitialPart), " walkers according to the CAS distribution."
            write(stdout, "(A,I15)") "Shift will remain fixed until the walker population reaches ", int(InitWalkers * nNodes)
            !PartFac is the number of walkers that should reside on the HF determinant
            PartFac = real(InitialPart, dp) / TotWeight
            tSinglePartPhase(:) = .true.
        end if

        !Now generate all excitations again, creating the required number of walkers on each one.
        DetIndex = 1
        NoatHF(:) = 0.0
        TotParts(:) = 0.0
        do i = 1, nCASDet
            if (tHPHF) then
                call EncodeBitDet(CASFullDets(:, i), iLutnJ)
                if (.not. IsAllowedHPHF(iLutnJ)) cycle
            end if
            iNode = DetermineDetNode(nel, CASFullDets(:, i), 0)
            if (iProcIndex == iNode) then
                !Number parts on this det = PartFac*Amplitude
                amp = CK(i, 1) * PartFac

                if (tRealCoeffByExcitLevel) ExcitLevel = FindBitExcitLevel(iLutnJ, iLutRef, nEl)
                if (tAllRealCoeff .or. &
                    & (tRealCoeffByExcitLevel .and. (ExcitLevel <= RealCoeffExcitThresh))) then
                    NoWalkers = amp
                else
                    NoWalkers = int(amp)
                    rat = amp - real(NoWalkers, dp)
                    r = genrand_real2_dSFMT()
                    if (abs(rat) > r) then
                        if (amp < 0.0_dp) then
                            NoWalkers = NoWalkers - 1
                        else
                            NoWalkers = NoWalkers + 1
                        end if
                    end if
                end if

                if (abs(NoWalkers) > 1.0e-12_dp) then
                    call EncodeBitDet(CASFullDets(:, i), iLutnJ)
                    if (DetBitEQ(iLutnJ, iLutRef(:, 1), nifd)) then
                        !Check if this determinant is reference determinant, so we can count number on hf.
                        do run = 1, inum_runs
                            NoatHF(run) = NoWalkers
                        end do
                    end if
                    call encode_det(CurrentDets(:, DetIndex), iLutnJ)
                    call clear_all_flags(CurrentDets(:, DetIndex))
                    do run = 1, inum_runs
                        temp_sign(run) = NoWalkers
                    end do
                    call encode_sign(CurrentDets(:, DetIndex), temp_sign)

                    ! Store the diagonal and off-diagonal matrix elements
                    HDiagTemp = get_diagonal_matel(CASFullDets(:, i), iLutnJ)
                    HOffDiagTemp = get_off_diagonal_matel(CASFullDets(:, i), iLutnJ)
                    call set_det_diagH(DetIndex, real(HDiagTemp, dp) - Hii)
                    call set_det_offdiagH(DetIndex, HOffDiagTemp)
                    call store_decoding(DetIndex, CASFullDets(:, i))

                    if (tTruncInitiator) then
                        !Set initiator flag if needed (always for HF)
                        call CalcParentFlag(DetIndex, iInit)
                    end if

                    DetHash = FindWalkerHash(CASFullDets(:, i), nWalkerHashes)
                    TempNode => HashIndex(DetHash)
                    ! If the first element in the list has not been used.
                    if (TempNode%Ind == 0) then
                        TempNode%Ind = DetIndex
                    else
                        do while (associated(TempNode%Next))
                            TempNode => TempNode%Next
                        end do
                        allocate(TempNode%Next)
                        nullify (TempNode%Next%Next)
                        TempNode%Next%Ind = DetIndex
                    end if
                    nullify (TempNode)

                    DetIndex = DetIndex + 1
                    do run = 1, inum_runs
                        TotParts(run) = TotParts(run) + abs(NoWalkers)
                    end do
                end if
            end if   !End if desired node
        end do

        TotWalkers = DetIndex - 1   !This is the number of occupied determinants on each node
        TotWalkersOld = TotWalkers

        !Set local&global variables
        TotPartsOld = TotParts
        call mpisumall(TotParts, AllTotParts)
        call mpisumall(NoatHF, AllNoatHF)
        call mpisumall(TotWalkers, AllTotWalkers)
        OldAllNoatHF = AllNoatHF
        do run = 1, inum_runs
            OldAllHFCyc(run) = AllNoatHF(run)
            OldAllAvWalkersCyc(run) = AllTotParts(run)
        end do
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        iter_data_fciqmc%tot_parts_old = AllTotPartsOld
        AllNoAbortedOld = 0.0_dp

        deallocate(CK, W, Hamil, CASBrr, CASDet, CASFullDets)

    end subroutine InitFCIMC_CAS


end module CAS_distribution_init
