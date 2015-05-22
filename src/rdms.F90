module rdms

! This file contains the routines used to calculate 1- and 2-RDM estimates.
! This is done on the fly to avoid having to histogram the full wavefunction
! which is extremely time and memory inefficient. In this way, these routines
! differ slightly from those in NatOrbsMod (which take a histogrammed 
! wavefunction usually truncated around double excitations and form the one
! electron RDM) but the basic formula is still the same.

! For example, he elements of the one electron reduced density matrix are given
! by:
!
! 1RDM_pq   = < Psi | a_p+ a_q | Psi > 
!           = < sum_i c_i D_i | a_p+ a_q | sum_j c_j D_j >
!
! where |Psi> is the full wavefunction, a_p+ is the creation operator and a_q is
! the annihilation operator. The elements 1RDM_pq therefore come from the sum
! of the contributions c_i*c_j from all pairs of determinants D_i and D_j which
! are related by a single excitation between p and q.

use Global_Utilities
use Parallel_neci
use bit_reps, only: NIfTot, NIfDBO, decode_bit_det, extract_bit_rep, &
                    encode_sign, extract_sign
use bit_rep_data, only: flag_deterministic, test_flag
use IntegralsData, only: UMAT
use UMatCache, only: UMatInd, GTID
use SystemData, only: NEl, nBasis, tStoreSpinOrbs, G1, BRR, lNoSymmetry, &
                      ARR, tUseMP2VarDenMat, Ecore, LMS, tHPHF, tFixLz, &
                      iMaxLz, tRef_Not_HF, tOddS_hphf, tROHF
use NatOrbsMod, only: NatOrbMat, NatOrbMatTag, Evalues, EvaluesTag, &
                      SetupNatOrbLabels
use CalcData, only: MemoryFacPart, NMCyc, InitiatorWalkNo
use OneEInts, only: TMAT2D
use FciMCData, only: MaxWalkersPart, MaxSpawned, Spawned_Parents, &
                     PreviousCycles, Spawned_Parents_Index, &
                     Spawned_ParentsTag, AccumRDMNorm_Inst, &
                     Spawned_Parents_IndexTag, Iter, AccumRDMNorm, &
                     InstNoatHF, tSinglePartPhase, AllAccumRDMNorm, iLutRef,&
                     HFDet_True, ilutHF_True, SpawnVec, &
                     SpawnVec2, SpawnVecTag, SpawnVec2Tag, SpawnedParts, &
                     SpawnedParts2, excit_gen_store_type, CurrentDets, &
                     IterRDMStart, ValidSpawnedList, &
                     TempSpawnedPartsInd, TempSpawnedParts, TotParts, &
                     TotWalkers, iLutHF, core_space, IterLastRDMFill, &
                     determ_sizes, determ_displs, &
                     full_determ_vecs_av, tFill_RDM, VaryShiftIter, &
                     IterRDM_HF, tFinalRDMEnergy
use LoggingData, only: RDMExcitLevel, tROFciDump, NoDumpTruncs, &
                   tExplicitAllRDM, tPrint1RDM, RDMEnergyIter, &
                   tDo_Not_Calc_RDMEnergy, tDiagRDM, tReadRDMs, &
                   tPopsfile, tNo_RDMs_to_read, twrite_RDMs_to_read, &
                   tWriteMultRDMs, tDumpForcesInfo, IterRDMonFly, &
                   tWrite_normalised_RDMs, IterWriteRDMs, tPrintRODump, &
                   tNoNOTransform, tTruncRODump, tRDMonfly, &
                   ThreshOccRDM, tThreshOccRDMDiag, tDipoles,&
                   tBrokenSymNOs, occ_numb_diff, tForceCauchySchwarz, &
                   tBreakSymNOs, RotNOs, tagRotNOs, local_cutoff, rottwo, &
                   rotthree, rotfour, tRDMInstEnergy,tFullHFAv, tWriteSpinFreeRDM
use RotateOrbsData, only: CoeffT1Tag, tTurnStoreSpinOff, NoFrozenVirt, &
                          SymLabelCounts2_rot,SymLabelList2_rot, &
                          SymLabelListInv_rot,NoOrbs, SpatOrbs, &
                          SymLabelCounts2_rotTag, SymLabelList2_rotTag, &
                          NoRotOrbs, SymLabelListInv_rotTag, NoOrbs
use DeterminantData, only: write_det
use hphf_integrals, only: hphf_sign
use HPHFRandExcitMod, only: FindExcitBitDetSym, FindDetSpinSym
use DetBitOps, only: TestClosedShellDet, FindBitExcitLevel, DetBitEQ, &
                     EncodeBitDet, DetBitLT, get_bit_excitmat
use DetCalcData, only: Det, FCIDets
use hist_data, only: AllHistogram, Histogram
use RotateOrbsMod, only: FourIndInts, FourIndIntsTag, PrintROFCIDUMP, &
                         PrintRepeatROFCIDUMP
use hist, only: find_hist_coeff_explicit
use sparse_arrays, only: sparse_core_ham, SparseCoreHamilTags, deallocate_sparse_ham, &
                        core_connections
use SymExcit3, only: GenExcitations3
use OneEInts, only: TMAT2D
use load_balance_calcnodes, only: DetermineDetNode
use global_det_data, only: get_iter_occ, get_av_sgn
use constants
use util_mod
use sort_mod

implicit none

integer, allocatable :: Sing_InitExcSlots(:),Sing_ExcList(:)
integer, allocatable :: Doub_InitExcSlots(:),Doub_ExcList(:)
integer(kind=n_int), allocatable :: Sing_ExcDjs(:,:),Sing_ExcDjs2(:,:)
integer(kind=n_int), allocatable :: Doub_ExcDjs(:,:),Doub_ExcDjs2(:,:)

integer :: Sing_ExcDjsTag,Sing_ExcDjs2Tag
integer :: Doub_ExcDjsTag,Doub_ExcDjs2Tag,UMATTempTag
integer :: Energies_unit, ActualStochSign_unit
integer :: NoSymLabelCounts, Rho_iiTag

real(dp), allocatable, target :: aaaa_RDM_inst(:,:), abab_RDM_inst(:,:), abba_RDM_inst(:,:)
real(dp), allocatable, target :: aaaa_RDM_full(:,:),abab_RDM_full(:,:), abba_RDM_full(:,:)
real(dp), allocatable, target :: bbbb_RDM_inst(:,:),baba_RDM_inst(:,:), baab_RDM_inst(:,:)
real(dp), allocatable, target :: bbbb_RDM_full(:,:),baba_RDM_full(:,:), baab_RDM_full(:,:)

integer :: aaaa_RDM_instTag,abab_RDM_instTag, abba_RDM_instTag
integer :: aaaa_RDM_fullTag,abab_RDM_fullTag, abba_RDM_fullTag
integer :: bbbb_RDM_instTag,baba_RDM_instTag, baab_RDM_instTag
integer :: bbbb_RDM_fullTag,baba_RDM_fullTag, baab_RDM_fullTag

real(dp), pointer :: aaaa_RDM(:,:) => null()
real(dp), pointer :: abab_RDM(:,:) => null()
real(dp), pointer :: abba_RDM(:,:) => null()
real(dp), pointer :: bbbb_RDM(:,:) => null()
real(dp), pointer :: baba_RDM(:,:) => null()
real(dp), pointer :: baab_RDM(:,:) => null()

real(dp), allocatable :: AllNodes_aaaa_RDM(:,:), AllNodes_abab_RDM(:,:), AllNodes_abba_RDM(:,:) 
real(dp), allocatable :: AllNodes_bbbb_RDM(:,:), AllNodes_baba_RDM(:,:), AllNodes_baab_RDM(:,:)  

real(dp), allocatable :: UMATTemp(:,:), Rho_ii(:)
real(dp), allocatable :: Lagrangian(:,:)

real(dp) :: OneEl_Gap,TwoEl_Gap, Normalisation,Trace_2RDM_Inst, Trace_2RDM, Trace_1RDM, norm
logical :: tCalc_RDMEnergy, tOpenShell
type(timer), save :: nElRDM_Time, FinaliseRDM_time, RDMEnergy_time
logical :: trotatedNOs = .false.

contains

    subroutine InitRDM()

        ! This routine initialises any of the arrays needed to calculate the
        ! reduced density matrix. It is used for both the explicit and
        ! stochastic RDMs.

        integer :: ierr,i, MemoryAlloc, MemoryAlloc_Root
        character(len=*), parameter :: this_routine='InitRDM'

        ! First thing is to check we're not trying to fill the RDMs in a way
        ! that is not compatible with the code (not every case has been
        ! accounted for yet).

#ifdef __CMPLX
        CAll Stop_All(this_routine,'Filling of reduced density matrices not working with &
                                    &complex walkers yet.')
#endif

        ! Only spatial orbitals for the 2-RDMs (and F12).
                
        if (tStoreSpinOrbs .and. (RDMExcitLevel .ne. 1)) &
            call stop_all(this_routine, '2-RDM calculations not set up for systems stored &
                                         &as spin orbitals.')

        if (tROHF .or. tStoreSpinOrbs) then
            tOpenShell = .true.
        else
            tOpenShell = .false.
        end if

        if (tExplicitAllRDM) then
            write(6,'(A)') " Explicitly calculating the reduced density matrices from the &
                                                        &FCIQMC wavefunction."
        else
            write(6,'(A)') " Stochastically calculating the reduced density matrices from the &
                            &FCIQMC wavefunction" 
            write(6,'(A)', advance='no') " incl. explicit connections to the following HF determinant:"
            call write_det (6, HFDet_True, .true.)
        end if

        if (RDMExcitLevel .eq. 1) then
            tCalc_RDMEnergy = .false.
        else
            ! If the RDMExcitLevel is 2 or 3 - and we're calculating the 2-RDM, 
            ! then we automatically calculate the energy unless we specifically
            ! say not to.
            if (tDo_Not_Calc_RDMEnergy) then
                tCalc_RDMEnergy = .false.            
            else
                tCalc_RDMEnergy = .true.
                write(6,'(A)') ' Calculating the energy from the reduced &
                &density matrix, this requires the 2 electron RDM from which the 1-RDM can also be constructed.'
            end if
        end if

        ! Have not got HPHF working with the explicit or truncated methods yet.
        ! Neither of these would be too difficult to implement.
        if (tHPHF .and. tExplicitAllRDM) call Stop_All('InitRDM',&
                'HPHF not set up with the explicit calculation of the RDM.')

        SpatOrbs = nBasis/2
        if (tOpenShell) then
            NoOrbs = nBasis
        else
            NoOrbs = SpatOrbs
        end if

        ! Here we're allocating arrays for the actual calculation of the RDM.
        MemoryAlloc = 0
        MemoryAlloc_Root = 0   ! Memory allocated in bytes.

        ! First for the storage of the actual 1- or 2-RMD.
        if (RDMExcitLevel .eq. 1) then

            ! This is the AllnElRDM, called NatOrbMat simply because we use the
            ! natural orbital routines to diagonalise etc. We don't have an
            ! instantaneous 1-RDM.
            allocate(NatOrbMat(NoOrbs, NoOrbs), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating NatOrbMat array,')
            call LogMemAlloc('NatOrbMat', NoOrbs**2, 8, this_routine, NatOrbMatTag, ierr)
            NatOrbMat(:,:) = 0.0_dp

            MemoryAlloc = MemoryAlloc + ( NoOrbs * NoOrbs * 8 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 ) 
        else
            ! If we're calculating the 2-RDM, the 1-RDM does not need to be
            ! calculated as well because all its info is in the 2-RDM anyway.

            ! The 2-RDM of the type alpha alpha alpha alpha (= beta beta beta beta).
            ! These *do not* include any 2-RDM(i,j,a,b) terms where i=j or a=b (if
            ! they're the same spin this can't happen).

            if (tRDMInstEnergy) then
                ! We will be filling up aaaa_RDM_inst as we go along, which need
                ! to be allocated per core.
                ! When calculating the energy, these will be summed over cores
                ! using an _inplace type command.
                ! To calculate the full energy of the RDM (i.e. over full accum.
                ! period), we need to allocate aaaa_RDM_full on the head nodes

                allocate(aaaa_RDM_inst(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating aaaa_RDM_inst array,')
                call LogMemAlloc('aaaa_RDM_inst',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,aaaa_RDM_instTag,ierr)
                aaaa_RDM_inst(:,:)=0.0_dp

                ! The 2-RDM of the type alpha beta beta alpha (= beta alpha alpha beta).
                ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b
                ! (these are the same as the abab elements).
                allocate(abba_RDM_inst(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating abba_RDM_inst array,')
                call LogMemAlloc('abba_RDM_inst', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, abba_RDM_instTag, ierr)
                abba_RDM_inst(:,:) = 0.0_dp

                MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 

                ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
                ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're
                ! different spin this is possible - hence the slightly different size to
                ! the aaaa array.
                allocate(abab_RDM_inst(((SpatOrbs*(SpatOrbs+1))/2), ((SpatOrbs*(SpatOrbs+1))/2)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating abab_RDM_inst array,')
                call LogMemAlloc('abab_RDM_inst', (((SpatOrbs*(SpatOrbs+1))/2)**2), 8, this_routine, abab_RDM_instTag, ierr)
                abab_RDM_inst(:,:) = 0.0_dp

                MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 

                ! Extra arrays for open shell systems.
                if (tOpenShell) then
                    allocate(bbbb_RDM_inst(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating bbbb_RDM_inst array,')
                    call LogMemAlloc('bbbb_RDM_inst', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, bbbb_RDM_instTag, ierr)
                    bbbb_RDM_inst(:,:) = 0.0_dp

                    allocate(baab_RDM_inst(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating baab_RDM_inst array,')
                    call LogMemAlloc('baab_RDM_inst',(((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, baab_RDM_instTag, ierr)
                    baab_RDM_inst(:,:) = 0.0_dp
                    MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 )
                    MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 )

                    allocate(baba_RDM_inst(((SpatOrbs*(SpatOrbs+1))/2), ((SpatOrbs*(SpatOrbs+1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating baba_RDM_inst array,')
                    call LogMemAlloc('baba_RDM_inst', (((SpatOrbs*(SpatOrbs+1))/2)**2), 8, this_routine, baba_RDM_instTag, ierr)
                    baba_RDM_inst(:,:) = 0.0_dp
                    MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 
                    MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 )
                end if

                if (iProcIndex .eq. 0) then
                    allocate(aaaa_RDM_full(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating aaaa_RDM_full array,')
                    call LogMemAlloc('aaaa_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, aaaa_RDM_fullTag, ierr)
                    aaaa_RDM_full(:,:)=0.0_dp

                    allocate(abab_RDM_full(((SpatOrbs*(SpatOrbs+1))/2), ((SpatOrbs*(SpatOrbs+1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating abab_RDM_full array,')
                    call LogMemAlloc('abab_RDM_full', (((SpatOrbs*(SpatOrbs+1))/2)**2), 8, this_routine, abab_RDM_fullTag, ierr)
                    abab_RDM_full(:,:) = 0.0_dp

                    allocate(abba_RDM_full(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat = ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating abba_RDM_full array,')
                    call LogMemAlloc('abba_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, abba_RDM_fullTag, ierr)
                    abba_RDM_full(:,:) = 0.0_dp

                    MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                    MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 
                    MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                    MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 

                    if (tOpenShell) then
                        allocate(bbbb_RDM_full(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                        if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating bbbb_RDM_full array,')
                        call LogMemAlloc('bbbb_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine,bbbb_RDM_fullTag,ierr)
                        bbbb_RDM_full(:,:) = 0.0_dp

                        allocate(baba_RDM_full(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
                        if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating baba_RDM_full array,')
                        call LogMemAlloc('baba_RDM_full',(((SpatOrbs*(SpatOrbs+1))/2)**2), 8,this_routine, baba_RDM_fullTag, ierr)
                        baba_RDM_full(:,:) = 0.0_dp

                        allocate(baab_RDM_full(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                        if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating baab_RDM_full array,')
                        call LogMemAlloc('baab_RDM_full',(((SpatOrbs*(SpatOrbs-1))/2)**2), 8,this_routine, baab_RDM_fullTag, ierr)
                        baab_RDM_full(:,:) = 0.0_dp

                        MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                        MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 
                        MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                        MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 

                    end if

                end if
                
                aaaa_RDM => aaaa_RDM_inst
                abba_RDM => abba_RDM_inst
                abab_RDM => abab_RDM_inst

                if (tOpenShell) then
                    bbbb_RDM => bbbb_RDM_inst
                    baab_RDM => baab_RDM_inst
                    baba_RDM => baba_RDM_inst
                end if

            else
                ! We're not calculating an instantaneous RDM energy.
                ! Put RDM contributions directly into 'full' arrays, which are
                ! now allocated every core
                allocate(aaaa_RDM_full(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating aaaa_RDM_full array,')
                call LogMemAlloc('aaaa_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine,aaaa_RDM_fullTag,ierr)
                aaaa_RDM_full(:,:) = 0.0_dp

                ! The 2-RDM of the type alpha beta beta alpha ( = beta alpha alpha beta).
                ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b
                ! (these are the same as the abab elements).
                allocate(abba_RDM_full(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating abba_RDM_full array,')
                call LogMemAlloc('abba_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, abba_RDM_fullTag, ierr)
                abba_RDM_full(:,:) = 0.0_dp

                MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 

                ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
                ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're
                ! different spin this is possible - hence the slightly different size
                ! to the aaaa array.
                allocate(abab_RDM_full(((SpatOrbs*(SpatOrbs+1))/2), ((SpatOrbs*(SpatOrbs+1))/2)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating abab_RDM_full array,')
                call LogMemAlloc('abab_RDM_full', (((SpatOrbs*(SpatOrbs+1))/2)**2), 8, this_routine, abab_RDM_fullTag, ierr)
                abab_RDM_full(:,:) = 0.0_dp

                MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 
                
                aaaa_RDM => aaaa_RDM_full
                abba_RDM => abba_RDM_full
                abab_RDM => abab_RDM_full

                if (tOpenShell) then
                    allocate(bbbb_RDM_full(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating bbbb_RDM_full array,')
                    call LogMemAlloc('bbbb_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, bbbb_RDM_fullTag, ierr)
                    bbbb_RDM_full(:,:) = 0.0_dp

                    allocate(baab_RDM_full(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating baab_RDM_full array,')
                    call LogMemAlloc('baab_RDM_full', (((SpatOrbs*(SpatOrbs-1))/2)**2), 8, this_routine, baab_RDM_fullTag, ierr)
                    baab_RDM_full(:,:) = 0.0_dp
                    MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                    MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 

                    allocate(baba_RDM_full(((SpatOrbs*(SpatOrbs+1))/2), ((SpatOrbs*(SpatOrbs+1))/2)), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating baba_RDM_full array,')
                    call LogMemAlloc('baba_RDM_full',(((SpatOrbs*(SpatOrbs+1))/2)**2), 8, this_routine, baba_RDM_fullTag, ierr)
                    baba_RDM_full(:,:) = 0.0_dp
                    MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 
                    MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 )* 8 ) 
      
                    bbbb_RDM => bbbb_RDM_full
                    baab_RDM => baab_RDM_full
                    baba_RDM => baba_RDM_full
                end if
            end if ! not instantaneous

            if (iProcindex .eq. 0) then
                if (tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
                    ! Still need to allocate 1-RDM to get nat orb occupation numbers.
                    allocate(NatOrbMat(NoOrbs, NoOrbs), stat=ierr)
                    if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating NatOrbMat array,')
                    call LogMemAlloc('NatOrbMat', NoOrbs**2,8, this_routine, NatOrbMatTag, ierr)
                    NatOrbMat(:,:) = 0.0_dp
                    MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 ) 
                end if
            end if

        end if            

        ! We then need to allocate the arrays for excitations etc when doing
        ! the explicit all calculation.
        if (tExplicitAllRDM) then

            ! We always calculate the single stuff - and if RDMExcitLevel is 1,
            ! this is all, otherwise calculate the double stuff too.

            ! This array actually contains the excitations in blocks of the
            ! processor they will be sent to. Only needed if the 1-RDM is the
            ! only thing being calculated.
            allocate(Sing_ExcDjs(0:NIfTot, nint((NEl*nBasis)*MemoryFacPart)), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Sing_ExcDjs array.')
            call LogMemAlloc('Sing_ExcDjs', nint(NEl*nBasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int, this_routine, Sing_ExcDjsTag, ierr)

            allocate(Sing_ExcDjs2(0:NIfTot, nint((NEl*nBasis)*MemoryFacPart)), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Sing_ExcDjs2 array.')
            call LogMemAlloc('Sing_ExcDjs2', nint(NEl*nBasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int, this_routine, Sing_ExcDjs2Tag, ierr)

            Sing_ExcDjs(:,:) = 0
            Sing_ExcDjs2(:,:) = 0

            MemoryAlloc = MemoryAlloc + ( (NIfTot + 1) * nint((NEl*nBasis)*MemoryFacPart) * size_n_int * 2 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 1) * nint((NEl*nBasis)*MemoryFacPart) * size_n_int * 2 ) 

            ! We need room to potentially generate N*M single excitations but
            ! these will be spread across each processor.

            OneEl_Gap=(real(NEl,dp)*real(nBasis,dp)*MemoryFacPart)/real(nProcessors,dp)

            ! This array contains the initial positions of the excitations
            ! for each processor.
            allocate(Sing_InitExcSlots(0:(nProcessors-1)), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Sing_InitExcSlots array,')
            do i = 0, nProcessors - 1
                Sing_InitExcSlots(i) = nint(OneEl_Gap*i) + 1
            end do

            ! This array contains the current position of the excitations as
            ! they're added.
            allocate(Sing_ExcList(0:(nProcessors-1)), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Sing_ExcList array,')
            Sing_ExcList(:) = Sing_InitExcSlots(:)

            if (RDMExcitLevel .ne. 1) then
                ! This array actually contains the excitations in blocks of
                ! the processor they will be sent to.
                allocate(Doub_ExcDjs(0:NIfTot,nint(((NEl*nBasis)**2)*MemoryFacPart)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating Doub_ExcDjs array.')
                call LogMemAlloc('Doub_ExcDjs', nint(((NEl*nBasis)**2)*MemoryFacPart)&
                                *(NIfTot+1), size_n_int, this_routine, Doub_ExcDjsTag, ierr)

                allocate(Doub_ExcDjs2(0:NIfTot, nint(((NEl*nBasis)**2)*MemoryFacPart)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Doub_ExcDjs2 array.')
                call LogMemAlloc('Doub_ExcDjs2',nint(((NEl*nBasis)**2)*MemoryFacPart)&
                                *(NIfTot+1), size_n_int, this_routine, Doub_ExcDjs2Tag, ierr)

                MemoryAlloc = MemoryAlloc + ( (NIfTot + 1) * nint(((NEl*nBasis)**2)*MemoryFacPart) * size_n_int * 2 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 1) * nint(((NEl*nBasis)**2)*MemoryFacPart) * size_n_int * 2 ) 

                ! We need room to potentially generate (N*M)^2 double excitations
                ! but these will be spread across each processor.        
                TwoEl_Gap = (((real(NEl,dp)*real(nBasis,dp))**2)*MemoryFacPart)/real(nProcessors,dp)

                Doub_ExcDjs(:,:) = 0
                Doub_ExcDjs2(:,:) = 0

                ! This array contains the initial positions of the excitations
                ! for each processor.
                allocate(Doub_InitExcSlots(0:(nProcessors-1)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Doub_InitExcSlots array,')
                do i = 0, nProcessors-1
                    Doub_InitExcSlots(i) = nint(TwoEl_Gap*i) + 1
                end do

                ! This array contains the current position of the excitations
                ! as they're added.
                allocate(Doub_ExcList(0:(nProcessors-1)), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Doub_ExcList array,')
                Doub_ExcList(:) = Doub_InitExcSlots(:)
            end if

        else

            ! Finally, we need to hold onto the parents of the spawned particles.
            ! This is not necessary if we're doing completely explicit calculations.
            allocate(Spawned_Parents(0:(NIfDBO+2), MaxSpawned), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating Spawned_Parents array,')
            call LogMemAlloc('Spawned_Parents', MaxSpawned*(NIfDBO+3), size_n_int,&
                                                this_routine,Spawned_ParentsTag, ierr)
            allocate(Spawned_Parents_Index(2,MaxSpawned),stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Spawned_Parents_Index array,')
            call LogMemAlloc('Spawned_Parents_Index', MaxSpawned*2,4, this_routine,&
                                                        Spawned_Parents_IndexTag, ierr)

            MemoryAlloc = MemoryAlloc + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 

            MemoryAlloc = MemoryAlloc + ( 2 * MaxSpawned * 4 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( 2 * MaxSpawned * 4 ) 

        end if

        if (iProcIndex .eq. 0) then
            write(6,"(A,F14.6,A,F14.6,A)") " Main RDM memory arrays consists of : ", &
                    & real(MemoryAlloc_Root,dp)/1048576.0_dp," Mb/Processor on the root, and ", &
                    & real(MemoryAlloc,dp)/1048576.0_dp," Mb/Processor on other processors."
        end if

        ! These parameters are set for the set up of the symmetry arrays, which
        ! are later used for the diagonalisation / rotation of the 1-RDMs.

        if (tOpenShell) then
            if (tFixLz) then
                NoSymLabelCounts = 16 * ( (2 * iMaxLz) + 1 )
            else
                NoSymLabelCounts = 16 
            end if
        else
            if (tFixLz) then
                NoSymLabelCounts = 8 * ( (2 * iMaxLz) + 1 )
            else
                NoSymLabelCounts = 8
            end if
        end if

        if ((RDMExcitLevel .eq. 1) .or. tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
            ! These arrays contain indexing systems to order the 1-RDM orbitals
            ! in terms of symmetry. This allows the diagonalisation of the RDMs
            ! to be done in symmetry blocks (a lot quicker/easier).
            ! The 2-RDM does not need to be reordered as it's never diagonalised. 

            allocate(SymLabelCounts2_rot(2,NoSymLabelCounts), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating SymLabelCounts2_rot array,')
            call LogMemAlloc('SymLabelCounts2_rot', 2*NoSymLabelCounts, 4, this_routine, SymLabelCounts2_rotTag, ierr)
            SymLabelCounts2_rot(:,:) = 0

            allocate(SymLabelList2_rot(NoOrbs), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating SymLabelList2_rot array,')
            call LogMemAlloc('SymLabelList2_rot', NoOrbs, 4, this_routine, SymLabelList2_rotTag, ierr)
            SymLabelList2_rot(:) = 0
     
            allocate(SymLabelListInv_rot(NoOrbs), stat=ierr)
            if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating SymLabelListInv_rot array,')
            call LogMemAlloc('SymLabelListInv_rot', NoOrbs, 4, this_routine, SymLabelListInv_rotTag, ierr)
            SymLabelListInv_rot(:) = 0   

            if ((iProcIndex .eq. 0) .and. tDiagRDM) then
                allocate(Evalues(NoOrbs), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Evalues array,')
                call LogMemAlloc('Evalues', NoOrbs, 8, this_routine, EvaluesTag, ierr)
                Evalues(:) = 0.0_dp

                allocate(Rho_ii(NoOrbs), stat=ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating Rho_ii array,')
                call LogMemAlloc('Rho_ii', NoOrbs, 8, this_routine, Rho_iiTag, ierr)
                Rho_ii(:) = 0.0_dp

            end if

            ! This routine actually sets up the symmetry labels for the 1-RDM.
            ! TODO : Merge this routine (and rotations later) with the NatOrbs file.
            call SetUpSymLabels_RDM() 

        end if            

        if (iProcIndex.eq.0) write(6,'(A)') " RDM memory allocation successful... "                    

        ! Open file to keep track of RDM Energies (if they're being calculated). 
        if ((iProcIndex.eq.0).and.tCalc_RDMEnergy) then
            Energies_unit = get_free_unit()
            open(Energies_unit, file='RDMEstimates', status='unknown', position='append')

            write(Energies_unit, "('#', 4X, 'Iteration', 6X, 'Energy numerator', 6X, 'Spin^2 numerator', 9X, 'Normalisation')")
        end if
        tFinalRDMEnergy = .false.

        Trace_2RDM = 0.0_dp
        Trace_2RDM_Inst = 0.0_dp
        AccumRDMNorm = 0.0_dp
        AccumRDMNorm_Inst = 0.0_dp
        ! AccumRDMNorm is the normalisation for the case where we're using a
        ! limited reference to calculate the RDM.

        ! Reads in the RDMs from a previous calculation, sets the accumulating
        ! normalisations, writes out the starting energy.
        if (tReadRDMs) then
            if (tSinglePartPhase(1) .or. tSinglePartPhase(inum_runs)) then
                write(6,'(A)') 'WARNING - Asking to read in the RDMs, but not varying shift from &
                                & the beginning of the calculation.'
                write(6,'(A)') 'Ignoring the request to read in the RDMs and starting again.'
                tReadRDMs = .false.
            else
                call Read_In_RDMs()
            end if
        end if

        ! By default, if we're writing out a popsfile (and doing an RDM
        ! calculation), we also write out the unnormalised RDMs that can be
        ! read in when restarting a calculation. If the NORDMSTOREAD option
        ! is on, these wont be printed.  
        if (TPopsFile .and. (.not. tno_RDMs_to_read)) twrite_RDMs_to_read = .true.

        nElRDM_Time%timer_name = 'nElRDMTime'
        FinaliseRDM_Time%timer_name = 'FinaliseRDMTime'
        RDMEnergy_Time%timer_name = 'RDMEnergyTime'

    end subroutine InitRDM

    subroutine zero_rdms()

        if (RDMExcitLevel.eq.1) then
            NatOrbMat(:,:) = 0.0_dp
        else
            aaaa_RDM(:,:) = 0.0_dp
            abab_RDM(:,:) = 0.0_dp
            abba_RDM(:,:) = 0.0_dp
            AllNodes_aaaa_RDM(:,:) = 0.0_dp
            AllNodes_abab_RDM(:,:) = 0.0_dp
            AllNodes_abba_RDM(:,:) = 0.0_dp

            if (tOpenShell) then
                bbbb_RDM(:,:) = 0.0_dp
                baba_RDM(:,:) = 0.0_dp
                baab_RDM(:,:) = 0.0_dp
                AllNodes_bbbb_RDM(:,:) = 0.0_dp
                AllNodes_baba_RDM(:,:) = 0.0_dp
                AllNodes_baab_RDM(:,:) = 0.0_dp
            end if

        end if

        Trace_2RDM = 0.0_dp
        Trace_2RDM_Inst = 0.0_dp
        AccumRDMNorm = 0.0_dp
        AccumRDMNorm_Inst = 0.0_dp

    end subroutine zero_rdms

    subroutine Read_In_RDMs()

        ! Reads in the arrays to restart the RDM calculation (and continue
        ! accumulating). These arrays are not normalised, so the trace is
        ! also calculated. The energy is then calculated (if required) from
        ! the RDMs read in only.

        logical :: exists_one
        logical :: exists_aaaa,exists_abab,exists_abba
        logical :: exists_bbbb,exists_baba,exists_baab
        integer :: RDM_unit, FileEnd
        integer :: i,j,a,b,Ind1,Ind2
        real(dp) :: Temp_RDM_Element, Norm_2RDM

        if (iProcIndex.eq.0) then 

            if (RDMExcitLevel.eq.1) then

                write(6,'(A)') ' Reading in the 1-RDM'

                ! The OneRDM will have been printed exactly as is. Without
                ! having been made hermitian, without being normalised, and in
                ! spatial orbitals if tStoreSpinOrbs is false.

                inquire(file='OneRDM_POPS', EXIST=exists_one)
                if (exists_one) then
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='OneRDM_POPS', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i,j,Temp_RDM_Element
                        if (FileEnd .gt. 0) call stop_all("Read_In_RDMs", "Error reading OneRDM_POPS")
                        if (FileEnd .lt. 0) exit

                        NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = Temp_RDM_Element
                    end do
                    close(RDM_unit)
                else
                    call stop_all('Read_In_RDMs','Attempting to read in the OneRDM, but &
                                    &the OneRDM_POPS file does not exist.')
                end if                                    

            else

                write(6,'(A)') ' Reading in the 2-RDMs'

                ! The TwoRDMs will have been printed exactly as they were.
                ! Without having been made hermitian, without being
                ! normalised, and in spatial orbitals. 

                ! Only read in the 2-RDMs (the 1-RDM becomes redundant).
                inquire(file='TwoRDM_POPS_aaaa', EXIST=exists_aaaa)
                inquire(file='TwoRDM_POPS_abab', EXIST=exists_abab)
                inquire(file='TwoRDM_POPS_abba', EXIST=exists_abba)

                if (tOpenShell)then
                    inquire(file='TwoRDM_POPS_bbbb', EXIST=exists_bbbb)
                    inquire(file='TwoRDM_POPS_baba', EXIST=exists_baba)
                    inquire(file='TwoRDM_POPS_baab', EXIST=exists_baab)
                end if

                if (exists_aaaa .and. exists_abab .and. exists_abba) then
                    ! All TOREAD RDM files are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='TwoRDM_POPS_aaaa', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all("Read_In_RDMs", "Error reading TwoRDM_POPS_aaaa")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        aaaa_RDM_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_abab', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all("Read_In_RDMs", "Error reading TwoRDM_POPS_abab")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        abab_RDM_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_abba', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i,j,a,b,Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all("Read_In_RDMs", "Error reading TwoRDM_POPS_abba")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        abba_RDM_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                else
                    write(6,*) 'exists_aaaa', exists_aaaa
                    write(6,*) 'exists_abab', exists_abab
                    write(6,*) 'exists_abba', exists_abba
                    call neci_flush(6)
                    call Stop_All('Read_in_RDMs',"Attempting to read in the RDMs, &
                                    &but at least one of the TwoRDM_a***_TOREAD files are missing.")
                end if

                if (tOpenShell .and. exists_bbbb .and. exists_baba .and. exists_baab) then
                    ! All TOREAD RDM files for open shell RDMs are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='TwoRDM_POPS_bbbb', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all("Read_In_RDMs", "Error reading TwoRDM_POPS_bbbb")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        bbbb_RDM_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_baba', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat = FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd.gt.0) call stop_all("Read_In_RDMs","Error reading TwoRDM_POPS_baba")
                        if (FileEnd.lt.0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        baba_RDM_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit,file='TwoRDM_POPS_baab',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat = FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all("Read_In_RDMs","Error reading TwoRDM_POPS_baab")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        baab_RDM_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                else if (tOpenShell) then
                    write(6,*) 'exists_bbbb', exists_bbbb
                    write(6,*) 'exists_baba', exists_baba
                    write(6,*) 'exists_baab', exists_baab
                    call neci_flush(6)
                    call Stop_All('Read_in_RDMs',"Attempting to read in the RDMs, &
                                    &but at least one of the TwoRDM_b***_TOREAD files are missing.")
                end if

            end if
        end if

        ! Calculate the energy for the matrices read in (if we're calculating more
        ! than the 1-RDM).
        if (tCalc_RDMEnergy) then
            call Calc_Energy_from_RDM(Norm_2RDM)
        end if

        ! Continue calculating the RDMs from the first iteration when the popsfiles
        ! (and RDMs) are read in. This overwrites the iteration number put in the input.
        IterRDMonFly = 0

    end subroutine Read_In_RDMs

    subroutine SetUpSymLabels_RDM() 

        ! This routine just sets up the symmetry labels so that the orbitals
        ! are ordered according to symmetry (all beta then all alpha if spin orbs).

        integer, allocatable :: SymOrbs_rot(:)
        integer :: LabOrbsTag, SymOrbs_rotTag, ierr, i, j, SpatSym, LzSym 
        integer :: lo, hi, Symi, SymCurr, Symi2, SymCurr2
        character(len=*), parameter :: this_routine = 'SetUpSymLabels_RDM'

        ! This is only allocated temporarily to be used to order the orbitals by.
        allocate(SymOrbs_rot(NoOrbs),stat=ierr)
        call LogMemAlloc('SymOrbs_rot',NoOrbs,4,this_routine,SymOrbs_rotTag,ierr)
        if (ierr .ne. 0) call Stop_All(this_routine,"Mem allocation for SymOrbs_rot failed.")

        ! Now we want to put the spatial orbital index, followed by the symmetry.
        SymLabelList2_rot(:) = 0
        SymOrbs_rot(:)=0

        ! *** STEP 1 *** Fill SymLabelList2_rot.
        ! Find the orbitals and order them in terms of symmetry.
        do i=1,SpatOrbs
            if (tOpenShell) then
                ! For open shell systems, all alpha are followed by all beta.
                SymLabelList2_rot(i) = BRR(2*i)
                SymLabelList2_rot(i+SpatOrbs) = BRR((2*i)-1)

                if (tFixLz) then
                    SpatSym = int(G1(BRR(2*i))%sym%S)
                    LzSym = int(G1(BRR(2*i))%Ml)
                    SymOrbs_rot(i) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )

                    SpatSym = int(G1(BRR((2*i)-1))%sym%S)
                    LzSym = int(G1(BRR((2*i)-1))%Ml)
                    SymOrbs_rot(i+SpatOrbs) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )
                else
                    SymOrbs_rot(i) = int(G1(BRR(2*i))%sym%S) 
                    SymOrbs_rot(i+SpatOrbs) = int(G1(BRR((2*i)-1))%sym%S) 
                end if
            else
                SymLabelList2_rot(i) = gtID(BRR(2*i))
                if (tFixLz) then
                    SpatSym = int(G1(BRR(2*i))%sym%S)
                    LzSym = int(G1(BRR(2*i))%Ml)
                    SymOrbs_rot(i) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )
                else
                    SymOrbs_rot(i) = int(G1(BRR(2*i))%sym%S)
                end if
                ! Orbital BRR(2*i) for i = 1 will be the beta orbital with the 
                ! second lowest energy - want the spatial orbital index to go with this.
                ! G1 also in spin orbitals - get symmetry of this beta orbital, will 
                ! be the same as the spatial orbital.
            end if
        end do

        call sort (SymOrbs_rot(1:SpatOrbs), SymLabelList2_rot(1:SpatOrbs))
        ! Sorts SymLabelList2_rot according to the order of SymOrbs_rot
        ! (i.e. in terms of symmetry).
        if (tOpenShell) &
            call sort (SymOrbs_rot(SpatOrbs+1:nBasis), SymLabelList2_rot(SpatOrbs+1:nBasis))
            ! Also do this for the beta set if spin orbitals.

        ! *** STEP 2 *** Fill SymLabelCounts2_rot_rot. This is like
        ! SymLabelCounts2_rot, but has a different ordering - BEWARE.
        ! SymLabelCounts(1,:) contains the position in SymLabelList2_rot
        ! where the symmetry index starts, SymLabelCounts(2,:) contains the
        ! number of orbitals in that symmetry index. Again if spin orbs, all
        ! alpha are followed by all beta - i.e. first 8 refer to alpha, second
        ! 8 to beta.

        if (lNoSymmetry) then
            ! If we are ignoring symmetry, all orbitals essentially have
            ! symmetry 0.
            SymLabelCounts2_rot(1,1) = 1
            SymLabelCounts2_rot(2,1) = SpatOrbs
            if (tOpenShell) then
                SymLabelCounts2_rot(1,9) = SpatOrbs+1
                SymLabelCounts2_rot(2,9) = SpatOrbs
            end if
        else 
            ! Otherwise we run through the occupied orbitals, counting the
            ! number with each symmetry (spatial and Lz) and noting where in
            ! SymLabelList2_rot each symmetry block starts.
            SymCurr = 0
            SymLabelCounts2_rot(1,1) = 1
            if (tOpenShell) then
                SymCurr2 = 0
                SymLabelCounts2_rot(1,9) = SpatOrbs + 1
            end if
            do i = 1,SpatOrbs
                if (tOpenShell) then
                    Symi = SymOrbs_rot(i)
                    Symi2 = SymOrbs_rot(i + SpatOrbs)
                else
                    Symi = SymOrbs_rot(i)
                end if
                SymLabelCounts2_rot(2,(Symi+1)) = SymLabelCounts2_rot(2,(Symi+1)) + 1
                if (Symi .ne. SymCurr) then
                    do j = SymCurr + 1, Symi
                        SymLabelCounts2_rot(1,(j+1)) = i
                    end do
                    SymCurr = Symi
                end if
                if (tOpenShell) then
                    SymLabelCounts2_rot(2,(Symi2+9)) = SymLabelCounts2_rot(2,(Symi2+9))+1
                    if (Symi2 .ne. SymCurr2) then
                        do j = SymCurr2 + 1, Symi2
                            SymLabelCounts2_rot(1,(j+9)) = i + SpatOrbs
                        end do
                        SymCurr2 = Symi2
                    end if
                end if
            end do
        end if

        ! Go through each symmetry group, making sure the orbitals are 
        ! ordered lowest to highest within each symmetry.
        do i = 1, NoSymLabelCounts
            if (SymLabelCounts2_rot(2,i) .ne. 0) then
                lo = SymLabelCounts2_rot(1, i)
                hi = lo + SymLabelCounts2_rot(2, i) - 1
                call sort (SymLabelList2_rot (lo:hi))
            end if
        end do

        ! Construct the inverse matrix. While SymLabelList2_rot takes a
        ! position and tells us what orbital is in it, we also might need to
        ! take an orbital and find out what position to put its contribution in.
        do i = 1, NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i)) = i
        end do

        ! Deallocate the arrays just used in this routine.
        deallocate(SymOrbs_rot)
        call LogMemDealloc(this_routine,SymOrbs_rotTag)

    end subroutine SetUpSymLabels_RDM

    subroutine DeAlloc_Alloc_SpawnedParts()

        ! When calculating the RDMs, we need to store the parent from which a
        ! child is spawned along with the children in the spawned array. This
        ! means a slightly larger array is communicated between processors,
        ! which there is no point in doing for the first part of the calculation.
        ! When we start calculating the RDMs this routine is called and the
        ! SpawnedParts array is made larger to accommodate the parents.

        integer :: ierr                               
        character(len=*), parameter :: this_routine='DeAlloc_Alloc_SpawnedParts'
        
        deallocate(SpawnVec)
        call LogMemDealloc(this_routine,SpawnVecTag)
        deallocate(SpawnVec2)
        call LogMemDealloc(this_routine,SpawnVec2Tag)
 
        allocate(SpawnVec(0:(NIftot+NIfDBO+2),MaxSpawned),stat=ierr)
        call LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+NIfDBO+3),size_n_int,this_routine,SpawnVecTag,ierr)
        allocate(SpawnVec2(0:(NIfTot+NIfDBO+2),MaxSpawned),stat=ierr)
        call LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+NIfDBO+3),size_n_int,this_routine,SpawnVec2Tag,ierr)

!        SpawnVec(:,:)=0
!        SpawnVec2(:,:)=0

        ! Point at correct spawning arrays
        SpawnedParts => SpawnVec
        SpawnedParts2 => SpawnVec2

        write(6,'(A54,F10.4,A4,F10.4,A13)') 'Memory requirement for spawned arrays increased from ',&
                                        real(((NIfTot+1)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp,' to ',&
                                        real(((NIfTot+NIfDBO+3)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp, ' Mb/Processor'

    end subroutine DeAlloc_Alloc_SpawnedParts

    subroutine extract_bit_rep_avsign_no_rdm(iLutnI, j, nI, SignI, FlagsI, IterRDMStartI, AvSignI, Store)

        ! This is just the standard extract_bit_rep routine for when we're not
        ! calculating the RDMs.    

        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: j
        integer, intent(out) :: nI(nel), FlagsI
        real(dp), dimension(lenof_sign), intent(out) :: SignI
        real(dp), dimension(lenof_sign), intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store
        
        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI, store)

        IterRDMStartI(:) = 0.0_dp
        AvSignI(:) = 0.0_dp

    end subroutine extract_bit_rep_avsign_no_rdm

    subroutine extract_bit_rep_avsign_norm(iLutnI, j, nI, SignI, FlagsI, IterRDMStartI, AvSignI, Store)

        ! The following extract_bit_rep_avsign routine extracts the bit
        ! representation of the current determinant, and calculate the average
        ! sign since this determinant became occupied. 

        ! In double run, we have to be particularly careful -- we need to start
        ! a new average when the determinant becomes newly occupied or
        ! unoccupied in either population (see CMO thesis). Additionally, we're
        ! also setting it up so that averages get restarted whenever we
        ! calculate the energy which saves a lot of faffing about, and storage
        ! of an extra set of RDMs, and is still unbiased. This is called for
        ! each determinant in the occupied list at the beginning of its FCIQMC
        ! cycle. It is used if we're calculating the RDMs with or without HPHF. 

        ! Input:    iLutnI (bit rep of current determinant).
        !           j - Which element in the CurrentDets array are we considering?
        ! Output:   nI, SignI, FlagsI after extract.                                              
        !           IterRDMStartI - new iteration the determinant became occupied (as a real).
        !           AvSignI - the new average walker population during this time (also real).

        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(out) :: nI(nel), FlagsI
        integer, intent(in) :: j
        real(dp), dimension(lenof_sign), intent(out) :: SignI
        real(dp), dimension(lenof_sign), intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store
        integer :: part_ind

        ! This is the iteration from which this determinant has been occupied.
        IterRDMStartI(1:lenof_sign) = get_iter_occ(j)
        
        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI)
            
        if (((Iter+PreviousCycles-IterRDMStart).gt.0) .and. &
            & (mod(((Iter-1)+PreviousCycles - IterRDMStart + 1),RDMEnergyIter) .eq. 0)) then 

            ! The previous iteration was one where we added in diagonal elements
            ! To keep things unbiased, we need to set up a new averaging block now.
            ! NB: if doing single run cutoff, note that doing things this way is now
            ! NOT the same as the technique described in CMO (and DMC's) thesis.
            ! Would expect diagonal elements to be slightly worse quality, improving
            ! as one calculates the RDM energy less frequently.  As this method is
            ! biased anyway, I'm not going to lose sleep over it.
            do part_ind = 1, lenof_sign
                AvSignI(part_ind) = SignI(part_ind)
                IterRDMStartI(part_ind) = real(Iter + PreviousCycles,dp)
            end do
        else
            ! Now let's consider other instances in which we need to start a new block:
            if (inum_runs.eq.2) then
#if defined(__DOUBLERUN) || defined(__PROG_NUMRUNS) || defined(__CMPLX)
                if ((SignI(1).eq.0).and.(IterRDMStartI(1).ne.0)) then
                    ! The population has just gone to zero on population 1.
                    ! Therefore, we need to start a new averaging block.
                    AvSignI(1) = 0
                    IterRDMStartI(1) = 0
                    AvSignI(2) = SignI(2)
                    IterRDMStartI(2) = real(Iter + PreviousCycles,dp)
                else if ((SignI(2) .eq. 0) .and. (IterRDMStartI(2) .ne. 0)) then
                    ! The population has just gone to zero on population 2.
                    ! Therefore, we need to start a new averaging block.
                    AvSignI(2) = 0
                    IterRDMStartI(2) = 0
                    AvSignI(1) = SignI(1)
                    IterRDMStartI(1) = real(Iter + PreviousCycles,dp)
                else if ((SignI(1) .ne. 0) .and. (IterRDMStartI(1) .eq. 0)) then
                    ! Population 1 has just become occupied.
                    IterRDMStartI(1) = real(Iter + PreviousCycles,dp)
                    IterRDMStartI(2) = real(Iter + PreviousCycles,dp)
                    AvSignI(1) = SignI(1)
                    AvSignI(2) = SignI(2)
                    if (SignI(2) .eq. 0) IterRDMStartI(2) = 0
                else if ((SignI(2) .ne. 0) .and. (IterRDMStartI(2) .eq. 0)) then
                    ! Population 2 has just become occupied.
                    IterRDMStartI(1) = real(Iter + PreviousCycles,dp)
                    IterRDMStartI(2) = real(Iter + PreviousCycles,dp)
                    AvSignI(1) = SignI(1)
                    AvSignI(2) = SignI(2)
                    if (SignI(1) .eq. 0) IterRDMStartI(1) = 0
                else
                    ! Nothing unusual has happened so update both populations
                    ! as normal.
                    do part_ind = 1, lenof_sign
                        ! Update the average population.
                        AvSignI(part_ind) = ( ((real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind)) * get_av_sgn(j, part_ind))&
                            + SignI(part_ind) ) / ( real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind) + 1.0_dp )
                    end do
                end if
#endif
            else
                do part_ind = 1, lenof_sign
                    ! If there is nothing stored there yet, the first iteration
                    ! the determinant became occupied is this one.
                    if (IterRDMStartI(part_ind) .eq. 0.0_dp) IterRDMStartI(part_ind) = real(Iter+PreviousCycles, dp)

                    ! Update the average population. This just comes out as the
                    ! current population (SignI) if this is the first  time the
                    ! determinant has become occupied.
                    AvSignI(part_ind) = ( ((real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind)) * get_av_sgn(j,part_ind)) &
                                    + SignI(part_ind) ) / ( real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind) + 1.0_dp )
                end do
            end if
        end if

    end subroutine extract_bit_rep_avsign_norm
    
    subroutine fill_rdm_diag_currdet_norm(iLutnI, nI, j, ExcitLevelI, tCoreSpaceDet)

        ! This routine calculates the diagonal RDM contribution, and explicit
        ! connections to the HF, from the current determinant. 

        ! j --> Which element of the main list CurrentDets are we considering?
        ! IterLastRDMFill is the number of iterations since the last time the
        ! RDM contributions were added in (often the frequency of the RDM
        ! energy calculation). 

        ! For the instantaneous RDMs we need to multiply the RDM contributions
        ! by either this, or the number of iterations the determinant has been
        ! occupied, which ever is fewer. For the full RDMs we need to multiply
        ! the RDM contributions by the number of iterations the determinant has
        ! been occupied.

        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: nI(nel), ExcitLevelI, j
        logical, intent(in), optional :: tCoreSpaceDet
        real(dp), dimension(lenof_sign) :: IterDetOcc
        integer(n_int) :: SpinCoupDet(0:nIfTot)
        integer :: nSpinCoup(nel), SignFac, HPHFExcitLevel, part_type
        real(dp) :: AvSignCurr(lenof_sign)
        integer :: IterLastRDMFill, AvSignIters, IterRDM
        
        ! This is the number of iterations this determinant has been occupied.
        IterDetOcc(1:lenof_sign) = real(Iter+PreviousCycles,dp) - get_iter_occ(j) + 1.0_dp
        AvSignIters = min(IterDetOcc(1), IterDetOcc(inum_runs))
        
        ! IterLastRDMFill is the number of iterations from the last time the
        ! energy was calculated.
        IterLastRDMFill = mod((Iter+PreviousCycles - IterRDMStart + 1),RDMEnergyIter)
        

        ! The number of iterations we want to weight this RDM contribution by is:
        if (IterLastRDMFill .gt. 0) then
            IterRDM = min(AvSignIters,IterLastRDMFill)
        else
            IterRDM = AvSignIters
        end if

        AvSignCurr = get_av_sgn(j)

        if (tHPHF) then
            if (.not.TestClosedShellDet(iLutnI)) then
                call Fill_Diag_RDM(nI, AvSignCurr/sqrt(2.0_dp), tCoreSpaceDet, IterRDM)

                
                ! C_X D_X = C_X / sqrt(2) [ D_I +/- D_I'] - for open shell dets,
                ! divide stored C_X by sqrt(2). 
                ! Add in I.
                call FindExcitBitDetSym(iLutnI, SpinCoupDet)
                call decode_bit_det (nSpinCoup, SpinCoupDet)
                ! Find out if it's + or - in the above expression.
                SignFac = hphf_sign(iLutnI)

                call Fill_Diag_RDM(nSpinCoup, real(SignFac,dp)*AvSignCurr/sqrt(2.0_dp), tCoreSpaceDet, IterRDM)

                ! For HPHF we're considering < D_I + D_I' | a_a+ a_b+ a_j a_i | D_I + D_I' >
                ! Not only do we have diagonal < D_I | a_a+ a_b+ a_j a_i | D_I > terms, but also cross terms
                ! < D_I | a_a+ a_b+ a_j a_i | D_I' > if D_I and D_I' can be connected by a single or double 
                ! excitation. Find excitation level between D_I and D_I' and add in the contribution if connected.
                HPHFExcitLevel = FindBitExcitLevel (iLutnI, SpinCoupDet, 2)
                if (HPHFExcitLevel.le.2) then 
                    call Add_RDM_From_IJ_Pair(nI, nSpinCoup, IterRDM*AvSignCurr(1)/sqrt(2.0_dp), &
                                            (real(SignFac,dp)*AvSignCurr(lenof_sign))/sqrt(2.0_dp), .true.)
                end if
            else

                ! HPHFs on, but determinant closed shell.
                call Fill_Diag_RDM(nI, AvSignCurr, tCoreSpaceDet, IterRDM)

            end if
            call Add_RDM_HFConnections_HPHF(iLutnI, nI, AvSignCurr, ExcitLevelI, IterRDM)   

        else
            ! Not using HPHFs.
            if (AvSignCurr(1)*AvSignCurr(lenof_sign) .ne. 0) call Fill_Diag_RDM(nI, AvSignCurr, tCoreSpaceDet, IterRDM)
            call Add_RDM_HFConnections_Norm(iLutnI, nI, AvSignCurr, ExcitLevelI, IterRDM)   

        end if

    end subroutine fill_rdm_diag_currdet_norm

    subroutine det_removed_fill_diag_rdm( iLutnI, j)

        ! This routine is called if a determinant is removed from the list of
        ! currently occupied. At this point we need to add in its diagonal
        ! contribution for the number of iterations it has been occupied (or
        ! since the contribution was last included). j --> which element of
        ! the main list CurrentDets are we considering.

        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: j
        integer :: nI(nel), ExcitLevel, IterLastRDMFill

        ! If the determinant is removed on an iteration that the diagonal RDM
        ! elements are  already being calculated, it will already have been
        ! counted.

        if (.not.((Iter.eq.NMCyc).or.(mod((Iter+PreviousCycles - IterRDMStart + 1),RDMEnergyIter).eq.0))) then
            ! The elements described above will have been already added in
            call decode_bit_det (nI, iLutnI)
            if (tRef_Not_HF) then
                ExcitLevel = FindBitExcitLevel (iLutHF_true, iLutnI, 2)
            else
                ExcitLevel = FindBitExcitLevel (iLutRef, iLutnI, 2)
            end if

            call fill_rdm_diag_currdet_norm(iLutnI, nI, j, ExcitLevel, .false.)

        end if

    end subroutine det_removed_fill_diag_rdm

! CMO up to here

    subroutine Add_RDM_HFConnections_Norm(iLutJ, nJ, AvSignJ, walkExcitLevel, IterRDM)

        ! This is called when we run over all TotWalkers in CurrentDets.    
        ! It is called for each CurrentDet which is a single or double of the HF.
        ! It explicitly adds in the HF - S/D connection, as if the HF were D_i and 
        ! the single or double D_j. This is the standard full space RDM calc (No HPHF).
        ! In this case the diagonal elements wll already be taken care of.

        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(NEl)
        real(dp), dimension(lenof_sign), intent(in) :: AvSignJ
        integer, intent(in) :: IterRDM
        integer, intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel, part_type

        ! Quick check that the HF population is being calculated correctly.
        if (.not.tFullHFAv) then
            ! If tFullHFAv, we continue the accumulation of AvNoAtHF even when
            ! InstNoAtHF is zero. Therefore, AvNoAtHF is allowed to be different
            ! to the AvSignJ stored in CurrentH for this det.
            if (walkExcitLevel.eq.0) then
                do part_type=1,lenof_sign
                    if (abs(AvSignJ(part_type)-AvNoatHF(part_type)) .gt. 1.0e-10_dp) then
                        write(6,*) 'HFDet_True', HFDet_True
                        write(6,*) 'nJ', nJ
                        write(6,*) 'iLutJ', iLutJ
                        write(6,*) 'AvSignJ', AvSignJ
                        write(6,*) 'AvNoatHF', AvNoatHF
                        write(6,*) "instnoathf", instnoathf
                        call Stop_All('Add_RDM_HFConnections_Norm','Incorrect average HF population.')
                    end if
                end do
            end if
        end if

        ! If we have a single or double, add in the connection to the HF,
        ! symmetrically.
        if ((walkExcitLevel .eq. 1) .or. (walkExcitLevel .eq. 2)) then
            call Add_RDM_From_IJ_Pair(HFDet_True, nJ, AvNoatHF(1), &
                                      (1.0_dp/real(lenof_sign,dp))*IterRDM*AvSignJ(lenof_sign), .true.)

            call Add_RDM_From_IJ_Pair(HFDet_True, nJ, AvNoatHF(lenof_sign), &
                                      (1.0_dp/real(lenof_sign,dp))*IterRDM*AvSignJ(1), .true.)
        end if

    end subroutine Add_RDM_HFConnections_Norm

    subroutine Add_RDM_HFConnections_HPHF(iLutJ,nJ,AvSignJ,walkExcitLevel,IterRDM)

        ! This is called when we run over all TotWalkers in CurrentDets.
        ! It is called for each CurrentDet which is a single or double of the HF.
        ! It adds in the HF - S/D connection. The diagonal elements will already
        ! have been taken care of by the extract routine.

        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer, intent(in) :: nJ(NEl)
        integer, intent(in) :: IterRDM
        real(dp), dimension(lenof_sign), intent(in) :: AvSignJ
        integer, intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel, part_type

        if (.not. tFullHFAv) then
            ! If tFullHFAv, we continue the accumulation of AvNoAtHF even
            ! when InstNoAtHF is zero. Therefore, AvNoAtHF is allowed to be
            ! different to the AvSignJ stored in CurrentH for this det.
            if (walkExcitLevel.eq.0) then
                do part_type = 1,lenof_sign
                    if (AvSignJ(part_type).ne.AvNoatHF(part_type)) then
                        write(6,*) 'AvSignJ',AvSignJ
                        write(6,*) 'AvNoatHF',AvNoatHF
                        call Stop_All('Add_RDM_HFConnections_HPHF','Incorrect average HF population.')
                    end if
                end do
            end if
        end if

        ! Now if the determinant is connected to the HF (i.e. single or double),
        ! add in the diagonal elements of this connection as well -
        ! symmetrically because no probabilities are involved.
        if ((walkExcitLevel .eq. 1) .or. (walkExcitLevel .eq. 2)) &
            call Fill_Spin_Coupled_RDM_v2(iLutHF_True, iLutJ, HFDet_True, nJ, AvNoatHF(1), IterRDM*AvSignJ(lenof_sign), .true.)

    end subroutine Add_RDM_HFConnections_HPHF

    subroutine calc_rdmbiasfac(p_spawn_rdmfac,p_gen,SignCurr,RDMBiasFacCurr)

        real(dp), intent(in) :: p_gen
        real(dp), intent(in) :: SignCurr
        real(dp), intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: p_spawn_rdmfac
        real(dp) :: p_notlist_rdmfac, p_spawn, p_not_spawn, p_max_walktospawn

        ! We eventually turn this real bias factor into an integer to be passed
        ! around with the spawned children and their parents - this only works
        ! with 64 bit at the moment.
        if (n_int .eq. 4) call Stop_All('attempt_create_normal', &
                        'the bias factor currently does not work with 32 bit integers.')

        ! Otherwise calculate the 'sign' of Di we are eventually going to add
        ! in as Di.Dj. Because we only add in Di.Dj when we successfully spawn
        ! from Di.Dj, we need to unbias (scale up) Di by the probability of this
        ! happening. We need the probability that the determinant i, with
        ! population n_i, will spawn on j. We only consider one instance of a
        ! pair Di,Dj, so just want the probability of any of the n_i walkers
        ! spawning at least once on Dj.

        ! P_successful_spawn(j | i)[n_i] =  1 - P_not_spawn(j | i)[n_i]
        ! P_not_spawn(j | i )[n_i] is the probability of none of the n_i walkers spawning on j from i.
        ! This requires either not generating j, or generating j and not succesfully spawning, n_i times.
        ! P_not_spawn(j | i )[n_i] = [(1 - P_gen(j | i)) + ( P_gen( j | i ) * (1 - P_spawn(j | i))]^n_i

        p_notlist_rdmfac = ( 1.0_dp - p_gen ) + ( p_gen * (1.0_dp - p_spawn_rdmfac) )

        ! The bias fac is now n_i / P_successful_spawn(j | i)[n_i].

        if (real(int(SignCurr),dp).ne.SignCurr) then
            ! There's a non-integer population on this determinant. We need to
            ! consider both possibilities - whether we attempted to spawn 
            ! int(SignCurr) times or int(SignCurr)+1 times.
            p_max_walktospawn = abs(SignCurr-real(int(SignCurr),dp))
            p_not_spawn = (1.0_dp - p_max_walktospawn)*(p_notlist_rdmfac**abs(int(SignCurr))) + &
                        p_max_walktospawn*(p_notlist_rdmfac**(abs(int(SignCurr))+1))

        else
            p_not_spawn = p_notlist_rdmfac**(abs(SignCurr))
        end if

        p_spawn = abs(1.0_dp - p_not_spawn)
        
        ! Always use instantaneous signs for stochastically sampled off-diag
        ! elements (see CMO thesis).
        RDMBiasFacCurr = SignCurr / p_spawn

    end subroutine calc_rdmbiasfac

    subroutine store_parent_with_spawned(RDMBiasFacCurr, WalkerNumber, iLutI, DetSpawningAttempts, iLutJ, procJ, part_type)

        ! We are spawning from iLutI to SpawnedParts(:,ValidSpawnedList(proc)).
        ! This routine stores the parent (D_i) with the spawned child (D_j) so
        ! that we can add in Ci.Cj to the RDM later on. The parent is NIfDBO
        ! integers long, and stored in the second part of the SpawnedParts array 
        ! from NIfTot+1 -> NIfTot+1 + NIfDBO.

        real(dp), intent(in) :: RDMBiasFacCurr
        integer, intent(in) :: WalkerNumber, procJ
        integer, intent(in) :: DetSpawningAttempts
        integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
        integer, intent(in) :: part_type
        logical :: tRDMStoreParent
        integer :: j

        if (RDMBiasFacCurr.eq.0.0_dp) then
            ! If RDMBiasFacCurr is exactly zero, any contribution from Ci.Cj will be zero 
            ! so it is not worth carrying on. 
            SpawnedParts(niftot+1:niftot+nifdbo+2, ValidSpawnedList(procJ)) = 0
        else

            ! First we want to check if this Di.Dj pair has already been accounted for.
            ! This means searching the Dj's that have already been spawned from this Di, to make sure 
            ! the new Di being spawned on here is not the same.
            ! The Dj children spawned by the current Di are being stored in the array TempSpawnedParts, 
            ! so that the reaccurance of a Di.Dj pair may be monitored.

            ! Store the Di parent with the spawned child, unless we find this Dj has already been spawned on.
            tRDMStoreParent = .true.

            ! Run through the Dj walkers that have already been spawned from this particular Di.
            ! If this is the first to be spawned from Di, TempSpawnedPartsInd will be zero, so we 
            ! just wont run over anything.

            do j = 1,TempSpawnedPartsInd
                if (DetBitEQ(iLutJ(0:NIfDBO),TempSpawnedParts(0:NIfDBO,j),NIfDBO)) then
                    ! If this Dj is found, we do not want to store the parent with this spawned walker.
                    tRDMStoreParent = .false.
                    exit
                end if
            end do

            if (tRDMStoreParent) then
                ! This is a new Dj that has been spawned from this Di.
                ! We want to store it in the temporary list of spawned parts which have come from this Di.
                if (WalkerNumber.ne.DetSpawningAttempts) then
                    ! Don't bother storing these if we're on the last walker, or if we only have one 
                    ! walker on Di.
                    TempSpawnedPartsInd = TempSpawnedPartsInd + 1
                    TempSpawnedParts(0:NIfDBO,TempSpawnedPartsInd) = iLutJ(0:NIfDBO)
                end if

                ! We also want to make sure the parent Di is stored with this Dj.
                SpawnedParts(niftot+1:niftot+nifdbo+1, ValidSpawnedList(procJ)) = iLutI(0:nifdbo) 

                ! We need to carry with the child (and the parent), the sign of the parent.
                ! In actual fact this is the sign of the parent divided by the probability of generating 
                ! that pair Di and Dj, to account for the 
                ! fact that Di and Dj are not always added to the RDM, but only when Di spawns on Dj.
                ! This RDMBiasFacCurr factor is turned into an integer to pass around to the relevant processors.
                SpawnedParts(niftot+nifdbo+2, ValidSpawnedList(procJ)) = &
                    transfer(RDMBiasFacCurr,SpawnedParts(niftot+nifdbo+2, ValidSpawnedList(procJ)))

            else
                ! This Di has already spawned on this Dj - don't store the Di parent with this child, 
                ! so that the pair is not double counted.  
                ! We are using the probability that Di spawns onto Dj *at least once*, so we don't want to 
                ! double count this pair.
                SpawnedParts(niftot+1:niftot+nifdbo+2, ValidSpawnedList(procJ)) = 0
            end if
        end if

    end subroutine store_parent_with_spawned

    subroutine check_fillRDM_DiDj(Spawned_No, iLutJ, realSignJ)

        ! The spawned parts contain the Dj's spawned by the Di's in CurrentDets.
        ! If the SpawnedPart is found in the CurrentDets list, it means that
        ! the Dj has a non-zero cj - and therefore the Di.Dj pair will have a
        ! non-zero ci.cj to contribute to the RDM. The index i tells us where
        ! to look in the parent array, for the Di's to go with this Dj.

        integer, intent(in) :: Spawned_No
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp), dimension(lenof_sign), intent(in) :: realSignJ
        integer :: ExcitLevel
    
        if (.not. DetBitEQ(iLutHF_True, iLutJ, NIfDBO)) then
                call DiDj_Found_FillRDM(Spawned_No, iLutJ, realSignJ)
        end if

    end subroutine check_fillRDM_DiDj
 
    subroutine DiDj_Found_FillRDM(Spawned_No, iLutJ, realSignJ)

        ! This routine is called when we have found a Di (or multiple Di's)
        ! spawning onto a Dj with sign /= 0 (i.e. occupied). We then want to
        ! run through all the Di, Dj pairs and add their coefficients 
        ! (with appropriate de-biasing factors) into the 1 and 2 electron RDM.

        use bit_rep_data, only: flag_deterministic, test_flag
        integer, intent(in) :: Spawned_No
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        real(dp),dimension(lenof_sign), intent(in) :: realSignJ
        integer :: i, j, nI(NEl), nJ(NEl), walkExcitLevel
        real(dp) :: part_realSignI
        integer :: dest_part_type, source_part_type
        logical :: tParity, tDetAdded

        ! Spawning from multiple parents, to iLutJ, which has SignJ.        

        ! We are at position Spawned_No in the SpawnedParts array.
        ! Spawned_Parents_Index(1,Spawned_No) is therefore the start position
        ! of the list of parents (Di's) which spawned on the Dj in
        ! SpawnedParts(Spawned_No). There are Spawned_Parents_Index(2,Spawned_No)
        ! of these parent Di's. Spawned_Parents(0:NIfDBO,x) is the determinant Di,
        ! Spawned_Parents(NIfDBO+1,x) is the un-biased ci.

        ! Run through all Di's.

        do i = Spawned_Parents_Index(1,Spawned_No), &
                Spawned_Parents_Index(1,Spawned_No) + Spawned_Parents_Index(2,Spawned_No) - 1 

            if (DetBitEQ(iLutHF_True,Spawned_Parents(0:NIfDBO,i),NIfDBO)) then
                ! We've already added HF - S, and HF - D symmetrically.
                ! Any connection with the HF has therefore already been added.
                cycle
            end if
            
            call decode_bit_det (nI, Spawned_Parents(0:NIfDBO,i))
            call decode_bit_det (nJ, iLutJ)


            part_realSignI = transfer( Spawned_Parents(NIfDBO+1,i), part_realSignI )

            ! The original spawning event (and the RealSignI) came from this population.
            source_part_type=Spawned_Parents(NIfDBO+2,i)

            !The sign contribution from J must come from the other population.
            if (source_part_type.eq.1) then
                dest_part_type=lenof_sign
            else
                dest_part_type=1
            end if

            ! Given the Di,Dj and Ci,Cj - find the orbitals involved in the
            ! excitation, and therefore the RDM elements we want to add the
            ! Ci.Cj to. We have to halve the contributions for DR as we're
            ! summing in pairs that originated from spawning events in both
            ! pop 1 and pop 2 -- i.e. doublecounted wrt diagonal elements
            if (tHPHF) then
                call Fill_Spin_Coupled_RDM_v2(Spawned_Parents(0:NIfDBO,i), iLutJ, nI, nJ, &
                           (1.0_dp/real(lenof_sign,dp))*part_realSignI, realSignJ(dest_part_type), .false.)
            else
                call Add_RDM_From_IJ_Pair(nI, nJ, (1.0_dp/real(lenof_sign,dp))*part_realSignI, realSignJ(dest_part_type), .false.)
            end if

        end do

    end subroutine DiDj_Found_FillRDM

    subroutine Fill_Spin_Coupled_RDM_v2(iLutnI,iLutnJ,nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)

        ! This routine does the same as Fill_Spin_Coupled_RDM, but hopefully
        ! more efficiently! It takes to HPHF functions, and calculate what
        ! needs to be summed into the RDMs

        integer(n_int), intent(in) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot)
        real(dp), intent(in) :: realSignI, realSignJ
        integer, intent(in) :: nI(NEl),nJ(NEl)
        logical, intent(in) :: tFill_CiCj_Symm
        integer(n_int) :: iLutnI2(0:NIfTot)
        integer :: nI2(NEl),nJ2(NEl)
        real(dp) :: NewSignJ,NewSignI,PermSignJ,PermSignI
        integer :: I_J_ExcLevel,ICoup_J_ExcLevel
        character(*), parameter :: t_r = 'Fill_Spin_Coupled_RDM_v2'

        if (TestClosedShellDet(iLutnI)) then
            if (tOddS_HPHF) then
                call stop_all(t_r,"Should not be any closed shell determinants in high S states")
            end if

            if (TestClosedShellDet(iLutnJ)) then
                ! Closed shell -> Closed shell - just as in determinant case
                call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
            else
                ! Closed shell -> open shell.
                call FindDetSpinSym(nJ,nJ2,NEl)
                NewSignJ = realSignJ/Root2
                call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,NewSignJ,tFill_CiCj_Symm)
                ! What is the permutation between Di and Dj'
                NewSignJ = NewSignJ * hphf_sign(iLutnJ)
                call Add_RDM_From_IJ_Pair(nI,nJ2,realSignI,NewSignJ,tFill_CiCj_Symm)
            end if

        else if (TestClosedShellDet(iLutnJ)) then
            ! Open shell -> closed shell
            call FindDetSpinSym(nI,nI2,NEl)
            NewSignI = realSignI/Root2
            call Add_RDM_From_IJ_Pair(nI,nJ,NewSignI,realSignJ,tFill_CiCj_Symm)
            ! What is the permutation between Di' and Dj?
            NewSignI = NewSignI * hphf_sign(iLutnI)
            call Add_RDM_From_IJ_Pair(nI2,nJ,NewSignI,realSignJ,tFill_CiCj_Symm)

        else
            ! Open shell -> open shell
            NewSignI = realSignI/Root2
            NewSignJ = realSignJ/Root2
            PermSignJ = NewSignJ * real(hphf_sign(iLutnJ),dp)
            PermSignI = NewSignI * real(hphf_sign(iLutnI),dp)
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            call FindDetSpinSym(nI,nI2,NEl)
            call FindDetSpinSym(nJ,nJ2,NEl)
            I_J_ExcLevel = FindBitExcitLevel(iLutnI, iLutnJ,2)
            ICoup_J_ExcLevel = FindBitExcitLevel(iLutnI2,iLutnJ,2)

            if (I_J_ExcLevel .le. 2) then
                call Add_RDM_From_IJ_Pair(nI,nJ,NewSignI,NewSignJ,tFill_CiCj_Symm) 
                ! Di -> Dj
                call Add_RDM_From_IJ_Pair(nI2,nJ2,PermSignI,PermSignJ,tFill_CiCj_Symm)   
                ! Di' -> Dj'  (both permuted sign)
            end if

            if (ICoup_J_ExcLevel .le. 2) then
                call Add_RDM_From_IJ_Pair(nI2,nJ,PermSignI,NewSignJ,tFill_CiCj_Symm)    
                ! Di' -> Dj  (i permuted sign)
                call Add_RDM_From_IJ_Pair(nI,nJ2,NewSignI,PermSignJ,tFill_CiCj_Symm)     
                ! Di  -> Dj'  (j permuted sign)
            end if
        end if

    end subroutine Fill_Spin_Coupled_RDM_v2

    subroutine Fill_Spin_Coupled_RDM(iLutnI,iLutnJ,nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)

        ! Above Fill_Spin_Coupled_RDM_v2 is more efficient version of this routine.
        ! If the two HPHF determinants we're considering consist of I + I' and J + J', 
        ! where X' is the spin coupled (all spins flipped) version of X,
        ! then we have already considered the I -> J excitation.
        ! And if I and J are connected by a double excitation, tDoubleConnection is
        ! true and we have also considered I' -> J'.
        ! But we need to also account for I -> J' and I' -> J.

        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot)
        integer, intent(in) :: nI(NEl), nJ(NEl)
        real(dp), intent(in) :: realSignI, realSignJ
        logical, intent(in) :: tFill_CiCj_Symm
        integer(kind=n_int) :: iLutnI2(0:NIfTot),iLutnJ2(0:NIfTot)
        integer :: Ex(2,2), SpinCoupI_J_ExcLevel, nI2(NEl), nJ2(NEl)
        integer :: SignFacI, SignFacJ, I_J_ExcLevel
        logical :: tParity
        real(dp) :: realSignFacI, realSignFacJ

        ! First we flip the spin of both determinants, and store I' and J'.
        ! Actually if I and J are related by a double excitation, we don't need J'.        

        ! First we flip the spin of I', and find out the excitation level between I' and J.
        ! If this is a double excitation, we don't actually need to find J' - we can just invert 
        ! the excitation matrix of the I' -> J transition.
        ! If this is anything above a double, we likewise don't need to find J', because I -> J' 
        ! will also have a 0 matrix element.

        I_J_ExcLevel = FindBitExcitLevel (iLutnI, iLutnJ, 2)

        if (.not. TestClosedShellDet(iLutnI)) then

            ! I is open shell, and so a spin coupled determinant I' exists.

            ! Find I'.
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            call decode_bit_det (nI2, iLutnI2)
            SignFacI = hphf_sign(iLutnI)
            realSignFacI = real(SignFacI,dp) / sqrt(2.0_dp)

            ! Find excitation level between I' and J - not necessarily the same as 
            ! that between I and J.
            SpinCoupI_J_ExcLevel = FindBitExcitLevel (iLutnI2, iLutnJ, 2)

            if ( (.not.(I_J_ExcLevel.le.2)) .and. (.not.(SpinCoupI_J_ExcLevel.le.2)) ) &
                call Stop_All('Fill_Spin_Coupled_RDM','No spin combination are connected.')
                
            if ( .not. TestClosedShellDet(iLutnJ) ) then
                
                ! Both I and J are open shell, need all 4 combinations.

                ! Find J'.
                call FindExcitBitDetSym(iLutnJ, iLutnJ2)
                call decode_bit_det (nJ2, iLutnJ2)
                SignFacJ = hphf_sign(iLutnJ)
                realSignFacJ = real(SignFacJ,dp) / sqrt(2.0_dp)

                if (I_J_ExcLevel.le.2) then

                    ! I -> J.
                    call Add_RDM_From_IJ_Pair(nI,nJ,(realSignI/sqrt(2.0_dp)),&
                                                (realSignJ/sqrt(2.0_dp)),tFill_CiCj_Symm)
 
                    ! I' -> J'.
                    call Add_RDM_From_IJ_Pair(nI2,nJ2,(realSignFacI*realSignI),&
                                              (realSignFacJ*realSignJ),tFill_CiCj_Symm)
                end if

                if (SpinCoupI_J_ExcLevel.le.2) then

                    ! I' -> J.
                    call Add_RDM_From_IJ_Pair(nI2,nJ,(realSignFacI*realSignI),&
                                                (realSignJ/sqrt(2.0_dp)),tFill_CiCj_Symm)

                    ! I -> J'.
                    call Add_RDM_From_IJ_Pair(nI, nJ2,(realSignI/sqrt(2.0_dp)),&
                                                 (realSignFacJ*realSignJ),tFill_CiCj_Symm)

                end if

            else
                
                ! I is open shell, but J is not.
                ! Need I -> J and I' -> J.

                ! I -> J.

                call Add_RDM_From_IJ_Pair(nI,nJ,(realSignI/sqrt(2.0_dp)),&
                                                        realSignJ,tFill_CiCj_Symm)

                ! I' -> J.
                call Add_RDM_From_IJ_Pair(nI2,nJ,(realSignFacI*realSignI),&
                                                        realSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I' -> J : ",nI2(:),nJ(:),realSignFacI*realSignI,realSignJ

            end if

        else if ( .not. TestClosedShellDet(iLutnJ) ) then
            ! This is the case where I is closed shell, but J is not.
            ! Need I -> J and I -> J'. 

            ! I -> J.
            if (I_J_ExcLevel.le.2) call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,&
                                               (realSignJ/sqrt(2.0_dp)),tFill_CiCj_Symm)

            ! Find J'.
            call FindExcitBitDetSym(iLutnJ, iLutnJ2)
            SignFacJ = hphf_sign(iLutnJ)
            realSignFacJ = real(SignFacJ,dp) / sqrt(2.0_dp)

            ! Find excitation level between I and J'.
            SpinCoupI_J_ExcLevel = FindBitExcitLevel (iLutnI, iLutnJ2, 2)

            if (SpinCoupI_J_ExcLevel.le.2) then
                call decode_bit_det (nJ2, iLutnJ2)
                
                ! I -> J'.
                call Add_RDM_From_IJ_Pair(nI,nJ2,realSignI,&
                                            (realSignFacJ*realSignJ),tFill_CiCj_Symm)

           end if

       else if (I_J_ExcLevel .le. 2) then

            ! I and J are both closed shell.

            ! Just I -> J.
            call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)

       end if

    end subroutine Fill_Spin_Coupled_RDM

    subroutine Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)

        ! This routine takes a pair of different determinants Di and Dj, and
        ! figures out which type of elements need to be added in to the RDM.

        integer, intent(in) :: nI(NEl), nJ(NEl)
        real(dp), intent(in) :: realSignI, realSignJ
        logical, intent(in) :: tFill_CiCj_Symm
        integer :: Ex(2,2),j
        logical :: tParity

        Ex(:,:) = 0
        Ex(1,1) = 2         ! Maximum excitation level - we know they are connected by
                            ! a double or single.
        tParity = .false.

        call GetExcitation(nI,nJ,NEl,Ex,tParity)
        ! Ex(1,:) comes out as the orbital(s) excited from, i.e. i,j
        ! Ex(2,:) comes out as the orbital(s) excited to, i.e. a,b.

        if (Ex(1,1).le.0) then
            ! Error.
            write(6,*) '*'
            write(6,*) 'nI',nI
            write(6,*) 'nJ',nJ
            write(6,*) 'Ex(:,:)',Ex(1,1),Ex(1,2),Ex(2,1),Ex(2,2)
            write(6,*) 'tParity',tParity
            write(6,*) 'realSignI',realSignI
            write(6,*) 'realSignJ',realSignJ
            write(6,*) '*'
            call neci_flush(6)
            call Stop_All('Add_RDM_From_IJ_Pair',&
                    'Excitation level between pair not 1 or 2 as it should be.')
        end if

        if ((Ex(1,2) .eq. 0) .and. (Ex(2,2) .eq. 0)) then
            
            ! Di and Dj are separated by a single excitation.
            ! Add in the contribution from this pair into the 1- and 2-RDM.
            
            call Fill_Sings_RDM(nI,Ex,tParity,realSignI,realSignJ,tFill_CiCj_Symm)
    
        else if (RDMExcitLevel .ne. 1) then

            ! Otherwise Di and Dj are connected by a double excitation.
            ! Add in this contribution to the 2-RDM (as long as we're
            ! calculating this obv).
            call Fill_Doubs_RDM(Ex,tParity,realSignI,realSignJ,tFill_CiCj_Symm)

        end if

    end subroutine Add_RDM_From_IJ_Pair

! =======================================================================================    
! EXPLICIT ROUTINES    
! =======================================================================================    

    subroutine Fill_ExplicitRDM_this_Iter(TotWalkers)
        integer(int64), intent(in) :: TotWalkers
        integer(kind=n_int) :: iLutnI(0:NIfTot)
        integer(int64) :: MaxTotWalkers,TotWalkIn(2),TotWalkOut(2)
        integer :: i,error
        real(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs
        logical :: blank_det
        integer, DIMENSION(lenof_sign) :: SignI, SignI2

        ! Run through the current determinants.
        ! Find the max number of determinants on a processor - all need to
        ! run through this number so that the communication can be done at
        ! all stages.

        TotWalkIn(1)=TotWalkers
        TotWalkIn(2)=iProcIndex

        call MPIAllReduceDatatype(TotWalkIn,1,MPI_MAXLOC,MPI_2integer,TotWalkOut)

        MaxTotWalkers=TotWalkOut(1)

        call set_timer(nElRDM_Time,30)

        do i = 1,int(MaxTotWalkers,sizeof_int)

            ! But if the actual number of determinants on this processor is
            ! less than the number  we're running through, feed in 0
            ! determinants and 0 sign.
            if (i .gt. TotWalkers) then
                iLutnI(:) = 0
                blank_det = .true.
            else
                iLutnI(:) = CurrentDets(:,i)
                blank_det = .false.
            end if

            call Add_ExplicitRDM_Contrib(iLutnI,blank_det)

        end do

        call halt_timer(nElRDM_Time)

    end subroutine Fill_ExplicitRDM_this_Iter

    subroutine Fill_Hist_ExplicitRDM_this_Iter(TotWalkers)

        integer(int64), intent(in) :: TotWalkers
        integer(kind=n_int) :: iLutnI(0:NIfTot)
        integer :: i,error
        real(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs,AllNode_norm
        logical :: blank_det
        real(dp), DIMENSION(lenof_sign) :: TempSign

        call set_timer(nElRDM_Time,30)

        call MPISumAll(Histogram, AllHistogram)

        norm = 0.0_dp
        if (iProcIndex .eq. 0) then
            do i = 1, Det
                norm = norm + AllHistogram(1,i)**2
            end do
            norm = sqrt(norm)
        end if

        call MPISumAll(norm, allNode_norm)
        norm = allNode_norm
 
        do i = 1, Det

            ! But if the actual number of determinants on this processor is
            ! less than the number we're running through, feed in 0
            ! determinants and 0 sign.
            if (Histogram(1,i) .eq. 0.0_dp) then
                iLutnI(:) = 0
                blank_det = .true.
            else
                iLutnI(:) = FCIDets(:,i)
                blank_det = .false.
            end if

            TempSign(1) = real(i,dp)
            call encode_sign(iLutnI,TempSign)

            call Add_Hist_ExplicitRDM_Contrib(iLutnI,blank_det)

        end do

        call halt_timer(nElRDM_Time)

    end subroutine Fill_Hist_ExplicitRDM_this_Iter

    subroutine Add_ExplicitRDM_Contrib(iLutnI,blank_det)

        ! This is the general routine for taking a particular determinant in
        ! the spawned list, D_i and adding it's contribution to the reduced
        ! density matrix.

        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        logical, intent(in) :: blank_det
        integer :: i
        
        ! Set up excitation arrays.
        ! These are blocked according to the processor the excitation would be
        ! on if occupied. In each block, the first entry is the sign of
        ! determinant D_i and the second the bit string of the determinant
        ! (these need to be sent along with the excitations). Each processor
        ! will have a different Di.

        Sing_ExcDjs(:,:)=0
        Sing_ExcList(:)=0
        Sing_ExcList(:) = Sing_InitExcSlots(:)

        do i=0,nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
            Sing_ExcList(i) = Sing_ExcList(i)+1
        end do
        if (RDMExcitLevel.ne.1) then
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            do i=0,nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                Doub_ExcList(i) = Doub_ExcList(i)+1
            end do
        end if

        if (.not.blank_det) call GenExcDjs(iLutnI)
        ! Out of here we will get a filled ExcDjs array with all the single or
        ! double excitations from Dj, this will be done for each proc. 

        ! We then need to send the excitations to the relevant processors.
        call SendProcExcDjs()
        ! This routine then calls SearchOccDets which takes each excitation
        ! and and binary searches the occupied determinants for this. If found,
        ! we re-find the orbitals and parity involved in the excitation, and
        ! add the c_i*c_j contributions to the corresponding matrix element.

    end subroutine Add_ExplicitRDM_Contrib

    subroutine Add_Hist_ExplicitRDM_Contrib(iLutnI,blank_det)

        ! This is the general routine for taking a particular determinant in
        ! the spawned list, D_i and adding it's contribution to the reduced
        ! density matrix.

        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        logical, intent(in) :: blank_det
        integer :: i
        
        ! Set up excitation arrays.
        ! These are blocked according to the processor the excitation would be
        ! on if occupied. In each block, the first entry is the sign of
        ! determinant D_i and the second the bit string of the determinant
        ! (these need to be sent along with the excitations). Each processor
        ! will have a different Di.

        Sing_ExcDjs(:,:)=0
        Sing_ExcList(:)=0
        Sing_ExcList(:) = Sing_InitExcSlots(:)

        do i=0,nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
            Sing_ExcList(i) = Sing_ExcList(i)+1
        end do
        if (RDMExcitLevel.ne.1) then
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            do i=0,nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                Doub_ExcList(i) = Doub_ExcList(i)+1
            end do
        end if

        if (.not.blank_det) call Gen_Hist_ExcDjs(iLutnI)
        ! Out of here we will get a filled ExcDjs array with all the single or
        ! double excitations  from Dj, this will be done for each proc. 

        ! We then need to send the excitations to the relevant processors.
        call Send_Hist_ProcExcDjs()
        ! This routine then calls SearchOccDets which takes each excitation
        ! and and binary searches the occupied determinants for this. If found,
        ! we re-find the orbitals and parity involved in the excitation, and add
        ! the c_i*c_j contributions to the corresponding matrix element.

    end subroutine Add_Hist_ExplicitRDM_Contrib

    subroutine GenExcDjs(iLutnI)

        ! This uses GenExcitations3 in symexcit3.F90 to generate all the
        ! possible either single or double excitations from D_i, finds the
        ! processor they would be on if occupied, and puts them in the
        ! SingExcDjs array according to that processor.

        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        integer(kind=n_int) :: iLutnJ(0:NIfTot)
        real(dp), dimension(lenof_sign) :: SignDi, SignDi2
        integer :: ExcitMat3(2,2), nI(NEl), nJ(NEl), Proc, FlagsDi
        integer :: a, b, CountTemp, exflag
        logical :: tAllExcitFound, tParity

        call extract_bit_rep (iLutnI, nI, SignDi, FlagsDi)
        ! Unfortunately uses the decoded determinant - might want to look at this.        

        call Fill_Diag_RDM(nI,SignDi,.false.)

!        CountTemp = 0

        ExcitMat3(:,:)=0
        ! Zeros in ExcitMat3 starts off at the first single excitation.
        tAllExcitFound=.false.
        ! This becomes true when all the excitations have been found.

        do while (.not.tAllExcitFound)
            exflag = 1
            call GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
                                                        tAllExcitFound,.true.)            
            ! Passed out of here is the singly excited determinant, nJ.
            ! Information such as the orbitals involved in the excitation and
            ! the parity is also found in this step, we are not currently
            ! storing this, and it is re-calculated later on (after the
            ! determinants are passed to the relevant processor) - but the
            ! speed of sending this information vs recalculating it will be
            ! tested. RDMExcitLevel is passed through, if this is 1, only
            ! singles are generated, if it is 2 only doubles are found.

            if (tAllExcitFound) exit

            iLutnJ(:) = 0
            call EncodeBitDet(nJ,iLutnJ)

            Proc = DetermineDetNode(nel,nJ,0)   
            ! This will return a value between 0 -> nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
            Sing_ExcList(Proc) = Sing_ExcList(Proc)+1
!            CountTemp = CountTemp + 1

            ! Want a quick test to see if arrays are getting full.
            if (Sing_ExcList(Proc).gt.nint(OneEl_Gap*(Proc+1))) then
                write(6,*) 'Proc',Proc
                write(6,*) 'Sing_ExcList',Sing_ExcList
                write(6,*) 'No. spaces for each proc',nint(OneEl_Gap)
                call Stop_All('GenExcDjs',&
                            'Too many excitations for space available.')
            end if
        end do

        if (RDMExcitLevel .ne. 1) then            

            ExcitMat3(:,:) = 0
            ! Zeros in ExcitMat3 starts off at the first single excitation.
            tAllExcitFound = .false.
            ! This becomes true when all the excitations have been found.

            do while (.not. tAllExcitFound)
                exflag = 2
                call GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
                                                            tAllExcitFound,.true.)            

                ! Passed out of here is the doubly excited determinant, nJ.
                ! Information such as the orbitals involved in the excitation
                ! and the parity is  also found in this step, we are not
                ! currently storing this, and it is re-calculated later on
                ! (after the determinants are passed to the relevant processor)
                ! - but the speed of sending this information vs recalculating
                ! it will be tested. RDMExcitLevel is passed through, if this
                ! is 1, only singles are generated, if it is 2 only doubles are
                ! found.

                if (tAllExcitFound) exit

                iLutnJ(:) = 0
                call EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nel,nJ,0)   
                !This will return a value between 0 -> nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being 
                ! stored in Doub_ExcDjs.

                Doub_ExcList(Proc) = Doub_ExcList(Proc)+1

                ! Want a quick test to see if arrays are getting full.
                if (Doub_ExcList(Proc).gt.nint(TwoEl_Gap*(Proc+1))) then
                    write(6,*) 'Proc',Proc
                    write(6,*) 'Doub_ExcList',Doub_ExcList
                    write(6,*) 'No. spaces for each proc',nint(TwoEl_Gap)
                    call Stop_All('GenExcDjs','Too many excitations for space available.')
                end if
            end do
        end if

    end subroutine GenExcDjs

    subroutine Gen_Hist_ExcDjs(iLutnI)

        ! This uses GenExcitations3 in symexcit3.F90 to generate all the
        ! possible either single or double excitations from D_i, finds the
        ! processor they would be on if occupied, and puts them in the
        ! SingExcDjs array according to that processor.

        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        integer(kind=n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: HistPos
        real(dp), dimension(lenof_sign) :: RealHistPos
        integer :: ExcitMat3(2,2), nI(NEl), nJ(NEl), Proc, FlagsDi
        integer :: a, b, CountTemp, exflag
        logical :: tAllExcitFound, tParity
        real(dp), dimension(lenof_sign) :: realSignDi

        call extract_bit_rep (iLutnI, nI, RealHistPos, FlagsDi)
        ! Unfortunately uses the decoded determinant - might want to look at this.
        HistPos=int(RealHistPos)
        
        realSignDi(1) = AllHistogram(1,HistPos(1))/norm
        realSignDi(lenof_sign) = AllHistogram(1,HistPos(1))/norm
        
        call Fill_Diag_RDM(nI,realSignDi,.false.)

!        CountTemp = 0

        ExcitMat3(:,:)= 0
        ! Zeros in ExcitMat3 starts off at the first single excitation.
        tAllExcitFound= .false.
        ! This becomes true when all the excitations have been found.        

        do while (.not. tAllExcitFound)
            exflag = 1
            call GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
                                                        tAllExcitFound,.true.)            

            ! Passed out of here is the singly excited determinant, nJ.
            ! Information such as the orbitals involved in the excitation and
            ! the parity is also found in this step, we are not currently
            ! storing this, and it is re-calculated later on (after the
            ! determinants are passed to the relevant processor) - but the speed
            ! of sending this information vs recalculating it will be tested.
            ! RDMExcitLevel is passed through, if this is 1, only singles are
            ! generated, if it is 2 only doubles are found.

            if (tAllExcitFound) exit

            iLutnJ(:) = 0
            call EncodeBitDet(nJ,iLutnJ)

            Proc = DetermineDetNode(nel,nJ,0)   
            ! This will return a value between 0 -> nProcessors-1.
            Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
            Sing_ExcList(Proc) = Sing_ExcList(Proc)+1
!            CountTemp = CountTemp + 1

            ! Want a quick test to see if arrays are getting full.
            if (Sing_ExcList(Proc) .gt. nint(OneEl_Gap*(Proc+1))) then
                write(6,*) 'Proc', Proc
                write(6,*) 'Sing_ExcList', Sing_ExcList
                write(6,*) 'No. spaces for each proc', nint(OneEl_Gap)
                call Stop_All('GenExcDjs',&
                            'Too many excitations for space available.')
            end if
        end do

        if (RDMExcitLevel .ne. 1) then            

            ExcitMat3(:,:) = 0
            ! Zeros in ExcitMat3 starts off at the first single excitation.
            tAllExcitFound=.false.
            ! This becomes true when all the excitations have been found.

            do while (.not. tAllExcitFound)
                exflag = 2
                call GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
                                                            tAllExcitFound,.true.)            

                ! Passed out of here is the doubly excited determinant, nJ.
                ! Information such as the orbitals involved in the excitation
                ! and the parity is also found in this step, we are not currently
                ! storing this, and it is re-calculated  later on (after the
                ! determinants are passed to the relevant processor) - but the
                ! speed of sending this information vs recalculating it will be
                ! tested. RDMExcitLevel is passed through, if this is 1, only
                ! singles are generated, if it is 2 only doubles are found.

                if (tAllExcitFound) exit

                iLutnJ(:) = 0
                call EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nel,nJ,0)   
                ! This will return a value between 0 -> nProcessors-1.
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being
                ! stored in Doub_ExcDjs.

                Doub_ExcList(Proc) = Doub_ExcList(Proc) + 1

                ! Want a quick test to see if arrays are getting full.            
                if (Doub_ExcList(Proc).gt.nint(TwoEl_Gap*(Proc+1))) then
                    write(6,*) 'Proc',Proc
                    write(6,*) 'Doub_ExcList',Doub_ExcList
                    write(6,*) 'No. spaces for each proc',nint(TwoEl_Gap)
                    call Stop_All('GenExcDjs','Too many excitations for space available.')
                end if
            end do
        end if

    end subroutine Gen_Hist_ExcDjs

    subroutine SendProcExcDjs()

        ! In this routine the excitations are sent to the relevant processors.
        ! Sent with them will be the Di they were excited from and its sign.
        ! Each processor will receive nProcessor number of lists with different
        ! Di determinants. The original Di's will (I think) still be in the
        ! original InitSingExcSlots positions. This follows the
        ! directannihilation algorithm closely.

        integer :: i,j
        integer(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        integer(MPIArg) :: sing_recvcounts(nProcessors)
        integer :: error,MaxSendIndex,MaxIndex
        integer(MPIArg) :: sing_recvdisps(nProcessors)
        integer(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i = 0, nProcessors-1
            sendcounts(i+1) = int(Sing_ExcList(i)-(nint(OneEl_Gap*i)+1),MPIArg)
            ! Sendcounts is the number of singly excited determinants we want
            ! to send for each processor (but goes from 1, not 0).
            disps(i+1) = nint(OneEl_Gap*i,MPIArg)
            ! and I think disps is the first slot for each processor - 1.            
        end do

        MaxSendIndex = Sing_ExcList(nProcessors-1)-1

        ! We now need to calculate the recvcounts and recvdisps -
        ! this is a job for AlltoAll
        sing_recvcounts(1:nProcessors)=0
        call MPIAlltoAll(sendcounts,1,sing_recvcounts,1,error)
        ! I think recvcounts(i) is the number of determinants sent from
        ! processor i.

        ! We can now get recvdisps from recvcounts, since we want the data to
        ! be contiguous after the move.
        sing_recvdisps(1) = 0
        do i = 2, nProcessors
            sing_recvdisps(i) = sing_recvdisps(i-1) + sing_recvcounts(i-1)
        end do

        MaxIndex = sing_recvdisps(nProcessors) + sing_recvcounts(nProcessors)
        ! But the actual number of integers we need to send is the calculated
        ! values * NIfTot+1.
        do i = 1, nProcessors
            sendcounts(i) = sendcounts(i)*(int(NIfTot+1,MPIArg))
            disps(i) = disps(i)*(int(NIfTot+1,MPIArg))
            sing_recvcounts(i) = sing_recvcounts(i)*(int(NIfTot+1,MPIArg))
            sing_recvdisps(i) = sing_recvdisps(i)*(int(NIfTot+1,MPIArg))
        end do

#ifdef PARALLEL
        call MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                            Sing_ExcDjs2,sing_recvcounts,sing_recvdisps,error)
#else
        Sing_ExcDjs2(0:NIfTot,1:MaxIndex)=Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

        call Sing_SearchOccDets(sing_recvcounts,sing_recvdisps)

        if (RDMExcitLevel .ne. 1) then
            do i = 0, nProcessors-1
                sendcounts(i+1) = int(Doub_ExcList(i)-(nint(TwoEl_Gap*i)+1),MPIArg)
                ! Sendcounts is the number of singly excited determinants we
                ! want to send for
                ! each processor (but goes from 1, not 0).
                disps(i+1) = nint(TwoEl_Gap*i,MPIArg)
                ! and I think disps is the first slot for each processor - 1.            
            end do

            MaxSendIndex = Doub_ExcList(nProcessors-1)-1

            ! We now need to calculate the recvcounts and recvdisps -
            ! this is a job for AlltoAll
            doub_recvcounts(1:nProcessors)=0
            call MPIAlltoAll(sendcounts,1,doub_recvcounts,1,error)
            ! I think recvcounts(i) is the number of determinants sent from
            ! processor i.

            ! We can now get recvdisps from recvcounts, since we want the data
            ! to be contiguous after the move.
            doub_recvdisps(1)=0
            do i=2,nProcessors
                doub_recvdisps(i)=doub_recvdisps(i-1)+doub_recvcounts(i-1)
            end do

            MaxIndex=doub_recvdisps(nProcessors)+doub_recvcounts(nProcessors)
            ! But the actual number of integers we need to send is the
            ! calculated values * NIfTot+1.
            do i=1,nProcessors
                sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
                disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
                doub_recvcounts(i)=doub_recvcounts(i)*(int(NIfTot+1,MPIArg))
                doub_recvdisps(i)=doub_recvdisps(i)*(int(NIfTot+1,MPIArg))
            end do

            ! This is the main send of all the single excitations to the
            ! corresponding processors.
#ifdef PARALLEL
            call MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                                    Doub_ExcDjs2,doub_recvcounts,doub_recvdisps,error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex)=Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            call Doub_SearchOccDets(doub_recvcounts,doub_recvdisps)

        end if
        
    end subroutine SendProcExcDjs

    subroutine Send_Hist_ProcExcDjs()

        ! In this routine the excitations are sent to the relevant processors.
        ! Sent with them will be the Di they were excited from and its sign.
        ! Each processor will receive nProcessor number of lists with different
        ! Di determinants. The original Di's will (I think) still be in the
        ! original InitSingExcSlots positions. This follows the
        ! directannihilation algorithm closely.

        integer :: i,j
        integer(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        integer(MPIArg) :: sing_recvcounts(nProcessors)
        integer :: error,MaxSendIndex,MaxIndex
        integer(MPIArg) :: sing_recvdisps(nProcessors)
        integer(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i = 0, nProcessors-1
            sendcounts(i+1)=int(Sing_ExcList(i)-(nint(OneEl_Gap*i)+1),MPIArg)
            ! Sendcounts is the number of singly excited determinants we want to send for
            ! each processor (but goes from 1, not 0).
            disps(i+1)=nint(OneEl_Gap*i,MPIArg)
            ! and I think disps is the first slot for each processor - 1.
        end do

        MaxSendIndex=Sing_ExcList(nProcessors-1)-1

        ! We now need to calculate the recvcounts and recvdisps -
        ! this is a job for AlltoAll
        sing_recvcounts(1:nProcessors)=0
        call MPIAlltoAll(sendcounts,1,sing_recvcounts,1,error)
        ! I think recvcounts(i) is the number of determinants sent from processor i.

        ! We can now get recvdisps from recvcounts, since we want the data to
        ! be contiguous after the move.
        sing_recvdisps(1) = 0
        do i = 2, nProcessors
            sing_recvdisps(i) = sing_recvdisps(i-1)+sing_recvcounts(i-1)
        end do

        MaxIndex=sing_recvdisps(nProcessors)+sing_recvcounts(nProcessors)
        ! But the actual number of integers we need to send is the calculated
        ! values * NIfTot+1.
        do i=1,nProcessors
            sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
            disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
            sing_recvcounts(i)=sing_recvcounts(i)*(int(NIfTot+1,MPIArg))
            sing_recvdisps(i)=sing_recvdisps(i)*(int(NIfTot+1,MPIArg))
        end do
#ifdef PARALLEL
        call MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                            Sing_ExcDjs2,sing_recvcounts,sing_recvdisps,error)
#else
        Sing_ExcDjs2(0:NIfTot,1:MaxIndex)=Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

        call Sing_Hist_SearchOccDets(sing_recvcounts,sing_recvdisps)


        if (RDMExcitLevel.ne.1) then            
            do i = 0, nProcessors-1
                sendcounts(i+1)=int(Doub_ExcList(i)-(nint(TwoEl_Gap*i)+1),MPIArg)
                ! Sendcounts is the number of singly excited determinants we
                ! want to send for each processor (but goes from 1, not 0).
                disps(i+1)=nint(TwoEl_Gap*i,MPIArg)
                ! and I think disps is the first slot for each processor - 1.            
            end do

            MaxSendIndex = Doub_ExcList(nProcessors-1)-1

            ! We now need to calculate the recvcounts and recvdisps -
            ! this is a job for AlltoAll
            doub_recvcounts(1:nProcessors)=0
            call MPIAlltoAll(sendcounts,1,doub_recvcounts,1,error)
            ! I think recvcounts(i) is the number of determinants sent from
            ! processor i.

            ! We can now get recvdisps from recvcounts, since we want the data
            ! to be contiguous after the move.
            doub_recvdisps(1) = 0
            do i = 2, nProcessors
                doub_recvdisps(i) = doub_recvdisps(i-1) + doub_recvcounts(i-1)
            end do

            MaxIndex= doub_recvdisps(nProcessors) + doub_recvcounts(nProcessors)
            ! But the actual number of integers we need to send is the calculated
            ! values * NIfTot+1.
            do i = 1, nProcessors
                sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
                disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
                doub_recvcounts(i)=doub_recvcounts(i)*(int(NIfTot+1,MPIArg))
                doub_recvdisps(i)=doub_recvdisps(i)*(int(NIfTot+1,MPIArg))
            end do

            ! This is the main send of all the single excitations to the
            ! corresponding processors.
#ifdef PARALLEL
            call MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                                    Doub_ExcDjs2,doub_recvcounts,doub_recvdisps,error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex)=Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            call Doub_Hist_SearchOccDets(doub_recvcounts,doub_recvdisps)

        end if

    end subroutine Send_Hist_ProcExcDjs

    subroutine Sing_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single
        ! excitations from each processor. These number sent from processor i
        ! is recvcounts(i), and the first 2 have information about the
        ! determinant Di from which the Dj's are single excitations (and it's sign).

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(kind=n_int) :: iLutnJ(0:NIfTot)
        real(dp), dimension(lenof_sign) :: SignDi,SignDj, SignDi2,SignDj2
        integer :: i, j, NoDets, StartDets, PartInd
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        real(dp) :: realSignDi, realSignDj

        ! Take each Dj, and binary search CurrentDets to see if it is occupied.

        do i=1,nProcessors
            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            if (NoDets .gt. 1) then
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                realSignDi = real(SignDi(1),dp)

                do j = StartDets+1, (NoDets+StartDets-1)

                    ! D_i is in the first spot - start from the second.                
                    iLutnJ(:)=Sing_ExcDjs2(:,j)
                    ! This binary searches CurrentDets between 1 and
                    ! TotWalkers for determinant iLutnJ. If found, tDetFound
                    ! will be true, and PartInd the index in CurrentDets where
                    ! the determinant is.
                    call BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)
                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and the parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations -
                        ! we know this is an excitation of level RDMExcitLevel. 
                        Ex(1,1) = 1
                        tParity = .false.

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

                        realSignDj = real(SignDj(1),dp)

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to.
                        call GetExcitation(nI,nJ,NEl,Ex,tParity)

                        if (Ex(1,1).le.0) call Stop_All('Sing_SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Sings_RDM(nI, Ex, tParity, realSignDi, realSignDj, .true.)

                        ! No normalisation factor just yet - possibly need to revise.
                    end if

                end do
            end if
        end do
      
    end subroutine Sing_SearchOccDets

    subroutine Doub_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single excitations
        ! from each processor. These number sent from processor i is
        ! recvcounts(i), and the first 2 have information  about the determinant
        ! Di from which the Dj's are single excitations (and it's sign).

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(kind=n_int) :: iLutnJ(0:NIfTot)
        real(dp), dimension(lenof_sign) :: SignDi,SignDj, SignDi2, SignDj2
        integer :: i, j, NoDets, StartDets, PartInd
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        real(dp) :: realSignDi,realSignDj

        ! Take each Dj, and binary search CurrentDets to see if it is occupied.
    
        do i = 1, nProcessors
            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            if (NoDets.gt.1) then
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                realSignDi = real(SignDi(1),dp)

                do j=StartDets+1,(NoDets+StartDets-1)
                    ! D_i is in the first spot - start from the second.
                    iLutnJ(:)=Doub_ExcDjs2(:,j)

                    ! This binary searches CurrentDets between 1 and TotWalkers
                    ! for determinant iLutnJ. If found, tDetFound will be true,
                    ! and PartInd the index in CurrentDets where the determinant is.
                    call BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and 
                        ! the parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations -
                        ! we know this is an excitation of level RDMExcitLevel. 
                        Ex(1,1) = 2
                        tParity = .false.

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

                        realSignDj = real(SignDj(1),dp)

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to. 
                        call GetExcitation(nI,nJ,NEl,Ex,tParity)

                        if (Ex(1,1).le.0) call Stop_All('SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Doubs_RDM(Ex,tParity,realSignDi,realSignDj,.true.)
                        
                        
                    end if
                end do
            end if
        end do
      
    end subroutine Doub_SearchOccDets

    subroutine Sing_Hist_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single
        ! excitations from each processor. These number sent from processor i
        ! is recvcounts(i), and the first 2 have information about the
        ! determinant Di from which the Dj's are single excitations (and it's sign).

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(kind=n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: HistPos
        real(dp), dimension(lenof_sign) :: RealHistPos

        integer :: i, j, NoDets, StartDets, PartInd, ExcitLevel
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        real(dp) :: realSignDi, realSignDj

        ! Take each Dj, and binary search CurrentDets to see if it is occupied.

        do i = 1, nProcessors

            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            if (NoDets.gt.1) then
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, RealHistPos, FlagsDi)

                HistPos=int(RealHistPos)

                realSignDi = AllHistogram(1,HistPos(1))/norm

                do j=StartDets+1,(NoDets+StartDets-1)

                    ! D_i is in the first spot - start from the second.                
                    iLutnJ(:)=Sing_ExcDjs2(:,j)
                    ! This binary searches CurrentDets between 1 and TotWalkers
                    ! for determinant iLutnJ. If found, tDetFound will be true,
                    ! and PartInd the index in CurrentDets where the determinant is.
                    call BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF_true, iLutnJ, NEl)
                    call find_hist_coeff_explicit (iLutnJ, ExcitLevel,PartInd,tDetFound)

                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and the
                        ! parity.
                        Ex(:,:) = 0
                        ! Ex(1,1) goes in as the max number of excitations - we
                        ! know this is an excitation of level RDMExcitLevel. 
                        Ex(1,1) = 1
                        tParity = .false.

                        call decode_bit_det(nJ,iLutnJ)

                        realSignDj = AllHistogram(1,PartInd)/norm

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to.    
                        call GetExcitation(nI,nJ,NEl,Ex,tParity)

                        if (Ex(1,1).le.0) call Stop_All('Sing_SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Sings_RDM(nI,Ex,tParity,realSignDi,realSignDj,.true.)

                        ! No normalisation factor just yet - possibly need to revise.                    
                    end if

                end do
            end if
        end do
      
    end subroutine Sing_Hist_SearchOccDets

    subroutine Doub_Hist_SearchOccDets(recvcounts,recvdisps)

        ! We now have arrays SingExcDjs2 which contain all the single excitations
        ! from each processor. These number sent from processor i is recvcounts(i),
        ! and the first 2 have information about the determinant Di from which
        ! the Dj's are single excitations (and it's sign).

        integer(MPIArg), intent(in) :: recvcounts(nProcessors),recvdisps(nProcessors)
        integer(kind=n_int) :: iLutnJ(0:NIfTot)
        integer, dimension(lenof_sign) :: HistPos
        real(dp), dimension(lenof_sign) :: RealHistPos
        integer :: i, j, NoDets, StartDets,PartInd, ExcitLevel
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        real(dp) :: realSignDi,realSignDj

        ! Take each Dj, and binary search the CurrentDets to see if it is occupied.
        
        do i=1,nProcessors

            ! Doing determinants from each processor separately because each
            ! has a different D_i it was excited from.
            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            if (NoDets.gt.1) then
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, RealHistPos, FlagsDi)

                HistPos=int(RealHistPos)

                realSignDi = AllHistogram(1,HistPos(1))/norm

                do j=StartDets+1,(NoDets+StartDets-1)

                    ! D_i is in the first spot - start from the second. 
                    iLutnJ(:)=Doub_ExcDjs2(:,j)

                    ! This binary searches CurrentDets between 1 and TotWalkers
                    ! for determinant iLutnJ. If found, tDetFound will be true,
                    ! and PartInd the index in CurrentDets where the  determinant is.
                    call BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF_True, iLutnJ, NEl)
                    call find_hist_coeff_explicit (iLutnJ, ExcitLevel, PartInd, tDetFound)

                    if (tDetFound) then
                        ! Determinant occupied; add c_i*c_j to the relevant
                        ! element of nElRDM. Need to first find the orbitals
                        ! involved in the excitation from D_i -> D_j and
                        ! the parity.
                        Ex(:,:)=0
                        ! Ex(1,1) goes in as the max number of excitations - we
                        ! know this is an excitation of level RDMExcitLevel. 
                        Ex(1,1)=2
                        tParity = .false.

                        call decode_bit_det(nJ,iLutnJ)
                        realSignDj = AllHistogram(1,PartInd)/norm

                        ! Ex(1,:) comes out as the orbital(s) excited from,
                        ! Ex(2,:) comes out as the orbital(s) excited to. 
                        call GetExcitation(nI,nJ,NEl,Ex,tParity)

                        if (Ex(1,1).le.0) call Stop_All('SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')
                        call Fill_Doubs_RDM(Ex,tParity,realSignDi,realSignDj,.true.)
                        
                        
                    end if
                end do
            end if
        end do
      
    end subroutine Doub_Hist_SearchOccDets

! =======================================================================================    
! THESE NEXT ROUTINES ARE GENERAL TO BOTH STOCHASTIC AND EXPLICIT    
! =======================================================================================    

    subroutine Fill_Diag_RDM(nI,realSignDi,tCoreSpaceDet, RDMItersIn)

        ! Fill diagonal elements of 1- and 2-RDM.
        ! These are < Di | a_i+ a_i | Di > and < Di | a_i+ a_j+ a_j a_i | Di >.

        integer, intent(in) :: nI(NEl)
        real(dp), dimension(lenof_sign), intent(in) :: realSignDi
        logical, intent(in), optional :: tCoreSpaceDet
        integer, intent(in), optional :: RDMItersIn
        integer :: i, j, iSpat, jSpat, Ind, iInd
        real(dp) :: ScaleContribFac
        integer :: RDMIters

        ! Need to add in the diagonal elements.
        
        ScaleContribFac=1.0
        
        if (.not.present(RDMItersIn)) then
            RDMIters=1.0_dp
        else
            RDMIters=RDMItersIn
        end if

        ! This is the single-run cutoff being applied (do not use in DR mode):
        if ((.not. tCoreSpaceDet) .or. .not.present(tCoreSpaceDet)) then
            ! Dets in the core space are never removed from main list, so
            ! strictly do not require corrections
            if (tThreshOccRDMDiag .and. (abs(RealSignDi(1)) .le. ThreshOccRDM)) ScaleContribFac=0.0_dp
        end if
        
        if (RDMExcitLevel.eq.1) then
            do i=1,NEl
                if (tOpenShell) then
                    iInd = SymLabelListInv_rot(nI(i))
                else 
                    ! SymLabelListInv_rot will be in spat orbitals too.
                    iInd = SymLabelListInv_rot(gtID(nI(i)))
                end if
                NatOrbMat(iInd,iInd) = NatOrbMat(iInd,iInd) &
                                          + ( realSignDi(1) * realSignDi(lenof_sign) * RDMIters )*ScaleContribFac 
            end do
        else
            ! Only calculating 2-RDM.
            ! nI(i) - spin orbital label. Odd=beta, even=alpha.
            do i=1,NEl - 1
                iSpat = gtID(nI(i))
                if (tOpenShell) iSpat = (nI(i)-1)/2 + 1

                ! Orbitals in nI ordered lowest to highest so nI(j) > nI(i),
                ! and jSpat >= iSpat (can only be equal if different spin).
                do j=i+1,NEl
                    jSpat = gtID(nI(j))
                     if (tOpenShell) jSpat = (nI(j)-1)/2 + 1 
               
                    ! either alpha alpha or beta beta -> aaaa/bbbb arrays.
                    if ( ((mod(nI(i),2).eq.1).and.(mod(nI(j),2).eq.1)) .or. &
                        ((mod(nI(i),2).eq.0).and.(mod(nI(j),2).eq.0)) ) then

                        ! Ind doesn't include diagonal terms (when iSpat == jSpat).
                        Ind=( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                        if (( mod(nI(i),2).eq.0) .or. (.not. tOpenShell))then
                            ! nI(i) is even --> aaaa.
                            aaaa_RDM( Ind, Ind ) = aaaa_RDM( Ind, Ind ) &
                                          + ( realSignDi(1) * realSignDi(lenof_sign) * RDMIters)*ScaleContribFac

                        else if ( mod(nI(i),2).eq.1)then
                            ! nI(i) is odd --> bbbb.
                            bbbb_RDM( Ind, Ind ) = bbbb_RDM( Ind, Ind ) &
                                          + ( realSignDi(1) * realSignDi(lenof_sign) * RDMIters)*ScaleContribFac

                        end if
                    ! either alpha beta or beta alpha -> abab/baba arrays.                                              
                    else

                        ! Ind does include diagonal terms (when iSpat == jSpat).
                        Ind=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat

                        if (jSpat.eq.iSpat)then
                                ! aSpat == bSpat == iSpat == jSpat terms are
                                ! saved in abab only.
                                abab_RDM( Ind, Ind ) = abab_RDM( Ind, Ind ) &
                                                + ( realSignDi(1) * realSignDi(lenof_sign) * RDMIters)*ScaleContribFac

                        else 

                            if ((mod(nI(i),2).eq.0) .or. (.not. tOpenShell))then
                                ! nI(i) is even ---> abab.
                                abab_RDM( Ind, Ind ) = abab_RDM( Ind, Ind ) &
                                          + ( realSignDi(1) * realSignDi(lenof_sign) * RDMIters)*ScaleContribFac
                            else if (mod(nI(i),2).eq.1)then
                                ! nI(i) is odd ---> baba.
                                baba_RDM( Ind, Ind ) = baba_RDM( Ind, Ind ) &
                                               + ( realSignDi(1) * realSignDi(lenof_sign) * RDMIters)*ScaleContribFac
                            end if

                       end if 

                    end if

                end do
            end do
        end if

    end subroutine Fill_Diag_RDM

    subroutine Fill_Sings_RDM(nI,Ex,tParity,realSignDi,realSignDj,tFill_CiCj_Symm)

        ! This routine adds in the contribution to the 1- and 2-RDM from
        ! determinants connected by a single excitation.

        integer, intent(in) :: nI(NEl), Ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: realSignDi, realSignDj
        logical, intent(in) :: tFill_CiCj_Symm
        integer :: k, Indik, Indak, iSpat, aSpat, kSpat, iInd, aInd
        real(dp) :: ParityFactor, ParityFactor2

        ParityFactor=1.0_dp
        if (tParity) ParityFactor=-1.0_dp

        if (RDMExcitLevel.eq.1) then

            ! SymLabelList2_rot(i) gives the orbital in position i
            ! SymLabelListInv_rot(i) gives the position orbital i should go in.
            if (tOpenShell) then
                iInd = Ex(1,1)
                aInd = Ex(2,1)
            else
                iInd = gtID(Ex(1,1))
                aInd = gtID(Ex(2,1))   ! These two must have the same spin.
            end if
            Indik = SymLabelListInv_rot(iInd)    ! Position of i 
            Indak = SymLabelListInv_rot(aInd)    ! Position of a.
            
            ! Adding to 1-RDM(i,a), ci.cj effectively.
            NatOrbMat( Indik, Indak ) = NatOrbMat( Indik, Indak ) + (ParityFactor * &
                                                             realSignDi * realSignDj )

            if (tFill_CiCj_Symm) then                                
                NatOrbMat( Indak, Indik ) = NatOrbMat( Indak, Indik ) + (ParityFactor * &
                                                             realSignDi * realSignDj )
            end if
        else
            ! Looking at elements of the type Gamma(i,k,a,k)

            ! The two determinants Di and Dj will have the same occupations
            ! except for the i and a. Any of the N-1 other electrons can be
            ! annihilated and created in the same orbital. So we run over
            ! all k = all N-1 other occupied orbitals.
            
            iSpat = gtID(Ex(1,1))
            aSpat = gtID(Ex(2,1))  ! These two must have the same spin.
            if (tOpenShell)then
                iSpat = (Ex(1,1)-1)/2 + 1
                aSpat = (Ex(2,1)-1)/2 + 1
            end if

            do k = 1, NEl                            

                kSpat = gtID(nI(k))
                if (tOpenShell) kSpat = (nI(k)-1)/2 + 1

                if (nI(k).ne.Ex(1,1)) then

                    if ((iSpat.eq.kSpat).or.(aSpat.eq.kSpat)) then
                        ! It is possible for i = k or a = k if they 
                        ! have different spins. the only arrays with
                        ! i = j or a = b are abab/baba.
                        ! -> abba/baab terms must be reordered to
                        ! become abab/baba

                        Indik=( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                        Indak=( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                        if ((iSpat .eq. aSpat) .or. (.not. tOpenShell) )then

                            ! The iSpat == jSpat == aSpat == bSpat term is
                            ! saved in the abab array only (not in baba)
                            abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                            if (tFill_CiCj_Symm) then
                                abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                            end if

                        else if (iSpat.eq.kSpat) then

                            ! We get the term k i -> k a, which is abab or baba.
                            ! If iSpat == kSpat and aSpat > kSpat, the indeces are
                            ! already ordered correctly. If they are not ordered
                            ! correctly (aSpat<kSpat), then we need to swap a and k
                            ! This results into an abba/baab term, but for equal
                            ! spatial orbitals, there is no abba/baab array, and
                            ! so we save the equivalent abab/baba term. Therefore
                            ! i and k have to be swapped as well, giving an abab/baba
                            ! term. Because there are none or two swaps, the parity does
                            ! not change!

                            if (aSpat .gt. kSpat) then

                                if ( (mod(Ex(2,1),2).eq.1) .or. (.not. tOpenShell) )then ! ki, ka -> last index beta
                                    abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if

                                else if (mod(Ex(2,1),2).eq.0)then ! ki, ka -> last index alpha
                                    baba_RDM( Indik, Indak ) = baba_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        baba_RDM( Indak, Indik ) = baba_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            else if  (aSpat .lt. kSpat) then

                                if ( (mod(Ex(2,1),2).eq.0) .or. (.not. tOpenShell) )then ! ik, ak -> third index alpha
                                    abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(2,1),2).eq.1)then ! ik, ak -> third index beta
                                    baba_RDM( Indik, Indak ) = baba_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        baba_RDM( Indak, Indik ) = baba_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            end if

                        else if (aSpat.eq.kSpat) then

                            if (iSpat .gt. kSpat) then 
                                if ( (mod(Ex(1,1),2) .eq. 1) .or. (.not. tOpenShell) )then ! ki, ka -> second index beta
                                    abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(1,1),2).eq.0)then ! ki, ka -> second index alpha
                                    baba_RDM( Indik, Indak ) = baba_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        baba_RDM( Indak, Indik ) = baba_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            else if  (iSpat .lt. kSpat) then

                                if ( (mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell) )then ! ik, ak -> first index alpha
                                    abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(1,1),2).eq.1)then ! ik, ak -> first index beta
                                    baba_RDM( Indik, Indak ) = baba_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        baba_RDM( Indak, Indik ) = baba_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                            realSignDi * realSignDj )
                                    end if
                                end if

                            end if
                        end if  ! a=k

                    else ! not (iSpat.eq.kSpat) .or. (aSpat.eq.kSpat))
                        ! Checking spins of i and k.
                        ! If same, i.e alpha alpha or beta beta -> aaaa array.
                        if ( ((mod(Ex(1,1),2).eq.1).and.(mod(nI(k),2).eq.1)) .or. &
                            ((mod(Ex(1,1),2).eq.0).and.(mod(nI(k),2).eq.0)) ) then

                            ! 2-RDM(i,j,a,b) is defined to have i < j and a < b, as that is how the unique 
                            ! indices are defined for i,j and a,b.
                            ! But the parity is defined so that the i -> a excitation is aligned.

                            ! I.e. we're adding these as nI(k),Ex(1,1) -> nI(k), Ex(2,1)
                            ! So if Ex(1,1) < nI(k), or Ex(2,1) < nI(k) then we need 
                            ! to switch the parity.
                            ParityFactor2 = ParityFactor
                            if ((Ex(1,1).lt.nI(k)).and.(Ex(2,1).gt.nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.0_dp)
                            if ((Ex(1,1).gt.nI(k)).and.(Ex(2,1).lt.nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.0_dp)

                            ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                            Indik=( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                            Indak=( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                            ! nI(k) even or odd. odd=bbbb, even=aaaa.                          
                            if ((mod(nI(k),2).eq.0) .or. (.not. tOpenShell)) then
                                aaaa_RDM( Indik, Indak ) = aaaa_RDM( Indik, Indak ) + ( ParityFactor2 * &
                                                                                    realSignDi * realSignDj )
                                if (tFill_CiCj_Symm) then
                                    aaaa_RDM( Indak, Indik ) = aaaa_RDM( Indak, Indik ) + ( ParityFactor2 * &
                                                                                    realSignDi * realSignDj )
                                end if

                            else if (mod(nI(k),2).eq.1)then
                                bbbb_RDM( Indik, Indak ) = bbbb_RDM( Indik, Indak ) + ( ParityFactor2 * &
                                                                                 realSignDi * realSignDj )
                                if (tFill_CiCj_Symm) then
                                        bbbb_RDM( Indak, Indik ) = bbbb_RDM( Indak, Indik ) + ( ParityFactor2 * &
                                                                                 realSignDi * realSignDj )
                                end if
                            end if

                        ! either abab/baba or abba/baab array. 
                        ! we distinguish between these because i<j and a<b.
                        else   ! abab/baba or abba/baab

                            if ( (Ex(1,1).lt.nI(k)).and.(Ex(2,1).lt.nI(k)) )then
                            !i k a k -> abab/baba 

                                ! It is possible for i = k or j = k if they are spat orbitals 
                                ! and have different spins.
                                Indik=( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                                Indak=( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                                !Ex(1,1) (i spin orb): first index even or odd 
                                if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell)) then
                                    abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                    end if
                                else if (mod(Ex(1,1),2).eq.1)then
                                    baba_RDM( Indik, Indak ) = baba_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                            baba_RDM( Indak, Indik ) = baba_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj ) 
                                    end if
                                end if

                            else if ( (Ex(1,1).gt.nI(k)).and.(Ex(2,1).gt.nI(k)) ) then
                            ! k i k a -> abab/baba 

                                Indik=( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                                Indak=( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                                !Ex(1,1) (i spin orb): second index even or odd 
                                if ((mod(Ex(1,1),2).eq.1) .or. (.not. tOpenShell)) then
                                    abab_RDM( Indik, Indak ) = abab_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                    if (tFill_CiCj_Symm) then
                                        abab_RDM( Indak, Indik ) = abab_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj ) 
                                    end if
                                else if (mod(Ex(1,1),2).eq.0)then
                                    baba_RDM( Indik, Indak ) = baba_RDM( Indik, Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                    if (tFill_CiCj_Symm) then
                                            baba_RDM( Indak, Indik ) = baba_RDM( Indak, Indik ) + ( ParityFactor * &
                                                                                    realSignDi * realSignDj )  
                                    end if
                                end if


                            else if ( (Ex(1,1).gt.nI(k)).and.(Ex(2,1).lt.nI(k)) ) then 
                            ! k i a k -> abba/baab

                                ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                                Indik=( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                                Indak=( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                                !i spin orb: second odd or even. 
                                if ((mod(Ex(1,1),2).eq.1) .or. (.not. tOpenShell))then
                                    abba_RDM( Indik, Indak ) = abba_RDM( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                    if (tFill_CiCj_Symm ) then
                                        if (.not. tOpenShell) then
                                            abba_RDM( Indak, Indik ) = abba_RDM( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        else
                                            baab_RDM( Indak, Indik ) = baab_RDM( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        end if
                                    end if

                                else if (mod(Ex(1,1),2).eq.0)then
                                    baab_RDM( Indik, Indak ) = baab_RDM( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abba_RDM( Indak, Indik ) = abba_RDM( Indak, Indik ) - ( ParityFactor * &
                                                                                     realSignDi * realSignDj )
                                    end if
                                end if

                            else if ( (Ex(1,1).lt.nI(k)).and.(Ex(2,1).gt.nI(k)) ) then 
                            ! i k k a -> abba/baab

                                ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                                Indik=( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                                Indak=( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                                !i spin orb: first index odd or even.
                                if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell))then
                                    abba_RDM( Indik, Indak ) = abba_RDM( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        if (.not. tOpenShell) then
                                            abba_RDM( Indak, Indik ) = abba_RDM( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        else
                                            baab_RDM( Indak, Indik ) = baab_RDM( Indak, Indik ) - ( ParityFactor * &
                                                                                    realSignDi * realSignDj )
                                        end if
                                    end if

                                else if (mod(Ex(1,1),2).eq.1)then
                                    baab_RDM( Indik, Indak ) = baab_RDM( Indik, Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    if (tFill_CiCj_Symm) then
                                        abba_RDM( Indak, Indik ) = abba_RDM( Indak, Indik ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                    end if
                                end if

                            end if ! order of k i k j

                        end if !abab/baba or abba/baab
                    end if  ! not (iSpat.eq.kSpat).or.(aSpat.eq.kSpat))

                end if 
            end do 

        end if  

    end subroutine Fill_Sings_RDM

    subroutine Fill_Doubs_RDM(Ex,tParity,realSignDi,realSignDj,tFill_CiCj_Symm)

        ! This routine adds in the contribution to the 2-RDM from determinants
        ! connected by a double excitation.

        integer, intent(in) :: Ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: realSignDi, realSignDj
        logical, intent(in) :: tFill_CiCj_Symm
        integer :: Indij, Indab, iSpat, jSpat, aSpat, bSpat
        real(dp) :: ParityFactor

        ! Adding to elements Gamma(i,j,a,b)

        ParityFactor=1.0_dp
        if (tParity) ParityFactor=-1.0_dp

        iSpat = gtID(Ex(1,1))
        jSpat = gtID(Ex(1,2))       ! Ex(1,1) < Ex(1,2)
        aSpat = gtID(Ex(2,1)) 
        bSpat = gtID(Ex(2,2))       ! Ex(2,1) < Ex(2,2)

        if (tOpenShell)then 
            iSpat = (Ex(1,1)-1)/2 + 1
            jSpat = (Ex(1,2)-1)/2 + 1
            aSpat = (Ex(2,1)-1)/2 + 1
            bSpat = (Ex(2,2)-1)/2 + 1
        end if         

        if ((iSpat.eq.jSpat).or.(aSpat.eq.bSpat)) then

            ! if i and a are different spin -> abba (but adding as abab - mult by -1).
            if ( ((mod(Ex(1,1),2).eq.0).and.(mod(Ex(2,1),2).eq.1)) .or. &
                ((mod(Ex(1,1),2).eq.1).and.(mod(Ex(2,1),2).eq.0)) ) &
                    ParityFactor = ParityFactor * (-1.0_dp)

            Indij=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
            Indab=( ( (bSpat-1) * bSpat ) / 2 ) + aSpat

            if ((iSpat.eq.jSpat).and.(aSpat.eq.bSpat))then
                ! abab and baba terms are equal and are saved in abab only.

                abab_RDM( Indij, Indab ) = abab_RDM( Indij, Indab ) + ( ParityFactor * &
                                                         realSignDi * realSignDj )
                if (tFill_CiCj_Symm) then
                    abab_RDM( Indab, Indij ) = abab_RDM( Indab, Indij ) + ( ParityFactor * &
                                                            realSignDi * realSignDj )
                end if

            else if (iSpat.eq.jSpat)then
               ! i and j may have to be swapped to get abab/baba term -> get
               ! spin from a (third index).
               if ((mod(Ex(2,1),2).eq.0) .or. (.not. tOpenShell))then
                    abab_RDM( Indij, Indab ) = abab_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        abab_RDM( Indab, Indij ) = abab_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               else if (mod(Ex(2,1),2).eq.1)then
                    baba_RDM( Indij, Indab ) = baba_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                  realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        baba_RDM( Indab, Indij ) = baba_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               end if

            else if (aSpat.eq.bSpat)then

               ! a and b may have to be swapped to get abab/baba term -> get
               ! spin from i (first index)
               if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell))then
                    abab_RDM( Indij, Indab ) = abab_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        abab_RDM( Indab, Indij ) = abab_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               else if (mod(Ex(1,1),2).eq.1)then
                    baba_RDM( Indij, Indab ) = baba_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                  realSignDi * realSignDj )
                   if (tFill_CiCj_Symm) then
                        baba_RDM( Indab, Indij ) = baba_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                   end if
               end if

            end if

          else ! not ((iSpat.eq.jSpat) .or. (aSpat.eq.bSpat))

            ! Checking spins of i and j (these must be same combination as a and b).
            ! If alpha alpha or beta beta -> aaaa array.
            if ( ((mod(Ex(1,1),2).eq.1).and.(mod(Ex(1,2),2).eq.1)) .or. &
                ((mod(Ex(1,1),2).eq.0).and.(mod(Ex(1,2),2).eq.0)) ) then

                ! Don't need to worry about diagonal terms, i can't equal j.
                ! jSpat > iSpat and bSpat > aSpat
                Indij = ( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                Indab = ( ( (bSpat-2) * (bSpat-1) ) / 2 ) + aSpat

                if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell))then
                    aaaa_RDM( Indij, Indab ) = aaaa_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                    if (tFill_CiCj_Symm) then
                            aaaa_RDM( Indab, Indij ) = aaaa_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                    end if

                else if (mod(Ex(1,1),2).eq.1)then
                    bbbb_RDM( Indij, Indab ) = bbbb_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                    if (tFill_CiCj_Symm) then
                          bbbb_RDM( Indab, Indij ) = bbbb_RDM( Indab, Indij ) + (ParityFactor * &
                                                                    realSignDi * realSignDj )
                    end if
                end if
                
            ! Either alpha beta or beta alpha -> abab array.
            else
                
                ! if when ordering i < j and a < b, is it abab or abba.

                ! i and a are the same spin -> abab/baba
                if ( ((mod(Ex(1,1),2).eq.0).and.(mod(Ex(2,1),2).eq.0)) .or. &
                    ((mod(Ex(1,1),2).eq.1).and.(mod(Ex(2,1),2).eq.1)) ) then

                    Indij=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
                    Indab=( ( (bSpat-1) * bSpat ) / 2 ) + aSpat

                    if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell)) then
                        abab_RDM( Indij, Indab ) = abab_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                            abab_RDM( Indab, Indij ) = abab_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                        end if
                    else if (mod(Ex(1,1),2).eq.1)then
                        baba_RDM( Indij, Indab ) = baba_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                                baba_RDM( Indab, Indij ) = baba_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                        end if
                    end if

                ! i and a are different spin -> abba
                ! the only double excitation case with Indij = Indab will go
                ! in here.
                else

                    ! Don't need to worry about diagonal terms, i can't equal j.
                    ! jSpat > iSpat and bSpat > aSpat
                    Indij=( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                    Indab=( ( (bSpat-2) * (bSpat-1) ) / 2 ) + aSpat

                    if ((mod(Ex(1,1),2).eq.0) .or. (.not. tOpenShell))then
                        abba_RDM( Indij, Indab ) = abba_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                            if (.not. tOpenShell)then
                                abba_RDM( Indab, Indij ) = abba_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                        realSignDi * realSignDj )
                            else
                                baab_RDM( Indab, Indij ) = baab_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                        realSignDi * realSignDj )
                            end if
                        end if

                    else if (mod(Ex(1,1),2).eq.1)then

                        baab_RDM( Indij, Indab ) = baab_RDM( Indij, Indab ) + ( ParityFactor * &
                                                                realSignDi * realSignDj )
                        if (tFill_CiCj_Symm) then
                                abba_RDM( Indab, Indij ) = abba_RDM( Indab, Indij ) + ( ParityFactor * &
                                                                    realSignDi * realSignDj )
                        end if

                    end if

                end if
            end if
        end if

    end subroutine Fill_Doubs_RDM

    subroutine FinaliseRDM()

        ! This routine finalises the one electron reduced density matrix stuff
        ! at the point of a softexit. This includes summing each of the
        ! individual matrices from each processor, and calling the
        ! diagonalisation routines if we want to get the occupation numbers.

        integer :: error
        real(dp) :: Norm_2RDM, Norm_2RDM_Inst
        real(dp) :: Norm_1RDM, Trace_1RDM, SumN_Rho_ii
        character(len=*), parameter :: this_routine='FinaliseRDM'

        call set_timer(FinaliseRDM_Time)

        write(6,*) ''

        if (tExplicitAllRDM) then
            write(6,*) '**** RDMs CALCULATED EXPLICITLY **** '
        else
            write(6,*) '**** RDMs CALCULATED STOCHASTIcallY **** '
        end if

        write(6,*) ''

        ! Combine the 1- or 2-RDM from all processors etc.

        if (RDMExcitLevel.eq.1) then
            call Finalise_1e_RDM(Norm_1RDM)  
        else
            ! We always want to calculate one final RDM energy, whether or not we're 
            ! calculating the energy throughout the calculation.
            ! Unless of course, only the 1-RDM is being calculated.

            ! Calculate the energy one last time - and write out everything we need.
            tFinalRDMEnergy = .true.

            !1RDM is contructed here (in calc_1RDM_energy)
            call Calc_Energy_from_RDM(Norm_2RDM)

            if (tPrint1RDM) then
                call Finalise_1e_RDM(Norm_1RDM)
            else if (tDiagRDM.and.(iProcIndex.eq.0)) then
                call calc_1e_norms(Trace_1RDM, Norm_1RDM, SumN_Rho_ii)
                write(6,*) ''
                write(6,'(A55,F30.20)') ' SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY &
                                            &HF ORBITALS: ',SumN_Rho_ii
            end if
            if (tDumpForcesInfo) then
                if (.not. tPrint1RDM) call Finalise_1e_RDM(Norm_1RDM)
                call Calc_Lagrangian_from_RDM(Norm_1RDM, Norm_2RDM)
                call convert_mats_Molpforces(Norm_1RDM, Norm_2RDM)
            end if

        end if

        call MPIBarrier(error)

        ! Call the routines from NatOrbs that diagonalise the one electron
        ! reduced density matrix.
        tRotatedNOs = .false. ! Needed for BrokenSymNo routine
        if (tDiagRDM) call find_nat_orb_occ_numbers()

        ! This is where we would likely call any further calculations of
        ! forces, etc.
        if (tDipoles) then
            if (.not. tPrint1RDM) call Finalise_1e_RDM(Norm_1RDM)
            call CalcDipoles(Norm_1RDM)
        end if

        ! After all the NO calculations are finished we'd like to do another
        ! rotation to obtain symmetry-broken natural orbitals
        if (tBrokenSymNOs) then
            call BrokenSymNO(occ_numb_diff)
        end if

        call halt_timer(FinaliseRDM_Time)
    
    end subroutine FinaliseRDM

    subroutine Finalise_1e_RDM(Norm_1RDM) 

        ! This routine takes the 1-RDM (NatOrbMat), normalises it, makes it 
        ! hermitian if required, and prints out the versions we're interested
        ! in. This is only ever called at the very end of a calculation.

        use Logging, only: twrite_RDMs_to_read, twrite_normalised_RDMs
                             
        integer :: i, ierr
        real(dp), intent(out) :: Norm_1RDM
        real(dp) :: Trace_1RDM, SumN_Rho_ii
        real(dp), allocatable :: AllNode_NatOrbMat(:,:)

        Norm_1RDM = 0.0_dp
        AllAccumRDMNorm = 0.0_dp

        if (RDMExcitLevel.eq.1) then

            allocate(AllNode_NatOrbMat(NoOrbs,NoOrbs),stat=ierr)
            
            call MPISumAll(NatOrbMat,AllNode_NatOrbMat)
            NatOrbMat=AllNode_NatOrbMat
            
            deallocate(AllNode_NatOrbMat)

        end if

        if (iProcIndex.eq.0) then 

            ! Find the normalisation.
            call calc_1e_norms(Trace_1RDM, Norm_1RDM, SumN_Rho_ii)

            ! Write out the unnormalised, non-hermitian OneRDM_POPS.
            if (twrite_RDMs_to_read) call Write_out_1RDM(Norm_1RDM,.false.)

            ! Enforce the hermiticity condition.  If the RDMExcitLevel is not 1, the 
            ! 1-RDM has been constructed from the hermitian 2-RDM, so this will not 
            ! be necessary.
            ! The HF_Ref and HF_S_D_Ref cases are not hermitian by definition.
            if (RDMExcitLevel.eq.1) then
                call make_1e_rdm_hermitian(Norm_1RDM)
                
                if (tForceCauchySchwarz)then
                    call Force_Cauchy_Schwarz(Norm_1RDM)
                end if

            end if
            
            ! Write out the final, normalised, hermitian OneRDM.                
            if (twrite_normalised_RDMs) call Write_out_1RDM(Norm_1RDM,.true.)

            write(6,'(A55,F30.20)') ' SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS: ',SumN_Rho_ii

        end if

    end subroutine Finalise_1e_RDM

    subroutine calc_1e_norms(Trace_1RDM, Norm_1RDM, SumN_Rho_ii)

        ! We want to 'normalise' the reduced density matrices. These are not
        ! even close to being normalised at the moment, because of the way
        ! they are calculated on the fly. They should be calculated from a
        ! normalised wavefunction. But we know that the trace of the one
        ! electron reduced density matrix must be equal to the number of the
        ! electrons. We can use this to find the factor we must divide the
        ! 1-RDM through by.

        real(dp), intent(out) :: Trace_1RDM, Norm_1RDM, SumN_Rho_ii
        integer :: i, HFDet_ID, BRR_ID

        Trace_1RDM = 0.0_dp
        Norm_1RDM = 0.0_dp

        do i = 1, NoOrbs
            Trace_1RDM = Trace_1RDM + NatOrbMat(i,i)
        end do

        Norm_1RDM = ( real(NEl,dp) / Trace_1RDM )
        
        ! Need to multiply each element of the 1 electron reduced density matrices 
        ! by NEl / Trace_1RDM,
        ! and then add it's contribution to the energy.
        
        ! Want to sum the diagonal elements of the 1-RDM for the HF orbitals.
        ! Given the HF orbitals, SymLabelListInv_rot tells us their position
        ! in the 1-RDM.
        SumN_Rho_ii = 0.0_dp
        do i = 1, NoOrbs

            ! Rho_ii is the diagonal elements of the 1-RDM. We want this
            ! ordered according to the energy of the orbitals. Brr has the
            ! orbital numbers in order of energy... i.e Brr(2) = the orbital
            ! index with the second lowest energy. Brr is always in spin
            ! orbitals. i gives the energy level, BRR gives the orbital,
            ! SymLabelListInv_rot gives the position of  this orbital in
            ! NatOrbMat.

            if (tDiagRDM) then
                if (tOpenShell) then
                    Rho_ii(i) = NatOrbMat(SymLabelListInv_rot(BRR(i)),SymLabelListInv_rot(BRR(i))) * Norm_1RDM
                else
                    BRR_ID = gtID(BRR(2*i))
                    Rho_ii(i) = NatOrbMat(SymLabelListInv_rot(BRR_ID),SymLabelListInv_rot(BRR_ID)) * Norm_1RDM
                end if
            end if
    
            if (i.le.NEl) then
                if (tOpenShell) then
                    SumN_Rho_ii = SumN_Rho_ii + &
                            ( NatOrbMat(SymLabelListInv_rot(HFDet_True(i)),SymLabelListInv_rot(HFDet_True(i))) &
                                * Norm_1RDM )
                else
                    HFDet_ID = gtID(HFDet_True(i))
                    SumN_Rho_ii = SumN_Rho_ii + &
                            ( NatOrbMat(SymLabelListInv_rot(HFDet_ID),SymLabelListInv_rot(HFDet_ID)) &
                                * Norm_1RDM ) / 2.0_dp
                end if
            end if
        end do

    end subroutine calc_1e_norms

    subroutine make_1e_rdm_hermitian(Norm_1RDM)

        ! Simply average the 1-RDM(i,j) and 1-RDM(j,i) elements which should
        ! be equal in a perfect world.

        real(dp), intent(in) :: Norm_1RDM
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp
        do i = 1, NoOrbs
            do j = i, NoOrbs
                if ((abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
                        (NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i))*Norm_1RDM))).gt.Max_Error_Hermiticity) &
                    Max_Error_Hermiticity = abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
                                                (NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i))*Norm_1RDM))

                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
                                            (NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i))*Norm_1RDM))

                Temp = (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) + &
                        NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)))/2.0_dp

                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = Temp
                NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = Temp
            end do
        end do

        ! Output the hermiticity errors.
        write(6,'(A29,F30.20)') ' MAX ABS ERROR IN 1RDM HERMITICITY', Max_Error_Hermiticity
        write(6,'(A29,F30.20)') ' SUM ABS ERROR IN 1RDM HERMITICITY', Sum_Error_Hermiticity

    end subroutine make_1e_rdm_hermitian

    subroutine Force_Cauchy_Schwarz(Norm_1RDM)

        real(dp), intent(in) :: Norm_1RDM
        integer :: i, j
        real(dp) :: UpperBound

        write(6,*) "Ensuring that Cauchy--Schwarz inequality holds."

        do i = 1, nBasis
            do j = 1, nBasis

                UpperBound = sqrt(NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(i))&
                    *NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(j)))

                if (abs(NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))) .gt. UpperBound)then

                    if (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) .lt. 0.0_dp)then
                        NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = -UpperBound
                    else if (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) .gt. 0.0_dp)then
                        NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = UpperBound
                    end if

                    write(6,*) "Changing element:", i, j
                else
                    cycle
                end if
            end do
        end do

    end subroutine Force_Cauchy_Schwarz

    subroutine Write_out_1RDM(Norm_1RDM,tNormalise)

        ! This routine writes out the OneRDM. If tNormalise is true, we are
        ! printing the normalised, hermitian matrix. Otherwise, Norm_1RDM is
        ! ignored and we print both 1-RDM(i,j) and 1-RDM(j,i) (in binary) 
        ! for the OneRDM_POPS file to be read in in a restart calculation.

        real(dp), intent(in) :: Norm_1RDM
        logical, intent(in) :: tNormalise
        integer :: i, j, iSpat, jSpat
        integer :: OneRDM_unit

        if (tNormalise) then
            ! Haven't got the capabilities to produce multiple 1-RDMs yet.
            write(6,*) 'Writing out the *normalised* 1 electron density matrix to file'
            call neci_flush(6)
            OneRDM_unit = get_free_unit()
            open(OneRDM_unit,file='OneRDM',status='unknown')
        else
            ! Only every write out 1 of these at the moment.
            write(6,*) 'Writing out the *unnormalised* 1 electron density matrix to file for reading in'
            call neci_flush(6)
            OneRDM_unit = get_free_unit()
            open(OneRDM_unit,file='OneRDM_POPS',status='unknown',form='unformatted')
        end if

        ! Currently always printing 1-RDM in spin orbitals.
        do i = 1, nBasis
            do j = 1, nBasis
                if (tOpenShell) then
                    if (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)).ne.0.0_dp) then 
                        if (tNormalise.and.(i.le.j)) then
                            write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) * Norm_1RDM
                        else if (.not.tNormalise) then
                            ! For the pops, we haven't made the 1-RDM hermitian yet, 
                            ! so print both the 1-RDM(i,j) and 1-RDM(j,i) elements.
                            ! This is written in binary.
                            write(OneRDM_unit) i,j,NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))
                        end if
                    end if
                else
                    iSpat = gtID(i)
                    jSpat = gtID(j)
                    if (NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)).ne.0.0_dp) then 
                        if (tNormalise.and.(i.le.j)) then
                            if (((mod(i,2).eq.0).and.(mod(j,2).eq.0)).or.&
                                ((mod(i,2).ne.0).and.(mod(j,2).ne.0))) then
                                write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                    ( NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)) &
                                                                    * Norm_1RDM ) / 2.0_dp
                            end if
                        else if (.not.tNormalise) then
                            ! The popsfile can be printed in spatial orbitals.
                            if ((mod(i,2).eq.0).and.(mod(j,2).eq.0)) then
                                write(OneRDM_unit) iSpat,jSpat, & 
                                    NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)) 
                            end if
                        end if
                    end if
                end if
            end do
        end do

        close(OneRDM_unit)

    end subroutine Write_out_1RDM

    subroutine Finalise_2e_RDM(Norm_2RDM_Inst, Norm_2RDM) 

        ! This routine sums, normalises, hermitian-ises, and prints the 2-RDMs.
        ! This may be called multiple times if we want to print multiple 2-RDMs.

        real(dp), intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        real(dp) :: AllAccumRDMNorm_Inst
        logical :: tmake_herm
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity
        integer :: ierr

        ! If Iter = 0, this means we have just read in the TwoRDM_POPS_a*** matrices into a***_RDM_full, and 
        ! just want to calculate the old energy.
        ! Don't need to do all this stuff here, because a***_RDM will be empty.

        if (((Iter+PreviousCycles).ne.0).and.((.not.tFinalRDMEnergy).or. &
            ((.not. tCalc_RDMEnergy).or.((Iter - VaryShiftIter(1)).le.IterRDMonFly) &
                      & .or.((Iter-VaryShiftIter(inum_runs)).le.IterRDMonFly) &
                      & .or. (mod((Iter+PreviousCycles-IterRDMStart)+1,RDMEnergyIter).ne.0)))) then


            allocate(AllNodes_aaaa_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            allocate(AllNodes_abba_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            allocate(AllNodes_abab_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
            
            ! The aaaa_RDM may be either inst or full, depending on whether we
            ! are calculating inst. energies or not
            call MPISumAll(aaaa_RDM(:,:), AllNodes_aaaa_RDM(:,:))
            call MPISumAll(abab_RDM(:,:), AllNodes_abab_RDM(:,:))
            call MPISumAll(abba_RDM(:,:), AllNodes_abba_RDM(:,:))
            
            aaaa_RDM(:,:) = AllNodes_aaaa_RDM(:,:)
            abab_RDM(:,:) = AllNodes_abab_RDM(:,:)
            abba_RDM(:,:) = AllNodes_abba_RDM(:,:)

            deallocate(AllNodes_aaaa_RDM)
            deallocate(AllNodes_abab_RDM)
            deallocate(AllNodes_abba_RDM)

            if (tOpenShell) then
                allocate(AllNodes_bbbb_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                allocate(AllNodes_baab_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                allocate(AllNodes_baba_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
         
                call MPISumAll(bbbb_RDM(:,:), AllNodes_bbbb_RDM(:,:))
                call MPISumAll(baba_RDM(:,:), AllNodes_baba_RDM(:,:))
                call MPISumAll(baab_RDM(:,:), AllNodes_baab_RDM(:,:))
                
                bbbb_RDM(:,:) = AllNodes_bbbb_RDM(:,:)
                baba_RDM(:,:) = AllNodes_baba_RDM(:,:)
                baab_RDM(:,:) = AllNodes_baab_RDM(:,:)

                deallocate(AllNodes_bbbb_RDM)
                deallocate(AllNodes_baba_RDM)
                deallocate(AllNodes_baab_RDM)
            end if


            ! The TwoElRDM on the root is now the sum of all 'instantaneous' RDMs
            ! (summed over the energy update cycle). Whereas TwoElRDM_full is
            ! accumulated over the entire run.
            if (tRDMInstEnergy .and. (iProcIndex.eq.0)) then
                aaaa_RDM_full(:,:)=aaaa_RDM_full(:,:)+aaaa_RDM(:,:)
                abab_RDM_full(:,:)=abab_RDM_full(:,:)+abab_RDM(:,:)
                abba_RDM_full(:,:)=abba_RDM_full(:,:)+abba_RDM(:,:)
                if (tOpenShell) then
                    bbbb_RDM_full(:,:)=bbbb_RDM_full(:,:)+bbbb_RDM(:,:)
                    baba_RDM_full(:,:)=baba_RDM_full(:,:)+baba_RDM(:,:)
                    baab_RDM_full(:,:)=baab_RDM_full(:,:)+baab_RDM(:,:)
                end if
            end if

            AllAccumRDMNorm_Inst = 0.0_dp
        end if
            
        if (iProcIndex.eq.0) then
            
            ! Calculate the normalisations.
            call calc_2e_norms(AllAccumRDMNorm_Inst, Norm_2RDM_Inst, Norm_2RDM)

            ! There's no need to explicitly make the RDM hermitian here, as the
            ! integrals are already hermitian -- when we calculate the energy,
            ! it comes out in the wash.
            
            ! Print out the relevant 2-RDMs.
            if ( tFinalRDMEnergy .or. &
                ( tWriteMultRDMs .and. (mod((Iter+PreviousCycles - IterRDMStart)+1,IterWriteRDMs).eq.0) ) ) then

                ! ********************************************
                ! SDS:
                ! WARNING: This variable has been set because otherwise
                !          conditional choices are made based on an
                !          uninitialised variable. This was set according to
                !          the current behaviour in the tests, but I have NO
                !          idea if that is correct.
                !          CMO: please advise.
                ! ********************************************
                tmake_herm = .true.

                if (tFinalRDMEnergy) then
                    ! Only ever want to print the POPS 2-RDMs (for reading in) at the end.
                    if (twrite_RDMs_to_read) call Write_out_2RDM(Norm_2RDM,.false.,.false.)
                end if

                ! This writes out the normalised, hermitian 2-RDMs.
                if (twrite_normalised_RDMs) call Write_out_2RDM(Norm_2RDM, .true., tmake_herm)

             end if
        end if

    end subroutine Finalise_2e_RDM

    subroutine calc_2e_norms(AllAccumRDMNorm_Inst, Norm_2RDM_Inst, Norm_2RDM)

        ! We want to 'normalise' the reduced density matrices. These are not
        ! even close to being normalised at the moment, because of the way
        ! they are calculated on the fly. They should be calculated from a
        ! normalised wavefunction.

        ! We also know that the trace of the two electron reduced density
        ! matrix must be equal to the number of electron pairs in the
        ! system = 1/2 N ( N - 1), so we can do the same for the 2RDM.

        real(dp), intent(in) :: AllAccumRDMNorm_Inst
        real(dp), intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        integer :: i

        ! Find the current, unnormalised trace of each matrix.
        ! TODO: This can be merged into the spin averaging when everything is working.

        Trace_2RDM_Inst = 0.0_dp
        Trace_2RDM = 0.0_dp

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            if (i.le.((SpatOrbs*(SpatOrbs-1))/2)) then
                if (tRDMInstEnergy) then
                    Trace_2RDM_Inst = Trace_2RDM_Inst + aaaa_RDM(i,i)
                    if (tOpenShell) Trace_2RDM_Inst = Trace_2RDM_Inst + bbbb_RDM(i,i)
                end if

                Trace_2RDM = Trace_2RDM + aaaa_RDM_full(i,i)
                 if (tOpenShell)   Trace_2RDM = Trace_2RDM + bbbb_RDM_full(i,i) 
            end if

            if (tRDMInstEnergy) then
                Trace_2RDM_Inst = Trace_2RDM_Inst + abab_RDM(i,i)
                if (tOpenShell) Trace_2RDM_Inst = Trace_2RDM_Inst + baba_RDM(i,i)
            end if

            Trace_2RDM = Trace_2RDM + abab_RDM_full(i,i)
            if (tOpenShell) Trace_2RDM = Trace_2RDM + baba_RDM_full(i,i)
        end do

        Norm_2RDM_Inst = 0.0_dp
        Norm_2RDM = 0.0_dp

        if (tRDMInstEnergy) Norm_2RDM_Inst = ( (0.50_dp * (real(NEl,dp) * (real(NEl,dp) - 1.0_dp))) / Trace_2RDM_Inst )
        Norm_2RDM = ( (0.50_dp * (real(NEl,dp) * (real(NEl,dp) - 1.0_dp))) / Trace_2RDM )

        ! Need to multiply each element of the 1 electron reduced density
        ! matrices by NEl / Trace_1RDM, and then add it's contribution to the energy.

    end subroutine calc_2e_norms

    subroutine make_2e_rdm_hermitian(Norm_2RDM, Max_Error_Hermiticity, Sum_Error_Hermiticity)

        ! This averages 2-RDM(i,j;a,b) and 2-RDM(a,b;i,j) or equivalently
        ! 2-RDM(Ind1,Ind2) and 2-RDM(Ind2,Ind1).

        real(dp), intent(in) :: Norm_2RDM
        real(dp), intent(out) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            do j = i+1, ((SpatOrbs*(SpatOrbs+1))/2)

                if ((i.le.((SpatOrbs*(SpatOrbs-1))/2)).and.(j.le.((SpatOrbs*(SpatOrbs-1))/2))) then

                    if ((abs((aaaa_RDM_full(i,j)*Norm_2RDM)-(aaaa_RDM_full(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((aaaa_RDM_full(i,j)*Norm_2RDM)-(aaaa_RDM_full(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((aaaa_RDM_full(i,j)*Norm_2RDM)-(aaaa_RDM_full(j,i)*Norm_2RDM))

                    Temp = (aaaa_RDM_full(i,j) + aaaa_RDM_full(j,i)) / 2.0_dp
                    aaaa_RDM_full(i,j) = Temp
                    aaaa_RDM_full(j,i) = Temp

                    if (tOpenShell)then
                        if ((abs((bbbb_RDM_full(i,j)*Norm_2RDM)-(bbbb_RDM_full(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                            Max_Error_Hermiticity = abs((bbbb_RDM_full(i,j)*Norm_2RDM)-(bbbb_RDM_full(j,i)*Norm_2RDM))

                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((bbbb_RDM_full(i,j)*Norm_2RDM)-(bbbb_RDM_full(j,i)*Norm_2RDM))

                        Temp = (bbbb_RDM_full(i,j) + bbbb_RDM_full(j,i)) / 2.0_dp
                        bbbb_RDM_full(i,j) = Temp
                        bbbb_RDM_full(j,i) = Temp
                    end if

                    if ((abs((abba_RDM_full(i,j)*Norm_2RDM)-(abba_RDM_full(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((abba_RDM_full(i,j)*Norm_2RDM)-(abba_RDM_full(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((abba_RDM_full(i,j)*Norm_2RDM)-(abba_RDM_full(j,i)*Norm_2RDM))

                    Temp = (abba_RDM_full(i,j) + abba_RDM_full(j,i)) / 2.0_dp
                    abba_RDM_full(i,j) = Temp
                    abba_RDM_full(j,i) = Temp

                    if (tOpenShell)then
                        if ((abs((baab_RDM_full(i,j)*Norm_2RDM)-(baab_RDM_full(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                            Max_Error_Hermiticity = abs((baab_RDM_full(i,j)*Norm_2RDM)-(baab_RDM_full(j,i)*Norm_2RDM))

                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((baab_RDM_full(i,j)*Norm_2RDM)-(baab_RDM_full(j,i)*Norm_2RDM))

                        Temp = (baab_RDM_full(i,j) + baab_RDM_full(j,i)) / 2.0_dp
                        baab_RDM_full(i,j) = Temp
                        baab_RDM_full(j,i) = Temp
                    end if

                end if
                
                if ((abs((abab_RDM_full(i,j)*Norm_2RDM)-(abab_RDM_full(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                    Max_Error_Hermiticity = abs((abab_RDM_full(i,j)*Norm_2RDM)-(abab_RDM_full(j,i)*Norm_2RDM))

                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((abab_RDM_full(i,j)*Norm_2RDM)-(abab_RDM_full(j,i)*Norm_2RDM))

                Temp = (abab_RDM_full(i,j) + abab_RDM_full(j,i)) / 2.0_dp
                abab_RDM_full(i,j) = Temp
                abab_RDM_full(j,i) = Temp

                if (tOpenShell)then
                    if ((abs((baba_RDM_full(i,j)*Norm_2RDM)-(baba_RDM_full(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((baba_RDM_full(i,j)*Norm_2RDM)-(baba_RDM_full(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((baba_RDM_full(i,j)*Norm_2RDM)-(baba_RDM_full(j,i)*Norm_2RDM))

                    Temp = (baba_RDM_full(i,j) + baba_RDM_full(j,i)) / 2.0_dp
                    baba_RDM_full(i,j) = Temp
                    baba_RDM_full(j,i) = Temp
                end if

            end do
        end do

    end subroutine make_2e_rdm_hermitian


    subroutine Write_out_2RDM(Norm_2RDM, tNormalise, tmake_herm)

        ! Writes out the 2-RDMs. If tNormalise is true, we print the normalised
        ! (hermitian) matrix. Otherwise we print the unnormalised 2-RDMs, and
        ! we print (in binary) both 2-RDM(Ind1,Ind2) and 2-RDM(Ind2,Ind1)
        ! because this matrix wont be hermitian.

        ! While, for instance, the TwoRDM_aaaa so far has actually been a sum
        ! of the aaaa elements and the bbbb elements.  We only want to print
        ! the aaaa elements.

        real(dp), intent(in) :: Norm_2RDM
        logical, intent(in) :: tNormalise, tmake_herm
        real(dp) :: ParityFactor,Divide_Factor 
        integer :: i, j, a, b, Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab, iunit_4
        integer :: No_Herm_Elements
        integer :: aaaa_RDM_unit, abab_RDM_unit, abba_RDM_unit
        integer :: bbbb_RDM_unit, baba_RDM_unit, baab_RDM_unit
        character(255) :: TwoRDM_aaaa_name, TwoRDM_abab_name, TwoRDM_abba_name
        character(255) :: TwoRDM_bbbb_name, TwoRDM_baba_name, TwoRDM_baab_name
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity, Sum_Herm_Percent 
        real(dp) :: Temp


        if (tNormalise) then
            write(6,*) 'Writing out the *normalised* 2 electron density matrix to file'
            call neci_flush(6)
            ! This takes the TwoRDM_aaaa, and if tWriteMultPops is true (and given that we've put 
            ! .true. in the 3rd position, it'll find the next unused TwoRDM_aaaa.X file name.
            call get_unique_filename('TwoRDM_aaaa',tWriteMultRDMs,.true.,1,TwoRDM_aaaa_name)
            aaaa_RDM_unit = get_free_unit()
            open(aaaa_RDM_unit,file=TwoRDM_aaaa_name,status='unknown')

            call get_unique_filename('TwoRDM_abab',tWriteMultRDMs,.true.,1,TwoRDM_abab_name)
            abab_RDM_unit = get_free_unit()
            open(abab_RDM_unit,file=TwoRDM_abab_name,status='unknown')

            call get_unique_filename('TwoRDM_abba',tWriteMultRDMs,.true.,1,TwoRDM_abba_name)
            abba_RDM_unit = get_free_unit()
            open(abba_RDM_unit,file=TwoRDM_abba_name,status='unknown')

            if (tOpenShell)then
                call get_unique_filename('TwoRDM_bbbb',tWriteMultRDMs,.true.,1,TwoRDM_bbbb_name)
                bbbb_RDM_unit = get_free_unit()
                open(bbbb_RDM_unit,file=TwoRDM_bbbb_name,status='unknown')

                call get_unique_filename('TwoRDM_baba',tWriteMultRDMs,.true.,1,TwoRDM_baba_name)
                baba_RDM_unit = get_free_unit()
                open(baba_RDM_unit,file=TwoRDM_baba_name,status='unknown')

                call get_unique_filename('TwoRDM_baab',tWriteMultRDMs,.true.,1,TwoRDM_baab_name)
                baab_RDM_unit = get_free_unit()
                open(baab_RDM_unit,file=TwoRDM_baab_name,status='unknown')
            end if

        else
            write(6,*) 'Writing out the *unnormalised* 2 electron density matrix to file for reading in'
            call neci_flush(6)
            aaaa_RDM_unit = get_free_unit()
            open(aaaa_RDM_unit,file='TwoRDM_POPS_aaaa',status='unknown',form='unformatted')
            abab_RDM_unit = get_free_unit()
            open(abab_RDM_unit,file='TwoRDM_POPS_abab',status='unknown',form='unformatted')
            abba_RDM_unit = get_free_unit()
            open(abba_RDM_unit,file='TwoRDM_POPS_abba',status='unknown',form='unformatted')

            if (tOpenShell)then
                bbbb_RDM_unit = get_free_unit()
                open(bbbb_RDM_unit,file='TwoRDM_POPS_bbbb',status='unknown',form='unformatted')
                baba_RDM_unit = get_free_unit()
                open(baba_RDM_unit,file='TwoRDM_POPS_baba',status='unknown',form='unformatted')
                baab_RDM_unit = get_free_unit()
                open(baab_RDM_unit,file='TwoRDM_POPS_baab',status='unknown',form='unformatted')
            end if
        end if
       
        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp
        Sum_Herm_Percent = 0.0_dp
        No_Herm_Elements = 0
        do i = 1, SpatOrbs

            do j = i, SpatOrbs

                Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                do a = 1, SpatOrbs

                    do b = a, SpatOrbs

                        Ind2_aa = ( ( (b-2) * (b-1) ) / 2 ) + a
                        Ind2_ab = ( ( (b-1) * b ) / 2 ) + a

                        ! usually each element will have two contributions (from aaaa and bbbb).
                        ! we then need to divide each by 2.
                        ! but in cases where i and j, and a and b, are in the same spatial 
                        ! orbital, there will be only one contribution.
                        if (((i.eq.j).and.(a.eq.b)) .or. tOpenShell) then
                            Divide_Factor = 1.0_dp
                        else
                            Divide_Factor = 2.0_dp
                        end if
                        
                        if ((i.ne.j).and.(a.ne.b)) then

                            if ( (aaaa_RDM_full(Ind1_aa,Ind2_aa).ne.0.0_dp).or.&
                                (aaaa_RDM_full(Ind2_aa,Ind1_aa).ne.0.0_dp) )then
                                ! If we're normalising (and have made the matrix hermitian) we only 
                                ! need to write out Ind1 < Ind2.
                                ! Otherwise we print out Ind1, Ind2 and Ind2, Ind1 so we can 
                                ! find the hermiticity error in the final matrix (after all runs).
                                if (tNormalise.and.(Ind1_aa.le.Ind2_aa)) then
                                    
                                    if ((abs((aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (aaaa_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                        Max_Error_Hermiticity = abs((aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (aaaa_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                            abs((aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (aaaa_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                        (abs((aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                - (aaaa_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                        (abs((aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                + (aaaa_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                    No_Herm_Elements = No_Herm_Elements + 1                                                      

                                    if (tmake_herm) then                                                            
                                        Temp = (aaaa_RDM_full(Ind1_aa,Ind2_aa) + aaaa_RDM_full(Ind2_aa,Ind1_aa)) / 2.0_dp

                                        aaaa_RDM_full(Ind1_aa,Ind2_aa) = Temp
                                        aaaa_RDM_full(Ind2_aa,Ind1_aa) = Temp
                                    end if

                                    if (tFinalRDMEnergy) then
                                        ! For the final calculation, the 2-RDMs will have been made hermitian.
                                        write(aaaa_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( aaaa_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                    else
                                        ! If we're printing the 2-RDMs early (using writeRDMSEVERY), the actual 
                                        ! matrix will not be hermitian, but we want to print a hermitian version.
                                        ! Average the values here.
                                        write(aaaa_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( ((aaaa_RDM_full(Ind1_aa,Ind2_aa) + aaaa_RDM_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                * Norm_2RDM ) / Divide_Factor
                                    end if
                                else if (.not.tNormalise) then
                                    ! for the POPS files, print everything to binary.
                                    ! no divide factor, we just read them in as is.
                                    write(aaaa_RDM_unit) i,j,a,b, &
                                            aaaa_RDM_full(Ind1_aa,Ind2_aa)
                                end if !tNormalise
                            end if !aaaa_RDM_full

                            if (tOpenShell) then

                                if ( (bbbb_RDM_full(Ind1_aa,Ind2_aa).ne.0.0_dp).or.&
                                    (bbbb_RDM_full(Ind2_aa,Ind1_aa).ne.0.0_dp) )then
                                    ! If we're normalising (and have made the matrix hermitian) we only 
                                    ! need to write out Ind1 < Ind2.
                                    ! Otherwise we print out Ind1, Ind2 and Ind2, Ind1 so we can 
                                    ! find the hermiticity error in the final matrix (after all runs).
                                    if (tNormalise.and.(Ind1_aa.le.Ind2_aa)) then

                                        if ((abs((bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (bbbb_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                            Max_Error_Hermiticity = abs((bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (bbbb_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                                abs((bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (bbbb_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                            (abs((bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (bbbb_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                            (abs((bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    + (bbbb_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1

                                        if (tmake_herm) then
                                            Temp = (bbbb_RDM_full(Ind1_aa,Ind2_aa) + bbbb_RDM_full(Ind2_aa,Ind1_aa)) / 2.0_dp

                                            bbbb_RDM_full(Ind1_aa,Ind2_aa) = Temp
                                            bbbb_RDM_full(Ind2_aa,Ind1_aa) = Temp
                                end if

                                        if (tFinalRDMEnergy) then
                                            ! For the final calculation, the 2-RDMs will have been made hermitian.
                                            write(bbbb_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                    ( bbbb_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            ! If we're printing the 2-RDMs early (using writeRDMSEVERY), the actual 
                                            ! matrix will not be hermitian, but we want to print a hermitian version.
                                            ! Average the values here.
                                            write(bbbb_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( ((bbbb_RDM_full(Ind1_aa,Ind2_aa) + bbbb_RDM_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                            end if
                                    else if (.not.tNormalise) then
                                        ! for the POPS files, print everything to binary.
                                        ! no divide factor, we just read them in as is.
                                        write(bbbb_RDM_unit) i,j,a,b, &
                                                bbbb_RDM_full(Ind1_aa,Ind2_aa)
                                    end if  ! tNormalise
                                end if  ! bbbb_RDM_full 
                            end if ! OpenShell

                            if (tOpenShell)then

                                if ( (abba_RDM_full(Ind1_aa,Ind2_aa) .ne. 0.0_dp) .or. &
                                     (baab_RDM_full(Ind2_aa,Ind1_aa) .ne. 0.0_dp) ) then

                                    if (tNormalise.and.(Ind1_aa .le. Ind2_aa)) then

                                        if ((abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (baab_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                            Max_Error_Hermiticity = abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (baab_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                                abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (baab_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                            (abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            - (baab_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                            (abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            + (baab_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1                                                  


                                        if (tmake_herm) then                                                            
                                            Temp = (abba_RDM_full(Ind1_aa,Ind2_aa) + baab_RDM_full(Ind2_aa,Ind1_aa)) / 2.0_dp
                                            abba_RDM_full(Ind1_aa,Ind2_aa) = Temp
                                            baab_RDM_full(Ind2_aa,Ind1_aa) = Temp
                                        end if

                                        if (tFinalRDMEnergy) then
                                            write(abba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( abba_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            write(abba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( ((abba_RDM_full(Ind1_aa,Ind2_aa) + abba_RDM_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                        * Norm_2RDM ) / Divide_Factor
                                        end if
                                    else if (.not.tNormalise) then
                                        write(abba_RDM_unit) i,j,a,b, &
                                            abba_RDM_full(Ind1_aa,Ind2_aa) 
                                    end if
                                end if ! abba_RDM_full baab_RDM_full

                                if ( (baab_RDM_full(Ind1_aa,Ind2_aa).ne.0.0_dp) .or. &
                                     (abba_RDM_full(Ind2_aa,Ind1_aa).ne.0.0_dp) ) then

                                if (tNormalise.and.(Ind1_aa.le.Ind2_aa)) then

                                        if ((abs((baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                            Max_Error_Hermiticity = abs((baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                                abs((baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                            (abs((baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                            (abs((baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            + (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1

                                        if (tmake_herm) then
                                            Temp = (baab_RDM_full(Ind1_aa,Ind2_aa) + abba_RDM_full(Ind2_aa,Ind1_aa)) / 2.0_dp
                                            baab_RDM_full(Ind1_aa,Ind2_aa) = Temp
                                            abba_RDM_full(Ind2_aa,Ind1_aa) = Temp  
                                        end if ! tmake_herm = .true.

                                        if (tFinalRDMEnergy) then
                                            write(baab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( baab_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            write(baab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( ((baab_RDM_full(Ind1_aa,Ind2_aa) + baab_RDM_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                        * Norm_2RDM ) / Divide_Factor
                                        end if  ! tFinalRDMEnergy = .true.

                                    else if (.not.tNormalise) then
                                        write(baab_RDM_unit) i,j,a,b, &
                                            baab_RDM_full(Ind1_aa,Ind2_aa) 
                                    end if ! tNormalise
                                end if ! baab_RDM_full abba_RDM_full

                            else ! not tOpenShell

                                if ( (abba_RDM_full(Ind1_aa,Ind2_aa).ne.0.0_dp).or.&
                                    (abba_RDM_full(Ind2_aa,Ind1_aa).ne.0.0_dp) ) then

                                    if (tNormalise.and.(Ind1_aa.le.Ind2_aa)) then
                                        if ((abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                        Max_Error_Hermiticity = abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                            abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                        (abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                        - (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                        (abs((abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                        + (abba_RDM_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1

                                        if (tmake_herm) then                                                            
                                            Temp = (abba_RDM_full(Ind1_aa,Ind2_aa) + abba_RDM_full(Ind2_aa,Ind1_aa)) / 2.0_dp
                                            abba_RDM_full(Ind1_aa,Ind2_aa) = Temp
                                            abba_RDM_full(Ind2_aa,Ind1_aa) = Temp
                                        end if

                                        if (tFinalRDMEnergy) then
                                            write(abba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( abba_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            write(abba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( ((abba_RDM_full(Ind1_aa,Ind2_aa) + abba_RDM_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                        end if
                                    else if (.not.tNormalise) then
                                        write(abba_RDM_unit) i,j,a,b, &
                                            abba_RDM_full(Ind1_aa,Ind2_aa) 
                                    end if
                                end if
                            end if 

                        end if  ! (i.ne.j) .and. (a.ne.b) 

                        if ( (abab_RDM_full(Ind1_ab,Ind2_ab).ne.0.0_dp).or.&
                            (abab_RDM_full(Ind2_ab,Ind1_ab).ne.0.0_dp) ) then

                            if ((abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                        - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                Max_Error_Hermiticity = abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                        - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                            Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                    abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                            Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                (abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                    - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                (abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                   + (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                            No_Herm_Elements = No_Herm_Elements + 1                                                        

                            if (tmake_herm) then                                                            
                                Temp = (abab_RDM_full(Ind1_ab,Ind2_ab) + abab_RDM_full(Ind2_ab,Ind1_ab)) / 2.0_dp
                                abab_RDM_full(Ind1_ab,Ind2_ab) = Temp
                                abab_RDM_full(Ind2_ab,Ind1_ab) = Temp
                            end if

                            if (tNormalise.and.(Ind1_ab.le.Ind2_ab)) then
                                if (tFinalRDMEnergy) then
                                    write(abab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                        ( abab_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                else
                                    write(abab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                        ( ((abab_RDM_full(Ind1_ab,Ind2_ab) + abab_RDM_full(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                * Norm_2RDM ) / Divide_Factor
                                end if
                            else if (.not.tNormalise) then
                                write(abab_RDM_unit) i,j,a,b, &
                                    abab_RDM_full(Ind1_ab,Ind2_ab) 
                            end if
                        end if

                        if (tOpenShell .and. ( (a .ne. b) .or. (i .ne. j) ))then

                            if ( (baba_RDM_full(Ind1_ab,Ind2_ab).ne.0.0_dp).or.&
                                (baba_RDM_full(Ind2_ab,Ind1_ab).ne.0.0_dp) ) then

                                if ((abs((baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (baba_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                    Max_Error_Hermiticity = abs((baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (baba_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                        abs((baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                            - (baba_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                    (abs((baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (baba_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                    (abs((baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        + (baba_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                                No_Herm_Elements = No_Herm_Elements + 1                                                        

                                if (tmake_herm) then                                                            
                                    Temp = (baba_RDM_full(Ind1_ab,Ind2_ab) + baba_RDM_full(Ind2_ab,Ind1_ab)) / 2.0_dp
                                    baba_RDM_full(Ind1_ab,Ind2_ab) = Temp
                                    baba_RDM_full(Ind2_ab,Ind1_ab) = Temp
                                end if

                                if (tNormalise.and.(Ind1_ab.le.Ind2_ab)) then
                                    if (tFinalRDMEnergy) then
                                        write(baba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( baba_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                    else
                                        write(baba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( ((baba_RDM_full(Ind1_ab,Ind2_ab) + baba_RDM_full(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                    end if
                                else if (.not.tNormalise) then
                                    write(baba_RDM_unit) i,j,a,b, &
                                        baba_RDM_full(Ind1_ab,Ind2_ab) 
                                end if
                            end if

                         else if (tOpenShell) then !a=b & i=j -> baba term saved in abab

                            if ( (abab_RDM_full(Ind1_ab,Ind2_ab).ne.0.0_dp).or.&
                                (abab_RDM_full(Ind2_ab,Ind1_ab).ne.0.0_dp) ) then

                                if ((abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                    Max_Error_Hermiticity = abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                        abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                            - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                    (abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                    (abs((abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        + (abab_RDM_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                                No_Herm_Elements = No_Herm_Elements + 1                                                        

                                if (tmake_herm) then                                                             
                                    Temp = (abab_RDM_full(Ind1_ab,Ind2_ab) + abab_RDM_full(Ind2_ab,Ind1_ab)) / 2.0_dp
                                    abab_RDM_full(Ind1_ab,Ind2_ab) = Temp
                                    abab_RDM_full(Ind2_ab,Ind1_ab) = Temp
                                end if

                                if (tNormalise.and.(Ind1_ab.le.Ind2_ab)) then
                                    if (tFinalRDMEnergy) then
                                        write(baba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( abab_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                    else
                                        write(baba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( ((abab_RDM_full(Ind1_ab,Ind2_ab) + abab_RDM_full(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                    end if
                                end if
                            end if

                         end if ! tOpenShell


                    end do  
                end do  

            end do  
        end do 

        close(aaaa_RDM_unit)
        close(abab_RDM_unit)
        close(abba_RDM_unit)
        if (tOpenShell)then
            close(bbbb_RDM_unit)
            close(baba_RDM_unit)
            close(baab_RDM_unit)
        end if

        if (tWriteSpinFreeRDM) call Write_spinfree_RDM(Norm_2RDM)

        if (tNormalise) then
            write(6,'(I15,F30.20,A20,A39)') Iter+PreviousCycles, Max_Error_Hermiticity, &
                                            '( Iteration,',' MAX ABS ERROR IN HERMITICITY )'
            write(6,'(I15,F30.20,A20,A39)') Iter+PreviousCycles, Sum_Error_Hermiticity, &
                                            '( Iteration,',' SUM ABS ERROR IN HERMITICITY )'
            write(6,'(I15,F30.20,A20,A51)') Iter+PreviousCycles, Sum_Herm_Percent/real(No_Herm_Elements,dp), &
                                            '( Iteration,',' AVERAGE ABS PERCENTAGE HERMITICITY ERROR )'
        end if

    end subroutine Write_out_2RDM

    subroutine Write_spinfree_RDM (Norm_2RDM)

        ! Write out the spinfree 2RDM for MPQC.

        real(dp), intent(in) :: Norm_2RDM
        integer :: i, j, a, b
        integer :: read_stat
        integer :: spinfree_RDM_unit

        write(6,*) "Writing out the spinfree RDM"
        spinfree_RDM_unit = get_free_unit()
        open(spinfree_RDM_unit, file="spinfree_TwoRDM", status="replace")
        
        do j = 1, SpatOrbs
            do b = 1, SpatOrbs
                do a = 1, SpatOrbs
                    do i = 1, SpatOrbs
                        if (abs(Find_Spatial_2RDM_Chem(i,j,a,b,Norm_2RDM)) > 1.0e-12) then

                            write(spinfree_RDM_unit,"(4I15,F30.20)") i, a, j, b, &
                                   Find_Spatial_2RDM_Chem(i,j,a,b,Norm_2RDM)
                        end if
                    end do
                end do
            end do
        end do

        write(spinfree_RDM_unit, "(4I15,F30.20)") -1, -1, -1, -1, -1.0_dp
        close(spinfree_RDM_unit)

    end subroutine Write_spinfree_RDM

    subroutine Calc_Energy_from_RDM(Norm_2RDM)

        ! This routine takes the 1 electron and 2 electron reduced density matrices 
        ! and calculated the energy they give.    
        ! The equation for the energy is as follows:
        !
        !   E = Tr(h1 1RDM) + 1/2 Tr(h2 2RDM)
        !
        ! where h1 are the 2 index integrals, and h2 the 4 index integrals.  The traces, Tr, 
        ! are given by:
        !   Tr(h1 1RDM) = Sum_i,j [ h1(i,j) 1RDM(j,i) ]
        !   Tr(h2 2RDM) = Sum_i,j;k,l [ h2(i,j;k,l) 2RDM(k,l;i,j) ]

        real(dp), intent(out) :: Norm_2RDM
        real(dp) :: Norm_2RDM_Inst
        integer :: i,j,a,b,Ind1_aa,Ind1_ab,Ind2_aa,Ind2_ab,ierr
        integer :: iSpin, jSpin, error
        real(dp) :: RDMEnergy_Inst, RDMEnergy, Coul, Exch, Parity_Factor
        real(dp) :: Coul_aaaa, Coul_bbbb, Coul_abba, Coul_baab, Coul_abab, Coul_baba
        real(dp) :: Exch_aaaa, Exch_bbbb, Exch_abba, Exch_baab, Exch_abab, Exch_baba
        real(dp) :: Trace_2RDM_New, RDMEnergy1, RDMEnergy2, spin_est

        call set_timer(RDMEnergy_Time,30)

        Trace_2RDM_New = 0.0_dp

        RDMEnergy_Inst = 0.0_dp
        RDMEnergy1 = 0.0_dp
        RDMEnergy2 = 0.0_dp
        RDMEnergy = 0.0_dp
    
        ! Normalise, make hermitian, print etc.
        call Finalise_2e_RDM(Norm_2RDM_Inst, Norm_2RDM)

        if (tFinalRDMEnergy) then
            write(6,*) ''
            write(6,*) 'Calculating the final RDM energy'
        end if

        if (iProcIndex.eq.0) then

            do i = 1, SpatOrbs
                iSpin = 2 * i

                do j = i, SpatOrbs
                    jSpin = 2 * j

                    Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                    Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                    do a = 1, SpatOrbs

                        ! Adding in contributions effectively from the 1-RDM (although these are calculated 
                        ! from the 2-RDM.
                        call calc_1RDM_energy(i,j,a,iSpin,jSpin, Norm_2RDM, Norm_2RDM_Inst, &
                                                    RDMEnergy_Inst, RDMEnergy1,RDmEnergy2)
                        
                        do b = a, SpatOrbs

                            Ind2_aa = ( ( (b-2) * (b-1) ) / 2 ) + a
                            Ind2_ab = ( ( (b-1) * b ) / 2 ) + a

                            ! UMAT in chemical notation.
                            ! In spin or spatial orbitals.
                            Coul = real(UMAT(UMatInd(i,j,a,b,0,0)),dp)
                            Exch = real(UMAT(UMatInd(i,j,b,a,0,0)),dp)
                            
                            if ((i.ne.j).and.(a.ne.b)) then
                                ! Cannot get i=j or a=b contributions in aaaa.

                                if (tStoreSpinOrbs)then
                                    Coul_aaaa=real(UMAT(UMatInd(2*i, 2*j, 2*a, 2*b, 0, 0)),dp)
                                    Coul_bbbb=real(UMAT(UMatInd(2*i-1, 2*j-1, 2*a-1, 2*b-1, 0, 0)),dp)
                                    Exch_aaaa=real(UMAT(UMatInd(2*i, 2*j, 2*b, 2*a, 0, 0)),dp)
                                    Exch_bbbb=real(UMAT(UMatInd(2*i-1, 2*j-1, 2*b-1, 2*a-1, 0, 0)),dp)     

                                    if (tRDMInstEnergy)then 
                                        RDMEnergy_Inst = RDMEnergy_Inst + (aaaa_RDM(Ind1_aa,Ind2_aa) &
                                                          *( Coul_aaaa - Exch_aaaa ) )
                                        RDMEnergy_Inst = RDMEnergy_Inst + (bbbb_RDM(Ind1_aa,Ind2_aa) &
                                                          *  ( Coul_bbbb - Exch_bbbb ) )
                                    end if

                                    RDMEnergy2 = RDMEnergy2 + ( aaaa_RDM_full(Ind1_aa,Ind2_aa) &
                                                          * Norm_2RDM * ( Coul_aaaa - Exch_aaaa ) )
                                    RDMEnergy2 = RDMEnergy2 + ( bbbb_RDM_full(Ind1_aa,Ind2_aa) &
                                                          * Norm_2RDM * ( Coul_bbbb - Exch_bbbb ) )
 
                                else 

                                    if (tRDMInstEnergy)then
                                        RDMEnergy_Inst = RDMEnergy_Inst + (aaaa_RDM(Ind1_aa,Ind2_aa) &
                                                                    *( Coul - Exch ) )
                                    end if

                                    RDMEnergy2 = RDMEnergy2 + ( aaaa_RDM_full(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * ( Coul - Exch ) ) 

                                    if (tOpenShell)then
                                        if (tRDMInstEnergy)then 
                                            RDMEnergy_Inst = RDMEnergy_Inst + (bbbb_RDM(Ind1_aa,Ind2_aa) &
                                                          *( Coul - Exch ) )
                                        end if

                                        RDMEnergy2 = RDMEnergy2 + ( bbbb_RDM_full(Ind1_aa,Ind2_aa) &
                                                          * Norm_2RDM * ( Coul - Exch ) )
                                    end if 

                                end if  


                                if (Ind1_aa.eq.Ind2_aa)then
                                    Trace_2RDM_New = Trace_2RDM_New + &
                                                        aaaa_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM
                                    if (tOpenShell) Trace_2RDM_New = Trace_2RDM_New + &
                                                        bbbb_RDM_full(Ind1_aa,Ind2_aa) * Norm_2RDM
                                end if
 
                                ! For abab cases, coul element will be non-zero, exchange zero.

                                if (tStoreSpinOrbs)then
                                    Coul_abab=real(UMAT(UMatInd(2*i, 2*j-1, 2*a, 2*b-1, 0, 0)),dp)
                                    Coul_baba=real(UMAT(UMatInd(2*i-1, 2*j, 2*a-1, 2*b, 0, 0)),dp)

                                    call neci_flush(6)

                                    if (tRDMInstEnergy)then
                                        RDMEnergy_Inst = RDMEnergy_Inst + ( abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        *  Coul_abab )
                                    end if
                                    RDMEnergy2 = RDMEnergy2 + ( abab_RDM_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Coul_abab ) 

                                    if (tRDMInstEnergy)then
                                        RDMEnergy_Inst = RDMEnergy_Inst + ( baba_RDM(Ind1_ab,Ind2_ab) &
                                                                        *  Coul_baba )
                                    end if
                                    RDMEnergy2 = RDMEnergy2 + ( baba_RDM_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Coul_baba ) 

                                else 
                                    if (tRDMInstEnergy)then
                                        RDMEnergy_Inst = RDMEnergy_Inst + ( abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        *  Coul )
                                        if (tOpenShell) RDMEnergy_Inst = RDMEnergy_Inst + ( baba_RDM(Ind1_ab,Ind2_ab) &
                                                                        *  Coul )
                                    end if

                                    RDMEnergy2 = RDMEnergy2 + ( abab_RDM_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Coul ) 
                                    if (tOpenShell) RDMEnergy2 = RDMEnergy2 + ( baba_RDM_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Coul) 

                                end if 

                                if (Ind1_ab.eq.Ind2_ab)then
                                    Trace_2RDM_New = Trace_2RDM_New + &
                                                        abab_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                    if (tOpenShell) Trace_2RDM_New = Trace_2RDM_New + &
                                                        baba_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                end if

                                ! For abba cases, coul element will be zero, exchange non-zero.

                                if (tStoreSpinOrbs)then
                                    Exch_abba=real(UMAT(UMatInd(2*i, 2*j-1, 2*b, 2*a-1, 0, 0)),dp)
                                    Exch_baab=real(UMAT(UMatInd(2*i-1, 2*j, 2*b-1, 2*a, 0, 0)),dp)

                                    if (tRDMInstEnergy)then 
                                        RDMEnergy_Inst = RDMEnergy_Inst - ( abba_RDM(Ind1_aa,Ind2_aa) &
                                                                        * Exch_abba )
                                        RDMEnergy_Inst = RDMEnergy_Inst - ( baab_RDM(Ind1_aa,Ind2_aa) &
                                                                        * Exch_baab )
                                    end if

                                    RDMEnergy2 = RDMEnergy2 - ( abba_RDM_full(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * Exch_abba ) 
                                    RDMEnergy2 = RDMEnergy2 - ( baab_RDM_full(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * Exch_baab ) 

                                else 

                                    if (tRDMInstEnergy)then 
                                        RDMEnergy_Inst = RDMEnergy_Inst - ( abba_RDM(Ind1_aa,Ind2_aa) &
                                                                        * Exch )
                                        if (tOpenShell) RDMEnergy_Inst = RDMEnergy_Inst - ( baab_RDM(Ind1_aa,Ind2_aa) &
                                                                        *  Exch )
                                    end if

                                    RDMEnergy2 = RDMEnergy2 - ( abba_RDM_full(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * Exch ) 
                                    if (tOpenShell) RDMEnergy2 = RDMEnergy2 - ( baab_RDM_full(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * Exch ) 

                               end if 

                           else if ( (i .eq. j) .or. (a .eq. b) )then
                                ! i = j or a = b
                                ! abab has both abab and abba elements in them effectively.
                                ! half will have non-zero coul, and half non-zero exchange.
                                ! For abab/baba Exch = 0, and for abba/baab Coul=0
                                ! abba/baab saved in abab/baba. Changes the sign. 
                                if (tStoreSpinOrbs)then
                                    Coul_abab=real(UMAT(UMatInd(2*i, 2*j-1, 2*a, 2*b-1, 0, 0)),dp)
                                    Coul_baba=real(UMAT(UMatInd(2*i-1, 2*j, 2*a-1, 2*b, 0, 0)),dp)

                                    if ( (i .eq. j) .and. (a .eq. b) ) then
                                        ! This term is saved in abab only
                                        Exch_abba=real(UMAT(UMatInd(2*i, 2*j-1, 2*b, 2*a-1, 0, 0)),dp)
    
                                        if (tRDMInstEnergy)then
                                            RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * abab_RDM(Ind1_ab,Ind2_ab) &
                                                                              * (Coul_abab+Exch_abba)
                                   
                                        end if 
    
                                        RDMEnergy2 = RDMEnergy2 + 0.5_dp * abab_RDM_full(Ind1_ab,Ind2_ab) &
                                                                  * Norm_2RDM * (Coul_abab+Exch_abba)

                                    else if (i .eq. j) then
                                        ! i = j : Swap first indeces to get abba/baab terms
                                        ! abba saved in baba, baab saved in abab (sign changes)
                                        Exch_abba=real(UMAT(UMatInd(2*j, 2*i-1, 2*b, 2*a-1, 0, 0)),dp) 
                                        Exch_baab=real(UMAT(UMatInd(2*j-1, 2*i, 2*b-1, 2*a, 0, 0)),dp) 
    
                                        if (tRDMInstEnergy)then
                                            RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * abab_RDM(Ind1_ab,Ind2_ab) &
                                                                              *  (Coul_abab+Exch_baab)
                                            RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * baba_RDM(Ind1_ab,Ind2_ab) &
                                                                              *  (Coul_baba+Exch_abba)
                                        end if 
    
                                        RDMEnergy2 = RDMEnergy2 +  0.5_dp * abab_RDM_full(Ind1_ab,Ind2_ab) &
                                                                   * Norm_2RDM * (Coul_abab+Exch_baab)


                                        RDMEnergy2 = RDMEnergy2 +  0.5_dp * baba_RDM_full(Ind1_ab,Ind2_ab) &
                                                                   * Norm_2RDM * (Coul_baba+Exch_abba)


                                    else if (a .eq. b) then
                                        ! a = b : Swap last indeces to get abba/baab terms
                                        ! abba saved in abab, baab saved in baba (sign changes)
                                        Exch_abba=real(UMAT(UMatInd(2*i, 2*j-1, 2*a, 2*b-1, 0, 0)),dp)
                                        Exch_baab=real(UMAT(UMatInd(2*i-1, 2*j, 2*a-1, 2*b, 0, 0)),dp)
    
                                        if (tRDMInstEnergy)then
                                            RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        *(Coul_abab+Exch_abba)
                                            RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * baba_RDM(Ind1_ab,Ind2_ab) &
                                                                        * (Coul_baba+Exch_baab)
                                        end if 
    
                                        RDMEnergy2 = RDMEnergy2 + 0.5_dp * abab_RDM_full(Ind1_ab,Ind2_ab) &
                                                                * Norm_2RDM * (Coul_abab+Exch_abba)


                                        RDMEnergy2 = RDMEnergy2 + 0.5_dp * baba_RDM_full(Ind1_ab,Ind2_ab) &
                                                                * Norm_2RDM * (Coul_baba+Exch_baab)
                                    end if

                                else

                                    if (tRDMInstEnergy)then 
                                        RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        * (Coul+Exch)
                                        if (tOpenShell) RDMEnergy_Inst = RDMEnergy_Inst  &
                                                                        + 0.5_dp *baba_RDM(Ind1_ab,Ind2_ab) &
                                                                        *  (Coul+Exch)
                                    end if 

                                    RDMEnergy2 = RDMEnergy2 + 0.5_dp * abab_RDM_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * (Coul+Exch)


                                    if (tOpenShell) RDMEnergy2 = RDMEnergy2 &
                                                            + 0.5_dp * baba_RDM_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * (Coul+Exch)

                                end if 

                                if (Ind1_ab.eq.Ind2_ab) then
                                    Trace_2RDM_New = Trace_2RDM_New + &
                                                        abab_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                    if (tOpenShell) then
                                        Trace_2RDM_New = Trace_2RDM_New + &
                                                        baba_RDM_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                    end if  
                                end if

                            end if

                       end do

                    end do
                end do
            end do

            if (tRDMInstEnergy) RDMEnergy_Inst = RDMEnergy_Inst + Ecore/Norm_2RDM_Inst
            RDMEnergy = RDMEnergy1 + RDMEnergy2 + Ecore 

            ! Calculate the instantaneous estimate of <S^2> using the 2RDM.
            spin_est = calc_2rdm_spin_estimate(Norm_2RDM_Inst)

            ! Obviously this 'instantaneous' energy is actually accumulated between energy 
            ! print outs.
            if (tRDMInstEnergy) then
                write(Energies_unit, "(1X,I13,3(2X,es20.13))") Iter+PreviousCycles, RDMEnergy_Inst, spin_est, &
                                                                1.0_dp/Norm_2RDM_Inst
            else
                write(Energies_unit, "(I31,F30.15)") Iter+PreviousCycles, RDMEnergy
            end if
            call neci_flush(Energies_unit)

            if (tFinalRDMEnergy) then
                write(6,*) 'Trace of 2-el-RDM before normalisation : ',Trace_2RDM
                write(6,*) 'Trace of 2-el-RDM after normalisation : ',Trace_2RDM_New
                write(6,*) 'Energy contribution from the 1-RDM: ',RDMEnergy1
                write(6,*) 'Energy contribution from the 2-RDM: ',RDMEnergy2
                write(6,'(A64,F30.20)') ' *TOTAL ENERGY* CALCULATED USING THE *REDUCED &
                                            &DENSITY MATRICES*:',RDMEnergy
                close(Energies_unit) 
            end if

        end if

        ! Zero all the 'instantaneous' stuff.
        if (tRDMInstEnergy) then
            aaaa_RDM(:,:) = 0.0_dp
            abab_RDM(:,:) = 0.0_dp
            abba_RDM(:,:) = 0.0_dp

            if (tOpenShell)then
                bbbb_RDM(:,:) = 0.0_dp
                baba_RDM(:,:) = 0.0_dp
                baab_RDM(:,:) = 0.0_dp
            end if

            AccumRDMNorm_Inst = 0.0_dp
            Trace_2RDM_Inst = 0.0_dp
        end if

        call halt_timer(RDMEnergy_Time)

    end subroutine Calc_Energy_from_RDM
    
    function calc_2rdm_spin_estimate(Norm_2RDM_Inst) result(spin_est)

        ! Return the (unnormalised) estimate of <S^2> from the instantaneous
        ! 2RDM estimates.

        real(dp), intent(in) :: Norm_2RDM_Inst
        integer :: i, j
        integer :: Ind1_aa, Ind1_ab
        real(dp) :: spin_est

        spin_est = 0.0_dp

        do i = 1, SpatOrbs

            do j = i+1, SpatOrbs

                Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                spin_est = spin_est + 2*aaaa_RDM(Ind1_aa, Ind1_aa) &
                                    - 2*abab_RDM(Ind1_ab, Ind1_ab) &
                                    + 4*abba_RDM(Ind1_aa, Ind1_aa)

                if (tOpenShell) then
                    spin_est = spin_est + 2*bbbb_RDM(Ind1_aa, Ind1_aa) &
                                        - 2*baba_RDM(Ind1_ab, Ind1_ab) &
                                        + 4*baab_RDM(Ind1_aa, Ind1_aa)
                end if

            end do

            ! i = j term.
            Ind1_ab = ( ( (i-1) * i ) / 2 ) + i

            spin_est = spin_est - 6*abab_RDM(Ind1_ab, Ind1_ab)
            if (tOpenShell) spin_est = spin_est - 6*baba_RDM(Ind1_ab, Ind1_ab)

        end do 

        spin_est = spin_est + 3*real(nel,dp)/Norm_2RDM_Inst

        spin_est = spin_est/4.0_dp

    end function calc_2rdm_spin_estimate

    subroutine Calc_Lagrangian_from_RDM(Norm_1RDM,Norm_2RDM)

        ! This routine takes the 1 electron and 2 electron reduced density
        ! matrices and calculated the Lagrangian term, X, required for the
        ! calculation of forces.

        ! The equation for X is as follows:
        !
        !   X_pq = Sum_r[h_pr 1RDM_qr] + 0.5*Sum_rst[(pr|st)[2RDM_qrst + 2RDM_rqst]]
        !
        !   where 2RDM is defined in chemical notation sense:
        !
        !   2RDM_ijkl = <Psi| a_i+ a_k+ a_l a_j|Psi>

        use IntegralsData, only: UMAT
        use UMatCache, only: UMatInd
        use RotateOrbsMod, only: SymLabelList2_rot
        use UMatCache, only: GTID
        use Logging, only: tDumpForcesInfo, tPrintLagrangian

        real(dp), intent(in) :: Norm_2RDM
        real(dp), intent(in) :: Norm_1RDM
        real(dp) :: Norm_2RDM_Inst
        integer :: p,q,r,s,t,ierr,stat
        integer :: pSpin, qSpin, rSpin, error
        real(dp) :: RDMEnergy_Inst, RDMEnergy, Coul, Exch, Parity_Factor 
        real(dp) :: Trace_2RDM_New, RDMEnergy1, RDMEnergy2
        real(dp) :: qrst, rqst
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity, Temp

        allocate(Lagrangian(SpatOrbs,SpatOrbs),stat=ierr)
        Lagrangian(:,:)=0.0_dp
        
        ! We will begin by calculating the Lagrangian in chemical notation - we will explicitely calculate
        ! both halves (X_pq and X_qp) in order to check if it is symmetric or not.
        !  - a symmetric Lagrangian is important in allowing us to use the derivative overlap matrix rather than the 
        !    coupled-perturbed coefficients when calculating the Forces later on 
        !    (see Sherrill, Analytic Gradients of CI Energies eq38)

        write(6,*) ''
        write(6,*) 'Calculating the Lagrangian X from the final density matrices'

        !Calculating the Lagrangian X in terms of spatial orbitals
        if (iProcIndex.eq.0) then

            do p = 1, SpatOrbs  !Run over spatial orbitals
                pSpin = 2*p    ! Picks out beta component
                do q = 1, SpatOrbs  !Want to calculate X(p,q) separately from X(q,p) for now to see if we're symmetric
                    do r = 1, SpatOrbs
                        rSpin=2*r
                        ! Adding in contributions effectively from the 1-RDM
                        ! We made sure earlier that the 1RDM is contructed, so we can call directly from this
                        if (tOpenShell) then
                            ! Include both aa and bb contributions 
                            Lagrangian(p,q)=Lagrangian(p,q)+(NatOrbMat(SymLabelListInv_rot(2*q),SymLabelListInv_rot(2*r)))*&
                                                                real(TMAT2D(pSpin,rSpin),8)*Norm_1RDM
                            Lagrangian(p,q)=Lagrangian(p,q)+(NatOrbMat(SymLabelListInv_rot(2*q-1),SymLabelListInv_rot(2*r-1)))*&
                                                                real(TMAT2D(pSpin-1,rSpin-1),8)*Norm_1RDM
                        else
                            !We will be here most often (?)
                            Lagrangian(p,q)=Lagrangian(p,q)+NatOrbMat(SymLabelListInv_rot(q),SymLabelListInv_rot(r))* &
                                                                real(TMAT2D(pSpin,rSpin),8)*Norm_1RDM
                        end if

                        do s = 1, SpatOrbs
                
                            ! Here we're looking to start adding on the contributions from the 2-RDM and the 2-el integrals
                            ! For X(p,q), these have the form 0.5*Sum_rst[(pr|st)[2RDM_qrst + 2RDM_rqst]]
                            
                            ! NOTE: In some notations, the factor of a half goes *inside* the 2RDM - 
                            ! consistent with Yamaguchi etc (see eq 11 in Sherrill, Analytic Gradients
                            ! of CI energies).  We will keep it outside for now, consistent with the storage
                            ! within neci

                            do t=1, SpatOrbs
                                
                                !Integral (pr|st) = <ps|rt>
                                !Give indices in PHYSICAL NOTATION
                                !NB, FCIDUMP is labelled in chemical notation
                                Coul = real(UMAT(UMatInd(p,s,r,t,0,0)),8)

                                qrst=Find_Spatial_2RDM_Chem(q,r,s,t, Norm_2RDM)
                                rqst=Find_Spatial_2RDM_Chem(r,q,s,t, Norm_2RDM)
                                
                                Lagrangian(p,q)=Lagrangian(p,q) + 0.5_dp*Coul*(qrst+rqst)
                            end do
                        end do
                    end do
                end do
            end do
       
            !! Now symmetrise (make hermitian, such that X_pq = X_qp) the Lagrangian X

            Max_Error_Hermiticity = 0.0_dp
            Sum_Error_Hermiticity = 0.0_dp
            do p = 1, SpatOrbs
                do q = p, SpatOrbs
                    if (abs(Lagrangian(p,q) - Lagrangian(q,p)).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs(Lagrangian(p,q)-Lagrangian(q,p))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity+abs(Lagrangian(p,q) - Lagrangian(q,p))

                    Temp = (Lagrangian(p,q) + Lagrangian(q,p))/2.0_dp
                    Lagrangian(p,q) = Temp
                    Lagrangian(q,p) = Temp
                end do
            end do

            ! Output the hermiticity errors.
            write(6,'(A40,F30.20)') ' MAX ABS ERROR IN Lagrangian HERMITICITY', Max_Error_Hermiticity
            write(6,'(A40,F30.20)') ' SUM ABS ERROR IN Lagrangian HERMITICITY', Sum_Error_Hermiticity

        end if

    end subroutine Calc_Lagrangian_from_RDM

    subroutine calc_1RDM_energy(i,j,a,iSpin,jSpin,Norm_2RDM,Norm_2RDM_Inst,&
                                                    RDMEnergy_Inst,RDMEnergy1,RDMEnergy2)

        ! This routine calculates the 1-RDM part of the RDM energy, and
        ! constructs the 1-RDM if required for diagonalisation or something.
        ! gamma(i,j) = [1/(NEl - 1)] * SUM_a Gamma(i,a,j,a) 
        ! want to calculate:    gamma(i,j) * h_ij
        ! h_ij => TMAT2D(iSpin,jSpin)
        ! iSpin = 2*i, jSpin = 2*j  -> alpha orbs

        use OneEInts, only: TMAT2D
        use Logging, only: tDiagRDM, tDumpForcesInfo, tDipoles

        integer, intent(in) :: i,j,a,iSpin,jSpin
        real(dp), intent(in) :: Norm_2RDM, Norm_2RDM_Inst
        real(dp), intent(inout) :: RDMEnergy_Inst, RDMEnergy1,RDmEnergy2
        real(dp) :: Parity_Factor, fac_doublecount
        integer :: Ind1_1e_ab, Ind2_1e_ab
        integer :: Ind1_1e_aa, Ind2_1e_aa
        integer :: iSpin_abab, iSpin_baba
        integer :: jSpin_abab, jSpin_baba
        integer :: iSpin_abba, iSpin_baab
        integer :: jSpin_abba, jSpin_baab
        logical :: t_abab_only, t_opposite_contri

        ! for i a -> j a excitation, when lined up as min max -> min max, 
        ! if a's are aligned, only a b a b arrays contain single excitations, 
        ! if a's not aligned, a b b a.
        ! all a a a a will contain single excitations.

        ! abab & baba terms
        if (((i.le.a).and.(j.le.a)).or.((i.ge.a).and.(j.ge.a))) then

            Ind1_1e_ab = ( ( (max(i,a)-1) * max(i,a) ) / 2 ) + min(i,a)
            Ind2_1e_ab = ( ( (max(j,a)-1) * max(j,a) ) / 2 ) + min(j,a)
            if (Ind1_1e_ab.ne.Ind2_1e_ab) then
            ! For Gamma elements corresponding to 1-RDMs ( Gamma(i,a,j,a) ), 
            ! we're only considering i =< j 
            ! therefore we need to sum in the opposite contribution too if i ne j.
                t_opposite_contri = .true.
            else  ! no opposite contribution from i = j term
                t_opposite_contri = .false.
            end if

            fac_doublecount = 1.0_dp
            if ( (i .lt. a) .or. (j .lt. a) )then
                ! i a j a  ->  i & j alpha for abab and beta for baba
                iSpin_abab = iSpin
                jSpin_abab = jSpin
                iSpin_baba = iSpin-1
                jSpin_baba = jSpin-1
                t_abab_only = .false. 
            else if ((i .gt. a) .or. (j .gt. a) )then
                ! a i a j->  i & j alpha for baba and beta for abab
                iSpin_abab = iSpin-1
                jSpin_abab = jSpin-1
                iSpin_baba = iSpin
                jSpin_baba = jSpin
                t_abab_only = .false. 
            else if ((i .eq. a) .and. (j .eq. a) )then
                ! a a a a -> i = a = j  abab and baba saved in abab array only! 
                ! -> count twice for close shell systems (fac_doublecount)
                iSpin_abab = iSpin
                jSpin_abab = jSpin
                iSpin_baba = iSpin-1
                jSpin_baba = jSpin-1
                t_abab_only = .true. 
                if (.not. tOpenShell) fac_doublecount=2.0_dp
            end if

            if (tRDMInstEnergy) then
                RDMEnergy_Inst = RDMEnergy_Inst + &
                           fac_doublecount*( (abab_RDM(Ind1_1e_ab,Ind2_1e_ab) ) &
                                               * real(TMAT2D(iSpin_abab,jSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )

                if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst + &
                           ( (abab_RDM(Ind2_1e_ab, Ind1_1e_ab) ) &
                                               * real(TMAT2D(jSpin_abab,iSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )
                if (tOpenShell) then !add baba terms
                    if (.not. t_abab_only) then
                        RDMEnergy_Inst = RDMEnergy_Inst + &
                                       ( (baba_RDM(Ind1_1e_ab,Ind2_1e_ab) ) &
                                       * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                       * (1.0_dp / real(NEl - 1,dp)) )

                        if (t_opposite_contri)  RDMEnergy_Inst = RDMEnergy_Inst + &
                                       ( (baba_RDM(Ind2_1e_ab,Ind1_1e_ab) ) &
                                       * real(TMAT2D(jSpin_baba,iSpin_baba),dp) &
                                       * (1.0_dp / real(NEl - 1,dp)) )
                    else   ! i=j=a -> baba_RDM saved in abab_RDM & t_opposite_contri = false
                        RDMEnergy_Inst = RDMEnergy_Inst + &
                                       ( (abab_RDM(Ind1_1e_ab,Ind2_1e_ab) ) &
                                       * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                       * (1.0_dp / real(NEl - 1,dp)) )
                    end if
                end if
            end if 

            RDMEnergy1 = RDMEnergy1 + fac_doublecount*( (abab_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                               * real(TMAT2D(iSpin_abab,jSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )

            if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + ( (abab_RDM_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM) &
                                               * real(TMAT2D(jSpin_abab,iSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )

            if (tOpenShell) then  ! add baba terms
                if (.not. t_abab_only) then
                    RDMEnergy1 = RDMEnergy1 + ( (baba_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                 * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                 * (1.0_dp / real(NEl - 1,dp)) )

                     if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + ( (baba_RDM_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM) &
                                 * real(TMAT2D(jSpin_baba,iSpin_baba),dp) &
                                 * (1.0_dp / real(NEl - 1,dp)) )
                else ! baba_RDM saved in abab_RDM  & t_opposite_contri = false
                    RDMEnergy1 = RDMEnergy1 + ( (abab_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                 * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                 * (1.0_dp / real(NEl - 1,dp)) )
                end if

            end if


            if ((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy ) then
               
                if (.not. tOpenShell) then
                 ! Spatial orbitals
                    NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                             NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                         + fac_doublecount*( abab_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                                 * (1.0_dp / real(NEl - 1,dp)) )

                    if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                             NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                         + ( abab_RDM_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM &
                                                 * (1.0_dp / real(NEl - 1,dp)) )

                else                      
                   NatOrbMat(SymLabelListInv_rot(iSpin_abab),SymLabelListInv_rot(jSpin_abab)) = &
                           NatOrbMat(SymLabelListInv_rot(iSpin_abab),SymLabelListInv_rot(jSpin_abab)) &
                                       + ( abab_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                       * (1.0_dp / real(NEl - 1,dp)) )

                   if (t_opposite_contri)  NatOrbMat(SymLabelListInv_rot(jSpin_abab),SymLabelListInv_rot(iSpin_abab)) = &
                           NatOrbMat(SymLabelListInv_rot(jSpin_abab),SymLabelListInv_rot(iSpin_abab)) &
                                       + ( abab_RDM_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM &
                                       * (1.0_dp / real(NEl - 1,dp)) ) 

                   if (.not. t_abab_only) then
                       NatOrbMat(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) = &
                            NatOrbMat(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) &
                                      + ( baba_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                      * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(jSpin_baba),SymLabelListInv_rot(iSpin_baba)) = &
                            NatOrbMat(SymLabelListInv_rot(jSpin_baba),SymLabelListInv_rot(iSpin_baba)) &
                                      + ( baba_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                      * (1.0_dp / real(NEl - 1,dp)) )
                   else ! i = j = a -> baba saved in abab & t_opposite_contri = false
                       NatOrbMat(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) = &
                          NatOrbMat(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) &
                                      + ( abab_RDM_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                      * (1.0_dp / real(NEl - 1,dp)) )
                   end if
                end if

           end if
 
       end if ! abab & baba terms

       !abba & baab & aaaa & bbbb terms
       if ((i.ne.a).and.(j.ne.a)) then   
           Ind1_1e_aa = ( ( (max(i,a)-2) * (max(i,a)-1) ) / 2 ) + min(i,a)
           Ind2_1e_aa = ( ( (max(j,a)-2) * (max(j,a)-1) ) / 2 ) + min(j,a)

           if (Ind1_1e_aa.ne.Ind2_1e_aa) then
           ! For Gamma elements corresponding to 1-RDMs (eg Gamma(i,a,a,j) ), 
           ! we're only considering i =< j 
           ! therefore we need to sum in the opposite contribution too.
               t_opposite_contri = .true.
           else  ! no opposite contribution from i = j term
               t_opposite_contri = .false.
           end if

           !abba & baab terms
           if ((i.ne.j).and.((i.lt.a).and.(j.gt.a)).or.((i.gt.a).and.(j.lt.a))) then 
               if ((i.lt.a).and.(j.gt.a))then
                   ! i a a j -> i & j alpha for abba and beta for baab
                   iSpin_abba = iSpin
                   jSpin_abba = jSpin
                   iSpin_baab = iSpin-1
                   jSpin_baab = jSpin-1
               else if ((i.gt.a).and.(j.lt.a))then
                   ! a i j a  -> i & j beta for abba and alpha for baab
                   iSpin_abba = iSpin-1
                   jSpin_abba = jSpin-1
                   iSpin_baab = iSpin
                   jSpin_baab = jSpin
               end if

               if (tRDMInstEnergy) then
                   RDMEnergy_Inst = RDMEnergy_Inst - &
                                                  ( (abba_RDM(Ind1_1e_aa,Ind2_1e_aa) ) &
                                                  * real(TMAT2D(iSpin_abba,jSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                   if (t_opposite_contri .and. (.not. tOpenShell) )then
                       RDMEnergy_Inst = RDMEnergy_Inst - &
                                                 ( (abba_RDM(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
                   end if

                   if (tOpenShell) then
                       ! abba becomes baab when i and j are swapped for the opposite contribution
                       if (t_opposite_contri ) RDMEnergy_Inst = RDMEnergy_Inst - &
                                                 ( (baab_RDM(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                       RDMEnergy_Inst = RDMEnergy_Inst - &
                                                  ( (baab_RDM(Ind1_1e_aa,Ind2_1e_aa) ) &
                                                  * real(TMAT2D(iSpin_baab,jSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst - &
                                                 ( (abba_RDM(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                  * real(TMAT2D(jSpin_baab,iSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
                   end if

               end if 

               RDMEnergy1 = RDMEnergy1 - ( (abba_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(iSpin_abba,jSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

               if (t_opposite_contri .and. (.not. tOpenShell) ) then
                   RDMEnergy1 = RDMEnergy1 - ( (abba_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
               end if

               if (tOpenShell) then  !add baab terms
                   ! abba becomes baab when i and j are swapped for the opposite contribution
                   if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 - ( (baab_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                   RDMEnergy1 = RDMEnergy1 - ( (baab_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(iSpin_baab,jSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                   if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 - ( (abba_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(jSpin_baab,iSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
               end if
                                                  
               if ((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy ) then
                   if (.not.tOpenShell) then
                   ! Spatial orbitals
                       NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                           NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                      - ( abba_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                              * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                          NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                      - ( abba_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                              * (1.0_dp / real(NEl - 1,dp)) )
                   else
                       NatOrbMat(SymLabelListInv_rot(iSpin_abba),SymLabelListInv_rot(jSpin_abba)) = &
                              NatOrbMat(SymLabelListInv_rot(iSpin_abba),SymLabelListInv_rot(jSpin_abba)) &
                                          - ( abba_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) ) 

                       if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(jSpin_abba),SymLabelListInv_rot(iSpin_abba)) = &
                              NatOrbMat(SymLabelListInv_rot(jSpin_abba),SymLabelListInv_rot(iSpin_abba)) &
                                          - ( baab_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) ) 

                       NatOrbMat(SymLabelListInv_rot(iSpin_baab),SymLabelListInv_rot(jSpin_baab)) = &
                              NatOrbMat(SymLabelListInv_rot(iSpin_baab),SymLabelListInv_rot(jSpin_baab)) &
                                          - ( baab_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(jSpin_baab),SymLabelListInv_rot(iSpin_baab)) = &
                              NatOrbMat(SymLabelListInv_rot(jSpin_baab),SymLabelListInv_rot(iSpin_baab)) &
                                          - ( abba_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
                   end if
               end if                   

           end if ! end of abba and baab terms

           ! aaaa & bbbb terms
           if (((i.lt.a).and.(j.lt.a)).or.((i.gt.a).and.(j.gt.a))) then
               Parity_Factor = 1.0_dp
           else
               Parity_Factor = -1.0_dp
           end if
          
           if (tRDmInstEnergy) then
               RDMEnergy_Inst = RDMEnergy_Inst + &
                          ( (aaaa_RDM(Ind1_1e_aa,Ind2_1e_aa) ) &
                                                * real(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
               if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst + &
                          ( (aaaa_RDM(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                * real(TMAT2D(jSpin,iSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
               if (tOpenShell) then
                   RDMEnergy_Inst = RDMEnergy_Inst +&
                          ( (bbbb_RDM(Ind1_1e_aa,Ind2_1e_aa)) &
                                                * real(TMAT2D(iSpin-1,jSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor)

                   if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst +&
                          ( (bbbb_RDM(Ind2_1e_aa,Ind1_1e_aa)) &
                                                * real(TMAT2D(jSpin-1,iSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor)
               end if

           end if 

           RDMEnergy1 = RDMEnergy1 + &
                           ( (aaaa_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

           if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + &
                           ( (aaaa_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(jSpin,iSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

           if (tOpenShell) then
               RDMEnergy1 = RDMEnergy1 + &
                           ( (bbbb_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(iSpin-1,jSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
               if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + &
                           ( (bbbb_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(jSpin-1,iSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
           end if

           if ((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then
               if (.not.tOpenShell)then
                   NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                            NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                        + ( aaaa_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                            NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                        + ( aaaa_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

               else

                   NatOrbMat(SymLabelListInv_rot(iSpin),SymLabelListInv_rot(jSpin)) = &
                            NatOrbMat(SymLabelListInv_rot(iSpin),SymLabelListInv_rot(jSpin)) &
                                        + ( aaaa_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(jSpin),SymLabelListInv_rot(iSpin)) = &
                            NatOrbMat(SymLabelListInv_rot(jSpin),SymLabelListInv_rot(iSpin)) &
                                        + ( aaaa_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   NatOrbMat(SymLabelListInv_rot(iSpin-1),SymLabelListInv_rot(jSpin-1)) = & 
                            NatOrbMat(SymLabelListInv_rot(iSpin-1),SymLabelListInv_rot(jSpin-1)) &
                                        + ( bbbb_RDM_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   if (t_opposite_contri) NatOrbMat(SymLabelListInv_rot(jSpin-1),SymLabelListInv_rot(iSpin-1)) = & 
                            NatOrbMat(SymLabelListInv_rot(jSpin-1),SymLabelListInv_rot(iSpin-1)) &
                                        + ( bbbb_RDM_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 
               end if
           end if

       end if 

    end subroutine calc_1RDM_energy

    subroutine find_nat_orb_occ_numbers()

        ! Diagonalises the 1-RDM (NatOrbMat), so that after this routine
        ! NatOrbMat is the eigenfunctions of the 1-RDM (the matrix transforming
        ! the MO's into the NOs). This also gets the NO occupation numbers
        ! (evaluse) and correlation entropy.

        integer :: ierr
        real(dp) :: SumDiag
        character(len=*), parameter :: this_routine='find_nat_orb_occ_numbers'

        if (iProcIndex .eq. 0) then
            
            ! Diagonalises the 1-RDM.  NatOrbMat goes in as the 1-RDM, comes out
            ! as the eigenvector of the 1-RDM (the matrix transforming the MO's
            ! into the NOs).
            call DiagRDM(SumDiag)

            ! Writes out the NO occupation numbers and evectors to files.
            call write_evales_and_transform_mat(SumDiag)

            if (tPrintRODump .and. tROHF) then               
                write(6,*) 'ROFCIDUMP not implemented for ROHF. Skip generation of ROFCIDUMP file.'
            else if (tPrintRODump) then
                allocate(FourIndInts(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat=ierr)
                call LogMemAlloc('FourIndInts',(NoOrbs**4), 8, this_routine, &
                                                        FourIndIntsTag, ierr)
                if (ierr .ne. 0) call Stop_All(this_routine, 'Problem allocating FourIndInts array,')

                ! Then, transform2ElInts.
                write(6,*) ''
                write(6,*) 'Transforming the four index integrals'
                call Transform2ElIntsMemSave_RDM()

                write(6,*) 'Re-calculating the fock matrix'
                call CalcFOCKMatrix_RDM()

                write(6,*) 'Refilling the UMAT and TMAT2D'

                ! The ROFCIDUMP is also printed out in here.
                call RefillUMATandTMAT2D_RDM()

                call neci_flush(6)

                call writebasis(6, G1, NoOrbs, ARR, BRR)

            end if
        end if

    end subroutine find_nat_orb_occ_numbers

    subroutine write_evales_and_transform_mat(SumDiag)

        real(dp), intent(in) :: SumDiag
        integer :: i, j, Evalues_unit, NatOrbs_unit, jSpat, jInd, NO_Number
        integer :: i_no,i_normal
        real(dp) :: Corr_Entropy, Norm_Evalues, SumN_NO_Occ
        logical :: tNegEvalue, tWrittenEvalue

        if (tOpenShell) then
            Norm_Evalues = SumDiag/real(NEl,dp)
        else
            Norm_Evalues = 2.0_dp*(SumDiag/real(NEl,dp))
        end if

        ! Write out normalised evalues to file and calculate the correlation
        ! entropy.
        Corr_Entropy = 0.0_dp
        Evalues_unit = get_free_unit()
        open(Evalues_unit,file='NO_OCC_NUMBERS',status='unknown')

        write(Evalues_unit,'(A)') '# NOs (natural orbitals) ordered by occupation number' 
        write(Evalues_unit,'(A)') '# MOs (HF orbitals) ordered by energy' 
        write(Evalues_unit,'(A1,A5,A30,A20,A30)') '#','NO','NO OCCUPATION NUMBER','MO','MO OCCUPATION NUMBER'
        tNegEvalue = .false.
        SumN_NO_Occ = 0.0_dp
        NO_Number = 1
        do i = 1, NoOrbs
            if (tOpenShell) then
                write(Evalues_unit,'(I6,G35.17,I15,G35.17)') i, Evalues(i)/Norm_Evalues, &
                                                                BRR(i), Rho_ii(i)
                if (Evalues(i) .gt. 0.0_dp) then
                    Corr_Entropy = Corr_Entropy - ( abs(Evalues(i)/ Norm_Evalues) &
                                                    * LOG(abs(Evalues(i)/ Norm_Evalues)) )
                else
                    tNegEvalue = .true.
                end if
                if (i .le. NEl) SumN_NO_Occ = SumN_NO_Occ + (Evalues(i)/Norm_Evalues)
            else
                write(Evalues_unit,'(I6,G35.17,I15,G35.17)') (2*i)-1,Evalues(i)/Norm_Evalues, &
                                                            BRR((2*i)-1), Rho_ii(i)/2.0_dp
                if (Evalues(i).gt.0.0_dp) then
                    Corr_Entropy = Corr_Entropy - (2.0_dp * ( abs(Evalues(i)/Norm_Evalues) &
                                                    * LOG(abs(Evalues(i)/Norm_Evalues)) ) )
                else
                    tNegEvalue = .true.
                end if
                write(Evalues_unit,'(I6,G35.17,I15,G35.17)') 2*i,Evalues(i)/Norm_Evalues, &
                                                            BRR(2*i), Rho_ii(i)/2.0_dp
                if (i .le. (NEl/2)) SumN_NO_Occ = SumN_NO_Occ + (2.0_dp * (Evalues(i)/Norm_Evalues))
            end if
        end do
        close(Evalues_unit)

        write(6,'(1X,A45,F30.20)') 'SUM OF THE N LARGEST NO OCCUPATION NUMBERS: ',SumN_NO_Occ
   
        write(6,'(1X,A20,F30.20)') 'CORRELATION ENTROPY', Corr_Entropy
        write(6,'(1X,A33,F30.20)') 'CORRELATION ENTROPY PER ELECTRON', Corr_Entropy / real(NEl,dp) 
        if (tNegEvalue) write(6,'(1X,"WARNING: Negative NO occupation numbers found.")')

        ! Write out the evectors to file.
        ! This is the matrix that transforms the molecular orbitals into the
        ! natural orbitals. Evalue(i) corresponds to Evector NatOrbsMat(1:nBasis,i)
        ! We just want the Evalues in the same order as above, but the
        ! 1:nBasis part (corresponding to the molecular orbitals), needs to
        ! refer to the actual orbital labels. Want these orbitals to preferably
        ! be in order, run through the orbital, need the position to find the
        ! corresponding NatOrbs element, use SymLabelListInv_rot
        if (.not. tNoNOTransform) then
            NatOrbs_unit = get_free_unit()
            open(NatOrbs_unit, file='NO_TRANSFORM', status='unknown')
            write(NatOrbs_unit,'(2A6,2A30)') '#   MO', 'NO', 'Transform Coeff', 'NO OCC NUMBER'
            ! write out in terms of spin orbitals, all alpha then all beta.
            NO_Number = 1
            do i_normal = 1, NoOrbs
                i_no = SymLabelListInv_rot(i_normal)
                tWrittenEvalue = .false.
                do j = 1, nBasis
                    ! Here i corresponds to the natural orbital, and j to the
                    ! molecular orbital. i is actually the spin orbital in
                    ! this case.
                    if (tOpenShell) then
                        jInd = j
                    else
                        if (mod(j,2).ne.0) then
                            jInd = gtID(j)
                        else
                            cycle
                        end if
                    end if

                    if (tWrittenEvalue) then
                        if (NatOrbMat(SymLabelListInv_rot(jInd),i_no) .ne. 0.0_dp) &
                            write(NatOrbs_unit,'(2I6,G35.17)') j,NO_Number,NatOrbMat(SymLabelListInv_rot(jInd),i_no)
                    else
                        if (NatOrbMat(SymLabelListInv_rot(jInd),i_no) .ne. 0.0_dp) then 
                            write(NatOrbs_unit,'(2I6,2G35.17)') j,NO_Number,NatOrbMat(SymLabelListInv_rot(jInd),i_no),&
                                                                                Evalues(i_no)/Norm_Evalues
                            tWrittenEvalue = .true.
                        end if
                    end if
                end do

                NO_Number = NO_Number + 1
                if (.not.tOpenShell) then
                    tWrittenEvalue = .false.
                    do j = 2, nBasis, 2
                        ! Here i corresponds to the natural orbital, and j to
                        ! the molecular orbital. i is actually the spin orbital
                        ! in this case.
                        jSpat = gtID(j)
                        if (tWrittenEvalue) then
                            if (NatOrbMat(SymLabelListInv_rot(jSpat),i_no) .ne. 0.0_dp) &
                                write(NatOrbs_unit,'(2I6,G35.17)') j,NO_Number,&
                                                                NatOrbMat(SymLabelListInv_rot(jSpat),i_no)
                        else
                            if (NatOrbMat(SymLabelListInv_rot(jSpat),i_no) .ne. 0.0_dp) then
                                write(NatOrbs_unit,'(2I6,2G35.17)') j,NO_Number,&
                                                        NatOrbMat(SymLabelListInv_rot(jSpat),i_no), &
                                                        Evalues(i_no)/Norm_Evalues
                                tWrittenEvalue = .true.
                            end if
                        end if

                    end do
                    NO_Number = NO_Number + 1
                end if
            end do
            close(NatOrbs_unit)
        end if

    end subroutine write_evales_and_transform_mat

    subroutine DiagRDM(SumTrace)

        ! The diagonalisation routine reorders the orbitals in such a way that
        ! the corresponding orbital labels are lost. In order to keep the spin
        ! and spatial symmetries, each symmetry must be fed into the
        ! diagonalisation routine separately. The best way to do this is to
        ! order the orbitals so that all the alpha orbitals follow all the beta
        ! orbitals, with the occupied orbitals first, in terms of symmetry, and
        ! the virtual second, also ordered by symmetry. This gives us
        ! flexibility w.r.t rotating only the occupied or only virtual and 
        ! looking at high spin states.

        real(dp), intent(out) :: SumTrace
        real(dp) :: SumDiagTrace
        real(dp), allocatable :: WORK2(:),EvaluesSym(:),NOMSym(:,:)
        integer :: ierr,i,j,spin,Sym,LWORK2,WORK2Tag,SymStartInd,NoSymBlock
        integer :: EvaluesSymTag,NOMSymTag,k,MaxSym
        logical :: tDiffSym, tDiffLzSym
        character(len=*), parameter :: this_routine='DiagRDM'

        ! Test that we're not breaking symmetry.
        ! And calculate the trace at the same time.
        SumTrace = 0.0_dp

        do i = 1, NoOrbs
            do j = 1, NoOrbs
                tDiffSym = .false.
                tDiffLzSym = .false.
                if (tOpenShell) then
                    if ((int(G1(SymLabelList2_rot(i))%sym%S).ne.&
                        int(G1(SymLabelList2_rot(j))%sym%S))) tDiffSym = .true.
                    if ((int(G1(SymLabelList2_rot(i))%Ml).ne.&
                        int(G1(SymLabelList2_rot(j))%Ml))) tDiffLzSym = .true.
                else
                    if ((int(G1(2*SymLabelList2_rot(i))%sym%S).ne.&
                        int(G1(2*SymLabelList2_rot(j))%sym%S))) tDiffSym = .true.
                    if ((int(G1(2*SymLabelList2_rot(i))%Ml).ne.&
                        int(G1(2*SymLabelList2_rot(j))%Ml))) tDiffLzSym = .true.
                end if
                if (tDiffSym) then
                    if (abs(NatOrbMat(i,j)).ge.1.0E-15_dp) then
                        write(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i',&
                                                                'Sym j','Matrix value'
                        if (tOpenShell) then                                                              
                            write(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(SymLabelList2_rot(i))%sym%S),&
                                    int(G1(SymLabelList2_rot(j))%sym%S),NatOrbMat(i,j)
                        else
                            write(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(2*SymLabelList2_rot(i))%sym%S),&
                                    int(G1(2*SymLabelList2_rot(j))%sym%S),NatOrbMat(i,j)
                        end if
                        if (tUseMP2VarDenMat) then
                            write(6,*) '**WARNING** - There is a non-zero NatOrbMat &
                            &value between orbitals of different symmetry.'
                            write(6,*) 'These elements will be ignored, and the symmetry &
                            &maintained in the final transformation matrix.'
                        else
                            write(6,*) 'k,SymLabelList2_rot(k),SymLabelListInv_rot(k)'
                            do k = 1,NoOrbs
                                write(6,*) k,SymLabelList2_rot(k),SymLabelListInv_rot(k)
                            end do
                            call neci_flush(6)
                            call Stop_All(this_routine,'Non-zero NatOrbMat value between &
                            &different symmetries.')
                        end if
                    end if
                    NatOrbMat(i,j) = 0.0_dp
                end if
                if (tDiffLzSym) then
                    if (abs(NatOrbMat(i,j)).ge.1.0E-15_dp) then
                        write(6,'(6A8,A40)') 'i','j','Label i','Label j','Lz i',&
                                                                'Lz j','Matrix value'
                        if (tOpenShell) then                                                              
                            write(6,'(6I8,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(SymLabelList2_rot(i))%Ml),&
                                    int(G1(SymLabelList2_rot(j))%Ml),NatOrbMat(i,j)
                        else
                            write(6,'(6I8,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(2*SymLabelList2_rot(i))%Ml),&
                                    int(G1(2*SymLabelList2_rot(j))%Ml),NatOrbMat(i,j)
                        end if
                        write(6,'(A)') ' **WARNING** - There is a non-zero NatOrbMat element &
                        &between orbitals of different Lz symmetry.'
                    end if
                end if
            end do
            SumTrace = SumTrace + NatOrbMat(i,i)
        end do

        write(6,*) ''
        write(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        call neci_flush(6)

        ! If we want to maintain the symmetry, we cannot have all the orbitals
        ! jumbled up when the  diagonaliser reorders the eigenvectors. Must
        ! instead feed each symmetry block in separately. This means that
        ! although the transformed orbitals are jumbled within the symmetry blocks, 
        ! the symmetry labels are all that are relevant and these are unaffected.
        Sym = 0
        LWORK2 = -1
        if (tOpenShell) then
            if (tFixLz) then
                MaxSym = (16 * ( ( 2 * iMaxLz ) + 1 ) ) - 1
            else
                MaxSym = 15
            end if
        else
            if (tFixLz) then
                MaxSym = (8 * ( ( 2 * iMaxLz ) + 1 ) ) - 1
            else
                MaxSym = 7
            end if
        end if
        do while (Sym .le. MaxSym)

            NoSymBlock = SymLabelCounts2_rot(2,Sym+1)

            SymStartInd = SymLabelCounts2_rot(1,Sym+1) - 1
            ! This is one less than the index that the symmetry starts, so that when we 
            ! run through i=1,..., we can start at SymStartInd+i.

            if (NoSymBlock.gt.1) then

                allocate(NOMSym(NoSymBlock,NoSymBlock),stat=ierr)
                call LogMemAlloc('NOMSym',NoSymBlock**2,8,this_routine,NOMSymTag,ierr)
                if (ierr .ne. 0) call Stop_All(this_routine,"Problem allocating NOMSym.")
                allocate(EvaluesSym(NoSymBlock),stat=ierr)
                call LogMemAlloc('EvaluesSym',NoSymBlock,8,this_routine,EvaluesSymTag,ierr)
                if (ierr .ne. 0) call Stop_All(this_routine,"Problem allocating EvaluesSym.")

                LWORK2=3*NoSymBlock+1
                allocate(WORK2(LWORK2),stat=ierr)
                call LogMemAlloc('WORK2',LWORK2,8,this_routine,WORK2Tag,ierr)
                if (ierr .ne. 0) call Stop_All(this_routine,"Problem allocating WORK2.")

                do j=1,NoSymBlock
                    do i=1,NoSymBlock
                        NOMSym(i,j)=NatOrbMat(SymStartInd+i,SymStartInd+j)
                    end do
                end do

                call dsyev('V', 'L', NoSymBlock, NOMSym, NoSymBlock, EvaluesSym, WORK2, LWORK2, ierr)
                ! NOMSym goes in as the original NOMSym, comes out as the 
                ! eigenvectors (Coefficients).
                ! EvaluesSym comes out as the eigenvalues in ascending order.

                do i=1,NoSymBlock
                    Evalues(SymStartInd+i)=EvaluesSym(NoSymBlock-i+1)
                end do

                ! CAREFUL if eigenvalues are put in ascending order, this may not be 
                ! correct, with the labelling system.
                ! may be better to just take coefficients and transform TMAT2DRot 
                ! in transform2elints.
                ! a check that comes out as diagonal is a check of this routine anyway.

                do j=1,NoSymBlock
                    do i=1,NoSymBlock
                        NatOrbMat(SymStartInd+i,SymStartInd+j)=NOMSym(i,NoSymBlock-j+1)
                    end do
                end do
                ! Directly fill the coefficient matrix with the eigenvectors from 
                ! the diagonalization.

                deallocate(WORK2)
                call LogMemDealloc(this_routine,WORK2Tag)

                deallocate(NOMSym)
                call LogMemDealloc(this_routine,NOMSymTag)

                deallocate(EvaluesSym)
                call LogMemDealloc(this_routine,EvaluesSymTag)

            else if (NoSymBlock .eq. 1) then
                ! The eigenvalue is the lone value, while the eigenvector is 1.

                Evalues(SymStartInd+1) = NatOrbMat(SymStartInd+1, SymStartInd+1)
                NatOrbMat(SymStartInd+1, SymStartInd+1) = 1.0_dp
            end if

            Sym = Sym + 1
        end do

        write(6,*) 'Matrix diagonalised'
        call neci_flush(6)

        SumDiagTrace = 0.0_dp
        do i = 1,NoOrbs
            SumDiagTrace=SumDiagTrace+Evalues(i)
        end do
        if ((abs(SumDiagTrace-SumTrace)).gt.1.0_dp) then
            write(6,*) 'Sum of diagonal NatOrbMat elements : ',SumTrace
            write(6,*) 'Sum of eigenvalues : ',SumDiagTrace
            write(6,*) 'WARNING : &
            &The trace of the 1RDM matrix before diagonalisation is '
            write(6,*) 'not equal to that after.'
        end if

        ! The MO's still correspond to SymLabelList2_rot.
        ! Although the NO's are slightly jumbled, they are only jumbled within
        ! their symmetry blocks. They still correspond to the symmetries of
        ! SymLabelList2_rot, which is the important part.

        ! But in order to look at the output, it is easier to consider them in
        ! terms of highest occupied to lowest occupied - i.e. in terms of the
        ! NO eigenvalues (occupation numbers).
        call OrderNatOrbMat()

    end subroutine DiagRDM

    subroutine OrderNatOrbMat()

        integer :: spin,i,j,ierr,StartSort,EndSort
        character(len=*), parameter :: this_routine = 'OrderRDM'
        integer, allocatable :: SymLabelList3_rot(:)
        real(dp), allocatable :: NatOrbMatTemp(:,:), EvaluesTemp(:)
        integer :: NatOrbMatTempTag, SymLabelList3_rotTag, EvaluesTempTag, Orb, New_Pos
        
        ! Here, if symmetry is kept, we are going to have to reorder the
        ! eigenvectors according to the size of the eigenvalues, while taking
        ! the orbital labels (and therefore symmetries) with them. This will
        ! be put back into MP2VDM from MP2VDMTemp.

        ! Want to reorder the eigenvalues from largest to smallest, taking the
        ! eigenvectors with them and the symmetry as well. If using spin
        ! orbitals, do this for the alpha spin and then the beta.

        allocate(NatOrbMatTemp(NoOrbs,NoOrbs), stat=ierr)
        call LogMemAlloc('NatOrbMatTemp', NoOrbs**2, 8,&
                            'OrderNatOrbMat', NatOrbMatTempTag, ierr)
        allocate(SymLabelList3_rot(NoOrbs), stat=ierr)
        call LogMemAlloc('SymLabelList3_rot', NoOrbs, 4,&
                            'OrderNatOrbMat', SymLabelList3_rotTag, ierr)
        allocate(EvaluesTemp(NoOrbs), stat=ierr)
        call LogMemAlloc('EvaluesTemp', NoOrbs, 4,&
                            'OrderNatOrbMat', EvaluesTempTag, ierr)

        ! Want to remember the original orbital ordering, as after the sort,
        ! the MO's will still have this ordering.
        SymLabelList3_rot(:) = SymLabelList2_rot(:)

        StartSort = 1
        EndSort = SpatOrbs

        ! Unfortunately this sort routine orders the orbitals in ascending
        ! order... which is not quite what we want.  Just remember this when
        ! printing out the Evalues.
        call sort (EValues(startSort:endSort), &
                   NatOrbMat(1:NoOrbs, startSort:endSort), &
                   SymLabelList2_rot(startSort:endSort))

        if (tOpenShell) then                  
            StartSort = SpatOrbs + 1
            EndSort = nBasis

            call sort (EValues(startSort:endSort), &
                       NatOrbMat(1:NoOrbs, startSort:endSort), &
                       SymLabelList2_rot(startSort:endSort))

        end if                       

        ! We now have the NO's ordered according to the size of their Evalues
        ! (occupation numbers).  This will have jumbled up their symmetries.
        ! Want to reorder the MO's to match this ordering (so that we only
        ! have one SymLabelList array).

        ! Need a new SymLabelListInv_rot too.        
        SymLabelListInv_rot(:) = 0   
        do i = 1, NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(NoOrbs-i+1)) = i
        end do

        NatOrbMatTemp(:,:) = NatOrbMat(:,:)
        NatOrbMat(:,:) = 0.0_dp

        do i = 1, NoOrbs
            do j = 1, NoOrbs

                ! In position j, the MO orbital Orb is currently there.
                Orb = SymLabelList3_rot(j)

                ! Want to move it to the position the NO's are in.
                New_Pos = SymLabelListInv_rot(Orb)

                ! But we also want to reverse the order of everything... 
                NatOrbMat(New_Pos,NoOrbs - i + 1)=NatOrbMatTemp(j,i)
            end do
        end do

        SymLabelList3_rot(:) = SymLabelList2_rot(:)
        EvaluesTemp(:) = Evalues(:)
        do i = 1, NoOrbs
            SymLabelList2_rot(i) = SymLabelList3_rot(NoOrbs - i + 1)
            Evalues(i) = EvaluesTemp(NoOrbs - i + 1)
        end do

        deallocate(NatOrbMatTemp)
        deallocate(SymLabelList3_rot)
        deallocate(EvaluesTemp)

    end subroutine OrderNatOrbMat

    subroutine Transform2ElIntsMemSave_RDM()

        ! This is an M^5 transform, which transforms all the two-electron
        ! integrals into the new basis described by the Coeff matrix.
        ! This is v memory inefficient and currently does not use any spatial 
        ! symmetry information.

        integer :: i,j,k,l,a,b,g,d,ierr,Temp4indintsTag,a2,d2,b2,g2
        real(dp), allocatable :: Temp4indints(:,:)
        character(len=*), parameter :: this_routine='Transform2ElIntsMemSave_RDM'
#ifdef __CMPLX
        call stop_all('Transform2ElIntsMemSave_RDM', &
                    'Rotating orbitals not implemented for complex orbitals.')
#endif
        
        ! Zero arrays from previous transform.

        allocate(Temp4indints(NoOrbs,NoOrbs),stat=ierr)
        call LogMemAlloc('Temp4indints',NoOrbs**2,8,&
                            'Transform2ElIntsMemSave_RDM',Temp4indintsTag,ierr)
        if (ierr .ne. 0) call Stop_All('Transform2ElIntsMemSave_RDM',&
                                    'Problem allocating memory to Temp4indints.')
 
        FourIndInts(:,:,:,:) = 0.0_dp

        ! Calculating the two-transformed, four index integrals.

        ! The untransformed <alpha beta | gamma delta> integrals are found from 
        ! UMAT(UMatInd(i,j,k,l,0,0)
        ! All our arrays are in spin orbitals - if tStoreSpinOrbs is false,
        ! UMAT will be in spatial orbitals - need to account for this.

        ! Running through 1,NoOrbs - the actual orbitals corresponding to that
        ! index are given by SymLabelList2_rot

        do b = 1, NoOrbs
            b2 = SymLabelList2_rot(b)
            do d = 1, NoOrbs
                d2 = SymLabelList2_rot(d)
                do a = 1, NoOrbs
                    a2 = SymLabelList2_rot(a)
                    do g = 1, NoOrbs
                        g2 = SymLabelList2_rot(g)

                        ! UMatInd in physical notation, but FourIndInts in
                        ! chemical (just to make it more clear in these
                        ! transformations). This means that here, a and g are
                        ! interchangable, and so are b and d.
                        FourIndInts(a,g,b,d) = real(UMAT(UMatInd(a2,b2,g2,d2,0,0)),dp)
                    end do
                end do

                Temp4indints(:,:) = 0.0_dp
                call dgemm('T','N',NoOrbs,NoOrbs,NoOrbs,1.0_dp,NatOrbMat(:,:),NoOrbs,&
                            FourIndInts(1:NoOrbs,1:NoOrbs,b,d),NoOrbs,0.0_dp,&
                            Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs)

                ! Temp4indints(i,g) comes out of here, so to transform g to k, 
                ! we need the transpose of this.

                call dgemm('T','T',NoOrbs,NoOrbs,NoOrbs,1.0_dp,NatOrbMat(:,:),NoOrbs,&
                            Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs,0.0_dp,&
                            FourIndInts(1:NoOrbs,1:NoOrbs,b,d),NoOrbs)
                ! Get Temp4indits02(i,k).
            end do
        end do
        
        ! Calculating the 3 transformed, 4 index integrals.
        ! 01=a untransformed,02=b,03=g,04=d
        do i = 1, NoOrbs
            do k = 1, NoOrbs

                Temp4indints(:,:) = 0.0_dp
                call dgemm('T','N',NoOrbs,NoOrbs,NoOrbs,1.0_dp,NatOrbMat(:,:),NoOrbs,&
                            FourIndInts(i,k,1:NoOrbs,1:NoOrbs),NoOrbs,0.0_dp,&
                            Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs)

                call dgemm('T','T',NoOrbs,NoOrbs,NoOrbs,1.0_dp,NatOrbMat(:,:),&
                            NoOrbs,Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs,0.0_dp,&
                            FourIndInts(i,k,1:NoOrbs,1:NoOrbs),NoOrbs)
            end do
        end do

        deallocate(Temp4indints)
        call LogMemDeAlloc('Transform2ElIntsMemSave_RDM',Temp4indintsTag)
 
    end subroutine Transform2ElIntsMemSave_RDM

    subroutine CalcFOCKMatrix_RDM()

        ! Calculate the fock matrix in the nat orb basis.

        integer :: i,j,k,l,a,b,ierr,ArrDiagNewTag
        real(dp) :: FOCKDiagSumHF,FOCKDiagSumNew
        character(len=*), parameter :: this_routine='CalcFOCKMatrix_RDM'
        real(dp), allocatable :: ArrDiagNew(:)

        ! This subroutine calculates and writes out the fock matrix for the
        ! transformed orbitals.
        ! ARR is originally the fock matrix in the HF basis.
        ! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

        ! When transforming the orbitals into approximate natural orbitals, we
        ! want to save memory, so don't bother calculating the whole matrix,
        ! just the diagonal elements that we actually need.

        allocate(ArrDiagNew(nBasis),stat=ierr)
        if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating ArrDiagNew array,')
        call LogMemAlloc('ArrDiagNew',nBasis,8,this_routine,ArrDiagNewTag,ierr)
        ArrDiagNew(:)=0.0_dp                     

        ! First calculate the sum of the diagonal elements, ARR.
        ! Check if this is already being done.
        FOCKDiagSumHF = 0.0_dp
        do a = 1, nBasis        
            FOCKDiagSumHF = FOCKDiagSumHF+Arr(a,2)
        end do

        ! Then calculate the fock matrix in the transformed basis, and the sum
        ! of the new diagonal elements.

        ! Our Arr in spin orbitals.
!        do j=1,NoOrbs
!            ArrNew(j,j)=Arr(2*j,2)
!        end do

        FOCKDiagSumNew = 0.0_dp
        do j = 1, NoOrbs
            l = SymLabelList2_rot(j)
            if (tStoreSpinOrbs) then
                ArrDiagNew(l) = 0.0_dp
            else
                ArrDiagNew(2*l) = 0.0_dp
                ArrDiagNew((2*l)-1) = 0.0_dp
            end if
            do a = 1, NoOrbs
                b = SymLabelList2_rot(a)
                if (tStoreSpinOrbs) then
                    ArrDiagNew(l)=ArrDiagNew(l)+(NatOrbMat(a,j)*ARR(b,2)*NatOrbMat(a,j))
                else
                    ArrDiagNew(2*l)=ArrDiagNew(2*l)+(NatOrbMat(a,j)*ARR(2*b,2)*NatOrbMat(a,j))
                    ArrDiagNew((2*l)-1)=ArrDiagNew((2*l)-1)+(NatOrbMat(a,j)*ARR((2*b)-1,2)*NatOrbMat(a,j))
                end if
            end do
            if (tStoreSpinOrbs) then
                FOCKDiagSumNew = FOCKDiagSumNew + (ArrDiagNew(l))
            else
                FOCKDiagSumNew = FOCKDiagSumNew + (ArrDiagNew(2*l))
                FOCKDiagSumNew = FOCKDiagSumNew + (ArrDiagNew((2*l)-1))
            end if
        end do
        ! If we are truncation the virtual space, only the unfrozen entries will 
        ! be transformed.

        ! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2)
        ! (ordered in terms of orbital number).
        ! ARR(:,2) needs to be ordered in terms of symmetry and then energy
        ! (like SymLabelList), so currently this ordering will not be correct
        ! when reading in qchem intDUMPS as the orbital number ordering is by energy.

        do j = 1,nBasis
            ARR(j,2) = ArrDiagNew(j)
            ARR(j,1) = ArrDiagNew(BRR(j))
        end do

        deallocate(ArrDiagNew)
        call LogMemDealloc(this_routine, ArrDiagNewTag)

    endsubroutine CalcFOCKMatrix_RDM

    subroutine RefillUMATandTMAT2D_RDM()

        ! UMat is in spin or spatial orbitals, TMAT2D only spin.
        ! This routine refills these to more easily write out the ROFCIDUMP,
        ! and originally to be able to continue a calculation (although I doubt
        ! this works at the moment).

        integer :: l,k,j,i,a,b,g,d,c,nBasis2,TMAT2DPartTag,ierr
        real(dp) :: NewTMAT
        real(dp), allocatable :: TMAT2DPart(:,:)
        character(len=*), parameter :: this_routine='RefillUMATandTMAT2D_RDM'
#ifdef __CMPLX
        call stop_all('RefillUMATandTMAT2D_RDM', &
                    'Rotating orbitals not implemented for complex orbitals.')
#endif

        ! TMAT2D is always in spin orbitals.
        allocate(TMAT2DPart(nBasis,nBasis),stat=ierr)
        if (ierr .ne. 0) call Stop_All(this_routine,'Problem allocating TMAT2DPart array,')
        call LogMemAlloc('TMAT2DPart',nBasis*nBasis,8,&
                                    'RefillUMAT_RDM',TMAT2DPartTag,ierr)
        TMAT2DPart(:,:) = 0.0_dp

        ! Make the UMAT elements the four index integrals.
        ! These are calculated by transforming the HF orbitals using the
        ! coefficients that have been found.
        do l = 1, NoOrbs
            d = SymLabelList2_rot(l)
            do k = 1, NoOrbs
                b = SymLabelList2_rot(k)
                do j = 1, NoOrbs
                    g = SymLabelList2_rot(j)
                    do i = 1, NoOrbs
                        a = SymLabelList2_rot(i)
                        ! The FourIndInts are in chemical notation, the UMatInd
                        ! in physical.
                        UMAT(UMatInd(a,b,g,d,0,0)) = FourIndInts(i,j,k,l)
                    end do
                end do
            end do
        end do

        ! Also calculate the 2 index integrals, and make these the elements
        ! of the TMAT2D matrix. TMAT2D is in spin orbitals.

        do a = 1,nBasis
            do k = 1,NoOrbs
                i = SymLabelList2_rot(k)
                NewTMAT = 0.0_dp
                do b = 1,NoOrbs
                    d = SymLabelList2_rot(b)
                    if (tStoreSpinOrbs) then
                        NewTMAT = NewTMAT + (NatOrbMat(b,k)*real(TMAT2D(d,a),dp))
                    else
                        NewTMAT = NewTMAT + (NatOrbMat(b,k)*real(TMAT2D(2*d,a),dp))
                    end if
                end do
                if (tStoreSpinOrbs) then
                    TMAT2DPart(i,a) = NewTMAT
                else
                    if (mod(a,2).eq.0) then
                        TMAT2DPart(2*i,a) = NewTMAT
                    else
                        TMAT2DPart((2*i)-1,a) = NewTMAT
                    end if
                end if
            end do
        end do

        do k = 1,nBasis
            do l = 1, NoOrbs
                j = SymLabelList2_rot(l)
                NewTMAT = 0.0_dp
                do a = 1, NoOrbs
                    c = SymLabelList2_rot(a)
                    if (tStoreSpinOrbs) then
                        NewTMAT = NewTMAT+(NatOrbMat(a,l)*TMAT2DPart(k,c))
                    else
                        NewTMAT = NewTMAT+(NatOrbMat(a,l)*TMAT2DPart(k,2*c))
                    end if
                end do
                if (tStoreSpinOrbs) then
                    TMAT2D(k,j) = NewTMAT
                else
                    if (mod(k,2) .eq. 0) then
                        TMAT2D(k,2*j) = NewTMAT
                    else
                        TMAT2D(k,(2*j)-1) = NewTMAT
                    end if
                end if
            end do
        end do

        deallocate(TMAT2DPart)
        call LogMemDeAlloc('RefillUMAT_RDM',TMAT2DPartTag)

        if (.not.trotatedNOs) then
            call PrintROFCIDUMP_RDM("ROFCIDUMP")
        end if

    endsubroutine RefillUMATandTMAT2D_RDM

    subroutine PrintROFCIDUMP_RDM(filename)

        ! This prints out a new FCIDUMP file in the same format as the old one.

        integer :: i,j,k,l,iunit, orb
        character(len=9) :: filename

!        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
!        call set_timer(PrintROFCIDUMP_Time,30)

        iunit = get_free_unit()
        open(iunit,file=filename,status='unknown') !'ROFCIDUMP',status='unknown')
        
        write(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NoOrbs,&
                                                ',NELEC=',NEl,',MS2=',LMS,','
        write(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,NoOrbs
            if (tStoreSpinOrbs) then
                write(iunit,'(I1,A1)',advance='no') int(G1(i)%sym%S)+1,','
            else
                if (tRotatedNOs.and.tBrokenSymNOs) then
                    write(iunit,'(I1,A1)',advance='no') 1,','
                else
                    write(iunit,'(I1,A1)',advance='no') int(G1(i*2)%sym%S)+1,','
                end if
            end if
        end do

        write(iunit,*) ''

        if (tStoreSpinOrbs) then
            write(iunit,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        else
            write(iunit,'(A7,I1,A12)') 'ISYM=',1,' UHF=.FALSE.'
        end if

        if (tFixLz) then
            write(iunit,'(A7)',advance='no') 'SYML='
            do i = 1, NoOrbs
                if (i .eq. NoOrbs) then
                    write(iunit,'(I3,A1)') -20,','
                else
                    write(iunit,'(I3,A1)',advance='no') -20,','
                end if
            end do
            write(iunit,'(A8)',advance='no') 'SYMLZ='
            do i = 1, NoOrbs
                orb = i
                if (.not. tStoreSpinOrbs) orb = 2 * orb
                write(iunit, '(i3,",")', advance='no') int(g1(orb)%ml)
            end do
            write(iunit,*)
        end if

        write(iunit,'(A5)') '&end'
       
        do i = 1, NoOrbs
            do j = 1, NoOrbs
                do l = 1, j
                    ! Potential to put symmetry in here, have currently taken it out, 
                    ! because when we're only printing non-zero values, it is kind 
                    ! of unnecessary - although it may be used to speed things up.
                    do k = 1, i
                        ! UMatInd is in physical notation <ij|kl>, but the indices
                        ! printed in the FCIDUMP are in chemical notation (ik|jl).
                        if ((abs(real(UMat(UMatInd(i,j,k,l,0,0)),dp))).ne.0.0_dp) &
                                write(iunit,'(F21.12,4I3)') &
                                real(UMat(UMatInd(i,j,k,l,0,0)),dp), i, k, j, l 
 
                    end do
                end do
           end do
        end do

        ! TMAT2D stored as spin orbitals.
        do i=1,NoOrbs
            ! Symmetry?
            do k=1,NoOrbs
                if (tStoreSpinOrbs) then
                    if ((real(TMAT2D(i,k),dp)).ne.0.0_dp) write(iunit,'(F21.12,4I3)') &
                                                        real(TMAT2D(i,k),dp),i,k,0,0
                else
                    if ((real(TMAT2D(2*i,2*k),dp)).ne.0.0_dp) write(iunit,'(F21.12,4I3)') &
                                                        real(TMAT2D(2*i,2*k),dp),i,k,0,0
                end if
            end do
        end do

        ! ARR has the energies of the orbitals (eigenvalues).
        ! ARR(:,2) has ordering we want.
        ! ARR is stored as spin orbitals.

        do k=1,NoOrbs
            if (tStoreSpinOrbs) then
                write(iunit,'(F21.12,4I3)') Arr(k,2), k, 0, 0, 0
            else
                write(iunit,'(F21.12,4I3)') Arr(2*k,2), k, 0, 0, 0
            end if
        end do

        write(iunit,'(F21.12,4I3)') ECore, 0, 0, 0, 0
        
        call neci_flush(iunit)

        close(iunit)

!        call halt_timer(PrintROFCIDUMP_Time)

    endsubroutine PrintROFCIDUMP_RDM

    subroutine BrokenSymNO(occ_numb_diff)

        ! This rouine finds natural orbitals (NOs) whose occupation
        ! numbers differ by a small relative threshold (occ_numb_diff) and
        ! rotates them by calling the Rotate2Orbs routine in order to
        ! break symmetry and maximally localise the NOs
        ! We'd like to keep both the original NOs (Natural Orbitals) 
        ! and the broken maximally localised NOs for checking
        ! This is not very (time-) efficient at the moment

        real(dp), intent(in) :: occ_numb_diff
        real(dp) :: diffnorm,SumDiag,sum_old,sum_new,selfint_old
        real(dp), allocatable :: no_store(:,:)
        integer, allocatable :: symlist_store(:)
        real(dp), allocatable :: trans_2orbs_coeffs(:,:)
        real(dp), allocatable :: selfint(:)
        integer, allocatable :: rotate_list(:,:),rotorbs(:,:)
        integer :: l1,l2,l3,l4,l5,l6,m,n
        integer :: iumat,jumat
        logical :: partnerfound,localdelocal

        allocate(trans_2orbs_coeffs(2,2))
        allocate(rotorbs(6,2))
        allocate(rotate_list((NoOrbs*(NoOrbs-1)),4))
        allocate(selfint(NoOrbs))
        allocate(no_store(NoOrbs,NoOrbs))
        allocate(symlist_store(NoOrbs))

        if (iProcIndex.eq.0) then

            ! Need to store NO coefficients since these are overwritten during 
            ! the orbital rotation.
            no_store = NatOrbMat
            symlist_store = SymLabelList2_rot

            ! For usage of other routines.
            do l1=1,NoOrbs
                SymLabelList2_rot(l1) = l1
            end do

            ! Normalisation.
            SumDiag = sum(Evalues)

            if (tStoreSpinOrbs) then
                diffnorm = SumDiag/dble(NEl)
            else
                diffnorm = 2.0_dp*(SumDiag/dble(NEl))
            end if
            
            trotatedNOs = .true.
     
            if (tStoreSpinorbs) then
                call Stop_all("BrokenSymNO","Broken symmetry NOs currently not implemented for UHF")
            end if

            write(6,*) '------------------------------------------------------------------------------'
            write(6,*) 'Localising NOs whose occupation numbers differ by less than threshold'
            write(6,*) '------------------------------------------------------------------------------'
            if (tBreakSymNOs) then
                write(6,*) 'Rotating specified NOs'
            else
                write(6,*) 'Threshold for orbitals to rotate:',occ_numb_diff
            end if

            ! Self-interactions.
            selfint(:) = 0.0_dp
            do l1=1,NoOrbs
                selfint(l1) = Umat(UmatInd(l1,l1,l1,l1,0,0))
            end do

            write(6,*) 'Self-interactions for NOs:'
            do l1=1,NoOrbs
                write(6,'(I3,3X,G25.12)') l1, selfint(l1)
            end do
            write(6,*) 'Sum of NO selfinteractions:',sum(selfint)
            selfint_old = sum(selfint)

            ! If the NOs to be rotated are specified in the input file.
            if (tBreakSymNOs) then
                rotate_list(:,:) = 0
                rotorbs(:,:) = 0
                m = 0
                do l1 = 2,(2*rottwo),2
                    m = m + 1
                    rotate_list(m,1) = RotNOs(l1-1)
                    rotate_list(m,2) = RotNOs(l1)
                end do
                do l1 = ((2*rottwo)+3),((2*rottwo)+(3*rotthree)),3
                    m = m + 1
                    rotate_list(m,1) = RotNOs(l1-2)
                    rotate_list(m,2) = RotNOs(l1-1)
                    rotate_list(m,3) = RotNOs(l1)
                end do
                do l1 = ((2*rottwo)+(3*rotthree)+4),((2*rottwo)+(3*rotthree)+(4*rotfour)),4
                    m = m + 1
                    rotate_list(m,1) = RotNOs(l1-3)
                    rotate_list(m,2) = RotNOs(l1-2)
                    rotate_list(m,3) = RotNOs(l1-1)
                    rotate_list(m,4) = RotNOs(l1)
                end do
            else
                ! If the threshold is used to generate a list of NOs to be
                ! rotated.

                ! Generate the list of orbitals which are rotated amongst each
                ! other.
                rotate_list(:,:) = 0
                rotorbs(:,:) = 0
                ! Need to account for spatial and spin orbital representations
                ! since orbitals of different spin cannot be mixed.
                ! List contains the NOs which are rotated.
                ! It can deal with a maximum of four NOs which are mixed.
                m = 0
                n = 1
                do l1 = 1, NoOrbs
                    if ((m .ne. 0) .and. (l1 .le. rotate_list(m,n))) cycle
                        partnerfound = .false.
                        n = 1
                        do l2 = (l1+1), NoOrbs

                            if ((abs((Evalues(l1)/diffnorm)-(Evalues(l2)/diffnorm))/abs((Evalues(l2)/diffnorm)))&
                                &.lt.occ_numb_diff) then
                            if (.not.partnerfound) then
                                m = m + 1
                                n = n + 1
                                rotate_list(m,1) = l1
                                rotate_list(m,2) = l2
                                partnerfound = .true.
                            else if (partnerfound) then
                                n = n + 1
                                ! this is for up to 2-fold degenearcy
!                                if (n.gt.2) then
!                                    n = 2
!                                    write(6,*) '***Warning***'
!                                    write(6,*) 'Threshold generated more than 2-fold degeneracy'
!                                    write(6,*) 'NOs around:',l2
!                                    cycle  ! don't want to rotate more than 2 orbitals
!                                end if
                                ! this is for up to four-fold degeneracy
                                if (n.gt.4) then
                                    n = 4
                                    write(6,*) '***Warning***'
                                    write(6,*) 'Threshold generated more than 4-fold degeneracy'
                                    write(6,*) 'NOs around:',l2
                                    cycle  ! don't want to rotate more than 4 orbitals
                                end if
                                rotate_list(m,n) = l2
                            end if
                        end if
                    end do
                end do
            end if

            write(6,*) 'The following pairs of orbitals will be rotated:'
            do l1=1,m
                write(6,'(I3,3X,4(I3))') l1,rotate_list(l1,:)
            end do

            NatOrbMat(:,:) = 0.0_dp
            do l1 = 1,NoOrbs
                NatOrbMat(l1,l1) = 1.0_dp
            end do

            ! Rotate two-fold degenerate pairs first.
            do l1 = 1,m
                ! If only two orbitals have the same occupation numbers.
                if (rotate_list(l1,3).eq.0) then
                    write(6,'(A20,4(I3))') 'Rotating NOs:',rotate_list(l1,:)
                    iumat = rotate_list(l1,1)
                    jumat = rotate_list(l1,2)
                    if (jumat.le.local_cutoff) then
                        localdelocal = .false.
                    else if (jumat.gt.local_cutoff) then
                        localdelocal = .true.
                    end if
                    call Rotate2Orbs(iumat,jumat,trans_2orbs_coeffs,selfint(iumat),&
                        &selfint(jumat),localdelocal)
                    ! The new NOs are 
                    ! phi_{i'} = cos a p_{i} + sin a p_{j}
                    ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                    NatorbMat(iumat,iumat) = trans_2orbs_coeffs(1,1)
                    NatorbMat(jumat,iumat) = trans_2orbs_coeffs(2,1)
                    NatorbMat(iumat,jumat) = trans_2orbs_coeffs(1,2)
                    NatorbMat(jumat,jumat) = trans_2orbs_coeffs(2,2)

                    write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)

                    selfint_old = sum(selfint)
                end if
           end do
           ! Transform integral corresponding to rotated NO.
           ! These are required for not printing out RPFCIDUMP or 
           ! BSFCIDUMP every time.
           trotatedNOs = .true.
           call Transform2ElIntsMemSave_RDM()
           call RefillUMATandTMAT2D_RDM()

           ! If three orbitals are degenerate.
           do l1 = 1, m
                 if ((rotate_list(l1,3).ne.0).and.(rotate_list(l1,4).eq.0)) then

                        sum_new = sum(selfint)
                        rotorbs(1,1) = 1
                        rotorbs(1,2) = 2
                        rotorbs(2,1) = 1
                        rotorbs(2,2) = 3
                        rotorbs(3,1) = 2
                        rotorbs(3,2) = 3

                        ! These have to be done self-consistently since all
                        ! three orbitals can intermix.
                        do
                            sum_old = sum_new
                            write(6,'(A20,4(I3))') 'Rotating NOs:', rotate_list(l1,:)
                            do l3 = 1,3
                                iumat = rotate_list(l1,rotorbs(l3,1))
                                jumat = rotate_list(l1,rotorbs(l3,2))
                                NatOrbMat(:,:) = 0.0_dp
                                do l4=1,NoOrbs
                                    if ((l4 .ne. iumat).and.(l4 .ne. jumat)) then
                                        NatOrbMat(l4,l4) = 1.0_dp
                                    end if
                                end do
                                if (jumat .le. local_cutoff) then
                                    localdelocal = .false.
                                else if (jumat .gt. local_cutoff) then
                                    localdelocal = .true.
                                end if
                                call Rotate2Orbs(iumat,jumat,&
                                    &trans_2orbs_coeffs,selfint(iumat),selfint(jumat)&
                                    &,localdelocal)
                                ! The new NOs are 
                                ! phi_{i'} = cos a p_{i} + sin a p_{j}
                                ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                                NatorbMat(iumat,iumat) = trans_2orbs_coeffs(1,1)
                                NatorbMat(jumat,iumat) = trans_2orbs_coeffs(2,1)
                                NatorbMat(iumat,jumat) = trans_2orbs_coeffs(1,2)
                                NatorbMat(jumat,jumat) = trans_2orbs_coeffs(2,2)
                                ! Transform integral corresponding to rotated NO.
                                ! These are required for not printing out
                                ! RPFCIDUMP or BSFCIDUMP every time.
                                trotatedNOs = .true.
                                call Transform2ElIntsMemSave_RDM()
                                call RefillUMATandTMAT2D_RDM()
                            end do
                            ! Check for convergence.
                            sum_new = sum(selfint)
                            write(6,'(A50,2G20.12)') 'Current and previous selfinteraction:',&
                                &sum_new,sum_old
                            if (abs(sum_new-sum_old).lt.1e-12_dp) then
                                exit
                            end if
                        end do

                    write(6,*) 'Sum of rotated NO self-interactions:', sum(selfint)

                    selfint_old = sum(selfint)

                else if ((rotate_list(l1,3).ne.0).and.(rotate_list(l1,4).ne.0)) then

                        sum_new = sum(selfint)
                        rotorbs(1,1) = 1
                        rotorbs(1,2) = 2
                        rotorbs(2,1) = 3
                        rotorbs(2,2) = 4
                        rotorbs(3,1) = 1
                        rotorbs(3,2) = 3
                        rotorbs(4,1) = 2
                        rotorbs(4,2) = 4
                        rotorbs(5,1) = 1
                        rotorbs(5,2) = 4
                        rotorbs(6,1) = 2
                        rotorbs(6,2) = 3

                        ! These have to be done self-consistently since all
                        ! three orbitals can intermix.
                        do
                            sum_old = sum_new
                            write(6,'(A20,4(I3))') 'Rotating NOs:',rotate_list(l1,:)
                            do l3 = 1, 3
                                NatOrbMat(:,:) = 0.0_dp
                                do l4 = 1, NoOrbs
                                    if ((l4 .ne. rotate_list(l1,1)) .and. (l4 .ne. rotate_list(l1,2))&
                                        &.and.(l4.ne.rotate_list(l1,3)) .and. (l4 .ne. rotate_list(l1,4))) then
                                            NatOrbMat(l4,l4) = 1.0_dp
                                    end if
                                end do
                                ! Rotate these two independently.
                                do l5 = 0, 1
                                    iumat = rotate_list(l1,rotorbs(((2*l3)-l5),1))
                                    jumat = rotate_list(l1,rotorbs(((2*l3)-l5),2))
                                    if (jumat.le.local_cutoff) then
                                        localdelocal = .false.
                                    else if (jumat.gt.local_cutoff) then
                                        localdelocal = .true.
                                    end if
                                    call Rotate2Orbs(iumat,jumat,&
                                        &trans_2orbs_coeffs,selfint(iumat),selfint(jumat),localdelocal)
                                    ! The new NOs are 
                                    ! phi_{i'} = cos a p_{i} + sin a p_{j}
                                    ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                                    NatorbMat(iumat,iumat) = trans_2orbs_coeffs(1,1)
                                    NatorbMat(jumat,iumat) = trans_2orbs_coeffs(2,1)
                                    NatorbMat(iumat,jumat) = trans_2orbs_coeffs(1,2)
                                    NatorbMat(jumat,jumat) = trans_2orbs_coeffs(2,2)
                                end do
                                ! Transform integral corresponding to rotated NO.
                                ! These are required for not printing out
                                ! RPFCIDUMP or BSFCIDUMP every time.
                                trotatedNOs = .true.
                                call Transform2ElIntsMemSave_RDM()
                                call RefillUMATandTMAT2D_RDM()
                            end do
                            ! Check for convergence.
                            sum_new = sum(selfint)

                            write(6,"(A50,2G20.12)") 'Current and previous selfinteractions:',&
                                &sum_new,sum_old

                            if (abs(sum_new-sum_old) .lt. 1e-12_dp) then
                                exit
                            end if
                        end do

                    write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)

                    selfint_old = sum(selfint)
                end if
           end do

            write(6,*) 'Final self-interactions for rotated NOs:'
            do l1=1,NoOrbs
                write(6,'(I3,3X,G25.12)') l1, selfint(l1)
            end do
            write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)

            write(6,*) '------------------------------------------------------'
            write(6,*) 'Writing out BSFCIDUMP...'
            write(6,*) '------------------------------------------------------'
            
            call PrintROFCIDUMP_RDM("BSFCIDUMP")
 
            ! Restore original NOs.
            NatOrbMat(:,:) = 0.0_dp
            NatOrbMat = no_store
            SymLabelList2_rot = symlist_store


        end if

        deallocate(trans_2orbs_coeffs)
        deallocate(rotate_list)
        deallocate(rotorbs)
        deallocate(selfint)
        deallocate(no_store)
        deallocate(symlist_store)
        
        if (tBreakSymNOs) then
            deallocate(RotNOs)
            call LogMemDealloc('BrokenSymNO',tagRotNOs)
        end if

    endsubroutine BrokenSymNO

    subroutine Rotate2Orbs(i,j,trans_2orbs_coeffs,selfintorb1,selfintorb2,localdelocal)

        ! This routine takes two orbitals i,j, and rotates them in order to
        ! maximally localise these. It employs an Edminston-Ruedenberg type
        ! localisation which maximises the self-interaction
        ! \sum_{i=1}^{2} \sum_{r,s,u,v} (c_{ir})*(c_{is})*c_{iu}c_{iv} <p_{i}p_{i}|u|p_{i}p_{i}>
        ! where p_{i} are the original NOs.
        ! The coefficients c are given by the following matrix:
        ! c =  cos a   sin a
        !      -sin a  cos a
        ! Then angle a is found by differentiating and setting it equal to 0
        ! which gives the following analytical expression of the form
        ! tan a = -x/y
        ! where x and y are sums of the original NO four index integrals.

        real(dp), allocatable, intent(inout) :: trans_2orbs_coeffs(:,:)
        real(dp), intent(inout) :: selfintorb1,selfintorb2
        real(dp) :: alpha2(17)
        !real(dp) :: secondderiv(2)
        real(dp) :: selfinteractions(17)
        real(dp) :: coeffcos,coeffsin,maxint
        integer :: maxangle(1)
        integer :: indicesij(2)
        integer, intent(in) :: i,j
        integer :: l1,l2,l3,l4,l5
        logical, intent(in) :: localdelocal

        indicesij(1) = i
        indicesij(2) = j
        trans_2orbs_coeffs(:,:) = 0.0_dp

        ! Umat(UMatInd(i,j,k,l,0,0)) contains the four-index integrals
        ! <ij|kl> (physical notation) in the NO basis

        coeffcos = Umat(UmatInd(i,i,i,j,0,0)) + Umat(UmatInd(i,i,j,i,0,0)) + Umat(UmatInd(i,j,i,i,0,0)) &
            & - Umat(UmatInd(i,j,j,j,0,0)) + Umat(UmatInd(j,i,i,i,0,0)) - Umat(UmatInd(j,i,j,j,0,0)) &
            & - Umat(UmatInd(j,j,i,j,0,0)) - Umat(UmatInd(j,j,j,i,0,0))

        coeffsin = -Umat(UmatInd(i,i,i,i,0,0)) + Umat(UmatInd(i,i,j,j,0,0)) + Umat(UmatInd(i,j,i,j,0,0)) &
            & + Umat(UmatInd(i,j,j,i,0,0)) + Umat(UmatInd(j,i,i,j,0,0)) + Umat(UmatInd(j,i,j,i,0,0)) &
            & + Umat(UmatInd(j,j,i,i,0,0)) - Umat(UmatInd(j,j,j,j,0,0))

        ! atan return a value in [-pi/2,pi/2]
        ! because of the 4*alpha in the equation there are 8 distinct solutions
        ! i.e. in the range 0,2*pi
        ! i.e. possible solutions are separated by (2*pi/8)=pi/4
        ! for safety 16 solutions are evaluated.
        alpha2(9) = atan((-coeffcos/coeffsin))
        alpha2(9) = alpha2(9)/4.0_dp
        do l1=8,1,-1
            alpha2(l1) = alpha2(l1+1) - (pi/4.0_dp)
        end do
        do l1=10,17
            alpha2(l1) = alpha2(l1-1) + (pi/4.0_dp)
        end do
        
        !! second derivatives to find maximum (necessary since the minimum, i.e. fully delocalised
        !! orbitals satisfy the same conditions
        !secondderiv(1) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(1))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(1)))
        !secondderiv(2) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(2))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(2)))

        ! Compute selfinteractions to check which one is largest.
        ! This is a better measure than the second derivatives.
        selfinteractions(:) = 0.0_dp 

        do l1=1,17
            trans_2orbs_coeffs(1,1) = cos(alpha2(l1))
            trans_2orbs_coeffs(2,1) = sin(alpha2(l1))
            trans_2orbs_coeffs(1,2) = -sin(alpha2(l1))
            trans_2orbs_coeffs(2,2) = cos(alpha2(l1))

            do l2=1,2
                do l3=1,2
                    do l4=1,2
                        do l5=1,2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                                &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                        end do
                    end do
                end do
            end do
            do l2=1,2
                do l3=1,2
                    do l4=1,2
                        do l5=1,2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,2)&
                                &*trans_2orbs_coeffs(l3,2)*&
                                &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                                &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                        end do
                    end do
                end do
            end do
        end do

        ! Choose the angle which maximises the self interactions.
        if (.not.localdelocal) then
            ! Maximally delocalised.
            maxangle = minloc(selfinteractions)
            maxint = minval(selfinteractions)
        else if (localdelocal) then
            ! Maximally localised.
            maxangle = maxloc(selfinteractions)
            maxint = maxval(selfinteractions)
        end if


        ! Return transformatin coefficients.
        trans_2orbs_coeffs(1,1) = cos(alpha2(maxangle(1)))
        trans_2orbs_coeffs(2,1) = sin(alpha2(maxangle(1)))
        trans_2orbs_coeffs(1,2) = -sin(alpha2(maxangle(1)))
        trans_2orbs_coeffs(2,2) = cos(alpha2(maxangle(1)))
 
        ! New sefl-interactions for transformed orbitals.
        selfintorb1 = 0.0_dp
        selfintorb2 = 0.0_dp
        do l2=1,2
            do l3=1,2
                do l4=1,2
                    do l5=1,2
                        selfintorb1 = selfintorb1 + trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*&
                            &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                            &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                        selfintorb2 = selfintorb2 + trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*&
                            &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                            &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                    end do
                end do
            end do
        end do

    end subroutine Rotate2Orbs

    subroutine DeallocateRDM()

        ! This routine just deallocates the arrays allocated in InitRDM.
        ! If the NECI calculation softexits before the RDMs start to fill,
        ! this is all that is called at the end.

        character(len=*), parameter :: this_routine='DeallocateRDM'

        if (tExplicitAllRDM) then

            ! This array contains the initial positions of the single
            ! excitations for each processor.
            deallocate(Sing_InitExcSlots)
 
            ! This array contains the current position of the single
            ! excitations as they're added.
            deallocate(Sing_ExcList)

            ! This array actually contains the single excitations in blocks of
            ! the processor they will be sent to.        
            deallocate(Sing_ExcDjs)
            call LogMemDeAlloc(this_routine,Sing_ExcDjsTag)

            deallocate(Sing_ExcDjs2)
            call LogMemDeAlloc(this_routine,Sing_ExcDjs2Tag)


            if (RDMExcitLevel .ne. 1) then
                ! This array contains the initial positions of the
                ! single excitations for each processor.
                deallocate(Doub_InitExcSlots)
 
                ! This array contains the current position of the single
                ! excitations as they're added.
                deallocate(Doub_ExcList)

                ! This array actually contains the single excitations in
                ! blocks of the  processor they will be sent to.        
                deallocate(Doub_ExcDjs)
                call LogMemDeAlloc(this_routine,Doub_ExcDjsTag)
     
                deallocate(Doub_ExcDjs2)
                call LogMemDeAlloc(this_routine,Doub_ExcDjs2Tag)
            end if

        else

            if (allocated(Spawned_Parents)) then
                deallocate(Spawned_Parents)
                call LogMemDeAlloc(this_routine,Spawned_ParentsTag)
            end if

            if (allocated(Spawned_Parents_Index)) then
                deallocate(Spawned_Parents_Index)
                call LogMemDeAlloc(this_routine,Spawned_Parents_IndexTag)
            end if

        end if

        if (allocated(NatOrbMat)) then
            deallocate(NatOrbMat)
            call LogMemDeAlloc(this_routine,NatOrbMatTag)
        end if

        if (allocated(Evalues)) then
            deallocate(Evalues)
            call LogMemDeAlloc(this_routine,EvaluesTag)
        end if

        if (allocated(Rho_ii)) then
            deallocate(Rho_ii)
            call LogMemDeAlloc(this_routine,Rho_iiTag)
        end if

        if (allocated(FourIndInts)) then
            deallocate(FourIndInts)
            call LogMemDeAlloc(this_routine,FourIndIntsTag)
        end if

        if (allocated(SymLabelCounts2_rot)) then
            deallocate(SymLabelCounts2_rot)
            call LogMemDeAlloc(this_routine,SymLabelCounts2_rotTag)
        end if

        if (allocated(SymLabelList2_rot)) then
            deallocate(SymLabelList2_rot)
            call LogMemDeAlloc(this_routine,SymLabelList2_rotTag)
        end if

        if (allocated(SymLabelListInv_rot)) then
            deallocate(SymLabelListInv_rot)
            call LogMemDeAlloc(this_routine,SymLabelListInv_rotTag)
        end if

        if (allocated(aaaa_RDM_inst)) then
            deallocate(aaaa_RDM_inst)
            call LogMemDeAlloc(this_routine,aaaa_RDM_instTag)
        end if

        if (allocated(abab_RDM_inst)) then
            deallocate(abab_RDM_inst)
            call LogMemDeAlloc(this_routine,abab_RDM_instTag)
        end if

        if (allocated(abba_RDM_inst)) then
            deallocate(abba_RDM_inst)
            call LogMemDeAlloc(this_routine,abba_RDM_instTag)
        end if

        if (allocated(bbbb_RDM_inst)) then
            deallocate(bbbb_RDM_inst)
            call LogMemDeAlloc(this_routine,bbbb_RDM_instTag)
        end if

        if (allocated(baba_RDM_inst)) then
            deallocate(baba_RDM_inst)
            call LogMemDeAlloc(this_routine,baba_RDM_instTag)
        end if

        if (allocated(baab_RDM_inst)) then
            deallocate(baab_RDM_inst)
            call LogMemDeAlloc(this_routine,baab_RDM_instTag)
        end if

        if (allocated(aaaa_RDM_full)) then
            deallocate(aaaa_RDM_full)
            call LogMemDeAlloc(this_routine,aaaa_RDM_fullTag)
        end if

        if (allocated(abab_RDM_full)) then
            deallocate(abab_RDM_full)
            call LogMemDeAlloc(this_routine,abab_RDM_fullTag)
        end if

        if (allocated(abba_RDM_full)) then
            deallocate(abba_RDM_full)
            call LogMemDeAlloc(this_routine,abba_RDM_fullTag)
        end if

        if (allocated(bbbb_RDM_full)) then
            deallocate(bbbb_RDM_full)
            call LogMemDeAlloc(this_routine,bbbb_RDM_fullTag)
        end if

        if (allocated(baba_RDM_full)) then
            deallocate(baba_RDM_full)
            call LogMemDeAlloc(this_routine,baba_RDM_fullTag)
        end if

        if (allocated(baab_RDM_full)) then
            deallocate(baab_RDM_full)
            call LogMemDeAlloc(this_routine,baab_RDM_fullTag)
        end if


    end subroutine DeallocateRDM

    subroutine convert_mats_Molpforces(Norm_1RDM, Norm_2RDM)

    use SystemData, only: nEl,nbasis,LMS,ECore
    use IntegralsData, only: nFrozen
    use NatOrbsMod, only: NatOrbMat
    use SymData, only: Sym_Psi, SymLabelCounts,nSymLabels
    use SymExcitDataMod, only: SpinOrbSymLabel,SymLabelCounts2
    use GenRandSymExcitNUMod, only: ClassCountInd, RandExcitSymLabelProd
    use sym_mod
    use UMatCache, only: GTID
    
    integer :: iblkq, iseccr, istat1, isyref, ms2
    integer :: posn1, posn2
    integer :: i, j, k, l
    integer :: myname, ifil, intrel, iout
    integer :: orb1, orb2, Sym_i, Sym_j, Sym_ij
    integer :: Sym_k, Sym_l, Sym_kl
    integer, dimension(8) :: iact, ldact !iact(:) # of active orbs per sym, 
                                         !ldact(:) - # Pairs of orbs that multiply to give given sym
    integer, dimension(8) :: icore, iclos 
    integer, dimension(nSymLabels):: blockstart1, blockstart2
    integer, dimension(nSymLabels) :: elements_assigned1, elements_assigned2
    integer :: FC_Lag_Len  !Length of the Frozen Core Lagrangian
    integer :: Len_1RDM, Len_2RDM, FCLag_Len
    real(dp) :: ijkl, jikl
    real(dp), intent(in) :: Norm_1RDM, Norm_2RDM
    real(dp), allocatable :: SymmetryPacked2RDM(:), SymmetryPacked1RDM(:) !SymPacked2RDM in CHEMICAL Notation
    real(dp), allocatable :: SymmetryPackedLagrangian(:)
    real(dp), allocatable :: FC_Lagrangian(:)
    integer :: iWfRecord, iWfSym
    real(dp) :: WfRecord
    character(len=*), parameter :: t_r='convert_mats_Molpforces'

#ifdef __int64
    intrel = 1
#else
    intrel = 2
#endif

    if (iProcIndex .eq. 0) then
        ! Calculating header information needed for the molpro dump routine.
        ! Considering all electron for now.
        icore(:) = 0  ! Number of frozen orbitals per symmetry
        iclos(:) = 0  ! icore(i) + Number of 'closed' orbitals per symmetry (iclos is the same as icore for us).
        call molpro_get_reference_info(iWfRecord, iWfSym)
        iblkq = iWfRecord ! record number for orbitals
        iseccr = 0          ! record number for core orbitals
        istat1 = 1        !ground state
        isyref = Sym_Psi + 1  !spatial symmetry of the wavefunction
#ifdef MOLPRO
        if (isyref .ne. iWfSym) call stop_all(t_r,"NECI and common/cref do not agree on irrep of wave function")
#endif
        ms2 = LMS  !2 * M_s
        myname = 5001 !Arbitrary file names
        ifil = 1
        ! ^- cgk: might want to use nexfre(igrsav) for these two.
        iout = molpro_get_iout()
        ! ^- thise one should come from common/tapes. 
        ldact(:) = 0
        iact(:) = 0
        Len_1RDM = 0
        Len_2RDM = 0
        blockstart1(:) = 0
        blockstart2(:) = 0
        elements_assigned1(:) = 0
        elements_assigned2(:) = 0
        
        
        ! Find out the number of orbital pairs that multiply to a given sym (ldact).
        do i=1,SpatOrbs  
            ! Run over spatial orbitals.
            do j = 1,  i    ! i .ge. j
                Sym_i=SpinOrbSymLabel(2*i)  ! Consider only alpha orbitals.
                Sym_j=SpinOrbSymLabel(2*j)
                Sym_ij=RandExcitSymLabelProd(Sym_i, Sym_j)
                ldact(Sym_ij+1)=ldact(Sym_ij+1)+1
            end do
        end do

        ! Calculate lengths of arrays, and where each sym block starts.
        do i = 0, nSymLabels-1

            ! CMO: Check if Sym_i goes 0-->7 or 1-->8.
            ! Find position of each symmetry block in sym-packed forms of RDMS 1 & 2.
            blockstart1(i+1)=Len_1RDM+1 ! N.B. Len_1RDM still being updated in this loop.
            blockstart2(i+1)=Len_2RDM+1 ! N.B. Len_2RDM still being updated in this loop.
            
            ! Count the number of active orbitals of the given symmetry.
            iact(i+1)=SymLabelCounts2(2,ClassCountInd(1,i,0)) !Count the number of active orbitals of the given symmetry.
            Len_1RDM=Len_1RDM+(iact(i+1)*(iact(i+1)+1)/2) ! add on # entries in sym-packed 1RDM for sym i.
            
            Len_2RDM = Len_2RDM + (ldact(i+1))**2 ! Assumes no frozen orbitals.
        end do

        FCLag_Len = SpatOrbs**2  ! Arbitrarily set this for now - we will not be printing it whilst nfrozen=0
        
        ! Allocate arrays accordingly.
        allocate(SymmetryPacked1RDM(Len_1RDM))
        allocate(SymmetryPackedLagrangian(Len_1RDM))
        allocate(SymmetryPacked2RDM(Len_2RDM))
        allocate(FC_Lagrangian(FCLag_Len))

        FC_Lagrangian(:) = 0  ! Frozen-core Lagrandian -- whilst we do all electron calcs.
        
        ! Constructing the Symmetry Packed arrays.
        ! We convert our 1RDM, Lagrangian and  2RDM into the required Molpro
        ! symmetry-packed format 2RDM is stored as a full square matrix,
        ! separated into symmetry blocks. We store D_ijkl (chemical notation,
        ! spatial orbs) where (i .ge. j) and (k .ge. l). For each symmetry X,
        ! there will be a block where (i,j) and (k,l) both have symmetry X
        ! making (ij,kl) totally symmetric, and D_ijkl (potentially) non-zero.
        ! 1RDM and Lagrangian are stored as upper triangles, separated by
        ! symmetry block

        SymmetryPacked2RDM(:) = 0.0_dp
        SymmetryPacked1RDM(:) = 0.0_dp
        SymmetryPackedLagrangian(:) = 0.0_dp

        do i = 1, SpatOrbs  !run over spatial orbitals, ALL ELECTRON ONLY
            do j = 1,  i    ! i .ge. j
                Sym_i = SpinOrbSymLabel(2*i)  !Consider only alpha orbitals
                Sym_j = SpinOrbSymLabel(2*j)
                Sym_ij = RandExcitSymLabelProd(Sym_i, Sym_j)
                if (Sym_ij .eq. 0) then
                    posn1 = blockstart1(Sym_i+1) + elements_assigned1(Sym_i+1)
                    ! Add pre-symmetrised contribution to the symmetry-packed 1-RDM.

                    if (tOpenShell) then
                        ! Include both aa and bb contributions.
                        SymmetryPacked1RDM(posn1)=&
                              (NatOrbMat(SymLabelListInv_rot(2*i),SymLabelListInv_rot(2*j))&
                             + NatOrbMat(SymLabelListInv_rot(2*i-1),SymLabelListInv_rot(2*j-1)))*Norm_1RDM
                    else
                        SymmetryPacked1RDM(posn1) = NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM
                    end if

                    ! Add in the symmetrised Lagrangian contribution to the
                    ! sym-packed Lagrangian
                    SymmetryPackedLagrangian(posn1) = Lagrangian(i,j)
                    elements_assigned1(Sym_i+1) = elements_assigned1(Sym_i+1) + 1
                end if
                do k = 1, SpatOrbs
                    do l = 1, k
                        Sym_k = SpinOrbSymLabel(2*k)  ! Consider only alpha orbitals
                        Sym_l = SpinOrbSymLabel(2*l)
                        Sym_kl = RandExcitSymLabelProd(Sym_k, Sym_l)

                        if (Sym_kl .eq. Sym_ij) then
                            posn2=blockstart2(Sym_ij+1) + elements_assigned2(Sym_ij+1)
                        
                            ! Extracts spat orb chemical notation term from spin
                            ! separated physical notation RDMs 
                            ijkl = Find_Spatial_2RDM_Chem(i, j, k, l, Norm_2RDM)
                            jikl = Find_Spatial_2RDM_Chem(j, i, k, l, Norm_2RDM)

                            SymmetryPacked2RDM(posn2) = 0.5*(ijkl+jikl)
                            elements_assigned2(Sym_ij+1) = elements_assigned2(Sym_ij+1) + 1
                        end if
                    end do
                end do
            end do
        end do
        
        call molpro_dump_mcscf_dens_for_grad(myname,ifil, &
            icore, iclos, iact, nEL, isyref, ms2, iblkq,&
            iseccr, istat1, SymmetryPacked1RDM, Len_1RDM, &
            SymmetryPacked2RDM, Len_2RDM, SymmetryPackedLagrangian, &
            Len_1RDM, FC_Lagrangian, FC_Lag_Len, &
            iout, intrel)
        
    end if
        
    end subroutine convert_mats_Molpforces

    function Find_Spatial_2RDM_Chem(p,q,r,s, Norm_2RDM) result(pqrs)

    ! This routine calculates the spatial orbital, chemical notation 2RDM component
    !                  D_pqrs = <Psi | p+ r+ s q | Psi>
    ! This is achieved by looking up the D_pr,qs component in the various spin-separated
    ! versions of the 2RDM that are currently stored with spatial orb numbering, in 
    ! physical notation (e.g. aaaa_RDM_full etc).

    ! To convert from spin to spatial orbitals, we need to apply the following:
    ! D_pr,qs = D_pr,qs(aaaa) + D_pr,qs(bbbb) + D_pr,qs(abab) + D_pr,qs(baba) (Eq. ***)

    ! We note now the following quirks of the aaaa_RDM_full-type arrays for the manner in
    ! which they store these components
    !     1. In most cases the current RDMs store the *sum* of the spin-inverted terms
    !          - ie, aaaa_RDM_full(pr,qs) contains the sum of the aaaa and bbbb contributions
    !     2. When p=r and q=s, there is only one contribution generated in NECI
    !          - ie, abab_RDM_full(pp,qq) contains only one of the two identical abab and baba contributions
    !          - Terms of this kind but be explicitly multiplied by two to satisfy Eq. *** above
    !          - This is stored in the "Mult_Factor"
    !     3. The existing 2RDMs only store terms with r>=p and s>=q
    !          - If we wish to look up a term with a different order to this, we must swap the
    !            order of the indices, considering the swapped spin and introducing appropriate signs
    !          - ie if p>r and s>q, D_pr,qs(abab) is found by looking up -D_rp,qs(abba)
     
    integer, intent(in) :: p,q,r,s
    real(dp), intent(in) :: Norm_2RDM
    real(dp) :: pqrs
    real(dp) :: Mult_Factor
    integer :: Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab 

    pqrs = 0.0_dp

    if ((p.eq.r).and.(q.eq.s)) then
        Mult_Factor = 2.0_dp
    else
        Mult_Factor = 1.0_dp
    end if
    
    if ((r.ge.p) .and. (s.ge.q)) then ! D_pr,qs correctly ordered.
        Ind1_aa = ( ( (r-2) * (r-1) ) / 2 ) + p
        Ind1_ab = ( ( (r-1) * r ) / 2 ) + p
        Ind2_aa = ( ( (s-2) * (s-1) ) / 2 ) + q
        Ind2_ab = ( ( (s-1) * s ) / 2 ) + q

        pqrs = pqrs + Mult_Factor*abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM
        if (tOpenShell) pqrs = pqrs + Mult_Factor*baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM

        if ((p .ne. r) .and. (q .ne. s)) then 
            pqrs=pqrs+Mult_Factor*aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs + Mult_Factor*bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
        end if

    else if ((p .gt. r) .and. (s .gt. q)) then ! Need to reorder D_pr,qs to -D_rp,qs.
        Ind1_aa = ( ( (p-2) * (p-1) ) / 2 ) + r
        Ind2_aa = ( ( (s-2) * (s-1) ) / 2 ) + q

        pqrs=pqrs-Mult_Factor*abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
        if (tOpenShell) pqrs = pqrs - Mult_Factor*baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM

        pqrs=pqrs-Mult_Factor*aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
        if (tOpenShell) pqrs = pqrs - Mult_Factor*bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM

    else if ((r .gt. p) .and. (q .gt. s)) then ! Need to reorder D_pr,qs to -D_pr,sq.
        Ind1_aa = ( ( (r-2) * (r-1) ) / 2 ) + p
        Ind2_aa = ( ( (q-2) * (q-1) ) / 2 ) + s

        pqrs = pqrs - Mult_Factor*abba_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
        if (tOpenShell) pqrs = pqrs - Mult_Factor*baab_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM

        pqrs = pqrs - Mult_Factor*aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
        if (tOpenShell) pqrs = pqrs - Mult_Factor*bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM

    else ! Must need to reorder D_pr,qs to D_rp,sq.
        Ind1_aa = ( ( (p-2) * (p-1) ) / 2 ) + r
        Ind1_ab = ( ( (p-1) * p ) / 2 ) + r
        Ind2_aa = ( ( (q-2) * (q-1) ) / 2 ) + s
        Ind2_ab = ( ( (q-1) * q ) / 2 ) + s

        pqrs = pqrs + Mult_Factor*abab_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM
        if (tOpenShell) pqrs = pqrs + Mult_Factor*baba_RDM_full(Ind1_ab,Ind2_ab)*Norm_2RDM

        if ((p .ne. r) .and. (q .ne. s)) then
            pqrs = pqrs + Mult_Factor*aaaa_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs + Mult_Factor*bbbb_RDM_full(Ind1_aa,Ind2_aa)*Norm_2RDM
        end if

    end if
        
    end function Find_Spatial_2RDM_Chem

    subroutine molpro_set_igrdsav(igrsav_)

        implicit double precision(a-h,o-z)
        implicit integer(i-n)
        character(len=*), parameter :: t_r='molpro_set_igrdsav'
#ifdef MOLPRO
        include "common/cwsave"
        ! tell gradient programs where gradient information is stored
        igrsav = igrsav_
#else
        call stop_all(t_r,'Should be being called from within MOLPRO')
#endif

    end subroutine molpro_set_igrdsav
    
    subroutine molpro_get_reference_info(iWfRecord, iWfSym)

        integer :: iWfSym,iWfRecord
#ifdef MOLPRO
        include "common/code"
        include "common/cref"
        ! wf(1) should contain the record number, as integer, of the last used orbital set.
        iWfRecord = wf(1)
        iWfSym = isyref
#else
        iWfRecord = 21402
#endif

    end subroutine molpro_get_reference_info
    
    function molpro_get_iout() result(iiout)

        integer :: iiout
#ifdef MOLPRO
        include "common/tapes"
        ! output i/o unit for main output (usually 6, but might be different
        ! if in logfile mode, e.g., during geometry optimization)
        iiout = iout
#else
        iiout = 6
#endif

    end function molpro_get_iout

    subroutine molpro_dump_mcscf_dens_for_grad(name,ifil, &
        icore,iclos,iact,nelec,isyref,ms2,iblkq,iseccr,istat1, &
        den1,lden1, &
        den2,lden2, &
        eps,leps, &
        epsc,lepsc,iout,intrel)

        implicit double precision(a-h,o-z)
    
        !  Use the following subroutine to dump the information needed for molpro
        !  forces calculations, such that it is in the exact format that molpro
        !  would have dumped from a MCSCF calculation.  When interfaced with Molpro,
        !  I understand that the molpro equivalent of this routine should be called. 
    
        integer,                            intent(inout) :: name
        integer,                            intent(inout) :: ifil
        integer, dimension(8),              intent(in)    :: icore
        integer, dimension(8),              intent(in)    :: iclos
        integer, dimension(8),              intent(in)    :: iact
        integer,                            intent(in)    :: nelec
        integer,                            intent(in)    :: isyref
        integer,                            intent(in)    :: ms2
        integer,                            intent(in)    :: iblkq
        integer,                            intent(in)    :: iseccr
        integer,                            intent(in)    :: istat1
        integer,                            intent(in)    :: lden1
        real(dp), dimension(lden1), intent(in)    :: den1
        integer,                            intent(in)    :: lden2
        real(dp), dimension(lden2), intent(in)    :: den2
        integer,                            intent(in)    :: leps
        real(dp), dimension(leps),  intent(in)    :: eps
        integer,                            intent(in)    :: lepsc
        real(dp), dimension(lepsc), intent(in)    :: epsc
        integer,                            intent(in)    :: iout
        integer,                            intent(in)    :: intrel
        integer                                           :: igrsav

        !> Write mcscf density in format needed for gradient program.
        !> \param[in,out] name record number to be written
        !> \param[in,out] ifil file number to be written
        !> \param[in] icore numbers of frozen core orbitals in each symmetry
        !> \param[in] iclos plus numbers of closed-shell orbitals in each symmetry
        !> \param[in] iact numbers of active orbitals in each symmetry
        !> \param[in] nelec number of active electrons -- DONE
        !> \param[in] isyref spatial symmetry of wavefunction -- DONE
        !> \param[in] ms2 spin quauntum number times 2 -- DONE
        !> \param[in] iblkq record number * 10 + file number for orbitals (typically 21402) -- DONE
        !> \param[in] iseccr record number * 10 + file number for frozen orbitals (typically 21002) -- DONE
        !> \param[in] istat1 state number (1 for ground state) -- DONE
        !> \param[in] den1 1-particle density matrix
        !> \param[in] lden1 size of den1
        !> \param[in] den2 2-particle density matrix
        !> \param[in] lden2 size of den2
        !> \param[in] eps Lagrangian
        !> \param[in] leps size of leps
        !> \param[in] epsc Lagrangian for frozen core
        !> \param[in] lepsc size of epsc
        !> \param[in] iout unit for output, eg. 6
        !> \param[in] intrel byte size ratio of integer to double precision, normally 1 or 2
        !> \param[out] igrsav where the record was written (needed in common/cwsave)

        integer, dimension(30)      :: header
        integer :: i, ncore
        character(len=6), parameter :: label='MCGRAD'
        integer :: lhead

#ifndef MOLPRO
          call stop_all(t_r,'Should not be here if not running through molpro')
#endif

        ncore = 0
        do i=1,8
            ncore = ncore + icore(i)
            header(1-1+i)=icore(i)
            header(1+7+i)=iclos(i)
            header(1+15+i)=iact(i)
        end do

        header(1+24)=nelec
        header(1+25)=isyref
        header(1+26)=ms2
        header(1+27)=iabs(iblkq)
        header(1+28)=iabs(iseccr)
        header(1+29)=iabs(istat1)
        lhead=30/intrel
        igrsav=10*(name-1)+ifil
      
        name=igrsav/10
        ifil=igrsav-10*name
#ifdef MOLPRO
        call reserv(lhead+lden1+lden2+leps+lepsc,ifil,name,-1)
        call writem(header,lhead,ifil,name,0,label)
        call writem(den1,lden1,ifil,name,lhead,label)
        call writem(den2,lden2,ifil,name,lhead+lden1,label)
        call writem(eps,leps,ifil,name,lhead+lden1+lden2,label)
        if (ncore.ne.0) call writem(epsc,lepsc,ifil,name,&
                      lhead+lden1+lden2+leps,label)
#endif
        write(iout,20) istat1,isyref,name,ifil
!       call setf(2,recmc,1)
      ! ^- cgk: not sure what this does.
20      format(/' Gradient information for state',i2,'.',i1,&
                  ' saved on record  ',i8,'.',i1)
        call molpro_set_igrdsav(igrsav)

        !igrsav=nexfre(igrsav)
        !name=igrsav/10
        !ifil=igrsav-10*name
        !call reserv(lhead+lden1+lden2+leps+lepsc,ifil,name,-1)
!        open (unit = ifil, file = "fciqmc_forces_info", form='UNFORMATTED', access='sequential')
!        write(ifil) header, den1, den2, eps
        !call writem(den2,lden2,ifil,name,lhead+lden1,label)
        !call writem(eps,leps,ifil,name,lhead+lden1+lden2,label)
        !if (ncore.ne.0) call writem(epsc,lepsc,ifil,name,
        !>                          lhead+lden1+lden2+leps,label)
!        write(iout,*) "istat1, isyref, name, ifil", istat1,isyref,name,ifil
!        write(iout,*) "header, den1, den2, eps", header, den1, den2, eps
!20      format(/' Gradient information for state',i2,'.',i1,
       ! >        ' saved on record  ',i8,'.',i1)
       ! return

      end subroutine molpro_dump_mcscf_dens_for_grad

    subroutine BinSearchParts_rdm(iLut,MinInd,MaxInd,PartInd,tSuccess)

        ! Do a binary search in CurrentDets, between the indices of MinInd and
        ! MaxInd. If successful, tSuccess will be true and  PartInd will be a
        ! coincident determinant. If there are multiple values, the chosen one
        ! may be any of them... If failure, then the index will be one less than
        ! the index that the particle would be in if it was present in the list.
        ! (or close enough!)

        integer(kind=n_int) :: iLut(0:NIfTot)
        integer :: MinInd,MaxInd,PartInd
        integer :: i,j,N,Comp
        logical :: tSuccess

        i=MinInd
        j=MaxInd
        if (i-j.eq.0) then
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO)
            if (Comp.eq.0) then
                tSuccess=.true.
                PartInd=MaxInd
                return
            else
                tSuccess=.false.
                PartInd=MinInd
            end if
        end if

        do while(j-i.gt.0)    ! End when the upper and lower bound are the same.
            N = (i+j)/2       ! Find the midpoint of the two indices.

            ! Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is
            ! more or 0 if they are the same
            Comp = DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO)

            if (Comp.eq.0) then
                ! Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                return
            else if ((Comp.eq.1).and.(i.ne.N)) then
                ! The value of the determinant at N is LESS than the determinant
                ! we're looking for. Therefore, move the lower bound of the search
                ! up to N. However, if the lower bound is already equal to N then
                ! the two bounds are consecutive and we have failed...
                i = N
            else if (i .eq. N) then


                if (i .eq. MaxInd-1) then
                    ! This deals with the case where we are interested in the
                    ! final/first entry in the list. Check the final entry of
                    ! the list and leave. We need to check the last index.
                    Comp = DetBitLT(CurrentDets(:,i+1), iLut(:), NIfDBO)
                    if (Comp .eq. 0) then
                        tSuccess = .true.
                        PartInd = i + 1
                        return
                    else if (Comp .eq. 1) then
                        ! final entry is less than the one we want.
                        tSuccess = .false.
                        PartInd = i + 1
                        return
                    else
                        tSuccess=.false.
                        PartInd=i
                        return
                    end if

                else if (i .eq. MinInd) then
                    tSuccess = .false.
                    PartInd = i
                    return
                else
                    i = j
                end if


            else if (Comp .eq. -1) then
                ! The value of the determinant at N is MORE than the determinant
                ! we're looking for. Move the upper bound of the search down to N.
                j = N
            else
                ! We have failed - exit loop.
                i = j
            end if

        end do

        ! If we have failed, then we want to find the index that is one less
        ! than where the particle would have been.
        tSuccess = .false.
        PartInd = max(MinInd,i-1)

    end subroutine BinSearchParts_rdm
    
    subroutine fill_RDM_offdiag_deterministic()

        use bit_rep_data, only: NIfD

        integer :: i, j
        integer :: SingEx(2,1), Ex(2,2)
        real(dp) :: InstSignI, InstSignJ
        real(dp) :: AvSignI, AvSignJ
        logical :: tParity
        integer(n_int) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer :: nI(nel), nJ(nel), IC
        integer :: IterRDM, connect_elem

        ! IterRDM will be the number of iterations that the contributions are
        ! ech weighted by.
        if (mod((iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0) then
            ! This numer of iterations is how regularly the energy is printed
            ! out.
            IterRDM = RDMEnergyIter
        else
            ! This must be the final iteration, as we've got tFill_RDM=.true.
            ! for an iteration where we wouldn't normally need the energy
            IterRDM = mod((Iter+PreviousCycles - IterRDMStart + 1), RDMEnergyIter)
        end if
        
        Ex(:,:) = 0

        ! Cycle over all core dets on this process.
        do i = 1, determ_sizes(iProcIndex)
            iLutI = core_space(:,determ_displs(iProcIndex)+i)
                         
            ! Connections to the HF are added in elsewhere, so skip them here.
            if (DetBitEq(iLutI, iLutHF_True, NifDBO)) cycle
           
            AvSignI = full_determ_vecs_av(1,determ_displs(iProcIndex)+i)

            call decode_bit_det(nI,iLutI)
 
            do j = 1, sparse_core_ham(i)%num_elements-1
                 ! Running over all non-zero off-diag matrix elements
                 ! Connections to whole space (1 row), excluding diagonal elements

                 ! Note: determ_displs holds sum(determ_sizes(0:proc-1))
                 ! Core space holds all the core determinants on every processor,
                 ! so we need to shuffle up to the range of indices corresponding
                 ! to this proc (using determ_displs) and then select the
                 ! correct one, i.
                 
                 iLutJ = core_space(:,core_connections(i)%positions(j))

                 ! Connections to the HF are added in elsewhere, so skip them here.
                 if (DetBitEq(iLutJ, iLutHF_True, NifDBO)) cycle
                 
                 AvSignJ = full_determ_vecs_av(inum_runs,core_connections(i)%positions(j))

                 connect_elem = core_connections(i)%elements(j)

                 IC = abs(connect_elem)

                 if (sign(1, connect_elem) .gt. 0) then
                     tParity = .false.
                 else
                     tParity = .true.
                 end if

                 if (tHPHF) then
                     call decode_bit_det(nJ, iLutJ)

                     call Fill_Spin_Coupled_RDM_v2(iLutI, iLutJ, nI, nJ, AvSignI*IterRDM, AvSignJ,.false.)
                 else
                     if (IC .eq. 1) then
                         ! Single excitation - contributes to 1- and 2-RDM
                         ! (if calculated).
                          
                         ! Note: get_bit_excitmat may be buggy (DetBitOps),
                         ! but will do for now as we need the Ex...
                         call get_bit_excitmat(iLutI(0:NIfD),iLutJ(0:NIfD), SingEx, IC)
                         Ex(:,1) = SingEx(:,1)
                        
                         ! No need to explicitly fill symmetrically as we'll
                         ! generate pairs of determinants both ways around using
                         ! the connectivity matrix.
                         call Fill_Sings_RDM(nI,Ex,tParity,AvSignI*IterRDM,AvSignJ,.false.)

                     else if ((IC .eq. 2) .and. (RDMExcitLevel .ne. 1)) then
                         
                         ! Note: get_bit_excitmat may be buggy (DetBitOps),
                         ! but will do for now as we need the Ex...
                         call get_bit_excitmat(iLutI(0:NIfD), iLutJ(0:NIfD), Ex, IC)
                         call Fill_Doubs_RDM(Ex,tParity,AvSignI*IterRDM,AvSignJ,.false.)
                     end if
                 end if
             end do
        end do

    end subroutine fill_RDM_offdiag_deterministic 

    ! To calculate the dipole moments, the casscf routine has to be called from molpro. i.e.

    ! {rhf;save,2103.2}
    ! {casscf,maxit=0;occ,10,4,4,1;closed,0,0,0,0;iprint,density,civector}
    ! gexpec,dm
    ! {fciqmc,iterations=10000,timestep=0.05,targetwalkers=10000,2RDMonFly,dipoles;core;orbital,2103.2}

    ! Notes: 
    !  o The orbitals used by the fciqmc module have to be the RHF orbitals, not the casscf natural ones.
    !  o The occ directive has to encompass the *whole* space
    !  o If core orbitals want to be frozen in the subsequent fciqmc calclation, then include these orbitals as 'closed'
    !        in the casscf call, and remove the 'core' command from fciqmc
    !  o If the system is too large to do a FCI casscf (which will hopefully normally be the case), then you can restrict
    !        the space by using the 'restict' directive directly after the 'wf' directive in the casscf, i.e.

    ! memory,128,m
    ! geometry={C;O,C,r};r=2.1316 bohr
    ! basis,VDZ;
    ! {rhf;save,2103.2}
    ! {casscf,maxit=0;occ,14,6,6,2;frozen,0,0,0,0;closed,2,0,0,0;
    ! wf,14,1,0;
    ! restrict,2,2,3.1,4.1;       !This defines a set of orbitals which must be doubly occupied in all configurations to cut down the
    ! restrict,2,2,1.2,1.3;       !size of the space, but ensure that the integrals are still calculated over the whole active space
    ! restrict,0,0,10.1,11.1,12.1,13.1,14.1;  !These are orbitals which must remain unoccupied to further 
                                              !reduce the size of the space
    ! restrict,0,0,3.2,4.2,5.2,6.2;
    ! restrict,0,0,3.3,4.3,5.3,6.3;
    ! iprint,density,civector}
    ! gexpec,dm
    ! {fciqmc,ITERATIONS=20,MAXATREF=50000,targetWALKERS=10000000;
    !   orbital,2103.2 }         
         !Note no 'core' directive included, since the core orbitals are 'closed' in the casscf and so will be ignored.

    subroutine CalcDipoles(Norm_1RDM)

#ifdef MOLPRO
        use outputResult
        use SymData, only: Sym_Psi,nSymLabels
        use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, ClassCountInd
        use SymExcitDataMod, only: SpinOrbSymLabel,SymLabelCounts2

        integer, dimension(nSymLabels) :: elements_assigned1, blockstart1
        real(dp) :: dipmom(3),znuc(3),zcor(3)
        real(dp), allocatable :: SymmetryPacked1RDM(:),zints(:,:)
        integer :: i, j, ipr, Sym_i, Sym_j, Sym_ij,
        integer :: posn1, isize, isyref, mxv, iout, ierr
        integer :: nt_frz(8), ntd_frz(8)
#endif
        ! The only thing needed is the 1RDM (normalized)
        real(dp), intent(in) :: Norm_1RDM
        character(len=*), parameter :: t_r='CalcDipoles'

#ifdef MOLPRO

        if (iProcIndex.eq.0) then
            iout = molpro_get_iout()
            ! We need to work out
            ! a) how molpro symmetry-packs UHF integrals (ROHF would be fine though)
            ! b) Ensure that the 1RDM is correctly calculated for UHF (It is always allocated as spatorbs)
            ! c) Modify this routine for contracting over spin-orbitals
            if (tOpenShell) call stop_all(t_r,'Not working for ROHF/UHF')

            isyref=Sym_Psi+1 ! Spatial symmetry of the wavefunction.

            ! Size of symmetry packed arrays (spatial).
            isize = 0
            blockstart1(:) = 0
            do i = 0,nSymLabels-1
                ! Find position of each symmetry block in sym-packed forms of RDMS 1 & 2.
                blockstart1(i+1)=isize+1 ! N.B. Len_1RDM still being updated in this loop.

                isize = isize + (SymLabelCounts2(2,ClassCountInd(1,i,0))*   &
                    (SymLabelCounts2(2,ClassCountInd(1,i,0))+1))/2 ! Counting alpha orbitals.
            end do

            nt_frz(:) = 0
            ntd_frz(:) = 0
            do i = 0,nSymLabels-1
                nt_frz(i+1) = SymLabelCounts2(2,ClassCountInd(1,i,0))
            end do

            do i = 2,8
                ntd_frz(i) = ntd_frz(i-1) + (nt_frz(i-1)*(nt_frz(i-1)+1))/2
            end do

            elements_assigned1(:) = 0
            allocate(SymmetryPacked1RDM(isize))
            SymmetryPacked1RDM(:) = 0.0_dp
            do i = 1,SpatOrbs ! Run over spatial orbitals, ALL ELECTRON ONLY.
                do j = 1, i ! i .ge. j
                    Sym_i = SpinOrbSymLabel(2*i)  ! Consider only alpha orbitals.
                    Sym_j = SpinOrbSymLabel(2*j)
                    Sym_ij = RandExcitSymLabelProd(Sym_i, Sym_j)
                    if (Sym_ij .eq. 0) then
                        if ((Sym_i+1) .gt. nSymLabels) call stop_all(t_r,'Error')
                        posn1 = blockstart1(Sym_i+1) + elements_assigned1(Sym_i+1)
                        if ((posn1 .gt. isize) .or. (posn1 .lt. 1)) then
                            call stop_all(t_r,'Error filling rdm')
                        end if

                        SymmetryPacked1RDM(posn1)=NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM
                        if (i .ne. j) then
                            ! Double the off-diagonal elements of the 1RDM, so
                            ! that when we contract over the symmetry packed
                            ! representation of the 1RDM, it is as if we are
                            ! also including the other half of the matrix
                            SymmetryPacked1RDM(posn1) = 2.0_dp*SymmetryPacked1RDM(posn1)
                        end if
                        elements_assigned1(Sym_i+1) = elements_assigned1(Sym_i+1) + 1
                    end if

                end do
            end do
                
            write(6,*) "Size of symmetry packed array: ", isize
            write(6,*) "Symmetry packed 1RDM: ", SymmetryPacked1RDM(:)
            dipmom(:) = 0.0_dp

            call clearvar('DMX')
            call clearvar('DMY')
            call clearvar('DMZ')

            allocate(zints(isize,3),stat=ierr)
            if (ierr .ne. 0) then
                write(6,*) "Alloc failed: ",ierr
                call stop_all(t_r,'Alloc failed')
            end if
            ! This now goes through an F77 wrapper file so that we can access the
            ! common blocks and check that the size of the symmetry packed arrays
            ! is correct.
            call GetDipMomInts(zints,isize,znuc,zcor,nt_frz,ntd_frz)

            do ipr=1,3
            
!                call pget(zints,ipr,znuc,zcor)

                ! Now, contract.
                do i = 1,isize
                    dipmom(ipr) = dipmom(ipr) - zints(i,ipr)*SymmetryPacked1RDM(i)
                end do
                dipmom(ipr) = dipmom(ipr) + znuc(ipr) - zcor(ipr)
            end do

            write(iout,"(A)") ""
            write(iout,"(A,3f15.8)") "DIPOLE MOMENT: ",dipmom(1:3)
            write(iout,"(A)") ""
            call output_result('FCIQMC','Dipole moment',dipmom(1:3),1,isyref,numberformat='3f15.8',debye=.TRUE.)
            mxv=1
            call setvar('DMX',dipmom(1),'AU',1,1,mxv,-1)
            call setvar('DMY',dipmom(2),'AU',1,1,mxv,-1)
            call setvar('DMZ',dipmom(3),'AU',1,1,mxv,-1)
            deallocate(zints,SymmetryPacked1RDM)
        end if

#else
        call warning_neci(t_r,'Cannot compute dipole moments if not running within molpro. Exiting...')
#endif

    end subroutine CalcDipoles

end module rdms
