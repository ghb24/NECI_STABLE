#include "macros.h"

! this is the initialization module for the real-time FCIQMC implementation

module real_time_init

    use real_time_data, only: t_real_time_fciqmc, gf_type, real_time_info, &
                              t_complex_ints, gf_overlap, temp_det_list, dyn_norm_red, &
                              temp_det_pointer, temp_det_hash, temp_freeslot, tOverpopulate, &
                              pert_norm, second_spawn_iter_data, DiagParts, stepsAlpha, &
                              DiagVec, normsize, valid_diag_spawns, tStabilizerShift, &
                              NoatHF_1, Annihilated_1, Acceptances_1, NoBorn_1, spawnBuf, &
                              SpawnFromSing_1, NoDied_1, NoAborted_1, NoRemoved_1, &
                              NoAddedInitiators_1, NoInitDets_1, NoNonInitDets_1, &
                              NoInitWalk_1, NoNonInitWalk_1, InitRemoved_1, tDynamicAlpha, &
                              AllNoatHF_1, AllNoatHF_1, AllGrowRate_1, AllGrowRateAbort_1, &
                              AllNoBorn_1, AllSpawnFromSing_1, AllNoDied_1, gf_count, &
                              AllAnnihilated_1, AllNoAborted_1, AllNoRemoved_1, allPopSnapshot, &
                              AllNoAddedInitiators_1, AllNoInitDets_1, AllNoNonInitDets_1, &
                              AllNoInitWalk_1, AllNoNonInitWalk_1, AllInitRemoved_1, &
                              AccRat_1, AllNoatDoubs_1, AllSumWalkersCyc_1, current_overlap, &
                              TotPartsStorage,  t_rotated_time, TotPartsPeak, asymptoticShift, &
                              tau_imag, tau_real, elapsedRealTime, elapsedImagTime, tNewOverlap, &
                              TotWalkers_orig, dyn_norm_psi, gs_energy, shift_damping, &
                              tStaticShift, MaxSpawnedDiag, tDynamicCoreSpace, overlap_states, &
                              overlap_real, overlap_imag, allGfs, tRealTimePopsfile, &
                              tLimitShift, nspawnMax, shiftLimit, numCycShiftExcess, &
                              TotPartsLastAlpha, alphaDamping, tDynamicDamping, iterInit, &
                              etaDamping, tStartVariation, rotThresh, stabilizerThresh, &
                              tInfInit,  popSnapshot, snapshotOrbs, phase_factors, tVerletSweep, &
			      numSnapshotOrbs, tLowerThreshold, t_kspace_operators, tVerletScheme, &
                              tLogTrajectory, tReadTrajectory, alphaCache, tauCache, trajFile, &
                              tGenerateCoreSpace, tGZero, wn_threshold, corespace_log_interval
    use real_time_procs, only: create_perturbed_ground, setup_temp_det_list, &
                               calc_norm, clean_overlap_states, openTauContourFile
    use verlet_aux, only: backup_initial_state, setup_delta_psi
    use constants, only: dp, n_int, int64, lenof_sign, inum_runs
    use Parallel_neci
    use ParallelHelper, only: iProcIndex, root, MPIbarrier, nNodes, MPI_SUM
    use util_mod, only: get_unique_filename, get_free_unit
    use Logging, only: tIncrementPops
    use kp_fciqmc_data_mod, only: tMultiplePopStart, tScalePopulation, &
                                  tOverlapPert, overlap_pert, scaling_factor
    use CalcData, only: tChangeProjEDet, tReadPops, tRestartHighPop, tFCIMC, tJumpShift, &
                        tStartSinglePart, tau, nmcyc, iPopsFileNoRead, tWritePopsNorm, &
                        tWalkContGrow, diagSft, pops_norm, InitWalkers, MemoryFacSpawn, &
                        StepsSft, tSemiStochastic, tTruncInitiator, tAddToInitiator
    use FciMCData, only: tSearchTau, alloc_popsfile_dets, pops_pert, tPopsAlreadyRead, &
                         tSinglePartPhase, iter_data_fciqmc, iter, PreviousCycles, &
                         AllGrowRate, spawn_ht, pDoubles, pSingles, TotParts, &
                         MaxSpawned, tSearchTauOption, TotWalkers, SumWalkersCyc, &
                         CurrentDets, popsfile_dets, MaxWalkersPart, WalkVecDets, &
                         SpawnedParts, determ_sizes
    use SystemData, only: lms, G1, nBasisMax, tHub, nel, tComplexWalkers_RealInts, nBasis, tReal
    use SymExcitDataMod, only: kTotal
    use sym_mod, only: MomPbcSym
    use perturbations, only: init_perturbation_annihilation, &
                             init_perturbation_creation
    use fcimc_initialisation, only: SetupParameters, InitFCIMCCalcPar, &
                                    init_fcimc_fn_pointers 
    use LoggingData, only: tZeroProjE, tFCIMCStats2
    use fcimc_output, only: write_fcimcstats2, WriteFciMCStatsHeader
    use replica_data, only: allocate_iter_data, set_initial_global_data
    use bit_rep_data, only: nifbcast, niftot, extract_sign, nifdbo
    use bit_reps, only: decode_bit_det
    use adi_references, only: setup_reference_space

    implicit none

    real(dp) :: benchmarkEnergy

contains

    subroutine init_real_time_calc_single()
      use real_time_procs, only: get_tot_parts
        ! this routine takes care of the correct setup of the real-time 
        ! calculation. like reading the popsfiles and preparing the start 
        ! of the calculation and setting certain global variables
        implicit none
        character(*), parameter :: this_routine = "init_real_time_calc_single"

        integer :: ierr

        print *, " Entering real-time FCIQMC initialisation "

        ! think about what variables have to be set for a succesful calc.

        ! also call the "normal" NECI setup routines to allow calculation
        call SetupParameters()

        ! for real-space hubbard, we conver to momentum operators if desired
        if(t_kspace_operators) call setup_momentum_operators()

        ! have to think about the the order about the above setup routines! 
        ! within this Init a readpops is called
        ! this function already produces the correctly perturbed ground state
        call InitFCIMCCalcPar()
        ! also init pointer here, and think about what options and defaults 
        ! should be set for a succsesfull init
        call init_fcimc_fn_pointers()

        ! then call the setup routine, which set all remaining needed quantities
        call setup_real_time_fciqmc()

        ! definetly read-in stored popsfile here. 
        ! need to store both <y(0)| and also create a_j y(0)> during read-in!
!         call read_popsfile_real_time()
        ! actually the InitFCIMCCalcPar should do that now correctly already

        ! do an MPIbarrier here.. although don't quite know why        
        call MPIBarrier(ierr)

        if(.not. tReadPops) call set_initial_global_data(TotWalkers, CurrentDets)
        
    end subroutine init_real_time_calc_single

    subroutine setup_real_time_fciqmc()
        ! this is the last setup routine, which depending on compilation,
        ! number of copies etc. sets up the final needed quantities to run 
        ! a simulation
      implicit none
        character(*), parameter :: this_routine = "setup_real_time_fciqmc"
        integer :: ierr,  run


        ! the new total momentum has to be constructed before the 
        ! time-evolved state is read in, as the latter deletes the 
        ! pops_pert, because perturbation and read-in are done in one
        ! function (dependencies...)
        if(tHub) then
           if(allocated(pops_pert)) then

              if(pops_pert(1)%nannihilate == 1) kTotal = kTotal &
                   - G1(pops_pert(1)%ann_orbs(1))%k
              if(pops_pert(1)%ncreate == 1) kTotal = kTotal &
                   + G1(pops_pert(1)%crtn_orbs(1))%k
              call MomPbcSym(kTotal,nBasisMax)
              print *, "New total momentum", kTotal
           endif
        endif

        ! allocate the according quantities! 
        ! n_time_steps have to be set here!
        print *, " Allocating greensfunction and wavefunction norm arrays!"
        ! allocate an additional slot for initial values
        if(numSnapshotOrbs>0) then 
           allocate(popSnapshot(numSnapshotOrbs),stat=ierr)
           allocate(allPopSnapshot(numSnapshotOrbs),stat=ierr)
           popSnapshot = 0.0_dp
           allPopSnapshot = 0.0_dp
        else
           allocate(popSnapshot(1),stat=ierr)
           allocate(allPopSnapshot(numSnapshotOrbs),stat=ierr)
           popSnapshot = 0.0_dp
           allPopSnapshot = 0.0_dp
        endif
        allocate(gs_energy(inum_runs),stat = ierr)
        allocate(temp_freeslot(MaxWalkersPart),stat = ierr)
        allocate(TotPartsPeak(inum_runs),stat = ierr)
        allocate(numCycShiftExcess(inum_runs), stat = ierr)
        numCycShiftExcess = 0
        ! allocate spawn buffer for verlet scheme
        if(tVerletScheme) allocate(spawnBuf(0:niftot,1:maxSpawned))

        TotPartsPeak = 0.0_dp
        gs_energy = benchmarkEnergy

        call init_overlap_buffers()

        if(tRealTimePopsfile) call readTimeEvolvedState()
        
        ! check for set lms.. i think that does not quite work yet 
        print *, "mz spin projection: ", lms

        print *, "tSinglePartPhase?:",tSinglePartPhase
        print *, "tWalkContGrow?", tWalkContGrow
        print *, "diagSft:", diagSft

        ! intialize the 2nd temporary determinant list needed in the 
        ! real-time fciqmc 

        ! also maybe use the spawn_ht hash table, so allocated it here! 
        call setup_temp_det_list()

        print *, "allocated(temp_det_list)?", allocated(temp_det_list)
        print *, "associated(temp_det_pointer)?", associated(temp_det_pointer)
        print *, "associated(temp_det_hash)?", associated(temp_det_hash)

        print *, "associated(spawn_ht)?", associated(spawn_ht)

        print *, "Allgrowrate: ", AllGrowRate
        ! print out the first infos on the calculation.. 
        ! although that definetly has to be changed for the real-time fciqm

        ! use new output format! 
        tFCIMCStats2 = .true.

        if (tFCIMCStats2) then
           call write_fcimcstats2(iter_data_fciqmc, initial = .true.)
        else
            call WriteFciMCStatsHeader()
        end if

        ! set the iter variable to 0 probably
        iter = 0
        
        ! and also the PreviousCycles var. since its essentially regarded as 
        ! a new calulcation
        PreviousCycles = 0

        ! for intermediate test_purposes turn off spawning to check if the 
        ! diagonal step works as intented 
!         pSingles = 0.0_dp
!         pDoubles = 0.0_dp

        ! also initialize the second_spawn_iter_data type
        call allocate_iter_data(second_spawn_iter_data)

        ! and also initialize the values: 
        second_spawn_iter_data%ndied = 0.0_dp
        second_spawn_iter_data%nborn = 0.0_dp
        second_spawn_iter_data%nannihil = 0.0_dp
        second_spawn_iter_data%naborted = 0.0_dp
        second_spawn_iter_data%nremoved = 0.0_dp
        second_spawn_iter_data%update_growth = 0.0_dp
        second_spawn_iter_data%update_growth_tot = 0.0_dp
        second_spawn_iter_data%tot_parts_old = TotParts
        second_spawn_iter_data%update_iters = 0

        TotPartsStorage = TotParts
        TotPartsLastAlpha = TotParts

        ! also intitialize the 2nd spawning array to deal with the 
        ! diagonal death step in the 2nd rt-fciqmc loop
        allocate(DiagVec(0:nifbcast, MaxWalkersPart), stat = ierr)

        DiagVec = 0

        DiagParts => DiagVec

        ! and the initial_spawn_slots equivalent
        ! although i think i can reuse the initialSpawnedSlots..
!         allocate(initial_diag_spawn_list(0:nNodes-1), stat = ierr) 

        valid_diag_spawns = 1

        do run = 1, inum_runs
           SumWalkersCyc(run) = SumWalkersCyc(run) + &
                sum(TotParts(min_part_type(run):max_part_type(run)))
        enddo

        ! also initialize all the relevant first RK step quantities.. 
        NoatHF_1 = 0.0_dp 
        Annihilated_1 = 0.0_dp
        Acceptances_1 = 0.0_dp
        NoBorn_1 = 0.0_dp
        SpawnFromSing_1 = 0 
        NoDied_1 = 0.0_dp
        NoAborted_1 = 0.0_dp
        NoRemoved_1 = 0.0_dp
        NoAddedInitiators_1 = 0
        NoInitDets_1 = 0
        NoNonInitDets_1 = 0
        NoInitWalk_1 = 0.0_dp
        NoNonInitWalk_1 = 0.0_dp
        InitRemoved_1 = 0
        
        ! also the global variables 
        AllNoatHF_1 = 0.0_dp
        AllNoatDoubs_1 = 0.0_dp
        AllGrowRate_1 = 0.0_dp
        AllGrowRateAbort_1 = 0
        AllNoBorn_1 = 0.0_dp
        AllSpawnFromSing_1 = 0
        AllNoDied_1 = 0.0_dp
        AllAnnihilated_1 = 0.0_dp
        AllNoAborted_1 = 0.0_dp
        AllNoRemoved_1 = 0.0_dp
        AllNoAddedInitiators_1 = 0
        AllNoInitDets_1 = 0
        AllNoNonInitDets_1 = 0
        AllNoInitWalk_1 = 0.0_dp
        AllNoNonInitWalk_1 = 0.0_dp
        AllInitRemoved_1 = 0
        AccRat_1 = 0.0_dp
        AllSumWalkersCyc_1 = 0.0_dp

        tVerletSweep = .false.
        if(tVerletScheme) then 
           call setup_delta_psi()
           call backup_initial_state()
           tau = tau/iterInit
        endif

        if(tGenerateCoreSpace) call initialize_corespace_construction()      

        if(tReadTrajectory) call read_in_trajectory()
        
        ! Set up the reference space for the adi-approach
        call setup_reference_space(tReadPops)

        call rotate_time()           

    end subroutine setup_real_time_fciqmc

    subroutine rotate_time()
      implicit none
      ! to avoid code multiplication
      if(t_rotated_time) then
         tau_imag = - sin(real_time_info%time_angle)*tau
         tau_real = cos(real_time_info%time_angle)*tau
      else
         tau_imag = 0.0_dp
         tau_real = tau
      endif
    end subroutine rotate_time

    subroutine init_overlap_buffers
      use CalcData, only: tSemiStochastic, ss_space_in
      use semi_stoch_gen, only: init_semi_stochastic
      use real_time_procs, only: reset_tot_parts
      implicit none
      ! this subroutine sets up everything required to compute green's functions
      integer :: ierr, j, i
      complex(dp), allocatable :: norm_buf(:)

      normsize = inum_runs**2      
      allocate(overlap_real(gf_count),overlap_imag(gf_count))
      allocate(gf_overlap(normsize,gf_count), stat = ierr)
      allocate(pert_norm(normsize,gf_count),stat = ierr)
      allocate(dyn_norm_psi(normsize),stat = ierr)
      dyn_norm_psi = 1.0_dp
      allocate(dyn_norm_red(normsize,gf_count), stat = ierr)
      allocate(current_overlap(normsize,gf_count),stat = ierr)

      gf_overlap = 0.0_dp

      ! also need to create the perturbed ground state to calculate the 
      ! overlaps to |y(t)> 
      call create_perturbed_ground()
      if(tSemiStochastic) call init_semi_stochastic(ss_space_in)
      ! If only the corespace time-evolution is to be taken, truncate the
      ! initial wavefunction to the corespace
      ! We currently do not truncate the overlap state too, but it might come
      ! later
      if(tGZero) then
         call truncate_initial_state()
         call truncate_overlap_states()
         call reset_tot_parts()
         TotWalkers = determ_sizes(iProcIndex)
      end if

      ! if a ground-state POPSFILE is used, we need to ensure coherence between the replicas
      if(.not. tRealTimePopsfile) call equalize_initial_phase()

      allocate(norm_buf(normsize),stat=ierr)
      ! to avoid dividing by 0 if not all entries get filled
      pert_norm = 1.0_dp
      norm_buf = calc_norm(CurrentDets,int(TotWalkers))
      call MPIReduce(norm_buf,MPI_SUM,dyn_norm_psi)
      do j = 1,gf_count
         ! calc. the norm of the perturbed ground states
         norm_buf = calc_norm(overlap_states(j)%dets,overlap_states(j)%nDets)
         print *, "CHECKSUM", TotWalkers, overlap_states(j)%nDets
         ! the norm (squared) can be obtained by reduction over all processes
         call MPIReduce(norm_buf,MPI_SUM,pert_norm(:,j))
         ! for diagonal green's functions, this is the same as pert_norm, but 
         ! in general, this general normalization is required.
         do i = 1, normsize
            dyn_norm_red(i,j) = sqrt(pert_norm(i,j)*dyn_norm_psi(i))
         end do
      enddo

      deallocate(norm_buf)

      if(.not. allocated(shift_damping)) then 
         allocate(shift_damping(inum_runs), stat = ierr)
         shift_damping = 0.0_dp
      endif

      if(tLogTrajectory) call openTauContourFile()

    end subroutine init_overlap_buffers

    ! write a routine which sets the defauls values, and maybe even turns off
    ! certain otherwise non-compatible variable to the appropriate values
    ! problem: this is called after most of the input is read
    subroutine set_real_time_defaults()
        implicit none

        ! by default, use runge-kutta instead of verlet
        tVerletScheme = .false.
        iterInit = 1

        ! todo: figure out quantities
        ! maximum number of spawn attempts per determinant if the number of
        ! total spawns exceeds some threshold
        nspawnMax = 1000

        ! energy offset
        benchmarkEnergy = 0.0_dp

        ! for the start definetly not change tau
        tSearchTau = .true.

        ! also set readpops to get the <y(0)| reference from the "normal"
        ! neci init routines
        tReadPops = .true.
        ! By default, apply a seperate perturbation to the overlap_states
        tNewOverlap = .true.

        ! startsinglepart does not work with popsfile and is not wanted too
        tStartSinglePart = .false.
        ! also, the shift is something infamous in real-time, so definitely no
        ! jumpshift
        tJumpShift = .false.

        ! but to ensure that the shift does not vary anymore, since there is 
        ! no such concept as the varying shift in the real-time fciqmc
        ! exception: when using rotated times, the shift still has to be considered
        tWalkContGrow = .true.

        ! probably not change reference anymore.. but check
        tChangeProjEDet = .true.

        ! and dont only restart on the highly populated dets
        tRestartHighPop = .false.

        ! nick also has some scale_population routine, think about 
        ! implementing this also for the real-time restart 
        tScalePopulation = .false.
        scaling_factor = 1.0_dp
        
        ! nick also has some multiple popsfile start routines..
        ! set the multiple popsstart with the number of replicas of mneci
        ! provided
        tMultiplePopStart = .false.
        
        ! by default, do not output any population of orbitals
        numSnapshotOrbs = 0

        ! from my way of outputting popsfiles i always do it in popsfile.n
        ! format -> so i probably have to set tIncrementPops and the count
        tIncrementPops = .true.
        iPopsFileNoRead = -1
        ! this will always result in *.0 being chosen as name, was rethought and
        ! decided to be good - this way, no files will be overwritten and 
        ! the read-in file is always the same

        ! have to set, that popsfile is not yet read in:
        tPopsAlreadyRead = .false.

        ! overwrite tfcimc flag to not enter the regular fcimc routine 
        tFCIMC = .false.

        ! usually only real-valued FCIDUMPs
        t_complex_ints = .false.

        ! probably should zero the projected energy, since its a total 
        ! different system 
        tZeroProjE = .false.

        ! setup_rotated_time: by default, pure real time is used
        t_rotated_time = .false.
        ! usually, the walker number shall be controlled
        tStaticShift = .false.
        asymptoticShift = 0.0_dp

        ! and in case of semistochastic approach, the core space shall be static
        tDynamicCoreSpace = .false.

        ! usually, systems with real integrals will be considered, but the walkers will
        ! always be complex
        tComplexWalkers_RealInts = .true.
        
        ! by default, the perturbation operators are in the same basis as the wavefunction
        t_kspace_operators = .false.

        ! if no gf kind is specified, only the overlap with the initial state will
        ! be considered -> only one overlap is obtained
        gf_count = 1
        allGfs = 0
        
        ! if starting a new calculation, we start at time 0 (arbitrary)
        elapsedRealTime = 0.0_dp
        elapsedImagTime = 0.0_dp

        ! by default, the initial state is taken from an ordinary popsfile
        ! if a time evolved state is desired, a second popsfile has to be supplied
        tRealTimePopsfile = .false.
        tStabilizerShift = .false.
        ! threshold population (in relation to the peak walker number) for activation of
        ! stabilizer shift
        stabilizerThresh = 0.8
        ! the merging of spawning events is done entirely automatically and therfore can not
        ! be switched on manually
        tLimitShift = .false.
        shiftLimit = 0.7_dp

        ! default values for dynamic rotation angle updating (it is not enabled by default)
        tDynamicAlpha = .false.
        stepsAlpha = 10
        alphaDamping = 0.05
        rotThresh = 0
	tLowerThreshold = .false.
        tOverpopulate = .false.
        tStartVariation = .false.

        ! the damping is constant by default
        tDynamicDamping = .false.
        etaDamping = 0.01

        ! and we allow normal usage of inititators
        tInfInit = .false.

        ! and we do not print out the contour separately
        tLogTrajectory = .false.
        tReadTrajectory = .false.

        tGenerateCoreSpace = .false.
        wn_threshold = 0.01
        corespace_log_interval = 300
        ! Get the full Green's function, not only in the corespace
        tGZero = .false.
    end subroutine set_real_time_defaults

    ! need a specific popsfile read function for the real-time calculation
    ! based upon georges previous implementation, but changed to handle 
    ! new code base
    ! need a routine which prepares the converged groundstates from an
    ! imaginary-time FCIQMC calculation 

    ! use code provided in the NECI GreensFuncs branch by George Booth
    ! in general reuse as much of the already provided functionality! 

    ! need a setup routine in the regular imag-time neci routine, which prints
    ! out the specified amount of GS wavefunctions
    ! and which sets the necessary flags and stuff

    ! i want through the CHANGEVARS facility start the print out of 
    ! the input specified print out of POPSFILES between certain intervals

    subroutine read_popsfile_real_time()
        use PopsfileMod, only: open_pops_head, FindPopsfileVersion, ReadPopsHeadv4
        implicit none

        integer :: iunit, popsversion, iPopLenof_Sign, iPopNel, iPopIter, &
                   PopNIfD, PopNIfY, PopNIfSgn, PopNIfFlag, PopNIfTot, &
                   PopBlockingIter, Popinum_runs, PopRandomHash(2056), &
                   read_nnodes, PopBalanceBlocks
        logical :: formpops, binpops, tPop64Bit, tPopHPHF, tPopLz
        integer(int64) :: iPopAllTotWalkers, read_walkers_on_nodes(0:nProcessors-1)
        real(dp) :: PopDiagSft(inum_runs), read_tau, PopSumNoatHF(lenof_sign), &
                    read_psingles, read_pparallel
        HElement_t(dp) :: PopAllSumENum(inum_runs)
        character(255) :: popsfile

        character(*), parameter :: this_routine = "read_popsfile_real_time"

        call open_pops_head(iunit,formpops,binpops)

        if (formpops) then
            ! this is the NON-binary read
            popsversion = FindPopsfileVersion(iunit)
            if (popsversion /= 4) then
                call stop_all(this_routine, "wrong POPSFILE version /= 4!")
            end if
            call ReadPopsHeadv4(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs, PopNIfFlag,PopNIfTot, &
                read_tau,PopBlockingIter, PopRandomHash, read_psingles, &
                read_pparallel, read_nnodes, read_walkers_on_nodes, PopBalanceBlocks)

        else 
            ! if popsfiles are stored in binary! there are seperate files for 
            ! the header and the actual population stats
            if (iProcIndex == root) then
                close(iunit)
                call get_unique_filename('POPSFILEBIN', tIncrementPops,.false.,&
                        iPopsFileNoRead,popsfile)
                open(iunit, file = popsfile, status = 'old', form = 'unformatted')
            end if
        end if

    end subroutine read_popsfile_real_time

    subroutine readTimeEvolvedState()
      use PopsfileMod, only : FindPopsfileVersion, ReadPopsHeadv4, InitFCIMC_pops, &
           open_pops_head
      use semi_stoch_gen, only: init_semi_stochastic
      use semi_stoch_procs, only: end_semistoch
      use real_time_procs, only: reset_core_space
      use CalcData, only: ss_space_in
      use FciMCData, only : TotImagTime
      implicit none
      
        integer :: iunit, popsversion, iPopLenof_Sign, iPopNel, iPopIter, &
                   PopNIfD, PopNIfY, PopNIfSgn, PopNIfFlag, PopNIfTot, &
                   PopBlockingIter, Popinum_runs, PopRandomHash(2056), &
                   read_nnodes, PopBalanceBlocks
        logical :: formpops, binpops, tPop64Bit, tPopHPHF, tPopLz
        integer(int64) :: iPopAllTotWalkers, read_walkers_on_nodes(0:nProcessors-1)
        real(dp) :: PopDiagSft(inum_runs), read_tau, PopSumNoatHF(lenof_sign), &
                    read_psingles, read_pparallel
        HElement_t(dp) :: PopAllSumENum(inum_runs)
        integer :: ierr

        character(255) :: rtPOPSFILE_name
        character(*), parameter :: this_routine = "readTimeEvolvedState"

        if(tSemiStochastic) then
           ! if semi-stochastic mode is enabled, it has to be disabled for read-in again
           ! as load balancing has to be performed
           call end_semistoch()
           call reset_core_space()
        endif

        binpops = .false.
        
        rtPOPSFILE_name = 'TIME_EVOLVED_POP'
  
        ! get the file containing the time evolved state
        call open_pops_head(iunit, formpops, binpops, rtPOPSFILE_name)
       
        popsversion = FindPopsfileVersion(iunit)
        if(popsversion /= 4) call stop_all(this_routine, "wrong popsfile version of TIME_EVOLVED_POP")
       
        call ReadPopsHeadv4(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
             iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
             PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs, PopNIfFlag,PopNIfTot, &
             read_tau,PopBlockingIter, PopRandomHash, read_psingles, &
             read_pparallel, read_nnodes, read_walkers_on_nodes, PopBalanceBlocks)

        ! at this point, we do not want to perturb the state and have no use for the 
        ! pops_pert variable anymore -> deallocate it

        ! read in the hacked shift_damping
        if(.not. allocated(shift_damping)) allocate(shift_damping(inum_runs),stat=ierr)
        shift_damping = PopSumNoatHF(1:inum_runs)

        call clear_pops_pert(pops_pert)
        
        ! read in the time evolved state and use it as initial state
        call InitFCIMC_pops(iPopAllTotWalkers, PopNIfSgn, iPopNel, read_nnodes, &
             read_walkers_on_nodes, pops_pert, &
             PopBalanceBLocks, PopDiagSft, rtPOPSFILE_name)

        call set_initial_times(read_tau, TotImagTime)

        ! if we disabled semi-stochastic mode temporarily, reenable it now
        if(tSemiStochastic) call init_semi_stochastic(ss_space_in)
      
    end subroutine readTimeEvolvedState

    subroutine set_initial_times(real_time, imag_time)
      implicit none
      real(dp), intent(in) :: real_time, imag_time
      
      ! with the inclusion of dynamic alpha, both, the real and the imaginary part
      ! have to be stored in the time-evolved popsfile
      ! the 
      elapsedImagTime = imag_time
      elapsedRealTime = real_time
    end subroutine set_initial_times

    subroutine equalize_initial_phase()
      use FciMCData, only: AllSumNoatHF
      use bit_rep_data, only: extract_sign
      use bit_reps, only: encode_sign
      implicit none
      
      integer :: i, signs(lenof_sign), iGf
      real(dp) :: tmp_sgn(lenof_sign)
      logical :: tSuccess

      signs = 1
      do i = 1, lenof_sign
         if(AllSumNoatHF(i)/AllSumNoatHF(1) < 0) then
            signs(i) = -1
         else
            signs(i) = 1
         endif
      enddo
      if(any(signs<0)) then
         do i = 1, TotWalkers
            call extract_sign(CurrentDets(:,i),tmp_sgn)
            tmp_sgn = tmp_sgn * signs
            call encode_sign(CurrentDets(:,i),tmp_sgn)
         end do 
         do iGf = 1, gf_count
            do i = 1, overlap_states(iGf)%nDets
               call extract_sign(overlap_states(iGf)%dets(:,i),tmp_sgn)
               tmp_sgn = tmp_sgn * signs
               call encode_sign(overlap_states(iGf)%dets(:,i),tmp_sgn)
            end do
         end do
      endif
         
    end subroutine equalize_initial_phase


    ! here, we convert the real-space annihilation/creation operators to momentum space
    subroutine setup_momentum_operators()
      use SystemData, only: G1
      implicit none
      
      integer :: k(3),state, momentum_state,spin, i, j
    
      if(gf_type == -1) then
         momentum_state = pops_pert(1)%ann_orbs(1)
      else if(gf_type == 1) then
         momentum_state = pops_pert(1)%crtn_orbs(1)
      else
         call stop_all("Equalize initial phase","Invalid momentum space green's function")
      endif
      allocate(phase_factors(nBasis/2))
      k = G1(momentum_state)%k
      spin = G1(momentum_state)%ms
      call clear_pops_pert(pops_pert)
      call clear_pops_pert(overlap_pert)
      allocate(pops_pert(nBasis/2))
      allocate(overlap_pert(nBasis/2))
      i = 0
      do state = 1, nBasis
         if(G1(state)%ms /= spin) cycle
         i = i + 1
         phase_factors(i) = 0
         do j = 1,3
            phase_factors(i) = phase_factors(i) + G1(state)%k(j) * k(j)
         end do
         call setup_single_perturbation(i,state,gf_type)
      end do
    end subroutine setup_momentum_operators

    subroutine setup_single_perturbation(pos,orbital,type)
      implicit none

      integer, intent(in) :: pos, orbital, type
      
      if(type == 1) then
         allocate(pops_pert(pos)%crtn_orbs(1)) 
         allocate(overlap_pert(pos)%crtn_orbs(1)) 
         overlap_pert(pos)%ncreate = 1
         pops_pert(pos)%ncreate = 1
         pops_pert(pos)%crtn_orbs(1) = orbital
         overlap_pert(pos)%crtn_orbs(1) = orbital
         call init_perturbation_creation(overlap_pert(pos))
         call init_perturbation_creation(pops_pert(pos))
      else
         allocate(pops_pert(pos)%ann_orbs(1)) 
         allocate(overlap_pert(pos)%ann_orbs(1)) 
         overlap_pert(pos)%nannihilate = 1
         pops_pert(pos)%nannihilate = 1
         pops_pert(pos)%ann_orbs(1) = orbital
         overlap_pert(pos)%ann_orbs(1) = orbital
         call init_perturbation_annihilation(overlap_pert(pos))
         call init_perturbation_annihilation(pops_pert(pos))
      end if
    end subroutine setup_single_perturbation

!------------------------------------------------------------------------------------------!

    subroutine dealloc_real_time_memory
      use replica_data, only: clean_iter_data
      implicit none
      
      integer :: ierr
      
      if(numSnapshotOrbs>0) deallocate(snapShotOrbs,stat=ierr)
      if(allocated(numCycShiftExcess)) deallocate(numCycShiftExcess, stat=ierr)
      deallocate(DiagVec,stat=ierr)
      call clean_iter_data(second_spawn_iter_data)
      deallocate(shift_damping, stat=ierr)
      deallocate(temp_freeslot, stat=ierr)
      deallocate(current_overlap,stat=ierr)
      deallocate(gs_energy,stat=ierr)
      deallocate(dyn_norm_red,stat=ierr)
      deallocate(dyn_norm_psi,stat=ierr)
      deallocate(pert_norm,stat=ierr)
      deallocate(gf_overlap,stat=ierr)
      deallocate(TotPartsPeak,stat=ierr)
      deallocate(temp_det_list,stat=ierr)
      call clean_overlap_states()
      call clear_pops_pert(pops_pert)
      call clear_pops_pert(overlap_pert)
      
      if(allocated(tauCache)) deallocate(tauCache)

    end subroutine dealloc_real_time_memory

!------------------------------------------------------------------------------------------!

    subroutine clear_pops_pert(perturbs)
      use FciMCData, only : perturbation
      implicit none
      
      type(perturbation), intent(inout), allocatable :: perturbs(:)
      integer :: j

      if(allocated(perturbs)) then
         do j = 1, gf_count
            if(allocated(perturbs(j)%crtn_orbs)) deallocate(perturbs(j)%crtn_orbs)
            if(allocated(perturbs(j)%ann_orbs)) deallocate(perturbs(j)%ann_orbs)
            deallocate(perturbs)
         end do
      endif
    end subroutine clear_pops_pert

!------------------------------------------------------------------------------------------!

    subroutine read_in_trajectory
      use CalcData, only: nmcyc
      use real_time_data, only: tauCache, alphaCache
      implicit none
      integer :: i, iunit
      logical :: checkTraj
      character(*), parameter :: this_routine = "read_in_trajectory"
      
      call get_unique_filename('tauContour',.false.,.false.,0,trajFile)
      iunit = get_free_unit()
      inquire(file=trajFile,exist=checkTraj)
      if(.not. checkTraj) call stop_all(this_routine,"No tauContour file detected.")
      

      ! We first need to read in the number of cycles
      open(iunit,file=trajFile,status='old',position='append')
      backspace(iunit)
      read(iunit,*) nmcyc
      close(iunit)

      ! Then, the cache for the values of alpha and tau is allocated
      allocate(tauCache(nmcyc))
      allocate(alphaCache(nmcyc))

      ! And then the values for alpha/tau
      open(iunit,file=trajFile,status='old')
      ! Here, we read in the timesteps and alphas that shall be used
      do i = 1,nmcyc
         read(iunit,*) tauCache(i), alphaCache(i)
      enddo
      close(iunit)
    end subroutine read_in_trajectory

!------------------------------------------------------------------------------------------!

    subroutine truncate_initial_state
      use semi_stoch_procs, only: check_determ_flag
      use hash, only: hash_table_lookup, remove_hash_table_entry
      use FciMCData, only: HashIndex
      use bit_reps, only: nullify_ilut, decode_bit_det
      implicit none
      integer :: i, nI(nel), DetHash, PartInd
      logical :: tSuccess
      
      do i = 1, TotWalkers
         if(.not. check_determ_flag(CurrentDets(:,i))) then
            call nullify_ilut(CurrentDets(:,i))
            call decode_bit_det(nI,CurrentDets(:,i))
            call hash_table_lookup(nI,CurrentDets(:,i),nifdbo,HashIndex,CurrentDets,&
                 PartInd,DetHash,tSuccess)
            if(tSuccess) call remove_hash_table_entry(HashIndex,nI,PartInd)
         endif
      enddo
    end subroutine truncate_initial_state

!------------------------------------------------------------------------------------------!

    subroutine truncate_overlap_states
      use semi_stoch_procs, only: check_determ_flag
      use FciMCData, only: HashIndex
      use bit_reps, only: nullify_ilut
      use hash, only: hash_table_lookup
      implicit none
      integer :: i, iGf, nI(nel), DetHash, PartInd
      logical :: tSuccess

      do iGf = 1, gf_count
         ! For each gf, truncate the corresponding overlap state
         do i = 1, overlap_states(iGf)%ndets
            call decode_bit_det(nI,overlap_states(iGf)%dets(:,i))
            call hash_table_lookup(nI,overlap_states(iGf)%dets(:,i),nifdbo,HashIndex,&
                 CurrentDets,PartInd,DetHash,tSuccess)
            if(tSuccess) then
               ! In principle, there are no non-core determinants left when this
               ! is called, but we do not want to depend on the order in which the
               ! truncation is carried out
               if(.not. check_determ_flag(CurrentDets(:,PartInd))) call nullify_ilut(&
                    overlap_states(iGf)%dets(:,i))
            else
               call nullify_ilut(overlap_states(iGf)%dets(:,i))
            endif
         enddo
      enddo
    end subroutine truncate_overlap_states

!------------------------------------------------------------------------------------------!

    subroutine set_deterministic_flag(i)
      use FciMCData, only: HashIndex, core_space
      use bit_rep_data, only: flag_deterministic
      use hash, only: hash_table_lookup
      use bit_reps, only: set_flag
      implicit none
      integer, intent(in) :: i
      integer ::  nI(nel), DetHash, PartInd
      logical :: tSuccess
      character(*), parameter :: this_routine = "set_deterministic_flag"

      ! For each core-state, we check if it is in the CurrentDets (which should
      ! always be the case in the initialization
      call decode_bit_det(nI,core_space(:,i))
      call hash_table_lookup(nI,core_space(:,i),nifdbo,HashIndex,CurrentDets,PartInd,&
           DetHash, tSuccess)
      if(tSuccess) then
         ! And then we set the deterministic flag
         call set_flag(CurrentDets(:,PartInd),flag_deterministic)
      else
         call stop_all(this_routine,"Deterministic state not present in CurrentDets")
      endif
    end subroutine set_deterministic_flag

!------------------------------------------------------------------------------------------!

    subroutine initialize_corespace_construction
      use hash, only: init_hash_table
      use FciMCData, only: nWalkerHashes, MaxWalkersPart
      use real_time_data, only: ssht, core_space_buf, csbuf_size
      implicit none
      integer :: ierr

      ! allocate the buffer for the corespace
      allocate(core_space_buf(0:niftot,MaxWalkersPart),stat = ierr)
      ! set the position to 0
      csbuf_size = 0
      ! setup the buffer hashtable
      allocate(ssht(nWalkerHashes),stat = ierr)
      call init_hash_table(ssht)

    end subroutine initialize_corespace_construction
    
end module real_time_init
