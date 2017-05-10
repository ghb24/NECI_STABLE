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
                              AllAnnihilated_1, AllNoAborted_1, AllNoRemoved_1, &
                              AllNoAddedInitiators_1, AllNoInitDets_1, AllNoNonInitDets_1, &
                              AllNoInitWalk_1, AllNoNonInitWalk_1, AllInitRemoved_1, &
                              AccRat_1, AllNoatDoubs_1, AllSumWalkersCyc_1, current_overlap, &
                              TotPartsStorage,  t_rotated_time, TotPartsPeak, asymptoticShift, &
                              tau_imag, tau_real, elapsedRealTime, elapsedImagTime, &
                              TotWalkers_orig, dyn_norm_psi, gs_energy, shift_damping, &
                              tStaticShift, MaxSpawnedDiag, tDynamicCoreSpace, overlap_states, &
                              overlap_real, overlap_imag, allGfs, tRealTimePopsfile, &
                              tLimitShift, nspawnMax, shiftLimit, numCycShiftExcess, &
                              TotPartsLastAlpha, alphaDamping, tDynamicDamping, iterInit, &
                              etaDamping, tStartVariation, rotThresh, stabilizerThresh, &
                              tInfInit,  popSnapshot, snapshotOrbs, phase_factors, tVerletSweep, &
			      numSnapshotOrbs, tLowerThreshold, t_kspace_operators, tVerletScheme
    use real_time_procs, only: create_perturbed_ground, setup_temp_det_list, &
                               calc_norm, clean_overlap_states
    use verlet_aux, only: backup_initial_state, setup_delta_psi
    use constants, only: dp, n_int, int64, lenof_sign, inum_runs
    use Parallel_neci
    use ParallelHelper, only: iProcIndex, root, MPIbarrier, nNodes, MPI_SUM
    use util_mod, only: get_unique_filename, get_free_unit
    use Logging, only: tIncrementPops
    use kp_fciqmc_data_mod, only: tMultiplePopStart, tScalePopulation, &
                                  tOverlapPert, overlap_pert, scaling_factor
    use CalcData, only: tChangeProjEDet, tReadPops, tRestartHighPop, tFCIMC, &
                        tStartSinglePart, tau, nmcyc, iPopsFileNoRead, tWritePopsNorm, &
                        tWalkContGrow, diagSft, pops_norm, InitWalkers, MemoryFacSpawn, &
                        StepsSft, tSemiStochastic, tTruncInitiator, tAddToInitiator
    use FciMCData, only: tSearchTau, alloc_popsfile_dets, pops_pert, tPopsAlreadyRead, &
                         tSinglePartPhase, iter_data_fciqmc, iter, PreviousCycles, &
                         AllGrowRate, spawn_ht, pDoubles, pSingles, TotParts, &
                         MaxSpawned, tSearchTauOption, TotWalkers, SumWalkersCyc, &
                         CurrentDets, popsfile_dets, MaxWalkersPart, WalkVecDets
    use SystemData, only: lms, G1, nBasisMax, tHub, nel, tComplexWalkers_RealInts, nBasis, tReal
    use SymExcitDataMod, only: kTotal
    use sym_mod, only: MomPbcSym
    use perturbations, only: init_perturbation_annihilation, &
                             init_perturbation_creation
    use fcimc_initialisation, only: SetupParameters, InitFCIMCCalcPar, &
                                    init_fcimc_fn_pointers 
    use kp_fciqmc_init, only: create_overlap_pert_vec
    use LoggingData, only: tZeroProjE, tFCIMCStats2
    use fcimc_output, only: write_fcimcstats2, WriteFciMCStatsHeader
    use replica_data, only: allocate_iter_data, set_initial_global_data
    use bit_rep_data, only: nifbcast, niftot
    use bit_reps, only: decode_bit_det

    implicit none

    logical :: tmemoryfacdiagset
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
        complex(dp), allocatable :: norm_buf(:)
        integer :: ierr,  j, run, indicator

        ! also need to create the perturbed ground state to calculate the 
        ! overlaps to |y(t)> 

        ! so have to call this routine before the InitFCIMCCalcPar, where the 
        ! time evolved y(t) will be stored in the CurrentDets array

        call create_perturbed_ground()

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

        if(tRealTimePopsfile) call readTimeEvolvedState()

        ! allocate the according quantities! 
        ! n_time_steps have to be set here!
        print *, " Allocating greensfunction and wavefunction norm arrays!"
        ! allocate an additional slot for initial values
        normsize = inum_runs**2
        if(numSnapshotOrbs>0) then 
           allocate(popSnapshot(numSnapshotOrbs),stat=ierr)
           popSnapshot = 0
        endif
        allocate(overlap_real(gf_count),overlap_imag(gf_count))
        allocate(gf_overlap(normsize,gf_count), stat = ierr)
        allocate(pert_norm(normsize,gf_count),stat = ierr)
        allocate(dyn_norm_psi(normsize),stat = ierr)
        dyn_norm_psi = 1.0_dp
        allocate(dyn_norm_red(normsize,gf_count), stat = ierr)
        allocate(gs_energy(inum_runs),stat = ierr)
        allocate(current_overlap(normsize,gf_count),stat = ierr)
        allocate(temp_freeslot(MaxWalkersPart),stat = ierr)
        allocate(TotPartsPeak(inum_runs),stat = ierr)
        if(.not. allocated(shift_damping)) then 
           allocate(shift_damping(inum_runs), stat = ierr)
           shift_damping = 0.0_dp
        endif
        if(tLimitShift) then 
           allocate(numCycShiftExcess(inum_runs), stat = ierr)
           numCycShiftExcess = 0
        endif
        ! allocate spawn buffer for verlet scheme
        if(tVerletScheme) allocate(spawnBuf(0:niftot,1:maxSpawned))

        gf_overlap = 0.0_dp
        TotPartsPeak = 0.0_dp
        gs_energy = benchmarkEnergy
        
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
           dyn_norm_red(:,j) = pert_norm(:,j)
        enddo

        deallocate(norm_buf)
        
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
        benchmarkEnergy = 0.0_dp
        numCycShiftExcess = 0

        tVerletSweep = .false.
        if(tVerletScheme) then 
           call setup_delta_psi()
           call backup_initial_state()
        endif

        call rotate_time()           

    end subroutine setup_real_time_fciqmc

    subroutine rotate_time()
      ! to avoid code multiplication
      if(t_rotated_time) then
         tau_imag = - sin(real_time_info%time_angle)*tau
         tau_real = cos(real_time_info%time_angle)*tau
      else
         tau_imag = 0.0_dp
         tau_real = tau
      endif
    end subroutine rotate_time

    ! need a real-time calc read_input routine to seperate that as much 
    ! from the rest of the code as possible! 
    subroutine real_time_read_input()
        use input_neci
        implicit none
        logical :: eof
        character(100) :: w
        character(*), parameter :: this_routine = "real_time_read_input"
        integer, parameter :: lesser = -1, greater = 1
        integer :: i,j
	integer, allocatable :: buffer(:)

        ! set the flag that this is a real time calculation
        t_real_time_fciqmc = .true.

        ! and set default values for the real-time calculation
        call set_real_time_defaults()

        real_time: do
            call read_line(eof)
            if (eof) then
                ! end of file reached
                exit
            end if
            call readu(w)

            select case (w)
            ! have to enter all the different input options here

            case("VERLET")
               ! using a verlet algorithm instead of the second order runge-kutta
               tVerletScheme = .true.
               if(item < nitems) call readi(iterInit)

            case ("DAMPING")
                ! to reduce the explosive spread of walkers through the 
                ! Hilbert space a small imaginery energy can be introduced in
                ! the Schroedinger equation id/dt y(t) = (H-E0-ie)y(t)
                call readf(real_time_info%damping)

             case("ROTATE-TIME")
                ! If the time is to be rotated by some angle time_angle to increase 
                ! stability, this can be set here
                t_rotated_time = .true.
                call readf(real_time_info%time_angle)

            ! use nicks perturbation & kp-fciqmc stuff here as much as 
            ! possible too

            ! just compute the time-evolution of a singly excited state (with
            ! reference to the ground state. This gives the contribution of
            ! this state to the spectrum
             case("SINGLE")
            ! deprecated, replace by MULTI
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.
                allocate(pops_pert(1))
                pops_pert%nannihilate = 1
                pops_pert%ncreate = 1
                allocate(pops_pert(1)%crtn_orbs(1))
                allocate(pops_pert(1)%ann_orbs(1))
                call readi(pops_pert(1)%ann_orbs(1))
                call readi(pops_pert(1)%crtn_orbs(1))
                call init_perturbation_annihilation(pops_pert(1))
                call init_perturbation_creation(pops_pert(1))
                allocate(overlap_pert(1))
                overlap_pert%nannihilate = 1
                overlap_pert%ncreate = 1
                allocate(overlap_pert(1)%crtn_orbs(1))
                allocate(overlap_pert(1)%ann_orbs(1))
                overlap_pert(1)%crtn_orbs(1)=pops_pert(1)%crtn_orbs(1)
                overlap_pert(1)%ann_orbs(1)=pops_pert(1)%ann_orbs(1)
                call init_perturbation_annihilation(overlap_pert(1))
                call init_perturbation_creation(overlap_pert(1))

            ! TODO: MAJOR: Add feature for applying multiple annihilation operators
            ! AND: Add feature for applying momentum annihilators to states in the
            ! real-space basis

             case("KSPACE")
                ! Apply the perturbations in kspace. This only does something for real
                ! space hubbard, else it is the default
                t_kspace_operators = .true.

             ! Arbitrary perturbation on the initial state, always get the overlap
             ! with the initial state
             case("MULTI")
                alloc_popsfile_dets = .true.
                tWritePopsNorm = .true.
                allocate(pops_pert(1))
                allocate(overlap_pert(1))
                allocate(buffer(nel))
                j=0
                ! first, read all orbitals to which particles shall be added
                do
                   if(item==nitems) exit
                   call readi(i)
                   ! -1 is the terminator for creation and indicates that all following 
                   ! orbitals are to be annihilated
                   if(i==-1) exit
                   j=j+1
                   ! as the size of pops_pert is unknown, use a buffer
                   buffer(j)=i
                enddo
                ! allocate creation operators
                if(j>0) then
                   pops_pert%ncreate = j
                   allocate(pops_pert(1)%crtn_orbs(j))
                ! we take the overlap with the initial state, so overlap_pert == pops_pert
                   overlap_pert%ncreate = j
                   allocate(overlap_pert(1)%crtn_orbs(j))
                   do i=1,j
                      pops_pert(1)%crtn_orbs(i)=buffer(i)
                      overlap_pert(1)%crtn_orbs(i)=pops_pert(1)%crtn_orbs(i)
                   end do
                endif
                j=0
                ! now, read in all orbitals from which particles shall be removed
                do
                   if(item==nitems) exit
                   j=j+1
                   call readi(buffer(j))
                enddo
                ! again, allocate annihilation operators
                if(j>0) then
                   pops_pert%nannihilate = j
                   overlap_pert%nannihilate = j
                   allocate(pops_pert(1)%ann_orbs(j))
                   allocate(overlap_pert(1)%ann_orbs(j))
                   do i=1,j
                      pops_pert(1)%ann_orbs(i)=buffer(i)
                      overlap_pert(1)%ann_orbs(i)=pops_pert(1)%ann_orbs(i)
                   end do
                endif
                call init_perturbation_annihilation(pops_pert(1))
                call init_perturbation_creation(pops_pert(1))
                call init_perturbation_annihilation(overlap_pert(1))
                call init_perturbation_creation(overlap_pert(1))
                deallocate(buffer)
                
            ! the most important info is if it is the photoemmission(lesser GF)
            ! or photoabsorption (greater GF) and the orbital we want the 
            ! corresponding operator apply on 
            ! the type of GF considered also changes the sign of the FT exponent

            ! decision for now: input a specific GF matrix element and the type 
            ! of the greensfunction to be calculated(lesser,greater) eg:
            ! lesser i j : <y(0)| a^+_i a_j |y(0)> 
            case ("LESSER")
               alloc_popsfile_dets = .true.
                ! lesser GF -> photo emission: apply a annihilation operator
                tOverlapPert = .true.
                tWritePopsNorm = .true.
                ! i probably also can use the overlap-perturbed routines 
                ! from nick
                ! but since applying <y(0)|a^+_i for all i is way cheaper 
                ! and should be done for all possible and allowed i. 
                ! and creating all those vectors should be done in the init
                ! step and stored, and then just calc. the overlap each time 
                ! step

                ! store the information of the type of greensfunction 
                gf_type = lesser
               

                ! probably have to loop over spin-orbitals dont i? yes!

                ! if no specific orbital is specified-> loop over all j! 
                ! but only do that later: input is a SPINORBITAL!
                if(item < nitems) then
                   allocate(pops_pert(1))
                   pops_pert%nannihilate = 1
                   allocate(pops_pert(1)%ann_orbs(1))
                   call readi(pops_pert(1)%ann_orbs(1))
                   call init_perturbation_annihilation(pops_pert(1)) 
                else 
                   call stop_all(this_routine, "Invalid input for Green's function")  
                endif 
                if (nitems == 3) then
                   gf_count = 1
                   !allocate the perturbation object

                   ! and also the lefthand perturbation object for overlap
                   allocate(overlap_pert(1))
                   overlap_pert%nannihilate = 1
                   allocate(overlap_pert(1)%ann_orbs(1))

                   ! read left hand operator first
                   call readi(overlap_pert(1)%ann_orbs(1))
                   call init_perturbation_annihilation(overlap_pert(1))

                else
                   if(nitems == 2) then
                      allGfs = 1
                   else
                      call stop_all(this_routine, "Invalid input for Green's function")   
                   endif
                endif
             case ("GREATER")
                ! greater GF -> photo absorption: apply a creation operator
                alloc_popsfile_dets = .true.
                tOverlapPert = .true.
                tWritePopsNorm = .true.

                ! i probably also can use the overlap-perturbed routines 
                ! from nick
                ! but since applying <y(0)|a_i for all i is way cheaper 
                ! and should be done for all possible and allowed i. 
                ! and creating all those vectors should be done in the init
                ! step and stored, and then just calc. the overlap each time 
                ! step

                ! store type of greensfunction
                gf_type = greater

                ! if no specific orbital is specified-> loop over all j! 
                ! but only do that later
                if(item < nitems) then
                    allocate(pops_pert(1))                   
                    pops_pert%ncreate = 1
                    allocate(pops_pert(1)%crtn_orbs(1))
                    call readi(pops_pert(1)%crtn_orbs(1))
                    call init_perturbation_creation(pops_pert(1)) 
                 else
                    call stop_all(this_routine, "Invalid input for Green's function")   
                endif
                if (nitems == 3) then
                    ! allocate the perturbation object
                    allocate(overlap_pert(1))
                    overlap_pert%ncreate = 1
                    allocate(overlap_pert(1)%crtn_orbs(1))
                    call readi(overlap_pert(1)%crtn_orbs(1))
                    call init_perturbation_creation(overlap_pert(1))
                else
                   if(nitems == 2) then
                      allGfs = 2
                   else
                      call stop_all(this_routine, "Invalid input for Green's function")   
                endif
             endif

            case ("SCALE-POPULATION")
                tScalePopulation = .true.
		
	    case("LOWER-THRESHOLD")
		tLowerThreshold = .true.

             case ("FULLY-ROTATED")
                ! for testing purposes, it is useful to do pure imaginary
                ! time evolution with the rotated time algorithm -> this is
                ! enabled by this keyword 
                ! in addition, this disables the usage of input POPSFILEs for
                ! more efficient ground state search (the real-time POPSFILE 
                ! read-in settings are not useful for ground state search)
                tStartSinglePart = .true.
                t_rotated_time = .true.
                tWalkContGrow = .true.
                real_time_info%time_angle = 2*atan(1.0_dp)
                
             case("PRINT-POP")
                ! include the time-dependent population of targeted orbitals into
                ! the output. This requires them to be evaluated on the fly
                numSnapShotOrbs = 0
                allocate(buffer(nitems+1))
                do
                   if(item < nitems) then
                      numSnapShotOrbs = numSnapShotOrbs + 1
                      ! nBasis is not defined at this point, so we cannot check if
                      ! there are too many items given - no serious input will contain
                      ! more arguments than basis states anyway
                      call readi(buffer(numSnapShotOrbs))
                   else
                      exit
                   endif
                end do
                allocate(snapShotOrbs(numSnapShotOrbs))
                snapShotOrbs(1:numSnapShotOrbs) = buffer(1:numSnapShotOrbs)
                deallocate(buffer)

             case("NOSHIFT")
                ! disabling the shift gives higher precision results as no
                ! renormalization of the norm by a dynamic factor is made
                ! note that the walker number will grow exponentially in this
                ! scenario, however
                asymptoticShift = 0.0_dp
                ! might want to set DiagSft = 0.0_dp but there migth also be some cases
                ! in which this is unwanted
                tStaticShift = .true.

             case("START-HF")
                ! do not read in an initial state from a POPSFILE and apply a perturbation
                ! but start right away in the HF as the initial state does not matter in 
                ! principle for the spectrum
                tReadPops = .false.
               
             case("STABILIZE-WALKERS")
                ! enabling this activates the dynamic shift as soon as the walker number drops
                ! below 80% of the peak value
                tStabilizerShift = .true.
                if(item < nitems) then
                   call readf(asymptoticShift)
                   tStaticShift = .true.
                endif

             case("ENERGY-BENCHMARK")
                ! one can specify an energy which shall be added as a global shift
                ! to the hamiltonian. Useful for getting transition energies
                call readf(benchmarkEnergy)

             case("DYNAMIC-CORE")
                tDynamicCoreSpace = .true.
                ! if dynamic core is set, the core space for semistochastic treatment is 
                ! updated every few hundred iterations according to the currently most
                ! occupied determinants
                
            case ("COMPLEX-INTEGRALS")
                ! in the real-time implementation, since we need the complex 
                ! functionality anyway, we have to additionally tell the 
                ! program, that the FCIDUMP input is complex
                ! the default is that they are real! 
                t_complex_ints = .true.
                
             case("NSPAWNMAX")
                ! specify a maximum number of spawn attempts per determinant in
                ! regulation mode (i.e. for large number of spawns)
                call readi(nspawnMax)
                
             case("RT-POPS")
                ! in addition to the 'normal' popsfile, a second one is supplied
                ! containing a time evolved state
                tRealTimePopsfile = .true.

             case("OVERPOPULATE")
                ! enabling sets the options for time-dependent shift and rotation 
                ! such that a positive shift will occur with a stable walker number
                t_rotated_time = .true.
                tDynamicAlpha = .true.
                ! this is done by pinning the shift to some positive value and
                ! then auto-adjusting the rotation
                tStaticShift = .true.
                ! here, rotation and shift variation have to start at the same point 
                ! (in principle, it is only required that the rotation does not start 
                ! before shift variation) to prevent the rotation from converging on
                ! its own, circumventing the overpopulation via positive shift
                tOverpopulate = .true.
                ! it is most efficient to turn on the shift after equilibration of the angle
                ! so this is done via the stabilize-walkers feature
                tStabilizerShift = .true.
                if(item < nitems) then
                   call readf(asymptoticShift)
                else
                   asymptoticShift = 2.0_dp
                endif
               
             case("DYNAMIC-ROTATION")
                ! this automatically adjusts the temporal rotation to find a minimal 
                ! alpha guaranteeing a fixed walker number
                tDynamicAlpha = .true.
                t_rotated_time = .true.
                if(item < nitems) call readf(alphaDamping)

             case("ROTATION-THRESHOLD")
                ! number of walkers at which the variation of rotation angle starts
                ! 0 by default
                call readi(rotThresh)

             case ("STEPSALPHA")
                ! length of the decay channel update cycle (in timesteps)
                ! i.e. angle of rotation and damping
                call readi(stepsAlpha)

             case("DYNAMIC-DAMPING")
                ! allow the damping to be time-dependent 
                ! optional: damping parameter for the adjustment of eta
                tDynamicDamping = .true.
                if(item < nitems) call readf(etaDamping)

             case("LIMIT-SHIFT")
                ! limits the shift to some maximum value. On short times, the threshold
                ! can be exceeded.
                tLimitShift = .true.
                ! optional argument: threshold value (absolute value!). Default is 3
                if(item < nitems) call readf(shiftLimit)

             case("INFINITE-INIT")
                ! use the initiator adaptiation without any inititators - works well
                ! in some real-time applications
                ! this is not equivalent to switching on initiators without the
                ! addtoinitiator keyword as infinite-init will also remove all
                ! existing inititators
                tInfInit = .true.
                tAddtoInitiator = .true.
                tTruncInitiator = .true.

            case ("ENDREALTIME")
                exit real_time

            case default
                ! copy from george: write possible options if a false input
                ! is given. todo

            end select 
        end do real_time

    end subroutine real_time_read_input

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

        ! for the start definetly not change tau
        tSearchTau = .true.

        ! also set readpops to get the <y(0)| reference from the "normal"
        ! neci init routines
        tReadPops = .true.

        ! startsinglepart does not work with popsfile and is not wanted too
        tStartSinglePart = .false.

        ! but to ensure that the shift does not vary anymore, since there is 
        ! no such concept as the varying shift in the real-time fciqmc
        ! exception: when using rotated times, the shift still has to be considered
        tWalkContGrow = .true.
        tSinglePartPhase = .true.

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

    end subroutine dealloc_real_time_memory

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
    
end module real_time_init
