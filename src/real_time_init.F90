#include "macros.h"

! this is the initialization module for the real-time FCIQMC implementation

module real_time_init

    use real_time_data, only: t_real_time_fciqmc, gf_type, real_time_info, &
                              t_complex_ints, gf_overlap, wf_norm, temp_det_list, &
                              temp_det_pointer, temp_det_hash, temp_freeslot, &
                              pert_norm, second_spawn_iter_data, DiagParts, &
                              DiagParts2, DiagVec, DiagVec2, valid_diag_spawn_list, &
                              NoatHF_1, Annihilated_1, Acceptances_1, NoBorn_1, &
                              SpawnFromSing_1, NoDied_1, NoAborted_1, NoRemoved_1, &
                              NoAddedInitiators_1, NoInitDets_1, NoNonInitDets_1, &
                              NoInitWalk_1, NoNonInitWalk_1, InitRemoved_1, &
                              AllNoatHF_1, AllNoatHF_1, AllGrowRate_1, AllGrowRateAbort_1, &
                              AllNoBorn_1, AllSpawnFromSing_1, AllNoDied_1, &
                              AllAnnihilated_1, AllNoAborted_1, AllNoRemoved_1, &
                              AllNoAddedInitiators_1, AllNoInitDets_1, AllNoNonInitDets_1, &
                              AllNoInitWalk_1, AllNoNonInitWalk_1, AllInitRemoved_1, &
                              AccRat_1, AllNoatDoubs_1, AllSumWalkersCyc_1, current_overlap, &
                              TotPartsStorage, TotWalkers_pert, t_rotated_time, time_angle, &
                              tau_imag, tau_real, elapsedRealTime, elapsedImagTime
    use real_time_procs, only: create_perturbed_ground, setup_temp_det_list, &
                               calc_norm
    use constants, only: dp, n_int, int64, lenof_sign, inum_runs
    use Parallel_neci
    use ParallelHelper, only: iProcIndex, root, MPIbarrier, nNodes, MPI_SUM
    use util_mod, only: get_unique_filename
    use Logging, only: tIncrementPops
    use kp_fciqmc_data_mod, only: scaling_factor, tMultiplePopStart, tScalePopulation, &
                                  tOverlapPert, overlap_pert, perturbed_ground
    use CalcData, only: tChangeProjEDet, tReadPops, tRestartHighPop, tFCIMC, &
                        tStartSinglePart, tau, nmcyc, iPopsFileNoRead, tWritePopsNorm, &
                        tWalkContGrow, diagSft, pops_norm
    use FciMCData, only: tSearchTau, alloc_popsfile_dets, pops_pert, tPopsAlreadyRead, &
                         tSinglePartPhase, iter_data_fciqmc, iter, PreviousCycles, &
                         AllGrowRate, spawn_ht, pDoubles, pSingles, TotParts, &
                         MaxSpawned, InitialSpawnedSlots, tSearchTauOption, TotWalkers, &
                         CurrentDets, popsfile_dets, MaxWalkersPart
    use SystemData, only: nBasis, lms, G1, nBasisMax, tHub, nel
    use SymExcitDataMod, only: kTotal
    use sym_mod, only: MomPbcSym
    use perturbations, only: init_perturbation_annihilation, &
                             init_perturbation_creation
    use fcimc_initialisation, only: SetupParameters, InitFCIMCCalcPar, &
                                    init_fcimc_fn_pointers 
    use kp_fciqmc_init, only: create_overlap_pert_vec
    use LoggingData, only: tZeroProjE, tFCIMCStats2
    use fcimc_output, only: write_fcimcstats2, WriteFciMCStatsHeader
    use replica_data, only: allocate_iter_data
    use bit_rep_data, only: nifbcast
    use bit_reps, only: decode_bit_det

    implicit none

contains

    subroutine init_real_time_calc_single()
      use real_time_procs, only: get_tot_parts
        ! this routine takes care of the correct setup of the real-time 
        ! calculation. like reading the popsfiles and preparing the start 
        ! of the calculation and setting certain global variables
        implicit none
        character(*), parameter :: this_routine = "init_real_time_calc_single"

        integer :: ierr
        integer :: nI(nel), kt(3)

        print *, " Entering real-time FCIQMC initialisation "

        ! think about what variables have to be set for a succesful calc.

        ! do a last test if all he input is compatible and correctly set
        ! and abort calculation if incompatible input is detected
        call check_input_real_time()

        ! also call the "normal" NECI setup routines to allow calculation
        call SetupParameters()

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
        if(tHub) then
           if(allocated(pops_pert)) then
              if(pops_pert(1)%nannihilate == 1) kTotal = kTotal &
                   - G1(pops_pert(1)%ann_orbs(1))%k
              if(pops_pert(1)%ncreate == 1) kTotal = kTotal &
                   + G1(pops_pert(1)%crtn_orbs(1))%k
              call MomPbcSym(kTotal,nBasisMax)
           endif
        endif
        print *, "New total momentum", kTotal
        call MPIBarrier(ierr)

    end subroutine init_real_time_calc_single

    subroutine setup_real_time_fciqmc()
        ! this is the last setup routine, which depending on compilation,
        ! number of copies etc. sets up the final needed quantities to run 
        ! a simulation
        character(*), parameter :: this_routine = "setup_real_time_fciqmc"
        real(dp) :: norm_buf(inum_runs), tzero_norm(inum_runs)
        integer :: ierr, run

        ! also need to create the perturbed ground state to calculate the 
        ! overlaps to |y(t)> 

        ! so have to call this routine before the InitFCIMCCalcPar, where the 
        ! time evolved y(t) will be stored in the CurrentDets array
        call create_perturbed_ground()

        ! change the flags dependent on the real-time input 
        if (real_time_info%t_equidistant_time) then
            print *, " The equdistant time step option is set for the real-time calculation"
            tSearchTau = .false.
            ! also have to turn off:
            tSearchTauOption = .false.
            if (real_time_info%time_step > 0.0_dp) then
                print *, " A specific time-step is chosen by input!"
                tau = real_time_info%time_step
            else
                print *, " No specific time-step is chosen by input! Use tau from Popsfile!"
                !tSearchTau = .false.
                real_time_info%time_step = tau
            end if

            print *, " time-step: ", real_time_info%time_step
        end if

        ! initialize the storage containers for the calculated overlaps
        ! depending on the input 
        if (real_time_info%n_time_steps < 0) then
            ! if no specific number if iterations is inputted
            print *, " No specific number of time-steps is chosen by input!"

            if (real_time_info%t_equidistant_time) then
                ! if the max-time is set in input calculate from that the number 
                ! of elements 
                if (real_time_info%max_time > 0.0_dp) then
                    print *, " A maximum time is chosen by input: ", real_time_info%max_time
                    real_time_info%n_time_steps = int(real_time_info%max_time / tau )
                    print *, " this leads to n_time_steps of: ", real_time_info%n_time_steps
                else
                    ! otherwise calculate it from the number of iterations
                    print *, " No maximum time is chosen by input! "
                    real_time_info%n_time_steps = nmcyc
                    print *, " Use maximum cycle number as n_time_steps: ", nmcyc
                    real_time_info%max_time = real(nmcyc,dp) * tau 
                    print *, " This leads to maximum time: ", real_time_info%max_time
                end if
            else
                ! if the time-step is not fixed, we can only limit the calculation
                ! by the maximum number of cycles. 
                ! we could essentially also define a maximum time, but we cannot 
                ! initialize the gf overlap and wf norm variables depending on it
                ! so we have to ensure that a nmcyc is defined if the equidistant 
                ! time flag is set to false! 
                print *, " Use maximum cycle number as n_time_steps: ", nmcyc
                real_time_info%n_time_steps = nmcyc
                if (real_time_info%max_time > 0.0_dp) then
                    print *, " Maximum time is chosen by input: ", real_time_info%max_time
                    print *, " Stop simulation at whatever is reached first: n_time_steps or max_time!"
                else
                    print *, " No maximum time is chosen by input! Stop simulation after n_time_steps! "
                end if
            end if
        else
            print *, " A specific number of time-steps is chosen by input! "
            print *, " n_time_steps: ", real_time_info%n_time_steps
            if (real_time_info%t_equidistant_time) then
                ! check if input is consistent.. although already do that 
                ! in check_input routine! 
                ! no do it here, since in check_input the popsfile tau is not 
                ! yet known! 
                if (real_time_info%max_time > 0.0_dp) then
                    if (real(real_time_info%n_time_steps,dp)*real_time_info%time_step &
                        /= real_time_info%max_time) then
                        print *, " number of timesteps and time step not congruent with max_time!"
                        real_time_info%max_time = real(real_time_info%n_time_steps,dp) &
                            * real_time_info%time_step
                        print *, " setting max_time to: ", real_time_info%max_time
                    end if
                else 
                    ! set max_time! 
                    print *, " no max_time inputted! set it to be congruent with time-steps"
                    real_time_info%max_time = real(real_time_info%n_time_steps,dp) &
                        * real_time_info%time_step
                    print *, " setting max_time to: ", real_time_info%max_time
                end if
            else
                print *, " non-equidistant time-step used!" 
                if (real_time_info%max_time > 0.0_dp) then
                    print *, " Stop simulation at whatever is reached first: n_time_steps or max_time!"
                    print *, " max_time: ", real_time_info%max_time
                else
                    print *, " stop simulation after n_time_steps is reached!"
                end if
            end if
        end if

        ! allocate the according quantities! 
        ! n_time_steps have to be set here!
        print *, " Allocating greensfunction and wavefunction norm arrays!"
        ! allocate an additional slot for initial values
        allocate(gf_overlap(lenof_sign,0:(real_time_info%n_time_steps+1)), stat = ierr)
        allocate(wf_norm(inum_runs,0:(real_time_info%n_time_steps+1)), stat = ierr)
        allocate(pert_norm(inum_runs),stat = ierr)
        allocate(current_overlap(lenof_sign),stat=ierr)
        allocate(temp_freeslot(MaxWalkersPart),stat=ierr)

        gf_overlap = 0.0_dp

        ! to avoid dividing by 0 if not all entries get filled
        wf_norm = 1.0_dp

        pert_norm = 1.0_dp
        ! calc. the norm of this perturbed ground-state
        norm_buf = calc_norm(perturbed_ground,int(TotWalkers_pert))
        ! the norm (squared) can be obtained by reduction over all processes
        call MPIReduce(norm_buf,MPI_SUM,pert_norm)
        
        ! and the same thing again for the initial state
        norm_buf = calc_norm(CurrentDets,TotWalkers)
        call MPIReduce(norm_buf,MPI_SUM,tzero_norm)

        ! set the fist value here for now
        do run = 1,inum_runs
           wf_norm(run,0) = sqrt(tzero_norm(run) * pert_norm(run))
        end do

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

        ! also intitialize the 2nd spawning array to deal with the 
        ! diagonal death step in the 2nd rt-fciqmc loop
        ! have to make this list as large as the original walker list? 
        ! hm, probably not, as due to the tau multiplicative factor 
        ! for not hihgly populated dets also a succesful death step is not 
        ! so likely.. make it as big as MaxSpawned for now, maybe increase 
        ! later if this turns out to be to small..
        allocate(DiagVec(0:nifbcast, MaxSpawned), DiagVec2(0:nifbcast, MaxSpawned),&
            stat = ierr)

        DiagVec = 0
        DiagVec2 = 0

        DiagParts => DiagVec
        DiagParts2 => DiagVec2

        ! also need to setup this valid spawn list and crap..
        allocate(valid_diag_spawn_list(0:nNodes-1), stat = ierr)

        ! and the initial_spawn_slots equivalent
        ! although i think i can reuse the initialSpawnedSlots..
!         allocate(initial_diag_spawn_list(0:nNodes-1), stat = ierr) 

        valid_diag_spawn_list(:) = InitialSpawnedSlots(:)

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
        elapsedRealTime = 0.0_dp
        elapsedImagTime = 0.0_dp

        call rotate_time()

    end subroutine setup_real_time_fciqmc

    subroutine rotate_time()
      ! to avoid code multiplication
      if(t_rotated_time) then
         tau_imag = - sin(time_angle)*tau
         tau_real = cos(time_angle)*tau
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

        integer :: i
        integer, parameter :: lesser = -1, greater = 1

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

            case ("EQUIDISTANT-TIME")
                ! use equidistant time: do this for the beginning 
                ! thats probably helpful is i calculate the GFs for 
                ! different orbitals, to calculate the spectrum from all 
                ! different contributions
                ! set tau-search to false later on dependent on that! 
                real_time_info%t_equidistant_time = .true.

                ! its possible to input the wanted time-step here too
                if (item < nitems) then
                    call readf(real_time_info%time_step)
                end if

            case ("MAX-TIME")
                ! input the targeted end time 
                ! a specified   
                call readf(real_time_info%max_time)

            case ("NUM-TIMESTEPS")
                ! also able to directly input the number of calculated 
                ! gf time-steps
                ! this option is compatible with and without tau-search
                ! but not directly with the max-time option
                call readi(real_time_info%n_time_steps)

            case ("DAMPING")
                ! to reduce the explosive spread of walkers through the 
                ! Hilbert space a small imaginery energy can be introduced in
                ! the Schroedinger equation id/dt y(t) = (H-E0-ie)y(t)
                call readf(real_time_info%damping)

             case("ROTATE-TIME")
                ! If the time is to be rotated by some angle time_angle to increase 
                ! stability, this can be set here
                t_rotated_time = .true.
                tWalkContGrow = .false.
                call readf(time_angle)

            ! use nicks perturbation & kp-fciqmc stuff here as much as 
            ! possible too

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
                if (item < nitems) then
                    !allocate the perturbation object
                    allocate(pops_pert(1))
                    ! and also the lefthand perturbation object for overlap
                    allocate(overlap_pert(1))

                    pops_pert%nannihilate = 1
                    overlap_pert%nannihilate = 1

                    allocate(pops_pert(1)%ann_orbs(1))
                    allocate(overlap_pert(1)%ann_orbs(1))

                    ! read left hand operator first
                    call readi(overlap_pert(1)%ann_orbs(1))
                    call readi(pops_pert(1)%ann_orbs(1))

                    call init_perturbation_annihilation(overlap_pert(1))
                    call init_perturbation_annihilation(pops_pert(1)) 

                else
                    ! otherwise the whole possible orbitals are to be applied

                    print *, "no specific annihilation orbitals specified! loop over all!"
                    call stop_all(this_routine, "not yet implemented!")

                    allocate(pops_pert(nBasis))

                    do i = 1, nBasis
                        pops_pert(i)%nannihilate = 1
                        allocate(pops_pert(i)%ann_orbs(1))

                        pops_pert(i)%ann_orbs(1) = i 

                        call init_perturbation_annihilation(pops_pert(i))

                    end do
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
                if (item < nitems) then
                    ! allocate the perturbation object
                    allocate(pops_pert(1))
                    allocate(overlap_pert(1))

                    pops_pert%ncreate = 1
                    overlap_pert%ncreate = 1

                    allocate(pops_pert(1)%crtn_orbs(1))
                    allocate(overlap_pert(1)%crtn_orbs(1))

                    call readi(overlap_pert(1)%crtn_orbs(1))
                    call readi(pops_pert(1)%crtn_orbs(1))


                    call init_perturbation_creation(overlap_pert(1))
                    call init_perturbation_creation(pops_pert(1)) 

                else
                    ! otherwise the whole possible orbitals are to be applied
                    ! input is a SPINORBITAL!
                    allocate(pops_pert(nBasis))

                    print *, "no specific creation orbital specified! loop over all!"
                    call stop_all(this_routine, "not yet implented!")

                    do i = 1, nBasis
                        pops_pert(i)%ncreate = 1
                        allocate(pops_pert(i)%crtn_orbs(1))

                        pops_pert(i)%crtn_orbs(1) = i 

                        call init_perturbation_creation(pops_pert(i))

                    end do
                end if


            case ("SCALE-POPULATION")
                tScalePopulation = .true.

             case ("FULLY-ROTATED")
                ! for testing purposes, it is useful to do pure imaginary
                ! time evolution with the rotated time algorithm -> this is
                ! enabled by this keyword 
                ! in addition, this disables the usage of input POPSFILEs for
                ! more efficient ground state search (the real-time POPSFILE 
                ! read-in settings are not useful for ground state search)
                tReadPops = .false.
                tStartSinglePart = .true.
                t_rotated_time = .true.
                tWalkContGrow = .false.
                time_angle = 4*atan(1.0_dp)/2.0_dp
                
            case ("COMPLEX-INTEGRALS")
                ! in the real-time implementation, since we need the complex 
                ! functionality anyway, we have to additionally tell the 
                ! program, that the FCIDUMP input is complex
                ! the default is that they are real! 
                t_complex_ints = .true.

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

        ! todo: figure out quantities

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
        ! tSinglePartPhase is not yet allocated during readinput!
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

        ! from my way of outputting popsfiles i always do it in popsfile.n
        ! format -> so i probably have to set tIncrementPops and the count
        tIncrementPops = .true.
        iPopsFileNoRead = -1

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

    end subroutine set_real_time_defaults

    subroutine check_input_real_time()
        ! routine to ensure all calculation parameter are set up correctly
        ! and abort otherwise 
        character(*), parameter :: this_routine = "check_input_real_time"

        ! if non-equidistant time-steps are used we have to either have 
        ! nmcyc or the n_time_steps quantity set to be able to allocate 
        ! the gf overlap and wf norm quantities
        if (.not. real_time_info%t_equidistant_time) then
            ASSERT(nmcyc > 0 .or. real_time_info%n_time_steps > 0)
        end if



    end subroutine check_input_real_time

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

    subroutine dealloc_real_time_memory
      use replica_data, only: clean_iter_data
      implicit none
      
      integer :: ierr
      
      deallocate(valid_diag_spawn_list,stat=ierr)
      deallocate(DiagVec,stat=ierr)
      call clean_iter_data(second_spawn_iter_data)
      deallocate(current_overlap,stat=ierr)
      deallocate(pert_norm,stat=ierr)
      deallocate(wf_norm,stat=ierr)
      deallocate(gf_overlap,stat=ierr)
      

    end subroutine dealloc_real_time_memory

end module real_time_init
