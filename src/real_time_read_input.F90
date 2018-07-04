module real_time_read_input_module
  use real_time_data
  use real_time_init, only: set_real_time_defaults, benchmarkenergy
  use FciMCData, only: alloc_popsfile_dets, pops_pert
  use CalcData, only: tAddToInitiator, tTruncInitiator, tWalkContGrow, tStartSinglePart, &
       tWritePopsNorm, tReadPops, ss_space_in, tSemiStochastic
  use perturbations, only: init_perturbation_creation, init_perturbation_annihilation
  use kp_fciqmc_data_mod, only: tOverlapPert, overlap_pert, tScalePopulation
  use SystemData, only: nel, tComplexWalkers_RealInts
  use bit_rep_data, only: tUseFlags
  use constants

  contains

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
               if(stepsAlpha .eq. 1) write(6,*) "Warning: STEPSALPHA is 1. Ignoring VERLET keyword"

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
            case("PROJECT-INITIAL-STATE")
               ! If we specify this, we do not create the overlap state specifically,
               ! but copy the currentdets. 
               ! This can greatly help to overcome memory problems, as it basically
               ! halves the memory required for initial state preparation
               tNewOverlap = .false.

            ! just compute the time-evolution of a singly excited state (with
            ! reference to the ground state. This gives the contribution of
            ! this state to the spectrum
             case("SINGLE")
                alloc_popsfile_dets = .true.
            ! deprecated, replace by MULTI
                tWritePopsNorm = .true.
                ! Now, overlap state and initial state are the same
                tNewOverlap = .false.
                allocate(pops_pert(1))
                pops_pert%nannihilate = 1
                pops_pert%ncreate = 1
                allocate(pops_pert(1)%crtn_orbs(1))
                allocate(pops_pert(1)%ann_orbs(1))
                call readi(pops_pert(1)%ann_orbs(1))
                call readi(pops_pert(1)%crtn_orbs(1))
                call init_perturbation_annihilation(pops_pert(1))
                call init_perturbation_creation(pops_pert(1))


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
                   
                   ! If the created and annihilated orbital are the same, we 
                   ! do not need to explicitly construct the projection state,
                   ! this might save a lot of memory
                   if(pops_pert(1)%ann_orbs(1) .eq. overlap_pert(1)%ann_orbs(1)) &
                        tNewOverlap = .false.


                else
                   if(nitems == 2) then
                      allGfs = 1
                   else
                      call stop_all(this_routine, "Invalid input for Green's function")   
                   endif
                endif
                
             case ("GREATER")
                alloc_popsfile_dets = .true.
                ! greater GF -> photo absorption: apply a creation operator
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

                   ! If the created and annihilated orbital are the same, we 
                   ! do not need to explicitly construct the projection state,
                   ! this might save a lot of memory
                   if(pops_pert(1)%crtn_orbs(1) .eq. overlap_pert(1)%crtn_orbs(1)) &
                        tNewOverlap = .false.

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
               ! indicates that the given rotation-threshold is not an upper
               ! but in fact a lower threshold, so the variation is switched
               ! on as soon as the walker number drops below
               ! this is a particularly useless thing in most cases, but for
               ! proving some stuff, it saves the day
		tLowerThreshold = .true.

             case ("FULLY-ROTATED")
                ! for testing purposes, it is useful to do pure imaginary
                ! time evolution with the rotated time algorithm -> this is
                ! enabled by this keyword 
                ! in addition, this disables the usage of input POPSFILEs for
                ! more efficient ground state search (the real-time POPSFILE 
                ! read-in settings are not useful for ground state search)
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
                tStartSinglePart = .true.
               
             case("STABILIZE-WALKERS")
                ! enabling this activates the dynamic shift as soon as the walker number drops
                ! below 80% of the peak value
                tStabilizerShift = .true.
                if(item < nitems) then
                   call readf(asymptoticShift)
                   tStaticShift = .true.
                endif

             case("UNCONSTRAINED-SHIFT")
                ! use an unconstrained shift mode that also allows
                ! negative shifts
                tOnlyPositiveShift = .false.
                write(iout,*) &
                     "WARNING: Using an unconstrained shift can lead to instabilities"

             case("HF-OVERLAP")
                ! take the overlap not with the initial state but with the perturbed
                ! reference
                tHFOverlap = .true.

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
                
             case("COMPLEXWALKERS-COMPLEXINTS")
                ! if we really use complex integrals, we have to tell as the
                ! default is using real integrals with complex walkers
                tComplexWalkers_RealInts = .false.
                
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
                if(stepsAlpha .eq. 1 .and. tVerletScheme) write(6,*) &
                     "Warning: STEPSALPHA is 1. Ignoring VERLET keyword"
                

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

             case("LOG-TRAJECTORY")
                ! This prints out the complex time trajectory in the form of alpha(iter)
                ! and tau(iter)
                tLogTrajectory = .true.
                
             case("GENERATE-CORESPACE")
                ! Now, we write out the most important determinants along the contour
                ! Also, the contour is logged
                tGenerateCoreSpace = .true.
                tLogTrajectory = .true.
                ! optionally, we can supply the number of states to log
                ss_space_in%tpops = .true.
                if(item < nitems) then 
                   call readi(ss_space_in%npops)
                else
                   ss_space_in%npops = 1000
                endif
                if(tSemiStochastic) call stop_all(this_routine, &
                     "GENERATE-CORESPACE NOT AVAILABLE IN SEMI-STOCHASTIC MODE")

             case("CORESPACE-THRESHOLD")
                ! Set the threshold from which on a determinant is in the corespace
                CALL readf(wn_threshold)
                
             case("CORESPACE-LOG-INTERVAL")
                ! Set the number of iterations after which we get the new candidates for the 
                ! corespace
                call readi(corespace_log_interval)
               
             case("READ-TRAJECTORY")
                ! This reads in a trajectory and performs the time-evolution along
                ! it
                tReadTrajectory = .true.
                
             case("LIVE-TRAJECTORY")
                ! Now we re-read the trajectory during runtime, this can be used to 
                ! use a trajectory that is currently being determined
                tReadTrajectory = .true.
                tLiveTrajectory = .true.

             case("CORESPACE-OVERLAP")
                ! Get the Green's function for the corespace only. This performs the 
                ! time-evolution only in the semistochastic space.
                tGZero = .true.
                ! If the corespace-greensfunction is to be obtained, semi-stochastic 
                ! has to be turned on
                if(.not. tSemiStochastic) call stop_all(this_routine, &
                     "CORESPACE-OVERLAP ONLY AVAILABLE IN SEMI-STOCHASTIC MODE")

            case ("ENDREALTIME")
                exit real_time

            case default
               call report("Keyword "//trim(w)//" not recognized in REALTIME block",.true.)

            end select 
        end do real_time

    end subroutine real_time_read_input

  end module real_time_read_input_module
