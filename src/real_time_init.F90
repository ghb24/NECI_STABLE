#include "macros.h"

! this is the initialization module for the real-time FCIQMC implementation

module real_time_init

    use real_time_data, only: t_real_time_fciqmc, gf_type
    use real_time_procs, only: create_perturbed_ground

    use constants, only: dp, n_int, int64, lenof_sign, inum_runs
    use Parallel_neci, only: nProcessors
    use ParallelHelper, only: iProcIndex, root
    use util_mod, only: get_unique_filename
    use Logging, only: tIncrementPops
    use CalcData, only: iPopsFileNoRead, tWritePopsNorm
    use kp_fciqmc_data_mod, only: scaling_factor, tMultiplePopStart, tScalePopulation, &
                                  tOverlapPert, overlap_pert
    use CalcData, only: tChangeProjEDet, tReadPops, tRestartHighPop, tFCIMC, &
                        tStartSinglePart
    use FciMCData, only: tSearchTau, alloc_popsfile_dets, pops_pert, tPopsAlreadyRead
    use SystemData, only: nBasis
    use perturbations, only: init_perturbation_annihilation, &
                             init_perturbation_creation
    use fcimc_initialisation, only: SetupParameters, InitFCIMCCalcPar, &
                                    init_fcimc_fn_pointers 

    use kp_fciqmc_init, only: create_overlap_pert_vec

    implicit none

contains

    subroutine init_real_time_calc_single()
        ! this routine takes care of the correct setup of the real-time 
        ! calculation. like reading the popsfiles and preparing the start 
        ! of the calculation and setting certain global variables
        implicit none
        character(*), parameter :: this_routine = "init_real_time_calc_single"

        print *, " Entering real-time FCIQMC initialisation "
        ! think about what variables have to be set for a succesful calc.

        ! do a last test if all he input is compatible and correctly set
        ! and abort calculation if incompatible input is detected
        call check_input_real_time()

        ! also call the "normal" NECI setup routines to allow calculation
        call SetupParameters()


        ! have to think about the the order about the above setup routines! 
        ! within this Init a readpops is called.. and 
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

    end subroutine init_real_time_calc_single

    subroutine setup_real_time_fciqmc()
        ! this is the last setup routine, which depending on compilation,
        ! number of copies etc. sets up the final needed quantities to run 
        ! a simulation
        character(*), parameter :: this_routine = "setup_real_time_fciqmc"

        ! also need to create the perturbed ground state to calculate the 
        ! overlaps to |y(t)> 

        ! so have to call this routine before the InitFCIMCCalcPar, where the 
        ! time evolved y(t) will be stored in the CurrentDets array
        call create_perturbed_ground()


    end subroutine setup_real_time_fciqmc

    ! need a real-time calc read_input routine to seperate that as much 
    ! from the rest of the code as possible! 
    subroutine real_time_read_input()
        use input_neci
        implicit none
        logical :: eof
        character(100) :: w
        character(*), parameter :: this_routine = "real_time_read_input"

        integer :: i

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
                ! lesser GF -> photo emission: apply a annihilation operator
                alloc_popsfile_dets = .true.
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
                gf_type = -1

                ! probably have to loop over spin-orbitals dont i? yes!

                ! if no specific orbital is specified-> loop over all j! 
                ! but only do that later: input is a SPINORBITAL!
                if (item < nitems) then
                    ! allocate the perturbation object
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
                end if

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
                gf_type = 1

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
    subroutine set_real_time_defaults()
        implicit none

        ! todo: figure out quantities

        ! for the start definetly not change tau
        tSearchTau = .false.

        ! also set readpops to get the <y(0)| reference from the "normal"
        ! neci init routines
        tReadPops = .true.

        ! startsinglepart does not work with popsfile and is not wanted too
        tStartSinglePart = .false.

        ! probably not change reference anymore.. but check
        tChangeProjEDet = .false.

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

    end subroutine set_real_time_defaults

    subroutine check_input_real_time()
        ! routine to ensure all calculation parameter are set up correctly
        ! and abort otherwise 
        character(*), parameter :: this_routine = "check_input_real_time"


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
                   PopBlockingIter, Popinum_runs, PopRandomHash(1024), &
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



end module real_time_init
