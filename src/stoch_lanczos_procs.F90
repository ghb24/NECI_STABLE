#include "macros.h"

module stoch_lanczos_procs

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use bit_rep_data, only: NIfTot
    use bit_reps, only: decode_bit_det
    use CalcData, only: tTruncInitiator, tStartSinglePart, InitialPart, InitWalkers
    use constants
    use dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData, only: ilutHF, HFDet, CurrentDets, SpawnedParts, SpawnedParts2, TotWalkers
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, HashIndex, nWalkerHashes
    use FciMCData, only: fcimc_iter_data, ll_node, MaxWalkersPart, tStartCoreGroundState
    use FciMCData, only: tPopsAlreadyRead, tHashWalkerList
    use FciMCParMod, only: create_particle, InitFCIMC_HF, SetupParameters, InitFCIMCCalcPar
    use FciMCParMod, only: init_fcimc_fn_pointers, WriteFciMCStats, WriteFciMCStatsHeader
    use FciMCParMod, only: rezero_iter_stats_each_iter
    use hash, only: FindWalkerHash, reset_hash_table, fill_in_hash_table
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias
    use Parallel_neci, only: MPIBarrier
    use procedure_pointers
    use SystemData, only: nel

    implicit none

    type stoch_lanczos_data
        ! If true then a ground-state calculation is being performed.
        logical :: tGround
        ! If true then a finite-temperature calculation is being performed.
        logical :: tFiniteTemp
        ! The number of different initial walker configurations to start
        ! Lanczos calculations from.
        integer :: nconfigs
        ! The number of Lanczos calculations to perform for each initial walker
        ! configuration.
        integer :: nrepeats
        ! The number of different Lanczos vectors to sample (the number of
        ! vectors which form the Krylov subspace at the end of a calculation).
        integer :: nkrylov_vecs
        ! The number of iterations to perform *between each Lanczos vector being
        ! sampled*.
        integer :: niters
    end type

    type(stoch_lanczos_data) :: lanczos

    integer :: TotWalkers_Lanc
    integer(n_int), allocatable :: init_lanczos_config(:,:)

contains

    subroutine stoch_lanczos_read_inp()

        use input_neci

        logical :: eof
        character(len=100) :: w

        lanczos%nconfigs = 1
        lanczos%nrepeats = 1

        read_inp: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)
            case("END-LANCZOS")
                exit read_inp
            case("GROUND-STATE")
                lanczos%tGround = .true.
            case("FINITE-TEMPERATURE")
                lanczos%tFiniteTemp = .true.
            case("NUM-INIT-CONFIGS")
                call geti(lanczos%nconfigs)
            case("NUM-RUNS-PER-CONFIG")
                call geti(lanczos%nrepeats)
            case("NUM-LANCZOS-VECS")
                call geti(lanczos%nkrylov_vecs)
            case("NUM-ITERS-PER-VEC")
                call geti(lanczos%niters)
            case default
                call report("Keyword "//trim(w)//" not recognized in stoch-lanczos block", .true.)
            end select
        end do read_inp

    end subroutine stoch_lanczos_read_inp

    subroutine init_stoch_lanczos(lanczos)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer :: ierr
        character (len=*), parameter :: t_r = "init_stoch_lanczos"

        if (.not. tHashWalkerList) call stop_all('t_r','Stochastic Lanczos can only be run using &
            &the linscalefcimcalgo option (the linear scaling algorithm).')

        tPopsAlreadyRead = .false.
        call SetupParameters()
        call InitFCIMCCalcPar()
        call init_fcimc_fn_pointers() 

        call WriteFciMCStatsHeader()
        call WriteFCIMCStats()

        if (n_int == 4) call stop_all('t_r', 'Use of RealCoefficients does not work with 32 bit &
             &integers due to the use of the transfer operation from dp reals to 64 bit integers.')

        ! If performing a finite-temperature calculation with more than one run for each initial
        ! configuration, we store this walker configuration so that we can restart from it later.
        if (lanczos%tFiniteTemp .and. lanczos%nrepeats > 1) then
            allocate(init_lanczos_config(0:NIfTot, MaxWalkersPart), stat=ierr)
        end if

        call MPIBarrier(ierr)

    end subroutine init_stoch_lanczos

    subroutine init_stoch_lanczos_iter(iter_data, determ_index)

        use FciMCData, only: fcimc_excit_gen_store, FreeSlot, iStartFreeSlot, iEndFreeSlot

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(out) :: determ_index

        ! Reset positions to spawn into in the spawning array.
        ValidSpawnedList = InitialSpawnedSlots

        ! Reset the array which holds empty slots in CurrentDets.
        FreeSlot(1:iEndFreeSlot) = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0

        ! Index for counting deterministic states.
        determ_index = 1

        call rezero_iter_stats_each_iter(iter_data)

    end subroutine init_stoch_lanczos_iter

    subroutine create_initial_config(lanczos, irun)

        type(stoch_lanczos_data), intent(in) :: lanczos
        integer, intent(in) :: irun
        integer :: DetHash, i, nwalkers, nwalkers_target

        call reset_hash_table(HashIndex)

        if (lanczos%tGround) then
            ! If starting from the core ground state then walkers will be added to the list in the
            ! semi-stochastic routines. Otherwise, call InitFCIMC_HF, which just adds the desired
            ! initial number of walkers to the Hartree-Fock determinant.
            if (.not. tStartCoreGroundState) call InitFCIMC_HF()

        else if (lanczos%tFiniteTemp) then
            if (irun == 1) then
                ! Convert the initial number of walkers to an integer.
                if (tStartSinglePart) then
                    nwalkers_target = int(InitialPart)
                else
                    nwalkers_target = int(InitWalkers)
                end if
                ! Finally, call the routine to create the walker distribution.
                call generate_init_config_basic(nwalkers_target, nwalkers)
                TotWalkers = int(nwalkers, int64)
                TotWalkers_Lanc = TotWalkers
                ! If starting from this configuration more than once, store it.
                if (lanczos%nrepeats > 1) init_lanczos_config(:, 1:TotWalkers) = CurrentDets(:, 1:TotWalkers)
            else if (irun > 1) then
                TotWalkers = TotWalkers_Lanc
                CurrentDets(:, 1:TotWalkers) = init_lanczos_config(:, 1:TotWalkers)
                call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, TotWalkers)
            end if
        end if

    end subroutine create_initial_config

    subroutine generate_init_config_basic(nwalkers, ndets)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        integer, intent(in) :: nwalkers
        integer, intent(out) :: ndets
        integer :: i, irun, excit, nattempts, DetHash
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_sign(lenof_sign)
        logical :: tInitiatorTemp
        type(fcimc_iter_data) :: temp_data
        integer(n_int), pointer :: PointTemp(:,:)

        ! Turn off the initiator method for the annihilation steps to be used here.
        tInitiatorTemp = tTruncInitiator
        tTruncInitiator = .false.

        ! Set the spawning slots to their starting positions.
        ValidSpawnedList = InitialSpawnedSlots

        ilut = 0

        do irun = 1, inum_runs
            walker_sign = 0.0_dp
            do i = 1, nwalkers
                ! Generate the determinant (output to ilut).
                call CreateRandomExcitLevDetUnbias(nel, HFDet, ilutHF, ilut, excit, nattempts)
                call decode_bit_det(nI, ilut)

                ! Choose whether the walker to be added has an amplitude of plus or minus one, with
                ! 0.5 chance of each.
                walker_sign(irun) = 1.0_dp
                r = genrand_real2_dSFMT()
                if (r < 0.5) walker_sign(irun) = -1.0_dp*walker_sign(irun)

                call create_particle(nI, ilut, walker_sign, 0, irun)

            end do
        end do

        ! Perform annihilation steps:
        ! Send the walkers to their correct processors. The resulting walkers will be stored in
        ! SpawnedParts2.
        call SendProcNewParts(ndets, tSingleProc = .false.)
        ! CompressSpawnedList works on SpawnedParts, not SpawnedParts2, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp
        call CompressSpawnedList(ndets, temp_data) 

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        CurrentDets = SpawnedParts

        ! Add the entries into the hash table.
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets)

        ValidSpawnedList = InitialSpawnedSlots
        SpawnedParts = 0

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

    end subroutine generate_init_config_basic

end module stoch_lanczos_procs
