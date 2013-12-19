#include  "macros.h"

module stochastic_lanczos

    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
    use bit_rep_data, only: NIfTot
    use bit_reps, only: decode_bit_det
    use CalcData, only: tTruncInitiator
    use constants
    use dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData, only: ilutHF, HFDet, CurrentDets, SpawnedParts, SpawnedParts2, TotWalkers
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots, HashIndex, nWalkerHashes
    use FciMCData, only: fcimc_iter_data, ll_node
    use FciMCParMod, only: create_particle
    use hash, only: FindWalkerHash
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias
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
        ! The number of separate Lanczos calculations to perform for each
        ! different initial walker configuration.
        integer :: nruns
        ! The number of different Lanczos vectors to sample (the number of
        ! vectors which form the Krylov subspace at the end of a calculation).
        integer :: nkrylov_vecs
        ! The number of iterations to perform *between each Lanczos vector is
        ! sampled*.
        integer :: niters
    end type

    type(stoch_lanczos_data) :: lanczos

contains

    subroutine stoch_lanczos_read_inp()

        use input_neci

        logical :: eof
        character(len=100) :: w

        read_inp: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)
            case("end-lanczos")
                exit read_inp
            case("ground-state")
                lanczos%tGround = .true.
            case("finite-temperature")
                lanczos%tFiniteTemp = .true.
            case("num-init-configs")
                call geti(lanczos%nconfigs)
            case("num-runs-per-config")
                call geti(lanczos%nruns)
            case("num-lanczos-vecs")
                call geti(lanczos%nkrylov_vecs)
            case("num-iters-per-vec")
                call geti(lanczos%niters)
            case default
                call report("Keyword "//trim(w)//" not recognized in stoch-lanczos block", .true.)
            end select
        end do read_inp

    end subroutine stoch_lanczos_read_inp

    subroutine generate_init_config_basic(nwalkers)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        integer, intent(in) :: nwalkers
        integer :: i, irun, excit, nattempts, DetHash, ndets
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_sign(lenof_sign)
        logical :: tInitiatorTemp
        type(ll_node), pointer :: temp_node
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

        TotWalkers = int(ndets, int64)

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        CurrentDets = SpawnedParts

        ! Add the entries into the hash table.
        do i = 1, ndets
            call decode_bit_det(nI, CurrentDets(:,i))
            DetHash = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(DetHash)
            ! If the first element in the list has not been used.
            if (temp_node%ind == 0) then
                temp_node%ind = i
            else
                do while (associated(temp_node%next))
                    temp_node => temp_node%next
                end do
                allocate(temp_node%next)
                nullify(temp_node%next%next)
                temp_node%next%ind = i
            end if
            nullify(temp_node)
        end do

        ValidSpawnedList = InitialSpawnedSlots
        SpawnedParts = 0

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

    end subroutine generate_init_config_basic

end module stochastic_lanczos
