#include "macros.h"
 
module kp_fciqmc_init
 
    use bit_rep_data
    use bit_reps, only: decode_bit_det, encode_sign
    use constants
    use Parallel_neci, only: iProcIndex, MPISum, MPISumAll, nProcessors
    use kp_fciqmc_data_mod

    implicit none

contains

    subroutine kp_fciqmc_read_inp(kp)

        use CalcData, only: tWritePopsNorm, iPopsFileNoRead
        use FciMCData, only: alloc_popsfile_dets
        use input_neci
        use LoggingData, only: tIncrementPops
        use perturbations, only: init_perturbation_creation, init_perturbation_annihilation
        use SystemData, only: tAllSymSectors, tSpn

        type(kp_fciqmc_data), intent(inout) :: kp
        logical :: eof
        character(len=100) :: w
        integer :: i, j, niters_temp, nvecs_temp, nreports_temp, npert, nexcit
        integer :: ras_size_1, ras_size_2, ras_size_3, ras_min_1, ras_max_3
        character (len=*), parameter :: t_r = "kp_fciqmc_read_inp"

        ! Default values.
        tExcitedStateKP = .false.
        kp%nconfigs = 1
        kp%nreports = 0
        kp%nrepeats = 1
        kp%nvecs = 0
        niters_temp = 1
        memory_factor_kp = 1.0_dp

        tFiniteTemp = .false.
        tExcitedInitState = .false.

        tCalcSpin = .false.

        nwalkers_per_site_init = 1.0_dp
        av_mc_excits_kp = 0.0_dp
        kp_hamil_exact_frac = 1.0_dp

        tMultiplePopStart = .false.
        tExactHamil = .false.
        tExactHamilSpawning = .false.
        tFullyStochasticHamil = .false.
        tInitCorrectNWalkers = .false.
        tOccDetermInit = .false.
        vary_niters = .false.
        tUseInitConfigSeeds = .false.
        tAllSymSectors = .false.
        tOutputAverageKPMatrices = .false.
        tOverlapPert = .false.
        tScalePopulation = .false.
        scaling_factor = 1.0_dp

        n_kp_pops = 0
        Occ_KP_CasOrbs = 0
        Virt_KP_CasOrbs = 0
        kp_mp1_ndets = 0

        read_inp: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)
            case("END-KP-FCIQMC")
                exit read_inp
            case("EXCITED-STATE-KP")
#ifdef __PROG_NUMRUNS
                tExcitedStateKP = .true.
                call geti(kp%nvecs)
#else
                call stop_all(t_r, "NECI must be compiled with multiple replicas (mneci) to use the excited-state-kp option.")
#endif
            case("FINITE-TEMPERATURE")
                tFiniteTemp = .true.

            case("CALC-SPIN")
                tCalcSpin = .true.

            case("MULTIPLE-POPS")
                tMultiplePopStart = .true.
                tIncrementPops = .true.
                iPopsFileNoRead = -1
            case("NUM-INIT-CONFIGS")
                call geti(kp%nconfigs)
            case("NUM-REPEATS-PER-INIT-CONFIG")
                call geti(kp%nrepeats)
            case("NUM-KRYLOV-VECS")
                call geti(nvecs_temp)
                if (kp%nvecs /= 0) then
                    if (kp%nvecs /= nvecs_temp) call stop_all(t_r, 'The number of values specified for the number of iterations &
                        &between Krylov vectors is not consistent with the number of Krylov vectors requested.')
                else
                    kp%nvecs = nvecs_temp
                    allocate(kp%niters(kp%nvecs))
                    kp%niters(kp%nvecs) = 0
                end if
            case("NREPORTS")
                call geti(nreports_temp)
                if (kp%nreports /= 0) then
                    if (kp%nreports /= nreports_temp) call stop_all(t_r, 'The number of values specified for the number of &
                        &iterations between reports is not consistent with the number of reports requested.')
                else
                    kp%nreports = nreports_temp
                    allocate(kp%niters(kp%nreports))
                    kp%niters(kp%nreports) = 0
                end if
            case("NUM-ITERS-BETWEEN-VECS")
                call geti(niters_temp)
            case("NUM-ITERS-BETWEEN-REPORTS")
                call geti(niters_temp)
            case("NUM-ITERS-BETWEEN-VECS-VARY")
                vary_niters = .true.
                if (kp%nvecs /= 0) then
                    if (kp%nvecs /= nitems) call stop_all(t_r, 'The number of values specified for the number of iterations &
                        &between Krylov vectors is not consistent with the number of Krylov vectors requested.')
                end if
                kp%nvecs = nitems
                if (.not. allocated(kp%niters)) allocate(kp%niters(kp%nvecs))
                do i = 1, kp%nvecs-1
                    call geti(kp%niters(i))
                end do
                kp%niters(kp%nvecs) = 0
            case("NUM-ITERS-BETWEEN-REPORTS-VARY")
                vary_niters = .true.
                if (kp%nreports /= 0) then
                    if (kp%nreports /= nitems) call stop_all(t_r, 'The number of values specified for the number of &
                        &iterations between reports is not consistent with the number of reports requested.')
                end if
                kp%nreports = nitems
                if (.not. allocated(kp%niters)) allocate(kp%niters(kp%nreports))
                do i = 1, kp%nreports-1
                    call geti(kp%niters(i))
                end do
                kp%niters(kp%nreports) = 0
            case("MEMORY-FACTOR")
                call getf(memory_factor_kp)
            case("NUM-WALKERS-PER-SITE-INIT")
                call getf(nwalkers_per_site_init)
            case("AVERAGEMCEXCITS-HAMIL")
                call getf(av_mc_excits_kp)
            case("EXACT-HAMIL-SPAWNING")
                tExactHamilSpawning = .true.
            case("EXACT-HAMIL-SPAWNING-FRAC")
                call getf(kp_hamil_exact_frac)
            case("EXACT-HAMIL")
                tExactHamil = .true.
            case("FULLY-STOCHASTIC-HAMIL")
                tFullyStochasticHamil = .true.
            case("INIT-CORRECT-WALKER-POP")
                tInitCorrectNWalkers = .true.
            case("OCCUPY-DETERM-INIT")
                tOccDetermInit = .true.
            case("INIT-CONFIG-SEEDS")
                if (kp%nconfigs == 0) call stop_all(t_r, 'Please use the num-init-configs option to enter the number of &
                                                          &initial configurations before using the init-vec-seeds option.')
                tUseInitConfigSeeds = .true.
                allocate(init_config_seeds(kp%nconfigs))
                do i = 1, kp%nconfigs
                    call geti(init_config_seeds(i))
                end do
            case("ALL-SYM-SECTORS")
                tAllSymSectors = .true.
            case("WRITE-AVERAGE-MATRICES")
                tOutputAverageKPMatrices = .true.
            case("SCALE-POPULATION")
                tScalePopulation = .true.
            case("EXCITED-INIT-STATE")
                tExcitedInitState = .true.
                ! Read in the number of excited states to be used.
                call readi(nexcit)
                allocate(kpfciqmc_ex_labels(nexcit))
                allocate(kpfciqmc_ex_weights(nexcit))

                ! Read in which excited states to use.
                call read_line(eof)
                do j = 1, nitems
                    call readi(kpfciqmc_ex_labels(j))
                end do

                ! Read in the relative weights for the trial excited states in
                ! the initial state.
                call read_line(eof)
                do j = 1, nitems
                    call getf(kpfciqmc_ex_weights(j))
                end do
                ! Normalise weights so that they sum to 1.
                kpfciqmc_ex_weights = kpfciqmc_ex_weights/sum(kpfciqmc_ex_weights)

            case("OVERLAP-PERTURB-ANNIHILATE")
                alloc_popsfile_dets = .true.
                tOverlapPert = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npert)
                if (.not. allocated(overlap_pert)) then
                    allocate(overlap_pert(npert))
                else
                    if (npert /= size(overlap_pert)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert
                    call read_line(eof)
                    overlap_pert(i)%nannihilate = nitems
                    allocate(overlap_pert(i)%ann_orbs(nitems))
                    do j = 1, nitems
                        call readi(overlap_pert(i)%ann_orbs(j))
                    end do
                    ! Create the rest of the annihilation-related
                    ! components of the pops_pert object.
                    call init_perturbation_annihilation(overlap_pert(i))
                end do

            case("OVERLAP-PERTURB-CREATION")
                alloc_popsfile_dets = .true.
                tOverlapPert = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                call readi(npert)
                if (.not. allocated(overlap_pert)) then
                    allocate(overlap_pert(npert))
                else
                    if (npert /= size(overlap_pert)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert
                    call read_line(eof)
                    overlap_pert(i)%ncreate = nitems
                    allocate(overlap_pert(i)%crtn_orbs(nitems))
                    do j = 1, nitems
                        call readi(overlap_pert(i)%crtn_orbs(j))
                    end do
                    ! Create the rest of the creation-related
                    ! components of the pops_pert object.
                    call init_perturbation_creation(overlap_pert(i))
                end do

            case("DOUBLES-TRIAL")
                kp_trial_space_in%tDoubles = .true.
            case("CAS-TRIAL")
                kp_trial_space_in%tCAS = .true.
                tSpn = .true.
                call geti(Occ_KP_CASOrbs)  ! Number of electrons in the CAS 
                call geti(Virt_KP_CASOrbs) ! Number of virtual spin-orbitals in the CAS
            case("RAS-TRIAL")
                kp_trial_space_in%tRAS = .true.
                call geti(ras_size_1)  ! Number of spatial orbitals in RAS1.
                call geti(ras_size_2)  ! Number of spatial orbitals in RAS2.
                call geti(ras_size_3)  ! Number of spatial orbitals in RAS3.
                call geti(ras_min_1)  ! Min number of electrons (alpha and beta) in RAS1 orbs. 
                call geti(ras_max_3)  ! Max number of electrons (alpha and beta) in RAS3 orbs.
                kp_ras%size_1 = int(ras_size_1,sp)
                kp_ras%size_2 = int(ras_size_2,sp)
                kp_ras%size_3 = int(ras_size_3,sp)
                kp_ras%min_1 = int(ras_min_1,sp)
                kp_ras%max_3 = int(ras_max_3,sp)
            case("MP1-TRIAL")
                kp_trial_space_in%tMP1 = .true.
                call geti(kp_mp1_ndets)
            case("HF-TRIAL")
                kp_trial_space_in%tHF = .true.
            case("POPS-TRIAL")
                kp_trial_space_in%tPops = .true.
                call geti(n_kp_pops)
            case("READ-TRIAL")
                kp_trial_space_in%tRead = .true.
            case("FCI-TRIAL")
                kp_trial_space_in%tFCI = .false.
            case default
                call report("Keyword "//trim(w)//" not recognized in kp-fciqmc block", .true.)
            end select
        end do read_inp

        if ( (.not. vary_niters)) then
            if (tExcitedStateKP) then
                if (.not. allocated(kp%niters)) allocate(kp%niters(kp%nreports))
                kp%niters = niters_temp
                kp%niters(kp%nreports) = 0
            else
                if (.not. allocated(kp%niters)) allocate(kp%niters(kp%nvecs))
                kp%niters = niters_temp
                kp%niters(kp%nvecs) = 0
            end if
        end if

    end subroutine kp_fciqmc_read_inp

    subroutine init_kp_fciqmc(kp)

        use CalcData, only: tSemiStochastic, tUseRealCoeffs, AvMCExcits
        use fcimc_initialisation, only: SetupParameters, InitFCIMCCalcPar, init_fcimc_fn_pointers
        use FciMCData, only: tPopsAlreadyRead, nWalkerHashes, SpawnVecKP
        use FciMCData, only: SpawnVecKP2, MaxSpawned, determ_space_size, determ_sizes
        use FciMCData, only: SpawnedPartsKP, SpawnedPartsKP2, MaxWalkersUncorrected
        use FciMCData, only: iter_data_fciqmc, spawn_ht, nhashes_spawn, tReplicaReferencesDiffer
        use FciMCParMod, only: WriteFciMCStatsHeader, write_fcimcstats2, tSinglePartPhase
        use hash, only: init_hash_table
        use LoggingData, only: tFCIMCStats2
        use Parallel_neci, only: MPIBarrier
        use SystemData, only: G1, tHeisenberg, tAllSymSectors
        use util_mod, only: int_fmt

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: ierr, krylov_vecs_mem, krylov_ht_mem, matrix_mem, spawn_ht_mem
        character (len=*), parameter :: t_r = "init_kp_fciqmc"

        ! Checks.
        if (.not. tUseRealCoeffs) call stop_all(t_r,'kp-fciqmc can only be run using &
            &real coefficients).')
        if (tExactHamil .and. nProcessors /= 1) call stop_all(t_r,'The exact-hamil &
            &option can only be used when running with one processor.')
        if (theisenberg .and. tAllSymSectors) call stop_all(t_r,'The option to use all &
            &symmetry sectors at once has not been implemented with the Heisenberg model.')
        if (n_int == 4) call stop_all(t_r, 'Use of RealCoefficients does not work with 32 bit &
             &integers due to the use of the transfer operation from dp reals to 64 bit integers.')
        if (tExcitedStateKP .and. inum_runs /= 2*kp%nvecs) call stop_all(t_r, 'When using the &
             &ExcitedStateKP option the number of replicas must be twice the number of states to &
             &be calculated.')

        ! Call external NECI initialisation routines.
        tPopsAlreadyRead = .false.
        call SetupParameters()
        call InitFCIMCCalcPar()
        call init_fcimc_fn_pointers() 

        write(6,'(/,12("="),1x,a9,1x,12("="))') "KP-FCIQMC"

        if (tSemiStochastic .and. (.not. tFullyStochasticHamil)) then
            tSemiStochasticKPHamil = .true.
        else
            tSemiStochasticKPHamil = .false.
        end if

        ! The number of elements required to store all replicas of all Krylov vectors.
        lenof_all_signs = lenof_sign_kp*kp%nvecs
        ! The total length of a bitstring containing all Krylov vectors.
        NIfTotKP = NIfDBO + lenof_all_signs + NIfFlag

        if (.not. tExcitedStateKP) then
            ! Allocate all of the KP arrays.
            nhashes_kp = nWalkerHashes
            TotWalkersKP = 0
            krylov_vecs_length = nint(MaxWalkersUncorrected*memory_factor_kp*kp%nvecs)
            nkrylov_amp_elems_tot = lenof_sign*kp%nvecs*krylov_vecs_length

            ! Allocate the krylov_vecs array.
            ! The number of MB of memory required to allocate krylov_vecs.
            krylov_vecs_mem = krylov_vecs_length*(NIfTotKP+1)*size_n_int/1000000
            write(6,'(a73,'//int_fmt(krylov_vecs_mem,1)//')') "About to allocate array to hold all Krylov vectors. &
                                           &Memory required (MB):", krylov_vecs_mem
            write(6,'(a13)',advance='no') "Allocating..."; call neci_flush(6)
            allocate(krylov_vecs(0:NIfTotKP, krylov_vecs_length), stat=ierr)
            if (ierr /= 0) then
                write(6,'(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating krylov_vecs array.")
            else
                write(6,'(1x,a5)') "Done."
            end if
            call neci_flush(6)
            krylov_vecs = 0_n_int

            ! Allocate the krylov_helems array.
            ! The number of MB of memory required to allocate krylov_helems.
            krylov_vecs_mem = krylov_vecs_length*size_n_int/1000000
            write(6,'(a103,'//int_fmt(krylov_vecs_mem,1)//')') "About to allocate array to hold diagonal Hamiltonian &
                                           &elements for Krylov vectors. Memory required (MB):", krylov_vecs_mem
            write(6,'(a13)',advance='no') "Allocating..."; call neci_flush(6)
            allocate(krylov_helems(krylov_vecs_length), stat=ierr)
            if (ierr /= 0) then
                write(6,'(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating krylov_helems array.")
            else
                write(6,'(1x,a5)') "Done."
            end if
            call neci_flush(6)
            krylov_helems = 0.0_dp

            ! Allocate the hash table to krylov_vecs.
            ! The number of MB of memory required to allocate krylov_vecs_ht.
            ! Each node requires 16 bytes.
            krylov_ht_mem = nhashes_kp*16/1000000
            write(6,'(a78,'//int_fmt(krylov_ht_mem,1)//')') "About to allocate hash table to the Krylov vector array. &
                                           &Memory required (MB):", krylov_ht_mem
            write(6,'(a13)',advance='no') "Allocating..."; call neci_flush(6)
            allocate(krylov_vecs_ht(nhashes_kp), stat=ierr)
            if (ierr /= 0) then
                write(6,'(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating krylov_vecs_ht array.")
            else
                write(6,'(1x,a5)') "Done."
                write(6,'(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                                  &increase as further nodes are added."
            end if

            call init_hash_table(krylov_vecs_ht)

            allocate(SpawnVecKP(0:NOffSgn+lenof_all_signs-1,MaxSpawned),stat=ierr)
            allocate(SpawnVecKP2(0:NOffSgn+lenof_all_signs-1,MaxSpawned),stat=ierr)
            SpawnVecKP(:,:) = 0_n_int
            SpawnVecKP2(:,:) = 0_n_int
            SpawnedPartsKP => SpawnVecKP
            SpawnedPartsKP2 => SpawnVecKP2

            if (tSemiStochastic) then
                allocate(partial_determ_vecs_kp(lenof_all_signs,determ_sizes(iProcIndex)), stat=ierr)
                allocate(full_determ_vecs_kp(lenof_all_signs,determ_space_size), stat=ierr)
                partial_determ_vecs_kp = 0.0_dp
                full_determ_vecs_kp = 0.0_dp
            end if

        end if

        ! Allocate the hash table to the spawning array.
        ! The number of MB of memory required to allocate spawn_ht.
        ! Each node requires 16 bytes.
        nhashes_spawn = 0.8*MaxSpawned
        spawn_ht_mem = nhashes_spawn*16/1000000
        write(6,'(a78,'//int_fmt(spawn_ht_mem,1)//')') "About to allocate hash table to the spawning array. &
                                       &Memory required (MB):", spawn_ht_mem
        write(6,'(a13)',advance='no') "Allocating..."; call neci_flush(6)
        allocate(spawn_ht(nhashes_spawn), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating spawn_ht array.")
        else
            write(6,'(1x,a5)') "Done."
            write(6,'(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                              &increase as further nodes are added."
        end if

        call init_hash_table(spawn_ht)

        ! (2*kp%nrepeats+16) arrays with (kp%nvecs**2) 8-byte elements each.
        matrix_mem = (2*kp%nrepeats+16)*(kp%nvecs**2)*8/1000000
        write(6,'(a66,'//int_fmt(matrix_mem,1)//')') "About to allocate various subspace matrices. &
                                       &Memory required (MB):", matrix_mem
        write(6,'(a13)',advance='no') "Allocating..."
        call neci_flush(6)

        if (tOverlapPert) then
            allocate(pert_overlaps(kp%nvecs))
            allocate(kp_all_pert_overlaps(kp%nvecs))
        end if

        allocate(kp_hamil_mean(kp%nvecs, kp%nvecs))
        allocate(kp_overlap_mean(kp%nvecs, kp%nvecs))
        allocate(kp_hamil_se(kp%nvecs, kp%nvecs))
        allocate(kp_overlap_se(kp%nvecs, kp%nvecs))

        allocate(kp_overlap_eigv(kp%nvecs))
        allocate(kp_init_overlaps(kp%nvecs))
        allocate(kp_overlap_eigenvecs(kp%nvecs, kp%nvecs))
        allocate(kp_transform_matrix(kp%nvecs, kp%nvecs))
        allocate(kp_inter_matrix(kp%nvecs, kp%nvecs))
        allocate(kp_eigenvecs_krylov(kp%nvecs, kp%nvecs))

        write(6,'(1x,a5)') "Done."

        ! If av_mc_excits_kp hasn't been set by the user, just use AvMCExcits.
        if (av_mc_excits_kp == 0.0_dp) av_mc_excits_kp = AvMCExcits

        MaxSpawnedEachProc = int(0.88_dp*real(MaxSpawned,dp)/nProcessors)

        if (tFCIMCStats2) then
            call write_fcimcstats2(iter_data_fciqmc, initial=.true.)
        else
            call WriteFciMCStatsHeader()
        end if
        call MPIBarrier(ierr)

        ! Store the initial state of tSinglePartPhase so that we can stop the
        ! shift from varying on subsequent repeats.
        tSinglePartPhaseKPInit = tSinglePartPhase

        ! Allow different replicas to have different references in this case.
        if (tExcitedStateKP) tReplicaReferencesDiffer = .true.

    end subroutine init_kp_fciqmc

    subroutine init_kp_fciqmc_repeat(iconfig, irepeat, nrepeats, nvecs)

        use CalcData, only: tStartSinglePart, InitialPart, InitWalkers, DiagSft, iPopsFileNoRead
        use FciMCData, only: iter, InputDiagSft, PreviousCycles, OldAllAvWalkersCyc, proje_iter
        use FciMCData, only: proje_iter_tot, AllGrowRate, SpawnedParts
        use FciMCParMod, only: tSinglePartPhase
        use hash, only: clear_hash_table
        use initial_trial_states
        use util_mod, only: int_fmt

        integer, intent(in) :: iconfig, irepeat, nrepeats, nvecs

        integer :: ndets_this_proc, nexcit
        real(dp), allocatable :: evecs_this_proc(:,:), init_vecs(:,:)

        write(6,'(1x,a22,'//int_fmt(irepeat,1)//')') "Starting repeat number", irepeat

        if (tExcitedStateKP) then
            nexcit = nvecs

            ! Create the trial excited states.
            call calc_trial_states(kp_trial_space_in, nexcit, ndets_this_proc, evecs_this_proc, SpawnedParts)
            ! Set the populations of these states to the requested value.
            call set_trial_populations(nexcit, ndets_this_proc, evecs_this_proc)
            ! Set the trial excited states as the FCIQMC wave functions.
            call set_trial_states(ndets_this_proc, evecs_this_proc, SpawnedParts)

            deallocate(evecs_this_proc)

        else if (tExcitedInitState) then
            nexcit = maxval(kpfciqmc_ex_labels)

            ! Create the trial excited states.
            call calc_trial_states(kp_trial_space_in, nexcit, ndets_this_proc, evecs_this_proc, SpawnedParts)
            ! Extract the desried initial excited states and average them.
            call create_init_excited_state(ndets_this_proc, evecs_this_proc, kpfciqmc_ex_labels, kpfciqmc_ex_weights, init_vecs)
            ! Set the populations of these states to the requested value.
            call set_trial_populations(1, ndets_this_proc, init_vecs)
            ! Set the trial excited states as the FCIQMC wave functions.
            call set_trial_states(ndets_this_proc, init_vecs, SpawnedParts)

            deallocate(evecs_this_proc, init_vecs)

        else
            ! If starting from multiple POPSFILEs then set this counter so that the
            ! correct POPSFILE is read in this time. To read in POPSFILE.x,
            ! iPopsFileNoRead needs to be set to -x-1. We want to read in POPSFILE
            ! numbers 0 to kp%nconfigs-1
            if (tMultiplePopStart) iPopsFileNoRead = -(iconfig-1)-1

            if (tOverlapPert .and. irepeat == 1) then
                pert_overlaps = 0.0_dp
                call create_overlap_pert_vec()
            end if

            call create_initial_config(iconfig, irepeat, nrepeats)

            call clear_hash_table(krylov_vecs_ht)
            krylov_vecs = 0_n_int
        end if

        ! Rezero all the necessary data.
        TotWalkersKP = 0
        nkrylov_amp_elems_used = 0
        iter = 0
        PreviousCycles = 0
        DiagSft = InputDiagSft
        if (tStartSinglePart) then
            OldAllAvWalkersCyc = InitialPart
        else
            OldAllAvWalkersCyc = InitWalkers*nProcessors
        end if
        proje_iter = 0.0_dp
        proje_iter_tot = 0.0_dp
        AllGrowRate = 0.0_dp

        ! Setting this variable to true stops the shift from varying instantly.
        tSinglePartPhase = tSinglePartPhaseKPInit

    end subroutine init_kp_fciqmc_repeat

    subroutine init_kp_fciqmc_iter(iter_data, determ_index)

        use FciMCData, only: FreeSlot, iStartFreeSlot, iEndFreeSlot, fcimc_iter_data, InitialSpawnedSlots
        use FciMCData, only: ValidSpawnedList, spawn_ht
        use FciMCParMod, only: rezero_iter_stats_each_iter
        use hash, only: clear_hash_table

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

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        call rezero_iter_stats_each_iter(iter_data)

    end subroutine init_kp_fciqmc_iter

    subroutine create_initial_config(iconfig, irepeat, nrepeats)

        use CalcData, only: tStartSinglePart, InitialPart, InitWalkers, tSemiStochastic, tReadPops
        use dSFMT_interface , only: dSFMT_init
        use fcimc_initialisation, only: InitFCIMC_HF
        use FciMCData, only: nWalkerHashes, HashIndex, pops_pert, SpawnedParts, TotWalkers, AllTotWalkers
        use FciMCData, only: TotParts, AllTotParts, TotPartsOld, AllTotPartsOld, kp_generate_time
        use FciMCData, only: tStartCoreGroundState, CurrentDets, HolesInList
        use hash, only: fill_in_hash_table, clear_hash_table
        use PopsfileMod, only: read_popsfile_wrapper
        use semi_stoch_procs, only: copy_core_dets_to_spawnedparts, fill_in_diag_helements
        use semi_stoch_procs, only: add_core_states_currentdet_hash, start_walkers_from_core_ground
        use semi_stoch_procs, only: check_determ_flag
        use timing_neci, only: get_total_time, set_timer, halt_timer

        integer, intent(in) :: iconfig, irepeat, nrepeats

        integer :: DetHash, nwalkers_int
        integer :: i
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign(lenof_sign_kp), TotPartsCheck(lenof_sign_kp)
        real(dp) :: nwalkers_target
        real(dp) :: norm, all_norm
        real(sp) :: total_time_before, total_time_after
        logical :: tCoreDet
        character(len=*), parameter :: t_r = "create_init_config"

        ! Clear everything from any previous repeats or starting configurations.
        call clear_hash_table(HashIndex)

        if (tStartSinglePart) then
            nwalkers_target = real(InitialPart,dp)
        else
            nwalkers_target = InitWalkers*nProcessors
        end if

        if (irepeat == 1) then
            if (.not. tFiniteTemp) then
                if (tReadPops) then
                    ! Call a wrapper function which will call the various functions
                    ! required to read in a popsfile.
                    if (allocated(pops_pert)) then
                        call read_popsfile_wrapper(pops_pert)
                    else
                        call read_popsfile_wrapper()
                    end if

                    if (tScalePopulation) then
                        call scale_population(CurrentDets, TotWalkers, nwalkers_target, TotPartsCheck, scaling_factor)
                        ! Update global data.
                        if (any(abs(TotPartsCheck-TotParts) > 1.0e-12_dp)) then
                            call stop_all(t_r, "Inconsistent values of TotParts calculated.")
                        end if
                        TotParts = TotParts*scaling_factor
                        TotPartsOld = TotParts
                        AllTotParts = AllTotParts*scaling_factor
                        AllTotPartsOld = AllTotParts
                    end if
                else
                    ! Put a walker on the Hartree-Fock, with the requested amplitude.
                    call InitFCIMC_HF()
                end if

                if (tSemiStochastic) then
                    ! core_space stores all core determinants from all processors. Move those on this
                    ! processor to SpawnedParts, which add_core_states_currentdet_hash uses.
                    call copy_core_dets_to_spawnedparts()
                    ! Any core space determinants which are not already in CurrentDets will be added
                    ! by this routine.
                    call add_core_states_currentdet_hash()
                    if (tStartCoreGroundState .and. (.not. tReadPops)) &
                        call start_walkers_from_core_ground(tPrintInfo = .false.)
                end if

            else if (tFiniteTemp) then
                ! Convert the initial number of walkers to an integer. Note that on multiple
                ! processors this may round up the requested number of walkers slightly.
                nwalkers_int = ceiling(nwalkers_target/real(nProcessors,dp))

                ! If requested, reset the random number generator with the requested seed
                ! before creating the random initial configuration.
                if (tUseInitConfigSeeds) call dSFMT_init((iProcIndex+1)*init_config_seeds(iconfig))

                write (6,'(a44)',advance='no') "# Generating initial walker configuration..."
                call set_timer(kp_generate_time)
                total_time_before = get_total_time(kp_generate_time)

                ! Create the random initial configuration. 
                if (tInitCorrectNWalkers) then
                    call generate_init_config_this_proc(nwalkers_int, nwalkers_per_site_init, tOccDetermInit)
                else
                    call generate_init_config_basic(nwalkers_int, nwalkers_per_site_init)
                end if

                call halt_timer(kp_generate_time)
                total_time_after = get_total_time(kp_generate_time)
                write(6,'(1x,a31,f9.3)') "Complete. Time taken (seconds):", total_time_after-total_time_before

            end if
        else if (irepeat > 1) then
            ! If repeating from a previsouly generated initial configuration, simpy reset the following
            ! data and copy the first Krylov vector (which is always the starting configuration) from
            ! the last run to CurrentDets.
            TotWalkers = TotWalkersInit
            AllTotWalkers = AllTotWalkersInit
            TotParts = TotPartsInit
            TotPartsOld = TotPartsInit
            AllTotParts = AllTotPartsInit
            AllTotPartsOld = AllTotPartsInit
            do i = 1, int(TotWalkers, sizeof_int)
                ! Copy across the bitstring encoding of the determinant and also the walker signs.
                CurrentDets(0:NOffSgn+lenof_sign_kp-1,i) = krylov_vecs(0:NOffSgn+lenof_sign_kp-1,i)
                ! Copy across the flags.
                CurrentDets(NIfTot,i) = krylov_vecs(NIfTotKP,i)
            end do
            call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, int(TotWalkers, sizeof_int), .true.)
        end if

        ! Calculate and store the diagonal element of the Hamiltonian for determinants in CurrentDets.
        call fill_in_diag_helements()

        ! If starting from this configuration more than once, store the relevant data for next time.
        if (nrepeats > 1 .and. irepeat == 1) then
            HolesInList = 0
            do i = 1, TotWalkers
                int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign_kp-1,i)
                call extract_sign (CurrentDets(:,i), real_sign)
                tCoreDet = check_determ_flag(CurrentDets(:,i))
                ! Don't add unoccpied determinants, unless they are core determinants.
                if (IsUnoccDet(real_sign) .and. (.not. tCoreDet)) HolesInList = HolesInList + 1
            end do
            TotWalkersInit = TotWalkers - HolesInList
            call MPISumAll(TotWalkersInit, AllTotWalkers)
            AllTotWalkersInit = AllTotWalkers
            TotPartsInit = TotParts
            AllTotPartsInit = AllTotParts
        end if

    end subroutine create_initial_config

    subroutine create_overlap_pert_vec()

        ! Read in the popsfile and apply perturbation operator overlap_pert.

        use FciMCData, only: TotWalkers, CurrentDets
        use PopsfileMod, only: read_popsfile_wrapper
        use util_mod, only: int_fmt

        integer :: mem_reqd, ierr
        character(len=*), parameter :: t_r = "create_overlap_pert_vec"

        if (allocated(perturbed_ground)) deallocate(perturbed_ground)

        ! Once this is finished, the vector that we want will be stored in
        ! CurrentDets. The total number of determinants will be TotWalkers.
        call read_popsfile_wrapper(overlap_pert)

        ! Print info about memory usage to the user.
        ! Memory required in MB.
        mem_reqd = TotWalkers*(NIfTotKP+1)*size_n_int/1000000

        write(6,'(a73,'//int_fmt(mem_reqd,1)//')') "About to allocate array to hold the perturbed &
                                           &ground state. Memory required (MB):", mem_reqd
        write(6,'(a13)',advance='no') "Allocating..."
        call neci_flush(6)
        allocate(perturbed_ground(0:NIfTot,TotWalkers), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array.")
        else
            write(6,'(1x,a5)') "Done."
        end if
        call neci_flush(6)

        perturbed_ground = CurrentDets(0:NIfTot,1:TotWalkers)

    end subroutine create_overlap_pert_vec

    subroutine generate_init_config_basic(nwalkers, nwalkers_per_site_init)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList
        use CalcData, only: tTruncInitiator, tSemiStochastic
        use dSFMT_interface, only: genrand_real2_dSFMT
        use fcimc_helper, only: create_particle
        use FciMCData, only: nWalkerHashes, fcimc_iter_data, HashIndex, SpawnedParts, SpawnedParts2
        use FciMCData, only: TotWalkers, CurrentDets, InitialSpawnedSlots, TotParts, TotPartsOld
        use FciMCData, only: AllTotParts, AllTotPartsOld, ValidSpawnedList, ilutHF, HFDet
        use hash, only: fill_in_hash_table
        use hilbert_space_size, only: CreateRandomExcitLevDetUnbias, create_rand_heisenberg_det
        use hilbert_space_size, only: create_rand_det_no_sym
        use replica_data, only: allocate_iter_data, clean_iter_data
        use semi_stoch_procs, only: copy_core_dets_to_spawnedparts, add_core_states_currentdet_hash
        use SystemData, only: nel, tHeisenberg, tAllSymSectors

        integer, intent(in) :: nwalkers
        real(dp) :: nwalkers_per_site_init
        integer :: i, ireplica, excit, nattempts, DetHash
        integer :: nspawns, ndets
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_amp, walker_sign(lenof_sign_kp)
        logical :: tInitiatorTemp
        type(fcimc_iter_data) :: unused_data
        integer(n_int), pointer :: PointTemp(:,:)

        call allocate_iter_data(unused_data)

        ! Turn off the initiator method for the annihilation steps to be used here.
        tInitiatorTemp = tTruncInitiator
        tTruncInitiator = .false.

        ! Set the spawning slots to their starting positions.
        ValidSpawnedList = InitialSpawnedSlots

        ilut = 0_n_int
        nspawns = ceiling(real(nwalkers,dp)/nwalkers_per_site_init)

        do i = 1, nspawns
            ! Generate a determinant (output to ilut).
            if (tAllSymSectors) then
                call create_rand_det_no_sym(ilut)
            else
                if (tHeisenberg) then
                    call create_rand_heisenberg_det(ilut)
                else
                    call CreateRandomExcitLevDetUnbias(nel, HFDet, ilutHF, ilut, excit, nattempts)
                end if
            end if
            call decode_bit_det(nI, ilut)

            ! Choose whether the walker should have a positive or negative amplitude, with
            ! 50% chance of each.
            walker_amp = nwalkers_per_site_init
            r = genrand_real2_dSFMT()
            if (r < 0.5_dp) walker_amp = -1.0_dp*walker_amp

            do ireplica = 1, inum_runs
                walker_sign = 0.0_dp
                walker_sign(ireplica) = walker_amp
                call create_particle(nI, ilut, walker_sign, ireplica)
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
        call CompressSpawnedList(ndets, unused_data) 

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        TotParts = 0.0_dp
        do i = 1, ndets
            CurrentDets(:,i) = SpawnedParts(:,i)
            walker_sign = transfer(CurrentDets(NOffSgn:NOffSgn+lenof_sign_kp-1, i), walker_sign)
            TotParts = TotParts + abs(walker_sign)
        end do
        TotPartsOld = TotParts

        ! Add the entries into the hash table.
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets, .true.)

        call MPISum(TotParts, AllTotParts)
        AllTotPartsOld = AllTotParts
        TotWalkers = int(ndets, int64)

        if (tSemiStochastic) then
            ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
            ! These routines will do this.
            call copy_core_dets_to_spawnedparts()
            call add_core_states_currentdet_hash()
        end if

        ValidSpawnedList = InitialSpawnedSlots

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

        call clean_iter_data(unused_data)

    end subroutine generate_init_config_basic

    subroutine generate_init_config_this_proc(nwalkers, nwalkers_per_site_init, tOccupyDetermSpace)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        use CalcData, only: tSemiStochastic
        use dSFMT_interface, only: genrand_real2_dSFMT
        use FciMCData, only: HashIndex, determ_sizes, determ_displs, TotWalkers, CurrentDets, HFDet
        use FciMCData, only: TotParts, TotPartsOld, AllTotParts, AllTotPartsOld, core_space, ilutHF
        use hash, only: DetermineDetNode, rm_unocc_dets_from_hash_table, hash_table_lookup
        use hash, only: add_hash_table_entry
        use hilbert_space_size, only: CreateRandomExcitLevDetUnbias, create_rand_heisenberg_det
        use hilbert_space_size, only: create_rand_det_no_sym
        use semi_stoch_procs, only: copy_core_dets_to_spawnedparts, add_core_states_currentdet_hash
        use SystemData, only: nel, tHeisenberg, tAllSymSectors

        integer, intent(in) :: nwalkers
        real(dp) :: nwalkers_per_site_init
        logical, intent(in) :: tOccupyDetermSpace
        integer :: proc, excit, nattempts, hash_val, det_ind, nI(nel)
        integer :: ideterm, ndets
        integer(n_int) :: ilut(0:NIfTot), int_sign(lenof_sign_kp)
        real(dp) :: real_sign_1(lenof_sign_kp), real_sign_2(lenof_sign_kp)
        real(dp) :: new_sign(lenof_sign_kp)
        real(dp) :: r
        logical :: tDetermAllOccupied, tDetFound

        ilut = 0_n_int
        ndets = 0
        TotParts = 0.0_dp

        ideterm = 0
        tDetermAllOccupied = .false.

        do
            ! If using the tOccupyDetermSpace option then we want to put walkers
            ! on states in the deterministic space first.
            ! If not, or if we have finished doing this, generate determinants
            ! randomly and uniformly.
            if (tOccupyDetermSpace .and. (.not. tDetermAllOccupied)) then
                ideterm = ideterm + 1
                ilut = core_space(:,ideterm + determ_displs(iProcIndex))
                call decode_bit_det(nI, ilut)
                ! If we have now occupied all deterministic states.
                if (ideterm == determ_sizes(iProcIndex)) tDetermAllOccupied = .true.
            else
                ! Generate a random determinant (returned in ilut).
                if (tAllSymSectors) then
                    call create_rand_det_no_sym(ilut)
                else
                    if (tHeisenberg) then
                        call create_rand_heisenberg_det(ilut)
                    else
                        call CreateRandomExcitLevDetUnbias(nel, HFDet, ilutHF, ilut, excit, nattempts)
                    end if
                end if

                ! If it doesn't belong to this processor, pick another.
                call decode_bit_det(nI, ilut)
                proc = DetermineDetNode(nel, nI, 0)
                if (proc /= iProcIndex) cycle
            end if

            ! Choose whether the walker should have a positive or negative amplitude, with
            ! 50% chance of each.
            real_sign_1 = nwalkers_per_site_init
            r = genrand_real2_dSFMT()
            if (r < 0.5_dp) real_sign_1 = -1.0_dp*real_sign_1
            int_sign = transfer(real_sign_1, int_sign)

            tDetFound = .false.
            ! Search the hash table to see if this determinant is in CurrentDets
            ! already.
            call hash_table_lookup(nI, ilut, NIfDBO, HashIndex, CurrentDets, det_ind, hash_val, tDetFound)
            if (tDetFound) then
                ! This determinant is already in CurrentDets. Just need to add
                ! the sign of this new walker on and update stats.
                int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign_kp-1, det_ind)
                real_sign_2 = transfer(int_sign, real_sign_2)
                new_sign = real_sign_1 + real_sign_2
                call encode_sign(CurrentDets(:, det_ind), new_sign)
                TotParts = TotParts - abs(real_sign_2) + abs(new_sign)
            else
                ! A new determiant needs to be added.
                ndets = ndets + 1
                det_ind = ndets
                ! Add the new determinant to the hash table.
                call add_hash_table_entry(HashIndex, det_ind, hash_val)
                ! Copy determinant data across.
                CurrentDets(0:NIfDBO, det_ind) = ilut(0:NIfDBO)
                CurrentDets(NOffSgn:NOffSgn+lenof_sign_kp-1, det_ind) = int_sign
                if (tUseFlags) CurrentDets(NOffFlag, det_ind) = 0_n_int
                TotParts = TotParts + abs(real_sign_1)
            end if

            ! If we've reached the total requested number of walkers, finish.
            if (TotParts(1) >= nwalkers) exit

        end do

        call MPISum(TotParts, AllTotParts)
        TotPartsOld = TotParts
        AllTotPartsOld = AllTotParts
        TotWalkers = int(ndets, int64)

        if (tSemiStochastic) then
            ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
            ! These routines will do this.
            call copy_core_dets_to_spawnedparts()
            call add_core_states_currentdet_hash()
        else
            ! Some determinants may have become occupied and then unoccupied in
            ! the course of the above. We need to remove the entries for these
            ! determinants from the hash table. For semi-stochastic calculations
            ! this is done in add_core_states_currentdet_hash.
            call rm_unocc_dets_from_hash_table(HashIndex, CurrentDets, ndets)
        end if

    end subroutine generate_init_config_this_proc

    subroutine create_init_excited_state(ndets_this_proc, trial_vecs, ex_state_labels, ex_state_weights, init_vec)

        integer, intent(in) :: ndets_this_proc
        real(dp), intent(in) :: trial_vecs(:,:)
        integer, intent(in) :: ex_state_labels(:)
        real(dp), intent(in) :: ex_state_weights(:)
        real(dp), allocatable, intent(out) :: init_vec(:,:)

        real(dp) :: real_sign(lenof_sign)
        integer :: i, j, ierr
        character(len=*), parameter :: t_r = "create_init_excited_state"

        allocate(init_vec(1,ndets_this_proc), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error in MPIScatterV call.")
         
        do i = 1, ndets_this_proc
            init_vec(1,i) = 0.0_dp
            do j = 1, size(ex_state_labels)
                init_vec(1,i) = init_vec(1,i) + ex_state_weights(j)*trial_vecs(ex_state_labels(j), i)
            end do
        end do

    end subroutine create_init_excited_state

    subroutine scale_population(walker_list, ndets, target_pop, input_pop, scaling_factor)

        ! Take an input list of walkers, find the total walker population in
        ! the list, and then multiply all the walker signs by some factor in
        ! order for the walker list to have the requested target population.

        integer(n_int), intent(inout) :: walker_list(:,:)
        integer(int64), intent(in) :: ndets
        real(dp), intent(in) :: target_pop
        real(dp), intent(out) :: input_pop(lenof_sign_kp)
        real(dp), intent(out) :: scaling_factor

        integer :: i
        real(dp) :: real_sign(lenof_sign_kp), all_input_pop(lenof_sign_kp)

        input_pop = 0.0_dp

        ! First, find the population of the walkers in walker_list.
        do i = 1, ndets
            call extract_sign(walker_list(:,i), real_sign)
            input_pop = input_pop + abs(real_sign)
        end do

        call MPISumAll(input_pop, all_input_pop)

        ! Just use the first particle type to determine the scaing factor.
        scaling_factor = target_pop/all_input_pop(1)

        ! Now multiply the walker signs by scaling_factor.
        do i = 1, ndets
            call extract_sign(walker_list(:,i), real_sign)
            real_sign = real_sign*scaling_factor
            call encode_sign(walker_list(:,i), real_sign)
        end do

    end subroutine scale_population

end module kp_fciqmc_init
