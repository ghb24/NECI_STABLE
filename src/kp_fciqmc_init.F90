#include "macros.h"

module kp_fciqmc_init

    use bit_rep_data
    use bit_reps, only: decode_bit_det, encode_sign
    use constants
    use Parallel_neci, only: iProcIndex, MPISum, MPISumAll, nProcessors
    use kp_fciqmc_data_mod
    use util_mod, only: near_zero
    use FciMCData, only: core_run
    use core_space_util, only: cs_replicas
    use input_parser_mod, only: FileReader_t, TokenIterator_t
    use fortran_strings, only: to_upper, to_lower, to_int, to_realdp
    implicit none

contains

    subroutine kp_fciqmc_read_inp(file_reader, kp)

        use CalcData, only: tWritePopsNorm, iPopsFileNoRead, tPairedReplicas
        use FciMCData, only: alloc_popsfile_dets
        use LoggingData, only: tIncrementPops
        use perturbations, only: init_perturbation_creation, init_perturbation_annihilation
        use SystemData, only: tAllSymSectors, tSpn

        type(kp_fciqmc_data), intent(inout) :: kp
        class(FileReader_t), intent(inout) :: file_reader
        logical :: eof
        type(TokenIterator_t) :: tokens
        character(len=100) :: w
        integer :: i, j, niters_temp, nvecs_temp, nreports_temp, npert, nexcit
        integer :: ras_size_1, ras_size_2, ras_size_3, ras_min_1, ras_max_3
        character(len=*), parameter :: t_r = "kp_fciqmc_read_inp", this_routine = t_r

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

        tPairedReplicas = .true.
#if defined(PROG_NUMRUNS_)
        nreplicas = 2
#endif
        tOrthogKPReplicas = .false.
        orthog_kp_iter = 0

        read_inp: do while (file_reader%nextline(tokens, skip_empty=.true.))
            w = to_upper(tokens%next())
            select case (w)
            case ("END-KP-FCIQMC")
                exit read_inp
            case ("EXCITED-STATE-KP")
#ifdef PROG_NUMRUNS_
                tExcitedStateKP = .true.
                kp%nvecs = to_int(tokens%next())
#else
                call stop_all(t_r, "NECI must be compiled with multiple replicas (mneci) to use the excited-state-kp option.")
#endif
            case ("FINITE-TEMPERATURE")
                tFiniteTemp = .true.

            case ("CALC-SPIN")
                tCalcSpin = .true.

            case ("MULTIPLE-POPS")
                tMultiplePopStart = .true.
                tIncrementPops = .true.
                iPopsFileNoRead = -1
            case ("NUM-INIT-CONFIGS")
                kp%nconfigs = to_int(tokens%next())
            case ("NUM-REPEATS-PER-INIT-CONFIG")
                kp%nrepeats = to_int(tokens%next())
            case ("NUM-KRYLOV-VECS")
                nvecs_temp = to_int(tokens%next())
                if (kp%nvecs /= 0) then
                    if (kp%nvecs /= nvecs_temp) call stop_all(t_r, 'The number of values specified for the number of iterations &
                        &between Krylov vectors is not consistent with the number of Krylov vectors requested.')
                else
                    kp%nvecs = nvecs_temp
                    allocate(kp%niters(kp%nvecs))
                    kp%niters(kp%nvecs) = 0
                end if
            case ("NREPORTS")
                nreports_temp = to_int(tokens%next())
                if (kp%nreports /= 0) then
                    if (kp%nreports /= nreports_temp) call stop_all(t_r, 'The number of values specified for the number of &
                        &iterations between reports is not consistent with the number of reports requested.')
                else
                    kp%nreports = nreports_temp
                    allocate(kp%niters(kp%nreports))
                    kp%niters(kp%nreports) = 0
                end if
            case ("NUM-ITERS-BETWEEN-VECS")
                niters_temp = to_int(tokens%next())
            case ("NUM-ITERS-BETWEEN-REPORTS")
                niters_temp = to_int(tokens%next())
            case ("NUM-ITERS-BETWEEN-VECS-VARY")
                vary_niters = .true.
                if (kp%nvecs /= 0) then
                    if (kp%nvecs /= tokens%size()) call stop_all(t_r, 'The number of values specified for the number of iterations &
                        &between Krylov vectors is not consistent with the number of Krylov vectors requested.')
                end if
                kp%nvecs = tokens%size()
                if (.not. allocated(kp%niters)) allocate(kp%niters(kp%nvecs))
                do i = 1, kp%nvecs - 1
                    kp%niters(i) = to_int(tokens%next())
                end do
                kp%niters(kp%nvecs) = 0
            case ("NUM-ITERS-BETWEEN-REPORTS-VARY")
                vary_niters = .true.
                if (kp%nreports /= 0) then
                    if (kp%nreports /= tokens%size()) call stop_all(t_r, 'The number of values specified for the number of &
                        &iterations between reports is not consistent with the number of reports requested.')
                end if
                kp%nreports = tokens%size()
                if (.not. allocated(kp%niters)) allocate(kp%niters(kp%nreports))
                do i = 1, kp%nreports - 1
                    kp%niters(i) = to_int(tokens%next())
                end do
                kp%niters(kp%nreports) = 0
            case ("MEMORY-FACTOR")
                memory_factor_kp = to_realdp(tokens%next())
            case ("NUM-WALKERS-PER-SITE-INIT")
                nwalkers_per_site_init = to_realdp(tokens%next())
            case ("AVERAGEMCEXCITS-HAMIL")
                av_mc_excits_kp = to_realdp(tokens%next())
            case ("EXACT-HAMIL-SPAWNING")
                tExactHamilSpawning = .true.
            case ("EXACT-HAMIL-SPAWNING-FRAC")
                kp_hamil_exact_frac = to_realdp(tokens%next())
            case ("EXACT-HAMIL")
                tExactHamil = .true.
            case ("FULLY-STOCHASTIC-HAMIL")
                tFullyStochasticHamil = .true.
            case ("INIT-CORRECT-WALKER-POP")
                tInitCorrectNWalkers = .true.
            case ("OCCUPY-DETERM-INIT")
                tOccDetermInit = .true.
            case ("INIT-CONFIG-SEEDS")
                if (kp%nconfigs == 0) call stop_all(t_r, 'Please use the num-init-configs option to enter the number of &
                                                          &initial configurations before using the init-vec-seeds option.')
                tUseInitConfigSeeds = .true.
                allocate(init_config_seeds(kp%nconfigs))
                do i = 1, kp%nconfigs
                    init_config_seeds(i) = to_int(tokens%next())
                end do
            case ("ALL-SYM-SECTORS")
                tAllSymSectors = .true.
            case ("WRITE-AVERAGE-MATRICES")
                tOutputAverageKPMatrices = .true.
            case ("SCALE-POPULATION")
                tScalePopulation = .true.
            case ("EXCITED-INIT-STATE")
                tExcitedInitState = .true.
                ! Read in the number of excited states to be used.
                nexcit = to_int(tokens%next())
                allocate(kpfciqmc_ex_labels(nexcit))
                allocate(kpfciqmc_ex_weights(nexcit))

                ! Read in which excited states to use.
                if (file_reader%nextline(tokens, skip_empty=.false.)) then
                    do j = 1, nexcit
                        kpfciqmc_ex_labels(j) = to_int(tokens%next())
                    end do
                else
                    call stop_all(this_routine, 'Unexpected EOF reached.')
                end if

                ! Read in the relative weights for the trial excited states in
                ! the initial state.
                if (file_reader%nextline(tokens, skip_empty=.false.)) then
                    do j = 1, nexcit
                        kpfciqmc_ex_weights(j) = to_realdp(tokens%next())
                    end do
                else
                    call stop_all(this_routine, 'Unexpected EOF reached.')
                end if

                ! Normalise weights so that they sum to 1.
                kpfciqmc_ex_weights = kpfciqmc_ex_weights / sum(kpfciqmc_ex_weights)

            case ("OVERLAP-PERTURB-ANNIHILATE")
                alloc_popsfile_dets = .true.
                tOverlapPert = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                npert = to_int(tokens%next())
                if (.not. allocated(overlap_pert)) then
                    allocate(overlap_pert(npert))
                else
                    if (npert /= size(overlap_pert)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert
                    if (file_reader%nextline(tokens, skip_empty=.false.)) then
                        overlap_pert(i)%nannihilate = tokens%remaining_items()
                        allocate(overlap_pert(i)%ann_orbs(tokens%remaining_items()))
                        do j = 1, size(overlap_pert(i)%ann_orbs)
                            overlap_pert(i)%ann_orbs(j) = to_int(tokens%next())
                        end do
                        ! Create the rest of the annihilation-related
                        ! components of the pops_pert object.
                        call init_perturbation_annihilation(overlap_pert(i))
                    else
                        call stop_all(this_routine, 'Unexpected EOF reached.')
                    end if
                end do

            case ("OVERLAP-PERTURB-CREATION")
                alloc_popsfile_dets = .true.
                tOverlapPert = .true.
                tWritePopsNorm = .true.

                ! Read in the number of perturbation operators which are about
                ! to be read in.
                npert = to_int(tokens%next())
                if (.not. allocated(overlap_pert)) then
                    allocate(overlap_pert(npert))
                else
                    if (npert /= size(overlap_pert)) then
                        call stop_all(t_r, "A different number of creation and annihilation perturbation have been requested.")
                    end if
                end if

                do i = 1, npert
                    if (file_reader%nextline(tokens, .false.)) then
                        overlap_pert(i)%ncreate = tokens%remaining_items()
                        allocate(overlap_pert(i)%crtn_orbs(tokens%remaining_items()))
                        do j = 1, size(overlap_pert(i)%crtn_orbs)
                            overlap_pert(i)%crtn_orbs(j) = to_int(tokens%next())
                        end do
                        ! Create the rest of the creation-related
                        ! components of the pops_pert object.
                        call init_perturbation_creation(overlap_pert(i))
                    else
                        call stop_all(this_routine, 'Unexpected EOF reached.')
                    end if
                end do

            case ("DOUBLES-TRIAL")
                kp_trial_space_in%tDoubles = .true.
            case ("CAS-TRIAL")
                kp_trial_space_in%tCAS = .true.
                tSpn = .true.
                kp_trial_space_in%occ_cas = to_int(tokens%next())  ! Number of electrons in the CAS
                kp_trial_space_in%virt_cas = to_int(tokens%next()) ! Number of virtual spin-orbitals in the CAS
            case ("RAS-TRIAL")
                kp_trial_space_in%tRAS = .true.
                ras_size_1 = to_int(tokens%next())  ! Number of spatial orbitals in RAS1.
                ras_size_2 = to_int(tokens%next())  ! Number of spatial orbitals in RAS2.
                ras_size_3 = to_int(tokens%next())  ! Number of spatial orbitals in RAS3.
                ras_min_1 = to_int(tokens%next()) ! Min number of electrons (alpha and beta) in RAS1 orbs.
                ras_max_3 = to_int(tokens%next()) ! Max number of electrons (alpha and beta) in RAS3 orbs.
                kp_trial_space_in%ras%size_1 = int(ras_size_1, sp)
                kp_trial_space_in%ras%size_2 = int(ras_size_2, sp)
                kp_trial_space_in%ras%size_3 = int(ras_size_3, sp)
                kp_trial_space_in%ras%min_1 = int(ras_min_1, sp)
                kp_trial_space_in%ras%max_3 = int(ras_max_3, sp)
            case ("MP1-TRIAL")
                kp_trial_space_in%tMP1 = .true.
                kp_trial_space_in%mp1_ndets = to_int(tokens%next())
            case ("HF-TRIAL")
                kp_trial_space_in%tHF = .true.
            case ("POPS-TRIAL")
                kp_trial_space_in%tPops = .true.
                kp_trial_space_in%npops = to_int(tokens%next())
            case ("READ-TRIAL")
                kp_trial_space_in%tRead = .true.
            case ("FCI-TRIAL")
                kp_trial_space_in%tFCI = .false.
            case ("UNPAIRED-REPLICAS")
                tPairedReplicas = .false.
#if defined(PROG_NUMRUNS_)
                nreplicas = 1
#endif
            case ("ORTHOGONALISE-REPLICAS")
                tOrthogKPReplicas = .true.
                if (tokens%remaining_items() > 0) then
                    orthog_kp_iter = to_int(tokens%next())
                end if
            case default
                call stop_all(this_routine, "Keyword "//trim(w)//" not recognized in kp-fciqmc block", .true.)
            end select
        end do read_inp

        if ((.not. vary_niters)) then
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

        use CalcData, only: tSemiStochastic, tUseRealCoeffs, AvMCExcits, tPairedReplicas
        use fcimc_initialisation, only: SetupParameters, InitFCIMCCalcPar, init_fcimc_fn_pointers
        use FciMCData, only: tPopsAlreadyRead, nWalkerHashes, SpawnVecKP
        use FciMCData, only: SpawnVecKP2, MaxSpawned
        use FciMCData, only: SpawnedPartsKP, SpawnedPartsKP2, MaxWalkersUncorrected
        use FciMCData, only: iter_data_fciqmc, spawn_ht, nhashes_spawn, tReplicaReferencesDiffer
        use FciMCParMod, only: WriteFciMCStatsHeader, write_fcimcstats2, tSinglePartPhase
        use hash, only: init_hash_table
        use LoggingData, only: tFCIMCStats2
        use Parallel_neci, only: MPIBarrier
        use SystemData, only: G1, tHeisenberg, tAllSymSectors
        use util_mod, only: int_fmt

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: ierr, krylov_vecs_mem, krylov_ht_mem, matrix_mem, spawn_ht_mem, i
        character(len=*), parameter :: t_r = "init_kp_fciqmc"

        ! Checks.
        if (.not. tUseRealCoeffs) call stop_all(t_r, 'kp-fciqmc can only be run using &
            &real coefficients).')
        if (tExactHamil .and. nProcessors /= 1) call stop_all(t_r, 'The exact-hamil &
            &option can only be used when running with one processor.')
        if (theisenberg .and. tAllSymSectors) call stop_all(t_r, 'The option to use all &
            &symmetry sectors at once has not been implemented with the Heisenberg model.')
        if (n_int == 4) call stop_all(t_r, 'Use of RealCoefficients does not work with 32 bit &
             &integers due to the use of the transfer operation from dp reals to 64 bit integers.')
        if (tExcitedStateKP) then
            if (tPairedReplicas .and. inum_runs /= 2 * kp%nvecs) then
                call stop_all(t_r, 'When using the ExcitedStateKP option with paired replicas, the &
                                   &number of replicas must be twice the number of states to be calculated.')
            else if ((.not. tPairedReplicas) .and. inum_runs /= kp%nvecs) then
                call stop_all(t_r, 'When using the ExcitedStateKP option without paired replicas, the &
                                   &number of replicas must be equal to the number of states to be calculated.')
            end if
        end if

        ! Call external NECI initialisation routines.
        tPopsAlreadyRead = .false.
        call SetupParameters()
        call InitFCIMCCalcPar()
        call init_fcimc_fn_pointers()

        write(stdout, '(/,12("="),1x,a9,1x,12("="))') "KP-FCIQMC"

        if (tSemiStochastic .and. (.not. tFullyStochasticHamil)) then
            tSemiStochasticKPHamil = .true.
        else
            tSemiStochasticKPHamil = .false.
        end if

#if defined(PROG_NUMRUNS_)
        ! The number of *repeated* replicas for a particular state.
        ! i.e. if tExcitedState = .true. then this is the number of repeated
        ! replicas for each excited state. If tExcitedState = .false. (doing
        ! the old KP-FCIQMC algorithm), this is the number of repeats for the
        ! KP-FCIQMC wave function.

#ifndef CMPLX_
        if (tPairedReplicas) then
            lenof_sign_kp = 2
        else
            lenof_sign_kp = 1
        end if
#else
        if (tPairedReplicas) then
            lenof_sign_kp = 4
        else
            lenof_sign_kp = 2
        end if
#endif
#endif

        ! The number of elements required to store all replicas of all Krylov vectors.
        lenof_all_signs = lenof_sign_kp * kp%nvecs
        ! The total length of a bitstring containing all Krylov vectors.
        ! +1 one is for the flag integer
        NIfTotKP = nifd + lenof_all_signs + 1

        ! Set up kp_ind_* arrays.
        allocate(kp_ind_1(kp%nvecs))
        allocate(kp_ind_2(kp%nvecs))

        if (tPairedReplicas) then
            do i = 1, kp%nvecs
                kp_ind_1(i) = 2 * i - 1
                kp_ind_2(i) = 2 * i
            end do
        else
            do i = 1, kp%nvecs
                kp_ind_1(i) = i
                kp_ind_2(i) = i
            end do
        end if

        if (.not. tExcitedStateKP) then
            ! Allocate all of the KP arrays.
            nhashes_kp = nWalkerHashes
            TotWalkersKP = 0
            krylov_vecs_length = nint(MaxWalkersUncorrected * memory_factor_kp * kp%nvecs)
            nkrylov_amp_elems_tot = lenof_sign * kp%nvecs * krylov_vecs_length

            ! Allocate the krylov_vecs array.
            ! The number of MB of memory required to allocate krylov_vecs.
            krylov_vecs_mem = krylov_vecs_length * (NIfTotKP + 1) * size_n_int / 1000000
            write(stdout, '(a73,'//int_fmt(krylov_vecs_mem, 1)//')') "About to allocate array to hold all Krylov vectors. &
                                           &Memory required (MB):", krylov_vecs_mem
            write(stdout, '(a13)', advance='no') "Allocating..."; call neci_flush(6)
            allocate(krylov_vecs(0:NIfTotKP, krylov_vecs_length), stat=ierr)
            if (ierr /= 0) then
                write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating krylov_vecs array.")
            else
                write(stdout, '(1x,a5)') "Done."
            end if
            call neci_flush(6)
            krylov_vecs = 0_n_int

            ! Allocate the krylov_helems array.
            ! The number of MB of memory required to allocate krylov_helems.
            krylov_vecs_mem = krylov_vecs_length * size_n_int / 1000000
            write(stdout, '(a103,'//int_fmt(krylov_vecs_mem, 1)//')') "About to allocate array to hold diagonal Hamiltonian &
                                           &elements for Krylov vectors. Memory required (MB):", krylov_vecs_mem
            write(stdout, '(a13)', advance='no') "Allocating..."; call neci_flush(6)
            allocate(krylov_helems(krylov_vecs_length), stat=ierr)
            if (ierr /= 0) then
                write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating krylov_helems array.")
            else
                write(stdout, '(1x,a5)') "Done."
            end if
            call neci_flush(6)
            krylov_helems = 0.0_dp

            ! Allocate the hash table to krylov_vecs.
            ! The number of MB of memory required to allocate krylov_vecs_ht.
            ! Each node requires 16 bytes.
            krylov_ht_mem = nhashes_kp * 16 / 1000000
            write(stdout, '(a78,'//int_fmt(krylov_ht_mem, 1)//')') "About to allocate hash table to the Krylov vector array. &
                                           &Memory required (MB):", krylov_ht_mem
            write(stdout, '(a13)', advance='no') "Allocating..."; call neci_flush(6)
            allocate(krylov_vecs_ht(nhashes_kp), stat=ierr)
            if (ierr /= 0) then
                write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
                call stop_all(t_r, "Error allocating krylov_vecs_ht array.")
            else
                write(stdout, '(1x,a5)') "Done."
                write(stdout, '(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                                  &increase as further nodes are added."
            end if

            call init_hash_table(krylov_vecs_ht)

            allocate(SpawnVecKP(0:IlutBits%ind_pop + lenof_all_signs - 1, MaxSpawned), stat=ierr)
            allocate(SpawnVecKP2(0:IlutBits%ind_pop + lenof_all_signs - 1, MaxSpawned), stat=ierr)
            SpawnVecKP(:, :) = 0_n_int
            SpawnVecKP2(:, :) = 0_n_int
            SpawnedPartsKP => SpawnVecKP
            SpawnedPartsKP2 => SpawnVecKP2

            if (tSemiStochastic) then
                associate(rep => cs_replicas(core_run))
                    allocate(partial_determ_vecs_kp(lenof_all_signs, rep%determ_sizes(iProcIndex)), stat=ierr)
                    allocate(full_determ_vecs_kp(lenof_all_signs, rep%determ_space_size), stat=ierr)
                end associate
                partial_determ_vecs_kp = 0.0_dp
                full_determ_vecs_kp = 0.0_dp
            end if

        end if

        ! Allocate the hash table to the spawning array.
        ! The number of MB of memory required to allocate spawn_ht.
        ! Each node requires 16 bytes.
        nhashes_spawn = int(0.8 * MaxSpawned)
        spawn_ht_mem = nhashes_spawn * 16 / 1000000
        write(stdout, '(a78,'//int_fmt(spawn_ht_mem, 1)//')') "About to allocate hash table to the spawning array. &
                                       &Memory required (MB):", spawn_ht_mem
        write(stdout, '(a13)', advance='no') "Allocating..."; call neci_flush(6)
        allocate(spawn_ht(nhashes_spawn), stat=ierr)
        if (ierr /= 0) then
            write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating spawn_ht array.")
        else
            write(stdout, '(1x,a5)') "Done."
            write(stdout, '(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                              &increase as further nodes are added."
        end if

        call init_hash_table(spawn_ht)

        ! (2*kp%nrepeats+16) arrays with (kp%nvecs**2) 8-byte elements each.
        matrix_mem = (2 * kp%nrepeats + 16) * (kp%nvecs**2) * 8 / 1000000
        write(stdout, '(a66,'//int_fmt(matrix_mem, 1)//')') "About to allocate various subspace matrices. &
                                       &Memory required (MB):", matrix_mem
        write(stdout, '(a13)', advance='no') "Allocating..."
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

        write(stdout, '(1x,a5)') "Done."

        ! If av_mc_excits_kp hasn't been set by the user, just use AvMCExcits.
        if (near_zero(av_mc_excits_kp)) av_mc_excits_kp = AvMCExcits

        ! Initialize
        kp_overlap_mean = 0.0_dp
        kp_hamil_mean = 0.0_dp

        MaxSpawnedEachProc = int(0.88_dp * real(MaxSpawned, dp) / nProcessors)

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

    subroutine init_kp_fciqmc_repeat(iconfig, irepeat, nrepeats, nvecs, iter_data)

        use CalcData, only: tStartSinglePart, InitialPart, InitWalkers, DiagSft, iPopsFileNoRead
        use CalcData, only: tPairedReplicas
        use FciMCData, only: iter, InputDiagSft, PreviousCycles, OldAllAvWalkersCyc, proje_iter
        use FciMCData, only: proje_iter_tot, AllGrowRate, SpawnedParts, fcimc_iter_data
        use FciMCParMod, only: tSinglePartPhase
        use fcimc_output, only: WriteFCIMCStats, write_fcimcstats2
        use hash, only: clear_hash_table
        use initial_trial_states
        use LoggingData, only: tFCIMCStats2, tPrintDataTables
        use util_mod, only: int_fmt

        integer, intent(in) :: iconfig, irepeat, nrepeats, nvecs

        integer :: ndets_this_proc, nexcit
        real(dp), allocatable :: evals(:)
        HElement_t(dp), allocatable :: evecs_this_proc(:, :), init_vecs(:, :)
        integer(MPIArg) :: space_sizes(0:nProcessors - 1), space_displs(0:nProcessors - 1)
        type(fcimc_iter_data), intent(in) :: iter_data

        write(stdout, '(1x,a22,'//int_fmt(irepeat, 1)//')') "Starting repeat number", irepeat

        if (tExcitedStateKP) then
            nexcit = nvecs
            allocate(evals(nexcit))

            ! Create the trial excited states.
            call calc_trial_states_lanczos(kp_trial_space_in, nexcit, ndets_this_proc, SpawnedParts, &
                                           evecs_this_proc, evals, space_sizes, space_displs)

            ! Set the populations of these states to the requested value.
            call set_trial_populations(nexcit, ndets_this_proc, evecs_this_proc)
            ! Set the trial excited states as the FCIQMC wave functions.
            call set_trial_states(ndets_this_proc, evecs_this_proc, SpawnedParts, .true., tPairedReplicas)

            deallocate(evecs_this_proc, evals)

        else if (tExcitedInitState) then
            nexcit = maxval(kpfciqmc_ex_labels)
            allocate(evals(nexcit))

            ! Create the trial excited states.
            call calc_trial_states_lanczos(kp_trial_space_in, nexcit, ndets_this_proc, SpawnedParts, &
                                           evecs_this_proc, evals, space_sizes, space_displs)

            ! Extract the desried initial excited states and average them.
            call create_init_excited_state(ndets_this_proc, evecs_this_proc, kpfciqmc_ex_labels, kpfciqmc_ex_weights, init_vecs)
            ! Set the populations of these states to the requested value.
            call set_trial_populations(1, ndets_this_proc, init_vecs)
            ! Set the trial excited states as the FCIQMC wave functions.
            call set_trial_states(ndets_this_proc, init_vecs, SpawnedParts, .true., tPairedReplicas)

            deallocate(evecs_this_proc, init_vecs, evals)

        else
            ! If starting from multiple POPSFILEs then set this counter so that the
            ! correct POPSFILE is read in this time. To read in POPSFILE.x,
            ! iPopsFileNoRead needs to be set to -x-1. We want to read in POPSFILE
            ! numbers 0 to kp%nconfigs-1
            if (tMultiplePopStart) iPopsFileNoRead = -(iconfig - 1) - 1

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
            OldAllAvWalkersCyc = InitWalkers * nProcessors
        end if
        proje_iter = 0.0_dp
        proje_iter_tot = 0.0_dp
        AllGrowRate = 0.0_dp

        ! Setting this variable to true stops the shift from varying instantly.
        tSinglePartPhase = tSinglePartPhaseKPInit

        ! Print out initial stats.
        if (tPrintDataTables) then
            if (tFCIMCStats2) then
                call write_fcimcstats2(iter_data)
            else
                call WriteFCIMCStats()
            end if
        end if

    end subroutine init_kp_fciqmc_repeat

    subroutine init_kp_fciqmc_iter(iter_data, determ_index)

        use FciMCData, only: FreeSlot, iStartFreeSlot, iEndFreeSlot, fcimc_iter_data, InitialSpawnedSlots
        use FciMCData, only: ValidSpawnedList, spawn_ht
        use FciMCParMod, only: rezero_iter_stats_each_iter
        use hash, only: clear_hash_table
        use rdm_data, only: rdm_definitions

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

        call rezero_iter_stats_each_iter(iter_data, rdm_definitions)

    end subroutine init_kp_fciqmc_iter

    subroutine create_initial_config(iconfig, irepeat, nrepeats)

        use CalcData, only: tStartSinglePart, InitialPart, InitWalkers, tSemiStochastic, tReadPops
        use dSFMT_interface, only: dSFMT_init
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
        integer(int64) :: i
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign(lenof_sign_kp), TotPartsCheck(lenof_sign_kp)
        real(dp) :: nwalkers_target
        real(dp) :: norm, all_norm
        real(dp) :: total_time_before, total_time_after
        logical :: tCoreDet
        character(len=*), parameter :: t_r = "create_init_config"

        ! Clear everything from any previous repeats or starting configurations.
        call clear_hash_table(HashIndex)

        if (tStartSinglePart) then
            nwalkers_target = real(InitialPart, dp)
        else
            nwalkers_target = InitWalkers * nProcessors
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
                        if (any(abs(TotPartsCheck - TotParts) > 1.0e-12_dp)) then
                            call stop_all(t_r, "Inconsistent values of TotParts calculated.")
                        end if
                        TotParts = TotParts * scaling_factor
                        TotPartsOld = TotParts
                        AllTotParts = AllTotParts * scaling_factor
                        AllTotPartsOld = AllTotParts
                    end if
                else
                    ! Put a walker on the Hartree-Fock, with the requested amplitude.
                    call InitFCIMC_HF()
                end if

                if (tSemiStochastic) then
                    ! core_space stores all core determinants from all processors. Move those on this
                    ! processor to SpawnedParts, which add_core_states_currentdet_hash uses.
                    call copy_core_dets_to_spawnedparts(cs_replicas(core_run))
                    ! Any core space determinants which are not already in CurrentDets will be added
                    ! by this routine.
                    call add_core_states_currentdet_hash(core_run)
                    if (tStartCoreGroundState .and. (.not. tReadPops)) &
                        call start_walkers_from_core_ground(tPrintInfo=.false., run=core_run)
                end if

            else if (tFiniteTemp) then
                ! Convert the initial number of walkers to an integer. Note that on multiple
                ! processors this may round up the requested number of walkers slightly.
                nwalkers_int = ceiling(nwalkers_target / real(nProcessors, dp))

                ! If requested, reset the random number generator with the requested seed
                ! before creating the random initial configuration.
                if (tUseInitConfigSeeds) call dSFMT_init((iProcIndex + 1) * init_config_seeds(iconfig))

                write(stdout, '(a44)', advance='no') "# Generating initial walker configuration..."
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
                write(stdout, '(1x,a31,f9.3)') "Complete. Time taken (seconds):", total_time_after - total_time_before

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
            do i = 1, int(TotWalkers)
                ! Copy across the bitstring encoding of the determinant and also the walker signs.
                CurrentDets(0:IlutBits%ind_pop + lenof_sign_kp - 1, i) = krylov_vecs(0:IlutBits%ind_pop + lenof_sign_kp - 1, i)
                ! Copy across the flags.
                CurrentDets(NIfTot, i) = krylov_vecs(NIfTotKP, i)
            end do
            call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, int(TotWalkers), .true.)
        end if

        ! Calculate and store the diagonal element of the Hamiltonian for determinants in CurrentDets.
        call fill_in_diag_helements()

        ! If starting from this configuration more than once, store the relevant data for next time.
        if (nrepeats > 1 .and. irepeat == 1) then
            HolesInList = 0
            do i = 1, TotWalkers
                int_sign = CurrentDets(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign_kp - 1, i)
                call extract_sign(CurrentDets(:, i), real_sign)
                tCoreDet = check_determ_flag(CurrentDets(:, i))
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
        mem_reqd = int(TotWalkers / 1000000_int64) * (NIfTotKP + 1) * size_n_int

        write(stdout, '(a73,'//int_fmt(mem_reqd, 1)//')') "About to allocate array to hold the perturbed &
                                           &ground state. Memory required (MB):", mem_reqd
        write(stdout, '(a13)', advance='no') "Allocating..."
        call neci_flush(6)
        allocate(perturbed_ground(0:NIfTot, TotWalkers), stat=ierr)
        if (ierr /= 0) then
            write(stdout, '(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating array.")
        else
            write(stdout, '(1x,a5)') "Done."
        end if
        call neci_flush(6)

        perturbed_ground = CurrentDets(0:NIfTot, 1:TotWalkers)

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
        integer :: nspawns, ndets, err
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_amp, walker_sign(lenof_sign_kp)
        logical :: tInitiatorTemp
        type(fcimc_iter_data) :: unused_data
        integer(n_int), pointer :: PointTemp(:, :)
        HElement_t(dp) :: hdiag_spawn

        call allocate_iter_data(unused_data)

        ! Turn off the initiator method for the annihilation steps to be used here.
        tInitiatorTemp = tTruncInitiator
        tTruncInitiator = .false.

        ! Set the spawning slots to their starting positions.
        ValidSpawnedList = InitialSpawnedSlots

        ilut = 0_n_int
        nspawns = ceiling(real(nwalkers, dp) / nwalkers_per_site_init)

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
            if (r < 0.5_dp) walker_amp = -1.0_dp * walker_amp

            do ireplica = 1, inum_runs
                walker_sign = 0.0_dp
                walker_sign(ireplica) = walker_amp
                call create_particle(nI, ilut, walker_sign, ireplica, hdiag_spawn, err)
            end do
        end do

        ! Perform annihilation steps:
        ! Send the walkers to their correct processors. The resulting walkers will be stored in
        ! SpawnedParts2.
        call SendProcNewParts(ndets, tSingleProc=.false.)
        ! CompressSpawnedList works on SpawnedParts, not SpawnedParts2, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp
        call CompressSpawnedList(ndets, unused_data)

        ! Finally, add the determinants in the spawned walker list to the main walker list.
        ! Copy the determinants themselves to CurrentDets.
        TotParts = 0.0_dp
        do i = 1, ndets
            CurrentDets(:, i) = SpawnedParts(0:NIfTot, i)
            walker_sign = transfer(CurrentDets(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign_kp - 1, i), walker_sign)
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
            call copy_core_dets_to_spawnedparts(cs_replicas(core_run))
            call add_core_states_currentdet_hash(core_run)
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
        use FciMCData, only: HashIndex, TotWalkers, CurrentDets, HFDet
        use FciMCData, only: TotParts, TotPartsOld, AllTotParts, AllTotPartsOld, ilutHF
        use core_space_util, only: cs_replicas
        use load_balance_calcnodes, only: DetermineDetNode
        use hash, only: rm_unocc_dets_from_hash_table, hash_table_lookup
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
        associate(rep => cs_replicas(1))
        do
            ! If using the tOccupyDetermSpace option then we want to put walkers
            ! on states in the deterministic space first.
            ! If not, or if we have finished doing this, generate determinants
            ! randomly and uniformly.
            if (tOccupyDetermSpace .and. (.not. tDetermAllOccupied)) then
                ideterm = ideterm + 1
                ilut = rep%core_space(:, ideterm + rep%determ_displs(iProcIndex))
                call decode_bit_det(nI, ilut)
                ! If we have now occupied all deterministic states.
                if (ideterm == rep%determ_sizes(iProcIndex)) tDetermAllOccupied = .true.
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
            if (r < 0.5_dp) real_sign_1 = -1.0_dp * real_sign_1
            int_sign = transfer(real_sign_1, int_sign)

            tDetFound = .false.
            ! Search the hash table to see if this determinant is in CurrentDets
            ! already.
            call hash_table_lookup(nI, ilut, nifd, HashIndex, CurrentDets, det_ind, hash_val, tDetFound)
            if (tDetFound) then
                ! This determinant is already in CurrentDets. Just need to add
                ! the sign of this new walker on and update stats.
                int_sign = CurrentDets(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign_kp - 1, det_ind)
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
                CurrentDets(0:nifd, det_ind) = ilut(0:nifd)
                CurrentDets(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign_kp - 1, det_ind) = int_sign
                CurrentDets(IlutBits%ind_flag, det_ind) = 0_n_int
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
            call copy_core_dets_to_spawnedparts(rep)
            call add_core_states_currentdet_hash(core_run)
        else
            ! Some determinants may have become occupied and then unoccupied in
            ! the course of the above. We need to remove the entries for these
            ! determinants from the hash table. For semi-stochastic calculations
            ! this is done in add_core_states_currentdet_hash.
            call rm_unocc_dets_from_hash_table(HashIndex, CurrentDets, ndets)
        end if
        end associate

    end subroutine generate_init_config_this_proc

    subroutine create_init_excited_state(ndets_this_proc, trial_vecs, ex_state_labels, ex_state_weights, init_vec)

        integer, intent(in) :: ndets_this_proc
        HElement_t(dp), intent(in) :: trial_vecs(:, :)
        integer, intent(in) :: ex_state_labels(:)
        real(dp), intent(in) :: ex_state_weights(:)
        HElement_t(dp), allocatable, intent(out) :: init_vec(:, :)

        real(dp) :: real_sign(lenof_sign)
        integer :: i, j, ierr
        character(len=*), parameter :: t_r = "create_init_excited_state"

        allocate(init_vec(1, ndets_this_proc), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error in MPIScatterV call.")

        do i = 1, ndets_this_proc
            init_vec(1, i) = h_cast(0.0_dp)
            do j = 1, size(ex_state_labels)
                init_vec(1, i) = init_vec(1, i) + ex_state_weights(j) * trial_vecs(ex_state_labels(j), i)
            end do
        end do

    end subroutine create_init_excited_state

    subroutine scale_population(walker_list, ndets, target_pop, input_pop, scaling_factor)

        ! Take an input list of walkers, find the total walker population in
        ! the list, and then multiply all the walker signs by some factor in
        ! order for the walker list to have the requested target population.

        integer(n_int), intent(inout) :: walker_list(:, :)
        integer(int64), intent(in) :: ndets
        real(dp), intent(in) :: target_pop
        real(dp), intent(out) :: input_pop(lenof_sign_kp)
        real(dp), intent(out) :: scaling_factor

        integer(int64) :: i
        real(dp) :: real_sign(lenof_sign_kp), all_input_pop(lenof_sign_kp)

        input_pop = 0.0_dp

        ! First, find the population of the walkers in walker_list.
        do i = 1, ndets
            call extract_sign(walker_list(:, i), real_sign)
            input_pop = input_pop + abs(real_sign)
        end do

        call MPISumAll(input_pop, all_input_pop)

        ! Just use the first particle type to determine the scaing factor.
        scaling_factor = target_pop / all_input_pop(1)

        ! Now multiply the walker signs by scaling_factor.
        do i = 1, ndets
            call extract_sign(walker_list(:, i), real_sign)
            real_sign = real_sign * scaling_factor
            call encode_sign(walker_list(:, i), real_sign)
        end do

    end subroutine scale_population

end module kp_fciqmc_init
