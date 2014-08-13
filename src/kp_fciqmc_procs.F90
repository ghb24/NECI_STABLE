#include "macros.h"
 
module kp_fciqmc_procs
 
    use AnnihilationMod, only: SendProcNewParts, CompressSpawnedList, RemoveDetHashIndex
    use bit_rep_data
    use bit_reps, only: decode_bit_det, encode_sign, flag_is_initiator
    use CalcData, only: tTruncInitiator, tStartSinglePart, InitialPart, InitWalkers
    use CalcData, only: tSemiStochastic, tReadPops, tUseRealCoeffs, tau, DiagSft
    use CalcData, only: AvMCExcits, tWritePopsNorm, iPopsFileNoRead, pops_norm
    use constants
    use DetBitOps, only: DetBitEq, EncodeBitDet, IsAllowedHPHF, FindBitExcitLevel
    use Determinants, only: get_helement
    use dSFMT_interface , only: dSFMT_init, genrand_real2_dSFMT
    use FciMCData
    use FciMCParMod, only: create_particle, InitFCIMC_HF, SetupParameters, InitFCIMCCalcPar
    use FciMCParMod, only: init_fcimc_fn_pointers, WriteFciMCStats, WriteFciMCStatsHeader
    use FciMCParMod, only: rezero_iter_stats_each_iter, tSinglePartPhase
    use gndts_mod, only: gndts
    use hash, only: FindWalkerHash, init_hash_table, reset_hash_table, fill_in_hash_table
    use hash, only: DetermineDetNode, remove_node
    use hilbert_space_size, only: CreateRandomExcitLevDetUnbias, create_rand_heisenberg_det
    use hilbert_space_size, only: create_rand_det_no_sym
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use LoggingData, only: tIncrementPops
    use Parallel_neci, only: MPIBarrier, iProcIndex, MPISum, MPIReduce, nProcessors, MPIAllReduce
    use ParallelHelper, only: root
    use perturbations, only: init_perturbation_creation, init_perturbation_annihilation
    use PopsfileMod, only: read_popsfile_wrapper
    use procedure_pointers
    use semi_stoch_procs, only: copy_core_dets_this_proc_to_spawnedparts, fill_in_CurrentH
    use semi_stoch_procs, only: add_core_states_currentdet_hash, start_walkers_from_core_ground
    use semi_stoch_procs, only: check_determ_flag
    use sym_mod, only: getsym
    use SystemData, only: nel, nbasis, BRR, nBasisMax, G1, tSpn, lms, tParity, SymRestrict
    use SystemData, only: BasisFn, tHeisenberg, tHPHF, tAllSymSectors
    use util_mod, only: get_free_unit

    implicit none

    type kp_fciqmc_data
        ! The number of different initial walker configurations to start
        ! calculations from.
        integer :: nconfigs
        ! The current configuration.
        integer, pointer :: iconfig
        ! The number of simulations to perform for each initial walker
        ! configuration.
        integer :: nrepeats
        ! The current repeat.
        integer, pointer :: irepeat
        ! The number of different Krylov vectors to sample (the number of
        ! vectors which form the Krylov subspace at the end of a calculation).
        integer :: nvecs
        ! The current Krylov vector.
        integer, pointer :: ivec
        ! The number of iterations to perform *between each Krylov vector being
        ! sampled*. niters(i) holds the number to be performed between the
        ! i-th and (i+1)th vectors. The final element is not set by the user,
        ! it is always set to 1 (because we don't want to go any further
        ! beyond the final Krylov vector).
        integer, allocatable :: niters(:)

        ! Stores of the overlap and projected Hamiltonian matrices.
        ! If tStoreKPMatrices is .true., then all matrices will be held for
        ! the current initial configuration, including all repeats. The third
        ! index will therefore run from 1 to nrepeats.
        ! If tKPMatrices is .false., then only the matrices from the current
        ! repeat will be held, and the third index will always be 1.
        real(dp), pointer :: overlap_matrices(:,:,:)
        real(dp), pointer :: hamil_matrices(:,:,:)
        ! Pointers to the matrices for the current repeat only.
        real(dp), pointer :: overlap_matrix(:,:)
        real(dp), pointer :: hamil_matrix(:,:)
    end type

    type(kp_fciqmc_data) :: kp

    ! Information for the krylov_vecs arrays, which holds all of the Krylov
    ! vectors together simultaneously.
    ! The number of hash values for the hash table used to access krylov_vecs.
    integer :: nhashes_kp
    ! The total number of slots in the krylov_vecs arrays (the number of
    ! different determinants which it can hold).
    integer :: krylov_vecs_length
    character(2) :: kp_length_fmt
    ! krylov_vecs_length is calculated as the total number of walkers requested
    ! (by the TotalWalkers option in the Calc block) multipled by this factor,
    ! and then multipled by the number of Krylov vectors (for now...).
    real(dp) :: memory_factor_kp
    ! The current number of different determinants held in krylov_vecs.
    integer :: TotWalkersKP
    integer(n_int), allocatable :: krylov_vecs(:,:)
    type(ll_node), pointer :: krylov_vecs_ht(:) 

    ! The number of elements in krylov_vecs which are used to store amplitudes.
    integer(int64) :: nkrylov_amp_elems_tot
    ! The nunber of amplitude elements in krylov_vecs which are non-zero.
    integer(int64) :: nkrylov_amp_elems_used

    ! The values of these variables for the current initial walker configuration.
    ! These are held in these variables so that, when we restart from this
    ! configuration, we can instantly reset the corresponding values.
    integer(int64) :: TotWalkersInit, AllTotWalkersInit
    real(dp) :: TotPartsInit(lenof_sign)
    real(dp) :: AllTotPartsInit(lenof_sign)

    ! If true then don't use a fixed number of iterations between each Krylov
    ! vector is taken, but vary the number.
    logical :: vary_niters

    ! The total sign length for all Krylov vectors together.
    integer :: lenof_sign_kp

    ! If true then calculate the projected Hamiltonian exactly (useful for
    ! testing only, in practice).
    logical :: tExactHamil
    ! If true, use the spawning from the main FCIQMC iterations to
    ! calculate the projected Hamiltonian.
    logical :: tHamilOnFly
    ! If true then use the semi-stochastic approach in the calculation
    ! of the projected Hamiltonian. This is on by default, if
    ! tSemiStochastic is true. If tFullyStochasticHamil if true then
    ! this logical will be false.
    logical :: tSemiStochasticKPHamil
    ! If true, don't use the semi-stochastic approach when calculating
    ! the projected Hamiltonian. False by default.
    logical :: tFullyStochasticHamil
    ! If true then a finite-temperature calculation is performed.
    logical :: tFiniteTemp
    ! If true then perform multiple kp-fciqmc calculations starting from
    ! several POPSFILEs in the running directory. POPSFILEs labelled
    ! between 0 and nconfigs-1 will be used.
    logical :: tMultiplePopStart
    ! For finite-temperature calculations, when creating the inital vector,
    ! this variable specifies how many walkers should be added to each
    ! chosen site.
    real(dp) :: nwalkers_per_site_init
    ! When estimating the projected Hamiltonian, this variable determines
    ! how many spawns (on average) each of the walkers in each of the
    ! Kyrlov vectors contribute.
    real(dp) :: av_mc_excits_kp
    ! If true then use generate_init_config_this_proc to generate the initial
    ! walker distribution for finite-temperature calculations. This will always
    ! generate the requested number of walkers (except for rounding when splitting
    ! this number between processors).
    logical :: tInitCorrectNWalkers
    ! If true then always occupy the deterministic space first when generating
    ! an initial configuration for finite-temperature calculations. Only after
    ! this is filled will other states be randomly occupied.
    logical :: tOccDetermInit
    ! In the projected Hamiltonian calculation, this holds the fraction of the
    ! spawning array to be filled before a communication round is performed.
    integer :: MaxSpawnedEachProc
    ! If true then the random number generator will be reset before genearting
    ! the next initial configuration (in a finite-temperature job) by using
    ! the values in init_config_seeds.
    logical :: tUseInitConfigSeeds
    ! See comments above.
    integer, allocatable :: init_config_seeds(:)
    ! If true, all repeats of the projected Hamiltonian and overlap matrices
    ! will be stored in memory for a given starting configuration until
    ! we move onto the next starting configuration. This means that we don't
    ! have to read all the matrices in again before averaging at the end
    ! of a tarting configuration, and so makes things a bit quicker. However,
    ! if this will use too much memory then this option can be turned off
    ! and the matrices will be read in before averaging instead.
    logical :: tStoreKPMatrices
    ! If true then the averaged projected Hamiltonian and overlap matrices
    ! will be output after completing all repeats of a given initial configuration.
    ! If more than one repeat has been performed, then the standard error
    ! will also be output.
    logical :: tOutputAverageKPMatrices

    ! If true then, after the perturbation operator has been applied (which
    ! will reduce the population), scale up the initial wave function to the
    ! requested population.
    logical :: tScalePopulation
    real(dp) :: scaling_factor

    ! If true then calculate the overlap of the final Hamiltonian eigenstates
    ! with a vector which is calculated by applying a perturbation operator
    ! (overlap_pert) to the popsfile wave function.
    logical :: tOverlapPert
    type(perturbation), save :: overlap_pert
    ! The result of overlap_pert applied to the wave function read in from
    ! a popsfile.
    integer(n_int), allocatable :: perturbed_ground(:,:)
    ! The overlaps of perturbed_ground with the various Krylov vectors.
    real(dp), allocatable :: pert_overlaps(:)
    ! The sum of pert_overlaps across all processes, calculated after all
    ! iterations in a repeat have been performed and only stored on the
    ! root process.
    real(dp), allocatable :: kp_all_pert_overlaps(:)

    ! Matrices used in the communication of the projected Hamiltonian and
    ! overlap matrices.
    real(dp), allocatable :: kp_mpi_matrices_in(:,:)
    real(dp), allocatable :: kp_mpi_matrices_out(:,:)

    ! After all repeats for a given initial configuration are complete, these
    ! arrays will hold the means and standard errors of the projected
    ! Hamiltonian and overlap matrices (unless only one sample was obtained,
    ! in which case the standard errors will be zero).
    real(dp), allocatable :: kp_hamil_mean(:,:)
    real(dp), allocatable :: kp_overlap_mean(:,:)
    real(dp), allocatable :: kp_hamil_se(:,:)
    real(dp), allocatable :: kp_overlap_se(:,:)

    ! These are used in the calculation of the final eigenvalues.
    ! Allocate these arrays once during initialisation and reuse them,
    ! as reallocating arrays over and over again isn't great to do...
    ! Arrays holding the eigenvalues of the Hamiltonian and overlap matrices.
    real(dp), allocatable :: kp_overlap_eigv(:)
    ! The overlap of the final vectors with the initial configuration.
    ! (the first Krylov vector).
    real(dp), allocatable :: kp_init_overlaps(:)
    ! The eigenvectors of the overlap matrix.
    real(dp), allocatable :: kp_overlap_eigenvecs(:,:)
    ! The matrix used to transform the Krylov vectors to a orthonormal basis.
    real(dp), allocatable :: kp_transform_matrix(:,:)
    ! Matrix used as temporary space during the transformation of the
    ! projected Hamiltonian into an orthonormal basis, in the Lowdin appraoch.
    real(dp), allocatable :: kp_inter_hamil(:,:)
    ! The final eigenvectors, but in the basis of Krylov vectors.
    ! This is used in the Lowdin approach: diagonalising kp_final_hamil
    ! will give the final eigenvectors in the orthonormal basis, then
    ! we do the inverse of the transformation procedure to get back to
    ! the Krylov basis.
    real(dp), allocatable :: kp_eigenvecs_krylov(:,:)

contains

    subroutine kp_fciqmc_read_inp()

        use input_neci

        logical :: eof
        character(len=100) :: w
        integer :: i, niters_temp, nvecs_temp
        character (len=*), parameter :: t_r = "kp_fciqmc_read_inp"

        ! Default values.
        kp%nconfigs = 1
        kp%nrepeats = 1
        kp%nvecs = 0
        niters_temp = 1
        memory_factor_kp = 1.0_dp
        nwalkers_per_site_init = 1.0_dp
        av_mc_excits_kp = 0.0_dp
        tFiniteTemp = .false.
        tMultiplePopStart = .false.
        tExactHamil = .false.
        tHamilOnFly = .false.
        tFullyStochasticHamil = .false.
        tInitCorrectNWalkers = .false.
        tOccDetermInit = .false.
        vary_niters = .false.
        tUseInitConfigSeeds = .false.
        tAllSymSectors = .false.
        tStoreKPMatrices = .true.
        tOutputAverageKPMatrices = .false.
        tOverlapPert = .false.
        tScalePopulation = .false.
        scaling_factor = 1.0_dp

        read_inp: do
            call read_line(eof)
            if (eof) then
                exit
            end if
            call readu(w)
            select case(w)
            case("END-KP-FCIQMC")
                exit read_inp
            case("FINITE-TEMPERATURE")
                tFiniteTemp = .true.
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
                end if
                kp%niters(kp%nvecs) = 0
            case("NUM-ITERS-BETWEEN-VECS")
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
            case("MEMORY-FACTOR")
                call getf(memory_factor_kp)
            case("NUM-WALKERS-PER-SITE-INIT")
                call getf(nwalkers_per_site_init)
            case("AVERAGEMCEXCITS-HAMIL")
                call getf(av_mc_excits_kp)
            case("EXACT-HAMIL")
                tExactHamil = .true.
            case("HAMIL-ON-FLY")
                tHamilOnFly = .true.
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
            case("WRITE-MATRICES-SEPARATE")
                tStoreKPMatrices = .false.
            case("WRITE-AVERAGE-MATRICES")
                tOutputAverageKPMatrices = .true.
            case("SCALE-POPULATION")
                tScalePopulation = .true.

            case("OVERLAP-PERTURB-ANNIHILATE")
                tOverlapPert = .true.
                tWritePopsNorm = .true.
                overlap_pert%nannihilate = nitems-1
                allocate(overlap_pert%ann_orbs(nitems-1))
                do i = 1, nitems-1
                    call readi(overlap_pert%ann_orbs(i))
                end do
                ! Create the rest of the annihilation-related
                ! components of the overlap_pert object.
                call init_perturbation_annihilation(overlap_pert)
            case("OVERLAP-PERTURB-CREATION")
                tOverlapPert = .true.
                tWritePopsNorm = .true.
                overlap_pert%ncreate = nitems-1
                allocate(overlap_pert%crtn_orbs(nitems-1))
                do i = 1, nitems-1
                    call readi(overlap_pert%crtn_orbs(i))
                end do
                ! Create the rest of the creation-related
                ! components of the overlap_pert object.
                call init_perturbation_creation(overlap_pert)

            case default
                call report("Keyword "//trim(w)//" not recognized in kp-fciqmc block", .true.)
            end select
        end do read_inp

        if (.not. vary_niters) then
            kp%niters = niters_temp
            kp%niters(kp%nvecs) = 0
        end if

    end subroutine kp_fciqmc_read_inp

    subroutine init_kp_fciqmc(kp)

        use SystemData, only: G1

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: ierr, krylov_vecs_memory, krylov_ht_memory, matrix_memory
        character(2) :: mem_fmt
        character (len=*), parameter :: t_r = "init_kp_fciqmc"

        ! Checks.
        if (.not. tHashWalkerList) call stop_all('t_r','kp-fciqmc can only be run using &
            &the linscalefcimcalgo option (the linear scaling algorithm).')
        if (.not. tUseRealCoeffs) call stop_all('t_r','kp-fciqmc can only be run using &
            &real coefficients).')
        if (tExactHamil .and. nProcessors /= 1) call stop_all('t_r','The exact-hamil &
            &option can only be used when running with one processor.')
        if (theisenberg .and. tAllSymSectors) call stop_all('t_r','The option to use all &
            &symmetry sectors at once has not been implemented with the Heisenberg model.')
        if (n_int == 4) call stop_all('t_r', 'Use of RealCoefficients does not work with 32 bit &
             &integers due to the use of the transfer operation from dp reals to 64 bit integers.')

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
        lenof_sign_kp = lenof_sign*kp%nvecs
        ! The total length of a bitstring containing all Krylov vectors.
        ! Note that the 1 is added because we store the diagonal Hamiltonian element in the array.
        NIfTotKP = NIfDBO + lenof_sign_kp + 1 + NIfFlag

        ! Allocate all of the KP arrays.
        nhashes_kp = nWalkerHashes
        TotWalkersKP = 0
        krylov_vecs_length = nint(MaxWalkersUncorrected*memory_factor_kp*kp%nvecs)
        write(kp_length_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(krylov_vecs_length)+1)))
        nkrylov_amp_elems_tot = lenof_sign*kp%nvecs*krylov_vecs_length

        ! Allocate the krylov_vecs array.
        ! The number of MB of memory required to allocate krylov_vecs.
        krylov_vecs_memory = krylov_vecs_length*(NIfTotKP+1)*size_n_int/1000000
        write(mem_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(krylov_vecs_memory)+1)))
        write(6,'(a73,1x,'//mem_fmt//')') "About to allocate array to hold all Krylov vectors. &
                                       &Memory required (MB):", krylov_vecs_memory
        write(6,'(a13)',advance='no') "Allocating..."
        call neci_flush(6)
        allocate(krylov_vecs(0:NIfTotKP, krylov_vecs_length), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating krylov_vecs array.")
        else
            write(6,'(1x,a5)') "Done."
        end if
        call neci_flush(6)
        krylov_vecs = 0_n_int

        ! Allocate the hash table to krylov_vecs.
        ! The number of MB of memory required to allocate krylov_vecs_ht.
        ! Each node requires 16 bytes.
        krylov_ht_memory = nhashes_kp*16/1000000
        write(mem_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(krylov_ht_memory)+1)))
        write(6,'(a78,1x,'//mem_fmt//')') "About to allocate hash table to the Krylov vector array. &
                                       &Memory required (MB):", krylov_ht_memory
        write(6,'(a13)',advance='no') "Allocating..."
        call neci_flush(6)
        allocate(krylov_vecs_ht(nhashes_kp), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(t_r, "Error allocating krylov_vecs array.")
        else
            write(6,'(1x,a5)') "Done."
            write(6,'(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                              &increase as further nodes are added."
        end if

        call init_hash_table(krylov_vecs_ht)

        if (tHamilOnFly) then
            allocate(SpawnVecKP(0:NIfTot,MaxSpawned),stat=ierr)
            SpawnVecKP(:,:) = 0_n_int
            SpawnedPartsKP => SpawnVecKP
            ! Do one extra iteration so that the Hamiltonian can be calculated for the final Krylov vector.
            kp%niters(kp%nvecs) = 1
        else
            allocate(SpawnVecKP(0:NOffSgn+lenof_sign_kp-1,MaxSpawned),stat=ierr)
            allocate(SpawnVecKP2(0:NOffSgn+lenof_sign_kp-1,MaxSpawned),stat=ierr)
            SpawnVecKP(:,:) = 0_n_int
            SpawnVecKP2(:,:) = 0_n_int
            SpawnedPartsKP => SpawnVecKP
            SpawnedPartsKP2 => SpawnVecKP2
            if (tSemiStochastic) then
                allocate(partial_determ_vecs_kp(lenof_sign_kp,determ_proc_sizes(iProcIndex)), stat=ierr)
                allocate(full_determ_vecs_kp(lenof_sign_kp,determ_space_size), stat=ierr)
                partial_determ_vecs_kp = 0.0_dp
                full_determ_vecs_kp = 0.0_dp
            end if
        end if

        if (tStoreKPMatrices) then
            ! (2*kp%nrepeats+16) arrays with (kp%nvecs**2) 8-byte elements each.
            matrix_memory = (2*kp%nrepeats+16)*(kp%nvecs**2)*8/1000000
        else
            ! 18 arrays with (kp%nvecs**2) 8-byte elements each.
            matrix_memory = 18*(kp%nvecs**2)*8/1000000
        end if
        write(mem_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(matrix_memory)+1)))
        write(6,'(a66,1x,'//mem_fmt//')') "About to allocate various subspace matrices. &
                                       &Memory required (MB):", matrix_memory
        write(6,'(a13)',advance='no') "Allocating..."
        call neci_flush(6)

        if (tOverlapPert) then
            allocate(pert_overlaps(kp%nvecs))
            allocate(kp_all_pert_overlaps(kp%nvecs))
        end if

        if (tStoreKPMatrices) then
            allocate(kp%overlap_matrices(kp%nvecs, kp%nvecs, kp%nrepeats), stat=ierr)
            allocate(kp%hamil_matrices(kp%nvecs, kp%nvecs, kp%nrepeats), stat=ierr)
        else
            allocate(kp%overlap_matrices(kp%nvecs, kp%nvecs, 1), stat=ierr)
            allocate(kp%hamil_matrices(kp%nvecs, kp%nvecs, 1), stat=ierr)
        end if

        allocate(kp_mpi_matrices_in(2*kp%nvecs, kp%nvecs))
        allocate(kp_mpi_matrices_out(2*kp%nvecs, kp%nvecs))

        allocate(kp_hamil_mean(kp%nvecs, kp%nvecs))
        allocate(kp_overlap_mean(kp%nvecs, kp%nvecs))
        allocate(kp_hamil_se(kp%nvecs, kp%nvecs))
        allocate(kp_overlap_se(kp%nvecs, kp%nvecs))

        allocate(kp_overlap_eigv(kp%nvecs))
        allocate(kp_init_overlaps(kp%nvecs))
        allocate(kp_overlap_eigenvecs(kp%nvecs, kp%nvecs))
        allocate(kp_transform_matrix(kp%nvecs, kp%nvecs))
        allocate(kp_inter_hamil(kp%nvecs, kp%nvecs))
        allocate(kp_eigenvecs_krylov(kp%nvecs, kp%nvecs))

        write(6,'(1x,a5)') "Done."

        ! If av_mc_excits_kp hasn't been set by the user, just use AvMCExcits.
        if (av_mc_excits_kp == 0.0_dp) av_mc_excits_kp = AvMCExcits

        MaxSpawnedEachProc = int(0.88*real(MaxSpawned,dp)/nProcessors)

        call WriteFciMCStatsHeader()
        call MPIBarrier(ierr)

    end subroutine init_kp_fciqmc

    subroutine init_kp_fciqmc_repeat(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        character(2) :: int_fmt

        write(int_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(kp%irepeat)+1)))
        write(6,'(1x,a22,1x,'//int_fmt//')') "Starting repeat number", kp%irepeat

        ! If starting from multiple POPSFILEs then set this counter so that the
        ! correct POPSFILE is read in this time. To read in POPSFILE.x,
        ! iPopsFileNoRead needs to be set to -x-1. We want to read in POPSFILE
        ! numbers 0 to kp%nconfigs-1
        if (tMultiplePopStart) iPopsFileNoRead = -(kp%iconfig-1)-1

        if (tOverlapPert .and. kp%irepeat == 1) then
            pert_overlaps = 0.0_dp
            call create_overlap_pert_vec()
        end if

        call create_initial_config(kp)

        call reset_hash_table(krylov_vecs_ht)
        krylov_vecs = 0_n_int

        ! Rezero all the necessary data.
        kp%overlap_matrix = 0.0_dp
        kp%hamil_matrix = 0.0_dp
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
        ! Setting this variable to true stops the shift from varying instantly.
        tSinglePartPhase = .true.

    end subroutine init_kp_fciqmc_repeat

    subroutine init_kp_fciqmc_iter(iter_data, determ_index)

        use FciMCData, only: FreeSlot, iStartFreeSlot, iEndFreeSlot

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

    end subroutine init_kp_fciqmc_iter

    subroutine create_initial_config(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: DetHash, nwalkers_int
        integer :: i
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign(lenof_sign), TotPartsCheck(lenof_sign)
        real(dp) :: nwalkers_target
        real(dp) :: norm, all_norm
        real(sp) :: total_time_before, total_time_after
        logical :: tCoreDet
        character(len=*), parameter :: t_r = "create_init_config"

        ! Clear everything from any previous repeats or starting configurations.
        call reset_hash_table(HashIndex)

        if (tStartSinglePart) then
            nwalkers_target = real(InitialPart,dp)
        else
            nwalkers_target = InitWalkers*nProcessors
        end if

        if (kp%irepeat == 1) then
            if (.not. tFiniteTemp) then
                if (tReadPops) then
                    ! Call a wrapper function which will call the various functions
                    ! required to read in a popsfile.
                    call read_popsfile_wrapper(pops_pert)

                    if (tScalePopulation) then
                        call scale_population(CurrentDets, TotWalkers, nwalkers_target, TotPartsCheck, scaling_factor)
                        ! Update global data.
                        if (any(abs(TotPartsCheck-TotParts) > 1.0e-12)) then
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
                    call copy_core_dets_this_proc_to_spawnedparts()
                    ! Any core space determinants which are not already in CurrentDets will be added
                    ! by this routine.
                    call add_core_states_currentdet_hash()
                    SpawnedParts = 0_n_int
                    if (tStartCoreGroundState .and. (.not. tReadPops)) &
                        call start_walkers_from_core_ground(tPrintInfo = .false.)
                end if

            else if (tFiniteTemp) then
                ! Convert the initial number of walkers to an integer. Note that on multiple
                ! processors this may round up the requested number of walkers slightly.
                nwalkers_target = ceiling(nwalkers_target/real(nProcessors,dp))

                ! If requested, reset the random number generator with the requested seed
                ! before creating the random initial configuration.
                if (tUseInitConfigSeeds) call dSFMT_init((iProcIndex+1)*init_config_seeds(kp%iconfig))

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
        else if (kp%irepeat > 1) then
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
                CurrentDets(0:NOffSgn+lenof_sign-1,i) = krylov_vecs(0:NOffSgn+lenof_sign-1,i)
                ! Copy across the flags.
                CurrentDets(NIfTot,i) = krylov_vecs(NIfTotKP,i)
            end do
            call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, int(TotWalkers, sizeof_int))
        end if

        ! Calculate and store the diagonal element of the Hamiltonian for determinants in CurrentDets.
        call fill_in_CurrentH()

        ! If starting from this configuration more than once, store the relevant data for next time.
        if (kp%nrepeats > 1 .and. kp%irepeat == 1) then
            HolesInList = 0
            do i = 1, TotWalkers
                int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,i)
                call extract_sign (CurrentDets(:,i), real_sign)
                tCoreDet = check_determ_flag(CurrentDets(:,i))
                ! Don't add unoccpied determinants, unless they are core determinants.
                if (IsUnoccDet(real_sign) .and. (.not. tCoreDet)) HolesInList = HolesInList + 1
            end do
            TotWalkersInit = TotWalkers - HolesInList
            call MPIAllReduce(TotWalkersInit, MPI_SUM, AllTotWalkers)
            AllTotWalkersInit = AllTotWalkers
            TotPartsInit = TotParts
            AllTotPartsInit = AllTotParts
        end if

    end subroutine create_initial_config

    subroutine create_overlap_pert_vec()

        ! Read in the popsfile and apply perturbation operator overlap_pert.

        integer :: mem_reqd, ierr
        character(2) :: mem_fmt
        character(len=*), parameter :: t_r = "create_overlap_pert_vec"

        if (allocated(perturbed_ground)) deallocate(perturbed_ground)

        ! Once this is finished, the vector that we want will be stored in
        ! CurrentDets. The total number of determinants will be TotWalkers.
        call read_popsfile_wrapper(overlap_pert)

        ! Print info about memory usage to the user.
        ! Memory required in MB.
        mem_reqd = TotWalkers*(NIfTotKP+1)*size_n_int/1000000
        write(mem_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(mem_reqd)+1)))

        write(6,'(a73,1x,'//mem_fmt//')') "About to allocate array to hold the perturbed &
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

        integer, intent(in) :: nwalkers
        real(dp) :: nwalkers_per_site_init
        integer :: i, ireplica, excit, nattempts, DetHash
        integer :: nspawns, ndets
        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        real(dp) :: r, walker_amp, walker_sign(lenof_sign)
        logical :: tInitiatorTemp
        type(fcimc_iter_data) :: unused_data
        integer(n_int), pointer :: PointTemp(:,:)

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
            if (r < 0.5) walker_amp = -1.0_dp*walker_amp

            do ireplica = 1, inum_runs
                walker_sign = 0.0_dp
                walker_sign(ireplica) = walker_amp
                call create_particle(nI, ilut, walker_sign, 0, ireplica)
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
            walker_sign = transfer(CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, i), walker_sign)
            TotParts = TotParts + abs(walker_sign)
        end do
        TotPartsOld = TotParts

        ! Add the entries into the hash table.
        call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, ndets)

        call MPIReduce(TotParts, MPI_SUM, AllTotParts)
        AllTotPartsOld = AllTotParts
        TotWalkers = int(ndets, int64)

        if (tSemiStochastic) then
            ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
            ! These routines will do this.
            call copy_core_dets_this_proc_to_spawnedparts()
            call add_core_states_currentdet_hash()
        end if

        ValidSpawnedList = InitialSpawnedSlots
        SpawnedParts = 0_n_int

        ! Turn the initiator method back on, if it was turned off at the start of this routine.
        tTruncInitiator = tInitiatorTemp

    end subroutine generate_init_config_basic

    subroutine generate_init_config_this_proc(nwalkers, nwalkers_per_site_init, tOccupyDetermSpace)

        ! This routine will distribute nwalkers walkers uniformly across all possible determinants.

        integer, intent(in) :: nwalkers
        real(dp) :: nwalkers_per_site_init
        logical, intent(in) :: tOccupyDetermSpace
        integer :: proc, excit, nattempts, idet, DetHash, det_ind, nI(nel)
        integer :: ideterm, ndets
        integer(n_int) :: ilut(0:NIfTot), int_sign(lenof_sign)
        type(ll_node), pointer :: prev, temp_node
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        real(dp) :: new_sign(lenof_sign), r
        logical :: tDetFound, tDetermAllOccupied
        character(len=*), parameter :: this_routine = "generate_init_config_this_proc"

        ilut = 0_n_int
        ndets = 0
        TotParts = 0.0_dp

        ideterm = 0
        tDetermAllOccupied = .false.

        do
            ! If using the tOccupyDetermSpace option then we want to put walkers
            ! on states in the deterministic space first.
            ! If not, or if we have finished doing this, generate determinants
            ! randoml and uniformly.
            if (tOccupyDetermSpace .and. (.not. tDetermAllOccupied)) then
                ideterm = ideterm + 1
                ilut = core_space(:,ideterm + determ_proc_indices(iProcIndex))
                call decode_bit_det(nI, ilut)
                ! If we have now occupied all deterministic states.
                if (ideterm == determ_proc_sizes(iProcIndex)) tDetermAllOccupied = .true.
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
            if (r < 0.5) real_sign_1 = -1.0_dp*real_sign_1
            int_sign = transfer(real_sign_1, int_sign)

            tDetFound = .false.
            DetHash = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(DetHash)
            prev => null()
            ! If the first element in the list for this hash value has been used.
            if (.not. temp_node%ind == 0) then
                ! Loop over all determinants with this hash value which are already in the list.
                do while (associated(temp_node))
                    if (DetBitEQ(CurrentDets(:,temp_node%ind), ilut, NIfDBO)) then
                        ! This determinant is already in the list.
                        tDetFound = .true.
                        int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, temp_node%ind)
                        real_sign_2 = transfer(int_sign, real_sign_2)
                        new_sign = real_sign_1 + real_sign_2
                        call encode_sign(CurrentDets(:, temp_node%ind), new_sign)
                        TotParts = TotParts - abs(real_sign_2) + abs(new_sign)
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    prev => temp_node
                    temp_node => temp_node%next
                end do

                if (.not. tDetFound) then
                    ! We need to add a new determinant in the next position in the list.
                    ! So create that next position!
                    allocate(prev%next)
                    temp_node => prev%next
                    nullify(temp_node%next)
                end if
            end if

            if (.not. tDetFound) then
                ! A new determiant needs to be added.
                ndets = ndets + 1
                TotParts = TotParts + abs(real_sign_1)
                det_ind = ndets
                temp_node%ind = det_ind
                ! Copy determinant data across.
                CurrentDets(0:NIfDBO, det_ind) = ilut(0:NIfDBO)
                CurrentDets(NOffSgn:NOffSgn+lenof_sign-1, det_ind) = int_sign
                if (tUseFlags) CurrentDets(NOffFlag, det_ind) = 0_n_int

                nullify(temp_node)
                nullify(prev)
            end if

            if (TotParts(1) >= nwalkers) exit

        end do

        call MPIReduce(TotParts, MPI_SUM, AllTotParts)
        TotPartsOld = TotParts
        AllTotPartsOld = AllTotParts
        TotWalkers = int(ndets, int64)

        if (tSemiStochastic) then
            ! Always need the core determinants to be at the top of CurrentDets, even when unoccupied.
            ! These routines will do this.
            call copy_core_dets_this_proc_to_spawnedparts()
            call add_core_states_currentdet_hash()
        else
            ! Remove the nodes of all unoccupied determinants from the hash table.
            ! For semi-stochastic calculations, this is done in add_core_states_currentdet_hash.
            do idet = 1, ndets
                call extract_sign(CurrentDets(:,idet),real_sign_1)
                if (.not. IsUnoccDet(real_sign_1)) cycle
                tDetFound = .false.
                call decode_bit_det(nI, CurrentDets(:, idet))
                DetHash = FindWalkerHash(nI, nWalkerHashes)
                temp_node => HashIndex(DetHash)
                prev => null()
                if (.not. temp_node%ind == 0) then
                    ! Loop over all determinants with this hash value.
                    do while (associated(temp_node))
                        if (temp_node%ind == idet) then
                            tDetFound = .true.
                            call remove_node(prev, temp_node)
                            exit
                        end if
                        ! Move on to the next determinant with this hash value.
                        prev => temp_node
                        temp_node => temp_node%next
                    end do
                end if
                ASSERT(tDetFound)
            end do
        end if

        SpawnedParts = 0_n_int

    end subroutine generate_init_config_this_proc

    subroutine scale_population(walker_list, ndets, target_pop, input_pop, scaling_factor)

        ! Take an input list of walkers, find the total walker population in
        ! the list, and then multiply all the walker signs by some factor in
        ! order for the walker list to have the requested target population.

        integer(n_int), intent(inout) :: walker_list(:,:)
        integer, intent(in) :: ndets
        real(dp), intent(in) :: target_pop
        real(dp), intent(out) :: input_pop(lenof_sign)
        real(dp), intent(out) :: scaling_factor

        integer :: i
        real(dp) :: real_sign(lenof_sign), all_input_pop(lenof_sign)

        input_pop = 0.0_dp

        ! First, find the population of the walkers in walker_list.
        do i = 1, ndets
            call extract_sign(walker_list(:,i), real_sign)
            input_pop = input_pop + abs(real_sign)
        end do

        call MPIAllReduce(input_pop, MPI_SUM, all_input_pop)

        ! Just use the first particle type to determine the scaing factor.
        scaling_factor = target_pop/all_input_pop(1)

        ! Now multiply the walker signs by scaling_factor.
        do i = 1, ndets
            call extract_sign(walker_list(:,i), real_sign)
            real_sign = real_sign*scaling_factor
            call encode_sign(walker_list(:,i), real_sign)
        end do

    end subroutine scale_population

    subroutine store_krylov_vec(kp)

        type(kp_fciqmc_data), intent(in) :: kp
        integer :: idet, iamp, sign_ind, hdiag_ind, flag_ind, DetHash, det_ind
        integer :: nI(nel)
        integer(n_int) :: temp, int_sign(lenof_sign)
        logical :: tDetFound, tCoreDet
        real(dp) :: amp_fraction, real_sign(lenof_sign)
        type(ll_node), pointer :: temp_node, prev
        character(2) :: int_fmt
        character(len=*), parameter :: t_r = "store_krylov_vec"

        write(6,'(a71)',advance='no') "# Adding the current walker configuration to the Krylov vector array..."
        call neci_flush(6)

        ! The index of the first element referring to the sign, for this ivec.
        sign_ind = NIfDBO + lenof_sign*(kp%ivec-1) + 1
        hdiag_ind = NIfDBO + lenof_sign_kp + 1
        if (tUseFlags) flag_ind = NIfDBO + lenof_sign_kp + 2

        ! Loop over all occupied determinants for this new Krylov vector.
        do idet = 1, TotWalkers
            int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,idet)
            call extract_sign (CurrentDets(:,idet), real_sign)
            tCoreDet = check_determ_flag(CurrentDets(:,idet))
            ! Don't add unoccpied determinants, unless they are core determinants.
            if (IsUnoccDet(real_sign) .and. (.not. tCoreDet)) cycle
            tDetFound = .false.
            call decode_bit_det(nI, CurrentDets(:,idet))
            DetHash = FindWalkerHash(nI, nhashes_kp)
            temp_node => krylov_vecs_ht(DetHash)

            ! If the first element in the list for this hash value has been used.
            if (.not. temp_node%ind == 0) then
                ! Loop over all determinants with this hash value which are already in the list.
                do while (associated(temp_node))
                    if (DetBitEQ(CurrentDets(:,idet), krylov_vecs(:,temp_node%ind), NIfDBO)) then
                        ! This determinant is already in the list.
                        det_ind = temp_node%ind
                        ! Add the amplitude for the new Krylov vector. The determinant and flag are
                        ! there already.
                        krylov_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = int_sign
                        tDetFound = .true.
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    prev => temp_node
                    temp_node => temp_node%next
                end do

                if (.not. tDetFound) then
                    ! We need to add a new determinant in the next position in the list.
                    ! So create that next position!
                    allocate(prev%next)
                    temp_node => prev%next
                    nullify(temp_node%next)
                end if
            end if

            if (.not. tDetFound) then
                ! A new determiant needs to be added.
                TotWalkersKP = TotWalkersKP + 1
                det_ind = TotWalkersKP
                temp_node%ind = det_ind

                if (TotWalkersKP > krylov_vecs_length) then
                    call stop_all(t_r, "There are no slots left in the krylov_vecs array for the next determinant. &
                                       &You can increase the size of this array using the memory-factor option in &
                                       &the kp-fciqmc block of the input file.")
                end if

                ! Copy determinant data across.
                krylov_vecs(0:NIfDBO,det_ind) = CurrentDets(0:NIfDBO,idet)
                krylov_vecs(sign_ind:sign_ind+lenof_sign-1,det_ind) = int_sign
                krylov_vecs(hdiag_ind,det_ind) = transfer(CurrentH(1,idet), temp)
                if (tUseFlags) krylov_vecs(flag_ind,det_ind) = CurrentDets(NOffFlag,idet)
            end if

            nullify(temp_node)
            nullify(prev)
            do iamp = 1, lenof_sign
                if (real_sign(iamp) /= 0_n_int) then
                    nkrylov_amp_elems_used = nkrylov_amp_elems_used + 1
                end if
            end do

        end do

        write(6,'(1x,a5)',advance='yes') "Done."
        write(int_fmt,'(a1,i1)') "i", ceiling(log10(real(abs(TotWalkersKP)+1)))
        write(6,'(a56,1x,'//int_fmt//',1x,a17,1x,'//kp_length_fmt//')') "# Number unique determinants in the Krylov &
                                             &vector array:", TotWalkersKP, "out of a possible", krylov_vecs_length
        amp_fraction = real(nkrylov_amp_elems_used,dp)/real(nkrylov_amp_elems_tot,dp)
        write(6,'(a69,1x,es10.4)') "# Fraction of the amplitude elements used in the Krylov vector array:", amp_fraction
        call neci_flush(6)

    end subroutine store_krylov_vec

    subroutine calc_overlap_matrix_elems(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, jvec, ind(kp%ivec)
        type(ll_node), pointer :: temp_node
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)

        associate(s_matrix => kp%overlap_matrix, ivec => kp%ivec)

            ! Just in case!
            s_matrix(1:ivec, ivec) = 0.0_dp
            s_matrix(ivec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in krylov_vecs, for each Krylov vector.
                ind(jvec) = NIfDBO + lenof_sign*(jvec-1) + 1
            end do

            ! Loop over all determinants in krylov_vecs.
            do idet = 1, TotWalkersKP
                int_sign = krylov_vecs(ind(ivec):ind(ivec)+lenof_sign-1, idet)
                real_sign_1 = transfer(int_sign, real_sign_1)
                if (IsUnoccDet(real_sign_1)) cycle
                ! Loop over all Krylov vectors currently stored.
                do jvec = 1, ivec
                    int_sign = krylov_vecs(ind(jvec):ind(jvec)+lenof_sign-1, idet)
                    real_sign_2 = transfer(int_sign, real_sign_1)
                    if (IsUnoccDet(real_sign_2)) cycle
#ifdef __DOUBLERUN
                    s_matrix(jvec,ivec) = s_matrix(jvec,ivec) + &
                        (real_sign_1(1)*real_sign_2(2) + real_sign_1(2)*real_sign_2(1))/2.0_dp
#else
                    s_matrix(jvec,ivec) = s_matrix(jvec,ivec) + real_sign_1(1)*real_sign_2(1)
#endif
                end do
            end do

            ! Fill in the lower-half of the overlap matrix.
            do jvec = 1, ivec
                s_matrix(ivec,jvec) = s_matrix(jvec,ivec)
            end do

        end associate

    end subroutine calc_overlap_matrix_elems

    subroutine calc_perturbation_overlap(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, sign_ind, DetHash, det_ind
        integer :: nI(nel)
        type(ll_node), pointer :: temp_node
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        real(dp) :: overlap
        logical :: tDetFound

        associate(s_matrix => kp%overlap_matrix, ivec => kp%ivec)

            overlap = 0.0_dp

            sign_ind = NIfDBO + lenof_sign*(ivec-1) + 1

            ! Loop over all determinants in perturbed_ground
            do idet = 1, size(perturbed_ground, dim=2)
                call extract_sign(perturbed_ground(:,idet), real_sign_1)
                if (IsUnoccDet(real_sign_1)) cycle

                ! Search to see if this determinant is in the Krylov vector.
                call decode_bit_det(nI, perturbed_ground(:,idet))
                DetHash = FindWalkerHash(nI, nhashes_kp)
                ! Point to the first node with this hash value in krylov_vecs.
                temp_node => krylov_vecs_ht(DetHash)
                if (temp_node%ind == 0) then
                    ! If there are no determinants at all with this hash value in krylov_vecs.
                    cycle
                else
                    tDetFound = .false.
                    do while (associated(temp_node))
                        if (DetBitEQ(perturbed_ground(:,idet), krylov_vecs(:,temp_node%ind), NIfDBO)) then
                            ! If this determinant has been found in krylov_vecs.
                            det_ind = temp_node%ind
                            tDetFound = .true.
                            exit
                        end if
                        ! Move on to the next determinant with this hash value.
                        temp_node => temp_node%next
                    end do
                    if (tDetFound) then
                        int_sign = krylov_vecs(sign_ind:sign_ind+lenof_sign-1, det_ind)
                        real_sign_2 = transfer(int_sign, real_sign_1)
#ifdef __DOUBLERUN
                        overlap = overlap + (real_sign_1(1)*real_sign_2(2) + real_sign_1(2)*real_sign_2(1))/2.0_dp
#else
                        overlap = overlap + real_sign_1(1)*real_sign_2(1)
#endif
                    end if
                end if
            end do

            pert_overlaps(ivec) = pert_overlaps(ivec) + overlap

        end associate

    end subroutine calc_perturbation_overlap

    subroutine calc_hamil_on_fly(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, jvec, ind(kp%ivec), nI(nel)
        integer :: det_ind, hdiag_ind, flag_ind, ideterm, DetHash
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        real(dp) :: temp
        type(ll_node), pointer :: temp_node
        logical :: tDetFound, tDeterm
        character(len=*), parameter :: t_r = "calc_hamil_elems_on_fly"

        associate(h_matrix => kp%hamil_matrix, ivec => kp%ivec)

            h_matrix(1:ivec, ivec) = 0.0_dp
            h_matrix(ivec, 1:ivec) = 0.0_dp

            do jvec = 1, ivec
                ! The first index of the sign in krylov_vecs, for each Krylov vector.
                ind(jvec) = NIfDBO + lenof_sign*(jvec-1) + 1
            end do
            hdiag_ind = NIfDBO + lenof_sign_kp + 1
            if (tUseFlags) flag_ind = NIfDBO + lenof_sign_kp + 2

            ideterm = 0

            ! Loop over all determinants in SpawnedPartsKP.
            do idet = 1, max_spawned_ind
                int_sign = SpawnedPartsKP(NOffSgn:NOffSgn+lenof_sign-1, idet)
                real_sign_1 = transfer(int_sign, real_sign_1)
                call decode_bit_det(nI, SpawnedPartsKP(:,idet))
                DetHash = FindWalkerHash(nI, nhashes_kp)
                ! Point to the first node with this hash value in krylov_vecs.
                temp_node => krylov_vecs_ht(DetHash)
                if (temp_node%ind == 0) then
                    ! If there are no determinants at all with this hash value in krylov_vecs.
                    cycle
                else
                    tDetFound = .false.
                    do while (associated(temp_node))
                        if (DetBitEQ(SpawnedPartsKP(:,idet), krylov_vecs(:,temp_node%ind), NIfDBO)) then
                            ! If this determinant has been found in krylov_vecs.
                            det_ind = temp_node%ind
                            tDetFound = .true.
                            exit
                        end if
                        ! Move on to the next determinant with this hash value.
                        temp_node => temp_node%next
                    end do
                    if (tDetFound) then
                        ! Add in the contribution to the projected Hamiltonian, for each Krylov vector.
                        do jvec = 1, ivec
                            int_sign = krylov_vecs(ind(jvec):ind(jvec)+lenof_sign-1, det_ind)
                            real_sign_2 = transfer(int_sign, real_sign_1)
                            if (IsUnoccDet(real_sign_2)) cycle
#ifdef __DOUBLERUN
                            h_matrix(jvec,ivec) = h_matrix(jvec,ivec) - &
                                (real_sign_1(1)*real_sign_2(2) + real_sign_1(2)*real_sign_2(1))/2.0_dp
#else
                            h_matrix(jvec,ivec) = h_matrix(jvec,ivec) - real_sign_1(1)*real_sign_2(1)
#endif
                        end do
                    end if
                end if
            end do

            ! Loop over all determinants in krylov_vecs.
            do idet = 1, TotWalkersKP

                real_sign_1 = 0.0_dp
                tDeterm = .false.
                if (tUseFlags) then
                    tDeterm = btest(krylov_vecs(flag_ind, idet), flag_deterministic + flag_bit_offset)
                end if

                if (tDeterm) then
                    ideterm  = ideterm + 1
                    real_sign_1 = - partial_determ_vector(:,ideterm) + &
                             (DiagSft+Hii) * tau * full_determ_vector(:, ideterm + determ_proc_indices(iProcIndex))
                else
                    int_sign = krylov_vecs(ind(ivec):ind(ivec)+lenof_sign-1, idet)
                    real_sign_1 = transfer(int_sign, real_sign_1)
                    real_sign_1 = tau * real_sign_1 * (transfer(krylov_vecs(hdiag_ind, idet), temp) + Hii)
                end if
                if (IsUnoccDet(real_sign_1)) cycle

                ! Loop over all Krylov vectors currently stored.
                do jvec = 1, ivec
                    int_sign = krylov_vecs(ind(jvec):ind(jvec)+lenof_sign-1, idet)
                    real_sign_2 = transfer(int_sign, real_sign_1)
                    if (IsUnoccDet(real_sign_2)) cycle
#ifdef __DOUBLERUN
                    h_matrix(jvec,ivec) = h_matrix(jvec,ivec) + &
                        (real_sign_1(1)*real_sign_2(2) + real_sign_1(2)*real_sign_2(1))/2.0_dp
#else
                    h_matrix(jvec,ivec) = h_matrix(jvec,ivec) + real_sign_1(1)*real_sign_2(1)
#endif
                end do
            end do

            do jvec = 1, ivec
                h_matrix(jvec,ivec) = h_matrix(jvec,ivec)/tau
                h_matrix(ivec,jvec) = h_matrix(jvec,ivec)
            end do

            if (tSemiStochastic) then
                if (ideterm /= determ_proc_sizes(iProcIndex)) then
                    write(6,*) "determ_proc_sizes(iProcIndex):", determ_proc_sizes(iProcIndex)
                    write(6,*) "ideterm:", ideterm
                    call neci_flush(6)
                    call stop_all(t_r, "An incorrect number of core determinants have been counted.")
                end if
            end if

        end associate

    end subroutine calc_hamil_on_fly

    subroutine calc_hamil_exact(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: i, j, idet, jdet, ic, hdiag_ind
        integer(n_int) :: ilut_1(0:NIfTot), ilut_2(0:NIfTot)
        integer(n_int) :: int_sign(lenof_sign_kp)
        integer :: nI(nel), nJ(nel)
        real(dp) :: real_sign_1(lenof_sign_kp), real_sign_2(lenof_sign_kp)
        real(dp) :: h_elem
        logical :: any_occ, occ_1, occ_2
        integer(4), allocatable :: occ_flags(:)

        hdiag_ind = NIfDBO + lenof_sign_kp + 1

        associate(h_matrix => kp%hamil_matrix)

            h_matrix = 0.0_dp

            allocate(occ_flags(TotWalkersKP))
            occ_flags = 0

            ilut_1 = 0_n_int
            ilut_2 = 0_n_int

            ! Check to see if there are any replica 1 or 2 walkers on this determinant.
            do idet = 1, TotWalkersKP
                int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, idet)

                any_occ = .false.
#ifdef __DOUBLERUN
                do i = 1, kp%nvecs
                    any_occ = any_occ .or. (int_sign(2*i-1) /= 0)
                end do
                if (any_occ) occ_flags = ibset(occ_flags(idet), 0)

                any_occ = .false.
                do i = 1, kp%nvecs
                    any_occ = any_occ .or. (int_sign(2*i) /= 0)
                end do
                if (any_occ) occ_flags = ibset(occ_flags(idet), 1)
#else
                do i = 1, kp%nvecs
                    any_occ = any_occ .or. (int_sign(i) /= 0)
                end do
                if (any_occ) occ_flags = ibset(occ_flags(idet), 0)
#endif
            end do

            ! Loop over all determinants in krylov_vecs.
            do idet = 1, TotWalkersKP
                ilut_1 = krylov_vecs(0:NIfDBO, idet)
                call decode_bit_det(nI, ilut_1)
                int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, idet)
                real_sign_1 = transfer(int_sign, real_sign_1)
                occ_1 = btest(occ_flags(idet),0)
                occ_2 = btest(occ_flags(idet),1)

                do jdet = idet, TotWalkersKP
#ifdef __DOUBLERUN
                    if (.not. ((occ_1 .and. btest(occ_flags(jdet),1)) .or. &
                        (occ_2 .and. btest(occ_flags(jdet),0)))) cycle
#else
                    if (.not. (occ_1 .and. btest(occ_flags(jdet),0)) ) cycle
#endif
                    ilut_2 = krylov_vecs(0:NIfDBO, jdet)
                    ic = FindBitExcitLevel(ilut_1, ilut_2)
                    if (ic > 2) cycle

                    call decode_bit_det(nJ, ilut_2)
                    int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, jdet)
                    real_sign_2 = transfer(int_sign, real_sign_1)
                    
                    if (idet == jdet) then
                        h_elem = transfer(krylov_vecs(hdiag_ind, idet), h_elem) + Hii
                    else
                        if (tHPHF) then
                            h_elem = hphf_off_diag_helement(nI, nJ, ilut_1, ilut_2)
                        else
                            h_elem = get_helement(nI, nJ, ic, ilut_1, ilut_2)
                        end if
                    end if

                    ! Finally, add in the contribution to all of the Hamiltonian elements.
                    do i = 1, kp%nvecs
                        do j = i, kp%nvecs
#ifdef __DOUBLERUN
                            if (idet == jdet) then
                                h_matrix(i,j) = h_matrix(i,j) + &
                                    h_elem*(real_sign_1(2*i-1)*real_sign_2(2*j) + &
                                    real_sign_1(2*i)*real_sign_2(2*j-1))/2
                            else
                                h_matrix(i,j) = h_matrix(i,j) + &
                                    h_elem*(real_sign_1(2*i-1)*real_sign_2(2*j) + &
                                    real_sign_1(2*i)*real_sign_2(2*j-1) + &
                                    real_sign_1(2*j-1)*real_sign_2(2*i) + &
                                    real_sign_1(2*j)*real_sign_2(2*i-1))/2
                            end if
#else
                            if (idet == jdet) then
                                h_matrix(i,j) = h_matrix(i,j) + &
                                    h_elem*real_sign_1(i)*real_sign_2(j)
                            else
                                h_matrix(i,j) = h_matrix(i,j) + &
                                    h_elem*(real_sign_1(i)*real_sign_2(j) + &
                                    real_sign_1(j)*real_sign_2(i))
                            end if
#endif
                        end do
                    end do

                end do
            end do

            do i = 1, kp%nvecs
                do j = 1, i-1
                    h_matrix(i,j) = h_matrix(j,i)
                end do
            end do

            deallocate(occ_flags)

        end associate

    end subroutine calc_hamil_exact

    subroutine communicate_kp_matrices(kp)

        ! Add all the overlap and projected Hamiltonian matrices together, with the result being
        ! held only on the root node.

        type(kp_fciqmc_data), intent(inout) :: kp

        kp_mpi_matrices_in(1:kp%nvecs, 1:kp%nvecs) = kp%overlap_matrix
        kp_mpi_matrices_in(kp%nvecs+1:2*kp%nvecs, 1:kp%nvecs) = kp%hamil_matrix

        call MPISum(kp_mpi_matrices_in, kp_mpi_matrices_out)

        if (iProcIndex == root) then
            kp%overlap_matrix = kp_mpi_matrices_out(1:kp%nvecs, 1:kp%nvecs)
            kp%hamil_matrix = kp_mpi_matrices_out(kp%nvecs+1:2*kp%nvecs, 1:kp%nvecs)
        end if

    end subroutine communicate_kp_matrices

    subroutine output_kp_matrices_wrapper(kp)

        type(kp_fciqmc_data), intent(in) :: kp

        if (iProcIndex == root) then
            call output_kp_matrices(kp, 'hamil  ', kp%hamil_matrices)
            call output_kp_matrices(kp, 'overlap', kp%overlap_matrices)
        end if

    end subroutine output_kp_matrices_wrapper

    subroutine output_kp_matrices(kp, stem, matrices)

        type(kp_fciqmc_data), intent(in) :: kp
        character(7), intent(in) :: stem
        real(dp), intent(in) :: matrices(:,:,:)
        character(2) :: ifmt, jfmt
        character(25) :: ind1, ind2, filename
        integer :: i, j, k, ilen, jlen, temp_unit, repeat_ind

        write(ind1,'(i15)') kp%iconfig
        if (tStoreKPMatrices) then
            filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
            repeat_ind = kp%irepeat
        else
            write(ind2,'(i15)') kp%irepeat
            filename = trim(trim(stem)//'.'//trim(adjustl(ind1))//'.'//trim(adjustl(ind2)))
            repeat_ind = 1
        end if
        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        ! Write all the components of the various estimates of the matrix, above and including the
        ! diagonal, one after another on separate lines.
        do i = 1, kp%nvecs
            ilen = ceiling(log10(real(abs(i)+1)))
            ! ifmt will hold the correct integer length so that there will be no spaces printed out.
            ! Note that this assumes that ilen < 10, which is very reasonable!
            write(ifmt,'(a1,i1)') "i", ilen
            do j = i, kp%nvecs
                jlen = ceiling(log10(real(abs(j)+1)))
                write(jfmt,'(a1,i1)') "i", jlen

                ! Write the index of the matrix element.
                write(temp_unit,'(a1,'//ifmt//',a1,'//jfmt//',a1)',advance='no') "(",i,",",j,")"
                do k = 1, repeat_ind
                    write(temp_unit,'(1x,es19.12)',advance='no') matrices(i,j,k)
                end do
                write(temp_unit,'()',advance='yes')
            end do
        end do

        close(temp_unit)

    end subroutine output_kp_matrices

    subroutine average_kp_matrices_wrapper(kp)

        type(kp_fciqmc_data), intent(in) :: kp

        call average_kp_matrices(kp%nrepeats, kp%hamil_matrices, kp_hamil_mean, kp_hamil_se)
        call average_kp_matrices(kp%nrepeats, kp%overlap_matrices, kp_overlap_mean, kp_overlap_se)
        if (tOutputAverageKPMatrices) then
            call output_average_kp_matrix(kp, 'av_hamil  ', kp_hamil_mean, kp_hamil_se)
            call output_average_kp_matrix(kp, 'av_overlap', kp_overlap_mean, kp_overlap_se)
        end if

    end subroutine average_kp_matrices_wrapper

    subroutine average_kp_matrices(nrepeats, matrices, mean, se)

        integer, intent(in) :: nrepeats
        real(dp), intent(in) :: matrices(:,:,:)
        real(dp), intent(inout) :: mean(:,:), se(:,:)
        integer :: irepeat

        mean = 0.0_dp
        se = 0.0_dp

        do irepeat = 1, nrepeats
            mean = mean + matrices(:,:,irepeat)
        end do
        mean = mean/nrepeats

        if (nrepeats > 1) then
            do irepeat = 1, nrepeats
                se = se + (matrices(:,:,irepeat)-mean)**2
            end do
            se = se/((nrepeats-1)*nrepeats)
            se = sqrt(se)
        end if

    end subroutine average_kp_matrices

    subroutine output_average_kp_matrix(kp, stem, mean, se)

        type(kp_fciqmc_data), intent(in) :: kp
        character(10), intent(in) :: stem
        real(dp), intent(inout) :: mean(:,:), se(:,:)
        integer :: irepeat
        character(2) :: ifmt, jfmt
        character(25) :: ind1, filename
        integer :: i, j, k, ilen, jlen, temp_unit, repeat_ind

        write(ind1,'(i15)') kp%iconfig
        filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        ! Write all the components of the various estimates of the matrix, above and including the
        ! diagonal, one after another on separate lines.
        do i = 1, kp%nvecs
            ilen = ceiling(log10(real(abs(i)+1)))
            ! ifmt will hold the correct integer length so that there will be no spaces printed out.
            ! Note that this assumes that ilen < 10, which is very reasonable!
            write(ifmt,'(a1,i1)') "i", ilen
            do j = i, kp%nvecs
                jlen = ceiling(log10(real(abs(j)+1)))
                write(jfmt,'(a1,i1)') "i", jlen

                ! Write the index of the matrix element.
                write(temp_unit,'(a1,'//ifmt//',a1,'//jfmt//',a1)',advance='no') "(",i,",",j,")"
                if (kp%nrepeats > 1) then
                    ! Write the mean and standard error.
                    write(temp_unit,'(1x,es19.12,1x,a3,es19.12)') mean(i,j), "+/-", se(i,j)
                else if (kp%nrepeats == 1) then
                    ! If we only have one sample then a standard error was not calculated, so
                    ! only output the mean.
                    write(temp_unit,'(1x,es19.12)') mean(i,j)
                end if
            end do
        end do

        close(temp_unit)

    end subroutine output_average_kp_matrix

    subroutine average_and_communicate_pert_overlaps(nrepeats)

        integer, intent(in) :: nrepeats

        pert_overlaps = pert_overlaps/nrepeats
        call MPISum(pert_overlaps, kp_all_pert_overlaps)
        write(6,*) "overlaps:", kp_all_pert_overlaps

    end subroutine average_and_communicate_pert_overlaps

    subroutine find_and_output_lowdin_eigv(kp)

        type(kp_fciqmc_data), intent(in) :: kp
        integer :: lwork, counter, i, nkeep, nkeep_len, temp_unit
        integer :: npositive, info, ierr
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: kp_final_hamil(:,:), kp_hamil_eigv(:)
        real(dp) :: kp_pert_energy_overlaps(kp%nvecs)
        character(2) :: nkeep_fmt
        character(7) :: string_fmt
        character(25) :: ind1, filename
        character(len=*), parameter :: stem = "lowdin"

        write(ind1,'(i15)') kp%iconfig
        filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')
    
        if (tWritePopsNorm) then
            write(temp_unit, '(4("-"),a41,25("-"))') "Norm of unperturbed initial wave function"
            write(temp_unit,'(1x,es19.12,/)') sqrt(pops_norm)
        end if

        ! Create the workspace for the diagonaliser.
        lwork = max(1,3*kp%nvecs-1)
        allocate(work(lwork), stat=ierr)

        kp_overlap_eigenvecs = kp_overlap_mean

        ! Now perform the diagonalisation.
        call dsyev('V', 'U', kp%nvecs, kp_overlap_eigenvecs, kp%nvecs, kp_overlap_eigv, work, lwork, info)

        npositive = 0
        write(temp_unit,'(4("-"),a26,40("-"))') "Overlap matrix eigenvalues"
        do i = 1, kp%nvecs
            write(temp_unit,'(1x,es19.12)') kp_overlap_eigv(i)
            if (kp_overlap_eigv(i) > 0.0_dp) npositive = npositive + 1
        end do

        if (tOverlapPert) then
            write(temp_unit,'(/,"# The overlap of each final Hamiltonian eigenstate with each &
                              &requested perturbed ground state will be printed. The first &
                              &printed overlap is that with the first Krylov vector. The second &
                              &printed overlap is that with the vector specified with the &
                              &OVERLAP-PERTURB-ANNIHILATE and OVERLAP-PERTURB-CREATION &
                              &options.")')
        end if

        do nkeep = 1, npositive

            allocate(kp_final_hamil(nkeep,nkeep))
            allocate(kp_hamil_eigv(nkeep))

            associate(transform_matrix => kp_transform_matrix(1:kp%nvecs, 1:nkeep), &
                      inter_hamil => kp_inter_hamil(1:kp%nvecs, 1:nkeep), &
                      eigenvecs_krylov => kp_eigenvecs_krylov(1:kp%nvecs, 1:nkeep), &
                      init_overlaps => kp_init_overlaps(1:nkeep))

                counter = 0
                do i = kp%nvecs-nkeep+1, kp%nvecs
                    counter = counter + 1
                    transform_matrix(:,counter) = kp_overlap_eigenvecs(:, i)/sqrt(kp_overlap_eigv(i))
                end do

                inter_hamil = matmul(kp_hamil_mean, transform_matrix)
                kp_final_hamil = matmul(transpose(transform_matrix), inter_hamil)

                call dsyev('V', 'U', nkeep, kp_final_hamil, nkeep, kp_hamil_eigv, work, lwork, info)

                eigenvecs_krylov = matmul(transform_matrix, kp_final_hamil)
                init_overlaps = matmul(kp_overlap_mean(1,:), eigenvecs_krylov)/scaling_factor

                if (tOverlapPert) kp_pert_energy_overlaps(1:nkeep) = matmul(kp_all_pert_overlaps, eigenvecs_krylov)

                nkeep_len = ceiling(log10(real(abs(nkeep)+1)))
                write(nkeep_fmt,'(a1,i1)') "i", nkeep_len
                write(string_fmt,'(i2,a5)') 15-nkeep_len, '("-")'
                write(temp_unit,'(/,4("-"),a37,1x,'//nkeep_fmt//',1x,a12,'//string_fmt//')') &
                    "Eigenvalues and overlaps when keeping", nkeep, "eigenvectors"
                do i = 1, nkeep
                    write(temp_unit,'(1x,es19.12,1x,es19.12)',advance='no') kp_hamil_eigv(i), init_overlaps(i)
                    if (tOverlapPert) write(temp_unit,'(1x,es19.12)',advance='no') kp_pert_energy_overlaps(i)
                    write(temp_unit,'()',advance='yes')
                end do

            end associate

            deallocate(kp_final_hamil)
            deallocate(kp_hamil_eigv)

        end do

        deallocate(work)

        close(temp_unit)

    end subroutine find_and_output_lowdin_eigv

    subroutine find_and_output_gram_schmidt_eigv(kp)

        type(kp_fciqmc_data), intent(in) :: kp
        integer :: lwork, counter, nkeep, nkeep_len, temp_unit
        integer :: npositive, info, ierr
        integer :: i, j, n, m
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: kp_final_hamil(:,:), kp_hamil_eigv(:)
        character(2) :: nkeep_fmt
        character(7) :: string_fmt
        character(25) :: ind1, filename
        character(len=*), parameter :: stem = "gram_schmidt"

        write(ind1,'(i15)') kp%iconfig
        filename = trim(trim(stem)//'.'//trim(adjustl(ind1)))
        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        ! Create the workspace for the diagonaliser.
        lwork = max(1,3*kp%nvecs-1)
        allocate(work(lwork), stat=ierr)

        ! Use the following allocated arrays as work space for the following routine. Not ideal, I know...
        allocate(kp_final_hamil(kp%nvecs, kp%nvecs))
        associate(S_tilde => kp_inter_hamil, N => kp_init_overlaps)
            call construct_gram_schmidt_transform_matrix(kp_overlap_mean, kp_transform_matrix, S_tilde, kp_final_hamil, N, kp%nvecs)
            npositive = 0
            do i = 1, kp%nvecs
                if (N(i) > 0.0_dp) then
                    npositive = npositive + 1
                else
                    exit
                end if
            end do

            write(temp_unit,'(4("-"),a21,45("-"))') "Normalisation factors"
            do i = 1, npositive+1
                write(temp_unit,'(1x,es19.12)') N(i)
            end do
        end associate
        deallocate(kp_final_hamil)

        do nkeep = 1, npositive

            allocate(kp_final_hamil(nkeep,nkeep))
            allocate(kp_hamil_eigv(nkeep))

            associate(S => kp_transform_matrix(1:kp%nvecs, 1:nkeep), &
                      init_overlaps => kp_init_overlaps(1:nkeep))

                kp_final_hamil = 0.0_dp
                do n = 1, nkeep
                    do m = 1, nkeep
                        do i = 1, n
                            do j = 1, m
                                kp_final_hamil(n,m) = kp_final_hamil(n,m) + S(i,n)*kp_hamil_mean(i,j)*S(j,m)
                            end do
                        end do
                    end do
                end do

                call dsyev('V', 'U', nkeep, kp_final_hamil, nkeep, kp_hamil_eigv, work, lwork, info)

                init_overlaps = kp_final_hamil(:,1)*sqrt(kp_overlap_mean(1,1))

                nkeep_len = ceiling(log10(real(abs(nkeep)+1)))
                write(nkeep_fmt,'(a1,i1)') "i", nkeep_len
                write(string_fmt,'(i2,a5)') 15-nkeep_len, '("-")'
                write(temp_unit,'(/,4("-"),a37,1x,'//nkeep_fmt//',1x,a12,'//string_fmt//')') &
                    "Eigenvalues and overlaps when keeping", nkeep, "eigenvectors"
                do i = 1, nkeep
                    write(temp_unit,'(1x,es19.12,1x,es19.12)') kp_hamil_eigv(i), init_overlaps(i)
                end do

            end associate

            deallocate(kp_final_hamil)
            deallocate(kp_hamil_eigv)

        end do

        deallocate(work)

        close(temp_unit)

    end subroutine find_and_output_gram_schmidt_eigv

    subroutine construct_gram_schmidt_transform_matrix(overlap, S, S_tilde, k, N, matrix_size)

        ! Construct the matrix, S, which transforms the Krylov vectors to a set of orthogonal vectors.

        ! We use the notation from the Appendix of Phys. Rev. B. 85, 205119 for the matrices and the
        ! indices, except we use 'm' instead of 'n' and 'l' instead of 'k' for the indices.

        ! Usage: overlap on input should contain the overlap matrix of the Krylov vectors that
        ! you want to produce a transformation matrix for. All other matrices input should be
        ! allocated on input, and should be the same size as overlap.

        real(dp), intent(inout) :: overlap(:,:), S(:,:), S_tilde(:,:), k(:,:), N(:)
        integer, intent(in) :: matrix_size
        integer :: m, i, l, p, q

        S = 0.0_dp
        S_tilde = 0.0_dp
        k = 0.0_dp
        N = 0.0_dp

        do m = 1, matrix_size
            S(m,m) = 1.0_dp
            S_tilde(m,m) = 1.0_dp
        end do

        N(1) = overlap(1,1)
        S(1,1) = 1/sqrt(N(1))

        do m = 2, matrix_size
            do i = 1, m-1
                do l = 1, i
                    k(i,m) = k(i,m) + S(l,i)*overlap(m,l)
                end do
            end do
            do p = 1, m-1
                do q = p, m-1
                    S_tilde(p,m) = S_tilde(p,m) - k(q,m)*S(p,q)
                end do
            end do
            do q = 1, m
                do p = 1, m
                    N(m) = N(m) + S_tilde(p,m)*S_tilde(q,m)*overlap(p,q)
                end do
            end do
            if (N(m) < 0.0_dp) return
            S(:,m) = S_tilde(:,m)/sqrt(N(m))
        end do

    end subroutine construct_gram_schmidt_transform_matrix

    subroutine print_populations_kp(kp)
    
        ! A useful test routine which will output the total walker population on both
        ! replicas, for each Krylov vector.

        type(kp_fciqmc_data), intent(in) :: kp
        integer :: ihash
        integer(n_int) :: int_sign(lenof_sign_kp)
        real(dp) :: real_sign(lenof_sign_kp), total_pop(lenof_sign_kp)
        type(ll_node), pointer :: temp_node

        int_sign = 0_n_int
        total_pop = 0.0_dp
        real_sign = 0.0_dp
        
        do ihash = 1, nhashes_kp
            temp_node => krylov_vecs_ht(ihash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    int_sign = krylov_vecs(NIfDBO+1:NIfDBO+lenof_sign_kp, temp_node%ind)
                    real_sign = transfer(int_sign, real_sign)
                    total_pop = total_pop + abs(real_sign)
                    temp_node => temp_node%next
                end do
            end if
        end do

        nullify(temp_node)

        write(6,*) "krylov_vecs populations:", total_pop

    end subroutine print_populations_kp

    subroutine print_amplitudes_kp(irepeat)

        ! A (*very* slow and memory intensive) test routine to print the current amplitudes (as stored
        ! in CurrentDets) of *all* determinants to a file. The amplitude of each replica will be printed
        ! one after the other. Since this is intended to be used with kp-fciqmc, irepeat is the number of
        ! the current repeat, but it will simply be used in naming the output file.

        ! Note that this routine will only work when using the tHashWalkerList option.

        integer, intent(in) :: irepeat
        integer, allocatable :: nI_list(:,:)
        integer :: temp(1,1), hf_ind, ndets
        integer :: i, j, ilen, counter, temp_unit, DetHash
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign(lenof_sign)
        type(ll_node), pointer :: temp_node
        type(BasisFn) :: iSym
        character(2) :: ifmt
        character(15) :: ind, filename

        ! Determine the total number of determinants.
        call gndts(nel, nbasis, BRR, nBasisMax, temp, .true., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)
        allocate(nI_list(nel, ndets))
        ! Generate the determinants and move them to nI_list.
        call gndts(nel, nbasis, BRR, nBasisMax, nI_list, .false., G1, tSpn, lms, tParity, SymRestrict, ndets, hf_ind)

        write(ind,'(i15)') irepeat
        filename = trim('amps.'//adjustl(ind))

        temp_unit = get_free_unit()
        open(temp_unit, file=trim(filename), status='replace')

        counter = 0

        do i = 1, ndets
            call getsym(nI_list(:,i), nel, G1, nBasisMax, iSym)
            ! Only carry on if the symmetry of this determinant is correct.
            if (iSym%Sym%S /= HFSym%Sym%S .or. iSym%Ms /= HFSym%Ms .or. iSym%Ml /= HFSym%Ml) cycle
            call EncodeBitDet(nI_list(:,i), ilut)
            if (.not. IsAllowedHPHF(ilut(0:NIfDBO))) cycle
            counter = counter + 1
            real_sign = 0.0_dp
            DetHash = FindWalkerHash(nI_list(:,i), nWalkerHashes)
            temp_node => HashIndex(DetHash)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if (DetBitEQ(ilut, CurrentDets(:,temp_node%ind), NIfDBO)) then
                        int_sign = CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,temp_node%ind)
                        real_sign = transfer(int_sign, real_sign)
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            ilen = ceiling(log10(real(abs(counter)+1)))
            ! ifmt will hold the correct integer length so that there will be no spaces printed out.
            ! Note that this assumes that ilen < 10, which is very reasonable!
            write(ifmt,'(a1,i1)') "i", ilen
            do j = 1, lenof_sign
                ! This assumes that lenof_sign < 10. Probably will always be 1 or 2.
                write(temp_unit,'(a1,'//ifmt//',a1,i1,a1,1x,es19.12)') "(", counter,",", j, ")", real_sign(j)
            end do
        end do

        close(temp_unit)

        deallocate(nI_list)

    end subroutine print_amplitudes_kp

end module kp_fciqmc_procs
