module kp_fciqmc_data_mod

    use CalcData, only: subspace_in
    use constants
    use FciMCData, only: ll_node, perturbation

    implicit none
    save

    type kp_fciqmc_data
        ! The number of different initial walker configurations to start
        ! calculations from.
        integer :: nconfigs
        ! This is only used in the excited-state KP algorithm.
        ! The number of times that data about the subspace Hamiltonian and
        ! overlap matrices are output during a single run.
        integer :: nreports
        ! The number of simulations to perform for each initial walker
        ! configuration.
        integer :: nrepeats
        ! The number of different Krylov vectors to sample (the number of
        ! vectors which form the Krylov subspace at the end of a calculation).
        integer :: nvecs
        ! The number of iterations to perform *between each Krylov vector being
        ! sampled*. niters(i) holds the number to be performed between the
        ! i-th and (i+1)th vectors. The final element is not set by the user,
        ! it is always set to 1 (because we don't want to go any further
        ! beyond the final Krylov vector).
        integer, allocatable :: niters(:)
    end type

    ! If true then we generate the initial Krylov subspace by using estimates
    ! of the lowest excited states, and project all of these trial states
    ! every iteration.
    logical :: tExcitedStateKP

    ! If true then calculate and output the spin squared operator in the
    ! Krylov basis.
    logical :: tCalcSpin

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
    ! The diagonal Hamiltonian elements of the vectors in krylov_vecs, in the
    ! same order.
    real(dp), allocatable :: krylov_helems(:)
    ! The hash is table used to access determinant data in krylov_vecs.
    type(ll_node), pointer :: krylov_vecs_ht(:) 

    ! These arrays are used if tExcitedStateKP = .false. in calc_projected_hamil
    ! in semi-stochasti calculations. They are used to store the deterministic
    ! vectors before and after the deterministic projection occurs. They have
    ! the same purpose as partial_determ_vecs and full_determ_vecs, except they
    ! are allocated to hold more vectors (all Krylov vectors), as is necessary
    ! if tExcitedStateKP = .false.
    HElement_t(dp), allocatable, dimension(:,:) :: partial_determ_vecs_kp
    HElement_t(dp), allocatable, dimension(:,:) :: full_determ_vecs_kp

    ! The number of elements in krylov_vecs which are used to store amplitudes.
    integer(int64) :: nkrylov_amp_elems_tot
    ! The nunber of amplitude elements in krylov_vecs which are non-zero.
    integer(int64) :: nkrylov_amp_elems_used

    ! The values of these variables for the current initial walker configuration.
    ! These are held in these variables so that, when we restart from this
    ! configuration, we can instantly reset the corresponding values.
    integer(int64) :: TotWalkersInit, AllTotWalkersInit
    real(dp), allocatable :: TotPartsInit(:)
    real(dp), allocatable :: AllTotPartsInit(:)

    ! If true then don't use a fixed number of iterations between each Krylov
    ! vector is taken, but vary the number.
    logical :: vary_niters

    ! The total sign length for all Krylov vectors together.
    integer :: lenof_all_signs

    ! If true then, in the stochastic determination of the projected
    ! Hamiltonian, perform spawning from more highly weighted determinants
    ! exactly.
    logical :: tExactHamilSpawning
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
    ! When estimating the projected Hamiltonian, this parameter will be
    ! used to decide how many determinants have *exact* spawning performed
    ! for them. Specifically, if the total population on a determinant
    ! (across all replicas) multiplied by kp_hamil_exact_frac is greater
    ! than or equal to the number of determinants connected to the
    ! Hartree-Fock, then exact spawning is performed.
    real(dp) :: kp_hamil_exact_frac
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

    ! This variable is used to store the state of tSinglePartPhase on input so
    ! that we can reset it to this value on each repeat.
    logical, allocatable :: tSinglePartPhaseKPInit(:)

    ! If true then calculate the overlap of the final Hamiltonian eigenstates
    ! with a vector which is calculated by applying a perturbation operator
    ! (overlap_pert) to the popsfile wave function.
    logical :: tOverlapPert
    type(perturbation), allocatable :: overlap_pert(:)
    ! The result of overlap_pert applied to the wave function read in from
    ! a popsfile.
    integer(n_int), allocatable :: perturbed_ground(:,:)
    ! The overlaps of perturbed_ground with the various Krylov vectors.
    real(dp), allocatable :: pert_overlaps(:)
    ! The sum of pert_overlaps across all processes, calculated after all
    ! iterations in a repeat have been performed and only stored on the
    ! root process.
    real(dp), allocatable :: kp_all_pert_overlaps(:)

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
    ! projected matrices into an orthonormal basis, in the Lowdin approach.
    real(dp), allocatable :: kp_inter_matrix(:,:)
    ! The final eigenvectors, but in the basis of Krylov vectors.
    ! This is used in the Lowdin approach: diagonalising kp_final_hamil
    ! will give the final eigenvectors in the orthonormal basis, then
    ! we do the inverse of the transformation procedure to get back to
    ! the Krylov basis.
    real(dp), allocatable :: kp_eigenvecs_krylov(:,:)

    ! If true then use an estimate of an excited state (or linear
    ! combination of excited states) to form the initial wave function.
    logical :: tExcitedInitState
    ! When using the excited-init-state option, this array is read in
    ! from the input file and specifies which excited states to take
    ! a weighted average of to form the initial state.
    integer, allocatable :: kpfciqmc_ex_labels(:)
    ! How much weight do we give to each trial state in the initial
    ! Krylov vector?
    real(dp), allocatable :: kpfciqmc_ex_weights(:)

    ! Type for the trial wave function space for excited-state calculations.
    type(subspace_in) :: kp_trial_space_in

    ! Arrays used to access the signs of different Krylov vector signs in the
    ! Krylov vector arrays. If tPairedReplicas is .false. then these arrays
    ! are identical, otherwise the first array gives access to the first set
    ! of vectors and the second array to the second set.
    integer, allocatable :: kp_ind_1(:), kp_ind_2(:)

    ! If true then perform an orthogonalisation step at the end of each
    ! iteration. This only applies to the CFQMC approach.
    logical :: tOrthogKPReplicas
    ! After which iteration should we start performing the orthogonalisation
    ! step?
    integer :: orthog_kp_iter
    
end module
