module rdm_data

    ! Module containing global data and derived types used for RDM calculation.

    ! The following is a description of some of the details of how RDMs are
    ! stored, including both technical aspects, and details on which elements
    ! are not stored due to symmetry. There is also a note on how variables
    ! are named, particularly orbital labels. These things are important for
    ! understanding various routines that perform operations on RDMs,
    ! particularly in rdm_finalising.

    ! Technical aspects of data structures
    ! ====================================

    ! 2-RDMs are stored in rdm_list_t objects. When new contributions to 2-RDMs
    ! are generated (through stochastic FCIQMC spawnings), they are added to
    ! a rdm_spawn_t object. two_rdm_spawn is the main global rdm_spawn_t object
    ! used for this purpose. At certain points (currently every iteration),
    ! these 2-RDM element 'spawnings' will be communicated to their correct
    ! processes. Another rdm_list_t object is used for this purpose (mainly
    ! two_rdm_recv). The resulting received RDM is then added into the existing
    ! main 2-RDM object through add_rdm_1_to_rdm_2, which simply performs
    ! numerical addition of these objects. In this way, the RDM (two_rdm_main)
    ! is accumulated over the entire period of RDM calculation, and the RDM
    ! is distributed.

    ! Currently, RDM elements are distributed across processes using the row
    ! label of the element. If there are n_proc processes and n_row rows, then
    ! each process holds n_row/n_proc rows, with the first process holding all
    ! the first n_row/n_proc rows, the second process holding the second such
    ! set, and so.

    ! Each rdm_list_t object holds not only the array of RDM elements, but also
    ! a hash table which refers to this list. This hash table is made up of
    ! linked lists, one for each hash value available (see hash.F90 for more
    ! details). With this, when adding a 2-RDM element, it can be checked
    ! relatively quickly if the element is already in the list, and the
    ! new element can be added directly into that existing position. This avoids
    ! having to perform an explicit sorting step, and also avoids extra memory
    ! being needed for repeated elements. Not only does the main RDM array have
    ! a hash_table for this purpose, but also the rdm_spawn_t objects (since
    ! they themselves hold an rdm_list_t object).

    ! Note on variable names
    ! ======================

    ! Throughout the RDM routines, we always try to use the variables i, j, k
    ! and l to refer to spin orbital labels, and p, q, r and s to refer to
    ! spatial orbitals labels. Please try and keep this, to avoid confusion!
    ! In some instances, particularly in rdm_data_utils routines, a routine
    ! can act on both spinned and spin-free 2-RDMs (for example, the routine
    ! add_rdm_1_to_rdm_2 will work on either types of RDM, each of which can
    ! be stored as a rdm_list_t object), in which case we usually use
    ! (i,j,k,l), but these could refer to either spin or spatial labels
    ! depending on the RDM.

    ! Similarly, when looping over all elements in an RDM array, we try to
    ! use ielem as the looping variable name. When looping over the various
    ! RDMs being sampled, we use irdm as the loop variable. Please try and
    ! avoid using i as a loop variable when there might be confusion with
    ! orbital labels, within RDM modules.

    ! Important points regarding which RDM elements are stored
    ! ========================================================

    ! For a 2-RDM element \Gamma_{ij,kl}, it is only stored if i<j and k<l.
    ! An element with i>j or k>l (or both) can be made into the above form by
    ! swapping i & j or k & l (or both), which introduces a minus sign. This
    ! is true even for RDMs formed from non-real wave functions.

    ! The above has a consequence which complicates some RDM processing
    ! routines: if i and j have the same spatial parts, and k and l *also*
    ! have the same spatial part, then we won't get both standard and
    ! spin-flipped versions of the element present in the RDM array - only one
    ! of them. Suppose 1 labels the first spatial orbital with beta spin, and
    ! 2 labels the first spatial orbital with alpha spin. Now consider the RDM
    ! element \Gamma_{12,12}. This can be generated and stored because i<j and
    ! k<l. But the spin-flipped version is \Gamma_{21,21}, which will never
    ! get generated or stored. This is not a problem because the latter is
    ! always equal to the former. However, this also leads to some edge cases
    ! which must be carefully considered when summing over spin. This is
    ! because all other types of elements in the spin-free 2-RDM end up getting
    ! counted twice when summing over all stored elements, because their
    ! spin-flipped partners *are* also present. See create_spinfree_2rdm for
    ! an example... Be careful!

    ! On the other hand, the symmetry \Gamma_{ij,kl} = \Gamma_{kl,ij}^* is
    ! *not* in use - 2-RDM elements are accumulated and stored both ways around.

    use constants, only: dp, n_int, int_rdm
    use FciMCData, only: ll_node
    use global_utilities, only: timer

    implicit none

    ! Data structures for RDM code.

    type one_rdm_t
        ! The 1-RDM object itself.
        real(dp), allocatable :: matrix(:,:)

        ! Eigenvalues of the 1-RDM.
        real(dp), allocatable :: evalues(:)
        ! Arrays to hold the diagonal of the 1-RDM, and the Lagrangian.
        real(dp), allocatable :: rho_ii(:)
        real(dp), allocatable :: lagrangian(:,:)

        integer :: rho_ii_tag, matrix_tag, evalues_tag

        ! In the 1-RDM matrix array, elements are not stored in the same order
        ! as the orbitals in the rest of the code. Instead, orbitals are
        ! reordered so that orbitals are first sorted by their symmetry label.
        ! sym_list_no(i) holds the position of the corresponding row/column
        ! in the RDM for orbital i. i.e. orbital i is mapped to sym_list_no(i).
        integer, allocatable :: sym_list_no(:)
        ! The inverse of sym_list_no, i.e. what orbital does row or column i
        ! in the 1-RDM correspond to?
        integer, allocatable :: sym_list_inv_no(:)
    end type one_rdm_t

    ! This data type is used for storing RDMs as 1D lists. A specific
    ! ordering of RDM elements in this list is not required - rather, a
    ! hash table can be used to access elements.
    type rdm_list_t
        ! The number of integers available to store signs, for each
        ! RDM entry in the elements array.
        integer :: sign_length = 0
        ! Array which holds the RDM elements.
        integer(int_rdm), allocatable :: elements(:,:)
        ! Hash table to the rdm array.
        type(ll_node), pointer :: hash_table(:)
        ! The allocated size of the elements array.
        integer :: max_nelements = 0
        ! The number of RDM elements currently entered into the elements array.
        integer :: nelements = 0
        ! Maximum number of unique hashes available in hash_table (not the
        ! number of currently unused ones, but the total number, i.e. the
        ! length of the hash_table array).
        integer :: nhashes = 0
    end type rdm_list_t

    ! This data type is used for accumulating contributions to an RDM, spawned
    ! on this processor. There are a collection of routines which then work
    ! specifically with this data structure, to perform communication.
    type rdm_spawn_t
        ! The number of rows in the RDM.
        integer :: nrows = 0

        ! This object holds the spawning array for the RDM elements before
        ! they are communicated, as well as relevant metadata which is stored
        ! with the RDM list (see rdm_list_t). Note that the nelements component
        ! won't be used or relevant, since the list won't be contiguous.
        type(rdm_list_t) :: rdm_send

        ! free_slots(i) holds the next available spawning slot in
        ! rdm_send%elements for processor i.
        integer, allocatable :: free_slots(:)
        ! init_free_slots(i) holds the index in rdm_send%elements where the
        ! very first RDM element to be sent to process i will be added.
        integer, allocatable :: init_free_slots(:)
    end type rdm_spawn_t

    type en_pert_t
        ! The number of integers available to store signs, for each
        ! RDM entry in the elements array.
        integer :: sign_length = 0
        ! Array which holds the RDM elements.
        integer(n_int), allocatable :: dets(:,:)
        ! Hash table to the rdm array.
        type(ll_node), pointer :: hash_table(:)
        ! The allocated size of the elements array.
        integer :: max_ndets = 0
        ! The number of determinants contributing to the perturbation
        ! currently stored in the dets array.
        integer :: ndets = 0
        ! The number of determinants contributing to the perturbation
        ! currently stored in the dets array, across all processors.
        integer :: ndets_all = 0
        ! Maximum number of unique hashes available in hash_table (not the
        ! number of currently unused ones, but the total number, i.e. the
        ! length of the hash_table array).
        integer :: nhashes = 0
    end type en_pert_t

    type rdm_estimates_t
        ! How many RDMs are being sampled.
        integer :: nrdms
        integer :: nrdms_standard
        integer :: nrdms_transition

        ! Unit of the separate file to which RDM estimates (such as energy and
        ! spin^2) are output.
        integer :: write_unit

        ! The following arrays have length nrdms - one estimate is held for
        ! each RDM.

        ! Arrays used to hold estimates from the *total* RDM (i.e. the array
        ! averaged over the whole RDM sampling period).
        real(dp), allocatable :: trace(:)
        real(dp), allocatable :: norm(:)
        real(dp), allocatable :: energy_1_num(:)
        real(dp), allocatable :: energy_2_num(:)
        real(dp), allocatable :: energy_num(:)
        real(dp), allocatable :: spin_num(:)
        real(dp), allocatable :: property(:,:)
        real(dp), allocatable :: energy_pert(:)
        real(dp), allocatable :: energy_pert_hf(:)

        ! Arrays used to hold estimates from the RDM over the *previous
        ! sampling block only*.
        real(dp), allocatable :: trace_inst(:)
        real(dp), allocatable :: norm_inst(:)
        real(dp), allocatable :: energy_1_num_inst(:)
        real(dp), allocatable :: energy_2_num_inst(:)
        real(dp), allocatable :: energy_num_inst(:)
        real(dp), allocatable :: spin_num_inst(:)
        real(dp), allocatable :: property_inst(:,:)
        real(dp), allocatable :: energy_pert_inst(:)
        real(dp), allocatable :: energy_pert_hf_inst(:)

        ! Hermiticity errors, i.e. \Gamma_{ij,kl} - \Gamma_{kl,ij}^*.
        ! The max_* array holds the maximum such error.
        ! The sum_* array holds the sum of all such errors.
        real(dp), allocatable :: max_error_herm(:)
        real(dp), allocatable :: sum_error_herm(:)

    end type rdm_estimates_t

    ! Data type used to define how many RDMs are being sampled and which states
    ! and FCIQMC simulations contribute to each of these RDMs. It also holds
    ! arrays which can be used to efficiently find all simulations that
    ! contribute to an RDM with a given other simulation, and a few other
    ! similar arrays.
    type rdm_definitions_t
        ! The total number of RDMs being calculated.
        ! Equal to nrdms_standard + nrdms_transition.
        integer :: nrdms = 0
        ! The number of 'standard' RDMs (i.e. non-transition RDMs) being
        ! calculated.
        integer :: nrdms_standard = 0
        ! The number of transition RDMs being calculated.
        integer :: nrdms_transition = 0

        ! state_labels(:,j) will store the labels of the *actual* wave functions
        ! (i.e., usually which excited state it is) contributing to the j'th RDM.
        integer, allocatable :: state_labels(:,:) ! (2, nrdms)
        ! sim_labels(:,j) will store the labels of the *FCIQMC* simulations
        ! (i.e. the 'replica' labels) which will be used to sample the j'th RDM
        ! being calculated.
        integer, allocatable :: sim_labels(:,:) ! (2, nrdms)

        ! For transition RDMs, with 2 replicas for each state, there will be 2
        ! copies of each transition RDM. This array simply specifies which of
        ! the 2 each RDM is - the first or second 'repeat'.
        integer, allocatable :: repeat_label(:) ! (nrdms)
        ! nrdms_per_sim(j) holds the number of different RDMS to which
        ! the FCIQMC simulation with label j contributes to.
        integer, allocatable :: nrdms_per_sim(:) ! (lenof_sign)
        ! sim_pairs(:,j) holds the list of the FCIQMC simulations labels
        ! which are paired with simulation j in contributing to RDMs.
        ! Elements which are not needed (due to a simulation not
        ! contributing to all RDMs) are set to 0.
        integer, allocatable :: sim_pairs(:,:) ! (nrdms, lenof_sign)
        ! rdm_labels(:,j) holds the list of RDM labels which simulation j
        ! contributes to. Elements which are not needed (due a simulation not
        ! contributing to all RDMs) are set to 0.
        integer, allocatable :: rdm_labels(:,:) ! (nrdms, lenof_sign)
        ! prefix for the names of the files to output to
        character(255) :: output_file_prefix
    end type rdm_definitions_t

    ! Global data.

    ! Factors which can be set by the user at input to modify the size of RDM
    ! arrays relative to their default sizes (as specified in the init_rdms
    ! routine).
    real(dp) :: rdm_main_size_fac = 1.0_dp
    real(dp) :: rdm_spawn_size_fac = 1.0_dp
    real(dp) :: rdm_recv_size_fac = 1.0_dp

    ! The primary global RDM objects.
    ! Arrays of objects, one for each 1-RDM being sampled.
    type(one_rdm_t), allocatable :: one_rdms(:) ! nrdms
    ! same for the initiator-only 1-RDMs
    type(one_rdm_t), allocatable :: inits_one_rdms(:) ! nrdms
    ! Object to hold spawnings to the 2-RDMs.
    type(rdm_spawn_t) :: two_rdm_spawn
    ! spawnings to the initiator space 2-RDMs
    type(rdm_spawn_t) :: two_rdm_inits_spawn
    ! Object to hold the main RDM itself, over the *entire* period of RDM
    ! sampling (note that this is not reset each sampling block).
    type(rdm_list_t) :: two_rdm_main
    ! Initiator-only RDMs
    type(rdm_list_t) :: two_rdm_inits
    ! Objects to hold the received RDM object, after communication of the
    ! spawned RDM list. This is then added into two_rdm_main.
    type(rdm_list_t) :: two_rdm_recv
    type(rdm_list_t) :: two_rdm_recv_2
    ! Object to hold RDM estimates.
    type(rdm_estimates_t) :: rdm_estimates
    type(rdm_estimates_t) :: inits_estimates
    ! Object which defines the states and FCIQMC simulations contributing
    ! to the various RDMs in the above arrays.
    type(rdm_definitions_t) :: rdm_definitions
    type(rdm_definitions_t) :: rdm_inits_defs

    logical :: tSetupInitsEst = .false.

    ! Object to hold information about the Epstein-Nesbet perturbation
    ! contributions.
    type(en_pert_t) :: en_pert_main

    ! The number of transition RDMs that the user asks for at input. This
    ! might be equal to half the number of tRDMs actually calculated, since
    ! when using the replica trick we calculate 2 of each tRDM.
    integer :: nrdms_transition_input
    ! This is the same as for state_labels in rdm_definition_t, but *only*
    ! deals with transition RDMs specifically. This array is used to hold the
    ! states specified by the user at input.
    integer, allocatable :: states_for_transition_rdm(:,:) ! (2, nrdms_transition_input)

    ! If true, then 2-RDM quantities will be output to a RDMEstimates file.
    logical :: print_2rdm_est

    ! Variable used in RDM calculations to specify that an open shell system
    ! is being studied.
    logical :: tOpenShell, tOpenSpatialOrbs

    ! Logical for natural orbital caluculation, to speficy whether orbitals
    ! have been rotated yet.
    logical :: tRotatedNOs = .false.

    ! Timers.
    type(timer), save :: nElRDM_Time, FinaliseRDMs_time, RDMEnergy_time

    ! ---- Data for using adaptive shift mode ------------------------!
    
    ! when using adaptive shift, the RDMs require a correction, namely
    ! the reference contribution

    real(dp) :: rdmCorrectionFactor, InstRDMCorrectionFactor, ThisRDMIter
    logical :: tApplyLC = .true.

    ! ---- Data for the explicit RDM code -----------------------------

    ! Arrays for when filling RDMs explicitly. See rdm_explicit for the
    ! relevant routines in that case.
    integer, allocatable :: Sing_InitExcSlots(:), Sing_ExcList(:)
    integer, allocatable :: Doub_InitExcSlots(:), Doub_ExcList(:)
    integer(n_int), allocatable :: Sing_ExcDjs(:,:), Sing_ExcDjs2(:,:)
    integer(n_int), allocatable :: Doub_ExcDjs(:,:), Doub_ExcDjs2(:,:)

    ! Tags for explicitly-filled RDMs.
    integer :: Sing_ExcDjsTag, Sing_ExcDjs2Tag
    integer :: Doub_ExcDjsTag, Doub_ExcDjs2Tag

    ! Normalisation factor used in explicit RDM code.
    real(dp) :: ExcNorm
    ! Variables related to the space in explicit RDM arrays above.
    real(dp) :: OneEl_Gap, TwoEl_Gap

    ! ---- End of data for the explicit RDM code ----------------------

end module rdm_data
