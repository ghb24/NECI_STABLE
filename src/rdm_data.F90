module rdm_data

    ! Module containing global data and derived types used for RDM calculation.

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

    type rdm_estimates_t
        ! How many RDMs are being sampled.
        integer :: nrdms

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

        ! Arrays used to hold estimates from the RDM over the *previous
        ! sampling block only*.
        real(dp), allocatable :: trace_inst(:)
        real(dp), allocatable :: norm_inst(:)
        real(dp), allocatable :: energy_1_num_inst(:)
        real(dp), allocatable :: energy_2_num_inst(:)
        real(dp), allocatable :: energy_num_inst(:)
        real(dp), allocatable :: spin_num_inst(:)

        ! Hermiticity errors, i.e. \Gamma_{ij,kl} - \Gamma_{kl,ij}^*.
        ! The max_* array holds the maximum such error.
        ! The sum_* array holds the sum of all such errors.
        real(dp), allocatable :: max_error_herm(:)
        real(dp), allocatable :: sum_error_herm(:)
    end type rdm_estimates_t

    ! Global data.

    ! The primary global RDM objects.
    ! Arrays of objects, one for each 1-RDM being sampled.
    type(one_rdm_t), allocatable :: one_rdms(:) ! nrdms
    ! Object to hold spawnings to the 2-RDMs.
    type(rdm_spawn_t) :: two_rdm_spawn
    ! Object to hold the main RDM itself, over the *entire* period of RDM
    ! sampling (note that this is not reset each sampling block).
    type(rdm_list_t) :: two_rdm_main
    ! Objects to hold the received RDM object, after communication of the
    ! spawned RDM list. This is then added into two_rdm_main.
    type(rdm_list_t) :: two_rdm_recv
    type(rdm_list_t) :: two_rdm_recv_2
    ! Object to hold RDM estimates.
    type(rdm_estimates_t) :: rdm_estimates

    ! The number of rdms being calculated in this simulation.
    integer :: nrdms = 0

    ! If true, then 2-RDM quantities will be output to a RDMEstimates file.
    logical :: print_2rdm_est

    ! Variable used in RDM calculations to specify that an open shell system
    ! is being studied.
    logical :: tOpenShell

    ! Logical for natural orbital caluculation, to speficy whether orbitals
    ! have been rotated yet.
    logical :: tRotatedNOs = .false.

    ! Timers.
    type(timer), save :: nElRDM_Time, FinaliseRDMs_time, RDMEnergy_time

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
