module rdm_data

    ! Module containing global data and derived types used for RDM calculation.

    use constants
    use FciMCData, only: ll_node
    use global_utilities, only: timer

    implicit none

    ! The number of rdms being calculated in this simulation.
    integer :: nrdms = 0

    ! If true, then the RDM energy will be calculated and output.
    ! WARNING: despite its name, this variable is also required to be true in
    ! WARNING: order to output 2-RDMs.
    logical :: tCalc_RDMEnergy

    ! Variable used in RDM calculations to specify that an open shell system
    ! is being studied.
    logical :: tOpenShell

    logical :: tRotatedNOs = .false.

    ! Unit of the separate file to which RDM estimates (such as energy and
    ! spin^2) are output.
    integer :: rdm_estimates_unit

    ! Arrays used as storage when summing RDMs over all processors.
    real(dp), allocatable :: AllNodes_RDM_small(:,:), AllNodes_RDM_large(:,:)

    ! Arrays for when filling arrays explicitly. See rdm_explicit for the
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

    ! Timers.
    type(timer), save :: nElRDM_Time, FinaliseRDMs_time, RDMEnergy_time

    ! Derived type to hold data for each RDM - other global data will be
    ! removed after purification work.

    type rdm_t
        ! Arrays to hold instantaneous estimates of 2-RDMs.
        ! The following three arrays are always used.
        real(dp), pointer :: aaaa_inst(:,:) => null(), abab_inst(:,:) => null(), abba_inst(:,:) => null()
        ! And the following three arrays are only used for open shell calculations
        ! (and possibly UHF systems in the future?).
        real(dp), pointer :: bbbb_inst(:,:) => null(), baba_inst(:,:) => null(), baab_inst(:,:) => null()

        ! Arrays to hold estimates of 2-RDMs, summed over the whole RDM calculation.
        real(dp), pointer :: aaaa_full(:,:) => null(), abab_full(:,:) => null(), abba_full(:,:) => null()
        real(dp), pointer :: bbbb_full(:,:) => null(), baba_full(:,:) => null(), baab_full(:,:) => null()

        ! Tags for the memory manager for the above RDM arrays.
        integer :: aaaa_instTag, abab_instTag, abba_instTag
        integer :: bbbb_instTag, baba_instTag, baab_instTag
        integer :: aaaa_fullTag, abab_fullTag, abba_fullTag
        integer :: bbbb_fullTag, baba_fullTag, baab_fullTag

        ! Pointers which can point to *_inst or *_full, depending on what the user
        ! has asked for.
        real(dp), pointer :: aaaa(:,:) => null()
        real(dp), pointer :: abab(:,:) => null()
        real(dp), pointer :: abba(:,:) => null()
        real(dp), pointer :: bbbb(:,:) => null()
        real(dp), pointer :: baba(:,:) => null()
        real(dp), pointer :: baab(:,:) => null()

        ! The 1-RDM.
        real(dp), allocatable :: matrix(:,:)

        ! Eigenvalues of the 1-RDM.
        real(dp), allocatable :: Evalues(:)
        ! Arrays to hold the diagonal of the 1-RDM, and the Lagrangian.
        real(dp), allocatable :: Rho_ii(:)
        real(dp), allocatable :: Lagrangian(:,:)

        integer :: Rho_iiTag, matrix_tag, EvaluesTag

        ! TODO: Comment.
        integer, allocatable :: sym_list_no(:)
        integer, allocatable :: sym_list_inv_no(:)

    end type rdm_t

    type one_rdm_t
        ! The 1-RDM.
        real(dp), allocatable :: matrix(:,:)

        ! Eigenvalues of the 1-RDM.
        real(dp), allocatable :: Evalues(:)
        ! Arrays to hold the diagonal of the 1-RDM, and the Lagrangian.
        real(dp), allocatable :: Rho_ii(:)
        real(dp), allocatable :: Lagrangian(:,:)

        integer :: Rho_iiTag, matrix_tag, EvaluesTag

        ! TODO: Comment.
        integer, allocatable :: sym_list_no(:)
        integer, allocatable :: sym_list_inv_no(:)

    end type one_rdm_t

    type rdm_estimates_old_t
        real(dp) :: RDMEnergy, RDMEnergy1, RDMEnergy2, RDMEnergy_Inst
        real(dp) :: Norm_2RDM, Norm_2RDM_Inst, Trace_2RDM, Trace_2RDM_normalised
        real(dp) :: spin_est
    end type rdm_estimates_old_t

    type rdm_estimates_t
        integer :: nrdms

        real(dp), allocatable :: trace(:)
        real(dp), allocatable :: norm(:)

        real(dp), allocatable :: energy_1_num(:)
        real(dp), allocatable :: energy_2_num(:)
        real(dp), allocatable :: energy_tot_num(:)

        real(dp), allocatable :: spin_num(:)

        real(dp), allocatable :: max_error_herm(:)
        real(dp), allocatable :: sum_error_herm(:)
    end type rdm_estimates_t

    ! Array of type rdm_t, for holding multiple different RDM instances.
    type(rdm_t), allocatable :: rdms(:)
    type(rdm_estimates_old_t), allocatable :: rdm_estimates_old(:)

    ! Data for parallel RDM implementation.

    ! This data type is used for storing RDMs as 1D lists. A specific
    ! ordering of RDM elements in this list is not required - rather, a
    ! hash table can be used to access elements.
    type rdm_list_t
        ! The number of integers available to store signs, for each
        ! RDM element in the elements array.
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
        ! with the RDM list. Note that the nelements component won't be used
        ! or relevant, since the list won't be contiguous.
        type(rdm_list_t) :: rdm_send

        ! free_slots(i) holds the next available spawning slot in
        ! rdm_send%elements for processor i.
        integer, allocatable :: free_slots(:)
        ! init_free_slots(i) holds the index in rdm_send%elements where the
        ! very first RDM element to be sent to process i will be added.
        integer, allocatable :: init_free_slots(:)
    end type rdm_spawn_t

    type(one_rdm_t), allocatable :: one_rdms(:)
    type(rdm_spawn_t) :: two_rdm_spawn
    type(rdm_list_t) :: rdm_main
    type(rdm_list_t) :: two_rdm_recv
    type(rdm_estimates_t) :: rdm_estimates

end module rdm_data
