module rdm_data

    ! Module containing global data and derived types used for RDM calculation.

    use constants
    use global_utilities, only: timer

    implicit none

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

    real(dp), allocatable :: Lagrangian(:,:)

    real(dp) :: Trace_2RDM, Trace_2RDM_Inst

    ! Timers.
    type(timer), save :: nElRDM_Time, FinaliseRDM_time, RDMEnergy_time

    ! Derived type to hold data for each RDM - other global data will be
    ! removed after purification work.

    type rdm_t

        ! Arrays to hold instantaneous estimates of 2-RDMs.
        ! The following three arrays are always used.
        real(dp), pointer :: aaaa_inst(:,:), abab_inst(:,:), abba_inst(:,:)
        ! And the following three arrays are only used for open shell calculations
        ! (and possibly UHF systems in the future?).
        real(dp), pointer :: bbbb_inst(:,:), baba_inst(:,:), baab_inst(:,:)

        ! Arrays to hold estimates of 2-RDMs, summed over the whole RDM calculation.
        real(dp), pointer :: aaaa_full(:,:), abab_full(:,:), abba_full(:,:)
        real(dp), pointer :: bbbb_full(:,:), baba_full(:,:), baab_full(:,:)

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

        ! Arrays to hold the diagonal of the 1-RDM, and the Lagrangian.
        real(dp), allocatable :: Rho_ii(:)
        real(dp), allocatable :: Lagrangian(:,:)
        integer :: Rho_iiTag

        ! RDM traces.
        real(dp) :: Trace_2RDM, Trace_2RDM_Inst

    end type rdm_t

    ! Array of type rdm_t, for holding multiple different RDM instances.
    type(rdm_t), allocatable :: rdms(:)

end module rdm_data
