module rdm_data

    ! Module containing global data and derived types used for RDM calculation.

    use constants
    use global_utilities, only: timer

    implicit none

    integer, allocatable :: Sing_InitExcSlots(:), Sing_ExcList(:)
    integer, allocatable :: Doub_InitExcSlots(:), Doub_ExcList(:)
    integer(kind=n_int), allocatable :: Sing_ExcDjs(:,:), Sing_ExcDjs2(:,:)
    integer(kind=n_int), allocatable :: Doub_ExcDjs(:,:), Doub_ExcDjs2(:,:)

    integer :: Sing_ExcDjsTag, Sing_ExcDjs2Tag
    integer :: Doub_ExcDjsTag, Doub_ExcDjs2Tag, UMATTempTag
    integer :: Energies_unit
    integer :: NoSymLabelCounts, Rho_iiTag

    real(dp), allocatable, target :: aaaa_RDM_inst(:,:), abab_RDM_inst(:,:), abba_RDM_inst(:,:)
    real(dp), allocatable, target :: aaaa_RDM_full(:,:), abab_RDM_full(:,:), abba_RDM_full(:,:)
    real(dp), allocatable, target :: bbbb_RDM_inst(:,:), baba_RDM_inst(:,:), baab_RDM_inst(:,:)
    real(dp), allocatable, target :: bbbb_RDM_full(:,:), baba_RDM_full(:,:), baab_RDM_full(:,:)

    integer :: aaaa_RDM_instTag, abab_RDM_instTag, abba_RDM_instTag
    integer :: aaaa_RDM_fullTag, abab_RDM_fullTag, abba_RDM_fullTag
    integer :: bbbb_RDM_instTag, baba_RDM_instTag, baab_RDM_instTag
    integer :: bbbb_RDM_fullTag, baba_RDM_fullTag, baab_RDM_fullTag

    real(dp), pointer :: aaaa_RDM(:,:) => null()
    real(dp), pointer :: abab_RDM(:,:) => null()
    real(dp), pointer :: abba_RDM(:,:) => null()
    real(dp), pointer :: bbbb_RDM(:,:) => null()
    real(dp), pointer :: baba_RDM(:,:) => null()
    real(dp), pointer :: baab_RDM(:,:) => null()

    real(dp), allocatable :: AllNodes_aaaa_RDM(:,:), AllNodes_abab_RDM(:,:), AllNodes_abba_RDM(:,:) 
    real(dp), allocatable :: AllNodes_bbbb_RDM(:,:), AllNodes_baba_RDM(:,:), AllNodes_baab_RDM(:,:)  

    real(dp), allocatable :: UMATTemp(:,:)
    real(dp), allocatable :: Rho_ii(:)
    real(dp), allocatable :: Lagrangian(:,:)

    logical :: tCalc_RDMEnergy

    real(dp) :: OneEl_Gap, TwoEl_Gap, Normalisation, Trace_2RDM_Inst, Trace_2RDM, Trace_1RDM, norm
    type(timer), save :: nElRDM_Time, FinaliseRDM_time, RDMEnergy_time

    logical :: tRotatedNOs = .false.

    logical :: tOpenShell

end module rdm_data
