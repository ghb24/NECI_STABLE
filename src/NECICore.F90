#include "macros.h"

! Because of optional arguments, an explicit interface for NECICore
! is required if you want to call it.
! You get it, with an ``#include "NECICore.h"``.
! One can't use modules because of required compiler independence.

Subroutine NECICore(iCacheFlag, tCPMD, tVASP, tMolpro_local, call_as_lib, &
                    int_name, filename_in, MemSize)
    != NECICore is the main outline of the NECI Program.
    != It provides a route for calling NECI when accessed as a library, rather
    != than as a standalone program.
    != In:
    !=    iCacheFlag: controls the behaviour of the 4-index integral cache.
    !=                Currently relevant only for CPMD and VASP calculations.
    !=                iCacheFlag=0: initialise and destroy the cache.
    !=                iCacheFlag=1: initialise but don't destroy the cache.
    !=                iCacheFlag=2: reuse and destroy the cache.
    !=                iCacheFlag=3: reuse and keep the cache.
    !=    tCPMD: True if doing a CPMD-based calculation.
    !=    tVASP: True if doing a VASP-based calculation.
    !=    call_as_lib: True if called as subroutine from external code.
    !=    int_name: is the name of the integral file to read in if necessary
    !=    filename: is the name of the input file to read in if necessary
    !=    MemSize: Memory limit in MB

    use ReadInput_neci, only: ReadInputMain
    use SystemData, only: tMolpro, tMolproMimic, MolproID, called_as_lib
    use MemoryManager

    ! main-level modules.
    use Calc, only: CalcDoCalc
    use CalcData, only: tUseProcsAsNodes
    use kp_fciqmc_procs, only: kp_fciqmc_data
    use Parallel_neci, only: MPINodes, iProcIndex, &
                             neci_MPIInit_called, neci_MPINodes_called
    use read_fci, only: FCIDUMP_name

    use ParallelHelper
    use UMatCache, only: UMat2D, tagUMat2D

    ! Utility modules.
    use global_utilities
    use constants
    use util_mod, only: get_free_unit

    USE MolproPlugin

    Implicit none
    integer, intent(in), optional :: iCacheFlag
    logical, intent(in), optional :: tCPMD, tVASP, tMolpro_local, call_as_lib
    character(*), intent(in), optional :: filename_in, int_name
    integer(int64), intent(in), optional :: MemSize
    type(timer), save :: proc_timer
    integer :: ios, iunit, iunit2, i, j, isfreeunit, iCacheFlag_
    character(*), parameter :: this_routine = 'NECICore'
    character(:), allocatable :: Filename
    logical :: toverride_input, tFCIDUMP_exist, tCPMD_, tVASP_
    type(kp_fciqmc_data) :: kp
    interface
        subroutine NECICodeInit(tCPMD, tVASP, called_as_lib, MemSize)
            import :: dp, int64
            implicit none
            logical, intent(in) :: tCPMD, tVASP, called_as_lib
            integer(int64), intent(in), optional :: MemSize
        end subroutine
    end interface

    def_default(iCacheFlag_, iCacheFlag, 0)
    def_default(tCPMD_, tCPMD, .false.)
    def_default(tVASP_, tVASP, .false.)
    def_default(tMolpro, tMolpro_local, .false.)
    def_default(called_as_lib, call_as_lib, .false.)
    def_default(FCIDUMP_name, int_name, 'FCIDUMP')
    def_default(filename, filename_in, '')

#ifdef SX
    call stop_all(this_routine, 'The NEC compiler does not produce a working &
                                &version of NECI.')
#endif

    neci_MPIInit_called = .false.
    neci_MPINodes_called = .false.

    ! Do the program initialisation.
    call NECICodeInit(tCPMD_, tVASP_, called_as_lib, MemSize)

    proc_timer%timer_name = 'NECICUBE  '
    call set_timer(proc_timer)

!   See ReadInputMain.  Causes the command line arguments to be checked for the input filename.
    toverride_input = .false.

    if (tMolpro) then
        IF (molpro_plugin) THEN
            FCIDUMP_name = TRIM(molpro_plugin_fcidumpname)
        ELSE
            FCIDUMP_name = adjustl(int_name)
        END IF
        inquire (file="FCIQMC_input_override", exist=toverride_input)
        if (toverride_input) then
            Filename = "FCIQMC_input_override"
        else
            IF (molpro_plugin) THEN
                filename = TRIM(molpro_plugin_datafilename)
            ELSE
                filename = filename_in
            end if
            MolproID = ''
            if (iProcIndex == Root) then
                !Now, extract the unique identifier for the input file that is read in.
#ifdef old_and_buggy
                i = 14
                j = 1
                do while (filename(i:i) /= ' ')
                    MolproID(j:j) = filename(i:i)
                    i = i + 1
                    j = j + 1
                end do
#else
                i = INDEX(filename, '/', .TRUE.) + 1
                j = INDEX(filename(i:), 'NECI'); IF (j /= 0) i = i + j - 1
                MolproID = filename(i:MIN(i + LEN(MolproID) - 1, LEN(filename)))
#endif
                write(iout, "(A,A)") "Molpro unique filename suffix: ", MolproID
            end if
        end if
    end if

    ios = 0
    if (.not. (tCPMD_ .or. tVASP_)) then
        ! CPMD and VASP calculations call the input parser *before* they call
        ! NECICore.  This is to allow the NECI input filename(s) to be specified
        ! easily from within the CPMD/VASP input files.
        call ReadInputMain(Filename, ios, toverride_input, kp)
        If (ios /= 0) call stop_all(this_routine, 'Error in Read')
    end if

    call MPINodes(tUseProcsAsNodes)  ! Setup MPI Node information - this is dependent upon knowing the job type configurations.

    call NECICalcInit(iCacheFlag_)

!   Actually do the calculations we're meant to.  :-)
    call CalcDoCalc(kp)

!   And all done: pick up after ourselves and settle down for a cup of tea.
    call NECICalcEnd(iCacheFlag_)

    call halt_timer(proc_timer)

    if (tMolpro .and. (.not. toverride_input) .and. (.not. tMolproMimic)) then
        !Delete the FCIDUMP unless we are overriding the input, or mimicing molpro run-time behaviour
        if (iProcIndex == 0) then
            inquire (file='FCIDUMP', exist=tFCIDUMP_exist)
            if (tFCIDUMP_exist) then
                iunit = get_free_unit()
                open(iunit, file=FCIDUMP_name, status='old', form='formatted')
                close(iunit)
                !close(iunit,status='delete')
            end if
        end if
    end if
    call NECICodeEnd(tCPMD_, tVASP_)

    return
End Subroutine NECICore

subroutine NECICodeInit(tCPMD, tVASP, call_as_lib, MemSize)
    use MolproPlugin
    != Initialise the NECI code.  This contains all the initialisation
    != procedures for the code (as opposed to the system in general): e.g. for
    != the memory handling scheme, the timing routines, any parallel routines etc.
    != In:
    !=    tCPMD: True if doing a CPMD-based calculation.
    !=    tVASP: True if doing a VASP-based calculation.

    ! Utility modules
    use MemoryManager, only: InitMemoryManager
    use timing_neci, only: time_at_all, init_timing
    use Parallel_neci, only: MPIInit
    use SystemData, only: tMolpro
    use CalcData, only: s_global_start
    use constants, only: dp, int64
    use util_mod, only: neci_etime

    implicit none
    logical, intent(in) :: tCPMD, tVASP, call_as_lib
    integer(int64), intent(in), optional :: MemSize
    real(dp) :: tend(2)

    time_at_all = .not. call_as_lib

    ! MPIInit contains dummy initialisation for serial jobs, e.g. so we
    ! can refer to the processor index being 0 for the parent processor.
    Call MPIInit(tCPMD .or. tVASP .or. tMolpro .or. call_as_lib) ! CPMD and VASP have their own MPI initialisation and termination routines.

    ! Measure when NECICore is called. We need to do this here, as molcas
    ! and molpro can call NECI part way through a run, so it is no use to time
    ! from when the _process_ began.
    ! As this can use MPI_WTIME, we can only call this after the MPIInit call
    s_global_start = neci_etime(tend)

    ! find out whether this is a Molpro plugin
    CALL MolproPluginInit(tMolpro)
    ! end find out whether this is a Molpro plugin
    ! If we use MPI_WTIME for timing, we have to call MPIInit first
    call init_timing()

    if (.not. TCPMD) then
        call InitMemoryManager(MemSize)
    end if
    call environment_report(tCPMD)

end subroutine NECICodeInit

subroutine NECICodeEnd(tCPMD, tVASP)
    != End the NECI code.  This contains all the termination
    != procedures for the code (as opposed to the system in general): e.g. for
    != the memory handling scheme, the timing routines, any parallel routines etc.
    != In:
    !=    tCPMD: True if doing a CPMD-based calculation.
    !=    tVASP: True if doing a VASP-based calculation.

    ! Utility modules
    use MemoryManager, only: LeaveMemoryManager
    use timing_neci, only: end_timing, print_timing_report
    use SystemData, only: tMolpro, tMolproMimic, called_as_lib, arr, brr, g1, &
                          tagArr, tagBrr, tagG1
    use Determinants, only: FDet, tagFDet
#ifdef PARALLEL
    use Parallel_neci, only: MPIEnd
    USE MolproPlugin
#endif

    implicit none
    logical, intent(in) :: tCPMD, tVASP
    INTEGER :: rank, ierr

!    CALL N_MEMORY_CHECK

    ! Cleanup any memory that hasn't been deallocated elsewhere, and isn't
    ! immediately obvious where to deallocate it...

    if (.not. tCPMD .and. .not. called_as_lib) call LeaveMemoryManager()
    call end_timing()
    call print_timing_report()

#ifdef PARALLEL
! Tell Molpro plugin server that we have finished
    CALL MolproPluginTerm(0)
! CPMD and VASP have their own MPI initialisation and termination routines.
    call MPIEnd(molpro_plugin .or. (tMolpro .and. (.not. tMolproMimic)) .or. tCPMD .or. tVASP .or. called_as_lib)
#endif

end subroutine NECICodeEnd

subroutine NECICalcInit(iCacheFlag)
    != Calculation specific initialisation: just a wrapper for the individual
    != initiialisation routines for each part of the calculation.
    != In:
    !=    iCacheFlag: controls the behaviour of the 4-index integral cache.
    !=                Currently relevant only for CPMD and VASP calculations.
    !=                iCacheFlag=0,1: initialise the cache.
    !=                iCacheFlag=2,3: reuse the cache from the previous NECI
    !=                                calculation from within the same CPMD/VASP
    !=                                calculation.

    use System, only: SysInit
    use SystemData, only: tRotateOrbs, tFindCINatOrbs, tGUGA, tUEG, &
                          t_ueg_transcorr, t_ueg_dump, tContact, t_mol_3_body
    use Integrals_neci, only: IntInit, IntFreeze, tPostFreezeHF, DumpFCIDUMP
    use IntegralsData, only: tDumpFCIDUMP
    use DetCalc, only: DetCalcInit, DoDetCalc
    use Determinants, only: DetPreFreezeInit, DetInit, DetPreFreezeInit_old
    use Calc, only: CalcInit
    use HFCalc, only: HFDoCalc
    use RotateOrbsMod, only: RotateOrbs
    use replica_data, only: init_replica_arrays
    use gen_coul_ueg_mod, only: GEN_Umat_TC, prep_ueg_dump, GEN_Umat_TC_Contact
    use LMat_mod, only: readLMat
    use guga_init, only: init_guga
    implicit none
    integer, intent(in) :: iCacheFlag

    ! Initialise the global data arrays, whose size depends on global values
    ! that used to be constant.
    ! These are essentially constant arrays available after the input has
    ! been read
    call init_replica_arrays()

!   Initlialize the system.  Sets up ...
!   Symmetry is a subset of the system
    call SysInit()

    if (tGUGA) call init_guga

!   Initialize the integrals.  This will read in integrals, as well as calculating
!   some relevant integrals if they are calculated
    call IntInit(iCacheFlag)

!   If we need to freeze some integrals, do so.  Previously we used to do an HF
!   calculation here.  Instead we relegate that to the main Calc section, and are
!   required to read in the relevant orbitals if necessary.

!   This will also call SysPostFreezeInit()
    if (tGUGA) then
        call DetPreFreezeInit_old()
    else
        call DetPreFreezeInit()
    end if

    !! we prepare the contribution of the 2 body transcorrelated operator

    If (tUEG .and. t_ueg_transcorr) then
!                  CALL GetUMatSize(nBasis,nEl,UMATINT)
!                  call shared_allocate("umat_TC3", umat_TC3, (/UMatInt/))
!                  !allocate(UMat(UMatInt), stat=ierr)
!                  LogAlloc(ierr, 'UMat_TC3', int(UMatInt),HElement_t_SizeB, tagUMat)
!                  UMat_TC3 = 0.0_dp
!                  write(6,*) "Size of UMat_TC3 is: ",UMATINT

        write(6, *) 'prepare the convolution part of the 2 body transcorrelated operator'

        If (tContact) then
            call GEN_Umat_TC_contact
        else
            call GEN_Umat_TC
        end if
        write(6, *) "The infinite sums for the transcorrelated approach is determined."

        if (t_ueg_dump) call prep_ueg_dump

    !!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    end if

    if (.not. tPostFreezeHF) then
        call HFDoCalc()
    end if

    call IntFreeze()
    if(t_mol_3_body) call readLMat()
    ! can i initialize the GUGA stuff here? after freezing? or otherwise
    ! it is incorrectly setup..
    ! try to init guga here..
    if (tGUGA) call init_guga

    if (tPostFreezeHF) call HFDoCalc()

    if (tDumpFCIDUMP) then
        !If wanted, write out a new FCIDUMP file (with the core frozen if necessary)
        call DumpFCIDUMP()
    end if

    call DetInit()

! Deal with the many-electron basis, setting up sym etc.
! If we are doing a determinant-enumerated calc, this will generate lists of determinants.
!    call DetFreeze()

!   Do any determinant-list based calculations if we have them.  (diagonalizations probably)
    call DetCalcInit()
    call DoDetCalc()

!   Do any initializations we need to do for calculations (e.g. ...?)

    call CalcInit()

    IF (tRotateOrbs .and. (.not. tFindCINatOrbs)) THEN
        CALL RotateOrbs()
    end if

end subroutine NECICalcInit

subroutine NECICalcEnd(iCacheFlag)
    != Calculation specific termination: just a wrapper for the individual
    != termination routines: deallocation etc.
    != In:
    !=    iCacheFlag: controls the behaviour of the 4-index integral cache.
    !=                Currently relevant only for CPMD and VASP calculations.
    !=                iCacheFlag=0,2: destroy the cache on exit.
    !=                iCacheFlag=1,3: keep the cache on exit: it will be re-used
    !=                                in subsequent calculations within the same
    !=                                call to CPMD/VASP.

    ! Main level modules.
    use System, only: SysCleanup
    use Integrals_neci, only: IntCleanup
    use Determinants, only: DetCleanup
    use Calc, only: CalcCleanup
    use replica_data, only: clean_replica_arrays
    use OneEInts, only: DestroyTMat, DestroyPropInts
    use Parallel_neci, only: clean_parallel
    use SymExcitDataMod, only: SpinOrbSymLabel, SymInvLabel
    use SystemData, only: arr, brr, g1, tagArr, tagBrr, tagG1
    use Determinants, only: FDet, tagFDet
    use MemoryManager
    use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots
    use LoggingData, only: tCalcPropEst

    implicit none
    integer, intent(in) :: iCacheFlag

!   Tidy up:
    call CalcCleanup()
    call DetCleanup()
    call IntCleanup(iCacheFlag)
    call DestroyTMAT(.true.)
    call DestroyTMAT(.false.)
    if (tCalcPropEst) call DestroyPropInts
    call SysCleanup()
    call clean_replica_arrays()
    call clean_parallel()

    if (allocated(SpinOrbSymLabel)) deallocate(SpinOrbSymLabel)
    if (allocated(SymInvLabel)) deallocate(SymInvLabel)
    if (allocated(ValidSpawnedList)) deallocate(ValidSpawnedList)
    if (allocated(InitialSpawnedSlots)) deallocate(InitialSpawnedSlots)
    if (associated(arr)) then
        deallocate(Arr)
        call LogMemDealloc('NECICore', tagArr)
    end if
    if (associated(brr)) then
        deallocate(Brr)
        call LogMemDealloc('NECICore', tagBrr)
    end if
    if (associated(G1)) then
        deallocate(G1)
        call LogMemDealloc('NECICore', tagG1)
    end if
    if (associated(FDet)) then
        deallocate(FDet)
        call LogMemDealloc('NECICore', tagFDet)
    end if

    return
end subroutine NECICalcEnd
