Subroutine  NECICore(iCacheFlag,tCPMD,tVASP)
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

    ! Main level modules.
    use System, only : SysInit, SysCleanup
    use Integrals, only : IntInit, IntFreeze, IntCleanup, tPostFreezeHF
    use DetCalc, only : DetCalcInit, DoDetCalc
    use Determinants, only : DetPreFreezeInit, DetInit, DetCleanup
    use Calc, only : CalcInit, CalcDoCalc, CalcCleanup
    use HFCalc, only: HFDoCalc

    Use ReadInput, only : ReadInputMain

    ! Utility modules.
#ifdef PARALLEL
    use Parallel, only: MPIInit, MPIEnd
#endif
    use MemoryManager, only: InitMemoryManager,LeaveMemoryManager
    use soft_exit
    use timing

    Implicit none
    Integer :: iCacheFlag
    type(timer), save :: proc_timer
    Logical :: tCPMD,tVASP
    integer :: ios
    character(255) :: Filename

    call init_timing()

!   See ReadInputMain.  Causes the command line arguments to be checked for the input filename.
    Filename="" 
    ios=0

#ifdef PARALLEL
    Call MPIInit(tCPMD.or.tVASP) ! CPMD and VASP have their own MPI initialisation and termination routines.
#endif

    if (.not.TCPMD) then
        call InitMemoryManager()
    end if
    call init_soft_exit()
    call environment_report(tCPMD)

    proc_timer%timer_name='NECICUBE  '
    call set_timer(proc_timer)

    if (.not.tCPMD.and..not.tVASP) then
        ! CPMD and VASP calculations call the input parser *before* they call
        ! NECICore.  This is to allow the NECI input filename(s) to be specified
        ! easily from within the CPMD/VASP input files.
        call ReadInputMain(Filename,ios)
        If (ios.ne.0) stop 'Error in Read'
    endif

!   Initlialize the system.  Sets up ...
!   Symmetry is a subset of the system
    call SysInit()

!   Initialize the integrals.  This will read in integrals, as well as calculating
!   some relevant integrals if they are calculated
    call IntInit(iCacheFlag)

!   If we need to freeze some integrals, do so.  Previously we used to do an HF
!   calculation here.  Instead we relegate that to the main Calc section, and are
!   required to read in the relevant orbitals if necessary.

!   This will also call SysPostFreezeInit()
    call DetPreFreezeInit()
    if (.not.tPostFreezeHF) call HFDoCalc()
    call IntFreeze()
    if (tPostFreezeHF) call HFDoCalc()
    call DetInit()

! Deal with the many-electron basis, setting up sym etc.
! If we are doing a determinant-enumerated calc, this will generate lists of determinants.
!    call DetFreeze()

!   Do any determinant-list based calculations if we have them.  (diagonalizations probably)
    call DetCalcInit()
    call DoDetCalc()

!   Do any initializations we need to do for calculations (e.g. ...?)
    call CalcInit()

!   Actually do the calculations we're meant to.  :-)
    call CalcDoCalc()

!   Tidy up: 
    call CalcCleanup()
    call DetCleanup()
    call IntCleanup(iCacheFlag)
    call SysCleanup()

#ifdef PARALLEL
    call MPIEnd(tCPMD.or.tVASP) ! CPMD and VASP have their own MPI initialisation and termination routines.
#endif

    CALL MEMORY_CHECK

    call halt_timer(proc_timer)

    call end_soft_exit()
    if (.not.tCPMD) call LeaveMemoryManager()
    call end_timing()
    call print_timing_report()

    return
End Subroutine NECICore
