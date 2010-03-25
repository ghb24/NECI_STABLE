Subroutine NECICore(iCacheFlag,tCPMD,tVASP)
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

    Use ReadInput, only : ReadInputMain

    ! main-level modules.
    use Calc, only: CalcDoCalc

    ! Utility modules.
    use global_utilities

    Implicit none
    integer,intent(in) :: iCacheFlag
    logical,intent(in) :: tCPMD,tVASP
    type(timer), save :: proc_timer
    integer :: ios
    character(255) :: Filename
    ! Do the program initialisation.
    call NECICodeInit(tCPMD,tVASP)

    proc_timer%timer_name='NECICUBE  '
    call set_timer(proc_timer)

!   See ReadInputMain.  Causes the command line arguments to be checked for the input filename.
    Filename="" 
    ios=0

    if (.not.tCPMD.and..not.tVASP) then
        ! CPMD and VASP calculations call the input parser *before* they call
        ! NECICore.  This is to allow the NECI input filename(s) to be specified
        ! easily from within the CPMD/VASP input files.
        call ReadInputMain(Filename,ios)
        If (ios.ne.0) stop 'Error in Read'
    endif

    call NECICalcInit(iCacheFlag)

!   Actually do the calculations we're meant to.  :-)
    call CalcDoCalc()

!   And all done: pick up after ourselves and settle down for a cup of tea.
    call NECICalcEnd(iCacheFlag)

    call halt_timer(proc_timer)

    call NECICodeEnd(tCPMD,tVASP)

    return
End Subroutine NECICore



subroutine NECICodeInit(tCPMD,tVASP)
    != Initialise the NECI code.  This contains all the initialisation
    != procedures for the code (as opposed to the system in general): e.g. for
    != the memory handling scheme, the timing routines, any parallel routines etc.
    != In:
    !=    tCPMD: True if doing a CPMD-based calculation.
    !=    tVASP: True if doing a VASP-based calculation.

    ! Utility modules
    use MemoryManager, only: InitMemoryManager
    use timing, only: init_timing
    use Parallel, only: MPIInit

    implicit none
    logical, intent(in) :: tCPMD,tVASP

    call init_timing()

    ! MPIInit contains dummy initialisation for serial jobs, e.g. so we
    ! can refer to the processor index being 0 for the parent processor.
    Call MPIInit(tCPMD.or.tVASP) ! CPMD and VASP have their own MPI initialisation and termination routines.

    if (.not.TCPMD) then
        call InitMemoryManager()
    end if
    call environment_report(tCPMD)

end subroutine NECICodeInit



subroutine NECICodeEnd(tCPMD,tVASP)
    != End the NECI code.  This contains all the termination 
    != procedures for the code (as opposed to the system in general): e.g. for
    != the memory handling scheme, the timing routines, any parallel routines etc.
    != In:
    !=    tCPMD: True if doing a CPMD-based calculation.
    !=    tVASP: True if doing a VASP-based calculation.

    ! Utility modules
    use MemoryManager, only: LeaveMemoryManager
    use timing, only: end_timing,print_timing_report
#ifdef PARALLEL
    use Parallel, only: MPIEnd
#endif

    implicit none
    logical, intent(in) :: tCPMD,tVASP

#ifdef PARALLEL
    call MPIEnd(tCPMD.or.tVASP) ! CPMD and VASP have their own MPI initialisation and termination routines.
#endif

    CALL N_MEMORY_CHECK

    if (.not.tCPMD) call LeaveMemoryManager()
    call end_timing()
    call print_timing_report()

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

    use System, only : SysInit
    use SystemData, only : tRotateOrbs,tFindCINatOrbs
    use Integrals, only : IntInit,IntFreeze,tPostFreezeHF
    use DetCalc, only : DetCalcInit,DoDetCalc
    use Determinants, only : DetPreFreezeInit,DetInit
    use Calc, only : CalcInit
    use HFCalc, only: HFDoCalc
    use RotateOrbsMod, only : RotateOrbs

    implicit none
    integer,intent(in) :: iCacheFlag
   

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
    
    IF(tRotateOrbs.and.(.not.tFindCINatOrbs)) THEN
        CALL RotateOrbs()
    ENDIF


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
    use Integrals, only: IntCleanup
    use Determinants, only: DetCleanup
    use Calc, only: CalcCleanup

    implicit none
    integer,intent(in) :: iCacheFlag

!   Tidy up: 
    call CalcCleanup()
    call DetCleanup()
    call IntCleanup(iCacheFlag)
    call SysCleanup()

    return
end subroutine NECICalcEnd
