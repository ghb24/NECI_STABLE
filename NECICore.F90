! NECICore is the main outline of the NECI Program.  It is called AFTER ReadInput has been called.

Subroutine  NECICore(iCacheFlag, tCPMD,tVASP)
    use MemoryManager, only: InitMemoryManager,LeaveMemoryManager
    use System, only : SysInit, SysCleanup
    use Integrals, only : IntInit, IntFreeze, IntCleanup, tPostFreezeHF
    use DetCalc, only : DetCalcInit, DoDetCalc
    use Determinants, only : DetPreFreezeInit, DetInit, DetCleanup
    use Calc, only : CalcInit, CalcDoCalc, CalcCleanup
    use HFCalc, only: HFDoCalc
    use timing, only: init_timing,end_timing
#ifdef PARALLEL
    use Parallel, only: MPIInit, MPIEnd
#endif
    Use ReadInput, only : ReadInputMain
    use timing
    Implicit none
!Set by CPMD to determine whether cache is saved
    Integer iCacheFlag
    INTEGER,SAVE :: iSub=0
    Logical tCPMD,tVASP
    integer ios
    character(255) Filename

    call init_timing()

    Filename="";
    ios=0

#ifdef PARALLEL
    Call MPIInit(tCPMD.or.tVASP)
#endif
    call environment_report(tCPMD)
    if (.not.TCPMD) then
        call InitMemoryManager()
    end if
    write (6,*) 'tVASP',tVASP
    call set_timer('NECICUBE  ',ISUB)
    if(.not.tCPMD.and..not.tVASP) THEN
        call ReadInputMain(Filename,ios)
        If (ios.ne.0) stop 'Error in Read'
    endif

! Initlialize the system.  Sets up ...
! Symmetry is a subset of the system
    call SysInit()
! Initialize the integrals.  This will read in integrals, as well as calculating
! some relevant integrals if they are calculated
    call IntInit(iCacheFlag)

! If we need to freeze some integrals, do so.  Previously we used to do an HF
! calculation here.  Instead we relegate that to the main Calc section, and are
! required to read in the relevant orbitals if necessary.

! This will also call SysPostFreezeInit()
    call DetPreFreezeInit()
    if(.not.tPostFreezeHF) call HFDoCalc()
    call IntFreeze()
    if(tPostFreezeHF) call HFDoCalc()
    call DetInit()

! Deal with the many-electron basis, setting up sym etc.
! If we are doing a determinant-enumerated calc, this will generate lists of determinants.
!    call DetFreeze()

! Do any determinant-list based calculations if we have them.  (diagonalizations probably)
    call DetCalcInit()
    call DoDetCalc()

! Do any initializations we need to do for calculations (e.g. ...?)
    call CalcInit()

! Actually do the calculations we're meant to
    call CalcDoCalc()

! tidy up
    call CalcCleanup()
    call DetCleanup()
    call IntCleanup(iCacheFlag)
    call SysCleanup()
#ifdef PARALLEL
    call MPIEnd(tCPMD)
#endif

    CALL MEMORY_CHECK

! ==-------------------------------------------------------------------==
    call halt_timer(ISUB)
! ==-------------------------------------------------------------------==

    if (.not.tCPMD) call LeaveMemoryManager
    call end_timing()
    call print_timing_report()

    return
End Subroutine NECICore
