! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
Subroutine NECICore(iCacheFlag,tCPMD,tVASP,tMolpro_local,int_name,filename_in)
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
    !=    int_name is the name of the integral file to read in if necessary
    !=    filename is the name of the input file to read in if necessary

    use ReadInput_neci, only : ReadInputMain
    use SystemData, only : tMolpro,tMolproMimic,MolproID

    ! main-level modules.
    use Calc, only: CalcDoCalc
    use CalcData, only: tUseProcsAsNodes
    use Parallel_neci, only: MPINodes, iProcIndex
    use read_fci, only: FCIDUMP_name
    use ParallelHelper, only: Root

    ! Utility modules.
    use global_utilities
    use util_mod, only: get_free_unit

    Implicit none
    integer,intent(in) :: iCacheFlag
    logical,intent(in) :: tCPMD,tVASP,tMolpro_local
    character(64), intent(in) :: filename_in,int_name
    type(timer), save :: proc_timer
    integer :: ios,iunit,i,j
    character(*), parameter :: this_routine = 'NECICore'
    character(64) :: Filename
    logical :: toverride_input,tFCIDUMP_exist
#ifdef MOLPRO
    include "common/tapes"
#else
    integer, parameter :: iout = 6
#endif
    
    tMolpro = tMolpro_local

#ifdef SX
    call stop_all(this_routine, 'The NEC compiler does not produce a working &
                                &version of NECI.')
#endif

    ! Do the program initialisation.
    call NECICodeInit(tCPMD,tVASP)

    proc_timer%timer_name='NECICUBE  '
    call set_timer(proc_timer)

!   See ReadInputMain.  Causes the command line arguments to be checked for the input filename.
    if(tMolpro) then
        FCIDUMP_name = adjustl(int_name)
        inquire(file="FCIQMC_input_override",exist=toverride_input)
        if(toverride_input) then
            Filename="FCIQMC_input_override"
        else
            filename=filename_in
            MolproID = ''
            if(iProcIndex.eq.Root) then
            !Now, extract the unique identifier for the input file that is read in.
                i=14
                j=1
                do while(filename(i:i).ne.' ')
                    MolproID(j:j) = filename(i:i)
                    i=i+1
                    j=j+1
                enddo
                write(iout,"(A,A)") "Molpro unique filename suffix: ",MolproID
            endif
        endif
    else
        Filename='' 
        FCIDUMP_name = 'FCIDUMP'
    endif
    ios=0

    if (.not.tCPMD.and..not.tVASP) then
        ! CPMD and VASP calculations call the input parser *before* they call
        ! NECICore.  This is to allow the NECI input filename(s) to be specified
        ! easily from within the CPMD/VASP input files.
        call ReadInputMain(Filename,ios,toverride_input)
        If (ios.ne.0) stop 'Error in Read'
    endif

    call MPINodes(tUseProcsAsNodes)  ! Setup MPI Node information - this is dependent upon knowing the job type configurations.

    call NECICalcInit(iCacheFlag)

!   Actually do the calculations we're meant to.  :-)
    call CalcDoCalc()

!   And all done: pick up after ourselves and settle down for a cup of tea.
    call NECICalcEnd(iCacheFlag)

    call halt_timer(proc_timer)
    
    if(tMolpro.and.(.not.toverride_input).and.(.not.tMolproMimic)) then
        !Delete the FCIDUMP unless we are overriding the input, or mimicing molpro run-time behaviour
        if(iProcIndex.eq.0) then
            inquire(file='FCIDUMP',exist=tFCIDUMP_exist)
            if(tFCIDUMP_exist) then
                iunit = get_free_unit()
                open(iunit,file=FCIDUMP_name,status='old',form='formatted')
                close(iunit,status='delete')
            endif
        endif
    endif


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
    use timing_neci, only: init_timing
    use Parallel_neci, only: MPIInit
    use SystemData, only : tMolpro

    implicit none
    logical, intent(in) :: tCPMD,tVASP

    call init_timing()

    ! MPIInit contains dummy initialisation for serial jobs, e.g. so we
    ! can refer to the processor index being 0 for the parent processor.
    Call MPIInit(tCPMD.or.tVASP.or.tMolpro) ! CPMD and VASP have their own MPI initialisation and termination routines.

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
    use timing_neci, only: end_timing,print_timing_report
    use SystemData, only : tMolpro,tMolproMimic
#ifdef PARALLEL
    use Parallel_neci, only: MPIEnd
#endif

    implicit none
    logical, intent(in) :: tCPMD,tVASP

#ifdef PARALLEL
! CPMD and VASP have their own MPI initialisation and termination routines.
    call MPIEnd((tMolpro.and.(.not.tMolproMimic)).or.tCPMD.or.tVASP) 
#endif

!    CALL N_MEMORY_CHECK

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
    use Integrals_neci, only : IntInit,IntFreeze,tPostFreezeHF
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
    use Integrals_neci, only: IntCleanup
    use Determinants, only: DetCleanup
    use Calc, only: CalcCleanup
    use shared_alloc, only: cleanup_shared_alloc

    implicit none
    integer,intent(in) :: iCacheFlag

!   Tidy up: 
    call CalcCleanup()
    call DetCleanup()
    call IntCleanup(iCacheFlag)
    call SysCleanup()
    call cleanup_shared_alloc ()

    return
end subroutine NECICalcEnd
