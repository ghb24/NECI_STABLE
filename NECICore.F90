! NECICore is the main outline of the NECI Program.  It is called AFTER ReadInput has been called.

Subroutine  NECICore(iCacheFlag, tCPMD)
    use MemoryManager, only: LeaveMemoryManager
    use System, only : SysInit, SysCleanup
    use Integrals, only : IntInit, IntFreeze, IntCleanup
    use DetCalc, only : DetCalcInit, DoDetCalc
    use Determinants, only : DetPreFreezeInit, DetInit, DetCleanup
    use Calc, only : CalcInit, CalcDoCalc, CalcCleanup
    use HFCalc, only: HFDoCalc
    Implicit none
!Set by CPMD to determine whether cache is saved
    Integer iCacheFlag
    INTEGER iSub
    Logical tCPMD

      call TimeTag()
      if (.not.TCPMD) call Envir()
      write (6,*)
      CALL TISET('NECICUBE  ',ISUB)



! Initlialize the system.  Sets up ...
!   Symmetry is a subset of the system
    call SysInit()
! Initialize the integrals.  This will read in integrals, as well as calculating some relevant integrals if they are calculated
    call IntInit(iCacheFlag)
! If we need to freeze some integrals, do so.
!  Previously we used to do an HF calculation here.  Instead we relegate that to the main Calc section, and are required to
!    read in the relevant orbitals if necessary.

!  This will also call SysPostFreezeInit()
    call DetPreFreezeInit()
    call HFDoCalc()
    call IntFreeze()
    call DetInit()
!  Deal with the many-electron basis, setting up sym etc.
!  If we are doing a determinant-enumerated calc, this will generate lists of determinants.
!    call DetFreeze()
!  Do any determinant-list based calculations if we have them.  (diagonalizations probably)
    call DetCalcInit()
    call DoDetCalc()
!  Do any initializations we need to do for calculations (e.g. ...?)
    call CalcInit()
!  Actually do the calculations we're meant to
    call CalcDoCalc()
! tidy up
    call CalcCleanup()
    call DetCleanup()
    call IntCleanup(iCacheFlag)
    call SysCleanup()
          CALL MEMORY_CHECK
! ==-------------------------------------------------------------------==
      CALL TIHALT('NECICUBE  ',ISUB)
! ==-------------------------------------------------------------------==
      CALL LeaveMemoryManager
      CALL TIPRI

End Subroutine NECICore
