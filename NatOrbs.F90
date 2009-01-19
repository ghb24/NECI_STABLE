!This file is primarily concerned with the creation of natural orbitals from a rotation of the previous orbitals.
!The 1-electron Reduced density matrix will be inputted, and the natural orbitals constructed. From there, the
!1 and 2 electron integrals will be transformed and replaced into UMat.


SUBROUTINE FindNatOrbs(OneRDM)

!First, diagonalize the 1-RDM...
    CALL Diag1RDM(OneRDM)

!Setup symmetry information needed...
!    CALL SetupSymNO()

!Find orbitals energy...
!    CALL FindOrbEnergies()

!Transform integrals...
!    CALL TransRotInts()

END SUBROUTINE FindNatOrbs

SUBROUTINE Diag1RDM(OneRDM)
    USE SystemData , only : NEl,nBasis
    USE Global_Utilities
    IMPLICIT NONE
    REAL*8 :: OneRDM(nBasis,nBasis)
    REAL*8 , ALLOCATABLE :: NOccNums(:),Work(:)
    INTEGER :: nOccNumsTag=0,iErr,WorkSize,WorkCheck,WorkTag=0,i
    CHARACTER(len=*), PARAMETER :: this_routine='Diag1RDM'

    ALLOCATE(NOccNums(nBasis),stat=ierr)
    CALL LogMemAlloc('NOccNums',nBasis,8,this_routine,NOccNumsTag,iErr)

!Find desired optimal workspace size
    WorkSize=-1
    CALL DSYEV('V','U',nBasis,OneRDM,nBasis,NOccNums,WorkCheck,WorkSize,iErr)
    WorkSize=WorkCheck

    WRITE(6,*) "Optimal size of scratch space for diagonalization = ",WorkCheck

    ALLOCATE(Work(WorkSize),stat=ierr)
    CALL LogMemAlloc('Work',WorkSize,8,this_routine,WorkTag,iErr)

    WRITE(6,"(A)",advance='no') "Diagonalizing 1-RDM to find approximate natural orbitals..."
    CALL FLUSH(6)

!Diagonalize... Matrix must be symmetric
    CALL DSYEV('V','U',nBasis,OneRDM,nBasis,NOccNums,Work,WorkSize,iErr)
    IF(iErr.eq.0) THEN
        WRITE(6,"(A)") "DONE!"
    ELSE
        WRITE(6,"(A)") "FAILED!"
        CALL Stop_All(this_routine,"Diagonalization of 1-RDM failed...")
    ENDIF

    DEALLOCATE(Work)
    CALL LogMemDealloc(this_routine,WorkTag)

    WRITE(6,*) "Occupation numbers of approximate natural spin-orbitals are: "
    do i=1,nBasis
        WRITE(6,*) i,NOccNums(i)    !Symmetry of orbital...?
    enddo
    CALL FLUSH(6)

END SUBROUTINE Diag1RDM
