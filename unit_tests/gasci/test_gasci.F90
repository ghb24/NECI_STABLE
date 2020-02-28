program test_gasci

  use disconnected_gasci
  use fruit
  use constants
  use read_fci, only: readfciint, initfromfcid, fcidump_name
  use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
  use IntegralsData, only: UMat, umat_win
  use Integrals_neci, only: IntInit, get_umat_el_normal
  use procedure_pointers, only: get_umat_el
  use SystemData, only: nel, nBasis, UMatEps, tStoreSpinOrbs, tReadFreeFormat, tCSF, &
       tReadInt, tGASSpinRecoupling
  use System, only: SysInit, SetSysDefaults
  use Parallel_neci, only: MPIInit
  use UMatCache, only: GetUMatSize, tTransGTID
  use OneEInts, only: Tmat2D
  use bit_rep_data, only: NIfTot, NIfDBO, NOffSgn, NIfSgn, extract_sign
  use bit_reps, only: encode_sign, decode_bit_det
  use DetBitOps, only: EncodeBitDet, DetBitEq
  use SymExcit3, only: countExcitations3, GenExcitations3
  use dSFMT_interface, only: dSFMT_init
  use FciMCData, only: pSingles, pDoubles, pParallel, excit_gen_store_type
  use Calc, only: CalcInit, SetCalcDefaults
  use Determinants, only: DetInit, DetPreFreezeInit, get_helement_normal

  implicit none

  call init_fruit()

  call gasci_test_driver()

  call fruit_summary()
  call fruit_finalize()

contains

  subroutine gasci_test_driver()
    implicit none
    integer :: nBasisMax(5,3), lms
    integer(int64) :: umatsize
    real(dp) :: ecore
    integer(MPIArg) :: ierr
    umatsize = 0
    nel = 5
    NIfDBO = 0
    NOffSgn = 1
    NIfSgn = 1
    NIfTot = 2

    fcidump_name = "FCIDUMP"
    UMatEps = 1.0e-8
    tStoreSpinOrbs = .false.
    tTransGTID = .false.
    tReadFreeFormat = .true.

    call MPIInit(.false.)

    call dSFMT_init(8)

    call SetCalcDefaults()
    call SetSysDefaults()
    tReadInt = .true.

    get_umat_el => get_umat_el_normal

    call initfromfcid(nel,nbasismax,nBasis,lms,.false.)

    call GetUMatSize(nBasis, umatsize)

    allocate(TMat2d(nBasis,nBasis))

    call shared_allocate_mpi(umat_win, umat, (/umatsize/))

    call readfciint(UMat,umat_win,nBasis,ecore,.false.)

    call SysInit()
    ! required: set up the spin info

    call DetInit()
    ! call SpinOrbSymSetup()

    call DetPreFreezeInit()

    call CalcInit()

    !call test_splitIndex()
    call loadGAS()
    call test_gasci_exgen()

    deallocate(TMat2D)
    call shared_deallocate_mpi(umat_win, UMat)

    call MPI_finalize(ierr)

  end subroutine gasci_test_driver

  subroutine test_gasci_exgen()
    implicit none
    integer :: nI(nel), nJ(nel)
    integer :: i, ex(2,maxExcit), exflag
    integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
    real(dp) :: pgen
    logical :: tPar, tAllExFound, tFound
    integer :: j, numEx, nSingles, nDoubles
    integer(n_int), allocatable :: allEx(:,:)
    integer, parameter :: sampleSize = 100000
    real(dp) :: pgenArr(lenof_sign), pTot, pNull
    HElement_t(dp) :: hel
    type(excit_gen_store_type), target :: store
    logical :: exDone(nel,nel,0:nBasis,0:nBasis)
    integer :: ic

    exDone = .false.

    ! some starting det
    do i = 1, nel
       nI(i) = i
    end do

    tCSF = .false.
    tGASSpinRecoupling = .false.
    call EncodeBitDet(nI,ilut)

    exflag = 3
    ! create a list of all singles and doubles for reference
    call CountExcitations3(nI,exflag,nSingles,nDoubles)
    allocate(allEx(0:(NIfTot+1),nSingles+nDoubles), source = 0_n_int)
    numEx = 0
    tAllExFound = .false.
    do
       call GenExcitations3(nI,ilut,nJ,exflag,ex,tPar,tAllExFound,.false.)
       if(tAllExFound) exit
       call encodeBitDet(nJ,ilutJ)
       numEx = numEx + 1
       allEx(0:NIfDBO,numEx) = ilutJ(0:NIfDBO)
    end do

    ! set the biases for excitation generation
    pParallel = 0.5_dp
    pSingles = 0.1_dp
    pDoubles = 0.9_dp

    pNull = 0.0_dp
    do i = 1, sampleSize
       call generate_nGAS_excitation(nI,ilut,nJ,ilutJ,exFlag,ic,&
            ex,tpar,pgen,hel,store,1)
       ! lookup the excitation
       tFound = .false.
       do j = 1, numEx
          if(DetBitEQ(ilutJ,allEx(:,j))) then
             if(allEx(NIfTot+1,j) > 1) then
                call extract_sign(allEx(:,j),pgenArr)
                if(abs(pgenArr(1)-pgen) > eps) then
                   write(iout,*) "ERROR: Conflicting pgen", pgenArr, pgen
                   write(iout,*) "Determinant", nJ
                   stop
                endif
             endif
             pgenArr = pgen
             call encode_sign(allEx(:,j),pgenArr)
             allEx(NIfTot+1,j) = allEx(NIfTot+1,j) + 1
             tFound = .true.
             exit
          end if
       end do
       if(.not. tFound .and. .not. nJ(1)==0) then
          call decode_bit_det(nJ,ilutJ)
          write(iout,*) "Error: Invalid excitation", nJ
          exit
       endif
       if(nJ(1) == 0) then
          if(.not. exDone(ex(1,1),ex(1,2),ex(2,1),ex(2,2)))then
             exDone(ex(1,1),ex(1,2),ex(2,1),ex(2,2)) = .true.
             pNull = pNull + pgen
          endif
       endif
    end do

    ! check that all excits have been generated and all pGens are right
    ! probability normalization
    pTot = pNull
    do i = 1, numEx
       call extract_sign(allEx(:,i),pgenArr)
       call decode_bit_det(nJ,allEx(:,i))
       write(iout,*) i, pgenArr(1), real(allEx(NIfTot+1,i))/real(sampleSize)
!       write(iout,*), nJ
       if(pgenArr(1) < eps .and. isValidExcit(allEx(:,i),ilut)) then
          if(abs(get_helement_normal(nJ,nI)) > eps) then
             write(iout,*) "Warning: missed det", nJ
             write(iout,*) "With matel", get_helement_normal(nI,nJ)
          endif
       endif
       pTot = pTot + pgenArr(1)
       if(pgenArr(1)<0) then
          write(iout,*) "Error: Negative pgen"
          exit
       endif
    end do
    write(iout,*) "Total prob. ", pTot
    write(iout,*) "pNull ", pNull

  end subroutine test_gasci_exgen
end program test_gasci
