program test_umat_hash

  use UMatHash
  use cachedExcitgen
  use fruit
  use constants
  use read_fci, only: readfciint, initfromfcid, fcidump_name
  use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
  use IntegralsData, only: UMat, umat_win
  use Integrals_neci, only: IntInit, get_umat_el_normal
  use procedure_pointers, only: get_umat_el
  use SystemData, only: nel, nBasis, UMatEps, tStoreSpinOrbs, tReadFreeFormat, tCSF, &
       tReadInt
  use System, only: SysInit, SetSysDefaults
  use Parallel_neci, only: MPIInit
  use UMatCache, only: GetUMatSize, tTransGTID
  use OneEInts, only: Tmat2D
  use bit_rep_data, only: NIfTot, NIfDBO, NOffSgn, NIfSgn, extract_sign
  use bit_reps, only: encode_sign, decode_bit_det
  use DetBitOps, only: EncodeBitDet, DetBitEq
  use SymExcit3, only: countExcitations3, GenExcitations3
  use dSFMT_interface, only: dSFMT_init
  use FciMCData, only: pSingles, pDoubles, pParallel
  use Calc, only: CalcInit, SetCalcDefaults
  use Determinants, only: DetInit, DetPreFreezeInit, get_helement_normal

  implicit none

  call init_fruit()

  call umat_hash_test_driver()

  call fruit_summary()
  call fruit_finalize()

contains

  subroutine umat_hash_test_driver()
    implicit none
    integer :: nBasisMax(5,3), lms
    integer(int64) :: umatsize
    real(dp) :: ecore
    integer :: ierr
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

    call GetUMatSize(nBasis, nel, umatsize)

    allocate(TMat2d(nBasis,nBasis))

    call shared_allocate_mpi(umat_win, umat, (/umatsize/))

    call readfciint(UMat,umat_win,nBasis,ecore,.false.)

    call SysInit()
    ! required: set up the spin info
    
    call DetInit()
    ! call SpinOrbSymSetup()
    
    call DetPreFreezeInit()

    call CalcInit()

    call initializeSparseUMat()

    !call test_splitIndex()
    call test_cached_exgen()

    deallocate(TMat2D)
    call shared_deallocate_mpi(umat_win, UMat)

    call MPI_finalize(ierr)

  end subroutine umat_hash_test_driver

  subroutine test_splitIndex()
    implicit none
    integer, parameter :: nInds = 50
    integer :: p,q,pq,nPQ
    logical :: check(nInds,nInds)
    integer :: backup

    backup = nStoreBasis
    nStoreBasis = nInds
    nPQ = fuseIndex(nInds,nInds)
    do pq = 1, nPQ
       call splitIndex(pq,p,q)
       if(check(p,q)) then
          write(iout,*) "Error - repeated index"
          stop
       endif
       check(p,q) = .true.
    end do

    do p = 1, nInds
       do q = p, nInds
          if(.not. check(p,q)) then
             write(iout,*) "Error - not all indices mapped", p,q
          endif
       enddo
    end do
    nStoreBasis = backup
  end subroutine test_splitIndex
  

  subroutine test_cached_exgen()
    implicit none
    integer :: nI(nel), nJ(nel)
    integer :: i, ex(2,2), exflag
    integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
    real(dp) :: pgen
    logical :: tPar, tAllExFound, tFound
    integer :: j, numEx, nSingles, nDoubles
    integer(n_int), allocatable :: allEx(:,:)
    integer, parameter :: sampleSize = 1000000
    real(dp) :: pgenArr(lenof_sign), pTot, pNull
    logical :: exDone(nel,nel,nBasis,nBasis)
    
    exDone = .false.

    ! some starting det
    do i = 1, nel
       nI(i) = i
    end do

    tCSF = .false.
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
       call generate_double_cached(nI,ilut,nJ,ilutJ,ex,pgen,tpar)
       ! lookup the excitation
       tFound = .false.
       do j = 1, numEx
          if(DetBitEQ(ilutJ,allEx(:,j))) then
             pgenArr = pgen
             call encode_sign(allEx(:,j),pgenArr)
             allEx(NIfTot+1,j) = allEx(NIfTot+1,j) + 1
             tFound = .true.
             exit
          end if
       end do
       if(.not. tFound .and. .not. all(nJ==0)) then
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
       if(pgenArr(1) < eps) &
           write(iout,*) "UMAT el ", get_helement_normal(nJ,(/1,2,3,4,5/))
       pTot = pTot + pgenArr(1)
    end do
    write(iout,*) "Total prob. ", pTot
    write(iout,*) "pNull ", pNull

  end subroutine test_cached_exgen
end program test_umat_hash
