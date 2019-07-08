program test_time_search

  use UMatHash
  use constants
  use read_fci, only: readfciint, initfromfcid, fcidump_name
  use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
  use IntegralsData, only: UMat, umat_win
  use Integrals_neci, only: IntInit, get_umat_el_normal
  use procedure_pointers, only: get_umat_el
  use SystemData, only: nel, nBasis, UMatEps, tStoreSpinOrbs, tReadFreeFormat, tCSF, &
       tReadInt
  use System, only: SysInit, SetSysDefaults
  use Parallel_neci, only: MPIInit, MPIEnd
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

  abstract interface
     function csum_search_t(arr, thresh)
       use constants
       real(dp), intent(in) :: arr(:)
       real(dp), intent(in) :: thresh
       integer :: csum_search_t
     end function csum_search_t
  end interface

  integer :: nBasisMax(5,3), lms
  integer(int64) :: umatsize
  real(dp) :: ecore
  procedure(csum_search_t), pointer :: csum_search
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

  write(iout,*) "Starting timing of linear search"
  csum_search => linearSearch
  call time_search(csum_search_wrapper)
  write(iout,*) "Starting timing of binary search"
  csum_search => binary_search_wrapper
  call time_search(csum_search_wrapper)
  write(iout,*) "Starting timing of false position search"
  csum_search => linearSearch_old
  call time_search(csum_search_wrapper)
  write(iout,*) "Starting timing of cache line search"
  csum_search => cacheLineSearch
  call time_search(csum_search_wrapper)
  write(iout,*) "Starting timing of alias method"
  call time_search(aliasSampling_wrapper)

  deallocate(TMat2D)
  call shared_deallocate_mpi(umat_win, UMat)

  call MPIEnd(.false.)

contains

  subroutine time_search(searchAlg)
    implicit none
    integer, parameter :: numTries = 1000000
    integer:: ind, pq, startPQ, endPQ, j
    real(dp) :: t1, t2
    integer(int64) :: maxInd, checksum

    interface 
       function searchAlg(startPQ, endPQ)
         use constants
         integer, intent(in) :: startPQ, endPQ
         integer :: searchAlg
       end function searchAlg
    end interface

    maxInd = fuseIndex(nStoreBasis,nStoreBasis)

    call cpu_time(t1)
    checksum = 0_int64
    do j = 1, numTries
       pq = int(genrand_real2_dSFMT()*maxInd)+1

       startPQ = posPQ(pq)
       endPQ = startPQ + nPQ(pq) - 1

       ind = searchAlg(startPQ,endPQ) + startPQ - 1
       checksum = mod(checksum + int(ind,int64),2_int64**31)
    end do
    call cpu_time(t2)
    write(iout,*) "Time taken: ", t2-t1
    write(iout,*) checksum
  end subroutine

  function aliasSampling_wrapper(startPQ, endPQ) result(ind)
    implicit none
    integer, intent(in) :: startPQ, endPQ
    integer :: ind

    ind = aliasSampling(biasTable(startPQ:endPQ), aliasTable(startPQ:endPQ))
    
  end function aliasSampling_wrapper
    
  function csum_search_wrapper(startPQ, endPQ) result(ind)
    implicit none
    integer, intent(in) :: startPQ, endPQ
    real(dp) :: thresh
    integer :: ind

    thresh = genrand_real2_dSFMT()*CumSparseUMat(endPQ)

    ind = csum_search(CumSparseUMat(startPQ:endPQ),thresh)
    
  end function csum_search_wrapper

  function binary_search_wrapper(arr, thresh) result(ind)
    implicit none
    real(dp), intent(in) :: arr(:)
    real(dp), intent(in) :: thresh
    integer :: ind

    ind = binary_search_first_ge(arr,thresh)
    
  end function binary_search_wrapper
end program
