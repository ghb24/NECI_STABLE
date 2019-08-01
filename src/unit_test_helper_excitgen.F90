module unit_test_helper_excitgen
  use constants
  use read_fci, only: readfciint, initfromfcid, fcidump_name
  use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
  use IntegralsData, only: UMat, umat_win
  use Integrals_neci, only: IntInit, get_umat_el_normal
  use procedure_pointers, only: get_umat_el, generate_excitation
  use SystemData, only: nel, nBasis, UMatEps, tStoreSpinOrbs, tReadFreeFormat, tCSF, &
       tReadInt
  use sort_mod
  use System, only: SysInit, SetSysDefaults
  use Parallel_neci, only: MPIInit, MPIEnd
  use UMatCache, only: GetUMatSize, tTransGTID
  use OneEInts, only: Tmat2D  
  use bit_rep_data, only: NIfTot, NIfDBO, NOffSgn, NIfSgn, extract_sign
  use bit_reps, only: encode_sign, decode_bit_det
  use DetBitOps, only: EncodeBitDet, DetBitEq
  use SymExcit3, only: countExcitations3, GenExcitations3
  use FciMCData, only: pSingles, pDoubles, pParallel, ilutRef, projEDet
  use SymExcitDataMod, only: excit_gen_store_type
  use Calc, only: CalcInit, SetCalcDefaults
  use dSFMT_interface, only: dSFMT_init, genrand_real2_dSFMT
  use Determinants, only: DetInit, DetPreFreezeInit, get_helement
  use util_mod, only: get_free_unit
  implicit none

  integer, parameter :: nelBase = 5
  integer, parameter :: nBasisBase = 12
  integer, parameter :: lmsBase = -1
  
contains

  subroutine test_excitation_generator(sampleSize, pTot, pNull)
    ! Test an excitation generator by creating a large number of excitations and
    ! compare the generated excits with a precomputed list of all excitations
    ! We thus make sure that
    !   a) all possible excitations are generated with some weight
    !   b) no invalid excitations are obtained
    implicit none
    integer, intent(in) :: sampleSize
    real(dp), intent(out) :: pTot, pNull
    integer :: nI(nel), nJ(nel)
    integer :: i, ex(2,2), exflag
    integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
    real(dp) :: pgen
    logical :: tPar, tAllExFound, tFound
    integer :: j, numEx, nSingles, nDoubles
    integer(n_int), allocatable :: allEx(:,:)
    real(dp) :: pgenArr(lenof_sign)
    real(dp) :: matel, matelN
    type(excit_gen_store_type) :: store
    logical :: exDoneDouble(0:nBasis,0:nBasis,0:nBasis,0:nBasis)
    logical :: exDoneSingle(0:nBasis,0:nBasis)
    integer :: ic, part, nFound, nullExcits
    HElement_t(dp) :: HEl
    exDoneDouble = .false.
    exDoneSingle = .false.

    ! some starting det - do NOT use the reference for the pcpp test, that would
    ! defeat the purpose
    do i = 1, nel
       if(2*i < nBasis) then
          nI(i) = 2*i-mod(i,2)
       else
          nI(i) = i
       endif
    end do
    call sort(nI)

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
    nullExcits = 0
    do i = 1, sampleSize
       store%tFilled = .false.
       call generate_excitation(nI,ilut,nJ,ilutJ,exFlag,ic,ex,tPar,pgen,HEl,store,part)
       ! lookup the excitation
       tFound = .false.
       do j = 1, numEx
          if(DetBitEQ(ilutJ,allEx(:,j))) then
             pgenArr = pgen
             call encode_sign(allEx(:,j),pgenArr)
             ! increase its counter by 1
             allEx(NIfTot+1,j) = allEx(NIfTot+1,j) + 1
             tFound = .true.
             exit
          end if
       end do
       ! if it was not found, and is not marked as invalid, we have a problem: this is not
       ! an excitaion
       if(.not. tFound .and. .not. nJ(1)==0) then
          call decode_bit_det(nJ,ilutJ)
          write(iout,*) "Error: Invalid excitation", nJ
          stop
       endif
       ! check if the generated excitation is invalid, if it is, mark this specific constellation
       ! so we do not double-count when calculating pNull
       if(nJ(1) == 0) then
          nullExcits = nullExcits + 1
          if(ic.eq.2) then
             if(.not. exDoneDouble(ex(1,1),ex(1,2),ex(2,1),ex(2,2)))then
                exDoneDouble(ex(1,1),ex(1,2),ex(2,1),ex(2,2)) = .true.
                pNull = pNull + pgen
             endif
          else if(ic.eq.1) then
             if(.not. exDoneSingle(ex(1,1),ex(1,2))) then
                exDoneSingle(ex(1,1),ex(1,2)) = .true.
                pNull = pNull + pgen
             endif
          endif
       endif
    end do

    ! check that all excits have been generated and all pGens are right
    ! probability normalization
    pTot = pNull
    matelN = 0.0
    do i = 1, numEx
       call decode_bit_det(nJ,allEx(:,i))
       matelN = matelN + abs(get_helement(nI,nJ))
    end do
    nFound = 0
    do i = 1, numEx
       call extract_sign(allEx(:,i),pgenArr)
       call decode_bit_det(nJ,allEx(:,i))
       if(pgenArr(1) > eps) then
          nFound = nFound + 1
          matel = get_helement(nI,nJ)
          write(iout,*) i, pgenArr(1), real(allEx(NIfTot+1,i))/real(sampleSize), &
               abs(matel)/(pgenArr(1)*matelN)
       endif
       pTot = pTot + pgenArr(1)
    end do
    write(iout,*) "Total prob. ", pTot
    write(iout,*) "pNull ", pNull
    write(iout,*) "Null ratio", nullExcits / real(sampleSize)
    write(iout,*) "In total", numEx, "excitations"
    write(iout,*) "Found", nFound, "excitations"

  end subroutine test_excitation_generator

  !------------------------------------------------------------------------------------------!

  subroutine init_excitgen_test()
    ! mimick the initialization of an FCIQMC calculation to the point where we can generate
    ! excitations with a weighted excitgen
    ! This requires setup of the basis, the symmetries and the integrals
    integer :: nBasisMax(5,3), lms
    integer(int64) :: umatsize
    real(dp) :: ecore
    umatsize = 0
    nel = nelBase
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

    call generate_random_integrals()    

    get_umat_el => get_umat_el_normal
    
    call initfromfcid(nel,nbasismax,nBasis,lms,.false.)
    lms = lmsBase
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
  end subroutine init_excitgen_test
  
  !------------------------------------------------------------------------------------------!

  subroutine finalize_excitgen_test()
    deallocate(TMat2D)
    call shared_deallocate_mpi(umat_win, UMat)
    call MPIEnd(.false.)
  end subroutine finalize_excitgen_test
  
  !------------------------------------------------------------------------------------------!

  ! generate an FCIDUMP file with random numbers with a given sparsity
  subroutine generate_random_integrals()
    real(dp), parameter :: sparse = 0.9
    real(dp), parameter :: sparseT = 0.1    
    integer :: i,j,k,l, iunit
    real(dp) :: r

    ! write the canonical FCIDUMP header
    iunit = get_free_unit()
    open(iunit,file="FCIDUMP")
    write(iunit,*) "&FCI NORB=",nBasisBase,",NELEC=",nelBase,"MS2=",lmsBase,","
    write(iunit,"(A)",advance="no") "ORBSYM="
    do i = 1, nBasisBase
       write(iunit,"(A)",advance="no") "1,"
    end do
    write(iunit,*)
    write(iunit,*) "ISYM=1,"
    write(iunit,*) "&END"
    ! generate random 4-index integrals with a given sparsity
    do i = 1, nBasisBase
       do j = 1, i
          do k = i, nBasisBase
             do l = 1, k
                r = genrand_real2_dSFMT()
                if(r < sparse) then
                   write(iunit, *) r*r, i,j,k,l
                endif
             end do
          end do
       end do
    end do
    ! generate random 2-index integrals with a given sparsity
    do i = 1, nBasisBase
       do j = 1, i
          r = genrand_real2_dSFMT()
          if(r < sparseT) then
             write(iunit,*) r, i,j
          endif
       end do
    end do

  end subroutine generate_random_integrals

  !------------------------------------------------------------------------------------------!  

  ! set the reference to the determinant with the first nel orbitals occupied
  subroutine set_ref()
    integer :: i
    allocate(projEDet(nel,1))
    do i = 1, nel
       projEDet(i,1) = i
    end do
    allocate(ilutRef(0:NifTot,1))
    call encodeBitDet(projEDet(:,1),ilutRef(:,1))
  end subroutine set_ref

  subroutine free_ref()
    deallocate(ilutRef)
    deallocate(projEDet)
  end subroutine free_ref
  
end module unit_test_helper_excitgen
