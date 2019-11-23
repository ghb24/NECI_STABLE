module LMat_calc
#ifdef USE_HDF5_
  use hdf5
#endif
  use hdf5_util
  use tc_three_body_data
  use LMat_Indexing, only: lMatIndSym, lMatIndSpin
  implicit none

  double precision, allocatable :: qwprod(:,:,:), ycoulomb(:,:,:,:), ycoulombAB(:,:,:,:)
#ifdef USE_HDF5_
  integer(hsize_t) :: nBasis, nGrid
#else
  integer :: nGrid
#endif
  logical :: ycoulombAB_exists
  integer(int64), allocatable :: lMatCalcHKeys(:)
  double precision, allocatable :: lMatCalcHVals(:)
  integer(int64) :: lMatIndMax
  integer(int64), allocatable :: lMatABCalcHKeys(:)
  double precision, allocatable :: lMatABCalcHVals(:)
  integer(int64) :: lMatABIndMax
  contains

  subroutine readLMatFactors(filename)

    character(*), intent(in) :: filename
    character(*), parameter :: nm_grp = "tcfactors", nm_nBasis = "nBasis", nm_nGrid = "nGrid", &
                               nm_weights="weights", nm_mo_vals="mo_vals", nm_ycoulomb="ycoulomb", nm_ycoulombAB="ycoulombAB"
    double precision, allocatable :: mo_vals(:,:), weights(:)
#ifdef USE_HDF5_
    integer(hid_t) :: err, file_id, grp_id, dataset, type_id
    integer(hsize_t) :: weights_dims(1), mo_vals_dims(2), ycoulomb_dims(4)
#endif
    integer i, a, b
    character(*), parameter :: this_routine = "readLMatFactors"
#ifdef USE_HDF5_
    write(iout,*) "**************** Reading TcFactors File ****************"
    call h5open_f(err)

    ! open the file
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err)

    call h5gopen_f(file_id, nm_grp, grp_id, err)

    call read_int64_attribute(grp_id, nm_nBasis, nBasis, required=.true.)
    write(iout,*) "Number of basis function: ", nBasis

    call read_int64_attribute(grp_id, nm_nGrid, nGrid, required=.true.)
    write(iout,*) "Number of grid points: ", nGrid

    ! load weights
    allocate(weights(nGrid))
    weights_dims = size(weights)
    call h5dopen_f(grp_id, nm_weights, dataset, err)
    call h5dget_type_f(dataset, type_id, err)
    call h5dread_f(dataset, type_id, weights, weights_dims, err)
    call h5tclose_f(type_id, err)
    call h5dclose_f(dataset, err)

    ! load mo_vals
    allocate(mo_vals(nGrid,nBasis))
    mo_vals_dims = size(mo_vals)
    call h5dopen_f(grp_id, nm_mo_vals, dataset, err)
    call h5dget_type_f(dataset, type_id, err)
    call h5dread_f(dataset, type_id, mo_vals, mo_vals_dims, err)
    call h5tclose_f(type_id, err)
    call h5dclose_f(dataset, err)

    ! load ycoulomb
    allocate(ycoulomb(3, nGrid, nBasis, nBasis))
    ycoulomb_dims = size(ycoulomb)
    call h5dopen_f(grp_id, nm_ycoulomb, dataset, err)
    call h5dget_type_f(dataset, type_id, err)
    call h5dread_f(dataset, type_id, ycoulomb, ycoulomb_dims, err)
    call h5tclose_f(type_id, err)
    call h5dclose_f(dataset, err)

    call h5lexists_f(grp_id, nm_ycoulombAB, ycoulombAB_exists, err)
    if (ycoulombAB_exists) then
        ! load ycoulombAB
        allocate(ycoulombAB(3, nGrid, nBasis, nBasis))
        call h5dopen_f(grp_id, nm_ycoulombAB, dataset, err)
        call h5dget_type_f(dataset, type_id, err)
        call h5dread_f(dataset, type_id, ycoulombAB, ycoulomb_dims, err)
        call h5tclose_f(type_id, err)
        call h5dclose_f(dataset, err)
    end if

    ! close the file, finalize hdf5
    call h5gclose_f(grp_id, err)
    call h5fclose_f(file_id, err)
    call h5close_f(err)

    if(.not. ycoulombAB_exists .and. tSpinCorrelator) then
        call stop_all(this_routine, "ycoulombAB is needed for Spin-Correlator but not found in TcFactors file.")
    end if

    !Combine weights and molecular oribtals
    allocate(qwprod(nGrid, nBasis, nBasis))
    do a=1,nBasis
        do b=1,nBasis
            qwprod(:,a, b) = mo_vals(:, a) * mo_vals(:, b) * weights
        end do
    end do
    deallocate(mo_vals)
    deallocate(weights)

    !ycoulom and ycoulom2 are symmetric in MO. However, only the upper part is filled
    !in the script producing the HDF5 file. Therefore, we fill the lower part here
    do a=1,nBasis
        do b=a+1,nBasis
            ycoulomb(:,:,a,b) = ycoulomb(:,:,b,a)
        end do
    end do

    if (ycoulombAB_exists) then
        do a=1,nBasis
            do b=a+1,nBasis
                ycoulombAB(:,:,a,b) = ycoulombAB(:,:,b,a)
            end do
        end do
    end if

    lMatIndMax = lMatIndSym(nBasis,nBasis,nBasis,nBasis,nBasis,nBasis)
    lMatCalcHSize = lMatCalcHFactor * lMatIndMax

    write(iout, *) "Total Size of LMat: ", lMatIndMax
    write(iout, *) "Size of LMatCalc Hash Table: ", lMatCalcHSize

    allocate(lMatCalcHKeys(lMatCalcHSize))
    allocate(lMatCalcHVals(lMatCalcHSize))
    do i=1,lMatCalcHSize
        lMatCalcHKeys(i) = -1
    end do
    lMatCalcHit = 0
    lMatCalcTot = 0
    lMatCalcHUsed = 0

    if (ycoulombAB_exists) then
        lMatABIndMax = lMatIndSpin(nBasis,nBasis,nBasis,nBasis,nBasis,nBasis)
        lMatABCalcHSize = lMatABCalcHFactor * lMatABIndMax

        write(iout, *) "Total Size of LMatAB: ", lMatABIndMax
        write(iout, *) "Size of LMatABCalc Hash Table: ", lMatABCalcHSize

        allocate(lMatABCalcHKeys(lMatABCalcHSize))
        allocate(lMatABCalcHVals(lMatABCalcHSize))
        do i=1,lMatABCalcHSize
            lMatABCalcHKeys(i) = -1
        end do
        lMatABCalcHit = 0
        lMatABCalcTot = 0
        lMatABCalcHUsed = 0
    endif

    write(iout,*) "********************************************************"
#else
    call stop_all(this_routine, 'HDF5 support not enabled at compile time')
#endif
  end subroutine readLMatFactors


  subroutine freeLMatFactors()

    deallocate(ycoulomb)
    deallocate(qwprod)
    deallocate(lMatCalcHKeys)
    deallocate(lMatCalcHVals)

    if (ycoulombAB_exists) then
        deallocate(ycoulombAB)
        deallocate(lMatABCalcHKeys)
        deallocate(lMatABCalcHVals)
    end if
  end subroutine freeLMatFactors


  function lMatCalc(i,k,m,j,l,n) result (matel)
    integer(int64), intent(in) :: i,j,k,l,m,n
    HElement_t(dp) :: matel
    integer(int64) :: hashKey, hashInd
    integer :: ii
    character(*), parameter :: this_routine = "lMatCalc"

    lMatCalcTot = lMatCalcTot + 1

    !First look up the value in the hash table
    hashKey = lMatIndSym(i,k,m,j,l,n)
    hashInd = mod(hashKey-1,lMatCalcHSize)+1 !Try whether other hash functions could imporve the hit rate

    if(hashKey==LMatCalcHKeys(hashInd))then
        lMatCalcHit = lMatCalcHit + 1
        matel = LMatCalcHVals(hashInd)
        return
    end if

    !It does not exist. So let's calculate it

    matel = 0.0_dp
    !Manual Optimization:
    !1-loop over first index is unrolled.
    !2-indpendent summations are done in seperate loops to enhance CPU-cahce utilization
    do ii = 1,nGrid
      matel = matel - qwprod(ii,m,n)*(&
            ycoulomb(1,ii,i,j)*ycoulomb(1,ii,k,l)+&
            ycoulomb(2,ii,i,j)*ycoulomb(2,ii,k,l)+&
            ycoulomb(3,ii,i,j)*ycoulomb(3,ii,k,l))
    end do
    do ii = 1,nGrid
      matel = matel - qwprod(ii,k,l)*(&
            ycoulomb(1,ii,i,j)*ycoulomb(1,ii,m,n)+&
            ycoulomb(2,ii,i,j)*ycoulomb(2,ii,m,n)+&
            ycoulomb(3,ii,i,j)*ycoulomb(3,ii,m,n))
    end do
    do ii = 1,nGrid
      matel = matel - qwprod(ii,i,j)*(&
            ycoulomb(1,ii,k,l)*ycoulomb(1,ii,m,n)+&
            ycoulomb(2,ii,k,l)*ycoulomb(2,ii,m,n)+&
            ycoulomb(3,ii,k,l)*ycoulomb(3,ii,m,n))
    end do

    !Store the new value in the hash table
    if(LMatCalcHKeys(hashInd)==-1)then
        LMatCalcHUsed = LMatCalcHUsed + 1
    end if
    LMatCalcHKeys(hashInd) = hashKey
    LMatCalcHVals(hashInd) = matel

  end function

  function lMatABCalc(i,k,m,j,l,n, spinMixture) result (matel)
    integer(int64), intent(in) :: i,j,k,l,m,n
    integer, intent(in) :: spinMixture
    HElement_t(dp) :: matel
    integer(int64) :: hashKey, hashInd
    integer :: ii
    character(*), parameter :: this_routine = "lMatCalc"

    lMatABCalcTot = lMatABCalcTot + 1

    !First look up the value in the hash table
    select case(spinMixture)
      case(1) ! ij has different spin
            hashKey = lMatIndSpin(i,k,m,j,l,n)
      case(2) !kl has different spin
            hashKey = lMatIndSpin(k,i,m,l,j,n)
      case(3) !mn has different spin
            hashKey = lMatIndSpin(m,k,i,n,l,j)
    end select
    hashInd = mod(hashKey-1,lMatABCalcHSize)+1 !Try whether other hash functions could imporve the hit rate

    if(hashKey==LMatABCalcHKeys(hashInd))then
        lMatABCalcHit = lMatABCalcHit + 1
        matel = LMatABCalcHVals(hashInd)
        return
    end if

    !It does not exist. So let's calculate it
    matel = 0.0_dp
    select case(spinMixture)
      case(1) ! ij has different spin
        do ii = 1,nGrid
          matel = matel - qwprod(ii,m,n)*(&
                ycoulombAB(1,ii,i,j)*ycoulomb(1,ii,k,l)+&
                ycoulombAB(2,ii,i,j)*ycoulomb(2,ii,k,l)+&
                ycoulombAB(3,ii,i,j)*ycoulomb(3,ii,k,l))
        end do
        do ii = 1,nGrid
          matel = matel - qwprod(ii,k,l)*(&
                ycoulombAB(1,ii,i,j)*ycoulomb(1,ii,m,n)+&
                ycoulombAB(2,ii,i,j)*ycoulomb(2,ii,m,n)+&
                ycoulombAB(3,ii,i,j)*ycoulomb(3,ii,m,n))
        end do
        do ii = 1,nGrid
          matel = matel - qwprod(ii,i,j)*(&
                ycoulombAB(1,ii,k,l)*ycoulombAB(1,ii,m,n)+&
                ycoulombAB(2,ii,k,l)*ycoulombAB(2,ii,m,n)+&
                ycoulombAB(3,ii,k,l)*ycoulombAB(3,ii,m,n))
        end do

      case(2) !kl has different spin
        do ii = 1,nGrid
          matel = matel - qwprod(ii,m,n)*(&
                ycoulomb(1,ii,i,j)*ycoulombAB(1,ii,k,l)+&
                ycoulomb(2,ii,i,j)*ycoulombAB(2,ii,k,l)+&
                ycoulomb(3,ii,i,j)*ycoulombAB(3,ii,k,l))
        end do
        do ii = 1,nGrid
          matel = matel - qwprod(ii,k,l)*(&
                ycoulombAB(1,ii,i,j)*ycoulombAB(1,ii,m,n)+&
                ycoulombAB(2,ii,i,j)*ycoulombAB(2,ii,m,n)+&
                ycoulombAB(3,ii,i,j)*ycoulombAB(3,ii,m,n))
        end do
        do ii = 1,nGrid
          matel = matel - qwprod(ii,i,j)*(&
                ycoulombAB(1,ii,k,l)*ycoulomb(1,ii,m,n)+&
                ycoulombAB(2,ii,k,l)*ycoulomb(2,ii,m,n)+&
                ycoulombAB(3,ii,k,l)*ycoulomb(3,ii,m,n))
        end do

      case(3) !mn has different spin
        do ii = 1,nGrid
          matel = matel - qwprod(ii,m,n)*(&
                ycoulombAB(1,ii,i,j)*ycoulombAB(1,ii,k,l)+&
                ycoulombAB(2,ii,i,j)*ycoulombAB(2,ii,k,l)+&
                ycoulombAB(3,ii,i,j)*ycoulombAB(3,ii,k,l))
        end do
        do ii = 1,nGrid
          matel = matel - qwprod(ii,k,l)*(&
                ycoulomb(1,ii,i,j)*ycoulombAB(1,ii,m,n)+&
                ycoulomb(2,ii,i,j)*ycoulombAB(2,ii,m,n)+&
                ycoulomb(3,ii,i,j)*ycoulombAB(3,ii,m,n))
        end do
        do ii = 1,nGrid
          matel = matel - qwprod(ii,i,j)*(&
                ycoulomb(1,ii,k,l)*ycoulombAB(1,ii,m,n)+&
                ycoulomb(2,ii,k,l)*ycoulombAB(2,ii,m,n)+&
                ycoulomb(3,ii,k,l)*ycoulombAB(3,ii,m,n))
        end do
      end select

    !Store the new value in the hash table
    if(LMatABCalcHKeys(hashInd)==-1)then
        LMatABCalcHUsed = LMatABCalcHUsed + 1
    end if
    LMatABCalcHKeys(hashInd) = hashKey
    LMatABCalcHVals(hashInd) = matel

  end function
end module LMat_calc
