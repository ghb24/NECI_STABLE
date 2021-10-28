#include "macros.h"

module LMat_calc
#ifdef USE_HDF5_
    use hdf5
#endif
    use hdf5_util
    use tc_three_body_data
    use LMat_Indexing, only: lMatIndSym, lMatIndSpin
    use util_mod, only: get_free_unit
    use mpi_wrapper, only: iProcIndex
    implicit none

    real(dp), allocatable :: qwprod(:, :, :), ycoulomb(:, :, :, :)
#ifdef USE_HDF5_
    integer(hsize_t) :: nBasis, nGrid
#else
    integer :: nGrid
#endif
    integer(int64), allocatable :: lMatCalcHKeys(:)
    real(dp), allocatable :: lMatCalcHVals(:)
    integer(int64) :: lMatIndMax
contains

    subroutine readLMatFactors()

        character(*), parameter :: filename = "tcfactors.h5"
        character(*), parameter :: nm_grp = "tcfactors", nm_nBasis = "nBasis", nm_nGrid = "nGrid", &
                                   nm_weights = "weights", nm_mo_vals = "mo_vals", nm_ycoulomb = "ycoulomb"
        real(dp), allocatable :: mo_vals(:, :), weights(:)
#ifdef USE_HDF5_
        integer(hid_t) :: err, file_id, grp_id, dataset, type_id
        integer(hsize_t) :: weights_dims(1), mo_vals_dims(2), ycoulomb_dims(4)
#endif
        integer i, a, b
        character(*), parameter :: this_routine = "readLMatFactors"
#ifdef USE_HDF5_
        write(stdout, *) "**************** Reading TcFactors File ****************"
        call h5open_f(err)

        ! open the file
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err)

        call h5gopen_f(file_id, nm_grp, grp_id, err)

        call read_int64_attribute(grp_id, nm_nBasis, nBasis, required=.true.)
        write(stdout, *) "Number of basis function: ", nBasis

        call read_int64_attribute(grp_id, nm_nGrid, nGrid, required=.true.)
        write(stdout, *) "Number of grid points: ", nGrid

        ! load weights
        allocate(weights(nGrid))
        weights_dims = size(weights)
        call h5dopen_f(grp_id, nm_weights, dataset, err)
        call h5dget_type_f(dataset, type_id, err)
        call h5dread_f(dataset, type_id, weights, weights_dims, err)
        call h5tclose_f(type_id, err)
        call h5dclose_f(dataset, err)

        ! load mo_vals
        allocate(mo_vals(nGrid, nBasis))
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

        ! close the file, finalize hdf5
        call h5gclose_f(grp_id, err)
        call h5fclose_f(file_id, err)
        call h5close_f(err)

        !Combine weights and molecular oribtals
        allocate(qwprod(nGrid, nBasis, nBasis))
        do a = 1, nBasis
            do b = 1, nBasis
                qwprod(:, a, b) = mo_vals(:, a) * mo_vals(:, b) * weights
            end do
        end do
        deallocate(mo_vals)
        deallocate(weights)

        !ycoulom and ycoulom2 are symmetric in MO. However, only the upper part is filled
        !in the script producing the HDF5 file. Therefore, we fill the lower part here
        do a = 1, nBasis
            do b = a + 1, nBasis
                ycoulomb(:, :, a, b) = ycoulomb(:, :, b, a)
            end do
        end do

        lMatIndMax = lMatIndSym(nBasis, nBasis, nBasis, nBasis, nBasis, nBasis)
        lMatCalcHSize = lMatCalcHFactor * lMatIndMax

        write(stdout, *) "Total Size of LMat: ", lMatIndMax
        write(stdout, *) "Size of LMatCalc Hash Table: ", lMatCalcHSize

        allocate(lMatCalcHKeys(lMatCalcHSize))
        allocate(lMatCalcHVals(lMatCalcHSize))
        do i = 1, lMatCalcHSize
            lMatCalcHKeys(i) = -1
        end do
        lMatCalcHit = 0
        lMatCalcTot = 0
        lMatCalcHUsed = 0

        write(stdout, *) "********************************************************"
#else
        call stop_all(this_routine, 'HDF5 support not enabled at compile time')
        unused_var(filename)
#endif
    end subroutine readLMatFactors

    subroutine freeLMatFactors()

        if (allocated(ycoulomb))      deallocate(ycoulomb)
        if (allocated(qwprod))        deallocate(qwprod)
        if (allocated(lMatCalcHKeys)) deallocate(lMatCalcHKeys)
        if (allocated(lMatCalcHVals)) deallocate(lMatCalcHVals)
    end subroutine freeLMatFactors

    function lMatCalc(i, k, m, j, l, n) result(matel)
        integer(int64), intent(in) :: i, j, k, l, m, n
        HElement_t(dp) :: matel
        integer(int64) :: hashKey, hashInd
        integer :: ii
        character(*), parameter :: this_routine = "lMatCalc"

        lMatCalcTot = lMatCalcTot + 1

        !First look up the value in the hash table
        hashKey = lMatIndSym(i, k, m, j, l, n)
        hashInd = mod(hashKey - 1, lMatCalcHSize) + 1 !Try whether other hash functions could imporve the hit rate

        if (hashKey == LMatCalcHKeys(hashInd)) then
            lMatCalcHit = lMatCalcHit + 1
            matel = LMatCalcHVals(hashInd)
            return
        end if

        !It does not exist. So let's calculate it

        matel = 0.0_dp
        !Manual Optimization:
        !1-loop over first index is unrolled.
        !2-indpendent summations are done in seperate loops to enhance CPU-cahce utilization
        do ii = 1, nGrid
            matel = matel - qwprod(ii, m, n) * ( &
                    ycoulomb(1, ii, i, j) * ycoulomb(1, ii, k, l) + &
                    ycoulomb(2, ii, i, j) * ycoulomb(2, ii, k, l) + &
                    ycoulomb(3, ii, i, j) * ycoulomb(3, ii, k, l))
        end do
        do ii = 1, nGrid
            matel = matel - qwprod(ii, k, l) * ( &
                    ycoulomb(1, ii, i, j) * ycoulomb(1, ii, m, n) + &
                    ycoulomb(2, ii, i, j) * ycoulomb(2, ii, m, n) + &
                    ycoulomb(3, ii, i, j) * ycoulomb(3, ii, m, n))
        end do
        do ii = 1, nGrid
            matel = matel - qwprod(ii, i, j) * ( &
                    ycoulomb(1, ii, k, l) * ycoulomb(1, ii, m, n) + &
                    ycoulomb(2, ii, k, l) * ycoulomb(2, ii, m, n) + &
                    ycoulomb(3, ii, k, l) * ycoulomb(3, ii, m, n))
        end do

        !Store the new value in the hash table
        if (LMatCalcHKeys(hashInd) == -1) then
            LMatCalcHUsed = LMatCalcHUsed + 1
        end if
        LMatCalcHKeys(hashInd) = hashKey
        LMatCalcHVals(hashInd) = matel

    end function

    subroutine read_rs_lmat_factors
        character(*), parameter :: filename_mo = "mos_in_r"
        character(*), parameter :: filename_ints = "x_w_ij_r"
        character(*), parameter :: filename_info = 'num_mos_grid'
        character(*), parameter :: this_routine = "read_rs_lmat_factors"

        integer :: iunit, ierr, i, j, k, l
        integer(int64) :: num_mos
        real(dp) :: integral
        real(dp), allocatable :: array_mos(:,:)

        root_print "Reading in range-separated TC factors"

        iunit = get_free_unit()
        ! Read the number of grid points and MOs
        open(iunit, file=filename_info, status='old')
        read(iunit, *, iostat=ierr) num_mos, ngrid
        close(iunit)

        root_print "with: ", ngrid, " grid points"
        root_print "and ", num_mos, " orbitals"

        ! read the MOs file
        allocate(array_mos(ngrid, num_mos), source=0.0_dp)
        iunit = get_free_unit()
        open(iunit, file=filename_mo, status='old')
        do
            read(iunit, *, iostat=ierr) i, j, array_mos(i, j)
            if (ierr < 0) then
                exit
            else if (ierr > 0) then
                call stop_all(this_routine, "error reading " // filename_mo)
            end if
        end do
        close(iunit)

        ! fill in the qwprod file:
        allocate(qwprod(ngrid, num_mos, num_mos), source=0.0_dp)
        do i = 1, num_mos
            do j = 1, num_mos
                qwprod(:, i, j) = array_mos(:, i) * array_mos(:, j)
            end do
        end do
        deallocate(array_mos)

        allocate(ycoulomb(3, ngrid, num_mos, num_mos), source=0.0_dp)
        iunit = get_free_unit()
        open(iunit, file=filename_ints, status='old')
        do
            read(iunit, *, iostat=ierr) i, j, k, l, integral
            if (ierr < 0) then
                exit
            else if (ierr > 0) then
                call stop_all(this_routine, "error reading " // filename_ints)
            end if
            ycoulomb(j, i, k, l) = integral
        end do
        close(iunit)


        lMatIndMax = lMatIndSym(num_mos, num_mos, num_mos, num_mos, num_mos, num_mos)
        lMatCalcHSize = lMatCalcHFactor * lMatIndMax

        root_print "Total Size of LMat: ", lMatIndMax
        root_print "Size of LMatCalc Hash Table: ", lMatCalcHSize

        allocate(lMatCalcHKeys(lMatCalcHSize))
        allocate(lMatCalcHVals(lMatCalcHSize))
        do i = 1, lMatCalcHSize
            lMatCalcHKeys(i) = -1
        end do
        lMatCalcHit = 0
        lMatCalcTot = 0
        lMatCalcHUsed = 0

        root_print "Done: reading in range-separated TC factors"

    end subroutine read_rs_lmat_factors

!     function rs_lmat_calc(i, j, m, k, l, n) result(matel)
!         integer(int64), intent(in) :: i, j, k, l, m, n
!         HElement_t(dp) :: matel
!         character(*), parameter :: this_routine = "rs_lmat_calc"
!
!         integer :: g, x
!
!         matel = 0.0_dp
!
!         do x = 1, 3
!             do g = 1, ngrid
!                 matel = matel + array_mos(g, i) * array_mos(g, k) &
!                     * array_ints(g, x, m, n) * array_ints(g, x, j, l)
!                 matel = matel + array_mos(g, j) * array_mos(g, l) &
!                     * array_ints(g, x, m, n) * array_ints(g, x, i, k)
!                 matel = matel + array_mos(g, m) * array_mos(g, n) &
!                     * array_ints(g, x, j, l) * array_ints(g, x, i, k)
!             end do
!         end do
!
!     end function rs_lmat_calc

end module LMat_calc
