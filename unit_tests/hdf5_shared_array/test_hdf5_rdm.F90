program test_hdf5_rdm
    use fruit
    use constants, only: dp
    use MPI_wrapper
    use Parallel_neci
#ifdef USE_HDF_
    use parallel_hdf5_utils, only: write_data_phdf5, read_data_phdf5
    use hdf5
#endif
    implicit none

    integer :: failed_count

    call MPIInit(.false.)
    call init_fruit()
#ifdef USE_HDF_
    call test_hdf5_ops()
#endif
    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)
    if (failed_count /= 0) call stop_all('test_hdf5_rdm', 'failed_tests')
    call MPIEnd(.false.)

contains

#ifdef USE_HDF_
    subroutine test_hdf5_ops()
      real(dp) :: rands(6,iProcIndex+1)
      integer :: wrt_buf(6,iProcIndex+1), rec_buf(6,(nProcessors*(nProcessors+1))/2)
      integer(hid_t) :: plist_id, file_id, grp_id
      integer(hdf_err) :: err
      integer :: mpi_err, len_arrays(nProcessors), low_bound, high_bound

      call h5open_f(err)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
      call h5pset_fapl_mpio_f(plist_id, CommGlobal, mpiInfoNull, err)
      call h5fcreate_f('unittest.h5', H5F_ACC_TRUNC_F, file_id, err, access_prp=plist_id)
      call h5pclose_f(plist_id, err)
      call h5gcreate_f(file_id, 'unit_test', grp_id, err)
      call MPIAllGather(size(wrt_buf, dim=2), len_arrays, mpi_err)

      call random_number(rands)
      wrt_buf(:,:) = floor(11 * rands)  ! random integers [0,10]

      call write_data_phdf5(wrt_buf, 'testdata', grp_id)
      call read_data_phdf5(rec_buf, 'testdata', grp_id)

      low_bound = sum(len_arrays(1:iProcIndex+1)) - size(wrt_buf, dim=2) + 1
      high_bound = sum(len_arrays(1:iProcIndex+1))
      call assert_true(all(rec_buf(:, low_bound:high_bound) == wrt_buf(:,:)))

      call h5gclose_f(grp_id, err)
      call h5fclose_f(file_id, err)
      call h5close_f(err)
      call h5garbage_collect_f(err)
    end subroutine test_hdf5_ops
#endif

end program
