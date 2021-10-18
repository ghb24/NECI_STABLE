! Unit test for gdata_io handlers
program test_gdata_io
    use constants
    use fruit
    use FciMCData, only: MaxWalkersPart
    use LoggingData, only: tAccumPops, tExplicitAllRDM, tRDMonFly
    use CalcData, only: tActivateLAS, tAutoAdaptiveShift, tContTimeFCIMC, tContTimeFull, &
        tReplicaEstimates, tPairedReplicas, tScaleBlooms, tSeniorInitiators, tStoredDets
    use gdata_io, only: gdata_io_t
    use global_det_data
    implicit none

    integer, parameter :: ndets = 1e4
    integer, parameter :: rd = 234
    real(dp), parameter :: acc_val = 1.234_dp
    real(dp), parameter :: tot_val = 2.48
    real(dp), parameter :: mr_val = 5.678_dp

    call init_fruit()
    call gdata_io_test_driver()
    call fruit_summary()
    call fruit_finalize()

contains

    !------------------------------------------------------------------------------------------!

    subroutine gdata_io_test_driver()
        call init_gdata_io_test()
        call test_sizes()
        call test_read()
        call test_write()
#ifdef USE_HDF_
        call test_write_hdf5()
#endif
    end subroutine gdata_io_test_driver

    !------------------------------------------------------------------------------------------!

    subroutine init_gdata_io_test()
        ! For (k)mneci and dneci, we need to set the number of runs
#if defined(PROG_NUMRUNS_) || defined(DOUBLERUN_)
        inum_runs = 4
#ifdef CMPLX_
        lenof_sign = 2*inum_runs
#else
        lenof_sign = inum_runs
#endif
#endif
        ! set all flags referred to in global_det_data initialization
        tSeniorInitiators = .false.
        tAutoAdaptiveShift = .true.
        tRDMonFly = .false.
        tExplicitAllRDM = .false.
        tContTimeFCIMC = .false.
        tContTimeFull = .false.
        tStoredDets = .false.
        tActivateLAS = .false.
        tScaleBlooms = .true.
        tPairedReplicas = .true.
        tReplicaEstimates = .false.
        MaxWalkersPart = ndets
        tAccumPops = .true.
        ! initialise the global det data array that we will read/write to
        call init_global_det_data(1,1)

    end subroutine init_gdata_io_test

    !------------------------------------------------------------------------------------------!

    subroutine test_sizes()
        ! Create a test object and check that its size matches the one dictated by input
        type(gdata_io_t) :: test_handler

        ! test the different options' sizes
        call test_handler%init_gdata_io(.true.,.true., .true., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == max_ratio_size + fvals_size + apvals_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.false.,.true., .true., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == max_ratio_size + apvals_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.true.,.false., .true., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == fvals_size + apvals_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.false.,.false., .true., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == apvals_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.true., .true., .false., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == max_ratio_size + fvals_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.false.,.true., .false., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == max_ratio_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.true.,.false., .false., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == fvals_size)
        call assert_true(test_handler%t_io())

        call test_handler%init_gdata_io(.false.,.false., .false., fvals_size, max_ratio_size, apvals_size)
        call assert_true(test_handler%entry_size() == 0)
        call assert_true(.not. test_handler%t_io())
    end subroutine test_sizes

    !------------------------------------------------------------------------------------------!

    subroutine test_read()
        ! Create a test object and read some data to global data
        type(gdata_io_t) :: test_handler
        real(dp), allocatable :: gdata_buf(:,:)

        call test_handler%init_gdata_io(.true., .true., .false., fvals_size, max_ratio_size, apvals_size)
        allocate(gdata_buf(test_handler%entry_size(),ndets))
        gdata_buf(1:inum_runs,:) = acc_val
        gdata_buf(inum_runs+1:2*inum_runs,:) = tot_val
        gdata_buf(fvals_size+1,:) = mr_val

        ! read the values from gdata_buf into global_det_data
        call test_handler%read_gdata(gdata_buf, ndets)

        call assert_equals(get_tot_spawns(rd,1), gdata_buf(1+inum_runs,rd))
        call assert_equals(get_acc_spawns(1,inum_runs), gdata_buf(inum_runs,1))
        call assert_equals(get_max_ratio(ndets), gdata_buf(fvals_size+1, ndets))

        deallocate(gdata_buf)
    end subroutine test_read

    !------------------------------------------------------------------------------------------!

    subroutine test_write()
        type(gdata_io_t) :: test_handler
        real(dp), allocatable :: gdata_buf(:,:)

        call set_global_det_data()

        call test_handler%init_gdata_io(.true., .true., .false., fvals_size, max_ratio_size, apvals_size)

        allocate(gdata_buf(test_handler%entry_size(),ndets))
        ! zero the buffer
        gdata_buf = 0.0_dp
        ! write the global_determinant_data to gdata_buf
        call test_handler%write_gdata(gdata_buf, ndets)

        ! the entries in the buffer have to match the original ones
        call assert_equals(get_tot_spawns(rd, 1), gdata_buf(1+inum_runs, rd))
        call assert_equals(get_acc_spawns(rd, 1), gdata_buf(1, rd))
        call assert_equals(get_max_ratio(rd), gdata_buf(fvals_size + 1, rd))

        deallocate(gdata_buf)
    end subroutine test_write

    !------------------------------------------------------------------------------------------!
#ifdef USE_HDF_
    subroutine test_write_hdf5()
        use hdf5, only: hsize_t
        type(gdata_io_t) :: test_handler
        integer(hsize_t), allocatable :: gdata_buf(:,:)

        call set_global_det_data()
        call test_handler%init_gdata_io(.true., .true., .false., fvals_size, max_ratio_size, apvals_size)

        allocate(gdata_buf(test_handler%entry_size(),ndets))
        gdata_buf = 0_hsize_t

        call test_handler%write_gdata_hdf5(gdata_buf, ndets)
        ! the entries in the buffer have to match the original ones
        call assert_equals(get_tot_spawns(rd, 1), transfer(gdata_buf(1, rd), acc_val))
        call assert_equals(get_acc_spawns(rd, 1), transfer(gdata_buf(inum_runs+1, rd), tot_val))
        call assert_equals(get_max_ratio(rd), transfer(gdata_buf(fvals_size + 1, rd), mr_val))
        deallocate(gdata_buf)
    end subroutine test_write_hdf5

#endif

    subroutine set_global_det_data()
        ! write some random values to some random part of global_determinant_data (otherwise 0)
        global_determinant_data = 0
        call set_max_ratio(mr_val, rd)
        global_determinant_data(pos_acc_spawns, rd) = acc_val
        global_determinant_data(pos_tot_spawns, rd) = tot_val
    end subroutine set_global_det_data

end program test_gdata_io
