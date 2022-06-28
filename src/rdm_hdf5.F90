#include "macros.h"

module rdm_hdf5

    ! Read and write RDMs using hdf5 as the data format
    !
    ! /archive/
    !     /rdms/
    !         /1100/
    !         /2200/
    !         /3300/
    !         /4400f/
    !         /norm/

    use MPI_wrapper
    use Parallel_neci
    use constants
    use hdf5_util
    use util_mod
    use MemoryManager, only: LogMemAlloc, LogMemDeAlloc
#ifdef USE_HDF_
    use hdf5
    use gdata_io, only: gdata_io_t, clone_signs, resize_attribute
#endif
    implicit none
    private
    public :: write_hdf5_rdms, read_hdf5_rdms

    character(*), parameter :: &
        nm_1rdm =         '/archive/rdms/1100', &
        nm_1rdm_values =  '/archive/rdms/1100/values', &
        nm_1rdm_indices = '/archive/rdms/1100/indices', &
        nm_2rdm =         '/archive/rdms/2200/', &
        nm_2rdm_values =  '/archive/rdms/2200/values', &
        nm_2rdm_indices = '/archive/rdms/2200/indices'

contains

    subroutine write_rdms_hdf5()
        character(*), parameter :: t_r = 'write_hdf5_rdms'
#ifdef USE_HDF_
        integer(hid_t) :: plist_id, file_id
        integer(hdf_err) :: err
        integer :: mpi_err
        character(255) :: filename

        ! TODO: make filename unique and increment if necessary
        filename = 'fciqmc.rdms.h5'

        write(stdout, *) "============== Writing HDF5 RDMs =============="
        write(stdout, *) "File name: ", filename

        call h5open_f(err)
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
        call h5pset_fapl_mpio_f(plist_id, CommGlobal, mpiInfoNull, err)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp=plist_id)
        call h5pclose_f(plist_id, err)

        ! TODO: write metadata to file
        ! write(stdout, *) "writing metadata"
        ! call write_rdm_metadata(file_id)

        call write_twordm_hdf5(file_id)
        ! later
        ! call write_onerdm_hdf5(file_id)
        ! call write_threerdm_hdf5(file_id)
        ! call write_contracted_fock_hdf5(file_id)

        call MPIBarrier(mpi_err)

        write(stdout, *) "closing RDM file"
        call h5fclose_f(file_id, err)
        call h5close_f(err)

        call h5garbage_collect_f(err)

        call MPIBarrier(mpi_err)
        write(stdout, *) "RDM file write successful"
#else
        call stop_all(t_r, 'HDF5 support not enabled at compile time')
#endif
    end subroutine write_rdmfile_hdf5


    subroutine write_2rdm_hdf5(parent)
        use ..., only: rdm_defs, rdm_recv,
        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: twordm_grp_id
        integer(hdf_err) :: err
        integer(hsize_t) :: arr_rdm%nelements(0:nProcessors - 1), totlength, &
                            write_offset(2)

        integer, allocatable :: indices(:,:)
        real(dp), allocatable :: values(:,:)
        character(*), parameter :: this_routine = "write_2rdm_hdf5"

        call h5gcreate_f(parent, nm_twordm, twordm_grp_id, err)

        ! how long are the arrays on each rank
        call MPIAllGather(rdm%nelements, arr_nelements, ierr)
        totlength = sum(arr_nelements)
        write_offset = [0_hsize_t, sum(arr_nelements(0:iProcIndex -1))]
        ! rdm_trace does not need to be synced?

        ! create the intermediate "indices" and "value" arrays
        do iproc = 0, (nProcessors - 1)
            if (iproc == iProcIndex) then  ! what does this mean?
                do irdm = 1, rdm_defs%nrdms
                    do ielem = 1, rdm%nelements
                        pqrs = rdm%elements(0, ielem)
                        call extract_sign_rdm(rdm%elements(:, ielem), rdm_sign)
                        rdm_sign = rdm_sign / rdm_trace
                        call extract_2_rdm_ind(pqrs, p, q, r, s, pq_, rs_)
                        pq = int(pq_); rs = int(rs_)
                        if (abs(rdm_sign(irdm)) > 1.e-12_dp) then
                            if (p >= q .and. pq >= rs .and. p >= r .and. p >= s) then
                                indices(1, ielem) = p; indices(2, ielem) = q
                                indices(3, ielem) = r; indices(4, ielem) = s
                                values(1, ielem) = rdm_sign(irdm)
                            end if
                        end if
                    end do
                end do
            end if
        end do

        ! write index arrays
        call write_2d_multi_arr_chunk_buff( &
            wfn_grp_id, 'indices', h5kind_to_type(int64, H5_INTEGER_KIND), &
            indices, &  ! array to write to file
            [4_hsize_t, int(printed_count, hsize_t)], &  ! 4 indices wide
            [0_hsize_t, 0_hsize_t], & ! offset
            [4_hsize_t, arr_rdm%nelements], & ! all dims
            [0_hsize_t, sum(arr_rdm%nelements(0:iProcIndex - 1))] & ! output offset
            )

        ! write value arrays
        call write_2d_multi_arr_chunk_buff( &
            wfn_grp_id, 'values', h5kind_to_type(dp, H5_REAL_KIND), &
            values, &  ! array to write to file
            ! what does nifd mean?
            [1_hsize_t, int(printed_count, hsize_t)], &  ! 1 index wide
            [0_hsize_t, 0_hsize_t], & ! offset
            [1_hsize_t, arr_rdm%nelements], & ! all dims
            [0_hsize_t, sum(arr_rdm%nelements(0:iProcIndex - 1))] & ! output offset
            )

        call h5gclose(twordm_grp_id, err)
    end subroutine write_2rdm_hdf5

end module rdm_hdf5
