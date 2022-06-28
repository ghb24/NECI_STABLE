#include "macros.h"

module rdm_hdf5

    use MPI_wrapper
    use Parallel_neci
    use constants
    use hdf5_util
    use util_mod
#ifdef USE_HDF_
    use hdf5
#endif
    use fortran_strings
    use rdm_data, only: rdm_definitions_t, rdm_list_t
    use guga_bitRepOps, only: extract_2_rdm_ind
    use rdm_data_utils, only: extract_sign_rdm, calc_separate_rdm_labels
    implicit none
    private
    public :: write_rdms_hdf5

    character(*), parameter :: nm_2rdm = '/archive/rdms/2200/'

contains

    subroutine write_rdms_hdf5(rdm_defs, rdm, rdm_trace)
        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        character(*), parameter :: t_r = 'write_hdf5_rdms'
#ifdef USE_HDF_
        integer(hid_t) :: plist_id, file_id
        integer(hdf_err) :: err
        integer :: mpi_err
        character(255) :: filename
        integer :: iroot

        do iroot = 1, rdm_defs%nrdms
            ! TODO: make filename unique and increment if necessary
            filename = 'fciqmc.rdms.' // str(iroot) // '.h5'

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

            call write_2rdm_hdf5(file_id, rdm_defs, rdm, rdm_trace, iroot)
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
        end do
#else
        call stop_all(t_r, 'HDF5 support not enabled at compile time')
#endif
    end subroutine write_rdms_hdf5


    subroutine write_2rdm_hdf5(parent, rdm_defs, rdm, rdm_trace, iroot)
        use SystemData, only: tGUGA
        integer(hid_t), intent(in) :: parent
        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        integer, intent(in) :: iroot
        integer(hid_t) :: twordm_grp_id
        integer(hdf_err) :: err
        integer(hsize_t) :: arr_nelements(0:nProcessors - 1), totlength, &
                            write_offset(2), n_rdm_el
        real(dp) :: rdm_sign(rdm%sign_length)
        integer :: ierr
        integer(n_int), allocatable :: indices(:,:)
        real(dp), allocatable :: values(:)
        integer :: iproc, p, q, r, s, pq, rs, ielem
        integer(n_int) :: pq_, rs_
        integer(int_rdm) :: pqrs

        call h5gcreate_f(parent, nm_2rdm, twordm_grp_id, err)
        ! rdm_nelements is a local copy of type hsize_t
        n_rdm_el = int(rdm%nelements, kind(hsize_t))
        call MPIAllGather(n_rdm_el, arr_nelements, ierr)
        totlength = sum(arr_nelements)
        write_offset = [0_hsize_t, sum(arr_nelements(0:iProcIndex - 1))]

        ! write these intermediate arrays to disk later
        ! TODO: lengths are wrong
        allocate(indices(4, rdm%nelements), source=0_n_int)
        allocate(values(rdm%nelements), source=0.0_dp)

        ! create the intermediate "indices" and "value" arrays
        do iproc = 0, (nProcessors - 1)
            if (iproc == iProcIndex) then  ! what does this mean?
                do ielem = 1, rdm%nelements
                    pqrs = rdm%elements(0, ielem)

                    call extract_sign_rdm(rdm%elements(:, ielem), rdm_sign)
                    ! Normalise.
                    rdm_sign = rdm_sign / rdm_trace
                    if (tGUGA) then
                        call extract_2_rdm_ind(pqrs, p, q, r, s, pq_, rs_)
                        pq = int(pq_); rs = int(rs_)
                    else
                        call calc_separate_rdm_labels(pqrs, pq, rs, p, s, q, r)
                    end if

                    if (abs(rdm_sign(iroot)) > 1.e-12_dp) then
                        if (p >= q .and. pq_ >= rs_ .and. p >= r .and. p >= s) then
                            indices(1, ielem) = p; indices(2, ielem) = q
                            indices(3, ielem) = r; indices(4, ielem) = s
                            values(ielem) = rdm_sign(iroot)
                        end if
                    end if
                end do
            end if
        end do

        ! write index arrays
        call write_2d_multi_arr_chunk_buff( &
            twordm_grp_id, &
            'indices', &
            h5kind_to_type(int64, H5_INTEGER_KIND), &
            indices, &  ! array to write to file
            [4_hsize_t, int(rdm%nelements, hsize_t)], &  ! 4 indices wide
            [0_hsize_t, 0_hsize_t], & ! offset
            [4_hsize_t, totlength], & ! all dims
            [0_hsize_t, sum(arr_nelements(0:iProcIndex - 1))] & ! output offset
        )
        ! write value arrays
        ! the function accepts only integer arrays
        call write_dp_1d_attribute(twordm_grp_id, 'values', values)

        call h5gclose_f(twordm_grp_id, err)
    end subroutine write_2rdm_hdf5

end module rdm_hdf5
