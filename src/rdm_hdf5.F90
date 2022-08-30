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

    implicit none
    private
    public :: write_rdms_hdf5


contains


    subroutine write_rdms_hdf5(rdm_defs, rdm, rdm_trace, one_rdms)
        use rdm_data, only: rdm_definitions_t, rdm_list_t, one_rdm_t
        use LoggingData, only: tWriteSpinFreeRDM, tPrint1RDM
        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        type(one_rdm_t), intent(inout), optional :: one_rdms(:)
        character(*), parameter :: this_routine = 'write_hdf5_rdms'
#ifdef USE_HDF_
        integer(hid_t) :: file_id, root_id, rdm_id, plist_id
        integer(hdf_err) :: err
        integer :: mpi_err
        character(255) :: filename
        integer :: iroot

        if (rdm_defs%nrdms_transition > 0) then
            call stop_all(this_routine, "Transition RDM support pending.")
        end if

        do iroot = 1, rdm_defs%nrdms
            if (iProcIndex == 0) then
                filename = 'fciqmc.rdms.' // str(iroot) // '.h5'
                write(stdout, *) "============== Writing HDF5 RDMs =============="
                write(stdout, *) "File name: ", trim(filename)
                write(stdout, *) "Regular RDMs saved in /archive/rdms/AA00/ where &
                                     &A denotes the number of fermion operators."
            end if
            call MPIBCast(filename)

            ! All HDF5 APIs that modify structural metadata are collective
            call h5open_f(err)
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
            ! setup file access property list for MPI-IO access.
            call h5pset_fapl_mpio_f(plist_id, CommGlobal, mpiInfoNull, err)
            ! create the file collectively
            ! H5F_ACC_TRUNC_F = overwrite existing file
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp=plist_id)
            call h5pclose_f(plist_id, err)
            call h5gcreate_f(file_id, 'archive', root_id, err)
            call h5gcreate_f(root_id, 'rdms', rdm_id, err)

            if (tPrint1RDM) then
                call write_1rdm_hdf5(rdm_id, one_rdms(iroot)%matrix)
                write(stdout, *) "1RDM written to file."
            end if
            if (tWriteSpinFreeRDM) then
                call write_2rdm_hdf5(rdm_id, rdm, rdm_trace, iroot)
                write(stdout, *) "2RDM written to file."
            end if
            ! if (tWrite3RDM) then call write_3rdm_hdf5(file_id)
            ! if (tWriteF4RDM) then call write_contracted_fock_hdf5(file_id)

            if (iProcIndex == 0) write(stdout, *) "closing RDM file."
            call h5gclose_f(rdm_id, err)
            call h5gclose_f(root_id, err)
            call h5fclose_f(file_id, err)
            call h5close_f(err)
            call h5garbage_collect_f(err)
            if (iProcIndex == 0) write(stdout, *) "RDM file write successful."
        end do
#else
        call stop_all(this_routine, 'HDF5 support not enabled at compile time.')
#endif
    end subroutine write_rdms_hdf5


    subroutine write_1rdm_hdf5(parent, one_rdm)
        use RotateOrbsData, only: ind => SymLabelListInv_rot
        integer(hid_t), intent(in) :: parent
        real(dp), intent(inout) :: one_rdm(:, :)
        integer(n_int), allocatable :: indices(:,:)
        real(dp), allocatable :: values(:)
        integer :: p, q, pq, dim_1rdm_vec
        integer(hid_t) :: onerdm_grp_id
        integer(hdf_err) :: err

        dim_1rdm_vec = size(one_rdm, dim=1) * (size(one_rdm, dim=1) + 1) / 2
        allocate(indices(2, dim_1rdm_vec), source=0_n_int)
        allocate(values(dim_1rdm_vec), source=0.0_dp)

        do p = 1, size(one_rdm, dim=1)
            do q = 1, size(one_rdm, dim=1)
                if (p >= q) then
                    pq = int(p * (p - 1)/2 + q)  ! Molcas style folding
                    indices(1, pq) = p; indices(2, pq) = q
                    if (abs(one_rdm(ind(p), ind(q))) > 1e-12_dp) then
                        values(pq) = one_rdm(ind(p), ind(q))  ! assumed to be normalised
                    else
                        values(pq) = 0.0_dp
                    end if
                end if
            end do
        end do

        call h5gcreate_f(parent, '1100', onerdm_grp_id, err)
        call write_int64_2d_dataset(onerdm_grp_id, 'indices', indices)
        call write_dp_1d_dataset(onerdm_grp_id, 'values', values)
        call h5gclose_f(onerdm_grp_id, err)
    end subroutine write_1rdm_hdf5


    subroutine write_2rdm_hdf5(parent, rdm, rdm_trace, iroot)
        use rdm_data, only: rdm_list_t
        use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm
        integer(hid_t), intent(in) :: parent
        type(rdm_list_t), intent(in) :: rdm
        integer(hid_t) :: plist_id
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        integer, intent(in) :: iroot
        integer(hid_t) :: grp_id, filespace, memspace, rank, dset_id
        integer(hdf_err) :: err
        real(dp) :: rdm_sign(rdm%sign_length)
        integer(int_rdm) :: pqrs
        integer :: iproc, ielem, &
                   p, q, r, s, pq, rs, tot_length_RDM, max_length_RDM
        integer, allocatable :: index(:), indices(:,:)
        real(dp), allocatable :: values(:)
        integer(hsize_t) :: dimsf(1), count(1), offset(1)
        integer :: mpirank, mpinum, mpierr

        ! create dynamic "indices" and "value" arrays to append to
        values = [real(dp) ::]
        ! index = [integer ::]
        do ielem = 1, rdm%nelements
            pqrs = rdm%elements(0, ielem)
            call extract_sign_rdm(rdm%elements(:, ielem), rdm_sign)
            rdm_sign = rdm_sign / rdm_trace
            call calc_separate_rdm_labels(pqrs, pq, rs, p, s, q, r)
            if (abs(rdm_sign(iroot)) > 1.e-12_dp) then
                if (p >= q .and. pq >= rs .and. p >= r .and. p >= s) then
                    ! index = [index, p]; index = [index, q]
                    ! index = [index, r]; index = [index, s]
                    values = [values, rdm_sign(iroot)]
                end if
            end if
        end do
        ! indices = reshape(index, [4, size(values)])
        write(stdout,*) 'length values: ', size(values)

        ! RDM elements are non-evenly distributed across ranks
        call MPIAllReduce(size(values), MPI_SUM, tot_length_RDM)
        dimsf = tot_length_RDM
        rank = 1

        ! create the 2RDM group
        call h5gcreate_f(parent, '2200', grp_id, err)
        ! rank and length of dataspace
        call h5screate_simple_f(rank, dimsf, filespace, err)
        ! create the corresponding dataset and get dset_id
        call h5dcreate_f(grp_id, 'values', H5T_NATIVE_DOUBLE, filespace, &
                         dset_id, err)

        ! writing from memory p-times a 1x1 block, see below
        count = 1
        call h5screate_simple_f(rank, count, memspace, err)

        ! set I/O pattern to "collective"
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

        call MPIAllReduce(size(values), MPI_MAX, max_length_RDM)
        call MPI_Comm_Rank(MPI_Comm_World, mpirank, mpierr)
        call MPI_Comm_Size(MPI_Comm_World, mpinum, mpierr)
        do p = 1, int(max_length_RDM)
            offset = mpirank + (p - 1) * mpinum
            if (size(values) < p) then
                call h5sselect_none_f(filespace, err)
            else
                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
                                           count, err)
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values(p), count, err, &
                            file_space_id=filespace, mem_space_id=memspace, &
                            xfer_prp=plist_id)
            ! call MPIBarrier(mpierr)
        end do

        call h5sclose_f(filespace, err)
        call h5sclose_f(memspace, err)
        call h5dclose_f(dset_id, err)
        call h5gclose_f(grp_id, err)
        call h5pclose_f(plist_id, err)

    end subroutine write_2rdm_hdf5


end module rdm_hdf5
