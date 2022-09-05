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

    !> Write all RDMs specified in the input to an HDF5 archive.
    subroutine write_rdms_hdf5(rdm_defs, rdm, rdm_trace, one_rdms)
      use rdm_data, only: rdm_definitions_t, rdm_list_t, one_rdm_t
      use LoggingData, only: tWriteSpinFreeRDM, tPrint1RDM
      !> Type carrying information such as number and type of RDM (regular/transition)
      type(rdm_definitions_t), intent(in) :: rdm_defs
      !> 2RDM data distributed over all MPI ranks
      type(rdm_list_t), intent(in) :: rdm
      !> Normalisation of the 2RDM
      real(dp), intent(in) :: rdm_trace(rdm%sign_length)
      !> Array carrying the 1RDM belonging to each root
      type(one_rdm_t), intent(inout), optional :: one_rdms(:)
      !> Required for the stop all to check for HDF5 library
      character(*), parameter :: this_routine = 'write_hdf5_rdms'
#ifdef USE_HDF_
      !> File and group handles
      integer(hid_t) :: file_id, root_id, rdm_id, plist_id
      !> HDF5 error code
      integer(hdf_err) :: err
      !> Name of the HDF5 archive
      character(255) :: filename
      !> The iroot'th eigenvector of Hamiltonian, e.g. 1 for ground state,
      !> 2 for first excited state, ...
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
        ! create the file collectively, H5F_ACC_TRUNC_F = overwrite existing file
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
        ! TODO: to be implemented later perhaps...
        ! if (tWrite3RDM) then
        !   call write_3rdm_hdf5(rdm_id, ...)
        !   write(stdout, *) "3RDM written to file."
        ! end if
        ! if (tWriteF4RDM) then
        !   call write_contracted_fock_hdf5(rdm_id, ...)
        !   write(stdout, *) "F.4RDM written to file."
        ! end if

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

    ! The steps for HDF5 I/O are always:
    ! 1. Obtain the dataset identifier.
    ! 2. Specify the memory datatype.
    ! 3. Specify the memory dataspace.
    ! 4. Specify the file dataspace.
    ! 5. Specify the transfer properties.
    ! 6. Perform the desired operation on the dataset.
    ! 7. Close the dataset.
    ! 8. Close the dataspace, datatype, and property list if
    ! necessary

    !> Write the 1RDM to an HDF5 archive.
    subroutine write_1rdm_hdf5(parent, one_rdm)
      use RotateOrbsData, only: ind => SymLabelListInv_rot
      !> HDF5 file handle of the parent directory.
      integer(hid_t), intent(in) :: parent
      !> 1RDM data redundant on each MPI rank
      real(dp), intent(inout) :: one_rdm(:, :)
      !> Used to accumulate the RDM indices before writing
      integer, allocatable :: index(:), indices(:,:)
      !> Accumulate the non-zero RDM values respectively
      real(dp), allocatable :: values(:)
      !> Spatial orbital indices
      integer :: p, q
      !> HDF5 file handle of 1RDM "1100" group
      integer(hid_t) :: grp_id
      !> HDF5 error code
      integer(hdf_err) :: err

      values = [real(dp) ::]
      index = [integer ::]
      ! All the 1RDM values are stored on all ranks redundantly anyway.
      ! Here the values and indices are only accumulated on the head node,
      ! to prevent writing the same data multiple times to the HDF5 archive.
      if (iProcIndex == 0) then
        do p = 1, size(one_rdm, dim=1)
          do q = 1, size(one_rdm, dim=1)
            if (p >= q) then
              if (abs(one_rdm(ind(p), ind(q))) > 1.e-12_dp) then
                index = [index, p]; index = [index, q]
                values = [values, one_rdm(ind(p), ind(q))]
              end if
            end if
          end do
        end do
      end if
      indices = reshape(index, [2, size(values)])

      call h5gcreate_f(parent, '1100', grp_id, err)
      call write_rdmvals_phdf5(values, 'values', grp_id)
      call write_rdmindices_phdf5(indices, 'indices', grp_id, rank_rdm=1)
      call h5gclose_f(grp_id, err)
    end subroutine write_1rdm_hdf5


    !> Write the 2RDM to an HDF5 archive.
    subroutine write_2rdm_hdf5(parent, rdm, rdm_trace, iroot)
      use rdm_data, only: rdm_list_t
      use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm
      !> HDF5 file handle of the parent directory.
      integer(hid_t), intent(in) :: parent
      !> 2RDM data distributed over all MPI ranks
      type(rdm_list_t), intent(in) :: rdm
      !> Normalisation of the 2RDM
      real(dp), intent(in) :: rdm_trace(rdm%sign_length)
      !> The 2RDM storage with symmetry introduces a sign that has to kept
      !> track of, refer to "rdm_data.F90" for further information
      real(dp) :: rdm_sign(rdm%sign_length)
      !> The iroot'th eigenvector of Hamiltonian, e.g. 1 for ground state,
      !> 2 for first excited state, ...
      integer, intent(in) :: iroot
      !> HDF5 file handle of 1RDM "1100" group
      integer(hid_t) :: grp_id
      !> HDF5 error code
      integer(hdf_err) :: err
      !> Folded 2RDM loop index, refer to stochastic GUGA-CASSCF paper 2021 SI
      integer(int_rdm) :: pqrs
      !> Spatial orbital and loop indices
      integer :: ielem, p, q, r, s, pq, rs
      !> Used to accumulate the RDM indices before writing
      integer, allocatable :: index(:), indices(:,:)
      !> Accumulate the non-zero RDM values respectively
      real(dp), allocatable :: values(:)

      ! create dynamic "indices" and "value" arrays and append to them in loop
      values = [real(dp) ::]
      index = [integer ::]
      do ielem = 1, rdm%nelements
          pqrs = rdm%elements(0, ielem)
          call extract_sign_rdm(rdm%elements(:, ielem), rdm_sign)
          rdm_sign = rdm_sign / rdm_trace
          call calc_separate_rdm_labels(pqrs, pq, rs, p, s, q, r)
          if (abs(rdm_sign(iroot)) > 1.e-12_dp) then
              if (p >= q .and. pq >= rs .and. p >= r .and. p >= s) then
                  index = [index, p]; index = [index, q]
                  index = [index, r]; index = [index, s]
                  values = [values, rdm_sign(iroot)]
              end if
          end if
      end do
      indices = reshape(index, [4, size(values)])

      call h5gcreate_f(parent, '2200', grp_id, err)  ! create the 2RDM group
      call write_rdmvals_phdf5(values, 'values', grp_id)
      call write_rdmindices_phdf5(indices, 'indices', grp_id, rank_rdm=2)
      call h5gclose_f(grp_id, err)
    end subroutine write_2rdm_hdf5


    !> Writes in parallel RDM values distributed over MPI ranks into an archive.
    subroutine write_rdmvals_phdf5(data, dataset_name, grp_id)
        !> Data to be written
        real(dp), intent(in) :: data(:)
        !> Name of dataset, either "values" or "indices"
        character(len=*), intent(in) :: dataset_name
        !> ID of group that dataset should be written into
        integer(hid_t), intent(in) :: grp_id
        !> Various filespace handles, rank of the tensor to be written
        integer(hid_t) :: filespace, memspace, dset_id, plist_id, rank
        !> dimension of dataset to be written, block size during writing and write offset
        integer(hsize_t) :: dimsf(1), count(1), offset(1)
        !> HDF5 error code
        integer(hdf_err) :: err
        !> Maximum length of data vector over all ranks,
        !> total length, i.e. sum of number of data on each rank
        integer :: max_len_data, tot_len_data
        !> List containing length of distributed data on each rank
        integer :: list_len_data(nProcessors)
        !> Generic loop index and MPI error code
        integer :: i, mpierr
        real(dp), allocatable :: rec_buf(:)
        logical :: file_exists

        call MPIAllReduce(size(data), MPI_MAX, max_len_data)
        call MPIAllGather(size(data), list_len_data, mpierr)
        call MPIAllReduce(size(data), MPI_SUM, tot_len_data)
        dimsf = tot_len_data
        rank = 1
        ! rank and length of dataspace
        call h5screate_simple_f(rank, dimsf, filespace, err)
        ! create the corresponding dataset and get dset_id
        call h5dcreate_f(grp_id, dataset_name, H5T_NATIVE_DOUBLE, filespace, &
                         dset_id, err)
        ! writing from memory p-times a 1x1 block, called "count", see below
        count = 1
        call h5screate_simple_f(rank, count, memspace, err)
        ! set I/O pattern to "collective"
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

        if (iProcIndex == 0) then
            offset = 0
        else
            offset = 0 + sum(list_len_data(1:iProcIndex))
        end if
        ! In HDF5 collective I/O mode all collective operations have to be
        ! called the same number of times by each MPI rank; however, this call
        ! may happen asynchronously. The write function has to be called twice
        ! here since on some rank the "data" array may be of length 0. In
        ! these cases a scalar is provided as input argument, but not written
        ! to disk!
        write(100+iProcIndex,*) "write from NECI"
        do i = 1, int(max_len_data)
            if (size(data) < i) then
                ! both memspace and filespace need to be selected as none,
                ! otherwise the write does not complete.
                call h5sselect_none_f(filespace, err)
                call h5sselect_none_f(memspace, err)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, 0.0_dp, count, err, &
                                file_space_id=filespace, mem_space_id=memspace, &
                                xfer_prp=plist_id)
            else
                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
                                           count, err)
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data(i), count, err, &
                                file_space_id=filespace, mem_space_id=memspace, &
                                xfer_prp=plist_id)
            write(100+iProcIndex,*) data(i)
            end if
            offset = offset + 1
            ! prevent going out of bounds
            if (offset(1) >= tot_len_data) offset = tot_len_data - 1
        end do

        ! only as test for now
        inquire(file=trim('fciqmc.rdms.1.h5'), exist=file_exists)
        if (file_exists) then
            call read_1ddata_phdf5(dset_id, rec_buf)
            write(100+iProcIndex,*) "read from HDF5"
            do i = 1, int(max_len_data)
                write(100+iProcIndex,*) rec_buf(i)
            end do
        end if

        call h5pclose_f(plist_id, err)
        call h5sclose_f(filespace, err)
        call h5sclose_f(memspace, err)
        call h5dclose_f(dset_id, err)
    end subroutine write_rdmvals_phdf5


    !> Writes in parallel RDM indices distributed over MPI ranks into an archive.
    subroutine write_rdmindices_phdf5(data, dataset_name, grp_id, rank_rdm)
        !> Data to be written
        integer, intent(in) :: data(:,:)
        !> Name of dataset
        character(len=*), intent(in) :: dataset_name
        !> ID of group that dataset should be written into
        integer(hid_t), intent(in) :: grp_id
        !> Rank of the N-body RDM, e.g. 1, 2, 3, 4, ...
        integer, intent(in) :: rank_rdm
        !> Various filespace handles, rank of the tensor to be written
        integer(hid_t) :: filespace, memspace, dset_id, plist_id, rank
        !> dimension of dataset to be written, block size during writing and write offset
        integer(hsize_t) :: dimsf(2), count(2), offset(2)
        !> HDF5 error code
        integer(hdf_err) :: err
        !> Maximum length of data vector over all ranks,
        !> total length, i.e. sum of number of data on each rank
        integer :: max_len_data, tot_len_data
        !> List containing length of distributed data on each rank
        integer :: list_len_data(nProcessors)
        !> Generic loop index and MPI error code
        integer :: i, mpierr

        call MPIAllReduce(size(data, dim=2), MPI_MAX, max_len_data)
        call MPIAllGather(size(data, dim=2), list_len_data, mpierr)
        call MPIAllReduce(size(data, dim=2), MPI_SUM, tot_len_data)
        dimsf = [2 * rank_rdm, tot_len_data]
        rank = 2
        ! rank and length of dataspace
        call h5screate_simple_f(rank, dimsf, filespace, err)
        ! create the corresponding dataset and get dset_id
        call h5dcreate_f(grp_id, dataset_name, H5T_NATIVE_INTEGER, filespace, &
                         dset_id, err)
        ! writing from memory p-times a 1 x (2 x rank) block, called "count", see below
        count = [2 * rank_rdm, 1]
        call h5screate_simple_f(rank, count, memspace, err)
        ! set I/O pattern to "collective"
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

        if (iProcIndex == 0) then
            offset = [0, 0]
        else
            offset = [0, sum(list_len_data(1:iProcIndex))]
        end if
        ! In HDF5 collective I/O mode all collective operations have to be
        ! called the same number of times by each MPI rank; however, this call
        ! may happen asynchronously. The write function has to be called twice
        ! here since on some rank the "data" array may be of length 0. In
        ! these cases a scalar is provided as input argument, but not written
        ! to disk!
        do i = 1, int(max_len_data)
            if (size(data, dim=2) < i) then
                ! both memspace and filespace need to be selected as none,
                ! otherwise the write does not complete.
                call h5sselect_none_f(filespace, err)
                call h5sselect_none_f(memspace, err)
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data(:,i), &
                                count, err, file_space_id=filespace, mem_space_id=memspace, &
                                xfer_prp=plist_id)
            else
                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
                                           count, err)
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data(:,i), count, err, &
                                file_space_id=filespace, mem_space_id=memspace, &
                                xfer_prp=plist_id)
            end if
            offset = [0_hsize_t, offset(2) + 1_hsize_t]
            ! prevent going out of bounds
            if (offset(2) >= tot_len_data) offset = [0_hsize_t, tot_len_data - 1_hsize_t]
        end do

        call h5pclose_f(plist_id, err)
        call h5sclose_f(filespace, err)
        call h5sclose_f(memspace, err)
        call h5dclose_f(dset_id, err)
    end subroutine write_rdmindices_phdf5


    !> Read RDM data from an HDF5 archive in parallel on all processors.
    !> Data has to be distributed afterwards.
    subroutine read_1ddata_phdf5(dset_id, rec_buf)
        !> ID of the data set
        integer(hid_t), intent(in) :: dset_id
        !> Receive buffer
        real(dp), allocatable, intent(out) :: rec_buf(:)
        !> Various filespace handles, rank of the tensor to be written
        integer(hid_t) :: filespace, memspace, plist_id, rank, type_id, native_type_id
        !> dimension of dataset to be written, block size during writing and write offset
        integer(hsize_t) :: dimsf(1), count(1), offset(1), maxdimsf(1), dimsm(1)
        !> HDF5 error code
        integer(hdf_err) :: err

        call h5dget_type_f(dset_id, type_id, err)
        ! HDF5 has an internal list which it searches either top to bottom
        ! or the other way around and uses the first match to assign a dataset type.
        ! The default search direction is ascending and for now I will leave it as is.
        call h5tget_native_type_f(type_id, H5T_DIR_ASCEND_F, native_type_id, err)
        call h5dget_space_f(dset_id, filespace, err)
        call h5sget_simple_extent_ndims_f(filespace, rank, err)
        call h5sget_simple_extent_dims_f(filespace, dimsf, maxdimsf, err)

        allocate(rec_buf(dimsf(1)))

        ! set I/O pattern to "collective"
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

        ! define memory dataspace and hyperslab
        offset = 0
        count = dimsf
        call h5screate_simple_f(rank, dimsm, memspace, err)
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, count, err)

        ! define filespace hyperslab
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, err)

        call h5dread_f(dset_id, native_type_id, rec_buf, dimsf, err, &
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=plist_id)

        call h5pclose_f(plist_id, err)
        call h5sclose_f(type_id, err)
        call h5sclose_f(native_type_id, err)
        call h5sclose_f(filespace, err)
        call h5sclose_f(memspace, err)
    end subroutine read_1ddata_phdf5


end module rdm_hdf5
