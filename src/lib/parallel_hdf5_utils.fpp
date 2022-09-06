#:include "../macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

#:set RANKS = range(1, 3)
#:set DATA_TYPES = ['integer', 'real']
#:def _rankdims(RANK)
$:'dimsf(1)' if RANK == 1 else _rankdims(RANK-1) + ', dimsf({num})'.format(num=RANK)
#:enddef _rankdims

module parallel_hdf5_utils

  use MPI_wrapper
  use Parallel_neci
#ifdef USE_HDF_
  use hdf5
#endif

  implicit none
  private
  public :: read_data_phdf5, write_data_phdf5

  interface write_data_phdf5
#:for DATA_TYPE in DATA_TYPES
  #:for RANK in RANKS
    module procedure write_${RANK}$d_${DATA_TYPE}$_data_phdf5
  #:endfor
#:endfor
  end interface

  interface read_data_phdf5
#:for DATA_TYPE in DATA_TYPES
  #:for RANK in RANKS
    module procedure read_${RANK}$d_${DATA_TYPE}$_data_phdf5
  #:endfor
#:endfor
  end interface

contains


#:for DATA_TYPE in DATA_TYPES
  #:for RANK in RANKS
    !> Writes in parallel contiguous up to 2D data distributed over MPI ranks into an archive.
    !> This routine assumes that the data varies across MPI ranks only in the final
    !> dimension
    subroutine write_${RANK}$d_${DATA_TYPE}$_data_phdf5(wrt_buf, dataset_name, grp_id)
      !> Data to be written
      #:if DATA_TYPE == 'real'
        ${DATA_TYPE}$(dp), intent(in) ${_ranksuffix(RANK)}$ :: wrt_buf
      #:elif DATA_TYPE == 'integer'
        ${DATA_TYPE}$, intent(in) ${_ranksuffix(RANK)}$ :: wrt_buf
      #:endif
      !> Name of dataset
      character(len=*), intent(in) :: dataset_name
      !> ID of group that dataset should be written into
      integer(hid_t), intent(in) :: grp_id
      !> Various filespace handles, rank of the tensor to be written
      integer(hid_t) :: filespace, memspace, dset_id, plist_id, rank
      !> dimension of dataset to be written, block size during writing and write offset
      integer(hsize_t) :: dimsf(${RANK}$), count(${RANK}$), offset(${RANK}$)
      !> HDF5 error code
      integer(hdf_err) :: err
      !> total length, i.e. sum of number of data on each MPI rank
      !> List containing length of distributed data on each MPI rank
      !> MPI error code
      integer :: tot_len_data, list_len_data(nProcessors), mpierr, data_shape(${RANK}$)

      data_shape = shape(wrt_buf)
      ! inquire about the size of the last dimension on each MPI rank and
      ! their sum which is required for mapping out the disk space.
      call MPIAllGather(data_shape(size(data_shape)), list_len_data, mpierr)
      call MPIAllReduce(data_shape(size(data_shape)), MPI_SUM, tot_len_data)
      #:if RANK == 1
        dimsf = [tot_len_data]
      #:elif RANK > 1
        dimsf = [data_shape(1:${RANK-1}$), tot_len_data]
      #:endif
      rank = ${RANK}$
      ! rank and length of dataspace
      call h5screate_simple_f(rank, dimsf, filespace, err)
      ! create the corresponding dataset and get dset_id
      #:if DATA_TYPE == 'real'
        call h5dcreate_f(grp_id, dataset_name, H5T_NATIVE_DOUBLE, &
                         filespace, dset_id, err)
      #:elif DATA_TYPE == 'integer'
        call h5dcreate_f(grp_id, dataset_name, H5T_NATIVE_INTEGER, &
                         filespace, dset_id, err)
      #:endif
      ! writing from memory a block called "count", mapping out the memspace
      count = shape(wrt_buf)
      call h5screate_simple_f(rank, count, memspace, err)

      ! Given two arrays of (4,2,2) and (4,2,20), the first one has offset
      ! (0,0,0) and the second (0,0,2)
      ! This also ensures that the ordering of different data sets can
      ! be correlated.
      if (iProcIndex == 0) then
          offset = 0 * data_shape
      else
      #:if RANK == 1
          offset = [sum(list_len_data(1:iProcIndex))]
      #:elif RANK > 1
          offset = [0 * data_shape(1:${RANK-1}$), sum(list_len_data(1:iProcIndex))]
      #:endif
      end if

      ! set I/O pattern to "collective"
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
      ! In HDF5 collective I/O mode all collective operations have to be
      ! called the same number of times by each MPI rank; however, this call
      ! may happen asynchronously. On some MPI ranks there might be no data,
      ! then file- and memspace have to be selected as none.
      if (data_shape(size(data_shape)) == 0) then
          call h5sselect_none_f(filespace, err)
          call h5sselect_none_f(memspace, err)
      else
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
                                     count, err)
      end if
      #:if DATA_TYPE == 'real'
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wrt_buf, count, err, &
      #:elif DATA_TYPE == 'integer'
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, wrt_buf, count, err, &
      #:endif
                      file_space_id=filespace, mem_space_id=memspace, &
                      xfer_prp=plist_id)

      call h5pclose_f(plist_id, err)
      call h5sclose_f(filespace, err)
      call h5sclose_f(memspace, err)
      call h5dclose_f(dset_id, err)
    end subroutine write_${RANK}$d_${DATA_TYPE}$_data_phdf5

    !> Read up to 2D data from an HDF5 archive in parallel on all processors.
    subroutine read_${RANK}$d_${DATA_TYPE}$_data_phdf5(rec_buf, dataset_name, grp_id)
      !> Receive buffer
      #:if DATA_TYPE == 'real'
        ${DATA_TYPE}$(dp), intent(inout) ${_ranksuffix(RANK)}$ :: rec_buf
      #:elif DATA_TYPE == 'integer'
        ${DATA_TYPE}$, intent(inout) ${_ranksuffix(RANK)}$ :: rec_buf
      #:endif
      !> Name of dataset
      character(len=*), intent(in) :: dataset_name
      !> ID of the group containing data set
      integer(hid_t), intent(in) :: grp_id
      !> Various filespace handles, rank of the tensor to be written
      integer(hid_t) :: dset_id, filespace, memspace, plist_id, rank
      !> dimension of dataset to be written, block size during writing and write offset
      integer(hsize_t) :: dimsf(${RANK}$), count(${RANK}$), offset(${RANK}$), &
                          maxdimsf(${RANK}$)
      !> HDF5 error code
      integer(hdf_err) :: err

      call h5dopen_f(grp_id, dataset_name, dset_id, err)
      call h5dget_space_f(dset_id, filespace, err)
      call h5sget_simple_extent_ndims_f(filespace, rank, err)
      call h5sget_simple_extent_dims_f(filespace, dimsf, maxdimsf, err)

      ! set I/O pattern to "collective"
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

      ! define memory dataspace and hyperslab
      offset = 0 * dimsf
      count = dimsf
      call h5screate_simple_f(rank, dimsf, memspace, err)
      call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, count, err)

      ! define filespace hyperslab
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, err)

      #:if DATA_TYPE == 'real'
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rec_buf, dimsf, err, &
      #:elif DATA_TYPE == 'integer'
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, rec_buf, dimsf, err, &
      #:endif
                       file_space_id=filespace, mem_space_id=memspace, &
                       xfer_prp=plist_id)

      call h5dclose_f(dset_id, err)
      call h5pclose_f(plist_id, err)
      call h5sclose_f(filespace, err)
      call h5sclose_f(memspace, err)
    end subroutine
  #:endfor
#:endfor

end module
