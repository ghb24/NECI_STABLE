#include "macros.h"

module hdf5_util

    ! This module contains helper functions for reading/writing elements to
    ! HDF5 files.
    !
    ! N.B. This does __not__ make them generic and interfaced.
    !
    ! It is important that we are explict about the data format in the hdf5
    ! file, as this is intended for clean interoperability. If we are to
    ! guarantee compatibility across programming languages, compilers and
    ! build configurations, then it helps to be explicit.

#ifdef __USE_HDF
    use iso_c_hack
    use constants
    use util_mod
    use hdf5
    use ParallelHelper
    implicit none

    interface read_int32_attribute
        module procedure read_int32_attribute_main
        module procedure read_int32_attribute_cast
    end interface

    interface write_int64_1d_dataset
        module procedure write_int64_1d_dataset_4
        module procedure write_int64_1d_dataset_8
    end interface

    interface read_int64_1d_dataset
        module procedure read_int64_1d_dataset_4
        module procedure read_int64_1d_dataset_8
    end interface

    interface write_log_scalar
        module procedure write_log_scalar_4
        module procedure write_log_scalar_8
    end interface

    interface read_log_scalar
        module procedure read_log_scalar_4
        module procedure read_log_scalar_8
    end interface

    interface write_int64_scalar
        module procedure write_int64_scalar_4
        module procedure write_int64_scalar_8
    end interface

    interface read_int64_scalar
        module procedure read_int64_scalar_4
        module procedure read_int64_scalar_8
    end interface

    integer :: tmp_lenof_sign

contains

    subroutine write_int32_attribute(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(in) :: val

        integer(hid_t) :: dataspace, attribute
        integer(hdf_err) :: err

        ! Create a scalar dataspace
        call h5screate_f(H5S_SCALAR_F, dataspace, err)

        ! Create the attribute with the correct type
        call h5acreate_f(parent, nm, H5T_NATIVE_INTEGER_4, dataspace, &
                         attribute, err)

        ! Write the data
        call h5awrite_f(attribute, H5T_NATIVE_INTEGER_4, val, [1_hsize_t], err)

        ! Close the residual handles
        call h5aclose_f(attribute, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine write_int64_attribute(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(in) :: val

        integer(hid_t) :: dataspace, attribute
        integer(hdf_err) :: err
        integer(int32), pointer :: ptr

        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5acreate_f(parent, nm, H5T_NATIVE_INTEGER_8, dataspace, &
                         attribute, err)
        call ptr_abuse_scalar(val, ptr)
        call h5awrite_f(attribute, H5T_NATIVE_INTEGER_8, ptr, [1_hsize_t], err)
        call h5aclose_f(attribute, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine read_int64_attribute(parent, nm, val, exists, default, &
                                    required)

        ! Read in a 64bit scalar attribute

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_int64_attribute'

        integer(hid_t) :: attribute
        integer(hdf_err) :: err
        logical(hdf_log) :: exists_
        integer(int32), pointer :: ptr

        call h5aexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5aopen_f(parent, nm, attribute, err)
            call ptr_abuse_scalar(val, ptr)
            call h5aread_f(attribute, H5T_NATIVE_INTEGER_8, ptr, &
                           [1_hsize_t], err)
            call h5aclose_f(attribute, err)
        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine read_int64_attribute

    subroutine write_string_attribute(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        character(*), intent(in) :: val

        integer(hid_t) :: dataspace, attribute, type_id
        integer(hdf_err) :: err

        ! Create an HDF type associated with a fortran string of _exactly_
        ! this length
        call h5tcopy_f(H5T_FORTRAN_S1, type_id, err)
        call h5tset_size_f(type_id, len(val, SIZE_T), err)

        ! Create a dataspace, and attribute, and write the string
        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5acreate_f(parent, nm, type_id, dataspace, attribute, err)
        call h5awrite_f(attribute, type_id, val, [1_hsize_t], err)
        call h5aclose_f(attribute, err)
        call h5sclose_f(dataspace, err)
        call h5tclose_f(type_id, err)

    end subroutine

    subroutine write_dp_1d_attribute(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(in) :: val(:)

        integer(hid_t) :: dataspace, attribute
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)

        dims = [size(val)]
        call h5screate_simple_f(1, dims, dataspace, err)
        call h5acreate_f(parent, nm, H5T_NATIVE_REAL_8, dataspace, attribute, &
                         err)
        call h5awrite_f(attribute, H5T_NATIVE_REAL_8, val, dims, err)
        call h5aclose_f(attribute, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine read_dp_1d_attribute(parent, nm, val, exists, default, required)

        ! Read in an array of real(dp)s

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(out) :: val(:)
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        real(dp), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_dp_1d_attribute'

        integer(hid_t) :: attribute
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)
        logical(hdf_log) :: exists_

        call h5aexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5aopen_f(parent, nm, attribute, err)
            dims = [size(val)]
            call check_attribute_params(attribute, nm, 8_hsize_t, H5T_FLOAT_F, dims)
            call h5aread_f(attribute, H5T_NATIVE_REAL_8, val, dims, err)
            call h5aclose_f(attribute, err)
        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6,*) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine

    subroutine write_dp_scalar(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(in) :: val

        integer(hid_t) :: dataspace, dataset
        integer(hdf_err) :: err

        ! Create a scalar dataspace, and dataset. Then write to it.
        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5dcreate_f(parent, nm, H5T_NATIVE_REAL_8, dataspace, &
                         dataset, err)
        if(iProcIndex.eq.0) call h5dwrite_f(dataset, H5T_NATIVE_REAL_8, val, [1_hsize_t], err)
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine read_dp_scalar(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        real(dp), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_dp_scalar'

        integer(hid_t) :: dataset
        integer(hdf_err) :: err
        logical(hdf_log) :: exists_

        ! Test if the relevant key exists
        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call h5dread_f(dataset, H5T_NATIVE_REAL_8, val, [1_hsize_t], err)
            if (err /= 0) &
                call stop_all(t_r, 'Read error')
            call h5dclose_f(dataset, err)
        endif

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine read_dp_scalar

    subroutine write_log_scalar_4(parent, nm, val)

        ! All logical values should be stored as 32bit integers.

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        logical(int32), intent(in) :: val

        integer(hid_t) :: dataspace, dataset
        integer(hdf_err) :: err
        integer(int32) :: tmp

        ! Store this as an integral value!
        if (val) then
            tmp = 1
        else
            tmp = 0
        end if

        ! Create a scalar dataspace, and dataset. Then write to it.
        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5dcreate_f(parent, nm, H5T_NATIVE_INTEGER_4, dataspace, &
                         dataset, err)
        if(iProcIndex.eq.0) call h5dwrite_f(dataset, H5T_NATIVE_INTEGER_4, tmp, [1_hsize_t], err)
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine write_log_scalar_8(parent, nm, val)

        ! Wrapper around the _4 routine

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        logical(int64), intent(in) :: val
        logical(int32) :: tmp

        tmp = val
        call write_log_scalar_4(parent, nm, tmp)

    end subroutine

    subroutine read_log_scalar_4(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        logical(int32), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        logical(int32), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_log_scalar_4'

        integer(hid_t) :: dataset
        integer(hdf_err) :: err
        logical(hdf_log) :: exists_
        integer(int32) :: tmp

        ! Test if the relevant key exists
        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call h5dread_f(dataset, H5T_NATIVE_INTEGER_4, tmp, [1_hsize_t], err)
            if (err /= 0) &
                call stop_all(t_r, 'Read error')
            if (tmp == 0) then
                val = .false.
            else
                val = .true.
            end if
            call h5dclose_f(dataset, err)
        endif

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine

    subroutine read_log_scalar_8(parent, nm, val, exists, default, required)

        ! Wrap the logical reading, such that it doesn't matter what type is
        ! being used internally by the code, the data format in the file
        ! should remain constant

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        logical(int64), intent(out) :: val
        logical, intent(out), optional :: exists
        logical(int32), intent(in), optional :: default
        logical, intent(in), optional :: required

        logical(int32) :: buf
        logical :: exists_

        call read_log_scalar_4(parent, nm, buf, exists_, default, required)

        if (exists_ .or. present(default)) &
            val = buf
        if (present(exists)) exists = exists_

    end subroutine

    subroutine write_int64_scalar_8(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(in) :: val

        integer(hid_t) :: dataspace, dataset
        integer(hdf_err) :: err
        integer(int32), pointer :: ptr

        ! Create a scalar dataspace, and dataset. Then write to it.
        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5dcreate_f(parent, nm, H5T_NATIVE_INTEGER_8, dataspace, &
                         dataset, err)
        call ptr_abuse_scalar(val, ptr)
        if(iProcIndex.eq.0) call h5dwrite_f(dataset, H5T_NATIVE_INTEGER_8, ptr, [1_hsize_t], err)
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine write_int64_scalar_4(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(in) :: val

        call write_int64_scalar_8(parent, nm, int(val, int64))

    end subroutine

    subroutine read_int64_scalar_8(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_int64_scalar'

        integer(hid_t) :: dataset
        integer(hdf_err) :: err
        logical(hdf_log) :: exists_
        integer(int32), pointer :: ptr

        ! Test if the relevant key exists
        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call ptr_abuse_scalar(val, ptr)
            call h5dread_f(dataset, H5T_NATIVE_INTEGER_8, ptr, [1_hsize_t], err)
            if (err /= 0) &
                call stop_all(t_r, 'Read error')
            call h5dclose_f(dataset, err)
        endif

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine read_int64_scalar_8

    subroutine read_int64_scalar_4(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default

        integer(int64) :: buf
        logical :: exists_

        call read_int64_scalar_8(parent, nm, buf, exists_, default, required)

        if (exists_ .or. present(default)) &
            val = int(buf, int32)
        if (present(exists)) exists = exists_

    end subroutine

    subroutine write_int64_1d_dataset_8(parent, nm, val)

        ! Write out a 64-bit array
        ! TODO: It may be worth refactoring multiple-types to remove the
        !       generic array writing code from the specific.

        integer(hid_t) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(in), target :: val(:)

        integer(hid_t) :: dataspace, dataset
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)
        integer(int32), pointer :: ptr(:)

        ! Create the appropriate dataspace
        dims = [size(val)]
        call h5screate_simple_f(1, dims, dataspace, err)

        ! Create the dataset with the correct type
        call h5dcreate_f(parent, nm, H5T_NATIVE_INTEGER_8, dataspace, &
                         dataset, err)

        ! write the data
        call ptr_abuse_1d(val, ptr)
        if(iProcIndex.eq.0) call h5dwrite_f(dataset, H5T_NATIVE_INTEGER_8, ptr, dims, err)

        ! Close the residual handles
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine write_int64_1d_dataset_4(parent, nm, val)

        ! Write out a 64-bit array

        integer(hid_t) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(in), target :: val(:)

        call write_int64_1d_dataset_8(parent, nm, int(val, int64))

    end subroutine

    subroutine write_dp_1d_dataset(parent, nm, val)

        integer(hid_t) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(in), target :: val(:)

        integer(hid_t) :: dataspace, dataset
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)

        ! Create the appropriate dataspace
        dims = [size(val)]
        call h5screate_simple_f(1, dims, dataspace, err)

        ! Create the dataset with the correct type
        call h5dcreate_f(parent, nm, H5T_NATIVE_REAL_8, dataspace, &
                         dataset, err)

        ! write the data
        if(iProcIndex.eq.0) call h5dwrite_f(dataset, H5T_NATIVE_REAL_8, val, dims, err)

        ! Close the residual handles
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine read_dp_1d_dataset(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(out), target :: val(:)
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        real(dp), intent(in), optional :: default(:)
        character(*), parameter :: t_r = 'read_dp_1d_dataset'

        integer(hid_t) :: dataset, type_id
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)
        logical(hdf_log) :: exists_
        real(dp), allocatable :: buf(:)

        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call h5dget_type_f(dataset, type_id, err)

            ! set up the read-buffer: We might need to add/remove replicas
            call setup_dp_1d_dataset_buffer(buf,val)

            ! Check dimensions and types.
            dims = [size(buf)]
            ! check versus the input, not the calculation's parameters
            call check_dataset_params(dataset, nm, 8_hsize_t, H5T_FLOAT_F, dims)

            ! And actually read the data.
            call h5dread_f(dataset, type_id, buf, dims, err)

            ! and move the data to val
            call move_dp_1d_dataset_buffer(val,buf)

            call h5tclose_f(type_id, err)
            call h5dclose_f(dataset, err)

        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine read_dp_1d_dataset

    subroutine setup_dp_1d_dataset_buffer(buf,val)
      ! allocate a buffer for reading dp_1d_datasets
      implicit none
      real(dp), allocatable, intent(out) :: buf(:)
      real(dp), target, intent(in) :: val(:)

      integer :: dims, ierr

      dims = size(val)
      ! if we change the number of replicas, we have to be careful
      if(lenof_sign /= tmp_lenof_sign) then
         if(dims .eq. lenof_sign) then
            allocate(buf(tmp_lenof_sign), stat = ierr)
         endif
      else
         ! else, the buffer is not very interesting
         allocate(buf(dims))
      endif
    end subroutine setup_dp_1d_dataset_buffer

    subroutine move_dp_1d_dataset_buffer(val,buf)
      ! moves the data from buf to val, eventually truncating/expanding it
      ! deallocates buf
      implicit none
      real(dp), allocatable, intent(inout) :: buf(:)
      real(dp), target, intent(inout) :: val(:)

      integer dimsVal, dimsBuf, ierr

      ! if buf is unallocated, this is not going anywhere
      if(.not. allocated(buf)) then
         write(iout,*) "WARNING: Trying to move data from empty buffer"
         return
      endif

      ! we need to check if the buffer can be copied 1:1
      dimsVal = size(val)
      dimsBuf = size(buf)
      if(dimsVal .eq. dimsBuf) then
         val(:) = buf(:)
      else
         ! now it depends if val is larger or smaller than buf
         if(dimsVal < dimsBuf) then
            ! depending, we either omit the last entries
            val(1:dimsVal) = buf(1:dimsBuf)
         else
            ! or copy the last one
            val(1:dimsBuf) = buf(1:dimsBuf)
            val(dimsBuf+1:dimsVal) = buf(dimsBuf)
         endif
      endif

      deallocate(buf)

    end subroutine move_dp_1d_dataset_buffer

    subroutine h5t_complex_t(dtype)

        ! Construct a datatype for complex numbers (dp)
        integer(hid_t), intent(out) :: dtype
        integer(hdf_err) :: err

        ! Complex numbers are two 8-byte floating points long
        call h5tcreate_f(H5T_COMPOUND_F, 16_hsize_t, dtype, err)
        call h5tinsert_f(dtype, "real", 0_hsize_t, H5T_NATIVE_REAL_8, err)
        call h5tinsert_f(dtype, "real", 0_hsize_t, H5T_NATIVE_REAL_8, err)
        call h5tinsert_f(dtype, "imag", 8_hsize_t, H5T_NATIVE_REAL_8, err)

    end subroutine

    subroutine write_cplx_1d_dataset(parent, nm, val)

        integer(hid_t) :: parent
        character(*), intent(in) :: nm
        complex(dp), intent(in), target :: val(:)

        integer(hid_t) :: dataspace, dataset, dtype
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)
        integer(int32), pointer :: ptr(:)

        ! Create the appropriate dataspace
        dims = [size(val)]
        call h5screate_simple_f(1, dims, dataspace, err)

        ! Create the dataset with the correct type
        call h5t_complex_t(dtype)
        call h5dcreate_f(parent, nm, dtype, dataspace, dataset, err)

        ! write the data
        call ptr_abuse_1d(val, ptr)
        if(iProcIndex.eq.0) call h5dwrite_f(dataset, dtype, ptr, dims, err)

        ! Close the residual handles
        call h5tclose_f(dtype, err)
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine read_cplx_1d_dataset(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        complex(dp), intent(out), target :: val(:)
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        complex(dp), intent(in), optional :: default(:)
        character(*), parameter :: t_r = 'read_dp_1d_dataset'

        integer(hid_t) :: dataset, type_id
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)
        integer(int32), pointer :: ptr(:)
        logical(hdf_log) :: exists_

        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call h5dget_type_f(dataset, type_id, err)

            ! Check dimensions and types.
            dims = [size(val)]
            call check_dataset_params(dataset, nm, 16_hsize_t, H5T_COMPOUND_F, &
                                      dims)

            ! And actually read the data.
            call ptr_abuse_1d(val, ptr)
            call h5dread_f(dataset, type_id, ptr, dims, err)

            call h5tclose_f(type_id, err)
            call h5dclose_f(dataset, err)
        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine read_cplx_1d_dataset


   subroutine write_2d_multi_arr_chunk_buff( &
                       parent, nm, itype, arr, arr_dims, mem_dims, mem_offset, &
                       dataspace_dims, dataspace_offset)

        ! Write a chunk of memory from each of the MPI processes into the
        ! specified place in the output file.
        !
        ! mem_dims   - the dimensions of the chunk of the memory array that we
        !              want to write
        ! mem_Offset - the offset of the aforementioned chunk from the start
        !              of the array
        ! dataspace_dims
        !            - the dimensions of the overall dataspace
        ! dataspace_offset
        !            - the offset at which we want to write this data.

        integer(hid_t), intent(in) :: parent, itype
        character(*), intent(in) :: nm
        integer(hsize_t) :: arr(1:,1:)
        integer(hsize_t), intent(in) :: arr_dims(2)
        integer(hsize_t), intent(in) :: mem_dims(2), mem_offset(2)
        integer(hsize_t), intent(in) :: dataspace_dims(2), dataspace_offset(2)

        integer(hsize_t) :: buff_dims(2)
        integer(hid_t) :: memspace, dataspace, dataset, plist_id
        integer(hdf_err) :: err
        integer(hsize_t), dimension(:,:), allocatable :: arr_buff
        integer(hsize_t) :: block_start, block_end, block_size, this_block_size
        type(c_ptr) :: cptr
        integer(int32), pointer :: ptr(:)
        integer(TagIntType) :: arr_buff_tag
        integer :: ierr

        ! Create an array in the target file with the size of the total amount
        ! to be written out, across the processors
        call h5screate_simple_f(2, dataspace_dims, dataspace, err)

        ! Create a property list to do multi-process writes
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)

        ! Create the dataset with the correct type
        call h5dcreate_f(parent, nm, itype, dataspace, dataset, err)

        !we use our own, contiguous buffers and independent MPI-IO to improve scaling
        !limit buffer size to 50 MB per task
        block_size=50000000/sizeof(arr(1,1))/mem_dims(1)
        block_size=min(block_size,mem_dims(2))

        buff_dims=[mem_dims(1), block_size]

        ! Create the source (memory) dataspace
        call h5screate_simple_f(2, buff_dims, memspace, err)

        allocate(arr_buff(mem_dims(1),block_size),stat = ierr)
        if(block_size.gt.0) &
             call LogMemAlloc('arr_buff',size(arr_buff),int(sizeof(arr_buff(1,1))),&
             'write_2d_multi',arr_buff_tag,ierr)
        block_start=1
        block_end=min(block_start+block_size-1,mem_dims(2))
        this_block_size=block_end-block_start+1
        do while (this_block_size.gt.0)
           arr_buff(:,1:this_block_size)=&
                arr(1+mem_offset(1):mem_offset(1)+mem_dims(1),block_start+mem_offset(2):block_end+mem_offset(2))

           !the last block might be smaller than block_size
           if (this_block_size.lt.block_size) call h5sselect_hyperslab_f(memspace, &
                H5S_SELECT_SET_F, [0_hsize_t,0_hsize_t], [buff_dims(1),this_block_size], err)

           call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                [dataspace_offset(1),dataspace_offset(2)+block_start-1], [buff_dims(1),this_block_size], err)

           ! Get access to a pointer to the array that will always be considered
           ! to have a valid type by the HDF5 library. For some reason the
           ! 64-bit, 2d array is not always given an interface...
           cptr=arr_2d_ptr(arr_buff)
           call c_f_pointer(cptr, ptr, [1])
           call h5dwrite_f(dataset, itype, ptr, buff_dims, err, memspace, &
                dataspace, xfer_prp=plist_id)
           block_start=block_start+block_size
           block_end=min(block_start+block_size-1,mem_dims(2))
           this_block_size=block_end-block_start+1
        end do

        deallocate(arr_buff)
        if(block_size.gt.0) call LogMemDealloc('write_2d_multi',arr_buff_tag)

        call h5dclose_f(dataset, err)
        call h5pclose_f(plist_id, err)
        call h5sclose_f(dataspace, err)
        call h5sclose_f(memspace, err)

    end subroutine write_2d_multi_arr_chunk_buff


    subroutine read_2d_multi_chunk(dataset, val, itype, dims, src_offset, &
                                   tgt_offset)

        ! Read in a specified chunk of the specified dataset.

        integer(hid_t), intent(in) :: dataset, itype
        integer(int64), intent(out) :: val(:,:)
        integer(hsize_t), intent(in) :: dims(2), src_offset(2), tgt_offset(2)

        integer(hid_t) :: plist_id, dataspace, memspace
        integer(hdf_err) :: err
        integer(hsize_t) :: mem_dims(2)
        integer(int32), pointer :: ptr(:,:)

        ! Create a property list to do multi-process reads
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        ! We use independent reads to contiguous buffers, which is the most
        ! reliable/scalable option. We use explicit buffering because there
        ! is a massive performance problem when writing directly to a
        ! (non-contiguous) hyperslab in SpawnedParts. Buffering by collective
        ! MPI-IO can create performance and memory problems for large task counts
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)

        ! Create the target (memory) dataspace, and select the appropriate
        ! hyperslab inside it.
        mem_dims = [size(val, 1), size(val, 2)]
        call h5screate_simple_f(2, mem_dims, memspace, err)
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, tgt_offset, &
                                   dims, err)

        ! Select the hyperslab in the array to be read
        call h5dget_space_f(dataset, dataspace, err)
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, src_offset, &
                                   dims, err)

        ! And do the actual read!
        call ptr_abuse_2d(val, ptr)
        call h5dread_f(dataset, itype, ptr, mem_dims, err, memspace, &
                       dataspace, xfer_prp=plist_id)

        call h5sclose_f(memspace, err)
        call h5sclose_f(dataspace, err)
        call h5pclose_f(plist_id, err)

    end subroutine

    subroutine read_int32_attribute_main(parent, nm, val, exists, default, &
                                         required)

        ! Read in a 32bit scalar attribute

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int32), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_int32_attribute_main'

        integer(hid_t) :: attribute
        integer(hdf_err) :: err
        logical(hdf_log) :: exists_

        call h5aexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5aopen_f(parent, nm, attribute, err)
            call h5aread_f(attribute, H5T_NATIVE_INTEGER_4, val, &
                           [1_hsize_t], err)
            call h5aclose_f(attribute, err)
        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine

    subroutine read_int32_attribute_cast(parent, nm, val, exists, default, &
                                         required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(out) :: val
        logical, intent(out), optional :: exists
        integer(int32), intent(in), optional :: default
        logical, intent(in), optional :: required
        integer(int32) :: val_tmp

        logical :: exists_

        call read_int32_attribute_main(parent, nm, val_tmp, exists_, default, &
                                       required)

        if (exists_ .or. present(default)) &
            val = int(val_tmp, int64)
        if (present(exists)) exists = exists_

    end subroutine

    subroutine read_string_attribute(parent, nm, val, exists, default, &
                                     required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        character(*), intent(out) :: val
        logical, intent(out), optional :: exists
        character(*), intent(in), optional :: default
        logical, intent(in), optional :: required
        character(*), parameter :: t_r = 'read_string_attribute'

        integer(hid_t) :: attribute, type_id
        integer(hdf_err) :: err
        integer(hsize_t) :: sz, buf_sz
        logical(hdf_log) :: exists_

        call h5aexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5aopen_f(parent, nm, attribute, err)
            call h5aget_type_f(attribute, type_id, err)
            ! TODO: test that this is a string type/class
            !call h5tget_class_f(type_id, class_id, err)
            call h5tget_size_f(type_id, sz, err)

            ! We can only read in if our buffer is big enough
            buf_sz = len(val)
            if (sz > buf_sz) then
                write(6,*) 'WARNING: Insufficient read buffer in routine ', t_r
                exists_ = .false.
            else
                ! Read in, and ensure that length/string termination is correct
                call h5aread_f(attribute, type_id, val, [1_hsize_t], err)
                val = val(1:sz)
            end if

            call h5tclose_f(type_id, err)
            call h5aclose_f(attribute, err)
        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = trim(default)

    end subroutine

    subroutine read_int64_1d_dataset_4( &
                                parent, nm, val, exists, default, required)

        ! Allow these values to be casted onto 32bit arrays if we are in a
        ! 32-bit build

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(out), target :: val(:)
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default(:)

        integer(int64) :: buf(size(val))

        call read_int64_1d_dataset_8(parent, nm, buf, exists, default, &
                                     required)
        val = buf

    end subroutine

    subroutine read_int64_1d_dataset_8( &
                                parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(out), target :: val(:)
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default(:)
        character(*), parameter :: t_r = 'read_int64_1d_dataset_8'

        integer(hid_t) :: dataset, type_id
        integer(hdf_err) :: err
        integer(hsize_t) :: dims(1)
        logical(hdf_log) :: exists_
        integer(int32), pointer :: ptr(:)

        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call h5dget_type_f(dataset, type_id, err)

            ! Check dimensions and types.
            dims = [size(val)]
            call check_dataset_params(dataset, nm, 8_hsize_t, H5T_INTEGER_F, dims)

            ! And actually read the data. Note that we manipulate the pointer
            ! to be a 32bit integer, as that is always present in the library.
            call ptr_abuse_1d(val, ptr)
            call h5dread_f(dataset, type_id, ptr, dims, err)

            call h5tclose_f(type_id, err)
            call h5dclose_f(dataset, err)
        end if

        if (present(required)) then
            if (required .and. .not. exists_) then
                write(6, *) nm
                call stop_all(t_r, "Required field does not exist")
            end if
        end if
        if (present(exists)) exists = exists_
        if (present(default) .and. .not. exists_) val = default

    end subroutine

    subroutine check_dataset_params(dataset, nm, sz, class_id, dims)

        ! Check that a dataset has the character that is required of it
        ! before reading in.

        integer(hid_t), intent(in) :: dataset
        integer(hdf_err), intent(in) :: class_id
        character(*), intent(in) :: nm
        integer(hsize_t), intent(in) :: sz
        integer(hsize_t), intent(in) :: dims(:)
        character(*), parameter :: t_r = 'check_dataset_params'

        integer(hid_t) :: rank, type_id

        integer(hid_t) :: dataspace
        integer(hdf_err) :: err, ds_class, ds_rank
        integer(hsize_t) :: ds_dims(size(dims)), ds_max_dims(size(dims)), ds_sz

        ! Get the type associated with the dataset. Check that it is an
        ! array with components that have the right number of bytes, and the
        ! correct base class type
        call h5dget_type_f(dataset, type_id, err)
        call h5tget_size_f(type_id, ds_sz, err)
        call h5tget_class_f(type_id, ds_class, err)
        call h5tclose_f(type_id, err)

        if (ds_sz /= sz .or. ds_class /= class_id) then
            write(6,*) 'Dataset name: ', nm
            call stop_all(t_r, "Invalid dataset type information found")
        end if

        ! Get the dataspace for the dataset. Check that the dataset has the
        ! requested dimensions
        call h5dget_space_f(dataset, dataspace, err)
        call h5sget_simple_extent_ndims_f(dataspace, ds_rank, err)

        rank = size(dims)
        if (rank /= ds_rank) then
            write(6,*) 'Dataset name: ', nm
            call stop_all(t_r, "Invalid dataset rank found")
        end if

        call h5sget_simple_extent_dims_f(dataspace, ds_dims, ds_max_dims, err)
        if (.not. all(dims == ds_dims)) then
            write(6,*) 'Dataset name: ', nm
            call stop_all(t_r, "Invalid dataset dimensions found")
        end if
        call h5sclose_f(dataspace, err)

    end subroutine check_dataset_params

    subroutine check_attribute_params(attribute, nm, sz, class_id, dims)

        ! Check that a attribute has the character that is required of it
        ! before reading in.

      integer(hid_t), intent(in) :: attribute
      integer(hdf_err), intent(in) :: class_id
      character(*), intent(in) :: nm
      integer(hsize_t), intent(in) :: sz
      integer(hsize_t), intent(in) :: dims(:)
      character(*), parameter :: t_r = 'check_attribute_params'

      integer(hid_t) :: rank, type_id

      integer(hid_t) :: dataspace
      integer(hdf_err) :: err, ds_class, ds_rank
      integer(hsize_t) :: ds_dims(size(dims)), ds_max_dims(size(dims)), ds_sz

      ! Get the type associated with the attribute. Check that it is an
      ! array with components that have the right number of bytes, and the
      ! correct base class type
      call h5aget_type_f(attribute, type_id, err)
      call h5tget_size_f(type_id, ds_sz, err)
      call h5tget_class_f(type_id, ds_class, err)
      call h5tclose_f(type_id, err)

      if (ds_sz /= sz .or. ds_class /= class_id) then
         write(6,*) 'Attribute name: ', nm
         call stop_all(t_r, "Invalid attribute type information found")
      end if

      ! Get the dataspace for the attribute. Check that the attribute has the
      ! requested dimensions
      call h5aget_space_f(attribute, dataspace, err)
      call h5sget_simple_extent_ndims_f(dataspace, ds_rank, err)

      rank = size(dims)
      if (rank /= ds_rank) then
         write(6,*) 'Attribute name: ', nm
         write(6,*) 'ranks', rank, ds_rank
         call stop_all(t_r, "Invalid attribute rank found")
      end if

      call h5sget_simple_extent_dims_f(dataspace, ds_dims, ds_max_dims, err)
      if (.not. all(dims == ds_dims)) then
         write(6,*) 'Attribute name: ', nm
         call stop_all(t_r, "Invalid attribute dimensions found")
      end if
      call h5sclose_f(dataspace, err)

    end subroutine check_attribute_params
#endif

end module
