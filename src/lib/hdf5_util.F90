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

    use iso_c_hack
    use constants
    use util_mod
    use hdf5
    implicit none

    interface read_int32_attribute
        module procedure read_int32_attribute_main
        module procedure read_int32_attribute_cast
    end interface

    interface read_int64_1d_dataset
        module procedure read_int64_1d_dataset_cast
        module procedure read_int64_1d_dataset_main
    end interface

    interface write_log_scalar
        module procedure write_log_scalar_4
        module procedure write_log_scalar_8
    end interface

    interface read_log_scalar
        module procedure read_log_scalar_4
        module procedure read_log_scalar_8
    end interface

contains

    subroutine write_int32_attribute(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int32), intent(in) :: val

        integer(hid_t) :: dataspace, attribute, err

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

        integer(hid_t) :: dataspace, attribute, err
        integer(int32), pointer :: ptr

        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5acreate_f(parent, nm, H5T_NATIVE_INTEGER_8, dataspace, &
                         attribute, err)
        call int32_pointer_abuse_scalar(val, ptr)
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

        integer(hid_t) :: attribute, err
        logical(hid_t) :: exists_
        integer(int32), pointer :: ptr

        call h5aexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5aopen_f(parent, nm, attribute, err)
            call int32_pointer_abuse_scalar(val, ptr)
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

        integer(hid_t) :: dataspace, attribute, type_id, err

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

    subroutine write_dp_scalar(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        real(dp), intent(in) :: val

        integer(hid_t) :: dataspace, dataset, err

        ! Create a scalar dataspace, and dataset. Then write to it.
        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5dcreate_f(parent, nm, H5T_NATIVE_REAL_8, dataspace, &
                         dataset, err)
        call h5dwrite_f(dataset, H5T_NATIVE_REAL_8, val, [1_hsize_t], err)
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

        integer(hid_t) :: dataset, err
        logical(hid_t) :: exists_

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

        integer(hid_t) :: dataspace, dataset, err
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
        call h5dwrite_f(dataset, H5T_NATIVE_INTEGER_4, tmp, [1_hsize_t], err)
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

        integer(hid_t) :: dataset, err
        logical(hid_t) :: exists_
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

    subroutine write_int64_scalar(parent, nm, val)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(in) :: val

        integer(hid_t) :: dataspace, dataset, err
        integer(int32), pointer :: ptr

        ! Create a scalar dataspace, and dataset. Then write to it.
        call h5screate_f(H5S_SCALAR_F, dataspace, err)
        call h5dcreate_f(parent, nm, H5T_NATIVE_INTEGER_8, dataspace, &
                         dataset, err)
        call int32_pointer_abuse_scalar(val, ptr)
        call h5dwrite_f(dataset, H5T_NATIVE_INTEGER_8, ptr, [1_hsize_t], err)
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine read_int64_scalar(parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(out) :: val
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default
        character(*), parameter :: t_r = 'read_int64_scalar'

        integer(hid_t) :: dataset, err
        logical(hid_t) :: exists_
        integer(int32), pointer :: ptr

        ! Test if the relevant key exists
        call h5lexists_f(parent, nm, exists_, err)
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call int32_pointer_abuse_scalar(val, ptr)
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

    end subroutine read_int64_scalar

    subroutine write_int64_1d_dataset(parent, nm, val)

        ! Write out a 64-bit array
        ! TODO: It may be worth refactoring multiple-types to remove the
        !       generic array writing code from the specific.

        integer(hid_t) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(in), target :: val(:)

        integer(hid_t) :: dataspace, dataset, err
        integer(hsize_t) :: dims(1)
        integer(int32), pointer :: ptr(:)

        ! Create the appropriate dataspace
        dims = [ubound(val) - lbound(val) + 1]
        call h5screate_simple_f(1_hid_t, dims, dataspace, err)

        ! Create the dataset with the correct type
        call h5dcreate_f(parent, nm, H5T_NATIVE_INTEGER_8, dataspace, &
                         dataset, err)

        ! write the data
        call int32_pointer_abuse(val, ptr)
        call h5dwrite_f(dataset, H5T_NATIVE_INTEGER_8, ptr, dims, err)

        ! Close the residual handles
        call h5dclose_f(dataset, err)
        call h5sclose_f(dataspace, err)

    end subroutine

    subroutine write_2d_multi_arr_chunk_offset( &
                       parent, nm, itype, cptr, arr_dims, mem_dims, mem_offset, &
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
        type(c_ptr), intent(in) :: cptr
        integer(hsize_t), intent(in) :: arr_dims(2)
        integer(hsize_t), intent(in) :: mem_dims(2), mem_offset(2)
        integer(hsize_t), intent(in) :: dataspace_dims(2), dataspace_offset(2)

        integer(hsize_t) :: dims(2)
        integer(hid_t) :: memspace, dataspace, dataset, plist_id, err
        integer(int32), pointer :: ptr(:)

        ! Create the source (memory) dataspace, and select the appropriate
        ! hyperslab inside it
        call h5screate_simple_f(2_hid_t, arr_dims, memspace, err)
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, mem_offset, &
                                   mem_dims, err)

        ! Create an array in the target file with the size of the total amount
        ! to be written out, across the processors
        call h5screate_simple_f(2_hid_t, dataspace_dims, dataspace, err)
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                   dataspace_offset, mem_dims, err)

        ! Create a property list to do multi-process collective writes
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

        ! Create the dataset with the correct type
        call h5dcreate_f(parent, nm, itype, dataspace, dataset, err)

        ! Get access to a pointer to the array that will always be considered
        ! to have a valid type by the HDF5 library. For some reason the
        ! 64-bit, 2d array is not always given an interface...
        call c_f_pointer(cptr, ptr, [1])

        call h5dwrite_f(dataset, itype, ptr, arr_dims, err, memspace, &
                        dataspace, xfer_prp=plist_id)

        call h5dclose_f(dataset, err)
        call h5pclose_f(plist_id, err)
        call h5sclose_f(dataspace, err)
        call h5sclose_f(memspace, err)

    end subroutine write_2d_multi_arr_chunk_offset

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

        integer(hid_t) :: attribute, err
        logical(hid_t) :: exists_

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

        integer(hid_t) :: attribute, err, type_id
        integer(hsize_t) :: sz, buf_sz
        logical(hid_t) :: exists_

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

    subroutine read_int64_1d_dataset_cast( &
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

        call read_int64_1d_dataset_main(parent, nm, buf, exists, default, &
                                        required)
        val = buf

    end subroutine

    subroutine read_int64_1d_dataset_main( &
                                parent, nm, val, exists, default, required)

        integer(hid_t), intent(in) :: parent
        character(*), intent(in) :: nm
        integer(int64), intent(out), target :: val(:)
        logical, intent(out), optional :: exists
        logical, intent(in), optional :: required
        integer(int64), intent(in), optional :: default(:)
        character(*), parameter :: t_r = 'read_int64_1d_dataset_main'

        integer(hid_t) :: dataset, err, type_id
        integer(hsize_t) :: dims(1)
        logical(hid_t) :: exists_
        integer(int32), pointer :: ptr(:)

        ! TODO: Fix this!!!
!        call h5dexists_f(parent, nm, exists_, err)
        exists_ = .true.
        if (exists_) then
            call h5dopen_f(parent, nm, dataset, err)
            call h5dget_type_f(dataset, type_id, err)

            ! Check dimensions and types.
            dims = [size(val)]
            call check_dataset_params(dataset, nm, 8, H5T_INTEGER_F, dims)

            ! And actually read the data. Note that we manipulate the pointer
            ! to be a 32bit integer, as that is always present in the library.
            call int32_pointer_abuse(val, ptr)
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

        integer(hid_t), intent(in) :: dataset, class_id
        character(*), intent(in) :: nm
        integer(hsize_t), intent(in) :: sz
        integer(hsize_t), intent(in) :: dims(:)
        character(*), parameter :: t_r = 'check_dataset_params'

        integer(hid_t) :: rank, type_id, ds_class
        
        integer(hid_t) :: ds_rank, dataspace, err
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

end module
