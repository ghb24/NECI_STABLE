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

end module
