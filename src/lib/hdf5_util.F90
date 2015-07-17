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

    use constants
    use hdf5
    implicit none

contains

    subroutine write_int32_attribute(parent, nm, val)

        integer(hid_t), intent(inout) :: parent
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

end module
