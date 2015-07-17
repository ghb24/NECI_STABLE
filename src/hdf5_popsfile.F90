#include "macros.h"

module hdf5_popsfile

    ! Read and write POPSFILES using hdf5 as the data format
    !
    ! Data structure:
    !
    !   --> Don't assume that only walker-data is going to be relevant
    !   --> Create groups appropriately.
    !
    ! /system/                - Details of the system that is being restarted
    !
    ! /calculation/           - Details of the calculation
    ! 
    ! /determinants/A         - Details of a determinental Hilbert space
    !     A: 

    use constants
    use hdf5_util
    use hdf5
    implicit none
    private

    ! Constants for naming various sections
    character(*), parameter :: &
            
            wfn_grp = 'wavefunction'

    public :: write_popsfile_hdf5, read_popsfile_hdf5

contains

    subroutine write_popsfile_hdf5()

        integer(hid_t) :: file_id, err

        ! Initialise the hdf5 fortran interface
        call h5open_f(err)

        ! TODO: Do sensible file handling here...
        call h5fcreate_f('popsfile.hdf5', H5F_ACC_TRUNC_F, file_id, err)

        call write_walkers(file_id)

        ! And we are done!
        call h5fclose_f(file_id, err)
        call h5close_f(err)

    end subroutine write_popsfile_hdf5


    subroutine read_popsfile_hdf5()

    end subroutine


    subroutine write_walkers(parent)

        ! Output the wavefunction information to the relevant groups in the
        ! wavefunction.

        integer(hid_t), intent(in) :: parent

        integer(hid_t) :: wfn_grp_id, err


        ! Firstly create the group for storing wavefunction info
        call h5gcreate_f(parent, wfn_grp, wfn_grp_id, err)


        ! And we are done
        call h5gclose_f(wfn_grp_id, err)
        


    end subroutine write_walkers

end module
