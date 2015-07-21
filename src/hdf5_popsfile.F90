#include "macros.h"

module hdf5_popsfile

    ! Read and write POPSFILES using hdf5 as the data format
    !
    ! Data structure:
    !
    !   --> Don't assume that only walker-data is going to be relevant
    !   --> Create groups appropriately.
    !
    ! A: vcs_ver             - The SHA ID of the git commit
    ! A: compiled_at         - The time of code compilation
    ! A: compiled_on         - The date of code compilation
    ! A: config              - The build configuration
    ! A: status              - (opt) indicates if local changes beyond SHA ID
    !
    ! /system/               - Details of the system that is being restarted
    !
    ! /calculation/          - Details of the calculation
    !     /random_hash/      - Random values used in the orbital mapping
    ! 
    ! /wavefunction/         - Details of a determinental Hilbert space
    !     A: width           - Width of the bit-rep in 64-bit integers
    !
    !     /ilut/             - The bit representations of the determinants
    !     /sgns/             - The occupation of the determinants

    use Parallel_neci
    use constants
    use hdf5_util
    use util_mod
    use hdf5
    implicit none
    private

    ! Constants for naming various sections
    character(*), parameter :: &

            nm_vcs_ver = 'vcs_ver', &
            nm_comp_time = 'compiled_at', &
            nm_comp_date = 'compiled_on', &
            nm_date = 'date', &
            nm_seq_no = 'seq_no', &
            nm_comp_config = 'config', &
            nm_comp_status = 'status', &

            nm_calc_grp = 'calculation', &
            nm_random_hash = 'random_hash', &
            
            nm_wfn_grp = 'wavefunction', &
            nm_rep_width = 'width', &
            nm_sgn_len = 'lenof_sign', &
            nm_ilut = 'ilut', &
            nm_sgns = 'sgns'

    public :: write_popsfile_hdf5, read_popsfile_hdf5

contains

    subroutine write_popsfile_hdf5()

        integer(hid_t) :: plist_id, file_id, err

        ! Initialise the hdf5 fortran interface
        call h5open_f(err)

        ! Set up a property list to ensure file handling across all nodes.
        ! TODO: Check if we should be using a more specific communicator
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULl, err)

        ! TODO: Do sensible file handling here...
        call h5fcreate_f('popsfile.hdf5', H5F_ACC_TRUNC_F, file_id, err, &
                         access_prp=plist_id)
        call h5pclose_f(plist_id, err)

        call write_metadata(file_id)
        call write_calc_data(file_id)
        call write_walkers(file_id)

        ! And we are done!
        call h5fclose_f(file_id, err)
        call h5close_f(err)

    end subroutine write_popsfile_hdf5


    subroutine read_popsfile_hdf5(dets)

        ! Read a popsfile in, prior to running a new calculation

        ! n.b. This reads into the specified array, to allow use of popsfiles
        !      for initialising perturbations, etc.

        integer(int64), intent(out) :: dets(:, :)
        integer(hid_t) :: err, file_id, plist_id
        integer :: tmp

        write(6,*)
        write(6,*) '=='
        write(6,*) "Reading in from HDF5 popsfile"
        write(6,*) '======================================'

        ! Initialise the hdf5 fortran interface
        call h5open_f(err)

        ! Set up a property list to ensure file handling across all nodes.
        ! TODO: Check if we should be using a more specific communicator
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULl, err)

        ! Open the popsfile
        call h5fopen_f('popsfile.hdf5', H5F_ACC_RDONLY_F, file_id, err, &
                       access_prp=plist_id)

        call read_metadata(file_id)
        call read_calc_data(file_id)

        ! And we are done
        call h5pclose_f(plist_id, err)
        call h5fclose_f(file_id, err)
        call h5close_f(err)

    end subroutine


    subroutine write_metadata(parent)

        use CalcData, only: calc_seq_no

        ! Output macroscopic metadata applicable to a restart file, which may
        ! be used for establishing providence of calculations, etc.

        integer(hid_t), intent(in) :: parent
        character(19) :: date_str
        integer :: date_values(8)

        ! TODO: Run by
        call write_string_attribute(parent, nm_vcs_ver, _VCS_VER)
        call write_string_attribute(parent, nm_comp_date, __DATE__)
        call write_string_attribute(parent, nm_comp_time, __TIME__)
        call write_string_attribute(parent, nm_comp_config, _CONFIG)
#ifdef _WORKING_DIR_CHANGES
        call write_string_attribute(parent, nm_comp_status, &
                                    "Working directory contains local changes")
#endif

        ! How many calculations have been run to get to this popsfile pt?
        call write_int32_attribute(parent, nm_seq_no, int(calc_seq_no, int32))

        ! When are we running this?
        call date_and_time(values=date_values)

        write(date_str, '(i4,"-",i2,"-",i2," ",i2,":",i2,":",i2)') &
            date_values(1:3), date_values(5:7)
        call write_string_attribute(parent, nm_date, date_str)

    end subroutine

    subroutine read_metadata(parent)

        use CalcData, only: calc_seq_no

        ! Read in the macroscopic metadata applicable to the restart file.

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: err, attribute

        logical :: exists
        character(100) :: str_buf

        write(6,*) 'Previous calculation'

        call read_string_attribute(parent, nm_date, str_buf, exists)
        if (exists) write(6,*) 'Date: ', trim(str_buf)
        call read_int32_attribute(parent, nm_seq_no,  calc_seq_no, &
                                  default=1_int32)
        write(6,*) 'Sequence no.:', calc_seq_no

        ! Output nice details for usability
        call read_string_attribute(parent, nm_vcs_ver, str_buf, exists)
        if (exists) write(6,*) 'VCS ver: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_config, str_buf, exists)
        if (exists) write(6,*) 'Build config: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_status, str_buf, exists)
        if (exists) write(6,*) 'Build status: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_date, str_buf, exists)
        if (exists) write(6,*) 'Build date: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_time, str_buf, exists)
        if (exists) write(6,*) 'Build time: ', trim(str_buf)
        write(6,*)

        ! Update values for the new calculation
        calc_seq_no = calc_seq_no + 1

    end subroutine

    subroutine write_calc_data(parent)

        use load_balance_calcnodes, only: RandomOrbIndex
        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: calc_grp, err

        ! Firstly create the group for storing calculation-related data
        call h5gcreate_f(parent, nm_calc_grp, calc_grp, err)

        ! Write out the random orbital mapping index
        call write_int64_1d_dataset(calc_grp, nm_random_hash, RandomOrbIndex)

        ! Clear stuff up
        call h5gclose_f(calc_grp, err)

    end subroutine

    subroutine read_calc_data(parent)

        use load_balance_calcnodes, only: RandomOrbIndex
        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: grp_id, err

        call h5gopen_f(parent, nm_calc_grp, grp_id, err)

        ! Read out the random orbital mapping index
        call read_int64_1d_dataset(grp_id, nm_random_hash, RandomOrbIndex, &
                                   required=.true.)

        call h5gclose_f(grp_id, err)

    end subroutine

    subroutine write_walkers(parent)

        use iso_c_hack
        use bit_rep_data, only: NIfD, NIfTot, NOffSgn
        use FciMCData, only: AllTotWalkers, CurrentDets, MaxWalkersPart, &
                             TotWalkers
        use CalcData, only: tUseRealCoeffs

        ! Output the wavefunction information to the relevant groups in the
        ! wavefunction.

        integer(hid_t), intent(in) :: parent
        type(c_ptr) :: cptr
        integer(int32), pointer :: ptr(:)
        integer(int32) :: boop

        character(*), parameter :: t_r = 'write_walkers'

        integer(hid_t) :: wfn_grp_id, dataspace, dataset, err, memspace
        integer(hid_t) :: plist_id

        integer(hsize_t) :: counts(0:nProcessors-1)
        integer(hsize_t) :: all_count

        integer(int32) :: bit_rep_width
        integer(hsize_t) :: mem_offset(2), write_offset(2)
        integer(hsize_t) :: dims(2), hyperdims(2)
        integer :: ierr

        ! TODO: Add a (slower) fallback routine for weird cases, odd HDF libs

        ! Firstly create the group for storing wavefunction info
        call h5gcreate_f(parent, nm_wfn_grp, wfn_grp_id, err)

        ! TODO: Refactor these chunks into their own little subroutines.
        ! We fix the format of the binary file. Thus if we are on a 32-bit
        ! build, we need to convert the data into 64-bit compatibile chunks.
        if (build_64bit) then
            bit_rep_width = NIfD + 1
        else
            bit_rep_width = 2 * (NIfD + 1)
            call stop_all(t_r, "Needs manual, careful, testing")
        end if

        ! Output the bit-representation data
        call write_int32_attribute(wfn_grp_id, nm_rep_width, bit_rep_width)
        call write_int32_attribute(wfn_grp_id, nm_sgn_len, &
                                   int(lenof_sign, int32))

        ! How many occuiped determinants are there on each of the processors
        call MPIAllGather(TotWalkers, counts, ierr)
        all_count = sum(counts)
        write_offset = [0_hsize_t, sum(counts(0:iProcIndex-1))]

        ! Write out the determinant bit-representations
        call write_2d_multi_arr_chunk_offset( &
                wfn_grp_id, nm_ilut, H5T_NATIVE_INTEGER_8, &
                arr_2d_ptr(CurrentDets), arr_2d_dims(CurrentDets), &
                [int(nifd+1, hsize_t), int(TotWalkers, hsize_t)], & ! dims
                [0_hsize_t, 0_hsize_t], & ! offset
                [int(nifd+1, hsize_t), all_count], & ! all dims
                [0_hsize_t, sum(counts(0:iProcIndex-1))] & ! output offset
        )

        ! Write out the sign values on each of the processors
        if (.not. tUseRealCoeffs) &
            call stop_all(t_r, "This could go badly...")

        call write_2d_multi_arr_chunk_offset( &
                wfn_grp_id, nm_sgns, H5T_NATIVE_REAL_8, &
                arr_2d_ptr(CurrentDets), arr_2d_dims(CurrentDets), &
                [int(lenof_sign, hsize_t), int(TotWalkers, hsize_t)], & ! dims
                [int(nOffSgn, hsize_t), 0_hsize_t], & ! offset
                [int(lenof_sign, hsize_t), all_count], & ! all dims
                [0_hsize_t, sum(counts(0:iProcIndex-1))] & ! output offset
        )

        ! And we are done
        call h5gclose_f(wfn_grp_id, err)
        


    end subroutine write_walkers

end module
