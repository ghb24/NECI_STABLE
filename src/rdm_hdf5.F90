#include "macros.h"

module rdm_hdf5

    use MPI_wrapper
    use Parallel_neci
    use constants
#ifdef USE_HDF_
    use parallel_hdf5_utils
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
      call write_data_phdf5(values, 'values', grp_id)
      call write_data_phdf5(indices, 'indices', grp_id)
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

      call h5gcreate_f(parent, '2200', grp_id, err)
      call write_data_phdf5(values, 'values', grp_id)
      call write_data_phdf5(indices, 'indices', grp_id)
      call h5gclose_f(grp_id, err)
    end subroutine write_2rdm_hdf5

end module rdm_hdf5
