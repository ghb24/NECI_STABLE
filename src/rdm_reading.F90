module rdm_reading

    ! Routines to read in various RDMS (1-RDM or 2-RDMs POPSFILES or spinfree
    ! 2-RDMs). Also, there are some routines to calculate and print 1-RDMs,
    ! designed to be used directly after reading in 2-RDMs, for cases where
    ! the user forgot to print them out.

    use constants

    implicit none

contains

    subroutine read_1rdm(rdm_defs, one_rdm, irdm)

        ! Read a 1-RDM POPSFILE file (with filename stem 'OneRDM_POPS.'), with
        ! label irdm. Copy it into the one_Rdm object.

        use rdm_data, only: one_rdm_t, rdm_definitions_t
        use RotateOrbsData, only: SymLabelListInv_rot
        use util_mod, only: get_free_unit, int_fmt

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(one_rdm_t), intent(inout) :: one_rdm
        integer, intent(in) :: irdm

        integer :: one_rdm_unit, file_end, i, j
        real(dp) :: rdm_sign
        logical :: file_exists
        character(20) :: filename
        character(len=*), parameter :: t_r = 'read_1rdm'

        write(6,'(1X,"Reading in the 1-RDMs...")')

        associate(state_labels => rdm_defs%state_labels, repeat_label => rdm_defs%repeat_label)
            if (state_labels(1,irdm) == state_labels(2,irdm)) then
                write(filename, '("OneRDM_POPS.",'//int_fmt(state_labels(1,irdm),0)//')') irdm
            else
                write(filename, '("OneRDM_POPS.",'//int_fmt(state_labels(1,irdm),0)//',"_",'&
                                                  //int_fmt(state_labels(2,irdm),0)//',".",i1)') &
                                    state_labels(1,irdm), state_labels(2,irdm), repeat_label(irdm)
            end if
        end associate

        inquire(file=trim(filename), exist=file_exists)

        if (.not. file_exists) then
            call stop_all(t_r, "Attempting to read in the 1-RDM from "//trim(filename)//", but this file does not exist.")
        end if

        one_rdm_unit = get_free_unit()
        open(one_rdm_unit, file=trim(filename), status='old', form='unformatted')

        do
            read(one_rdm_unit, iostat=file_end) i, j, rdm_sign
            if (file_end > 0) call stop_all(t_r, "Error reading "//trim(filename)//".")
            ! file_end < 0 => end of file reached.
            if (file_end < 0) exit

            one_rdm%matrix(SymLabelListInv_rot(i), SymLabelListInv_rot(j)) = rdm_sign
        end do

        close(one_rdm_unit)

    end subroutine read_1rdm

    subroutine read_2rdm_popsfile(rdm, spawn)

        ! Read 2-RDMs into the rdm object, from an RDM_POPSFILE file.
        ! The spawn object is used to perform parallel communication of
        ! elements to their correct process. Already is done by the root
        ! process. The spawn object is cleared at the end of this routine.

        use hash, only: clear_hash_table, FindWalkerHash, add_hash_table_entry
        use Parallel_neci, only: iProcIndex, nProcessors
        use rdm_data, only: rdm_list_t, rdm_spawn_t
        use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm, add_to_rdm_spawn_t
        use rdm_data_utils, only: communicate_rdm_spawn_t_wrapper, annihilate_rdm_list
        use util_mod, only: get_free_unit

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn

        integer(int_rdm) :: ijkl, rdm_entry(0:rdm%sign_length)
        integer :: ij, kl, i, j, k, l, ij_proc_row, sign_length_old
        integer :: ielem, pops_unit, file_end, hash_val
        real(dp) :: rdm_sign(rdm%sign_length)
        logical :: file_exists, nearly_full, finished, all_finished
        character(len=*), parameter :: t_r = 'read_2rdm_popsfile'

        write(6,'(1X,"Reading in the 2-RDMs...")')

        ! If we're about to fill up the spawn list, perform a communication.
        nearly_full = .false.
        ! Have we finished adding RDM elements to the spawned list?
        finished = .false.

        inquire(file='RDM_POPSFILE', exist=file_exists)

        if (.not. file_exists) then
            call stop_all(t_r, "Attempting to read in the 2-RDM from RDM_POPSFILE, but this file does not exist.")
        end if

        ! Only let the root processor do the reading in.
        if (iProcIndex == 0) then
            pops_unit = get_free_unit()
            open(pops_unit, file='RDM_POPSFILE', status='old', form='unformatted')

            ! Read in the first line, which holds sign_length for the
            ! printed RDM - check it is consistent.
            read(pops_unit, iostat=file_end) sign_length_old
            if (sign_length_old /= rdm%sign_length) then
                call stop_all(t_r, "Error reading RDM_POPSFILE - the number of RDMs printed in the &
                                   &popsfile is different to the number to be sampled.")
            else if (file_end > 0) then
                call stop_all(t_r, "Error reading the first line of RDM_POPSFILE.")
            end if

            do
                read(pops_unit, iostat=file_end) rdm_entry
                if (file_end > 0) call stop_all(t_r, "Error reading RDM_POPSFILE.")
                ! file_end < 0 => end of file reached.
                if (file_end < 0) exit

                ! If the spawned list is nearly full, perform a communication.
                if (nearly_full) then
                    call communicate_rdm_spawn_t_wrapper(spawn, rdm, finished, all_finished)
                    nearly_full = .false.
                end if

                call calc_separate_rdm_labels(rdm_entry(0), ij, kl, i, j, k, l)
                call extract_sign_rdm(rdm_entry, rdm_sign)

                call add_to_rdm_spawn_t(spawn, i, j, k, l, rdm_sign, .false., nearly_full)
            end do

            close(pops_unit)
        end if

        finished = .true.
        ! Keep performing communications until all RDM spawnings on every
        ! processor have been communicated.
        do
            call communicate_rdm_spawn_t_wrapper(spawn, rdm, finished, all_finished)
            if (all_finished) exit
        end do

        call annihilate_rdm_list(rdm)

        ! Fill in the hash table to the RDM. Clear it first, just in case.
        call clear_hash_table(rdm%hash_table)
        do ielem = 1, rdm%nelements
            call calc_separate_rdm_labels(rdm%elements(0,ielem), ij, kl, i, j, k, l)
            hash_val = FindWalkerHash((/i,j,k,l/), size(rdm%hash_table))
            call add_hash_table_entry(rdm%hash_table, ielem, hash_val)
        end do

        ! Clear the spawn object.
        spawn%free_slots = spawn%init_free_slots(0:nProcessors-1)
        call clear_hash_table(spawn%rdm_send%hash_table)

    end subroutine read_2rdm_popsfile

    subroutine read_spinfree_2rdm_files(rdm_defs, rdm, spawn)

        ! Read spinfree 2-RDMs from standard NECI spinfree_TwoRDM files.

        ! The resulting object will be in the form of a spinfree 2-RDM, *not*
        ! of a standard NECI RDM, so care should be taken - one can not expect
        ! most routines to act correctly on this resulting object. For example,
        ! only spatial orbital labels or stored, symmetry is not taken account
        ! of (so i<j and k<l is not enforced) and orbital labels are in a
        ! different order to that of traditional NECI 2-RDMs.

        use hash, only: clear_hash_table, FindWalkerHash, add_hash_table_entry
        use Parallel_neci, only: iProcIndex, nProcessors
        use rdm_data, only: rdm_list_t, rdm_spawn_t
        use rdm_data, only: rdm_definitions_t
        use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm, add_to_rdm_spawn_t
        use rdm_data_utils, only: communicate_rdm_spawn_t_wrapper, annihilate_rdm_list
        use util_mod, only: get_free_unit, int_fmt

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn

        integer(int_rdm) :: ijkl, rdm_entry(0:rdm%sign_length)
        integer :: ij, kl, i, j, k, l, ij_proc_row, sign_length_old
        integer :: irdm, ielem, iunit(rdm%sign_length), ierr, hash_val
        real(dp) :: rdm_sign(rdm%sign_length)
        logical :: file_exists, nearly_full, finished, all_finished
        logical :: file_end(rdm%sign_length)
        character(30) :: rdm_filename(rdm%sign_length)
        character(len=*), parameter :: t_r = 'read_2rdm_popsfile'

        write(6,'(1X,"Reading in the spinfree 2-RDMs...")')

        ! If we're about to fill up the spawn list, perform a communication.
        nearly_full = .false.
        ! Have we finished adding RDM elements to the spawned list?
        finished = .false.
        ! Have we reached the end of the various files?
        file_end = .false.

        associate(state_labels => rdm_defs%state_labels, repeat_label => rdm_defs%repeat_label)
            do irdm = 1, rdm%sign_length
                if (state_labels(1,irdm) == state_labels(2,irdm)) then
                    write(rdm_filename(irdm), '("spinfree_TwoRDM.",'//int_fmt(state_labels(1,irdm),0)//')') irdm
                else
                    write(rdm_filename(irdm), '("spinfree_TwoRDM.",'//int_fmt(state_labels(1,irdm),0)//',"_",'&
                                                                    //int_fmt(state_labels(2,irdm),0)//',".",i1)') &
                                        state_labels(1,irdm), state_labels(2,irdm), repeat_label(irdm)
                end if
                inquire(file=trim(rdm_filename(irdm)), exist=file_exists)

                if (.not. file_exists) then
                    call stop_all(t_r, "Attempting to read in a spinfree 2-RDM from "//trim(rdm_filename(irdm))//", &
                                       &but this file does not exist.")
                end if
            end do
        end associate

        ! Only let the root processor do the reading in.
        if (iProcIndex == 0) then
            do irdm = 1, rdm%sign_length
                iunit(irdm) = get_free_unit()
                open(iunit(irdm), file=trim(rdm_filename(irdm)), status='old')
            end do

            do
                do irdm = 1, rdm%sign_length
                    if (.not. file_end(irdm)) then
                        rdm_sign = 0.0_dp
                        read(iunit(irdm), "(4I15, F30.20)", iostat=ierr) i, j, k, l, rdm_sign(irdm)
                        if (ierr > 0) call stop_all(t_r, "Error reading "//trim(rdm_filename(irdm))//".")
                        ! The final line is printed as follows, by convention.
                        if (i == -1 .and. j == -1 .and. k == -1 .and. l == -1) then
                            file_end(irdm) = .true.
                            exit
                        end if

                        ! If the spawned list is nearly full, perform a communication.
                        if (nearly_full) then
                            call communicate_rdm_spawn_t_wrapper(spawn, rdm, finished, all_finished)
                            nearly_full = .false.
                        end if

                        call add_to_rdm_spawn_t(spawn, i, j, k, l, rdm_sign, .true., nearly_full)
                    end if
                end do
                if (all(file_end)) exit
            end do

            do irdm = 1, rdm%sign_length
                close(iunit(irdm))
            end do
        end if

        finished = .true.
        ! Keep performing communications until all RDM spawnings on every
        ! processor have been communicated.
        do
            call communicate_rdm_spawn_t_wrapper(spawn, rdm, finished, all_finished)
            if (all_finished) exit
        end do

        call annihilate_rdm_list(rdm)

        ! Fill in the hash table to the RDM. Clear it first, just in case.
        call clear_hash_table(rdm%hash_table)
        do ielem = 1, rdm%nelements
            call calc_separate_rdm_labels(rdm%elements(0,ielem), ij, kl, i, j, k, l)
            hash_val = FindWalkerHash((/i,j,k,l/), size(rdm%hash_table))
            call add_hash_table_entry(rdm%hash_table, ielem, hash_val)
        end do

        ! Clear the spawn object.
        spawn%free_slots = spawn%init_free_slots(0:nProcessors-1)
        call clear_hash_table(spawn%rdm_send%hash_table)

    end subroutine read_spinfree_2rdm_files

    ! Routines to calculate and print 1-RDMs directly after reading in 2-RDMs.

    subroutine print_1rdms_from_2rdms_wrapper(rdm_defs, one_rdms, two_rdms, open_shell)

        ! Wrapper function to calculate and print 1-RDMs from 2-RDMs, as might
        ! be useful if the user forgot to print 1-RDMs in a calculation.

        use Parallel_neci, only: MPISumAll, iProcIndex
        use rdm_data, only: rdm_definitions_t, one_rdm_t, rdm_list_t
        use rdm_estimators, only: calc_rdm_trace
        use rdm_finalising, only: calc_1rdms_from_2rdms, write_1rdm, finalise_1e_rdm
        use SystemData, only: nel

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        type(rdm_list_t), intent(in) :: two_rdms
        logical, intent(in) :: open_shell

        integer :: irdm
        real(dp) :: rdm_trace(two_rdms%sign_length), rdm_trace_all(two_rdms%sign_length)
        real(dp) :: rdm_norm_all(two_rdms%sign_length), norm_1rdm(size(one_rdms))
        character(len=*), parameter :: t_r = 'print_1rdms_from_2rdms_wrapper'

        if (size(one_rdms) == 0) then
            call stop_all(t_r, "You have asked to print 1-RDMs but they are not allocated. Make sure &
                                &that you have asked for both 1-RDMs and 2-RDMs to be calculated by &
                                &seting the RDMExcitLevel to 3, i.e. 'CALCRDMONFLY 3 ...' in input options.")
        end if

        !call calc_rdm_trace(two_rdms, rdm_trace)
        !call MPISumAll(rdm_trace, rdm_trace_all)
        ! RDMs are normalised so that their trace is nel*(nel-1)/2.
        !rdm_norm_all = rdm_trace_all*2.0_dp/(nel*(nel-1))

        ! The read-in spinfree 2-RDMs should already be normalised.
        rdm_norm_all = 1.0_dp
        norm_1rdm = 1.0_dp

        call calc_1rdms_from_2rdms(rdm_defs, one_rdms, two_rdms, rdm_norm_all, open_shell)

        ! we need to finalise the 1rdms
        call finalise_1e_rdm(rdm_defs, one_rdms, norm_1rdm, .false.)

        if(iProcIndex==0) then
           do irdm = 1, size(one_rdms)
              call write_1rdm(rdm_defs, one_rdms(irdm)%matrix, irdm, norm_1rdm(irdm), .true., .false.)
           end do
        endif

    end subroutine print_1rdms_from_2rdms_wrapper

    subroutine print_1rdms_from_sf2rdms_wrapper(rdm_defs, one_rdms, two_rdms)

        ! Wrapper routine to calculate and print spinfree 1-RDMs from
        ! spinfree 2-RDMs, as read in directly from a spinfree_TwoRDM file
        ! which was output by NECI. This routine will *not* work on standard
        ! format NECI RDM objects, nor is it capable of outputting spinned
        ! 1-RDMs in open shell systems, only spinfree 1-RDMs.

        use Parallel_neci, only: MPISumAll
        use rdm_data, only: rdm_definitions_t, one_rdm_t, rdm_list_t
        use rdm_finalising, only: calc_1rdms_from_spinfree_2rdms, write_1rdm

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        type(rdm_list_t), intent(in) :: two_rdms

        integer :: irdm
        real(dp) :: rdm_norm_all(two_rdms%sign_length), norm_1rdm
        character(len=*), parameter :: t_r = 'print_1rdms_from_sf2rdms_wrapper'

        if (size(one_rdms) == 0) then
            call stop_all(t_r, "You have asked to print 1-RDMs but they are not allocated. Make sure &
                                &that you have asked for both 1-RDMs and 2-RDMs to be calculated by &
                                &seting the RDMExcitLevel to 3, i.e. 'CALCRDMONFLY 3 ...' in input options.")
        end if

        ! The read-in spinfree 2-RDMs should already be normalised.
        rdm_norm_all = 1.0_dp

        call calc_1rdms_from_spinfree_2rdms(one_rdms, two_rdms, rdm_norm_all)

        do irdm = 1, size(one_rdms)
            call write_1rdm(rdm_defs, one_rdms(irdm)%matrix, irdm, 1.0_dp, .true., .false.)
        end do

    end subroutine print_1rdms_from_sf2rdms_wrapper

end module rdm_reading
