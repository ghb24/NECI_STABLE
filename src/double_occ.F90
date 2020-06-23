#include "macros.h"

! make on module file for now which stores all the relevant stuff
! for the double occupancy measurements
module double_occ_mod

    use SystemData, only: nel, nbasis
    use bit_rep_data, only: nifd, niftot
    use constants, only: n_int, lenof_sign, write_state_t, dp, int_rdm, inum_runs
    use ParallelHelper, only: iProcIndex, root
    use CalcData, only: tReadPops, StepsSft
    use LoggingData, only: tMCOutput, t_calc_double_occ_av, t_spin_measurements
    use util_mod
    use FciMCData, only: iter, PreviousCycles, norm_psi, totwalkers, all_norm_psi_squared
    use rdm_data, only: rdm_list_t
    use Parallel_neci, only: iProcIndex, nProcessors, MPISumAll, MPIAllreduce
    use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm, &
                              calc_combined_rdm_label
    use UMatCache, only: spatial
    use sort_mod, only: sort
    use FciMCData, only: fcimc_iter_data

    implicit none

    ! data storage part:
    ! i need local and global storage containers for <n_u n_d>
    ! and maybe also instantaneous, averaged, summed over quantities

    ! instantaneous quantities: they are just used to see the progress of
    ! the double occupancy, they are reset every iteration
    real(dp) :: inst_double_occ = 0.0_dp, all_inst_double_occ = 0.0_dp

    ! summed over quantities: after the shift has set in this quantitiy
    ! is used to calculate the double occupancy of the converged WF
    ! i could also print this at every step over the summed norm to see the
    ! convergence progress of this quantitiy
    real(dp) :: sum_double_occ = 0.0_dp, all_sum_double_occ = 0.0_dp

    ! it seems there is no record of the averaged norm hm..
    real(dp) :: sum_norm_psi_squared = 0.0_dp

    real(dp), allocatable :: spin_up_occ(:), spin_down_occ(:), spin_diff(:), &
                             double_occ_vec(:)

    real(dp), allocatable :: inst_spin_diff(:), all_inst_spin_diff(:)
    real(dp), allocatable :: inst_spatial_doub_occ(:), all_inst_spatial_doub_occ(:), &
                             sum_double_occ_vec(:), sum_spin_diff(:)

contains

    subroutine init_spin_measurements()
        ! routine to initialize spin measurement vectors
        character(*), parameter :: this_routine = "init_spin_measurements"
        integer :: ierr

        if (allocated(spin_up_occ))               deallocate(spin_up_occ)
        if (allocated(spin_down_occ))             deallocate(spin_down_occ)
        if (allocated(spin_diff))                 deallocate(spin_diff)
        if (allocated(double_occ_vec))            deallocate(double_occ_vec)
        if (allocated(inst_spin_diff))            deallocate(inst_spin_diff)
        if (allocated(all_inst_spin_diff))        deallocate(all_inst_spin_diff)
        if (allocated(inst_spatial_doub_occ))     deallocate(inst_spatial_doub_occ)
        if (allocated(all_inst_spatial_doub_occ)) deallocate(all_inst_spatial_doub_occ)
        if (allocated(sum_double_occ_vec))        deallocate(sum_double_occ_vec)
        if (allocated(sum_spin_diff))             deallocate(sum_spin_diff)

        allocate(spin_up_occ(nbasis/2))
        allocate(spin_down_occ(nBasis/2))
        allocate(spin_diff(nBasis/2))
        allocate(double_occ_vec(nBasis/2))
        allocate(inst_spin_diff(nBasis/2))
        allocate(all_inst_spin_diff(nBasis/2))
        allocate(inst_spatial_doub_occ(nBasis/2))
        allocate(all_inst_spatial_doub_occ(nBasis/2))
        allocate(sum_double_occ_vec(nBasis/2))
        allocate(sum_spin_diff(nBasis/2))

        spin_up_occ = 0.0_dp
        spin_down_occ = 0.0_dp
        spin_diff = 0.0_dp
        double_occ_vec = 0.0_dp
        inst_spin_diff = 0.0_dp
        all_inst_spin_diff = 0.0_dp
        inst_spatial_doub_occ = 0.0_dp
        all_inst_spatial_doub_occ = 0.0_dp
        sum_double_occ_vec = 0.0_dp
        sum_spin_diff = 0.0_dp

    end subroutine init_spin_measurements

    subroutine deallocate_spin_measurements()
        ! deallocate the spin occupation measurement vectors
        character(*), parameter :: this_routine = "deallocate_spin_measurements"

        if (allocated(spin_up_occ))        deallocate(spin_up_occ)
        if (allocated(spin_down_occ))      deallocate(spin_down_occ)
        if (allocated(spin_diff))          deallocate(spin_diff)
        if (allocated(double_occ_vec))     deallocate(double_occ_vec)
        if (allocated(inst_spin_diff))     deallocate(inst_spin_diff)
        if (allocated(all_inst_spin_diff)) deallocate(all_inst_spin_diff)
        if (allocated(inst_spatial_doub_occ)) deallocate(inst_spatial_doub_occ)
        if (allocated(all_inst_spatial_doub_occ)) deallocate(all_inst_spatial_doub_occ)
        if (allocated(sum_double_occ_vec)) deallocate(sum_double_occ_vec)

    end subroutine deallocate_spin_measurements

    subroutine measure_double_occ_and_spin_diff(ilut, nI, real_sgn)
        ! routine to measure double occupancy and spin difference for each
        ! orbital
        use UMatCache, only: gtid
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: ni(nel)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        character(*), parameter :: this_routine = "measure_double_occ_and_spin_diff"

        integer :: i, spin_orb, spat_orb(nel)
        real(dp) :: contrib

        ASSERT(allocated(spin_diff))
        ASSERT(allocated(double_occ_vec))

        ! i can calculate the C_i^2 at the beginning already since it is
        ! always the same and i guess i have atleast on contribution atleast
        ! for each occupied orbital
#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
        ! i do not want to deal with complex runs for now..
        call stop_all(this_routine, &
            "complex double occupancy measurement not yet implemented!")
#else
        contrib = real_sgn(1) * real_sgn(2)
#endif
#else
        contrib = abs(real_sgn(1))**2
#endif

        spat_orb = gtid(nI)

        i = 1
        do while (i < nel + 1)
            spin_orb = nI(i)
            if (is_beta(spin_orb)) then
                ! check if it is a doubly occupied orb
                if (IsDoub(ilut,spin_orb)) then
                    ! then we want to add to the double_occ vector

                    inst_spatial_doub_occ(spat_orb(i)) = &
                        inst_spatial_doub_occ(spat_orb(i)) + contrib

                    ! and we can skip the even alpha orbital in nI
                    i = i + 2
                else
                    ! beta spin contributes negatively!
                    inst_spin_diff(spat_orb(i)) = inst_spin_diff(spat_orb(i)) &
                                                  - contrib

                    i = i + 1
                end if

            else
                ! the way i plan to set it up, we check beta spins in the
                ! same orbital first.. so it can't be doubly occupied at
                ! this point!
                inst_spin_diff(spat_orb(i)) = inst_spin_diff(spat_orb(i)) &
                                              + contrib

                i = i + 1
            end if
        end do

        ! this should be it or?
    end subroutine measure_double_occ_and_spin_diff

    subroutine finalize_double_occ_and_spin_diff()
        ! routine to communicate all the data from all nodes and also
        ! print out the information i guess..
        character(*), parameter :: this_routine = "finalize_double_occ_and_spin_diff"

        real(dp), allocatable :: all_double_occ_vec(:), all_spin_diff(:), &
                                 all_sum_double_occ_vec(:), all_sum_spin_diff(:)


        allocate(all_double_occ_vec(nBasis/2))
        allocate(all_spin_diff(nBasis/2))
        allocate(all_sum_double_occ_vec(nBasis/2))
        allocate(all_sum_spin_diff(nBasis/2))

        all_double_occ_vec = 0.0_dp
        all_spin_diff = 0.0_dp
        all_sum_double_occ_vec = 0.0_dp
        all_sum_spin_diff = 0.0_dp

        call MPIAllreduce(spin_diff, MPI_SUM, all_spin_diff)
        call MPIAllreduce(double_occ_vec, MPI_SUM, all_double_occ_vec)
        call MPIAllreduce(sum_double_occ_vec, MPI_SUM, all_sum_double_occ_vec)
        call MPIAllreduce(sum_spin_diff, MPI_SUM, all_sum_spin_diff)

        if (iProcIndex == Root) then
            ! first ouput the double occupancy and spin difference for each
            ! spatial orbital
            print *, "double occupancy for each orbital: "
            print *, sum_double_occ_vec / (sum_norm_psi_squared * real(StepsSft,dp))

            print *, "spin difference for each orbital: "
            print *, sum_spin_diff / (sum_norm_psi_squared * real(StepsSft,dp))

            ! and also calculate the summed value over all orbitals to
            ! compare it with the other method and two_rdms

            print *, "total double double occupancy: "
            print *, 2.0_dp * sum(sum_double_occ_vec) / &
                (sum_norm_psi_squared * real(nBasis,dp) * real(StepsSft,dp))

            print *, "total spin difference (should have to do with Ms!)"
            print *, 2.0_dp * sum(sum_spin_diff) / &
                (sum_norm_psi_squared * real(nbasis,dp) * real(StepsSft,dp))

        end if

        call deallocate_spin_measurements()

        deallocate(all_double_occ_vec)
        deallocate(all_spin_diff)
        deallocate(all_sum_double_occ_vec)

    end subroutine finalize_double_occ_and_spin_diff

    subroutine measure_spin_occupation_ilut(ilut, real_sgn)
        ! routine to measure the spin occupation for each spatial orbital
        ! individually
        ! for now it is done in a really slow fashion, if this ends up
        ! to be interesting it should be optimized
        ! this is done for ilut inputs..
        use DetBitOps, only: MaskBeta, MaskAlpha
        integer(n_int), intent(in) :: ilut(0:niftot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        character(*), parameter :: this_routine = "measure_spin_occupation_ilut"

        integer(n_int) :: alpha(0:niftot), beta(0:niftot)
        integer :: i, x, y

        ! first check the alpha and beta occupations:
        alpha = iand(ilut, MaskAlpha)
        beta = iand(ilut, MaskBeta)

        ! i have to take care about how many integers are used to encode
        ! the orbitals
        do i = 1, nBasis/2
            ! i have to determine how to convert (i) to the correct
            ! index in the ilut..
            ! x = which integer?
            ! y = which position in integer?
            if (btest(alpha(x),y)) then
#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
            call stop_all(this_routine, &
                "complex double occupancy measurement not yet implemented!")
#else
                spin_up_occ(i) = spin_up_occ(i) + real_sgn(1) * real_sgn(2)
#endif
#else
                spin_up_occ(i) = spin_up_occ(i) + abs(real_sgn(1))**2
#endif
            end if

            if (btest(beta(x),y)) then
#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
            call stop_all(this_routine, &
                "complex double occupancy measurement not yet implemented!")

#else
                spin_down_occ(i) = spin_down_occ(i) + real_sgn(1) * real_sgn(2)
#endif
#else
                spin_down_occ(i) = spin_down_occ(i) + abs(real_sgn(1))**2
#endif
            end if
        end do

    end subroutine measure_spin_occupation_ilut

    subroutine measure_spin_occupation_nI(nI, real_sgn)
        ! routine same as above, but for nI input of the determinant
        use UMatCache, only: gtid
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        character(*), parameter :: this_routine = "measure_spin_occupation_nI"
        integer :: i, spin_orb, spat_orb

        ! in this routine i loop over the nI and check if and which spins
        ! are occupied

        do i = 1, nel
            spin_orb = nI(i)
            spat_orb = gtid(spin_orb)
            if (is_alpha(spin_orb)) then
#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
            call stop_all(this_routine, &
                "complex double occupancy measurement not yet implemented!")
#else
                spin_up_occ(spat_orb) = spin_up_occ(spat_orb) + real_sgn(1) * real_sgn(2)
#endif
#else
                spin_up_occ(spat_orb) = spin_up_occ(spat_orb) + abs(real_sgn(1))**2
#endif
            end if
            if (is_beta(spin_orb)) then
#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
            call stop_all(this_routine, &
                "complex double occupancy measurement not yet implemented!")
#else
                spin_down_occ(spat_orb) = spin_down_occ(spat_orb) + real_sgn(1) * real_sgn(2)
#endif
#else
                spin_down_occ(spat_orb) = spin_up_occ(spat_orb) + abs(real_sgn(1))**2
#endif
            end if
        end do

    end subroutine measure_spin_occupation_nI


    function count_double_orbs(ilut) result(n_double_orbs)
        ! returns the number of doubly occupied orbitals for a given
        ! determinant.
        use DetBitOps, only: count_open_orbs
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer :: n_double_orbs
        character(*), parameter :: this_routine = "count_double_occupancy"

        integer :: n_open_orbs

        n_open_orbs = count_open_orbs(ilut(0:nifd))

        ! the doubly occupied orbitals are just the number of electrons
        ! minus the number of open orbitals divided by 2
        n_double_orbs = (nel - n_open_orbs)/2

    end function count_double_orbs


    function get_double_occupancy(ilut, real_sgn) result(double_occ)
        ! function to get the contribution to the double occupancy for a
        ! given determinant, by calculating the
        ! |C_I|^2 *\sum_i <n_iu n_id>/n_spatial_orbs
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        real(dp) :: double_occ
        character(*), parameter :: this_routine = "get_double_occupancy"

        integer :: n_double_orbs
        real(dp) :: frac_double_orbs
        integer(n_int) :: sgn(lenof_sign)

#ifdef CMPLX_
        complex(dp) :: complex_sgn
#endif

        ! here i need to be careful, if it is a double or single run and
        ! stuff
        n_double_orbs = count_double_orbs(ilut(0:nifd))

        ! i only need the fraction of doubly occupied orbitals out of all
        ! spatial orbitals
        frac_double_orbs = 2.0_dp * real(n_double_orbs,dp) / real(nbasis,dp)

        ! now i want to sum in the walkers from the runs..
        ! todo: have to figure out how to access the different runs

        ! extract the walker occupation

        ! do i want to do that for complex walkers also?? i guess so..
        ! to get it running do it only for  single run for now!
        ! do not do the division here, but only in the output!
#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
        call stop_all(this_routine, &
            "complex double occupancy measurement not yet implemented!")
        ! i in the case of complex runs i guess that the double_occ
        ! can be complex too.. atleast in the case of doubleruns..
#else
        ! i essentially only need two runs!
        double_occ = real_sgn(1) * real_sgn(2) * frac_double_orbs
#endif
#else
#ifdef CMPLX_
        call stop_all(this_routine, &
            "complex double occupancy measurement not yet implemented!")
#endif
        double_occ = real_sgn(1)**2 * frac_double_orbs
#endif

    end function get_double_occupancy

    subroutine rezero_double_occ_stats
        ! at the end of each cycle i should rezero the non-summed and
        ! non-averaged quantities to 0
        character(*), parameter :: this_routine = "rezero_double_occ_stats"

        ! there is more stuff to come..
        inst_double_occ = 0.0_dp

        if (t_spin_measurements) then
            inst_spatial_doub_occ = 0.0_dp
        end if

    end subroutine rezero_double_occ_stats

    subroutine rezero_spin_diff()
        character(*), parameter :: this_routine = "rezero_spin_diff"

        inst_spin_diff = 0.0_dp

    end subroutine rezero_spin_diff

    subroutine write_spin_diff_stats(initial)
        ! routine to print out the instant spin-difference to a file
        ! routine to print out the instant spin-difference to a file
        logical, intent(in), optional :: initial
        character(*), parameter :: this_routine = "write_spin_diff_stats"

        type(write_state_t), save :: state
        logical, save :: inited = .false.
        integer :: i
        character(12) :: num

        def_default(state%init,initial,.false.)

        if (iProcIndex == root .and. .not. inited) then
            state%funit = get_free_unit()
            call init_spin_diff_output(state%funit)
            inited = .true.
        end if

        if (iProcIndex == root) then
            if (state%init .or. state%prepend) then
                write(state%funit, '("#")', advance = 'no')
                state%prepend = state%init

            else if (.not. state%prepend) then
                write(state%funit, '(" ")', advance = 'no')

            end if

            state%cols = 0
            state%cols_mc = 0
            state%mc_out = tMCOutput

            call stats_out(state, .false., iter + PreviousCycles, 'Iter.')

            do i = 1, nBasis/2
                write(num, '(i12)') i

                call stats_out(state,.false., all_inst_spin_diff(i)/ &
                    (sum(all_norm_psi_squared) / real(inum_runs,dp) * real(StepsSft,dp)), &
                    'orbital ' // trim(adjustl(num)))

            end do

            write(state%funit, *)
            call neci_flush(state%funit)

        end if

    end subroutine write_spin_diff_stats

    subroutine init_spin_diff_output(funit)
        ! routine to initialize the instant spin-diff output
        integer, intent(in) :: funit
        character(*), parameter :: this_routine = "init_spin_diff_output"
        character(30) :: filename
        character(43) :: filename2
        character(12) :: num
        logical :: exists
        integer :: i, ierr

        filename = "spin_diff_stats"

        if (tReadPops) then
            open(funit, file = filename, status = 'unknown', position = 'append')
        else
            inquire(file = filename, exist = exists)

            if (exists) then
                i = 1
                do while(exists)
                    write(num, '(i12)') i
                    filename2 = trim(adjustl(filename)) // "." // &
                        trim(adjustl(num))

                    inquire(file = filename2, exist = exists)
                    if (i > 10000) call stop_all(this_routine, &
                        "error finding free spin_diff_stats")

                    i = i + 1
                end do

                call rename(filename, filename2)
            end if

            open(funit, file = filename, status = 'unknown', iostat = ierr)

        end if

    end subroutine init_spin_diff_output

    subroutine init_spat_doub_occ_stats(funit)
        ! i need a routine to initialize the additional output, which I
        ! think should go into a seperate file for now!
        integer, intent(in) :: funit
        character(*), parameter :: this_routine = "init_spat_doub_occ_stats"
        character(30) :: filename
        character(43) :: filename2
        character(12) :: num
        logical :: exists
        integer :: i, ierr

        filename = "spatial_double_occupancy_stats"

        if (tReadPops) then
            open(funit, file = filename, status = 'unknown', position = 'append')

        else

            inquire(file=filename, exist = exists)

            ! rename the existing file an create a new one
            if (exists) then

                i = 1
                do while(exists)
                    write(num, '(i12)') i
                    filename2 = trim(adjustl(filename)) // "." // &
                                trim(adjustl(num))

                    inquire(file=filename2, exist = exists)
                    if (i > 10000) call stop_all(this_routine, &
                            "error finding free spatial_double_occupancy_stats")

                    i = i + 1
                end do

                ! i am not sure where this routine is defined:
                call rename(filename, filename2)
            end if

            open(funit, file=filename, status='unknown', iostat=ierr)

        end if


    end subroutine init_spat_doub_occ_stats

    subroutine write_spat_doub_occ_stats(initial)
        ! routine to write out the double occupancy data
        logical, intent(in), optional :: initial
        character(*), parameter :: this_routine = "write_spat_doub_occ_stats"

        type(write_state_t), save :: state
        logical, save :: inited = .false.
        integer :: i
        character(12) :: num

        def_default(state%init,initial,.false.)

        ! If the output file hasn't been opened yet, then create it.
        if (iProcIndex == Root .and. .not. inited) then
            state%funit = get_free_unit()
            call init_spat_doub_occ_stats(state%funit)
            inited = .true.
        end if

        if (iProcIndex == root) then
            if (state%init .or. state%prepend) then
                write(state%funit, '("#")', advance = 'no')
                state%prepend = state%init

            else if (.not. state%prepend) then
                write(state%funit, '(" ")', advance = 'no')

            end if

            state%cols = 0
            state%cols_mc = 0
            state%mc_out = tMCOutput

            call stats_out(state,.false., iter + PreviousCycles, 'Iter.')

            do i = 1, nBasis/2
                write(num, '(i12)') i

                call stats_out(state,.false., all_inst_spatial_doub_occ(i)/ &
                    (sum(all_norm_psi_squared) / real(inum_runs,dp) * real(StepsSft,dp)), &
                    'orbital ' // trim(adjustl(num)))

            end do

            ! And we are done
            write(state%funit, *)
            call neci_flush(state%funit)

        end if

    end subroutine write_spat_doub_occ_stats

    subroutine write_double_occ_stats(initial)
        ! routine to write out the double occupancy data
        logical, intent(in), optional :: initial
        character(*), parameter :: this_routine = "write_double_occ_stats"

        type(write_state_t), save :: state
        logical, save :: inited = .false.

        ! Provide default 'initial' option
        if (present(initial)) then
            state%init = initial
        else
            state%init = .false.
        end if

        ! If the output file hasn't been opened yet, then create it.
        if (iProcIndex == Root .and. .not. inited) then
            state%funit = get_free_unit()
            call init_double_occ_output(state%funit)
            inited = .true.
        end if

        if (iProcIndex == root) then
            if (state%init .or. state%prepend) then
                write(state%funit, '("#")', advance = 'no')
                state%prepend = state%init

            else if (.not. state%prepend) then
                write(state%funit, '(" ")', advance = 'no')

            end if

            state%cols = 0
            state%cols_mc = 0
            state%mc_out = tMCOutput

            call stats_out(state,.false., iter + PreviousCycles, 'Iter.')
            call stats_out(state,.false., all_inst_double_occ / &
                (real(StepsSft,dp) * sum(all_norm_psi_squared) / real(inum_runs, dp)), 'Double Occ.')
            if (t_calc_double_occ_av) then
               if(.not. near_zero(sum_norm_psi_squared)) then
                  call stats_out(state,.false., sum_double_occ / &
                      (real(StepsSft,dp) * sum_norm_psi_squared), 'DoubOcc Av')
               else
                  call stats_out(state,.false.,0.0_dp,'DoubOcc Av')
               endif
            else
                call stats_out(state,.false.,0.0_dp, 'DoubOcc Av')
            end if

            ! And we are done
            write(state%funit, *)
            call neci_flush(state%funit)

        end if

    end subroutine write_double_occ_stats

    subroutine init_double_occ_output(funit)
        ! i need a routine to initialize the additional output, which I
        ! think should go into a seperate file for now!
        integer, intent(in) :: funit
        character(*), parameter :: this_routine = "init_double_occ_output"
        character(30) :: filename
        character(43) :: filename2
        character(12) :: num
        logical :: exists
        integer :: i, ierr

        filename = "double_occupancy_stats"

        if (tReadPops) then
            open(funit, file = filename, status = 'unknown', position = 'append')

        else

            inquire(file=filename, exist = exists)

            ! rename the existing file an create a new one
            if (exists) then

                i = 1
                do while(exists)
                    write(num, '(i12)') i
                    filename2 = trim(adjustl(filename)) // "." // &
                                trim(adjustl(num))

                    inquire(file=filename2, exist = exists)
                    if (i > 10000) call stop_all(this_routine, &
                            "error finding free double_occupancy_stats")

                    i = i + 1
                end do

                ! i am not sure where this routine is defined:
                call rename(filename, filename2)
            end if

            open(funit, file=filename, status='unknown', iostat=ierr)

        end if

    end subroutine init_double_occ_output

    subroutine calc_double_occ_from_rdm(rdm, rdm_trace, inst_occ)
        ! also write a routine which calculates the double occupancy from the
        ! 2-rdm, if it has been calculated!
        type(rdm_list_t), intent(inout) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        real(dp), intent(out), optional :: inst_occ
        character(*), parameter :: this_routine = "calc_double_occ_from_rdm"

        integer :: ielem, ij, kl, i, j, k, l, p, q, r, s, iproc, irdm, ierr, hash_val
        integer(int_rdm) :: ijkl
        real(dp) :: rdm_sign(rdm%sign_length)
        real(dp) :: double_occ(rdm%sign_length), all_double_occ(rdm%sign_length)
        real(dp), allocatable :: spatial_double_occ(:), all_spatial_double_occ(:)
        logical :: tSuccess

        double_occ = 0.0_dp
        ! just a quick addition to calculate spatially resolved double
        ! occupancy
        if (t_spin_measurements) then
            allocate(spatial_double_occ(nBasis/2))
            allocate(all_spatial_double_occ(nBasis/2))
            spatial_double_occ = 0.0_dp
            all_spatial_double_occ = 0.0_dp
        end if
        ! todo: find out about the flags to ensure the rdm was actually
        ! calculated!
        ! seperately on all processors do this summation and then
        ! communicate the results in the end..
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
            call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

            ! normalise
            rdm_sign = rdm_sign / rdm_trace

            ! convert to spatial orbitals:
            p = spatial(i)
            q = spatial(j)
            r = spatial(k)
            s = spatial(l)

            ! only consider the diagonal elements!
            if (p == q .and. p == r .and. p == s) then
                ASSERT(.not. same_spin(i,j))
                ASSERT(.not. same_spin(k,l))

                ! add up all the diagonal contributions
                double_occ = double_occ + rdm_sign

                if (t_spin_measurements) spatial_double_occ(p) = rdm_sign(1)

            end if
        end do

        ! at the end average over the spatial orbitals
        double_occ = 2.0_dp * double_occ / real(nbasis, dp)

        ! MPI communicate:
        call MPISumAll(double_occ, all_double_occ)

        if (present(inst_occ)) then
            inst_occ = double_occ(1)
        end if

        if (t_spin_measurements) then
            call MPIAllreduce(spatial_double_occ, MPI_SUM, all_spatial_double_occ)
        end if

        if (iProcIndex == root) then
            print *, "======"
            print *, "Double occupancy from RDM: ", all_double_occ
            print *, "======"

            if (t_spin_measurements) then

                print *, "======"
                print *, "spatially resolved double occupancy from RDM: "
                print *, all_spatial_double_occ
                print *, "======"
            end if
        end if

    end subroutine calc_double_occ_from_rdm

end module double_occ_mod
