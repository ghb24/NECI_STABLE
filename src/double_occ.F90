#include "macros.h"

! make on module file for now which stores all the relevant stuff 
! for the double occupancy measurements
module double_occ_mod

    use SystemData, only: nel, nbasis
    use bit_rep_data, only: nifd, nOffSgn, niftot
    use constants, only: n_int, lenof_sign, write_state_t, dp
    use ParallelHelper, only: iProcIndex, root
    use CalcData, only: tReadPops
    use LoggingData, only: tMCOutput
    use util_mod
    use FciMCData, only: iter, PreviousCycles, norm_psi, totwalkers
    implicit none

    ! data storage part: 
    ! i need local and global storage containers for <n_u n_d> 
    ! and maybe also instantaneous, averaged, summed over quantities
    real(dp) :: inst_double_occ = 0.0_dp

contains 

    function count_double_orbs(ilut) result(n_double_orbs) 
        ! returns the number of doubly occupied orbitals for a given 
        ! determinant. 
        use DetBitOps, only: count_open_orbs
        integer(n_int), intent(in) :: ilut(0:nifd) 
        integer :: n_double_orbs 
        character(*), parameter :: this_routine = "count_double_occupancy"

        integer :: n_open_orbs

        n_open_orbs = count_open_orbs(ilut)

        ! the doubly occupied orbitals are just the number of electrons 
        ! minus the number of open orbitals divided by 2
        n_double_orbs = (nel - n_open_orbs)/2
    end function count_double_orbs

    
    function get_double_occupancy(ilut, real_sgn) result(double_occ)
        ! function to get the contribution to the double occupancy for a 
        ! given determinant, by calculating the 
        ! |C_I|^2 *\sum_i <n_iu n_id>/n_spatial_orbs
        integer(n_int), intent(in) :: ilut(NIfTot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        real(dp) :: double_occ
        character(*), parameter :: this_routine = "get_double_occupancy"

        integer :: n_double_orbs
        real(dp) :: frac_double_orbs
        integer(n_int) :: sgn(lenof_sign)

#ifdef __CMPLX
        complex(dp) :: complex_sgn
#endif

        ! here i need to be careful, if it is a double or single run and 
        ! stuff
        n_double_orbs = count_double_orbs(ilut(0:nifd))

        ! i only need the fraction of doubly occupied orbitals out of all 
        ! spatial orbitals
        frac_double_orbs = real(n_double_orbs,dp) / (real(nbasis,dp) / 2.0_dp)

        ! now i want to sum in the walkers from the runs..
        ! todo: have to figure out how to access the different runs

        ! extract the walker occupation
!         sgn = ilut(nOffSgn:nOffSgn + lenof_sign -1)
!         real_sgn = transfer(sgn, real_sgn)

        ! do i want to do that for complex walkers also?? i guess so..
        ! to get it running do it only for  single run for now! 
        double_occ = abs(real_sgn(1))**2 * frac_double_orbs / norm_psi(1)

    end function get_double_occupancy

    subroutine rezero_double_occ_stats
        ! at the end of each cycle i should rezero the non-summed and 
        ! non-averaged quantities to 0
        character(*), parameter :: this_routine = "rezero_double_occ_stats"

        ! there is more stuff to come..
        inst_double_occ = 0.0_dp

    end subroutine rezero_double_occ_stats

    subroutine write_double_occ_stats(iter_data, initial)
        ! routine to write out the double occupancy data
        type(fcimc_iter_data), intent(in) :: iter_data
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
        
            call stats_out(state,.true., iter + PreviousCycles, 'Iter.')
            call stats_out(state, .true., inst_double_occ/totwalkers, 'Double Occ.')

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
        
end module double_occ_mod
