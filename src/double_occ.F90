#include "macros.h"

! make on module file for now which stores all the relevant stuff 
! for the double occupancy measurements
module double_occ_mod

    use SystemData, only: nel, nbasis
    use bit_rep_data, only: nifd, nOffSgn, niftot
    use constants, only: n_int, lenof_sign, write_state_t, dp, inum_runs
    use ParallelHelper, only: iProcIndex, root
    use CalcData, only: tReadPops
    use LoggingData, only: tMCOutput, t_calc_double_occ_av
    use util_mod
    use FciMCData, only: iter, PreviousCycles, norm_psi, totwalkers, all_norm_psi_squared
    use Parallel_neci, only: iProcIndex, nProcessors, MPISumAll
    use sort_mod, only: sort

    ! only after merge with master: 
!     use constants, only: int_rdm
!     use rdm_data, only: rdm_list_t
!     use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm
!     use UMatCache, only: spatial
#ifdef __DEBUG
!     use Determinants, only: writedetbit
#endif

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

contains 

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
#ifdef __DEBUG
        print *, "double occupied orbs for determinant:", n_double_orbs
        call writedetbit(6, ilut, .true.)
#endif

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
        ! do not do the division here, but only in the output! 
#ifdef __PROG_NUMRUNS
#ifdef __CMPLX 
        ! todo
#else
        ! i essentially only need two runs! 
        double_occ = real_sgn(1) * real_sgn(2) * frac_double_orbs
#endif
#else
        double_occ = abs(real_sgn(1))**2 * frac_double_orbs 
#endif

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
       
#ifdef __DEBUG
            print *, "double_occ: ", all_inst_double_occ
            print *, "norm: ", all_norm_psi_squared
            print *, "inum_runs: ", real(inum_runs,dp)
            print *, "summed double occ: ", sum_double_occ
            print *, "summed norm: ", sum_norm_psi_squared
#endif

            call stats_out(state,.false., iter + PreviousCycles, 'Iter.')
            call stats_out(state,.false., all_inst_double_occ / & 
                (sum(all_norm_psi_squared) / real(inum_runs, dp)), 'Double Occ.')
!             call stats_out(state, .true., inst_double_occ / norm_psi(1), 'Double Occ.')
            if (t_calc_double_occ_av) then
                call stats_out(state,.false., sum_double_occ / sum_norm_psi_squared, 'DoubOcc Av')
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
        
    ! only usable if i have rdms in GUGA:
!     subroutine calc_double_occ_from_rdm(rdm, rdm_trace, nrdms_to_print)
!         ! also write a routine which calculates the double occupancy from the 
!         ! 2-rdm, if it has been calculated! 
!         type(rdm_list_t), intent(inout) :: rdm
!         real(dp), intent(in) :: rdm_trace(rdm%sign_length)
!         integer, intent(in) :: nrdms_to_print
!         character(*), parameter :: this_routine = "calc_double_occ_from_rdm"
! 
!         integer :: ielem, ij, kl, i, j, k, l, p, q, r, s, iproc, irdm, ierr
!         integer(int_rdm) :: ijkl 
!         real(dp) :: rdm_sign(rdm%sign_length)
!         real(dp) :: double_occ(rdm%sign_length), all_double_occ(rdm%sign_length)
! 
!         double_occ = 0.0_dp
!         ! todo: find out about the flags to ensure the rdm was actually 
!         ! calculated! 
! 
!         call sort(rdm%elements(:,1:rdm%nelements))
! 
!         ! i have to do that over all processors i guess since the rdms are 
!         ! stored in a distributed way! 
!         ! although i could do that with an MPI communication
!         ! seperately on all processors do this summation and then 
!         ! communicate the results in the end.. 
!         do ielem = 1, rdm%nelements
!             ijkl = rdm%elements(0,ielem)
!             call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
!             call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)
! 
!             ! normalise 
!             rdm_sign = rdm_sign / rdm_trace 
! 
!             ! convert to spatial orbitals:
!             p = spatial(i)
!             q = spatial(j) 
!             r = spatial(k)
!             s = spatial(l)
! 
!             ! only consider the diagonal elements! 
!             if (p == q .and. p == r .and. p == s) then
!                 ASSERT(is_alpha(i) .and. is_beta(j) .and. is_alpha(k) .and. is_beta(l))
! 
!                 ! add up all the diagonal contributions
!                 double_occ = double_occ + rdm_sign
! 
!             end if
!         end do
! 
!         ! at the end average over the spatial orbitals 
!         double_occ = double_occ / (real(nbasis, dp) / 2.0_dp)
! 
!         ! MPI communicate: 
!         call MPISumAll(double_occ, all_double_occ)
! 
!         if (iProcIndex == root) then
!             print *, "======"
!             print *, "Double occupancy from RDM: ", all_double_occ
!             print *, "======"
!         end if
!         ! and i guess i should write it to a file too
! !         open(iunit, file = 'double_occ_from_rdm', status = 'unknown', iostat = ierr)
! !         write(
! 
! 
!     end subroutine calc_double_occ_from_rdm


end module double_occ_mod
