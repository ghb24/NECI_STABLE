module rdm_estimators

    use bit_rep_data, only: NIfTot
    use constants

    implicit none

contains

    subroutine init_rdm_estimates_t(est, nrdms, open_output_file)

        ! Initialise an rdm_estimates_t object. Allocate arrays to be large
        ! enough to hold estimates for nrdms RDMs.

        ! Also, if open_output_file is true, and if this is the processor with
        ! label 0, then open an RDMEstimates file, and write the file's header.

        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_estimates_t
        use util_mod, only: get_free_unit

        type(rdm_estimates_t), intent(out) :: est
        integer, intent(in) :: nrdms
        logical, intent(in) :: open_output_file

        integer :: ierr

        ! Store the number of RDMs.
        est%nrdms = nrdms

        ! Estimates over the entire RDM sampling period.
        allocate(est%trace(nrdms), stat=ierr)
        allocate(est%norm(nrdms), stat=ierr)
        allocate(est%energy_1_num(nrdms), stat=ierr)
        allocate(est%energy_2_num(nrdms), stat=ierr)
        allocate(est%energy_num(nrdms), stat=ierr)
        allocate(est%spin_num(nrdms), stat=ierr)

        ! "Instantaneous" estimates over the previous sampling block.
        allocate(est%trace_inst(nrdms), stat=ierr)
        allocate(est%norm_inst(nrdms), stat=ierr)
        allocate(est%energy_1_num_inst(nrdms), stat=ierr)
        allocate(est%energy_2_num_inst(nrdms), stat=ierr)
        allocate(est%energy_num_inst(nrdms), stat=ierr)
        allocate(est%spin_num_inst(nrdms), stat=ierr)

        ! Hermiticity errors, for the final RDMs.
        allocate(est%max_error_herm(nrdms), stat=ierr)
        allocate(est%sum_error_herm(nrdms), stat=ierr)

        est%trace = 0.0_dp
        est%norm = 0.0_dp
        est%energy_1_num = 0.0_dp
        est%energy_2_num = 0.0_dp
        est%energy_num = 0.0_dp
        est%spin_num = 0.0_dp

        est%trace_inst = 0.0_dp
        est%norm_inst = 0.0_dp
        est%energy_1_num_inst = 0.0_dp
        est%energy_2_num_inst = 0.0_dp
        est%energy_num_inst = 0.0_dp
        est%spin_num_inst = 0.0_dp

        est%max_error_herm = 0.0_dp
        est%sum_error_herm = 0.0_dp

        ! If appropriate, create a new RDMEstimates file.
        if (iProcIndex == 0 .and. open_output_file) then
            est%write_unit = get_free_unit()
            call write_rdm_est_file_header(est%write_unit, nrdms)
        end if

    end subroutine init_rdm_estimates_t

    subroutine write_rdm_est_file_header(write_unit, nrdms)

        ! Open a new RDMEstimates file (overwriting any existing file), and
        ! write a header to it, appropriate for when we are sampling nrdms RDMs.

        integer, intent(in) :: write_unit, nrdms

        integer :: irdm

        open(write_unit, file='RDMEstimates', status='unknown', position='append')

        write(write_unit, '("#", 4X, "Iteration")', advance='no')

        do irdm = 1, nrdms
            write(write_unit, '(4x,"Energy numerator",1x,i2)', advance='no') irdm
            write(write_unit, '(4x,"Spin^2 numerator",1x,i2)', advance='no') irdm
            write(write_unit, '(7x,"Normalisation",1x,i2)', advance='no') irdm
        end do

        write(write_unit,'()')

    end subroutine write_rdm_est_file_header

    subroutine calc_rdm_estimates_wrapper(est, rdm)

        ! Calculate the estimates for the RDM stored in rdm. The full estimates
        ! are stored using this object, and also instantaneous estimates. The
        ! instantaneous estimates are calculated by subtracting the previous
        ! stored values from the newly calculated ones. This works so long as
        ! the estimates are linear functions of the RDMs, which they will be
        ! for any observable.

        use Parallel_neci, only: MPISumAll
        use rdm_data, only: rdm_estimates_t, rdm_list_t
        use rdm_parallel, only: calc_rdm_trace, calc_rdm_spin, calc_rdm_energy
        use SystemData, only: nel, ecore

        type(rdm_estimates_t), intent(inout) :: est
        type(rdm_list_t), intent(in) :: rdm

        real(dp) :: rdm_trace(est%nrdms), rdm_norm(est%nrdms)
        real(dp) :: rdm_energy_1(est%nrdms), rdm_energy_2(est%nrdms)
        real(dp) :: rdm_spin(est%nrdms)

        ! Use the _inst variables as temporary variables to store the current
        ! total values. These are updated at the end of this routine.
        est%trace_inst = est%trace
        est%norm_inst = est%norm
        est%energy_1_num_inst = est%energy_1_num
        est%energy_2_num_inst = est%energy_2_num
        est%energy_num_inst = est%energy_num
        est%spin_num_inst = est%spin_num

        ! Calculate the new total values.

        call calc_rdm_trace(rdm, rdm_trace)
        call MPISumAll(rdm_trace, est%trace)

        ! RDMs are normalised so that their trace is nel*(nel-1)/2.
        rdm_norm = rdm_trace*2.0_dp/(nel*(nel-1))
        est%norm = est%trace*2.0_dp/(nel*(nel-1))

        ! The 1- and 2- electron operator contributions to the RDM energy.
        call calc_rdm_energy(rdm, rdm_energy_1, rdm_energy_2)
        call MPISumAll(rdm_energy_1, est%energy_1_num)
        call MPISumAll(rdm_energy_2, est%energy_2_num)
        ! The *total* energy, including the core contribution.
        est%energy_num = est%energy_1_num + est%energy_2_num + ecore*est%norm

        ! Estimate of the expectation value of the spin squared operator
        ! (equal to S(S+1) for spin quantum number S).
        call calc_rdm_spin(rdm, rdm_norm, rdm_spin)
        call MPISumAll(rdm_spin, est%spin_num)

        ! Calculate the instantaneous values by subtracting the old total
        ! values from the new total ones.
        est%trace_inst = est%trace - est%trace_inst
        est%norm_inst = est%norm - est%norm_inst
        est%energy_1_num_inst = est%energy_1_num - est%energy_1_num_inst
        est%energy_2_num_inst = est%energy_2_num - est%energy_2_num_inst
        est%energy_num_inst = est%energy_num - est%energy_num_inst
        est%spin_num_inst = est%spin_num - est%spin_num_inst

    end subroutine calc_rdm_estimates_wrapper

    subroutine rdm_output_wrapper(est, rdm, rdm_recv, spawn)

        ! Call routines to output RDMs in all requested forms.

        ! We also calculate the Hermitian errors here, since this is something
        ! we typically want to do at the same point (the very end of a
        ! simulation usually), and this requires large parallel
        ! communications, as does the printing.

        use LoggingData, only: tWrite_normalised_RDMs, tWriteSpinFreeRDM
        use rdm_data, only: rdm_estimates_t, rdm_list_t, rdm_spawn_t, tOpenShell
        use rdm_parallel, only: print_rdms_spin_sym_wrapper, print_spinfree_2rdm_wrapper
        use rdm_parallel, only: calc_hermitian_errors

        type(rdm_estimates_t), intent(inout) :: est
        ! IMPORTANT: rdm is not actually modified by this routine, despite
        ! needing inout status.
        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_list_t), intent(inout) :: rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn

        call calc_hermitian_errors(rdm, rdm_recv, spawn, est%norm, est%max_error_herm, est%sum_error_herm)

        if (tWriteSpinFreeRDM) call print_spinfree_2rdm_wrapper(rdm, rdm_recv, spawn, est%norm)
        if (tWrite_Normalised_RDMs) call print_rdms_spin_sym_wrapper(rdm, rdm_recv, spawn, est%norm, tOpenShell)

    end subroutine rdm_output_wrapper

    subroutine write_rdm_estimates(est, final_output)

        ! Write RDM estimates to the RDMEstimates file. Specifically, the
        ! numerator of the energy and spin^2 estimators are output, as is the
        ! trace (the denominator of these estimator).

        ! Also, if final_output is true, then output the final total estimates
        ! to standard output, and close the RDMEstimates unit.

        use FciMCData, only: Iter, PreviousCycles
        use LoggingData, only: tRDMInstEnergy
        use rdm_data, only: rdm_estimates_t
        use util_mod, only: int_fmt

        type(rdm_estimates_t), intent(in) :: est
        logical, intent(in) :: final_output

        integer :: irdm

        if (tRDMInstEnergy) then
            write(est%write_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
            do irdm = 1, est%nrdms
                write(est%write_unit, '(3(3x,es20.13))', advance='no') &
                    est%energy_num_inst(irdm), est%spin_num_inst(irdm), est%norm_inst(irdm)
            end do
            write(est%write_unit,'()')
        else
            write(est%write_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
            do irdm = 1, est%nrdms
                write(est%write_unit, '(3(3x,es20.13))', advance='no') &
                    est%energy_num(irdm), est%spin_num(irdm), est%norm(irdm)
            end do
            write(est%write_unit, '()')
        end if

        call neci_flush(est%write_unit)

        if (final_output) then
            do irdm = 1, est%nrdms
                write(6,'(1x,"FINAL ESTIMATES FOR RDM",1X,'//int_fmt(irdm)//',":",)') irdm
                write(6,'(1x,"Trace of 2-el-RDM before normalisation:",1x,es17.10)') est%trace(irdm)
                write(6,'(1x,"Trace of 2-el-RDM after normalisation:",1x,es17.10)') est%trace(irdm)/est%norm(irdm)
                write(6,'(1x,"Energy contribution from the 1-RDM:",1x,es17.10)') est%energy_1_num(irdm)/est%norm(irdm)
                write(6,'(1x,"Energy contribution from the 2-RDM:",1x,es17.10)') est%energy_2_num(irdm)/est%norm(irdm)
                write(6,'(1x,"*TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*:",1x,es20.13,/)') &
                    est%energy_num(irdm)/est%norm(irdm)
            end do

            close(est%write_unit)
        end if

    end subroutine write_rdm_estimates

end module rdm_estimators
