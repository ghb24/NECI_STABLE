module rdm_estimators

    use bit_rep_data, only: NIfTot
    use constants

    implicit none

contains

    subroutine init_rdm_estimates_t(est, nrdms)

        use rdm_data, only: rdm_estimates_t

        type(rdm_estimates_t), intent(out) :: est
        integer, intent(in) :: nrdms

        integer :: ierr

        est%nrdms = nrdms

        allocate(est%trace(nrdms), stat=ierr)
        allocate(est%norm(nrdms), stat=ierr)

        allocate(est%energy_1_num(nrdms), stat=ierr)
        allocate(est%energy_2_num(nrdms), stat=ierr)
        allocate(est%energy_tot_num(nrdms), stat=ierr)

        allocate(est%spin_num(nrdms), stat=ierr)

        allocate(est%max_error_herm(nrdms), stat=ierr)
        allocate(est%sum_error_herm(nrdms), stat=ierr)

        est%trace = 0.0_dp
        est%norm = 0.0_dp

        est%energy_1_num = 0.0_dp
        est%energy_2_num = 0.0_dp
        est%energy_tot_num = 0.0_dp

        est%spin_num = 0.0_dp

        est%max_error_herm = 0.0_dp
        est%sum_error_herm = 0.0_dp

    end subroutine init_rdm_estimates_t

    subroutine rdm_output_wrapper(est, rdm, rdm_recv, spawn)

        use FciMCData, only: tFinalRDMEnergy
        use LoggingData, only: tWrite_normalised_RDMs, tWriteSpinFreeRDM
        use Parallel_neci, only: MPISumAll
        use rdm_data, only: rdm_estimates_t, rdm_list_t, rdm_spawn_t, tOpenShell
        use rdm_parallel, only: calc_rdm_trace, calc_rdm_spin, calc_rdm_energy, print_rdms_with_spin
        use rdm_parallel, only: print_rdms_spin_sym_wrapper, print_spinfree_2rdm_wrapper
        use rdm_parallel, only: calc_hermitian_errors
        use SystemData, only: nel, ecore

        type(rdm_estimates_t), intent(inout) :: est
        ! IMPORTANT: rdm is not actually modified by this routine, despite
        ! needing inout status.
        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_list_t), intent(inout) :: rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn

        real(dp) :: rdm_trace(est%nrdms), rdm_norm(est%nrdms)
        real(dp) :: rdm_energy_1(est%nrdms), rdm_energy_2(est%nrdms)
        real(dp) :: rdm_spin(est%nrdms)

        call calc_rdm_trace(rdm, rdm_trace)
        call MPISumAll(rdm_trace, est%trace)

        rdm_norm = rdm_trace*2.0_dp/(nel*(nel-1))
        est%norm = est%trace*2.0_dp/(nel*(nel-1))

        call calc_rdm_energy(rdm, rdm_energy_1, rdm_energy_2)
        call MPISumAll(rdm_energy_1, est%energy_1_num)
        call MPISumAll(rdm_energy_2, est%energy_2_num)
        est%energy_tot_num = est%energy_1_num + est%energy_2_num + ecore*est%norm

        call calc_rdm_spin(rdm, rdm_norm, rdm_spin)
        call MPISumAll(rdm_spin, est%spin_num)

        call calc_hermitian_errors(rdm, rdm_recv, spawn, est%norm, est%max_error_herm, est%sum_error_herm)

        if (tWriteSpinFreeRDM .and. tFinalRDMEnergy) call print_spinfree_2rdm_wrapper(rdm, rdm_recv, spawn, est%norm)
        if (tFinalRDMEnergy .and. tWrite_Normalised_RDMs) then
            call print_rdms_spin_sym_wrapper(rdm, rdm_recv, spawn, est%norm, tOpenShell)
        end if

    end subroutine rdm_output_wrapper

    subroutine write_rdm_estimates(est, est_old)

        use FciMCData, only: tFinalRDMEnergy, Iter, PreviousCycles
        use LoggingData, only: tRDMInstEnergy
        use rdm_data, only: rdm_estimates_t, rdm_estimates_unit
        use rdm_data, only: rdm_estimates_old_t
        use util_mod, only: int_fmt

        type(rdm_estimates_t), intent(in) :: est
        type(rdm_estimates_old_t), intent(in) :: est_old(:)

        integer :: i

        if (tRDMInstEnergy) then
            write(rdm_estimates_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
            do i = 1, est%nrdms
                write(rdm_estimates_unit, '(6(3x,es20.13))', advance='no') &
                    est_old(i)%RDMEnergy_Inst, est_old(i)%spin_est, 1.0_dp/est_old(i)%Norm_2RDM_Inst, &
                    est%energy_tot_num(i), est%spin_num(i), est%norm(i)
            end do
            write(rdm_estimates_unit,'()')

        else

            write(rdm_estimates_unit, '(1x,i13)') Iter+PreviousCycles
            do i = 1, est%nrdms
                write(rdm_estimates_unit, '(3(3x,es20.13))', advance='no') &
                    est_old(i)%RDMEnergy, est_old(i)%spin_est, 1.0_dp/est_old(i)%Norm_2RDM_Inst
            end do
            write(rdm_estimates_unit, '()')
        end if

        call neci_flush(rdm_estimates_unit)

        if (tFinalRDMEnergy) then
            do i = 1, est%nrdms
                write(6,'(1x,"FINAL ESTIMATES FOR RDM",1X,'//int_fmt(i)//',":",)') i
                write(6,'(1x,"Trace of 2-el-RDM before normalisation:",1x,es17.10)') est%trace(i)
                write(6,'(1x,"Trace of 2-el-RDM after normalisation:",1x,es17.10)') est%trace(i)/est%norm(i)
                write(6,'(1x,"Energy contribution from the 1-RDM:",1x,es17.10)') est%energy_1_num(i)/est%norm(i)
                write(6,'(1x,"Energy contribution from the 2-RDM:",1x,es17.10)') est%energy_2_num(i)/est%norm(i)
                write(6,'(1x,"*TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*:",1x,es20.13,/)') &
                    est%energy_tot_num(i)/est%norm(i)
            end do
            close(rdm_estimates_unit)
        end if

    end subroutine write_rdm_estimates

end module rdm_estimators
