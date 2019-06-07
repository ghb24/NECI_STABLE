#include "macros.h"

module rdm_estimators

    use CalcData, only: tAdaptiveShift
    use bit_rep_data, only: NIfTot
    use constants
    use rdm_data, only: rdm_list_t, rdm_spawn_t
    use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm, calc_rdm_trace
    use SystemData, only: tGUGA

    implicit none

contains

    subroutine init_rdm_estimates_t(est, nrdms_standard, nrdms_transition, open_output_file, &
         filename)

        ! Initialise an rdm_estimates_t object. Allocate arrays to be large
        ! enough to hold estimates for nrdms_srandard+nrdms_transition RDMs.

        ! Also, if open_output_file is true, and if this is the processor with
        ! label 0, then open an RDMEstimates file, and write the file's header.

        use CalcData, only: tEN2
        use LoggingData, only: tCalcPropEst, iNumPropToEst
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_estimates_t
        use util_mod, only: get_free_unit

        type(rdm_estimates_t), intent(out) :: est
        integer, intent(in) :: nrdms_standard, nrdms_transition
        logical, intent(in) :: open_output_file
        character(*), intent(in), optional :: filename

        integer :: nrdms, ierr
        character(255) :: rdm_filename

        nrdms = nrdms_standard + nrdms_transition

        ! Store the number of RDMs.
        est%nrdms = nrdms
        est%nrdms_standard = nrdms_standard
        est%nrdms_transition = nrdms_transition

        ! Estimates over the entire RDM sampling period.
        allocate(est%trace(nrdms), stat=ierr)
        allocate(est%norm(nrdms), stat=ierr)
        allocate(est%energy_1_num(nrdms), stat=ierr)
        allocate(est%energy_2_num(nrdms), stat=ierr)
        allocate(est%energy_num(nrdms), stat=ierr)
        allocate(est%spin_num(nrdms), stat=ierr)
        if (tCalcPropEst) allocate(est%property(iNumPropToEst,nrdms), stat=ierr)
        if (tEN2) allocate(est%energy_pert(nrdms_standard), stat=ierr)

        ! "Instantaneous" estimates over the previous sampling block.
        allocate(est%trace_inst(nrdms), stat=ierr)
        allocate(est%norm_inst(nrdms), stat=ierr)
        allocate(est%energy_1_num_inst(nrdms), stat=ierr)
        allocate(est%energy_2_num_inst(nrdms), stat=ierr)
        allocate(est%energy_num_inst(nrdms), stat=ierr)
        allocate(est%spin_num_inst(nrdms), stat=ierr)
        if (tCalcPropEst) allocate(est%property_inst(iNumPropToEst,nrdms), stat=ierr)
        if (tEN2) allocate(est%energy_pert_inst(nrdms_standard), stat=ierr)

        ! Hermiticity errors, for the final RDMs.
        allocate(est%max_error_herm(nrdms), stat=ierr)
        allocate(est%sum_error_herm(nrdms), stat=ierr)

        est%trace = 0.0_dp
        est%norm = 0.0_dp
        est%energy_1_num = 0.0_dp
        est%energy_2_num = 0.0_dp
        est%energy_num = 0.0_dp
        est%spin_num = 0.0_dp
        if (tCalcPropEst) est%property = 0.0_dp
        if (tEN2) est%energy_pert = 0.0_dp

        est%trace_inst = 0.0_dp
        est%norm_inst = 0.0_dp
        est%energy_1_num_inst = 0.0_dp
        est%energy_2_num_inst = 0.0_dp
        est%energy_num_inst = 0.0_dp
        est%spin_num_inst = 0.0_dp
        if (tCalcPropEst) est%property_inst = 0.0_dp
        if (tEN2) est%energy_pert_inst = 0.0_dp

        est%max_error_herm = 0.0_dp
        est%sum_error_herm = 0.0_dp

        ! If appropriate, create a new RDMEstimates file.
        if (iProcIndex == 0 .and. open_output_file) then
            est%write_unit = get_free_unit()
            if(present(filename)) then
               rdm_filename = filename
            else
               rdm_filename = "RDMEstimates"
            endif
            call write_rdm_est_file_header(est%write_unit, nrdms_standard, nrdms_transition, &
                 rdm_filename)
        else
            ! If we don't have a file open with this unit, set it to something
            ! unique, so we can easily check, and which will cause an obvious
            ! error if we don't.
            est%write_unit = huge(est%write_unit)
        end if

    end subroutine init_rdm_estimates_t

    subroutine dealloc_rdm_estimates_t(est)

        ! Initialise an rdm_estimates_t object. Allocate arrays to be large
        ! enough to hold estimates for nrdms RDMs.

        ! Also, if open_output_file is true, and if this is the processor with
        ! label 0, then open an RDMEstimates file, and write the file's header.

        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_estimates_t
        use util_mod, only: get_free_unit

        type(rdm_estimates_t), intent(inout) :: est

        integer :: ierr

        if (allocated(est%trace)) deallocate(est%trace, stat=ierr)
        if (allocated(est%norm)) deallocate(est%norm, stat=ierr)
        if (allocated(est%energy_1_num)) deallocate(est%energy_1_num, stat=ierr)
        if (allocated(est%energy_2_num)) deallocate(est%energy_2_num, stat=ierr)
        if (allocated(est%energy_num)) deallocate(est%energy_num, stat=ierr)
        if (allocated(est%spin_num)) deallocate(est%spin_num, stat=ierr)
        if (allocated(est%property)) deallocate(est%property, stat=ierr)
        if (allocated(est%trace_inst)) deallocate(est%trace_inst, stat=ierr)
        if (allocated(est%norm_inst)) deallocate(est%norm_inst, stat=ierr)
        if (allocated(est%energy_1_num_inst)) deallocate(est%energy_1_num_inst, stat=ierr)
        if (allocated(est%energy_2_num_inst)) deallocate(est%energy_2_num_inst, stat=ierr)
        if (allocated(est%energy_num_inst)) deallocate(est%energy_num_inst, stat=ierr)
        if (allocated(est%spin_num_inst)) deallocate(est%spin_num_inst, stat=ierr)
        if (allocated(est%property_inst)) deallocate(est%property_inst, stat=ierr)
        if (allocated(est%max_error_herm)) deallocate(est%max_error_herm, stat=ierr)
        if (allocated(est%sum_error_herm)) deallocate(est%sum_error_herm, stat=ierr)
        if (allocated(est%energy_pert)) deallocate(est%energy_pert, stat=ierr)
        if (allocated(est%energy_pert_inst)) deallocate(est%energy_pert_inst, stat=ierr)

        ! Close the RDMEstimates unit, if it was opened on this processor.
        ! The following was what it was set to if it was not opened in
        ! init_rdm_estimates_t, so don't attempt a close in that case.
        if (est%write_unit /= huge(est%write_unit)) close(est%write_unit)

    end subroutine dealloc_rdm_estimates_t

    subroutine write_rdm_est_file_header(write_unit, nrdms_standard, nrdms_transition, &
         filename)

        ! Open a new RDMEstimates file (overwriting any existing file), and
        ! write a header to it, appropriate for when we are sampling nrdms RDMs.

        use CalcData, only: tEN2
        use LoggingData, only: tCalcPropEst, iNumPropToEst

        integer, intent(in) :: write_unit, nrdms_standard, nrdms_transition
        character(255), intent(in) :: filename

        integer :: irdm, iprop

        open(write_unit, file=trim(filename), status='unknown', position='append')

        write(write_unit, '("#", 4X, "Iteration")', advance='no')

        do irdm = 1, nrdms_standard
            write(write_unit, '(4x,"Energy numerator",1x,i2)', advance='no') irdm
            if (.not. tGUGA) then
                write(write_unit, '(4x,"Spin^2 numerator",1x,i2)', advance='no') irdm
            end if
            if (tEN2) then
                write(write_unit, '(7x,"EN2 numerator",1x,i2)', advance='no') irdm
                write(write_unit, '(3x,"Var+EN2 numerator",1x,i2)', advance='no') irdm
            end if
            if (tCalcPropEst) then
                do iprop = 1,iNumPropToEst
                    write(write_unit, '(4x,"Property(",i2,")",1x,i2)', advance='no') iprop, irdm
                end do
            end if
            write(write_unit, '(7x,"Normalisation",1x,i2)', advance='no') irdm
        end do

        do irdm = nrdms_standard+1, nrdms_standard+nrdms_transition
            if (tCalcPropEst) then
                do iprop = 1,iNumPropToEst
                    write(write_unit, '(4x,"Property(",i2,")",1x,i2)', advance='no') iprop, irdm
                end do
            end if
        end do

        write(write_unit,'()')

    end subroutine write_rdm_est_file_header

    subroutine calc_2rdm_estimates_wrapper(rdm_defs, est, rdm, en_pert)

        ! Calculate the estimates for the 2-RDM stored in rdm. The full estimates
        ! are stored using this object, and also instantaneous estimates. The
        ! instantaneous estimates are calculated by subtracting the previous
        ! stored values from the newly calculated ones. This works so long as
        ! the estimates are linear functions of the RDMs, which they will be
        ! for any observable.

        use CalcData, only: tEN2
        use Parallel_neci, only: MPISumAll
        use rdm_data, only: rdm_estimates_t, rdm_list_t, rdm_definitions_t
        use rdm_data, only: en_pert_t
        use SystemData, only: nel, ecore
        use OneEInts, only: PropCore
        use LoggingData, only: iNumPropToEst, tCalcPropEst

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_estimates_t), intent(inout) :: est
        type(rdm_list_t), intent(in) :: rdm
        type(en_pert_t), intent(in) :: en_pert

        integer :: irdm, iprop
        real(dp) :: rdm_trace(est%nrdms), rdm_norm(est%nrdms)
        real(dp) :: rdm_energy_1(est%nrdms), rdm_energy_2(est%nrdms)
        real(dp) :: rdm_spin(est%nrdms)
        real(dp) :: rdm_prop(iNumProptoEst,est%nrdms)
        real(dp) :: energy_pert(est%nrdms_standard)

        ! Use the _inst variables as temporary variables to store the current
        ! total values. These are updated at the end of this routine.
        est%trace_inst = est%trace
        est%norm_inst = est%norm
        est%energy_1_num_inst = est%energy_1_num
        est%energy_2_num_inst = est%energy_2_num
        est%energy_num_inst = est%energy_num
        est%spin_num_inst = est%spin_num
        if (tCalcPropEst) est%property_inst = est%property

        ! Calculate the new total values.

        call calc_rdm_trace(rdm, rdm_trace)
        call MPISumAll(rdm_trace, est%trace)

        ! RDMs are normalised so that their trace is nel*(nel-1)/2.
        rdm_norm = rdm_trace*2.0_dp/(nel*(nel-1))
        est%norm = est%trace*2.0_dp/(nel*(nel-1))

        ! For the transition RDMs, we want to calculate the norms using the
        ! non-transition RDMs.
        do irdm = est%nrdms_standard+1, est%nrdms
            est%norm(irdm) = sqrt( est%norm(rdm_defs%state_labels(1,irdm)) * est%norm(rdm_defs%state_labels(2,irdm)) )
        end do

        ! The 1- and 2- electron operator contributions to the RDM energy.
        call calc_rdm_energy(rdm, rdm_energy_1, rdm_energy_2)
        call MPISumAll(rdm_energy_1, est%energy_1_num)
        call MPISumAll(rdm_energy_2, est%energy_2_num)
        ! The *total* energy, including the core contribution.
        est%energy_num = est%energy_1_num + est%energy_2_num + ecore*est%norm

        ! Estimate of the expectation value of the spin squared operator
        ! (equal to S(S+1) for spin quantum number S).
        if (.not. tGUGA) then
            call calc_rdm_spin(rdm, rdm_norm, rdm_spin)
            call MPISumAll(rdm_spin, est%spin_num)
        end if

        if (tCalcPropEst) then 
            ! Estimate of the properties using different property integrals
            ! and all the standard and transition rdms. 
            call calc_rdm_prop(rdm, rdm_prop)
            call MPISumAll(rdm_prop, est%property)
            ! Add the contribution from the core (zero body) part of the perturbation. 
            do iprop = 1, iNumPropToEst
                est%property(iprop,:) = est%property(iprop,:) + est%norm*PropCore(iprop)
            end do
        end if

        ! Calculate the instantaneous values by subtracting the old total
        ! values from the new total ones.
        est%trace_inst = est%trace - est%trace_inst
        est%norm_inst = est%norm - est%norm_inst
        est%energy_1_num_inst = est%energy_1_num - est%energy_1_num_inst
        est%energy_2_num_inst = est%energy_2_num - est%energy_2_num_inst
        est%energy_num_inst = est%energy_num - est%energy_num_inst
        est%spin_num_inst = est%spin_num - est%spin_num_inst
        if(tCalcPropEst) est%property_inst = est%property - est%property_inst

        ! For the EN Perturbation terms, we clear them at the start of
        ! every RDM averaging cycle, so they're treated a bit differently.
        if (tEN2) then
            call calc_en_pert_energy(en_pert, est%energy_num_inst(1:en_pert%sign_length), &
                                     est%norm_inst(1:en_pert%sign_length), energy_pert)

            call MPISumAll(energy_pert, est%energy_pert_inst)

            est%energy_pert = est%energy_pert + est%energy_pert_inst
        end if

    end subroutine calc_2rdm_estimates_wrapper

    subroutine write_rdm_estimates(rdm_defs, est, final_output, write_to_separate_file, &
         tInitsRDM)

        ! Write RDM estimates to the RDMEstimates file. Specifically, the
        ! numerator of the energy and spin^2 estimators are output, as is the
        ! trace (the denominator of these estimator).

        ! Also, if final_output is true, then output the final total estimates
        ! to standard output, and close the RDMEstimates unit.

        use CalcData, only: tEN2
        use FciMCData, only: Iter, PreviousCycles
        use LoggingData, only: tRDMInstEnergy, tCalcPropEst, iNumPropToEst
        use rdm_data, only: rdm_estimates_t, rdm_definitions_t
        use util_mod, only: int_fmt

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_estimates_t), intent(in) :: est
        logical, intent(in) :: final_output, write_to_separate_file, tInitsRDM

        integer :: irdm, iprop

        write(6,*) "Writing RDMs to file at iteration ", iter

        if (write_to_separate_file) then
            if (tRDMInstEnergy) then
                write(est%write_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
                do irdm = 1, est%nrdms_standard
                    write(est%write_unit, '(3x,es20.13)', advance='no') &
                        est%energy_num_inst(irdm)
                    
                    if (.not. tGUGA) then
                        write(est%write_unit, '(3x,es20.13)', advance='no') &
                            est%spin_num_inst(irdm)
                    end if
                    if (tEN2) then
                        write(est%write_unit,'(2(3x,es20.13))', advance='no') &
                            est%energy_pert_inst(irdm), est%energy_pert_inst(irdm) + est%energy_num_inst(irdm)
                    end if
                    if (tCalcPropEst) then
                        do iprop=1,iNumPropToEst
                            write(est%write_unit,'(3x,es20.13)', advance='no') &
                                est%property_inst(iprop,irdm)
                        end do 
                    end if
                    write(est%write_unit, '(3x,es20.13)', advance='no') &
                        est%norm_inst(irdm)
                end do
                do irdm = est%nrdms_standard+1, est%nrdms_standard+est%nrdms_transition
                    if(tCalcPropEst) then
                        do iprop = 1, iNumPropToEst
                            write(est%write_unit,'(3x,es20.13)', advance='no') &
                                est%property_inst(iprop,irdm)
                        enddo 
                    endif
                end do
                write(est%write_unit,'()')
            else
                write(est%write_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
                do irdm = 1, est%nrdms_standard
                    write(est%write_unit, '(3x,es20.13)', advance='no') &
                        est%energy_num(irdm)
                    if (.not. tGUGA) then
                        write(est%write_unit, '(3x,es20.13)', advance='no') &
                            est%spin_num(irdm)
                    end if
                    if (tEN2) then
                        write(est%write_unit,'(2(3x,es20.13))', advance='no') &
                            est%energy_pert(irdm), est%energy_pert(irdm) + est%energy_num(irdm)
                    end if
                    if (tCalcPropEst) then
                        do iprop = 1, iNumPropToEst
                            write(est%write_unit,'(3x,es20.13)', advance='no') &
                                est%property(iprop,irdm)
                        end do 
                    end if
                    write(est%write_unit, '(3x,es20.13)', advance='no') &
                        est%norm(irdm)
                end do
                do irdm = est%nrdms_standard+1, est%nrdms_standard+est%nrdms_transition
                    if(tCalcPropEst) then
                        do iprop=1,iNumPropToEst
                            write(est%write_unit,'(3x,es20.13)',advance='no') &
                                est%property(iprop,irdm)
                        enddo 
                    endif
                end do
                write(est%write_unit, '()')
            end if

            call neci_flush(est%write_unit)
        end if

        if (final_output) then
            ! Banner for the start of the 2-RDM section in output.
           if(tInitsRDM) then
              write(6,'(1x,2("="),1x,"INFORMATION FOR FINAL 2-RDMS (Initiators)",1x,57("="),/)')
           else if(tAdaptiveShift) then
              write(6,'(1x,2("="),1x,"INFORMATION FOR FINAL 2-RDMS (Lagrangian)",1x,57("="),/)')
           else
              write(6,'(1x,2("="),1x,"INFORMATION FOR FINAL 2-RDMS",1x,57("="),/)')
           endif

            do irdm = 1, est%nrdms_standard
                write(6,'(1x,"2-RDM ESTIMATES FOR STATE",1x,'//int_fmt(irdm)//',":",)') irdm

                write(6,'(1x,"Trace of 2-el-RDM before normalisation:",1x,es17.10)') est%trace(irdm)
                write(6,'(1x,"Trace of 2-el-RDM after normalisation:",1x,es17.10)') est%trace(irdm)/est%norm(irdm)

                write(6,'(1x,"Energy contribution from the 1-RDM:",1x,es17.10)') est%energy_1_num(irdm)/est%norm(irdm)
                write(6,'(1x,"Energy contribution from the 2-RDM:",1x,es17.10)') est%energy_2_num(irdm)/est%norm(irdm)
                write(6,'(1x,"*TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*:",1x,es20.13,/)') &
                    est%energy_num(irdm)/est%norm(irdm)

                if (tEN2) then
                    write(6,'(1x,"EN2 corrections are below. Note that these may have a much &
                                  &larger error bar than the",/," variational energy above. Please do a &
                                  &blocking analysis rather than just using the energies below.")')
                    write(6,'(1x,"EN2 energy correction:",1x,es17.10)') est%energy_pert(irdm)/est%norm(irdm)
                    write(6,'(1x,"*TOTAL ENERGY* including the EN2 correction:",1x,es17.10,/)') &
                        (est%energy_num(irdm) + est%energy_pert(irdm))/est%norm(irdm)
                end if

                ! Hermiticity error measures.
                write(6,'(1x,"Hermiticty error estimates:")')
                write(6,'(1x,i15,f30.20,5x,a41)') Iter+PreviousCycles, est%max_error_herm(irdm), &
                                                '(Iteration, MAX ABS ERROR IN HERMITICITY)'
                write(6,'(1x,i15,f30.20,5x,a41,/)') Iter+PreviousCycles, est%sum_error_herm(irdm), &
                                                '(Iteration, SUM ABS ERROR IN HERMITICITY)'
            end do
            do irdm = est%nrdms_standard+1, est%nrdms
                associate(state_labels => rdm_defs%state_labels, repeat_label => rdm_defs%repeat_label)
                    write(6,'(1x,"2-RDM ESTIMATES FOR TRANSITION",1x,'//int_fmt(state_labels(2,irdm))//'," -> ",&
                              &'//int_fmt(state_labels(1,irdm))//',1x,"(",i1,")",)') &
                              state_labels(2,irdm), state_labels(1,irdm), repeat_label(irdm)
                end associate

                write(6,'(1x,"Trace of 2-el-RDM before normalisation:",1x,es17.10)') est%trace(irdm)
                write(6,'(1x,"Trace of 2-el-RDM after normalisation:",1x,es17.10,/)') est%trace(irdm)/est%norm(irdm)

                ! Hermiticity difference measures - these shouldn't be zero for
                ! transition RDMs, but it useful to give them to the test
                ! suite, to make sure somebody doesn't change something to
                ! start adding an RDM element on the wrong side of the diagonal.
                write(6,'(1x,"Hermiticty difference estimates, for test suite:")')
                write(6,'(1x,i15,f30.20,5x,a41)') Iter+PreviousCycles, est%max_error_herm(irdm), &
                                                '(Iteration, MAX ABS DIFF IN HERMITICITY)'
                write(6,'(1x,i15,f30.20,5x,a41,/)') Iter+PreviousCycles, est%sum_error_herm(irdm), &
                                                '(Iteration, SUM ABS DIFF IN HERMITICITY)'
            end do

            ! Banner for the end of the 2-RDM section in output.
            write(6,'(1x,89("="))')
        end if

    end subroutine write_rdm_estimates

    subroutine calc_rdm_energy(rdm, rdm_energy_1, rdm_energy_2)

        ! Calculate both the 1- and 2-electron contributions of the
        ! (unnormalised) energy from the 2-RDM object in rdm, and output them
        ! to rdm_energy_1 and rdm_energy_2.

        use rdm_data, only: rdm_list_t
        use rdm_integral_fns, only: one_elec_int, two_elec_int
        use SystemData, only: nel

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(out) :: rdm_energy_1(rdm%sign_length)
        real(dp), intent(out) :: rdm_energy_2(rdm%sign_length)

        integer(int_rdm) :: ijkl
        integer :: ielem, ij, kl, i, j, k, l
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_energy_1 = 0.0_dp
        rdm_energy_2 = 0.0_dp

        ! Loop over all elements in the 2-RDM.
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

            ! Decode pqrs label into p, q, r and s labels.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

            ! The 2-RDM contribution to the energy:
            rdm_energy_2 = rdm_energy_2 + rdm_sign*two_elec_int(i,j,k,l)
            ! The 1-RDM contribution to the energy:
            if (i == k) rdm_energy_1 = rdm_energy_1 + rdm_sign*one_elec_int(j,l)/(nel-1)
            if (j == l) rdm_energy_1 = rdm_energy_1 + rdm_sign*one_elec_int(i,k)/(nel-1)
            if (i == l) rdm_energy_1 = rdm_energy_1 - rdm_sign*one_elec_int(j,k)/(nel-1)
            if (j == k) rdm_energy_1 = rdm_energy_1 - rdm_sign*one_elec_int(i,l)/(nel-1)
        end do

    end subroutine calc_rdm_energy

    subroutine calc_rdm_prop(rdm, rdm_prop)

        use rdm_data, only: rdm_list_t
        use rdm_integral_fns, only: GetPropInts
        use SystemData, only: nel
        use LoggingData, only: iNumPropToEst

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(out) :: rdm_prop(iNumPropToEst,rdm%sign_length)

        integer(int_rdm) :: ijkl
        integer :: ielem, iprop, ij, kl, i, j, k, l
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_prop = 0.0_dp

        ! Loop over all elements in the 2-RDM.
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

            ! Decode pqrs label into p, q, r and s labels.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

            do iprop=1, iNumPropToEst
            ! We get dot product of the 1-RDM and one-electron property integrals:
                if (i == k) rdm_prop(iprop,:) = rdm_prop(iprop,:) + rdm_sign*GetPropInts(j,l,iprop)/(nel-1)
                if (j == l) rdm_prop(iprop,:) = rdm_prop(iprop,:) + rdm_sign*GetPropInts(i,k,iprop)/(nel-1)
                if (i == l) rdm_prop(iprop,:) = rdm_prop(iprop,:) - rdm_sign*GetPropInts(j,k,iprop)/(nel-1)
                if (j == k) rdm_prop(iprop,:) = rdm_prop(iprop,:) - rdm_sign*GetPropInts(i,l,iprop)/(nel-1)
            end do
        end do
    end subroutine calc_rdm_prop

    subroutine calc_rdm_spin(rdm, rdm_norm, rdm_spin)

        ! Return the (unnormalised) estimate of <S^2> from the instantaneous
        ! 2-RDM estimates.

        use rdm_data, only: rdm_spawn_t
        use SystemData, only: nel
        use UMatCache, only: spatial

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_norm(rdm%sign_length)
        real(dp), intent(out) :: rdm_spin(rdm%sign_length)

        integer(int_rdm) :: ijkl
        integer :: ielem
        integer :: ij, kl, i, j, k, l ! spin orbitals
        integer :: p, q, r, s ! spatial orbitals
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_spin = 0.0_dp

        ! Loop over all RDM elements.
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            ! Obtain spin orbital labels.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
            ! Obtain spatial orbital labels.
            p = spatial(i); q = spatial(j);
            r = spatial(k); s = spatial(l);

            ! Note to the reader for the following code: if mod(i,2) == 1 then
            ! i is a beta (b) orbital, if mod(i,2) == 0 then it is an alpha (a)
            ! obrital.

            ! The following if-statement allows IJIJ spatial combinations.
            if (p == r .and. q == s) then
                ! If we get to this point then we definitely have a contribution
                ! to add in, so extract the sign.
                call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

                ! If all labels have the same spatial part (IIII):
                if (p == q) then
                    if (is_beta(i) .and. is_alpha(j) .and. is_beta(k) .and. is_alpha(l)) then
                        rdm_spin = rdm_spin - 6.0_dp*rdm_sign
                    end if

                else
                    ! We only get here if the spatial parts obey IJIJ, for I /= J:

                    ! The following if-statement allows the following spin combinations:
                    ! aaaa, bbbb, abab and baba.
                    if (mod(i,2) == mod(k,2) .and. mod(j,2) == mod(l,2)) then

                        if (mod(i,2) == mod(j,2)) then
                            ! aaaa and bbbb.
                            rdm_spin = rdm_spin + 2.0_dp*rdm_sign
                        else
                            ! abab and baba.
                            rdm_spin = rdm_spin - 2.0_dp*rdm_sign
                        end if
                    else
                        ! We only get here if the spin parts are abba or baab.
                        rdm_spin = rdm_spin + 4.0_dp*rdm_sign
                    end if

                end if
            end if
        end do

        rdm_spin = rdm_spin + 3.0_dp*real(nel,dp)*rdm_norm
        rdm_spin = rdm_spin/4.0_dp

    end subroutine calc_rdm_spin

    subroutine calc_en_pert_energy(en_pert, rdm_energy_num, rdm_norm, &
                                   energy_pert)

        ! Calculate the Epstein-Nesbet perturbation energy.

        ! This is defined as
        !
        ! E_{PN} = \sum_{a} \frac{ ( \sum_i H_{ai} \psi_i )^2 }{ E_{RDM} - E_{aa} }.
        !
        ! The values ( \sum_i H_{ai} \psi_i )^2 are what is stored as the signs
        ! in the en_pert%dets objects.

        use bit_rep_data, only: nIfDBO, NIfTot
        use bit_reps, only: decode_bit_det
        use determinants, only: get_helement
        use FciMCData, only: Hii
        use hphf_integrals, only: hphf_diag_helement
        use rdm_data, only: en_pert_t
        use rdm_data_utils, only: extract_sign_EN
        use SystemData, only: nel, tHPHF

        type(en_pert_t), intent(in) :: en_pert
        real(dp), intent(in) :: rdm_energy_num(en_pert%sign_length)
        real(dp), intent(in) :: rdm_norm(en_pert%sign_length)
        real(dp), intent(out) :: energy_pert(en_pert%sign_length)

        integer(n_int) :: ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: idet, istate
        real(dp) :: h_aa
        real(dp) :: contrib(en_pert%sign_length)
        real(dp) :: contrib_rdm(en_pert%sign_length)

        energy_pert = 0.0_dp
        ilut = 0_n_int

        ! Loop over all determinants.
        do idet = 1, en_pert%ndets
            ilut(0:NIfDBO) = en_pert%dets(0:NIfDBO,idet)
            call decode_bit_det(nI, ilut)

            if (tHPHF) then
                h_aa = hphf_diag_helement(nI, ilut)
            else
                h_aa = get_helement(nI, nI, 0)
            end if

            call extract_sign_EN(en_pert%sign_length, en_pert%dets(:,idet), contrib)

            do istate = 1, en_pert%sign_length
                contrib_rdm(istate) = contrib(istate)/( (rdm_energy_num(istate)/rdm_norm(istate)) - h_aa )
            end do

            energy_pert = energy_pert + contrib_rdm
        end do

    end subroutine calc_en_pert_energy

    subroutine calc_hermitian_errors(rdm, rdm_recv, spawn, rdm_norm, max_error_herm_all, sum_error_herm_all)

        ! Calculate the hermiticity errors, i.e.
        ! 
        ! \Gamma_{ij,kl} - \Gamma_{kl,ij}^*,
        ! 
        ! which should be equal to zero for an exact 2-RDM.

        ! Specifically we return the largest such error, and the sum of all
        ! such errors, to max_error_herm_all and sum_error_herm_all.

        ! Use the rdm_recv and spawn objects as space for performing a
        ! communication of a 2-RDM. The hermiticity error 2-RDM which is
        ! communicated is defined as
        !
        ! \Gamma^{error}_{ij,kl} = \Gamma_{ij,kl} - \Gamma_{kl,ij}^*,
        !
        ! for ij < kl.

        use hash, only: clear_hash_table
        use Parallel_neci, only: MPIAllReduce, nProcessors
        use ParallelHelper, only: MPI_SUM, MPI_MAX
        use rdm_data_utils, only: add_to_rdm_spawn_t, communicate_rdm_spawn_t_wrapper, annihilate_rdm_list

        type(rdm_list_t), intent(in) :: rdm
        type(rdm_list_t), intent(inout) :: rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn
        real(dp), intent(in) :: rdm_norm(rdm%sign_length)
        real(dp), intent(out) :: max_error_herm_all(rdm%sign_length)
        real(dp), intent(out) :: sum_error_herm_all(rdm%sign_length)

        integer(int_rdm) :: ijkl
        integer :: ielem
        integer :: ij, kl, i, j, k, l ! spin orbitals
        real(dp) :: rdm_sign(rdm%sign_length)
        real(dp) :: max_error_herm(rdm%sign_length), sum_error_herm(rdm%sign_length)
        logical :: nearly_full, finished, all_finished

        ! If we're about to fill up the spawn list, perform a communication.
        nearly_full = .false.
        ! Have we finished adding RDM elements to the spawned list?
        finished = .false.
        rdm_recv%nelements = 0

        ! Clear the spawn object before we use it.
        spawn%free_slots = spawn%init_free_slots(0:nProcessors-1)
        call clear_hash_table(spawn%rdm_send%hash_table)

        ! Loop over all RDM elements.
        do ielem = 1, rdm%nelements
            ! If the spawned list is nearly full, perform a communication.
            if (nearly_full) then
                call communicate_rdm_spawn_t_wrapper(spawn, rdm_recv, finished, all_finished)
                nearly_full = .false.
            end if

            ijkl = rdm%elements(0,ielem)
            ! Obtain spin orbital labels and the RDM element sign.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
            call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

            ! If in the lower half of the RDM, reflect to the upper half and
            ! include with a minus sign.
            if (tGUGA) then 
                call add_to_rdm_spawn_t(spawn, i, j, k, l, rdm_sign, .false., nearly_full)
            else
                if (ij > kl) then
                    call add_to_rdm_spawn_t(spawn, k, l, i, j, -rdm_sign, .false., nearly_full)
                else if (ij < kl) then
                    call add_to_rdm_spawn_t(spawn, i, j, k, l, rdm_sign, .false., nearly_full)
                end if
            end if
        end do

        finished = .true.
        ! Keep performing communications until all RDM spawnings on every
        ! processor have been communicated.
        do
            call communicate_rdm_spawn_t_wrapper(spawn, rdm_recv, finished, all_finished)
            if (all_finished) exit
        end do

        call annihilate_rdm_list(rdm_recv)

        max_error_herm = 0.0_dp
        sum_error_herm = 0.0_dp

        ! Find the largest error and sum of errrors on this processor.
        do ielem = 1, rdm_recv%nelements
            call extract_sign_rdm(rdm_recv%elements(:,ielem), rdm_sign)
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
            if (ij == kl) cycle
            rdm_sign = abs(rdm_sign)
            max_error_herm = max(max_error_herm, rdm_sign)
            sum_error_herm = sum_error_herm + rdm_sign
        end do

        ! The input 2-RDM wasn't normalised, so need to normalise these.
        max_error_herm = max_error_herm/rdm_norm 
        sum_error_herm = sum_error_herm/rdm_norm

       ! Find the largest error and sum of errors across all processors.
       call MPIAllReduce(max_error_herm, MPI_MAX, max_error_herm_all)
       call MPIAllReduce(sum_error_herm, MPI_SUM, sum_error_herm_all)

    end subroutine calc_hermitian_errors

end module rdm_estimators
