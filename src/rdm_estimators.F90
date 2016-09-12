#include "macros.h"

module rdm_estimators

    use bit_rep_data, only: NIfTot
    use constants
    use rdm_data, only: rdm_list_t, rdm_spawn_t
    use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm

    implicit none

contains

    subroutine init_rdm_estimates_t(est, nrdms_standard, nrdms_transition, open_output_file)

        ! Initialise an rdm_estimates_t object. Allocate arrays to be large
        ! enough to hold estimates for nrdms_srandard+nrdms_transition RDMs.

        ! Also, if open_output_file is true, and if this is the processor with
        ! label 0, then open an RDMEstimates file, and write the file's header.

        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_estimates_t
        use util_mod, only: get_free_unit

        type(rdm_estimates_t), intent(out) :: est
        integer, intent(in) :: nrdms_standard, nrdms_transition
        logical, intent(in) :: open_output_file

        integer :: nrdms, ierr

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
            call write_rdm_est_file_header(est%write_unit, nrdms_standard)
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
        if (allocated(est%trace_inst)) deallocate(est%trace_inst, stat=ierr)
        if (allocated(est%norm_inst)) deallocate(est%norm_inst, stat=ierr)
        if (allocated(est%energy_1_num_inst)) deallocate(est%energy_1_num_inst, stat=ierr)
        if (allocated(est%energy_2_num_inst)) deallocate(est%energy_2_num_inst, stat=ierr)
        if (allocated(est%energy_num_inst)) deallocate(est%energy_num_inst, stat=ierr)
        if (allocated(est%spin_num_inst)) deallocate(est%spin_num_inst, stat=ierr)
        if (allocated(est%max_error_herm)) deallocate(est%max_error_herm, stat=ierr)
        if (allocated(est%sum_error_herm)) deallocate(est%sum_error_herm, stat=ierr)

        ! Close the RDMEstimates unit, if it was opened on this processor.
        ! The following was what it was set to if it was not opened in
        ! init_rdm_estimates_t, so don't attempt a close in that case.
        if (est%write_unit /= huge(est%write_unit)) close(est%write_unit)

    end subroutine dealloc_rdm_estimates_t

    subroutine write_rdm_est_file_header(write_unit, nrdms_standard)

        ! Open a new RDMEstimates file (overwriting any existing file), and
        ! write a header to it, appropriate for when we are sampling nrdms RDMs.

        integer, intent(in) :: write_unit, nrdms_standard

        integer :: irdm

        open(write_unit, file='RDMEstimates', status='unknown', position='append')

        write(write_unit, '("#", 4X, "Iteration")', advance='no')

        do irdm = 1, nrdms_standard
            write(write_unit, '(4x,"Energy numerator",1x,i2)', advance='no') irdm
            write(write_unit, '(4x,"Spin^2 numerator",1x,i2)', advance='no') irdm
            write(write_unit, '(7x,"Normalisation",1x,i2)', advance='no') irdm
        end do

        write(write_unit,'()')

    end subroutine write_rdm_est_file_header

    subroutine calc_2rdm_estimates_wrapper(rdm_defs, est, rdm)

        ! Calculate the estimates for the 2-RDM stored in rdm. The full estimates
        ! are stored using this object, and also instantaneous estimates. The
        ! instantaneous estimates are calculated by subtracting the previous
        ! stored values from the newly calculated ones. This works so long as
        ! the estimates are linear functions of the RDMs, which they will be
        ! for any observable.

        use Parallel_neci, only: MPISumAll
        use rdm_data, only: rdm_estimates_t, rdm_list_t, rdm_definitions_t
        use SystemData, only: nel, ecore

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_estimates_t), intent(inout) :: est
        type(rdm_list_t), intent(in) :: rdm

        integer :: irdm
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

    end subroutine calc_2rdm_estimates_wrapper

    subroutine write_rdm_estimates(rdm_defs, est, final_output, write_to_separate_file)

        ! Write RDM estimates to the RDMEstimates file. Specifically, the
        ! numerator of the energy and spin^2 estimators are output, as is the
        ! trace (the denominator of these estimator).

        ! Also, if final_output is true, then output the final total estimates
        ! to standard output, and close the RDMEstimates unit.

        use FciMCData, only: Iter, PreviousCycles
        use LoggingData, only: tRDMInstEnergy
        use rdm_data, only: rdm_estimates_t, rdm_definitions_t
        use util_mod, only: int_fmt

        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_estimates_t), intent(in) :: est
        logical, intent(in) :: final_output, write_to_separate_file

        integer :: irdm

        if (write_to_separate_file) then
            if (tRDMInstEnergy) then
                write(est%write_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
                do irdm = 1, est%nrdms_standard
                    write(est%write_unit, '(3(3x,es20.13))', advance='no') &
                        est%energy_num_inst(irdm), est%spin_num_inst(irdm), est%norm_inst(irdm)
                end do
                write(est%write_unit,'()')
            else
                write(est%write_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
                do irdm = 1, est%nrdms_standard
                    write(est%write_unit, '(3(3x,es20.13))', advance='no') &
                        est%energy_num(irdm), est%spin_num(irdm), est%norm(irdm)
                end do
                write(est%write_unit, '()')
            end if

            call neci_flush(est%write_unit)
        end if

        if (final_output) then
            ! Banner for the start of the 2-RDM section in output.
            write(6,'(1x,2("="),1x,"INFORMATION FOR FINAL 2-RDMS",1x,57("="),/)')

            do irdm = 1, est%nrdms_standard
                write(6,'(1x,"2-RDM ESTIMATES FOR STATE",1x,'//int_fmt(irdm)//',":",)') irdm

                write(6,'(1x,"Trace of 2-el-RDM before normalisation:",1x,es17.10)') est%trace(irdm)
                write(6,'(1x,"Trace of 2-el-RDM after normalisation:",1x,es17.10)') est%trace(irdm)/est%norm(irdm)

                write(6,'(1x,"Energy contribution from the 1-RDM:",1x,es17.10)') est%energy_1_num(irdm)/est%norm(irdm)
                write(6,'(1x,"Energy contribution from the 2-RDM:",1x,es17.10)') est%energy_2_num(irdm)/est%norm(irdm)
                write(6,'(1x,"*TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*:",1x,es20.13,/)') &
                    est%energy_num(irdm)/est%norm(irdm)

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

    subroutine calc_rdm_trace(rdm, rdm_trace)

        ! Calculate trace of the 2-RDM in the rdm object, and output it to
        ! rdm_trace.

        ! This trace is defined as
        !
        ! rdm_trace = \sum_{ij} \Gamma_{ij,ij},
        !
        ! where \Gamma_{ij,kl} is the 2-RDM stored in rdm, and i and j are
        ! spin orbital labels.

        use rdm_data, only: rdm_spawn_t

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(out) :: rdm_trace(rdm%sign_length)

        integer(int_rdm) :: ijkl
        integer :: ielem
        integer :: ij, kl, i, j, k, l ! spin orbitals
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_trace = 0.0_dp

        ! Loop over all RDM elements.
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

            ! If this is a diagonal element, add the element to the trace.
            if (ij == kl) then
                call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)
                rdm_trace = rdm_trace + rdm_sign
            end if
        end do

    end subroutine calc_rdm_trace

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
            if (ij > kl) then
                call add_to_rdm_spawn_t(spawn, k, l, i, j, -rdm_sign, .false., nearly_full)
            else if (ij < kl) then
                call add_to_rdm_spawn_t(spawn, i, j, k, l, rdm_sign, .false., nearly_full)
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
