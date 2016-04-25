#include "macros.h"

! Routines to perform final operations on RDMs, which includes calculating
! 1-RDMs from 2-RDMs, apply symmetries, and writing RDMs out.

module rdm_finalising

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use rdm_data, only: rdm_list_t, rdm_spawn_t, one_rdm_t
    use rdm_data_utils, only: calc_separate_rdm_labels, extract_sign_rdm, add_to_rdm_spawn_t
    use rdm_data_utils, only: communicate_rdm_spawn_t
    use util_mod

    implicit none

contains

    subroutine finalise_rdms(one_rdms, two_rdms, rdm_recv, spawn, rdm_estimates)

        ! Wrapper routine, called at the end of a simulation, which in turn
        ! calls all required finalisation routines.

#ifdef _MOLCAS_
        use EN2MOLCAS, only: NECI_E
#endif
        use LoggingData, only: tBrokenSymNOs, occ_numb_diff, RDMExcitLevel, tExplicitAllRDM
        use LoggingData, only: tPrint1RDM, tDiagRDM, tDumpForcesInfo, tDipoles
        use Parallel_neci, only: iProcIndex, MPIBarrier, MPIBCast
        use rdm_data, only: tRotatedNos, FinaliseRDMs_Time, tOpenShell, print_2rdm_est
        use rdm_data, only: rdm_list_t, rdm_spawn_t, one_rdm_t, rdm_estimates_t
        use rdm_estimators, only: calc_2rdm_estimates_wrapper, write_rdm_estimates
        use rdm_nat_orbs, only: find_nat_orb_occ_numbers, BrokenSymNo
        use util_mod, only: set_timer, halt_timer

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        type(rdm_list_t), intent(inout) :: two_rdms, rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn
        type(rdm_estimates_t), intent(inout) :: rdm_estimates

        integer :: i, ierr
        real(dp) :: norm_1rdm, trace_1rdm, SumN_Rho_ii

        call set_timer(FinaliseRDMs_Time)

        if (tExplicitAllRDM) then
            write(6,'(/,"**** RDMs CALCULATED EXPLICITLY ****",1X,/)')
        else
            write(6,'(/,"**** RDMs CALCULATED STOCHASTICALLY ****",1X,/)')
        end if

        ! Combine the 1- or 2-RDM from all processors, etc.

        ! Stuff using the 2-RDMs:
        if (RDMExcitLevel /= 1) then
            ! We always want to calculate one final RDM energy, whether or not we're
            ! calculating the energy throughout the calculation.
            ! Unless of course, only the 1-RDM is being calculated.

            ! Calculate the RDM estmimates from the final few iterations,
            ! since it was last calculated. But only do this if they're
            ! actually to be printed.
            if (print_2rdm_est) then
                call calc_2rdm_estimates_wrapper(rdm_estimates, two_rdms)
            end if

            ! Output the final 2-RDMs themselves, in all forms desired.
            call output_2rdm_wrapper(rdm_estimates, two_rdms, rdm_recv, spawn)

            ! Calculate the 1-RDMs from the 2-RDMS, if required.
            if (tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
                call calc_1rdms_from_2rdms(one_rdms, two_rdms, rdm_estimates%norm, tOpenShell)
            end if
        end if

        ! Stuff using the 1-RDMs:
        do i = 1, size(one_rdms)
            if (RDMExcitLevel == 1) then
                call Finalise_1e_RDM(one_rdms(i)%matrix, one_rdms(i)%rho_ii, i, norm_1rdm, .false.)
            else
                if (tPrint1RDM) then
                    call Finalise_1e_RDM(one_rdms(i)%matrix, one_rdms(i)%rho_ii, i, norm_1rdm, .false.)
                else if (tDiagRDM .and. iProcIndex == 0) then
                    call calc_1e_norms(one_rdms(i)%matrix, one_rdms(i)%rho_ii, trace_1rdm, norm_1rdm, SumN_Rho_ii)
                    write(6,'(/,1X,"SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS:",1X,F20.13)') SumN_Rho_ii
                end if
            end if

            call MPIBarrier(ierr)

            ! Call the routines from NatOrbs that diagonalise the one electron
            ! reduced density matrix.
            tRotatedNOs = .false. ! Needed for BrokenSymNo routine
            if (tDiagRDM) call find_nat_orb_occ_numbers(one_rdms(i), i)

            ! After all the NO calculations are finished we'd like to do another
            ! rotation to obtain symmetry-broken natural orbitals
            if (tBrokenSymNOs) then
                call BrokenSymNO(one_rdms(i)%evalues, occ_numb_diff)
            end if
        end do

        ! Write the final instantaneous RDM estimates, and also the final
        ! report of the total RDM estimates.
        if (iProcIndex == 0) call write_rdm_estimates(rdm_estimates, .true., print_2rdm_est)
#ifdef _MOLCAS_
        if (print_2rdm_est) then
            NECI_E = rdm_estimates%rdm_energy_tot_accum(1)
            call MPIBarrier(ierr)
            call MPIBCast(NECI_E)
            write(6,*) 'NECI_E at rdm_general.f90 ', NECI_E
        end if
#endif

        call halt_timer(FinaliseRDMs_Time)

    end subroutine finalise_rdms

    subroutine output_2rdm_wrapper(est, rdm, rdm_recv, spawn)

        ! Call routines to output RDMs in all requested forms.

        ! We also calculate the Hermitian errors here, since this is something
        ! we typically want to do at the same point (the very end of a
        ! simulation usually), and this requires large parallel
        ! communications, as does the printing.

        use LoggingData, only: tWrite_normalised_RDMs, tWriteSpinFreeRDM, tWrite_RDMs_to_read
        use rdm_data, only: rdm_estimates_t, rdm_list_t, rdm_spawn_t, tOpenShell
        use rdm_estimators, only: calc_hermitian_errors

        type(rdm_estimates_t), intent(inout) :: est
        ! IMPORTANT: rdm is not actually modified by this routine, despite
        ! needing inout status.
        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_list_t), intent(inout) :: rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn

        call calc_hermitian_errors(rdm, rdm_recv, spawn, est%norm, est%max_error_herm, est%sum_error_herm)

        if (tWriteSpinFreeRDM) call print_spinfree_2rdm_wrapper(rdm, rdm_recv, spawn, est%norm)
        if (tWrite_Normalised_RDMs) call print_rdms_spin_sym_wrapper(rdm, rdm_recv, spawn, est%norm, tOpenShell)
        if (tWrite_RDMs_to_read) call print_rdm_popsfile(rdm)

        !call read_rdm_popsfile(rdm, spawn)

    end subroutine output_2rdm_wrapper

    subroutine calc_1rdms_from_spinfree_2rdms(one_rdms, two_rdms, rdm_trace)

        ! For each 2-RDM in two_rdms, calculate the corresponding spin-free
        ! 1-RDM:
        !
        ! \gamma^{spinfree}_{p,q} = \frac{1}{N-1} \sum_a \Gamma^{spinfree}_{pa,qa}
        !
        ! Here, p, q and a are spatial labels. N is the number of electrons.

        ! The spinfree 1-RDM is defined in terms of the spinned 1-RDM by:
        !
        ! \gamma^{spinfree}_{p,q} = \gamma_{p\alpha,q\alpha) + \gamma_{p\beta,q\beta)

        ! IMPORTANT: This routine should *only* be used by taking *spin-free*
        ! 2-RDMs as the input. Specifically, it takes spin-free RDMs as returned
        ! by the create_spinfree_2rdm routine, which does *not* restrict the
        ! labels allowed. Inputting 2-RDMs in other will give incorrect results!

        ! The output 1-RDM elements are sorted in the standard form: elements
        ! are indexed using the SymLabelListInv_rot array, so that the 1-RDMs
        ! will be in block-diagonal form, with elements within each symmetry
        ! block stored together.

        use Parallel_neci, only: MPISumAll
        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: nel
        use UMatCache, only: spatial

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        type(rdm_list_t), intent(in) :: two_rdms
        real(dp), intent(in) :: rdm_trace(:)

        integer(int_rdm) :: pqrs
        integer :: ielem, irdm, ierr
        integer :: pq, rs, p, q, r, s
        real(dp) :: rdm_sign(two_rdms%sign_length)
        real(dp), allocatable :: temp_rdm(:,:)

        do irdm = 1, size(one_rdms)
            one_rdms(irdm)%matrix = 0.0_dp
        end do

        ! Loop over all elements of the 2-RDM, \Gamma_{pq,rs}, where p, q, r
        ! and s are spatial labels. If at least two spatial indices are the
        ! same then we have a contribution to the 1-RDM.
        do ielem = 1, two_rdms%nelements
            pqrs = two_rdms%elements(0,ielem)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(pqrs, pq, rs, r, s, q, p)

            call extract_sign_rdm(two_rdms%elements(:,ielem), rdm_sign)

            associate(ind => SymLabelListInv_rot)
                ! An element of the form \Gamma_{pa,ra}.
                if (q == s) then
                    do irdm = 1, size(one_rdms)
                        one_rdms(irdm)%matrix(ind(p), ind(r)) = one_rdms(irdm)%matrix(ind(p), ind(r)) + rdm_sign(irdm)
                    end do
                end if
            end associate
        end do

        ! Allocate a temporary array in which to receive the MPI communication.
        allocate(temp_rdm(size(one_rdms(1)%matrix,1), size(one_rdms(1)%matrix,2)), stat=ierr)

        ! Perform a sum over all processes, for each 1-RDM being sampled.
        do irdm = 1, size(one_rdms)
            call MPISumAll(one_rdms(irdm)%matrix, temp_rdm)
            ! Copy summed RDM back to the main array, and normalise.
            one_rdms(irdm)%matrix = temp_rdm / (rdm_trace(irdm)*real(nel-1,dp))
        end do

        deallocate(temp_rdm, stat=ierr)

    end subroutine calc_1rdms_from_spinfree_2rdms

    subroutine calc_1rdms_from_2rdms(one_rdms, two_rdms, rdm_trace, open_shell)

        ! For each 2-RDM in two_rdms, if open_shell is true then calculate the
        ! full spinned 1-RDM, otherwise calculate the spinfree 1-RDM. The
        ! former case is defined by:
        !
        ! \gamma_{i,j} = \frac{1}{N-1} \sum_k \Gamma_{ik,jk}
        !
        ! Here, i, j and k are spatial labels. N is the number of electrons.
        !
        ! The spinfree case is then a contraction over the spin labels of the
        ! spinned 1-RDM:
        !
        ! \gamma^{spinfree}_{p,q} = \gamma_{p\alpha,q\alpha) + \gamma_{p\beta,q\beta)
        !
        ! where p and q are spatial labels.

        ! The output 1-RDM elements are sorted in the standard form: elements
        ! are indexed using the SymLabelListInv_rot array, so that the 1-RDMs
        ! will be in block-diagonal form, with elements within each symmetry
        ! block stored together.

        use Parallel_neci, only: MPISumAll
        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: nel
        use UMatCache, only: spatial

        type(one_rdm_t), intent(inout) :: one_rdms(:)
        type(rdm_list_t), intent(in) :: two_rdms
        real(dp), intent(in) :: rdm_trace(:)
        logical, intent(in) :: open_shell

        integer(int_rdm) :: ijkl
        integer :: ielem, irdm, ierr
        integer :: ij, kl, i, j, k, l
        integer :: p, q, r, s
        real(dp) :: rdm_sign(two_rdms%sign_length)
        real(dp), allocatable :: temp_rdm(:,:)

        do irdm = 1, size(one_rdms)
            one_rdms(irdm)%matrix = 0.0_dp
        end do

        ! Loop over all elements of the 2-RDM, \Gamma_{pq,rs}, where p, q, r
        ! and s are spatial labels. If at least two spatial indices are the
        ! same then we have a contribution to the 1-RDM.
        do ielem = 1, two_rdms%nelements
            ijkl = two_rdms%elements(0,ielem)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

            ! For closed shell systems we work with spatial orbitals, to
            ! calculate spin-free 1RDMs.
            if (open_shell) then
                p = i; q = j;
                r = k; s = l;
            else
                p = spatial(i); q = spatial(j);
                r = spatial(k); s = spatial(l);
            end if

            call extract_sign_rdm(two_rdms%elements(:,ielem), rdm_sign)

            ! If abba or baab term - swap last two indices and sign.
            if (.not. same_spin(i,k)) then
                if (open_shell) then
                    r = l; s = k;
                else
                    r = spatial(l); s = spatial(k);
                end if
                rdm_sign = -rdm_sign
            end if

            associate(ind => SymLabelListInv_rot)

                ! An element of the form \Gamma_{aq,as}.
                if (p == r) then
                    do irdm = 1, size(one_rdms)
                        one_rdms(irdm)%matrix(ind(q), ind(s)) = one_rdms(irdm)%matrix(ind(q), ind(s)) + rdm_sign(irdm)
                    end do
                end if
                ! An element of the form \Gamma_{pa,ra}.
                if (q == s) then
                    do irdm = 1, size(one_rdms)
                        one_rdms(irdm)%matrix(ind(p), ind(r)) = one_rdms(irdm)%matrix(ind(p), ind(r)) + rdm_sign(irdm)
                    end do
                end if

                ! The below cases give contributions by swapping one pair of
                ! indices. Only include these contributions if we have aaaa or
                ! bbbb terms. This because if we had a term with spin signature
                ! abab (for example), then swapping as below would give abba
                ! or baab terms, which don't contribute to the 1-RDM.
                if (same_spin(k,l)) then
                    ! An element of the form \Gamma_{pa,as}.
                    if (p == s) then
                        do irdm = 1, size(one_rdms)
                            one_rdms(irdm)%matrix(ind(q), ind(r)) = one_rdms(irdm)%matrix(ind(q), ind(r)) - rdm_sign(irdm)
                        end do
                    end if
                    ! An element of the form \Gamma_{aq,ra}.
                    if (q == r) then
                        do irdm = 1, size(one_rdms)
                            one_rdms(irdm)%matrix(ind(p), ind(s)) = one_rdms(irdm)%matrix(ind(p), ind(s)) - rdm_sign(irdm)
                        end do
                    end if
                end if

            end associate
        end do

        ! Allocate a temporary RDM array.
        allocate(temp_rdm(size(one_rdms(1)%matrix,1), size(one_rdms(1)%matrix,2)), stat=ierr)

        ! Make every RDM symmetric. This could have been done when adding
        ! contribution in above, but hopefully the code will be clearer if
        ! done here.
        do irdm = 1, size(one_rdms)
            ! Use temp_rdm as temporary space for the transpose, to (hopefully)
            ! prevent a temporary array being created in the sum below.
            temp_rdm = transpose(one_rdms(irdm)%matrix)
            one_rdms(irdm)%matrix = (one_rdms(irdm)%matrix + temp_rdm)/2.0_dp
        end do

        ! Perform a sum over all processes, for each 1-RDM being sampled.
        do irdm = 1, size(one_rdms)
            call MPISumAll(one_rdms(irdm)%matrix, temp_rdm)
            ! Copy summed RDM back to the main array, and normalise.
            one_rdms(irdm)%matrix = temp_rdm / (rdm_trace(irdm)*real(nel-1,dp))
        end do

        deallocate(temp_rdm, stat=ierr)

    end subroutine calc_1rdms_from_2rdms

    subroutine make_hermitian_rdm(rdm, spawn, rdm_recv)

        ! Take the RDM in the rdm object, and output a new RDM which is the
        ! same but with Hermiticy applied to it, i.e., the elements above and
        ! below the diagonal are averaged appropriately.

        ! If rdm_recv is input then the new RDM will be output to this object.
        ! If not, then the RDM in the rdm object will be overwritten. However,
        ! the hash table in these objects will *not* be updated.

        use rdm_data_utils, only: annihilate_rdm_list

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn
        type(rdm_list_t), optional, intent(inout) :: rdm_recv

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        integer :: p_temp, q_temp
        real(dp) :: rdm_sign(rdm%sign_length)

        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
            call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

            ! Factor of a half to account for prevent double-counting, and
            ! instead average elements from above and below the diagonal.
            if (pq /= rs) rdm_sign = 0.5_dp*rdm_sign

            ! If in the lower half of the RDM, reflect to the upper half.
            if (pq > rs) then
                call add_to_rdm_spawn_t(spawn, r, s, p, q, rdm_sign, .false.)
            else
                call add_to_rdm_spawn_t(spawn, p, q, r, s, rdm_sign, .false.)
            end if
        end do

        if (present(rdm_recv)) then
            call communicate_rdm_spawn_t(spawn, rdm_recv)
            call annihilate_rdm_list(rdm_recv)
        else
            call communicate_rdm_spawn_t(spawn, rdm)
            call annihilate_rdm_list(rdm)
        end if

    end subroutine make_hermitian_rdm

    subroutine apply_symmetries_for_output(rdm, spawn, open_shell, rdm_recv)

        ! This routine will take in rdm, and output a new rdm which will have
        ! all appropriate symmetries applied so that the latter RDM can be
        ! passed to the routine to write RDMs.

        ! If rdm_recv is input then the new RDM will be output to this object.
        ! If not, then the RDM in the rdm object will be overwritten. However,
        ! the hash table in these objects will *not* be updated.

        ! The input RDM should already have hermiticy symmetry applied to it.

        ! WARNING: To clarify potential confusion, we point out that this
        ! routine also applies hermiticy again, but *only* for the spatial
        ! labels, not the full spin labels. It does so specifically using a
        ! particular legacy ordering. This is not a mistake - for open shell
        ! systems we do need to apply full hermiticy, but also need to apply spatial
        ! hermiticy because spatial labels are written within each output file.

        use rdm_data_utils, only: annihilate_rdm_list

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn
        logical, intent(in) :: open_shell
        type(rdm_list_t), optional, intent(inout) :: rdm_recv

        integer(int_rdm) :: ijkl
        integer :: ielem, ij, kl, i, j, k, l
        integer :: p, q, r, s
        integer :: pq_legacy, rs_legacy
        real(dp) :: rdm_sign(rdm%sign_length)

        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0,ielem)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)
            call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)

            ! When there are two elements which are guaranteed to be exactly the
            ! same, we usually only want to print one of them (and to average
            ! over the equal terms). This function returns the labels we want.
            call apply_legacy_output_ordering(i, j, k, l, rdm_sign, pq_legacy, rs_legacy)

            ! For closed shell systems, want bbbb -> aaaa, baba -> abab,
            ! baab -> abba, by flipping all spins, which is a symmetry for such
            ! systems. If the first label has beta spin then we definitely want
            ! to move this RDM element to the equivalent flipped term so go ahead
            ! and flip all the spins.
            if (.not. open_shell) then
                if (is_beta(i)) then
                    ! The ab_pair macro swaps alpha and beta spins of a label
                    ! while keeping the spatial orbital unchanged.
                    i = ab_pair(i)
                    j = ab_pair(j)
                    k = ab_pair(k)
                    l = ab_pair(l)
                end if
                ! If the spatial parts of i and j are the same, and the spatial
                ! parts of k and l are *also* the same, then the RDM element won't
                ! have been added into both equivalent spin-flipped arrays
                ! because i<j and k<l is enforced), so we don't count twice.
                if (.not. (is_in_pair(i,j) .and. is_in_pair(k,l))) then
                    ! Also, if (i,j) and (k,l) have the same spatial parts, but
                    ! different spin parts ((alpha,beta) and (beta,alpha), or
                    ! vice versa) then they only occur once, again because we
                    ! enforce i<j and k<l for all stored RDM elements.
                    if (.not. (pq_legacy == rs_legacy .and. (.not. ij == kl))) then
                        rdm_sign = rdm_sign*0.5_dp
                    end if
                end if
            end if

            call add_to_rdm_spawn_t(spawn, i, j, k, l, rdm_sign, .false.)

            if (open_shell) then
                ! For open shell systems, if i and j have the same spatial parts,
                ! and k and l do too, then we only have baba spin signature,
                ! (because we enforce i<j, k<l) but we'd like to print out abab too.
                if (is_in_pair(i,j) .and. is_in_pair(k,l)) then
                    call add_to_rdm_spawn_t(spawn, j, i, l, k, rdm_sign, .false.)
                end if

                ! Because we enforce hermiticy symmetry in the output, we would
                ! only print the following terms with baab. We want to print it
                ! with abba too here, so do that.
                if (pq_legacy == rs_legacy .and. (.not. ij == kl)) then
                    call add_to_rdm_spawn_t(spawn, k, l, i, j, rdm_sign, .false.)
                end if
            end if
        end do

        if (present(rdm_recv)) then
            call communicate_rdm_spawn_t(spawn, rdm_recv)
            call annihilate_rdm_list(rdm_recv)
        else
            call communicate_rdm_spawn_t(spawn, rdm)
            call annihilate_rdm_list(rdm)
        end if

    end subroutine apply_symmetries_for_output

    pure subroutine apply_legacy_output_ordering(i, j, k, l, rdm_sign, pq_legacy, rs_legacy)

        ! Enforce the symmetries of RDMs to only keep certain combinations of
        ! i, j, k, l spin labels, where a redundancy exists. Whenever we have
        ! an unused combination, flip/swap labels (and the sign if necessary).

        ! For example the 2-RDM is hermitian, so if ij /= kl, then we only need
        ! to print either \Gamma_{ij,kl} or \Gamma{kl,ij}, but not both. Which
        ! combinations we decide to print is decided below, which is purely a
        ! legacy decision (as far as I know!). See comments below for defintions
        ! of what we keep.

        use SystemData, only: nbasis
        use UMatCache, only: spatial

        integer, intent(inout) :: i, j, k, l
        real(dp), intent(inout) :: rdm_sign(:)
        integer, intent(out) :: pq_legacy, rs_legacy

        integer :: p, q, r, s
        integer :: i_temp, j_temp

        ! RDMs are output in files labelled by their spin signatures:
        ! aaaa, abab, abba, bbbb, baba or baab.
        ! Within each file, therefore, only spatial orbital labels are printed.
        ! Thus, we need to use spatial orbitals to determine which RDM elements
        ! are  to kept, and which transformed.
        p = spatial(i); q = spatial(j);
        r = spatial(k); s = spatial(l);

        ! When we calculate the combined labels, pq and rs, we would
        ! usually have p and q swapped below, and similarly with r and s.
        ! However, the old RDM files prints only RDM elements with pq < rs,
        ! where pq and rs are defined as follows.
        pq_legacy = (q-1)*nbasis + p
        rs_legacy = (s-1)*nbasis + r

        ! Apply symmetry (for *real* RDMs), to only print elements from one
        ! half of the RDM, using the legacy ordering.
        if (pq_legacy > rs_legacy) then
            i_temp = i; j_temp = j;
            i = k; j = l;
            k = i_temp; l = j_temp;
        end if

        ! If either i and j have the same spatial part, of k and l have the
        ! same spatial part, and we have a spin signature with 2 alphas and
        ! 2 betas, then the convention is to output it as either abab or
        ! baba, but *not* as abba or baab. If we have abba or baab in this
        ! case then we have to swap two indices and introduce a minus sign.
        ! Because we enforce i<j and k<l in all RDM elements, there are only
        ! two possibilities to consider:
        if (is_in_pair(i,j) .and. is_beta(i) .and. is_alpha(j) .and. &
                is_alpha(k) .and. is_beta(l)) then
            i = ab_pair(i)
            j = ab_pair(j)
            rdm_sign = -rdm_sign
        else if (is_in_pair(k,l) .and. is_alpha(i) .and. is_beta(j) .and. &
                 is_beta(k) .and. is_alpha(l)) then
            k = ab_pair(k)
            l = ab_pair(l)
            rdm_sign = -rdm_sign
        end if

    end subroutine apply_legacy_output_ordering

    subroutine print_rdms_spin_sym_wrapper(rdm, rdm_recv, spawn, rdm_trace, open_shell)

        ! Compress the full spinned-RDMs by summing over spin-equivalent terms
        ! (i.e. aaaa and bbbb rdms), and also applying symmetry of (*real*)
        ! RDMs. The result will be stored in rdm_recv. Then, print it out to a
        ! file.

        ! IMPORTANT: Although the rdm object has inout status, it will *not*
        ! be modified. The inout status is to allow for the optional possibility
        ! of updating the first argument of make_hermitian_rdm, which is not
        ! used here.

        use hash, only: clear_hash_table

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_list_t), intent(inout) :: rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        logical, intent(in) :: open_shell

        spawn%free_slots = spawn%init_free_slots
        call clear_hash_table(spawn%rdm_send%hash_table)

        call make_hermitian_rdm(rdm, spawn, rdm_recv)

        call apply_symmetries_for_output(rdm_recv, spawn, open_shell)
        call print_rdms_with_spin(rdm_recv, rdm_trace)

    end subroutine print_rdms_spin_sym_wrapper

    subroutine create_spinfree_2rdm(rdm, spawn, rdm_recv)

        ! Take an standard (spinned) 2-RDM, stored in rdm, and output the
        ! spinfree version of it to the spawn%rdm_recv object.

        ! The input RDM has elements equal to:
        !
        ! \Gamma_{ij,kl} = < a^+_i a^+_j a_l a_k >
        !
        ! where i, j, k and l are spin orbital labels, and the output spinfree
        ! RDM has elements equal to:
        !
        ! \Gamma^{spinfree}_{pq,rs} = \sum_{x,y} < a^+_{p,x} a^+_{q,y} a_{s,y} a_{r,x} >
        !
        ! where p, q, r, s are spatial orbital labels, and x and y are spin
        ! labels (alpha or beta) which are summed over.

        ! Thus, all terms with spin signature aaaa, abab, bbbb or baba are
        ! summed together. Terms with spin signature abba or baab have their
        ! final two spin orbital labels swapped (introducing a minus sign), so
        ! that they give a contribution to the resulting spinfree RDM element.

        use rdm_data_utils, only: annihilate_rdm_list
        use SystemData, only: nbasis
        use UMatCache, only: spatial

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn
        type(rdm_list_t), optional, intent(inout) :: rdm_recv

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        integer :: pq_spat, rs_spat
        integer :: p_spat, q_spat, r_spat, s_spat
        integer :: r_orig, s_orig
        real(dp) :: rdm_sign(rdm%sign_length)

        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
            call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

            ! Store the original labels, before we possibly swap them.
            r_orig = r; s_orig = s;

            ! If this term is abba or baab then we can make it abab or baba by
            ! swapping the last two indices, which introduces a minus sign.
            ! It will then contribute to a spinfree 2-RDM element.
            if (.not. same_spin(p,r)) then
                s = r_orig
                r = s_orig
                rdm_sign = -rdm_sign
            end if

            ! Get the spatial orbital labels from the spin orbital ones.
            p_spat = spatial(p); q_spat = spatial(q);
            r_spat = spatial(r); s_spat = spatial(s);
            ! The 'combined' labels.
            pq_spat = (p_spat-1)*nbasis + q_spat
            rs_spat = (r_spat-1)*nbasis + s_spat

            ! If the RDM is not symmetrised then the same term will be added
            ! from both above below the diagonal, so in this case we want a
            ! factor of a half to average and not double count.
            if (pq_spat /= rs_spat) rdm_sign = rdm_sign*0.5_dp

            ! Due to the fact that RDM elements are only stored with p < q and
            ! r < s, the following terms are only stored with baba spin, never
            ! with abab. Double this term to make up for it.
            if (p_spat == q_spat .and. r_spat == s_spat) rdm_sign = 2.0_dp*rdm_sign

            ! Add all spinfree 2-RDM elements corresponding to these labels.
            call add_rdm_elements(p_spat, q_spat, r_spat, s_spat, rdm_sign, spawn)

            ! If this is an aaaa or bbbb term then *minus* this RDM element will
            ! be equal to the equivalent RDM element with the last two labels
            ! swapped. So, add this contribution into that RDM element. We
            ! don't have to do this, but doing so applies some extra averaging.
            ! Want to apply all the averaging possible over equivalent elements.
            if (same_spin(p, q)) then
                ! Re-extract sign in case it has been modified.
                call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

                ! Swap the spatial labels.
                r_spat = spatial(s_orig); s_spat = spatial(r_orig);
                rs_spat = (r_spat-1)*nbasis + s_spat

                if (pq_spat /= rs_spat) rdm_sign = rdm_sign*0.5_dp
                rdm_sign = -rdm_sign

                call add_rdm_elements(p_spat, q_spat, r_spat, s_spat, rdm_sign, spawn)
            end if

        end do

        if (present(rdm_recv)) then
            call communicate_rdm_spawn_t(spawn, rdm_recv)
            call annihilate_rdm_list(rdm_recv)
        else
            call communicate_rdm_spawn_t(spawn, rdm)
            call annihilate_rdm_list(rdm)
        end if

    contains

        subroutine add_rdm_elements(p_spat, q_spat, r_spat, s_spat, rdm_sign, spawn)

            ! Add in the single contribution rdm_sign to the following elements
            ! of the spinfree 2-RDM:
            !
            ! \Gamma^{spinfree}_{pq,rs} = \sum_{x,y} < a^+_{p,x} a^+_{q,y} a_{s,y} a_{r,x} >
            ! \Gamma^{spinfree}_{qp,sr} = \sum_{x,y} < a^+_{q,x} a^+_{p,y} a_{r,y} a_{s,x} >
            ! \Gamma^{spinfree}_{rs,pq} = \sum_{x,y} < a^+_{r,x} a^+_{s,y} a_{q,y} a_{p,x} >
            ! \Gamma^{spinfree}_{sr,qp} = \sum_{x,y} < a^+_{s,x} a^+_{r,y} a_{p,y} a_{q,x} >
            !
            ! where x and y are spin labels which are summed over in the final
            ! result.
            !
            ! For a *REAL* spinfree 2-RDM, all of these elements are rigorously
            ! equal, so it is appropriate that we add all contributions in
            ! together like this.
            !
            ! The if-statements in here prevent adding to the same RDM element
            ! twice.

            integer, intent(in) :: p_spat, q_spat, r_spat, s_spat
            real(dp), intent(in) :: rdm_sign(:)
            type(rdm_spawn_t), intent(inout) :: spawn

            ! RDM element \Gamma_{pq,rs}.
            call add_to_rdm_spawn_t(spawn, p_spat, q_spat, r_spat, s_spat, rdm_sign, .true.)

            ! RDM element \Gamma_{qp,sr}.
            if (.not. (p_spat == q_spat .and. r_spat == s_spat)) then
                call add_to_rdm_spawn_t(spawn, q_spat, p_spat, s_spat, r_spat, rdm_sign, .true.)
            end if

            if (pq_spat /= rs_spat) then
                ! RDM element \Gamma_{rs,pq}.
                call add_to_rdm_spawn_t(spawn, r_spat, s_spat, p_spat, q_spat, rdm_sign, .true.)

                ! RDM element \Gamma_{sr,qp}.
                if (.not. (p_spat == q_spat .and. r_spat == s_spat)) then
                    call add_to_rdm_spawn_t(spawn, s_spat, r_spat, q_spat, p_spat, rdm_sign, .true.)
                end if
            end if

        end subroutine add_rdm_elements

    end subroutine create_spinfree_2rdm

    subroutine print_spinfree_2rdm_wrapper(rdm, rdm_recv, spawn, rdm_trace)

        ! IMPORTANT: Although the rdm object has inout status, it will *not*
        ! be modified. The inout status is to allow for the optional possibility
        ! of updating the first argument of create_spinfree_2rdm, which is not
        ! used here.

        use hash, only: clear_hash_table

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_list_t), intent(inout) :: rdm_recv
        type(rdm_spawn_t), intent(inout) :: spawn
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)

        spawn%free_slots = spawn%init_free_slots
        call clear_hash_table(spawn%rdm_send%hash_table)

        call create_spinfree_2rdm(rdm, spawn, rdm_recv)
        call print_spinfree_2rdm(rdm_recv, rdm_trace)

    end subroutine print_spinfree_2rdm_wrapper

    subroutine print_rdm_popsfile(rdm)

        use Parallel_neci, only: MPIBarrier
        use util_mod, only: get_free_unit

        type(rdm_list_t), intent(inout) :: rdm

        integer :: ielem, iproc, ierr, pops_unit

        do iproc = 0, nProcessors-1

            if (iproc == iProcIndex) then
                ! Let the first processor clear the file to start with.
                if (iproc == 0) then
                    pops_unit = get_free_unit()
                    open(pops_unit, file='RDM_POPSFILE', status='replace', form='unformatted')
                    ! Let the first processor start by printing the number of
                    ! RDMs being sampled.
                    write(pops_unit) rdm%sign_length
                else
                    pops_unit = get_free_unit()
                    open(pops_unit, file='RDM_POPSFILE', status='old', position='append', form='unformatted')
                end if

                do ielem = 1, rdm%nelements
                    write(pops_unit) rdm%elements(:,ielem)
                end do

                close(pops_unit)
            end if

            ! Wait for the current processor to finish printing its RDM elements.
            call MPIBarrier(ierr)
        end do

    end subroutine print_rdm_popsfile

    subroutine read_rdm_popsfile(rdm, spawn)

        use hash, only: clear_hash_table, fill_in_hash_table
        use hash, only: FindWalkerHash, add_hash_table_entry
        use Parallel_neci, only: MPIBarrier
        use SystemData, only: nbasis
        use util_mod, only: get_free_unit

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(in) :: spawn

        integer(int_rdm) :: ijkl, rdm_entry(0:rdm%sign_length)
        integer :: ij, kl, i, j, k, l, ij_proc_row, sign_length_old
        integer :: iproc, elem_proc, pops_unit, file_end, hash_val, ierr
        character(len=*), parameter :: t_r = 'read_rdm_popsfile'

        ! Make sure that the RDM is empty first.
        rdm%nelements = 0
        rdm%elements = 0_int_rdm
        call clear_hash_table(rdm%hash_table)

        do iproc = 0, nProcessors-1

            if (iproc == iProcIndex) then
                pops_unit = get_free_unit()
                open(pops_unit, file='RDM_POPSFILE', status='old', form='unformatted')

                ! Read in the first line, which holds sign_length for the
                ! printed RDM - check it is consistent.
                read(pops_unit, iostat=file_end) sign_length_old
                if (sign_length_old /= rdm%sign_length) then
                    call stop_all(t_r, "Error reading RDM_POPSFILE - the number of RDMs printed in the &
                                       &popsfile is different to the number to be sampled.")
                else if (file_end > 0) then
                    call stop_all(t_r, "Error reading first line of RDM_POPSFILE.")
                end if

                do
                    read(pops_unit, iostat=file_end) rdm_entry
                    if (file_end > 0) call stop_all(t_r, "Error reading RDM_POPSFILE.")
                    ! file_end < 0 => end of file reached.
                    if (file_end < 0) exit

                    call calc_separate_rdm_labels(rdm_entry(0), ij, kl, i, j, k, l)
                    ij_proc_row = nbasis*(i-1) - i*(i-1)/2 + j - i
                    ! Calculate the process for the element.
                    elem_proc = (ij_proc_row-1)*nProcessors/spawn%nrows

                    if (elem_proc == iProcIndex) then
                        ! Add the element to the RDM list.
                        rdm%nelements = rdm%nelements + 1
                        rdm%elements(:, rdm%nelements) = rdm_entry
                        ! Add in the entry to the hash table.
                        hash_val = FindWalkerHash((/i,j,k,l/), size(rdm%hash_table))
                        call add_hash_table_entry(rdm%hash_table, rdm%nelements, hash_val)
                    else if (elem_proc > iProcIndex) then
                        ! If we've reached the end of the section of the
                        ! popsfile for this processor, then finish reading.
                        exit
                    end if
                end do

                close(pops_unit)
            end if

            ! Wait for the current processor to finish reading its RDM elements.
            call MPIBarrier(ierr)
        end do

    end subroutine read_rdm_popsfile

    subroutine print_rdms_with_spin(rdm, rdm_trace)

        ! Print the RDM stored in rdm to files, normalised by rdm_trace.

        ! This routine will print out *all* the spin cobminations separately,
        ! including both aaaa and bbbb arrays, and all other combinations.

        ! The files are called 'TwoRDM_aaaa', 'TwoRDM_abab', 'TwoRDM_abba', etc...

        use Parallel_neci, only: MPIBarrier
        use sort_mod, only: sort
        use UMatCache, only: spatial
        use util_mod, only: get_free_unit

        type(rdm_list_t), intent(inout) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)

        integer(int_rdm) :: pqrs
        integer :: i, irdm, ierr, iproc, write_unit
        integer :: iunit_aaaa, iunit_abab, iunit_abba
        integer :: iunit_bbbb, iunit_baba, iunit_baab
        integer :: pq, rs, p, q, r, s
        integer :: p_spat, q_spat, r_spat, s_spat
        real(dp) :: rdm_sign(rdm%sign_length)
        character(3) :: sgn_len, suffix

        ! Store rdm%sign_length as a string, for the formatting string.
        write(sgn_len,'(i3)') rdm%sign_length

        call sort(rdm%elements(:,1:rdm%nelements))

        do iproc = 0, nProcessors-1
            do irdm = 1, rdm%sign_length
                write(suffix, '('//int_fmt(irdm,0)//')') irdm

                if (iproc == iProcIndex) then

                    ! Open all the files to be written to:
                    ! Let the first processor clear all the files to start with.
                    if (iproc == 0) then
                        iunit_aaaa = get_free_unit()
                        open(iunit_aaaa, file='TwoRDM_aaaa.'//trim(suffix), status='replace')
                        iunit_abab = get_free_unit()
                        open(iunit_abab, file='TwoRDM_abab.'//trim(suffix), status='replace')
                        iunit_abba = get_free_unit()
                        open(iunit_abba, file='TwoRDM_abba.'//trim(suffix), status='replace')
                        iunit_bbbb = get_free_unit()
                        open(iunit_bbbb, file='TwoRDM_bbbb.'//trim(suffix), status='replace')
                        iunit_baba = get_free_unit()
                        open(iunit_baba, file='TwoRDM_baba.'//trim(suffix), status='replace')
                        iunit_baab = get_free_unit()
                        open(iunit_baab, file='TwoRDM_baab.'//trim(suffix), status='replace')
                    else
                        iunit_aaaa = get_free_unit()
                        open(iunit_aaaa, file='TwoRDM_aaaa.'//trim(suffix), status='old', position='append')
                        iunit_abab = get_free_unit()
                        open(iunit_abab, file='TwoRDM_abab.'//trim(suffix), status='old', position='append')
                        iunit_abba = get_free_unit()
                        open(iunit_abba, file='TwoRDM_abba.'//trim(suffix), status='old', position='append')
                        iunit_bbbb = get_free_unit()
                        open(iunit_bbbb, file='TwoRDM_bbbb.'//trim(suffix), status='old', position='append')
                        iunit_baba = get_free_unit()
                        open(iunit_baba, file='TwoRDM_baba.'//trim(suffix), status='old', position='append')
                        iunit_baab = get_free_unit()
                        open(iunit_baab, file='TwoRDM_baab.'//trim(suffix), status='old', position='append')
                    end if

                    do i = 1, rdm%nelements
                        pqrs = rdm%elements(0,i)
                        ! Obtain spin orbital labels.
                        call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
                        call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

                        ! Normalise.
                        rdm_sign = rdm_sign/rdm_trace

                        p_spat = spatial(p); q_spat = spatial(q);
                        r_spat = spatial(r); s_spat = spatial(s);

                        ! Find out what the spin labels are, and print the RDM
                        ! element to the appropriate file.
                        if (is_alpha(p) .and. is_alpha(q) .and. is_alpha(r) .and. is_alpha(s)) then
                            write_unit = iunit_aaaa
                        else if (is_alpha(p) .and. is_beta(q) .and. is_alpha(r) .and. is_beta(s)) then
                            write_unit = iunit_abab
                        else if (is_alpha(p) .and. is_beta(q) .and. is_beta(r) .and. is_alpha(s)) then
                            write_unit = iunit_abba
                        else if (is_beta(p) .and. is_beta(q) .and. is_beta(r) .and. is_beta(s)) then
                            write_unit = iunit_bbbb
                        else if (is_beta(p) .and. is_alpha(q) .and. is_beta(r) .and. is_alpha(s)) then
                            write_unit = iunit_baba
                        else if (is_beta(p) .and. is_alpha(q) .and. is_alpha(r) .and. is_beta(s)) then
                            write_unit = iunit_baab
                        end if

                        if (abs(rdm_sign(irdm)) > 1.e-12_dp) then
                            write(write_unit,'(4i6,'//trim(sgn_len)//'g25.17)') p_spat, q_spat, r_spat, s_spat, rdm_sign(irdm)
                        end if
                    end do

                    close(iunit_aaaa); close(iunit_abab); close(iunit_abba); close(iunit_abab);
                    close(iunit_bbbb); close(iunit_baba); close(iunit_baab); close(iunit_abba);
                end if
            end do

            ! Wait for the current processor to finish printing its RDM elements.
            call MPIBarrier(ierr)
        end do

    end subroutine print_rdms_with_spin

    subroutine print_spinfree_2rdm(rdm, rdm_trace)

        ! Print all the RDM elements stored in rdm to a single file (for each
        ! state being sampled).

        ! The stem of the filenames is "spinfree_TwoRDM", and a final line
        ! required by MPQC to read spinfree 2-RDMs is also printed. This
        ! routine also assumes that the RDM element labels are already in
        ! spatial form, performing no transformation from spin to spatial
        ! form. This routine is therefore appropriate for printing spinfree
        ! 2-RDMs.

        use Parallel_neci, only: MPIBarrier
        use ParallelHelper, only: root
        use sort_mod, only: sort
        use util_mod, only: get_free_unit

        type(rdm_list_t), intent(inout) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)

        integer(int_rdm) :: pqrs
        integer :: i, irdm, iunit, iproc, ierr
        integer :: pq_spat, rs_spat
        integer :: p_spat, q_spat, r_spat, s_spat
        real(dp) :: rdm_sign(rdm%sign_length)
        character(30) :: rdm_filename

        call sort(rdm%elements(:,1:rdm%nelements))

        do iproc = 0, nProcessors-1
            if (iproc == iProcIndex) then

                ! Loop over all RDMs beings sampled.
                do irdm = 1, rdm%sign_length
                    write(rdm_filename, '("spinfree_TwoRDM.",'//int_fmt(irdm,0)//')') irdm
                    ! Open the file to be written to.
                    iunit = get_free_unit()
                    ! Let the first process clear the file, if it already exist.
                    if (iproc == 0) then
                        open(iunit, file=rdm_filename, status='replace')
                    else
                        open(iunit, file=rdm_filename, status='old', position='append')
                    end if

                    do i = 1, rdm%nelements
                        pqrs = rdm%elements(0,i)
                        ! Obtain spin orbital labels.
                        call calc_separate_rdm_labels(pqrs, pq_spat, rs_spat, r_spat, s_spat, q_spat, p_spat)
                        call extract_sign_rdm(rdm%elements(:,i), rdm_sign)
                        ! Normalise.
                        rdm_sign = rdm_sign/rdm_trace

                        if (abs(rdm_sign(irdm)) > 1.e-12_dp) then
                            write(iunit,"(4I15, F30.20)") p_spat, q_spat, r_spat, s_spat, rdm_sign(irdm)
                        end if
                    end do

                    ! The following final line is required by (I assume!) MPQC.
                    ! Let the last process print it.
                    if (iProcIndex == nProcessors-1) then
                        write(iunit, "(4I15, F30.20)") -1, -1, -1, -1, -1.0_dp
                    end if

                    close(iunit)
                end do
            end if

            ! Wait for the current processor to finish printing its RDM elements.
            call MPIBarrier(ierr)
        end do

    end subroutine print_spinfree_2rdm


    ! ------- Routines for finalising 1-RDMs ---------------------------------

    subroutine Finalise_1e_RDM(matrix, matrix_diag, irdm, norm_1rdm, tOldRDMs)

        ! This routine takes the 1-RDM (matrix), normalises it, makes it
        ! hermitian if required, and prints out the versions we're interested
        ! in. This is only ever called at the very end of a calculation.

        use LoggingData, only: twrite_RDMs_to_read, twrite_normalised_RDMs, tForceCauchySchwarz
        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: iProcIndex, MPISumAll
        use RotateOrbsData, only: NoOrbs

        real(dp), intent(inout) :: matrix(:,:)
        real(dp), intent(inout) :: matrix_diag(:)
        integer, intent(in) :: irdm
        real(dp), intent(out) :: norm_1rdm
        logical, intent(in) :: tOldRDMs

        integer :: ierr
        real(dp) :: trace_1rdm, SumN_Rho_ii
        real(dp), allocatable :: AllNode_one_rdm(:,:)

        norm_1rdm = 0.0_dp

        if (RDMExcitLevel == 1) then
            allocate(AllNode_one_rdm(NoOrbs, NoOrbs), stat=ierr)

            call MPISumAll(matrix, AllNode_one_rdm)
            matrix = AllNode_one_rdm

            deallocate(AllNode_one_rdm)
        end if

        if (iProcIndex == 0) then
            ! Find the normalisation.
            call calc_1e_norms(matrix, matrix_diag, trace_1rdm, norm_1rdm, SumN_Rho_ii)

            ! Write out the unnormalised, non-hermitian OneRDM_POPS.
            if (twrite_RDMs_to_read) call write_1rdm(matrix, irdm, norm_1rdm, .false., tOldRDMs)

            ! Enforce the hermiticity condition.  If the RDMExcitLevel is not 1, the
            ! 1-RDM has been constructed from the hermitian 2-RDM, so this will not
            ! be necessary.
            ! The HF_Ref and HF_S_D_Ref cases are not hermitian by definition.
            if (RDMExcitLevel == 1) then
                call make_1e_rdm_hermitian(matrix, norm_1rdm)

                if (tForceCauchySchwarz) then
                    call Force_Cauchy_Schwarz(matrix)
                end if
            end if

            ! Write out the final, normalised, hermitian OneRDM.
            if (tWrite_normalised_RDMs) call write_1rdm(matrix, irdm, norm_1rdm, .true., tOldRDMs)

            write(6,'(1X,"SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS:",1X,F20.13)') SumN_Rho_ii
        end if

    end subroutine Finalise_1e_RDM

    subroutine Force_Cauchy_Schwarz(matrix)

        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: nbasis

        real(dp), intent(inout) :: matrix(:,:)

        integer :: i, j
        real(dp) :: UpperBound

        write(6,'("Ensuring that Cauchy--Schwarz inequality holds.")')

        associate(ind => SymLabelListInv_rot)
            do i = 1, nbasis
                do j = 1, nbasis
                    UpperBound = sqrt(matrix(ind(i),ind(i)) * matrix(ind(j),ind(j)))

                    if (abs(matrix(ind(i), ind(j))) > UpperBound) then

                        if (matrix(ind(i), ind(j)) < 0.0_dp) then
                            matrix(ind(i), ind(j)) = -UpperBound
                        else if (matrix(ind(i), ind(j)) > 0.0_dp) then
                            matrix(ind(i), ind(j)) = UpperBound
                        end if

                        write(6,'("Changing element:")') i, j
                    else
                        cycle
                    end if
                end do
            end do
        end associate

    end subroutine Force_Cauchy_Schwarz

    subroutine calc_1e_norms(matrix, matrix_diag, trace_1rdm, norm_1rdm, SumN_Rho_ii)

        ! We want to 'normalise' the reduced density matrices. These are not
        ! even close to being normalised at the moment, because of the way
        ! they are calculated on the fly. They should be calculated from a
        ! normalised wavefunction. But we know that the trace of the one
        ! electron reduced density matrix must be equal to the number of the
        ! electrons. We can use this to find the factor we must divide the
        ! 1-RDM through by.

        use FciMCData, only: HFDet_True
        use LoggingData, only: tDiagRDM
        use rdm_data, only: tOpenShell
        use RotateOrbsData, only: SymLabelListInv_rot, NoOrbs
        use SystemData, only: BRR, nel
        use UMatCache, only: gtID

        real(dp), intent(in) :: matrix(:,:)
        real(dp), intent(inout) :: matrix_diag(:)
        real(dp), intent(out) :: trace_1rdm, norm_1rdm, SumN_Rho_ii

        integer :: i, HFDet_ID, BRR_ID

        trace_1rdm = 0.0_dp
        norm_1rdm = 0.0_dp

        do i = 1, NoOrbs
            trace_1rdm = trace_1rdm + matrix(i,i)
        end do

        norm_1rdm = real(nel, dp) / trace_1rdm

        ! Need to multiply each element of the 1 electron reduced density matrices
        ! by nel / trace_1rdm,
        ! and then add it's contribution to the energy.

        ! Want to sum the diagonal elements of the 1-RDM for the HF orbitals.
        ! Given the HF orbitals, SymLabelListInv_rot tells us their position
        ! in the 1-RDM.
        SumN_Rho_ii = 0.0_dp

        do i = 1, NoOrbs
            ! Rho_ii is the diagonal elements of the 1-RDM. We want this
            ! ordered according to the energy of the orbitals. Brr has the
            ! orbital numbers in order of energy... i.e Brr(2) = the orbital
            ! index with the second lowest energy. Brr is always in spin
            ! orbitals. i gives the energy level, BRR gives the orbital,
            ! SymLabelListInv_rot gives the position of  this orbital in
            ! one_rdm.

            associate(ind => SymLabelListInv_rot)
                if (tDiagRDM) then
                    if (tOpenShell) then
                        matrix_diag(i) = matrix(ind(BRR(i)), ind(BRR(i))) * norm_1rdm
                    else
                        BRR_ID = gtID(BRR(2*i))
                        matrix_diag(i) = matrix(ind(BRR_ID), ind(BRR_ID)) * norm_1rdm
                    end if
                end if

                if (i <= nel) then
                    if (tOpenShell) then
                        SumN_Rho_ii = SumN_Rho_ii + ( matrix(ind(HFDet_True(i)), ind(HFDet_True(i))) * norm_1rdm )
                    else
                        HFDet_ID = gtID(HFDet_True(i))
                        SumN_Rho_ii = SumN_Rho_ii + ( matrix(ind(HFDet_ID), ind(HFDet_ID)) * norm_1rdm ) / 2.0_dp
                    end if
                end if
            end associate
        end do

    end subroutine calc_1e_norms

    subroutine make_1e_rdm_hermitian(matrix, norm_1rdm)

        ! Simply average the 1-RDM(i,j) and 1-RDM(j,i) elements which should
        ! be equal in a perfect world.

        use RotateOrbsData, only: SymLabelListInv_rot, NoOrbs

        real(dp), intent(inout) :: matrix(:,:)
        real(dp), intent(in) :: norm_1rdm

        real(dp) :: max_error_herm, sum_error_herm
        integer :: i, j
        real(dp) :: temp

        max_error_herm = 0.0_dp
        sum_error_herm = 0.0_dp

        associate(ind => SymLabelListInv_rot)
            do i = 1, NoOrbs
                do j = i, NoOrbs
                    if ((abs((matrix(ind(i),ind(j))*norm_1rdm) - (matrix(ind(j),ind(i))*norm_1rdm))) > max_error_herm) then
                        max_error_herm = abs(matrix(ind(i),ind(j))*norm_1rdm - matrix(ind(j), ind(i))*norm_1rdm)
                    end if

                    sum_error_herm = sum_error_herm + abs(matrix(ind(i),ind(j))*norm_1rdm - matrix(ind(j),ind(i))*norm_1rdm)

                    temp = (matrix(ind(i),ind(j)) + matrix(ind(j),ind(i)))/2.0_dp
                    matrix(ind(i), ind(j)) = temp
                    matrix(ind(j), ind(i)) = temp
                end do
            end do
        end associate

        ! Output the hermiticity errors.
        write(6,'(1X,"MAX ABS ERROR IN 1RDM HERMITICITY",F20.13)') max_error_herm
        write(6,'(1X,"MAX ABS ERROR IN 1RDM HERMITICITY",F20.13)') sum_error_herm

    end subroutine make_1e_rdm_hermitian

    subroutine write_1rdm(one_rdm, irdm, norm_1rdm, tNormalise, tOldRDMs)

        ! This routine writes out the OneRDM. If tNormalise is true, we are
        ! printing the normalised, hermitian matrix. Otherwise, norm_1rdm is
        ! ignored and we print both 1-RDM(i,j) and 1-RDM(j,i) (in binary)
        ! for the OneRDM_POPS file to be read in in a restart calculation.

        use rdm_data, only: tOpenShell
        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: nbasis
        use UMatCache, only: gtID
        use util_mod, only: get_free_unit, int_fmt

        real(dp), intent(in) :: one_rdm(:,:)
        integer, intent(in) :: irdm
        real(dp), intent(in) :: norm_1rdm
        logical, intent(in) :: tNormalise
        logical, intent(in) :: tOldRDMs

        integer :: i, j, iSpat, jSpat
        integer :: one_rdm_unit
        character(20) :: filename

        if (tNormalise) then
            ! Haven't got the capabilities to produce multiple 1-RDMs yet.
            write(6,'(1X,"Writing out the *normalised* 1 electron density matrix to file")')
            call neci_flush(6)
            one_rdm_unit = get_free_unit()
#ifdef _MOLCAS_
            call molcas_open(one_rdm_unit, "ONERDM")
#else
            if (tOldRDMs) then
                write(filename, '("OneRDM_old.",'//int_fmt(irdm,0)//')') irdm
                open(one_rdm_unit, file=trim(filename), status='unknown')
            else
                write(filename, '("OneRDM.",'//int_fmt(irdm,0)//')') irdm
                open(one_rdm_unit, file=trim(filename), status='unknown')
            end if
#endif
        else
            ! Only every write out 1 of these at the moment.
            write(6,'(1X,"Writing out the *unnormalised* 1 electron density matrix to file for reading in")')
            call neci_flush(6)
            one_rdm_unit = get_free_unit()
            if (tOldRDMs) then
                write(filename, '("OneRDM_POPS_old.",'//int_fmt(irdm,0)//')') irdm
                open(one_rdm_unit, file=trim(filename), status='unknown', form='unformatted')
            else
                write(filename, '("OneRDM_POPS.",'//int_fmt(irdm,0)//')') irdm
                open(one_rdm_unit, file=trim(filename), status='unknown', form='unformatted')
            end if
        end if

        ! Currently always printing 1-RDM in spin orbitals.
        associate(ind => SymLabelListInv_rot)
            do i = 1, nbasis
                do j = 1, nbasis
                    if (tOpenShell) then
                        if (abs(one_rdm(ind(i), ind(j))) > 1.0e-12_dp) then
                            if (tNormalise .and. (i <= j)) then
                                write(one_rdm_unit,"(2I6,G25.17)") i, j, one_rdm(ind(i), ind(j)) * norm_1rdm
                            else if (.not. tNormalise) then
                                ! For the pops, we haven't made the 1-RDM hermitian yet,
                                ! so print both the 1-RDM(i,j) and 1-RDM(j,i) elements.
                                ! This is written in binary.
                                write(one_rdm_unit) i, j, one_rdm(ind(i), ind(j))
                            end if
                        end if
                    else
                        iSpat = gtID(i)
                        jSpat = gtID(j)
                        if (abs(one_rdm(ind(iSpat), ind(jSpat))) > 1.0e-12_dp) then
                            if (tNormalise .and. (i <= j)) then
                                if ((mod(i,2) == 0 .and. mod(j,2) == 0) .or. &
                                    (mod(i,2) /= 0 .and. mod(j,2) .ne. 0)) then
                                    write(one_rdm_unit,"(2I6,G25.17)") i, j, &
                                        ( one_rdm(ind(iSpat),ind(jSpat)) * norm_1rdm ) / 2.0_dp
                                end if
                            else if (.not. tNormalise) then
                                ! The popsfile can be printed in spatial orbitals.
                                if (mod(i,2) == 0 .and. mod(j,2) == 0) then
                                    write(one_rdm_unit) iSpat, jSpat, one_rdm(ind(iSpat), ind(jSpat))
                                end if
                            end if
                        end if
                    end if
                end do
            end do
        end associate

        close(one_rdm_unit)

    end subroutine write_1rdm

end module rdm_finalising
