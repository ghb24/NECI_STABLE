#include "macros.h"

module rdm_parallel

    ! This module contains all routines used for the calculations of reduced
    ! density matrices, when distributed across all MPI processes.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use rdm_data, only: nrdms, rdm_list_t, rdm_spawn_t, one_rdm_t
    use util_mod

    implicit none

contains

    subroutine init_rdm_list_t(rdm, sign_length, max_nelements, nhashes)

        ! Initialise an rdm_list_t object.

        ! Out: rdm - rdm_list_t object to be initialised.
        ! In:  sign_length - the number of signs which can be stored for each element.
        ! In:  max_nelements - the length of the rdm%elements array.
        ! In:  nhashes - the number of unique hashes for indexing the hash table.

        use hash, only: init_hash_table

        type(rdm_list_t), intent(out) :: rdm
        integer, intent(in) :: sign_length, max_nelements, nhashes

        integer :: ierr

        rdm%sign_length = sign_length
        rdm%max_nelements = max_nelements
        rdm%nhashes = nhashes

        allocate(rdm%elements(0:sign_length, max_nelements))
        rdm%elements = 0_int_rdm

        allocate(rdm%hash_table(nhashes), stat=ierr)
        call init_hash_table(rdm%hash_table)

    end subroutine init_rdm_list_t

    subroutine init_rdm_spawn_t(spawn, nrows, sign_length, max_nelements_send, max_nelements_recv, nhashes)

        ! Initialise an rdm_spawn_t object.

        ! * IMPORTANT * - This routine will not allocate the hash table for the
        ! received RDM, spawn%rdm_recv%hash_table, since it usually won't be
        ! used. If it is needed for a particular application then it must be
        ! allocated and initialised separately.

        ! Out: spawn - rdm_spawn_t object to be initialised.
        ! In:  nrows - the number of rows in the RDM.
        ! In:  sign_length - the number of signs which can be stored for each element.
        ! In:  max_nelements_send - the length of the spawn%rdm_send%elements array.
        ! In:  max_nelements_recv - the length of the spawn%rdm_recv%elements array.
        ! In:  nhashes - the number of unique hashes for indexing spawn%rdm_send%hash_table.

        type(rdm_spawn_t), intent(out) :: spawn
        integer, intent(in) :: nrows, sign_length, nhashes
        integer, intent(in) :: max_nelements_send, max_nelements_recv
        integer :: i, ierr
        real(dp) :: slots_per_proc

        spawn%nrows = nrows

        call init_rdm_list_t(spawn%rdm_send, sign_length, max_nelements_send, nhashes)
        ! Don't need the hash table for the received list, so pass 0 for nhashes.
        call init_rdm_list_t(spawn%rdm_recv, sign_length, max_nelements_recv, 0)

        allocate(spawn%free_slots(0:nProcessors-1), stat=ierr)
        allocate(spawn%init_free_slots(0:nProcessors-1), stat=ierr)

        ! Equally divide RDM rows across all processors.
        slots_per_proc = real(max_nelements_send, dp)/real(nProcessors, dp)
        do i = 0, nProcessors-1
            spawn%init_free_slots(i) = nint(slots_per_proc*i)+1
        end do

        ! Set the free slots array to its initial value.
        spawn%free_slots = spawn%init_free_slots

    end subroutine init_rdm_spawn_t

    pure subroutine calc_combined_rdm_label(p, q, r, s, pqrs)

        ! Combine the four 2-RDM spin orbital labels into unique integers.
        ! p and q are combined into one number, pq. r and s are combined into
        ! one number, rs. Both of these are then combined into one single
        ! number, pqrs. The largest value for pqrs is M^4, where M is the
        ! number of spin orbitals.
        
        ! The compression defined in this routine will not give a fully
        ! compressed RDM index labelling, because it allows a separate pq
        ! integer if p and q are equal, which will always give a zero
        ! RDM element, and the same for r and s. It also doesn't take
        ! spatial symmetry into account. But this is fine if one just
        ! seeks a unique combined label for each combination of individual
        ! spin orbital labels.

        ! In: p, q, r, s - spin orbitals of the RDM contribution.
        ! Out: pqrs - Label combining p, q, r and s.

        use SystemData, only: nbasis

        integer, intent(in) :: p, q, r, s
        integer(int_rdm), intent(out) :: pqrs
        integer :: pq, rs

        pq = (p-1)*nbasis + q
        rs = (r-1)*nbasis + s
        pqrs = (pq-1)*(nbasis**2) + rs

    end subroutine calc_combined_rdm_label

    pure subroutine calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)

        use SystemData, only: nbasis

        integer(int_rdm), intent(in) :: pqrs
        integer, intent(out) :: pq, rs, p, q, r, s

        rs = mod(pqrs-1, nbasis**2) + 1
        pq = (pqrs - rs)/(nbasis**2) + 1

        q = mod(pq-1, nbasis) + 1
        p = (pq - q)/nbasis + 1

        s = mod(rs-1, nbasis) + 1
        r = (rs - s)/nbasis + 1

    end subroutine calc_separate_rdm_labels

    pure subroutine extract_sign_rdm(rdm_entry, real_sign)

        integer(int_rdm), intent(in) :: rdm_entry(0:)
        real(dp), intent(out) :: real_sign(size(rdm_entry)-1)
        integer(int_rdm) :: int_sign(size(rdm_entry)-1)

        int_sign = rdm_entry(1:size(int_sign))
        real_sign = transfer(int_sign, real_sign)

    end subroutine extract_sign_rdm

    pure subroutine encode_sign_rdm(rdm_entry, real_sign)

        integer(int_rdm), intent(inout) :: rdm_entry(0:)
        real(dp), intent(in) :: real_sign(1:size(rdm_entry)-1)
        integer(int_rdm) :: int_sign(1:size(rdm_entry)-1)

        int_sign = transfer(real_sign, int_sign)
        rdm_entry(1:size(int_sign)) = int_sign

    end subroutine encode_sign_rdm

    subroutine add_to_rdm_spawn_t(spawn, p, q, r, s, contrib_sign, spinfree)

        ! In/Out: rdm_spawn - the rdm_spawn_t object to which contributions will be added.
        ! In: p, q, r, s - spin orbitals of the RDM contribution, with p<q, r<s.
        ! In: contrib_sign - the sign (amplitude) of the contribution to be added.

        use hash, only: hash_table_lookup, add_hash_table_entry
        use SystemData, only: nbasis

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: p, q, r, s
        real(dp), intent(in) :: contrib_sign(spawn%rdm_send%sign_length)
        logical, intent(in) :: spinfree

        integer :: pq_compressed, proc, ind, hash_val
        integer(int_rdm) :: pqrs
        real(dp) :: real_sign_old(spawn%rdm_send%sign_length), real_sign_new(spawn%rdm_send%sign_length)
        logical :: tSuccess, list_full
        character(*), parameter :: this_routine = 'add_to_rdm_spawn_t'

        associate(rdm => spawn%rdm_send)

            ! Calculate combined RDM labels. A different ordering is used for
            ! outputting RDMs with and without spin. The following definitions
            ! will aid ordering via a sort operation later.
            if (.not. spinfree) then
                call calc_combined_rdm_label(p, q, r, s, pqrs)
            else
                call calc_combined_rdm_label(r, s, q, p, pqrs)
            end if

            ! Search to see if this RDM element is already in the RDM array.
            ! If it, tSuccess will be true and ind will hold the position of the
            ! entry in rdm%elements.
            call hash_table_lookup((/p,q,r,s/), (/pqrs/), 0, rdm%hash_table, rdm%elements, ind, hash_val, tSuccess)

            if (tSuccess) then
                ! Extract the existing sign.
                call extract_sign_rdm(rdm%elements(:,ind), real_sign_old)
                ! Update the total sign.
                real_sign_new = real_sign_old + contrib_sign
                ! Encode the new sign.
                call encode_sign_rdm(rdm%elements(:,ind), real_sign_new)
            else
                ! Determine the process label.
                if (.not. spinfree) then
                    ! The following maps (p,q), with p<q, to single integers
                    ! with no gaps. It is benefical to have no gaps here, for
                    ! good load balancing. The final integers are ordered so
                    ! that p is dominant over q.
                    pq_compressed = nbasis*(p-1) - p*(p-1)/2 + q - p
                    ! Calculate the process for the element.
                    proc = (pq_compressed-1)*nProcessors/spawn%nrows
                else
                    ! For spin-free case, we halve the number of labels. Also,
                    ! the last two labels are dominant in the ordering, so use
                    ! these instead, to allow writing out in the correct order.
                    pq_compressed = (nbasis/2)*(r-1) + s
                    proc = (pq_compressed-1)*nProcessors/(nbasis**2/4)
                end if

                ! Check that there is enough memory for the new spawned RDM entry.
                list_full = .false.
                if (proc == nProcessors - 1) then
                    if (spawn%free_slots(proc) > rdm%max_nelements) list_full = .true.
                else
                    if (spawn%free_slots(proc) > spawn%init_free_slots(proc+1)) list_full = .true.
                end if
                if (list_full) then
                    write(6,'("Attempting to add an RDM contribution to the spawned list on processor:",&
                               &1X,'//int_fmt(proc,0)//')') proc
                    write(6,'("No memory slots available for this spawn.")')
                    call stop_all(this_routine, "Out of memory for spawned RDM contributions.")
                end if

                rdm%elements(0, spawn%free_slots(proc)) = pqrs
                call encode_sign_rdm(rdm%elements(:, spawn%free_slots(proc)), contrib_sign)

                call add_hash_table_entry(rdm%hash_table, spawn%free_slots(proc), hash_val)

                spawn%free_slots(proc) = spawn%free_slots(proc) + 1
            end if

            !if (p > q .or. r > s) then
            !    write(6,'("p, q, r, s:", 1X,'//int_fmt(p,0)//', 1X,'//int_fmt(q,0)//', &
            !               &1X,'//int_fmt(r,0)//', 1X,'//int_fmt(s,0)//')') p, q, r, s
            !    call stop_all(this_routine,"Incorrect ordering of RDM orbitals passed to RDM spawning routine.")
            !end if

        end associate

    end subroutine add_to_rdm_spawn_t

    subroutine communicate_rdm_spawn_t(spawn)

        use hash, only: clear_hash_table
        use Parallel_neci, only: MPIAlltoAll, MPIAlltoAllv

        type(rdm_spawn_t), intent(inout) :: spawn

        integer :: i, ierr
        integer(MPIArg) :: send_sizes(0:nProcessors-1), recv_sizes(0:nProcessors-1)
        integer(MPIArg) :: send_displs(0:nProcessors-1), recv_displs(0:nProcessors-1)
        character(*), parameter :: this_routine = 'communicate_rdm_spawn_t'

        ! How many rows of data to send to each processor.
        do i = 0, nProcessors-1
            send_sizes(i) = int(spawn%free_slots(i) - spawn%init_free_slots(i), MPIArg)
        end do

        ! The displacement of the beginning of each processor's section of the
        ! free_slots array, relative to the first element of this array.
        send_displs = int(spawn%init_free_slots-1, MPIArg)

        call MPIAlltoAll(send_sizes, 1, recv_sizes, 1, ierr)

        recv_displs(0) = 0
        do i = 1, nProcessors-1
            recv_displs(i) = recv_displs(i-1) + recv_sizes(i-1)
        end do

        ! The total number of RDM elements sent to this processor.
        spawn%rdm_recv%nelements = int(sum(recv_sizes), sizeof_int)

        if (spawn%rdm_recv%nelements > spawn%rdm_recv%max_nelements) then
            write(6,'("Attempting to receive RDM elements on processor:",&
                       &1X,'//int_fmt(iProcIndex,0)//')') iProcIndex
            write(6,'("Insufficient space in the receiving RDM array.")')
            call stop_all(this_routine, "Not enough memory to communicate RDM elements.")
        end if

        send_sizes = send_sizes*size(spawn%rdm_send%elements,1)
        recv_sizes = recv_sizes*size(spawn%rdm_send%elements,1)
        send_displs = send_displs*size(spawn%rdm_send%elements,1)
        recv_displs = recv_displs*size(spawn%rdm_send%elements,1)

        ! Perform the communication.
        call MPIAlltoAllv(spawn%rdm_send%elements, send_sizes, send_displs, &
                          spawn%rdm_recv%elements, recv_sizes, recv_displs, ierr)

        ! Now we can reset the free_slots array and reset the hash table.
        spawn%free_slots = spawn%init_free_slots
        call clear_hash_table(spawn%rdm_send%hash_table)

    end subroutine communicate_rdm_spawn_t

    subroutine add_rdm_1_to_rdm_2(rdm_1, rdm_2)

        use hash, only: hash_table_lookup, add_hash_table_entry

        type(rdm_list_t), intent(in) :: rdm_1
        type(rdm_list_t), intent(inout) :: rdm_2

        integer :: i, pq, rs, p, q, r, s
        integer(int_rdm) :: pqrs
        integer :: ind, hash_val
        real(dp) :: real_sign_old(rdm_2%sign_length), real_sign_new(rdm_2%sign_length)
        real(dp) :: spawn_sign(rdm_2%sign_length)
        logical :: tSuccess
        character(*), parameter :: this_routine = 'add_rdm_1_to_rdm_2'

        do i = 1, rdm_1%nelements
            ! Decode the compressed RDM labels.
            pqrs = rdm_1%elements(0,i)
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)

            ! Extract the spawned sign.
            call extract_sign_rdm(rdm_1%elements(:,i), spawn_sign)

            ! Search to see if this RDM element is already in the RDM 2.
            ! If it, tSuccess will be true and ind will hold the position of the
            ! element in rdm.
            call hash_table_lookup((/p,q,r,s/), (/pqrs/), 0, rdm_2%hash_table, rdm_2%elements, ind, hash_val, tSuccess)

            if (tSuccess) then
                ! Extract the existing sign.
                call extract_sign_rdm(rdm_2%elements(:,ind), real_sign_old)

                ! Update the total sign.
                real_sign_new = real_sign_old + spawn_sign
                ! Encode the new sign.
                call encode_sign_rdm(rdm_2%elements(:,ind), real_sign_new)
            else
                ! Check that there is enough memory for the new RDM element.
                if (rdm_2%nelements+1 > rdm_2%max_nelements) then
                    write(6,'("Ran out of memory while adding new elements to the RDM array.")')
                    call stop_all(this_routine, "Out of memory for RDM elements.")
                end if

                ! Update the rdm array, and its hash table, and the number of
                ! RDM elements.
                rdm_2%nelements = rdm_2%nelements + 1
                rdm_2%elements(0, rdm_2%nelements) = pqrs
                call encode_sign_rdm(rdm_2%elements(:, rdm_2%nelements), spawn_sign)
                call add_hash_table_entry(rdm_2%hash_table, rdm_2%nelements, hash_val)
            end if
        end do

    end subroutine add_rdm_1_to_rdm_2

    subroutine calc_rdm_energy(rdm, rdm_energy)

        use rdm_data, only: rdm_list_t
        use rdm_integral_fns, only: one_elec_int, two_elec_int
        use SystemData, only: nel

        type(rdm_list_t), intent(inout) :: rdm
        real(dp), intent(out) :: rdm_energy(rdm%sign_length)

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_energy = 0.0_dp

        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

            ! Decode pqrs label into p, q, r and s labels.
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)

            ! The 2-RDM contribution to the energy:
            rdm_energy = rdm_energy + rdm_sign*two_elec_int(p,q,r,s)
            ! The 1-RDM contribution to the energy:
            if (p == r) rdm_energy = rdm_energy + rdm_sign*one_elec_int(q,s)/(nel-1)
            if (q == s) rdm_energy = rdm_energy + rdm_sign*one_elec_int(p,r)/(nel-1)
            if (p == s) rdm_energy = rdm_energy - rdm_sign*one_elec_int(q,r)/(nel-1)
            if (q == r) rdm_energy = rdm_energy - rdm_sign*one_elec_int(p,s)/(nel-1)
        end do

    end subroutine calc_rdm_energy

    subroutine calc_rdm_trace(rdm, rdm_trace)

        use rdm_data, only: rdm_spawn_t
        use rdm_integral_fns, only: one_elec_int, two_elec_int
        use SystemData, only: nel

        type(rdm_list_t), intent(inout) :: rdm
        real(dp), intent(out) :: rdm_trace(rdm%sign_length)

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_trace = 0.0_dp

        ! Just working with one processor and the spawned list for now,
        ! as a check - this will be corrected.
        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)

            if (pq == rs) then
                call extract_sign_rdm(rdm%elements(:,i), rdm_sign)
                rdm_trace = rdm_trace + rdm_sign
            end if
        end do

        ! When dividing the RDM by the output trace, we want the new
        ! normalisation to be (nel*(nel-1))/2.
        rdm_trace = rdm_trace * 2.0_dp/(nel*(nel-1))

    end subroutine calc_rdm_trace

    subroutine calc_rdm_spin(rdm, rdm_trace, rdm_spin)

        ! Return the (unnormalised) estimate of <S^2> from the instantaneous
        ! 2RDM estimates.

        use rdm_data, only: rdm_spawn_t
        use SystemData, only: nel
        use UMatCache, only: spatial

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        real(dp), intent(out) :: rdm_spin(rdm%sign_length)

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        integer :: p_spat, q_spat, r_spat, s_spat
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_spin = 0.0_dp

        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            ! Obtain spin orbital labels.
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
            ! Obtain spatial orbital labels.
            p_spat = spatial(p); q_spat = spatial(q);
            r_spat = spatial(r); s_spat = spatial(s);

            ! Note to the reader for the following code: if mod(p,2) == 1 then
            ! p is a beta (b) orbital, if mod(p,2) == 0 then it is an alpha (a) obrital.

            ! The following if statement allows IJIJ spatial combinations.
            if (p_spat == r_spat .and. q_spat == s_spat) then
                ! If we get to this point then we definitely have a contribution
                ! to add in, so extract the sign.
                call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

                ! If all labels have the same spatial part (IIII):
                if (p_spat == q_spat) then
                    if (is_beta(p) .and. is_alpha(q) .and. is_beta(r) .and. is_alpha(s)) then
                        rdm_spin = rdm_spin - 6.0_dp*rdm_sign
                    end if

                else
                    ! We only get here if the spatial parts obey IJIJ, for I /= J:

                    ! The following if statement allows the following spin combinations:
                    ! aaaa, bbbb, abab and baba.
                    if (mod(p,2) == mod(r,2) .and. mod(q,2) == mod(s,2)) then

                        if (mod(p,2) == mod(q,2)) then
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

        rdm_spin = rdm_spin + 3.0_dp*real(nel,dp)*rdm_trace

        rdm_spin = rdm_spin/4.0_dp

    end subroutine calc_rdm_spin

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

    subroutine make_hermitian_rdm(rdm, spawn)

        type(rdm_list_t), intent(in) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        integer :: p_temp, q_temp
        real(dp) :: rdm_sign(rdm%sign_length)

        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
            call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

            ! If in the lower half of the RDM, reflect to the upper half.
            if (pq > rs) then
                p_temp = p; q_temp = q;
                p = r; q = s;
                r = p_temp; s = q_temp;
            end if

            call add_to_rdm_spawn_t(spawn, p, q, r, s, rdm_sign, .false.)
        end do

    end subroutine make_hermitian_rdm

    subroutine apply_spin_sym_to_rdm(rdm, spawn)

        ! Loop through all elements of the RDM, and whenever we have a bbbb
        ! term, move it to the equivalent aaaa element, and similarly for baba
        ! to abab and baab and abba (here a means alpha spin, b means beta spin).
        ! Also, only include terms above the diagonal - so average off-diagonal
        ! terms of the RDM. The ordering conforms with that of the printed
        ! RDMs from NECI prior to distribution of RDMs.

        use SystemData, only: nbasis
        use UMatCache, only: spatial

        type(rdm_list_t), intent(in) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn

        integer(int_rdm) :: pqrs
        integer :: i, pq, rs, p, q, r, s
        integer :: pq_legacy, rs_legacy
        integer :: p_temp, q_temp, r_temp
        integer :: p_spat, q_spat, r_spat, s_spat
        real(dp) :: rdm_sign(rdm%sign_length)

        do i = 1, rdm%nelements
            pqrs = rdm%elements(0,i)
            ! Obtain spin orbital labels and the RDM element.
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
            call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

            ! Get the spatial orbital labels from the spin orbital ones.
            p_spat = spatial(p); q_spat = spatial(q);
            r_spat = spatial(r); s_spat = spatial(s);

            ! When we calculate the combined labels, pq and rs, we would
            ! usually have p and q swapped below, and similarly with r and s.
            ! However, to be consistent with the old RDM code's printing
            ! system, we use the following. These will be used in deciding
            ! which of the equivalent elements we print when applying the
            ! symmetry of the RDM about the diagonal.
            ! Note we use spatial labels here, too.
            pq_legacy = (q_spat-1)*nbasis + p_spat
            rs_legacy = (s_spat-1)*nbasis + r_spat

            ! Apply symmetry (for *real* RDMs), to only print elements from
            ! one half of the RDM, using the legacy ordering.
            if (pq_legacy > rs_legacy) then
                p_temp = p; q_temp = q;
                p = r; q = s;
                r = p_temp; s = q_temp;
            end if

            ! If the first label has beta spin then we definitely want to move
            ! this RDM element to the equivalent alpha term, so go ahead and
            ! flip all the spins. Otherwise, move straight on.
            if (is_beta(p)) then
                ! The ab_pair macro swaps alpha and beta spins of a label while
                ! keeping the spatial orbital unchanged.
                p = ab_pair(p)
                q = ab_pair(q)
                r = ab_pair(r)
                s = ab_pair(s)
            end if

            ! If either p and q or r and s have the same spatial labels
            ! then we enforce the ordering abab by swapping orbitals and
            ! introducing minus signs as necessary.
            if (p_spat == q_spat .or. r_spat == s_spat) then
                ! If the above condition is met then we must have either abab
                ! or abba at this point. In the latter case, just swap the
                ! spins, and introduce a minus sign.
                ! * This logic is only valid if we are working in a fixed Ms
                ! spin sector *
                if (is_beta(r)) then
                    r = ab_pair(r)
                    s = ab_pair(s)
                    rdm_sign = -rdm_sign
                end if
            end if

            ! If the spatial parts of p and q are the same, and the spatial
            ! parts of r and s are *also* the same, then the RDM element won't
            ! have been added into both equivalent spin-flipped arrays, so we
            ! don't count twice. But otherwise, include a factor a half to
            ! account for double counting.
            if (.not. (is_in_pair(p,q) .and. is_in_pair(r,s))) rdm_sign = rdm_sign*0.5_dp

            ! If we're on a spatial off-diagonal element then the element is
            ! included twice from both above and below the diagonal. So divide
            ! by two to take the average.
            ! (Note that the element can also be above or below the diagonal
            ! by having the same spatial part and different spins, but this is
            ! taken care of when summing over spin-flipped configurations - no
            ! extra double counting occurs, so it is correct to consider
            ! spatial orbitals here, not spin orbitals).
            if (pq_legacy /= rs_legacy) rdm_sign = rdm_sign*0.5_dp

            call add_to_rdm_spawn_t(spawn, p, q, r, s, rdm_sign, .false.)
        end do

        call communicate_rdm_spawn_t(spawn)
        call annihilate_rdm_list(spawn%rdm_recv)

    end subroutine apply_spin_sym_to_rdm

    subroutine print_rdms_spin_sym_wrapper(rdm, spawn, rdm_trace)

        ! Compress the full spinned-RDMs by summing over spin-equivalent terms
        ! (i.e. aaaa and bbbb rdms), and also applying symmetry of (*real*)
        ! RDMs. The result will be stored in spawn%rdm_recv. Then, print with
        ! out to a file.

        use hash, only: clear_hash_table

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)

        spawn%free_slots = spawn%init_free_slots
        call clear_hash_table(spawn%rdm_send%hash_table)

        call apply_spin_sym_to_rdm(rdm, spawn)
        call print_rdms_with_spin(spawn%rdm_recv, rdm_trace)

    end subroutine print_rdms_spin_sym_wrapper

    subroutine create_spinfree_2rdm(rdm, spawn)

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

        use SystemData, only: nbasis
        use UMatCache, only: spatial

        type(rdm_list_t), intent(in) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn

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
            ! r < s, the following terms are only store with baba spin, never
            ! with abab. Double this term to make up for it.
            if (p_spat == q_spat .and. r_spat == s_spat) rdm_sign = 2.0_dp*rdm_sign

            ! Add all spinfree 2-RDM elements corresponding to these labels.
            call add_rdm_elements(p_spat, q_spat, r_spat, s_spat, rdm_sign, spawn)

            ! If this is an aaaa or bbbb term then *minus* this RDM element will
            ! be equal to the equivalent RDM element with the last two labels
            ! swapped. So, add this contribution into that RDM element. We
            ! don't have to do this, but doing so applies some extra averaging.
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

        call communicate_rdm_spawn_t(spawn)
        call annihilate_rdm_list(spawn%rdm_recv)

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
            ! esult.
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

    subroutine print_spinfree_2rdm_wrapper(rdm, spawn, rdm_trace)

        use hash, only: clear_hash_table

        type(rdm_list_t), intent(inout) :: rdm
        type(rdm_spawn_t), intent(inout) :: spawn
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)

        spawn%free_slots = spawn%init_free_slots
        call clear_hash_table(spawn%rdm_send%hash_table)

        call create_spinfree_2rdm(rdm, spawn)
        call print_spinfree_2rdm(spawn%rdm_recv, rdm_trace)

    end subroutine print_spinfree_2rdm_wrapper

    subroutine annihilate_rdm_list(rdm)

        ! Perform annihilation on RDM elements in rdm%elements.

        ! * WARNING * rdm%hash_table is not used by this routine, and will not
        ! be updated by it after reordering and compressing rdm%elements.

        use sort_mod, only: sort

        type(rdm_list_t), intent(inout) :: rdm

        integer :: i, counter
        integer(int_rdm) :: pqrs
        real(dp) :: rdm_sign(rdm%sign_length), summed_rdm_sign(rdm%sign_length)

        if (rdm%nelements > 0) then
            call sort(rdm%elements(:,1:rdm%nelements))

            summed_rdm_sign = 0.0_dp
            counter = 0

            do i = 1, rdm%nelements-1
                call extract_sign_rdm(rdm%elements(:,i), rdm_sign)
                summed_rdm_sign = summed_rdm_sign + rdm_sign

                ! Is the next RDM element the same as this one? If so then
                ! don't keep this element yet, but wait until all signs on
                ! this RDM element have been summed (annihilated).
                if (.not. (rdm%elements(0,i) == rdm%elements(0, i+1))) then
                    ! If this element is zero for all RDMs, then don't keep it.
                    if (any(abs(summed_rdm_sign) > 1.0e-12_dp)) then
                        counter = counter + 1
                        rdm%elements(0, counter) = rdm%elements(0,i)
                        call encode_sign_rdm(rdm%elements(:, counter), summed_rdm_sign)
                    end if
                    summed_rdm_sign = 0.0_dp
                end if
            end do

            ! Account for the final element separately.
            call extract_sign_rdm(rdm%elements(:,rdm%nelements), rdm_sign)
            summed_rdm_sign = summed_rdm_sign + rdm_sign
            if (any(abs(summed_rdm_sign) > 1.0-12_dp)) then
                counter = counter + 1
                rdm%elements(0, counter) = rdm%elements(0,rdm%nelements)
                call encode_sign_rdm(rdm%elements(:, counter), summed_rdm_sign)
            end if

            rdm%nelements = counter
        end if

    end subroutine annihilate_rdm_list

    subroutine print_rdms_with_spin(rdm, rdm_trace)

        ! Print the RDM stored in rdm to files, normalised by rdm_trace.

        ! This routine will print out *all* the spin cobminations separately,
        ! including both aaaa and bbbb arrays, and all other combinations.

        ! The files are called 'rdm_aaaa', 'rdm_abab', 'rdm_abba', etc...

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
                        open(iunit_aaaa, file='rdm_aaaa.'//trim(suffix), status='replace')
                        iunit_abab = get_free_unit()
                        open(iunit_abab, file='rdm_abab.'//trim(suffix), status='replace')
                        iunit_abba = get_free_unit()
                        open(iunit_abba, file='rdm_abba.'//trim(suffix), status='replace')
                        iunit_bbbb = get_free_unit()
                        open(iunit_bbbb, file='rdm_bbbb.'//trim(suffix), status='replace')
                        iunit_baba = get_free_unit()
                        open(iunit_baba, file='rdm_baba.'//trim(suffix), status='replace')
                        iunit_baab = get_free_unit()
                        open(iunit_baab, file='rdm_baab.'//trim(suffix), status='replace')
                    else
                        iunit_aaaa = get_free_unit()
                        open(iunit_aaaa, file='rdm_aaaa.'//trim(suffix), status='old', position='append')
                        iunit_abab = get_free_unit()
                        open(iunit_abab, file='rdm_abab.'//trim(suffix), status='old', position='append')
                        iunit_abba = get_free_unit()
                        open(iunit_abba, file='rdm_abba.'//trim(suffix), status='old', position='append')
                        iunit_bbbb = get_free_unit()
                        open(iunit_bbbb, file='rdm_bbbb.'//trim(suffix), status='old', position='append')
                        iunit_baba = get_free_unit()
                        open(iunit_baba, file='rdm_baba.'//trim(suffix), status='old', position='append')
                        iunit_baab = get_free_unit()
                        open(iunit_baab, file='rdm_baab.'//trim(suffix), status='old', position='append')
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
                    write(rdm_filename, '("new_spinfree_TwoRDM.",'//int_fmt(irdm,0)//')') irdm
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

end module rdm_parallel
