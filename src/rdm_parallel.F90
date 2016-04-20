#include "macros.h"

module rdm_parallel

    ! This module contains all routines used for the calculations of reduced
    ! density matrices, when distributed across all MPI processes.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use rdm_data, only: rdm_list_t, rdm_spawn_t, one_rdm_t
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

    subroutine init_rdm_spawn_t(spawn, nrows, sign_length, max_nelements_send, nhashes)

        ! Initialise an rdm_spawn_t object.

        ! Out: spawn - rdm_spawn_t object to be initialised.
        ! In:  nrows - the number of rows in the RDM.
        ! In:  sign_length - the number of signs which can be stored for each element.
        ! In:  max_nelements_send - the length of the spawn%rdm_send%elements array.
        ! In:  nhashes - the number of unique hashes for indexing spawn%rdm_send%hash_table.

        type(rdm_spawn_t), intent(out) :: spawn
        integer, intent(in) :: nrows, sign_length, nhashes
        integer, intent(in) :: max_nelements_send
        integer :: i, ierr
        real(dp) :: slots_per_proc

        spawn%nrows = nrows

        call init_rdm_list_t(spawn%rdm_send, sign_length, max_nelements_send, nhashes)

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

    subroutine communicate_rdm_spawn_t(spawn, rdm_recv)

        use hash, only: clear_hash_table
        use Parallel_neci, only: MPIAlltoAll, MPIAlltoAllv

        type(rdm_spawn_t), intent(inout) :: spawn
        type(rdm_list_t), intent(inout) :: rdm_recv

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
        rdm_recv%nelements = int(sum(recv_sizes), sizeof_int)

        if (rdm_recv%nelements > rdm_recv%max_nelements) then
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
                          rdm_recv%elements, recv_sizes, recv_displs, ierr)

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

end module rdm_parallel
