#include "macros.h"

module rdm_parallel

    ! This module contains all routines used for the calculations of reduced
    ! density matrices, when distributed across all MPI processes.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use rdm_data, only: rdm_list_t, rdm_spawn_t
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

    subroutine init_rdm_spawn_t(spawn, nrows, sign_length, max_nelements, max_nelements_recv, nhashes)

        ! Initialise an rdm_spawn_t object.

        ! Out: spawn - rdm_spawn_t object to be initialised.
        ! In:  nrows - the number of rows in the RDM.
        ! In:  sign_length - the number of signs which can be stored for each element.
        ! In:  max_nelements - the length of the spawn%rdm%elements array.
        ! In:  max_nelements_recv - the length of the spawn%rdm_recv array.
        ! In:  nhashes - the number of unique hashes for indexing spawn%rdm%hash_table.

        type(rdm_spawn_t), intent(out) :: spawn
        integer, intent(in) :: nrows, sign_length, nhashes
        integer, intent(in) :: max_nelements, max_nelements_recv
        integer :: i, ierr
        real(dp) :: slots_per_proc

        spawn%nrows = nrows
        spawn%max_nelements_recv = max_nelements_recv

        call init_rdm_list_t(spawn%rdm, sign_length, max_nelements, nhashes)

        allocate(spawn%rdm_recv(0:sign_length, max_nelements_recv))
        spawn%rdm_recv = 0_int_rdm

        allocate(spawn%free_slots(0:nProcessors-1), stat=ierr)
        allocate(spawn%init_free_slots(0:nProcessors-1), stat=ierr)

        ! Equally divide RDM rows across all processors.
        slots_per_proc = real(max_nelements, dp)/real(nProcessors, dp)
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

    subroutine add_to_rdm_spawn_t(spawn, p, q, r, s, contrib_sign)

        ! In/Out: rdm_spawn - the rdm_spawn_t object to which contributions will be added.
        ! In: p, q, r, s - spin orbitals of the RDM contribution, with p<q, r<s.
        ! In: contrib_sign - the sign (amplitude) of the contribution to be added.

        use hash, only: hash_table_lookup, add_hash_table_entry
        use SystemData, only: nbasis

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: p, q, r, s
        real(dp), intent(in) :: contrib_sign(spawn%rdm%sign_length)

        integer :: pq_compressed, proc, ind, hash_val
        integer(int_rdm) :: pqrs
        real(dp) :: real_sign_old(spawn%rdm%sign_length), real_sign_new(spawn%rdm%sign_length)
        logical :: tSuccess, list_full
        character(*), parameter :: this_routine = 'add_to_rdm_spawn_t'

        ! Calculate combined RDM labels.
        call calc_combined_rdm_label(p, q, r, s, pqrs)

        ! Search to see if this RDM element is already in the RDM array.
        ! If it, tSuccess will be true and ind will hold the position of the
        ! entry in spawn%rdm%elements.
        call hash_table_lookup((/p,q,r,s/), (/pqrs/), 0, spawn%rdm%hash_table, spawn%rdm%elements, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! Extract the existing sign.
            call extract_sign_rdm(spawn%rdm%elements(:,ind), real_sign_old)
            ! Update the total sign.
            real_sign_new = real_sign_old + contrib_sign
            ! Encode the new sign.
            call encode_sign_rdm(spawn%rdm%elements(:,ind), real_sign_new)
        else
            ! The following maps (p,q), with p<q, to single integers with no gaps.
            ! It is benefical to have no gaps here, for good load balancing.
            ! The final integers are ordered so that p is dominant over q.
            pq_compressed = nbasis*(p-1) - p*(p-1)/2 + q - p
            ! Calculate processor for the element.
            proc = (pq_compressed-1)*nProcessors/spawn%nrows

            ! Check that there is enough memory for the new spawned RDM entry.
            list_full = .false.
            if (proc == nProcessors - 1) then
                if (spawn%free_slots(proc) > spawn%rdm%max_nelements) list_full = .true.
            else
                if (spawn%free_slots(proc) > spawn%init_free_slots(proc+1)) list_full = .true.
            end if
            if (list_full) then
                write(6,'("Attempting to add an RDM contribution to the spawned list on processor:",&
                           &1X,'//int_fmt(proc,0)//')') proc
                write(6,'("No memory slots available for this spawn.")')
                call stop_all(this_routine, "Out of memory for spawned RDM contributions.")
            end if

            spawn%rdm%elements(0, spawn%free_slots(proc)) = pqrs
            call encode_sign_rdm(spawn%rdm%elements(:, spawn%free_slots(proc)), contrib_sign)

            call add_hash_table_entry(spawn%rdm%hash_table, spawn%free_slots(proc), hash_val)

            spawn%free_slots(proc) = spawn%free_slots(proc) + 1
        end if

        if (p > q .or. r > s) then
            write(6,'("p, q, r, s:", 1X,'//int_fmt(p,0)//', 1X,'//int_fmt(q,0)//', &
                       &1X,'//int_fmt(r,0)//', 1X,'//int_fmt(s,0)//')') p, q, r, s
            call stop_all(this_routine,"Incorrect ordering of RDM orbitals passed to RDM spawning routine.")
        end if

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
        spawn%nelements_recv = int(sum(recv_sizes), sizeof_int)

        if (spawn%nelements_recv > spawn%max_nelements_recv) then
            write(6,'("Attempting to receive RDM elements on processor:",&
                       &1X,'//int_fmt(iProcIndex,0)//')') iProcIndex
            write(6,'("Insufficient space in the receiving RDM array.")')
            call stop_all(this_routine, "Not enough memory to communicate RDM elements.")
        end if

        send_sizes = send_sizes*size(spawn%rdm%elements,1)
        recv_sizes = recv_sizes*size(spawn%rdm%elements,1)
        send_displs = send_displs*size(spawn%rdm%elements,1)
        recv_displs = recv_displs*size(spawn%rdm%elements,1)

        ! Perform the communication.
        call MPIAlltoAllv(spawn%rdm%elements, send_sizes, send_displs, &
                          spawn%rdm_recv, recv_sizes, recv_displs, ierr)

        ! Now we can reset the spawn%rdm%elements and free_slots array.
        spawn%free_slots = spawn%init_free_slots
        call clear_hash_table(spawn%rdm%hash_table)

    end subroutine communicate_rdm_spawn_t

    subroutine add_rdm_1_to_rdm_2(rdm_1_elems, num_rdm_1_elems, rdm_2)

        use hash, only: hash_table_lookup, add_hash_table_entry

        integer(int_rdm), intent(in) :: rdm_1_elems(0:,:)
        integer, intent(in) :: num_rdm_1_elems
        type(rdm_list_t), intent(inout) :: rdm_2

        integer :: i, pq, rs, p, q, r, s
        integer(int_rdm) :: pqrs
        integer :: ind, hash_val
        real(dp) :: real_sign_old(rdm_2%sign_length), real_sign_new(rdm_2%sign_length)
        real(dp) :: spawn_sign(rdm_2%sign_length)
        logical :: tSuccess
        character(*), parameter :: this_routine = 'add_rdm_1_to_rdm_2'

        do i = 1, num_rdm_1_elems
            ! Decode the compressed RDM labels.
            pqrs = rdm_1_elems(0,i)
            call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)

            ! Extract the spawned sign.
            call extract_sign_rdm(rdm_1_elems(:,i), spawn_sign)

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
            p_spat = (p-1)/2 + 1
            q_spat = (q-1)/2 + 1
            r_spat = (r-1)/2 + 1
            s_spat = (s-1)/2 + 1

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

    subroutine print_rdms_with_spin(rdm, rdm_trace)

        use Parallel_neci, only: MPIBarrier
        use sort_mod, only: sort
        use util_mod, only: get_free_unit

        type(rdm_list_t), intent(inout) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)

        integer(int_rdm) :: pqrs
        integer :: i, ierr, iproc, write_unit
        integer :: iunit_aaaa, iunit_abab, iunit_abba
        integer :: iunit_bbbb, iunit_baba, iunit_baab
        integer :: pq, rs, p, q, r, s
        integer :: p_spat, q_spat, r_spat, s_spat
        real(dp) :: rdm_sign(rdm%sign_length)
        character(3) :: sgn_len

        ! Store rdm%sign_length as a string, for the formatting string.
        write(sgn_len,'(i3)') rdm%sign_length

        call sort(rdm%elements(:,1:rdm%nelements))

        do iproc = 0, nProcessors-1
            if (iproc == iProcIndex) then

                ! Open all the files to be written to:
                ! Let the first processor clear all the files to start with.
                if (iproc == 0) then
                    iunit_aaaa = get_free_unit()
                    open(iunit_aaaa, file='test_aaaa', status='replace')
                    iunit_abab = get_free_unit()
                    open(iunit_abab, file='test_abab', status='replace')
                    iunit_abba = get_free_unit()
                    open(iunit_abba, file='test_abba', status='replace')
                    iunit_bbbb = get_free_unit()
                    open(iunit_bbbb, file='test_bbbb', status='replace')
                    iunit_baba = get_free_unit()
                    open(iunit_baba, file='test_baba', status='replace')
                    iunit_baab = get_free_unit()
                    open(iunit_baab, file='test_baab', status='replace')
                else
                    iunit_aaaa = get_free_unit()
                    open(iunit_aaaa, file='test_aaaa', status='old', position='append')
                    iunit_abab = get_free_unit()
                    open(iunit_abab, file='test_abab', status='old', position='append')
                    iunit_abba = get_free_unit()
                    open(iunit_abba, file='test_abba', status='old', position='append')
                    iunit_bbbb = get_free_unit()
                    open(iunit_bbbb, file='test_bbbb', status='old', position='append')
                    iunit_baba = get_free_unit()
                    open(iunit_baba, file='test_baba', status='old', position='append')
                    iunit_baab = get_free_unit()
                    open(iunit_baab, file='test_baab', status='old', position='append')
                end if

                do i = 1, rdm%nelements
                    pqrs = rdm%elements(0,i)
                    ! Obtain spin orbital labels.
                    call calc_separate_rdm_labels(pqrs, pq, rs, p, q, r, s)
                    call extract_sign_rdm(rdm%elements(:,i), rdm_sign)

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

                    p_spat = (p-1)/2 + 1
                    q_spat = (q-1)/2 + 1
                    r_spat = (r-1)/2 + 1
                    s_spat = (s-1)/2 + 1

                    write(write_unit,'(4i6,'//trim(sgn_len)//'g25.17)') p_spat, q_spat, r_spat, s_spat, rdm_sign/rdm_trace
                end do

                close(iunit_aaaa)
                close(iunit_abab)
                close(iunit_abba)
                close(iunit_bbbb)
                close(iunit_baba)
                close(iunit_baab)
                close(iunit_abab)
                close(iunit_abba)
            end if

            ! Wait for the current processor to finish printing its RDM elements.
            call MPIBarrier(ierr)
        end do

    end subroutine print_rdms_with_spin

end module rdm_parallel
