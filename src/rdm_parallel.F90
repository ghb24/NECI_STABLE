#include "macros.h"

module rdm_parallel

    ! This module contains all routines used for the calculations of reduced
    ! density matrices, when distributed across all MPI processes.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use rdm_data, only: rdm_spawn_t
    use util_mod

    implicit none

contains

    subroutine init_rdm_spawn_t(spawn, nrdms, nrows, contribs_length, nhashes_rdm)

        ! Initialise an rdm_spawn_t object.

        ! Out: spawn - rdm_spawn_t object to be initialised.
        ! In:  nrdms - the number of RDMs to be stored in the array.
        ! In:  nrows - the number of rows in the RDM.
        ! In:  contribs_length - the length of the RDM spawning arrays.
        ! In:  nhashes_rdm - the number of unique hashes for indexing the hash table.

        use hash, only: init_hash_table

        type(rdm_spawn_t), intent(out) :: spawn
        integer, intent(in) :: nrdms, nrows, contribs_length, nhashes_rdm
        integer :: i, ierr
        real(dp) :: slots_per_proc

        spawn%nrdms = nrdms
        spawn%nrows = nrows
        spawn%contribs_length = contribs_length
        spawn%nhashes_rdm = nhashes_rdm

        allocate(spawn%contribs(0:nrdms, contribs_length))
        allocate(spawn%contribs_recv(0:nrdms, contribs_length))

        allocate(spawn%hash_table(nhashes_rdm), stat=ierr)
        call init_hash_table(spawn%hash_table)

        allocate(spawn%free_slots(0:nProcessors-1), stat=ierr)
        allocate(spawn%init_free_slots(0:nProcessors-1), stat=ierr)

        ! Equally divide contribs across all processors.
        slots_per_proc = real(contribs_length, dp)/real(nProcessors, dp)
        do i = 0, nProcessors-1
            spawn%init_free_slots(i) = nint(slots_per_proc*i)+1
        end do

        ! Set the free slots array to its initial value.
        spawn%free_slots = spawn%init_free_slots

    end subroutine init_rdm_spawn_t

    pure subroutine calc_combined_rdm_labels(i, j, k, l, ijkl)

        ! Combine the four 2-RDM spin orbital labels into unique integers.
        ! i and j are combined into one number, ij. k and l are combined into
        ! one number, kl. Both of these are then combined into one single
        ! number, ijkl. The largest value for ijkl is M^4, where M is the
        ! number of spin orbitals.
        
        ! The compression defined in this routine will not give a fully
        ! compressed RDM index labelling, because it allows a separate ij
        ! integer if i and j are equal, which will always give a zero
        ! RDM element, and the same for k and l. It also doesn't take
        ! spatial symmetry into account. But this is fine if one just
        ! seeks a unique combined label for each combination of individual
        ! spin orbital labels.

        ! In: i, j, k l - spin orbitals of the RDM contribution.
        ! Out: ijkl - Label combining i, j, k and l.

        use SystemData, only: nbasis

        integer, intent(in) :: i, j, k, l
        integer(int_rdm), intent(out) :: ijkl
        integer :: ij, kl

        ij = (i-1)*nbasis + j
        kl = (k-1)*nbasis + l
        ijkl = (ij-1)*(nbasis**2) + kl

    end subroutine calc_combined_rdm_labels

    pure subroutine calc_separate_rdm_labels(ijkl, i, j, k, l)

        use SystemData, only: nbasis

        integer(int_rdm), intent(in) :: ijkl
        integer, intent(out) :: i, j, k, l
        integer :: ij, kl

        kl = mod(ijkl, nbasis**2)
        ij = (ijkl - kl)/(nbasis**2) + 1

        j = mod(ij, nbasis)
        i = (ij - i)/nbasis + 1

        l = mod(kl, nbasis)
        k = (kl - i)/nbasis + 1

    end subroutine calc_separate_rdm_labels

    pure subroutine extract_sign_rdm(rdm_entry, real_sign)

        integer(n_int), intent(in) :: rdm_entry(0:)
        real(dp), intent(out) :: real_sign(size(rdm_entry)-1)
        integer(n_int) :: int_sign(size(rdm_entry)-1)

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

        ! In/Out: rdm_spawn - the rdm_spawn_t type to which contributions will be added.
        ! In: p, q, r, s - spin orbitals of the RDM contribution, with p<q, r<s.
        ! In: contrib_sign - the sign (amplitude) of the contribution to be added.

        use hash, only: hash_table_lookup, add_hash_table_entry

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: p, q, r, s
        real(dp), intent(in) :: contrib_sign(spawn%nrdms)

        integer :: pq_compressed, proc, ind, hash_val
        integer(int_rdm) :: pqrs
        real(dp) :: real_sign_old(spawn%nrdms), real_sign_new(spawn%nrdms)
        logical :: tSuccess, list_full
        character(*), parameter :: this_routine = 'add_to_rdm_spawn_t'

        ! Calculate combined RDM labels.
        call calc_combined_rdm_labels(p, q, r, s, pqrs)

        ! Search to see if this RDM element is already in the contribs array.
        ! If it, tSuccess will be true and ind will hold the position of the
        ! entry in contribs.
        call hash_table_lookup((/p,q,r,s/), (/pqrs/), 0, spawn%hash_table, spawn%contribs, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! Extract the existing sign.
            call extract_sign_rdm(spawn%contribs(:,ind), real_sign_old)
            ! Update the total sign.
            real_sign_new = real_sign_old + contrib_sign
            ! Encode the new sign.
            call encode_sign_rdm(spawn%contribs(:,ind), real_sign_new)
        else
            ! The following maps (p,q), with p<q, to single integers with no gaps.
            ! It is benefical to have no gaps here, for good load balancing.
            pq_compressed = (q-1)*(q-2)/2 + p
            ! Calculate processor for the element.
            proc = (pq_compressed-1)*nProcessors/spawn%nrows

            ! Check that there is enough memory for the new spawned RDM entry.
            list_full = .false.
            if (proc == nProcessors - 1) then
                if (spawn%free_slots(proc) > spawn%contribs_length) list_full = .true.
            else
                if (spawn%free_slots(proc) > spawn%init_free_slots(proc+1)) list_full = .true.
            end if
            if (list_full) then
                write(6,'("Attempting to add an RDM contribution to the spawned list on processor:",&
                           &1X,'//int_fmt(proc,0)//')') proc
                write(6,'("No memory slots available for this spawn.")')
                call stop_all(this_routine, "Out of memory for spawned RDM contributions.")
            end if

            spawn%contribs(0, spawn%free_slots(proc)) = pqrs
            call encode_sign_rdm(spawn%contribs(:, spawn%free_slots(proc)), contrib_sign)

            call add_hash_table_entry(spawn%hash_table, spawn%free_slots(proc), hash_val)

            spawn%free_slots(proc) = spawn%free_slots(proc) + 1
            write(6,*) "num_entries:", spawn%free_slots(proc), "on proc", proc; flush(6)
        end if

        if (p > q .or. r > s) call stop_all("Here","Here")

    end subroutine add_to_rdm_spawn_t

    subroutine calc_rdm_energy(spawn)

        use rdm_data, only: rdm_spawn_t
        use rdm_integral_fns, only: one_elec_int, two_elec_int
        use SystemData, only: nel

        type(rdm_spawn_t), intent(inout) :: spawn

        integer :: i, pqrs, p, q, r, s
        real(dp) :: rdm_sign(spawn%nrdms), rdm_energy(spawn%nrdms)

        rdm_energy = 0.0_dp

        ! Just working with one processor and the spawned list for now,
        ! as a check - this will be corrected.
        do i = 1, spawn%free_slots(0)-1
            pqrs = spawn%contribs(0,i)
            call calc_separate_rdm_labels(pqrs, p, q, r, s)
            call extract_sign_rdm(spawn%contribs(:,i), rdm_sign)

            ! The 2-RDM contribution to the energy:
            rdm_energy = rdm_energy + rdm_sign*two_elec_int(p,q,r,s)
            ! The 1-RDM contribution to the energy:
            if (p == r) rdm_energy = rdm_energy + rdm_sign*one_elec_int(q,s)/(nel-1)
            if (q == s) rdm_energy = rdm_energy + rdm_sign*one_elec_int(p,r)/(nel-1)
            if (p == s) rdm_energy = rdm_energy - rdm_sign*one_elec_int(q,r)/(nel-1)
            if (q == r) rdm_energy = rdm_energy - rdm_sign*one_elec_int(p,s)/(nel-1)
        end do

    end subroutine calc_rdm_energy

end module rdm_parallel
