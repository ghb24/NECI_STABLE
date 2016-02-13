#include "macros.h"

module rdm_parallel

    ! This module contains all routines used for the calculations of reduced
    ! density matrices, when distributed across all MPI processes.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use SystemData, only: nel, nbasis
    use util_mod

    implicit none

contains

    subroutine init_parallel_rdms(nrdms)

        ! This routine initialises all arrays and data used for the parallel
        ! implementation of RDMs.

        use hash, only: init_hash_table
        use rdm_data, only: rdm_arr , rdm_arr_ht, rdm_spawn, rdm_spawn_ht
        use rdm_data, only: rdm_spawn_recv, rdm_arr_length, nhashes_rdm, rdm_nrows

        integer, intent(in) :: nrdms
        integer :: ierr

        rdm_nrows = nbasis*(nbasis-1)/2

        ! For now, create RDM arrays big enough so that *all* RDM elements on
        ! a particular processor can be stored, using the usual approximations
        ! to take symmetry into account. Include a factor of 1.2 to account for
        ! factors such as imperfect load balancing (which affects the spawned
        ! array).
        rdm_arr_length = 1.2*rdm_nrows/(32*nProcessors)
        nhashes_rdm = 0.8*rdm_arr_length

        allocate(rdm_spawn(0:nrdms, rdm_arr_length))
        allocate(rdm_spawn_recv(0:nrdms, rdm_arr_length))
        allocate(rdm_arr(0:nrdms, rdm_arr_length))

        ! Allocate and initialise hash tables.
        allocate(rdm_spawn_ht(nhashes_rdm), stat=ierr)
        call init_hash_table(rdm_spawn_ht)

        allocate(rdm_spawn_ht(nhashes_rdm), stat=ierr)
        call init_hash_table(rdm_arr_ht)

    end subroutine init_parallel_rdms

    pure subroutine calc_combined_rdm_labels(i, j, k, l, ij, kl, ijkl)

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
        ! Out: ij - Label combining i and j.
        ! Out: kl - Label combining k and l.
        ! Out: ijkl - Label combining i, j, k and l.

        integer, intent(in) :: i, j, k, l
        integer, intent(out) :: ij, kl
        integer(int_rdm), intent(out) :: ijkl

        ij = (i-1)*nbasis + j
        kl = (k-1)*nbasis + l
        ijkl = (ij-1)*(nbasis**2) + kl

    end subroutine calc_combined_rdm_labels

    subroutine calc_separate_rdm_labels(ijkl, i, j, k, l)

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

    subroutine add_to_rdm_spawn(i, j, k, l, nrdms, contrib_sign, spawn_arr, spawn_arr_ht, spawn_slots, init_spawn_slots)

        ! In: i, j, k l - spin orbitals of the RDM contribution.
        ! In: contrib_sign - the amplitude of the contribution to the RDM.
        ! In/Out: spawn_arr - the RDM spawning array to add to.
        ! In/Out: spawn_arr_ht - the hash table to spawn_arr.
        ! In/Out: spawn_slots - array holding the next free slots in
        !             spawn_arr for each processor.
        ! In: init_spawn_slots - array holding the positions of the first
        !         entry in spawn_arr for each processor.

        use hash, only: hash_table_lookup, add_hash_table_entry
        use rdm_data, only: rdm_nrows

        integer, intent(in) :: i, j, k, l, nrdms
        real(dp), intent(in) :: contrib_sign(nrdms)
        integer(int_rdm), intent(inout) :: spawn_arr(0:,:)
        type(ll_node), pointer, intent(inout) :: spawn_arr_ht(:)
        integer, intent(inout) :: spawn_slots(:)
        integer, intent(in) :: init_spawn_slots(:)

        integer :: ij, kl, proc, ind, hash_val
        integer(int_rdm) :: ijkl
        real(dp) :: real_sign_old(nrdms), real_sign_new(nrdms)
        logical :: tSuccess, list_full
        character(*), parameter :: this_routine = 'add_to_rdm_spawn'

        ! Calculate combined RDM labels.
        call calc_combined_rdm_labels(i, j, k, l, ij, kl, ijkl)

        ! Search to see if this RDM element is already in the spawn_arr array.
        ! If it, tSuccess will be true and ind will hold the position of the
        ! entry in spawn_arr.
        call hash_table_lookup((/i,j,k,l/), (/ijkl/), 0, spawn_arr_ht, spawn_arr, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! Extract the existing sign.
            call extract_sign_rdm(spawn_arr(:,ind), real_sign_old)
            ! Update the total sign.
            real_sign_new = real_sign_old + contrib_sign
            ! Encode the new sign.
            call encode_sign_rdm(spawn_arr(:,ind), real_sign_new)
        else
            ! Calculate processor for the element.
            proc = ij*nProcessors/rdm_nrows

            ! Check that there is enough memory for the new spawned RDM entry.
            list_full = .false.
            if (proc == nProcessors - 1) then
                if (spawn_slots(proc) > size(spawn_arr,2)) list_full = .true.
            else
                if (spawn_slots(proc) > init_spawn_slots(proc+1)) list_full = .true.
            end if
            if (list_full) then
                write(6,'("Attempting to add an RDM contribution to the spawned list on processor:",&
                           &1X,'//int_fmt(proc,0)//')') proc
                write(6,'("No memory slots available for this spawn.")')
                call stop_all(this_routine, "Out of memory for spawned RDM contributions.")
            end if

            spawn_arr(0, spawn_slots(proc)) = ijkl
            call encode_sign_rdm(spawn_arr(:,spawn_slots(proc)), contrib_sign)

            call add_hash_table_entry(spawn_arr_ht, spawn_slots(proc), hash_val)

            spawn_slots(proc) = spawn_slots(proc) + 1
        end if

    end subroutine add_to_rdm_spawn

end module rdm_parallel
