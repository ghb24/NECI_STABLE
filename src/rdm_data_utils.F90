#include "macros.h"

module rdm_data_utils

    ! This module contains routines for performing common and general
    ! operations with the RDM data structures defined in rdm_data. These
    ! data structures are one_rdm_t, rdm_list_t and rdm_spawn_t.

    ! For example, initialisation routines for these objects are held here,
    ! as are routines to add RDMs together, perform annihilation within an
    ! rdm_list_t object, and adding to and communicating with rdm_spawn_t
    ! objects.

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use Parallel_neci, only: iProcIndex, nProcessors
    use rdm_data, only: rdm_list_t, rdm_spawn_t, one_rdm_t, en_pert_t
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

        integer :: iproc, ierr
        real(dp) :: slots_per_proc

        spawn%nrows = nrows

        call init_rdm_list_t(spawn%rdm_send, sign_length, max_nelements_send, nhashes)

        allocate(spawn%free_slots(0:nProcessors-1), stat=ierr)
        ! init_free_slots has one extra element compared to free_slots. This
        ! is set equal to the total number of elements + 1, which allows us to
        ! avoid an extra if-statement for an edge case in add_to_rdm_spawn_t.
        allocate(spawn%init_free_slots(0:nProcessors), stat=ierr)

        ! Equally divide RDM rows across all processors.
        slots_per_proc = real(max_nelements_send, dp)/real(nProcessors, dp)
        do iproc = 0, nProcessors-1
            spawn%init_free_slots(iproc) = nint(slots_per_proc*iproc)+1
        end do
        ! For edge cases - see comment above.
        spawn%init_free_slots(nProcessors) = max_nelements_send + 1

        ! Set the free slots array to its initial value.
        spawn%free_slots = spawn%init_free_slots(0:nProcessors-1)

    end subroutine init_rdm_spawn_t

    subroutine init_en_pert_t(en_pert, sign_length, max_ndets, nhashes)

        ! Initialise an en_pert_t object.

        ! Out: en_pert - en_pert_t object to be initialised.
        ! In:  sign_length - the number of signs which can be stored for each element.
        ! In:  max_ndets - the length of the en_pert%dets array.
        ! In:  nhashes - the number of unique hashes for indexing the hash table.

        use hash, only: init_hash_table

        type(en_pert_t), intent(out) :: en_pert
        integer, intent(in) :: sign_length, max_ndets, nhashes

        integer :: ierr

        en_pert%sign_length = sign_length
        en_pert%max_ndets = max_ndets
        en_pert%nhashes = nhashes
        en_pert%ndets = 0

        allocate(en_pert%dets(0:sign_length+NIfDBO, max_ndets))
        en_pert%dets = 0_n_int

        allocate(en_pert%hash_table(nhashes), stat=ierr)
        call init_hash_table(en_pert%hash_table)

    end subroutine init_en_pert_t

    subroutine init_one_rdm_t(one_rdm, norbs)

        ! Initialise a one_rdm_t object.

        ! Out: one_rdm - one_rdm_t object to be initialised.
        ! In:  norbs - the number of orbitals to be indexed in the 1-RDM.

        use rdm_data, only: one_rdm_t

        type(one_rdm_t), intent(out) :: one_rdm
        integer, intent(in) :: norbs

        integer :: ierr
        character(*), parameter :: t_r = 'init_one_rdm_t'

        allocate(one_rdm%matrix(norbs, norbs), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, 'Problem allocating 1-RDM array.')
        call LogMemAlloc('one_rdm%matrix', norbs**2, 8, t_r, one_rdm%matrix_tag, ierr)
        one_rdm%matrix = 0.0_dp

        allocate(one_rdm%evalues(norbs), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, 'Problem allocating evalues array,')
        call LogMemAlloc('one_rdm%evalues', norbs, 8, t_r, one_rdm%evalues_tag, ierr)
        one_rdm%evalues = 0.0_dp

        allocate(one_rdm%rho_ii(norbs), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, 'Problem allocating 1-RDM diagonal array (rho_ii).')
        call LogMemAlloc('one_rdm%rho_ii', norbs, 8, t_r, one_rdm%rho_ii_tag, ierr)
        one_rdm%rho_ii = 0.0_dp

        allocate(one_rdm%sym_list_no(norbs), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, 'Problem allocating sym_list_no array.')
        allocate(one_rdm%sym_list_inv_no(norbs), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, 'Problem allocating sym_list_inv_no array.')

    end subroutine init_one_rdm_t

    subroutine init_rdm_definitions_t(rdm_defs, nrdms_standard, nrdms_transition, states_for_transition_rdm,filename)

        ! Set up arrays which define which RDMs and tRDMs are being sampled by
        ! specifying which states and FCIQMC simualtions are involved in those
        ! RDMs.

        ! For example, state_labels specifies which states are used to
        ! construct each RDM. sim_labels specifies which FCIQMC simulations
        ! are used to construct each RDM. See the rdm_data module for more
        ! description on the various arrays set up here.

        use rdm_data, only: rdm_definitions_t

        type(rdm_definitions_t), intent(out) :: rdm_defs
        integer, intent(in) :: nrdms_standard, nrdms_transition
        integer, intent(in) :: states_for_transition_rdm(:,:) ! (2, nrdms_transition/nreplicas)
        character(*), intent(in), optional :: filename

        integer :: nrdms, irdm, counter

        nrdms = nrdms_standard + nrdms_transition
        rdm_defs%nrdms = nrdms
        rdm_defs%nrdms_standard = nrdms_standard
        rdm_defs%nrdms_transition = nrdms_transition

        ! state_labels(:,j) will store the labels of the *actual* wave
        ! functions (i.e., which eigenstate is being sampled) contributing to
        ! the j'th RDM.
        allocate(rdm_defs%state_labels(2, nrdms))
        ! The 'standard' (non-transition) RDMs.
        do irdm = 1, nrdms_standard
            rdm_defs%state_labels(1,irdm) = irdm
            rdm_defs%state_labels(2,irdm) = irdm
        end do
        ! The transition RDMs - these were passed in in
        ! states_for_transition_rdm, which in practice is read in from the
        ! input and then passed to this routine.
        if (nrdms_transition > 0) then
            if (nreplicas == 1) then
                rdm_defs%state_labels(:,nrdms_standard+1:nrdms_standard+nrdms_transition) = states_for_transition_rdm
            else if (nreplicas == 2) then
                do irdm = 2, nrdms_transition, 2
                    ! In this case, there are two transition RDMs sampled for
                    ! each one the user requested, because there are two
                    ! combinations of replicas which can be used.
                    rdm_defs%state_labels(:,nrdms_standard+irdm-1) = states_for_transition_rdm(:,irdm/2)
                    rdm_defs%state_labels(:,nrdms_standard+irdm)   = states_for_transition_rdm(:,irdm/2)
                end do
            end if
        end if

        ! For transition RDMs, with 2 replicas for each state, there will be 2
        ! copies of each transition RDM. This array simply holds which of the
        ! 2 each RDM is - the first or second repeat.
        allocate(rdm_defs%repeat_label(nrdms))
        do irdm = 1, nrdms_standard
            rdm_defs%repeat_label(irdm) = 1
        end do
        do irdm = 1, nrdms_transition
            rdm_defs%repeat_label(irdm+nrdms_standard) = nreplicas - mod(irdm,nreplicas)
        end do

        ! sim_labels(:,j) will store the labels of the *FCIQMC* wave functions
        ! (i.e. the 'replica' labels) which will be used to sample the j'th RDM
        ! being calculated.
        allocate(rdm_defs%sim_labels(2, nrdms))
        do irdm = 1, nrdms
            if (nreplicas == 1) then
                rdm_defs%sim_labels(1,irdm) = rdm_defs%state_labels(1,irdm)
                rdm_defs%sim_labels(2,irdm) = rdm_defs%state_labels(2,irdm)
            else if (nreplicas == 2) then
                rdm_defs%sim_labels(1,irdm) = rdm_defs%state_labels(1,irdm)*nreplicas-mod(rdm_defs%repeat_label(irdm),2)
                rdm_defs%sim_labels(2,irdm) = rdm_defs%state_labels(2,irdm)*nreplicas-mod(rdm_defs%repeat_label(irdm)+1,2)
            end if
        end do

        allocate(rdm_defs%nrdms_per_sim(lenof_sign))
        allocate(rdm_defs%sim_pairs(nrdms, lenof_sign))
        allocate(rdm_defs%rdm_labels(nrdms, lenof_sign))
        rdm_defs%nrdms_per_sim = 0
        rdm_defs%sim_pairs = 0
        rdm_defs%rdm_labels = 0
        do irdm = 1, nrdms
            ! Count the number of times each replica label is contributing to
            ! an RDM.
            rdm_defs%nrdms_per_sim(rdm_defs%sim_labels(2,irdm))  = rdm_defs%nrdms_per_sim(rdm_defs%sim_labels(2,irdm)) + 1
            ! The number of RDMs which we have currently counted, for this replica.
            counter = rdm_defs%nrdms_per_sim(rdm_defs%sim_labels(2,irdm))
            ! Store which replica is paired with this, for this particular RDM
            ! sampling.
            rdm_defs%sim_pairs(counter, rdm_defs%sim_labels(2,irdm)) = rdm_defs%sim_labels(1,irdm)
            ! Store which RDM this pair of signs contributes to.
            rdm_defs%rdm_labels(counter, rdm_defs%sim_labels(2,irdm)) = irdm

            ! Do the same as above, but now for cases when spawning from the
            ! 'second' replica in the pair, *but only if it is different*.
            if (rdm_defs%sim_labels(1,irdm) /= rdm_defs%sim_labels(2,irdm)) then
                rdm_defs%nrdms_per_sim(rdm_defs%sim_labels(1,irdm))  = &
                    rdm_defs%nrdms_per_sim(rdm_defs%sim_labels(1,irdm)) + 1
                counter = rdm_defs%nrdms_per_sim(rdm_defs%sim_labels(1,irdm))
                rdm_defs%sim_pairs(counter, rdm_defs%sim_labels(1,irdm)) = rdm_defs%sim_labels(2,irdm)
                rdm_defs%rdm_labels(counter, rdm_defs%sim_labels(1,irdm)) = irdm
            end if
        end do

        if(present(filename)) then
           rdm_defs%output_file_prefix = filename
        else
           rdm_defs%output_file_prefix = 'TwoRDM'
        endif

    end subroutine init_rdm_definitions_t

    subroutine dealloc_rdm_list_t(rdm)

        ! Deallocate an rdm_list_t object.

        ! In/Out: rdm - rdm_list_t object to be deallocated.

        use hash, only: clear_hash_table

        type(rdm_list_t), intent(inout) :: rdm

        integer :: ierr

        if (allocated(rdm%elements)) deallocate(rdm%elements, stat=ierr)

        call clear_hash_table(rdm%hash_table)
        deallocate(rdm%hash_table, stat=ierr)
        nullify(rdm%hash_table)

    end subroutine dealloc_rdm_list_t

    subroutine dealloc_rdm_spawn_t(spawn)

        ! Deallocate an rdm_spawn_t object.

        ! In/Out: spawn - rdm_spawn_t object to be deallocated.

        use hash, only: clear_hash_table

        type(rdm_spawn_t), intent(inout) :: spawn

        integer :: ierr

        call dealloc_rdm_list_t(spawn%rdm_send)

        if (allocated(spawn%free_slots)) deallocate(spawn%free_slots, stat=ierr)
        if (allocated(spawn%init_free_slots)) deallocate(spawn%init_free_slots, stat=ierr)

    end subroutine dealloc_rdm_spawn_t

    subroutine dealloc_en_pert_t(en_pert)

        ! Deallocate an en_pert_t object.

        ! In/Out: en_pert - en_pert_t object to be deallocated.

        use hash, only: clear_hash_table

        type(en_pert_t), intent(inout) :: en_pert

        integer :: ierr

        if (allocated(en_pert%dets)) deallocate(en_pert%dets, stat=ierr)

        call clear_hash_table(en_pert%hash_table)
        deallocate(en_pert%hash_table, stat=ierr)
        nullify(en_pert%hash_table)

    end subroutine dealloc_en_pert_t

    subroutine dealloc_one_rdm_t(one_rdm)

        ! Deallocate a one_rdm_t object.

        ! In/Out: one_rdm - one_rdm_t object to be deallocated.

        use rdm_data, only: one_rdm_t

        type(one_rdm_t), intent(inout) :: one_rdm

        integer :: ierr
        character(*), parameter :: t_r = 'dealloc_one_rdm_t'

        if (allocated(one_rdm%matrix)) then
            deallocate(one_rdm%matrix, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem deallocating 1-RDM array.')
            call LogMemDeAlloc(t_r, one_rdm%matrix_tag)
        end if
        if (allocated(one_rdm%evalues)) then
            deallocate(one_rdm%evalues, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem deallocating evalues array,')
            call LogMemDeAlloc(t_r, one_rdm%evalues_tag)
        end if
        if (allocated(one_rdm%rho_ii)) then
            deallocate(one_rdm%rho_ii, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem deallocating 1-RDM diagonal array (rho_ii).')
            call LogMemDeAlloc(t_r, one_rdm%rho_ii_tag)
        end if
        if (allocated(one_rdm%sym_list_no)) then
            deallocate(one_rdm%sym_list_no, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem deallocating sym_list_no array.')
        end if
        if (allocated(one_rdm%sym_list_inv_no)) then
            deallocate(one_rdm%sym_list_inv_no, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem deallocating sym_list_inv_no array.')
        end if

    end subroutine dealloc_one_rdm_t

    pure subroutine clear_rdm_list_t(rdm)

        use hash, only: clear_hash_table

        type(rdm_list_t), intent(inout) :: rdm

        rdm%nelements = 0
        rdm%elements = 0_int_rdm
        call clear_hash_table(rdm%hash_table)

    end subroutine clear_rdm_list_t

    pure subroutine clear_one_rdms(one_rdms)

        ! Clear all the arrays for the one_rdm_t objects passed in (except for
        ! sym_list arrays, which we'd expect to be constant throughout a
        ! simulation, except for special cases like basis rotations).

        type(one_rdm_t), intent(inout) :: one_rdms(:)

        integer :: irdm

        do irdm = 1, size(one_rdms)
            if (allocated(one_rdms(irdm)%matrix)) one_rdms(irdm)%matrix = 0.0_dp
            if (allocated(one_rdms(irdm)%evalues)) one_rdms(irdm)%evalues = 0.0_dp
            if (allocated(one_rdms(irdm)%rho_ii)) one_rdms(irdm)%rho_ii = 0.0_dp
            if (allocated(one_rdms(irdm)%lagrangian)) one_rdms(irdm)%lagrangian = 0.0_dp
        end do

    end subroutine clear_one_rdms

    pure subroutine calc_combined_rdm_label(i, j, k, l, ijkl)

        ! Combine the four 2-RDM orbital labels into unique integers. i and j
        ! are combined into one number, ij. k and l are combined into one
        ! number, kl. Both of these are then combined into one single
        ! number, ijkl. Assuming (i,j,k,l) are *spin* orbitals labels (which
        ! they usually will be but not necessarily), the largest value for ijkl
        ! is M^4, where M is the number of spin orbitals.

        ! The compression defined in this routine will not give a fully
        ! compressed RDM index labelling, because it allows a separate ij
        ! integer if i and j are equal, even though this RDM element is never
        ! accessed, and the same for k and l. It also doesn't take spatial
        ! symmetry into account. But this is fine if one just seeks a unique
        ! combined label for each combination of individual orbital labels.

        ! In: i, j, k, l - orbital labels for the RDM contribution.
        ! Out: ijkl - Label combining i, j, k, l.

        use SystemData, only: nbasis

        integer, intent(in) :: i, j, k, l ! spin or spatial orbitals
        integer(int_rdm), intent(out) :: ijkl
        integer :: ij, kl

        ij = (i-1)*nbasis + j
        kl = (k-1)*nbasis + l
        ijkl = (ij-1)*(int(nbasis, int_rdm)**2) + kl

    end subroutine calc_combined_rdm_label

    pure subroutine calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

        ! Decode the four orbital labels stored in the input ijkl, i.e. do the
        ! opposite of calc_combined_rdm_label.

        use SystemData, only: nbasis

        integer(int_rdm), intent(in) :: ijkl
        integer, intent(out) :: ij, kl, i, j, k, l ! spin or spatial orbitals

        kl = int(mod(ijkl - 1, int(nbasis, int_rdm)**2)) + 1
        ij = int((ijkl - int(kl, int_rdm)) / int(nbasis, int_rdm)**2) + 1

        j = mod(ij-1, nbasis) + 1
        i = (ij - j)/nbasis + 1

        l = mod(kl-1, nbasis) + 1
        k = (kl - l)/nbasis + 1

    end subroutine calc_separate_rdm_labels

    pure subroutine extract_sign_rdm(rdm_entry, real_sign)

        ! Extract and decode the RDM sign stored in rdm_entry. This input
        ! entry has kind int_rdm, and contains the encoded spin orbital labels
        ! in the 0'th element. We want just the other elements, and in their
        ! real(dp) form.

        integer(int_rdm), intent(in) :: rdm_entry(0:)
        real(dp), intent(out) :: real_sign(size(rdm_entry)-1)
        integer(int_rdm) :: int_sign(size(rdm_entry)-1)

        int_sign = rdm_entry(1:size(int_sign))
        real_sign = transfer(int_sign, real_sign)

    end subroutine extract_sign_rdm

    pure subroutine encode_sign_rdm(rdm_entry, real_sign)

        ! Take the real RDM elements stored in real_sign, and encode them
        ! as integers with kind int_rdm, and then place these encoded signs in
        ! the positive elements of the rdm_entry array (the 0'th element
        ! contains the encoded spin orbital labels).

        integer(int_rdm), intent(inout) :: rdm_entry(0:)
        real(dp), intent(in) :: real_sign(1:size(rdm_entry)-1)
        integer(int_rdm) :: int_sign(1:size(rdm_entry)-1)

        int_sign = transfer(real_sign, int_sign)
        rdm_entry(1:size(int_sign)) = int_sign

    end subroutine encode_sign_rdm

    pure subroutine extract_sign_EN(sign_length, ilut, real_sign)

        integer, intent(in) :: sign_length
        integer(n_int), intent(in) :: ilut(0:)
        real(dp), intent(out) :: real_sign(sign_length)
        integer(n_int) :: sign(sign_length)

        sign = ilut(NIfDBO+1:NIfDBO+sign_length)
        real_sign = transfer(sign, real_sign)

    end subroutine extract_sign_EN

    pure subroutine encode_sign_EN(sign_length, ilut, real_sign)

        integer, intent(in) :: sign_length
        integer(n_int), intent(inout) :: ilut(0:)
        real(dp), intent(in) :: real_sign(sign_length)
        integer(n_int) :: sign(sign_length)

        sign = transfer(real_sign, sign)
        ilut(NIfDBO+1:NIfDBO+sign_length) = sign

    end subroutine encode_sign_EN

    subroutine add_to_rdm_spawn_t(spawn, i, j, k, l, contrib_sign, spinfree, nearly_full)

        ! In/Out: rdm_spawn - the rdm_spawn_t object to which contributions will be added.
        ! In: i, j, k, l - orbitals labels for the RDM contribution, with i<j, k<l.
        ! In: contrib_sign - the sign (amplitude) of the contribution to be added.
        ! In: spinfree - is the RDM being created to be output directly in spinfree form?
        ! In/Out: nearly_full - make this logical true if we come close to filling a
        !     processes section of the spawning list.

        use hash, only: hash_table_lookup, add_hash_table_entry
        use SystemData, only: nbasis

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: i, j, k, l ! spin or spatial orbital
        real(dp), intent(in) :: contrib_sign(spawn%rdm_send%sign_length)
        logical, intent(in) :: spinfree
        logical, intent(inout), optional :: nearly_full

        integer(int_rdm) :: ijkl
        integer :: ij_compressed, proc, ind, hash_val, slots_left
        real(dp) :: real_sign_old(spawn%rdm_send%sign_length), real_sign_new(spawn%rdm_send%sign_length)
        logical :: tSuccess
        character(*), parameter :: t_r = 'add_to_rdm_spawn_t'

        associate(rdm => spawn%rdm_send)

            ! Calculate combined RDM labels. A different ordering is used for
            ! outputting RDMs with and without spin. The following definitions
            ! will aid ordering via a sort operation later.
            if (.not. spinfree) then
                call calc_combined_rdm_label(i, j, k, l, ijkl)
            else
                call calc_combined_rdm_label(k, l, j, i, ijkl)
            end if

            ! Search to see if this RDM element is already in the RDM array.
            ! If it, tSuccess will be true and ind will hold the position of the
            ! entry in rdm%elements.
            call hash_table_lookup((/i,j,k,l/), (/ijkl/), 0, rdm%hash_table, rdm%elements, ind, hash_val, tSuccess)

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
                    ij_compressed = nbasis*(i-1) - i*(i-1)/2 + j - i
                    ! Calculate the process for the element.
                    proc = (ij_compressed-1)*nProcessors/spawn%nrows
                else
                    ! For spin-free case, we halve the number of labels. Also,
                    ! the last two labels are dominant in the ordering, so use
                    ! these instead, to allow writing out in the correct order.
                    ij_compressed = (nbasis/2)*(k-1) + l
                    proc = (ij_compressed-1)*nProcessors/(nbasis**2/4)
                end if

                ! Check that there is enough memory for the new spawned RDM entry.
                slots_left = spawn%init_free_slots(proc+1) - spawn%free_slots(proc)

                if (slots_left <= 0) call try_rdm_spawn_realloc(spawn, proc, spinfree)

                if (present(nearly_full)) then
                    ! 10 chosen somewhat arbitrarily, although there are times
                    ! when we call this routine 8 times in a row, so best not
                    ! to make any smaller.
                    if (slots_left <= 10) nearly_full = .true.
                end if

                rdm%elements(0, spawn%free_slots(proc)) = ijkl
                call encode_sign_rdm(rdm%elements(:, spawn%free_slots(proc)), contrib_sign)

                call add_hash_table_entry(rdm%hash_table, spawn%free_slots(proc), hash_val)

                spawn%free_slots(proc) = spawn%free_slots(proc) + 1
            end if

        end associate

    end subroutine add_to_rdm_spawn_t

    subroutine communicate_rdm_spawn_t(spawn, rdm_recv)

        ! Perform communication of RDM elements in the spawn object, to the
        ! rdm_recv object. The hash table and free_slots array for the spawn
        ! object are then reset at the end of this routine. The newly
        ! received spawnings will be added to rdm_recv *without* overwriting
        ! elements currently in the list (which is useful for situations where
        ! multiple communications are required, due to limited space in the
        ! spawning array).

        use hash, only: clear_hash_table
        use Parallel_neci, only: MPIAlltoAll, MPIAlltoAllv

        type(rdm_spawn_t), intent(inout) :: spawn
        type(rdm_list_t), intent(inout) :: rdm_recv

        integer :: iproc, nelements_old, new_nelements, ierr, i
        integer(MPIArg) :: send_sizes(0:nProcessors-1), recv_sizes(0:nProcessors-1)
        integer(MPIArg) :: send_displs(0:nProcessors-1), recv_displs(0:nProcessors-1)

        nelements_old = rdm_recv%nelements

        ! How many rows of data to send to each processor.
        do iproc = 0, nProcessors-1
            send_sizes(iproc) = int(spawn%free_slots(iproc) - spawn%init_free_slots(iproc), MPIArg)
            ! The displacement of the beginning of each processor's section of the
            ! free_slots array, relative to the first element of this array.
            send_displs(iproc) = int(spawn%init_free_slots(iproc)-1,MPIArg)
        end do

        ! this does not work with some compilers
        !send_displs = int(spawn%init_free_slots(0:nProcessors-1) - 1, MPIArg)

        call MPIAlltoAll(send_sizes, 1, recv_sizes, 1, ierr)

        recv_displs(0) = 0
        do iproc = 1, nProcessors-1
            recv_displs(iproc) = recv_displs(iproc-1) + recv_sizes(iproc-1)
        end do

        ! The total number of RDM elements in the list after the receive.
        new_nelements = rdm_recv%nelements + int(sum(recv_sizes), sizeof_int)

        ! If we don't have enough memory in the receiving list, try
        ! reallocating it to be big enough.
        if (new_nelements > rdm_recv%max_nelements) then
            call try_rdm_list_realloc(rdm_recv, new_nelements, .true.)
        end if

        ! Update the number of valid RDM elements in the received list.
        rdm_recv%nelements = new_nelements

        send_sizes = send_sizes*size(spawn%rdm_send%elements,1)
        recv_sizes = recv_sizes*size(spawn%rdm_send%elements,1)
        send_displs = send_displs*size(spawn%rdm_send%elements,1)
        recv_displs = recv_displs*size(spawn%rdm_send%elements,1)

        rdm_recv%elements(:,nelements_old+1:) = 0
        ! Perform the communication.
        call MPIAlltoAllv(spawn%rdm_send%elements, send_sizes, send_displs, &
                          rdm_recv%elements(:,nelements_old+1:), &
                          recv_sizes, recv_displs, ierr)

        ! Now we can reset the free_slots array and reset the hash table.
        do iproc = 0, nProcessors - 1
           spawn%free_slots(iproc) = spawn%init_free_slots(iproc)
        end do
        call clear_hash_table(spawn%rdm_send%hash_table)

    end subroutine communicate_rdm_spawn_t

    subroutine communicate_rdm_spawn_t_wrapper(spawn, rdm_recv, finished, all_finished)

        ! This is a wrapper function around communicate_rdm_spawn_t, which is
        ! useful for certain routines in rdm_finalising and rdm_estimators.
        ! Upon calling this routine, each process will input the logical
        ! finished as either .true. or .false., with the former indicating
        ! that this will be the final call to communicate_rdm_spawn_t by the
        ! routine. This routine then checks if ever routine is on its final
        ! required call, and if so, the calling functions can know to not
        ! perform any further communications.

        use Parallel_neci, only: MPIAllGather

        type(rdm_spawn_t), intent(inout) :: spawn
        type(rdm_list_t), optional, intent(inout) :: rdm_recv
        logical, intent(in) :: finished
        logical, intent(inout) :: all_finished

        logical :: finished_array(nProcessors)
        integer :: ierr

        call communicate_rdm_spawn_t(spawn, rdm_recv)

        ! Find if all processes have finished their communication.
        call MPIAllGather(finished, finished_array, ierr)
        all_finished = all(finished_array)

    end subroutine communicate_rdm_spawn_t_wrapper

    subroutine try_rdm_list_realloc(rdm_recv, new_nelements, recv_list)

        ! For cases where the receiving RDM array is not big enough for a
        ! communication, try and reallocate it to be big enough. This also
        ! requires a temporary array to be allocated, to store the current
        ! state of the receive list.

        ! recv_list should be input as true if reallocating an receiving RDM
        ! object. It should be false if reallocating the main array in the
        ! subroutine add_rdm_1_to_rdm_2. The only difference this makes is
        ! in the message output.

        type(rdm_list_t), intent(inout) :: rdm_recv
        integer, intent(in) :: new_nelements
        logical, intent(in) :: recv_list

        integer :: old_nelements, memory_old, memory_new, ierr
        integer(int_rdm), allocatable :: temp_elements(:,:)
        character(*), parameter :: t_r = 'try_rdm_list_realloc'

        ! The number of elements currently filled in the RDM array.
        old_nelements = rdm_recv%nelements

        if (recv_list) then
            write(6,'("WARNING: There is not enough space in the current RDM array to receive all of the &
                      &communicated RDM elements. We will now try and reallocate this array to be large &
                      &enough. If there is not sufficient memory then the program may crash.")'); call neci_flush(6)
        else
            write(6,'("WARNING: There is not enough space in the current RDM array to add the received &
                      &RDM elements to the main RDM array. We will now try and reallocate this array to be 1.5 &
                      &times larger. If there is not sufficient memory then the program may crash.")'); call neci_flush(6)
        end if

        ! Memory of the old and new arrays, in bytes.
        memory_old = rdm_recv%max_nelements*(rdm_recv%sign_length+1)*size_int_rdm
        memory_new = new_nelements*(rdm_recv%sign_length+1)*size_int_rdm

        write(6,'("Old RDM array had the following size (MB):", f14.6)') real(memory_old,dp)/1048576.0_dp
        write(6,'("Required new RDM array must have the following size (MB):", f14.6)') real(memory_new,dp)/1048576.0_dp

        if (old_nelements > 0) then
            ! Allocate a temporary array to copy the old RDM list to, while we
            ! reallocate that array.
            allocate(temp_elements(0:rdm_recv%sign_length, old_nelements), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error while allocating temporary array to hold existing &
                                              &RDM receive array.")
            temp_elements = rdm_recv%elements(:,1:old_nelements)
        end if

        deallocate(rdm_recv%elements, stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error while deallocating existing RDM receive array.")

        allocate(rdm_recv%elements(0:rdm_recv%sign_length, new_nelements), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error while allocating RDM receive array to the new larger size.")

        ! Update the maximum number of elements for the rdm_recv object.
        rdm_recv%max_nelements = new_nelements

        if (old_nelements > 0) then
            ! Copy the existing elements back, and deallocate the temorary array.
            rdm_recv%elements(:,1:old_nelements) = temp_elements

            deallocate(temp_elements, stat=ierr)
            if (ierr /= 0) call stop_all(t_r, "Error while deallocating temporary RDM array.")
        end if

    end subroutine try_rdm_list_realloc

    subroutine try_rdm_spawn_realloc(spawn, proc, spinfree)

        use hash, only: update_hash_table_ind

        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: proc
        logical, intent(in) :: spinfree

        real(dp) :: slots_per_proc_new
        integer :: old_max_length, new_max_length, nstored
        integer :: memory_old, memory_new, pos_diff, iproc, ierr
        integer :: new_init_slots(0:nProcessors)
        integer :: i, j, k, l, ij, kl, idet, hash_val
        integer(int_rdm), allocatable :: temp_elements(:,:), ijkl
        character(*), parameter :: t_r = 'try_rdm_spawn_realloc'

        associate(rdm => spawn%rdm_send)

        old_max_length = rdm%max_nelements
        new_max_length = 2*rdm%max_nelements

        slots_per_proc_new = real(new_max_length, dp)/real(nProcessors, dp)

        ! Create new init_free_slots array.
        do iproc = 0, nProcessors-1
            new_init_slots(iproc) = nint(slots_per_proc_new*iproc)+1
        end do
        new_init_slots(nProcessors) = new_max_length + 1

        write(6,'("WARNING: There is not enough space in the current RDM spawning array to store the &
                  &RDM elements to be sent to process",'//int_fmt(proc,1)//',". We will now try and &
                  &reallocate the entire RDM spawning array to be twice its current size. If there is &
                  &not sufficient memory then the program may crash. This also requires recreating the &
                  &hash table to some of this object, which may take some time.")') proc; call neci_flush(6)

        ! Memory of the old and new arrays, in bytes.
        memory_old = old_max_length*(rdm%sign_length+1)*size_int_rdm
        memory_new = new_max_length*(rdm%sign_length+1)*size_int_rdm

        write(6,'("Old RDM spawning array had the following size (MB):", f14.6)') real(memory_old,dp)/1048576.0_dp
        write(6,'("Required new array must have the following size (MB):", f14.6)') real(memory_new,dp)/1048576.0_dp

        ! Allocate a temporary array to copy the old RDM list to, while we
        ! reallocate that array.
        allocate(temp_elements(0:rdm%sign_length, old_max_length), stat=ierr)

        if (ierr /= 0) call stop_all(t_r, "Error while allocating temporary array to hold existing &
                                          &RDM spawning array.")
        temp_elements(:, 1:old_max_length) = rdm%elements

        deallocate(rdm%elements, stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error while deallocating existing RDM spawning array.")

        allocate(rdm%elements(0:rdm%sign_length, new_max_length), stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error while allocating RDM spawning array to the new larger size.")
        ! Update the maximum number of elements for the spawning array.
        rdm%max_nelements = new_max_length

        ! Loop over all processes, copy RDM spawns into the new list in the
        ! correct new positons, and update the hash table as necessary.
        do iproc = 0, nProcessors-1
            ! The number of RDM elements actually filled in in this processor's
            ! section of the spawning array.
            nstored = spawn%free_slots(iproc) - spawn%init_free_slots(iproc)

            ! Copy RDM elements back across from the temporary array.
            rdm%elements(:, new_init_slots(iproc):new_init_slots(iproc)+nstored-1) = &
                temp_elements(:, spawn%init_free_slots(iproc):spawn%init_free_slots(iproc)+nstored-1)

            ! Update the free_slots array as necessary.
            spawn%free_slots(iproc) = new_init_slots(iproc) + nstored

            ! How much the beginning of the sections for this processor have changed,
            pos_diff = new_init_slots(iproc) - spawn%init_free_slots(iproc)

            ! Update hash table. Don't need to update proc 1 section, since this
            ! does not get moved.
            if (iproc > 0) then
                if (spinfree) then
                    do idet = new_init_slots(iproc), new_init_slots(iproc)+nstored-1
                        call calc_separate_rdm_labels(rdm%elements(0,idet), ij, kl, k, l, j, i)
                        call update_hash_table_ind(rdm%hash_table, (/i,j,k,l/), idet-pos_diff, idet)
                    end do
                else
                    do idet = new_init_slots(iproc), new_init_slots(iproc)+nstored-1
                        call calc_separate_rdm_labels(rdm%elements(0,idet), ij, kl, i, j, k, l)
                        call update_hash_table_ind(rdm%hash_table, (/i,j,k,l/), idet-pos_diff, idet)
                    end do
                end if
            end if
        end do

        spawn%init_free_slots = new_init_slots

        deallocate(temp_elements, stat=ierr)
        if (ierr /= 0) call stop_all(t_r, "Error while deallocating temporary array.")

        end associate

    end subroutine try_rdm_spawn_realloc

    subroutine add_rdm_1_to_rdm_2(rdm_1, rdm_2, scale_factor)

        ! Take the RDM elements in the rdm_1 object, and add them to the rdm_2
        ! object. The has table for rdm_2 will also be updated. This literally
        ! performs the numerical addition of the two RDM objects.
        use SystemData, only: nel, nBasis
        use hash, only: hash_table_lookup, add_hash_table_entry
        use Parallel_neci, only: MPISumAll

        type(rdm_list_t), intent(in) :: rdm_1
        type(rdm_list_t), intent(inout) :: rdm_2
        real(dp), intent(in), optional :: scale_factor

        integer :: ielem
        integer(int_rdm) :: ijkl
        integer :: ij, kl, i, j, k, l ! spin or spatial orbitals
        integer :: ind, hash_val
        real(dp) :: real_sign_old(rdm_2%sign_length), real_sign_new(rdm_2%sign_length)
        real(dp) :: spawn_sign(rdm_2%sign_length)
        real(dp) :: internal_scale_factor(rdm_1%sign_length), rdm_trace(rdm_1%sign_length)
        logical :: tSuccess
        character(*), parameter :: t_r = 'add_rdm_1_to_rdm_2'

        if(present(scale_factor)) then
           ! normalize and rescale the rdm_1 if requested here
           call calc_rdm_trace(rdm_1, rdm_trace)
           call MPISumAll(rdm_trace, internal_scale_factor)
           internal_scale_factor = scale_factor*(nel*(nel-1))/(2*internal_scale_factor)
        else
           internal_scale_factor = 1.0_dp
        endif

        if(rdm_1%sign_length .ne. rdm_2%sign_length) call stop_all(t_r,"nrdms mismatch")

        do ielem = 1, rdm_1%nelements
            ! Decode the compressed RDM labels.
            ijkl = rdm_1%elements(0,ielem)
            call calc_separate_rdm_labels(ijkl, ij, kl, i, j, k, l)

            ! Extract the spawned sign.
            call extract_sign_rdm(rdm_1%elements(:,ielem), spawn_sign)

            ! Search to see if this RDM element is already in the RDM 2.
            ! If it, tSuccess will be true and ind will hold the position of the
            ! element in rdm.
            if(any((/i,j,k,l/)<=0) .or. any((/i,j,k,l/)>nBasis)) then
               write(iout,*) "Invalid rdm element", i,j,k,l,ijkl, ielem, rdm_1%nelements
               call stop_all(t_r,"Erroneous indices")
            endif
            call hash_table_lookup((/i,j,k,l/), (/ijkl/), 0, rdm_2%hash_table, rdm_2%elements, ind, hash_val, tSuccess)

            if (tSuccess) then
                ! Extract the existing sign.
                call extract_sign_rdm(rdm_2%elements(:,ind), real_sign_old)

                ! Update the total sign.
                real_sign_new = real_sign_old + spawn_sign*internal_scale_factor
                ! Encode the new sign.
                call encode_sign_rdm(rdm_2%elements(:,ind), real_sign_new)
            else
                ! If we don't have enough memory in rdm_2, try increasing its
                ! size to be 1.5 times bigger.
                if (rdm_2%nelements+1 > rdm_2%max_nelements) then
                    call try_rdm_list_realloc(rdm_2, int(1.5*rdm_2%max_nelements), .false.)
                end if

                ! Update the rdm array, and its hash table, and the number of
                ! RDM elements.
                rdm_2%nelements = rdm_2%nelements + 1
                rdm_2%elements(0, rdm_2%nelements) = ijkl
                call encode_sign_rdm(rdm_2%elements(:, rdm_2%nelements), spawn_sign)
                call add_hash_table_entry(rdm_2%hash_table, rdm_2%nelements, hash_val)
            end if
        end do

    end subroutine add_rdm_1_to_rdm_2

    subroutine scale_rdm(rdm, scale_factor)
      ! rescale all entries of one rdm by a scale factor
      ! required in the adaptive shift correction, where a weighted sum
      ! of two rdms is taken
      implicit none
      type(rdm_list_t), intent(inout) :: rdm
      real(dp) :: scale_factor(rdm%sign_length)

      integer :: i, j
      real(dp) :: tmp_sign(rdm%sign_length)

      do i = 1, rdm%nelements
         call extract_sign_rdm(rdm%elements(:,i), tmp_sign)
         tmp_sign = tmp_sign * scale_factor
         call encode_sign_rdm(rdm%elements(:,i), tmp_sign)
      end do
    end subroutine scale_rdm

    subroutine annihilate_rdm_list(rdm)

        ! Perform annihilation on RDM elements in rdm%elements.

        ! * WARNING * rdm%hash_table is not used by this routine, and will not
        ! be updated by it after reordering and compressing rdm%elements.

        use sort_mod, only: sort

        type(rdm_list_t), intent(inout) :: rdm

        integer :: ielem, counter
        real(dp) :: rdm_sign(rdm%sign_length), summed_rdm_sign(rdm%sign_length)

        if (rdm%nelements > 0) then
            call sort(rdm%elements(:,1:rdm%nelements))

            summed_rdm_sign = 0.0_dp
            counter = 0

            do ielem = 1, rdm%nelements-1
                call extract_sign_rdm(rdm%elements(:,ielem), rdm_sign)
                summed_rdm_sign = summed_rdm_sign + rdm_sign

                ! Is the next RDM element the same as this one? If so then
                ! don't keep this element yet, but wait until all signs on
                ! this RDM element have been summed (annihilated).
                if (.not. (rdm%elements(0,ielem) == rdm%elements(0, ielem+1))) then
                    ! If this element is zero for all RDMs, then don't keep it.
                    if (any(abs(summed_rdm_sign) > 1.0e-12_dp)) then
                        counter = counter + 1
                        rdm%elements(0, counter) = rdm%elements(0,ielem)
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

    subroutine add_to_en_pert_t(en_pert, nI, ilut, contrib_sign)

        ! In/Out: en_pert - the en_pert_t object to which contributions will be added.
        ! In: nI - A list of the occupied orbitals in the determinant.
        ! In: ilut - The determinant in a bitstring form.
        ! In: contrib_sign - the sign (amplitude) of the contribution to be added.

        use hash, only: hash_table_lookup, add_hash_table_entry
        use SystemData, only: nel

        type(en_pert_t), intent(inout) :: en_pert
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), intent(in) :: contrib_sign(en_pert%sign_length)

        integer :: ind, hash_val, slots_left
        real(dp) :: real_sign_old(en_pert%sign_length), real_sign_new(en_pert%sign_length)
        logical :: tSuccess
        character(*), parameter :: t_r = 'add_to_en_pert_t'

        ! Search to see if this determinant is already in the dets array.
        ! If it, tSuccess will be true and ind will hold the position of the
        ! entry in en_pert%dets.
        call hash_table_lookup(nI, ilut, NIfDBO, en_pert%hash_table, en_pert%dets, ind, hash_val, tSuccess)

        if (tSuccess) then
            ! Extract the existing sign.
            call extract_sign_EN(en_pert%sign_length, en_pert%dets(:,ind), real_sign_old)
            ! Update the total sign.
            real_sign_new = real_sign_old + contrib_sign
            ! Encode the new sign.
            call encode_sign_EN(en_pert%sign_length, en_pert%dets(:,ind), real_sign_new)
        else
            en_pert%ndets = en_pert%ndets + 1

            ! Check that there is enough memory for the new determinant.
            slots_left = en_pert%max_ndets - en_pert%ndets

            if (slots_left < 0) then
                write(6,'("ERROR: No space left in the EN2 array. Aborting to prevent incorrect results...")')
                call neci_flush(6)
                call stop_all(t_r, 'No space left in the EN2 array. Please increase memoryfacspawn.')
            else if (slots_left < 20) then
                write(6,'("WARNING: Less than 20 slots left in EN2 array. The program will abort &
                           &when there are no slots remaining.")'); call neci_flush(6)
            end if

            en_pert%dets(0:NIfDBO, en_pert%ndets) = ilut(0:NIfDBO)
            call encode_sign_EN(en_pert%sign_length, en_pert%dets(:, en_pert%ndets), contrib_sign)

            call add_hash_table_entry(en_pert%hash_table, en_pert%ndets, hash_val)
        end if

    end subroutine add_to_en_pert_t

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

end module rdm_data_utils
