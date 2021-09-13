#include "macros.h"

! This module contains a type and routines for defining and creating a sparse Hamiltonian.
! The type created to store this information is sparse_matrix_real. For an N-by-N matrix,
! one creates a 1d array of type sparse_matrix_real. Each element of this array stores the
! information on one single row of the matrix - the positions of the non-zero elements
! and their values, in the same order. This completley defines the matrix.

! Note that in some parts of NECI, the sparse nature is used in a different form (i.e. in
! the Lanczos code).

module sparse_arrays

    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: decode_bit_det, nifguga
    use CalcData, only: tReadPops
    use constants
    use DetBitOps, only: DetBitEq, CountBits, TestClosedShellDet
    use Determinants, only: get_helement
    use FciMCData, only: SpawnedParts, Hii
    use core_space_util, only: cs_replicas, sparse_matrix_real, sparse_matrix_int, &
                               core_space_t
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement, &
                              hphf_off_diag_helement_opt
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier, MPIAllGatherV
    use SystemData, only: tHPHF, nel
    use global_det_data, only: set_det_diagH
    use shared_rhash, only: shared_rhash_t
    use SystemData, only: tGUGA
    use guga_excitations, only: actHamiltonian, &
                                calc_guga_matrix_element
    use guga_bitRepOps, only: convert_ilut_toGUGA, extract_h_element, &
                              init_csf_information
    use util_mod, only: binary_search, near_zero
    use guga_data, only: tag_excitations, ExcitationInformation_t
    use guga_matrixElements, only: calcDiagMatEleGuga_nI

    implicit none

    type trial_hashtable
        ! All the states with this hash value.
        integer(n_int), allocatable, dimension(:, :) :: states
        ! The number of clashes for ths hash value.
        integer :: nclash
    end type trial_hashtable

    type core_hashtable
        ! The indices of states with this hash table.
        integer, allocatable, dimension(:) :: ind
        ! The number of clashes for this hash value.
        integer :: nclash
    end type core_hashtable

    type(sparse_matrix_real), allocatable, dimension(:) :: sparse_ham
    integer(TagIntType), allocatable, dimension(:, :) :: SparseHamilTags

    ! For quick access it is often useful to have just the diagonal elements. Note,
    ! however, that they *are* stored in sparse_ham too.
    real(dp), allocatable, dimension(:) :: hamil_diag
    integer(TagIntType) :: HDiagTag

    type(trial_hashtable), allocatable, dimension(:) :: trial_ht
    type(trial_hashtable), allocatable, dimension(:) :: con_ht

    ! --- For when using the determ-proj-approx-hamil option -----
    type(shared_rhash_t) :: var_ht
    type(sparse_matrix_real), allocatable, dimension(:) :: approx_ham

contains

    subroutine calculate_sparse_hamiltonian_non_hermitian(num_states, ilut_list)
        ! same routine as below, but for non-hermitian Hamiltonians in the
        ! case of transcorrelation
        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character(len=*), parameter :: t_r = "calculate_sparse_hamiltonian_non_hermitian"

        allocate(sparse_ham(num_states))
        allocate(SparseHamilTags(2, num_states))
        allocate(hamiltonian_row(num_states), stat=ierr)
        call LogMemAlloc('hamiltonian_row', num_states, 8, t_r, HRTag, ierr)
        allocate(hamil_diag(num_states), stat=ierr)
        call LogMemAlloc('hamil_diag', num_states, 8, t_r, HDiagTag, ierr)
        allocate(sparse_row_sizes(num_states), stat=ierr)
        call LogMemAlloc('sparse_row_sizes', num_states, bytes_int, t_r, SRTag, ierr)

        ! Set each element to one to count the diagonal elements straight away.
        sparse_row_sizes = 1

        if (tGUGA) then
            call stop_all(t_r, "modify get_helement for GUGA!")
        end if

        do i = 1, num_states

            hamiltonian_row = 0.0_dp

            call decode_bit_det(nI, ilut_list(:, i))

            ! we have to loop over everything in case on non-hermiticity
            do j = 1, num_states

                call decode_bit_det(nJ, ilut_list(:, j))
                if (i == j) then
                    if (tHPHF) then
                        hamiltonian_row(i) = hphf_diag_helement(nI, ilut_list(:, i))
                    else
                        hamiltonian_row(i) = get_helement(nI, nI, 0)
                    end if
                    hamil_diag(i) = hamiltonian_row(i)
                else
                    if (tHPHF) then
                        !TODO: do i need <I|H|J> or <J|H|I>?
                        hamiltonian_row(j) = hphf_off_diag_helement( &
                                             nI, nJ, ilut_list(:, i), ilut_list(:, j))
                    else
                        hamiltonian_row(j) = get_helement( &
                                             nI, nJ, ilut_list(:, i), ilut_list(:, j))
                    end if
                    if (abs(hamiltonian_row(j)) > EPS) then
                        ! i think in the non-hermitian i only need to update
                        ! one of the counters..
                        sparse_row_sizes(i) = sparse_row_sizes(i) + 1
                    end if
                end if
            end do

            call allocate_sparse_ham_row(sparse_ham, i, sparse_row_sizes(i), &
                                         "sparse_ham", SparseHamilTags(:, i))

            sparse_ham(i)%elements = 0.0_dp
            sparse_ham(i)%positions = 0
            sparse_ham(i)%num_elements = sparse_row_sizes(i)

            ! now fill in the elements, all of them
            counter = 1
            do j = 1, num_states
                if (abs(hamiltonian_row(j)) > EPS) then
                    sparse_ham(i)%positions(counter) = j
                    sparse_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
            end do

        end do

        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(sparse_row_sizes, stat=ierr)
        call LogMemDealloc(t_r, SRTag, ierr)

    end subroutine calculate_sparse_hamiltonian_non_hermitian

    subroutine calculate_sparse_hamiltonian(num_states, ilut_list)

        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character(len=*), parameter :: t_r = "calculate_sparse_hamiltonian"

        integer :: pos, nexcits
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:, :)
        type(ExcitationInformation_t) :: excitInfo

        allocate(sparse_ham(num_states))
        allocate(SparseHamilTags(2, num_states))
        allocate(hamiltonian_row(num_states), stat=ierr)
        call LogMemAlloc('hamiltonian_row', num_states, 8, t_r, HRTag, ierr)
        allocate(hamil_diag(num_states), stat=ierr)
        call LogMemAlloc('hamil_diag', num_states, 8, t_r, HDiagTag, ierr)
        allocate(sparse_row_sizes(num_states), stat=ierr)
        call LogMemAlloc('sparse_row_sizes', num_states, bytes_int, t_r, SRTag, ierr)
        allocate(sparse_diag_positions(num_states), stat=ierr)
        call LogMemAlloc('sparse_diag_positions', num_states, bytes_int, t_r, SDTag, ierr)
        allocate(indices(num_states), stat=ierr)
        call LogMemAlloc('indices', num_states, bytes_int, t_r, ITag, ierr)

        ! Set each element to one to count the diagonal elements straight away.
        sparse_row_sizes = 1

        ! also need to change this routine here for the guga implementation

        do i = 1, num_states

            hamiltonian_row = 0.0_dp

            call decode_bit_det(nI, ilut_list(:, i))

            ! sparse_diag_positions(i) stores the number of non-zero elements in row i of the
            ! Hamiltonian, up to and including the diagonal element.
            ! sparse_row_sizes(i) stores this number currently as all non-zero elements before
            ! the diagonal have been counted (as the Hamiltonian is symmetric).
            sparse_diag_positions(i) = sparse_row_sizes(i)

            if (tGUGA) call init_csf_information(ilut_list(0:nifd, i))

            do j = i, num_states

                call decode_bit_det(nJ, ilut_list(:, j))

                ! If on the diagonal of the Hamiltonian.
                if (i == j) then
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_diag_helement(nI, ilut_list(:, i))
                    else if (tGUGA) then
                        hamiltonian_row(j) = calcDiagMatEleGuga_nI(nI)
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, 0)
                    end if
                    hamil_diag(j) = hamiltonian_row(j)
                else
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, ilut_list(:, i), &
                                                                    ilut_list(:, j))
                    else if (tGUGA) then
                        call calc_guga_matrix_element(ilut_list(:, i), ilut_list(:, j), &
                                                      excitInfo, hamiltonian_row(j), .true., 1)
#ifdef CMPLX_
                        hamiltonian_row(j) = conjg(hamiltonian_row(j))
#endif
                        ! call calc_guga_matrix_element(ilut_list(:,j), ilut_list(:,i), &
                        !         excitInfo, hamiltonian_row(j), .true., 2)
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, ilut_list(:, i), &
                                                          ilut_list(:, j))
                    end if
                    if (abs(hamiltonian_row(j)) > 0.0_dp) then
                        ! If element is nonzero, update the following sizes.
                        sparse_row_sizes(i) = sparse_row_sizes(i) + 1
                        sparse_row_sizes(j) = sparse_row_sizes(j) + 1
                    end if
                end if
            end do
            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_ham, i, sparse_row_sizes(i), "sparse_ham", SparseHamilTags(:, i))

            sparse_ham(i)%elements = 0.0_dp
            sparse_ham(i)%positions = 0
            sparse_ham(i)%num_elements = sparse_row_sizes(i)

            ! Now fill in all matrix elements beyond and including the diagonal, as these are
            ! stored in hamiltonian_row.
            counter = sparse_diag_positions(i)
            do j = i, num_states
                if (abs(hamiltonian_row(j)) > 0.0_dp) then
                    sparse_ham(i)%positions(counter) = j
                    sparse_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
            end do
        end do

        ! At this point, sparse_ham has been allocated with the correct sizes, but only the
        ! matrix elements above and including the diagonal have been filled in. Now we must
        ! fill in the other elements. To do this, cycle through every element above the
        ! diagonal and fill in every corresponding element below it:
        indices = 1
        do i = 1, num_states
            do j = sparse_diag_positions(i) + 1, sparse_row_sizes(i)

                sparse_ham(sparse_ham(i)%positions(j))%&
                    &positions(indices(sparse_ham(i)%positions(j))) = i
                sparse_ham(sparse_ham(i)%positions(j))%&
                    &elements(indices(sparse_ham(i)%positions(j))) = h_conjg(sparse_ham(i)%elements(j))

                indices(sparse_ham(i)%positions(j)) = indices(sparse_ham(i)%positions(j)) + 1

            end do
        end do

        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(sparse_row_sizes, stat=ierr)
        call LogMemDealloc(t_r, SRTag, ierr)
        deallocate(sparse_diag_positions, stat=ierr)
        call LogMemDealloc(t_r, SDTag, ierr)
        deallocate(indices, stat=ierr)
        call LogMemDealloc(t_r, ITag, ierr)

    end subroutine calculate_sparse_hamiltonian

    subroutine calculate_sparse_ham_par(num_states, ilut_list, tPrintInfo)

        integer(MPIArg), intent(in) :: num_states(0:nProcessors - 1)
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states(iProcIndex))
        logical, intent(in) :: tPrintInfo
        integer(MPIArg) :: disps(0:nProcessors - 1)
        integer :: i, j, row_size, counter, num_states_tot, ierr, bytes_required
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:, :) :: temp_store
        integer(TagIntType) :: TempStoreTag, HRTag, SDTag
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calculate_sparse_ham_par"

        integer :: pos, nexcits
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:, :)
        type(ExcitationInformation_t) :: excitInfo

        num_states_tot = int(sum(num_states), sizeof_int)
        disps(0) = 0
        do i = 1, nProcessors - 1
            disps(i) = disps(i - 1) + num_states(i - 1)
        end do

        safe_realloc_e(sparse_ham, (num_states(iProcIndex)), ierr)
        safe_realloc_e(SparseHamilTags, (2, num_states(iProcIndex)), ierr)
        safe_realloc_e(hamiltonian_row, (num_states_tot), ierr)
        call LogMemAlloc('hamiltonian_row', num_states_tot, 8, t_r, HRTag, ierr)
        safe_realloc_e(hamil_diag, (num_states(iProcIndex)), ierr)
        call LogMemAlloc('hamil_diag', int(num_states(iProcIndex), sizeof_int), 8, t_r, HDiagTag, ierr)
        safe_realloc_e(temp_store, (0:NIfTot, num_states_tot), ierr)
        call LogMemAlloc('temp_store', num_states_tot * (NIfTot + 1), 8, t_r, TempStoreTag, ierr)

        ! Stick together the determinants from all processors, on all processors.
        call MPIAllGatherV(ilut_list(:, 1:num_states(iProcIndex)), temp_store, num_states, disps)

        ! Loop over all determinants on this processor.
        do i = 1, num_states(iProcIndex)

            call decode_bit_det(nI, ilut_list(:, i))

            row_size = 0
            hamiltonian_row = 0.0_dp
            ! Loop over all determinants on all processors.

            if (tGUGA) call init_csf_information(ilut_list(0:nifd, i))

            do j = 1, num_states_tot

                call decode_bit_det(nJ, temp_store(:, j))

                ! If on the diagonal of the Hamiltonian.
                if (DetBitEq(ilut_list(:, i), temp_store(:, j), nifd)) then
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_diag_helement(nI, ilut_list(:, i))
                    else if (tGUGA) then
                        hamiltonian_row(j) = calcDiagMatEleGuga_nI(nI)
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, 0)
                    end if
                    hamil_diag(i) = hamiltonian_row(j)
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, ilut_list(:, i), temp_store(:, j))
                    else if (tGUGA) then
                        call calc_guga_matrix_element(ilut_list(:, i), temp_store(:, j), &
                                                      excitInfo, hamiltonian_row(j), .true., 1)
#ifdef CMPLX_
                        hamiltonian_row(j) = conjg(hamiltonian_row(j))
#endif
                        ! call calc_guga_matrix_element(temp_store(:,j), ilut_list(:,i), &
                        !     excitInfo, hamiltonian_row(j), .true., 2)
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, ilut_list(:, i), temp_store(:, j))
                    end if
                    if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                end if
            end do

            if (tPrintInfo) then
                if (i == 1) then
                    bytes_required = row_size * (8 + bytes_int)
                    write(stdout, '(1x,a43)') "About to allocate first row of Hamiltonian."
                    write(stdout, '(1x,a40,1x,i8)') "The memory (bytes) required for this is:", bytes_required
                    write(stdout, '(1x,a71,1x,i7)') "The total number of determinants (and hence rows) on this processor is:", &
                        num_states(iProcIndex)
                    write(stdout, '(1x,a58,1x,i7)') "The total number of determinants across all processors is:", num_states_tot
                    write(stdout, '(1x,a77,1x,i7)') "It is therefore expected that the total memory (MB) required will be roughly:", &
                        num_states_tot * bytes_required / 1000000
                else if (mod(i, 1000) == 0) then
                    write(stdout, '(1x,a23,1x,i7)') "Finished computing row:", i
                end if
            end if

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_ham, i, row_size, "sparse_ham", SparseHamilTags(:, i))

            sparse_ham(i)%elements = 0.0_dp
            sparse_ham(i)%positions = 0
            sparse_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, num_states_tot
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + disps(iProcIndex))) then
                    sparse_ham(i)%positions(counter) = j
                    sparse_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        call MPIBarrier(ierr)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)
        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)

    end subroutine calculate_sparse_ham_par

    subroutine calc_determ_hamil_sparse(rep)

        use SystemData, only: t_3_body_excits, t_mol_3_body, t_ueg_transcorr
        type(core_space_t), intent(inout) :: rep
        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:, :) :: temp_store
        integer(TagIntType) :: HRTag, TempStoreTag
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row

        integer :: pos, nExcit
        integer(n_int) :: ilutG(0:nifguga)
        integer(n_int), pointer :: excitations(:, :)
        type(ExcitationInformation_t) :: excitInfo

        character(len=*), parameter :: this_routine = "calc_determ_hamil_sparse"

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        HElement_t(dp) :: tmp_mat, tmp_mat_2

        allocate(rep%sparse_core_ham(rep%determ_sizes(iProcIndex)), stat=ierr)
        allocate(rep%SparseCoreHamilTags(2, rep%determ_sizes(iProcIndex)))
        allocate(hamiltonian_row(rep%determ_space_size), stat=ierr)
        call LogMemAlloc('hamiltonian_row', int(rep%determ_space_size, sizeof_int), 8, this_routine, HRTag, ierr)
        allocate(rep%core_ham_diag(rep%determ_sizes(iProcIndex)), stat=ierr)
        allocate(temp_store(0:NIfTot, rep%determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', rep%determ_space_size * (NIfTot + 1), 8, this_routine, TempStoreTag, ierr)

        ! Stick together the deterministic states from all processors, on
        ! all processors.
        ! n.b. Explicitly use 0:NIfTot, as NIfTot may not equal NIfBCast
        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:rep%determ_sizes(iProcIndex)), &
                           temp_store(0:niftot, 1:), rep%determ_sizes, rep%determ_displs)


        ! Loop over all deterministic states on this processor.
        do i = 1, rep%determ_sizes(iProcIndex)

            ilutI = SpawnedParts(0:niftot, i)
            call decode_bit_det(nI, IlutI)

            row_size = 0
            hamiltonian_row = 0.0_dp

            if (tGUGA) call init_csf_information(ilutI(0:nifd))

            ! Loop over all deterministic states.
            do j = 1, rep%determ_space_size

                ilutJ = temp_store(:, j)
                call decode_bit_det(nJ, ilutJ)

                ! If on the diagonal of the Hamiltonian.
                if (DetBitEq(IlutI, ilutJ, nifd)) then
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_diag_helement(nI, IlutI) - Hii
                    else
                        ! for guga: the diagonal is fine, since i overwrite
                        ! that within get_helement
                        hamiltonian_row(j) = get_helement(nI, nJ, 0) - Hii
                    end if
                    rep%core_ham_diag(i) = hamiltonian_row(j)
                    ! We calculate and store the diagonal matrix element at
                    ! this point for later access.
                    if (.not. tReadPops) &
                        call set_det_diagH(i, Real(hamiltonian_row(j), dp))
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, IlutI, IlutJ)
                    else if (tGUGA) then
                        ! for the off-diagonal elements i have to call the GUGA
                        ! specific function
                        ! but this is a waste.. i do not have to do that for
                        ! every nJ i could just check the list generated
                        ! by H|nI>..
                        call calc_guga_matrix_element(IlutI, IlutJ, &
                                                      excitInfo, tmp_mat, .true., 1)
#ifdef DEBUG_
                        call calc_guga_matrix_element(IlutI, IlutJ, &
                                                      excitInfo, tmp_mat_2, .true., 2)
                        if (.not. near_zero(tmp_mat - tmp_mat_2)) then
                            call stop_all(this_routine, "type 1 and 2 do not agree!")
                        end if
                        call calc_guga_matrix_element(IlutJ, IlutI, &
                                                      excitInfo, tmp_mat_2, .true., 2)
                        if (.not. near_zero(tmp_mat - tmp_mat_2)) then
                            call stop_all(this_routine, "not hermititan!")
                        end if
#endif

#ifdef CMPLX_
                        hamiltonian_row(j) = conjg(tmp_mat)
#else
                        hamiltonian_row(j) = tmp_mat
#endif
                    else

                        tmp = ieor(IlutI(0:NIfD), IlutJ(0:NIfD))
                        tmp = iand(IlutI(0:NIfD), tmp)
                        IC = CountBits(tmp, NIfD)

                        if (IC <= maxExcit) then
                            hamiltonian_row(j) = get_helement(nI, nJ, IC, ilutI, IlutJ)
                        end if
                    end if
                    if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                end if
            end do
            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(rep%sparse_core_ham, i, row_size, "sparse_core_ham", rep%SparseCoreHamilTags(:, i))

            rep%sparse_core_ham(i)%elements = 0.0_dp
            rep%sparse_core_ham(i)%positions = 0
            rep%sparse_core_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, rep%determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + rep%determ_displs(iProcIndex))) then
                    rep%sparse_core_ham(i)%positions(counter) = j
                    rep%sparse_core_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        ! Don't time the mpi_barrier call, because if we did then we wouldn't
        ! be able separate out some of the core Hamiltonian creation time from
        ! the MPIBarrier calls in the main loop.
        call MPIBarrier(ierr, tTimeIn=.false.)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(this_routine, TempStoreTag, ierr)
        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(this_routine, HRTag, ierr)
!         deallocate(temp_store_nI, stat=ierr)

    end subroutine calc_determ_hamil_sparse

    subroutine calc_determ_hamil_sparse_hphf(rep)
        type(core_space_t), intent(inout) :: rep
        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:, :) :: temp_store
        integer, allocatable :: temp_store_nI(:, :)
        integer(TagIntType) :: HRTag, TempStoreTag
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calc_determ_hamil_sparse_hphf"

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC
        logical :: CS_I
        logical, allocatable :: cs(:)

        allocate(rep%sparse_core_ham(rep%determ_sizes(iProcIndex)), stat=ierr)
        allocate(rep%SparseCoreHamilTags(2, rep%determ_sizes(iProcIndex)))
        allocate(hamiltonian_row(rep%determ_space_size), stat=ierr)
        call LogMemAlloc('hamiltonian_row', int(rep%determ_space_size, sizeof_int), 8, t_r, HRTag, ierr)
        allocate(rep%core_ham_diag(rep%determ_sizes(iProcIndex)), stat=ierr)
        allocate(temp_store(0:NIfTot, rep%determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', rep%determ_space_size * (NIfTot + 1), 8, t_r, TempStoreTag, ierr)
        allocate(temp_store_nI(nel, rep%determ_space_size), stat=ierr)
        allocate(cs(rep%determ_space_size), stat=ierr)

        ! Stick together the deterministic states from all processors, on
        ! all processors.
        ! n.b. Explicitly use 0:NIfTot, as NIfTot may not equal NIfBCast
        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:rep%determ_sizes(iProcIndex)), &
                           temp_store, rep%determ_sizes, rep%determ_displs)

        do i = 1, rep%determ_space_size
            call decode_bit_det(temp_store_nI(:, i), temp_store(:, i))
            cs(i) = TestClosedShellDet(temp_store(:, i))
        end do

        ! Loop over all deterministic states on this processor.
        do i = 1, rep%determ_sizes(iProcIndex)

            !call decode_bit_det(nI, SpawnedParts(:, i))
            nI = temp_store_nI(:, i + rep%determ_displs(iProcIndex))

            row_size = 0
            hamiltonian_row = 0.0_dp

            CS_I = cs(i + rep%determ_displs(iProcIndex))

            ! Loop over all deterministic states.
            do j = 1, rep%determ_space_size

                !call decode_bit_det(nJ, temp_store(:,j))
                nJ = temp_store_nI(:, j)

                ! If on the diagonal of the Hamiltonian.
                if (j == i + rep%determ_displs(iProcIndex)) then
                    hamiltonian_row(j) = hphf_diag_helement(nI, SpawnedParts(:, i)) - Hii
                    rep%core_ham_diag(i) = hamiltonian_row(j)
                    ! We calculate and store the diagonal matrix element at
                    ! this point for later access.
                    if (.not. tReadPops) &
                        call set_det_diagH(i, Real(hamiltonian_row(j), dp))
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    tmp = ieor(SpawnedParts(0:NIfD, i), temp_store(0:NIfD, j))
                    tmp = iand(SpawnedParts(0:NIfD, i), tmp)
                    IC = CountBits(tmp, NIfD)

                    if (IC <= maxExcit .or. ((.not. CS_I) .and. (.not. cs(j)))) then
                        hamiltonian_row(j) = hphf_off_diag_helement_opt(nI, SpawnedParts(:, i), temp_store(:, j), IC, CS_I, cs(j))
                        if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                    end if
                end if

            end do

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(rep%sparse_core_ham, i, row_size, "sparse_core_ham", rep%SparseCoreHamilTags(:, i))

            rep%sparse_core_ham(i)%elements = 0.0_dp
            rep%sparse_core_ham(i)%positions = 0
            rep%sparse_core_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, rep%determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + rep%determ_displs(iProcIndex))) then
                    rep%sparse_core_ham(i)%positions(counter) = j
                    rep%sparse_core_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        ! Don't time the mpi_barrier call, because if we did then we wouldn't
        ! be able separate out some of the core Hamiltonian creation time from
        ! the MPIBarrier calls in the main loop.
        call MPIBarrier(ierr, tTimeIn=.false.)

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)
        deallocate(hamiltonian_row, stat=ierr)
        call LogMemDealloc(t_r, HRTag, ierr)
        deallocate(temp_store_nI, stat=ierr)
        deallocate(cs, stat=ierr)

    end subroutine calc_determ_hamil_sparse_hphf

    subroutine calc_approx_hamil_sparse_hphf(rep)
        type(core_space_t), intent(in) :: rep
        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer, allocatable :: temp_store_nI(:, :)
        HElement_t(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calc_approx_hamil_sparse_hphf"

        integer(n_int) :: tmp(0:NIfD)
        integer :: IC
        logical :: CS_I, var_state_i, var_state_j
        logical, allocatable :: cs(:)

        allocate(approx_ham(rep%determ_sizes(iProcIndex)), stat=ierr)
        allocate(hamiltonian_row(rep%determ_space_size), stat=ierr)
        allocate(temp_store_nI(nel, rep%determ_space_size), stat=ierr)
        allocate(cs(rep%determ_space_size), stat=ierr)

        do i = 1, rep%determ_space_size
            call decode_bit_det(temp_store_nI(:, i), rep%core_space(:, i))
            cs(i) = TestClosedShellDet(rep%core_space(:, i))
        end do

        ! Loop over all deterministic states on this processor.
        do i = 1, rep%determ_sizes(iProcIndex)

            !call decode_bit_det(nI, SpawnedParts(:, i))
            nI = temp_store_nI(:, i + rep%determ_displs(iProcIndex))

            var_state_i = is_var_state(rep%core_space(:, i + rep%determ_displs(iProcIndex)), nI)

            row_size = 0
            hamiltonian_row = 0.0_dp

            CS_I = cs(i + rep%determ_displs(iProcIndex))

            ! Loop over all deterministic states.
            do j = 1, rep%determ_space_size

                !call decode_bit_det(nJ, core_space(:,j))
                nJ = temp_store_nI(:, j)

                ! If on the diagonal of the Hamiltonian.
                if (j == i + rep%determ_displs(iProcIndex)) then
                    hamiltonian_row(j) = hphf_diag_helement(nI, rep%core_space(:, i + rep%determ_displs(iProcIndex))) - Hii
                    !core_ham_diag(i) = hamiltonian_row(j)
                    ! We calculate and store the diagonal matrix element at
                    ! this point for later access.
                    if (.not. tReadPops) &
                        call set_det_diagH(i, Real(hamiltonian_row(j), dp))
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    var_state_j = is_var_state(rep%core_space(:, j), nJ)

                    ! Only add a matrix element if both states are variational
                    ! states, or if one of them is (the rectangular portion of
                    ! H connected var_space to the space of connections to it).
                    if (var_state_i .or. var_state_j) then
                        tmp = ieor(rep%core_space(0:NIfD, i + rep%determ_displs(iProcIndex)), rep%core_space(0:NIfD, j))
                        tmp = iand(rep%core_space(0:NIfD, i + rep%determ_displs(iProcIndex)), tmp)
                        IC = CountBits(tmp, NIfD)

                        if (IC <= maxExcit .or. ((.not. CS_I) .and. (.not. cs(j)))) then

                            hamiltonian_row(j) = hphf_off_diag_helement_opt(nI, &
                                rep%core_space(:, i + rep%determ_displs(iProcIndex)), &
                                rep%core_space(:, j), IC, CS_I, cs(j))

                            if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                        end if
                    end if
                end if

            end do

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            allocate(approx_ham(i)%elements(row_size), stat=ierr)
            allocate(approx_ham(i)%positions(row_size), stat=ierr)

            approx_ham(i)%elements = 0.0_dp
            approx_ham(i)%positions = 0
            approx_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, rep%determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + rep%determ_displs(iProcIndex))) then
                    approx_ham(i)%positions(counter) = j
                    approx_ham(i)%elements(counter) = hamiltonian_row(j)
                    counter = counter + 1
                end if
                if (counter == row_size + 1) exit
            end do

        end do

        ! Don't time the mpi_barrier call, because if we did then we wouldn't
        ! be able separate out some of the core Hamiltonian creation time from
        ! the MPIBarrier calls in the main loop.
        call MPIBarrier(ierr, tTimeIn=.false.)

        deallocate(hamiltonian_row, stat=ierr)
        deallocate(temp_store_nI, stat=ierr)
        deallocate(cs, stat=ierr)

    end subroutine calc_approx_hamil_sparse_hphf

    subroutine allocate_sparse_ham_row(sparse_matrix, row, sparse_row_size, sparse_matrix_name, sparse_tags)

        ! Allocate a single row and add it to the memory manager.

        type(sparse_matrix_real), intent(inout) :: sparse_matrix(:)
        integer, intent(in) :: row, sparse_row_size
        character(len=*), intent(in) :: sparse_matrix_name
        integer(TagIntType), intent(inout) :: sparse_tags(2)
        integer :: ierr
        character(len=1024) :: string_row
        character(len=1024) :: var_name
        character(len=*), parameter :: t_r = "allocate_sparse_ham_row"

        write(string_row, '(I10)') row

        var_name = trim(sparse_matrix_name)//"_"//trim(string_row)//"_elements"
        allocate(sparse_matrix(row)%elements(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, 8, t_r, sparse_tags(1), ierr)

        var_name = trim(sparse_matrix_name)//"_"//trim(string_row)//"_positions"
        allocate(sparse_matrix(row)%positions(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, bytes_int, t_r, sparse_tags(2), ierr)

    end subroutine allocate_sparse_ham_row

    subroutine deallocate_core_hashtable(ht)

        type(core_hashtable), intent(inout), allocatable :: ht(:)

        integer :: i, ierr

        if (allocated(ht)) then
            do i = 1, size(ht)
                if (allocated(ht(i)%ind)) then
                    deallocate(ht(i)%ind, stat=ierr)
                    if (ierr /= 0) write(stdout, '("Error when deallocating core hashtable ind array:",1X,i8)') ierr
                end if
            end do

            deallocate(ht, stat=ierr)
            if (ierr /= 0) write(stdout, '("Error when deallocating core hashtable:",1X,i8)') ierr
        end if

    end subroutine deallocate_core_hashtable

    subroutine deallocate_trial_hashtable(ht)

        type(trial_hashtable), intent(inout), allocatable :: ht(:)

        integer :: i, ierr

        if (allocated(ht)) then
            do i = 1, size(ht)
                if (allocated(ht(i)%states)) then
                    deallocate(ht(i)%states, stat=ierr)
                    if (ierr /= 0) write(stdout, '("Error when deallocating trial hashtable states array:",1X,i8)') ierr
                end if
            end do

            deallocate(ht, stat=ierr)
            if (ierr /= 0) write(stdout, '("Error when deallocating core hashtable:",1X,i8)') ierr
        end if

    end subroutine deallocate_trial_hashtable

    function is_var_state(ilut, nI) result(var_state)

        use FciMCData, only: var_space, var_space_size_int
        use hash, only: FindWalkerHash

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer(int64) :: hash_val
        integer(int64) :: pos
        logical :: var_state

        hash_val = FindWalkerHash(nI, var_space_size_int)

        call var_ht%callback_lookup(hash_val, pos, var_state, loc_verify)
    contains

        function loc_verify(ind) result(match)
            integer(int64), intent(in) :: ind
            logical :: match

            match = all(ilut(0:nifd) == var_space(0:nifd, ind))

        end function loc_verify
    end function is_var_state

end module sparse_arrays
