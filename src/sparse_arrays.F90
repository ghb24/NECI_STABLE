#include  "macros.h"

! This module contains a type and routines for defining and creating a sparse Hamiltonian.
! The type created to store this information is sparse_matrix_real. For an N-by-N matrix,
! one creates a 1d array of type sparse_matrix_real. Each element of this array stores the
! information on one single row of the matrix - the positions of the non-zero elements
! and their values, in the same order. This completley defines the matrix.

! Note that in some parts of NECI, the sparse nature is used in a different form (i.e. in
! the Lanczos code).

module sparse_arrays

    use CalcData, only: tRegenDiagHEls, tReadPops
    use constants
    use bit_rep_data, only: NIfTot, NIfDBO
    use bit_reps, only: decode_bit_det
    use Determinants, only: get_helement
    use FciMCData, only: determ_space_size, determ_proc_sizes, determ_proc_indices, &
                         SpawnedParts, CurrentH, Hii, core_ham_diag
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBarrier, MPIAllGatherV
    use SystemData, only: tHPHF, nel

    implicit none

    type sparse_matrix_real
        real(dp), allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_real

    type sparse_matrix_int
        integer, allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_int

    type trial_hashtable
        ! All the states with this hash value.
        integer(n_int), allocatable, dimension(:,:) :: states
        ! The number of clashes for ths hash value.
        integer :: nclash
    end type trial_hashtable

    type(sparse_matrix_real), allocatable, dimension(:) :: sparse_ham
    integer(TagIntType), allocatable, dimension(:,:) :: SparseHamilTags

    ! For quick access it is often useful to have just the diagonal elements. Note,
    ! however, that they *are* stored in sparse_ham too.
    real(dp), allocatable, dimension(:) :: hamil_diag
    integer(TagIntType) :: HDiagTag

    ! The core Hamiltonian for semi-stochastiic simulations.
    type(sparse_matrix_real), allocatable, dimension(:) :: sparse_core_ham
    integer(TagIntType), allocatable, dimension(:,:) :: SparseCoreHamilTags

    ! Stores the parities for all connected pairs of states in the core space.
    type(sparse_matrix_int), allocatable, dimension(:) :: core_connections

    type(trial_hashtable), allocatable, dimension(:) :: trial_ht
    type(trial_hashtable), allocatable, dimension(:) :: con_ht

contains

    subroutine calculate_sparse_hamiltonian(num_states, ilut_list)

        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        real(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character(len=*), parameter :: t_r = "calculate_sparse_hamiltonian"

        allocate(sparse_ham(num_states))
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

        allocate(SparseHamilTags(2, num_states))

        ! Set each element to one to count the diagonal elements straight away.
        sparse_row_sizes = 1

        do i = 1, num_states

            hamiltonian_row = 0.0_dp

            call decode_bit_det(nI, ilut_list(:, i))

            ! sparse_diag_positions(i) stores the number of non-zero elements in row i of the
            ! Hamiltonian, up to and including the diagonal element.
            ! sparse_row_sizes(i) stores this number currently as all non-zero elements before
            ! the diagonal have been counted (as the Hamiltonian is symmetric).
            sparse_diag_positions(i) = sparse_row_sizes(i)

            do j = i, num_states

                call decode_bit_det(nJ, ilut_list(:, j))

                ! If on the diagonal of the Hamiltonian.
                if (i == j) then
                    if (.not. tHPHF) then
                        hamiltonian_row(j) = get_helement(nI, nJ, 0)
                    else
                        hamiltonian_row(j) = hphf_diag_helement(nI, ilut_list(:, i))
                    end if
                    hamil_diag(j) = hamiltonian_row(j)
                else
                    if (.not. tHPHF) then
                        hamiltonian_row(j) = get_helement(nI, nJ, ilut_list(:, i), &
                                                                 ilut_list(:, j))
                    else
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, ilut_list(:, i), &
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
            call allocate_sparse_ham_row(sparse_ham, i, sparse_row_sizes(i), "sparse_ham", SparseHamilTags(:,i)) 

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
                    &elements(indices(sparse_ham(i)%positions(j))) = sparse_ham(i)%elements(j)

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

    subroutine calc_determ_hamil_sparse()

        integer :: i, j, row_size, counter, ierr
        integer :: nI(nel), nJ(nel)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer(TagIntType) :: HRTag, TempStoreTag
        real(dp), allocatable, dimension(:) :: hamiltonian_row
        character(len=*), parameter :: t_r = "calculate_det_hamiltonian_sparse"

        allocate(sparse_core_ham(determ_proc_sizes(iProcIndex)), stat=ierr)
        allocate(SparseCoreHamilTags(2, determ_proc_sizes(iProcIndex)))
        allocate(hamiltonian_row(determ_space_size), stat=ierr)
        call LogMemAlloc('hamiltonian_row', int(determ_space_size,sizeof_int), 8, t_r, HRTag, ierr)
        allocate(core_ham_diag(determ_proc_sizes(iProcIndex)), stat=ierr)
        allocate(temp_store(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', maxval(determ_proc_sizes)*(NIfTot+1), 8, t_r, &
                         TempStoreTag, ierr)

        ! Stick together the deterministic states from all processors, on all processors.
        call MPIAllGatherV(SpawnedParts(:,1:determ_proc_sizes(iProcIndex)), temp_store, &
                       determ_proc_sizes, determ_proc_indices)

        ! Loop over all deterministic states on this processor.
        do i = 1, determ_proc_sizes(iProcIndex)

            call decode_bit_det(nI, SpawnedParts(:, i))

            row_size = 0
            hamiltonian_row = 0.0_dp

            ! Loop over all deterministic states.
            do j = 1, determ_space_size

                call decode_bit_det(nJ, temp_store(:,j))

                ! If on the diagonal of the Hamiltonian.
                if (all( SpawnedParts(0:NIfDBO, i) == temp_store(0:NIfDBO, j) )) then
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_diag_helement(nI, SpawnedParts(:,i)) - Hii
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, 0) - Hii
                    end if
                    core_ham_diag(i) = hamiltonian_row(j)
                    ! We calculate and store CurrentH at this point for ease.
                    if ((.not. tRegenDiagHEls) .and. (.not. tReadPops)) CurrentH(1,i) = hamiltonian_row(j)
                    ! Always include the diagonal elements.
                    row_size = row_size + 1
                else
                    if (tHPHF) then
                        hamiltonian_row(j) = hphf_off_diag_helement(nI, nJ, SpawnedParts(:,i), temp_store(:,j))
                    else
                        hamiltonian_row(j) = get_helement(nI, nJ, SpawnedParts(:, i), temp_store(:, j))
                    end if
                    if (abs(hamiltonian_row(j)) > 0.0_dp) row_size = row_size + 1
                end if

            end do

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so allocate it.
            call allocate_sparse_ham_row(sparse_core_ham, i, row_size, "sparse_core_ham", SparseCoreHamilTags(:,i))

            sparse_core_ham(i)%elements = 0.0_dp
            sparse_core_ham(i)%positions = 0
            sparse_core_ham(i)%num_elements = row_size

            counter = 1
            do j = 1, determ_space_size
                ! If non-zero or a diagonal element.
                if (abs(hamiltonian_row(j)) > 0.0_dp .or. (j == i + determ_proc_indices(iProcIndex)) ) then
                    sparse_core_ham(i)%positions(counter) = j
                    sparse_core_ham(i)%elements(counter) = hamiltonian_row(j)
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

    end subroutine calc_determ_hamil_sparse

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

        write (string_row, '(I10)') row

        var_name = trim(sparse_matrix_name)//"_"//trim(string_row)//"_elements"
        allocate(sparse_matrix(row)%elements(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, 8, t_r, sparse_tags(1), ierr)

        var_name = trim(sparse_matrix_name)//"_"//trim(string_row)//"_positions"
        allocate(sparse_matrix(row)%positions(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, bytes_int, t_r, sparse_tags(2), ierr)

    end subroutine allocate_sparse_ham_row

    subroutine deallocate_sparse_ham(sparse_matrix, sparse_matrix_name, sparse_tags)

        ! Deallocate the whole array, and remove all rows from the memory manager.

        type(sparse_matrix_real), intent(inout), allocatable :: sparse_matrix(:)
        character(len=*), intent(in) :: sparse_matrix_name
        integer(TagIntType), intent(inout), allocatable :: sparse_tags(:,:)
        integer :: sparse_matrix_size, i, ierr
        character(len=*), parameter :: t_r = "deallocate_sparse_ham"

        sparse_matrix_size = size(sparse_matrix)

        do i = sparse_matrix_size, 1, -1

            deallocate(sparse_matrix(i)%elements, stat=ierr)
            call LogMemDealloc(t_r, sparse_tags(1,i), ierr)

            deallocate(sparse_matrix(i)%positions, stat=ierr)
            call LogMemDealloc(t_r, sparse_tags(2,i), ierr)

        end do

        deallocate(sparse_tags)
        deallocate(sparse_matrix)

    end subroutine deallocate_sparse_ham

end module sparse_arrays
