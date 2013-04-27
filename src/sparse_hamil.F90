#include  "macros.h"

! This module contains a type and routines for defining and creating a sparse Hamiltonian.
! The type created to store this information is sparse_matrix_info. For an N-by-N matrix,
! one creates a 1d array of type sparse_matrix_info. Each element of this array stores the
! information on one single row of the matrix - the positions of the non-zero elements
! and their values, in the same order. This completley defines the matrix.

! Note that in some parts of NECI, the sparse nature is used in a different form (i.e. in
! the Lanczos code).

module sparse_hamil

    use constants
    use bit_rep_data, only: NIfTot
    use bit_reps, only: decode_bit_det
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use SystemData, only: tHPHF, nel

    implicit none

    type sparse_matrix_info
        real(dp), allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_info

    type(sparse_matrix_info), allocatable, dimension(:) :: sparse_ham

    ! For quick access it is often useful to have just the diagonal elements. Note,
    ! however, that they *are* stored in sparse_ham too.
    real(dp), allocatable, dimension(:) :: hamil_diag

    integer(TagIntType) :: HDiagTag
    integer(TagIntType), allocatable, dimension(:,:) :: SparseHamilTags

contains

    subroutine calculate_sparse_hamiltonian(num_states, ilut_list)

        integer, intent(in) :: num_states
        integer(n_int), intent(in) :: ilut_list(0:NIfTot, num_states)
        integer :: i, j, counter, ierr
        integer :: nI(nel), nJ(nel)
        real(dp), allocatable, dimension(:) :: hamiltonian_row
        integer, allocatable, dimension(:) :: sparse_diag_positions, sparse_row_sizes, indices
        integer(TagIntType) :: HRTag, SRTag, SDTag, ITag
        character (len=*), parameter :: this_routine = "calculate_sparse_hamiltonian"

        allocate(sparse_ham(num_states))
        allocate(hamiltonian_row(num_states), stat=ierr)
        call LogMemAlloc('hamiltonian_row', num_states, 8, this_routine, HRTag, ierr)
        allocate(hamil_diag(num_states), stat=ierr)
        call LogMemAlloc('hamil_diag', num_states, 8, this_routine, HDiagTag, ierr)
        allocate(sparse_row_sizes(num_states), stat=ierr)
        call LogMemAlloc('sparse_row_sizes', num_states, bytes_int, this_routine, SRTag, ierr)
        allocate(sparse_diag_positions(num_states), stat=ierr)
        call LogMemAlloc('sparse_diag_positions', num_states, bytes_int, this_routine, SDTag, ierr)
        allocate(indices(num_states), stat=ierr)
        call LogMemAlloc('indices', num_states, bytes_int, this_routine, ITag, ierr)

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

            ! Now we know the number of non-zero elements in this row of the Hamiltonian, so
            ! allocate it.
            call allocate_sparse_ham_row(i, sparse_row_sizes(i))

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
        call LogMemDealloc(this_routine, HRTag, ierr)
        deallocate(sparse_row_sizes, stat=ierr)
        call LogMemDealloc(this_routine, SRTag, ierr)
        deallocate(sparse_diag_positions, stat=ierr)
        call LogMemDealloc(this_routine, SDTag, ierr)
        deallocate(indices, stat=ierr)
        call LogMemDealloc(this_routine, ITag, ierr)

    end subroutine calculate_sparse_hamiltonian

    subroutine allocate_sparse_ham_row(row, sparse_row_size)

        ! Allocate a single row and add it to the memory manager.

        integer, intent(in) :: row, sparse_row_size
        integer :: ierr
        character (len=1024) :: string_row
        character (len=1024) :: var_name
        character (len=*), parameter :: this_routine = "allocate_sparse_ham_row"

        write (string_row, '(I10)') row

        var_name = "sparse_ham_"//trim(string_row)//"_elements"
        allocate(sparse_ham(row)%elements(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, 8, this_routine, &
                         SparseHamilTags(1,row), ierr)

        var_name = "sparse_ham_"//trim(string_row)//"_positions"
        allocate(sparse_ham(row)%positions(sparse_row_size), stat=ierr)
        call LogMemAlloc(var_name, sparse_row_size, bytes_int, this_routine, &
                         SparseHamilTags(2,row), ierr)

    end subroutine allocate_sparse_ham_row

    subroutine deallocate_sparse_ham()

        ! Deallocate the whole array, and remove all rows from the memory manager.

        integer :: sparse_ham_size, i, ierr
        character (len=*), parameter :: this_routine = "allocate_sparse_ham"

        sparse_ham_size = size(sparse_ham)

        do i = sparse_ham_size, 1, -1

            deallocate(sparse_ham(i)%elements, stat=ierr)
            call LogMemDealloc(this_routine, SparseHamilTags(1,i), ierr)

            deallocate(sparse_ham(i)%positions, stat=ierr)
            call LogMemDealloc(this_routine, SparseHamilTags(2,i), ierr)
        end do

        deallocate(SparseHamilTags)
        deallocate(sparse_ham)

    end subroutine deallocate_sparse_ham

end module sparse_hamil
