#include "macros.h"

module core_space_util
    use mpi
    use constants
    use shared_rhash, only: shared_rhash_t
    use shared_memory_mpi
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use ParallelHelper
    implicit none

    private
    public :: core_space_t, cs_replicas, sparse_matrix_real, sparse_matrix_int, &
        deallocate_sparse_ham, min_pt, max_pt

    integer, parameter :: min_pt = 1
    integer, parameter :: max_pt = rep_size

    type sparse_matrix_real
        HElement_t(dp), allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_real

    type sparse_matrix_int
        integer, allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: positions
        integer :: num_elements
    end type sparse_matrix_int
    
    ! Info type containing the characterization of the core space on a given replica
    type :: core_space_t
        ! determ_sizes(i) holds the core space size on processor i.
        integer(MPIArg), allocatable, dimension(:) :: determ_sizes
        ! determ_displs(i) holds sum(determ_sizes(i-1)), that is, the
        ! total number of core states on all processors up to processor i.
        ! (determ_displs(1) == 0).
        integer(MPIArg), allocatable, dimension(:) :: determ_displs
        ! determ_last(i) holds the final index belonging process i.
        integer(MPIArg), allocatable, dimension(:) :: determ_last
        ! The total size of the core space on all processors.
        integer(MPIArg) :: determ_space_size
        ! determ_space_size_int is identical to determ_space_size, but converted
        ! to the default integer kind.
        integer :: determ_space_size_int
        ! This vector will store the indicies of the deterministic states in CurrentDets. This is worked out in the main loop.
        integer, allocatable, dimension(:) :: indices_of_determ_states

        ! This stores the entire core space from all processes, on each process.
        integer(n_int), pointer, dimension(:,:) :: core_space => null()
        integer(n_int), pointer, dimension(:,:) :: core_space_direct => null()
        integer(MPIArg) :: core_space_win

        type(shared_rhash_t) :: core_ht

        ! The core Hamiltonian for semi-stochastiic simulations.
        type(sparse_matrix_real), allocatable, dimension(:) :: sparse_core_ham
        integer(TagIntType), allocatable, dimension(:,:) :: SparseCoreHamilTags
        ! This stores all the amplitudes of the walkers in the deterministic space. This vector has the size of the part
        ! of the deterministic space stored on *this* processor only. It is therefore used to store the deterministic vector
        ! on this processor, before it is combined to give the whole vector, which is stored in full_determ_vecs.
        ! Later in the iteration, it is also used to store the result of the multiplication by the core Hamiltonian on
        ! full_determ_vecs.
        real(dp), allocatable, dimension(:,:) :: partial_determ_vecs
        real(dp), allocatable, dimension(:,:) :: full_determ_vecs
        real(dp), allocatable, dimension(:,:) :: full_determ_vecs_av
        integer(TagIntType) :: FDetermTag, FDetermAvTag, PDetermTag, IDetermTag, CoreSpaceTag
        ! Stores the parities for all connected pairs of states in the core space.
        type(sparse_matrix_int), allocatable, dimension(:) :: core_connections        

        ! The diagonal elements of the core-space Hamiltonian (with Hii taken away).
        real(dp), allocatable, dimension(:) :: core_ham_diag

        ! ilut Sign range in which this core space operates
        integer :: max_run, min_run

        ! Is this a global core space?
        logical :: t_global
    contains        
        procedure :: alloc_wf
        procedure :: associate_run
        procedure :: dealloc

        procedure :: max_part
        procedure :: min_part

        procedure :: first_run
        procedure :: last_run
    end type core_space_t

    ! Each replica has its own core space information
    type(core_space_t), allocatable :: cs_replicas(:)
    
contains

    ! Interface function giving access to the core space replicas: This serves to
    ! switch between models with one core space per replica or one shared core space
    ! efficiently
    pure function max_part(this) result(ir)
        class(core_space_t), intent(in) :: this
        integer :: ir

        unused_var(this)
        ir = max_part_type(this%max_run)
    end function max_part

    pure function min_part(this) result(ir)
        class(core_space_t), intent(in) :: this
        integer :: ir

        unused_var(this)
        ir = min_part_type(this%min_run)
    end function min_part

!------------------------------------------------------------------------------------------!    

    pure function last_run(this) result(ir)
        class(core_space_t), intent(in) :: this
        integer :: ir

        ir = this%max_run
    end function last_run

    pure function first_run(this) result(ir)
        class(core_space_t), intent(in) :: this
        integer :: ir

        ir = this%min_run
    end function first_run

    !------------------------------------------------------------------------------------------!

    ! Associate this core space with a replica
    subroutine associate_run(this, t_global, run)
        class(core_space_t), intent(inout) :: this
        logical, intent(in) :: t_global
        integer, intent(in) :: run

        this%t_global = t_global
        if(this%t_global) then
            this%min_run = 1
            this%max_run = inum_runs
        else
            this%min_run = run
            this%max_run = run
        end if        
    end subroutine associate_run

    !------------------------------------------------------------------------------------------!

    
    subroutine alloc_wf(this)
        class(core_space_t), intent(inout) :: this        
        integer :: vec_size
        integer :: ierr
        character(*), parameter :: t_r = "core_space_t%alloc_wf"

        ! Store the operating range of this core space        
        if(this%t_global) then
            vec_size = lenof_sign
        else
            vec_size = rep_size
        end if

        ! Allocate the vectors to store the walker amplitudes and the deterministic Hamiltonian.       
        allocate(this%full_determ_vecs(vec_size,this%determ_space_size), stat=ierr)
        call LogMemAlloc('full_determ_vecs', this%determ_space_size_int*lenof_sign, &
            8, t_r, this%FDetermTag, ierr)
        allocate(this%full_determ_vecs_av(vec_size,this%determ_space_size), stat=ierr)
        call LogMemAlloc('full_determ_vecs_av', this%determ_space_size_int*lenof_sign, &
            8, t_r, this%FDetermAvTag, ierr)
        allocate(this%partial_determ_vecs(vec_size,this%determ_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('partial_determ_vecs', int(this%determ_sizes(iProcIndex), &
            sizeof_int)*lenof_sign, 8, t_r, this%PDetermTag, ierr)
        
        this%full_determ_vecs = 0.0_dp
        this%full_determ_vecs_av = 0.0_dp
        this%partial_determ_vecs = 0.0_dp
        
    end subroutine alloc_wf

!------------------------------------------------------------------------------------------!    

    subroutine dealloc(this)
        class(core_space_t), intent(inout) :: this
        character(*), parameter :: t_r = "core_space_t%dealloc"
        integer :: ierr
        call this%core_ht%dealloc()
        if (allocated(this%indices_of_determ_states)) then
            deallocate(this%indices_of_determ_states, stat=ierr)
            call LogMemDealloc(t_r, this%IDetermTag, ierr)
        end if
        if (associated(this%core_space)) then
            this%core_space => null()
            call shared_deallocate_mpi(this%core_space_win, this%core_space_direct)
            call LogMemDealloc(t_r, this%CoreSpaceTag, ierr)
        end if

        if (allocated(this%determ_sizes)) then
            deallocate(this%determ_sizes, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating determ_sizes:",1X,i8)') ierr
        end if
        if (allocated(this%determ_displs)) then
            deallocate(this%determ_displs, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating determ_displs:",1X,i8)') ierr
        end if
        if (allocated(this%determ_last)) then
            deallocate(this%determ_last, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating determ_last:",1X,i8)') ierr
        end if

        if (allocated(this%partial_determ_vecs)) then
            deallocate(this%partial_determ_vecs, stat=ierr)
            call LogMemDealloc(t_r, this%PDetermTag, ierr)
        end if
        if (allocated(this%full_determ_vecs)) then
            deallocate(this%full_determ_vecs, stat=ierr)
            call LogMemDealloc(t_r, this%FDetermTag, ierr)
        end if
        if (allocated(this%full_determ_vecs_av)) then
            deallocate(this%full_determ_vecs_av, stat=ierr)
            call LogMemDealloc(t_r, this%FDetermAvTag, ierr)
        end if

        call deallocate_sparse_ham(this%sparse_core_ham, this%SparseCoreHamilTags)

        call deallocate_sparse_matrix_int(this%core_connections)
        if (allocated(this%core_ham_diag)) then
            deallocate(this%core_ham_diag, stat=ierr)
        end if        

    end subroutine dealloc

    subroutine deallocate_sparse_ham(sparse_matrix, sparse_tags)

        ! Deallocate the whole array, and remove all rows from the memory manager.

        type(sparse_matrix_real), intent(inout), allocatable :: sparse_matrix(:)
        integer(TagIntType), intent(inout), allocatable :: sparse_tags(:,:)
        integer :: sparse_matrix_size, i, ierr
        character(len=*), parameter :: t_r = "deallocate_sparse_ham"

        sparse_matrix_size = size(sparse_matrix)

        do i = sparse_matrix_size, 1, -1

            deallocate(sparse_matrix(i)%elements, stat=ierr)
            !call LogMemDealloc(t_r, sparse_tags(1,i), ierr)

            deallocate(sparse_matrix(i)%positions, stat=ierr)
            !call LogMemDealloc(t_r, sparse_tags(2,i), ierr)

        end do

        if (allocated(sparse_tags)) deallocate(sparse_tags)
        if (allocated(sparse_matrix)) deallocate(sparse_matrix)

    end subroutine deallocate_sparse_ham

    subroutine deallocate_sparse_matrix_int(sparse_mat)

        type(sparse_matrix_int), intent(inout), allocatable :: sparse_mat(:)

        integer :: i, ierr

        if (allocated(sparse_mat)) then
            do i = 1, size(sparse_mat)
                if (allocated(sparse_mat(i)%elements)) then
                    deallocate(sparse_mat(i)%elements, stat=ierr)
                    if (ierr /= 0) write(6,'("Error when deallocating sparse matrix elements array:",1X,i8)') ierr
                end if
                if (allocated(sparse_mat(i)%positions)) then
                    deallocate(sparse_mat(i)%positions, stat=ierr)
                    if (ierr /= 0) write(6,'("Error when deallocating sparse matrix positions array:",1X,i8)') ierr
                end if
            end do

            deallocate(sparse_mat, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating sparse matrix array:",1X,i8)') ierr
        end if

    end subroutine deallocate_sparse_matrix_int    

end module core_space_util
