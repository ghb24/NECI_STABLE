#include "macros.h"
module hamiltonian_linalg
    !   robert.anderson@kcl.ac.uk
    !   June 2016
    !
    !   Linear algebra routines for CI hamiltonians
    !
    use constants
    use ras_data
    use FciMCData, only : hamiltonian
    use sparse_arrays, only: sparse_ham, hamil_diag, HDiagTag
    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg, MPIBarrier
    use Parallel_neci, only: MPIBCast, MPIGatherV, MPIAllGather
    use ParallelHelper, only: root
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc

    implicit none
    ! The value of hamil_type specifies what form the Hamiltonian is stored in.
    ! The following options are currently available:
    integer, parameter :: full_hamil_type = 1
    integer, parameter :: sparse_hamil_type = 2
    integer, parameter :: parallel_sparse_hamil_type = 3
    integer, parameter :: direct_ci_type = 4

    type HamiltonianCalcType
        integer :: hamil_type
        ! The dimension of the vector space we are working in, as determined by the number
        ! of rows and columns in the Hamiltonian matrix.
        integer :: space_size 
        ! For parallel calculations, store all spaces sizes on each processor
        integer(MPIArg), allocatable, dimension(:) :: space_sizes
        ! For parallel calculations, this vector is the size of the space on this processor. This
        ! vector is used to store the output of H |ket> on this processor.
        ! (formerly partial_davidson_vector)
        real(dp), allocatable, dimension(:) :: partial_H_ket
        ! displacements necessary for communication of H |ket>
        ! (formerly davidson_disps)
        integer(MPIArg), allocatable, dimension(:) :: partial_H_ket_disps
        ! All algorithms for solving large eigenproblems involve a unitary rotation of the 
        ! Hamiltonian into a smaller basis. basis_vectors(:,i) is the ith such unit vector
        real(dp), allocatable, dimension(:,:) :: basis_vectors
        ! the location in the basis is the HF determinant
        integer :: hfindex
        ! the hamiltonian projected into basis_vectors
        real(dp), allocatable, dimension(:,:) :: projected_hamil
        ! we'll usually need some working space for diagonalisation of H in the small basis
        real(dp), allocatable, dimension(:,:) :: projected_hamil_work
        ! For parallel calculations, only the processor with label root performs the main
        real(dp), allocatable, dimension(:) :: temp_in, temp_out
        logical :: skip_calc
        logical :: t_store_subspace_basis
        integer :: max_subspace_size
    end type

    ! leaving these as global data, since direct CI currently only works for Davidson
    type(ras_vector), allocatable, dimension(:,:,:) :: direct_ci_inp, direct_ci_out

    contains

    subroutine InitHamiltonianCalc(this, print_info, hamil_type, max_subspace_size, t_store_subspace_basis)
        use direct_ci, only: create_ham_diag_direct_ci
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings
        use ras, only: find_ras_size
        use util_mod, only: int_fmt

        type(HamiltonianCalcType), intent(inout) :: this
        logical, intent(in) :: print_info
        integer :: i, mem_reqd, residual_mem_reqd, ierr, space_size
        integer, intent(in) :: hamil_type, max_subspace_size
        logical, intent(in) :: t_store_subspace_basis
        character(*), parameter :: t_r = "InitHamiltonianCalc"
        integer(MPIArg) :: mpi_temp
        real(dp), allocatable :: hamil_diag_temp(:)
       
        this%hamil_type = hamil_type
        this%t_store_subspace_basis = t_store_subspace_basis
        this%max_subspace_size = max_subspace_size

        associate(&
            space_size => this%space_size,&
            hfindex => this%hfindex&
        )

        ! allocate and define the hamiltonian diagonal, if not done so already.
        if (.not. allocated(hamil_diag)) then
            if (hamil_type == direct_ci_type) then
                allocate(hamil_diag(space_size), stat=ierr)
                call logmemalloc("hamil_diag", space_size, 8, t_r, hdiagtag, ierr)
                call create_ham_diag_direct_ci(davidson_ras, davidson_classes, davidson_strings, hamil_diag)
            else if (hamil_type == full_hamil_type) then
                space_size = size(hamiltonian,1)
                allocate(hamil_diag(space_size), stat=ierr)
                call logmemalloc("hamil_diag", space_size, 8, t_r, hdiagtag, ierr)
                do i = 1, space_size
                    hamil_diag(i) = hamiltonian(i,i)
                end do
            else if (hamil_type == parallel_sparse_hamil_type .or. hamil_type == sparse_hamil_type) then
                    ! in the case of sparse implementations, the diagonal should be
                    ! created when the hamiltonian itself is.
                    call stop_all("t_r", "the diagonal of the hamiltonian has not been allocated. &
                                         &cannot perform Davidson or Lanczos calculation.")
            end if
        end if

        space_size = size(hamil_diag)

        if (hamil_type == parallel_sparse_hamil_type) then

            allocate(this%space_sizes(0:nProcessors-1))
            allocate(this%partial_H_ket_disps(0:nProcessors-1))
            allocate(this%partial_H_ket(space_size))

            mpi_temp = int(space_size, MPIArg)
            call MPIAllGather(mpi_temp, this%space_sizes, ierr)
            ! The total space size across all processors.
            space_size = int(sum(this%space_sizes), sizeof_int)
            allocate(hamil_diag_temp(space_size))
            this%partial_H_ket_disps(0) = 0
            do i = 1, nProcessors-1
                this%partial_H_ket_disps(i) = sum(this%space_sizes(:i-1))
            end do
            call MPIGatherV(hamil_diag, hamil_diag_temp, this%space_sizes, this%partial_H_ket_disps, ierr)

            if (iProcIndex == root) then
                allocate(hamil_diag(space_size))
                hamil_diag = hamil_diag_temp
            end if
            deallocate(hamil_diag_temp)
        end if


        if (print_info) then
            write(6,'(1x,"number of determinants on this process:",'//int_fmt(space_size,1)//')') space_size; call neci_flush(6)
        end if

        if (iprocindex == root) then
            hfindex = maxloc((-hamil_diag),1)

            ! the memory required to allocate each of basis_vectors and
            ! multipied_basis_vectors, in mb.
            mem_reqd = max_subspace_size*space_size*8/1000000

            ! allocate the necessary arrays:
            if (t_store_subspace_basis) then
                safe_calloc(this%basis_vectors, (space_size, max_subspace_size), 0.0_dp)
                if (print_info) then
                    write(6,'(1x,"allocating array to hold subspace vectors (",'//int_fmt(mem_reqd,0)//',1x,"mb).")') mem_reqd
                    call neci_flush(6)
                endif
            endif
            safe_calloc(this%projected_hamil, (max_subspace_size,max_subspace_size), 0.0_dp)
            safe_calloc(this%projected_hamil_work, (max_subspace_size,max_subspace_size), 0.0_dp)
        else
            safe_malloc(this%temp_in, (space_size))
            safe_malloc(this%temp_out, (space_size))
        end if

        this%skip_calc = .false.
        if (print_info) then
            write(6,'(1x,"Hamiltonian calculation setup complete.",/)')
        endif
        call neci_flush(6)

        end associate

    end subroutine InitHamiltonianCalc

    subroutine DestroyHamiltonianCalc(this)
        type(HamiltonianCalcType), intent(inout) :: this
        safe_free(this%space_sizes)
        safe_free(this%partial_H_ket)
        safe_free(this%partial_H_ket_disps)
        safe_free(this%basis_vectors)
        safe_free(this%projected_hamil)
        safe_free(this%projected_hamil_work)
        safe_free(this%temp_in)
        safe_free(this%temp_out)
    end subroutine DestroyHamiltonianCalc

    subroutine multiply_hamil_and_vector(this, input_vector, output_vector)

        type(HamiltonianCalcType), intent(inout) :: this
        HElement_t(dp), intent(in) :: input_vector(:)
        HElement_t(dp), intent(out) :: output_vector(:)

        associate(hamil_type => this%hamil_type)

        if (hamil_type == full_hamil_type) then
            call multiply_hamil_and_vector_full(this, input_vector, output_vector)
        else if (hamil_type == sparse_hamil_type) then
            call mult_hamil_vector_sparse(this, input_vector, output_vector)
        else if (hamil_type == parallel_sparse_hamil_type) then
            call mult_hamil_vector_par_sparse(this, input_vector, output_vector)
        else if (hamil_type == direct_ci_type) then
            call mult_hamil_vector_direct_ci(input_vector, output_vector)
        end if

        end associate

    end subroutine multiply_hamil_and_vector

    subroutine multiply_hamil_and_vector_full(this, input_vector, output_vector)

        type(HamiltonianCalcType), intent(inout) :: this
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)

        ! This function performs y := alpha*a*x + beta*y
        ! N specifies not to use the transpose of a.
        ! space_size is the number of rows in a.
        ! space_size is the number of columns of a.
        ! alpha = 1.0_dp.
        ! a = hamiltonian.
        ! space_size is the first dimension of a.
        ! input x = input_vector.
        ! 1 is the increment of the elements of x.
        ! beta = 0.0_dp.
        ! output y = output_vector.
        ! 1 is the incremenet of the elements of y.

        associate(space_size => this%space_size)

        call dgemv('N', &
                   space_size, &
                   space_size, &
                   1.0_dp, &
                   hamiltonian, &
                   space_size, &
                   input_vector, &
                   1, &
                   0.0_dp, &
                   output_vector, &
                   1)

       end associate

    end subroutine multiply_hamil_and_vector_full

    subroutine mult_hamil_vector_sparse(this, input_vector, output_vector)

        type(HamiltonianCalcType), intent(inout) :: this
        HElement_t(dp), intent(in) :: input_vector(:)
        HElement_t(dp), intent(out) :: output_vector(:)
        integer :: i, j

        associate(space_size => this%space_size)
        output_vector = 0.0_dp

        do i = 1, space_size
            do j = 1, sparse_ham(i)%num_elements
                output_vector(i) = output_vector(i) + sparse_ham(i)%elements(j)*input_vector(sparse_ham(i)%positions(j))
            end do
        end do

        end associate

    end subroutine mult_hamil_vector_sparse

    subroutine mult_hamil_vector_par_sparse(this, input_vector, output_vector)

        type(HamiltonianCalcType), intent(inout) :: this
        HElement_t(dp), intent(in) :: input_vector(:)
        HElement_t(dp), intent(out) :: output_vector(:)
        integer :: i, j, ierr

        ! Use output_vector as temporary space.
        output_vector = input_vector

        call MPIBarrier(ierr, tTimeIn=.false.)

        call MPIBCast(output_vector)

        this%partial_H_ket = 0.0_dp

        do i = 1, this%space_sizes(iProcIndex)
            do j = 1, sparse_ham(i)%num_elements
                this%partial_H_ket(i) = this%partial_H_ket(i) + &
                    sparse_ham(i)%elements(j)*output_vector(sparse_ham(i)%positions(j))
            end do
        end do

        call MPIGatherV(this%partial_H_ket, output_vector, this%space_sizes, this%partial_H_ket_disps, ierr)

    end subroutine mult_hamil_vector_par_sparse

    subroutine mult_hamil_vector_direct_ci(input_vector, output_vector)

        use direct_ci, only: perform_multiplication, transfer_from_block_form, transfer_to_block_form
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits
        use SystemData, only: ecore

        HElement_t(dp), intent(in) :: input_vector(:)
        HElement_t(dp), intent(out) :: output_vector(:)

        ! The davidson code uses a single vector to store amplitudes. However, the direct CI code
        ! works in terms of alpha and beta strings and so uses block matrices. This routine will
        ! transfer the vector to block form.

        call transfer_to_block_form(davidson_ras, davidson_classes, input_vector, direct_ci_inp)

        call perform_multiplication(davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits, &
                direct_ci_inp, direct_ci_out)

        call transfer_from_block_form(davidson_ras, davidson_classes, output_vector, direct_ci_out)

        ! The above multiplication does not include the nuclear-nuclear energy, so add this
        ! contribution now.
        output_vector = output_vector + ecore*input_vector

    end subroutine mult_hamil_vector_direct_ci


end module hamiltonian_linalg
