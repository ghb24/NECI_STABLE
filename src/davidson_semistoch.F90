#include "macros.h"
module davidson_semistoch

    ! This module performs the Davidson method to find the ground state of a diagonally-
    ! dominant matrix. For details of the theory behind the method, see i.e:
    ! http://web.mit.edu/bolin/www/Project-Report-18.335J.pdf

    use constants
    use FciMCData, only: core_ham_diag, DavidsonTag
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg, MPIBarrier
    use Parallel_neci, only: MPIBCast, MPIGatherV, MPIAllGather
    use ParallelHelper, only: root
    use ras_data
    use sparse_arrays, only: sparse_core_ham, HDiagTag

    implicit none

    integer, parameter :: max_num_davidson_iters = 20
    real(dp), parameter :: residual_norm_target = 0.0000001_dp

    ! To cut down on the amount of global data, introduce a derived type to hold a Davidson session
    type davidson_ss
        ! Total space size.
        integer :: space_size
        ! Space size on this process.
        integer :: space_size_this_proc
        ! Sizes on each process.
        integer(MPIArg), allocatable :: sizes(:)
        ! Displacements of each section of the vector across processes
        integer(MPIArg), allocatable :: displs(:)
        ! All algorithms for solving large eigenproblems involve a unitary rotation of the 
        ! Hamiltonian into a smaller basis. basis_vectors(:,i) is the ith such unit vector
        HElement_t(dp), allocatable, dimension(:,:) :: basis_vectors
        ! This array stores the basis vectors multiplied by H in its columns, i.e.
        ! multiplied_basis_vectors(:,1) = H*basis_vector(:,1).
        real(dp), allocatable :: multiplied_basis_vectors(:,:)
        ! By diagonalising the projected Hamiltonian we get an estimate at the ground state in
        ! the basis of those basis vectors stored in the basis_vectors array. davidson_eigenvector
        ! stores this same state, but in the *original* basis set. It therefore has a dimension
        ! the same size as the vector space.
        real(dp), allocatable :: davidson_eigenvector(:)
        ! This array holds the components of davidson_eigenvector in the basis of Krylov vectors.
        real(dp), allocatable :: eigenvector_proj(:)
        ! The residual is defined as r = H*v - E*v, where H is the Hamiltonian matrix, v is the
        ! ground state estimate (stored in davidson_eigenvector) and E is the corresponding
        ! energy eigenvalue. If v is an exact eigenstate then all the components of the residual
        ! are zero.
        real(dp), allocatable :: residual(:)
        ! As noted above, if davidson_eigenvector holds an exact eigenstate then the residual
        ! will have all zero components and this norm (the standard Euclidean norm) will be zero.
        ! Hence it is a measure of how converged the solution is.
        real(dp) :: residual_norm
        real(dp) :: davidson_eigenvalue
        ! temp vectors for real matrix-vector calculations even when compiling in complex mode
        real(dp), allocatable :: temp_in(:), temp_out(:)
        ! the hamiltonian projected into basis_vectors
        real(dp), allocatable :: projected_hamil(:,:)
        ! we'll usually need some working space for diagonalisation of H in the small basis
        real(dp), allocatable :: projected_hamil_work(:,:)
        ! For parallel calculations, this vector is the size of the space on this processor. This
        ! vector is used to store the output of H |ket> on this processor.
        HElement_t(dp), allocatable :: partial(:)
    end type davidson_ss

    interface multiply_hamil_and_vector_ss
        module procedure mult_ham_vector_real_ss
        module procedure mult_ham_vector_complex_ss
    end interface

    contains

    subroutine perform_davidson_ss(this, print_info_in)

        logical, intent(in) :: print_info_in
        logical :: print_info
        integer :: i
        real(dp) :: start_time, end_time
        type(davidson_ss), intent(inout) :: this

        ! Only let the root processor print information.
        print_info = print_info_in .and. (iProcIndex == root)
        
        call init_davidson_ss(this, print_info)

        if (print_info) write(6,'(1X,"Iteration",4X,"Residual norm",12X,"Energy",7X,"Time")'); call neci_flush(6)

        do i = 2, max_num_davidson_iters

            start_time = MPI_WTIME()

            if (iProcIndex == root) call subspace_expansion_ss(this, i)

            call project_hamiltonian_ss(this, i)

            if (iProcIndex == root) call subspace_extraction_ss(this, i)

            call calculate_residual_ss(this, i)

            call calculate_residual_norm_ss(this)

            end_time = MPI_WTIME()

            if (print_info) write(6,'(8X,i2,3X,f14.9,2x,f16.10,2x,f9.3)') i-1, this%residual_norm, &
                this%davidson_eigenvalue, end_time-start_time; call neci_flush(6)

            if (this%residual_norm < residual_norm_target) exit

        end do

        if (print_info) write(6,'(/,1x,"Final calculated energy:",1X,f16.10)') this%davidson_eigenvalue

        call free_davidson_ss(this)

    end subroutine perform_davidson_ss

    subroutine init_davidson_ss(this, print_info)
    
        ! This subroutine initialises the Davdison method by allocating the necessary arrays,
        ! defining the initial basis vector and projected Hamiltonian, and setting an initial
        ! guess at the ground state eigenvalue. It also calculates the corresponding residual
        ! which is needed to expand the space.

        use util_mod, only: int_fmt
        
        type(davidson_ss), intent(inout) :: this

        logical, intent(in) :: print_info

        integer :: i, hfindex, mem_reqd, residual_mem_reqd, ierr
        real(dp), allocatable :: hamil_diag_temp(:)
        logical :: skip_calc
        integer(MPIArg) :: mpi_temp
        character (len=*), parameter :: t_r = "init_davidson_ss"

        associate( &
            davidson_eigenvalue => this%davidson_eigenvalue, &
            space_size => this%space_size, &
            space_size_this_proc => this%space_size_this_proc &
        )

        space_size_this_proc = size(core_ham_diag)

        allocate(this%displs(0:nProcessors-1))
        allocate(this%sizes(0:nProcessors-1))
        allocate(this%partial(space_size_this_proc))

        mpi_temp = int(space_size_this_proc, MPIArg)
        call MPIAllGather(mpi_temp, this%sizes, ierr)
        ! The total space size across all processors.
        space_size = int(sum(this%sizes), sizeof_int)

        this%displs(0) = 0
        do i = 1, nProcessors-1
            this%displs(i) = sum(this%displs(:i-1))
        end do

        ! if a davidson calculation has already been performed, this array might still be
        ! allocated, so check!
        if (allocated(this%davidson_eigenvector)) then
            deallocate(this%davidson_eigenvector, stat=ierr)
            call logmemdealloc(t_r, davidsontag, ierr)
        end if
        safe_calloc_e(this%davidson_eigenvector, (space_size), 0.0_dp, ierr)
        call logmemalloc("davidson_eigenvector", space_size, 8, t_r, davidsontag, ierr)

        allocate(hamil_diag_temp(space_size))
        call MPIGatherV(core_ham_diag, hamil_diag_temp, this%sizes, this%displs, ierr)

        ! if there is only one state in the space being diagonalised:
        if (space_size == 1) then
            this%davidson_eigenvector(1) = 1.0_dp
            if (iprocindex == root) davidson_eigenvalue = hamil_diag_temp(1)
            call mpibcast(davidson_eigenvalue)
            skip_calc = .true.
            return
        end if

        safe_calloc(this%projected_hamil, (max_num_davidson_iters, max_num_davidson_iters), 0.0_dp)
        safe_calloc(this%projected_hamil_work, (max_num_davidson_iters, max_num_davidson_iters), 0.0_dp)

        if (iprocindex == root) then
            hfindex = maxloc((-hamil_diag_temp),1)
            ! the memory required to allocate each of basis_vectors and
            ! multipied_basis_vectors, in mb.
            mem_reqd = (max_num_davidson_iters*space_size*8)/1000000
            ! the memory required to allocate residual.
            residual_mem_reqd = space_size*8/1000000

            ! allocate the necessary arrays:
            if (print_info) then
                write(6,'(1x,"allocating array to hold subspace vectors (",'//int_fmt(mem_reqd,0)//',1x,"mb).")') mem_reqd
                call neci_flush(6)
            endif
            safe_calloc(this%basis_vectors, (space_size, max_num_davidson_iters), 0.0_dp)

            if (print_info) then
                write(6,'(1x,"allocating array to hold multiplied krylov vectors (",'&
                //int_fmt(mem_reqd,0)//',1x,"mb).")') mem_reqd
                call neci_flush(6)
            end if
            safe_calloc(this%multiplied_basis_vectors, (space_size, max_num_davidson_iters), 0.0_dp)

            if (print_info) then
                write(6,'(1x,"allocating array to hold the residual vector (",'&
                //int_fmt(residual_mem_reqd,0)//',1x,"mb).",/)') residual_mem_reqd
                call neci_flush(6)
            end if
            safe_calloc(this%residual, (space_size), 0.0_dp)

            safe_calloc(this%eigenvector_proj, (max_num_davidson_iters), 0.0_dp)

            ! for the initial basis vector, choose the hartree-fock state:
            !this%super%basis_vectors(hfindex, 1) = 1.0_dp
            this%basis_vectors(hfindex, 1) = 1.0_dp
            ! choose the hartree-fock state as the initial guess at the ground state, too.
            this%eigenvector_proj(1) = 1.0_dp
            this%davidson_eigenvector(hfindex) = 1.0_dp

            ! fill in the projected hamiltonian so far.
            !this%super%projected_hamil(1,1) = hamil_diag(hfindex)
            this%projected_hamil(1,1) = hamil_diag_temp(hfindex)
            ! take the initial eigenvalue to be the hartree-fock energy minus some small
            ! amount. this value cannot be exactly the hartree-fock energy, as this will
            ! result in dividing by zero in the subspace expansion step.
            davidson_eigenvalue = hamil_diag_temp(hfindex) - 0.001_dp
        else
            safe_malloc(this%temp_in, (space_size))
            safe_malloc(this%temp_out, (space_size))
        end if

        if (print_info) write(6,'(1x,"calculating the initial residual vector...")',advance='no'); call neci_flush(6)

        ! check that multiplying the initial vector by the hamiltonian doesn't give back
        ! the same vector. if it does then the initial vector (the hf determinant) is
        ! the ground state, so just keep that and exit the calculation.
        ! also, the result of the multiplied basis vector is used to calculate the
        ! initial residual vector, if the above condition is not true.
        skip_calc = .false.
        if (iprocindex == root) then
            call multiply_hamil_and_vector_ss(this%davidson_eigenvector, this%multiplied_basis_vectors(:,1), &
                                              this%partial, this%sizes, this%displs)
        else
            call multiply_hamil_and_vector_ss(this%davidson_eigenvector, this%temp_out, this%partial, this%sizes, this%displs)
        end if

        if (iprocindex == root) then
            if (all(abs(this%multiplied_basis_vectors(:,1)-core_ham_diag(hfindex)*this%davidson_eigenvector) < 1.0e-12_dp)) then
                skip_calc = .true.
                davidson_eigenvalue = core_ham_diag(hfindex)
            end if
        end if

        call mpibcast(skip_calc)
        if (skip_calc) return

        ! calculate the intial residual vector.
        call calculate_residual_ss(this, 1)
        call calculate_residual_norm_ss(this)

        if (print_info) write(6,'(1x,"done.",/)'); call neci_flush(6)

        end associate

    end subroutine init_davidson_ss

    subroutine subspace_expansion_ss(this, basis_index)

        type(davidson_ss), intent(inout) :: this
        integer, intent(in) :: basis_index
        integer :: i
        real(dp) :: dot_prod, norm

        ! Create the new basis state from the residual. This step performs
        ! t = (D - EI)^(-1) r,
        ! where D is the diagonal of the Hamiltonian matrix, E is the eigenvalue previously
        ! calculated, I is the identity matrix and r is the residual.

        do i = 1, this%space_size
            this%basis_vectors(i, basis_index) = this%residual(i)/(core_ham_diag(i) - this%davidson_eigenvalue)
        end do

        ! This step then maskes the new basis vector orthogonal to all other basis vectors, by doing
        ! t <- t - (t,v)v
        ! for each basis vector v, where (t,v) denotes the dot product.
        do i = 1, basis_index - 1
            dot_prod = dot_product(this%basis_vectors(:,basis_index), this%basis_vectors(:,i))
            this%basis_vectors(:, basis_index) = &
                this%basis_vectors(:, basis_index) - dot_prod*this%basis_vectors(:,i)
        end do

        ! Finally we calculate the norm of the new basis vector and then normalise it to have a norm of 1.
        ! The new basis vector is stored in the next available column in the basis_vectors array.
        norm = dot_product(this%basis_vectors(:,basis_index), this%basis_vectors(:,basis_index))
        norm = sqrt(norm)
        this%basis_vectors(:,basis_index) = this%basis_vectors(:,basis_index)/norm

    end subroutine subspace_expansion_ss


    subroutine subspace_extraction_ss(this, basis_index)

        type(davidson_ss), intent(inout) :: this
        integer, intent(in) :: basis_index
        integer :: lwork, info
        real(dp), allocatable, dimension(:) :: work
        real(dp) :: eigenvalue_list(basis_index)

        ! Scrap space for the diagonaliser.
        lwork = max(1,3*basis_index-1)
        allocate(work(lwork))

        ! This routine diagonalises a symmetric matrix, A.
        ! V tells the routine to calculate eigenvalues *and* eigenvectors.
        ! U tells the routine to get the upper half of A (it is symmetric).
        ! basis_index is the number of rows and columns in A.
        ! A = projected_hamil_work. This matrix stores the eigenvectors in its columns on output.
        ! basis_index is the leading dimension of A.
        ! eigenvalue_list stores the eigenvalues on output.
        ! work is scrap space.
        ! lwork is the length of the work array.
        ! info = 0 on output is diagonalisation is successful.
        call dsyev(&
            'V', &
            'U', &
            basis_index, &
            this%projected_hamil_work(1:basis_index,1:basis_index), &
            basis_index, &
            eigenvalue_list, &
            work, &
            lwork, &
            info &
        )

        this%davidson_eigenvalue = eigenvalue_list(1)
        ! The first column stores the ground state.
        this%eigenvector_proj(1:basis_index) = this%projected_hamil_work(1:basis_index,1)

        deallocate(work)

        ! eigenvector_proj stores the eigenstate in the basis of vectors stored in the array
        ! basis_vectors. We now want it in terms of the original basis. To get this, multiply
        ! eigenvector_proj by basis_vectors(:, 1:basis_index):

        ! This function performs y := alpha*a*x + beta*y
        ! N specifies not to use the transpose of a.
        ! space_size is the number of rows in a.
        ! basis_index is the number of columns of a.
        ! alpha = 1.0_dp.
        ! a = basis_vectors(:,1:basis_index).
        ! space_size is the first dimension of a.
        ! input x = eigenvector_proj(1:basis_index).
        ! 1 is the increment of the elements of x.
        ! beta = 0.0_dp.
        ! output y = davidson_eigenvector.
        ! 1 is the incremenet of the elements of y.
        call dgemv('N', &
                   this%space_size, &
                   basis_index, &
                   1.0_dp, &
                   this%basis_vectors(:,1:basis_index), &
                   this%space_size, &
                   this%eigenvector_proj(1:basis_index), &
                   1, &
                   0.0_dp, &
                   this%davidson_eigenvector, &
                   1)
        
    end subroutine subspace_extraction_ss

    subroutine project_hamiltonian_ss(this, basis_index)

        type(davidson_ss), intent(inout) :: this
        integer, intent(in) :: basis_index
        integer :: i

        if (iProcIndex == root) then
            ! Multiply the new basis_vector by the hamiltonian and store the result in
            ! multiplied_basis_vectors.
            call multiply_hamil_and_vector_ss(real(this%basis_vectors(:,basis_index), dp), &
                this%multiplied_basis_vectors(:,basis_index), this%partial, this%sizes, this%displs)

            ! Now multiply U^T by (H U) to find projected_hamil. The projected Hamiltonian will
            ! only differ in the new final column and row. Also, projected_hamil is symmetric.
            ! Hence, we only need to calculate the final column, and use this to update the final
            ! row also.
            do i = 1, basis_index
                this%projected_hamil(i, basis_index) = &
                    dot_product(this%basis_vectors(:, i), this%multiplied_basis_vectors(:,basis_index))
                this%projected_hamil(basis_index, i) = this%projected_hamil(i, basis_index)
            end do

            ! We will use the scrap Hamiltonian to pass into the diagonaliser later, since it
            ! overwrites this input matrix with the eigenvectors. Hence, make sure the scrap space
            ! stores the updated projected Hamiltonian.
            this%projected_hamil_work = this%projected_hamil
        else
            call multiply_hamil_and_vector_ss(this%temp_in, this%temp_out, this%partial, this%sizes, this%displs)
        end if

    end subroutine project_hamiltonian_ss

    subroutine mult_ham_vector_real_ss(input_vector, output_vector, partial, sizes, displs)

        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)
        real(dp), intent(out) :: partial(:)
        integer(MPIArg), intent(in) :: sizes(0:), displs(0:)

        integer :: i, j, ierr

        ! Use output_vector as temporary space.
        output_vector = input_vector

        call MPIBarrier(ierr, tTimeIn=.false.)

        call MPIBCast(output_vector)

        partial = 0.0_dp

        do i = 1, sizes(iProcIndex)
            do j = 1, sparse_core_ham(i)%num_elements
                partial(i) = partial(i) + sparse_core_ham(i)%elements(j)*output_vector(sparse_core_ham(i)%positions(j))
            end do
        end do

        ! this templated routine forbids mixing real and complex arrays
        call MPIGatherV(real(partial, dp), output_vector, sizes, displs, ierr)

    end subroutine mult_ham_vector_real_ss

    subroutine mult_ham_vector_complex_ss(input_vector, output_vector, partial, sizes, displs)

        complex(dp), intent(in) :: input_vector(:)
        complex(dp), intent(out) :: output_vector(:)
        HElement_t(dp), intent(out) :: partial(:)
        integer(MPIArg), intent(in) :: sizes(0:), displs(0:)

        integer :: i, j, ierr

        ! Use output_vector as temporary space.
        output_vector = input_vector

        call MPIBarrier(ierr, tTimeIn=.false.)

        call MPIBCast(output_vector)

        partial = 0.0_dp

        do i = 1, sizes(iProcIndex)
            do j = 1, sparse_core_ham(i)%num_elements
                partial(i) = partial(i) + sparse_core_ham(i)%elements(j)*output_vector(sparse_core_ham(i)%positions(j))
            end do
        end do
#ifdef __CMPLX
        call MPIGatherV(partial, output_vector, sizes, displs, ierr)
#else
        call MPIGatherV(cmplx(partial, 0.0_dp,dp), output_vector, sizes, displs, ierr)
#endif

    end subroutine mult_ham_vector_complex_ss

    subroutine calculate_residual_ss(this, basis_index)

        type(davidson_ss), intent(inout) :: this
        integer, intent(in) :: basis_index

        ! This routine calculates the residual, r, corresponding to the new estimate of the
        ! ground state, stored in davidson_eigenvector. This is defined as
        ! r = Hv - Ev,
        ! where H is the Hamiltonian, v is the ground state vector estimate and E is the 
        ! ground state energy estimate.

        if (iProcIndex == root) then 
            ! Calculate r = Hv - Ev:
            ! Note that, here, eigenvector_proj holds the components of v in the Krylov basis,
            ! and multiplied_basis_vectors holds the Krylov vectors multiplied by H, hence
            ! the matmul below does indeed retturn Hv.
            this%residual = matmul(this%multiplied_basis_vectors(:,1:basis_index), this%eigenvector_proj(1:basis_index))
            this%residual = this%residual - this%davidson_eigenvalue*this%davidson_eigenvector
        end if

    end subroutine calculate_residual_ss

    subroutine calculate_residual_norm_ss(this)

        type(davidson_ss), intent(inout) :: this

        ! This subroutine calculates the Euclidean norm of the reisudal vector, r:
        ! residual_norm^2 = \sum_i r_i^2
        if (iProcIndex == root) then
            this%residual_norm = sqrt(dot_product(this%residual, this%residual))
        end if

        !if (this%super%hamil_type == parallel_sparse_hamil_type) call MPIBCast(this%residual_norm)
        call MPIBCast(this%residual_norm)

    end subroutine calculate_residual_norm_ss

    subroutine free_davidson_ss(this)

        type(davidson_ss), intent(inout) :: this

        ! we are now done with these arrays
        safe_free(this%basis_vectors)
        safe_free(this%multiplied_basis_vectors)
        safe_free(this%residual)
        safe_free(this%eigenvector_proj)
        safe_free(this%temp_in)
        safe_free(this%temp_out)
        ! but keep the davidson eigenvector

    end subroutine free_davidson_ss

    subroutine destroy_davidson_ss(this)

        type(davidson_ss), intent(inout) :: this

        ! deallocate the davidson vector as well
        call free_davidson_ss(this)
        safe_free(this%davidson_eigenvector)

    end subroutine destroy_davidson_ss

end module davidson_semistoch
