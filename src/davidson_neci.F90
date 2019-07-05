#include "macros.h"
module davidson_neci

    ! This module performs the Davidson method to find the ground state of a diagonally-
    ! dominant matrix. For details of the theory behind the method, see i.e:
    ! http://web.mit.edu/bolin/www/Project-Report-18.335J.pdf

    use constants
    use FciMCData, only: hamiltonian, DavidsonTag
    use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg, MPIBarrier
    use Parallel_neci, only: MPIBCast, MPIGatherV, MPIAllGather
    use ParallelHelper, only: root
    use ras_data
    use sparse_arrays, only: sparse_ham, hamil_diag, HDiagTag
    use hamiltonian_linalg, only: &
        full_hamil_type, &
        sparse_hamil_type, &
        parallel_sparse_hamil_type, &
        direct_ci_type, &
        HamiltonianCalcType, &
        initHamiltonianCalc, &
        multiply_hamil_and_vector, &
        direct_ci_inp, &
        direct_ci_out

    implicit none

    integer, parameter :: max_num_davidson_iters = 25
    real(dp), parameter :: residual_norm_target = 0.0000001_dp

    ! To cut down on the amount of global data, introduce a derived type to hold a Davidson session
    type DavidsonCalcType
        ! "super type"
        type(HamiltonianCalcType) :: super
        ! This array stores the basis vectors multiplied by H in its columns, i.e.
        ! multiplied_basis_vectors(:,1) = H*basis_vector(:,1).
        real(dp), allocatable, dimension(:,:) :: multiplied_basis_vectors
        ! By diagonalising the projected Hamiltonian we get an estimate at the ground state in
        ! the basis of those basis vectors stored in the basis_vectors array. davidson_eigenvector
        ! stores this same state, but in the *original* basis set. It therefore has a dimension
        ! the same size as the vector space.
        real(dp), allocatable, dimension(:) :: davidson_eigenvector
        ! This array holds the components of davidson_eigenvector in the basis of Krylov vectors.
        real(dp), allocatable, dimension(:) :: eigenvector_proj
        ! The residual is defined as r = H*v - E*v, where H is the Hamiltonian matrix, v is the
        ! ground state estimate (stored in davidson_eigenvector) and E is the corresponding
        ! energy eigenvalue. If v is an exact eigenstate then all the components of the residual
        ! are zero.
        real(dp), allocatable, dimension(:) :: residual
        ! As noted above, if davidson_eigenvector holds an exact eigenstate then the residual
        ! will have all zero components and this norm (the standard Euclidean norm) will be zero.
        ! Hence it is a measure of how converged the solution is.
        real(dp) :: residual_norm
        real(dp) :: davidson_eigenvalue
        ! temp vectors for real matrix-vector calculations even when compiling in complex mode
        real(dp), allocatable, dimension(:) :: temp_in, temp_out
    end type


    contains

    subroutine perform_davidson(this, hamil_type_in, print_info_in)
      use mpi
        integer, intent(in) :: hamil_type_in
        logical, intent(in) :: print_info_in
        logical :: print_info
        integer :: i
        real(dp) :: start_time, end_time
        type(DavidsonCalcType), intent(inout) :: this

        ! Only let the root processor print information.
        print_info = print_info_in .and. (iProcIndex == root)
        
        call InitDavidsonCalc(this, print_info, hamil_type_in)

        if (print_info) write(6,'(1X,"Iteration",4X,"Residual norm",12X,"Energy",7X,"Time")'); call neci_flush(6)

        do i = 2, max_num_davidson_iters

            if (this%super%skip_calc) exit

            start_time = MPI_WTIME()

            if (iProcIndex == root) call subspace_expansion(this, i)

            call project_hamiltonian(this, i)

            if (iProcIndex == root) call subspace_extraction(this, i)

            call calculate_residual(this, i)

            call calculate_residual_norm(this)

            end_time = MPI_WTIME()

            if (print_info) write(6,'(8X,i2,3X,f14.9,2x,f16.10,2x,f9.3)') i-1, this%residual_norm, &
                this%davidson_eigenvalue, end_time-start_time; call neci_flush(6)

            if (this%residual_norm < residual_norm_target) exit

        end do

        if (print_info) write(6,'(/,1x,"Final calculated energy:",1X,f16.10)') this%davidson_eigenvalue

        call FreeDavidsonCalc(this)

    end subroutine perform_davidson

    subroutine InitDavidsonCalc(this, print_info, hamil_type)
    
        ! This subroutine initialises the Davdison method by allocating the necessary arrays,
        ! defining the initial basis vector and projected Hamiltonian, and setting an initial
        ! guess at the ground state eigenvalue. It also calculates the corresponding residual
        ! which is needed to expand the space.

        use direct_ci, only: create_ham_diag_direct_ci
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings
        use ras, only: find_ras_size
        use util_mod, only: int_fmt
        
        type(DavidsonCalcType), intent(inout) :: this

        logical, intent(in) :: print_info
        integer, intent(in) :: hamil_type

        logical :: skip_calc
        integer :: i, mem_reqd, residual_mem_reqd, ierr
        integer(mpiarg) :: mpi_temp
        real(dp), allocatable :: hamil_diag_temp(:)
        character (len=*), parameter :: t_r = "init_davidson"

        call InitHamiltonianCalc(this%super, print_info, hamil_type, max_num_davidson_iters, .true., .false.)

        associate( &
            davidson_eigenvalue => this%davidson_eigenvalue, &
            space_size => this%super%space_size, &
            hfindex => this%super%hfindex &
        )

        ! if a davidson calculation has already been performed, this array might still be
        ! allocated, so check!
        if (allocated(this%davidson_eigenvector)) then
            deallocate(this%davidson_eigenvector, stat=ierr)
            call logmemdealloc(t_r, davidsontag, ierr)
        end if
        safe_calloc_e(this%davidson_eigenvector, (space_size), 0.0_dp, ierr)
        call logmemalloc("davidson_eigenvector", space_size, 8, t_r, davidsontag, ierr)

        ! if there is only one state in the space being diagonalised:
        if (space_size == 1) then
            this%davidson_eigenvector(1) = 1.0_dp
            if (iprocindex == root) davidson_eigenvalue = hamil_diag(1)
            call mpibcast(davidson_eigenvalue)
            skip_calc = .true.
            call mpibcast(skip_calc)
            this%super%skip_calc = skip_calc
            return
        end if

        if (iprocindex == root) then
            hfindex = maxloc((-hamil_diag),1)
            ! the memory required to allocate each of basis_vectors and
            ! multipied_basis_vectors, in mb.
            mem_reqd = (max_num_davidson_iters*space_size*8)/1000000
            ! the memory required to allocate residual.
            residual_mem_reqd = space_size*8/1000000

            ! allocate the necessary arrays:
            if (print_info) write(6,'(1x,"allocating array to hold multiplied krylov vectors (",'&
                //int_fmt(mem_reqd,0)//',1x,"mb).")') mem_reqd; call neci_flush(6)
            safe_calloc(this%multiplied_basis_vectors, (space_size, max_num_davidson_iters), 0.0_dp)
            safe_calloc(this%eigenvector_proj, (max_num_davidson_iters), 0.0_dp)
            if (print_info) write(6,'(1x,"allocating array to hold the residual vector (",'&
                //int_fmt(residual_mem_reqd,0)//',1x,"mb).",/)') residual_mem_reqd; call neci_flush(6)
            safe_calloc(this%residual, (space_size), 0.0_dp)

            ! for the initial basis vector, choose the hartree-fock state:
            this%super%basis_vectors(hfindex, 1) = 1.0_dp
            ! choose the hartree-fock state as the initial guess at the ground state, too.
            this%eigenvector_proj(1) = 1.0_dp
            this%davidson_eigenvector(hfindex) = 1.0_dp

            ! fill in the projected hamiltonian so far.
            this%super%projected_hamil(1,1) = hamil_diag(hfindex)
            ! take the initial eigenvalue to be the hartree-fock energy minus some small
            ! amount. this value cannot be exactly the hartree-fock energy, as this will
            ! result in dividing by zero in the subspace expansion step.
            davidson_eigenvalue = hamil_diag(hfindex) - 0.001_dp
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
            call multiply_hamil_and_vector(this%super, this%davidson_eigenvector, this%multiplied_basis_vectors(:,1))
        else
            call multiply_hamil_and_vector(this%super, this%davidson_eigenvector, this%temp_out)
        end if
        if (iprocindex == root) then
            if (all(abs(this%multiplied_basis_vectors(:,1)-hamil_diag(hfindex)*this%davidson_eigenvector) < 1.0e-12_dp)) then
                skip_calc = .true.
                davidson_eigenvalue = hamil_diag(hfindex)
            end if
        end if
        if (hamil_type == parallel_sparse_hamil_type) call mpibcast(skip_calc)
        this%super%skip_calc = skip_calc
        if (skip_calc) return
        ! calculate the intial residual vector.
        call calculate_residual(this, 1)
        call calculate_residual_norm(this)

        if (print_info) write(6,'(1x,"done.",/)'); call neci_flush(6)

        end associate

    end subroutine InitDavidsonCalc

    subroutine subspace_expansion(this, basis_index)
        type(DavidsonCalcType), intent(inout) :: this
        integer, intent(in) :: basis_index
        integer :: i
        real(dp) :: dot_prod, norm

        ! Create the new basis state from the residual. This step performs
        ! t = (D - EI)^(-1) r,
        ! where D is the diagonal of the Hamiltonian matrix, E is the eigenvalue previously
        ! calculated, I is the identity matrix and r is the residual.

        do i = 1, this%super%space_size
            this%super%basis_vectors(i, basis_index) = this%residual(i)/(hamil_diag(i) - this%davidson_eigenvalue)
        end do

        ! This step then maskes the new basis vector orthogonal to all other basis vectors, by doing
        ! t <- t - (t,v)v
        ! for each basis vector v, where (t,v) denotes the dot product.
        do i = 1, basis_index - 1
            dot_prod = dot_product(this%super%basis_vectors(:,basis_index), this%super%basis_vectors(:,i))
            this%super%basis_vectors(:, basis_index) = &
                this%super%basis_vectors(:, basis_index) - dot_prod*this%super%basis_vectors(:,i)
        end do

        ! Finally we calculate the norm of the new basis vector and then normalise it to have a norm of 1.
        ! The new basis vector is stored in the next available column in the basis_vectors array.
        norm = dot_product(this%super%basis_vectors(:,basis_index), this%super%basis_vectors(:,basis_index))
        norm = sqrt(norm)
        this%super%basis_vectors(:,basis_index) = this%super%basis_vectors(:,basis_index)/norm

    end subroutine subspace_expansion


    subroutine subspace_extraction(this, basis_index)

        type(DavidsonCalcType), intent(inout) :: this
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
            this%super%projected_hamil_work(1:basis_index,1:basis_index), &
            basis_index, &
            eigenvalue_list, &
            work, &
            lwork, &
            info &
        )

        this%davidson_eigenvalue = eigenvalue_list(1)
        ! The first column stores the ground state.
        this%eigenvector_proj(1:basis_index) = this%super%projected_hamil_work(1:basis_index,1)

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
                   this%super%space_size, &
                   basis_index, &
                   1.0_dp, &
                   this%super%basis_vectors(:,1:basis_index), &
                   this%super%space_size, &
                   this%eigenvector_proj(1:basis_index), &
                   1, &
                   0.0_dp, &
                   this%davidson_eigenvector, &
                   1)
        
    end subroutine subspace_extraction

    subroutine project_hamiltonian(this, basis_index)

        type(DavidsonCalcType), intent(inout) :: this
        integer, intent(in) :: basis_index
        integer :: i

        if (iProcIndex == root) then
            ! Multiply the new basis_vector by the hamiltonian and store the result in
            ! multiplied_basis_vectors.
            call multiply_hamil_and_vector(this%super, real(this%super%basis_vectors(:,basis_index), dp), &
                this%multiplied_basis_vectors(:,basis_index))

            ! Now multiply U^T by (H U) to find projected_hamil. The projected Hamiltonian will
            ! only differ in the new final column and row. Also, projected_hamil is symmetric.
            ! Hence, we only need to calculate the final column, and use this to update the final
            ! row also.
            do i = 1, basis_index
                this%super%projected_hamil(i, basis_index) = &
                    dot_product(this%super%basis_vectors(:, i), this%multiplied_basis_vectors(:,basis_index))
                this%super%projected_hamil(basis_index, i) = this%super%projected_hamil(i, basis_index)
            end do

            ! We will use the scrap Hamiltonian to pass into the diagonaliser later, since it
            ! overwrites this input matrix with the eigenvectors. Hence, make sure the scrap space
            ! stores the updated projected Hamiltonian.
            this%super%projected_hamil_work = this%super%projected_hamil
        else
            call multiply_hamil_and_vector(this%super, this%temp_in, this%temp_out)
        end if

    end subroutine project_hamiltonian

    subroutine calculate_residual(this, basis_index)

        type(DavidsonCalcType), intent(inout) :: this
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

    end subroutine calculate_residual

    subroutine calculate_residual_norm(this)

        type(DavidsonCalcType), intent(inout) :: this
        ! This subroutine calculates the Euclidean norm of the reisudal vector, r:
        ! residual_norm^2 = \sum_i r_i^2
        if (iProcIndex == root) then
            this%residual_norm = sqrt(dot_product(this%residual, this%residual))
        end if

        if (this%super%hamil_type == parallel_sparse_hamil_type) call MPIBCast(this%residual_norm)

    end subroutine calculate_residual_norm

    subroutine FreeDavidsonCalc(this)
        use hamiltonian_linalg, only: DestroyHamiltonianCalc
        type(DavidsonCalcType), intent(inout) :: this
        ! destroy the super type instance
        call DestroyHamiltonianCalc(this%super)
        ! we are now done with these arrays
        safe_free(this%multiplied_basis_vectors)
        safe_free(this%residual)
        safe_free(this%eigenvector_proj)
        safe_free(this%temp_in)
        safe_free(this%temp_out)
        ! but keep the davidson eigenvector
    end subroutine FreeDavidsonCalc

    subroutine DestroyDavidsonCalc(this)
        type(DavidsonCalcType), intent(inout) :: this
        ! deallocate the davidson vector as well
        call FreeDavidsonCalc(this)
        safe_free(this%davidson_eigenvector)
    end subroutine DestroyDavidsonCalc

    function davidson_direct_ci_init(print_info_in) result (this)
        use bit_rep_data, only: NIfD
        use direct_ci, only: create_direct_ci_arrays
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits
        use ras, only: initialise_ras_space, find_ras_size

        type(DavidsonCalcType) :: this
        logical, intent(in) :: print_info_in
        integer :: class_i, class_j, j, sym_i, sym_j

        write(6,'(/,1X,"Beginning Direct CI Davidson calculation.",/)'); call neci_flush(6)


        call initialise_ras_space(davidson_ras, davidson_classes)
        ! The total hilbert space dimension of calculation to be performed.
        call find_ras_size(davidson_ras, davidson_classes, this%super%space_size)

        allocate(davidson_strings(-1:tot_nelec, davidson_ras%num_strings))
        allocate(davidson_iluts(0:NIfD, davidson_ras%num_strings))
        allocate(davidson_excits(davidson_ras%num_strings))

        ! Create the arrays used by the direct CI multiplication.
        call create_direct_ci_arrays(davidson_ras, davidson_classes, davidson_strings, &
                davidson_iluts, davidson_excits)

        ! Allocate input and output direct CI vectors.
        allocate(direct_ci_inp(size(davidson_classes),size(davidson_classes),0:7))
        allocate(direct_ci_out(size(davidson_classes),size(davidson_classes),0:7))
        do class_i = 1, size(davidson_classes)
            do j = 1, davidson_classes(class_i)%num_comb
                class_j = davidson_classes(class_i)%allowed_combns(j)
                do sym_i = 0, 7
                    sym_j = ieor(int(HFSym_ras,sizeof_int), sym_i)
                    if (davidson_classes(class_i)%num_sym(sym_i) == 0) cycle
                    if (davidson_classes(class_j)%num_sym(sym_j) == 0) cycle
                    allocate(direct_ci_inp(class_i,class_j,sym_i)%&
                        elements(1:davidson_classes(class_i)%num_sym(sym_i),1:davidson_classes(class_j)%num_sym(sym_j)))
                    allocate(direct_ci_out(class_i,class_j,sym_i)%&
                        elements(1:davidson_classes(class_i)%num_sym(sym_i),1:davidson_classes(class_j)%num_sym(sym_j)))
                end do
            end do
        end do
    end function davidson_direct_ci_init

    subroutine davidson_direct_ci_end(this)
        type(DavidsonCalcType), intent(inout) :: this
        integer :: ierr
        if (allocated(hamil_diag)) then
            deallocate(hamil_diag, stat=ierr)
            call LogMemDealloc("davidson_direct_ci_end", HDiagTag, ierr)
        end if
        if (allocated(this%davidson_eigenvector)) then
            deallocate(this%davidson_eigenvector, stat=ierr)
            call LogMemDealloc("davidson_direct_ci_end", DavidsonTag, ierr)
        end if

        write(6,'(/,1X,"Direct CI Davidson calculation complete.",/)'); call neci_flush(6)

        write(6,"(1X,a10,f16.10)") "GROUND E =", this%davidson_eigenvalue; call neci_flush(6)
    end subroutine davidson_direct_ci_end

end module davidson_neci
