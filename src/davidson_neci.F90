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

implicit none

integer :: max_num_davidson_iters = 25
real(dp) :: residual_norm_target = 0.00000001_dp

integer :: hamil_type
! The value of hamil_type specifies what form the Hamiltonian is stored in.
! The following options are currently available:
integer :: full_hamil_type = 1
integer :: sparse_hamil_type = 2
integer :: parallel_sparse_hamil_type = 3
integer :: direct_ci_type = 4

! The dimension of the vector space we are working in, as determined by the number
! of rows and columns in the Hamiltonian matrix.
integer :: space_size
! This array stores the basis vectors used in its columns, i.e. basis_vector(:,1) stores
! the components of the first basis vector.
real(dp), allocatable, dimension(:,:) :: basis_vectors
! The projected Hamiltonian is H_p = U^T H U, where U is the array of basis vectors.
real(dp), allocatable, dimension(:,:) :: projected_hamil
! This is used as scrap space for the projected Hamiltonian.
real(dp), allocatable, dimension(:,:) :: projected_hamil_scrap
! By diagonalising the projected Hamiltonian we get an estimate at the ground state in
! the basis of those basis vectors stored in the basis_vectors array. davidson_eigenvector
! stores this same state, but in the *original* basis set. It therefore has a dimension
! the same size as the vector space.
real(dp), allocatable, dimension(:) :: davidson_eigenvector
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

! For parallel calculations, this vector is the size of the space on this processor. This
! vector is used to store the output of multiplication by the Hamiltonian on this processor.
real(dp), allocatable, dimension(:) :: partial_davidson_vector
! For parallel calculations, only the processor with label root performs the main
! davidson calculation. These vectors are used as temporary space for the other processors.
real(dp), allocatable, dimension(:) :: temp_in, temp_out

! For parallel calculations, store all spaces sizes on each processor, and the
! displacements necessary for communication.
integer(MPIArg), allocatable, dimension(:) :: space_sizes, davidson_disps

type(ras_vector), allocatable, dimension(:,:,:) :: direct_ci_inp, direct_ci_out

integer(TagIntType) :: ResidualTag

    contains

    subroutine perform_davidson(input_hamil_type, print_info)

        integer, intent(in) :: input_hamil_type
        logical, intent(in) :: print_info
        logical :: tSkipCalc
        integer :: i

        hamil_type = input_hamil_type

        call init_davidson(tSkipCalc)

        do i = 2, max_num_davidson_iters

            if (tSkipCalc) exit

            if (iProcIndex == root) call subspace_expansion(i)

            call project_hamiltonian(i)

            if (iProcIndex == root) call subspace_extraction(i)

            call calculate_residual()

            call calculate_residual_norm()

            if (print_info) write(6,'(a10,1X,i2,5X,a14,1X,f12.10,5x,a7,1x,f15.10)') &
                "Iteration:", i-1, "residual norm:", residual_norm, "energy:", davidson_eigenvalue

            if (residual_norm < residual_norm_target) exit
            
        end do

        if (print_info) write(6,'(a24,1X,f14.10)') "Final calculated energy:", davidson_eigenvalue

        call end_davidson()

    end subroutine perform_davidson

    subroutine init_davidson(tSkipCalc)
    
        ! This subroutine initialises the Davdison method by allocating the necessary arrays,
        ! defining the initial basis vector and projected Hamiltonian, and setting an initial
        ! guess at the ground state eigenvalue. It also calculates the corresponding residual
        ! which is needed to expand the space.

        use direct_ci, only: create_ham_diag_direct_ci
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings
        use ras, only: find_ras_size

        logical, intent(out) :: tSkipCalc

        integer :: i, HFindex, ierr
        integer(MPIArg) :: mpi_temp
        real(dp), allocatable, dimension(:) :: hamil_diag_temp, test_vector(:)
        character (len=*), parameter :: t_r = "init_davidson"

        ! Allocate and define the Hamiltonian diagonal, if not done so already.
        if (.not. allocated(hamil_diag)) then
            if (hamil_type == direct_ci_type) then
                allocate(hamil_diag(space_size), stat=ierr)
                call LogMemAlloc("hamil_diag", space_size, 8, t_r, HDiagTag, ierr)
                call create_ham_diag_direct_ci(davidson_ras, davidson_classes, davidson_strings, hamil_diag)
            else if (hamil_type == full_hamil_type) then
                space_size = size(hamiltonian,1)
                allocate(hamil_diag(space_size), stat=ierr)
                call LogMemAlloc("hamil_diag", space_size, 8, t_r, HDiagTag, ierr)
                do i = 1, space_size
                    hamil_diag(i) = hamiltonian(i,i)
                end do
            else if (hamil_type == parallel_sparse_hamil_type .or. &
                     hamil_type == sparse_hamil_type) then
                    ! In the case of sparse implementations, the diagonal should be
                    ! created when the Hamiltonian itself is.
                    call stop_all("t_r", "The diagonal of the Hamiltonian has not been allocated. Cannot perform &
                                         &Davidson calculation.")
            end if
        end if

        space_size = size(hamil_diag)

        if (hamil_type == parallel_sparse_hamil_type) then
            allocate(partial_davidson_vector(space_size))
            allocate(space_sizes(0:nProcessors-1))
            allocate(davidson_disps(0:nProcessors-1))
            mpi_temp = int(space_size, MPIArg)
            call MPIAllGather(mpi_temp, space_sizes, ierr)
            ! The total space size across all processors.
            space_size = int(sum(space_sizes), sizeof_int)
            allocate(hamil_diag_temp(space_size))
            davidson_disps(0) = 0
            do i = 1, nProcessors-1
                davidson_disps(i) = sum(space_sizes(:i-1))
            end do
            call MPIGatherV(hamil_diag, hamil_diag_temp, space_sizes, davidson_disps, ierr)

            if (iProcIndex == root) then
                deallocate(hamil_diag)
                allocate(hamil_diag(space_size))
                hamil_diag = hamil_diag_temp
            end if
            deallocate(hamil_diag_temp)
        end if

        ! If a davidson calculation has already been performed, this array might still be
        ! allocated, so check!
        if (allocated(davidson_eigenvector)) then
            deallocate(davidson_eigenvector, stat=ierr)
            call LogMemDealloc(t_r, DavidsonTag, ierr)
        end if
        allocate(davidson_eigenvector(space_size), stat=ierr)
        call LogMemAlloc("davidson_eigenvector", space_size, 8, t_r, DavidsonTag, ierr)
        davidson_eigenvector = 0.0_dp

        if (iProcIndex == root) then
            HFindex = maxloc((-hamil_diag),1)

            ! Allocate the necessary arrays:
            allocate(basis_vectors(space_size, max_num_davidson_iters))
            allocate(projected_hamil(max_num_davidson_iters,max_num_davidson_iters))
            allocate(projected_hamil_scrap(max_num_davidson_iters,max_num_davidson_iters))
            allocate(residual(space_size))
            call LogMemAlloc("residual", space_size, 8, t_r, ResidualTag, ierr)
            projected_hamil = 0.0_dp
            basis_vectors = 0.0_dp
            residual = 0.0_dp

            ! If there is only one state in the space being diagonalised:
            if (space_size == 1) then
                davidson_eigenvector(1) = 1.0_dp
                davidson_eigenvalue = hamil_diag(1)
                max_num_davidson_iters = 0
                return
            end if

            ! For the initial basis vector, choose the Hartree-Fock state:
            basis_vectors(HFindex, 1) = 1.0_dp
            ! Choose the Hartree-Fock state as the initial guess at the ground state, too.
            davidson_eigenvector(HFindex) = 1.0_dp

            ! Fill in the projected Hamiltonian so far.
            projected_hamil(1,1) = hamil_diag(HFindex)
            ! Take the initial eigenvalue to be the Hartree-Fock energy minus some small
            ! amount. This value cannot be exactly the Hartree-Fock energy, as this will
            ! result in dividing by zero in the subspace expansion step.
            davidson_eigenvalue = hamil_diag(HFindex) - 0.001_dp

        else
            allocate(temp_in(space_size))
            allocate(temp_out(space_size))
        end if

        ! Check that multiplying the initial vector by the Hamiltonian doesn't give back
        ! the same vector. If it does then the initial vector (the HF determinant) is
        ! the ground state, so just keep that and exit the calculation.
        tSkipCalc = .false.
        allocate(test_vector(space_size))
        call multiply_hamil_and_vector(davidson_eigenvector, test_vector)
        if (iProcIndex == root) then
            if (all(abs(test_vector-hamil_diag(HFindex)*davidson_eigenvector) < 1.0e-12)) then
                tSkipCalc = .true.
                davidson_eigenvalue = hamil_diag(HFindex)
            end if
        end if
        deallocate(test_vector)
        if (hamil_type == parallel_sparse_hamil_type) call MPIBCast(tSkipCalc)
        if (tSkipCalc) return

        ! Calculate the corresponding residual.
        call calculate_residual()

    end subroutine init_davidson

    subroutine subspace_expansion(basis_index)

        integer, intent(in) :: basis_index
        integer :: i
        real(dp) :: dot_prod, ddot, norm

        ! Create the new basis state from the residual. This step performs
        ! t = (D - EI)^(-1) r,
        ! where D is the diagonal of the Hamiltonian matrix, E is the eigenvalue previously
        ! calculated, I is the identity matrix and r is the residual.
        do i = 1, space_size
            basis_vectors(i, basis_index) = residual(i)/(hamil_diag(i) - davidson_eigenvalue)
        end do

        ! This step then maskes the new basis vector orthogonal to all other basis vectors, by doing
        ! t <- t - (t,v)v
        ! for each basis vector v, where (t,v) denotes the dot product.
        do i = 1, basis_index - 1
            dot_prod = ddot(space_size, basis_vectors(:,basis_index), 1, basis_vectors(:,i), 1)
            basis_vectors(:, basis_index) = basis_vectors(:, basis_index) - dot_prod*basis_vectors(:,i)
        end do

        ! Finally we calculate the norm of the new basis vector and then normalise it to have a norm of 1.
        ! The new basis vector is stored in the next available column in the basis_vectors array.
        norm = ddot(space_size, basis_vectors(:,basis_index), 1, basis_vectors(:,basis_index), 1)
        norm = sqrt(norm)
        basis_vectors(:,basis_index) = basis_vectors(:,basis_index)/norm

    end subroutine subspace_expansion

    subroutine project_hamiltonian(basis_index)

        integer, intent(in) :: basis_index
        integer :: i
        real(dp) :: multiplied_basis_vector(space_size), ddot

        if (iProcIndex == root) then
            ! Multiply the new basis_vector by the hamiltonian and store the result in
            ! multiplied_basis_vector.
            call multiply_hamil_and_vector(basis_vectors(:,basis_index), &
                multiplied_basis_vector)

            ! Now multiply U^T by (H U) to find projected_hamil. The projected Hamiltonian will
            ! only differ in the new final column and row. Also, projected_hamil is symmetric.
            ! Hence, we only need to calculate the final column, and use this to update the final
            ! row also.
            do i = 1, basis_index
                projected_hamil(i, basis_index) = ddot(space_size, basis_vectors(:, i), 1, multiplied_basis_vector, 1)
                projected_hamil(basis_index, i) = projected_hamil(i, basis_index)
            end do

            ! We will use the scrap Hamiltonian to pass into the diagonaliser later, since it
            ! overwrites this input matrix with the eigenvectors. Hence, make sure the scrap space
            ! stores the updated projected Hamiltonian.
            projected_hamil_scrap = projected_hamil
        else
            call multiply_hamil_and_vector(temp_in, temp_out)
        end if

    end subroutine project_hamiltonian

    subroutine subspace_extraction(basis_index)

        integer, intent(in) :: basis_index
        integer :: lwork, info
        real(dp), allocatable, dimension(:) :: work
        real(dp) :: eigenvalue_list(basis_index)
        real(dp) :: eigenvector_proj(basis_index)

        ! Scrap space for the diagonaliser.
        lwork = max(1,3*basis_index-1)
        allocate(work(lwork))

        ! This routine diagonalises a symmetric matrix, A.
        ! V tells the routine to calculate eigenvalues *and* eigenvectors.
        ! U tells the routine to get the upper half of A (it is symmetric).
        ! basis_index is the number of rows and columns in A.
        ! A = projected_hamil_scrap. This matrix stores the eigenvectors in its columns on output.
        ! basis_index is the leading dimension of A.
        ! eigenvalue_list stores the eigenvalues on output.
        ! work is scrap space.
        ! lwork is the length of the work array.
        ! info = 0 on output is diagonalisation is successful.
        call dsyev('V', 'U', basis_index, projected_hamil_scrap(1:basis_index,1:basis_index), basis_index, &
                       eigenvalue_list, work, lwork, info)

        davidson_eigenvalue = eigenvalue_list(1)
        ! The first column stores the ground state.
        eigenvector_proj = projected_hamil_scrap(1:basis_index,1)

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
        ! input x = eigenvector_proj.
        ! 1 is the increment of the elements of x.
        ! beta = 0.0_dp.
        ! output y = davidson_eigenvector.
        ! 1 is the incremenet of the elements of y.
        call dgemv('N', &
                   space_size, &
                   basis_index, &
                   1.0_dp, &
                   basis_vectors(:,1:basis_index), &
                   space_size, &
                   eigenvector_proj, &
                   1, &
                   0.0_dp, &
                   davidson_eigenvector, &
                   1)
        
    end subroutine subspace_extraction

    subroutine calculate_residual()

        ! This routine calculates the residual, r, corresponding to the new estimate of the
        ! ground state, stored in davidson_eigenvector. This is defined as
        ! r = Hv - Ev,
        ! where H is the Hamiltonian, v is the ground state vector estimate and E is the 
        ! ground state energy estimate.

        ! First multiply davidson_eigenvector by the Hamiltonian.
        if (iProcIndex == root) then 
            call multiply_hamil_and_vector(davidson_eigenvector, residual)
            ! Then simply take of the eigenvalue multiplied by the davidson_eigenvector...
            residual = residual - davidson_eigenvalue*davidson_eigenvector
        else
            call multiply_hamil_and_vector(temp_in, temp_out)
        end if

    end subroutine calculate_residual

    subroutine calculate_residual_norm()

        ! This subroutine calculates the Euclidean norm of the reisudal vector, r:
        ! residual_norm^2 = \sum_i r_i^2

        real(dp) :: ddot

        if (iProcIndex == root) then
            residual_norm = ddot(space_size, residual, 1, residual, 1)
            residual_norm = sqrt(residual_norm)
        end if

        if (hamil_type == parallel_sparse_hamil_type) call MPIBCast(residual_norm)

    end subroutine calculate_residual_norm

    subroutine multiply_hamil_and_vector(input_vector, output_vector)

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)

        if (hamil_type == full_hamil_type) then
            call multiply_hamil_and_vector_full(input_vector, output_vector)
        else if (hamil_type == sparse_hamil_type) then
            call mult_hamil_vector_sparse(input_vector, output_vector)
        else if (hamil_type == parallel_sparse_hamil_type) then
            call mult_hamil_vector_par_sparse(input_vector, output_vector)
        else if (hamil_type == direct_ci_type) then
            call mult_hamil_vector_direct_ci(input_vector, output_vector)
        end if

    end subroutine multiply_hamil_and_vector

    subroutine multiply_hamil_and_vector_full(input_vector, output_vector)

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)

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

    end subroutine multiply_hamil_and_vector_full

    subroutine mult_hamil_vector_sparse(input_vector, output_vector)

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)
        integer :: i, j

        output_vector = 0.0_dp

        do i = 1, space_size
            do j = 1, sparse_ham(i)%num_elements
                output_vector(i) = output_vector(i) + sparse_ham(i)%elements(j)*input_vector(sparse_ham(i)%positions(j))
            end do
        end do

    end subroutine mult_hamil_vector_sparse

    subroutine mult_hamil_vector_par_sparse(input_vector, output_vector)

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)
        integer :: i, j, ierr

        ! Use output_vector as temporary space.
        output_vector = input_vector

        call MPIBarrier(ierr)

        call MPIBCast(output_vector)

        partial_davidson_vector = 0.0_dp

        do i = 1, space_sizes(iProcIndex)
            do j = 1, sparse_ham(i)%num_elements
                partial_davidson_vector(i) = partial_davidson_vector(i) + &
                    sparse_ham(i)%elements(j)*output_vector(sparse_ham(i)%positions(j))
            end do
        end do

        call MPIGatherV(partial_davidson_vector, output_vector, space_sizes, davidson_disps, ierr)

    end subroutine mult_hamil_vector_par_sparse

    subroutine mult_hamil_vector_direct_ci(input_vector, output_vector)

        use direct_ci, only: perform_multiplication, transfer_from_block_form, transfer_to_block_form
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits
        use SystemData, only: ecore

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)

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

    subroutine end_davidson()

        integer :: ierr

        ! Deallocate all Davidson arrays. Note that the eigenvector is not deallocated,
        ! so that it can be used later.
        if (allocated(basis_vectors)) deallocate(basis_vectors)
        if (allocated(projected_hamil)) deallocate(projected_hamil)
        if (allocated(projected_hamil_scrap)) deallocate(projected_hamil_scrap)
        if (allocated(partial_davidson_vector)) deallocate(partial_davidson_vector)
        if (allocated(temp_in)) deallocate(temp_in)
        if (allocated(temp_out)) deallocate(temp_out)
        if (allocated(residual)) then
            deallocate(residual, stat=ierr)
            call LogMemDealloc("end_davidson", ResidualTag, ierr)
        end if

    end subroutine end_davidson

    subroutine davidson_direct_ci_init()

        use bit_rep_data, only: NIfD
        use direct_ci, only: create_direct_ci_arrays
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits
        use ras, only: initialise_ras_space, find_ras_size

        integer :: class_i, class_j, j, sym_i, sym_j

        call initialise_ras_space(davidson_ras, davidson_classes)
        ! The total hilbert space dimension of calculation to be performed.
        call find_ras_size(davidson_ras, davidson_classes, space_size)

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

    end subroutine davidson_direct_ci_init

    subroutine davidson_direct_ci_end()

        integer :: ierr

        if (allocated(hamil_diag)) then
            deallocate(hamil_diag, stat=ierr)
            call LogMemDealloc("davidson_direct_ci_end", HDiagTag, ierr)
        end if
        if (allocated(davidson_eigenvector)) then
            deallocate(davidson_eigenvector, stat=ierr)
            call LogMemDealloc("davidson_direct_ci_end", DavidsonTag, ierr)
        end if

        write(6,"(/,a10,f19.9)") "GROUND E =", davidson_eigenvalue

    end subroutine davidson_direct_ci_end

end module davidson_neci
