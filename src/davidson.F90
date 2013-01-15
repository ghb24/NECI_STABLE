module davidson

! This module performs the Davidson method to find the ground state of a diagonally-
! dominant matrix. For details of the theory behind the method, see i.e:
! http://web.mit.edu/bolin/www/Project-Report-18.335J.pdf

use constants
use FciMCData, only: hamiltonian, sparse_matrix_info, sparse_hamil, hamil_diag

implicit none

integer :: max_num_davidson_iters = 50
real(dp) :: residual_norm_target = 0.0000001

! If this is true then the Hamiltonian is must be stored in a sparse form, using
! the sparse_matrix_info type (see the FciMCData module). This is then assumed
! in all the Hamiltonian matrix calls.
logical :: sparse_multiply = .false.

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

    contains

    subroutine perform_davidson(sparse)

        ! This subroutine performs the main loop for the algorithm, from which each of the
        ! induvidual steps are called.

        logical, intent(in) :: sparse
        integer :: i

        if (sparse) sparse_multiply = .true.

        call init_davidson()

        do i = 2, max_num_davidson_iters

            call subspace_expansion(i)

            call project_hamiltonian(i)

            call subspace_extraction(i)

            call calculate_residual()

            call calculate_residual_norm()

            print *, "iteration:", i, "residual norm:", residual_norm
            ! If the solution is sufficiently converged then exit the loop.
            if (residual_norm < residual_norm_target) exit
            
        end do

        call end_davidson()

    end subroutine perform_davidson

    subroutine init_davidson()
    
        use Determinants, only: get_helement
        use FciMCData, only: HFDet

        ! This subroutine initialises the Davdison method by allocating the necessary arrays,
        ! defining the initial basis vector and projected Hamiltonian, and setting an initial
        ! guess at the ground state eigenvalue. It also calculates the corresponding residual
        ! which is needed to expand the space.

        integer :: i, HFindex

        ! If not using the sparse storage of the Hamiltonian, allocate the Hamiltonian diagonal.
        if (.not. allocated(hamil_diag)) then
            allocate(hamil_diag(size(hamiltonian,1)))
            do i = 1, size(hamiltonian,1)
                hamil_diag(i) = hamiltonian(i,i)
            end do
        end if

        HFindex = maxloc(abs(hamil_diag),1)

        ! The space size is determined by the size of the Hamiltonian matrix:
        space_size = size(hamil_diag)

        ! Allocate the necessary arrays:
        allocate(basis_vectors(space_size, max_num_davidson_iters))
        allocate(projected_hamil(1,1))
        allocate(projected_hamil_scrap(1,1))
        ! If a davidson calculation has already been performed, this array might still be
        ! allocated, so check!
        if (allocated(davidson_eigenvector)) deallocate(davidson_eigenvector)
        allocate(davidson_eigenvector(space_size))
        allocate(residual(space_size))
        projected_hamil = 0.0_dp
        basis_vectors = 0.0_dp
        davidson_eigenvector = 0.0_dp
        residual = 0.0_dp

        ! For the initial basis vector, choose the Hartree-Fock state:
        basis_vectors(HFindex, 1) = 1.0_dp
        ! Choose the Hartree-Fock state as the initial guess at the ground state, too.
        davidson_eigenvector(HFindex) = 1.0_dp

        ! Fill in the projected Hamiltonian so far.
        projected_hamil(1,1) = hamil_diag(HFindex)
        ! Take the initial eigenvalue to be the Hartree-Fock energy minus some small
        ! amount. This value cannot be exactly the Hartree-Fock energy, as this will
        ! result in dividing by zero in the subspace expansion step.
        davidson_eigenvalue = hamil_diag(1) - 0.001

        ! Calculate the corresponding residual:
        call calculate_residual()

    end subroutine init_davidson

    subroutine subspace_expansion(basis_index)

        integer, intent(in) :: basis_index
        integer :: i
        real(dp) :: dot_prod, norm

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
            dot_prod = dot_product(basis_vectors(:,basis_index), basis_vectors(:,i))
            basis_vectors(:, basis_index) = basis_vectors(:, basis_index) - dot_prod*basis_vectors(:,i)
        end do

        ! Finally we calculate the norm of the new basis vector and then normalise it to have a norm of 1.
        norm = 0.0_dp
        do i = 1, space_size
            norm = norm + basis_vectors(i, basis_index)*basis_vectors(i, basis_index)
        end do
        norm = sqrt(norm)

        ! The new basis vector is stored in the next avaliable column in the basis_vectors array.
        basis_vectors(:, basis_index) = basis_vectors(:, basis_index)/norm

    end subroutine subspace_expansion

    subroutine project_hamiltonian(basis_index)

        integer, intent(in) :: basis_index
        integer :: i
        real(dp) :: multiplied_basis_vector(space_size)

        ! Slightly hacky: Store projected_hamil in the scrap space and then reallocate
        ! projected_hamil to be of a larger size. Then reallocate the scrap space too.
        projected_hamil_scrap = projected_hamil
        deallocate(projected_hamil)
        allocate(projected_hamil(basis_index, basis_index))
        projected_hamil = 0.0_dp
        projected_hamil(1:(basis_index-1), 1:(basis_index-1)) = &
            projected_hamil_scrap(1:(basis_index-1), 1:(basis_index-1))
        deallocate(projected_hamil_scrap)
        allocate(projected_hamil_scrap(basis_index, basis_index))

        ! Multiply the new basis_vector by the hamiltonian and store the result in
        ! multiplied_basis_vector.
        call multiply_hamil_and_vector(basis_vectors(:,basis_index), &
            multiplied_basis_vector)

        ! Now multiply U^T by (H U) to find projected_hamil. The projected Hamiltonian will
        ! only differ in the new final column and row. Also, projected_hamil is symmetric.
        ! Hence, we only need to calculate the final column, and use this to update the final
        ! row also.
        do i = 1, basis_index
            projected_hamil(i, basis_index) = projected_hamil(i, basis_index) + &
                dot_product(basis_vectors(:, i), multiplied_basis_vector)
            projected_hamil(basis_index, i) = projected_hamil(i, basis_index)
        end do

        ! We will use the scrap Hamiltonian to pass into the diagonaliser later, since it
        ! overwrites this input matrix with the eigenvectors. Hence, make sure the scrap space
        ! stores the updated projected Hamiltonian.
        projected_hamil_scrap = projected_hamil

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
        ! V specifies tells the routine to calculate eigenvalues *and* eigenvectors.
        ! U tells the routine to get the upper half of A (it is symmetric).
        ! basis_index is the number of rows and columns in A.
        ! A = projected_hamil_scrap. This matrix stores the eigenvectors in its columns on output.
        ! basis_index is the leading dimension of A.
        ! eigenvalue_list stores the eigenvalues on output.
        ! work is scrap space.
        ! lwork is the length of the work array.
        ! info = 0 on output is diagonalisation is successful.
        call dsyev('V', 'U', basis_index, projected_hamil_scrap, basis_index, &
                       eigenvalue_list, work, lwork, info)

        davidson_eigenvalue = eigenvalue_list(1)
        print *, davidson_eigenvalue
        ! The first column stores the ground state.
        eigenvector_proj = projected_hamil_scrap(:,1)

        deallocate(work)

        ! eigenvector_proj stores the eigenstate in the basis of vectors stored in the array
        ! basis_vectors. We now want it in terms of the original basis. To get this, multiply
        ! eigenvector_proj by basis_vectors(:, 1:Basis_index):

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
        call multiply_hamil_and_vector(davidson_eigenvector, residual)

        ! Then simply take of the eigenvalue multiplied by the davidson_eigenvector...
        residual = residual - davidson_eigenvalue*davidson_eigenvector

    end subroutine calculate_residual

    subroutine calculate_residual_norm()

        ! This subroutine calculates the Euclidean norm of the reisudal vector, r:
        ! residual_norm^2 = \sum_i r_i^2

        integer :: i

        residual_norm = 0.0_dp

        do i = 1, space_size
            residual_norm = residual_norm + residual(i)*residual(i)
        end do

        residual_norm = sqrt(residual_norm)

    end subroutine calculate_residual_norm

    subroutine multiply_hamil_and_vector(input_vector, output_vector)

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)

        if (sparse_multiply) then
            call multiply_hamil_and_vector_sparse(input_vector, output_vector)
        else if (.not. sparse_multiply) then
            call multiply_hamil_and_vector_naive(input_vector, output_vector)
        end if

    end subroutine multiply_hamil_and_vector

    subroutine multiply_hamil_and_vector_naive(input_vector, output_vector)

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

    end subroutine multiply_hamil_and_vector_naive

    subroutine multiply_hamil_and_vector_sparse(input_vector, output_vector)

        real(dp), intent(in) :: input_vector(space_size)
        real(dp), intent(out) :: output_vector(space_size)
        integer :: i, j

        output_vector = 0.0_dp

        do i = 1, space_size
            do j = 1, sparse_hamil(i)%num_elements
                output_vector(i) = output_vector(i) + &
                    sparse_hamil(i)%elements(j) * &
                    input_vector(sparse_hamil(i)%positions(j))
            end do
        end do

    end subroutine multiply_hamil_and_vector_sparse

    subroutine end_davidson()

        ! Deallocate all Davidson arrays. Note that the eigenvector is not deallocated
        ! so that it can be used later.
        deallocate(basis_vectors)
        deallocate(projected_hamil)
        deallocate(projected_hamil_scrap)
        deallocate(residual)

    end subroutine end_davidson

end module davidson
