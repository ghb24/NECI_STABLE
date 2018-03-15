#include "macros.h" 

! a small module with functions needed for the unit-tests
module unit_test_helpers

    use constants, only: dp, EPS

    use lattice_mod, only: get_helement_lattice

    use Determinants, only: get_helement

    use SystemData, only: t_lattice_model, nOccAlpha, nOccBeta, &
                          trans_corr_param_2body, omega

    implicit none

contains

    function calc_eigenvalues(matrix) result(e_values)
        real(dp), intent(in) :: matrix(:,:) 
        real(dp) :: e_values(size(matrix,1))

        integer :: n, info, work(3*size(matrix,1))
        real(dp) :: tmp_matrix(size(matrix,1),size(matrix,2)),dummy_val(size(matrix,1))
        real(dp) :: dummy_vec_1(1,size(matrix,1)), dummy_vec_2(1,size(matrix,1))

        n = size(matrix,1)

        tmp_matrix = matrix 

        call dgeev('N','N', n, tmp_matrix, n, e_values, &
            dummy_val, dummy_vec_1,1,dummy_vec_2,1,work,3*n,info)

    end function calc_eigenvalues

    subroutine eig(matrix, e_values, e_vectors) 
        ! for very restricted matrices do a diag routine here! 
        real(dp), intent(in) :: matrix(:,:) 
        real(dp), intent(out) :: e_values(size(matrix,1))
        real(dp), intent(out), optional :: e_vectors(size(matrix,1),size(matrix,1))

        ! get the specifics for the eigenvectors still.. 
        ! i think i need a bigger work, and maybe also a flag for how many 
        ! eigenvectors i want.. maybe also the number of eigenvalues.. 
        integer :: n, info 
        real(dp) :: work(4*size(matrix,1)), tmp_matrix(size(matrix,1),size(matrix,2))
        real(dp) :: dummy_vec(1,size(matrix,1)), dummy(size(matrix,1))


        ! and convention is: we only want the right eigenvectors!!
        ! and always assume real-only eigenvalues
        if (present(e_vectors)) then 

            n = size(matrix,1)

            tmp_matrix = matrix 

            call dgeev('N','V',n,tmp_matrix, n, e_values, dummy, dummy_vec, & 
                1, e_vectors, n, work, 4*n, info)

        else
            e_values = calc_eigenvalues(matrix)
        end if
    end subroutine eig

    function create_hamiltonian(list_nI) result(hamil) 
        ! quite specific hamiltonian creation for my tests.. 
        integer, intent(in) :: list_nI(:,:) 
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j 

        do i = 1, size(list_nI,2)
            do j = 1, size(list_nI,2)
                hamil(i,j) = get_helement_lattice(list_nI(:,j),list_nI(:,i))
            end do
        end do

    end function create_hamiltonian

    function create_hamiltonian_old(list_nI) result(hamil) 
        ! try to also create the hamiltonian with the old implementation.. 
        ! although i think there needs to be more setup done.. 
        integer, intent(in) :: list_nI(:,:)
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j 

        t_lattice_model = .false.
        do i = 1, size(list_nI,2)
            do j = 1, size(list_nI,2)
                hamil(i,j) = get_helement(list_nI(:,j),list_nI(:,i))
            end do
        end do
        t_lattice_model = .true.

    end function create_hamiltonian_old

    subroutine print_matrix(matrix, iunit) 
        ! print a 2-D real matrix 
        class(*), intent(in) :: matrix(:,:)
        integer, intent(in), optional :: iunit
        
        integer :: i 

        select type (matrix)
        type is (integer)
            if (present(iunit)) then 
                do i = 1, size(matrix,1)
                    write(iunit,*) matrix(i,:)
                end do
            else
                do i = 1, size(matrix,1)
                    print *, matrix(i,:)
                end do
            end if
        type is (real(dp))
            if (present(iunit)) then 
                do i = 1, size(matrix,1)
                    write(iunit,*) matrix(i,:)
                end do
            else
                do i = 1, size(matrix,1)
                    print *, matrix(i,:)
                end do
            end if
        end select


    end subroutine print_matrix

    function similarity_transform(H, t_mat_opt) result(trans_H)
        HElement_t(dp), intent(in) :: H(:,:)
        real(dp), intent(in), optional :: t_mat_opt(:,:)
        real(dp) :: trans_H(size(H,1),size(H,2))

        real(dp) :: t_mat(size(H,1),size(H,2))

        if (present(t_mat_opt)) then 
            t_mat = t_mat_opt
        else
            ! otherwise assume the on-site correlation factor is used
            t_mat = get_tranformation_matrix(H, nOccAlpha*nOccBeta) 
        end if

        trans_H = blas_matmul(blas_matmul(matrix_exponential(-t_mat), H), matrix_exponential(t_mat))

    end function similarity_transform

    function get_tranformation_matrix(hamil, n_pairs) result(t_matrix)
        ! n_pairs is actually also a global system dependent quantitiy.. 
        ! which actually might be helpful.. but input it here! 
        HElement_t(dp), intent(in) :: hamil(:,:)
        integer, intent(in) :: n_pairs
        real(dp) :: t_matrix(size(hamil,1),size(hamil,2))

        integer :: i, j

        t_matrix = 0.0_dp 

        do i = 1, size(hamil,1)
            do j = 1, size(hamil,1)
                if (i == j) then 
                    t_matrix(i,i) = n_pairs 
                else 
                    if (abs(hamil(i,j)) > EPS) then 
                        t_matrix(i,j) = sign(1.0_dp, hamil(i,j))
                    end if
                end if
            end do
        end do

        t_matrix = trans_corr_param_2body/omega * t_matrix

    end function get_tranformation_matrix

    function blas_matmul(A, B) result(C)
        ! a basic wrapper to the most fundamental matrix mult with blas 
        real(dp), intent(in) :: A(:,:), B(:,:) 
        real(dp) :: C(size(A,1),size(A,2))

        integer :: n 

        n = size(A,1) 

        call dgemm('N', 'N', n, n, n, 1.0_dp, A, n, B, n, 0.0_dp, C, n)

    end function blas_matmul

    function linspace(start_val, end_val, n_opt) result(vec) 
        real(dp), intent(in) :: start_val, end_val 
        integer, intent(in), optional :: n_opt
        real(dp), allocatable :: vec(:) 

        integer :: n, i 
        real(dp) :: dist 

        ! set default 
        if (present(n_opt)) then 
            n = n_opt
        else 
            n = 100
        end if

        dist = (end_val - start_val) / real(n - 1, dp) 

        allocate(vec(n)) 

        vec = [ ( start_val + i * dist, i = 0,n-1)]

    end function linspace 

    function matrix_exponential(matrix) result(exp_matrix)
        ! calculate the matrix exponential of a real, symmetric 2-D matrix with lapack 
        ! routines
        ! i need A = UDU^-1
        ! e^A = Ue^DU^-1
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: exp_matrix(size(matrix,1),size(matrix,2))

        ! maybe i need to allocate this stuff:
        real(dp) :: vectors(size(matrix,1),size(matrix,2)), values(size(matrix,1))
        real(dp) :: work(3*size(matrix,1)-1), inverse(size(matrix,1),size(matrix,2))
        real(dp) :: exp_diag(size(matrix,1),size(matrix,2))
        integer :: info, n

        n = size(matrix,1)

        ! first i need to diagonalise the matrix and calculate the 
        ! eigenvectors 
        vectors = matrix

        call dsyev('V', 'U', n, vectors, n, values, work, 3*n-1,info)

        ! now i have the eigenvectors, which i need the inverse of 
        ! it is rotation only or? so i would just need a transpose or?
!         inverse = matrix_inverse(vectors) 
        inverse = transpose(vectors)

        ! i need to construct exp(eigenvalues) as a diagonal matrix! 
        exp_diag = matrix_diag(exp(values))

        exp_matrix = blas_matmul(blas_matmul(vectors,exp_diag),inverse)
!         exp_matrix = matmul(matmul(vectors,exp_diag),inverse)
        
    end function matrix_exponential

    function matrix_diag(vector) result(diag) 
        ! constructs a diagonal matrix with the vector on the diagonal 
        real(dp), intent(in) :: vector(:)
        real(dp) :: diag(size(vector),size(vector))

        integer :: i 

        diag = 0.0_dp

        do i = 1, size(vector)
            diag(i,i) = vector(i)
        end do

    end function matrix_diag

    function matrix_inverse(matrix) result(inverse)
        ! from fortran-wiki! search there for "matrix+inversion"
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: inverse(size(matrix,1),size(matrix,2))
        character(*), parameter :: this_routine = "matrix_inverse"

        real(dp) :: work(size(matrix,1))
        integer :: ipiv(size(matrix,1))
        integer :: n, info

        inverse = matrix 
        n = size(matrix,1)

        call dgetrf(n,n,inverse,n,ipiv,info)

        if (info /= 0) call stop_all(this_routine, "matrix singular!")

        call dgetri(n, inverse, n, ipiv, work, n, info)

        if (info /= 0) call stop_all(this_routine, "matrix inversion failed!")

    end function matrix_inverse


end module unit_test_helpers


