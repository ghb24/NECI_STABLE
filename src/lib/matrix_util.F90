#include "macros.h"

module matrix_util
    use constants, only: sp, dp, EPS
    use util_mod, only: near_zero, get_free_unit
    use sort_mod, only: sort
    implicit none
    private
    public :: eig, print_matrix, matrix_exponential, det, blas_matmul, linspace, norm, &
        calc_eigenvalues, check_symmetric, find_degeneracies, eig_sym, norm_cmplx, &
        store_hf_coeff, my_minloc, my_minval, matrix_inverse, print_vec

    interface linspace
        module procedure linspace_sp
        module procedure linspace_dp
    end interface linspace

    interface eig
        module procedure eig_cmplx
        module procedure eig_real
    end interface eig

    interface calc_eigenvalues
        module procedure calc_eigenvalues_real
        module procedure calc_eigenvalues_cmplx
    end interface calc_eigenvalues

contains

    subroutine print_vec(vec, filename, t_index, t_zero)
        class(*), intent(in) :: vec(:)
        character(*), intent(in), optional :: filename
        logical, intent(in), optional :: t_index, t_zero

        logical :: t_index_, t_zero_
        integer :: iunit, i
        def_default(t_index_, t_index, .false.)
        def_default(t_zero_, t_zero, .false.)

        select type(vec)
        type is (integer)
            if (present(filename)) then
                iunit = get_free_unit()
                open(iunit, file = filename, status = 'replace', action = 'write')

                if (t_zero_) then
                    if (t_index_) then
                        write(iunit, *) 0, 0.0_dp
                    else
                        write(iunit, *) 0.0_dp
                    end if
                end if


                if (t_index_) then
                    do i = 1, size(vec,1)
                        write(iunit, *) i, vec(i)
                    end do
                else
                    do i = 1, size(vec,1)
                        write(iunit, *) vec(i)
                    end do
                end if

                close(iunit)
            else
                do i = 1, size(vec,1)
                    print *, vec(i)
                end do
            end if
        type is (real(dp))
            if (present(filename)) then
                iunit = get_free_unit()
                open(iunit, file = filename, status = 'replace', action = 'write')

                if (t_zero_) then
                    if (t_index_) then
                        write(iunit, *) 0, 0.0_dp
                    else
                        write(iunit, *) 0.0_dp
                    end if
                end if


                if (t_index_) then
                    do i = 1, size(vec,1)
                        write(iunit, *) i, vec(i)
                    end do
                else
                    do i = 1, size(vec,1)
                        write(iunit, *) vec(i)
                    end do
                end if

                close(iunit)
            else
                do i = 1, size(vec,1)
                    print *, vec(i)
                end do
            end if

        end select

    end subroutine print_vec

    subroutine eig_real(matrix, e_values, e_vectors, t_left_ev)
        ! for very restricted matrices do a diag routine here!
        real(dp), intent(in) :: matrix(:,:)
        real(dp), intent(out) :: e_values(size(matrix,1))
        real(dp), intent(out), optional :: e_vectors(size(matrix,1),size(matrix,1))
        logical, intent(in), optional :: t_left_ev

        ! get the specifics for the eigenvectors still..
        ! i think i need a bigger work, and maybe also a flag for how many
        ! eigenvectors i want.. maybe also the number of eigenvalues..
        integer :: n, info
        real(dp) :: work(4*size(matrix,1)), tmp_matrix(size(matrix,1),size(matrix,2))
        real(dp) :: left_ev(size(matrix,1),size(matrix,1)), dummy_eval(size(matrix,1))
        real(dp) :: right_ev(size(matrix,1),size(matrix,1))
        integer :: sort_ind(size(matrix,1))
        character :: left, right


        ! and convention is: we only want the right eigenvectors!!
        ! and always assume real-only eigenvalues
        if (present(e_vectors)) then

            if (present(t_left_ev)) then
                if (t_left_ev) then
                    left = 'V'
                    right = 'N'
                else
                    left = 'N'
                    right = 'V'
                end if
            else
                left = 'N'
                right = 'V'
            end if

            n = size(matrix,1)

            tmp_matrix = matrix

            left = 'V'
            right = 'V'

            call dgeev(&
                left, &
                right, &
                n, &
                tmp_matrix, &
                n, &
                e_values, &
                dummy_eval, &
                left_ev, &
                n, &
                right_ev, &
                n, &
                work, &
                4*n, &
                info)

            sort_ind = [(n, n = 1, size(matrix,1))]

            call sort(e_values, sort_ind)

            if (present(t_left_ev)) then
                if (t_left_ev) then
                    e_vectors = left_ev(:,sort_ind)
                else
                    e_vectors = right_ev(:,sort_ind)
                end if
            else
                e_vectors = right_ev(:,sort_ind)
            end if

        else
            e_values = calc_eigenvalues(matrix)
        end if

    end subroutine eig_real


    subroutine eig_cmplx(matrix, e_values, e_vectors, t_left_ev)
        ! for very restricted matrices do a diag routine here!
        complex(dp), intent(in) :: matrix(:,:)
        real(dp), intent(out) :: e_values(size(matrix,1))
        complex(dp), intent(out), optional :: e_vectors(size(matrix,1),size(matrix,1))
        logical, intent(in), optional :: t_left_ev
        character(*), parameter :: this_routine = 'eig_cmplx'

        ! get the specifics for the eigenvectors still..
        ! i think i need a bigger work, and maybe also a flag for how many
        ! eigenvectors i want.. maybe also the number of eigenvalues..
        integer :: n, info
        complex(dp) :: work(4*size(matrix,1)), tmp_matrix(size(matrix,1),size(matrix,2))
        complex(dp) :: left_ev(size(matrix,1),size(matrix,1))
        real(dp), allocatable :: rwork(:)
        real(dp) :: right_ev(size(matrix,1),size(matrix,1))
        integer :: sort_ind(size(matrix,1))
        character(len=1) :: left, right


        ! and convention is: we only want the right eigenvectors!!
        ! and always assume real-only eigenvalues
        if (present(e_vectors)) then

            if (present(t_left_ev)) then
                if (t_left_ev) then
                    left = 'V'
                    right = 'N'
                else
                    left = 'N'
                    right = 'V'
                end if
            else
                left = 'N'
                right = 'V'
            end if

            n = size(matrix,1)

            tmp_matrix = matrix

            left = 'V'
            right = 'V'


            allocate(rwork(max(1,3*n-2)))
            call zheev(&
                 left, &
                 right, &
                 n, &
                 tmp_matrix, &
                 n, &
                 e_values, &
                 work, &
                 4*n, &
                 rwork, &
                 info)
            if (info /= 0) call stop_all(this_routine, 'Failed in BLAS call.')
            deallocate(rwork)

            sort_ind = [(n, n = 1, size(matrix,1))]

            call sort(e_values, sort_ind)

            if (present(t_left_ev)) then
                if (t_left_ev) then
                    e_vectors = left_ev(:,sort_ind)
                else
                    e_vectors = right_ev(:,sort_ind)
                end if
            else
                e_vectors = right_ev(:,sort_ind)
            end if

        else
            e_values = calc_eigenvalues(matrix)
        end if

    end subroutine eig_cmplx

    function calc_eigenvalues_real(matrix) result(e_values)
        real(dp), intent(in) :: matrix(:,:)
        real(dp) :: e_values(size(matrix,1))
        character(*), parameter :: this_routine = 'calc_eigenvalues_real'

        integer :: n, info
        real(dp) :: work(3*size(matrix,1))
        real(dp) :: tmp_matrix(size(matrix,1),size(matrix,2)),dummy_val(size(matrix,1))
        real(dp) :: dummy_vec_1(1,size(matrix,1)), dummy_vec_2(1,size(matrix,1))

        n = size(matrix,1)

        tmp_matrix = matrix
        call dgeev('N','N', n, tmp_matrix, n, e_values, &
            dummy_val, dummy_vec_1,1,dummy_vec_2,1,work,3*n,info)
        if (info /= 0) call stop_all(this_routine, 'Failed in BLAS call.')
        call sort(e_values)

    end function calc_eigenvalues_real

    function calc_eigenvalues_cmplx(matrix) result(e_values)
        complex(dp), intent(in) :: matrix(:,:)
        real(dp) :: e_values(size(matrix,1))
        character(*), parameter :: this_routine = 'calc_eigenvalues_cmplx'

        integer :: n, info
        complex(dp) :: work(3*size(matrix,1))
        complex(dp) :: tmp_matrix(size(matrix,1),size(matrix,2))
        real(dp), allocatable :: rwork(:)

        n = size(matrix,1)

        tmp_matrix = matrix
        allocate(rwork(max(1, 3*n - 2)))
        call zheev('N','N', n, tmp_matrix, n, e_values, work, 3*n, rwork, info)
        if (info /= 0) call stop_all(this_routine, 'Failed in BLAS call.')
        deallocate(rwork)
        call sort(e_values)

    end function calc_eigenvalues_cmplx

    subroutine eig_sym(matrix, e_values, e_vectors)
        real(dp), intent(in) :: matrix(:,:)
        real(dp), intent(out) :: e_values(size(matrix,1))
        real(dp), intent(out), optional :: e_vectors(size(matrix,1),size(matrix,2))
        character(*), parameter :: this_routine = 'eig_sym'
        integer :: n, info, lwork
        character(1) :: jobz
        real(dp) :: tmp_matrix(size(matrix,1),size(matrix,2))
        real(dp) :: work(3*size(matrix,1))

        n = size(matrix,1)
        lwork = 3*n

        tmp_matrix = matrix

        if (present(e_vectors)) then
            jobz = 'V'
        else
            jobz = 'N'
        end if

        call dsyev(&
            jobz, &
            'U', &
            n, &
            tmp_matrix, &
            n, &
            e_values, &
            work, &
            lwork, &
            info)
        if (info /= 0) call stop_all(this_routine, 'Failed in BLAS call.')

        if (present(e_vectors)) then
            e_vectors = tmp_matrix
        end if

    end subroutine eig_sym

    logical function check_symmetric(matrix)
        ! function to check if a given matrix is symmetric
        ! for a square matrix!
        real(dp), intent(in) :: matrix(:,:)

        real(dp) :: diff(size(matrix,1),size(matrix,2))

        diff = matrix - transpose(matrix)

        check_symmetric = .false.

        if (near_zero(sum(abs(diff)))) check_symmetric = .true.

    end function check_symmetric

    subroutine print_matrix(matrix, iunit)
        ! print a 2-D real matrix
        class(*), intent(in) :: matrix(:,:)
        integer, intent(in), optional :: iunit

        integer :: i, j, tmp_unit

        select type (matrix)
        type is (integer)
            if (present(iunit)) then
                do i = lbound(matrix,1), ubound(matrix,1)
                    write(iunit,*) matrix(i,:)
                end do
            else
                do i = lbound(matrix,1), ubound(matrix,1)
                    print *, matrix(i,:)
                end do
            end if
        type is (real(dp))
            if (present(iunit)) then
                do i = lbound(matrix,1),ubound(matrix,1)
                    do j = lbound(matrix,2), ubound(matrix,2) - 1
                        write(iunit,'(G25.17)', advance = 'no') matrix(i,j)
                    end do
                    write(iunit,'(G25.17)', advance = 'yes') matrix(i,j)
                end do
            else
                do i = lbound(matrix,1),ubound(matrix,1)
                    print *, matrix(i,:)
                end do
            end if
        type is (complex(dp))
            if (present(iunit)) then
                tmp_unit = iunit
            else
                tmp_unit = 6
            end if
            do i = lbound(matrix,1),ubound(matrix,1)
                do j = lbound(matrix,2), ubound(matrix,2)
                    if (j < ubound(matrix,2)) then
                        write(tmp_unit,fmt = '(F10.8,SP,F10.8,"i",1x)', advance = 'no') matrix(i,j)
                    else
                        write(tmp_unit,fmt = '(F10.8,SP,F10.8,"i")', advance = 'yes') matrix(i,j)
                    end if
                end do
            end do
        end select


    end subroutine print_matrix


    real(dp) function det(matrix)
        real(dp), intent(in) :: matrix(:,:)

        integer :: n, i, info
        integer, allocatable :: ipiv(:)
        real(dp) :: sgn
        real(dp), allocatable :: tmp_matrix(:,:)

        n = size(matrix,1)
        allocate(tmp_matrix(n,n), source = matrix)
        allocate(ipiv(n))

        ipiv = 0

        call dgetrf(n,n,tmp_matrix,n,ipiv,info)

        det = 1.0_dp

        do i = 1, N
            det = det * tmp_matrix(i,i)
        end do

        sgn = 1.0_dp
        do i = 1, n
            if (ipiv(i) /= i) then
                sgn = -sgn
            end if
        end do

        det = sgn * det

    end function det


    function blas_matmul(A, B) result(C)
        ! a basic wrapper to the most fundamental matrix mult with blas
        HElement_t(dp), intent(in) :: A(:,:), B(:,:)
        HElement_t(dp) :: C(size(A,1),size(A,2))

        integer :: n

        n = size(A,1)
#ifdef CMPLX_
        call zgemm('N','N', n, n, n, cmplx(1.0_dp,0.0_dp,kind=dp), A, n, B, n, &
            cmplx(1.0_dp,0.0_dp,kind=dp), C, n)
#else
        call dgemm('N', 'N', n, n, n, 1.0_dp, A, n, B, n, 0.0_dp, C, n)
#endif
    end function blas_matmul

    function linspace_sp(start_val, end_val, n_opt) result(vec)
        real(sp), intent(in) :: start_val, end_val
        integer, intent(in), optional :: n_opt
        real(sp), allocatable :: vec(:)

        integer :: n, i
        real(sp) :: dist

        ! set default
        if (present(n_opt)) then
            n = n_opt
        else
            n = 100
        end if

        dist = (end_val - start_val) / real(n - 1, sp)

        allocate(vec(n))

        vec = [ ( start_val + i * dist, i = 0,n-1)]

    end function linspace_sp


    function linspace_dp(start_val, end_val, n_opt) result(vec)
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

    end function linspace_dp

    function matrix_exponential(matrix) result(exp_matrix)
        ! calculate the matrix exponential of a real, symmetric 2-D matrix with lapack
        ! routines
        ! i need A = UDU^-1
        ! e^A = Ue^DU^-1
        HElement_t(dp), intent(in) :: matrix(:,:)
        HElement_t(dp) :: exp_matrix(size(matrix,1),size(matrix,2))

        ! maybe i need to allocate this stuff:
        HElement_t(dp) :: vectors(size(matrix,1),size(matrix,2))
        real(dp) :: values(size(matrix,1))
        HElement_t(dp) :: work(3*size(matrix,1)-1)
        HElement_t(dp) :: inverse(size(matrix,1),size(matrix,2))
        HElement_t(dp) :: exp_diag(size(matrix,1),size(matrix,2))
        integer :: info, n

        n = size(matrix,1)

        ! first i need to diagonalise the matrix and calculate the
        ! eigenvectors
        vectors = matrix
#ifdef CMPLX_
        block
            real(dp), allocatable :: rwork(:)
            allocate(rwork(max(1, 3*n - 2)))
            call zheev('V', 'U', n, vectors, n, values, work, 3*n-1, rwork, info)
            deallocate(rwork)
        end block
#else
        call dsyev('V', 'U', n, vectors, n, values, work, 3*n-1,info)
#endif
        ! now i have the eigenvectors, which i need the inverse of
        ! it is rotation only or? so i would just need a transpose or?
        inverse = transpose(vectors)

        ! i need to construct exp(eigenvalues) as a diagonal matrix!
        exp_diag = matrix_diag(exp(values))

        exp_matrix = blas_matmul(blas_matmul(vectors,exp_diag), inverse)

    end function matrix_exponential

    function matrix_diag(vector) result(diag)
        ! constructs a diagonal matrix with the vector on the diagonal
        real(dp), intent(in) :: vector(:)
        HElement_t(dp) :: diag(size(vector),size(vector))

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

    real(dp) function norm(vec,p_in)
        ! function to calculate the Lp norm of a given vector
        ! if p_in = -1 this indicates the p_inf norm
        real(dp), intent(in) :: vec(:)
        integer, intent(in), optional :: p_in
        integer :: p, i

        if (present(p_in)) then
!             ASSERT(p_in == -1 .or. p_in >= 0)
            p = p_in
        else
            ! default is the L2 norm
            p = 2
        end if

        if (p == -1) then
            norm = maxval(abs(vec))
        else
            norm = 0.0_dp
            do i = 1, size(vec)
                norm = norm + abs(vec(i))**p
            end do

            norm = norm**(1.0_dp/real(p,dp))
        end if

    end function norm

    subroutine store_hf_coeff(e_values, e_vecs, target_state, hf_coeff, hf_ind, gs_ind)
        real(dp), intent(in) :: e_values(:), e_vecs(:,:)
        integer, intent(in), optional :: target_state
        real(dp), intent(out) :: hf_coeff
        integer, intent(out) :: hf_ind, gs_ind

        real(dp) :: gs_vec(size(e_values))
        integer :: target_state_
        def_default(target_state_,target_state,1)

        gs_ind = my_minloc(e_values, target_state)

        gs_vec = abs(e_vecs(:,gs_ind))

        hf_ind = maxloc(gs_vec,1)
        hf_coeff = gs_vec(hf_ind)

    end subroutine store_hf_coeff

    pure real(dp) function my_minval(vec, target_state)
        real(dp), intent(in) :: vec(:)
        integer, intent(in), optional :: target_state

        if (present(target_state)) then
            my_minval = vec(my_minloc(vec,target_state))
        else
            my_minval = minval(vec)
        end if

    end function my_minval


    pure integer function my_minloc(vec, target_state)
        real(dp), intent(in) :: vec(:)
        integer, intent(in), optional :: target_state

        logical :: flag(size(vec))
        integer :: i


        if (present(target_state)) then
            flag = .true.
            do i = 1, target_state
                my_minloc = minloc(vec, dim = 1, mask = flag)
                flag(my_minloc) = .false.
            end do
        else
            my_minloc = minloc(vec, 1)
        end if

    end function my_minloc

    real(dp) function norm_cmplx(vec,p_in)
        ! function to calculate the Lp norm of a given vector
        ! if p_in = -1 this indicates the p_inf norm
        complex(dp), intent(in) :: vec(:)
        integer, intent(in), optional :: p_in
        integer :: p, i

        if (present(p_in)) then
            p = p_in
        else
            ! default is the L2 norm
            p = 2
        end if

        if (p == -1) then
            norm_cmplx = maxval(abs(vec))
        else
            norm_cmplx = 0.0_dp
            do i = 1, size(vec)
                norm_cmplx = norm_cmplx + abs(vec(i))**p
            end do

            norm_cmplx = norm_cmplx**(1.0_dp/real(p,dp))
        end if

    end function norm_cmplx


    subroutine find_degeneracies(e_values, ind, pairs)
        ! find the indices of degenerate eigenvalues
        ! ind will have as many rows as degenerate eigenvalues exist
        ! and the columns will be the maximum number of degeneracy + 1
        ! since in the first column the number of degenerate eigenvalues are
        ! stored!
        ! it assumes the eigenvalues are sorted!!
        ! in pairs the paired indices of the degenerate eigenvalue are stored!
        real(dp), intent(in) :: e_values(:)
        integer, intent(out), allocatable :: ind(:,:), pairs(:,:)

        integer :: i,j,tmp_ind(size(e_values),size(e_values)+1), e_ind
        integer :: max_val

        tmp_ind = 0
        e_ind = 1
        i = 1

        do while (i < size(e_values) .and. e_ind < size(e_values))
            j = 0
            do while(e_ind + j <= size(e_values))
                if (abs(e_values(e_ind) - e_values(e_ind+j)) < 10.e-8) then
                    tmp_ind(i,j+2) = e_ind+j
                    j = j + 1
                else
                    exit
                end if
            end do
            tmp_ind(i,1) = j
            i = i + 1
            e_ind = e_ind + j
        end do

        ! deal with end-value specifically
        if (e_ind == size(e_values)) then
            tmp_ind(i,1) = 1
            tmp_ind(i,2) = e_ind
        end if

        max_val = maxval(tmp_ind(:,1))+1
        allocate(ind(i-1,max_val), source = tmp_ind(1:i-1,1:max_val))

        if (max_val == 2) then
            ! if no degeneracies
            allocate(pairs(size(e_values),1))
            pairs = 0
            return
        end if
        allocate(pairs(size(e_values),max_val-2))
        pairs = 0
        ! do it in a stupid way and reuse the created array ind
        do i = 1, size(ind,1)
            if (ind(i,1) > 1) then
                do j = 2, ind(i,1) + 1
                    pairs(ind(i,j),:) = pack(ind(i,2:ind(i,1)+1), &
                        ind(i,2:ind(i,1)+1) /= ind(i,j))
                end do
            end if
        end do

    end subroutine find_degeneracies


end module
