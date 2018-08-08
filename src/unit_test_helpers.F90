#include "macros.h" 

! a small module with functions needed for the unit-tests
module unit_test_helpers

    use constants, only: dp, EPS, n_int

    use lattice_mod, only: get_helement_lattice, lattice

    use Determinants, only: get_helement

    use SystemData, only: t_lattice_model, nOccAlpha, nOccBeta, &
                          trans_corr_param_2body, omega, nel, nBasis, &
                          arr, brr, nBasis, bhub, tGUGA

    use fcimcdata, only: excit_gen_store_type

    use util_mod, only: binary_search, choose

    use bit_reps, only: decode_bit_det
    
    use bit_rep_data, only: niftot, nifd

    use DetBitOps, only: ilut_lt, ilut_gt, EncodeBitDet, findbitexcitlevel, &
                         count_open_orbs

    use sort_mod, only: sort

    use ras, only: sort_orbitals

    implicit none

    abstract interface
        subroutine generate_all_excits_t(nI, n_excits, det_list) 
            use SystemData, only: nel 
            use constants, only: n_int
            integer, intent(in) :: nI(nel) 
            integer, intent(out) :: n_excits
            integer(n_int), intent(out), allocatable :: det_list(:,:)
        end subroutine generate_all_excits_t
    end interface

contains

    subroutine setup_arr_brr(in_lat) 
        class(lattice), intent(in) :: in_lat

        integer :: i 

        if (associated(arr)) deallocate(arr) 
        allocate(arr(nBasis,2))
        if (associated(brr)) deallocate(brr) 
        allocate(brr(nBasis))

        brr = [(i, i = 1, nBasis)]
        arr = 0.0_dp 

        do i = 1, nbasis 
            arr(i,:) = bhub * in_lat%dispersion_rel_orb(get_spatial(i))
        end do

        call sort(arr(1:nBasis,1), brr(1:nBasis), nskip = 2)
        call sort(arr(2:nBasis,1), brr(2:nBasis), nskip = 2)
! 
!         print *, "arr: " 
!         do i = 1, nBasis 
!             print *, arr(i,:) 
!         end do
!         print *, "brr: " 
!         do i = 1, nBasis
!             print *, brr(i)
!         end do

    end subroutine setup_arr_brr

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

        call sort(e_values)

    end function calc_eigenvalues

    subroutine eig(matrix, e_values, e_vectors, t_left_ev) 
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

!             call print_matrix(tmp_matrix)

!             print *, "bounds left_ev: ", lbound(left_ev,1),ubound(left_ev,1), &
!                 lbound(left_ev,2),ubound(left_ev,2)
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

!             print *, "bounds left_ev: ", lbound(left_ev,1),ubound(left_ev,1), &
!                 lbound(left_ev,2),ubound(left_ev,2)

!             print *, "Re(eval):", e_values
!             print *, "Im(eval):", dummy_eval(sort_ind)

!             print *, "left evectors: "
!             call print_matrix(left_ev(:,sort_ind))
!             print *, "right evectors: "
!             call print_matrix(right_ev(:,sort_ind))

!             print *, "tmp matrix:"
!             call print_matrix(tmp_matrix)

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

    end subroutine eig

    subroutine eig_sym(matrix, e_values, e_vectors) 
        real(dp), intent(in) :: matrix(:,:)
        real(dp), intent(out) :: e_values(size(matrix,1))
        real(dp), intent(out), optional :: e_vectors(size(matrix,1),size(matrix,2))
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

        if (sum(diff) < EPS) check_symmetric = .true.

    end function check_symmetric

    function create_hamiltonian(list_nI) result(hamil) 
        ! quite specific hamiltonian creation for my tests.. 
        integer, intent(in) :: list_nI(:,:) 
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j 

        hamil = h_cast(0.0_dp)

        do i = 1, size(list_nI,2)
            do j = 1, size(list_nI,2)
                hamil(i,j) = get_helement_lattice(list_nI(:,j),list_nI(:,i))
            end do
        end do

    end function create_hamiltonian

    function create_spin_dependent_hopping(list_nI, spin_opt) result(hamil)
        ! function to create a spin-dependent hopping matrix for the 
        ! exact tests of the spin-dependent hoppint transcorrelation
        ! is no spin (1 alpha, -1 beta) is given then alpha hopping is the 
        ! default
        integer, intent(in) :: list_nI(:,:)
        integer, intent(in), optional :: spin_opt
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j, spin, ex(2,2), ic
        integer(n_int) :: ilutI(0:NifTot), ilutJ(0:niftot)
        logical :: tpar

        hamil = h_cast(0.0_dp)

        if (present(spin_opt)) then 
            spin = spin_opt
        else 
            spin = 1
        end if

        do i = 1, size(list_nI,2)
            call EncodeBitDet(list_nI(:,i), ilutI)
            do j = 1, size(list_nI,2)
                call EncodeBitDet(list_nI(:,j), ilutJ)

                ic = findbitexcitlevel(ilutI,ilutJ)

                if (ic /= 1) cycle

                ex(1,1) = 1
                call GetBitExcitation(ilutI,ilutJ,ex,tpar)

                if (.not. same_spin(ex(1,1),ex(2,1))) cycle

                if (get_spin_pn(ex(1,1)) == spin) then
                    hamil(i,j) = get_helement_lattice(list_nI(:,j),list_nI(:,i))
                end if

            end do
        end do

    end function create_spin_dependent_hopping

    function create_hamiltonian_old(list_nI) result(hamil) 
        ! try to also create the hamiltonian with the old implementation.. 
        ! although i think there needs to be more setup done.. 
        integer, intent(in) :: list_nI(:,:)
        HElement_t(dp) :: hamil(size(list_nI,2),size(list_nI,2))

        integer :: i, j 

        t_lattice_model = .false.
        if (tGUGA) then
            call stop_all("create_hamiltonian_old", &
                "modify get_helement for GUGA")
        end if
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
                    write(iunit,*) matrix(i,:)
                end do
            else
                do i = lbound(matrix,1),ubound(matrix,1)
                    print *, matrix(i,:)
                end do
            end if
        end select


    end subroutine print_matrix

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

    function create_all_spin_flips(nI_in) result(spin_flips)
        ! takes a given spin-configuration in nI representation and 
        ! creates all possible states with flipped spin and same ms
        integer, intent(in) :: nI_in(:)
        integer, allocatable :: spin_flips(:,:)

        integer :: nI(size(nI_in))
        integer :: num, n_open, ms, n_states
        integer(n_int) :: ilutI(0:niftot)

        nI = nI_in

        num = size(nI)

        ! it is not necessary to have nI sorted at input, so do that now
        call sort_orbitals(nI)

        call EncodeBitDet(nI, ilutI)

        n_open = count_open_orbs(ilutI)
        ms = sum(get_spin_pn(nI))

        ! the number of possible spin distributions:
        n_states = int(choose(n_open, n_open/2 + ms))

        allocate(spin_flips(num,n_states))
        spin_flips = 0
        ! the first det will be the original one
        spin_flips(:,1) = nI

        call stop_all("here", "not yet implemented!")


    end function create_all_spin_flips

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

    subroutine run_excit_gen_tester(excit_gen, excit_gen_name, opt_nI, opt_n_iters, & 
            gen_all_excits) 
        use procedure_pointers, only: generate_excitation_t
        use util_mod, only: binary_search

        procedure(generate_excitation_t) :: excit_gen 
        integer, intent(in), optional :: opt_nI(nel), opt_n_iters
        procedure(generate_all_excits_t), optional :: gen_all_excits
        character(*), intent(in) :: excit_gen_name
        character(*), parameter :: this_routine = "run_excit_gen_tester"

        integer :: i, nI(nel), n_iters 
        integer :: default_n_iters = 100000
        integer :: default_n_dets = 1

        integer(n_int) :: ilut(0:niftot), tgt_ilut(0:niftot)
        integer :: nJ(nel), n_excits, ex(2,2), ic, ex_flag, i_unused = 0
        type(excit_gen_store_type) :: store 
        logical :: tPar, found_all
        real(dp) :: pgen, contrib
        HElement_t(dp) :: hel 
        integer(n_int), allocatable :: det_list(:,:)
        real(dp), allocatable :: contrib_list(:)
        logical, allocatable :: generated_list(:) 
        integer :: n_dets, n_generated, pos

        ASSERT(nel > 0)
        ! and also nbasis and stuff.. 
        ASSERT(nbasis > 0) 
        ASSERT(nel <= nbasis) 

        ! use some default values if not provided: 
        ! nel must be set! 
        if (present(opt_nI)) then 
            nI = opt_nI
        else 
            ! use HF-det as default
            nI = [(i, i = 1, nel)]
        end if

        if (present(opt_n_iters)) then 
            n_iters = opt_n_iters
        else 
            n_iters = default_n_iters
        end if

        ! i have to rewrite this routine, to a part which 
        ! creates all excitations and another who runs on possibly 
        ! multiple excitations! 
!         if (present(opt_n_dets)) then 
!             n_dets = opt_n_dets
!         else 
            n_dets = default_n_dets
!         end if

        ! the problem here is now.. we want to calulate all the possible 
        ! excitations.. i would need a general routine, which does that 
        ! with only the hamiltonian knowledge.. this is not really there.. 
        ! i should have a routine: 
        ! for this special setup, which is tested.. 

        ! for some special systems we should provide a routine to 
        ! calculate all the excitations (hubbard, UEG eg.) 
        if (present(gen_all_excits)) then 
            call gen_all_excits(nI, n_excits, det_list)
        else 
            call gen_all_excits_default(nI, n_excits, det_list)
        end if

        print *, "total possible excitations: ", n_excits
        do i = 1, n_excits 
            call writebitdet(6, det_list(:,i),.true.)
        end do

        ! call this below now for the number of specified determinants 
        ! (also use excitations of the first inputted, to be really 
        !   consistent) 

        n_dets = min(n_dets, n_excits) 

        print *, "---------------------------------"
        print *, "testing: ", excit_gen_name
        print *, "for ", n_dets, " determinants" 
        print *, " and ", n_iters, " iterations "

        call EncodeBitDet(nI, ilut) 

        print *, "for starting determinant: ", nI 

        ! Lists to keep track of things
        allocate(generated_list(n_excits))
        allocate(contrib_list(n_excits))
        generated_list = .false.
        contrib_list = 0

        n_generated = 0
        contrib = 0.0_dp

        do i = 1, n_iters
            if (mod(i, 1000) == 0) then 
                print *, i, "/" ,n_iters, " - ", contrib / real(n_excits*i,dp) 
            end if
            call excit_gen(nI, ilut, nJ, tgt_ilut, ex_flag, ic, ex, tpar, pgen, & 
                        hel, store) 

            if (nJ(1) == 0) cycle 
            call EncodeBitDet(nJ, tgt_ilut) 
            pos = binary_search(det_list, tgt_ilut, nifd+1)
            if (pos < 0) then 
                print *, "nJ: ", nJ 
                print *, "ilutJ:", tgt_ilut
                call stop_all(this_routine, 'Unexpected determinant generated')
            else 
                generated_list(pos) = .true. 
                n_generated = n_generated + 1 

                contrib = contrib + 1.0_dp / pgen 
                contrib_list(pos) = contrib_list(pos) + 1.0_dp / pgen 
            end if 
        end do

        print *, n_generated, " dets generated in ", n_iters, " iterations " 
        print *, 100.0_dp * (n_iters - n_generated) / real(n_iters,dp), "% abortion rate" 
        print *, "Averaged contribution: ", contrib / real(n_excits*n_iters,dp)

        ! check all dets are generated: 
        ASSERT(all(generated_list))

        print *, "=================================="
        print *, "Contribution List: "
        do i = 1, n_excits 
            ic = findbitexcitlevel(ilut, det_list(:,i))
            call writebitdet(6, det_list(:,i), .false.)
            print *, contrib_list(i)/real(n_iters,dp), "|", ic
        end do
        ! and check the uniformity of the excitation generation
!         ASSERT(all(abs(contrib_list / n_iters - 1.0_dp) < 0.01_dp))

    end subroutine run_excit_gen_tester

    subroutine create_hilbert_space(nI, n_states, state_list_ni, state_list_ilut, & 
            gen_all_excits_opt)
        ! a really basic routine, which creates the whole hilbert space based 
        ! on a input determinant and other quantities, like symmetry sectors, 
        ! already set outside the routine. for now this is specifically 
        ! implemented for the k-space hubbard model, since i still need to 
        ! test the transcorrelated approach there! 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_states
        integer, intent(out), allocatable :: state_list_ni(:,:) 
        integer(n_int), intent(out), allocatable :: state_list_ilut(:,:) 
        procedure(generate_all_excits_t), optional :: gen_all_excits_opt
        character(*), parameter :: this_routine = "create_hilbert_space"

        procedure(generate_all_excits_t), pointer :: gen_all_excits
        integer(n_int), allocatable :: excit_list(:,:), temp_list_ilut(:,:) 
        integer, allocatable :: temp_list_ni(:,:) 
        integer :: n_excits, n_total, tmp_n_states, cnt, i, j, pos
        integer(n_int) :: ilutI(0:niftot) 

        ! determine the type of system by the gen_all_excits routine 

        if (present(gen_all_excits_opt)) then 
            gen_all_excits => gen_all_excits_opt
        else
            gen_all_excits => gen_all_excits_default
        end if

        ! estimate the total number of excitations 
        n_total = int(choose(nBasis/2, nOccAlpha) * choose(nBasis/2,nOccBeta))

        n_states = 1 
        allocate(temp_list_ilut(0:niftot, n_total)) 
        allocate(temp_list_ni(nel, n_total)) 

        call EncodeBitDet(nI, ilutI)

        temp_list_ilut(:,1) = ilutI 
        temp_list_ni(:,1) = nI

        tmp_n_states = 1
        ! thats a really inefficient way to do it: 
        ! think of smth better at some point! 
        do while (.true.) 

            ! and i need a counter, which counts the number of added 
            ! excitations to the whole list.. i guess if this is 0 
            ! the whole hilbert space is reached.. 
            cnt = 0

            ! i need a temporary loop variable
            do i = 1, tmp_n_states 
                call gen_all_excits(temp_list_ni(:,i), n_excits, excit_list) 

                ! now i have to check if those states are already in the list
                do j = 1, n_excits 

                    pos = binary_search(temp_list_ilut(:,1:(tmp_n_states + cnt)), & 
                        excit_list(:,j), nifd+1)

                    ! if not yet found: 
                    if (pos < 0) then 
                        ! insert it at the correct place
                        ! does - pos give me the correct place then? 
                        pos = -pos
                        ! lets try.. and temp_list is always big enough i think..
                        ! first move
                        temp_list_ilut(:,(pos+1):tmp_n_states+cnt+1) = & 
                            temp_list_ilut(:,pos:(tmp_n_states+cnt))

                        temp_list_ni(:,(pos+1):(tmp_n_states+cnt+1)) = & 
                            temp_list_ni(:,pos:(tmp_n_states+cnt)) 
                        ! then insert 
                        temp_list_ilut(:,pos) = excit_list(:,j)
                        
                        call decode_bit_det(temp_list_ni(:,pos), excit_list(:,j))

                        ! and increase the number of state counter 
                        cnt = cnt + 1
                    else 
                        ! if already found i do not need to do anything i 
                        ! guess.. 
                    end if
                end do

            end do
            tmp_n_states = tmp_n_states + cnt 

            ! and somehow i need an exit criteria, if we found all the states.. 
            if (cnt == 0) exit
        end do

        n_states = tmp_n_states

        ! it should be already sorted or?? i think so.. 
        ! or does binary_search not indicate the position
        allocate(state_list_ni(nel,n_states), source = temp_list_ni(:,1:n_states))
        allocate(state_list_ilut(0:niftot,n_states), source = temp_list_ilut(:,1:n_states))

    end subroutine create_hilbert_space

    subroutine gen_all_excits_default(nI, n_excits, det_list) 
        use SymExcit3, only: CountExcitations3, GenExcitations3
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits 
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        character(*), parameter :: this_routine = "gen_all_excits_default"

        integer :: n_singles, n_doubles, n_dets, ex(2,2), ex_flag
        integer :: nJ(nel) 
        logical :: tpar, found_all 
        integer(n_int) :: ilut(0:niftot)

        n_excits = -1

        call EncodeBitDet(nI, ilut) 

        ! for reference in the "normal" case it looks like that: 
        call CountExcitations3(nI, 2, n_singles, n_doubles) 

        n_excits = n_singles + n_doubles 

        print *, "n_singles: ", n_singles
        print *, "n_doubles: ", n_doubles
        
        allocate(det_list(0:niftot,n_excits)) 
        n_dets = 0
        found_all = .false. 
        ex = 0
        ex_flag = 2 
        call GenExcitations3 (nI, ilut, nJ, ex_flag, ex, tpar, found_all, &
                              .false.)

        do while (.not. found_all)
            n_dets = n_dets + 1
            call EncodeBitDet (nJ, det_list(:,n_dets))

            call GenExcitations3 (nI, ilut, nJ, ex_flag, ex, tpar, &
                                  found_all, .false.)
        end do

        if (n_dets /= n_excits) then
            print *, "expected number of excitations: ", n_excits
            print *, "actual calculated ones: ", n_dets
            call stop_all(this_routine,"Incorrect number of excitations found")
        end if

        ! Sort the dets, so they are easy to find by binary searching
        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_excits_default


end module unit_test_helpers


