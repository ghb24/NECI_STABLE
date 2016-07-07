#include "macros.h"
module lanczos_general
    !   robert.anderson@kcl.ac.uk
    !   June 2016
    !
    !   Implementation of the algorithm documented in:
    !   http://dx.doi.org/10.1016/0024-3795(80)90167-6
    !
    use constants
    use FciMCData, only: hamiltonian, LanczosTag
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


    ! maximum size of subspace
    integer, parameter :: max_lanczos_vecs_default = 100
    ! the maximum difference between successive eigenvalues for each state
    ! required to trigger a loop exit before max_lanczos_vecs has been reached.
    real(dp), parameter :: convergence_error = 1.0e-10
    ! the algorithm is terminated when the norm of the last Lanczos vector
    ! reaches zero
    real(dp), parameter :: beta_epsilon = 1.0e-16
    ! If the overlap between any two Ritz vectors exceeds this value, a stop_all
    ! will be thrown
    real(dp), parameter :: orthog_tolerance = 5.0e-4

    type LanczosCalcType
        ! "super type"
        type(HamiltonianCalcType) :: super
        ! if we're not storing the subspace basis, we'll need:
        HElement_t(dp), allocatable :: v_1(:), v_k(:), v_k_minus_1(:)
        ! number of states to find
        integer :: n_states
        ! The beta elements which can't fit in T:
        HElement_t(dp) :: beta_0, beta_1
        ! work array for calculating each lanczos vector in turn
        HElement_t(dp), allocatable :: lanczos_vector(:)
        ! last estimate of each eigenvalue
        real(dp), allocatable :: T_eigenvalues_old(:)
        ! ritz values: exact eigenvalues of the tridiagonal matrix in super%projected_hamil
        real(dp), allocatable :: T_eigenvalues(:)
        ! exact eigenvectors of the tridiagonal matrix in super%projected_hamil
        HElement_t(dp), allocatable :: T_eigenvectors(:,:)
        ! approximate eivenvectors of the hamiltonian
        HElement_t(dp), allocatable :: ritz_vectors(:,:)
    end type

    contains

    ! We're storing alpha and beta arrays in T directly, so let's use getters 
    ! and setters to make this easier
    pure function getAlpha(this, i) result (val)
        type(LanczosCalcType), intent(in) :: this
        integer, intent(in) :: i
        HElement_t(dp) :: val
        val = this%super%projected_hamil(i,i)
    end function getAlpha

    subroutine setAlpha(this, i, val)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: i
        HElement_t(dp), intent(in) :: val
        this%super%projected_hamil(i,i) = val
    end subroutine setAlpha

    pure function getBeta(this, i) result (val)
        type(LanczosCalcType), intent(in) :: this
        integer, intent(in) :: i
        HElement_t(dp) :: val
        if (i>1) then
            val = this%super%projected_hamil(i,i-1)
            return
        elseif (i==0) then
            val = this%beta_0 
        else
            val = this%beta_1
        endif
    end function getBeta

    subroutine setBeta(this, i, val)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: i
        ! subdiagonal and superdiagonal elements of T are norms => always real
        real(dp), intent(in) :: val
        if (i>1) then
            this%super%projected_hamil(i-1,i) = cmplx(val)
            this%super%projected_hamil(i,i-1) = cmplx(val)
            return
        elseif (i==0) then
            this%beta_0 = val
        else
            this%beta_1 = val
        endif
    end subroutine setBeta

    subroutine addVec(this, k, vec)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: k
        HElement_t(dp), intent(in) :: vec(:)
        if (this%super%t_store_subspace_basis) then
            this%super%basis_vectors(:,k) = vec
        else
            this%v_k_minus_1 = this%v_k
            if (k==2) this%v_1 = this%v_k
            this%v_k = vec
        endif
    end subroutine addVec

    subroutine normaliseVec(this, k)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: k
        if (this%super%t_store_subspace_basis) then
            !associate(v=>this%super%basis_vectors(:,k))
            !v = v /sqrt(dble(dot_product(h_conjg(v), v)))
            this%super%basis_vectors(:,k) = this%super%basis_vectors(:,k) / &
                sqrt(dble(dot_product(h_conjg(this%super%basis_vectors(:,k)), this%super%basis_vectors(:,k))))
!            end associate
        else
            associate(v=>this%v_k)
            v = v /sqrt(dble(dot_product(h_conjg(v), v)))
            end associate
        endif
    end subroutine normaliseVec

    subroutine InitLanczosCalc(this, print_info, hamil_type, n_states, max_lanczos_vecs, t_store_subspace_basis)
        type(LanczosCalcType), intent(out) :: this
        integer, intent(in) :: hamil_type, n_states, max_lanczos_vecs
        logical, intent(in) :: t_store_subspace_basis
        character (len=*), parameter :: t_r = "init_lanczos"
        logical, intent(in) :: print_info

        integer :: i, HFindex, mem_reqd, residual_mem_reqd, ierr
        integer(MPIArg) :: mpi_temp
        real(dp), allocatable :: hamil_diag_temp(:)
        real(dp), allocatable :: lowest_energies(:)
        integer, allocatable :: lowest_energy_det_indices(:)

        call InitHamiltonianCalc(this%super, print_info, hamil_type, max_lanczos_vecs, t_store_subspace_basis)
        associate(&
            space_size => this%super%space_size, &
            max_subspace_size => this%super%max_subspace_size &
        )

        if (n_states > space_size) then
            call stop_all(t_r, "Not enough determinants in the space to produce the required number of approximate eigenpairs")
        endif

        this%n_states = n_states
        ! these will be larger than neccessary until the final iteration:
        safe_malloc(this%T_eigenvalues, (max_subspace_size))
        safe_malloc(this%T_eigenvalues_old, (n_states))
        safe_malloc(this%T_eigenvectors, (max_subspace_size, max_subspace_size))

        safe_calloc(this%lanczos_vector, (space_size), 0.0_dp)

        if (.not.t_store_subspace_basis) then
            safe_malloc(this%v_1, (max_subspace_size))
            safe_malloc(this%v_k, (max_subspace_size))
            safe_malloc(this%v_k_minus_1, (max_subspace_size))
        endif

        ! check for previous allocation of the eigenvector estimates
        if (allocated(this%ritz_vectors)) then
            deallocate(this%ritz_vectors, stat=ierr)
            call logmemdealloc(t_r, lanczosTag, ierr)
        end if
        ! nstates columns, holding approx eigenvectors of length space_size
        safe_calloc_e(this%ritz_vectors, (n_states, space_size), 0.0_dp, ierr)
        call logmemalloc("ritz_vectors", space_size*n_states, 8, t_r, lanczosTag, ierr)

        safe_malloc(hamil_diag_temp, (size(hamil_diag)))
        hamil_diag_temp = -hamil_diag
        safe_calloc(lowest_energies, (n_states), 0.0_dp)
        safe_calloc(lowest_energy_det_indices, (n_states), 0)
        lowest_energy_det_indices(1) = maxloc(hamil_diag_temp,1)
        lowest_energies(1) = hamil_diag(lowest_energy_det_indices(1))
        do i = 2, n_states
            lowest_energy_det_indices(i) = maxloc(hamil_diag_temp,1,mask=(hamil_diag_temp)<-lowest_energies(i-1))
            lowest_energies(i) = hamil_diag(lowest_energy_det_indices(i))
        enddo
        safe_free(hamil_diag_temp)

        ! If there is only one determinant per state in the space being diagonalised:
        if (space_size == n_states) then
            if (iProcIndex == root) then
                this%T_eigenvalues(:) = lowest_energies(:)
                do i = 1, n_states
                    this%T_eigenvectors(i,lowest_energy_det_indices(i)) = 1.0_dp
                enddo
            endif
            call MPIBCast(this%T_eigenvalues)
            this%super%skip_calc = .true.
            return
        end if

        ! with the memory management done, now we can start setting up the pre-iteration
        ! numerics. Set the first basis vector to be an equal, normalised superposition
        ! of the nstates lowest lying determinants
        do i = 1, n_states
            this%lanczos_vector(lowest_energy_det_indices(i)) = 1.0_dp/sqrt(dble(n_states))
        enddo

        debug_line(8,this%lanczos_vector)


        call addVec(this, 1, this%lanczos_vector)
        call normaliseVec(this, 1)

        call project_hamiltonian_lanczos(this, 1)   

        safe_free(lowest_energies)
        safe_free(lowest_energy_det_indices)
        end associate

    end subroutine InitlanczosCalc

    subroutine FreeLanczosCalc(this)
        ! deallocate all working arrays, but leave results
        use hamiltonian_linalg, only: DestroyHamiltonianCalc
        type(LanczosCalcType), intent(inout) :: this
        ! destroy the super type instance
        safe_free(this%v_1)
        safe_free(this%v_k)
        safe_free(this%v_k_minus_1)
        safe_free(this%T_eigenvalues_old)
        safe_free(this%T_eigenvectors)
        call DestroyHamiltonianCalc(this%super)
        safe_free(this%lanczos_vector)
    end subroutine

    subroutine DestroyLanczosCalc(this)
        type(LanczosCalcType), intent(inout) :: this
        call FreeLanczosCalc(this)
        ! deallocate results as well
        safe_free(this%T_eigenvalues)
        safe_free(this%ritz_vectors)
    end subroutine DestroyLanczosCalc

    subroutine perform_orthogonality_test(this)
        type(LanczosCalcType), intent(inout) :: this
        integer :: i, j
        HElement_t(dp) :: overlap, normi, normj
        character(*), parameter :: t_r = "perform_orthogonality_test"
        if (this%n_states>1) then
            do i = 2, this%n_states
                do j = 1, i-1
                    overlap = dot_product(h_conjg(this%ritz_vectors(i,:)), this%ritz_vectors(j,:))
                    if (abs(overlap) > orthog_tolerance) then
                        call stop_all(t_r, "Ritz vector overlap exceeds tolerance. n_states possibly too large for space_size") 
                    endif
                enddo
            enddo
            write(6, '(" Ritz vectors are mutually orthogonal to a tolerance of "5ES16.7)') orthog_tolerance
        endif
    end subroutine perform_orthogonality_test

    subroutine project_hamiltonian_lanczos(this, basis_index)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: basis_index

        if (this%super%t_store_subspace_basis) then
            if (iprocindex==root) then
                call multiply_hamil_and_vector(this%super, this%super%basis_vectors(:,basis_index), this%lanczos_vector)
            else
                call multiply_hamil_and_vector(this%super, this%super%basis_vectors(:,basis_index), this%super%temp_out)
            endif
        else
            if (iprocindex==root) then
                call multiply_hamil_and_vector(this%super, this%v_k, this%lanczos_vector)
            else
                call multiply_hamil_and_vector(this%super, this%v_k, this%super%temp_out)
            endif
        endif
    end subroutine project_hamiltonian_lanczos

    subroutine diagonalise_tridiagonal(this, N, t_calc_eigenvectors)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: N
        logical, intent(in) :: t_calc_eigenvectors
        integer :: lwork, info
        HElement_t(dp), allocatable :: cwork(:)
        real(dp), allocatable :: rwork(:)
        character :: jobz

        lwork = max(1,3*N+1)
        safe_malloc(rwork, (lwork))
        if (t_calc_eigenvectors) then
            jobz = 'V'
        else
            jobz = 'N'
        endif

        this%T_eigenvalues_old(1:this%n_states) = this%T_eigenvalues(1:this%n_states)
        this%super%projected_hamil_work(1:N,1:N) = this%super%projected_hamil(1:N,1:N)
#ifdef __CMPLX
        call pretty_print(8, this%super%projected_hamil_work(1:N,1:N))
        safe_malloc(cwork, (lwork))
        call zheev(&
             jobz, &
            'U', &
             N , &
             this%super%projected_hamil_work(1:N,1:N), &
             N, &
             this%T_eigenvalues(1:N), &
             cwork, &
             lwork, &
             rwork, &
             info &
         )
        safe_free(cwork)
#else
        call dsyev(&
             jobz, &
            'U', &
             N , &
             this%super%projected_hamil_work(1:N,1:N), &
             N, &
             this%T_eigenvalues(1:N), &
             rwork, &
             lwork, &
             info &
         )
#endif
        safe_free(rwork)
        if (t_calc_eigenvectors) then
            ! move the eigenvectors out of the working array
            this%T_eigenvectors(1:this%n_states,1:N) = &
                this%super%projected_hamil_work(1:this%n_states,1:N)
        endif

    end subroutine diagonalise_tridiagonal

    function check_deltas(this) result (t)
        type(LanczosCalcType), intent(inout) :: this
        logical :: t
        t = .not.any(abs(this%T_eigenvalues(1:this%n_states) - this%T_eigenvalues_old(1:this%n_states))>convergence_error)
    end function check_deltas

    subroutine compute_ritz_vectors(this, k)
        ! now we project the eigenvalues of T into the full space to obtain
        ! estimates for the state vectors
        !   rows, columns:
        ! ritz_vectors(space_size, n_states) = basis_vectors(space_size, k) x T_eigenvectors(k, n_states)
        !   we want to do (column-major):
        ! ritz_vectors(n_states, space_size) = basis_vectors(k, space_size) x T_eigenvectors(n_states, k)
        !   but we have:
        ! basis_vectors(space_size, k), and T_eigenvectors(n_states, k)
        !   so only transpose basis_vectors
        
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: k
        
        associate(space_size => this%super%space_size)

        if (this%super%t_store_subspace_basis) then
#ifdef __CMPLX
            call zgemm(&
                'T', &
                'N', &
                space_size, &
                this%n_states, &
                k, &
                1.0_dp, &
                this%super%basis_vectors(1:space_size,1:k), &
                k, &
                this%T_eigenvectors(1:this%n_states, 1:k), &
                k, &
                0.0_dp, &
                this%ritz_vectors(1:this%n_states, 1:space_size), &
                space_size &
            )
#else
            call dgemm(&
                'T', &
                'N', &
                space_size, &     ! rows in basis_vectors**T (columns in basis_vectors) => rows in ritz_vectors
                this%n_states, &  ! columns in T_eigenvectors => columns in ritz_vectors
                k, &              ! columns in basis_vectors, and rows in T_eigenvectors
                1.0_dp, &
                this%super%basis_vectors(1:space_size,1:k), &           
                k, &
                this%T_eigenvectors(1:this%n_states, 1:k), &
                k, &
                0.0_dp, &
                this%ritz_vectors(1:this%n_states, 1:space_size), &
                space_size &
            )
#endif
        else
            ! TODO: must reconstitute the basis vectors one at a time
        endif
        end associate

    end subroutine compute_ritz_vectors

    subroutine perform_lanczos(this, n_states, t_store_subspace_basis, hamil_type, print_info_in)
        use hamiltonian_linalg, only : orthogonalise_against_previous_basis_vectors
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: n_states, hamil_type
        logical, intent(in) :: t_store_subspace_basis, print_info_in
        logical :: print_info
        integer :: k, i, ierr
        real(sp) :: start_time, end_time
        HElement_t(dp) :: overlap

        ! Only let the root processor print information.
        print_info = print_info_in .and. (iProcIndex == root)

        if (.not. allocated(this%lanczos_vector)) then
            call InitLanczosCalc(this, print_info, hamil_type, n_states, max_lanczos_vecs_default, t_store_subspace_basis)
        endif

        if (print_info) write(6,'(1X,"Perfoming a Lanczos Diagonalisation of the trial space")')
        if (print_info) write(6,'(/,1X,"Iteration",4x,"State",12X,"Energy",7X,"Time")'); call neci_flush(6)
        associate(&
            basis_vectors => this%super%basis_vectors, &
            v_k => this%v_k, &
            v_k_minus_1 => this%v_k_minus_1 &
        )

        do k = 1, this%super%max_subspace_size-1


            if (this%super%skip_calc) exit
            call cpu_time(start_time)

            if (iprocindex==root) then
                if (t_store_subspace_basis) then
                    !debug_line(8, this%lanczos_vector)
                    !debug_line(8, basis_vectors(:,k))
                    overlap = dot_product(h_conjg(basis_vectors(:,k)), this%lanczos_vector)
                    debug_line(8, overlap)
                    call setAlpha(this, k, cmplx(real(overlap)))
                    debug_line(8, getAlpha(this, k))
                    call addVec(this, k+1, this%lanczos_vector-getAlpha(this, k)*basis_vectors(:,k))
                    debug_line(8, h_conjg(basis_vectors(:,k+1)))
                    debug_line(8, basis_vectors(:,k+1))

                    call setBeta(this, k+1, sqrt(dble(dot_product(h_conjg(basis_vectors(:,k+1)), basis_vectors(:,k+1)))))
                    if (abs(getBeta(this, k+1))<beta_epsilon) then
                        write(6, *) " Lanczos vector norm reached zero"
                        exit
                    endif
                    debug_line(8, (dot_product(h_conjg(basis_vectors(:,k+1)), basis_vectors(:,k+1))))
                    basis_vectors(:,k+1) = basis_vectors(:,k+1) / getBeta(this, k+1)
                    ! Although Lanczos vectors computed in exact arithmetic are orthogonal,
                    ! those generated on a machine will stray increasingly from orthogonal with
                    ! each iteration: so do GS procedure here
                    write(8,*)
                    debug_line(8, basis_vectors(:,k+1))
                    call orthogonalise_against_previous_basis_vectors(this%super, k+1)
                    debug_line(8, basis_vectors(:,k+1))
                    write(8,*)
                else
                    overlap = dot_product(h_conjg(this%v_k), this%lanczos_vector)
                    call setAlpha(this, k, overlap)
                    call addVec(this, k+1, this%lanczos_vector-getAlpha(this, k)*v_k)
                    call setBeta(this, k+1, sqrt(dble(dot_product(h_conjg(v_k), v_k))))
                    v_k = v_k / getBeta(this, k+1)
                endif
            endif
            call MPIBarrier(ierr)
            call project_hamiltonian_lanczos(this, k+1)

            call cpu_time(end_time)

            if (print_info) then
                do i=1, this%n_states
                    write(6,'(8X,i4,3X,i2,2x,f16.10,2x,f9.3)') k, i, this%T_eigenvalues(i), end_time-start_time
                    call neci_flush(6)
                enddo
                write(6,*)
            endif

            if (iprocindex==root) then
                if (t_store_subspace_basis) then
                    this%lanczos_vector = this%lanczos_vector - getBeta(this, k+1)*basis_vectors(:,k)
                else
                    this%lanczos_vector = this%lanczos_vector - getBeta(this, k+1)*v_k_minus_1
                endif

                call diagonalise_tridiagonal(this, k, .false.)
                if (check_deltas(this)) then
                    if (print_info) then
                        write(6, '(i2" eigenvalues(s) were successfully converged to within "5ES16.7)') this%n_states, convergence_error
                        call neci_flush(6)
                    endif
                    exit
                endif
            endif
        end do
        
        if (iprocindex==root) then
            call diagonalise_tridiagonal(this, k-1, .true.)
            call compute_ritz_vectors(this, k)
            call perform_orthogonality_test(this)
        endif

        if (print_info) then
            write(6,'(/,1x,"Final calculated energies:")')
            do i=1, this%n_states
                write(6,'(1x,"State",1X,i2,4X,f16.10)') i, this%T_eigenvalues(i)
                call neci_flush(6)
            enddo
        endif

        ! wait for final diagonalisation before freeing any memory
        call MPIBarrier(ierr)
        call FreeLanczosCalc(this)
        end associate
    end subroutine perform_lanczos

    subroutine pretty_print(out_unit, mat)
        integer :: out_unit, shp(2), i
        HElement_t(dp) :: mat(:,:)
        shp = shape(mat)
        do i = 1, shp(2)
            write(out_unit, *) real(mat(:,i))
        enddo
        write(out_unit, *)
    end subroutine

end module lanczos_general
