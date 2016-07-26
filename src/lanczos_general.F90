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
        direct_ci_out, &
        inner_product, &
        euclidean_norm_square, &
        euclidean_norm


    ! maximum size of subspace
    integer, parameter :: max_lanczos_vecs_default = 40
    ! maximum number of main loop restarts
    integer, parameter :: max_restarts = 10
    ! the maximum difference between successive eigenvalues for each state
    ! required to trigger a loop exit before max_lanczos_vecs has been reached.
    real(dp), parameter :: convergence_error = 1.0e-10
    ! the algorithm is terminated when the norm of the last Lanczos vector
    ! reaches zero
    real(dp), parameter :: beta_threshold = 1.0e-7
    ! If the overlap between any two Ritz vectors exceeds this value, a stop_all
    ! will be thrown
    real(dp), parameter :: orthog_tolerance = 1.0e-5

    type LanczosCalcType
        ! "super type" for common hamiltonian data
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
        real(dp), allocatable :: ritz_values_old(:)
        ! ritz values: exact eigenvalues of the tridiagonal matrix in super%projected_hamil
        real(dp), allocatable :: ritz_values(:)
        ! exact eigenvectors of the tridiagonal matrix in super%projected_hamil
        HElement_t(dp), allocatable :: T_eigenvectors(:,:)
        ! approximate eigenvectors of the hamiltonian
        HElement_t(dp), allocatable :: ritz_vectors(:,:)
        ! final eigenvalue estimates
        real(dp), allocatable :: eigenvalues(:)
        ! final eigenvector estimates
        HElement_t(dp), allocatable :: eigenvectors(:,:)
        ! keep track of the states which are already converged
        logical, allocatable :: t_states_converged(:)
    end type

    contains

    ! We're storing alpha and beta arrays in T directly, so let's use getters 
    ! and setters to make this easier
    pure function getAlpha(this, i) result (val)
        type(LanczosCalcType), intent(in) :: this
        integer, intent(in) :: i
        real(dp) :: val
        val = this%super%projected_hamil(i,i)
    end function getAlpha

    pure subroutine setAlpha(this, i, val)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: i
        real(dp), intent(in) :: val
        this%super%projected_hamil(i,i) = val
    end subroutine setAlpha

    pure function getBeta(this, i) result (val)
        type(LanczosCalcType), intent(in) :: this
        integer, intent(in) :: i
        real(dp) :: val
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
            this%super%projected_hamil(i-1,i) = val
            this%super%projected_hamil(i,i-1) = val
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

    subroutine InitLanczosCalc(this, print_info, hamil_type, n_states, max_lanczos_vecs, t_store_subspace_basis, t_orthogonalise)
        type(LanczosCalcType), intent(out) :: this
        integer, intent(in) :: hamil_type, n_states, max_lanczos_vecs
        logical, intent(in) :: t_store_subspace_basis, t_orthogonalise
        character (len=*), parameter :: t_r = "init_lanczos"
        logical, intent(in) :: print_info

        integer :: i, HFindex, mem_reqd, residual_mem_reqd, ierr
        integer(MPIArg) :: mpi_temp
        real(dp), allocatable :: hamil_diag_temp(:)
        real(dp), allocatable :: lowest_energies(:)
        integer, allocatable :: lowest_energy_det_indices(:)

        call InitHamiltonianCalc(this%super, print_info, hamil_type, max_lanczos_vecs, t_store_subspace_basis, t_orthogonalise)
        associate(&
            space_size => this%super%space_size, &
            max_subspace_size => this%super%max_subspace_size &
        )

        if (n_states > space_size) then
            call stop_all(t_r, "Not enough determinants in the space to produce the required number of approximate eigenpairs")
        endif

        this%n_states = n_states
        ! these will be larger than neccessary until the final iteration:
        safe_malloc(this%ritz_values, (max_subspace_size))
        safe_calloc(this%eigenvalues, (n_states), 0.0_dp)
        safe_calloc(this%t_states_converged, (n_states), .false.)
        safe_malloc(this%ritz_values_old, (n_states))
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
        safe_calloc_e(this%ritz_vectors, (space_size, n_states), 0.0_dp, ierr)
        call logmemalloc("ritz_vectors", space_size*n_states, 8, t_r, lanczosTag, ierr)

        safe_calloc_e(this%eigenvectors, (space_size, n_states), 0.0_dp, ierr)
        call logmemalloc("eigenvectors", space_size*n_states, 8, t_r, lanczosTag, ierr)

        if (iProcIndex==root) then

            safe_malloc(hamil_diag_temp, (size(hamil_diag)))
            hamil_diag_temp(:) = -hamil_diag(:)
            safe_calloc(lowest_energies, (n_states), 0.0_dp)
            safe_calloc(lowest_energy_det_indices, (n_states), 0)
            lowest_energy_det_indices(1) = maxloc(hamil_diag_temp,dim=1)
            lowest_energies(1) = hamil_diag(lowest_energy_det_indices(1))
            do i = 2, n_states
                lowest_energy_det_indices(i) = maxloc(hamil_diag_temp,dim=1,mask=hamil_diag_temp<-lowest_energies(i-1))
                lowest_energies(i) = hamil_diag(lowest_energy_det_indices(i))
            enddo
            safe_free(hamil_diag_temp)
        endif

        ! If there is only one determinant per state in the space being diagonalised:
        if (space_size == n_states) then
            if (iProcIndex == root) then
                this%eigenvalues(:) = lowest_energies(:)
                do i = 1, n_states
                    this%eigenvectors(lowest_energy_det_indices(i), i) = 1.0_dp
                enddo
            endif
            this%super%skip_calc = .true.
            return
        end if

        if (iProcIndex==root) then
            ! with the memory management done, now we can start setting up the pre-iteration
            ! numerics. Set the first basis vector to be an equal, normalised superposition
            ! of the nstates lowest lying determinants

            do i = 1, n_states
                this%lanczos_vector(lowest_energy_det_indices(i)) = 1.0_dp/sqrt(dble(n_states))
            enddo

            safe_free(lowest_energies)
            safe_free(lowest_energy_det_indices)
        endif
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
        safe_free(this%ritz_values_old)
        safe_free(this%T_eigenvectors)
        safe_free(this%ritz_values)
        safe_free(this%ritz_vectors)
        safe_free(this%t_states_converged)
        call DestroyHamiltonianCalc(this%super)
        safe_free(this%lanczos_vector)
    end subroutine

    subroutine DestroyLanczosCalc(this)
        type(LanczosCalcType), intent(inout) :: this
        call FreeLanczosCalc(this)
        ! deallocate results as well
        safe_free(this%eigenvalues)
        safe_free(this%eigenvectors)
    end subroutine DestroyLanczosCalc

    subroutine perform_orthogonality_test(this)
        type(LanczosCalcType), intent(inout) :: this
        integer :: i, j
        real(dp) :: overlap, largest_overlap
        character(*), parameter :: t_r = "perform_orthogonality_test"
        largest_overlap = 0.0_dp
        if (this%n_states>1) then
            do i = 2, this%n_states
                do j = 1, i-1
                    overlap = inner_product(this%ritz_vectors(:,i), this%ritz_vectors(:,j))

                    if (abs(overlap) > largest_overlap) then
                        largest_overlap = abs(overlap)
                    endif
                    if (abs(overlap) > orthog_tolerance) then
                        write(6, '(" Largest Ritz vector overlap: "5ES16.7)') largest_overlap
                        call stop_all(t_r, "Ritz vector overlap exceeds tolerance. n_states possibly too large for space_size") 
                    endif
                enddo
            enddo
            write(6, '(" Largest Ritz vector overlap: "5ES16.7)') largest_overlap
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
        real(dp), allocatable :: rwork(:)
        character :: jobz

        lwork = max(1,3*N+1)
        safe_malloc(rwork, (lwork))
        if (t_calc_eigenvectors) then
            jobz = 'V'
        else
            jobz = 'N'
        endif

        if (.not. t_calc_eigenvectors) then
            this%ritz_values_old(1:this%n_states) = this%ritz_values(1:this%n_states)
            this%super%projected_hamil_work(1:N,1:N) = this%super%projected_hamil(1:N,1:N)
        endif

        call dsyev(&
             jobz, &
            'U', &
             N , &
             this%super%projected_hamil_work(1:N,1:N), &
             N, &
             this%ritz_values(1:N), &
             rwork, &
             lwork, &
             info &
        )
        safe_free(rwork)

        if (t_calc_eigenvectors) then
            ! move the eigenvectors out of the working array
            this%T_eigenvectors(1:N,1:this%n_states) = &
                cmplx(this%super%projected_hamil_work(1:N, 1:this%n_states))
        endif
    end subroutine diagonalise_tridiagonal

    function check_deltas(this, k) result (t)
        type(LanczosCalcType), intent(inout) :: this
        integer :: k
        logical :: t
        if (k<=this%n_states) then
            t = .false.
            return
        endif
        t=.not.any(abs(this%ritz_values(1:this%n_states)-this%ritz_values_old(1:this%n_states))>convergence_error)
    end function check_deltas

    function check_delta(this, k, i_state) result (t)
        type(LanczosCalcType), intent(inout) :: this
        integer :: k, i_state
        logical :: t
        if (k<=this%n_states) then
            t = .false.
            return
        endif
        t = (abs(this%ritz_values(i_state)-this%ritz_values_old(i_state))<convergence_error)
    end function check_delta

    subroutine compute_ritz_vectors(this, k)
        ! now we project the eigenvalues of T into the full space to obtain
        ! estimates for the state vectors
        !   rows, columns:                                       M                                N
        ! ritz_vectors(n_states, space_size) = T_eigenvectors(n_states, k) x basis_vectors(k, space_size)
        !   we want to do (column-major):
        ! ritz_vectors(space_size, n_states) = T_eigenvectors(k, n_states) x basis_vectors(space_size, k)
        
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: k
        
        associate(space_size => this%super%space_size)

        if (this%super%t_store_subspace_basis) then
            ! for each state, find the ritz vector
            ! TODO: do this with BLAS
            do i = 1, this%n_states
                do j = 1, space_size
                    this%ritz_vectors(j,i) = inner_product(this%super%basis_vectors(j,1:k), this%T_eigenvectors(1:k,i))
                enddo
            enddo
        else
            ! TODO: must reconstitute the basis vectors one at a time
        endif
        end associate
    end subroutine compute_ritz_vectors

    function get_rayleigh_quotient(this, i_state) result (exp_val)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: i_state
        integer :: ierr
        HElement_t(dp), allocatable :: H_ket(:)
        real(dp) :: exp_val

        safe_malloc(H_ket, (this%super%space_size))

        call MPIBCast(this%eigenvectors)
        if (iprocindex==root) then
            call multiply_hamil_and_vector(this%super, this%eigenvectors(1:this%super%space_size, i_state), H_ket)
        else
            call multiply_hamil_and_vector(this%super, this%eigenvectors(1:this%super%space_size, i_state), this%super%temp_out)
        endif
        if (iprocindex==root) then
            exp_val = real(inner_product(this%eigenvectors(1:this%super%space_size, i_state), H_ket), dp)
            exp_val = exp_val/euclidean_norm_square(this%eigenvectors(1:this%super%space_size, i_state))
        else
            exp_val=0.0_dp
        endif
        safe_free(H_ket)
    end function get_rayleigh_quotient

    function compute_residual_norm(this, i_state) result (norm)
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: i_state
        integer :: ierr
        HElement_t(dp), allocatable :: H_ket(:)
        real(dp) :: norm

        safe_malloc(H_ket, (this%super%space_size))

        call MPIBCast(this%eigenvectors)
        call MPIBCast(this%eigenvalues)
        if (iprocindex==root) then
            call multiply_hamil_and_vector(this%super, this%eigenvectors(1:this%super%space_size, i_state), H_ket)
        else
            call multiply_hamil_and_vector(this%super, this%eigenvectors(1:this%super%space_size, i_state), this%super%temp_out)
        endif
        norm = euclidean_norm(H_ket-this%eigenvalues(i_state)*this%eigenvectors(1:this%super%space_size, i_state))
        safe_free(H_ket)
    end function

    subroutine perform_lanczos(this, n_states, t_store_subspace_basis, hamil_type, print_info_in, t_orthogonalise)
        use hamiltonian_linalg, only : orthogonalise_against_previous_basis_vectors
        type(LanczosCalcType), intent(inout) :: this
        integer, intent(in) :: n_states, hamil_type
        logical, intent(in) :: t_store_subspace_basis, print_info_in, t_orthogonalise
        logical :: print_info, t_deltas_pass
        integer :: k, i, ierr, nplaces, delta_check_outcome, restart_counter
        real(sp) :: start_time, end_time
        real(dp) :: exp_val, overlap
        character(40) :: main_output_fmt, final_output_fmt
        
        ! format specifier formatting!
        ! ensure that enough places are displayed to show convergence to the desired accuracy
        nplaces = -int(log10(convergence_error))
        write(main_output_fmt, '("(8X,i4,3X,i2,2x,f",i2,".",i2,",2x,f9.3)")') nplaces+7, nplaces
        write(final_output_fmt, '("(1x,",A7,",i2,4X,f",i2,".",i2,")")') '"State"', nplaces+7, nplaces

        ! Only let the root processor print information.
        print_info = print_info_in .and. (iProcIndex == root)

        if (.not. allocated(this%lanczos_vector)) then
            call InitLanczosCalc(this, print_info, hamil_type, n_states, &
                max_lanczos_vecs_default, t_store_subspace_basis, t_orthogonalise)
        endif

        if (print_info) then
            write(6,'(1X,"Perfoming a Lanczos Diagonalisation of the trial space")')
            if (this%super%t_orthogonalise) then
                write(6,'(/,1X,"Orthogonalising Lanczos vectors")')
            else
                write(6,'(/,1X,"Not orthogonalising Lanczos vectors")')
            endif
        endif
        associate(&
            basis_vectors => this%super%basis_vectors, &
            v_k => this%v_k, &
            v_k_minus_1 => this%v_k_minus_1 &
        )
    
        ! start Lanczos restart loop
        restart_counter = 0
        do
            if (iprocindex==root) then
                call addVec(this, 1, this%lanczos_vector)
            endif
            call MPIBCast(this%lanczos_vector)
            call MPIBCast(this%super%basis_vectors(:,1))
            ! set lanczos vector 1 to H*v_1
            call project_hamiltonian_lanczos(this, 1)   

            if (print_info) then
                write(6,'(/,1X,"Iteration",4x,"State",12X,"Energy",7X,"Time")'); call neci_flush(6)
            endif
            ! start Lanczos main loop
            do k = 1, this%super%max_subspace_size-1
                if (this%super%skip_calc) exit

                call cpu_time(start_time)

                if (iprocindex==root) then
                    if (t_store_subspace_basis) then
                        overlap = real(inner_product(basis_vectors(:,k), this%lanczos_vector),dp)
                        call setAlpha(this, k, overlap)
                        call addVec(this, k+1, this%lanczos_vector-getAlpha(this, k)*basis_vectors(:,k))
                        call setBeta(this, k+1, euclidean_norm(basis_vectors(:,k+1)))
                        if (abs(getBeta(this, k+1))<beta_threshold) then
                            write(6, *) " Lanczos vector norm reached zero"
                            exit
                        endif
                        basis_vectors(:,k+1) = basis_vectors(:,k+1) / getBeta(this, k+1)
                    else
                        overlap = real(dot_product(this%v_k, this%lanczos_vector),dp)
                        call setAlpha(this, k, overlap)
                        call addVec(this, k+1, this%lanczos_vector-getAlpha(this, k)*v_k)
                        call setBeta(this, k+1, sqrt(dble(dot_product(v_k, v_k))))
                        v_k = v_k / getBeta(this, k+1)
                    endif
                endif
                call MPIBarrier(ierr)
                call MPIBCast(this%super%basis_vectors(:,k+1))
                if (this%super%t_orthogonalise) then
                    ! Although Lanczos vectors computed in exact arithmetic are orthogonal,
                    ! those generated on a machine will stray increasingly from orthogonal with
                    ! each iteration: so do GS procedure here
                    call orthogonalise_against_previous_basis_vectors(this%super, k+1)
                endif
                call project_hamiltonian_lanczos(this, k+1)

                call cpu_time(end_time)

                if (iprocindex==root) then
                    if (t_store_subspace_basis) then
                        this%lanczos_vector = this%lanczos_vector - getBeta(this, k+1)*basis_vectors(:,k)
                    else
                        this%lanczos_vector = this%lanczos_vector - getBeta(this, k+1)*v_k_minus_1
                    endif

                    call diagonalise_tridiagonal(this, k, .false.)
                    ! if all states have converged we can end the main loop prematurely
                    t_deltas_pass = .true.
                    do i=1, this%n_states
                        t_deltas_pass = t_deltas_pass .and. check_delta(this, k, i)
                        if (.not.this%t_states_converged(i)) then
                            if (print_info) then
                                write(6,trim(main_output_fmt)) k, i, this%ritz_values(i), end_time-start_time
                                call neci_flush(6)
                            endif
                        endif
                    enddo
                    write(6,*)
                endif

                call MPIBCast(t_deltas_pass)
                if (t_deltas_pass) then
                    exit
                endif
                call MPIBCast(this%lanczos_vector)
            ! end Lanczos main loop
            end do
            
            if (iprocindex==root) then
                call diagonalise_tridiagonal(this, k, .true.)
                call compute_ritz_vectors(this, k)
            endif

            ! we now evaluate the quality of the convergence achieved. If a state is sufficiently 
            ! well converged, save it and its energy otherwise, add the current best ritz vector
            ! for the state to the restarting Lanczos vector such that we restart with an equal
            ! superposition of the best guesses for the unconverged states.
            if (iProcIndex==root) then
                this%lanczos_vector(:) = h_cast(0.0_dp)
                do i = 1, this%n_states
                    if (this%t_states_converged(i)) then
                        cycle
                    endif
                    if (.not. check_delta(this, k, i)) then
                        this%lanczos_vector = this%lanczos_vector + this%ritz_vectors(:,i)
                    else
                        this%eigenvectors(:,i) = this%ritz_vectors(:,i)
                        this%eigenvalues(i) = this%ritz_values(i)
                        this%t_states_converged(i) = .true.
                    endif
                enddo
                this%lanczos_vector = this%lanczos_vector / euclidean_norm(this%lanczos_vector)
            endif
            call MPIBCast(this%t_states_converged)

            if (.not.any(.not.this%t_states_converged)) then
                exit
            endif

            restart_counter = restart_counter + 1
            if (restart_counter<max_restarts) then
                if(print_info) then
                    write(6,'(/,1x,"Maximum iteration number reached, restarting Lanczos algorithm (",i3,")")') &
                        restart_counter
                endif
            else
                if(print_info) then
                    write(6,'(/,1x,"Maximum restart number reached. Some states may not be converged.")')
                endif
                if (iprocindex==root) then
                    do i = 1, this%n_states
                        if (.not. this%t_states_converged(i)) then
                            this%eigenvectors(:,i) = this%ritz_vectors(:,i)
                            this%eigenvalues(i) = this%ritz_values(i)
                        endif
                    enddo
                endif
                exit
            endif
        ! end Lanczos restart loop
        enddo

        if (check_deltas(this, k)) then
            if (print_info) then
                write(6, '(i2" eigenvalues(s) were successfully converged to within "5ES16.7)') &
                    this%n_states, convergence_error
                call neci_flush(6)
            endif
        endif

        if (iprocindex==root) then
            call perform_orthogonality_test(this)
        endif

        if (print_info) then
            write(6,'(/,1x,"Final calculated energies:")')
            do i=1, this%n_states
                write(6,final_output_fmt) i, this%eigenvalues(i)
                call neci_flush(6)
            enddo
        endif

        ! how good are the ritz vectors?
        write(6,'(/,1x,"Ritz vector expectation energies:")')
        do i=1, this%n_states
            exp_val = get_rayleigh_quotient(this, i)
            if (print_info) then
                write(6,final_output_fmt) i, exp_val
                call neci_flush(6)
            endif
        enddo

        write(6,'(/,1x,"Ritz vector residual norms:")')
        do i=1, this%n_states
            exp_val = compute_residual_norm(this, i)
            if (print_info) then
                write(6,final_output_fmt) i, exp_val
                call neci_flush(6)
            endif
        enddo
        write(6,'(/,1x,"End of Lanczos procedure.",/)')

        ! wait for final diagonalisation before freeing any memory
        call MPIBarrier(ierr)
        call FreeLanczosCalc(this)
        end associate
    end subroutine perform_lanczos

end module lanczos_general
