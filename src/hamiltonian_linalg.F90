#include "macros.h"
module hamiltonian_linalg
    !   robert.anderson@kcl.ac.uk
    !   June 2016
    !
    !   (Parallel) linear algebra routines for hermitian matrices
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
    !Use the determinat with lowest energy as HF, otherwise use the value in hfindex
    !This is mainly added so that HF determinant in FCI-Davidson can be specified.
    logical :: tCalcHFIndex = .True.

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
        HElement_t(dp), allocatable, dimension(:) :: partial_H_ket
        ! displacements necessary for communication of H |ket>
        ! (formerly davidson_disps)
        integer(MPIArg), allocatable, dimension(:) :: partial_H_ket_disps
        ! All algorithms for solving large eigenproblems involve a unitary rotation of the
        ! Hamiltonian into a smaller basis. basis_vectors(:,i) is the ith such unit vector
        HElement_t(dp), allocatable, dimension(:,:) :: basis_vectors
        ! the location in the basis is the HF determinant
        integer :: hfindex
        ! the hamiltonian projected into basis_vectors
        real(dp), allocatable, dimension(:,:) :: projected_hamil
        ! we'll usually need some working space for diagonalisation of H in the small basis
        real(dp), allocatable, dimension(:,:) :: projected_hamil_work
        ! For parallel calculations, only the processor with label root performs the main
        HElement_t(dp), allocatable, dimension(:) :: temp_in, temp_out
        logical :: skip_calc
        logical :: t_store_subspace_basis
        ! should we orthogonalise each new basis vector against previous basis vectors?
        logical :: t_orthogonalise
        integer :: max_subspace_size
    end type

    ! leaving these as global data, since direct CI currently only works for Davidson
    type(ras_vector), allocatable, dimension(:,:,:) :: direct_ci_inp, direct_ci_out


    interface inner_product
        module procedure inner_product_real
        module procedure inner_product_complex
    end interface

    interface euclidean_norm
        module procedure euclidean_norm_real
        module procedure euclidean_norm_complex
    end interface

    interface euclidean_norm_square
        module procedure euclidean_norm_square_real
        module procedure euclidean_norm_square_complex
    end interface

    interface multiply_hamil_and_vector
        module procedure multiply_hamil_and_vector_real
        module procedure multiply_hamil_and_vector_complex
    end interface

    interface multiply_hamil_and_vector_full
        module procedure multiply_hamil_and_vector_full_real
        module procedure multiply_hamil_and_vector_full_complex
    end interface

    interface mult_hamil_vector_sparse
        module procedure mult_hamil_vector_sparse_real
        module procedure mult_hamil_vector_sparse_complex
    end interface

    interface mult_hamil_vector_par_sparse
        module procedure mult_hamil_vector_par_sparse_real
        module procedure mult_hamil_vector_par_sparse_complex
    end interface

    interface mult_hamil_vector_direct_ci
        module procedure mult_hamil_vector_direct_ci_real
        module procedure mult_hamil_vector_direct_ci_complex
    end interface

    contains

    subroutine InitHamiltonianCalc(this, print_info, hamil_type, max_subspace_size, t_store_subspace_basis, t_orthogonalise)
        use direct_ci, only: create_ham_diag_direct_ci
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings
        use ras, only: find_ras_size
        use util_mod, only: int_fmt

        type(HamiltonianCalcType), intent(inout) :: this
        logical, intent(in) :: print_info
        integer :: i, mem_reqd, residual_mem_reqd, ierr, space_size
        integer, intent(in) :: hamil_type, max_subspace_size
        logical, intent(in) :: t_store_subspace_basis, t_orthogonalise
        character(*), parameter :: t_r = "InitHamiltonianCalc"
        integer(MPIArg) :: mpi_temp
        real(dp), allocatable :: hamil_diag_temp(:)

        this%hamil_type = hamil_type
        this%t_store_subspace_basis = t_store_subspace_basis
        this%t_orthogonalise = t_orthogonalise

        associate(&
            space_size => this%space_size, &
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
        ! may not exceed size of space)
        this%max_subspace_size = max_subspace_size!min(max_subspace_size, space_size)

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
                this%partial_H_ket_disps(i) = this%partial_H_ket_disps(i-1) &
                                            + this%space_sizes(i-1)
            end do
            call MPIGatherV(hamil_diag, hamil_diag_temp, this%space_sizes, this%partial_H_ket_disps, ierr)

            if (iProcIndex == root) then
                safe_realloc(hamil_diag, (space_size))
                hamil_diag = hamil_diag_temp
            end if
            deallocate(hamil_diag_temp)
        end if


        if (print_info) then
            write(6,'(1x,"number of determinants in total:",'//int_fmt(space_size,1)//')') space_size; call neci_flush(6)
        end if

        if(tCalcHFIndex)then
            hfindex = maxloc((-hamil_diag),1)
        end if

        ! the memory required to allocate each of basis_vectors and
        ! multipied_basis_vectors, in mb.
        mem_reqd = max_subspace_size*this%space_size*8/1000000

        ! allocate the necessary arrays:
        if (t_store_subspace_basis) then
            safe_calloc(this%basis_vectors, (this%space_size, max_subspace_size), 0.0_dp)
            if (print_info) then
                write(6,'(1x,"allocating array to hold subspace vectors (",'//int_fmt(mem_reqd,0)//',1x,"mb).")') mem_reqd
                call neci_flush(6)
            endif
        endif
        safe_calloc(this%projected_hamil, (max_subspace_size,max_subspace_size), 0.0_dp)
        safe_calloc(this%projected_hamil_work, (max_subspace_size,max_subspace_size), 0.0_dp)
        if (iprocindex /= root) then
            safe_malloc(this%temp_in, (this%space_size))
            safe_malloc(this%temp_out, (this%space_size))
        end if

        this%skip_calc = .false.
        if (print_info) then
            write(6,'(1x,"Hamiltonian calculation setup complete.",/)')
        endif
        call neci_flush(6)

        end associate

    end subroutine InitHamiltonianCalc


    subroutine orthogonalise_against_previous_basis_vectors(this, basis_index)
        type(HamiltonianCalcType), intent(inout) :: this
        integer, intent(in) :: basis_index
        integer :: i

        associate(v => this%basis_vectors(:,basis_index))
        do i = 1, basis_index-1
            v = v - this%basis_vectors(:,i)*inner_product(v,this%basis_vectors(:,i)) &
                /inner_product(this%basis_vectors(:,i), this%basis_vectors(:,i))
        enddo
        !v = v/euclidean_norm(v)
        end associate

    end subroutine orthogonalise_against_previous_basis_vectors

    subroutine build_full_hamiltonian_from_sparse(this, full_ham)
        type(HamiltonianCalcType), intent(inout) :: this
        HElement_t(dp), intent(out) :: full_ham(:,:)
        HElement_t(dp), allocatable :: vec(:)
        integer :: i, ierr
        write(6,*) "converting sparse hamiltonian to full"
        safe_malloc(vec, (this%space_size))
        do i = 1, this%space_size
            call MPIBarrier(ierr)
            vec = h_cast(0.0_dp)
            vec(i) = h_cast(1.0_dp)
            if (iprocindex==root) then
                call multiply_hamil_and_vector(this, vec, full_ham(i,1:this%space_size))
            else
                call multiply_hamil_and_vector(this, vec, this%temp_out)
            endif
            write(6,*) i,"hamiltonian columns converted"
            call MPIBarrier(ierr)
        enddo
        safe_free(vec)
    end subroutine build_full_hamiltonian_from_sparse

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

    ! static BLAS wrappers
    function inner_product_real(u, v) result (inner_product)
        real(dp), intent(in) :: u(:), v(:)
        real(dp) :: inner_product, ddot
        character(*), parameter :: this_routine = "inner_product_real"
        ASSERT(size(u)==size(v))
        inner_product = ddot(size(u), u, 1, v, 1)
    end function inner_product_real

    function inner_product_complex(u, v) result (inner_product)
        complex(dp), intent(in) :: u(:), v(:)
        complex(dp) :: inner_product, zdotc
        character(*), parameter :: this_routine = "inner_product_complex"
        ASSERT(size(u)==size(v))
        inner_product = dot_product(u,v)
    end function inner_product_complex

    function euclidean_norm_square_real(u) result (euclidean_norm_square)
        real(dp), intent(in) :: u(:)
        real(dp) :: euclidean_norm_square
        euclidean_norm_square = inner_product(u, u)
    end function euclidean_norm_square_real

    function euclidean_norm_square_complex(u) result (euclidean_norm_square)
        complex(dp), intent(in) :: u(:)
        real(dp) :: euclidean_norm_square
        euclidean_norm_square = real(inner_product(u, u), dp)
    end function euclidean_norm_square_complex

    function euclidean_norm_real(u) result (euclidean_norm)
        real(dp), intent(in) :: u(:)
        real(dp) :: euclidean_norm
        euclidean_norm = sqrt(euclidean_norm_square(u))
    end function euclidean_norm_real

    function euclidean_norm_complex(u) result (euclidean_norm)
        complex(dp), intent(in) :: u(:)
        real(dp) :: euclidean_norm
        euclidean_norm = sqrt(euclidean_norm_square(u))
    end function euclidean_norm_complex

    subroutine pretty_print(out_unit, mat)
        integer :: out_unit, shp(2), i
        real(dp) :: mat(:,:)
        shp = shape(mat)
        do i = 1, shp(2)
            write(out_unit, *) real(mat(:,i))
        enddo
        write(out_unit, *)
    end subroutine

    ! Hamiltonian-ket multiplication
    subroutine multiply_hamil_and_vector_real(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)

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
    end subroutine multiply_hamil_and_vector_real

    subroutine multiply_hamil_and_vector_complex(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        complex(dp), intent(in) :: input_vector(:)
        complex(dp), intent(out) :: output_vector(:)

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
    end subroutine multiply_hamil_and_vector_complex

    subroutine multiply_hamil_and_vector_full_real(this, input_vector, output_vector)
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
        ! 1 is the increment of the elements of y.

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
    end subroutine multiply_hamil_and_vector_full_real

    subroutine multiply_hamil_and_vector_full_complex(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        complex(dp), intent(in) :: input_vector(:)
        complex(dp), intent(out) :: output_vector(:)

        associate(space_size => this%space_size)

        call zgemv('N', &
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
    end subroutine multiply_hamil_and_vector_full_complex

    subroutine mult_hamil_vector_sparse_real(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)
        integer :: i, j

        associate(space_size => this%space_size)
        output_vector = 0.0_dp

        do i = 1, space_size
            do j = 1, sparse_ham(i)%num_elements
                output_vector(i) = output_vector(i) + real(sparse_ham(i)%elements(j)*input_vector(sparse_ham(i)%positions(j)), dp)
            end do
        end do

        end associate
    end subroutine mult_hamil_vector_sparse_real

    subroutine mult_hamil_vector_sparse_complex(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        complex(dp), intent(in) :: input_vector(:)
        complex(dp), intent(out) :: output_vector(:)
        integer :: i, j

        associate(space_size => this%space_size)
        output_vector = 0.0_dp

        do i = 1, space_size
            do j = 1, sparse_ham(i)%num_elements
                output_vector(i) = output_vector(i) + sparse_ham(i)%elements(j)*input_vector(sparse_ham(i)%positions(j))
            end do
        end do

        end associate
    end subroutine mult_hamil_vector_sparse_complex

    subroutine mult_hamil_vector_par_sparse_real(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)
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
        ! this templated routine forbids mixing real and complex arrays
        call MPIGatherV(real(this%partial_H_ket, dp), output_vector, this%space_sizes, this%partial_H_ket_disps, ierr)
    end subroutine mult_hamil_vector_par_sparse_real

    subroutine mult_hamil_vector_par_sparse_complex(this, input_vector, output_vector)
        type(HamiltonianCalcType), intent(inout) :: this
        complex(dp), intent(in) :: input_vector(:)
        complex(dp), intent(out) :: output_vector(:)
        integer :: i, j, ierr

        ! Use output_vector as temporary space.
        output_vector = input_vector

        call MPIBarrier(ierr, tTimeIn=.false.)

        call MPIBCast(output_vector)

        this%partial_H_ket = 0.0_dp

        do i = 1, this%space_sizes(iProcIndex)
            do j = 1, sparse_ham(i)%num_elements
#ifdef __CMPLX
                this%partial_H_ket(i) = this%partial_H_ket(i) + &
                    sparse_ham(i)%elements(j)*output_vector(sparse_ham(i)%positions(j))
#endif
            end do
        end do
#ifdef __CMPLX
        call MPIGatherV(this%partial_H_ket, output_vector, this%space_sizes, this%partial_H_ket_disps, ierr)
#else
        call MPIGatherV(cmplx(this%partial_H_ket,0.0_dp,dp), output_vector, this%space_sizes, this%partial_H_ket_disps, ierr)
#endif
        !call MPIGatherV(this%partial_H_ket, output_vector, this%space_sizes, this%partial_H_ket_disps, ierr)
    end subroutine mult_hamil_vector_par_sparse_complex

    subroutine mult_hamil_vector_direct_ci_real(input_vector, output_vector)
        use direct_ci, only: perform_multiplication, transfer_from_block_form, transfer_to_block_form
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits
        use SystemData, only: ecore

        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)

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
    end subroutine mult_hamil_vector_direct_ci_real


    subroutine mult_hamil_vector_direct_ci_complex(input_vector, output_vector)
        use direct_ci, only: perform_multiplication, transfer_from_block_form, transfer_to_block_form
        use FciMCData, only: davidson_ras, davidson_classes, davidson_strings, davidson_iluts, davidson_excits
        use SystemData, only: ecore

        complex(dp), intent(in) :: input_vector(:)
        complex(dp), intent(out) :: output_vector(:)
        character(*), parameter :: t_r = "mult_hamil_vector_direct_ci_complex"

        unused_var(input_vector)
        output_vector = 0.0d0
        call stop_all(t_r, "not yet implemented for complex CI coefficients")

    end subroutine mult_hamil_vector_direct_ci_complex

end module hamiltonian_linalg
