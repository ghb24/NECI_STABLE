!   robert.anderson@kcl.ac.uk
!   June 2016
!
!   Linear algebra routines for FCI hamiltonians
!
module hamiltonian_linalg
    use FciMCData, only : hamiltonian
    use sparse_arrays, only: sparse_ham, hamil_diag

    implicit none
    ! The value of hamil_type specifies what form the Hamiltonian is stored in.
    ! The following options are currently available:
    integer, parameter :: full_hamil_type = 1
    integer, parameter :: sparse_hamil_type = 2
    integer, parameter :: parallel_sparse_hamil_type = 3
    integer, parameter :: direct_ci_type = 4

    type HamiltonianCalcType
        integer :: hamil_type
        ! The dimension of the vector space we are working in, as determined by the number
        ! of rows and columns in the Hamiltonian matrix.
        integer :: space_size 
        ! For parallel calculations, store all spaces sizes on each processor
        integer(MPIArg), allocatable, dimension(:) :: space_sizes
        ! For parallel calculations, this vector is the size of the space on this processor. This
        ! vector is used to store the output of multiplication by the Hamiltonian on this processor.
        ! (formerly partial_davidson_vector)
        real(dp), allocatable, dimension(:) :: partial_working_vector
        logical :: skip_calc
    end type

    contains

    function initHamiltonianCalc(hamil_type, space_size, skip_calc) result (hamiltonianCalc)
        integer, intent(in) :: hamil_type, space_size
        logical, intent(in) :: skip_calc
        type(HamiltonianCalcType) :: hamiltonianCalc
        hamiltonianCalc%hamil_type = hamil_type
        hamiltonianCalc%space_size = space_size
        hamiltonianCalc%skip_calc = skip_calc
    end function

    subroutine multiply_hamil_and_vector(thisCalc, input_vector, output_vector)

        type(HamiltonianCalcType), intent(in) :: thisCalc
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)

        associate(thisCalc%hamil_type => hamil_type)

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

    subroutine multiply_hamil_and_vector_full(thisCalc, input_vector, output_vector)

        type(HamiltonianCalcType), intent(in) :: thisCalc
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
        ! 1 is the incremenet of the elements of y.

        associate(space_size => thisCalc%space_size)

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

    subroutine mult_hamil_vector_sparse(thisCalc, input_vector, output_vector)

        type(HamiltonianCalcType), intent(in) :: thisCalc
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)
        integer :: i, j

        associate(space_size => thisCalc%space_size)
        output_vector = 0.0_dp

        do i = 1, space_size
            do j = 1, sparse_ham(i)%num_elements
                output_vector(i) = output_vector(i) + sparse_ham(i)%elements(j)*input_vector(sparse_ham(i)%positions(j))
            end do
        end do

    end subroutine mult_hamil_vector_sparse

    subroutine mult_hamil_vector_par_sparse(thisCalc, input_vector, output_vector)

        type(HamiltonianCalcType), intent(in) :: thisCalc
        real(dp), intent(in) :: input_vector(:)
        real(dp), intent(out) :: output_vector(:)
        integer :: i, j, ierr

        ! Use output_vector as temporary space.
        associate(space_sizes => thisCalc%space_sizes)
        output_vector = input_vector

        call MPIBarrier(ierr, tTimeIn=.false.)

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

end module hamiltonian_linalg
