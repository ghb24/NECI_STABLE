#include "macros.h"

module semi_stochastic

    use bit_rep_data, only: flag_deterministic
    use bit_reps, only: NIfD, NIfTot, decode_bit_det_bitwise, encode_bit_rep
    use CalcData, only: tRegenDiagHEls
    use constants
    use detBitOps, only: ilut_lt, ilut_gt
    use Determinants, only: get_helement
    use FciMCData, only: HFDet, ilutHF, iHFProc, Hii, CurrentDets, CurrentH, &
                         deterministic_proc_sizes, deterministic_proc_indices, &
                         det_vector, result_det_vector, core_hamiltonian
    use GenRandSymExcitNUMod, only: enumerate_all_single_excitations, &
                                    enumerate_all_double_excitations, spat_excit_t
    use hash , only : DetermineDetNode
    use hphf_integrals, only: hphf_diag_helement
    use Parallel_neci, only: iProcIndex, nProcessors, MPIBCast, MPIBarrier
    use sort_mod, only: sort
    use SystemData, only: nel, tHPHF

    implicit none

contains

    subroutine init_semi_stochastic()

        ! Initialise the semi-stochastic information. This includes enumerating a list
        ! of all determinants or CSFs in the deterministic space and calculating and
        ! storing the resulting Hamiltonian matrix elements. The lists which will store
        ! the psip vectors in the deterministic space are also allocated.

        type(spat_excit_t), target :: gen_store
        integer :: i, j, col_index, iproc, IC
        integer :: nI(nel), nJ(nel)
        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer :: det_space_size, ierr
        integer :: no_zero, no_not_zero
        real :: fraction_zero

        no_zero = 0
        no_not_zero = 0

        allocate(deterministic_proc_sizes(0:nProcessors-1))
        allocate(deterministic_proc_indices(0:nProcessors-1))
        deterministic_proc_sizes = 0
        deterministic_proc_indices = 0

        ! Count the Hartree-Fock determinant (The Hartree-Fock determinant has already been
        ! added to the main list as usual in the subroutine InitFCIMCCalcPar).
        deterministic_proc_sizes(iHFProc) = 1

        ! This condition tells the enumerating subroutines to initialise the loop.
        ilut(0) = -1

        ! Find the first single excitation.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        ! Subroutine to add this state to the spawned list if on this processor.
        call add_basis_state_to_list(ilut, nI)

        ! When no more basis functions are found, this value if returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)
            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)

            ! If a determinant is returned (if we did not find the final one last time.)
            if (ilut(0) /= -1) call add_basis_state_to_list(ilut, nI)
        end do

        ! Now generate the double excitations...

        ilut(0) = -1
        ! The first double excitation.
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_basis_state_to_list(ilut, nI)

        do while(ilut(0) /= -1)
            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)

            if (ilut(0) /= -1) call add_basis_state_to_list(ilut, nI)
        end do

        ! All excitations have now been generated.

        ! We now know the size of the deterministic space on each processor, so now find
        ! the total size of the space and also allocate vectors to store psip amplitudes
        ! and the deterministic Hamiltonian.
        det_space_size = sum(deterministic_proc_sizes)
        allocate(det_vector(deterministic_proc_sizes(iProcIndex)))
        allocate(result_det_vector(det_space_size))
        allocate(core_hamiltonian(deterministic_proc_sizes(iProcIndex), det_space_size))
        det_vector = 0.0_dp
        result_det_vector = 0.0_dp

        ! Calculate the indices in the full vector at which the various processors take
        ! over, relative to the first index position in the vector (ie, the first value
        ! in this vector will be 0, *not* 1 - this is required for mpiallgatherv later).
        do i = 0, nProcessors-1
            if (i == 0) then
                deterministic_proc_indices(i) = 0
            else
                deterministic_proc_indices(i) = deterministic_proc_indices(i-1) + &
                                                  deterministic_proc_sizes(i-1)
            end if
        end do

        !write(6,*) deterministic_proc_indices

        ! Sort the list of basis states so that it is correctly ordered in the predefined
        ! order which is always kept throughout the simulation.
        call sort(CurrentDets(:,1:deterministic_proc_sizes(iProcIndex)), ilut_lt, ilut_gt)

        ! temp_store is storage space for bitstrings so that the Hamiltonian matrix
        ! elements can be calculated.
        allocate(temp_store(0:NIfTot, maxval(deterministic_proc_sizes)))

        !write (6,*) iHFProc
        !write (6,*) "size = ", deterministic_proc_sizes(iProcIndex)
        !write (6,*) "proc = ", iProcIndex
        !write (6,*) "no procs = ", nProcessors

        ! Calcuation of the Hamiltonian matrix elements to be stored. Loop over each
        ! processor in order and broadcast the sorted vector of bitstrings from the
        ! processor to all other processors so that all matrix elements can be found.
        do iproc = 0, nProcessors-1

            ! If we are broadcasting from this processor next, transfer the bitstrings
            ! to the array temp_store.
            if (iproc == iProcIndex) temp_store(:,1:deterministic_proc_sizes(iproc)) = &
                                      CurrentDets(:,1:deterministic_proc_sizes(iproc))

            ! Perform the broadcasting to other all other processors.
            call MPIBCast(temp_store, size(temp_store), iproc)

            !write (6,*) "iproc = ", iproc, "temp_store = ", temp_store(2,:)
            !write (6,*) ""

            ! The index of the column before the first column of the block of the
            ! Hamiltonian currently being calculated.
            col_index = deterministic_proc_indices(iproc)

            call MPIBarrier(ierr)

            ! Loop over all the elements in the block of the Hamiltonian corresponding
            ! to these two prcoessors.
            do i = 1, deterministic_proc_sizes(iProcIndex)
                do j = 1, deterministic_proc_sizes(iproc)
                    !write (6,*) iprocIndex, iproc, i, j
                    !call neci_flush(6)

                    ! Get nI form of the basis functions.
                    call decode_bit_det_bitwise(nI, CurrentDets(:, i))
                    call decode_bit_det_bitwise(nJ, temp_store(:, j))

                    ! If on the diagonal of the Hamiltonian.
                    if ((iprocIndex == iproc) .and. (i == j)) then
                        core_hamiltonian(i, col_index + j) = &
                                             get_helement(nI, nJ, 0) - Hii
                        ! We calculate and store CurrentH at this point for ease.
                        if (.not.tRegenDiagHEls) CurrentH(1,i) = &
                                             core_hamiltonian(i, col_index + j)
                    else
                        core_hamiltonian(i, col_index + j) = &
                                         get_helement(nI, nJ, CurrentDets(:, i), &
                                                            temp_store(:, j))
                    end if

                    !if (abs(core_hamiltonian(i, col_index + j)) > 0.0_dp) then
                    !    no_not_zero = no_not_zero + 1
                    !else
                    !    no_zero = no_zero + 1
                    !end if

                end do
            end do

        end do

        !fraction_zero = real(no_zero)/real(no_zero + no_not_zero)

        call MPIBarrier(ierr)

        !write (6,*) "Fraction zero = ", fraction_zero
        !call neci_flush(6)

        deallocate(temp_store)

        !do i = 1, deterministic_proc_sizes(iProcIndex)
        !    write(6,*) core_hamiltonian(i,:)
        !    write(6,*) ""
        !end do

        !call stop_all("here","here")

    end subroutine init_semi_stochastic

    subroutine add_basis_state_to_list(ilut, nI)

        ! This subroutine, called from init_semi_stochastic, takes a bitstring, finds
        ! the nI representation, decides if the basis state lives on this processor
        ! and, if so, adds it to the spawned list. It also increases the vector sizes
        ! being calculated.

        ! In: ilut - The determinant in a bitstring form.
        ! Out: nI - The list of occupied orbitals in an array.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(out) :: nI(nel)
        integer(n_int) :: flags
        integer :: proc, sgn(lenof_sign)

        sgn = 0
        ! Flag to specify that these basis states are in the deterministic space.
        flags = 0
        flags = ibset(flags,flag_deterministic)

        ! Find the nI representation of determinant.
        call decode_bit_det_bitwise(nI, ilut)
        ! Find the processor which this state belongs to.
        proc = DetermineDetNode(nI,0)
        ! Increase the size of the deterministic space on the correct processor.
        deterministic_proc_sizes(proc) = deterministic_proc_sizes(proc) + 1
        ! If this determinant belongs to this processor, add it to the main list.
        if (proc == iProcIndex) call encode_bit_rep(CurrentDets(:, &
                                 deterministic_proc_sizes(iProcIndex)), ilut, sgn, flags)

    end subroutine add_basis_state_to_list

end module semi_stochastic
