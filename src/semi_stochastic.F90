#include "macros.h"

module semi_stochastic

    use bit_reps, only: NIfTot, decode_bit_det_bitwise
    use constants
    use FciMCData, only: HFDet, ilutHF, iHFProc
    use GenRandSymExcitNUMod, only: enumerate_all_single_excitations, &
                                    enumerate_all_double_excitations, spat_excit_t
    use hash , only : DetermineDetNode
    use Parallel_neci, only: iProcIndex
    use SystemData, only: nel, det_vector, result_det_vector

    implicit none

contains

    subroutine init_semi_stochastic()

        ! Initialise the semi-stochastic information. This includes enumerating a list
        ! of all determinants or CSFs in the deterministic space and calculating and
        ! storing the resulting Hamiltonian matrix elements. The lists which will store
        ! the psip vectors in the deterministic space are also allocated.

        type(spat_excit_t), target :: gen_store
        integer :: i, nI(nel)
        integer(n_int) :: ilut(0:NIfTot)
        integer :: det_vector_size, result_det_vector_size

        result_det_vector_size = 0
        ! Initialise as 1 to count the Hartree-Fock determinant.
        det_vector_size = 1

        ! If the Hartree-Fock determinant is on this processor.
        if (iProcIndex == iHFProc) result_det_vector_size = 1 

        ! This condition tells the enumerating subroutines to initialise the loop.
        ilut(0) = -1

        ! Find the first single excitation.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        ! Subroutine to add this state to the spawned list if on this processor.
        call add_basis_fn_to_list(ilut, det_vector_size, result_det_vector_size, nI)

        ! When no more basis functions are found, this value if returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)

            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, &
                                                   gen_store)

            ! If a determinant is returned (if we did not find the last one last time.)
            if (ilut(0) /= -1) call add_basis_fn_to_list(ilut, det_vector_size, &
                                                         result_det_vector_size, nI)
        end do

        ! Now generate double excitations.

        ilut(0) = -1
        ! The first double excitation.
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        call add_basis_fn_to_list(ilut, det_vector_size, result_det_vector_size, nI)

        do while(ilut(0) /= -1)

            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, &
                                                   gen_store)

            if (ilut(0) /= -1) call add_basis_fn_to_list(ilut, det_vector_size, &
                                                         result_det_vector_size, nI)
        end do

        ! Now that we know the total size of det_vector and result_det_vector, allocate
        ! them.
        allocate(det_vector(det_vector_size))
        allocate(result_det_vector(result_det_vector_size))

    end subroutine init_semi_stochastic

    subroutine add_basis_fn_to_list(ilut, vector_size, result_vector_size, nI)

        ! This subroutine, called from init_semi_stochastic, takes a bitstring, finds
        ! the nI representation, decides if the basis state lives on this processor
        ! and, if so, adds it to the spawned list. It also increases the two vector
        ! sizes being calculated.

        ! In: ilut - The determinant in a bitstring form.
        ! Inout : vector_size - This will store the size of the deterministic space.
        ! Inout : result_vector_size - The size of the vector storing the deterministic
        !          amplitudes on this processor after multiplication by core_hamilonian.
        ! Out: nI - The list of occupied orbitals in an array.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(inout) :: vector_size, result_vector_size
        integer, intent(out) :: nI(nel)
        integer :: proc

        ! Find the nI representation of determinant.
        call decode_bit_det_bitwise(nI, ilut)
        ! Find the processor which this state belongs to.
        proc = DetermineDetNode(nI,0)
        ! If this determinant belongs to this processor, increase the vector size and
        ! add it to the spawned list.        
        vector_size = vector_size + 1
        if (proc == iProcIndex) result_vector_size = result_vector_size + 1

    end subroutine add_basis_fn_to_list

end module semi_stochastic
