#include "macros.h"

module semi_stochastic

    use bit_reps, only: NIfTot, decode_bit_det_bitwise
    use CalcData, only: tRegenDiagHEls
    use constants
    use detBitOps, only: ilut_lt, ilut_gt
    use Determinants, only: get_helement
    use FciMCData, only: HFDet, ilutHF, iHFProc, Hii, CurrentDets
    use GenRandSymExcitNUMod, only: enumerate_all_single_excitations, &
                                    enumerate_all_double_excitations, spat_excit_t
    use hash , only : DetermineDetNode
    use hphf_integrals, only: hphf_diag_helement
    use Parallel_neci, only: iProcIndex
    use SystemData, only: nel, det_vector, result_det_vector, tHPHF

    implicit none

contains

    subroutine init_semi_stochastic()

        ! Initialise the semi-stochastic information. This includes enumerating a list
        ! of all determinants or CSFs in the deterministic space and calculating and
        ! storing the resulting Hamiltonian matrix elements. The lists which will store
        ! the psip vectors in the deterministic space are also allocated.

        type(spat_excit_t), target :: gen_store
        integer :: i, nI(nel)
        integer(n_int) :: ilut(0:NIfTot), flags
        integer :: det_vector_size, result_det_vector_size
        Element_t :: HDiag

        result_det_vector_size = 0
        ! Initialise as 1 to count the Hartree-Fock determinant (Note the Hartree-Fock
        ! determinant has already been added to the main list as usual in the subroutine
        ! InitFCIMCCalcPar).
        use constants , on
        det_vector_size = 1

        ! Create the flag to specify that the basis function is in the deterministic 
        ! space.
        flags = 0
        flags = ibset(flags,flag_deterministic)

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

        ! Sort the list of basis states so that it is correctly ordered in an predefined
        ! order which is always kept throughout the simulation.
        call sort(CurrentDets(:,1:result_vector_size), ilut_lt, ilut_gt)

        ! Now that we know the total size of det_vector and result_det_vector, allocate
        ! them.
        allocate(det_vector(det_vector_size))
        allocate(result_det_vector(result_det_vector_size))

        det_vector = 0.0_dp
        result_det_vector = 0.0_dp

        ! Calculate the Hamiltonian diagonal matrix elements if we are storing them.
        if(.not.tRegenDiagHEls) then
            do i = 1, result_vector_size
                if (tHPHF) then
                    HDiag = hphf_diag_helement (nI,ilut)
                else
                    HDiag = get_helement (nI, nI, 0)
                endif
                CurrentH(1,i)=real(HDiag,dp)-Hii
            end do
        end if

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
        if (proc == iProcIndex) then
            result_vector_size = result_vector_size + 1
            call encode_bit_rep(CurrentDets(:, result_vector_size), iLutJ, 0, flags)
        end if

    end subroutine add_basis_fn_to_list

end module semi_stochastic
