#include "macros.h"

module semi_stochastic

    use bit_reps, only: NIfTot, decode_bit_det_bitwise
    use constants
    use DeterminantData, only: write_det
    use FciMCData, only: HFDet, ilutHF
    use GenRandSymExcitNUMod, only: enumerate_all_single_excitations, &
                                    enumerate_all_double_excitations, spat_excit_t
    use SystemData, only: nel

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

        call write_det(6, HFDet, .true.)

        ! This condition tells the enumerating subroutines to initialise the loops.
        ilut(0) = -1
        ! Find the first determinant.
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        ! When no more basis functions are found, this value if returned and the loop
        ! is exited.
        do while(ilut(0) /= -1)

            call decode_bit_det_bitwise(nI, ilut)
            call write_det(6, nI, .true.)

            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, &
                                                   gen_store)
        end do

        write(6,*) ""

        ilut(0) = -1
        call enumerate_all_double_excitations (ilutHF, HFDet, ilut, gen_store)
        do while(ilut(0) /= -1)

            call decode_bit_det_bitwise(nI, ilut)
            call write_det(6, nI, .true.)

            call enumerate_all_double_excitations (ilutHF, HFDet, ilut, &
                                                   gen_store)
        end do

    end subroutine init_semi_stochastic

end module semi_stochastic
