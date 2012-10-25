#include "macros.h"

module semi_stochastic

    use bit_reps, only: NIfTot, decode_bit_det_bitwise
    use constants
    use DeterminantData, only: write_det
    use FciMCData, only: HFDet, ilutHF
    use GenRandSymExcitNUMod, only: enumerate_all_single_excitations, &
                                    spat_excit_t
    use SystemData, only: nel

    implicit none

contains

    subroutine init_semi_stochastic()

        type(spat_excit_t), target :: gen_store
        integer :: i, nI(nel)
        logical :: first_det
        integer(n_int) :: ilut(0:NIfTot)

        call write_det(6, HFDet, .true.)

        ilut(0) = -1
        call enumerate_all_single_excitations (ilutHF, HFDet, ilut, gen_store)
        do while(ilut(0) /= -1)

            call decode_bit_det_bitwise(nI, ilut)
            call write_det(6, nI, .true.)

            call enumerate_all_single_excitations (ilutHF, HFDet, ilut, &
                                                   gen_store)
        end do

    end subroutine init_semi_stochastic

end module semi_stochastic
