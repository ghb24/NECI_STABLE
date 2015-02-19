#include "macros.h"
module orthogonalise

    use FciMCData, only: TotWalkers, CurrentDets, all_norm_psi_squared
    use bit_reps, only: extract_sign, encode_sign
    use Parallel_neci
    use constants
    implicit none

contains

    subroutine orthogonalise_replicas ()

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.
        !
        ! |psi_2'> = |psi_2> - (|psi_1><psi_1|psi_2>)/(<psi_1|psi_1>)

        integer :: j
        real(dp) :: sgn(lenof_sign), delta, scal_prod
        real(dp) :: psi_squared(lenof_sign), all_psi_squared(lenof_sign)
        character(*), parameter :: this_routine = 'orthogonalise_replicas'

        ASSERT(inum_runs == 2)
        ASSERT(lenof_sign == 2)

        ! We need the norm of the wavefunction to do anything. Don't trust
        ! the global values here, as they aren't valid until after the
        ! communication --> not set yet.
        !
        ! We want them to be valid for the new psi...
        psi_squared = 0
        scal_prod = 0
        do j = 1, int(TotWalkers, sizeof_int)
            call extract_sign(CurrentDets(:, j), sgn)
            if (IsUnoccDet(sgn)) cycle
            psi_squared = psi_squared + sgn**2
            scal_prod = scal_prod + sgn(1) * sgn(2)
        end do
        call MPISumAll(psi_squared, all_psi_squared)


#ifndef __PROG_NUMRUNS
        call stop_all(this_routine, "orthogonalise replicas requires mneci.x")
#else
        ! Calculate the change
        do j = 1, int(TotWalkers, sizeof_int)
            
            ! Adjust the wavefunctions
            call extract_sign(CurrentDets(:,j), sgn)
            delta = - sgn(1) * sgn(2) / psi_squared(1)
            sgn(2) = sgn(2) + delta

            ! And stochastically round, so that the minimum particle sign
            ! is maintained in an unbiased way.
            ! 
            ! ... TODO: Don't worry about this for now...

            call encode_sign(CurrentDets(:,j), sgn)

        end do
#endif

    end subroutine

end module
