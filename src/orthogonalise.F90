#include "macros.h"
module orthogonalise

    use FciMCData, only: TotWalkers, CurrentDets, all_norm_psi_squared, &
                         NoBorn, NoDied, fcimc_iter_data
    use bit_reps, only: extract_sign, encode_sign
    use dSFMT_interface, only: genrand_real2_dSFMT
    use Parallel_neci
    use constants
    implicit none

contains

    subroutine orthogonalise_replicas (iter_data)

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.
        !
        ! |psi_2'> = |psi_2> - (|psi_1><psi_1|psi_2>)/(<psi_1|psi_1>)

        integer :: j
        real(dp) :: sgn(lenof_sign), delta, scal_prod, all_scal_prod, r
        real(dp) :: sgn_orig
        real(dp) :: psi_squared(lenof_sign), all_psi_squared(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'orthogonalise_replicas'

        ASSERT(inum_runs == 2)
        ASSERT(lenof_sign == 2)

#ifndef __PROG_NUMRUNS
        call stop_all(this_routine, "orthogonalise replicas requires mneci.x")
#else

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
        call MPISumAll(scal_prod, all_scal_prod)

        ! Calculate the change
        do j = 1, int(TotWalkers, sizeof_int)
            
            ! Adjust the wavefunctions
            call extract_sign(CurrentDets(:,j), sgn)
            delta = - sgn(1) * all_scal_prod / all_psi_squared(1)
            sgn_orig = sgn(2)
            sgn(2) = sgn(2) + delta

            ! And stochastically round, so that the minimum particle sign
            ! is maintained in an unbiased way.
            if (abs(sgn(2)) < 1.0_dp) then
                r = genrand_real2_dSFMT()
                if (r > abs(sgn(2))) then
                    sgn(2) = sign(1.0_dp, sgn(2))
                else
                    sgn(2) = 0.0_dp
                endif
            end if

            call encode_sign(CurrentDets(:,j), sgn)

            ! Now we need to do our accounting to make sure that the NoBorn/
            ! Died/etc. counters remain reasonable.
            !
            ! n.b. we don't worry about the delta=0 case, as adding 0 to the
            !      accumulators doesn't cause errors.
            if (sgn(2) >= 0 .eqv. sgn_orig >= 0) then
                NoDied(2) = NoDied(2) + abs(sgn_orig)
                iter_data%ndied(2) = iter_data%ndied(2) + abs(sgn_orig)
                NoBorn(2) = NoBorn(2) + abs(sgn(2))
                iter_data%nborn(2) = iter_data%nborn(2) + abs(sgn(2))
            else if (abs(sgn(2)) >= abs(sgn_orig)) then
                NoDied(2) = NoBorn(2) + abs(sgn(2) - sgn_orig)
                iter_data%ndied(2) = iter_data%ndied(2) + abs(sgn(2) -sgn_orig)
            else
                NoBorn(2) = NoBorn(2) + abs(sgn(2) - sgn_orig)
                iter_data%nborn(2) = iter_data%nborn(2) + abs(sgn(2) -sgn_orig)
            end if


        end do

        ! Note that we shouldn't be able to kill all particles on a site, as
        ! we can only change run 2 if run 1 is occupied, and run 1 doesn't
        ! change...
        ASSERT(.not. IsUnoccDet(sgn))
#endif

    end subroutine

end module
