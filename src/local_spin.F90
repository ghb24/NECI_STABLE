#include "macros.h"

module local_spin

      use constants, only: dp, n_int, lenof_sign
      use bit_rep_data, only: IlutBits
      use LoggingData, only: num_local_spin_orbs


      implicit none

      real(dp) :: inst_local_spin = 0.0_dp, all_local_spin = 0.0_dp
      real(dp) :: sum_local_spin = 0.0_dp, all_sum_local_spin = 0.0_dp

      ! it seems there is no record of the averaged norm hm..
      real(dp) :: sum_norm_psi_squared = 0.0_dp

contains

    subroutine measure_local_spin(ilut, nI, real_sgn)
        debug_function_name("measure_local_spin")
        integer(n_int), intent(in) :: ilut(0:IlutBits%len_tot)
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: real_sgn(lenof_sign)

        real(dp) :: coeff

#if defined PROG_NUMRUNS_ || defined DOUBLERUN_
#ifdef CMPLX_
        ! i do not want to deal with complex runs for now..
        call stop_all(this_routine, &
                      "complex double occupancy measurement not yet implemented!")
#else
        coeff = real_sgn(1) * real_sgn(2)
#endif
#else
        coeff = abs(real_sgn(1))**2
#endif
        ! the current b vector should be fine to get the total spin
        loc_spin = 2*currentB_ilut(num_local_spin_orbs) * &
            (2*currentB_ilut(num_local_spin_orbs) + 1.0_dp)

        inst_local_spin = inst_local_spin + coeff * loc_spin

    end subroutine measure_local_spin

    subroutine finalize_local_spin_measurement()
        debug_function_name("finalize_local_spin_measurement")

        call MPIAllreduce(local_spin, MPI_SUM, all_local_spin)

        if (iProcIndex == root) then
            print *, "local spin up to orbital: ", num_local_spin_orbs

        end if

    end subroutine finalize_local_spin_measurement

    subroutine rezero_local_spin_stats()
        inst_local_spin = 0.0_dp
    end subroutine rezero_local_spin_stats
end module local_spin
