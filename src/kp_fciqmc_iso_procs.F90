#include "macros.h"

module kp_fciqmc_iso_procs

    use bit_rep_data
    use constants
    use kp_fciqmc_procs

    use FciMCData, only: TotWalkers, CurrentDets

contains

    subroutine calc_overlap_matrix(kp)

        type(kp_fciqmc_data), intent(inout) :: kp
        integer :: idet, jvec, ind(kp%ivec)
        integer(n_int) :: int_sign(lenof_sign)
        real(dp) :: real_sign_1(lenof_sign), real_sign_2(lenof_sign)

        associate(s_matrix => kp%overlap_matrix)

            ! Just in case!
            s_matrix = 0.0_dp

            do jvec = 1, kp%nvecs
                ! The first index of the sign in CurrentDets, for each vector.
                ind(jvec) = NIfDBO + lenof_sign*(jvec-1) + 1
            end do

            ! Loop over all determinants in CurrentDets.
            do idet = 1, TotWalkers
                do ivec = 1, kp%nvecs
                    int_sign = CurrentDets(ind(ivec):ind(ivec)+1, idet)
                    real_sign_1 = transfer(int_sign, real_sign_1)
                    if (IsUnoccDet(real_sign_1)) cycle
                    ! Loop over all Krylov vectors currently stored.
                    do jvec = 1, ivec
                        int_sign = CurrentDets(ind(jvec):ind(jvec)+1, idet)
                        real_sign_2 = transfer(int_sign, real_sign_1)
                        if (IsUnoccDet(real_sign_2)) cycle

                        s_matrix(jvec,ivec) = s_matrix(jvec,ivec) + &
                            (real_sign_1(1)*real_sign_2(2) + real_sign_1(2)*real_sign_2(1))/2.0_dp
                    end do
                end do
            end do

            ! Fill in the lower-half of the overlap matrix.
            do jvec = 1, ivec
                s_matrix(ivec,jvec) = s_matrix(jvec,ivec)
            end do

        end associate

    end subroutine calc_overlap_matrix

end module kp_fciqmc_iso_procs
