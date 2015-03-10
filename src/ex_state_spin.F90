#include "macros.h"

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use kp_fciqmc_data_mod, only: lenof_all_signs

module ex_state_spin

    subroutine calc_projected_spin(nvecs, krylov_array, krylov_ht, ndets, spin_matrix)

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        type(ll_node), pointer, intent(in) :: krylov_ht(:)
        integer, intent(in) :: ndets
        real(dp), intent(out) :: spin_matrix(:,:)

        integer :: idet, iel, jel, iorb, jorb, i, j
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        integer :: nI_parent(nel)
        logical :: ialpha, jalpha, ibeta, jbeta
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: real_sign(lenof_all_signs)

        ilut_parent = 0_n_int

        ! Loop over all states and add all contributions.
        do idet = 1, ndets

            ! The 'parent' determinant to consider 'spawning' from.
            ilut_parent(0:NIfDBO) = krylov_array(0:NIfDBO,idet)

            call decode_bit_det(nI_parent, ilut_parent)

            do iel = 1, nel
                iorb = nI(i)
                ibeta = is_beta(iorb)
                ialpha = .not. ibeta

                do jel = 1, nel
                    jorb = nI(j)
                    jbeta = is_beta(jorb)
                    jalpha = .not. jbeta

                    ! If both are beta orbitals, or both are alpha orbitals,
                    ! then add the contributions to spin matrix.
                    if ((ibeta .and. jbeta) .or. (ialpha .and. jalpha)) then
                        int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, idet)
                        real_sign = transfer(int_sign, real_sign)

                        do i = 1, nvecs
                            do j = 1, nvecs
                                spin_matrix(i,j) = 0.25_dp*(real_sign(2*i-1)*real_sign(2*j) + &
                                                            real_sign(2*i)*real_sign(2*j-1))/2.0_dp
                            end do
                        end do
                    end if
                end do
            end do
        end do

    end subroutine calc_projected_spin

end module ex_state_spin
