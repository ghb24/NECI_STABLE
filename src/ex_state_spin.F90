#include "macros.h"

module ex_state_spin

    use bit_rep_data, only: NIfTot, NIfDBO, NOffSgn
    use constants
    use kp_fciqmc_data_mod, only: lenof_all_signs

    implicit none

contains

    subroutine calc_projected_spin(nvecs, krylov_array, krylov_ht, ndets, spin_matrix)

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: EncodeBitDet
        use FciMCData, only: ll_node
        use get_excit, only: make_double
        use hash, only: FindWalkerHash
        use SystemData, only: nel

        integer, intent(in) :: nvecs
        integer(n_int), intent(in) :: krylov_array(0:,:)
        type(ll_node), pointer, intent(in) :: krylov_ht(:)
        integer, intent(in) :: ndets
        real(dp), intent(out) :: spin_matrix(:,:)

        integer :: idet, iel, jel, iorb, jorb, i, j
        integer(n_int) :: ilut_child(0:NIfTot), ilut_parent(0:NIfTot)
        integer :: nI_parent(nel), nI_child(nel)
        logical :: ialpha, jalpha, ibeta, jbeta
        integer(n_int) :: int_sign(lenof_all_signs)
        real(dp) :: parent_sign(lenof_all_signs), child_sign(lenof_all_signs)
        real(dp) :: parity_fac
        integer :: ex(2,2), hash_val, det_ind
        logical :: tParity, tDetFound
        type(ll_node), pointer :: temp_node

        ilut_parent = 0_n_int

        ! Loop over all states and add all contributions.
        do idet = 1, ndets

            ! The 'parent' determinant to consider 'spawning' from.
            ilut_parent(0:NIfDBO) = krylov_array(0:NIfDBO,idet)

            ! Get the real walker amplitudes.
            int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, idet)
            parent_sign = transfer(int_sign, parent_sign)

            call decode_bit_det(nI_parent, ilut_parent)

            do iel = 1, nel
                iorb = nI_parent(iel)
                ibeta = is_beta(iorb)
                ialpha = .not. ibeta

                do jel = iel+1, nel

                    jorb = nI_parent(jel)
                    jbeta = is_beta(jorb)
                    jalpha = .not. jbeta

                    ! If both are beta orbitals, or both are alpha orbitals,
                    ! then add the contributions to spin matrix.
                    if ((ibeta .and. jbeta) .or. (ialpha .and. jalpha)) then
                        do i = 1, nvecs
                            do j = i, nvecs
                                spin_matrix(i,j) = spin_matrix(i,j) + &
                                    0.5_dp*(parent_sign(2*i-1)*parent_sign(2*j) + &
                                            parent_sign(2*i)*parent_sign(2*j-1))/2.0_dp
                            end do
                        end do
                    else
                        ! One orbital has alpha spin, the other has beta spin.

                        ! First add the elements corresponding to 'spawning'
                        ! from and to the same determinant.
                        do i = 1, nvecs
                            do j = i, nvecs
                                spin_matrix(i,j) = spin_matrix(i,j) - &
                                    0.5_dp*(parent_sign(2*i-1)*parent_sign(2*j) + &
                                            parent_sign(2*i)*parent_sign(2*j-1))/2.0_dp
                            end do
                        end do

                        ! Now for the more tricky bit - 'spawning' to
                        ! different determinants.
                        
                        ! First check if the two electrons are in the same
                        ! spatial orbital. If they are then this will actually
                        ! 'spawn' to the same determinant, so just add the
                        ! contribution in.
                        if (is_in_pair(iorb, jorb)) then
                            ! Only take each contribution once.
                            do i = 1, nvecs
                                do j = i, nvecs
                                    spin_matrix(i,j) = spin_matrix(i,j) - &
                                        (parent_sign(2*i-1)*parent_sign(2*j) + &
                                         parent_sign(2*i)*parent_sign(2*j-1))/2.0_dp
                                end do
                            end do
                        else
                            ! If one of the new orbitals is already occupied
                            ! then there will be no contribution.
                            if (IsOcc(ilut_parent, ab_pair(iorb)) .or. IsOcc(ilut_parent, ab_pair(jorb))) cycle
                            ! Otherwise, we need to spawn to a different
                            ! determinant. Create it and find the parity.
                            call make_double(nI_parent, nI_child, iel, jel, ab_pair(iorb), ab_pair(jorb), ex, tParity)
                            call EncodeBitDet(nI_child, ilut_child)

                            if (tParity) then
                                parity_fac = 1.0_dp
                            else
                                parity_fac = -1.0_dp
                            end if

                            hash_val = FindWalkerHash(nI_child, size(krylov_ht))
                            temp_node => krylov_ht(hash_val)

                            if (temp_node%ind == 0) then
                                ! If there are no determinants at all with this hash value in krylov_array.
                                cycle
                            else
                                tDetFound = .false.
                                do while (associated(temp_node))
                                    if ( all(ilut_child(0:NIfDBO) == krylov_array(0:NIfDBO,temp_node%ind)) ) then
                                        ! If this CurrentDets determinant has been found in krylov_array.
                                        det_ind = temp_node%ind
                                        tDetFound = .true.
                                        exit
                                    end if
                                    ! Move on to the next determinant with this hash value.
                                    temp_node => temp_node%next
                                end do
                                if (tDetFound) then
                                    int_sign = krylov_array(NOffSgn:NOffSgn+lenof_all_signs-1, det_ind)
                                    child_sign = transfer(int_sign, parent_sign)
                                    if (IsUnoccDet(child_sign)) cycle

                                    ! Finally, add in the contribution to the projected Hamiltonian
                                    ! for each pair of Krylov vectors.
                                    do i = 1, nvecs
                                        do j = i, nvecs
                                            spin_matrix(i,j) = spin_matrix(i,j) - &
                                                1.0_dp*parity_fac*(parent_sign(2*i-1)*child_sign(2*j) + &
                                                                   parent_sign(2*i)*child_sign(2*j-1))/2.0_dp
                                        end do
                                    end do
                                end if
                            end if

                        end if
                    end if

                end do ! Over all occupied electrons > iel (jel).
            end do ! Over all occupied electrons in this determiant (iel).

        end do ! Over all occupied determinants.

        ! Symmetrise the spin matrix.
        do i = 1, nvecs
            do j = 1, i-1
                spin_matrix(i,j) = spin_matrix(j,i)
            end do
        end do

    end subroutine calc_projected_spin

end module ex_state_spin
