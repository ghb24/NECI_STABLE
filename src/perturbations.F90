module perturbations

    use constants
    use FciMCData, only: perturbation

contains

    subroutine init_perturbation_annihilation(pert)

        ! Create the 'ann_elems' and 'ann_bits' components of a perturbation object.
        ! The 'ann_orbs' components should be allocated and filled in already, as
        ! should the nannihilate component.

        type(perturbation), intent(inout) :: pert
        integer :: i

        if (.not. allocated(pert%ann_elems)) allocate(pert%ann_elems(pert%nannihilate))
        if (.not. allocated(pert%ann_bits)) allocate(pert%ann_bits(pert%nannihilate))
        do i = 1, pert%nannihilate
            pert%ann_elems(i) = (pert%ann_orbs(i)-1)/bits_n_int
            pert%ann_bits(i) = mod(pert%ann_orbs(i)-1, bits_n_int)
        end do

    end subroutine init_perturbation_annihilation

    subroutine init_perturbation_creation(pert)

        ! Create the 'crtn_elems' and 'crtn_bits' components of a perturbation object.
        ! The 'crtn_orbs' components should be allocated and filled in already, as
        ! should the ncreate component.

        type(perturbation), intent(inout) :: pert
        integer :: i

        if (.not. allocated(pert%crtn_elems)) allocate(pert%crtn_elems(pert%ncreate))
        if (.not. allocated(pert%crtn_bits)) allocate(pert%crtn_bits(pert%ncreate))
        do i = 1, pert%ncreate
            pert%crtn_elems(i) = (pert%crtn_orbs(i)-1)/bits_n_int
            pert%crtn_bits(i) = mod(pert%crtn_orbs(i)-1, bits_n_int)
        end do

    end subroutine init_perturbation_creation

    subroutine apply_perturbation_array(perturbs, ndets, dets_in, dets_out)

        use bit_rep_data, only: NIfTot
        use bit_reps, only: add_ilut_lists

        type(perturbation), intent(in) :: perturbs(:)
        integer, intent(inout) :: ndets
        integer(n_int), intent(in) :: dets_in(0:,:) ! First dimension must be 0:NIfTot
        integer(n_int), intent(out) :: dets_out(0:,:) ! First dimension must be 0:NIfTot

        integer :: i, ndets_init, ndets_pert_1, ndets_pert_2
        integer(n_int), allocatable :: temp_dets_1(:,:), temp_dets_2(:,:)

        ndets_init = ndets

        if (size(perturbs) > 1) then
            ! Just allocate these arrays to be the same size as the CurrentDets
            ! array.
            allocate(temp_dets_1(0:NIfTot, MaxWalkersPart), stat=ierr) 
            allocate(temp_dets_2(0:NIfTot, MaxWalkersPart), stat=ierr) 

            ! Apply the first perturbation.
            ndets_pert_1 = ndets_init
            ! Note that ndets_pert_1 is altered by this routine.
            call apply_perturbation(perturbs(1), ndets_pert_1, dets_in, temp_dets_1)

            do i = 2, size(perturbs)
                ndets_pert_2 = ndets_init
                ! Note that ndets_pert_2 is altered by this routine.
                call apply_perturbation(perturbs(i), ndets_pert_2, dets_in, temp_dets_2)
                call add_ilut_lists(ndets_pert_1, ndets_pert_2, .false., temp_dets_1, temp_dets_2, dets_out, ndets)
                ! If we still have more perturbations to apply, copy the result
                ! to temp_dets_1. Else, exit the final result in dets_out.
                if (i /= size(perturbs)) then
                    ndets_pert_1 = ndets
                    temp_dets_1(:,1:ndets) = dets_out(:,1:ndets)
                end if
            end do

            deallocate(temp_dets_1)
            deallocate(temp_dets_2)
        else if (size(perturbs) == 1) then
            ! Simply apply the single perturbation using the final input and
            ! output arrays.
            call apply_perturbation(perturbs(1), ndets, dets_in, dets_out)
        end if

    end subroutine apply_perturbation_array

    subroutine apply_perturbation(perturb, ndets, dets_in, dets_out)

        ! Take in a list of determinants (dets_in) and apply a pertubation
        ! to each determinant. As we go we shuffle down determinants to fill in
        ! gaps opened up by removed determinants.

        ! WARNING: This routine uses SpanwedParts and the annihilation routine
        ! SendProcNewParts to perform communication of the perturbed
        ! determinants. It will overwrite whatever is in SpawnedParts.

        use AnnihilationMod, only: SendProcNewParts
        use bit_rep_data, only: NIfTot, NIfDBO, extract_sign
        use bit_reps, only: encode_sign, decode_bit_det
        use DetBitOps, only: ilut_lt, ilut_gt
        use hash, only: reset_hash_table, fill_in_hash_table, DetermineDetNode
        use FciMCData, only: tHashWalkerList, HashIndex, SpawnedParts, SpawnedParts2
        use FciMCData, only: ValidSpawnedList, InitialSpawnedSlots
        use sort_mod, only: sort
        use SystemData, only: nel

        type(perturbation), intent(in) :: perturb
        integer, intent(inout) :: ndets
        integer(n_int), intent(in) :: dets_in(0:,:) ! First dimension must be 0:NIfTot
        integer(n_int), intent(out) :: dets_out(0:,:) ! First dimension must be 0:NIfTot

        integer(n_int) :: ilut(0:NIfTot)
        integer(n_int), pointer :: PointTemp(:,:)
        integer :: i, nremoved, proc
        integer :: nI(nel)

        ! If the perturbation is the identity operator then just return.
        if (perturb%nannihilate == 0 .and. perturb%ncreate == 0) return

        nremoved = 0
        ! Reset the spawning slot positions in SpawnedParts.
        ValidSpawnedList = InitialSpawnedSlots

        do i = 1, ndets
            ! Copy the ilut so that we don't alter the input list.
            ilut = dets_in(:,i)
            call perturb_det(ilut, perturb)
            
            if (all(ilut(0:NIfDBO) == 0_n_int)) then
                nremoved = nremoved + 1
            else
                call decode_bit_det(nI, ilut)
                proc = DetermineDetNode(nel,nI,0)
                SpawnedParts(:, ValidSpawnedList(proc)) = ilut
                ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
            end if
        end do

        ndets = ndets - nremoved

        ! Send perturbed determinants to their new processors.
        call SendProcNewParts(ndets, tSingleProc=.false.)
 
        ! The result of SendProcNewParts is now stored in an array pointed to
        ! by SpawnedParts2. We want it in an array pointed to by SpawnedParts,
        ! so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp
        nullify(PointTemp)

        ! Now move the contents of SpawnedParts to ilut_list.
        do i = 1, ndets
            dets_out(:,i) = SpawnedParts(:,i)
        end do

        call sort(dets_out(:, 1:ndets), ilut_lt, ilut_gt)

    end subroutine apply_perturbation

    subroutine perturb_det(ilut, perturb)

        ! This routine takes a determinant encoded in ilut and applies a
        ! collection of creation and annihilation operators, which together is
        ! referred to as the perturbation operator.

        ! ***IMPORTANT*** The ordering of the annihilation and creation within
        ! the perturbation operator is assumed to be as follows. All
        ! annihilation operators are applied before all creation operators.
        ! Annihilation and creation operators are themselves ordered so that
        ! operators for orbitals with the lowest orbital number are applied
        ! last (leftmost), and higher orbitals numbers applied first
        ! (rightmost). In NECI, the ordering with respect to alpha and beta is
        ! (beta, alpha, beta, alpha...). Note that the order of the creation and
        ! annihilation operators in the input arrays is unimportant.

        ! The a_bits(i)'th bit of ilut(a_elems(i)) should be the bit of the
        ! orbital for the i'th annihilation operator, and similarly for
        ! c_bits and c_elems for creation operators.

        use bit_rep_data, only: NIfTot, NIfDBO, NIfD, extract_sign
        use bit_reps, only: encode_sign
        use DetBitOps, only: CountBits

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        type(perturbation), intent(in) :: perturb

        integer(n_int) :: combined_mask(0:NIfD), new_mask(0:NIfD), ones
        integer(n_int) :: original_ilut(0:NIfTot), masked_ilut(0:NIfD)
        integer :: i, j, num_minus_signs
        real(dp) :: new_sign_factor, real_sign(lenof_sign)

        ! If the perturbation is the identity operator then just return.
        if (perturb%nannihilate == 0 .and. perturb%ncreate == 0) return

        original_ilut = ilut

        associate(a_elems => perturb%ann_elems, a_bits => perturb%ann_bits, &
                  c_elems => perturb%crtn_elems, c_bits => perturb%crtn_bits)

            ! First, loop over all creation and annihilation operators and see if they
            ! destroy the determinant encoded in ilut.

            do i = 1, perturb%nannihilate
                if ( .not. btest(ilut(a_elems(i)), a_bits(i)) ) then
                    ilut(0:NIfDBO) = 0_n_int
                    return
                end if
                ilut(a_elems(i)) = ibclr(ilut(a_elems(i)), a_bits(i))
            end do

            do i = 1, perturb%ncreate
                if ( btest(ilut(c_elems(i)), c_bits(i)) ) then
                    ilut(0:NIfDBO) = 0_n_int
                    return
                end if
                ilut(c_elems(i)) = ibset(ilut(c_elems(i)), c_bits(i))
            end do

            ! If we get to this point then the determinant was not destroyed by the
            ! creation and annihilation operators, else we would have returned.

            ! Now, determine whether or not the application of these operators
            ! introduces a minus sign. We consider the application of the operators
            ! when they are already ordered by the basis ordering used in NECI.
            ! All annihilation operators are applied before creation operators.

            ! We need to consider moving the creation and annihilation operators
            ! in the perturbing operator through the creation operators which
            ! create the original determinant (by acting on the vacuum).

            ! For a creation or annihilation operator for orbital p in the
            ! perturbing operator, we want to consider moving it through creation
            ! operators in the original determinant for orbitals less than p.
            ! Suppose there are N such operators. Then the ordering introduces a
            ! factor of (-1)^N. So we just need to count the number of orbitals
            ! before orbital p in the determinant. We can do this by creating a bit
            ! mask with 1's for all orbitals *before* (and not including) p in the
            ! original determinant, and then performing an and operation between
            ! this mask and the ilut, and then counting the number of bits set in
            ! the resulting bitstring.

            ! However, in general we want to consider several annihilation and
            ! creation operators in the perturbing operator. Instead of doing the
            ! above for each operator, which is slow, we can do it all in one go.
            ! If a creation operator in the original ilut is 'passed' an odd number
            ! of times then we want to count it, otherwise we don't. So we want to
            ! create a bit mask that has 1's for all orbitals that would be passed
            ! an odd number of times by the applied creation and annihilation
            ! operators in the perturbation operator, and 0's for all other
            ! orbitals. This can be done by performing an ieor on the bit masks
            ! for induvidual operators. Once an exclusive or has been performed for
            ! each annihilation and creation operator in the perturbing operator,
            ! we can finally perform an and with the original ilut,and then count
            ! the number of set bits. If this is odd, include a minus sign.

            ones = not(0_n_int)
            combined_mask = 0_n_int

            do i = 1, perturb%nannihilate
                do j = 0, NIfD
                    ! Create a mask which has 0's in all orbitals which might be
                    ! passed in the ordering, and 1's elsewhere. 
                    if (j < a_elems(i)) then
                        new_mask(j) = 0_n_int
                    else if (j == a_elems(i)) then
                        new_mask(j) = ishft(ones, a_bits(i))
                    else
                        new_mask(j) = ones
                    end if
                end do

                ! Now make it so that it has 1's for all orbitals that might be
                ! passed and 0's elsewhere.
                new_mask = not(new_mask)

                ! Now perform an ieor to combine this mask with the combined masks
                ! for all previously considered orbitals.
                combined_mask = ieor(combined_mask, new_mask)
            end do

            ! Now do exactly the same for the creation operators.
            do i = 1, perturb%ncreate
                do j = 0, NIfD
                    if (j < c_elems(i)) then
                        new_mask(j) = 0_n_int
                    else if (j == c_elems(i)) then
                        new_mask(j) = ishft(ones, c_bits(i))
                    else
                        new_mask(j) = ones
                    end if
                end do
                new_mask = not(new_mask)
                combined_mask = ieor(combined_mask, new_mask)
            end do

            ! Finally apply the mask and count the resulting number of set bits...
            masked_ilut = iand(combined_mask, original_ilut(0:NIfD))
            num_minus_signs = CountBits(masked_ilut, NIfD)
            ! ...and apply the new sign.
            !new_sign_factor = (-1_n_int)**(num_minus_signs)
            new_sign_factor = (-1.0_dp)**num_minus_signs
            call extract_sign(ilut, real_sign)
            real_sign = real_sign*new_sign_factor
            call encode_sign(ilut, real_sign)

            call extract_sign(ilut, real_sign)

        end associate

    end subroutine perturb_det

end module perturbations
