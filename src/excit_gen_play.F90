module excit_gens

    use SystemData, only: nel
    use SymExcit3, only: CountExcitations3, GenExcitations3
    use procedure_pointers, only: excit_gen_store_type, spawned_info_t
    use dSFMT_interface, only: genrand_real2_dSFMT
    use Determinants, only: get_helement, write_det
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use bit_rep_data, only: NIfTot
    use bit_reps, only: decode_bit_det
    use constants
    use sort_mod
    implicit none

contains

    subroutine gen_excit_hel_weighted (nI, ilutI, nJ, ilutJ, exFlag, store, &
                                       spawned_info)

        ! A really laborious, slow, explicit and brute force method to
        ! generating all excitations in proportion to their connection
        ! strength. This demonstrates the maximum possible value of tau that
        ! can be used.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in), target :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        type(excit_gen_store_type), intent(inout), target :: store
        type(spawned_info_t), intent(inout), target :: spawned_info
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_hel_weighted'

        integer :: nsing, ndoub, nexcit

        ! Count how many singles and doubles there are!
        call CountExcitations3 (nI, 3, nsing, ndoub)
        nexcit = nsing + ndoub

        call gen_excit_hel_local (nI, ilutI, nJ, ilutJ, exFlag, store, &
                                  spawned_info, nexcit)

    end subroutine


    subroutine gen_excit_hel_local (nI, ilutI, nJ, ilutJ, exFlag, store, &
                                    spawned_info, nexcit)

        ! A really laborious, slow, explicit and brute force method to
        ! generating all excitations in proportion to their connection
        ! strength. This demonstrates the maximum possible value of tau that
        ! can be used.

        integer, intent(in) :: nI(nel), exFlag, nexcit
        integer(n_int), intent(in), target :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel)
        type(excit_gen_store_type), intent(inout), target :: store
        type(spawned_info_t), intent(inout), target :: spawned_info
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        character(*), parameter :: this_routine = 'gen_excit_hel_weighted'

        integer(n_int) :: iluts(0:NIfTot, nexcit)
        HElement_t :: hels(nexcit), hel_sum, hel_cum
        integer :: excit_count, ex(2,2), i, flag
        logical :: found_all, par

        ! Generate two lists. One with all of the available excitations, and
        ! one with their HElement values
        excit_count = 0
        found_all = .false.
        hel_sum = 0
        nJ = 0
        flag = 3
        ex = 0
        call GenExcitations3(nI, ilutI, nJ, flag, ex, par, found_all, &
                             .false.)
        do while (.not. found_all)
            excit_count = excit_count + 1
            call EncodeBitDet(nJ, iluts(:, excit_count))
            hels(excit_count) = abs(get_helement(nI, nJ))
            hel_sum = hel_sum + hels(excit_count)
            call GenExcitations3(nI, ilutI, nJ, flag, ex, par, found_all, &
                                 .false.)
        end do

        if (excit_count /= nexcit) &
            call stop_all(this_routine,"Incorrect number of excitations found")

        ! Sort the lists!!!
        call sort(hels, iluts)

        ! Pick a random cumulative helement value
        hel_cum = hel_sum * genrand_real2_dSFMT()

        ! Work from the end, and stop when we have got where we need to be!
        do i = nexcit, 1, -1
            hel_cum = hel_cum - hels(i)
            if (hel_cum <= 0) exit
        end do

        ! Just in case we get shafted by rounding errors
        if (i < 1) i = 1

        ! Now we know what we are returning
        ilutJ = iluts(:, i)
        call decode_bit_det (nJ, ilutJ)
        spawned_info%ex(1,1) = 2
        call GetBitExcitation (ilutI, ilutJ, spawned_info%ex, &
                               spawned_info%tParity)
        spawned_info%ic = FindBitExcitLevel(ilutI, ilutJ, 2)
        spawned_info%pgen = hels(i) / hel_sum

    end subroutine

end module
