#include "macros.h"

module spatial_initiator

    use constants
    use bit_reps, only: NIfDBO, NIfTot, NIfD, extract_part_sign, &
                        encode_part_sign, encode_bit_rep
    use DetBitOps, only: MaskAlpha, MaskBeta, spatial_bit_det
    use FciMCData, only: CurrentInits, no_spatial_init_dets, max_inits
    use util_mod, only: binary_search

    implicit none

contains

    subroutine add_initiator_list (Ilut)
        
        ! Add an initiator to the spatial initiator list.
        !
        ! In: IlutI - The spin orbital bit representation of the determinant
        !             to add to the spatial initiator list.

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer(n_int), dimension(0:NIfD) :: alpha, beta, a_sft, b_sft
        integer(n_int) :: spat(0:NIfTot)
        real(dp) :: sgn(lenof_sign)
        integer :: flag, pos
        
        ! Obtain standard spatial representation
        spat = spatial_bit_det (ilut)

        if (no_spatial_init_dets == 0) then
            ! If list is empty, start it.
            sgn(1) = 1
            call encode_bit_rep (CurrentInits(:,1), spat, sgn, flag)
            no_spatial_init_dets = 1
        else
            pos = binary_search(CurrentInits(:, 1:no_spatial_init_dets), &
                                spat(0:NIfD), NIfD+1)
            if (pos < 0) then
                ! If det not in the list, add it.
                sgn(1) = 1
                CurrentInits(:, (-pos)+1:no_spatial_init_dets+1) = &
                    CurrentInits(:, (-pos):no_spatial_init_dets)
                call encode_bit_rep (CurrentInits(:,-pos), spat, sgn, flag)
                no_spatial_init_dets = no_spatial_init_dets + 1
            else
                ! If we have found the det already, then increment its count.
                sgn(1) = extract_part_sign (CurrentInits(:, pos), 1) + 1
                call encode_part_sign (CurrentInits(:, pos), sgn(1), 1)
            endif
        endif

    end subroutine

    subroutine rm_initiator_list (ilut)

        ! Remove an initiator from the spatial initiator list.
        !
        ! In: Ilut - The spin orbital bit representation of the determinant
        !             to remove from the spatial initiator list.

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer(n_int), dimension(0:NIfD) :: alpha, beta, a_sft, b_sft
        integer(n_int) :: spat(0:NIfTot)
        real(dp) :: sgn
        integer :: pos
        character(*), parameter :: this_routine = 'rm_initiator_list'

        ! Obtain standard spatial representation
        spat = spatial_bit_det (ilut)

        ! Find the spatial initiator in the list. If it is not there, then
        ! something is very wrong.
        pos = binary_search(CurrentInits(:, 1:no_spatial_init_dets), &
                            spat(0:NIfD), NIfD+1)
        if (pos < 0) &
            call stop_all (this_routine, "Spatial initiator to remove from &
                                         &list not found.")
        
        sgn = extract_part_sign (CurrentInits(:,pos), 1)
        ASSERT(sgn > 0)
        if (sgn == 1) then
            ! Remove this determinant from the list.
            if (no_spatial_init_dets > 1) then
                CurrentInits(:, pos:no_spatial_init_dets - 1) = &
                    CurrentInits(:, pos+1:no_spatial_init_dets)
            endif
            no_spatial_init_dets = no_spatial_init_dets - 1
        else
            ! If more than one determinant with this spatial configuration is
            ! in the list, just decrement the count.
            call encode_part_sign (CurrentInits(:,pos), sgn - 1, 1)
        endif

    end subroutine

    function is_spatial_init (ilut) result (bInit)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer(n_int), dimension(0:NIfD) :: alpha, beta, a_sft, b_sft
        integer(n_int) :: spat(0:NIfTot)
        integer :: pos
        logical :: bInit
        
        ! By default, we assume it isn't a spatial initiator
        bInit = .false.
        if (no_spatial_init_dets == 0) return

        ! Obtain standard spatial representation
        spat = spatial_bit_det (ilut)

        ! If we find the spatial det in the list, then this is a spatial
        ! initiator
        pos = binary_search(CurrentInits(:, 1:no_spatial_init_dets), &
                            spat(0:NIfD), NIfD+1)
        if (pos > 0) bInit = .true.

    end function

end module

