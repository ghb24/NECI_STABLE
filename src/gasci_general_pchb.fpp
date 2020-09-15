#include "macros.h"
#:include "macros.fpph"
#:include "algorithms.fpph"

module gasci_general_pchb
    use constants, only: n_int, dp, int64, maxExcit
    use util_mod, only: choose, factrl, stop_all
    use FciMCData, only: excit_gen_store_type
    use SystemData, only: nEl
    use bit_rep_data, only: NIfTot

    use gasci, only: GAS_specification, GASSpec_t

    implicit none

    private

    public :: gen_general_GASCI_pchb, get_partitions

contains

    !>  @brief
    !>  The disconnected_GAS_PCHB excitation generator subroutine.
    !>
    !>  @details
    !>  This is a wrapper around `disconnected_GAS_PCHB%gen_excit`
    !>  to match the function pointer interface.
    !>  The interface is common to all excitation generators, see proc_ptrs.F90
    subroutine gen_general_GASCI_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex_mat, tParity, pGen, hel, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex_mat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'gen_GASCI_pchb'

        @:unused_var(exFlag, part_type, store)
        @:ASSERT(GAS_specification%contains(nI))

    end subroutine gen_general_GASCI_pchb

    !> @brief
    !> Get the ordered partitions of n into k summands.
    !>
    !> @details
    !> Get all possible solutions for the k dimensional hypersurface.
    !> \f[x_1 + ... + x_k = n  \f]
    !> by taking into account the order.
    !> \f[ 1 + 0 = 1 \f] is different from
    !> \f[ 0 + 1 = 1 \f].
    !> The German wikipedia has a nice article
    !> https://de.wikipedia.org/wiki/Partitionsfunktion#Geordnete_Zahlpartitionen
    pure function get_partitions(k, n) result(res)
        integer, intent(in) :: k, n
        integer(int64) :: res(k, n_partitions(k, n))
        integer :: idx_part, j

        idx_part = 1
        res(:, idx_part) = 0_int64
        res(1, idx_part) = n

        do idx_part = 2, size(res, 2) - 1
            res(:, idx_part) = res(:, idx_part - 1)

            do j = size(res, 1), 2, -1
                if (res(j - 1, idx_part) >= res(j, idx_part) .and. res(j - 1, idx_part) > 0) exit
            end do

            ! Transfer 1 from left and everything from right to res(j)
            res(j, idx_part) = res(j, idx_part) + 1_int64 + sum(res(j + 1 :, idx_part))
            res(j + 1 :, idx_part) = 0_int64
            res(j - 1, idx_part) = res(j - 1, idx_part) - 1_int64
        end do

        res(:, idx_part) = 0_int64
        res(size(res, 1), idx_part) = n
    end function

    pure elemental function n_partitions(k, n) result(res)
        integer, intent(in) :: k, n
        integer :: res
        res = choose(n + k - 1, k - 1)
    end function

end module gasci_general_pchb
