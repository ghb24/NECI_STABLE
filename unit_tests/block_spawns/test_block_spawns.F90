#include "macros.h"

module test_block_spawns_mod
    use fruit, only: run_test_case, assert_true
    use constants, only: n_int
    use block_spawns, only: init_block_list, block_list, is_excitation_aborted,&
                            blocklist_hash_table
    use SystemData, only: nEL, nbasis
    use bit_reps, only: decode_bit_det
    use bit_rep_data, only: nifd
    use util_mod, only: get_free_unit, delete_file
    use hash, only: RandomHash2
    use DetBitOps, only: EncodeBitDet

    better_implicit_none

    private

    public :: block_spawns_test_driver

contains

    subroutine block_spawns_test_driver()
        nEL = 6
        nbasis = 12
        ! Hard coded hash initialization
        allocate(RandomHash2(12), source=0)
        RandomHash2 = [6313, 2856, 5705, 7890, 2242, 503, &
                       5423, 1637, 104, 1121, 3750, 4292]
        call run_test_case(test_init_block_list, "test_init_block_list")
        call run_test_case(test_is_excitation_aborted, &
                           "test_is_excitation_aborted")

    end subroutine block_spawns_test_driver

    subroutine test_init_block_list
        character(len=255) :: blocklist_name, blocklist_content
        integer :: nI(nel), i, iunit
        integer, allocatable :: ref_list(:, :)

        blocklist_name = "blocklist_file"

        iunit = get_free_unit()

        open(newunit=iunit, file=blocklist_name, status='replace', &
             action='write')
        blocklist_content = "1 2 3 4-6" // new_line('A') // &
                            "1 4 5 8 9 12" // new_line('A') // &
                            "2-5 10 11" // new_line('A') // &
                            "4-9"
        write(iunit, '(A)') blocklist_content
        close(iunit)

        call init_block_list(blocklist_name, block_list, blocklist_hash_table)

        open (newunit=iunit, file=blocklist_name, status='old')
        close (iunit, status='delete')

        allocate(ref_list(nel, size(block_list, 2)), source=0)
        ref_list(:, 1) = [1, 2, 3, 4, 5, 6]
        ref_list(:, 2) = [1, 4, 5, 8, 9, 12]
        ref_list(:, 3) = [2, 3, 4, 5, 10, 11]
        ref_list(:, 4) = [4, 5, 6, 7, 8, 9]

        do i = 1, size(block_list, 2)
            call decode_bit_det(nI, block_list(:, i))
            call assert_true(all(ref_list(:, i) == nI))
        end do

    end subroutine test_init_block_list

    subroutine test_is_excitation_aborted()
        integer :: nI(nEL)
        integer(kind=n_int) :: ilut(0:nifd)

        nI = [1, 2, 3, 4, 5, 6]
        call EncodeBitDet(nI, ilut)
        call assert_true(is_excitation_aborted(nI, ilut))
        nI = [1, 4, 5, 8, 9, 12]
        call EncodeBitDet(nI, ilut)
        call assert_true(is_excitation_aborted(nI, ilut))
        nI = [2, 3, 4, 5, 10, 11]
        call EncodeBitDet(nI, ilut)
        call assert_true(is_excitation_aborted(nI, ilut))
        nI = [4, 5, 6, 7, 8, 9]
        call EncodeBitDet(nI, ilut)
        call assert_true(is_excitation_aborted(nI, ilut))
        nI = [2, 5, 6, 7, 8, 9]
        call EncodeBitDet(nI, ilut)
        call assert_true(.not. is_excitation_aborted(nI, ilut))
        nI = [1, 4, 5, 8, 11, 12]
        call EncodeBitDet(nI, ilut)
        call assert_true(.not. is_excitation_aborted(nI, ilut))

    end subroutine test_is_excitation_aborted

end module test_block_spawns_mod

!------------------------------------------------------------------------------!

program test_block_spawns
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, get_failed_count
    use test_block_spawns_mod, only: block_spawns_test_driver

    better_implicit_none

    integer :: failed_count

    call init_fruit()

    call block_spawns_test_driver()

    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
   if (failed_count /= 0) stop -1
end program test_block_spawns
