module test_gasci_mod
    use fruit
    use gasci, only: GAS_specification_t, iGAS_from_spinorb, iGAS_from_spatorb, &
        contains_nI, get_nGAS
    implicit none
    private
    public :: test_igas_from_spatorb, test_igas_from_spinorb, test_contains_nI

contains

        ! integer, parameter :: nEl = 4, nDets = 5
        ! integer, parameter :: &
        !   nI(nEl, nDets) = reshape([1, 1, 3, 3, &
        !                             1, 2, 3, 4, &
        !                             1, 1, 2, 2, &
        !                             1, 1, 2, 3, &
        !                             1, 3, 3, 4], [nEl, nDets])

    subroutine test_igas_from_spatorb()
        type(GAS_specification_t) :: GAS_spec

        GAS_spec%n_orbs = [2, 4]
        GAS_spec%n_min = [2, 4]
        GAS_spec%n_max = [2, 4]

        call assert_equals([1, 1, 2, 2], f([1, 2, 3, 4]), 4)

        contains

        elemental function f(orbidx) result(iGAS)
            integer, intent(in) :: orbidx
            integer :: iGAS
            iGAS = iGAS_from_spatorb(GAS_spec, orbidx)
        end function
    end subroutine

    subroutine test_igas_from_spinorb()
        type(GAS_specification_t) :: GAS_spec

        GAS_spec%n_orbs = [2, 4]
        GAS_spec%n_min = [2, 4]
        GAS_spec%n_max = [2, 4]

        call assert_equals([1, 1, 1, 1, 2, 2, 2, 2], &
                          f([1, 2, 3, 4, 5, 6, 7, 8]), 8)

        contains

        elemental function f(orbidx) result(iGAS)
            integer, intent(in) :: orbidx
            integer :: iGAS
            iGAS = iGAS_from_spinorb(GAS_spec, orbidx)
        end function
    end subroutine

    subroutine test_contains_nI()
        type(GAS_specification_t) :: GAS_spec

        GAS_spec%n_orbs = [2, 4]
        GAS_spec%n_min = [2, 4]
        GAS_spec%n_max = [2, 4]

        call assert_equals(.true., contains_nI(GAS_spec, [1, 2, 5, 6]))
        call assert_equals(.true., contains_nI(GAS_spec, [1, 3, 5, 6]))

        call assert_equals(.false., contains_nI(GAS_spec, [1, 2, 3, 4]))
        call assert_equals(.false., contains_nI(GAS_spec, [1, 2, 3, 5]))

        call assert_equals(.false., contains_nI(GAS_spec, [5, 6, 7, 8]))
        call assert_equals(.false., contains_nI(GAS_spec, [1, 6, 7, 8]))
    end subroutine

end module test_gasci_mod

program test_gasci_program

    use mpi
    use fruit
    use test_gasci_mod, only: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_nI

    implicit none
    integer :: failed_count, err

    call mpi_init(err)

    call init_fruit()

    call test_gasci_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

    call mpi_finalize(err)

contains

    subroutine test_gasci_driver()
        call run_test_case(test_igas_from_spatorb, "test_igas_from_spatorb")
        call run_test_case(test_igas_from_spinorb, "test_igas_from_spinorb")
        call run_test_case(test_contains_nI, "test_contains_nI")
    end subroutine
end program test_gasci_program
