module test_orb_idx_mod
    use fruit
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, SpinProj_t, size, &
        calc_spin, alpha, beta, operator(==)
    use excitation_types, only: NoExc_t, SingleExc_t, DoubleExc_t, excite
    implicit none
    private
    public :: test_calc_spin, test_conversion, test_excite


contains

    subroutine test_calc_spin()
        type(SpinOrbIdx_t) :: orbs
        type(SpinProj_t), allocatable :: expected(:), calculated(:)

        orbs = SpinOrbIdx_t([1, 3, 4, 5, 7])
        expected = [beta, beta, alpha, beta, beta]
        calculated = calc_spin(orbs)
        call assert_equals(expected%val, calculated%val, size(orbs))

    end subroutine

    subroutine test_conversion()
        type(SpatOrbIdx_t) :: orbs
        type(SpinOrbIdx_t) :: expected, calculated
        orbs = SpatOrbIdx_t([1, 3, 4])

        expected = SpinOrbIdx_t([1, 2, 5, 6, 7, 8])
        calculated = SpinOrbIdx_t(orbs)
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([1, 5, 7])
        calculated = SpinOrbIdx_t(orbs, m_s=beta)
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([2, 6, 8])
        calculated = SpinOrbIdx_t(orbs, m_s=alpha)
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([2, 4, 6])
        calculated = SpinOrbIdx_t([1, 2, 3, 4, 5, 6], m_s=alpha)
        call assert_true(all(expected == calculated))
    end subroutine

    subroutine test_excite()
        type(SpinOrbIdx_t) :: reference
        reference = SpinOrbIdx_t([1, 2, 3])

        call assert_true(all(reference == excite(reference, NoExc_t())))
        call assert_true(all(SpinOrbIdx_t([2, 3, 5]) == excite(reference, SingleExc_t(1, 5))))
        call assert_true(all(SpinOrbIdx_t([3, 4, 5]) == excite(reference, DoubleExc_t(1, 5, 2, 4))))

        reference = SpinOrbIdx_t([1, 2, 3, 11, 12, 14])
        call assert_true(all(SpinOrbIdx_t([2, 3, 5, 11, 12, 14]) == excite(reference, SingleExc_t(1, 5))))

    end subroutine

end module test_orb_idx_mod

program test_orb_idx_and_exc_types

    use mpi
    use fruit
    use test_orb_idx_mod, only: test_calc_spin, test_conversion, &
        test_excite



    implicit none
    integer :: failed_count, err

    integer :: n

    call mpi_init(err)

    call init_fruit()

    call test_orb_idx_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_orb_idx_and_exc_types_program', 'failed_tests')

    call mpi_finalize(err)

contains

    subroutine test_orb_idx_driver()
        call run_test_case(test_calc_spin, "test_calc_spin")
        call run_test_case(test_conversion, "test_conversion")
        call run_test_case(test_excite, "test_excite")
    end subroutine
end program test_orb_idx_and_exc_types
