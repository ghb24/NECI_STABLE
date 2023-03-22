module test_orb_idx_mod
    use fruit, only: assert_true, assert_false
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, SpinProj_t, size, &
        calc_spin, alpha, beta, operator(==)
    use excitation_types, only: Excite_0_t, Excite_1_t, Excite_2_t, excite, is_canonical, &
        canonicalize
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
        call assert_true(all(expected == calculated))
        call assert_true(size(expected) == size(calculated))

    end subroutine

    subroutine test_conversion()
        type(SpatOrbIdx_t) :: orbs
        type(SpinOrbIdx_t) :: expected, calculated
        orbs = SpatOrbIdx_t([1, 3, 4])

        expected = SpinOrbIdx_t([1, 2, 5, 6, 7, 8])
        calculated = SpinOrbIdx_t(orbs)
        call assert_true(all(expected == calculated))

        expected = spinorbidx_t([1, 5, 7])
        calculated = spinorbidx_t(orbs, m_s=beta)
        call assert_true(all(expected == calculated))

        expected = spinorbidx_t([2, 6, 8])
        calculated = spinorbidx_t(orbs, m_s=alpha)
        call assert_true(all(expected == calculated))

        expected = spinorbidx_t([2, 4, 6])
        calculated = spinorbidx_t([1, 2, 3, 4, 5, 6], m_s=alpha)
        call assert_true(all(expected == calculated))
    end subroutine

    subroutine test_excite()
        type(SpinOrbIdx_t) :: reference
        reference = SpinOrbIdx_t([1, 2, 3])

        call assert_true(all(reference == excite(reference, Excite_0_t())))
        call assert_true(all(SpinOrbIdx_t([2, 3, 5]) == excite(reference, Excite_1_t(1, 5))))
        call assert_true(all(SpinOrbIdx_t([3, 4, 5]) == excite(reference, canonicalize(Excite_2_t(1, 5, 2, 4)))))

        call assert_false(is_canonical(Excite_2_t(1, 5, 2)))
        call assert_false(is_canonical(Excite_2_t(1, 5, 2, 4)))
        call assert_true(is_canonical(Excite_2_t(1, 5, 2, 6)))

        reference = SpinOrbIdx_t([1, 2, 3, 11, 12, 14])
        call assert_true(all(SpinOrbIdx_t([2, 3, 5, 11, 12, 14]) == excite(reference, Excite_1_t(1, 5))))

    end subroutine

end module test_orb_idx_mod

program test_orb_idx_and_exc_types

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use test_orb_idx_mod, only: test_calc_spin, test_conversion, &
        test_excite
    use Parallel_neci, only: MPIInit, MPIEnd



    implicit none
    integer :: failed_count
    logical :: err

    call MPIInit(err)

    call init_fruit()

    call test_orb_idx_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_orb_idx_and_exc_types_program', 'failed_tests')

    call MPIEnd(err)

contains

    subroutine test_orb_idx_driver()
        call run_test_case(test_calc_spin, "test_calc_spin")
        call run_test_case(test_conversion, "test_conversion")
        call run_test_case(test_excite, "test_excite")
    end subroutine
end program test_orb_idx_and_exc_types
