module test_orb_idx_mod
    use fruit
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, Spin_t, size, &
        calc_spin, spin => spin_values, operator(==)
    implicit none
    private
    public :: test_calc_spin, test_conversion


contains

    subroutine test_calc_spin()
        type(SpinOrbIdx_t) :: orbs
        type(Spin_t) :: expected
        orbs = SpinOrbIdx_t([1, 3, 4, 5, 7])
        expected = Spin_t([spin%beta, spin%beta, spin%alpha, spin%beta, spin%beta])
! TODO(Oskar, Kai): AAAAAAHHHHH
!         expected = Spin_t([1, 2, 3, 4, 5])
!         expected%m_s = [1, 2, 3, 4, 5]

        associate(calculated => calc_spin(orbs))
            write(*, *) expected%m_s
            write(*, *) calculated%m_s
            call assert_equals(expected%m_s, calculated%m_s, size(orbs))
        end associate

    end subroutine

    subroutine test_conversion()
        type(SpatOrbIdx_t) :: orbs
        orbs = SpatOrbIdx_t([1, 3, 4])

        associate(expected => SpinOrbIdx_t([1, 2, 5, 6, 7, 8]), &
                  calculated => SpinOrbIdx_t(orbs))
            call assert_equals(expected%idx, calculated%idx, size(calculated))
        end associate

        associate(expected => SpinOrbIdx_t([1, 5, 7]), &
                  calculated => SpinOrbIdx_t(orbs, spins=Spin_t([spin%beta])))
            call assert_equals(expected%idx, calculated%idx, size(calculated))
        end associate

        associate(expected => SpinOrbIdx_t([2, 6, 8]), &
                  calculated => SpinOrbIdx_t(orbs, spins=Spin_t([spin%alpha])))
            call assert_equals(expected%idx, calculated%idx, size(calculated))
        end associate

        associate(expected => SpinOrbIdx_t([1, 6, 7]), &
                  calculated => &
                        SpinOrbIdx_t(orbs, &
                            spins=Spin_t([spin%beta, spin%alpha, spin%beta])))
            call assert_equals(expected%idx, calculated%idx, size(calculated))
        end associate

        associate(expected => SpinOrbIdx_t([2, 4, 6]), &
                  calculated => &
                        SpinOrbIdx_t([1, 2, 3, 4, 5, 6], Spin_t([spin%alpha])))
            call assert_true(all(expected == calculated))
        end associate
    end subroutine

end module test_orb_idx_mod

program test_orb_idx

    use mpi
    use fruit
    use test_orb_idx_mod, only: test_calc_spin, test_conversion


    implicit none
    integer :: failed_count, err

    integer :: n

    call mpi_init(err)

    call init_fruit()

    call test_orb_idx_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_orb_idx_program', 'failed_tests')

    call mpi_finalize(err)

contains

    subroutine test_orb_idx_driver()
        call run_test_case(test_calc_spin, "test_calc_spin")
!         call run_test_case(test_conversion, "test_conversion")
    end subroutine
end program test_orb_idx
