#include "macros.h"

module test_index_convention_mod
    use fruit, only: assert_equals, assert_true, run_test_case
    use constants, only: dp, n_int
    use SD_spin_purification_mod, only: S2_expval, spin_momentum, &
        get_open_shell, spin_q_num, S2_expval_exc, dyn_S2_expval_exc, &
        nI_invariant_S2_expval_exc
    use sets_mod, only: operator(.in.)
    use orb_idx_mod, only: sigma => calc_spin_raw, alpha, beta, get_spat
    use excitation_types, only: Excitation_t, Excite_0_t, Excite_1_t, Excite_2_t, Excite_3_t, &
        is_canonical
    use util_mod, only: operator(.isclose.)
    implicit none
    private
    public :: test_driver

    abstract interface
        real(dp) pure function h_t(P, Q)
            import :: dp
            implicit none
            integer, intent(in) :: P, Q
        end function

        real(dp) pure function g_t(P, Q, R, S)
            import :: dp
            implicit none
            integer, intent(in) :: P, Q, R, S
        end function
    end interface

contains

    subroutine test_driver()
        call run_test_case(test_spin_from_integral, "test_spin_from_integral")
    end subroutine

    subroutine test_spin_from_integral
        call assert_true(h_spin(2, 2) .isclose. h_spin(4, 4))

        block
            integer :: I, J, K, L
            integer, parameter :: n_spin_orb = 4
            type(Excite_2_t) :: exc

            do J = 1, n_spin_orb
                do I = 1, J - 1
                    write(*, *) g_spin(I, I, J, J), g_spin(I, J, J, I), I, J, J, I
                    do L = 1, n_spin_orb
                        if (L .in. [I, J]) cycle
                        do K = 1, L - 1
                            if (K .in. [I, J]) cycle
                            exc = Excite_2_t(K, I, L, J)
                            if (is_canonical(exc)) then
                                if (.not. (nI_invariant_S2_expval_exc(exc) .isclose. eval_Excite_2_t(exc, g_spin))) then
                                    write(*, *)
                                    write(*, *) exc%val(1, :)
                                    write(*, *) exc%val(2, :)
                                    write(*, *) nI_invariant_S2_expval_exc(exc)
                                    write(*, *) eval_Excite_2_t(exc, g_spin)
                                    write(*, *)
                                end if
                                call assert_true(nI_invariant_S2_expval_exc(exc) .isclose. eval_Excite_2_t(exc, g_spin))
                            end if
                        end do
                    end do
                end do
            end do
        end block
    end subroutine

    real(dp) pure function eval_Excite_2_t(exc, g)
        type(Excite_2_t), intent(in) :: exc
        procedure(g_t) :: g
        associate(I => exc%val(2, 1), J => exc%val(2, 2), K => exc%val(1, 1), L => exc%val(1, 2))
            eval_Excite_2_t = g(I, K, J, L) - g(I, L, J, K)
        end associate
    end function
    !
    ! real(dp) pure function eval_Excite_1_t(nI, exc, h, g)
    !     type(Excite_2_0), intent(in) :: exc
    !     procedure(g_t) :: g
    !     associate(I => exc%val(2, 1), J => exc%val(2, 2), K => exc%val(1, 1), L => exc%val(1, 2))
    !         eval_Excite_2_t = g(I, K, J, L) - g(I, L, J, K)
    !     end associate
    ! end function

    real(dp) pure function h_spin(P, Q)
        integer, intent(in) :: P, Q
        if (P == Q .and. sigma(P) == alpha) then
            h_spin = 1._dp
        else
            h_spin = 0._dp
        end if
    end function

    real(dp) pure function g_spin(P, Q, R, S)
        integer, intent(in) :: P, Q, R, S
        ! if (sigma(P) == sigma(Q) &
        !         .and. sigma(R) == sigma(S) &
        !         .and. get_spat(P) == get_spat(S) &
        !         .and. get_spat(R) == get_spat(Q) &
        !     ) then
        !     if (.true.) then
        !         g_spin = -1._dp
        !     else
        !         g_spin = -1._dp
        !     end if
        !
        if (sigma(P) == sigma(Q) .and. sigma(S) == sigma(R) &
            .and. get_spat(P) == get_spat(S) .and. get_spat(R) == get_spat(Q) &
            .and. get_spat(P) /= get_spat(R)) then
            g_spin = -1._dp
        else
            g_spin = 0._dp
        end if
    end function

end module

program test_index_convention_prog

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd

    use test_index_convention_mod, only: test_driver

    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_index_convention_prog', 'failed_tests')

        call MPIEnd(.false.)
    end block
end program
