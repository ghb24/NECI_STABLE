#include "macros.h"

! Here we evaluate the spin operator using the 2nd quantisation expression
! with integrals.
module test_index_convention_mod
    use fruit, only: assert_equals, assert_true, run_test_case
    use constants, only: dp, n_int
    use SD_spin_purification_mod, only: S2_expval, spin_momentum, &
        get_open_shell, spin_q_num, S2_expval_exc, dyn_S2_expval_exc, &
        nI_invariant_S2_expval_exc
    use sets_mod, only: operator(.in.), set, operator(.U.)
    use orb_idx_mod, only: sigma => calc_spin_raw, alpha, beta, get_spat, sum, &
        SpinProj_t
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
                do I = 1, n_spin_orb
                    if (I == J) cycle
                    do L = 1, n_spin_orb
                        if (L .in. set([I, J])) cycle
                        do K = 1, n_spin_orb
                            if ((K .in. set([I, J])) .or. K == L) cycle
                            exc = Excite_2_t(K, I, L, J)
                            if (sum(sigma(exc%val(1, :))) == sum(sigma(exc%val(2, :)))) then
                                call assert_true(nI_invariant_S2_expval_exc(exc) .isclose. eval_Excite_2_t(exc, g_spin))
                            end if
                        end do
                    end do
                end do
            end do

            call assert_true(S2_expval_exc([2, 4], Excite_0_t()) .isclose. eval_Excite_0_t([2, 4], h_spin, g_spin))
            call assert_true(S2_expval_exc([1, 3], Excite_0_t()) .isclose. eval_Excite_0_t([1, 3], h_spin, g_spin))

            call assert_true(S2_expval_exc([1, 2, 3], Excite_0_t()) .isclose. eval_Excite_0_t([1, 2, 3], h_spin, g_spin))
            call assert_true(S2_expval_exc([1, 2, 3, 4, 6], Excite_0_t()) .isclose. eval_Excite_0_t([1, 2, 3, 4, 6], h_spin, g_spin))
            call assert_true(S2_expval_exc([1, 2, 3, 4], Excite_0_t()) .isclose. eval_Excite_0_t([1, 2, 3, 4], h_spin, g_spin))
            call assert_true(S2_expval_exc([integer::], Excite_0_t()) .isclose. eval_Excite_0_t([integer::], h_spin, g_spin))


            call assert_true(S2_expval_exc([2, 4], Excite_1_t(2, 6)) .isclose. eval_Excite_1_t([2, 4], Excite_1_t(2, 6), h_spin, g_spin))
            call assert_true(S2_expval_exc([1, 3], Excite_1_t(3, 5)) .isclose. eval_Excite_1_t([1, 3], Excite_1_t(3, 5), h_spin, g_spin))
        end block
    end subroutine

    pure function eval_Excite_0_t(nI, h, g) result(res)
        integer, intent(in) :: nI(:)
        procedure(h_t) :: h
        procedure(g_t) :: g
        real(dp) :: res

        integer :: i, j

        res = 0._dp
        do i = 1, size(nI)
            res = res + h(nI(i), nI(i))
        end do

        do i = 1, size(nI)
            do j = 1, size(nI)
                if (i == j) cycle
                res = res + &
                    (g(nI(i), nI(i), nI(j), nI(j)) - g(nI(i), nI(j), nI(j), nI(i))) / 2._dp
            end do
        end do
    end function


    pure function eval_Excite_1_t(nI, exc, h, g) result(res)
        integer, intent(in) :: nI(:)
        type(Excite_1_t), intent(in) :: exc
        procedure(h_t) :: h
        procedure(g_t) :: g
        real(dp) :: res
        integer :: i

        associate(src => exc%val(1, 1), tgt => exc%val(2, 1))
            res = h(src, tgt)
            do i = 1, size(nI)
                associate(R => nI(i))
                    res = res + (g(src, tgt, R, R) - g(src, R, R, tgt))
                end associate
            end do
        end associate
    end function


    real(dp) pure function eval_Excite_2_t(exc, g)
        type(Excite_2_t), intent(in) :: exc
        procedure(g_t) :: g
        associate(I => exc%val(2, 1), J => exc%val(2, 2), K => exc%val(1, 1), L => exc%val(1, 2))
            eval_Excite_2_t = g(I, K, J, L) - g(I, L, J, K)
        end associate
    end function

    real(dp) pure function h_spin(P, Q)
        integer, intent(in) :: P, Q
        if (P == Q) then
            h_spin = 3._dp / 4._dp
        else
            h_spin = 0._dp
        end if
    end function

    real(dp) pure function g_spin(P, Q, R, S)
        integer, intent(in) :: P, Q, R, S
        g_spin = 0._dp
        if (get_spat(P) == get_spat(Q) .and. get_spat(R) == get_spat(S)) then
            if (sigma(P) == sigma(S) .and. sigma(Q) == sigma(R) .and. sigma(P) /= sigma(Q)) then ! abba, baab, i.e. Exchange
                g_spin = 1._dp
            else if (sigma(P) == sigma(Q) .and. sigma(Q) == sigma(R) .and. sigma(R) == sigma(S)) then ! aaaa, bbbb
                g_spin = 1._dp / 2._dp
            else if (sigma(P) == sigma(Q) .and. sigma(R) == sigma(S) .and. sigma(P) /= sigma(R)) then ! aabb, bbaa
                g_spin = -1._dp / 2._dp
            end if
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
