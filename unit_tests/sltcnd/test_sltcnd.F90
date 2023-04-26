module test_sltcnd_mod
    use fruit, only: assert_true, assert_false, run_test_case
    use constants, only: dp, n_int
    use sets_mod, only: operator(.complement.)

    use bit_reps, only: decode_bit_det
    use unit_test_helper_excitgen, only: &
        init_excitgen_test, finalize_excitgen_test, &
        generate_random_integrals, RandomFcidumpWriter_t
    use symexcit3, only: gen_excits
    use excitation_types, only: get_excitation, Excitation_t, Excite_0_t, Excite_1_t, Excite_2_t
    use sltcnd_mod, only: sltcnd_excit, diagH_after_exc
    use util_mod, only: operator(.isclose.)

    use sltcnd_mod, only: get_gamma
    implicit none
    private
    public :: test_sltcnd_driver

    integer :: i
    integer, parameter :: system_size = 30
    integer, parameter :: nI(system_size) = [(i, i = 1, system_size)], n_spat_orb = system_size

contains

    subroutine test_diagH_from_exc()
        integer :: i

        integer(n_int), allocatable :: det_list(:, :)
        class(Excitation_t), allocatable :: exc
        integer :: nJ(size(nI)), n_excits, ic
        logical :: tParity

        call init_excitgen_test(nI, &
            RandomFcidumpWriter_t(&
                n_spat_orb, nI, sparse=0.7_dp, sparseT=0.7_dp) &
        )

        associate(E_0 => sltcnd_excit(nI, Excite_0_t()))
            do ic = 1, 2
                call gen_excits(nI, n_excits, det_list, ic)
                do i = 1, size(det_list, 2)
                    call decode_bit_det(nJ, det_list(:, i))
                    call get_excitation(nI, nJ, ic, exc, tParity)

                    select type (exc)
                    type is (Excite_1_t)
                        call assert_true(sltcnd_excit(nJ, Excite_0_t()) .isclose. diagH_after_exc(nI, E_0, exc))
                    type is (Excite_2_t)
                        call assert_true(sltcnd_excit(nJ, Excite_0_t()) .isclose. diagH_after_exc(nI, E_0, exc))
                    end select
                end do
            end do
        end associate

        call finalize_excitgen_test()
    end subroutine

    subroutine test_timing()
        use timing_neci, only: timer, set_timer, halt_timer, get_total_time
        integer :: i
        type(timer) :: time_direct, time_from_exc

        type(Excite_2_t), allocatable :: exc(:)
        integer, allocatable :: nJs(:, :)

        call init_excitgen_test(nI, &
            RandomFcidumpWriter_t(&
                n_spat_orb, nI, sparse=1.0_dp, sparseT=1.0_dp) &
        )

        block
            integer(n_int), allocatable :: det_list(:, :)
            integer :: n_excits

            call gen_excits(nI, n_excits, det_list, 2)
            allocate(nJs(size(nI), n_excits), exc(n_excits))
            do i = 1, size(det_list, 2)
                call decode_bit_det(nJs(:, i), det_list(:, i))
                exc(i)%val(1, :) = nI .complement. nJs(:, i)
                exc(i)%val(2, :) = nJs(:, i) .complement. nI
            end do
        end block

        block
            real(dp) :: E_accum_direct
            E_accum_direct = 0._dp

            call set_timer(time_direct)
            do i = 1, size(nJs, 2)
                E_accum_direct = E_accum_direct + sltcnd_excit(nJs(:, i), Excite_0_t())
            end do
            call halt_timer(time_direct)
        end block

        associate(E_0 => sltcnd_excit(nI, Excite_0_t()))
        block
            real(dp) :: E_accum_from_exc
            E_accum_from_exc = 0._dp

            call set_timer(time_from_exc)
            associate(E_0 => sltcnd_excit(nI, Excite_0_t()))
                do i = 1, size(exc)
                    E_accum_from_exc = E_accum_from_exc + diagH_after_exc(nI, E_0, exc(i))
                end do
            end associate
            call halt_timer(time_from_exc)
        end block
        end associate

        write(*, *) get_total_time(time_direct) / get_total_time(time_from_exc)


        write(*, *) get_gamma(1, 2)
        write(*, *) get_gamma(1, 2)


        call finalize_excitgen_test()

    end subroutine


    subroutine test_sltcnd_driver()
        call run_test_case(test_diagH_from_exc, "test_diagH_from_exc")
        call run_test_case(test_timing, "test_timing")
    end subroutine

end module test_sltcnd_mod

program test_sltcnd

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count
    use util_mod, only: stop_all
    use test_sltcnd_mod, only: test_sltcnd_driver
    use Parallel_neci, only: MPIInit, MPIEnd



    implicit none
    integer :: failed_count
    logical :: err

    call MPIInit(err)

    call init_fruit()

    call test_sltcnd_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_sltcnd', 'failed_tests')

    call MPIEnd(err)

end program
