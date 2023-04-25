module test_sltcnd_mod
    use fruit, only: assert_true, assert_false, run_test_case
    use constants, only: dp, n_int

    use bit_reps, only: decode_bit_det
    use unit_test_helper_excitgen, only: &
        init_excitgen_test, finalize_excitgen_test, &
        generate_random_integrals, RandomFcidumpWriter_t
    use symexcit3, only: gen_excits
    use excitation_types, only: get_excitation, Excitation_t, Excite_0_t, Excite_1_t, Excite_2_t
    use sltcnd_mod, only: sltcnd_excit, diagH_after_exc
    use util_mod, only: operator(.isclose.)
    implicit none
    private
    public :: test_sltcnd_driver


contains

    subroutine test_diagH_from_exc()
        integer, parameter :: nI(4) = [1, 4, 7, 10], n_spat_orb = 10

        integer(n_int), allocatable :: det_list(:, :)
        class(Excitation_t), allocatable :: exc
        integer :: nJ(size(nI)), n_excits, i, ic
        logical :: tParity

        call init_excitgen_test(nI, &
            RandomFcidumpWriter_t(&
                n_spat_orb, nI, sparse=0.7_dp, sparseT=0.7_dp) &
        )

        associate(E_0 => sltcnd_excit(nI, Excite_0_t()))
            do ic = 1, 1
                call gen_excits(nI, n_excits, det_list, ic)
                do i = 1, size(det_list, 2)
                    call decode_bit_det(nJ, det_list(:, i))
                    call get_excitation(nI, nJ, ic, exc, tParity)

                    select type (exc)
                    type is (Excite_1_t)
                        if (.not. (sltcnd_excit(nJ, Excite_0_t()) .isclose. diagH_after_exc(E_0, exc))) then
                            write(*, *)
                            write(*, *) nI
                            write(*, *) exc%val(1, :)
                            write(*, *) exc%val(2, :)
                            write(*, *) nJ
                            write(*, *) sltcnd_excit(nJ, Excite_0_t()), diagH_after_exc(E_0, exc)
                        end if

                    type is (Excite_2_t)
                        write(*, *) exc%val(1, :)
                        write(*, *) exc%val(2, :)
                        write(*, *) sltcnd_excit(nJ, Excite_0_t()), diagH_after_exc(E_0, exc)
                    end select
                end do
            end do

        end associate



        call finalize_excitgen_test()
    end subroutine


    subroutine test_sltcnd_driver()
        call run_test_case(test_diagH_from_exc, "test_diagH_from_exc")
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
