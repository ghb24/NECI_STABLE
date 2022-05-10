module test_cases
    use fruit
    use sets_mod, only: is_sorted, disjoint, subset, operator(.U.), operator(.cap.), &
        operator(.complement.), suc => special_union_complement, set, is_set
    implicit none
    private
    public :: test_sets_mod_driver


contains

    subroutine test_is_sorted()
        call assert_true(is_sorted([1, 2, 3]))
        call assert_true(is_sorted([1., 2., 3.]))
        call assert_true(is_sorted([integer::]))
        call assert_false(is_sorted([1., 4., 3.]))
        call assert_false(is_sorted([1., 4., 3.]))
    end subroutine


    subroutine test_create_set()
        call assert_true(all([1, 2, 3] == set([3, 3, 1, 2, 2])))
        call assert_true(all([integer::] == set([integer::])))
        call assert_true(all([1] == set([1, 1, 1])))
        call assert_true(all([1, 2] == set([2, 1, 1, 1])))

        call assert_true(is_set([1, 2, 3]))
        call assert_true(is_set([integer::]))
        call assert_false(is_set([1, 2, 3, 4, 4]))
        call assert_false(is_set([3, 2, 2]))
    end subroutine


    subroutine test_disjoint()
        call assert_true(disjoint([1, 2, 3], [4, 5, 6]))
        call assert_true(disjoint([1, 2, 3], [5, 6]))
        call assert_true(disjoint([1, 2, 3], [integer::]))

        call assert_false(disjoint([1, 2, 3], [1, 2, 3]))
        call assert_false(disjoint([1, 2, 3], [3]))
        call assert_false(disjoint([1, 2, 3], [1, 3]))
    end subroutine

    subroutine test_subset()
        call assert_true(subset([1, 2, 3], [1, 2, 3]))
        call assert_true(subset([1, 2], [1, 2, 3]))
        call assert_true(subset([integer::], [1, 2, 3]))

        call assert_false(subset([1, 2, 3], [1, 2]))
        call assert_false(subset([1, 2, 3], [integer::]))
    end subroutine

    subroutine test_union()
        associate(expected => [1, 2, 3], calculated => [1, 2] .U. [3])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [1, 2, 3], calculated => [1, 2, 3] .U. [integer::])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [1, 2, 3], calculated => [1, 2, 3] .U. [2, 3])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [1, 2, 3, 4], calculated => [1, 2, 3] .U. [2, 3, 4])
            call assert_equals(expected, calculated, size(calculated))
        end associate
    end subroutine

    subroutine test_intersect()
        associate(expected => [1, 2, 3], calculated => [1, 2, 3] .cap. [1, 2, 3])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [integer::], calculated => [1, 2, 3] .cap. [integer::])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [2, 3], calculated => [1, 2, 3] .cap. [2, 3])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [2, 3], calculated => [1, 2, 3] .cap. [2, 3, 4])
            call assert_equals(expected, calculated, size(calculated))
        end associate
    end subroutine

    subroutine test_complement()
        associate(expected => [integer::], calculated => [1, 2, 3] .complement. [1, 2, 3])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [1, 2, 3], calculated => [1, 2, 3] .complement. [integer::])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [1], calculated => [1, 2, 3] .complement. [2, 3])
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [1], calculated => [1, 2, 3] .complement. [2, 3, 4])
            call assert_equals(expected, calculated, size(calculated))
        end associate
    end subroutine

    subroutine test_special_union_complement()
        associate(expected => [1, 2, 4], calculated => suc([1, 2, 3], [4], [3]))
            call assert_equals(expected, calculated, size(calculated))
        end associate
        associate(expected => [0, 2, 4], &
                  calculated => suc([1, 2, 3], [0, 4], [1, 3]))
            call assert_equals(expected, calculated, size(calculated))
        end associate

        associate(expected => [0, 4, 5], &
                  calculated => suc([1, 2, 3], [0, 4, 5], [1, 2, 3]))
            call assert_equals(expected, calculated, size(calculated))
        end associate

        associate(expected => [1, 2, 3], &
                  calculated => suc([1, 2, 3], [integer::], [integer::]))
            call assert_equals(expected, calculated, size(calculated))
        end associate
    end subroutine

    subroutine test_sets_mod_driver()
        call run_test_case(test_is_sorted, "test_is_sorted")
        call run_test_case(test_create_set, "test_create_set")
        call run_test_case(test_disjoint, "test_disjoint")
        call run_test_case(test_subset, "test_subset")
        call run_test_case(test_union, "test_union")
        call run_test_case(test_intersect, "test_intersect")
        call run_test_case(test_complement, "test_complement")
        call run_test_case(test_special_union_complement, "test_special_union_complement")
    end subroutine

end module test_cases

program test_sets_mod

    use mpi
    use fruit
    use test_cases, only: test_sets_mod_driver

    implicit none
    integer :: failed_count, err

    call mpi_init(err)

    call init_fruit()

    call test_sets_mod_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_sets_mod', 'failed_tests')

    call mpi_finalize(err)

contains
end program test_sets_mod
