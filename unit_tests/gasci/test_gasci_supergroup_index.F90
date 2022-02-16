module test_gasci_supergroup_index_mod
    use fruit
    use constants, only: dp, int64, n_int
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero, choose
    use util_mod, only: factrl, intswap, cumsum

    use gasci, only: LocalGASSpec_t, CumulGASSpec_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t, composition_idx, get_compositions, &
        get_lowest_supergroup, get_highest_supergroup

    implicit none
    private
    public :: test_gasci_driver

contains

    subroutine test_gasci_driver()
        call run_test_case(test_compositioning, "test_compositioning")
        call run_test_case(test_supergroup_indexer_class, "test_supergroup_indexer_class")
        call run_test_case(test_lowest_supergroup, "test_lowest_supergroup")
        call run_test_case(test_highest_supergroup, "test_highest_supergroup")
        call run_test_case(test_next_supergroup, "test_next_supergroup")
        call run_test_case(test_count_supergroups, "test_count_supergroups")
    end subroutine

    subroutine test_compositioning()
        integer, allocatable :: compositions(:, :)
        integer :: i, n, k
        logical :: correct

        correct = .true.
        do n = 1, 10
            do k = 1, 10
                compositions = get_compositions(k, n)
                do i = 1, size(compositions, 2)
                    if (i /= composition_idx(compositions(:, i))) correct = .false.
                end do
            end do
        end do
        call assert_true(correct)
    end subroutine

    subroutine test_supergroup_indexer_class()
        block
            type(LocalGASSpec_t) :: GAS_spec
            type(SuperGroupIndexer_t) :: indexer
            integer :: i, j
            integer, allocatable :: supergroups(:, :)
            logical :: correct

            GAS_spec = LocalGASSpec_t(&
                n_min=[5,  4,  4,  4,  5], &
                n_max=[7,  8,  8,  8,  7], &
                spat_GAS_orbs = [([(j, i = 1, 6)], j = 1, 5)])
            call assert_true(GAS_spec%is_valid())

            indexer = SuperGroupIndexer_t(GAS_spec, 30)
            supergroups = indexer%get_supergroups()

            correct = .true.
            do i = 1, size(supergroups, 2)
                if (i /= indexer%idx_supergroup(supergroups(:, i))) correct = .false.
            end do
            call assert_true(correct)
        end block

        block
            type(LocalGASSpec_t) :: GAS_spec
            type(SuperGroupIndexer_t) :: indexer

            GAS_spec = LocalGASSpec_t(&
                n_min=[0, 0, 1], &
                n_max=[2, 2, 2], &
                spat_GAS_orbs = [1, 1, 2, 2, 3, 3])
            call assert_true(GAS_spec%is_valid())

            indexer = SuperGroupIndexer_t(GAS_spec, 3)

            call assert_true(5 == indexer%n_supergroups())

            call assert_true(1 == indexer%idx_nI([1, 2, 9]))
            call assert_true(1 == indexer%idx_nI([1, 3, 9]))
            call assert_true(1 == indexer%idx_nI([1, 3, 10]))

            call assert_true(2 == indexer%idx_nI([1, 5, 10]))
            call assert_true(2 == indexer%idx_nI([1, 7, 10]))
            call assert_true(2 == indexer%idx_nI([1, 5, 10]))

            call assert_true(3 == indexer%idx_nI([1, 10, 11]))
            call assert_true(3 == indexer%idx_nI([3, 10, 11]))

            call assert_true(4 == indexer%idx_nI([5, 7, 11]))
            call assert_true(4 == indexer%idx_nI([6, 7, 11]))

            call assert_true(5 == indexer%idx_nI([7, 10, 11]))
        end block
    end subroutine


    subroutine test_count_supergroups()
        type(SuperGroupIndexer_t), allocatable :: sg_indexer(:)
        integer, allocatable :: expected(:)
        integer :: i

        sg_indexer = [(get_sg_indexer(i), i = 1, 5)]
        expected = [1, 5, 25, 125, 625]

        do i = 1, size(sg_indexer)
            call assert_equals(expected(i), sg_indexer(i)%n_supergroups())
        end do

    contains

        elemental function get_benzene_GAS_spec(n_benz, n_exc) result(res)
            !! Create the GAS specification for [n_benz * (6, 6)] active spaces.
            integer, intent(in) :: n_benz, n_exc
            type(CumulGASSpec_t) :: res
            integer :: i, j
            res = CumulGASSpec_t(&
                cn_min=[(6 * i - n_exc, i = 1, n_benz - 1), n_benz * 6], &
                cn_max=[(6 * i + n_exc, i = 1, n_benz - 1), n_benz * 6], &
                spat_GAS_orbs = [([(j, i = 1, 6)], j = 1, n_benz)])
        end function

        function get_sg_indexer(n_benz) result(res)
            integer, intent(in) :: n_benz
            integer, parameter :: n_exc = 2
            type(SuperGroupIndexer_t) :: res
            res = SuperGroupIndexer_t(get_benzene_GAS_spec(n_benz, n_exc), n_benz * 6)
        end function
    end subroutine

    subroutine test_lowest_supergroup
        integer, allocatable :: calculated(:), expected(:)
        integer :: N
        block
            type(LocalGASSpec_t) :: GAS_spec

            GAS_spec = LocalGASSpec_t(n_min=[0, 1, 2], n_max=[1, 2, 3], &
                                      spat_GAS_orbs=[1, 1, 1, 2, 2, 2, 3, 3, 3])
            N = 3
            expected = [0, 1, 2]
            calculated = get_lowest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))

            N = 5
            expected = [1, 2, 2]
            calculated = get_lowest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))


            GAS_spec = LocalGASSpec_t(n_min=[0, 1, 1], n_max=[3, 2, 1], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 5
            expected = [2, 2, 1]
            calculated = get_lowest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))

            GAS_spec = LocalGASSpec_t(n_min=[0, 1, 2], n_max=[3, 2, 2], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 5
            expected = [2, 1, 2]
            calculated = get_lowest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))
        end block

        block
            type(CumulGASSpec_t) :: GAS_spec

            GAS_spec = CumulGASSpec_t(cn_min=[0, 1, 4], cn_max=[1, 2, 4], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 4
            expected = [1, 1, 2]
            calculated = get_lowest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))


            GAS_spec = CumulGASSpec_t(cn_min=[1, 2, 5], cn_max=[4, 4, 5], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 5
            expected = [2, 2, 1]
            calculated = get_lowest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))
        end block
    end subroutine

    subroutine test_highest_supergroup
        integer, allocatable :: calculated(:), expected(:)
        integer :: N
        block
            type(LocalGASSpec_t) :: GAS_spec

            GAS_spec = LocalGASSpec_t(n_min=[0, 1, 2], n_max=[1, 2, 3], &
                                      spat_GAS_orbs=[1, 1, 1, 2, 2, 2, 3, 3, 3])
            N = 3
            expected = [0, 1, 2]
            calculated = get_highest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))

            N = 5
            expected = [0, 2, 3]
            calculated = get_highest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))


            GAS_spec = LocalGASSpec_t(n_min=[0, 1, 1], n_max=[3, 2, 1], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 5
            expected = [2, 2, 1]
            calculated = get_highest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))

            GAS_spec = LocalGASSpec_t(n_min=[0, 1, 2], n_max=[3, 2, 2], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 5
            expected = [1, 2, 2]
            calculated = get_highest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))
        end block

        block
            type(CumulGASSpec_t) :: GAS_spec

            GAS_spec = CumulGASSpec_t(cn_min=[0, 1, 1], &
                                      cn_max=[2, 3, 3], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 3
            expected = [0, 1, 2]
            calculated = get_highest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))


            GAS_spec = CumulGASSpec_t(cn_min=[1, 2, 5], cn_max=[4, 4, 5], &
                                      spat_GAS_orbs=[1, 2, 3])
            N = 5
            expected = [1, 2, 2]
            calculated = get_highest_supergroup(GAS_spec, N)
            call assert_equals(size(calculated), size(expected))
            call assert_equals(calculated, expected, size(calculated))
        end block


    end subroutine

    subroutine test_next_supergroup

    end subroutine

end module


program test_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_supergroup_index_mod, only: test_gasci_driver


    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        block
            integer, allocatable :: M(:)

            M = [1, 2, 3]
            write(*, *) M
            write(*, *) eoshift(M, shift=-1)
            write(*, *) eoshift(M, shift=+1)

        end block

        call init_fruit()

        call test_gasci_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

end program test_gasci_program
