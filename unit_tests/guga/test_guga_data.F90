module test_guga_data_mod
    use fruit, only: assert_true, assert_false, assert_equals, run_test_case
    use guga_data, only: ijkl_Index_t, ijab_Index_t, &
              DistinctDoubleType_t, distinct_doubles, DistinctDouble_t, &
              excit_type
    implicit none
    private
    public :: test_driver

contains

    subroutine test_driver()
        call run_test_case(test_indices, "test_indices")
        call run_test_case(test_distinct_double_type, "test_distinct_double_type")
        call run_test_case(test_distinct_double, "test_distinct_double_type")
        ! call run_test_case(test_fuse, "test_fuse")
    end subroutine

    subroutine test_indices()
        !! Here we test the explicit ijkl vs. ijab indexing
        !!
        !! This includes also conversion between the two and
        !! normalization.
        type(ijkl_Index_t) :: ijkl_idx

        ijkl_idx = ijkl_Index_t([3, 4, 1, 2])
        call assert_false(ijkl_idx%is_normalized())
        block
            integer, parameter :: expected(4) = [3, 4, 1, 2]
            integer :: i, j, k, l
            call assert_true(all(ijkl_idx%val == expected))

            call ijkl_idx%extract_ijkl(i, j, k, l)
            call assert_true(all([i, j, k, l] == expected))
        end block


        associate(ijab_idx => ijkl_idx%convert())
            block
                integer, parameter :: expected(4) = [4, 2, 3, 1]
                integer :: i, j, a, b
                call assert_true(all(ijab_idx%val == expected))

                call ijab_idx%extract_ijab(i, j, a, b)
                call assert_true(all([i, j, a, b] == expected))
            end block
        end associate

        ijkl_idx = ijkl_idx%normalize()
        call assert_true(ijkl_idx%is_normalized())
        block
            integer, parameter :: expected(4) = [1, 2, 3, 4]
            integer :: i, j, k, l
            call assert_true(all(ijkl_idx%val == expected))

            call ijkl_idx%extract_ijkl(i, j, k, l)
            call assert_true(all([i, j, k, l] == expected))
        end block


        associate(ijab_idx => ijkl_idx%convert())
            block
                integer, parameter :: expected(4) = [2, 4, 1, 3]
                integer :: i, j, a, b
                call assert_true(all(ijab_idx%val == expected))

                call ijab_idx%extract_ijab(i, j, a, b)
                call assert_true(all([i, j, a, b] == expected))
            end block
        end associate
    end subroutine

    subroutine test_distinct_double_type
        !! Here we test if the type of distinct double excitation is correctly found
        associate(dd => distinct_doubles)
            associate(typ => DistinctDoubleType_t(ijkl_Index_t([2, 1, 2, 3])))
                call assert_true(typ == dd%jijk)
                call assert_true(typ%repr == excit_type%single_overlap_L_to_R)
                call assert_true(typ%other == excit_type%single_overlap_L_to_R)
            end associate
            associate(typ => DistinctDoubleType_t(ijkl_Index_t([1, 2, 3, 2])))
                call assert_true(typ == dd%ijkj)
                call assert_true(typ%repr == excit_type%single_overlap_R_to_L)
                call assert_true(typ%other == excit_type%single_overlap_R_to_L)
            end associate
            associate(typ => DistinctDoubleType_t(ijkl_Index_t([1, 3, 2, 3])))
                call assert_true(typ == dd%ikjk)
                call assert_true(typ%repr == excit_type%fullstop_raising)
                call assert_true(typ%other == excit_type%fullstop_raising)
            end associate
            associate(typ => DistinctDoubleType_t(ijkl_Index_t([3, 1, 3, 2])))
                call assert_true(typ == dd%kikj)
                call assert_true(typ%repr == excit_type%fullstop_lowering)
                call assert_true(typ%other == excit_type%fullstop_lowering)
            end associate
            associate(typ => DistinctDoubleType_t(ijkl_Index_t([2, 1, 3, 3])))
                call assert_true(typ == dd%jikk)
                call assert_true(typ%repr == excit_type%fullstop_L_to_R)
                call assert_true(typ%other == excit_type%non_overlap)
            end associate
            associate(typ => DistinctDoubleType_t(ijkl_Index_t([1, 2, 3, 4])))
                call assert_true(typ == dd%ijkl)
                call assert_true(typ%repr == excit_type%double_R_to_L_to_R)
                call assert_true(typ%other == excit_type%non_overlap)
            end associate
        end associate
    end subroutine


    subroutine test_distinct_double
        !! Here we test if the type of distinct double excitation is correctly found
        !! and the corresponding index for representation is correctly obtained.
        associate(dd => distinct_doubles)
            associate(exc => DistinctDouble_t(ijkl_Index_t([2, 1, 2, 3])))
                call assert_true(exc%typ == dd%jijk)
                call assert_true(exc%typ%repr == excit_type%single_overlap_L_to_R)
                call assert_true(exc%typ%other == excit_type%single_overlap_L_to_R)
                associate(ijkl_idx => exc%get_repr_index())
                    ! The excitations without non-overlap are boring
                    call assert_true(all(ijkl_idx%val == [2, 1, 2, 3]))
                end associate
            end associate
            associate(exc => DistinctDouble_t(ijkl_Index_t([1, 2, 3, 2])))
                call assert_true(exc%typ == dd%ijkj)
                call assert_true(exc%typ%repr == excit_type%single_overlap_R_to_L)
                call assert_true(exc%typ%other == excit_type%single_overlap_R_to_L)
            end associate
            associate(exc => DistinctDouble_t(ijkl_Index_t([1, 3, 2, 3])))
                call assert_true(exc%typ == dd%ikjk)
                call assert_true(exc%typ%repr == excit_type%fullstop_raising)
                call assert_true(exc%typ%other == excit_type%fullstop_raising)
            end associate
            associate(exc => DistinctDouble_t(ijkl_Index_t([3, 1, 3, 2])))
                call assert_true(exc%typ == dd%kikj)
                call assert_true(exc%typ%repr == excit_type%fullstop_lowering)
                call assert_true(exc%typ%other == excit_type%fullstop_lowering)
            end associate
            associate(exc => DistinctDouble_t(ijkl_Index_t([2, 1, 3, 3])))
                call assert_true(exc%typ == dd%jikk)
                call assert_true(exc%typ%repr == excit_type%fullstop_L_to_R)
                call assert_true(exc%typ%other == excit_type%non_overlap)
                associate(ijkl_idx => exc%get_repr_index())
                    call assert_true(all(ijkl_idx%val == [2, 3, 3, 1]))
                end associate
            end associate
            associate(exc => DistinctDouble_t(ijkl_Index_t([1, 20, 300, 4000])))
                call assert_true(exc%typ == dd%ijkl)
                call assert_true(exc%typ%repr == excit_type%double_R_to_L_to_R)
                call assert_true(exc%typ%other == excit_type%non_overlap)
                associate(ijkl_idx => exc%get_repr_index())
                    call assert_true(all(exc%typ%idx_for_repr == [1, 4, 3, 2]))
                    call assert_true(all(ijkl_idx%val == [1, 4000, 300, 20]))
                end associate
            end associate
        end associate
    end subroutine

end module

program test_guga_data_prog

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_guga_data_mod, only: test_driver
    use util_mod, only: stop_all

    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_guga_data', 'failed_tests')

        call MPIEnd(.false.)
    end block

end program
