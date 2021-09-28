! replace _template by whatever one needs
program test_cc_amplitudes

    use fruit
    use cc_amplitudes
    use systemdata, only: nel, nbasis
    use FciMCData, only: projedet, ilutref
    use detbitops, only: encodebitdet
    use bit_reps, only: init_bit_rep

    implicit none

    integer :: failed_count

    call init_fruit()
    call cc_amplitudes_test_driver
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

    subroutine cc_amplitudes_test_driver

        call run_test_case(get_ind_test, "get_ind_test")
        call run_test_case(binomial_test, "binomial_test")
        call run_test_case(calc_number_of_excitations_test, "calc_number_of_excitations_test")
        call run_test_case(calc_n_parallel_excitations_test, "calc_n_parallel_excitations_test")
        call run_test_case(get_ex_test, "get_ex_test")

    end subroutine cc_amplitudes_test_driver

    subroutine get_ex_test


        type(cc_amplitude) :: cc(2)

        integer :: test(2,1), test2(2,2)

        cc(1)%order = 1
        cc(2)%order = 2

        print *, ""
        print *, "testing: get_ex "

        nel = 2
        nbasis = 4

        test(1,1) = 1; test(2,1) = 3
        call assert_equals(test, cc(1)%get_ex(1), 2, 1)
        test(2,1) = 4
        call assert_equals(test, cc(1)%get_ex(2), 2, 1)
        test(1,1) = 2
        test(2,1) = 3
        call assert_equals(test, cc(1)%get_ex(3), 2, 1)
        test(2,1) = 4
        call assert_equals(test, cc(1)%get_ex(4), 2, 1)

        test2(1,:) = [1,2]
        test2(2,:) = [3,4]

        nel = 4
        nbasis = 8

        test(1,1) = 1
        test(2,1) = 5
        call assert_equals(test, cc(1)%get_ex(1), 2, 1)
        test(2,1) = 6
        call assert_equals(test, cc(1)%get_ex(2), 2, 1)
        test(2,1) = 7
        call assert_equals(test, cc(1)%get_ex(3), 2, 1)
        test(2,1) = 8
        call assert_equals(test, cc(1)%get_ex(4), 2, 1)

        test(1,1) = 2
        test(2,1) = 5
        call assert_equals(test, cc(1)%get_ex(5), 2, 1)
        test(2,1) = 6
        call assert_equals(test, cc(1)%get_ex(6), 2, 1)
        test(2,1) = 7
        call assert_equals(test, cc(1)%get_ex(7), 2, 1)
        test(2,1) = 8
        call assert_equals(test, cc(1)%get_ex(8), 2, 1)

        test(1,1) = 4
        test(2,1) = 7
        call assert_equals(test, cc(1)%get_ex(nel*(nbasis-nel)-1), 2, 1)
        test(2,1) = 8
        call assert_equals(test, cc(1)%get_ex(nel*(nbasis-nel)), 2, 1)

        test2(2,:) = [5,6]

        test2(1,:) = [1,4]

        nel = 6
        test(1,1) = 1
        test(2,1) = 7

        call assert_equals(test, cc(1)%get_ex(1), 2, 1)
        test(2,1) = 8
        call assert_equals(test, cc(1)%get_ex(2), 2, 1)
        test(1,1) = 2
        test(2,1) = 7
        call assert_equals(test, cc(1)%get_ex(3), 2, 1)

        test2(1,:) = [1,2]
        test2(2,:) = [7,8]

        nel = 2

        test(1,1) = 1
        test(2,1) = 3
        call assert_equals(test, cc(1)%get_ex(1), 2, 1)
        test(2,1) = 4
        call assert_equals(test, cc(1)%get_ex(2), 2, 1)
        test(2,1) = 5
        call assert_equals(test, cc(1)%get_ex(3), 2, 1)
        test(2,1) = 6
        call assert_equals(test, cc(1)%get_ex(4), 2, 1)
        test(2,1) = 7
        call assert_equals(test, cc(1)%get_ex(5), 2, 1)
        test(2,1) = 8
        call assert_equals(test, cc(1)%get_ex(6), 2, 1)
        test(1,1) = 2
        test(2,1) = 3
        call assert_equals(test, cc(1)%get_ex(7), 2, 1)
        test(2,1) = 4
        call assert_equals(test, cc(1)%get_ex(8), 2, 1)

        nel = -1
        nbasis = -1

    end subroutine get_ex_test

    subroutine get_ind_test

        type(cc_amplitude) :: cc(2)

        print *, ""
        print *, "testing: get_ind"



        cc(1)%order = 1
        cc(2)%order = 2

        ! assume a closed shell reference with the first (nel) spin-orbitals
        ! occupied.. or implement the function in such a way that it always
        ! behaves like it would be like that!
        nel = 1
        nbasis = 2
        call init_bit_rep()

        allocate(projedet(1,1), source = 1)
        allocate(ilutRef(0:niftot,1), source = 0_n_int)
        call encodebitdet([1], ilutRef(:,1))

        call setup_ind_matrix_singles()

        call assert_equals(0, cc(1)%get_ind([1],[2]))
        call assert_equals(0, cc(1)%get_ind([1],[1]))
        call assert_equals(0, cc(1)%get_ind([2],[1]))
        call assert_equals(0, cc(1)%get_ind([2],[2]))

        nel = 2
        nbasis = 4

        call encodebitdet([1,2], ilutRef(:,1))

        call setup_ind_matrix_singles()

        call assert_equals(1, cc(1)%get_ind([1],[3]))
        call assert_equals(0, cc(1)%get_ind([1],[4]))
        call assert_equals(0, cc(1)%get_ind([2],[3]))
        call assert_equals(2, cc(1)%get_ind([2],[4]))

        ! and now doubles:

        nel = 4
        nbasis = 8

        call encodebitdet([1,2,3,4], ilutRef(:,1))
        call setup_ind_matrix_singles()

        call assert_equals(0, cc(1)%get_ind([1],[2]))
        call assert_equals(0, cc(1)%get_ind([1],[3]))
        call assert_equals(1, cc(1)%get_ind([1],[5]))
        call assert_equals(0, cc(1)%get_ind([1],[6]))
        call assert_equals(2, cc(1)%get_ind([1],[7]))
        call assert_equals(0, cc(1)%get_ind([1],[8]))

        call assert_equals(0, cc(1)%get_ind([2],[5]))
        call assert_equals(3, cc(1)%get_ind([2],[6]))
        call assert_equals(0, cc(1)%get_ind([2],[7]))
        call assert_equals(4, cc(1)%get_ind([2],[8]))

        call assert_equals(8, cc(1)%get_ind([4],[8]))
        call assert_equals(0, cc(1)%get_ind([4],[7]))
        call assert_equals(6, cc(1)%get_ind([3],[7]))

        ! and the double excitations are a bit more tricky or?
        ! what if elec_ind > orb_ind?? error? or can that happen if
        ! reference is not ordered or closed shell??

        nel = -1
        nbasis = -1
        deallocate(ilutRef)
        deallocate(projedet)

    end subroutine get_ind_test

    subroutine binomial_test
        print *, ""
        print *, "testing: binomial"

        call assert_equals(1, binomial(1,1))
        call assert_equals(1, binomial(2,0))
        call assert_equals(1, binomial(3,0))
        call assert_equals(1, binomial(4,0))
        call assert_equals(1, binomial(2,0))
        call assert_equals(1, binomial(3,0))
        call assert_equals(1, binomial(4,0))

        call assert_equals(2, binomial(2,1))
        call assert_equals(3, binomial(3,1))
        call assert_equals(4, binomial(4,1))

        call assert_equals(10, binomial(5,2))
        call assert_equals(10, binomial(5,3))

        call assert_equals(20, binomial(6,3))
        call assert_equals(35, binomial(7,4))

        call assert_equals(1, binomial(0,0))
        call assert_equals(0, binomial(0,1))
        call assert_equals(0, binomial(0,100))
        call assert_equals(0, binomial(0,2))

    end subroutine binomial_test

    subroutine calc_number_of_excitations_test

        print *, ""
        print *, "calc_number_of_excitations"

        call assert_equals([8,18], calc_number_of_excitations(2,2,2,4),2)
        call assert_equals([8,18,8], calc_number_of_excitations(2,2,3,4),2)
        call assert_equals([8,18,8,1], calc_number_of_excitations(2,2,4,4),2)

        call assert_equals([7,13,3,0], calc_number_of_excitations(1,2,4,4), 2)
        call assert_equals([7,13,3,0], calc_number_of_excitations(2,1,4,4), 2)

        call assert_equals([9,9,1,0], calc_number_of_excitations(0,3,4,6), 2)

        call assert_equals([10,27,12,0], calc_number_of_excitations(1,3,4,5), 2)
        call assert_equals([10,27,12,0], calc_number_of_excitations(3,1,4,5), 2)

        call assert_equals([12,42,36,9,0], calc_number_of_excitations(2,3,5,5), 2)

    end subroutine calc_number_of_excitations_test

    subroutine calc_n_parallel_excitations_test

        print *, ""
        print *, "testing: calc_n_parallel_excitations "

        call assert_equals(4, calc_n_parallel_excitations(2,4,1))
        call assert_equals(6, calc_n_parallel_excitations(2,5,1))
        call assert_equals(2, calc_n_parallel_excitations(2,3,1))
        call assert_equals(1, calc_n_parallel_excitations(1,2,1))
        call assert_equals(3, calc_n_parallel_excitations(3,4,1))

        call assert_equals(1, calc_n_parallel_excitations(2,4,2))
        call assert_equals(3, calc_n_parallel_excitations(2,5,2))

    end subroutine calc_n_parallel_excitations_test

end program test_cc_amplitudes
