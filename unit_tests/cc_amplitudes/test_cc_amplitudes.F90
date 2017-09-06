! replace _template by whatever one needs 
program test_cc_amplitudes

    use fruit
    use cc_amplitudes

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

        call run_test_case(dongxia_amplitudes_test, "dongxia_amplitudes_test")
        call run_test_case(get_ind_test, "get_ind_test")
        call run_test_case(get_ex_test, "get_ex_test")

    end subroutine cc_amplitudes_test_driver

    subroutine get_ex_test

        use systemdata, only: nel, nbasis 

        type(cc_amplitude) :: cc(2)

        integer :: test(2,1), test2(2,2)

        cc(1)%order = 1 
        cc(2)%order = 2 

        print *, ""
        print *, "testing: get_ex "

        nel = 1 
        nbasis = 2 
        test(1,1) = 1 
        test(2,1) = 2

        print *, ""
        print *, "nel, nbasis: ", nel, nbasis
        call assert_equals(test, cc(1)%get_ex(1), 2, 1)

        nel = 2
        nbasis = 4 

        print *, ""
        print *, "nel, nbasis: ", nel, nbasis

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

!         call assert_equals(test2, cc(2)%get_ex(1), 2, 2)

        nel = 4
        nbasis = 8 

        print *, ""
        print *, "nel, nbasis: ", nel, nbasis
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

        test2(1,:) = [1,2]
        test2(2,:) = [5,6]
! 
!         call assert_equals(test2, cc(2)%get_ex(1), 2, 2)
!         test2(2,2) = 7
!         call assert_equals(test2, cc(2)%get_ex(2), 2, 2)
!         test2(2,2) = 8
!         call assert_equals(test2, cc(2)%get_ex(3), 2, 2)
!         test2(2,:) = [6,7]
!         call assert_equals(test2, cc(2)%get_ex(4), 2, 2)
!         test2(2,2) = 8
!         call assert_equals(test2, cc(2)%get_ex(5), 2, 2)
!         test2(2,1) = 7
!         call assert_equals(test2, cc(2)%get_ex(6), 2, 2)
! 
!         test2(1,2) = 3
!         test2(2,:) = [5,6]
!         call assert_equals(test2, cc(2)%get_ex(7), 2, 2)

        test2(1,:) = [1,4]
!         test2(2,:) = [5,6]
!         call assert_equals(test2, cc(2)%get_ex(13), 2, 2)
!         test2(1,:) = [2,3]
!         test2(2,:) = [5,6]
!         call assert_equals(test2, cc(2)%get_ex(19), 2, 2)
!         test2(1,:) = [3,4]
!         test2(2,:) = [7,8]
!         call assert_equals(test2, cc(2)%get_ex((nel*(nel-1))**2/4), 2, 2)

        nel = 6 
        test(1,1) = 1
        test(2,1) = 7

        print *, ""
        print *, "nel, nbasis: ", nel, nbasis
        call assert_equals(test, cc(1)%get_ex(1), 2, 1)
        test(2,1) = 8
        call assert_equals(test, cc(1)%get_ex(2), 2, 1)
        test(1,1) = 2
        test(2,1) = 7
        call assert_equals(test, cc(1)%get_ex(3), 2, 1)

        test2(1,:) = [1,2]
        test2(2,:) = [7,8]

!         call assert_equals(test2, cc(2)%get_ex(1), 2, 2)
!         test2(1,:) = [1,3]
! 
!         call assert_equals(test2, cc(2)%get_ex(2), 2, 2)
!         print *, cc(2)%get_ex(2)

        nel = 2 
        print *, ""
        print *, "nel, nbasis: ", nel, nbasis

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
        use systemdata, only: nel, nbasis

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
        call assert_equals(1, cc(1)%get_ind([1],[2]))

        nel = 2
        nbasis = 4
        
        call assert_equals(1, cc(1)%get_ind([1],[3]))
        call assert_equals(2, cc(1)%get_ind([1],[4]))
        call assert_equals(3, cc(1)%get_ind([2],[3]))
        call assert_equals(4, cc(1)%get_ind([2],[4]))

        ! and now doubles: 
!         call assert_equals(1, cc(2)%get_ind([1,2],[3,4]))

        nel = 4 
        nbasis = 8 
        call assert_equals(1, cc(1)%get_ind([1],[5]))
        call assert_equals(2, cc(1)%get_ind([1],[6]))
        call assert_equals(3, cc(1)%get_ind([1],[7]))
        call assert_equals(4, cc(1)%get_ind([1],[8]))

        call assert_equals(5, cc(1)%get_ind([2],[5]))
        call assert_equals(6, cc(1)%get_ind([2],[6]))
        call assert_equals(7, cc(1)%get_ind([2],[7]))
        call assert_equals(8, cc(1)%get_ind([2],[8]))

        call assert_equals(nel*(nbasis-nel), cc(1)%get_ind([4],[8]))
        call assert_equals(nel*(nbasis-nel)-1, cc(1)%get_ind([4],[7]))

        ! and now doubles: 
!         call assert_equals(1, cc(2)%get_ind([1,2],[5,6]))
!         call assert_equals(2, cc(2)%get_ind([1,2],[5,7]))
!         call assert_equals(3, cc(2)%get_ind([1,2],[5,8]))
!         call assert_equals(4, cc(2)%get_ind([1,2],[6,7]))
!         call assert_equals(5, cc(2)%get_ind([1,2],[6,8]))
!         call assert_equals(6, cc(2)%get_ind([1,2],[7,8]))
! 
!         call assert_equals(7, cc(2)%get_ind([1,3],[5,6]))
!         call assert_equals(8, cc(2)%get_ind([1,3],[5,7]))
!         call assert_equals(9, cc(2)%get_ind([1,3],[5,8]))
!         call assert_equals(10, cc(2)%get_ind([1,3],[6,7]))
! 
!         call assert_equals(13, cc(2)%get_ind([1,4],[5,6]))
!         call assert_equals(19, cc(2)%get_ind([2,3],[5,6]))

!         call assert_equals((nel*(nel-1))**2/4, cc(2)%get_ind([3,4],[7,8]))
        ! and the double excitations are a bit more tricky or? 
        ! what if elec_ind > orb_ind?? error? or can that happen if 
        ! reference is not ordered or closed shell??

        nel = 6 
        call assert_equals(1, cc(1)%get_ind([1],[7]))
        call assert_equals(2, cc(1)%get_ind([1],[8]))
        call assert_equals(3, cc(1)%get_ind([2],[7]))
        call assert_equals(4, cc(1)%get_ind([2],[8]))
        
!         call assert_equals(1, cc(2)%get_ind([1,2],[7,8]))
!         call assert_equals(2, cc(2)%get_ind([1,3],[7,8]))
!         call assert_equals(3, cc(2)%get_ind([1,4],[7,8]))
! 
!         call assert_equals(6, cc(2)%get_ind([1,6],[7,8]))
! 
!         call assert_equals(7, cc(2)%get_ind([2,3],[7,8]))
!         call assert_equals(8, cc(2)%get_ind([2,4],[7,8]))

        nel = 2 
        call assert_equals(1, cc(1)%get_ind([1],[3]))
        call assert_equals(2, cc(1)%get_ind([1],[4]))
        call assert_equals(6, cc(1)%get_ind([1],[8]))
        call assert_equals(7, cc(1)%get_ind([2],[3]))

!         call assert_equals(1, cc(2)%get_ind([1,2],[3,4]))
!         call assert_equals(2, cc(2)%get_ind([1,2],[3,5]))
!         call assert_equals(3, cc(2)%get_ind([1,2],[3,6]))
!         call assert_equals(4, cc(2)%get_ind([1,2],[3,7]))
!         call assert_equals(5, cc(2)%get_ind([1,2],[3,8]))
! 
!         call assert_equals(6, cc(2)%get_ind([1,2],[4,5]))

        nel = -1 
        nbasis = -1 

    end subroutine get_ind_test

    subroutine dongxia_amplitudes_test

        print *, ""
        print *, "test dongxia_amplitudes "

        call assert_true(.false.)

    end subroutine dongxia_amplitudes_test

end program test_cc_amplitudes
