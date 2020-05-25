! replace _template by whatever one needs
program test_guga_pchb_excitgen

    use fruit
    ! use template

    implicit none

    integer :: failed_count

    call init_fruit()
    call guga_pchb_test_driver()
    call fruit_summary()
    call fruit_finalize()

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

    subroutine my_run_test_case( routine, testing_name, tested_name)
        interface
            subroutine routine()
            end subroutine
        end interface
        character(*), intent(in) :: testing_name, tested_name

        print *, ""
        print *, " ========================================================= "
        print *, "  testing: ", tested_name,"()"
        print *, " ========================================================= "


        call run_test_case(routine, testing_name)

        print *, ""
        print *, " ========================================================= "
        print *, "  testing: ", tested_name, "()  DONE! "
        print *, " ========================================================= "
        print *, ""

    end subroutine my_run_test_case

    subroutine guga_pchb_test_driver()

        call my_run_test_case(pick_uniform_spatial_hole_test, &
            "pick_uniform_spatial_hole_test", "pick_uniform_spatial_hole")

    end subroutine guga_pchb_test_driver

    subroutine pick_uniform_spatial_hole_test

        call assert_true(.false.)

    end subroutine pick_uniform_spatial_hole_test

end program test_guga_pchb_excitgen

