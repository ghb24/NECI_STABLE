module fruit_extensions

    use fruit, only: run_test_case
    implicit none
    private
    public :: my_run_test_case

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

end module

