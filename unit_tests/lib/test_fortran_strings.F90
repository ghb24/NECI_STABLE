module test_fortran_strings_mod
    use fruit
    use constants, only: dp
    use fortran_strings, only: split, Token_t, can_be_real, can_be_int, str, &
        to_upper, to_lower
    implicit none
    private
    public :: test_driver



contains

    subroutine test_split()

        block
            type(Token_t), allocatable :: tokens(:)
            character(:), allocatable :: physicists


            physicists = "Curie Schroedinger Newton"
            tokens = split(physicists, ' ')
            call assert_true(tokens(1)%str == 'Curie')
            call assert_true(tokens(2)%str == 'Schroedinger')
            call assert_true(tokens(3)%str == 'Newton')


            physicists = "  Curie    Schroedinger   Newton  "
            tokens = split(physicists, ' ')
            call assert_true(tokens(1)%str == 'Curie')
            call assert_true(tokens(2)%str == 'Schroedinger')
            call assert_true(tokens(3)%str == 'Newton')
        end block

        block
            type(Token_t), allocatable :: tokens(:)
            character(:), allocatable :: expr

            expr = "5*3"
            tokens = split(expr, '*')
            call assert_true(tokens(1)%str == '5')
            call assert_true(tokens(2)%str == '3')

            expr = "*5**3*"
            tokens = split(expr, '*')
            call assert_true(tokens(1)%str == '5')
            call assert_true(tokens(2)%str == '3')
        end block
    end subroutine

    subroutine test_if_numbers()

        call assert_true(can_be_real("5"))
        call assert_true(can_be_int("5"))

        call assert_true(can_be_real("5."))
        call assert_false(can_be_int("5.2"))

        call assert_false(can_be_real("asdf"))
        call assert_false(can_be_int("asdf"))
        call assert_false(can_be_real(""))
        call assert_false(can_be_int(""))
        call assert_false(can_be_real("5..3"))
        call assert_false(can_be_int("5..3"))

    end subroutine


    subroutine test_conversion()

        call assert_true(str(5) == ("5"))

        call assert_true(str(5., 1) == "0.5E+01")

        call assert_true(str(5.312e8, 1) == "0.5E+09")
        call assert_true(str(5.312e8, 2) == "0.53E+09")
        call assert_true(str(5.312e8, 3) == "0.531E+09")

    end subroutine

    subroutine test_driver()
        call run_test_case(test_split, "test_split")
        call run_test_case(test_if_numbers, "test_if_numbers")
        call run_test_case(test_conversion, "test_conversion")
    end subroutine

end module test_fortran_strings_mod

program test_fortran_strings_program

    use mpi
    use fruit
    use test_fortran_strings_mod, only: test_driver

    implicit none
    integer :: failed_count

    block

        call init_fruit()

        call test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_fortran_strings_program', 'failed_tests')

    end block

contains

end program
