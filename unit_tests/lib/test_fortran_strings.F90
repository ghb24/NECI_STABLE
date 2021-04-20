module test_fortran_strings_mod
    use fruit
    use constants, only: dp
    use fortran_strings, only: split, Token_t
    implicit none
    private
    public :: test_split



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

end module test_fortran_strings_mod

program test_fortran_strings_program

    use mpi
    use fruit
    use test_fortran_strings_mod, only: test_split

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

    subroutine test_driver()
        call run_test_case(test_split, "test_split")
    end subroutine
end program
