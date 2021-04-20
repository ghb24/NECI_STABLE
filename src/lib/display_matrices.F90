#include "macros.h"

module display_matrices
    use constants, only: dp, iout
    use util_mod, only: get_free_unit
    implicit none
    private
    public :: write_matrix

    interface write_matrix
        module procedure write_matrix_1D, write_matrix_2D
    end interface

contains

    subroutine write_matrix_2D(M, dec_places, unit_id, advance)
        use fortran_strings, only: str
        implicit none
        real(dp), intent(in) :: M(:, :)
        integer, intent(in), optional :: dec_places, unit_id
        logical, intent(in), optional :: advance
        integer :: unit_id_, i
        logical :: advance_

        def_default(unit_id_, unit_id, iout)
        def_default(advance_, advance, .true.)

        write(unit_id_, '(A)') '['
        do i = 1, size(M, 1) - 1
            write(unit_id_, '(A)', advance='no') ' '
            call write_matrix(M(i, :), dec_places, unit_id_, advance=.false.)
            write(unit_id_, '(A)') ','
        end do
        write(unit_id_, '(A)', advance='no') ' '
        call write_matrix(M(i, :), dec_places, unit_id_, advance=.true.)

        write(unit_id_, '(A)', advance=merge('yes', 'no ', advance_)) ']'
    end subroutine

    subroutine write_matrix_1D(M, dec_places, unit_id, advance)
        use fortran_strings, only: str
        implicit none
        real(dp), intent(in) :: M(:)
        integer, intent(in), optional :: dec_places, unit_id
        logical, intent(in), optional :: advance
        integer :: dec_places_, unit_id_
        logical :: advance_

        character(:), allocatable :: fmter, fmter_no_comma
        integer :: i

        def_default(dec_places_, dec_places, 5)
        def_default(unit_id_, unit_id, iout)
        def_default(advance_, advance, .true.)

        fmter = '(E'//str(dec_places_ + 6)//'.'//str(dec_places_)//', A2)'
        fmter_no_comma = '(E'//str(dec_places_ + 6)//'.'//str(dec_places_)//')'

        write(unit_id_, '(A)', advance='no') '['
        do i = 1, size(M) - 1
            write(unit_id_, fmter, advance='no') M(i), ', '
        end do
        write(unit_id_, fmter_no_comma, advance='no') M(size(M))
        write(unit_id_, '(A)', advance=merge('yes', 'no ', advance_)) ']'
    end subroutine

end module
