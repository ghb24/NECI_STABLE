#include "macros.h"

module display_matrices
    use constants, only: dp, iout
    use util_mod, only: get_free_unit
    implicit none
    private
    public :: write_matrix

contains

    subroutine write_matrix(M, dec_places, unit_id)
        use fortran_strings, only: str
        implicit none
        real(dp), intent(in) :: M(:, :)
        integer, intent(in), optional :: dec_places, unit_id
        integer :: dec_places_, unit_id_

        character(:), allocatable :: fmter, fmter_no_comma
        integer :: i, j

        def_default(dec_places_, dec_places, 5)
        def_default(unit_id_, unit_id, iout)

        fmter = '(E'//str(dec_places_ + 6)//'.'//str(dec_places_)//', A2)'
        fmter_no_comma = '(E'//str(dec_places_ + 6)//'.'//str(dec_places_)//')'

        do i = 1, size(M, 1)
            do j = 1, size(M, 2) - 1
                write(unit_id_, fmter, advance='no') M(i, j), ', '
            end do
            write(unit_id_, fmter_no_comma, advance='no') M(i, size(M, 2))
            write(unit_id_, *)
        end do
    end subroutine write_matrix

end module
