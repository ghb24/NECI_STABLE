#include "macros.h"

module DeterminantData
    use bit_rep_data, only: NIfTot
    use constants
    use SystemData, only: nel, nbasis
    use MemoryManager, only: TagIntType
    implicit none
    private
    public :: write_det, write_det_len, tagFDet, Fdet, calculated_ms, &
        get_lexicographic

    integer, pointer :: FDet(:)
    integer(TagIntType) :: tagFDet

    integer :: calculated_ms

    type lexicographic_store
        integer, pointer :: dorder(:) => null()
        integer, pointer :: open_orbs(:) => null()
        integer, pointer :: open_indices(:) => null()
        integer :: nopen, nup
    end type

contains

    subroutine write_det(nunit, nI, lTerm)

        ! Write the specified determinant (or CSF) to the output unit nunit.
        ! Terminate the line (with '\n') if lTerm = .true.
        !
        ! In: nunit - Output (file) unit
        !     nI    - Determinant to output
        !     lTerm - Terminate line with newline character?

        integer, intent(in) :: nunit, nI(nel)
        logical, intent(in) :: lTerm

        call write_det_len(nunit, nI, nel, lterm)
    end subroutine write_det

    subroutine write_det_len(nunit, nI, nlen, lterm)

        ! Worker function for the above. Can be accessed to print an unusual
        ! lengthed determinant.

        integer, intent(in) :: nunit, nlen, nI(nlen)
        logical, intent(in) :: lTerm
        integer :: i, elec

        ! Start with a bracket, and loop over all the electrons
        write(nunit, '("(")', advance='no')
        do i = 1, nlen
            elec = nI(i)

            ! Write out the orbital number
            write(nunit, '(i5)', advance='no') elec
            if (i /= nlen) write(nunit, '(",")', advance='no')
        end do

        ! Close the written determinant off
        write(nunit, '(")")', advance='no')
        if (lTerm) write(nunit, *)
    end subroutine write_det_len

    subroutine get_lexicographic(dorder, nopen, nup)

        ! Unlike the csf version, this uses 1 == alpha, 0 = beta.

        integer, intent(in) :: nopen, nup
        integer, intent(inout) :: dorder(nopen)
        integer :: comb(nup)
        integer :: i, j
        logical :: bInc

        ! Initialise
        if (dorder(1) == -1) then
            dorder(1:nup) = 0
            dorder(nup + 1:nopen) = 1
        else
            ! Get the list of positions of the beta electrons
            j = 0
            do i = 1, nopen
                if (dorder(i) == 0) then
                    j = j + 1
                    comb(j) = i

                    ! Have we reached the last possibility?
                    if (j == 1 .and. i == nopen - nup + 1) then
                        dorder(1) = -1
                        return
                    end if

                    if (j == nup) exit
                end if
            end do

            do i = 1, nup
                bInc = .false.
                if (i == nup) then
                    bInc = .true.
                else if (i < nup) then
                    if (comb(i + 1) /= comb(i) + 1) bInc = .true.
                end if

                if (bInc) then
                    comb(i) = comb(i) + 1
                    exit
                else
                    comb(i) = i
                end if
            end do

            dorder = 1
            dorder(comb) = 0
        end if
    end subroutine

end module
