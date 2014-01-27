#include "macros.h"

module DeterminantData
    use SystemData, only: nel, tCSF
    use csf_data, only: iscsf, csf_orbital_mask, csf_yama_bit
    use MemoryManager, only: TagIntType
    implicit none

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

    subroutine write_det (nunit, nI, lTerm)

        ! Write the specified determinant (or CSF) to the output unit nunit.
        ! Terminate the line (with '\n') if lTerm = .true.
        !
        ! In: nunit - Output (file) unit
        !     nI    - Determinant to output
        !     lTerm - Terminate line with newline character?

        integer, intent(in) :: nunit, nI(nel)
        logical, intent(in) :: lTerm

        call write_det_len (nunit, nI, nel, lterm)
    end subroutine write_det

    subroutine write_det_len (nunit, nI, nlen, lterm)

        ! Worker function for the above. Can be accessed to print an unusual
        ! lengthed determinant.

        integer, intent(in) :: nunit, nlen, nI(nlen)
        logical, intent(in) :: lTerm
        integer :: i, elec
        logical open_shell, bCSF

        if(nlen>0) then
           ! Is this a csf?
           bCSF = tCSF .and. iscsf(nI)
           open_shell = .false.
        endif

        ! Start with a bracket, and loop over all the electrons
        write(nunit,'("(")',advance='no')
        do i=1,nlen
            ! If this is a csf, extract the orbital number, and test
            ! if we have passed all the closed shell electrons
            if (bCSF) then
                if ((.not.open_shell) .and. &
                    btest(nI(i), csf_yama_bit)) open_shell = .true.

                elec = iand(nI(i), csf_orbital_mask)
            else
                elec = nI(i)
            endif

            ! Write out the orbital number, and +/- for open shell csf e-
            if (open_shell) then
                write(nunit,'(i4)',advance='no') elec
                if (btest(nI(i), csf_yama_bit)) then
                    write(nunit,'("+")',advance='no')
                else
                    write(nunit,'("-")',advance='no')
                endif
            else
                write(nunit,'(i5)',advance='no') elec
            endif
            if (i /= nlen) write(nunit,'(",")',advance='no')
        enddo

        ! Close the written determinant off
        write(nunit,'(")")',advance='no')
        if (lTerm) write(nunit,*)
    end subroutine write_det_len

    subroutine get_lexicographic (dorder, nopen, nup)

        ! Unlike the csf version, this uses 1 == alpha, 0 = beta.

        integer, intent(in) :: nopen, nup
        integer, intent(inout) :: dorder(nopen)
        integer :: comb(nup)
        integer :: i, j
        logical :: bInc

        ! Initialise
        if (dorder(1) == -1) then
            dorder(1:nup) = 0
            dorder(nup+1:nopen) = 1
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
                    endif

                    if (j == nup) exit
                endif
            enddo

            do i = 1, nup
                bInc = .false.
                if (i == nup) then
                    bInc = .true.
                else if (i < nup) then
                    if (comb(i+1) /= comb(i) + 1) bInc = .true.
                endif

                if (bInc) then
                    comb(i) = comb(i) + 1
                    exit
                else
                    comb(i) = i
                endif
            enddo

            dorder = 1
            dorder(comb) = 0
        endif
    end subroutine



end module
