module get_excit

    use constants
    use SystemData, only: nel
    use bit_rep_data, only: NIfTot
    use DeterminantData, only: write_det
    implicit none


contains

#ifndef __DEBUG
    pure &
#endif
    subroutine make_single (nI, nJ, elec, tgt, ex, tParity)

        integer, intent(in) :: nI(nel), elec, tgt
        integer, intent(out) :: ex(2,2), nJ(nel)
        logical, intent(out) :: tParity
        character(*), parameter :: this_routine = 'make_excit'
        integer :: i, src

        ! Initialise return values, Returned excitation matrix includes the
        ! orbitals of the source, rather than their position.
        nJ = nI
        ex(1,1) = nI(elec)
        ex(2,1) = tgt

        ! Get the parity
        src = nI(elec)
        if (src < tgt) then

            ! How far do we have to move the unaffected orbitals?
            do i = elec+1, nel
                if (tgt < nJ(i)) then
                    nJ(i-1) = tgt
                    exit
                else
                    nJ(i-1) = nJ(i)
                end if
            end do
            if (i == nel+1) nJ(nel) = tgt

        else

            ! How far do we have to move the unaffected orbitals?
            do i = elec-1, 1, -1
                if (tgt > nJ(i)) then
                    nJ(i+1) = tgt
                    exit
                else
                    nJ(i+1) = nJ(i)
                end if
            end do
            if (i == 0) nJ(1) = tgt

        end if

        ! Magic! Avoids conditional tests.
        tParity = .not. btest(elec - i, 0)

#ifdef __DEBUG
        ! This is a useful (but O[N]) check to test the generated determinant.
        call IsSymAllowedExcit(nI, nJ, 2, ExcitMat)
#endif

    end subroutine


#ifndef __DEBUG
    pure &
#endif
    subroutine make_double (nI, nJ, elec1, elec2, tgt1, tgt2, ex, tParity)

        integer, intent(in) :: nI(nel), elec1, elec2, tgt1, tgt2
        integer, intent(out) :: ex(2,2), nJ(nel)
        logical, intent(out) :: tParity
        character(*), parameter :: this_routine = 'make_excit'
        integer :: i, k, elecs(2), srcs(2), tgts(2), pos_moved

        ! Get the source/target terms in order!
        elecs = (/min(elec1, elec2), max(elec2, elec2)/)
        tgts = (/min(tgt1, tgt2), max(tgt1, tgt2)/)

        ! Fill in the excitation matrix
        srcs = nI(elecs)
        ex(1,1:2) = srcs
        ex(2,1:2) = tgts

        ! Initialise return value
        nJ = nI

        ! As we move these around we need to do some playing!
        if (srcs(1) < tgts(1) .and. srcs(2) < tgts(1)) then
            elecs(2) = elecs(2) - 1
        end if

        ! Count how far we have moved normal orbitals
        pos_moved = 0
        do k = 1, 2

            ! If we need to search up or down depends on the relative sizes
            ! of src/tgt
            if (srcs(k) < tgts(k)) then

                ! How far do we have to move the unaffected orbitals?
                if (elecs(k) == nel) then
                    i = nel+1
                    nJ(nel) = tgts(k)
                else
                    do i = elecs(k)+1, nel
                        if (tgts(k) < nJ(i)) then
                            nJ(i-1) = tgts(k)
                            exit
                        else
                            nJ(i-1) = nJ(i)
                        end if
                    end do
                    if (i == nel+1) nJ(nel) = tgts(k)

                end if

            else

                ! How far do we have to move the unaffected orbitals?
                if (elecs(k) == 1) then
                    i = 0
                    nJ(1) = tgts(k)
                else
                    do i = elecs(k)-1, 1, -1
                        if (tgts(k) > nJ(i)) then
                            nJ(i+1) = tgts(k)
                            exit
                        else
                            nJ(i+1) = nJ(i)
                        end if
                    end do
                    if (i == 0) nJ(1) = tgts(k)
                end if

            end if

            pos_moved = pos_moved + elecs(k) - i + 1

        end do

        tParity = btest(pos_moved, 0)
!        parity = 1 - 2 * modulo(pos_moved, 2)

#ifdef __DEBUG
        ! This is a useful (but O[N]) check to test the generated determinant.
        call IsSymAllowedExcit(nI, nJ, 2, ExcitMat)
#endif

    end subroutine


end module
