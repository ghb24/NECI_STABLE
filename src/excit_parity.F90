module get_excit

    use constants
    use SystemData, only: nel
    use bit_rep_data, only: NIfTot
    use DeterminantData, only: write_det
    implicit none


contains

    pure subroutine make_single (nI, ilutI, nJ, ilutJ, elec, tgt, ex, parity)

        integer, intent(in) :: nI(nel), elec, tgt
        integer, intent(out) :: ex(2,2), parity, nJ(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
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
        parity = 2 * modulo(elec - i, 2) - 1

    end subroutine


    subroutine make_double (nI, ilutI, nJ, ilutJ, elec1, elec2, tgt1, &
                                 tgt2, ex, parity)

        integer, intent(in) :: nI(nel), elec1, elec2, tgt1, tgt2
        integer, intent(out) :: ex(2,2), parity, nJ(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
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
                do i = elecs(k)+1, nel
                    if (tgts(k) < nJ(i)) then
                        nJ(i-1) = tgts(k)
                        exit
                    else
                        nJ(i-1) = nJ(i)
                    end if
                end do
                if (i == nel+1) nJ(nel) = tgts(k)

            else

                ! How far do we have to move the unaffected orbitals?
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

            pos_moved = pos_moved + elecs(k) - i + 1

        end do

        if (btest(pos_moved, 0)) then
            parity = -1
        else
            parity = 1
        end if
!        parity = 1 - 2 * modulo(pos_moved, 2)

    end subroutine


end module
