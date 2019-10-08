#include "macros.h"
module get_excit

    use constants
    use SystemData, only: nel
    use bit_rep_data, only: NIfTot
    use DeterminantData, only: write_det
    use sym_general_mod, only: SymAllowedExcit
    implicit none


contains

    subroutine make_single (nI, nJ, elec, tgt, ex, tParity)

        integer, intent(in) :: nI(nel), elec, tgt
        integer, intent(out) :: ex(2,2), nJ(nel)
        logical, intent(out) :: tParity
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'make_single'
#endif
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
        tParity = mod(elec - i, 2) == 0

#ifdef __DEBUG
        ! This is a useful (but O[N]) check to test the generated determinant.
        if (.not. SymAllowedExcit(nI, nJ, 1, ex)) &
            call stop_all(this_routine, 'Generating invalid excitation')
#endif

    end subroutine


    subroutine make_double (nI, nJ, elec1, elec2, tgt1, tgt2, ex, tParity)

        integer, intent(in) :: nI(:), elec1, elec2, tgt1, tgt2
        integer, intent(inout) :: ex(2,2), nJ(:)
        logical, intent(out) :: tParity
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'make_double'
#endif

        integer :: i, k, elecs(2), srcs(2), tgts(2), pos_moved

        ! Get the source/target terms in order!
        elecs = (/min(elec1, elec2), max(elec1, elec2)/)
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
        if (.not. SymAllowedExcit(nI, nJ, 2, ex)) then
            print *, "nI: ", nI
            print *, "nJ: ", nJ
            print *, "elecs: ", ex(1,:)
            print *, "orbs: ", ex(2,:)
            call stop_all(this_routine, 'Generated invalid excitation')
        end if
#endif

    end subroutine

    function exciteIlut(ilut,src,orbs) result(ilutJ)
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(in) :: src(2), orbs(2)
      integer(n_int) :: ilutJ(0:NIfTot)

      ilutJ = ilut
      clr_orb (ilutJ, src(1))
      clr_orb (ilutJ, src(2))
      set_orb (ilutJ, orbs(1))
      set_orb (ilutJ, orbs(2))

    end function exciteIlut



end module
