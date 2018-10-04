module hdiag_from_excit

    use bit_reps, only: NIfTot
    use constants, only: dp, n_int
    use DetBitOps, only: CalcOpenOrbs, TestClosedShellDet, FindBitExcitLevel
    use Integrals_neci, only: get_umat_el
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use OneEInts, only: GetTMatEl
    use sltcnd_mod, only: sltcnd
    use SystemData, only: nel, tExch, G1, tReltvy, tHPHF, tOddS_HPHF
    use UMatCache, only: gtID

    implicit none

    contains

    function get_hdiag_from_excit(nI, nJ, iLutnJ, IC, ex, hel_old) result(hel)

        integer, intent(in) :: nI(nel), nJ(nel), IC, ex(2,2)
        integer(n_int), intent(in) :: iLutnJ(0:NIfTot)
        HElement_t(dp), intent(in) :: hel_old

        HElement_t(dp) :: hel
        
        if (IC == 1) then
            hel = get_hdiag_from_sing_excit(nI, ex(:,1), hel_old)
        else if (IC == 2) then
            hel = get_hdiag_from_doub_excit(nI, ex, hel_old)
        else
            hel = 0.0_dp
        end if

        if (tHPHF) call correct_hdiag_hphf(nJ, iLutnJ, hel)

    end function get_hdiag_from_excit

    function get_hdiag_from_sing_excit(nI, ex, hel_old) result(hel)

        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).

        integer, intent(in) :: nI(nel), ex(2)
        HElement_t(dp), intent(in) :: hel_old

        HElement_t(dp) :: hel
        integer :: id(nel), id_ex(2), i

        ! Correct one electron integral contribution
        hel = hel_old - GetTMATEl(ex(1), ex(1)) + GetTMATEl(ex(2), ex(2))

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)
        id_ex = gtID(ex)

        do i = 1, nel
            if (nI(i) /= ex(1)) then
                hel = hel - get_umat_el (id_ex(1), id(i), id_ex(1), id(i)) &
                          + get_umat_el (id_ex(2), id(i), id_ex(2), id(i))

                if ( tReltvy .or. (G1(ex(1))%Ms == G1(nI(i))%Ms) ) then
                    hel = hel + get_umat_el (id_ex(1), id(i), id(i), id_ex(1)) &
                              - get_umat_el (id_ex(2), id(i), id(i), id_ex(2))
                end if
            end if
        end do

    end function get_hdiag_from_sing_excit

    function get_hdiag_from_doub_excit(nI, ex, hel_old) result(hel)

        ! Calculate the  by the SlaterCondon Rules when the two
        ! determinants are the same (so we only need to specify one).

        integer, intent(in) :: nI(nel), ex(2,2)
        HElement_t(dp), intent(in) :: hel_old

        HElement_t(dp) :: hel
        integer :: id(nel), ex_ordered(2,2), id_ex(2,2), i

        ! Correct one electron integral contribution
        hel = hel_old - GetTMATEl(ex(1,1), ex(1,1)) + GetTMATEl(ex(2,1), ex(2,1)) &
                      - GetTMATEl(ex(1,2), ex(1,2)) + GetTMATEl(ex(2,2), ex(2,2))

        ex_ordered = ex
        if (G1(ex(1,1))%Ms /= G1(ex(2,1))%Ms) then
            ex_ordered(2,1) = ex(2,2)
            ex_ordered(2,2) = ex(2,1)
        end if

        ! Obtain the spatial rather than spin indices if required
        id = gtID(nI)
        id_ex = gtID(ex_ordered)

        do i = 1, nel
            if (nI(i) /= ex_ordered(1,1)) then
                hel = hel - get_umat_el (id_ex(1,1), id(i), id_ex(1,1), id(i)) &
                          + get_umat_el (id_ex(2,1), id(i), id_ex(2,1), id(i))

                if ( tReltvy .or. (G1(ex_ordered(1,1))%Ms == G1(nI(i))%Ms) ) then
                    hel = hel + get_umat_el (id_ex(1,1), id(i), id(i), id_ex(1,1)) &
                              - get_umat_el (id_ex(2,1), id(i), id(i), id_ex(2,1))
                end if
            end if
        end do

        do i = 1, nel
            if (nI(i) /= ex_ordered(1,1) .and. nI(i) /= ex_ordered(1,2)) then
                hel = hel - get_umat_el (id_ex(1,2), id(i), id_ex(1,2), id(i)) &
                          + get_umat_el (id_ex(2,2), id(i), id_ex(2,2), id(i))

                if ( tReltvy .or. (G1(ex_ordered(1,2))%Ms == G1(nI(i))%Ms) ) then
                    hel = hel + get_umat_el (id_ex(1,2), id(i), id(i), id_ex(1,2)) &
                              - get_umat_el (id_ex(2,2), id(i), id(i), id_ex(2,2))
                end if
            end if
        end do

        hel = hel - get_umat_el (id_ex(1,2), id_ex(2,1), id_ex(1,2), id_ex(2,1)) &
                  + get_umat_el (id_ex(2,1), id_ex(2,2), id_ex(2,1), id_ex(2,2))

        if (G1(ex_ordered(1,1))%Ms == G1(ex_ordered(1,2))%Ms) then
            hel = hel + get_umat_el (id_ex(1,2), id_ex(2,1), id_ex(2,1), id_ex(1,2)) &
                      - get_umat_el (id_ex(2,1), id_ex(2,2), id_ex(2,2), id_ex(2,1))
        end if

    end function get_hdiag_from_doub_excit

    function get_hdiag_bare_hphf(nI, iLutnI, hel_old) result(hel)

        ! <X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are
        ! the two spin-coupled dets which make up X.

        integer, intent(in) :: nI(nel) 
        integer(n_int), intent(in) :: iLutnI(0:NIfTot)
        HElement_t(dp), intent(in) :: hel_old

        HElement_t(dp) :: hel

        integer(n_int) :: iLutnI2(0:NIfTot)
        integer :: ExcitLevel, OpenOrbs
        HElement_t(dp) :: MatEl2

        hel = hel_old

        if (.not. TestClosedShellDet(iLutnI)) then
            ! See if there is a cross-term. If there is, then remove this
            ! from hel_old to get the desired value
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            ExcitLevel = FindBitExcitLevel(iLutnI, iLutnI2, 2)
            if (ExcitLevel <= 2) then
                call CalcOpenOrbs (iLutnI, OpenOrbs)
                MatEl2 = sltcnd (nI, iLutnI, iLutnI2)

                if (tOddS_HPHF) then
                    if (mod(OpenOrbs,2) == 1) then
                        ! Subtract cross terms if determinant is antisymmetric.
                        hel = hel - MatEl2
                    else
                        hel = hel + MatEl2
                    end if
                else
                    if (mod(OpenOrbs,2) == 1) then
                        ! Subtract cross terms if determinant is antisymmetric.
                        hel = hel + MatEl2
                    else
                        hel = hel - MatEl2
                    end if
                end if
            end if
        end if

    end function get_hdiag_bare_hphf

    subroutine correct_hdiag_hphf(nJ, iLutnJ, hel)

        ! <X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are
        ! the two spin-coupled dets which make up X.

        integer, intent(in) :: nJ(nel)
        integer(n_int), intent(in) :: iLutnJ(0:NIfTot)
        HElement_t(dp), intent(inout) :: hel

        integer(n_int) :: iLutnJ2(0:NIfTot)
        integer :: ExcitLevel, OpenOrbs
        HElement_t(dp) :: MatEl2

        if (.not. TestClosedShellDet(iLutnJ)) then
            ! See if there is a cross-term. If there is, then remove this
            ! from hel_old to get the desired value
            call FindExcitBitDetSym(iLutnJ, iLutnJ2)
            ExcitLevel = FindBitExcitLevel(iLutnJ, iLutnJ2, 2)
            if (ExcitLevel <= 2) then
                call CalcOpenOrbs (iLutnJ, OpenOrbs)
                MatEl2 = sltcnd (nJ, iLutnJ, iLutnJ2)

                if (tOddS_HPHF) then
                    if (mod(OpenOrbs,2) == 1) then
                        ! Subtract cross terms if determinant is antisymmetric.
                        hel = hel + MatEl2
                    else
                        hel = hel - MatEl2
                    end if
                else
                    if (mod(OpenOrbs,2) == 1) then
                        ! Subtract cross terms if determinant is antisymmetric.
                        hel = hel - MatEl2
                    else
                        hel = hel + MatEl2
                    end if
                end if
            end if
        end if

    end subroutine correct_hdiag_hphf

end module hdiag_from_excit
