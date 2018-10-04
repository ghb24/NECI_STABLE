module hdiag_from_excit

    use constants, only: dp
    use Integrals_neci, only: get_umat_el
    use OneEInts, only: GetTMatEl
    use SystemData, only: nel, tExch, G1, tReltvy
    use UMatCache, only: gtID

    implicit none

    contains

    function get_hdiag_from_excit(nI, IC, ex, hel_old) result(hel)

        integer, intent(in) :: nI(nel), IC, ex(2,2)
        HElement_t(dp), intent(in) :: hel_old

        HElement_t(dp) :: hel
        
        if (IC == 1) then
            hel = get_hdiag_from_sing_excit(nI, ex(:,1), hel_old)
        else if (IC == 2) then
            hel = get_hdiag_from_doub_excit(nI, ex, hel_old)
        else
            hel = 0.0_dp
        end if

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

end module hdiag_from_excit
