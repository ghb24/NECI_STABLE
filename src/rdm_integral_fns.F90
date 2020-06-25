module rdm_integral_fns

    use constants
    use excitation_types, only: DoubleExc_t

    implicit none

contains

    pure function one_elec_int(i, j) result(integral)

        use OneEInts, only: GetTMatEl

        integer, intent(in) :: i, j
        real(dp) :: integral

        integral = GetTMatEl(i, j)

    end function one_elec_int

    function two_elec_int(i, j, k, l) result(integral)

        use sltcnd_mod, only: sltcnd_2_kernel
        use UMatCache, only: gtID

        integer, intent(in) :: i, j, k, l
        real(dp) :: integral

        integral = real(sltcnd_2_kernel(DoubleExc_t(j, l, i, k)), dp)

    end function two_elec_int

    function GetPropInts(i, j, iprop) result(integral)

        use OneEInts, only: OneEPropInts
        integer, intent(in) :: i, j, iprop
        real(dp) :: integral

        integral = OneEPropInts(i, j, iprop)

    end function GetPropInts

end module rdm_integral_fns
