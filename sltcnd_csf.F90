module sltcnd_csf_mod
    use SystemData, only: nel
    use HElem
    use timing

contains
    type(HElement) function sltcnd_csf (nI, nJ, IC)
        
        ! Use the Slater-Condon Rules to evaluate the H matrix element between
        ! two determinants. Assume CSF ordering of orbitals (closed pairs
        ! followed by open shell electrons). However, this is NOT to be passed
        ! CSFS - it is to evaluate the component determinants.
        !
        ! In:  nI(nel) nJ(nel) - The determinants to evaluate
        !      IC              - The number of orbitals nI, nJ differ by
        ! Ret: sltcnd_csf      - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: IC
        
        ! For debugging purposes, measure the time contribution.
        type(timer), save :: slt_time
        slt_time%timer_name = 'sltcnd_csf'
        call set_timer(slt_time)
        
        select case (IC)
        case (0)
            ! The determinants are exactly the same
            sltcnd_csf = sltcnd_csf_0 (nI, nJ)

        case (1)
            ! The determinants differ by only one orbital
            sltcnd_csf = sltcnd_csf_1 (nI, nJ)

        case (2)
            ! The determinants differ by two orbitals
            sltcnd_csf = sltcnd_csf_2 (nI, nJ)

        case default
            ! The determinants differ by more than two orbitals
            sltcnd_csf%v = 0
        endselect

        call halt_timer(slt_time)
    end function


    type (HElement) function sltcnd_csf_0 (nI, nJ)
        integer, intent(in) :: nI(nel), nJ(nel)
    end function sltcnd_csf_0

    type (HElement) function sltcnd_csf_1 (nI, nJ)
        integer, intent(in) :: nI(nel), nJ(nel)
    end function sltcnd_csf_1
    
    type (HElement) function sltcnd_csf_2 (nI, nJ)
        integer, intent(in) :: nI(nel), nJ(nel)
    end function sltcnd_csf_2

        


end module
