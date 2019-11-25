#include "macros.h"

module cepa_shifts 

    use constants, only: dp, inum_runs
    use replica_data, only: diagsft 
    use SystemData, only: nel, nOccAlpha, nBasis
    use cc_amplitudes, only: t_cc_amplitudes, cc_singles_factor, & 
                        cc_doubles_factor, cc_triples_factor, cc_quads_factor

    implicit none 

    ! maybe i should start and store all the necessary flags and variables 
    ! in the files the functionality belongs to.. is this less messy? 

    logical :: t_cepa_shift = .false. 
    character(4) :: cepa_method
    real(dp) :: aqcc_factor = 0.0_dp
    logical :: t_apply_full_cepa = .false.

    ! there are different sorts of shift i can apply.. which and also there 
    ! is a choice on which excitation level i want to apply those 
    ! in  0 = <ijab|H - E_0 - D_ij|Psi>
    ! and 0 =   <ia|H - E_0 - D_i |Psi>
    ! method        D_ij                            D_i
    ! -------------------------------------------------
    ! cepa(0)       0                               0
    ! cepa(1)       1/2\sum_k(e_ik + e_jk)          \sum_k e_ik
    ! cepa(3)   -e_ij + \sum_k(e_ik + e_jk)     -e_ii + 2\sum_k e_ik
    ! acpf          2 E_c/N                          2 E_c / N
    ! aqcc      [1 - (N-3)(N-2)/(N(N-1))]E_c    [1 - (N-3)(N-2)/(N(N-1))]E_c

    ! so the best thing to do would be to make the input in the Calc Block: 
    ! CEPA 'method' 

    ! and depending on the method string determine the shifts D_ij and D_i
    abstract interface 
        function cepa_shift_t(run)
            use constants, only: dp 
            integer, intent(in) :: run
            real(dp) :: cepa_shift_t
        end function cepa_shift_t

        function cepa_shift_ex_level_t(run, ex_level)
            use constants, only: dp 
            integer, intent(in) :: run, ex_level 
            real(dp) :: cepa_shift_ex_level_t
        end function cepa_shift_ex_level_t

    end interface

    procedure(cepa_shift_t), pointer :: cepa_shift_single
    procedure(cepa_shift_t), pointer :: cepa_shift_double

    procedure(cepa_shift_ex_level_t), pointer :: cepa_shift
    
contains 

    subroutine init_cepa_shifts() 
        use DetBitOps, only: TestClosedShellDet
        use FciMCData, only: ilutref

        character(*), parameter :: this_routine = "init_cepa_shifts"

        integer :: i
        print *, "init cepa shifts "

        if (.not. allocated(ilutref)) then 
            call stop_all(this_routine, "reference det not yet set!")
        end if

        do i = 1, inum_runs
            if (.not. TestClosedShellDet(ilutref(:,i))) then 
                call stop_all(this_routine, "Cepa shifts only for closed shell reference!")
            end if
        end do

        select case(trim(adjustl(cepa_method))) 

        case ('0') 
            ! here the shift has to cancel the correlation energy, but i can't 
            ! point to the shift.. hm.. i guess i can't do that so nicely.. 
            cepa_shift_single => cepa_0
            cepa_shift_double => cepa_0

        case ('acpf') 

            ! here it gets tricky.. it would be nice if we have the shift here.
            ! i actually should us procedure pointers i guess.. 
            cepa_shift_single => cepa_acpf
            cepa_shift_double => cepa_acpf
        
        case ('aqcc') 

            ! is it orbital or electrons here? 
            if (nel <= 3) then 
                call stop_all(this_routine, "not enough electrons for aqcc shift!")
            end if

            aqcc_factor = (1.0_dp - real((nel - 3)*(nel - 2),dp)/real(nel*(nel - 1), dp))

            cepa_shift_single => cepa_aqcc
            cepa_shift_double => cepa_aqcc
!             
        case default 

            call stop_all(this_routine, "not recognised cepa shift!")

        end select 

        ! i have to point to the cc-version or to the CISD version too.. 
        if (t_cc_amplitudes) then 
            cepa_shift => cepa_shift_cc
        else
            cepa_shift => cepa_shift_cisd
        end if

    end subroutine init_cepa_shifts

    real(dp) function cepa_shift_cc(run, ex_level) 
        ! this routine should get used when we want to adapt the cepa-shift
        ! with biasing through the occupation of the higher order 
        ! excitations 
        integer, intent(in) :: run, ex_level
        if (ex_level == 1) then
            cepa_shift_cc = cepa_shift_single(run) * cc_singles_factor()
        else if (ex_level == 2) then 
            cepa_shift_cc = cepa_shift_double(run) * cc_doubles_factor()
        else 
            if (t_apply_full_cepa) then 
                if (ex_level == 3) then
                    ! i guess the shift could be the same for triples, just the 
                    ! factor changes
                    cepa_shift_cc = cepa_shift_double(run) * cc_triples_factor()
                else if (ex_level == 4) then 
                    cepa_shift_cc = cepa_shift_double(run) * cc_quads_factor()
                else 
                    ! the farther one goes out the less occupied the excited space 
                    ! becomes.. so maybe it is fair to assume to apply the 
                    ! full shift.. 
                    cepa_shift_cc = cepa_shift_double(run) 
                end if
            else 
                cepa_shift_cc = 0.0_dp
            end if
        end if

    end function cepa_shift_cc

    real(dp) function cepa_shift_cisd(run, ex_level) 
        ! this is the bare cepa shift, which would be used in a CISD calculation
        ! with no higher excitations then doubles
        integer, intent(in) :: run, ex_level
        if (ex_level == 1) then 
            cepa_shift_cisd = cepa_shift_single(run) 
        else if (ex_level == 2) then 
            cepa_shift_cisd = cepa_shift_double(run) 
        else
            if (t_apply_full_cepa) then 
                ! in this approach we apply the same shift, especially to 
                ! the higher order of excitations.. 
                ! but maybe also with an correction.. tbd
                cepa_shift_cisd = cepa_shift_double(run) 
            else
                ! with the change in the death-step to - (S - D) i then have to 
                ! set this to 0 here, if i want to apply the full shift
!                 cepa_shift = diagsft(run)
                cepa_shift_cisd = 0.0_dp
            end if
        end if
    end function cepa_shift_cisd

    real(dp) function cepa_0(run)
        integer, intent(in) :: run
        ! change this back so it cancels the actual shift 
        ! so we can adjust this outside with the cc-amplitudes estimate
!         cepa_0 = 0.0_dp
        cepa_0 = diagsft(run)
    end function cepa_0

    real(dp) function cepa_acpf(run) 
        integer, intent(in) :: run 
        ! do i use the shift or the projected energy here?? tbd
!         cepa_acpf = 2.0_dp * diagsft(run) / real(nel, dp)
        ! change the implementation, so that it actually gives 
        ! S - D = 2S/N -> D = S(1 - 2/N)
        cepa_acpf = diagsft(run)*(1.0_dp - 2.0_dp/real(nel,dp))
    
    end function cepa_acpf

    real(dp) function cepa_aqcc(run) 
        integer, intent(in) :: run 
!         cepa_aqcc = aqcc_factor * diagsft(run)
        cepa_aqcc = diagsft(run)*(1.0_dp - aqcc_factor)
    end function cepa_aqcc

end module cepa_shifts
