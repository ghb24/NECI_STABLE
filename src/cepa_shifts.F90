#include "macros.h"

module cepa_shifts 

    use constants, only: dp, inum_runs
    use replica_data, only: diagsft 
    use SystemData, only: nel

    implicit none 

    ! maybe i should start and store all the necessary flags and variables 
    ! in the files the functionality belongs to.. is this less messy? 

    logical :: t_cepa_shift = .false. 
    character(4) :: cepa_method
    real(dp) :: aqcc_factor = 0.0_dp

!     real(dp), allocatable :: cepa_shift_single(:), cepa_shift_double(:)

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

!         function cepa_shift_ex_level_t(run, ex_level)
!             use constants, only: dp 
!             integer, intent(in) :: run, ex_level 
!             real(dp) :: cepa_shift_ex_level_t
!         end function cepa_shift_ex_level_t

    end interface

    procedure(cepa_shift_t), pointer :: cepa_shift_single
    procedure(cepa_shift_t), pointer :: cepa_shift_double

!     procedure(cepa_shift_ex_level_t), pointer :: cepa_shift
    
contains 

    subroutine init_cepa_shifts() 

        character(*), parameter :: this_routine = "init_cepa_shifts"

        ! i have to allocate it for each replica
!         allocate(cepa_shift_single(inum_runs))
!         allocate(cepa_shift_double(inum_runs))

        select case(trim(adjustl(cepa_method))) 

        case ('0') 
            ! here the shift has to cancel the correlation energy, but i can't 
            ! point to the shift.. hm.. i guess i can't do that so nicely.. 
            cepa_shift_single => cepa_0
            cepa_shift_double => cepa_0

        case ('1') 
            
            cepa_shift_single => cepa_1_single
            cepa_shift_double => cepa_1_double
            call stop_all(this_routine, "cepa(1) not yet implemented!")

        case ('3') 

            cepa_shift_single => cepa_3_single
            cepa_shift_double => cepa_3_double

            call stop_all(this_routine, "cepa(3) not yet implemented!")

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

    end subroutine init_cepa_shifts

    real(dp) function cepa_shift(run, ex_level) 
        integer, intent(in) :: run, ex_level
        if (ex_level == 1) then 
            cepa_shift = cepa_shift_single(run) 
        else if (ex_level == 2) then 
            cepa_shift = cepa_shift_double(run) 
        else
            cepa_shift = diagsft(run)
        end if
    end function cepa_shift

    real(dp) function cepa_0(run)
        integer, intent(in) :: run
        cepa_0 = 0.0_dp
    end function cepa_0

    real(dp) function cepa_1_single(run)
        integer, intent(in) :: run
        ! todo
    end function cepa_1_single

    real(dp) function cepa_1_double(run)
        integer, intent(in) :: run 
        ! todo
    end function cepa_1_double

    real(dp) function cepa_3_single(run) 
        integer, intent(in) :: run 
        ! todo 
    end function cepa_3_single

    real(dp) function cepa_3_double(run) 
        integer, intent(in) :: run 
        ! todo 
    end function cepa_3_double

    real(dp) function cepa_acpf(run) 
        integer, intent(in) :: run 
        ! do i use the shift or the projected energy here?? tbd
        cepa_acpf = 2.0_dp * diagsft(run) / real(nel, dp)
    end function cepa_acpf

    real(dp) function cepa_aqcc(run) 
        integer, intent(in) :: run 
        cepa_aqcc = aqcc_factor * diagsft(run)
    end function cepa_aqcc

end module cepa_shifts
