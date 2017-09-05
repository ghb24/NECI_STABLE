#include "macros.h"

module cepa_shifts 

    use constants, only: dp, inum_runs
    use replica_data, only: diagsft 
    use SystemData, only: nel, nOccAlpha, nBasis

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
        use DetBitOps, only: TestClosedShellDet
        use FciMCData, only: projedet

        character(*), parameter :: this_routine = "init_cepa_shifts"

        integer :: i
        ! i have to allocate it for each replica
!         allocate(cepa_shift_single(inum_runs))
!         allocate(cepa_shift_double(inum_runs))

        print *, "init cepa shifts "
        if (.not. allocated(projedet)) then 
            call stop_all(this_routine, "reference det not yet set!")
        end if

        do i = 1, inum_runs
            if (.not. TestClosedShellDet(projedet(:,i))) then 
                call stop_all(this_routine, "Cepa shifts only for closed shell reference!")
            end if
        end do

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

    function calc_number_of_excitations(n_alpha, n_beta, max_excit, n_orbs) &
            result(n_excits)
        ! i might want a routine which calculates the number of all possible 
        ! excitations for each excitation level 
        integer, intent(in) :: n_alpha, n_beta, max_excit, n_orbs
        integer :: n_excits(max_excit)

        integer :: n_parallel(max_excit,2), i, j, k

        ! the number of all possible single excitations, per spin is: 
        
        n_parallel = 0
        n_excits = 0

        do i = 1, max_excit
            n_parallel(i,1) = calc_n_parallel_excitations(n_alpha, n_orbs, i)
            n_parallel(i,2) = calc_n_parallel_excitations(n_beta, n_orbs, i)
        end do

        ! i can set this up in a general way.. 

        do i = 1, max_excit

            ! the 2 parallel spin species always are added 
            n_excits(i) = sum(n_parallel(i,:)) 

            ! and i have to zig-zag loop over the other possible combinations 
            ! to lead to this level of excitation.. 
            j = i - 1 
            k = 1 

            do while (j >= k) 

                n_excits(i) = n_excits(i) + n_parallel(j,1) * n_parallel(k,2) 

                ! and if j /= k do it the other way around too 
                if (j /= k) then 
                    n_excits(i) = n_excits(i) + n_parallel(j,2) * n_parallel(k,1)

                end if
                ! and in/decrease the counters 
                j = j - 1 
                k = k + 1

                ! in this way i calculate eg.: 
                ! n_ij^ab = n_ij_ab(u + d) + n_i^a(u) * n_i^a(d) 
                ! and 
                ! n_ijk^abc = n_ijk^abc(u + d) + n_ij^ab(u) * n_i^a(d) and spin-flipped..

            end do
        end do

    end function calc_number_of_excitations

    function calc_n_single_excits(n_elecs, n_orbs) result(n_single_excits)
        ! this calculates the number of possible single excitations for a 
        ! given spin species of electrons and the number of spatial orbitals!
        integer, intent(in) :: n_elecs, n_orbs
        integer :: n_single_excits

        ! it is the number of electrons which can be excited and the number 
        ! of available orbitals for this spin!
!         n_single_excits = n_elecs * (n_orbs - n_elecs)

        ! with binomial it is just: 
        n_single_excits = binomial(n_elecs, 1) * binomial(n_orbs - n_elecs, 1)

    end function calc_n_single_excits

    function calc_n_parallel_excitations(n_elecs, n_orbs, ic) result(n_parallel)
        ! this function determines the number of parallel spin excitaiton 
        ! with each electrons having the same spin for a given excitation 
        ! level and number of electrons and available orbitals for this spin 
        integer, intent(in) :: n_elecs, n_orbs, ic 
        integer :: n_parallel

        n_parallel = binomial(n_elecs, ic) * binomial(n_orbs - n_elecs, ic) 

    end function calc_n_parallel_excitations

    integer function binomial(n, k)
        ! write a new binomial function using the fortran2008 standard 
        ! gamma function 
        integer, intent(in) :: n, k 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "binomial"
#endif

        ASSERT(n >= 0)
        ASSERT(k >= 0)

        if (k > n) then 
            binomial = 0
        else
            binomial = nint(gamma(real(n) + 1) / (gamma(real(n - k) + 1) * gamma(real(k) + 1)))
        end if

    end function binomial 

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
