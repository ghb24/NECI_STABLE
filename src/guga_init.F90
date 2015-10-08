#include "macros.h"
! guga module containing all necessary functionality needed to initialize 
! a guga simulation
module guga_init
    ! module use statements
    use SystemData, only: tCSF, tSPN, tHPHF, lNoSymmetry, STOT, nEl, nReps, &
        nBasis, tGUGA, tNoBrillouin, tExactSizeSpace, tUHF
    use hist_data, only: tHistSpawn
    use LoggingData, only: tCalcFCIMCPsi, tPrintOrbOcc
    use spin_project, only: tSpinProject
    use bit_rep_data, only: tUseFlags, nIfD
    use guga_data, only: init_guga
    ! variable declaration
    implicit none


contains

    subroutine checkInputGUGA()
        ! routine to check if all the input parameters given are consistent
        ! and otherwise stops the excecution
        ! is called inf checkinput() in file readinput.F90
        character(*), parameter :: routineName = 'checkInputGUGA'

        ! check in certain system options have conflicts:
        if (tCSF) then
            call stop_all(routineName, &
                "Cannot use two CSF implementations tGUGA and tCSF at the same time!")
        end if

        if (tSPN) then
            call stop_all(routineName, &
                "GUGA not yet implemented with spin restriction SPIN-RESTRICT!")
        end if

        if (tHPHF) then
            call stop_all(routineName, &
                "GUGA not compatible with HPHF option!")
        end if

        if (tSpinProject) then
            call stop_all(routineName, &
                "GUGA not compatible with tSpinProject!")
        end if

        if (.not. lNoSymmetry) then
            call stop_all(routineName, &
                "GUGA not yet implemented with k-point symmetry!")
        end if
        
        if (tExactSizeSpace) then
            call stop_all(routineName, &
                "calculation of exact Hilbert space size not yet implemented with GUGA!")
        end if

        ! probably also have to assert against all the hist and exact 
        ! calculation flags.. also rdms... and certain excitation generators
        if (tHistSpawn) then
            call stop_all(routineName, &
                "HISTSPAWN not yet compatible with GUGA!")
        end if

        if (tCalcFCIMCPsi) then
            call stop_all(routineName, &
                "PRINTFCIMCPSI not yet compatible with GUGA!")
        end if

        if (tPrintOrbOcc) then
            call stop_all(routineName, &
                "PRINTORBOCCS not yet implemented with GUGA!")
        end if


        if (.not. tNoBrillouin) then
            call stop_all(routineName, &
                "Brillouin theorem not valid for GUGA approach!(I think atleast...)")
        end if

        ! also check if provided input values match:
        ! CONVENTION: STOT in units of h/2!
        if (STOT .gt. nEl) then
            call stop_all(routineName, &
                "total spin S in units of h/2 cannot be higher than number of electrons!")
        end if

        if (mod(STOT,2) /= mod(nEl,2)) then
            call stop_all(routineName, &
                "number of electrons nEl and total spin S in units of h/2 must have same parity!")
        end if

        ! maybe more to come...
        ! UHF basis is also not compatible with guga? or not... or atleast 
        ! i am not yet implementing it in such a way so it can work
        if (tUHF) then
            call stop_all(routineName, &
                "GUGA approach and UHF basis not yet (or never?) compatible!")
        end if

        ! assert that tUseFlags is set, to be able to encode deltaB values
        ! in the ilut representation for excitation generation
        !if (.not.tUseFlags) then
        !    call stop_all(routineName, &
        !        "tUseFlags has to be .true. to encode deltaB values in ilut!")
        !end if
        ! cannot do assertion here, since flag is set afterwards, but changed 
        ! corresponding code, so flag is set.
        
    end subroutine checkInputGUGA

    subroutine SysInitGUGA()
        ! in general this just prints out information of the state of the 
        ! guga calculation and to avoid having repeating code elements, does 
        ! not check input options and values again
        ! this routine is called in SysInit() of System_neci.F90
        write(6,*) ' ************ Using the GUGA-CSF implementation **********'
        write(6,*) ' Restricting the total spin of the system, tGUGA : ', tGUGA
        write(6,'(A,I5)') '  Restricting total spin S in units of h/2 to ', STOT
        write(6,*) ' So eg. S = 1 corresponds to one unpaired electron '
        write(6,*) ' not quite sure yet how to deal with extensively used m_s quantum number..'
        write(6,*) ' NOTE: for now, although SPIN-RESTRICT is set off, internally m_s(LMS) '
        write(6,*) ' is set to STOT, to make use of reference determinant creations already implemented'
        write(6,*) ' Since NECI always seems to take the beta orbitals first for open shell or '
        write(6,*) ' spin restricted systems, associate those to positively coupled +h/2 orbitals '
        write(6,*) ' to always ensure a S >= 0 value!' 
        write(6,*) ' *********************************************************'
        ! maybe have to do more variable setting here in the future...
        ! probably should set the nReps quantity here
        ! maybe after all I do not need nReps quantity.. 
!         nReps = nBasis / 2


        ! have to ensure to use flags to be able to encode deltaB values
        ! in the ilut bit representation for the excitation generation
        tUseFlags = .true.


        ! initialize the pointers...
        ! todo: not sure if here..
        call init_guga

    end subroutine SysInitGUGA
! 
!     subroutine init_bit_rep_guga()
!         ! initialize special guga determinant structure for excitation 
!         ! generation
!         character(*), parameter :: this_routine = "init_bit_rep_guga"
! 
!         ! Structure of a bit representation:
! 
!         ! | 0-NIfD: Det | x0 | x1 | deltaB |
!         !
!         ! -------
!         ! (NIfD + 1) * 64-bits              Orbital rep.
!         !  1         * 64-bits              x0 matrix element
!         !  1         * 64-bits              x1 matrix element
!         !  1         * 32-bits              deltaB value
! 
!         ! number of integers needed for orbital structure
!         ! essentially could take this from the already initialized "normal"
!         ! determinant
! !         nIfD = int(basis / bit_n_int)
! 
!         ! integers like nIfY are zero anyway in the guga run
! 
!         ! we are storing real coefficients for the x0 and x1 matrix elements
!         ! so 64 bit integers have to be used and we cant store flags and 
!         ! matrix elements in the same integer
! 
!         ! according to the other initiatior function flags MUST be last
!         ! or else the sort function dont work anymore!
!         ! how does this affect my implementation? since i want to use 
!         ! add_ilut_lists which call ilut_lt and ilut_gt? 
!         ! probably have to use flags like normal and encode deltaB into those
! 
!         ! the only thing is that i need the "flags" to be at the last entry
!         ! but doesnt matter what it contains so i need nIfD + 3 entries
! 
!         ! so essentially only call this function after the init_bit_rep()
!         ! call because it just needs the initialised nIfD value
! 
!         ! and i have to rewrite set- and getDelta B and write new custom
!         ! extract and encode part sign functions to correctly adress the 
!         ! new ilut arrays
!         nIfGUGA = nIfD + 3
! 
!     end subroutine init_bit_rep_guga
! 
! 
end module guga_init





         
