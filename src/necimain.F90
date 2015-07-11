! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module EN2Molcas
   use constants
   implicit NONE
   save
   real(dp) :: NECI_E
end module EN2Molcas

#if defined(_MOLCAS_)

      subroutine NECImain(NECIen)
        USE EN2MOLCAS, only : NECI_E
        character(64) :: dummy1,dummy2
        real(8), intent (out) :: NECIen
        write(6,*) "STARTING NECI from Molcas"
        dummy1=' '
        dummy2=' '
       ! Indicate not called by CPMD, VASP, Molpro
        call NECICore(0,.False.,.False.,.false.,.true.,dummy1,dummy2)
       ! Once we got excited states energies we will add them here to ENER array.
!        write(6,*) 'NECI_E in necimain', NECI_E
        NECIen=NECI_E
      end subroutine NECImain

#elif defined(CBINDMPI)
    ! If calling from C, then we need to have an available fortran calling
    ! point available to the C-start point
    
    subroutine neci_main_c () bind(c)
        implicit none
        character(64) :: dummy1,dummy2

#ifdef __DEBUG
        write(6,*) 'STARTING NECI'
#endif
        dummy1=' '
        dummy2=' '
        ! Indicate not called by CPMD, VASP, Molpro
        call NECICore (0, .false., .false., .false., .false., dummy1, dummy2)

    end subroutine

#else

    ! necimain is the entry point for a standalone NECI.  It reads in an 
    ! input, and then runs the NECI Core
    program NECI
        implicit none
        character(64) :: dummy1,dummy2

#ifdef __DEBUG
        write(6,*) "STARTING NECI"
#endif

        dummy1=' '
        dummy2=' '
        ! Indicate not called by CPMD, VASP, Molpro
        call NECICore(0, .False., .False., .false., .false., dummy1, dummy2)

    end program NECI

#endif
