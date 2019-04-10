#if defined(CBINDMPI)
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
