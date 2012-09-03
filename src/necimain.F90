

#ifdef CBINDMPI

    ! If calling from C, then we need to have an available fortran calling
    ! point available to the C-start point
    
    subroutine neci_main_c () bind(c)

        write(6,*) 'STARTING NECI'
        ! Indicate not called by CPMD, VASP, Molpro
        call NECICore (0, .false., .false., .false.,' ',' ')

    end subroutine

#else

    ! necimain is the entry point for a standalone NECI.  It reads in an 
    ! input, and then runs the NECI Core
    program NECI

        implicit none
        write(6,*) "STARTING NECI"
        ! Indicate not called by CPMD, VASP, Molpro
        call NECICore(0,.False.,.False.,.false.,' ',' ')

    end program NECI

#endif
