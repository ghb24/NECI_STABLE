    ! necimain is the entry point for a standalone NECI.  It reads in an
    ! input, and then runs the NECI Core
    program NECI
        implicit none

#include "NECICore.h"

#ifdef DEBUG_
        write(stdout, *) "STARTING NECI"
#endif

        ! Indicate not called by CPMD, VASP, Molpro
        call NECICore()

    end program NECI
