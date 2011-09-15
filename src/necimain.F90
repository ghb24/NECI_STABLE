!  necimain is the entry point for a standalone NECI.  It reads in an input, and then runs the NECI Core
#ifndef MOLPRO 
Program NECI

    implicit none
    write(6,*) "STARTING NECI"
    call NECICore(0,.False.,.False.,.false.) !Indicate not called by CPMD/VASP/Molpro

End Program NECI
#endif
