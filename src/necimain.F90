!  necimain is the entry point for a standalone NECI.  It reads in an input, and then runs the NECI Core
Program NECI

    implicit none
    write(6,*) "STARTING NECI"
    call NECICore(0,.False.,.False.) !Indicate not called by CPMD/VASP

End Program NECI
