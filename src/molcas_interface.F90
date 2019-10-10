subroutine NECImain(fcidmp, input_name, NECIen)

    use constants, only : dp, iout
    use rdm_finalising, only : RDM_energy
    implicit none
#include "NECICore.h"
    character(*), intent(in) :: fcidmp, input_name
    real(dp), intent (out) :: NECIen

    write(iout, *) "STARTING NECI from Molcas"
    call NECICore(call_as_lib=.true., int_name=fcidmp, filename_in=input_name)
    ! Once we got excited states energies we will add them here to ENER array.
    NECIen = RDM_energy

end subroutine NECImain
