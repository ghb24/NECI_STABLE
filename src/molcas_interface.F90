subroutine NECImain(fcidmp, input_name, MemSize, NECIen)

    use constants, only: dp, stdout, int64
    use rdm_finalising, only: RDM_energy
    implicit none
#include "NECICore.h"
    character(*), intent(in) :: fcidmp, input_name
    integer(int64), intent(in) :: MemSize
    real(dp), intent(out) :: NECIen

    write(stdout, *) "STARTING NECI from Molcas"
    call NECICore(call_as_lib=.true., int_name=fcidmp, filename_in=input_name, &
                  MemSize=MemSize)
    ! Once we got excited states energies we will add them here to ENER array.
    NECIen = RDM_energy

end subroutine NECImain
