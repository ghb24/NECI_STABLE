#ifndef __CMPLX

#include "macros.h"

module guga_rdm
    ! RDM module specifically for the GUGA spin-adapted implementation 

contains 


    subroutine calc_explicit_1_rdm()
        ! routine to calculate the one-RDM explicitly in the GUGA formalism
        ! the one RDM is given by rho_ij = <Psi|E_ij|Psi> 
        ! so, similar to the actHamiltonian routine I need to loop over all 
        ! the orbital indices (i,j) and calculate the action of a given E_ij 
        ! on the sampled wavefunction |Psi> and the resulting overlap with <Psi|

    end subroutine calc_explicit_1_rdm

end module guga_rdm

#endif
