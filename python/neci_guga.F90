!> Python interface for initializing the GUGA plugin
subroutine init_guga(S)
!    use guga_plugin, only: init_guga_plugin
    integer, intent(in) :: S
    !f2py intent(in) :: S
!    call init_guga_plugin(stot_ = S)
end subroutine init_guga

!> Python interface for getting a guga matrix element (diagonal)
subroutine diag_matel(nel, nI, matel)
    use guga_plugin, only: guga_matel
    integer, intent(in) :: nel, nI(nel)
    double precision :: matel
    !f2py intent(in) :: nI
    !f2py intent(hide), depend(nI) :: nel = shape(nI, dim = 0)
    !f2py intent(out) :: matel

    matel = guga_matel(nI,nI)
end subroutine diag_matel

!> Run a full neci calculation with a given permuation of orbitals
!! and return the reference weight
subroutine run_neci(norb, perm, weight)
    use FciMCData, only: fciqmc_run_ref_weight
    integer, intent(in) :: norb, perm(norb)
#include "NECICore.h"
    
    call NECICore(0, .false., .false., .true.)
    weight = fciqmc_run_ref_weight    
end subroutine run_neci
