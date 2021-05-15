!> Python interface for initializing the GUGA plugin
subroutine init_guga(S, nel)
    use guga_plugin, only: init_guga_plugin
    integer, intent(in) :: S
    integer, intent(in) :: nel
    !f2py intent(in) :: S
    !f2py integer :: nel = -1
    call init_guga_plugin(stot_ = S, nel_ = nel)
end subroutine init_guga

!> Python interface for finalizing the GUGA plugin
subroutine clear_guga(err)
    use guga_plugin, only: finalize_guga_plugin
    integer, intent(out) :: err
    !f2py intent(out) :: err

    call finalize_guga_plugin()
    err = 0
end subroutine clear_guga

!> Python interface for getting a guga matrix element (diagonal)
subroutine diag_matel(nel, nI, matel)
    use guga_plugin, only: guga_matel
    integer, intent(in) :: nel, nI(nel)
    double precision, intent(out) :: matel
    !f2py intent(in) :: nI
    !f2py intent(hide), depend(nI) :: nel = shape(nI, dim = 0)
    !f2py intent(out) :: matel

    matel = guga_matel(nI,nI)
end subroutine diag_matel

!> Run a full neci calculation with a given permuation of orbitals
!! and return the reference weight
subroutine run_neci(norb, perm, weight)    
    use FciMCData, only: fciqmc_run_ref_weight
    use read_fci, only: load_orb_perm
    integer, intent(in) :: norb, perm(norb)
    double precision, intent(out) :: weight
    !f2py intent(in) :: perm
    !f2py intent(hide), depend(perm) :: norb = shape(perm, dim = 0)
    !f2py intent(out) :: weight
#include "NECICore.h"
    call load_orb_perm(perm)
    call NECICore(0, .false., .false., .true., filename_in = "neci.inp")
    weight = fciqmc_run_ref_weight    
end subroutine run_neci
