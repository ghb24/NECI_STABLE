!> Python interface for initializing the GUGA plugin
subroutine init_guga(fcidump_path, S, nel, norb)
    use guga_plugin, only: init_guga_plugin
    character(*), intent(in) :: fcidump_path
    integer, intent(in) :: S, nel, norb

    !f2py intent(in) :: fcidump_path
    !f2py intent(in) :: S
    !f2py intent(in) :: nel
    !f2py intent(in) :: norb

    call init_guga_plugin(&
        fcidump_path, stot_ = S, nel_ = nel, nbasis_ = norb*2, nSpatOrbs_ = norb)

end subroutine init_guga

!> Python interface for finalizing the GUGA plugin
subroutine clear_guga(err)
    use guga_plugin, only: finalize_guga_plugin

    integer, intent(out) :: err

    !f2py intent(out) :: err

    call finalize_guga_plugin()
    err = 0

end subroutine clear_guga

!> Python interface for getting a guga matrix element
subroutine csf_matel(nIlen, nJlen, nI, nJ, matel)
    use guga_plugin, only: guga_matel
    use SystemData, only: nel
    use util_mod, only: stop_all

    integer, intent(in) :: nIlen, nI(nIlen), nJlen, nJ(nJlen)
    double precision, intent(out) :: matel
    character(25), parameter :: t_r = 'pylib/csf_matel'

    !f2py intent(in) :: nI, nJ
    !f2py intent(hide), depend(nI) :: nIlen = shape(nI, dim = 0)
    !f2py intent(hide), depend(nJ) :: nJlen = shape(nJ, dim = 0)
    !f2py intent(out) :: matel

    if (nIlen /= nel) then
        call stop_all(t_r, 'The number of elements of nI is not the same as '&
            // 'initialized nel')
    else if (nJlen /= nel) then
        call stop_all(t_r, 'The number of elements of nJ is not the same as '&
            // 'initialized nel')
    end if

    matel = guga_matel(nI, nJ)

end subroutine csf_matel

!> Run a full neci calculation with a given permuation of orbitals
!! and return the reference weight
subroutine run_neci(norb, perm, weight)
    use FciMCData, only: fciqmc_run_ref_weight
    use read_fci, only: load_orb_perm
    use SystemData, only: nSpatOrbs
    use util_mod, only: stop_all

    integer, intent(in) :: norb, perm(norb)
    double precision, intent(out) :: weight
    character(25), parameter :: t_r = 'pylib/run_neci'

    !f2py intent(in) :: perm
    !f2py intent(hide), depend(perm) :: norb = shape(perm, dim = 0)
    !f2py intent(out) :: weight

#include "NECICore.h"

    if (norb /= nSpatOrbs) then
        call stop_all(&
            t_r, 'The length of perm is not the same as initialized norb')
    end if

    call load_orb_perm(perm)

    call NECICore(0, .false., .false., .true., filename_in = "neci.inp")

    weight = fciqmc_run_ref_weight

end subroutine run_neci
