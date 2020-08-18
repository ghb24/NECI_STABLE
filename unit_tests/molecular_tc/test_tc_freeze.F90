program test_tc_freeze
  use constants
  use Parallel_neci, only: MPIInit, MPIEnd, iProcIndex_intra
  use fruit
  use unit_test_helper_excitgen
  use util_mod, only: operator(.isclose.), operator(.div.)
  use sltcnd_mod, only: initSltCndPtr, dyn_sltcnd_excit
  use SystemData, only: t_mol_3_body, nBasis, ECore
  use IntegralsData, only: nFrozen
  use OneEInts, only: TMat2D
  use excitation_types, only: NoExc_t, SingleExc_t, DoubleExc_t, Excitation_t
  use orb_idx_mod, only: SpinProj_t
  use LMat_mod, only: readLMat, freeLMat

  implicit none
  integer, parameter :: n_con = 10
  integer, parameter :: nI_ref(n_con) = (/1,2,3,4,5,6,7,8,9,10/)

  call MPIInit(.false.)
  call init_fruit()
  call tc_freeze_test_driver()
  call fruit_summary()
  call fruit_finalize()
  call MPIEnd(.false.)

contains

    subroutine tc_freeze_test_driver()
        type(NoExc_t) :: exc_0
        type(SingleExc_t) :: exc_1
        
        ! Initialize the matrix element calculation
        call init_excitgen_test(&
            n_el=n_con, fcidump_writer=FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
        t_mol_3_body = .true.
        call initSltCndPtr()

        ! Test the diagonal elements
        call test_freeze((/1,1,2,2,1,1/), exc_0)
        call test_freeze((/1,1,2,1,1,2/), exc_0)
        call test_freeze((/1,2,2,1,1,1/), exc_0)
        call test_freeze((/1,1,3,1,1,3/), exc_0)
        call test_freeze((/1,1,3,1,3,1/), exc_0)
        call test_freeze((/1,3,3,1,1,1/), exc_0)
        call test_freeze((/1,1,1,1,3,3/), exc_0)                
        call test_freeze((/1,2,3,1,2,3/), exc_0)
        call test_freeze((/1,2,3,1,2,3/), exc_0)
        call test_freeze((/1,2,3,2,1,3/), exc_0)
        call test_freeze((/2,2,3,1,1,3/), exc_0)
        call test_freeze((/2,3,2,1,1,3/), exc_0)
        call test_freeze((/1,3,4,1,3,4/), exc_0)
        call test_freeze((/1,3,4,3,1,4/), exc_0)
        call test_freeze((/1,3,4,4,1,3/), exc_0)
        call test_freeze((/3,3,4,1,1,4/), exc_0)
        call test_freeze((/1,3,3,4,1,4/), exc_0)

        ! Test the single excitation matrix elements
        exc_1 = SingleExc_t((/3,11/))
        call test_freeze((/1,1,3,1,1,11/), exc_1)
        call test_freeze((/1,2,3,1,2,11/), exc_1)
        call test_freeze((/1,2,3,2,1,11/), exc_1)
        call test_freeze((/1,2,3,1,11,2/), exc_1)
        call test_freeze((/1,2,3,11,1,2/), exc_1)
        call test_freeze((/1,1,3,11,2,2/), exc_1)
        call test_freeze((/11,2,3,1,1,2/), exc_1)          
        
    end subroutine tc_freeze_test_driver

    subroutine test_freeze(inds, exc)
        integer, intent(in) :: inds(6)
        class(Excitation_t), intent(in) :: exc
        logical :: t_par
        real(dp) :: e_ref, e_freeze

        call reset_ints()
        call write_single_tcdump(inds)
        nFrozen = 0
        call readLMat()
        e_ref = dyn_sltcnd_excit(nI_ref, exc, t_par)
        call freeLMat()
        nFrozen = 4
        call readLMat()
        e_freeze = dyn_sltcnd_excit(nI_ref, exc, t_par) + ECore
        call freeLMat()

        call assert_true(e_freeze .isclose. e_ref)
        write(iout, *) e_freeze, e_ref
    end subroutine test_freeze

    subroutine write_single_tcdump(inds)
        integer, intent(in) :: inds(6)
        integer :: iunit

        iunit = get_free_unit()
        open(iunit, file = "TCDUMP", status = 'replace')

        write(iunit, *), 1.0, inds
    end subroutine write_single_tcdump

    subroutine random_fcidump(iunit)
        integer, intent(in) :: iunit
        call generate_random_integrals(&
            iunit, n_el=n_con, n_spat_orb=14, sparse=0.9_dp, sparseT=0.1_dp, total_ms=SpinProj_t(0))
    end subroutine random_fcidump

    subroutine reset_ints()
        integer :: i
        if(iProcIndex_intra == 0) then
            UMat = 0.0_dp
        end if
        TMat2D = 0.0_dp
        ECore = 0.0_dp
        
    end subroutine reset_ints
  
end program  test_tc_freeze
