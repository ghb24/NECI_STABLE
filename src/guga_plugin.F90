#include "macros.h"
module guga_plugin
  use constants
  use DetBitOps
  use SystemData
  use guga_excitations
  use guga_matrixElements
  use guga_data
  use guga_init
  use read_fci
  use FciMCData
  use OneEInts, only: TMat2d
  use UMatCache, only: tTransGTID, GetUMatSize, tumat2d, umat2d, tdeferred_umat2d
  use Calc, only: SetCalcDefaults, CalcInit
  use System, only: SetSysDefaults, SysInit
  use OneEInts, only: TMat2d
  use shared_memory_mpi, only: shared_allocate_mpi
  use IntegralsData, only: umat_win, umat
  use DetCalc, only: DetCalcInit
  use LoggingData, only: tRDMonfly, tExplicitAllRDM  
  use shared_memory_mpi, only: shared_allocate_mpi
  use IntegralsData, only: umat_win, umat
  use Parallel_neci, only: MPIInit
  use LoggingData, only: tRDMonfly, tExplicitAllRDM
  use Integrals_neci, only: get_umat_el_normal
  use procedure_pointers, only: get_umat_el
  
  implicit none
  private
  public :: init_guga_plugin, guga_matel

contains

  subroutine init_guga_plugin(stot_, t_testmode_, nel_, nbasis_, &
    nSpatOrbs_)
    integer, intent(in), optional :: stot_
    logical, intent(in), optional :: t_testmode_
    integer, intent(in), optional :: nel_
    integer, intent(in), optional :: nbasis_
    integer, intent(in), optional :: nSpatOrbs_
    logical :: t_testmode
    integer(int64) :: umatsize
    def_default(t_testmode, t_testmode_, .false.)
    def_default(nel, nel_, 0)
    def_default(nbasis, nbasis_, 0)
    def_default(nSpatOrbs, nSpatOrbs_, 0)
    def_default(stot, stot_, 0)    
    call MPIInit(.false.)
    umatsize = 0
    lms = 0
    tGUGA = .true.

    tRDMonfly = .true.
    tFillingStochRDMOnFly = .true.
    call init_bit_rep()

    tGen_sym_guga_mol = .true.
    tgen_guga_weighted = .true.
    tdeferred_umat2d = .true.
    tumat2d = .false.
    ! set this to false before the init to setup all the ilut variables    
    tExplicitAllRDM = .false.

    call init_guga()

    fcidump_name = "FCIDUMP"
    UMatEps = 1.0e-8
    tStoreSpinOrbs = .false.
    tTransGTID = .false.
    tReadFreeFormat = .true.

    call MPIInit(.false.)

    call dSFMT_init(1)

    call SetCalcDefaults()
    call SetSysDefaults()
    tReadInt = .true.

    if(t_testmode) call generate_uniform_integrals()
    
    get_umat_el => get_umat_el_normal

    call initfromfcid(nel,nbasismax,nBasis,lms,.false.)

    call GetUMatSize(nBasis, umatsize)

    allocate(TMat2d(nBasis,nBasis))

    call shared_allocate_mpi(umat_win, umat, (/umatsize/))

    call readfciint(UMat,umat_win,nBasis,ecore,.false.)
    call SysInit()
    ! required: set up the spin info

    call DetInit()
    call DetCalcInit()
    ! call SpinOrbSymSetup()

    call DetPreFreezeInit()

    call CalcInit()
    
  end subroutine init_guga_plugin

  function guga_matel(nI, nJ) result(matel)
    integer, intent(in) :: nI(nel), nJ(nel)
    HElement_t(dp) :: matel
    integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
    type(ExcitationInformation_t) :: excitInfo

    if(all(nI == nJ)) then
      matel = calcDiagMatEleGuga_nI(nI)
    else
      call EncodeBitDet(nI, ilutI)
      call EncodeBitDet(nJ, ilutJ)
      call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, matel, &
      .true., 2)
    end if
    
  end function guga_matel
  
end module guga_plugin
