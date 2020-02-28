#include "macros.h"
MODULE System

    use SystemData
    use CalcData, only: TAU, tTruncInitiator, InitiatorWalkNo, &
                        occCASorbs, virtCASorbs, tPairedReplicas, tInitializeCSF, &
                        S2Init, tDynamicAvMcEx
    use FciMCData, only: tGenMatHEl, t_initialized_roi
    use sort_mod
    use SymExcitDataMod, only: tBuildOccVirtList, tBuildSpinSepLists
    use constants
    use iso_c_hack
    use read_fci, only: FCIDUMP_name
    use util_mod, only: error_function, error_function_c,&
        near_zero, operator(.isclose.)
    use lattice_mod, only: lattice, lat
    use k_space_hubbard, only: setup_symmetry_table
    use breathing_Hub, only: setupMomIndexTable, setupBreathingCont
    use tc_three_body_data, only: LMatEps, tSparseLMat
    use ParallelHelper, only: iprocindex, root

    IMPLICIT NONE

    contains

    subroutine SetSysDefaults()
      != Set defaults for Calc data items.

      use default_sets
      USE SymData, only: tAbelianFastExcitGen
      USE SymData, only: tStoreStateList
      use OneEInts, only: tOneElecDiag
      implicit none

      ! Default from SymExcitDataMod
      tBuildOccVirtList = .false.
!     SYSTEM defaults - leave these as the default defaults
!     Any further addition of defaults should change these after via
!     specifying a new set of DEFAULTS.
      tReltvy = .false.
      tComplexOrbs_RealInts = .false.
      tComplexWalkers_RealInts = .false.
      tReadFreeFormat=.true.
      tMolproMimic=.false.
      tNoSingExcits=.false.
      tOneElecDiag=.false.
      tMCSizeTruncSpace=.false.
      iMCCalcTruncLev=0
      tOddS_HPHF=.false.
      tRotatedOrbsReal=.false.  !This is set if compiled in real, but reading in a complex FCIDUMP.
      tISKFuncs=.false.       !This is for kpoint symmetry with inversion so that determinants can be combined.
      tKPntSym=.false.        !This is for k-point symmetry with the symrandexcit2 excitation generators.
      tNoSinglesPossible = .false.
      t_mol_3_body = .false.
      tMCSizeSpace=.false.
      CalcDetPrint=1000
      CalcDetCycles=10000
      tFixLz=.false.
      tStoreSpinOrbs=.false.    !by default we store/lookup integrals as spatial integrals
      tNoBrillouin=.true.
      tBrillouinsDefault=.true.
      tROHF=.false.
      tCacheFCIDUMPInts=.false.
      tHPHFInts=.false.
      tHPHF=.false.
      tMaxHLGap=.false.
      UMatEps = 1.0e-8
      LMatEps = 1.0e-10
      tExactSizeSpace=.false.
      iRanLuxLev=3      !This is the default level of quality for the random number generator.
      tNoSymGenRandExcits=.false.
      tNonUniRandExcits=.true.
      tCycleOrbs=.false.
      TSTARBIN=.false.
      TREADINT=.false.
      THFORDER=.false.
      TDFREAD=.false.
      tRIIntegrals=.false.
      TPBC=.false.
      TCPMD=.false.
      tVASP=.false.
      THUB=.false.
      TUEG=.false.
      tUEG2=.false.
      tHeisenberg = .false.
      tLatticeGens =.false.
      tNoFailAb=.false.
      LMS=0
      TSPN=.false.
      TCSF=.false.
      tInitializeCSF = .false.
      TCSFOLD = .false.
      csf_trunc_level = 0
      tTruncateCSF = .false.
      STOT=0
      TPARITY = .false.
      IParity(:)=0
      dimen = 3
      NMAXX = 0
      NMAXY = 0
      NMAXZ = 0
      NMSH=32
      BOX=1.0_dp
      BOA=1.0_dp
      COA=1.0_dp
      TUSEBRILLOUIN=.false.
      tUHF=.false.
      FUEGRS=0.0_dp
      iPeriodicDampingType=0
      fRc=0.0_dp
      TEXCH=.true.
      UHUB = 4
      BHUB = -1
      btHub = 0.0_dp
      TREAL = .false.
      tUEGTrueEnergies = .false.
      tUEGSpecifyMomentum = .false.
      tUEGOffset = .false.
      TTILT = .false.
      TALPHA = .false.
      ALPHA = 0.0_dp
      ISTATE = 1
      OrbECutoff=1e20_dp
      tOrbECutoff=.false.
      gCutoff=1e20_dp ! This shouldn't be used
      tgCutoff=.false.
      tMadelung=.false.
      Madelung=0.0_dp
      tUEGFreeze=.false.
      FreezeCutoff=1e20_dp ! This shouldn't be used
      tMP2UEGRestrict=.false.
      tStoreAsExcitations=.false.
      TBIN=.false.
      tAbelianFastExcitGen=.true.
      TNoRenormRandExcits=.false.
      tStoreStateList=.false.       !This will be turned to true by default if not in abelian symmetry
      tAssumeSizeExcitgen=.false.
      tLagrange=.false.
      tShake=.false.
      lNoSymmetry=.false.
      tRotateOrbs=.false.
      TimeStep=0.01_dp
      ConvergedForce=0.001_dp
      ShakeConverged=0.001_dp
      tShakeApprox=.false.
      tShakeDelay=.false.
      ShakeStart=1
      tROIteration=.false.
      ROIterMax=5000
      tShakeIter=.false.
      ShakeIterMax=5
      tSeparateOccVirt=.false.
      tOffDiagSqrdMin=.false.
      tOffDiagSqrdMax=.false.
      tOffDiagMin=.false.
      tOffDiagMax=.false.
      tDoubExcMin=.false.
      tHijSqrdMin=.false.
      tOneElIntMax=.false.
      tOnePartOrbEnMax=.false.
      tHFSingDoubExcMax=.false.
      tERLocalization=.false.
      tVirtCoulombMax=.false.
      tVirtExchangeMin=.false.
      tRotateOccOnly=.false.
      tRotateVirtOnly=.false.
      tSpinOrbs=.false.
      tReadInCoeff=.false.
      tUseMP2VarDenMat=.false.
      tUseHFOrbs=.false.
      tFindCINatOrbs=.false.
      DiagWeight=1.0_dp
      OffDiagWeight=1.0_dp
      OneElWeight=1.0_dp
      OrbEnMaxAlpha=1.0_dp
      MaxMinFac=1
      DiagMaxMinFac=-1
      OneElMaxMinFac=1
      tDiagonalizehij=.false.
      tHFNoOrder=.false.
      tSymIgnoreEnergies=.false.
      tPickVirtUniform = .false.
      modk_offdiag = .false.
      tAllSymSectors = .false.
      tGenHelWeighted = .false.
      tGen_4ind_weighted = .false.
      tGen_4ind_part_exact = .false.
      tGen_4ind_lin_exact = .false.
      tGen_4ind_reverse = .false.
      tUEGNewGenerator = .false.
      tGen_4ind_2 = .false.
      tGen_4ind_2_symmetric = .false.
      tmodHub = .false.
      t_uniform_excits = .false.
      t_mol_3_body = .false.
      t_ueg_3_body = .false.
      t_ueg_transcorr = .false.
      t_trcorr_gausscutoff = .false.
      t_ueg_dump = .false.
      t_exclude_3_body_excits = .false.
      t_pcpp_excitgen = .false.
      t_pchb_excitgen = .false.
      ! use weighted singles for the pchb excitgen?
      t_pchb_weighted_singles = .false.
      tMultiReplicas = .false.
      tGiovannisBrokenInit = .false.
      ! GAS options
      tGASSpinRecoupling = .false.
      tGAS = .false.

#ifdef PROG_NUMRUNS_
      ! by default, excitation generation already creates matrix elements
      tGenMatHEl = .true.
      tInfSumTCCalc= .false.
      tInfSumTCPrint= .false.
      tInfSumTCRead= .false.
      PotentialStrength=1.0_dp
      TranscorrCutoff=0
      TranscorrIntCutoff=0
      TranscorrGaussCutoff=1.d0
      TContact=.false.
      TUnitary=.false.
      Tperiodicinmom=.false.
      t12FoldSym = .false.
      t_initialized_roi = .false.

      inum_runs = 1
#ifdef CMPLX_
      lenof_sign = 2
#else
      lenof_sign = 1
#endif
#endif

!Feb08 defaults:
      IF(Feb08) THEN
          !...add defaults...
      ENDIF

! Coulomb damping function currently removed.
      FCOULDAMPBETA=-1.0_dp
      COULDAMPORB=0
      k_offset = 0.0_dp

    end subroutine SetSysDefaults

    SUBROUTINE SysReadInput()
      USE input_neci
      USE SymData, only: tAbelianFastExcitGen
      USE SymData, only: tStoreStateList
      use OneEInts, only: tOneElecDiag
      IMPLICIT NONE
      LOGICAL eof
      CHARACTER (LEN=100) w
      INTEGER I,Odd_EvenHPHF,Odd_EvenMI
      integer :: ras_size_1, ras_size_2, ras_size_3, ras_min_1, ras_max_3, itmp
      character(*), parameter :: t_r = 'SysReadInput'

      ! The system block is specified with at least one keyword on the same
      ! line, giving the system type being used.
      call readu(w)
      select case(w)

      case("DFREAD")             ! Instead, specify DensityFitted within the system block.
          TREADINT = .true.
          TDFREAD = .true.
          call readu(w)
          select case(w)
          case("ORDER")
              THFORDER = .true.
          case("NOORDER")
              THFNOORDER = .true.
          end select
      case("BINREAD")            ! Instead, specify Binary within the system block.
          TREADINT=.true.
          TBIN=.true.
          call readu(w)
          select case(w)
          case("ORDER")
              THFORDER=.true.
          case("NOORDER")
              THFNOORDER = .true.
          end select
      case("READ","GENERIC")
          TREADINT = .true.
          call readu(w)
          select case(w)
          case("ORDER")
              THFORDER = .true.
          case("NOORDER")
              THFNOORDER = .true.
          end select

      case("HUBBARD")
          THUB = .true.
          TPBC=.true.

          if (item < nitems) then
              ! this indicates the new hubbard implementation
              ! for consistency turn off the old hubbard indication
              ! and for now this is only done for the real-space hubbard
              ! not for the k-space implementation todo
              ! and do i need to turn of tpbc also? try
              ! use the already provided setup routine and just modify the
              ! necessary stuff, like excitation generators!
              t_new_hubbard = .true.
              call readl(w)
              select case (w)
              case ('real-space','real')
                  treal = .true.
                  t_new_real_space_hubbard = .true.
                  t_lattice_model = .true.

                  ! if no further input is given a provided fcidump is
                  ! assumed! but this still needs to be implemented
                  ! this fcidump gives the lattice structure!
                  if (item < nitems) then
                      call readl(lattice_type)
                  else
                      lattice_type = 'read'
                  end if
              case ('momentum-space','k-space','momentum')
                  ! reuse most of the old initialisation for the k-space
                  ! hubbard. one has to be really careful to initialize all
                  ! the correct stuff especially for the matrix element
                  ! calculation with the HPHF option turned on!
                  t_k_space_hubbard = .true.
                  t_lattice_model = .true.
                  tKPntSym = .true.

              case default
                  print *, w
                  call Stop_All(t_r, "not recognised keyword!")
              end select
          end if

      case ('TJ','TJ-MODEL')
          t_tJ_model = .true.
          t_lattice_model = .true.
          ! misuse the hubbard initialisation
          tHub = .true.
          tpbc = .true.
          treal = .true.

      case ('HEISENBERG')
          ! should i misuse the already provided setup for the hubbard
          ! model again? .. maybe..
          ! maybe i should use a general flag like t_lattice_model
          ! especially for the matrix element evaluation and stuff
          t_heisenberg_model = .true.
          t_lattice_model = .true.
          thub = .true.
          tpbc = .true.
          treal = .true.

      case("RIINTEGRALS")
          tRIIntegrals = .true.
          tReadInt=.true.
      case("UEG")
          TUEG = .true.
          tOneElecDiag=.true.   !One electron integrals diagonal
      case("VASP")
          tVASP= .true.
      case("CPMD")
          TCPMD = .true.
          call readu(w)
          select case(w)
          case("ORDER")
              THFORDER = .true.
          end select
      case("BOX")
          tOneElecDiag=.true.   !One electron integrals diagonal
      case default
          call report("System type "//trim(w)//" not valid",.true.)
      end select

      ! Now parse the rest of the system block.
system: do
        call read_line(eof)
        if (eof) then
            call report("Incomplete input file",.true.)
        end if
        call readu(w)
        select case(w)

        ! Options for molecular (READ) systems: control how the integral file
        ! is read in.
        case("BINARY")
            tBin=.true.
        case("DENSITYFITTED")
            tDFRead=.true.
        ! General options.
        case("RIINTEGRALS")
           tRIIntegrals = .true.
        case("READCACHEINTS")
           tCacheFCIDUMPInts=.true.
        case("ELECTRONS","NEL")
            call geti(NEL)
        case("SPIN-RESTRICT")
            if(item.lt.nitems) then
               call geti(LMS)
            else
               LMS=0
            endif
            TSPN = .true.
        case("INITIAL-SPIN")
           call getf(S2Init)
           tInitializeCSF = .true.
        case("CSF")
            if(item.lt.nitems) then
               call geti(STOT)
            else
               STOT=0
            endif
            TCSF = .true.
        case("TRUNCATE-CSF")
            if (item < nitems) then
                call geti(csf_trunc_level)
            else
                csf_trunc_level = 0
            endif
            tTruncateCSF = .true.
        case("CSF-OLD")
            if(item.lt.nitems) then
               call geti(STOT)
            else
               STOT=0
            endif
            TCSFOLD = .true.
        case("SYMIGNOREENERGIES")
            tSymIgnoreEnergies=.true.
        case("NOSYMMETRY")
            lNoSymmetry=.true.
            IF(tHub) THEN
                CALL Stop_All("SysReadInput","Cannot turn off symmetry with the hubbard model.")
            ENDIF

        case("FREEFORMAT")
            ! Relax the formatting requirements for reading in FCIDUMP files.
            !
            ! For historical reasons, QChem uses a very fixed format for
            ! outputting FCIDUMP files. As a result the columns of orbital
            ! indices will merge whenever there are more than 99 spatial
            ! orbitals. To correctly read these files a FIXED format is
            ! required for reading. Obviously, this is non-ideal when reading
            ! FCIDUMP formats from elsewhere.
            !
            ! The QChem behaviour used to be default, but this has been
            ! deprecated. To obtain the fixed behaviour use
            ! "FREEFORMAT OFF" or "FREEFORMAT FALSE"

            tReadFreeFormat = .true.
            if (item < nitems) then
                call readu(w)
                select case(w)
                case("OFF", "FALSE")
                    tReadFreeFormat = .false.
                case default
                end select
            end if

        case("SYM")
            TPARITY = .true.
            do I = 1,4
              call geti(IPARITY(I))
            end do
! the last number is the symmetry specification - and is placed in position 5
            IPARITY(5)=IPARITY(4)
            IPARITY(4)=0
            tSymSet = .true.
        case("USEBRILLOUINTHEOREM")
          TUSEBRILLOUIN=.TRUE.
          tNoBrillouin=.false.
          tBrillouinsDefault=.false.
        case("NOBRILLOUINTHEOREM")
            tNoBrillouin=.true.
            tBrillouinsDefault=.false.
        case("UHF")
! This keyword is required if we are doing an open shell calculation
! but do not want to include singles in the energy calculations.
            tUHF=.true.
        case("RS")
            call getf(FUEGRS)
        case("EXCHANGE-CUTOFF")
            iPeriodicDampingType=2
            if(item.lt.nitems) then
               call getf(fRc)
            endif
        case("EXCHANGE-ATTENUATE")
            iPeriodicDampingType=1
            if(item.lt.nitems) then
               call getf(fRc)
            endif
        case("EXCHANGE")
            call readu(w)
            select case(w)
               case("ON")
                  TEXCH=.TRUE.
               case("OFF")
                  TEXCH=.FALSE.
               case default
                  call report("EXCHANGE "//trim(w)//" not valid",.true.)
            end select
        case("COULOMB")
            call report("Coulomb feature removed",.true.)
!            call getf(FCOUL)
        case("COULOMB-DAMPING")
            call report("Coulomb damping feature removed",.true.)
!            call readu(w)
!            select case(w)
!            case("ENERGY")
!               call getf(FCOULDAMPMU)
!               call getf(FCOULDAMPBETA)
!            case("ORBITAL")
!               call geti(COULDAMPORB)
!               call getf(FCOULDAMPBETA)
!            end select
        case("ENERGY-CUTOFF")
          tOrbECutoff=.true.
          call getf(OrbECutoff)
        case("G-CUTOFF")
          tgCutoff=.true.
          call getf(gCutoff)
        case("FREEZE-CUTOFF")
            tUEGFreeze=.true.
            call getf(FreezeCutoff)
        case("MADELUNG")
          tMadelung=.true.
          call getf(Madelung)
        case("UEG2")
           tUEG2 = .true.
        case("STORE-AS-EXCITATIONS")
           tStoreAsExcitations=.true.
        case("MP2-UEG-RESTRICT")
           tMP2UEGRestrict=.true.
           call geti(kiRestrict(1))
           call geti(kiRestrict(2))
           call geti(kiRestrict(3))
           call geti(kiMsRestrict)
           call geti(kjRestrict(1))
           call geti(kjRestrict(2))
           call geti(kjRestrict(3))
           call geti(kjMsRestrict)

        ! Options for model systems (electrons in a box/Hubbard).
        case("CELL")
            call geti(NMAXX)
            call geti(NMAXY)
            call geti(NMAXZ)

            ! misuse the cell keyword to set this up to also have the
            ! hubbard setup already provided..
!             if (t_new_real_space_hubbard) then
!                length_x = NMAXX
!                length_y = NMAXY
!            end if

       case('SPIN-TRANSCORR')
           ! make a spin-dependent transcorrelation factor
           t_spin_dependent_transcorr = .true.
           if (item < nitems) then
               call getf(trans_corr_param)
           else
               trans_corr_param = 0.1_dp
           end if
           t_non_hermitian = .true.

        case("NONHERMITIAN")
           ! just use a non-hermitian Hamiltonian, no additional tweaks
           t_non_hermitian = .true.

        case('MOLECULAR-TRANSCORR')
            t_non_hermitian = .true.
            ! optionally supply the three-body integrals of the TC Hamiltonian
            t_3_body_excits = .true.
            if(item < nitems) then
               call readu(w)
               select case(w)
               case("3-BODY")
                  t_mol_3_body = .true.
                  max_ex_level = 3
                  ! this uses a uniform excitation generator, switch off matrix
                  ! element computation for HPHF
                  tGenMatHEl = .false.
               case("UEG")
                  t_mol_3_body = .true.
                  tGenMatHEl = .false.
                  t12FoldSym = .true.
                  tNoSinglesPossible = .true.
               case default
                   t_mol_3_body = .true.
                   tGenMatHEl = .false.
               end select
            endif
            if(t_mol_3_body) max_ex_level = 3

       case('UEG-TRANSCORR')
           t_ueg_transcorr = .true.
           t_non_hermitian = .true.
            do while(item < nitems)
               call readu(w)
               select case(w)
               case("3-BODY")
                  tTrcorrExgen = .false.
                  tTrCorrRandExgen = .true.
                  t_ueg_3_body = .true.
                  tGenMatHEl = .false.
                  max_ex_level = 3
                  tRPA_tc= .false.

               case("TRCORR-EXCITGEN")
                  tTrcorrExgen = .true.
                  tTrCorrRandExgen = .false.

               case("RAND-EXCITGEN")
                  tTrCorrRandExgen = .true.

!              case default
!                 t_ueg_3_body = .false.
!                 tTrcorrExgen = .true.
!                 tTrCorrRandExgen = .false.


               end select
!               write(6,*) tTrcorrExgen, tTrCorrRandExgen, t_ueg_3_body
            enddo
!               call stop_all('debug stop')

       case('UEG-DUMP')
           t_ueg_dump = .true.

       case('EXCLUDE-3-BODY-EX')
          ! Do not generate 3-body excitations, even in the molecular-transcorr mode
          t_exclude_3_body_excits = .true.

       case ('TRANSCORRELATED', 'TRANSCORR', 'TRANS-CORR')
           ! activate the transcorrelated Hamiltonian idea from hongjun for
           ! the real-space hubbard model
           t_trans_corr = .true.
           t_non_hermitian = .true.

           if (item < nitems) then
               call getf(trans_corr_param)
           else
               ! defaul value 1 for now, since i have no clue how this behaves
               trans_corr_param = 1.0_dp
           end if

        case ("TRANSCORR-NEW")
            t_trans_corr = .true.
            t_trans_corr_new = .true.
           t_non_hermitian = .true.

           if (item < nitems) then
               call getf(trans_corr_param)
           else
               ! defaul value 1 for now, since i have no clue how this behaves
               trans_corr_param = 1.0_dp
           end if

        case ('2-BODY-TRANSCORR', '2-BODY-TRANS-CORR', '2-BODY-TRANSCORRELATED','TRANSCORR-2BODY')
           ! for the tJ model there are 2 choices of the transcorrelation
           ! indicate that here!
           t_trans_corr_2body = .true.
           t_non_hermitian = .true.

           if (item < nitems) then
               call getf(trans_corr_param_2body)

           else
               trans_corr_param_2body = 0.25_dp
           end if

           ! if it is the k-space hubbard also activate 3-body excitations here
           if (t_k_space_hubbard) then
               t_3_body_excits = .true.
               max_ex_level = 3
           end if

        case ('NEIGHBOR-TRANSCORR','TRANSCORR-NEIGHBOR','N-TRANSCORR')
            t_trans_corr_2body = .true.
            t_non_hermitian = .true.

            if (item < nitems) then
                call getf(trans_corr_param_2body)
            else
                trans_corr_param_2body = 0.25_dp
            end if

           ! if it is the k-space hubbard also activate 3-body excitations here
           if (t_k_space_hubbard) then
               t_3_body_excits = .true.
               max_ex_level = 3
           end if

        case ("TRANSCORR-HOP","HOP-TRANSCORR")
            t_trans_corr_hop = .true.
            t_non_hermitian = .true.

            if (item < nitems) then
                call getf(trans_corr_param)
            else
                trans_corr_param = 0.5_dp
            end if


       ! Options for the type of the reciprocal lattice (eg sc, fcc, bcc)
        case("REAL_LATTICE_TYPE")
            call readl(real_lattice_type)
       ! Options for the dimension (1, 2, or 3)
         case("DIMENSION")
            call geti(dimen)

        ! Options for transcorrelated method (only: UEG 2D 3D, Homogeneous 1D 3D
        ! gas with contact interaction)
         case("TRANSCORRCUTOFF")
           if(item < nitems) then
                call readu(w)
                select case(w)
                  case("GAUSS")
                    if(dimen.ne.1) stop 'Gauss cutoff is developed only for 1D!'
                    call getf(TranscorrGaussCutoff)
                    t_trcorr_gausscutoff = .true.

                  case("STEP")
                    t_trcorr_gausscutoff = .false.
                    call geti(TranscorrCutoff)
                    if(tContact.and.dimen.eq.3) then
                       tInfSumTCCalc=.true.
                       call geti(TranscorrIntCutoff)
                    endif

                end select

           endif

       !Options for turn off the RPA term(only: transcorrelated homogeneous 1D
       !gas with contact interaction), tRPA_tc is set to true for two particles,
       ! but it is turned off, if 3-body interactions are used
         case("NORPATC")
           tRPA_tc=.false.

         case("PERIODICINMOMSPACE")
           Tperiodicinmom=.true.

       ! Contact interaction for homogenous one dimensional Fermi gas is applied
         case("CONTACTINTERACTION")
           tContact=.true.
           call getf(PotentialStrength)
           if(dimen.ne.1) &
                stop 'Contact interaction only for 1D!'

       ! Contact interaction forthe three dimensional Fermi gas in the unitary
       ! interaction
         case("CONTACTINTERACTIONUNITARY")
           tContact=.true.
           tUnitary=.true.
           if(dimen.ne.3) stop 'Unitary regime only for 3D!'

       ! Option for infinite summation at transcorrelated method for homogeneous
       ! 3D gas with contact interaction
         case("TRANSCORRINFSUM")
           if(item < nitems) then
                call readu(w)
                select case(w)
                        case("CALC")
                                tInfSumTCCalc=.true.

                        case("CALCANDPRINT")
                                tInfSumTCCalc=.true.
                                tInfSumTCPrint=.true.

                        case("READ")
                                tInfSumTCCalc=.false.
                                tInfSumTCRead=.true.
                end select
           endif



        ! This means that no a is generated when b would be made and rejected
        ! O(N^2) loop makes this a poor choice for larger systems.
        case("UEGNOFAIL")
            tNoFailAb = .true.
        ! These are the new lattice excitation generators that conserve momentum
        ! during excitation generation for efficiency
        case("LATTICE-EXCITGEN")
            tLatticeGens =.true.
        ! use the simplified random excitation generator for k-space hubbard that
        ! does not use a cumulative list (it is much more efficient)
        case("UNIFORM-EXCITGEN")
            t_uniform_excits = .true.
        case("MESH")
            call geti(NMSH)
        case("BOXSIZE")
            call getf(BOX)
            if(item.lt.nitems) then
               call getf(BOA)
               call getf(COA)
            else
               BOA=1.0_dp
               COA=1.0_dp
            endif
        case("U")
            call getf(UHUB)
        case("B")
            call getf(BHUB)

        case ("J")
            ! specify the tJ exchange here, the default is 1.0
            ! this could also be used for the heisenberg model..
            call getf(exchange_j)
         case("C")
            call getf(btHub)
            tmodHub = .true.

        case ("NEXT-NEAREST-HOPPING")
            call getf(nn_bhub)

        case("REAL")
            TREAL = .true.
            ! in case of the real-space lattice also turn off symmetries
            lNoSymmetry = .true.
        case("APERIODIC")
            TPBC = .false.


        case("TWISTED-BC")
            ! use of twisted boundary conditions for the cubic and tilted
            ! hubbard lattice model
            t_twisted_bc = .true.
            call getf(twisted_bc(1))
            if (item < nitems) then
                call getf(twisted_bc(2))
            else
                ! if only one input apply same twist in x and y direction
                twisted_bc(2) = twisted_bc(1)
            endif

        case ("OPEN-BC")
            ! open boundary implementation for the real-space hubbard
            ! model
            ! only applicable for 1 and 2 dimensional lattices!
            ! for the square lattice on can specify mixed boundary conditions
            ! (periodic + open) but not in the tilted where always full
            ! open boundary conditions are used if this keyword is present!

            if (item < nitems) then
                call readu(w)

                select case (w)
                case ("X")
                    ! onl open in x-direction
                    t_open_bc_x = .true.

                case ("Y")
                    ! only open in y-direction
                    t_open_bc_y = .true.

                case ("XY")
                    ! open in both directions
                    t_open_bc_x = .true.
                    t_open_bc_y = .true.

                end select
            else
                t_open_bc_x = .true.
                t_open_bc_y = .true.
            end if

        case("LATTICE")
            ! new hubbard implementation
            ! but maybe think of a better way to init that..
            ! the input has to be like:
            ! lattice [type] [len_1] [*len_2]
            ! where length to is optional if it is necessary to input it.
!             tHub = .false.
!             treal = .false.
!             lNoSymmetry = .true.
            ! this treal is not true.. now we also have k-space hubbard lattice
            ! support
!             treal = .true.
!             t_new_real_space_hubbard = .true.

            ! set some defaults:
            lattice_type = "read"

            length_x = -1
            length_y = -1

!             tPBC = .false.

            if (item < nitems) then
               ! use only new hubbard flags in this case
               call readl(lattice_type)
            end if

            if (item < nitems) then
                call geti(length_x)
            end if

            if (item < nitems) then
                call geti(length_y)
            end if

            if (item < nitems) then
                call geti(length_z)
            end if

            if (t_k_space_hubbard) then
                lat => lattice(lattice_type, length_x, length_y, length_z, &
                    .not. t_open_bc_x, .not. t_open_bc_y, .not. t_open_bc_z,'k-space')
            else if (t_new_real_space_hubbard) then
                lat => lattice(lattice_type, length_x, length_y, length_z, &
                    .not. t_open_bc_x, .not. t_open_bc_y, .not. t_open_bc_z,'real-space')
            else
                lat => lattice(lattice_type, length_x, length_y, length_z, &
                    .not. t_open_bc_x, .not. t_open_bc_y, .not. t_open_bc_z)

            end if

            ! maybe i have to reuse the cell input functionality or set it
            ! here also, so that the setup is not messed up
            ! i should set those quantities here again..
            nmaxx = length_x
            nmaxy = length_y
            nmaxz = 1

        case("UEG-OFFSET")
            tUEGOffset=.true.
            call getf(k_offset(1))
            call getf(k_offset(2))
            call getf(k_offset(3))
        case("UEG-SCALED-ENERGIES")
            tUEGTrueEnergies=.true.
        case("UEG-MOMENTUM")
            tUEGSpecifyMomentum=.true.
            call geti(k_momentum(1))
            call geti(k_momentum(2))
            call geti(k_momentum(3))
        case("TILT")
            TTILT = .true.
            call geti(ITILTX)
            call geti(ITILTY)
        case("ALPHA")
            TALPHA = .true.
            call getf(ALPHA)
        case("STATE")
            call geti(ISTATE)
            if ( ISTATE /= 1 ) then
                call report("Require ISTATE to be left set as 1",.true.)
            end if
        case("MODK-OFFDIAG")
            modk_offdiag = .true.
        case("FAST-EXCITGEN")
            tAbelianFastExcitGen=.true.
    ! tAbelianFastExcitGen is a temporary flag.
    !  It is used to indicate that if an Abelian symmetry group is presents
    !   the excitation generators should use optimized routines
    !   to take this into account.  Not all excitation generator functions
    !   currently work with this.  USE WITH CARE
            if (item.lt.nitems) then
                call readu(w)
                select case(w)
                case("OFF")
                   tAbelianFastExcitGen=.false.
                end select
            end if
        case("NORENORMRANDEXCITS")
!    Since we have already calculated the number of excitations possible for each symmetry type, there
!    no need to renormalise all excitations with weight 1. As long as pairs of allowed occupied and
!    virtual orbitals can be chosen without any bias, then we can generate random excitations in O[1] time.
!    This is default off since it will change previous results, however it is strongly recommended to be
!    on for virtually all unweighted MC calculations, since it should speed up generation, especially in
!    low symmetry and/or large systems.
            tNoRenormRandExcits=.true.
        case("STORESTATELIST")
!This flag indicates that we want to store the full list of symmetry state pairs.
!This is done by default in non-abelian symmetries, but in abelian symmetries, will
!only be done if specified. This may be wanted, since it means that random excitations
!will be created quicker currently, since NoRenormRandExcits can only work with it on in
!abelian symmetry.
            tStoreStateList=.true.
        case("ASSUMESIZEEXCITGEN")
!This flag indicates that we want to setup the excitations faster, by returning a maximum size of the
!generators instantly in the first setup, and then we do not need to run through the excitations twice.
!However, certain bits of the generator are not stored, namely nAllowPPS and SymProds, as well as the
!iterator information. This means that we can only randomly generate excitations, and any attempt to
!use these assumed size excitation generators to generate the whole list of excitations, will result
!in bad, bad times.
            tAssumeSizeExcitgen=.true.
        case("MOMINVSYM")
            call stop_all(t_r,'Deprecated function. Look in defunct_code folder if you want to see it')
        case("HPHF")
            tHPHF=.true.
            if(item.lt.nitems) then
                call geti(Odd_EvenHPHF)
                if(Odd_EvenHPHF.eq.1) then
                    !Want to converge onto an Odd S State
                    tOddS_HPHF=.true.
                elseif(Odd_EvenHPHF.eq.0) then
                    !Want to converge onto an Even S State
                    !tOddS_HPHF should be false by default.
                else
                    call stop_all("SysReadInput","Invalid variable given to HPHF option: 0 = Even S; 1 = Odd S")
                endif
            endif
        case("ROTATEORBS")
! The ROTATEORBS calculation initiates a routine which takes the HF orbitals
! and finds the optimal set of transformation coefficients to fit a particular criteria specified below.
! This new set of orbitals can then used to produce a ROFCIDUMP file and perform the FCIMC calculation.
            tRotateOrbs=.true.
            if (item.lt.nitems) then
                call Getf(TimeStep)
                call Getf(ConvergedForce)
            end if
! The SHAKE orthonormalisation algorithm is automatically turned on with a default of 5 iterations.
            tShake=.true.
            tShakeIter=.true.

        case("DIAGONALIZEHIJ")
! Instead of doing an orbital rotation according to the P.E's below, this keyword sets the rotate orbs routine to
! take the input <i|h|j>, and diagonalise it to give the rotation coefficients.
            tDiagonalizehij=.true.

        case("HFSINGDOUBEXCMAX")
! This maximises the square of the four index integrals corresponding to single and double excitations from
! the HF.
            tHFSingDoubExcMax=.true.
            MaxMinFac=-1

        case("OFFDIAGSQRDMIN")
            tOffDiagSqrdMin=.true.
            MaxMinFac=1
            IF(item.lt.nitems) THEN
                call Getf(OffDiagWeight)
            ELSE
                OffDiagWeight=1.0
            ENDIF
! This sets the orbital rotation to find the coefficients which minimise the sum of the |<ij|kl>|^2 elements.
! The following real value sets the importance of the minimisation relative to the ER maximisation, providing
! the ERLOCALIZATION keyword is also present.  This is the same for all OFFDIAG keywords below.

        case("OFFDIAGSQRDMAX")
            tOffDiagSqrdMax=.true.
            MaxMinFac=-1
            IF(item.lt.nitems) THEN
                call Getf(OffDiagWeight)
            ELSE
                OffDiagWeight=1.0
            ENDIF
! This sets the orbital rotation to find the coefficients which maximise the sum of the |<ij|kl>|^2 elements.

        case("OFFDIAGMIN")
            tOffDiagMin=.true.
            MaxMinFac=1
            IF(item.lt.nitems) THEN
                call Getf(OffDiagWeight)
            ELSE
                OffDiagWeight=1.0
            ENDIF
! This sets the orbital rotation to find the coefficients which minimise the sum of the <ij|kl> elements.

        case("OFFDIAGMAX")
            tOffDiagMax=.true.
            MaxMinFac=-1
            IF(item.lt.nitems) THEN
                call Getf(OffDiagWeight)
            ELSE
                OffDiagWeight=1.0
            ENDIF
! This sets the orbital rotation to find the coefficients which maximise the sum of the <ij|kl> elements.

        case("DOUBEXCITEMIN")
            tDoubExcMin=.true.
            MaxMinFac=1
            IF(item.lt.nitems) THEN
                call Getf(OffDiagWeight)
            ELSE
                OffDiagWeight=1.0
            ENDIF
! This sets the orbital rotation to find the coefficients which minimise the double excitation hamiltonian elements.

        case("HIJSQRDMIN")
            tHijSqrdMin=.true.
            OneElMaxMinFac=1
            tRotateVirtOnly=.true.
            IF(item.lt.nitems) THEN
                call Getf(OneElWeight)
            ELSE
                OneElWeight=1.0
            ENDIF
! This sets the orbital rotation to find the coefficients which minimis the square of the one electron integrals, <i|h|j>.
! i can be occupied or virtual, j only virtual, i=<j.

        case("ONEELINTMAX")
            tOneElIntMax=.true.
            MaxMinFac=-1
            tRotateVirtOnly=.true.
! This sets the orbital rotation to find the coefficients which maximise the one electron integrals, <i|h|i>.

        case("ONEPARTORBENMAX")
            tOnePartOrbEnMax=.true.
            MaxMinFac=-1
            tRotateVirtOnly=.true.
            IF(item.lt.nitems) THEN
                call Getf(OrbEnMaxAlpha)
            ELSE
                OrbEnMaxAlpha=1.0_dp
            ENDIF
! This sets the orbital rotation to find the coefficients which maximise the one particle orbital energies.
! i.e maximise the fock energies, epsilon_i = <i|h|i> + sum_j [<ij||ij>].

        case("ERLOCALIZATION")
            tERLocalization=.true.
            DiagMaxMinFac=-1
            IF(item.lt.nitems) THEN
                call Getf(DiagWeight)
            ELSE
                DiagWeight=1.0
            ENDIF
! This sets the orbital rotation to an Edmiston-Reudenberg localisation.  This maximises the self repulsion energy, i.e
! maximises the sum of the <ii|ii> terms.

        case("VIRTCOULOMBMAX")
            tVirtCoulombMax=.true.
            MaxMinFac=-1
            tRotateVirtOnly=.true.
! This sets the orbital rotation to maximise the sum_i_j <ij|ij> terms where both i and j are virtual spatial orbitals.

        case("VIRTEXCHANGEMIN")
            tVirtExchangeMin=.true.
            MaxMinFac=1
            tRotateVirtOnly=.true.
! This sets the orbital rotation to minimise the sum_i_j <ij|ji> terms where both i and j are virtual spatial orbitals.

        case("SHAKE")
! This will use the shake algorithm to iteratively enforce orthonormalisation
!  on the rotation coefficients calculated in the ROTATEORBS
! routine.  It finds a force matrix which moves the coefficients at a tangent
! to the constraint surface, from here, a normalisation
! will require just a small adjustment to ensure complete orthonormalisation,
! but not majorly affecting the new coefficients.
            tShake=.true.
            IF(item.lt.nitems) THEN
               call getf(ShakeConverged)
            ENDIF

        case("SHAKEAPPROX")
! This turns on the shake approximation algorithm.
! To be used if the matrix inversion required in the full shake algorithm cannot
! be performed.
! The approximation applies the iterative scheme to find lambda,
! to each constraint in succession, rather than simultaneously.
            tShake=.false.
            tShakeApprox=.true.

        case("SHAKEITER")
! Much like 'ROIteration', this overwrites the convergence criteria for the iterations of the shake constraints
! and instead performs only the number of iterations specified on this line.
            tShakeIter=.true.
            call Geti(ShakeIterMax)

        case("SHAKEDELAY")
! This option sets the shake orthonomalisation algorithm to only kick in after a certain number of rotatation iterations.  This
! potentially allows a large shift in the coefficients away from their starting point, followed by orthonormalisation to a
! significantly different position.
            tShakeDelay=.true.
            call Geti(ShakeStart)
        case("LAGRANGE")
! This will use a non-iterative lagrange multiplier for each component of each rotated vector
! in the rotateorbs routines in order to
! attempt to maintain orthogonality. This currently does not seem to work too well!
            tShake=.false.
            tLagrange=.true.

        case("SEPARATEOCCVIRT")
! This option applies to the orbital rotation.
! If present, the virtual and occuppied orbitals are localised separately.
! This has the  advantage of keeping the HF reference determinant the same.
            tSeparateOccVirt=.true.

        case("ROTATEOCCONLY")
! This option applies to orbital rotation.  It separates the orbitals into occupied and virtual,
! and rotates the occupied while keeping the virtual as the HF.
            tSeparateOccVirt=.true.
            tRotateOccOnly=.true.

        case("ROTATEVIRTONLY")
! This option rotates the virtual orbitals while keeping the occupied as the HF.
            tSeparateOccVirt=.true.
            tRotateVirtOnly=.true.

        case("ROITERATION")
! Specifying this keyword overwrites the convergence limit from the ROTATEORBS line,
! and instead runs the orbital rotation for as many
! iterations as chosen on this line.
            tROIteration=.true.
            call Geti(ROIterMax)

        case("MAXHLGAP")
            tMaxHLGap=.true.

        case("ROTATEDORBS")
! This is a flag which is to be included if the FCIDUMP being read in is one containing rotated orbitals.
! The only difference is that the
! H elements for single excitations are no longer 0 (as for HF),
! and walkers on singly excited determinants must be included in the energy
! calculations.
            tRotatedOrbs=.true.
            tNoBrillouin=.true.
            tBrillouinsDefault=.false.

        case("SPINORBS")
! This flag simply uses spin orbitals to perform the rotation rather than spatial orbitals.
! IF UHF=.true. is present in the FCIDUMP file, this will happen automatically, otherwise this keyword is required.
            tSpinOrbs=.true.

        case("READINTRANSMAT")
! This sets the rotation routine to read in a transformation matrix (coefft1) given by the file TRANSFORMMAT,
! and use the transformation
! routine and print rofcidump routines to transform the orbitals and print out a new dump file.
            tReadInCoeff=.true.

        case("USEMP2VDM")
! Once in the rotation routine, the MP2 variational density matrix is calculated,
! and this is used to transform the orbitals and print and
! new FCIDUMP file.
            tUseMP2VarDenMat=.true.
            tSeparateOccVirt=.true.
            tRotateVirtOnly=.true.
            tShake=.false.

        case("USECINATORBS")
! This rotation option is slightly different, it first requires a spawning calculation
! from which the amplitudes of the wavefunction are
! obtained.  From here, the one electron reduced density matrix is calculated, and the eigenvectors
!  of this are used to rotate the HF orbitals.
! A new ROFCIDUMP file is then produced in the new natural orbitals.
            tFindCINatOrbs=.true.
            tShake=.false.

        case("USEHFORBS")
            tUseHFOrbs=.true.
            tShake=.false.
            tSeparateOccVirt=.true.

        case("RANLUXLEV")
!This is the level of quality for the random number generator. Values go from 1 -> 4. 3 is default.
            call readi(iRanLuxLev)

        case("CALCEXACTSIZESPACE")
!This option will calculate the exact size of the symmetry allowed space of determinants. Will scale badly.
            tExactSizeSpace=.true.

        case("CALCMCSIZETRUNCSPACE")
!This option will approximate the exact size of the symmetry allowed truncated space of determinants by MC.
!The variance on the value will decrease as 1/N_steps
            tMCSizeTruncSpace=.true.
            CALL Geti(iMCCalcTruncLev)
            CALL GetiLong(CalcDetCycles)
            CALL GetiLong(CalcDetPrint)

        case("CALCMCSIZESPACE")
!This option will approximate the exact size of the symmetry allowed space of determinants by MC.
! The variance on the value will decrease as 1/N_steps
            tMCSizeSpace=.true.
            CALL GetiLong(CalcDetCycles)
            CALL GetiLong(CalcDetPrint)

        case("NONUNIFORMRANDEXCITS")
!This indicates that the new, non-uniform O[N] random excitation generators are to be used.
!CYCLETHRUORBS can be useful if we have small basis sets or v high restrictive symmetry and will eliminate
!large numbers of unsuccessful random draws of orbitals by calculating how many allowed orbitals there are
!and cycling through them until the allowed one is drawn, rather than randomly drawing and redrawing until
!an allowed orbital is found. For large basis sets, the chance of drawing a forbidden orbital is small
!enough that this should be an unneccesary expense.
            tNonUniRandExcits=.true.
            do while (item.lt.nitems)
                call readu(w)
                select case(w)
                    case("CYCLETHRUORBS")
                        tCycleOrbs=.true.
                    case("NOSYMGEN")
!Do not generate symmetry-allowed excitations only, but all excitations. Spin-symmetry is still taken into account.
                        tNoSymGenRandExcits=.true.
                    case("IMPORTANCESAMPLE")
!Importance sample the excitations for FCIMCPar
                        CALL Stop_All("ReadSysInp","IMPORTANCESAMPLE option depreciated")
!                        tImportanceSample=.true.
                    case("PICK-VIRT-UNIFORM")
                        ! Pick virtual orbitals randomly and uniformly in the
                        ! 3rd generation of random excitation generators
                        ! (symrandexcit3.F90)
                        tPickVirtUniform = .true.
                        tBuildOccVirtList = .true.
                    case("PICK-VIRT-UNIFORM-MAG")
                        ! Pick virtual orbitals randomly and uniformly in the
                        ! 3rd generation of random excitation generators
                        ! (symrandexcit3.F90)
                        tPickVirtUniform = .true.
                        tBuildOccVirtList = .true.
                        tBuildSpinSepLists = .true.
                    case("HEL-WEIGHTED-SLOW")
                        ! Pick excitations from any site with a generation
                        ! probability proportional to the connectiong HElement
                        ! --> Lots of enumeration
                        ! --> Very slow
                        ! --> Maximal possible values of tau.
                        tGenHelWeighted = .true.
                    case("4IND-WEIGHTED")
                        ! Weight excitations based on the magnitude of the
                        ! 4-index integrals (ij|ij)
                        tGen_4ind_weighted = .true.
                    case("4IND-REVERSE")
                        ! Weight excitations based on the magnitude of the
                        ! actual hamiltonian matrix elements (at least for
                        ! doubles). This is effectively the "reverse" of
                        ! 4IND-WEIGHTED as above.
                        tGen_4ind_reverse = .true.
                    case("4IND-WEIGHTED-PART-EXACT")
                        ! Weight excitations as in 4IND-WEIGHTED, except for
                        ! double excitations with the same spin which are
                        ! weighted according to:
                        ! sqrt(((ii|aa) + (jj|aa))(<ij|ab>-<ij|ba>))
                        tGen_4ind_weighted = .true.
                        tGen_4ind_part_exact = .true.
                    case("4IND-WEIGHTED-LIN-EXACT")
                        ! Weight excitations as in 4IND-WEIGHTED, except for
                        ! double excitations with the same spin which are
                        ! weighted according to:
                        ! (1/M)(<ij|ab> - <ij|ba>)
                        ! (The second half of this only affecting the choice
                        ! of electron b)
                        tGen_4ind_weighted = .true.
                        tGen_4ind_lin_exact = .true.
                    case("4IND-WEIGHTED-2")
                        ! Second attempt at 4ind weighted generator
                        !
                        ! This version disables symmetric selection of the
                        ! orbitals. I.e. for non-parallel electron choice, the
                        ! orbitals are chosen one spin first, then the other
                        !
                        ! --> Very slightly better Hij/pgen ratio
                        ! --> Much faster runtime
                        tGen_4ind_2 = .true.
                        tGen_4ind_part_exact = .true.
                        tGen_4ind_2_symmetric = .false.

                    case("4IND-WEIGHTED-UNBOUND")
                        ! [Werner Dobrautz 26.4.2017:]
                        ! new tests for optimizations for the latest
                        ! excitation generator implementation, without using
                        ! the artificial lower bounds for the energy
                        tGen_4ind_2 = .true.
                        tGen_4ind_part_exact = .true.
                        tGen_4ind_2_symmetric = .false.
                        tGen_4ind_unbound = .true.

                        ! make a few small tests for the frequency histograms
                        if (item < nitems) then
                            call readu(w)

                            select case (w)
                            case ("IIAA")
                                ! weight with <ii|aa> instead of <ia|ia>
                                t_iiaa = .true.

                            case ("RATIO")
                                ! weigh with the ratio <ia|ia>/<ja|ja>
                                t_ratio = .true.

                            case ("IIAA-RATIO", "RATIO-IIAA")
                                t_iiaa = .true.
                                t_ratio = .true.

                            end select
                        end if


                    case("4IND-WEIGHTED-2-SYMMETRIC")
                        ! The other version of this generator. This permits
                        ! selecting orbitals in both directions
                        tGen_4ind_2 = .true.
                        tGen_4ind_part_exact = .true.
                        tGen_4ind_2_symmetric = .true.

                    case("PCPP")
                       ! the precomputed power-pitzer excitation generator
                       t_pcpp_excitgen = .true.

                    case("PCHB")
                       ! the precomputed heat-bath excitation generator (uniform singles)
                       t_pchb_excitgen = .true.
                    case("UEG")
                        ! Use the new UEG excitation generator.
                        ! TODO: This probably isn't the best way to do this
                        tUEGNewGenerator = .true.
                        tLatticeGens = .true.

                    case default
                        call Stop_All("ReadSysInp",trim(w)//" not a valid keyword")
                end select
            enddo

        case("PCHB-WEIGHTED-SINGLES")
            ! Enable using weighted single excitations with the pchb excitation generator
            t_pchb_weighted_singles = .true.

        case("SPAWNLISTDETS")
!This option will mean that a file called SpawnOnlyDets will be read in,
! and only these determinants will be allowed to be spawned at.
            CALL Stop_All("ReadSysInp","SPAWNLISTDETS option depreciated")
!            tListDets=.true.

        case("UMATEPSILON")

            ! For systems that read in from an FCIDUMP file, any two-electron
            ! integrals are screened against a threshold parameter. Below
            ! this they are ignored.
            !
            ! By default, this parameter is 10-e8, but it can be changed here.
            call readf(UMatEps)

         case("LMATEPSILON")
            ! Six-index integrals are screened, too, with the default being 1e-10
            call readf(LMatEps)

        case("NOSINGEXCITS")
!This will mean that no single excitations are ever attempted to be generated.
            tNoSingExcits=.true.
        case("DIAGONALTMAT")
!Implies that the orbital basis are eigenfunctions of the KE operator, so TMAT can be stored as a diagonal matrix to save space.
          tOneElecDiag=.true.   !One electron integrals diagonal
        case("ROHF")
!This is an option for open-shell systems to specify that the integrals are *restricted* open-shell integrals.
!This will save memory (around a factor of 16) for the integral storage,
!  but the FCIDUMP file should be the same as before (ie in UHF form).
            tROHF=.true.
            tNoBrillouin=.true.
            tBrillouinsDefault=.false.
            IF(tFindCINatOrbs) CALL Stop_All("ReadSysInp","For orbital rotations of open shell "&
            &//"systems, UMAT must be stored in spin &
            & orbitals - cannot be compressed using ROHF.")

        case("LZTOT")
            tFixLz=.true.
            call readi(LzTot)
        case("KPOINTS")
            tKPntSym=.true.
        case("MOLPROMIMIC")
            !Mimic the run-time behaviour of molpros NECI implementation
            tMolpro=.true.
            tMolproMimic=.true.
        case("READ_ROFCIDUMP")
            ! Overwrite current FCIDUMP name, and instead look for a file
            ! called "ROFCIDUMP".
            FCIDUMP_name = 'ROFCIDUMP'
        case("COMPLEXORBS_REALINTS")
            !We have complex orbitals, but real integrals. This means that we only have 4x permutational symmetry,
            !so we need to check the (momentum) symmetry before we look up any integrals
            tComplexOrbs_RealInts = .true.

         case("COMPLEXWALKERS-REALINTS")
            ! In case complex walkers shall be used but not complex basis functions,
            ! such that the integrals are real and have full symmetry
            tComplexWalkers_RealInts = .true.

        case("SYSTEM-REPLICAS")
            ! How many copies of the simulation do we want to run in parallel?
            ! This can only be done using mneci.x, where the size of the
            ! representation (i.e. lenof_sign) is permitted to vary at runtime
#ifdef PROG_NUMRUNS_
            call readi(inum_runs)
            tMultiReplicas = .true.
#ifdef CMPLX_
            lenof_sign = 2*inum_runs
#else
            lenof_sign = inum_runs
#endif
            if (inum_runs > inum_runs_max) then
                write(6,*) 'Maximum SYSTEM-REPLICAS: ', inum_runs_max
                call stop_all(t_r, 'SYSTEM-REPLICAS is greater than maximum &
                                   &permitted value')
            end if
#else
            call readi(itmp)
#ifdef DOUBLERUN_
            if (itmp /= 2) then
#else
            if (itmp /= 1) then
#endif
                call stop_all(t_r, "mneci.x build must be used for running &
                                   &with multiple simultaneous replicas")
            endif
#endif

        case("HEISENBERG")
            tHeisenberg = .true.

        case("GIOVANNIS-BROKEN-INIT")
            ! Giovanni's scheme for initialising determinants with the correct
            ! spin an symmetry properties in a wider range of cases than
            ! currently supported.

            ! Looks nice, but it currently breaks lots of other stuff!
            tGiovannisBrokenInit = .true.

         case("SPIN-CONSERVING-GAS")
            tGAS = .true.
            tGASSpinRecoupling = .true.

         case("PART-CONSERVING-GAS")
            tGAS = .true.
            tGASSpinRecoupling = .false.

        case("ENDSYS")
            exit system
        case default
            call report("Keyword "                                    &
   &          //trim(w)//" not recognized in SYSTEM block",.true.)
        end select
      end do system

      if(NEL.eq.0)                                                    &
   &     call report("Number of electrons cannot be zero.",.true.)

      if (.not. tUEG2) then
          if(THUB.OR.TUEG.OR..NOT.(TREADINT.OR.TCPMD.or.tVASP)) then
            if(NMAXX.EQ.0)                                               &
            &        call report("Must specify CELL "                          &
            &        //"- the number of basis functions in each dim.",         &
            &        .true.)
            if(.NOT.THUB.AND. near_zero(BOX))                                &
            &        call report("Must specify BOX size.",.true.)
            if(TTILT.AND..NOT.THUB)                                      &
            &        call report("TILT can only be specified with HUBBARD.",.true.)
          endif
      end if

    END SUBROUTINE SysReadInput



    Subroutine SysInit
      Use global_utilities
      use SymData, only: tAbelian,TwoCycleSymGens,nSymLabels
      use constants, only: Pi, Pi2, THIRD
      use legacy_data, only: CSF_NBSTART
      use read_fci
      use sym_mod
      use SymExcitDataMod, only: kPointToBasisFn
      implicit none
      character(*), parameter :: this_routine='SysInit'
      integer ierr

      CHARACTER(len=1) CPAR(3)
      CHARACTER(len=3) CPARITY
! For init of mom
      TYPE(BasisFN) G

!  For the UEG
      real(dp) FKF,Rs

! General variables
      INTEGER i,j,k,l,iG,k2_max_2d
      INTEGER len
      type(timer), save :: proc_timer
      real(dp) SUM
! Called functions
      TYPE(BasisFN) FrzSym
      logical kallowed
      real(dp) dUnscaledE
      real(dp), allocatable :: arr_tmp(:,:)
      integer, allocatable :: brr_tmp(:)
!     UEG2
      integer :: AllocateStatus
      real(dp), parameter :: EulersConst = 0.5772156649015328606065120900824024_dp


!      write (6,*)
!      call TimeTag()
!      if (.not.TCPMD) call Envir()
!      write (6,*)

      ECORE=0.0_dp

! //AJWT TBR
!      IFDET=0
!      TRHOIJND=.false.
      proc_timer%timer_name='SysInit   '
      call set_timer(proc_timer)

      ! if a total spin was set, set the spin projection if unspecified
      if(tInitializeCSF .and. .not. TSPN) then
         ! S is physical spin, LMS is an integer (-2*Ms)
         LMS = int(2 * S2Init)
         TSPN = .true.
         write(iout,*) 'Spin projection unspecified, assuming Ms=S'
      end if

!C ==-------------------------------------------------------------------==
!C..Input parameters
      WRITE(6,'(A)') '======== SYSTEM =========='
      WRITE(6,'(A,I5)') '  NUMBER OF ELECTRONS : ' , NEL
      IF(TSPN) THEN
          WRITE(6,*) ' Restricting the spin state of the system, TSPN : ' , TSPN
      ELSE
          WRITE(6,*) ' No restriction on the spin state of the system, TSPN : ' , TSPN
      ENDIF
      NBASISMAX(1:5,1:7)=0
      TSPINPOLAR=.FALSE.
      DO I=1,3
         SymRestrict%k(I)=IPARITY(I)
      ENDDO
      SymRestrict%Ms=IPARITY(4)
      SymRestrict%Sym%s=IPARITY(5)


      IF(TSPN) THEN
!C.. If we're doing monte carlo without generating a list of
!C.. determinants, we cannot as yet force spin or parity, except by
!C.. restricting the basis set.  This will only work for Ms=NEL/2
!C         IF((.NOT.TBLOCK).AND.(.NOT.TCALCHMAT.OR.NTAY.LE.0)) THEN
!C            WRITE(6,*) 'TSPN set to TRUE.  Determinant list not being',
!C     &         ' used for MC.  Forcing MS=Nel/2'
!C            IF(MOD(NEL,2).EQ.0) THEN
!C               LMS=NEL/2
!C            ELSE
!C               LMS=NEL
!C            ENDIF
!C            TSPINPOLAR=.TRUE.
!C         ENDIF
          IF(MOD(LMS+NEL*2,2).NE.MOD(NEL,2)) THEN
            WRITE(6,*) 'LMS=',LMS,' not achievable with',NEL,' electrons'
            WRITE(6,*) 'Resetting LMS'
            LMS=MOD(NEL,2)
          ENDIF
          LMS2=LMS
      ELSE
          if (tUEG2) then
              TSPN =.TRUE.
              if (dimen==1) then
!                 1D calculations are always spin polarised
                  LMS = NEL
                  !TSPINPOLAR = .TRUE.
              else
                  LMS = 0
              end if
               LMS2=LMS
          end if
      ENDIF
      WRITE(6,*) ' GLOBAL MS : ' , LMS

      IF((NBASISMAX(2,3).eq.1).or.tROHF) THEN
!If we are dealing with an open shell system, the calculation of symreps will sometimes fail.
!This will have consequences for the rest of the program, so in a slightly hacky way, we can simply
!not reorder the orbitals by energy, so that they remain in symmetries.
!The reason it fails it that it looks for a complete set of orbitals which are degenerate
!and ignores those in the symmetry classification. Unfortunately in ROHF and UHF there aren't such things.
          WRITE(6,'(A)') "  Open shell system - SYMIGNOREENERGIES set.  "
!          tHFNOORDER=.true.
          tSymIgnoreEnergies=.true.
      ENDIF

      if(tUseBrillouin) THEN
         WRITE(6,'(A)') "  Using Brillouin's Theorem to ignore single excitations"
      endif
      if(tStoreAsExcitations) THEN
         write(6,'(A)') "  Storing determinants as excitations from the HF determinant.  WARNING this may not work!"
         IF(nEL.lt.8) call stop_all(this_routine, '  tStoreAsExcitations requires nEl>=8.')
      endif

      ! Conditions for using old CSF routines
      if (tCSFOld) then
          write (6, '(A)') "************ Using Old CSF routines ************"
          if (tSPN) then
              write(6, '(A,I3)') "  Restricting total spin*2 to ", STOT
              if (LMS > STOT) &
                  call stop_all (this_routine, "Cannot have LMS>STOT")

              ! Encode the symmetry for the total spin in LMS
              LMS2 = LMS + STOT*CSF_NBSTART
          endif
          nBasisMax(4,7) = 1
      endif

      ! Conditions for using CSFs
      if (tCSF) then
          write (6, '(A)') "************ Using CSFs for calculation **********"

          if (LMS > STOT) call stop_all (this_routine, "Cannot have LMS>STOT")

          if (tHPHF) then
              call stop_all (this_routine, "CSFs not compatible with HPHF")
          endif

      endif

      if (tTruncateCSF .and. (.not. tCSF)) then
          call stop_all (this_routine, "CSFs required to use truncate-csf")
      endif

      TwoCycleSymGens=.false.
      IF(TCPMD) THEN
         WRITE(6,'(A)') '  *** GENERIC SYSTEM USING KOHN-SHAM ORBITALS ***  '
         CALL CPMDSYSTEMINIT(LEN)
         IF(TPARITY) THEN
            WRITE(6,"(A)",advance='no') '  SYMMETRIES : '
            CALL WRITEALLSYM(5,SymRestrict,.true.)
         ENDIF
         IF(THFORDER) WRITE(6,'(A)')      "  Ordering according to 1-electron energies.  "
      ELSEIF(tVASP) THEN
         WRITE(6,'(A)') '  *** GENERIC SYSTEM USING HF ORBITALS PRODUCED BY VASP ***  '
         CALL VASPSystemInit(LEN)
      ELSEIF(TREADINT) THEN
!C.. we read in the integrals from FCIDUMP and ignore most config
!C..
         WRITE(6,'(A)') '  *** GENERIC SYSTEM ***  '
         IF(THUB) THEN
            THUB=.FALSE.
            WRITE(6,'(A)') "  Setting THUB=.FALSE.  "
         ENDIF
         TwoCycleSymGens=.true.
         IF(TDFREAD) THEN
            WRITE(6,'(A)') "  Reading Density fitted integrals.  "
            LMSBASIS=LMS
            CALL InitDFBasis(nBasisMax,Len)
         ELSEIF(tRIIntegrals.or.tCacheFCIDUMPInts) THEN
!tCacheFCIDUMPInts means that we read in all the integrals from the FCIDUMP integral file, but store them contiguously in the cache
            LMSBASIS=LMS
            IF(tRIIntegrals) THEN
                WRITE(6,'(A)') "  Reading RI integrals.  "
                CALL InitRIBasis(nBasisMax,Len)
                LMSBASIS=LMS
            ELSE
                WRITE(6,'(A)') "  Reading in all integrals from FCIDUMP file, but storing them in a cache...  "
                tAbelian=.true.
            ENDIF
            CALL INITFROMFCID(NEL,NBASISMAX,LEN,LMSBASIS,TBIN)
            nBasisMax(3,3)=1    !No momentum conservation
            nBasisMax(4,1)=-1   !God knows...
            nBasisMax(4,2)=1    !Ditto
!.. Correspond to ISS=0
!            NBASISMAX(2,3)=-1  !indicate we generate on the fly
            iSpinSkip=-1
!NBASISMAX(2,3)  !indicate we generate on the fly
!C.. say we're a UHF det so all singles are 0
!            IF(LMS.EQ.0) THEN
!               NBASISMAX(4,5)=1
!            ELSE
!               NBASISMAX(4,5)=2
!            ENDIF
            IF(NBASISMAX(2,3).EQ.1) then
                WRITE(6,'(A)') " Unrestricted calculation.  Cave Arthropodia.  "
            ELSEIF(tROHF) then
                WRITE(6,'(A)') "  High-spin restricted calculation. Seriously Cave Arthropodia.  "
            ENDIF
         ELSE
            LMSBASIS=LMS
!            WRITE(6,*) "TBIN:",tBin
            CALL INITFROMFCID(NEL,NBASISMAX,LEN,LMSBASIS,TBIN)
!C.. say we're a UHF det so all singles are 0
!            IF(LMS.EQ.0) THEN
!               NBASISMAX(4,5)=1
!            ELSE
!               NBASISMAX(4,5)=2
!            ENDIF
            IF(tStoreSpinOrbs) then
                WRITE(6,'(A)') "  Unrestricted calculation.  Cave Arthropodia.  "
            ELSEIF(tROHF) then
                WRITE(6,'(A)') "  High-spin restricted calculation. Seriously Cave Arthropodia.  "
            ENDIF
         ENDIF
      ELSE

    ! ======================================================
      if (tUEG2) then
          !determines the recip lattice, the real volume, R_s and the Fermi vector
          call LatticeInit(RS, FKF)

          ! engergy cutoff not given -> set to 2* Fermi vector
          if(.not. torbEcutoff) then
              orbEcutoff =(2.0_dp*FKF/k_lattice_constant)**2.0_dp
              torbEcutoff = .true.
          end if
          ! if cell size not given
          if(NMAXX == 0 .and.  NMAXY == 0 .and. NMAXZ == 0) then
              call CalcCell
          end if

          TTILT=.FALSE.
          ALAT=0.0_dp   !shouldn't be used in the UEG2 part...
          NMAX=MAX(NMAXX,NMAXY,NMAXZ)
          NNR=NMSH*NMSH*NMSH
          WRITE(6,'(A)') '  *** In UEG2 ***  '
          WRITE(6,'(A)') '  *** UNIFORM ELECTRON GAS CALCULATION ***  '
          WRITE(6,'(A,F20.16)') '  Electron Gas Rs set to ',FUEGRS
          IF(TPARITY) THEN
              WRITE(6,*) ' MOMENTUM : ',(IPARITY(I),I=1,3)
          ENDIF

           WRITE(6,'(A,I5)') '  Dimension : ' , Dimen
           WRITE(6,*) '  Reciprocal lattice constant : ' ,  k_lattice_constant
           WRITE(6,'(A,I5)') '  NMAXX : ' , NMAXX
           WRITE(6,'(A,I5)') '  NMAXY : ' , NMAXY
           WRITE(6,'(A,I5)') '  NMAXZ : ' , NMAXZ
           WRITE(6,'(A,I5)') '  NMSH : ' , NMSH
           WRITE(6,*) " Wigner-Seitz radius Rs=",RS
           WRITE(6,*) " Fermi vector kF^2=",FKF**2
           WRITE(6,*) " Fermi Energy EF=",FKF*FKF/2
           write(6,*) " Unscaled fermi vector kF=", FKF/k_lattice_constant
           WRITE(6,*) " Unscaled Fermi Energy nmax**2=",(FKF*FKF/2)/(0.5*(2*PI/ALAT(5))**2)
           IF(.not. (OrbECutoff .isclose. 1e-20_dp)) WRITE(6,*) " Orbital Energy Cutoff:",OrbECutoff
           WRITE(6,'(1X,A,F19.5)') '  VOLUME : ' , OMEGA
           WRITE(6,*) ' TALPHA : ' , TALPHA
           WRITE(6,'(1X,A,F19.5)') '  ALPHA : ' , ALPHA

          ALPHA=(OMEGA)**THIRD*ALPHA
           WRITE(6,'(1X,A,F19.5)') '  SCALED ALPHA : ' , ALPHA
           WRITE(6,*) 'Madelung constant: ',  Madelung

!C..
!C..Calculate number of basis functions
!C.. UEG allows from -NMAX->NMAX
          IF(TSPINPOLAR) THEN
             NBASISMAX(4,1)=1
!C.. spinskip
             NBASISMAX(2,3)=1
          ELSE
             NBASISMAX(4,1)=-1
!C.. spinskip
!  If spinskip is unset
             IF(NBASISMAX(2,3).EQ.0) NBASISMAX(2,3)=2
          ENDIF

          NBASISMAX(4,2)=1
          NBASISMAX(1,1)=-NMAXX
          NBASISMAX(1,2)=NMAXX
          NBASISMAX(2,1)=-NMAXY
          NBASISMAX(2,2)=NMAXY
          NBASISMAX(3,1)=-NMAXZ
          NBASISMAX(3,2)=NMAXZ
          NBASISMAX(1,3)=-1
          LEN=(2*NMAXX+1)*(2*NMAXY+1)*(2*NMAXZ+1)*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
!C.. UEG
          NBASISMAX(3,3)=-1

!C.. we actually store twice as much in arr as we need.
!C.. the ARR(1:LEN,1) are the energies of the orbitals ordered according to
!C.. BRR.  ARR(1:LEN,2) are the energies of the orbitals with default
!C.. ordering.
!C.. ARR is reallocated in IntFreezeBasis if orbitals are frozen so that it
!C.. has the correct size and shape to contain the eigenvalues of the active
!C.. basis.
          WRITE(6,'(A,I5)') "  NUMBER OF SPIN ORBITALS IN BASIS : ", Len
          Allocate(Arr(LEN,2),STAT=ierr)
          LogAlloc(ierr,'Arr',2*LEN,8,tagArr)
          ! // TBR
          !      IP_ARRSTORE=IP_ARR
          ARR=0.0_dp
          Allocate(Brr(LEN),STAT=ierr)
          LogAlloc(ierr,'Brr',LEN,4,tagBrr)
          BRR(1:LEN)=0
          Allocate(G1(Len),STAT=ierr)
          LogAlloc(ierr,'G1',LEN,BasisFNSizeB,tagG1)
          G1(1:LEN)=NullBasisFn


          IF(TCPMD) THEN
              WRITE(6,'(A)') '*** INITIALIZING BASIS FNs FROM CPMD ***'
              CALL CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
              NBASIS=LEN
              iSpinSkip=NBasisMax(2,3)
          ELSEIF(tVASP) THEN
              WRITE(6,'(A)') '*** INITIALIZING BASIS FNs FROM VASP ***'
              CALL VASPBasisInit(ARR,BRR,G1,LEN) ! This also modifies nBasisMax
              NBASIS=LEN
              iSpinSkip=NBasisMax(2,3)
          ELSEIF(TREADINT.AND.TDFREAD) THEN
              WRITE(6,'(A)') '*** Creating Basis Fns from Dalton output ***'
              call InitDaltonBasis(Arr,Brr,G1,Len)
              nBasis=Len
              call GenMolpSymTable(1,G1,nBasis)
          ELSEIF(TREADINT) THEN
          !This is also called for tRiIntegrals and tCacheFCIDUMPInts
              WRITE(6,'(A)') '*** CREATING BASIS FNs FROM FCIDUMP ***'
              CALL GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN,TBIN)
              NBASIS=LEN

  !C.. we're reading in integrals and have a molpro symmetry table
              IF(lNoSymmetry) THEN
                  WRITE(6,*) "Turning Symmetry off"
                  DO I=1,nBasis
                      G1(I)%Sym%s=0
                  ENDDO
                  CALL GENMOLPSYMTABLE(1,G1,NBASIS)
                  DO I=1,nBasis
                      G1(I)%Sym%s=0
                  ENDDO
              ELSE
                  CALL GENMOLPSYMTABLE(NBASISMAX(5,2)+1,G1,NBASIS)
              ENDIF

          ELSE
!C.. Create plane wave basis functions

              WRITE(6,*) "Creating plane wave basis."

              IG=0
              DO I=NBASISMAX(1,1),NBASISMAX(1,2)
                  DO J=NBASISMAX(2,1),NBASISMAX(2,2)
                      DO K=NBASISMAX(3,1),NBASISMAX(3,2)
                          DO L=NBASISMAX(4,1),NBASISMAX(4,2),2
                              !if (dimen ==1 .AND. L ==1) cycle
                              G%k(1)=I
                              G%k(2)=J
                              G%k(3)=K
                              G%Ms=L

                              IF(KALLOWED(G,NBASISMAX)) THEN
                                  CALL GetUEGKE(I,J,K,ALAT,tUEGTrueEnergies,tUEGOffset,k_offset,SUM,dUnscaledE)
                                  IF(dUnscaledE.gt.OrbECutoff) CYCLE
                                  IG=IG+1
                                  ARR(IG,1)=SUM
                                  ARR(IG,2)=SUM
                                  BRR(IG)=IG
                                  !C..These are the quantum numbers: n,l,m and sigma
                                  G1(IG)%K(1)=I
                                  G1(IG)%K(2)=J
                                  G1(IG)%K(3)=K
                                  G1(IG)%MS=L
                                  G1(IG)%Sym=TotSymRep()
                              ENDIF

                          ENDDO
                      ENDDO
                  ENDDO
              ENDDO
              if(real_lattice_type == 'sc' .AND. maxval(G1%K(1)) .ge. NMAXX)  WRITE(6,*) 'ERROR IN CELL CALCULATION! '

!C..Check to see if all's well
              WRITE(6,*) ' NUMBER OF BASIS FUNCTIONS : ' , IG
              NBASIS=IG

              ! calculate k-vectors in cartesian coordinates
              allocate(kvec(NBASIS, 3), STAT = AllocateStatus)
              IG=0
              DO I=1, NBASIS
                  IG=IG+1
                  kvec(IG, 1)=  k_lattice_vectors(1,1)*G1(IG)%K(1) &
                                      +k_lattice_vectors(2,1)*G1(IG)%K(2)&
                                      +k_lattice_vectors(3,1)*G1(IG)%K(3)
                  kvec(IG, 2)=  k_lattice_vectors(1,2)*G1(IG)%K(1)&
                                      +k_lattice_vectors(2,2)*G1(IG)%K(2)&
                                      +k_lattice_vectors(3,2)*G1(IG)%K(3)
                  kvec(IG, 3)=  k_lattice_vectors(1,3)*G1(IG)%K(1)&
                                      +k_lattice_vectors(2,3)*G1(IG)%K(2)&
                                      +k_lattice_vectors(3,3)*G1(IG)%K(3)
              ENDDO

              IF(LEN.NE.IG) THEN
                  if(OrbECutoff.gt.-1e20_dp) then
                      write(6,*) " Have removed ", LEN-IG, " high energy orbitals "
                      ! Resize arr and brr.
                      allocate(arr_tmp(nbasis,2),brr_tmp(nbasis),stat=ierr)
                      arr_tmp = arr(1:nbasis,:)
                      brr_tmp = brr(1:nbasis)
                      deallocate(arr,brr,stat=ierr)
                      LogDealloc(tagarr)
                      LogDealloc(tagbrr)
                      allocate(arr(nbasis,2),brr(nbasis),stat=ierr)
                      LogAlloc(ierr,'Arr',2*nbasis,8,tagArr)
                      LogAlloc(ierr,'Brr',nbasis,4,tagBrr)
                      arr = arr_tmp
                      brr = brr_tmp
                      deallocate(arr_tmp, brr_tmp, stat=ierr)
                  else
                      WRITE(6,*) " LEN=",LEN,"IG=",IG
                      call stop_all(this_routine, ' LEN NE IG ')
                  endif
              ENDIF

              if (.not. tHub) CALL GENMOLPSYMTABLE(1,G1,NBASIS)
          ENDIF

          IF(tFixLz) THEN
              WRITE(6,'(A)') "****** USING Lz SYMMETRY *******"
              WRITE(6,'(A,I5)') "Pure spherical harmonics with complex orbitals used to constrain Lz to: ",LzTot
              WRITE(6,*) "Due to the breaking of the Ml degeneracy, the fock energies are slightly wrong, "&
              &//"on order of 1.0e-4_dp - do not use for MP2!"
              if(nsymlabels.gt.4) then
                  call stop_all(this_routine,"D2h point group detected. Incompatable with Lz symmetry conserving "&
                  &//"orbitals. Have you transformed these orbitals into spherical harmonics correctly?!")
              endif
          ENDIF

!C..        (.NOT.TREADINT)
!C.. Set the initial symmetry to be totally symmetric
          FrzSym=NullBasisFn
          FrzSym%Sym=TotSymRep()
          CALL SetupFreezeSym(FrzSym)
!C..Now we sort them using SORT2 and then SORT

!C.. This sorts ARR and BRR into order of ARR [AJWT]
          IF(.NOT.THFNOORDER) THEN
              CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
          ELSE
              !.. copy the default ordered energies.
              CALL DCOPY(NBASIS,ARR(1,1),1,ARR(1,2),1)
          ENDIF
!      WRITE(6,*) THFNOORDER, " THFNOORDER"

          if(.not.tMolpro) then
          !If we are calling from molpro, we write the basis later (after reordering)
              CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)
          endif

          IF(NEL.GT.NBASIS) call stop_all(this_routine, 'MORE ELECTRONS THAN BASIS FUNCTIONS')
          CALL neci_flush(6)

          NOCC=NEl/2
          IF(TREADINT) THEN
      !C.. we're reading in integrals and have a molpro symmetry table
              IF(lNoSymmetry) THEN
                  WRITE(6,*) "Turning Symmetry off"
                  CALL GENMOLPSYMREPS()
              ELSE
                  CALL GENMOLPSYMREPS()
              ENDIF
          ELSEIF(TCPMD) THEN
      !C.. If TCPMD, then we've generated the symmetry table earlier,
      !C.. but we still need the sym reps table.
              call stop_all(this_routine,'CPMD interface depricated')
!              CALL GENCPMDSYMREPS(G1,NBASIS,ARR,1.e-5_dp)
          ELSEIF(tVASP) THEN
      !C.. If VASP-based calculation, then we've generated the symmetry table earlier,
      !C.. but we still need the sym reps table. DEGENTOL=1.0e-6_dp. CHECK w/AJWT.
              CALL GENSYMREPS(G1,NBASIS,ARR,1.e-6_dp)
          ELSEIF(THUB.AND..NOT.TREAL) THEN
              CALL GenHubMomIrrepsSymTable(G1,nBasis,nBasisMax)
              CALL GENHUBSYMREPS(NBASIS,ARR,BRR)
              CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)
          ELSE
      !C.. no symmetry, so a simple sym table
              CALL GENMOLPSYMREPS()
          ENDIF

      !// TBR
      !      WRITE(6,*) ' ETRIAL : ',ETRIAL
      !      IF(FCOUL.NE.1.0_dp)  WRITE(6,*) "WARNING: FCOUL is not 1.0_dp. FCOUL=",FCOUL
          IF(FCOULDAMPBETA.GT.0) WRITE(6,*) "FCOUL Damping.  Beta ",FCOULDAMPBETA," Mu ",FCOULDAMPMU
          call halt_timer(proc_timer)

          !calculate tau if not given
          if (TAU .lt. 0.0_dp) then
              call  CalcTau
          end if

          if(tMadelung .AND. near_zero(Madelung) .AND. dimen == 3) then
             Madelung=calc_madelung()
          else if (tMadelung .AND. near_zero(Madelung) .AND. dimen /= 3) then
              call stop_all (this_routine, "Calculation of Madelung constant works in 3D only!")
          end if

      return
      endif  !UEG2
! ======================================================

          IF(TUEG) THEN
             WRITE(6,'(A)') '  *** UNIFORM ELECTRON GAS CALCULATION ***  '
             if (iPeriodicDampingType /= 0) then
                 ! We are using either a screened or an attenuated Coulomb
                 ! potential for calculating the exchange integrals.
                 ! This means that we need to be able to distinguish between
                 ! exchange integrals and normal Coulomb integrals and hence we
                 ! should refer to spin-orbitals throughout.
                 nBasisMax(2,3) = 1
                 tStoreSpinOrbs = .true.
             end if

             if(t_ueg_transcorr)WRITE(6,*) 'Using Transcorrelated method on UEG'

             IF (.not. near_zero(FUEGRS)) THEN
                WRITE(6,'(A,F20.16)') '  Electron Gas Rs set to ',FUEGRS
                WRITE(6,'(A,I4)') '  DIMENSION set to ',dimen
                if(dimen==3) then
                 OMEGA=BOX*BOX*BOX*BOA*COA
!C.. required density is (3/(4 pi rs^3))
!C.. need omega to be (NEL* 4 pi rs^3 / 3)
!C.. need box to be (NEL*4 pi/(3 BOA COA))^(1/3) rs
                 BOX=(NEL*4.0_dp*PI/(3.0_dp*BOA*COA))**(1.0_dp/3.0_dp)
                 BOX=BOX*FUEGRS
                 WRITE(6,'(A, F20.16)') "  Resetting box size to ", BOX
                else if(dimen==2) then
                 OMEGA=BOX*BOX*BOA
!C.. required density is (1/( pi rs^2))
!C.. need omega to be (NEL* pi rs^2)
!C.. need box to be (NEL* pi/ BOA)^(1/2) rs
                 BOX=(NEL*PI/BOA)**(1.0_dp/2.0_dp)
                 BOX=BOX*FUEGRS
                 WRITE(6,'(A, F20.16)') "  Resetting box size to ", BOX
                else
                 WRITE(6,'(A, I4)') " Dimension problem  ", dimen
                stop
                end if

             ENDIF
          ENDIF
          IF(THUB) WRITE(6,'(A)') '  *** HUBBARD MODEL ***  '
!C..
          IF(.NOT.THUB.AND..NOT.TUEG) THEN
             ! Just have even/odd symmetry so it's a two cycle symmetry
             ! generation issue.
             TwoCycleSymGens = .true.
             WRITE(6,'(A)') "  Electron in cubic box.  "
             IF(TPARITY) THEN
                WRITE(6,'(A)') '  *******************************  '
                WRITE(6,*) ' PARITY IS ON '
                DO I=1,3
                   IF(IPARITY(I).EQ.1) THEN
                      CPAR(I)='G'
                    ELSEIF(IPARITY(I).EQ.-1) THEN
                      CPAR(I)='U'
                    ELSE
                      call stop_all(this_routine, ' !!! PROBLEM WITH PARITY !!! ')
                    ENDIF
                ENDDO
                CPARITY=CPAR(1)//CPAR(2)//CPAR(3)
                WRITE(6,*) ' PARITY : ' , CPARITY
             ELSE
                WRITE(6,*) ' PARITY IS OFF '
             ENDIF
             WRITE(6,*) ' ******************************* '

!  //TBR
!         IF((.NOT.TBLOCK).AND.(.NOT.TCALCHMAT.OR.NTAY.LT.0)) STOP 'CANNOT USE PARITY WITHOUT LIST OF DETS'
          ELSE
             IF(TPARITY) THEN
                WRITE(6,*) ' MOMENTUM : ',(IPARITY(I),I=1,3)
             ENDIF
          ENDIF
!C..

         ! W.D: are those variable ever used actually?
          NMAX=MAX(NMAXX,NMAXY,NMAXZ)
          NNR=NMSH*NMSH*NMSH
          WRITE(6,'(A,I5)') '  NMAXX : ' , NMAXX
          WRITE(6,'(A,I5)') '  NMAXY : ' , NMAXY
          WRITE(6,'(A,I5)') '  NMAXZ : ' , NMAXZ
          WRITE(6,'(A,I5)') '  NMSH : ' , NMSH
!C.. 2D check
          IF(NMAXZ.EQ.0) THEN
             WRITE(6,'(A)') ' NMAXZ=0.  2D calculation using C/A=1/A  '
             COA=1/BOX
          ENDIF

!C..
          IF(THUB) THEN
             WRITE(6,'(1X,A,F19.5)') '  HUBBARD T : ' , BHUB
             WRITE(6,'(1X,A,F19.5)') '  HUBBARD U : ' , UHUB
             if (abs(nn_bhub) > EPS) then
                 WRITE(6,'(1X,A,F19.5)') '  HUBBARD T* : ' , nn_bhub
                 print *, "Also next-nearest neighbor hopping!"
             end if
             IF(TTILT) WRITE(6,*) ' TILTED LATTICE: ',ITILTX, ",",ITILTY
             IF(TTILT.AND.ITILTX.GT.ITILTY) call stop_all(this_routine, 'ERROR: ITILTX>ITILTY')
             if (t_new_hubbard) then
                 if (iprocindex == root) then
                     print *, "New Hubbard Implementation! "
                     print *, "lattice used: "
                     call lat%print_lat()
                 end if
             end if
          ELSE
             WRITE(6,'(1X,A,F19.5)') '  BOX LENGTH : ' , BOX
             WRITE(6,'(1X,A,F19.5)') '  B/A : ' , BOA
             WRITE(6,'(1X,A,F19.5)') '  C/A : ' , COA
             TTILT=.FALSE.
          ENDIF
          ALAT(1)=BOX
          ALAT(2)=BOX*BOA
          ALAT(3)=BOX*COA
          IF(near_zero(fRc) .AND. iPeriodicDampingType /= 0) THEN
             ALAT(4)=BOX*((BOA*COA)/(4*PI/3))**THIRD
          ELSE
             ALAT(4)=fRc
          ENDIF
!      ALAT(4)=2*BOX*(BOA*COA)**(1/3.0_dp)

          IF(THUB) THEN
              if (t_new_hubbard) then
                 omega = real(lat%get_nsites(), dp)
                 if (iprocindex == root) then
                     print *, " periodic boundary conditions: ", lat%is_periodic()
                     print *, "Real space basis: ", t_new_real_space_hubbard
                 end if
             else
                 WRITE(6,*) ' X-LENGTH OF HUBBARD CHAIN:', NMAXX
                 WRITE(6,*) ' Y-LENGTH OF HUBBARD CHAIN:', NMAXY
                 WRITE(6,*) ' Z-LENGTH OF HUBBARD CHAIN:', NMAXZ
                 WRITE(6,*) ' Periodic Boundary Conditions:',TPBC
                 WRITE(6,*) ' Real space basis:',TREAL

             IF(TTILT.AND.THUB) THEN
                OMEGA=real(NMAXX,dp)*NMAXY*(ITILTX*ITILTX+ITILTY*ITILTY)
             ELSE
                OMEGA=real(NMAXX,dp)*(NMAXY)*(NMAXZ)
             ENDIF
             end if
             RS=1.0_dp
          ELSE
            if(dimen==3)then


             OMEGA=ALAT(1)*ALAT(2)*ALAT(3)
             RS=(3.0_dp*OMEGA/(4.0_dp*PI*NEL))**THIRD
             ALAT(5)=RS
             IF(iPeriodicDampingType.NE.0) THEN
                IF(iPeriodicDampingType.EQ.1) THEN
                   WRITE(6,*) " Using attenuated Coulomb potential for exchange interactions."
                ELSEIF(iPeriodicDampingType.EQ.2) THEN
                   WRITE(6,*) " Using cut-off Coulomb potential for exchange interactions."
                ENDIF

                WRITE(6,*) " Rc cutoff: ",ALAT(4)
             ENDIF
             WRITE(6,*) " Wigner-Seitz radius Rs=",RS
             FKF=(9*PI/4)**THIRD/RS
             WRITE(6,*) " Fermi vector kF=",FKF
             WRITE(6,*) " Fermi Energy EF=",FKF*FKF/2
             WRITE(6,*) " Unscaled Fermi Energy nmax**2=",(FKF*FKF/2)/(0.5*(2*PI/ALAT(5))**2)
            else if(dimen==2)then
             OMEGA=ALAT(1)*ALAT(2)
             RS=dsqrt(OMEGA/(PI*NEL))
             ALAT(5)=RS
             IF(iPeriodicDampingType.NE.0) THEN
                IF(iPeriodicDampingType.EQ.1) THEN
                   WRITE(6,*) " Using attenuated Coulomb potential for exchange interactions."
                ELSEIF(iPeriodicDampingType.EQ.2) THEN
                   WRITE(6,*) " Using cut-off Coulomb potential for exchange interactions."
                ENDIF

                WRITE(6,*) " Rc cutoff: ",ALAT(4)
             ENDIF
             WRITE(6,*) " Wigner-Seitz radius Rs=",RS
             FKF=dsqrt(2.0_dp)/RS
             WRITE(6,*) " Fermi vector kF=",FKF
             WRITE(6,*) " Fermi Energy EF=",FKF*FKF/2
             WRITE(6,*) " Unscaled Fermi Energy nmax**2=",(FKF*FKF/2)/(0.5*(2*PI/ALAT(5))**2) !??????????

            else if (dimen==1) then
                OMEGA=ALAT(1)
                RS=OMEGA/(2*NEL)
                ALAT(5)=RS
                write(6,*)'Be cautios, the 1D rs and kF values have not been checked thorougHly!'
                IF(iPeriodicDampingType.NE.0) THEN
                IF(iPeriodicDampingType.EQ.1) THEN
                   WRITE(6,*) " Using attenuated Coulomb potential for exchange interactions."
                ELSEIF(iPeriodicDampingType.EQ.2) THEN
                   WRITE(6,*) " Using cut-off Coulomb potential for exchange interactions."
                ENDIF

                WRITE(6,*) " Rc cutoff: ",ALAT(4)
             ENDIF
             WRITE(6,*) " Wigner-Seitz radius Rs=",RS
             FKF=PI*RS
             WRITE(6,*) " Fermi vector kF=",FKF
             WRITE(6,*) " Fermi Energy EF=",FKF*FKF/2
             WRITE(6,*) " Unscaled Fermi Energy nmax**2=",(FKF*FKF/2)/(0.5*(2*PI/ALAT(5))**2)

            else
                write(6,*)'Dimension problem'
                stop
            end if


          ENDIF
          IF(.not. (OrbECutoff .isclose. 1e-20_dp)) WRITE(6,*) " Orbital Energy Cutoff:",OrbECutoff
          WRITE(6,'(1X,A,F19.5)') '  VOLUME : ' , OMEGA
          WRITE(6,*) ' TALPHA : ' , TALPHA
          WRITE(6,'(1X,A,F19.5)') '  ALPHA : ' , ALPHA
          ALPHA=MIN(ALAT(1),ALAT(2),ALAT(3))*ALPHA
          WRITE(6,'(1X,A,F19.5)') '  SCALED ALPHA : ' , ALPHA

!C..
!C..Calculate number of basis functions
!C.. UEG allows from -NMAX->NMAX
          IF(TSPINPOLAR) THEN
             NBASISMAX(4,1)=1
!C.. spinskip
             NBASISMAX(2,3)=1
          ELSE
             NBASISMAX(4,1)=-1
!C.. spinskip
!  If spinskip is unset
             IF(NBASISMAX(2,3).EQ.0) NBASISMAX(2,3)=2
          ENDIF
          NBASISMAX(4,2)=1
          IF(THUB) THEN
              if (t_new_hubbard) then
                ! i need a new setup routine for this for the new generic
                ! hubbard setup! essentialy this just sets up NBASISMAX ..
                ! what do i need from this legacy variable??
                len = 2*lat%get_nsites()
                if (t_k_space_hubbard) then
                    ! this indicates pbc and k-space.
                    ! especially for the addelecsym function!
                    NBASISMAX(1,3) = 0

                else if (t_new_real_space_hubbard) then
                    NBASISMAX(1,3) = 4
                    NBASISMAX(3,3) = 0
                end if

              else
                 IF(TTILT) THEN
                    CALL SETBASISLIM_HUBTILT(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,TPBC,ITILTX,ITILTY)
                    ! is supported now!
    !                 IF(TREAL) call stop_all(this_routine, 'REAL TILTED HUBBARD NOT SUPPORTED')
                  ELSE
                    CALL SETBASISLIM_HUB(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,TPBC,TREAL)
                 ENDIF
             end if
          ELSEIF(TUEG) THEN
             NBASISMAX(1,1)=-NMAXX
             NBASISMAX(1,2)=NMAXX
             NBASISMAX(2,1)=-NMAXY
             NBASISMAX(2,2)=NMAXY
             NBASISMAX(3,1)=-NMAXZ
             NBASISMAX(3,2)=NMAXZ
             NBASISMAX(1,3)=-1
             if(dimen==3)then
              LEN=(2*NMAXX+1)*(2*NMAXY+1)*(2*NMAXZ+1)*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
             else if (dimen==2)then
              LEN=(2*NMAXX+1)*(2*NMAXY+1)*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
             else if (dimen ==1) then
                LEN=(2*NMAXX+1)*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
             else
              WRITE(6,'(A, I4)') " Dimension problem  ", dimen
              stop
             end if
!C.. UEG
             NBASISMAX(3,3)=-1

          ELSE
             NBASISMAX(1,1)=1
             NBASISMAX(1,2)=NMAXX
             NBASISMAX(2,1)=1
             NBASISMAX(2,2)=NMAXY
             NBASISMAX(3,1)=1
             NBASISMAX(3,2)=NMAXZ
             NBASISMAX(1,3)=0
             LEN=NMAXX*NMAXY*NMAXZ*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
             NBASISMAX(1,3)=0
!C.. particle in box
             NBASISMAX(3,3)=-2
             tAbelian=.true.
          ENDIF
      ENDIF
!C..         (.NOT.TREADINT)

     if (t_new_hubbard) then
         ! [W.D. 25.1.2018]
         ! ignore the old thub keyword and try to set everything up
         ! in a standalone fashion for the new hubbard implementation,
         ! since otherwise this causes a lot of conflict with other
         ! assumptions in the old implementation!
         ! todo!

     end if

!C.. we actually store twice as much in arr as we need.
!C.. the ARR(1:LEN,1) are the energies of the orbitals ordered according to
!C.. BRR.  ARR(1:LEN,2) are the energies of the orbitals with default
!C.. ordering.
!C.. ARR is reallocated in IntFreezeBasis if orbitals are frozen so that it
!C.. has the correct size and shape to contain the eigenvalues of the active
!C.. basis.
      WRITE(6,'(A,I5)') "  NUMBER OF SPIN ORBITALS IN BASIS : ", Len
      Allocate(Arr(LEN,2),STAT=ierr)
      LogAlloc(ierr,'Arr',2*LEN,8,tagArr)
! // TBR
!      IP_ARRSTORE=IP_ARR
      ARR=0.0_dp
      Allocate(Brr(LEN),STAT=ierr)
      LogAlloc(ierr,'Brr',LEN,4,tagBrr)
      BRR(1:LEN)=0
      Allocate(G1(Len),STAT=ierr)
      LogAlloc(ierr,'G1',LEN,BasisFNSizeB,tagG1)
      G1(1:LEN)=NullBasisFn
      IF(TCPMD) THEN
         WRITE(6,'(A)') '*** INITIALIZING BASIS FNs FROM CPMD ***'
         CALL CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
         NBASIS=LEN
         iSpinSkip=NBasisMax(2,3)
      ELSEIF(tVASP) THEN
         WRITE(6,'(A)') '*** INITIALIZING BASIS FNs FROM VASP ***'
         CALL VASPBasisInit(ARR,BRR,G1,LEN) ! This also modifies nBasisMax
         NBASIS=LEN
         iSpinSkip=NBasisMax(2,3)
      ELSEIF(TREADINT.AND.TDFREAD) THEN
         WRITE(6,'(A)') '*** Creating Basis Fns from Dalton output ***'
         call InitDaltonBasis(Arr,Brr,G1,Len)
         nBasis=Len
         call GenMolpSymTable(1,G1,nBasis)
      ELSEIF(TREADINT) THEN
!This is also called for tRiIntegrals and tCacheFCIDUMPInts
         WRITE(6,'(A)') '*** CREATING BASIS FNs FROM FCIDUMP ***'
         CALL GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN,TBIN)
         NBASIS=LEN
!C.. we're reading in integrals and have a molpro symmetry table
         IF(lNoSymmetry) THEN
            WRITE(6,*) "Turning Symmetry off"
            DO I=1,nBasis
               G1(I)%Sym%s=0
            ENDDO
            CALL GENMOLPSYMTABLE(1,G1,NBASIS)
            DO I=1,nBasis
               G1(I)%Sym%s=0
            ENDDO
         ELSE
            CALL GENMOLPSYMTABLE(NBASISMAX(5,2)+1,G1,NBASIS)
         ENDIF
      ELSE
!C.. Create plane wave basis functions
        if (treal) then
         WRITE(6,*) "Creating real-space basis."
        else
         WRITE(6,*) "Creating plane wave basis."
        end if
        if (t_new_hubbard) then
            BRR = [(i, i = 1, 2*lat%get_nsites())]
            IG = 2*lat%get_nsites()

            if (t_new_real_space_hubbard) then
                ! i have to do everything what is done below with my new
                ! lattice class:
                ! something like: (loop over spin-orbital!)
                ARR = 0.0_dp
            else
                ARR(:,1) = [(bhub*lat%dispersion_rel_spin_orb(i), i = 1, IG)]


                ARR(:,2) = [(bhub*lat%dispersion_rel_spin_orb(i), i = 1, IG)]
            end if

            do i = 1, lat%get_nsites()
                ! todo: not sure if i really should set up the k-vectors
                ! for the real-space! maybe the convention actually is to
                ! set them all to 0! check that!
                ! do i still want to save the "real-space positions" in the
                ! k-vectors? i could.. but do i need it?
                G1(2*i-1)%k = lat%get_k_vec(i)
                G1(2*i-1)%ms = -1
                G1(2*i-1)%Sym = TotSymRep()

                G1(2*i)%k = lat%get_k_vec(i)
                G1(2*i)%ms = 1
                G1(2*i)%Sym = TotSymRep()
                ! and in the k-space i still need to
            end do
        else if(dimen==3) then

         IG=0
         DO I=NBASISMAX(1,1),NBASISMAX(1,2)
           DO J=NBASISMAX(2,1),NBASISMAX(2,2)
             DO K=NBASISMAX(3,1),NBASISMAX(3,2)
               DO L=NBASISMAX(4,1),NBASISMAX(4,2),2
                  G%k(1)=I
                  G%k(2)=J
                  G%k(3)=K
                  G%Ms=L
                  ! change to implement the tilted real-space
                  if ((treal .and. .not. ttilt) .or. KALLOWED(G,nBasisMax))then
!                   IF((THUB.AND.(TREAL.OR..NOT.TPBC)).OR.KALLOWED(G,NBASISMAX)) THEN
                    IF(THUB) THEN
!C..Note for the Hubbard model, the t is defined by ALAT(1)!
                       call setupMomIndexTable()
                       call setupBreathingCont(2*btHub/OMEGA)
                       IF(TPBC) THEN
                       CALL HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ELSE
                      CALL HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ENDIF
                    ELSEIF(TUEG) THEN
                       CALL GetUEGKE(I,J,K,ALAT,tUEGTrueEnergies,tUEGOffset,k_offset,SUM,dUnscaledE)
!                      Bug for  non-trueEnergy
!                      IF(dUnscaledE.gt.OrbECutoff) CYCLE
                       IF(tUEGTrueEnergies)then
                        IF(dUnscaledE.gt.OrbECutoff) CYCLE
                       ELSE
                        IF(SUM.gt.OrbECutoff) CYCLE
                       END IF

                    ELSE
                       SUM=(BOX**2)*((I*I/ALAT(1)**2)+(J*J/ALAT(2)**2)+(K*K/ALAT(3)**2))
                    ENDIF
                    IF(.NOT.TUEG.AND.SUM.GT.OrbECutoff) CYCLE
                    IG=IG+1
                    ARR(IG,1)=SUM
                    ARR(IG,2)=SUM
                    BRR(IG)=IG
!C..These are the quantum numbers: n,l,m and sigma
                    G1(IG)%K(1)=I
                    G1(IG)%K(2)=J
                    G1(IG)%K(3)=K
                    G1(IG)%MS=L
                    G1(IG)%Sym=TotSymRep()

                  ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
        else if(dimen==2) then
          IG=0
          k2_max_2d=int(OrbECutoff)+1
         do  k=0, k2_max_2d
         DO I=NBASISMAX(1,1),NBASISMAX(1,2)
           DO J=NBASISMAX(2,1),NBASISMAX(2,2)
               DO L=NBASISMAX(4,1),NBASISMAX(4,2),2
                  G%k(1)=I
                  G%k(2)=J
                  G%k(3)=0
                  G%Ms=L
                  ! change to implement the tilted real-space
                  if ((treal .and. .not. ttilt) .or. KALLOWED(G,nBasisMax))then
!                   IF((THUB.AND.(TREAL.OR..NOT.TPBC)).OR.KALLOWED(G,NBASISMAX)) THEN
                    IF(THUB) THEN
!C..Note for the Hubbard model, the t is defined by ALAT(1)!
                       IF(TPBC) THEN
                       CALL HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ELSE
                      CALL HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ENDIF
                    ELSEIF(TUEG) THEN
                       CALL GetUEGKE(I,J,0,ALAT,tUEGTrueEnergies,tUEGOffset,k_offset,SUM,dUnscaledE)
!                      Bug for  non-trueEnergy
!                      IF(dUnscaledE.gt.OrbECutoff) CYCLE
                       IF(tUEGTrueEnergies)then
                        IF(dUnscaledE.gt.OrbECutoff) CYCLE
                       ELSE
                        if((sum.gt.(k*1.d0+1.d-20)).or.(sum.lt.(k*1.d0-1.d-20)))CYCLE
                        IF(SUM.gt.OrbECutoff) CYCLE
                       END IF

                    ELSE
                       SUM=(BOX**2)*((I*I/ALAT(1)**2)+(J*J/ALAT(2)**2))
                    ENDIF
                    IF(.NOT.TUEG.AND.SUM.GT.OrbECutoff) CYCLE
                    IG=IG+1
                    ARR(IG,1)=SUM
                    ARR(IG,2)=SUM
                    BRR(IG)=IG
!C..These are the quantum numbers: n,l,m and sigma
                    G1(IG)%K(1)=I
                    G1(IG)%K(2)=J
                    G1(IG)%K(3)=0
                    G1(IG)%MS=L
                    G1(IG)%Sym=TotSymRep()
                  ENDIF
               ENDDO
           ENDDO
         ENDDO
         end do
!C..Check to see if all's well
        else if (dimen==1) then
         IG=0
         DO I=NBASISMAX(1,1),NBASISMAX(1,2)
           DO J=NBASISMAX(2,1),NBASISMAX(2,2)
             DO K=NBASISMAX(3,1),NBASISMAX(3,2)
               DO L=NBASISMAX(4,1),NBASISMAX(4,2),2
                  G%k(1)=I
                  G%k(2)=J
                  G%k(3)=K
                  G%Ms=L
                  ! change to implement the tilted real-space
                  if ((treal .and. .not. ttilt) .or. KALLOWED(G,nBasisMax))then
!                   IF((THUB.AND.(TREAL.OR..NOT.TPBC)).OR.KALLOWED(G,NBASISMAX))
!                   THEN
                    IF(THUB) THEN
!C..Note for the Hubbard model, the t is defined by ALAT(1)!
                       call setupMomIndexTable()
                       call setupBreathingCont(2*btHub/OMEGA)
                       IF(TPBC) THEN
                       CALL HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ELSE
                      CALL HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
                       ENDIF
                    ELSEIF(TUEG) THEN
                       CALL GetUEGKE(I,J,K,ALAT,tUEGTrueEnergies,tUEGOffset,k_offset,SUM,dUnscaledE)
!                      Bug for  non-trueEnergy
!                      IF(dUnscaledE.gt.OrbECutoff) CYCLE
                       IF(tUEGTrueEnergies)then
                        IF(dUnscaledE.gt.OrbECutoff) CYCLE
                       ELSE
                        IF(SUM.gt.OrbECutoff) CYCLE
                       END IF

                    ELSE
                       SUM=(BOX**2)*((I*I/ALAT(1)**2)+(J*J/ALAT(2)**2)+(K*K/ALAT(3)**2))
                    ENDIF
                    IF(.NOT.TUEG.AND.SUM.GT.OrbECutoff) CYCLE
                    IG=IG+1
                    ARR(IG,1)=SUM
                    ARR(IG,2)=SUM
                    BRR(IG)=IG
!C..These are the quantum numbers: n,l,m and sigma
                    G1(IG)%K(1)=I
                    G1(IG)%K(2)=J
                    G1(IG)%K(3)=K
                    G1(IG)%MS=L
                    G1(IG)%Sym=TotSymRep()

                  ENDIF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
        else
              WRITE(6,'(A, I4)') " Dimension problem  ", dimen
              stop
        end if
!C..Check to see if all's well
         WRITE(6,*) ' NUMBER OF BASIS FUNCTIONS : ' , IG
         NBASIS=IG
         IF(LEN.NE.IG) THEN
            IF(OrbECutoff.gt.-1e20_dp) then
               write(6,*) " Have removed ", LEN-IG, " high energy orbitals "
               ! Resize arr and brr.
               allocate(arr_tmp(nbasis,2),brr_tmp(nbasis),stat=ierr)
               arr_tmp = arr(1:nbasis,:)
               brr_tmp = brr(1:nbasis)
               deallocate(arr,brr,stat=ierr)
               LogDealloc(tagarr)
               LogDealloc(tagbrr)
               allocate(arr(nbasis,2),brr(nbasis),stat=ierr)
               LogAlloc(ierr,'Arr',2*nbasis,8,tagArr)
               LogAlloc(ierr,'Brr',nbasis,4,tagBrr)
               arr = arr_tmp
               brr = brr_tmp
               deallocate(arr_tmp, brr_tmp, stat=ierr)
            else
               WRITE(6,*) " LEN=",LEN,"IG=",IG
               call stop_all(this_routine, ' LEN NE IG ')
            endif
         ENDIF
         ! also turn off symmetries for the real-space basis
         if (thub .and. treal) then
             WRITE(6,*) "Turning Symmetry off for the real-space Hubbard"
             DO I = 1,nBasis
                 G1(I)%Sym%s = 0
             end do
             call GENMOLPSYMTABLE(1,G1,NBASIS)
             DO I = 1, nBasis
                 G1(I)%Sym%s = 0
             end do
         end if
         if (.not. tHub) CALL GENMOLPSYMTABLE(1,G1,NBASIS)
      ENDIF

      IF(tFixLz) THEN
          WRITE(6,'(A)') "****** USING Lz SYMMETRY *******"
          WRITE(6,'(A,I5)') "Pure spherical harmonics with complex orbitals used to constrain Lz to: ",LzTot
          WRITE(6,*) "Due to the breaking of the Ml degeneracy, the fock energies are slightly wrong, "&
          &//"on order of 1.0e-4_dp - do not use for MP2!"
          if(nsymlabels.gt.4) then
              call stop_all(this_routine,"D2h point group detected. Incompatable with Lz symmetry conserving "&
              &//"orbitals. Have you transformed these orbitals into spherical harmonics correctly?!")
          endif
      ENDIF

!C..        (.NOT.TREADINT)
!C.. Set the initial symmetry to be totally symmetric
      FrzSym=NullBasisFn
      FrzSym%Sym=TotSymRep()
      CALL SetupFreezeSym(FrzSym)
!C..Now we sort them using SORT2 and then SORT

!C.. This sorts ARR and BRR into order of ARR [AJWT]
      IF(.NOT.THFNOORDER) THEN
          CALL ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
      ELSE
!.. copy the default ordered energies.
          CALL DCOPY(NBASIS,ARR(1,1),1,ARR(1,2),1)
      ENDIF
!      WRITE(6,*) THFNOORDER, " THFNOORDER"
      if(.not.tMolpro) then
          !If we are calling from molpro, we write the basis later (after reordering)
          CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)
      endif
      IF(NEL.GT.NBASIS) &
          call stop_all(this_routine, 'MORE ELECTRONS THAN BASIS FUNCTIONS')
      CALL neci_flush(6)
      IF(TREAL.AND.THUB) THEN
!C.. we need to allow integrals between different spins
         NBASISMAX(2,3)=1
      ENDIF

      NOCC=NEl/2
      IF(TREADINT) THEN
!C.. we're reading in integrals and have a molpro symmetry table
         IF(lNoSymmetry) THEN
            WRITE(6,*) "Turning Symmetry off"
            CALL GENMOLPSYMREPS()
         ELSE
            CALL GENMOLPSYMREPS()
         ENDIF
      ELSEIF(TCPMD) THEN
!C.. If TCPMD, then we've generated the symmetry table earlier,
!C.. but we still need the sym reps table.
         CALL GENCPMDSYMREPS(G1,NBASIS,ARR)
      ELSEIF(tVASP) THEN
!C.. If VASP-based calculation, then we've generated the symmetry table earlier,
!C.. but we still need the sym reps table. DEGENTOL=1.0e-6_dp. CHECK w/AJWT.
         CALL GENSYMREPS(G1,NBASIS,ARR,1.e-6_dp)
      ELSEIF(THUB.AND..NOT.TREAL) THEN
          ! [W.D. 25.1.2018:]
          ! this also has to be changed for the new hubbard implementation
          if (t_k_space_hubbard) then
              ! or i change the function below to account for the new
              ! implementation
              call setup_symmetry_table()
!               call gen_symreps()


          else
             CALL GenHubMomIrrepsSymTable(G1,nBasis,nBasisMax)
         end if
         ! this function does not make sense..
         CALL GENHUBSYMREPS(NBASIS,ARR,BRR)
         CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)
      ELSE
!C.. no symmetry, so a simple sym table
         CALL GENMOLPSYMREPS()
      ENDIF

!C..
!// TBR
!      WRITE(6,*) ' TREAD : ' , TREAD


!// TBR
!      WRITE(6,*) ' ETRIAL : ',ETRIAL
!      IF(FCOUL.NE.1.0_dp)  WRITE(6,*) "WARNING: FCOUL is not 1.0_dp. FCOUL=",FCOUL
      IF(FCOULDAMPBETA.GT.0) WRITE(6,*) "FCOUL Damping.  Beta ",FCOULDAMPBETA," Mu ",FCOULDAMPMU
      call halt_timer(proc_timer)
    End Subroutine SysInit



    Subroutine SysCleanup()

      use sym_mod, only: EndSym

      CALL ENDSYM()

    End Subroutine SysCleanup



    logical function AreSameSpatialOrb(i,j)
      ! Test whether spin orbitals i and j are from the same spatial orbital.
      ! Returns true if i and j are the *same* spin orbital or are the alpha
      ! and beta spin orbitals of the same spatial orbital.
      integer :: i,j
      integer :: a,b
      AreSameSpatialOrb=.false.
      if (i.eq.j) then
          AreSameSpatialOrb=.true.
      else
          a=min(i,j)
          b=max(i,j)
          if (G1(a)%Ms.eq.-1.and.b-a.eq.1) then
              ! a is the alpha and b is the beta of the same spatial orbital.
              AreSameSpatialOrb=.true.
          end if
      end if
    end function AreSameSpatialOrb


!================================================================
!         UEG2 Subroutines
!================================================================
SUBROUTINE LatticeInit(RS, FKF)
 !  initiates  the reciprocal lattice, real volume and the Fermi vector

    real(dp), intent(out) :: RS, FKF

    !   check dimension
    if(dimen==3) then ! 3D
        OMEGA=4.0_dp/3.0_dp*PI*FUEGRS**3*NEL
        RS=(3.0_dp*OMEGA/(4.0_dp*PI*NEL))**THIRD
        FKF=(9*PI/4)**THIRD/RS
    ! define  lattice vectors and lattice constant in reciprocal space
        if (real_lattice_type == "sc") then
            k_lattice_constant = 2.0_dp*PI/OMEGA**THIRD
            k_lattice_vectors(1,1:3) = (/1, 0, 0 /)
            k_lattice_vectors(2,1:3) = (/0, 1, 0 /)
            k_lattice_vectors(3,1:3) = (/0, 0, 1 /)
        else if (real_lattice_type == "bcc") then
            k_lattice_constant = 2.0_dp*PI/(2.0_dp*OMEGA)**THIRD
            k_lattice_vectors(1,1:3) = (/0, 1, 1 /)
            k_lattice_vectors(2,1:3) = (/1, 0, 1 /)
            k_lattice_vectors(3,1:3) = (/1, 1, 0 /)
        else if (real_lattice_type == "fcc") then
            k_lattice_constant =2.0_dp*PI/(4.0_dp*OMEGA)**THIRD
            k_lattice_vectors(1,1:3) = (/-1, 1, 1 /)
            k_lattice_vectors(2,1:3) = (/1, -1, 1 /)
            k_lattice_vectors(3,1:3) = (/1, 1, -1 /)
        else
            write(6,'(A)')  'lattice type not valid'
        end if
    else if (dimen==2) then !2D
        write(6,'(A)') ' NMAXZ=0 : 2D calculation'
        OMEGA=PI*FUEGRS**2*NEL
        RS=(OMEGA/(PI*NEL))**(1.0_dp/2.0_dp)
        FKF=sqrt(2.0_dp)/RS
        ! define  lattice vectors and lattice constant in reciprocal space
        k_lattice_constant = 2.0_dp*PI/OMEGA**(1.0_dp/2.0_dp)
        k_lattice_vectors(1,1:3) = (/1, 0, 0 /)
        k_lattice_vectors(2,1:3) = (/0, 1, 0 /)
        k_lattice_vectors(3,1:3) = (/0, 0, 0 /)
    else if (dimen==1) then !1D
        write(6,'(A)') ' NMAXZ=0,  NMAXY=0 : 1D calculation'
        OMEGA=2.0_dp*FUEGRS*NEL
        RS=OMEGA/(2.0_dp*NEL)
        FKF=(PI/2.0_dp)/RS  !for spin polarised simulation
        ! define  lattice vectors and lattice constant in reciprocal space
        k_lattice_constant = 2.0_dp*PI/OMEGA
        k_lattice_vectors(1,1:3) = (/1, 0, 0 /)
        k_lattice_vectors(2,1:3) = (/0, 0, 0 /)
        k_lattice_vectors(3,1:3) = (/0, 0, 0 /)
    else
        write(6,'(A)') 'Problem with dimension! '
    endif
    return

END SUBROUTINE LatticeInit


SUBROUTINE CalcCell
    !Detemines the cell size for a given cutoff and lattice type

    integer :: ii, jj, kk, EE
    logical :: under_cutoff

    if (real_lattice_type == "sc" .OR. dimen .lt. 3) then
        NMAXX=int(sqrt(orbEcutoff))+1
        if(dimen .gt. 1) NMAXY=int(sqrt(orbEcutoff))+1
        if(dimen .gt. 2) NMAXZ=int(sqrt(orbEcutoff))+1
    else if (real_lattice_type == "bcc" .or. real_lattice_type == "fcc") then
        ! calculate needed cell size
        ii = 0  ! ii is always positiv. jj varies from -ii to ii, kk from -|jj| to |jj|
        under_cutoff = .true.
        do while (ii .le. int(orbEcutoff) .and. under_cutoff) !until  no E < cutoff was found
            under_cutoff = .false.
            jj =-ii
            do while (abs(jj) .le. abs(ii) .and. .not. under_cutoff) !until E < cutoff is found or jj =ii
                kk = -abs(jj)
                do while (abs(kk) .le. abs(jj) .and. .not. under_cutoff)!until E < cutoff is found or kk=jj
                    !calculate unscaled energy for ii, jj, kk
                    EE =(k_lattice_vectors(1,1)*ii+k_lattice_vectors(2,1)*jj+k_lattice_vectors(3,1)*kk)**2
                    EE =EE +(k_lattice_vectors(1,2)*ii+k_lattice_vectors(2,2)*jj+k_lattice_vectors(3,2)*kk)**2
                    EE =EE +(k_lattice_vectors(1,3)*ii+k_lattice_vectors(2,3)*jj+k_lattice_vectors(3,3)*kk)**2
                    if ( EE .le. orbEcutoff) under_cutoff = .true.
                    kk=kk+1
                end do
                  jj = jj+1
            end do
            ii = ii+1
        end do
        NMAXX=ii!-1
        NMAXY=ii!-1
        NMAXZ=ii!-1
    end if ! lattice type
    return

END SUBROUTINE CalcCell


SUBROUTINE CalcTau
    !Detemines tau for a given lattice type

    if(dimen == 3) then ! 3D
        TAU = (k_lattice_constant**2* OMEGA) / (4.0_dp*PI) !Hij_min**-1
        if (tTruncInitiator) TAU = TAU*InitiatorWalkNo
        if (tHPHF) TAU = TAU /sqrt(2.0_dp)
        TAU = 0.9_dp*TAU*4.0_dp/(NEL*(NEL-1))/(NBASIS-NEL)
        if (TAU .gt. k_lattice_constant**(-2)/OrbEcutoff) then
            TAU= 1.0_dp/(k_lattice_constant**(2)*OrbEcutoff)  !using Hii
            write(6,*) '***************** Tau set by using Hii *******************************'
            !write(6,*) 1.0_dp/((2.0_dp*PI/Omega**third)**2*orbEcutoff)
        else
            write(6,*) 'Tau set by using Hji'
        end if

    else if (dimen ==2) then !2D
        TAU = (k_lattice_constant * OMEGA)/(2.0_dp*PI)  !Hij_min**-1
        if (tTruncInitiator) TAU = TAU*InitiatorWalkNo
        if (tHPHF) TAU = TAU /sqrt(2.0_dp)
        TAU = 0.9_dp*TAU*4.0_dp/(NEL*(NEL-1))/(NBASIS-NEL)
        if (TAU .gt. k_lattice_constant**(-2)/OrbEcutoff) then
            !!!!!!!! NOT WORKING YET!!!!!!!
            TAU= 1.0_dp/(k_lattice_constant**(2)*OrbEcutoff)  !using Hii
            write(6,*) '***************** Tau set by using Hii *******************************'
        else
            write(6,*) 'Tau set by using Hji'
        end if

    else if (dimen ==1) then !1D
        TAU = OMEGA/ (-2.0_dp*log(1.0_dp/(2.0_dp*sqrt(orbEcutoff))))
        if (tTruncInitiator) TAU = TAU*InitiatorWalkNo
        if (tHPHF) TAU = TAU /sqrt(2.0_dp)
        TAU = 0.9_dp*TAU*4.0_dp/(NEL*(NEL-1))/(NBASIS-NEL)
        if (TAU .gt. 0.9_dp* 1.0_dp/(0.5_dp*(k_lattice_constant)**2*NEL*OrbEcutoff))  then
            TAU=0.9_dp* 1.0_dp/(0.5_dp*(k_lattice_constant)**2*NEL*OrbEcutoff)   !using Hii
            write(6,*) '***************** Tau set by using Hii *******************************'
        else
            write(6,*) 'Tau set by using Hji'
        end if

    endif !dimension
    write(6, *) 'Tau set to: ', TAU
    return
END SUBROUTINE CalcTau

  function calc_madelung() result(res)

    real(dp)  :: kappa
    integer :: i1, i2, i3, i4
    integer :: n2
    real(dp)  :: k2,ek2,recipsum2
    real(dp) :: t1, modr,er2,realsum2
    real(dp)  :: term2,term4
    integer :: r_lattice_vectors(3, 3)
    real(dp) :: r_lattice_constant
    integer:: kvecX, kvecY, kvecZ
    integer :: tvecX, tvecY, tvecZ
    real(dp) ::  inacc_madelung, temp_sum
    real(dp) :: k2max
    integer :: jj
    real(dp) :: xx, yy, zz
    real(dp) :: res

    inacc_madelung = 1.0d-15 ! Cannot be set lower than 1.0d-15

    !-------------------------------   3D lattice   -----------------------------------------

    kappa=2.8_dp/OMEGA**(1.0_dp/3.0_dp)
    if (real_lattice_type == "sc") then
        r_lattice_constant = OMEGA**THIRD
        r_lattice_vectors(1,1:3) = (/1, 0, 0 /)
        r_lattice_vectors(2,1:3) = (/0, 1, 0 /)
        r_lattice_vectors(3,1:3) = (/0, 0, 1 /)
    else if (real_lattice_type == "fcc") then
        r_lattice_constant =0.5_dp*(2.0_dp*OMEGA)**THIRD
        r_lattice_vectors(1,1:3) = (/-1, 1, 1 /)
        r_lattice_vectors(2,1:3) = (/1, -1, 1 /)
        r_lattice_vectors(3,1:3) = (/1, 1, -1 /)
    else if (real_lattice_type == "bcc") then
        r_lattice_constant =0.5_dp*(4.0_dp*OMEGA)**THIRD
        r_lattice_vectors(1,1:3) = (/0, 1, 1 /)
        r_lattice_vectors(2,1:3) = (/1, 0, 1 /)
        r_lattice_vectors(3,1:3) = (/1, 1, 0 /)
    else
        write(6,'(A)')  'lattice type not valid'
    end if

    term2=-pi/(kappa**2.0_dp*OMEGA)
!     write(6,*) term2, "term2"
    term4=-2.0_dp*kappa/sqrt(pi)
!         write(6,*) term4, "term4"

    recipsum2=0.0_dp
    temp_sum = 0.0_dp
    k2max =0.0_dp
    i4 =1
    do while(temp_sum*(1.0_dp+inacc_madelung) .le. recipsum2)
        temp_sum = recipsum2
        recipsum2 =0.0_dp
        do i1=-i4,i4
            do i2=-i4,i4
                do i3=-i4, i4
                    kvecX=k_lattice_vectors(1,1)*i1+k_lattice_vectors(2,1)*i2+k_lattice_vectors(3,1)*i3
                    kvecY=k_lattice_vectors(1,2)*i1+k_lattice_vectors(2,2)*i2+k_lattice_vectors(3,2)*i3
                    kvecZ=k_lattice_vectors(1,3)*i1+k_lattice_vectors(2,3)*i2+k_lattice_vectors(3,3)*i3
                    n2=kvecX**2+kvecY**2+kvecZ**2
                    k2=(k_lattice_constant/2.0_dp/PI)**2.0_dp*real(n2, dp)
                    ek2=(1.0_dp/OMEGA)*(1.0_dp/(pi*k2))*exp(-pi**2.0_dp*k2/kappa**2.0_dp)
                    if( n2 .gt. k2max) k2max = n2
                    if (n2.ne.0) then
!                         write(6,*) k2,ek2 ! for testing
                        recipsum2=recipsum2+ek2
                    endif
                enddo
            enddo
        enddo
        i4 = i4+1
!             write(6,*) "i4 kmax^2", i4-1, k2max*k_lattice_constant**2
    enddo

    realsum2=0.0_dp
    temp_sum = 0.0_dp
    i4 =1
    do while (temp_sum*(1.0_dp+inacc_madelung) .le. realsum2)
        temp_sum = realsum2
        realsum2=0.0_dp
        do i1=-i4,i4
            do i2=-i4,i4
                do i3=-i4,i4
                    tvecX=r_lattice_vectors(1,1)*i1+r_lattice_vectors(2,1)*i2+r_lattice_vectors(3,1)*i3
                    tvecY=r_lattice_vectors(1,2)*i1+r_lattice_vectors(2,2)*i2+r_lattice_vectors(3,2)*i3
                    tvecZ=r_lattice_vectors(1,3)*i1+r_lattice_vectors(2,3)*i2+r_lattice_vectors(3,3)*i3
                    n2=tvecX**2+tvecY**2+tvecZ**2
                    t1 =r_lattice_constant*sqrt(real(n2, dp))
                    if (.not. near_zero(t1)) then
                        er2=error_function_c(kappa*t1)/t1
                        realsum2=realsum2+er2
                    endif
                enddo
            enddo
        enddo
!             write(6,*) 'i4, realsum2' , i4, realsum2
        i4 = i4+1
    enddo
  ! write(6,*) "real space", realsum2


    !-------------------------------   output   -----------------------------------------
    res =realsum2+recipsum2+term2+term4

    write(6,*) "Calculating Madelung Constant - Fraser et al. PRB 53 4 1814"
    write(6,*) "kappa taken from CASINO manual to be", kappa
    write(6,*) "k2_max", k2max, k2max*k_lattice_constant**2
!     write(6,*) "omega", OMEGA
!     write(6, *) "klatticeconstant", k_lattice_constant
!     write(6, *) "rlatticeconstant", r_lattice_constant
    write(6,*) "Madelung constant", res

    return
end function


END MODULE System

SUBROUTINE WRITEBASIS(NUNIT,G1,NHG,ARR,BRR)
  ! Write out the current basis to unit nUnit
  use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB, k_lattice_vectors
  use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB, nel, tUEG2
  use DeterminantData, only: fdet
  use sym_mod, only: writesym
  use constants, only: dp
  IMPLICIT NONE
  INTEGER NUNIT,NHG,BRR(NHG),I
  integer :: pos
  TYPE(BASISFN) G1(NHG)
  real(dp) ARR(NHG,2), unscaled_energy, kvecX, kvecY, kvecZ

  ! nb. Cannot use EncodeBitDet as would be easy, as nifd, niftot etc are not
  !     filled in yet. --> track pos.
  if (.not. associated(fdet)) &
      write(nunit,'("HF determinant not yet defined.")')
  pos = 1
!=============================================
  if (tUEG2) then

      DO I=1,NHG
!     kvectors in cartesian coordinates
          kvecX=  k_lattice_vectors(1,1)*G1(BRR(I))%K(1) &
                      +k_lattice_vectors(2,1)* G1(BRR(I))%K(2) &
                      +k_lattice_vectors(3,1)*G1(BRR(I))%K(3)
          kvecY=  k_lattice_vectors(1,2)*G1(BRR(I))%K(1) &
                      +k_lattice_vectors(2,2)* G1(BRR(I))%K(2) &
                      +k_lattice_vectors(3,2)*G1(BRR(I))%K(3)
          kvecZ=  k_lattice_vectors(1,3)*G1(BRR(I))%K(1) &
                      +k_lattice_vectors(2,3)* G1(BRR(I))%K(2) &
                      +k_lattice_vectors(3,3)*G1(BRR(I))%K(3)

          unscaled_energy=((kvecX)**2+(kvecY)**2+(kvecZ)**2)

          WRITE(NUNIT,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
          CALL WRITESYM(NUNIT,G1(BRR(I))%SYM,.FALSE.)
          WRITE(NUNIT,'(I4)',advance='no') G1(BRR(I))%Ml
          WRITE(NUNIT,'(3F19.9)', advance='no')  ARR(I,1), ARR(BRR(I),2), unscaled_energy

          if (associated(fdet)) then
              pos=1
              do while (pos < nel .and. fdet(pos) < brr(i))
                  pos = pos + 1
              enddo
              if (brr(i) == fdet(pos)) write (nunit, '(" #")', advance='no')
          endif
          write (nunit,*)
      ENDDO
      RETURN
  end if !UEG2
!=============================================
  DO I=1,NHG
      WRITE(NUNIT,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
      CALL WRITESYM(NUNIT,G1(BRR(I))%SYM,.FALSE.)
      WRITE(NUNIT,'(I4)',advance='no') G1(BRR(I))%Ml
      WRITE(NUNIT,'(2F19.9)', advance='no')  ARR(I,1),ARR(BRR(I),2)
      if (associated(fdet)) then
          pos=1
          do while (pos < nel .and. fdet(pos) < brr(i))
              pos = pos + 1
          enddo
          if (brr(i) == fdet(pos)) write (nunit, '(" #")', advance='no')
      endif
      write (nunit,*)
  ENDDO
  RETURN
END SUBROUTINE WRITEBASIS



SUBROUTINE ORDERBASIS(NBASIS,ARR,BRR,ORBORDER,NBASISMAX,G1)
  use SystemData, only: BasisFN
  use sort_mod
  use util_mod, only: NECI_ICOPY
  use constants, only: dp
  implicit none
  INTEGER NBASIS,BRR(NBASIS),ORBORDER(8,2),nBasisMax(5,*)
  INTEGER BRR2(NBASIS)
  TYPE(BASISFN) G1(NBASIS)
  real(dp) ARR(NBASIS,2),ARR2(NBASIS,2)
  INTEGER IDONE,I,J,IBFN,ITOT,ITYPE,ISPIN
  real(dp) OEN
  character(*), parameter :: this_routine = 'ORDERBASIS'
  IDONE=0
  ITOT=0
!.. copy the default ordered energies.
  CALL DCOPY(NBASIS,ARR(1,1),1,ARR(1,2),1)
  CALL DCOPY(NBASIS,ARR(1,1),1,ARR2(1,2),1)
  WRITE(6,*) ''
  WRITE(6,"(A,8I4)") "Ordering Basis (Closed): ", (ORBORDER(I,1),I=1,8)
  WRITE(6,"(A,8I4)") "Ordering Basis (Open  ): ", (ORBORDER(I,2),I=1,8)
  IF(NBASISMAX(3,3).EQ.1) THEN
!.. we use the symmetries of the spatial orbitals
     DO ITYPE=1,2
        IBFN=1
        DO I=1,8
           DO J=1,ORBORDER(I,ITYPE)
              DO WHILE(IBFN.LE.NBASIS.AND.(G1(IBFN)%SYM%s.LT.I-1.OR.BRR(IBFN).EQ.0))
                 IBFN=IBFN+1
              ENDDO
              IF(IBFN.GT.NBASIS) THEN
                 call stop_all(this_routine, "Cannot find enough basis fns of correct symmetry")
              ENDIF
              IDONE=IDONE+1
              BRR2(IDONE)=IBFN
              BRR(IBFN)=0
              ARR2(IDONE,1)=ARR(IBFN,1)
              IBFN=IBFN+1
           ENDDO
        ENDDO
        ! Beta sort
        call sort (arr2(itot+1:idone,1), brr2(itot+1:idone), nskip=2)
        ! Alpha sort
        call sort (arr2(itot+2:idone,1), brr2(itot+2:idone), nskip=2)
        ITOT=IDONE
     ENDDO
     DO I=1,NBASIS
        IF(BRR(I).NE.0) THEN
           ITOT=ITOT+1
           BRR2(ITOT)=BRR(I)
           ARR2(ITOT,1)=ARR(I,1)
        ENDIF
     ENDDO
     CALL NECI_ICOPY(NBASIS,BRR2,1,BRR,1)
     CALL DCOPY(NBASIS,ARR2(1,1),1,ARR(1,1),1)
  ENDIF
! beta sort
  call sort (arr(idone+1:nbasis,1), brr(idone+1:nbasis), nskip=2)
! alpha sort
  call sort (arr(idone+2:nbasis,1), brr(idone+2:nbasis), nskip=2)
!.. We need to now go through each set of degenerate orbitals, and make
!.. the correct ones are paired together in BRR otherwise bad things
!.. happen in FREEZEBASIS
!.. We do this by ensuring that within a degenerate set, the BRR are in
!.. ascending order
!         IF(NBASISMAX(3,3).EQ.1) G1(3,BRR(1))=J
  DO ISPIN=0,1
     OEN=ARR(1+ISPIN,1)
     J=1+ISPIN
     ITOT=2
     DO I=3+ISPIN,NBASIS,2
        IF(ABS(ARR(I,1)-OEN).GT.1.0e-4_dp) THEN
!.. We don't have degenerate orbitals
!.. First deal with the last set of degenerate orbitals
!.. We sort them into order of BRR
           call sort (brr(i-itot:i-1), arr(i-itot:i-1,1), nskip=2)
!.. now setup the new degenerate set.
           J=J+2
           ITOT=2
        ELSE
           ITOT=ITOT+2
        ENDIF
        OEN=ARR(I,1)
        IF(NBASISMAX(3,3).EQ.1) THEN
!.. If we've got a generic spatial sym or hf we mark degeneracies
!               G(3,BRR(I))=J
        ENDIF
     ENDDO
! i is now nBasis+2
     call sort (brr(i-itot:i-2), arr(i-itot:i-2,1), nskip=2)
  ENDDO
END subroutine ORDERBASIS




!dUnscaledEnergy gives the energy without reference to box size and without any offset.
SUBROUTINE GetUEGKE(I,J,K,ALAT,tUEGTrueEnergies,tUEGOffset,k_offset,Energy,dUnscaledEnergy)

   use SystemData, only: tUEG2, k_lattice_vectors, k_lattice_constant
   use constants, only: Pi, Pi2, THIRD
   use constants, only: dp
   IMPLICIT NONE
   INTEGER I,J,K
   real(dp) ALat(3),k_offset(3),Energy,E
   LOGICAL tUEGOffset, tUEGTrueEnergies
   real(dp) ::  dUnscaledEnergy
   integer :: kvecX, kvecY, kvecZ
   !==================================
   ! initialize unscaled energy for the case of not using tUEGTrueEnergies
   dunscaledEnergy = 0.0_dp
   if (tUEG2) then
      ! kvectors in cartesian coordinates
      kvecX=k_lattice_vectors(1,1)*I+k_lattice_vectors(2,1)*J+k_lattice_vectors(3,1)*K
      kvecY=k_lattice_vectors(1,2)*I+k_lattice_vectors(2,2)*J+k_lattice_vectors(3,2)*K
      kvecZ=k_lattice_vectors(1,3)*I+k_lattice_vectors(2,3)*J+k_lattice_vectors(3,3)*K

       IF(tUEGTrueEnergies) then
           if(tUEGOffset) then
              E=(kvecX+k_offset(1))**2+(kvecY+k_offset(2))**2+(kvecZ+k_offset(3))**2
           else
              E=(kvecX)**2+(kvecY)**2+(kvecZ)**2
           endif
           Energy=0.5_dp*E*k_lattice_constant**2
           dUnscaledEnergy=((kvecX)**2+(kvecY)**2+(kvecZ)**2)
       ELSE
           Energy=((kvecX)**2+(kvecY)**2+(kvecZ)**2)
       ENDIF

       return
   endif
   !==================================
   IF(tUEGTrueEnergies) then
       IF(tUEGOffset) then
          E=((I+k_offset(1))**2/ALAT(1)**2)
          E=E+((J+k_offset(2))**2/ALAT(2)**2)
          E=E+((K+k_offset(3))**2/ALAT(3)**2)
       else
          E=(I*I/ALAT(1)**2)
          E=E+(J*J/ALAT(2)**2)
          E=E+(K*K/ALAT(3)**2)
       endif
       Energy=0.5*4*PI*PI*E
       dUnscaledEnergy=(I*I)
       dUnscaledEnergy=dUnscaledEnergy+(J*J)
       dUnscaledEnergy=dUnscaledEnergy+(K*K)
   ELSE
       E=(I*I)
       E=E+(J*J)
       E=E+(K*K)
       Energy=E
    ENDIF
END SUBROUTINE GetUEGKE
