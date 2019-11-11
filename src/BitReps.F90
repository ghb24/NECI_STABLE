#include "macros.h"

module bit_reps
    use FciMCData, only: CurrentDets, WalkVecDets, MaxWalkersPart, tLogNumSpawns
    use SystemData, only: nel, tCSF, tTruncateCSF, nbasis, csf_trunc_level
    use CalcData, only: tTruncInitiator, tUseRealCoeffs, tSemiStochastic, &
                        tCSFCore, tTrialWavefunction, semistoch_shift_iter, &
                        tStartTrialLater, tPreCond, tReplicaEstimates, tStoredDets

    use csf_data, only: csf_yama_bit, csf_test_bit
    use constants, only: lenof_sign, end_n_int, bits_n_int, n_int, dp,sizeof_int
    use DetBitOps, only: count_open_orbs, CountBits
    use bit_rep_data
    use SymExcitDataMod, only: excit_gen_store_type, tBuildOccVirtList, &
                               tBuildSpinSepLists, &
                               OrbClassCount, ScratchSize, SymLabelList2, &
                               SymLabelCounts2
    use sym_general_mod, only: ClassCountInd
    use global_det_data, only: get_determinant
    implicit none

    ! Structure of a bit representation:

    ! | 0-NIfD: Det | Yamanouchi | Sign(Re) | Sign(Im) | Flags |
    !
    ! -------
    ! (NIfD + 1) * 64-bits              Orbital rep.
    !  NIfY      * 32-bits              Yamanouchi symbol
    !  1         * 32-bits              Signs (Re)
    ! (1         * 32-bits if needed)   Signs (Im)
    ! (1         * 32-bits if needed)   Flags


    interface set_flag
        module procedure set_flag_single
        module procedure set_flag_general
    end interface

    ! Which decoding function do we want to use?
    interface decode_bit_det
!        module procedure decode_bit_det_bitwise
        module procedure decode_bit_det_chunks
        module procedure decode_bit_det_lists
    end interface

    integer, parameter :: l1(1:33)=(/0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,2,1,2,0,0,0/)
    integer, parameter :: l2(1:33)=(/0,0,0,1,3,0,0,0,0,0,0,0,2,1,3,0,0,0,0,0,0,2,2,3,0,0,0,0,0,0,3,1,2/)
    integer, parameter :: l3(1:33)=(/3,0,0,0,0,0,1,4,0,0,0,0,0,0,0,2,1,4,0,0,0,0,0,0,2,2,4,0,0,0,0,0,0/)
    integer, parameter :: l4(1:33)=(/3,1,2,4,0,0,0,0,0,2,3,4,0,0,0,0,0,0,3,1,3,4,0,0,0,0,0,3,2,3,4,0,0/)
    integer, parameter :: l5(1:33)=(/0,0,0,4,1,2,3,4,0,0,0,0,1,5,0,0,0,0,0,0,0,2,1,5,0,0,0,0,0,0,2,2,5/)
    integer, parameter :: l6(1:33)=(/0,0,0,0,0,0,3,1,2,5,0,0,0,0,0,2,3,5,0,0,0,0,0,0,3,1,3,5,0,0,0,0,0/)
    integer, parameter :: l7(1:33)=(/3,2,3,5,0,0,0,0,0,4,1,2,3,5,0,0,0,0,2,4,5,0,0,0,0,0,0,3,1,4,5,0,0/)
    integer, parameter :: l8(1:33)=(/0,0,0,3,2,4,5,0,0,0,0,0,4,1,2,4,5,0,0,0,0,3,3,4,5,0,0,0,0,0,4,1,3/)
    integer, parameter :: l9(1:33)=(/4,5,0,0,0,0,4,2,3,4,5,0,0,0,0,5,1,2,3,4,5,0,0,0,1,6,0,0,0,0,0,0,0/)
    integer, parameter :: l10(1:33)=(/2,1,6,0,0,0,0,0,0,2,2,6,0,0,0,0,0,0,3,1,2,6,0,0,0,0,0,2,3,6,0,0,0/)
    integer, parameter :: l11(1:33)=(/0,0,0,3,1,3,6,0,0,0,0,0,3,2,3,6,0,0,0,0,0,4,1,2,3,6,0,0,0,0,2,4,6/)
    integer, parameter :: l12(1:33)=(/0,0,0,0,0,0,3,1,4,6,0,0,0,0,0,3,2,4,6,0,0,0,0,0,4,1,2,4,6,0,0,0,0/)
    integer, parameter :: l13(1:33)=(/3,3,4,6,0,0,0,0,0,4,1,3,4,6,0,0,0,0,4,2,3,4,6,0,0,0,0,5,1,2,3,4,6/)
    integer, parameter :: l14(1:33)=(/0,0,0,2,5,6,0,0,0,0,0,0,3,1,5,6,0,0,0,0,0,3,2,5,6,0,0,0,0,0,4,1,2/)
    integer, parameter :: l15(1:33)=(/5,6,0,0,0,0,3,3,5,6,0,0,0,0,0,4,1,3,5,6,0,0,0,0,4,2,3,5,6,0,0,0,0/)
    integer, parameter :: l16(1:33)=(/5,1,2,3,5,6,0,0,0,3,4,5,6,0,0,0,0,0,4,1,4,5,6,0,0,0,0,4,2,4,5,6,0/)
    integer, parameter :: l17(1:33)=(/0,0,0,5,1,2,4,5,6,0,0,0,4,3,4,5,6,0,0,0,0,5,1,3,4,5,6,0,0,0,5,2,3/)
    integer, parameter :: l18(1:33)=(/4,5,6,0,0,0,6,1,2,3,4,5,6,0,0,1,7,0,0,0,0,0,0,0,2,1,7,0,0,0,0,0,0/)
    integer, parameter :: l19(1:33)=(/2,2,7,0,0,0,0,0,0,3,1,2,7,0,0,0,0,0,2,3,7,0,0,0,0,0,0,3,1,3,7,0,0/)
    integer, parameter :: l20(1:33)=(/0,0,0,3,2,3,7,0,0,0,0,0,4,1,2,3,7,0,0,0,0,2,4,7,0,0,0,0,0,0,3,1,4/)
    integer, parameter :: l21(1:33)=(/7,0,0,0,0,0,3,2,4,7,0,0,0,0,0,4,1,2,4,7,0,0,0,0,3,3,4,7,0,0,0,0,0/)
    integer, parameter :: l22(1:33)=(/4,1,3,4,7,0,0,0,0,4,2,3,4,7,0,0,0,0,5,1,2,3,4,7,0,0,0,2,5,7,0,0,0/)
    integer, parameter :: l23(1:33)=(/0,0,0,3,1,5,7,0,0,0,0,0,3,2,5,7,0,0,0,0,0,4,1,2,5,7,0,0,0,0,3,3,5/)
    integer, parameter :: l24(1:33)=(/7,0,0,0,0,0,4,1,3,5,7,0,0,0,0,4,2,3,5,7,0,0,0,0,5,1,2,3,5,7,0,0,0/)
    integer, parameter :: l25(1:33)=(/3,4,5,7,0,0,0,0,0,4,1,4,5,7,0,0,0,0,4,2,4,5,7,0,0,0,0,5,1,2,4,5,7/)
    integer, parameter :: l26(1:33)=(/0,0,0,4,3,4,5,7,0,0,0,0,5,1,3,4,5,7,0,0,0,5,2,3,4,5,7,0,0,0,6,1,2/)
    integer, parameter :: l27(1:33)=(/3,4,5,7,0,0,2,6,7,0,0,0,0,0,0,3,1,6,7,0,0,0,0,0,3,2,6,7,0,0,0,0,0/)
    integer, parameter :: l28(1:33)=(/4,1,2,6,7,0,0,0,0,3,3,6,7,0,0,0,0,0,4,1,3,6,7,0,0,0,0,4,2,3,6,7,0/)
    integer, parameter :: l29(1:33)=(/0,0,0,5,1,2,3,6,7,0,0,0,3,4,6,7,0,0,0,0,0,4,1,4,6,7,0,0,0,0,4,2,4/)
    integer, parameter :: l30(1:33)=(/6,7,0,0,0,0,5,1,2,4,6,7,0,0,0,4,3,4,6,7,0,0,0,0,5,1,3,4,6,7,0,0,0/)
    integer, parameter :: l31(1:33)=(/5,2,3,4,6,7,0,0,0,6,1,2,3,4,6,7,0,0,3,5,6,7,0,0,0,0,0,4,1,5,6,7,0/)
    integer, parameter :: l32(1:33)=(/0,0,0,4,2,5,6,7,0,0,0,0,5,1,2,5,6,7,0,0,0,4,3,5,6,7,0,0,0,0,5,1,3/)
    integer, parameter :: l33(1:33)=(/5,6,7,0,0,0,5,2,3,5,6,7,0,0,0,6,1,2,3,5,6,7,0,0,4,4,5,6,7,0,0,0,0/)
    integer, parameter :: l34(1:33)=(/5,1,4,5,6,7,0,0,0,5,2,4,5,6,7,0,0,0,6,1,2,4,5,6,7,0,0,5,3,4,5,6,7/)
    integer, parameter :: l35(1:33)=(/0,0,0,6,1,3,4,5,6,7,0,0,6,2,3,4,5,6,7,0,0,7,1,2,3,4,5,6,7,0,1,8,0/)
    integer, parameter :: l36(1:33)=(/0,0,0,0,0,0,2,1,8,0,0,0,0,0,0,2,2,8,0,0,0,0,0,0,3,1,2,8,0,0,0,0,0/)
    integer, parameter :: l37(1:33)=(/2,3,8,0,0,0,0,0,0,3,1,3,8,0,0,0,0,0,3,2,3,8,0,0,0,0,0,4,1,2,3,8,0/)
    integer, parameter :: l38(1:33)=(/0,0,0,2,4,8,0,0,0,0,0,0,3,1,4,8,0,0,0,0,0,3,2,4,8,0,0,0,0,0,4,1,2/)
    integer, parameter :: l39(1:33)=(/4,8,0,0,0,0,3,3,4,8,0,0,0,0,0,4,1,3,4,8,0,0,0,0,4,2,3,4,8,0,0,0,0/)
    integer, parameter :: l40(1:33)=(/5,1,2,3,4,8,0,0,0,2,5,8,0,0,0,0,0,0,3,1,5,8,0,0,0,0,0,3,2,5,8,0,0/)
    integer, parameter :: l41(1:33)=(/0,0,0,4,1,2,5,8,0,0,0,0,3,3,5,8,0,0,0,0,0,4,1,3,5,8,0,0,0,0,4,2,3/)
    integer, parameter :: l42(1:33)=(/5,8,0,0,0,0,5,1,2,3,5,8,0,0,0,3,4,5,8,0,0,0,0,0,4,1,4,5,8,0,0,0,0/)
    integer, parameter :: l43(1:33)=(/4,2,4,5,8,0,0,0,0,5,1,2,4,5,8,0,0,0,4,3,4,5,8,0,0,0,0,5,1,3,4,5,8/)
    integer, parameter :: l44(1:33)=(/0,0,0,5,2,3,4,5,8,0,0,0,6,1,2,3,4,5,8,0,0,2,6,8,0,0,0,0,0,0,3,1,6/)
    integer, parameter :: l45(1:33)=(/8,0,0,0,0,0,3,2,6,8,0,0,0,0,0,4,1,2,6,8,0,0,0,0,3,3,6,8,0,0,0,0,0/)
    integer, parameter :: l46(1:33)=(/4,1,3,6,8,0,0,0,0,4,2,3,6,8,0,0,0,0,5,1,2,3,6,8,0,0,0,3,4,6,8,0,0/)
    integer, parameter :: l47(1:33)=(/0,0,0,4,1,4,6,8,0,0,0,0,4,2,4,6,8,0,0,0,0,5,1,2,4,6,8,0,0,0,4,3,4/)
    integer, parameter :: l48(1:33)=(/6,8,0,0,0,0,5,1,3,4,6,8,0,0,0,5,2,3,4,6,8,0,0,0,6,1,2,3,4,6,8,0,0/)
    integer, parameter :: l49(1:33)=(/3,5,6,8,0,0,0,0,0,4,1,5,6,8,0,0,0,0,4,2,5,6,8,0,0,0,0,5,1,2,5,6,8/)
    integer, parameter :: l50(1:33)=(/0,0,0,4,3,5,6,8,0,0,0,0,5,1,3,5,6,8,0,0,0,5,2,3,5,6,8,0,0,0,6,1,2/)
    integer, parameter :: l51(1:33)=(/3,5,6,8,0,0,4,4,5,6,8,0,0,0,0,5,1,4,5,6,8,0,0,0,5,2,4,5,6,8,0,0,0/)
    integer, parameter :: l52(1:33)=(/6,1,2,4,5,6,8,0,0,5,3,4,5,6,8,0,0,0,6,1,3,4,5,6,8,0,0,6,2,3,4,5,6/)
    integer, parameter :: l53(1:33)=(/8,0,0,7,1,2,3,4,5,6,8,0,2,7,8,0,0,0,0,0,0,3,1,7,8,0,0,0,0,0,3,2,7/)
    integer, parameter :: l54(1:33)=(/8,0,0,0,0,0,4,1,2,7,8,0,0,0,0,3,3,7,8,0,0,0,0,0,4,1,3,7,8,0,0,0,0/)
    integer, parameter :: l55(1:33)=(/4,2,3,7,8,0,0,0,0,5,1,2,3,7,8,0,0,0,3,4,7,8,0,0,0,0,0,4,1,4,7,8,0/)
    integer, parameter :: l56(1:33)=(/0,0,0,4,2,4,7,8,0,0,0,0,5,1,2,4,7,8,0,0,0,4,3,4,7,8,0,0,0,0,5,1,3/)
    integer, parameter :: l57(1:33)=(/4,7,8,0,0,0,5,2,3,4,7,8,0,0,0,6,1,2,3,4,7,8,0,0,3,5,7,8,0,0,0,0,0/)
    integer, parameter :: l58(1:33)=(/4,1,5,7,8,0,0,0,0,4,2,5,7,8,0,0,0,0,5,1,2,5,7,8,0,0,0,4,3,5,7,8,0/)
    integer, parameter :: l59(1:33)=(/0,0,0,5,1,3,5,7,8,0,0,0,5,2,3,5,7,8,0,0,0,6,1,2,3,5,7,8,0,0,4,4,5/)
    integer, parameter :: l60(1:33)=(/7,8,0,0,0,0,5,1,4,5,7,8,0,0,0,5,2,4,5,7,8,0,0,0,6,1,2,4,5,7,8,0,0/)
    integer, parameter :: l61(1:33)=(/5,3,4,5,7,8,0,0,0,6,1,3,4,5,7,8,0,0,6,2,3,4,5,7,8,0,0,7,1,2,3,4,5/)
    integer, parameter :: l62(1:33)=(/7,8,0,3,6,7,8,0,0,0,0,0,4,1,6,7,8,0,0,0,0,4,2,6,7,8,0,0,0,0,5,1,2/)
    integer, parameter :: l63(1:33)=(/6,7,8,0,0,0,4,3,6,7,8,0,0,0,0,5,1,3,6,7,8,0,0,0,5,2,3,6,7,8,0,0,0/)
    integer, parameter :: l64(1:33)=(/6,1,2,3,6,7,8,0,0,4,4,6,7,8,0,0,0,0,5,1,4,6,7,8,0,0,0,5,2,4,6,7,8/)
    integer, parameter :: l65(1:33)=(/0,0,0,6,1,2,4,6,7,8,0,0,5,3,4,6,7,8,0,0,0,6,1,3,4,6,7,8,0,0,6,2,3/)
    integer, parameter :: l66(1:33)=(/4,6,7,8,0,0,7,1,2,3,4,6,7,8,0,4,5,6,7,8,0,0,0,0,5,1,5,6,7,8,0,0,0/)
    integer, parameter :: l67(1:33)=(/5,2,5,6,7,8,0,0,0,6,1,2,5,6,7,8,0,0,5,3,5,6,7,8,0,0,0,6,1,3,5,6,7/)
    integer, parameter :: l68(1:33)=(/8,0,0,6,2,3,5,6,7,8,0,0,7,1,2,3,5,6,7,8,0,5,4,5,6,7,8,0,0,0,6,1,4/)
    integer, parameter :: l69(1:33)=(/5,6,7,8,0,0,6,2,4,5,6,7,8,0,0,7,1,2,4,5,6,7,8,0,6,3,4,5,6,7,8,0,0/)
    integer, parameter :: l70(1:27)=(/7,1,3,4,5,6,7,8,0,7,2,3,4,5,6,7,8,0,8,1,2,3,4,5,6,7,8/)

    ! Some (rather nasty) data for the chunkwise decoding
    integer, parameter :: decode_map_arr(0:8,0:255) = reshape(&
        (/l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,&
        l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38,l39,l40,&
        l41,l42,l43,l44,l45,l46,l47,l48,l49,l50,l51,l52,l53,l54,l55,l56,l57,l58,l59,l60,&
        l61,l62,l63,l64,l65,l66,l67,l68,l69,l70/),(/9,256/) )

    private :: l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26
    private :: l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38,l39,l40,l41,l42,l43,l44,l45,l46,l47,l48,l49,l50
    private :: l51,l52,l53,l54,l55,l56,l57,l58,l59,l60,l61,l62,l63,l64,l65,l66,l67,l68,l69,l70

contains

    subroutine allocate_currentdets ()

        ! Allocate memory of the correct size for the currentdets array.

        integer :: ierr
        character(*), parameter :: this_routine = 'allocate_currentdets'

        allocate (WalkVecDets(0:NIfTot, MaxWalkersPart), stat=ierr)
        if (ierr /= 0) &
            call stop_all (this_routine, "Allocation failed for WalkVecDets")

    end subroutine allocate_currentdets

    subroutine init_bit_rep ()

        ! Set the values of nifd etc.

        character(*), parameter :: this_routine = 'init_bit_rep'

        ! This indicates the upper-bound for the determinants when expressed
        ! in bit-form. This will equal int(nBasis/32).
        ! The actual total length for a determinant in bit form will be
        ! NoIntforDet+1 + nIfY (which is the size of the Yamanouchi Symbol
        nIfD = int(nbasis / bits_n_int)

        ! Could use only 32-bits for this, except that it makes it very
        ! tricky to do do anything like sorting, as the latter 32-bits of the
        ! integer would contain random junk.
        NOffY = NIfD + 1
        if (tCSF) then
            if (tTruncateCSF) then
                NIfY = int(csf_trunc_level / bits_n_int) + 1
            else
                NIfY = int(nel / bits_n_int) + 1
            endif
        else
            NIfY = 0
        endif
        if (NIfY > 1) &
            call stop_all (this_routine, "CSFs with more than bits_n_int &
                          &open-shell electrons are not supported, and are &
                          &probably not a good idea.")

        ! The signs array
        NOffSgn = NOffY + NIfY
        NIfSgn = lenof_sign
#ifdef __PROG_NUMRUNS
        write(6,*) 'Calculation supports multiple parallel runs'
#elif defined(__DOUBLERUN)
        WRITE(6,*) "Double run in use."
#endif
#if defined(__CMPLX)
        WRITE(6,*) "Complex walkers in use."
#endif
        write(6,*) 'Number of simultaneous walker distributions: ',inum_runs
        write(6,*) 'Number of sign components in bit representation of determinant: ', NIfSgn

        ! The number of integers used for sorting / other bit manipulations
        NIfDBO = NIfD + NIfY

#ifdef __PROG_NUMRUNS
        if (lenof_sign_max /= 20) then
            call stop_all(this_routine, "Invalid build configuration. Update &
                         &flags to account for new lenof_sign_max, then &
                         &update this message")
        end if
#endif

! If we are using programattic lenofsign, then we also need to use separate
! integers for the flags, as the number of initiator/parent flags increases
! dramatically!

        ! K.G. 24.08.18
        ! Flags are being used in basically every calculation,
        ! considering recent developments, the possibility not to
        ! use flags is obsolete
        NIfFlag = 1

        NOffFlag = NOffSgn + NIfSgn

        ! N.B. Flags MUST be last!!!!!
        !      If we change this bit, then we need to adjust ilut_lt and
        !      ilut_gt.

        ! The total number of bits_n_int-bit integers used - 1
        NIfTot = NIfD + NIfY + NIfSgn + NIfFlag

        WRITE(6,"(A,I6)") "Setting integer length of determinants as bit-strings to: ", NIfTot + 1
        WRITE(6,"(A,I6)") "Setting integer bit-length of determinants as bit-strings to: ", bits_n_int

        ! By default we DO NOT initialise RDM parts of the bit rep now
        bit_rdm_init = .false.

        !
        ! The broadcasted information, used in annihilation, may require more
        ! information to be used.
        ! TODO: We may not always need the flags array. Test that...

        ! Create space for broadcasting the parent particle coefficient?
        ! ghb removed this ability on 14/4/16
        nOffParentCoeff = NIfTot + 1
        nIfParentCoeff = 0

        NIfBCast = NIfTot + nIfParentCoeff

        ! sometimes, we also need to store the number of spawn events
        ! in this iteration
        NSpawnOffset = NIfTot + 1
        if(tLogNumSpawns) then
            ! then, there is an extra integer in spawnedparts just behind
            ! the ilut noting the number of spawn events
            nOffParentCoeff = nOffParentCoeff + 1
            NIfBCast = NIfBCast + 1
        end if

        ! If we need to communicate the diagonal Hamiltonian element
        ! for the spawning
        if (tPreCond .or. tReplicaEstimates) then
            NOffSpawnHDiag = NIfBCast + 1
            NIfBCast = NIfBCast + 1
        end if

    end subroutine

    subroutine extract_bit_rep (ilut, nI, real_sgn, flags, j, store)

        ! Extract useful terms out of the bit-representation of a walker

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(out) :: nI(nel), flags
        integer, intent(in), optional :: j
        type(excit_gen_store_type), intent(inout), optional :: store
        real(dp), intent(out) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

        if(tStoredDets .and. present(j)) then
           nI = get_determinant(j)
        else if (tBuildOccVirtList .and. present(store)) then
            if(tBuildSpinSepLists) then
                call decode_bit_det_spinsep (nI, ilut, store)
            else
                call decode_bit_det_lists (nI, ilut, store)
            endif
        else
            call decode_bit_det (nI, ilut)
        endif

        sgn = iLut(NOffSgn:NOffSgn+lenof_sign-1)
        real_sgn = transfer(sgn, real_sgn)

        flags = int(iLut(NOffFlag), sizeof_int)

    end subroutine extract_bit_rep

    !Extract all flags as a single integer
    function extract_flags (iLut) result(flags)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer :: flags

        flags = int(ilut(NOffFlag), sizeof_int)

    end function extract_flags

    !Extract the sign (as a real_dp) for a particular element in the lenof_sign "array"
    pure function extract_part_sign (ilut, part_type) result(real_sgn)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: part_type
        real(dp) :: real_sgn
        real_sgn = transfer( ilut(nOffSgn + part_type - 1), real_sgn)
    end function

    pure function extract_run_sign(ilut, run) result(sgn)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: run
        HElement_t(dp) :: sgn

        ! Strange bug in compiler
        unused_var(run)
#ifdef __CMPLX
        sgn = cmplx(extract_part_sign(ilut, min_part_type(run)), extract_part_sign(ilut, max_part_type(run)))
#else
        sgn = extract_part_sign(ilut, min_part_type(run))
#endif
    end function


    !From the determinants, array of signs, and flag integer, create the
    !"packaged walker" representation
    pure subroutine encode_bit_rep (ilut, Det, real_sgn, flag)
        integer(n_int), intent(out) :: ilut(0:nIfTot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        integer(n_int), intent(in) :: Det(0:NIfDBO)
        integer, intent(in) :: flag
        integer(n_int) :: sgn(lenof_sign)

        iLut(0:NIfDBO) = Det

        sgn = transfer(real_sgn, sgn)
        iLut(NOffSgn:NOffSgn+NIfSgn-1) = sgn

        ilut(NOffFlag) = int(flag,n_int)

    end subroutine encode_bit_rep

    subroutine encode_flags (ilut, flag)

        ! Add new flag information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flag

        iLut(NOffFlag) = int(flag, n_int)

    end subroutine encode_flags

    pure function get_initiator_flag(sgn_index) result (flag)
        integer, intent(in) :: sgn_index
        integer :: flag
        ! Strange bug in compiler
        unused_var(sgn_index)
        ! map 1->1, 2->1, 3->3, 4->3, 5->5, 6->5 for complex,
        ! as the initiator flag is stored in the "real" bit
        ! of each run
        flag = flag_initiator(min_part_type(part_type_to_run(sgn_index)))
    end function get_initiator_flag

    pure function get_initiator_flag_by_run(run) result (flag)
        integer, intent(in) :: run
        integer :: flag
        ! Strange bug in compiler
        unused_var(run)
        ! map 1->1, 2->3, 3->5, 4->7 for complex
        flag = flag_initiator(min_part_type(run))
    end function get_initiator_flag_by_run

    pure function any_run_is_initiator(ilut) result (t)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer :: run
        logical :: t
        t = .false.
        do run = 1, inum_runs
            if (test_flag(ilut, get_initiator_flag_by_run(run))) then
                t = .true.
                return
            endif
        enddo
    end function any_run_is_initiator

    pure function all_runs_are_initiator(ilut) result(t)
      integer(n_int), intent(in) :: ilut(0:niftot)
      integer :: run
      logical :: t
      t = .true.
      do run = 1, inum_runs
        if(.not. test_flag(ilut, get_initiator_flag_by_run(run))) then
           t = .false.
           return
        endif
      end do
    end function all_runs_are_initiator

    subroutine clear_all_flags (ilut)

        ! Clear all of the flags

        integer(n_int), intent(inout) :: ilut(0:niftot)

        ilut(NOffFlag) = 0_n_int

    end subroutine clear_all_flags



    subroutine encode_sign (ilut, real_sgn)

        ! Add new sign information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

        sgn = transfer(real_sgn, sgn)
        iLut(NOffSgn:NOffSgn+NIfSgn-1) = sgn

    end subroutine encode_sign

    subroutine encode_run_sign (ilut, real_sgn, imag_sgn, run)

        ! Encode only the real AND imaginary component of the sign for the
        ! walker. Sign argument is now a scalar.
        !
        ! In:    real_sgn  - The new sign component
        !        imag_sgn  - The new imaginary sign component
        !        run - Update given run. 1 ==> inum_runs
        ! InOut:  ilut     - The bit representation to update
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(in) :: run
        real(dp), intent(in) :: real_sgn, imag_sgn
        character(*), parameter :: this_routine='encode_run_sign'

        ASSERT(run<=inum_runs)
        call encode_part_sign(ilut, real_sgn, min_part_type(run))
#ifdef __CMPLX
        call encode_part_sign(ilut, imag_sgn, max_part_type(run))
#else
        unused_var(imag_sgn)
#endif
    end subroutine encode_run_sign


    subroutine encode_part_sign (ilut, real_sgn, part_type)

        ! Encode only the real OR imaginary component of the sign for the
        ! walker. Sign argument is now a scalar.
        !
        ! In:    real_sgn  - The new sign component
        !        part_type - Update real/imaginary part. 1 ==> Re, 2 ==> Im
        ! InOut:  ilut     - The bit representation to update

        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(in) :: part_type
        real(dp), intent(in) :: real_sgn
        integer(n_int) :: sgn

        sgn = transfer(real_sgn, sgn)
        iLut(NOffSgn+part_type-1) = sgn

    end subroutine encode_part_sign

    subroutine nullify_ilut (ilut)

        ! Sets the sign of a determinant to equal zero.
        integer(n_int), intent(inout) :: ilut(0:NIfTot)

        iLut(NOffSgn:NOffSgn+NIfSgn-1) = transfer(0.0_dp, 0_n_int)

    end subroutine

    subroutine nullify_ilut_part (ilut, part_type)

        ! Sets the sign of the walkers of the specified particle type on
        ! a determinant to equal zero.
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(in) :: part_type

        iLut(NOffSgn+part_type-1) = transfer(0.0_dp, 0_n_int)

    end subroutine


    subroutine set_flag_general (ilut, flg, state)

        ! Set or clear the specified flag (0 indexed) according to
        ! the value in state.
        !
        ! In:    flg   - Integer index of flag to set
        !        state - Flag will be set if state is true.
        ! InOut: ilut  - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
        logical, intent(in) :: state

        if (state) then
            call set_flag_single (ilut, flg)
        else
            call clr_flag (ilut, flg)
        endif
    end subroutine set_flag_general

    subroutine set_flag_single (ilut, flg)

        ! Set the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to set
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
!        integer :: off, ind

!        ind = NOffFlag + flg / bits_n_int
!        off = mod(flg, bits_n_int)
!        ilut(ind) = ibset(ilut(ind), off)

        ! This now assumes that we do not have more flags than bits in an
        ! integer.
        ilut(NOffFlag) = ibset(ilut(NOffFlag), flg)

    end subroutine set_flag_single

    subroutine copy_flag (ilut_src, ilut_dest, flg)

        ! Copy the selected flag between iluts

        integer(n_int), intent(in) :: ilut_src(0:niftot)
        integer(n_int), intent(inout) :: ilut_dest(0:niftot)
        integer, intent(in) :: flg
        logical :: state

        state = test_flag (ilut_src, flg)
        call set_flag_general (ilut_dest, flg, state)

    end subroutine


    subroutine clr_flag (ilut, flg)

        ! Clear the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to clear
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
!        integer :: off, ind

!        ind = NOffFlag + flg / bits_n_int
!        off = mod(flg, bits_n_int)
!        ilut(ind) = ibclr(ilut(ind), off)

!This now assumes that we do not have more flags than bits in an integer.
        ilut(NOffFlag) = ibclr(ilut(NOffFlag), flg)

    end subroutine clr_flag

    function bit_parent_zero(ilut) result(zero)

        ! Used by the RDM functions
        ! Is the communicated parent zero?

        integer(n_int), intent(in) :: ilut(0:nIfBCast)
        logical :: zero
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'bit_parent_zero'
#endif

        ASSERT(bit_rdm_init)

        zero = all(ilut(NOffParent:NOffParent + NIfDBO) == 0)

    end function

    subroutine extract_parent(ilut, parent_ilut)

        integer(n_int), intent(in) :: ilut(0:nIfBCast)
        integer(n_int), intent(out) :: parent_ilut(0:NIfDBO)
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'extract_parent'
#endif

        ASSERT(bit_rdm_init)

        parent_ilut = ilut(nOffParent:nOffParent + NIfDBO)

    end subroutine

    subroutine encode_parent(ilut, ilut_parent, RDMBiasFacCurr)

        integer(n_int), intent(inout) :: ilut(0:NIfBCast)
        integer(n_int), intent(in) :: ilut_parent(0:NIfTot)
        real(dp), intent(in) :: RDMBiasFacCurr
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'encode_parent'
#endif

        ASSERT(bit_rdm_init)

        ilut(nOffParent:nOffParent + nIfDBO) = ilut_parent(0:NIfDBO)

        ilut(nOffParent + nIfDBO + 1) = &
            transfer(RDMBiasFacCurr, ilut(nOffParent + nIfDBO + 1))
        ! store the flag
        ilut(nOffParent + nIfDBO + 2) = ilut_parent(NIfTot)

    end subroutine

    subroutine zero_parent(ilut)

        integer(n_int), intent(inout) :: ilut(0:nIfBCast)
#ifdef __DEBUG
        character(*), parameter :: this_routine = 'zero_parent'
#endif

        ASSERT(bit_rdm_init)

        ilut(nOffParent:nOffParent+nIfDBO+1) = 0

    end subroutine

    subroutine set_parent_coeff(ilut, coeff)

        ! Store the coefficient of the parent walker of a spawn for more
        ! complex initiator logic. This option was removed by ghb on 14/4/16.
        ! so this routine should never be called

        integer(n_int), intent(inout) :: ilut(0:nIfBCast)
        real(dp), intent(in) :: coeff
        character(*), parameter :: this_routine = 'set_parent_coeff'

        call stop_all(this_routine,'Routine deprecated')

        ASSERT(nIfParentCoeff == 1)
        ilut(nOffParentCoeff) = transfer(coeff, ilut(nOffParentCoeff))

    end subroutine

    function extract_parent_coeff(ilut) result(coeff)

        ! Obtain the coefficient of the parent walker of a spawn for more
        ! complex initiator logic. This option was deprecated by ghb on 14/4/16.

        integer(n_int), intent(in) :: ilut(0:nIfBCast)
        real(dp) :: coeff
        character(*), parameter :: this_routine = 'extract_parent_coeff'

        call stop_all(this_routine,'Routine deprecated by ghb on 14/4/16')

        ASSERT(nIfParentCoeff == 1)

        coeff = transfer(ilut(nOffParentCoeff), coeff)

    end function

    subroutine encode_spawn_hdiag(ilut, hel)

        integer(n_int), intent(inout) :: ilut(0:NIfBCast)
        HElement_t(dp), intent(in) :: hel

        ilut(nOffSpawnHDiag) = transfer(hel, ilut(nOffSpawnHDiag))

    end subroutine encode_spawn_hdiag

    function extract_spawn_hdiag(ilut) result(hel)

        integer(n_int), intent(in) :: ilut(0:nIfBCast)

        HElement_t(dp) :: hel

        hel = transfer(ilut(nOffSpawnHDiag), hel)

    end function extract_spawn_hdiag

    subroutine log_spawn(ilut)

      ! set the spawn counter to 1
      implicit none
      integer(n_int), intent(inout) :: ilut(0:NIfBCast)

      ilut(NSpawnOffset) = 1
    end subroutine log_spawn

    subroutine increase_spawn_counter(ilut)
      ! increase the spawn counter by 1
      implicit none
      integer(n_int), intent(inout) :: ilut(0:NIfBCast)

      ilut(NSPawnOffset) = ilut(NSpawnOffset) + 1

    end subroutine increase_spawn_counter

    function get_num_spawns(ilut) result(nSpawn)
      ! read the number of spawns to this det so far
      implicit none
      integer(n_int), intent(inout) :: ilut(0:NIfBCast)
      integer :: nSpawn

      nSpawn = int(ilut(nSpawnOffset))

    end function get_num_spawns



    ! function test_flag is in bit_rep_data
    ! This avoids a circular dependence with DetBitOps.

    subroutine encode_det (ilut, Det)

        ! Add new det information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer(n_int), intent(in) :: Det(0:NIfDBO)

        iLut(0:NIfDBO) = Det

    end subroutine encode_det

    subroutine decode_bit_det_lists (nI, iLut, store)

        ! This routine decodes a determinant in bit form and constructs
        ! the natural ordered integer representation of the det.
        !
        ! It also constructs lists of the occupied and unoccupied orbitals
        ! within a symmetry.

        integer(n_int), intent(in) :: iLut(0:niftot)
        integer, intent(out) :: nI(:)
        type(excit_gen_store_type), intent(inout) :: store
        integer :: i, j, elec, orb, ind, virt(ScratchSize)
        integer :: nel_loc

        nel_loc = size(nI)

        ! Initialise the class counts
        store%ClassCountOcc = 0
        virt = 0

        elec = 0
        do i = 0, NIfD
            do j = 0, end_n_int
                orb = (i * bits_n_int) + (j + 1)
                ind = ClassCountInd(orb)
                if (btest(iLut(i), j)) then
                    !An electron is at this orbital
                    elec = elec + 1
                    nI(elec) = orb

                    ! Update class counts
                    store%ClassCountOcc(ind) = store%ClassCountOcc(ind) + 1

                    ! Store orbital INDEX in list of occ. orbs.
                    store%occ_list(store%ClassCountOcc(ind), ind) = elec

                    if (elec == nel_loc) exit
                else
                    ! Update count
                    virt(ind) = virt(ind) + 1

                    ! Store orbital in list of unocc. orbs.
                    store%virt_list(virt(ind), ind) = orb
                endif
            enddo
            if (elec == nel_loc) exit
        enddo

        ! Give final class count
        store%ClassCountUnocc = OrbClassCount - store%ClassCountOcc
        store%tFilled = .true.
        store%scratch3(1) = -1

        ! Fill in the remaineder of the virtuals list
        forall (ind = 1:ScratchSize)
            !if (virt(ind) /= store%ClassCountUnocc(ind)) then
                store%virt_list ( &
                    virt(ind) + 1 : &
                    store%ClassCountUnocc(ind), ind) = &
                SymLabelList2 (&
                    SymLabelCounts2(1, ind) + virt(ind) + &
                        store%ClassCountOcc(ind) : &
                    SymLabelCounts2(1, ind) + OrbClassCount(ind) - 1)
            !endif
        endforall

    end subroutine

    pure function getExcitationType(ExMat, IC) result(exTypeFlag)
        integer, intent(in) :: ExMat(2,2), IC
        integer :: exTypeFlag

        if (IC==1) then
            if (is_beta(ExMat(2,1)) .neqv. is_beta(ExMat(1,1))) then
                exTypeFlag = 3
                return
            else
                exTypeFlag = 1
            endif

        elseif (IC==2) then
            if (is_beta(ExMat(1,1)) .and. is_beta(ExMat(1,2))) then
                ! elec orbs are both beta
                if (is_beta(ExMat(2,1)) .and. is_beta(ExMat(2,2))) then
                    ! virt orbs are both beta
                    exTypeFlag = 2
                    return
                elseif (is_alpha(ExMat(2,1)) .and. is_alpha(ExMat(2,2))) then
                    ! virt orbs are both alpha
                    exTypeFlag = 5
                    return
                else
                    ! one of the spins changes
                    exTypeFlag = 4
                    return
                endif
            elseif (is_alpha(ExMat(1,1)) .and. is_alpha(ExMat(1,2))) then
                ! elec orbs are both alpha
                if (is_alpha(ExMat(2,1)) .and. is_alpha(ExMat(2,2))) then
                    ! virt orbs are both alpha
                    exTypeFlag = 2
                    return
                elseif (is_beta(ExMat(2,1)) .and. is_beta(ExMat(2,2))) then
                    ! virt orbs are both beta
                    exTypeFlag = 5
                    return
                else
                    ! one of the spins changes
                    exTypeFlag = 4
                    return
                endif
            else
                ! elec orb spins are different
                if (is_beta(ExMat(2,1)) .neqv. is_beta(ExMat(2,2))) then
                    ! virt orbs are of opposite spin
                    exTypeFlag = 2
                    return
                else
                    ! virt orbs are of the same spin
                    exTypeFlag = 4
                    return
                endif
            endif
        endif

    end function

    subroutine decode_bit_det_spinsep (nI, iLut, store)

        ! This routine decodes a determinant in bit form and constructs
        ! the natural ordered integer representation of the det.
        !
        ! It also constructs lists of the occupied and unoccupied orbitals
        ! within a symmetry.

        integer(n_int), intent(in) :: iLut(0:niftot)
        integer, intent(out) :: nI(:)
        type(excit_gen_store_type), intent(inout) :: store
        integer :: i, j, elec, orb, ind, virt(ScratchSize)
        integer :: nel_loc

        nel_loc = size(nI)

        ! Initialise the class counts
        store%ClassCountOcc = 0
        virt = 0

        elec = 0
        store%nel_alpha = 0

        do i = 0, NIfD
            do j = 0, end_n_int
                orb = (i * bits_n_int) + (j + 1)
                ind = ClassCountInd(orb)
                if (btest(iLut(i), j)) then
                    !An electron is at this orbital
                    elec = elec + 1
                    nI(elec) = orb

                    ! is the orbital spin alpha or beta?
                    if (mod(ind,2)==1) then
                        ! alpha
                        store%nel_alpha = store%nel_alpha+1
                        store%nI_alpha(store%nel_alpha) = orb
                        store%nI_alpha_inds(store%nel_alpha) = elec
                    else
                        store%nI_beta(elec-store%nel_alpha) = orb
                        store%nI_beta_inds(elec-store%nel_alpha) = elec
                    endif

                    ! Update class counts
                    store%ClassCountOcc(ind) = store%ClassCountOcc(ind) + 1

                    ! Store orbital INDEX in list of occ. orbs.
                    store%occ_list(store%ClassCountOcc(ind), ind) = elec

                    if (elec == nel_loc) exit
                else
                    ! Update count
                    virt(ind) = virt(ind) + 1
        !            write(iout,*) "filling virt"
                    ! Store orbital in list of unocc. orbs.
                    store%virt_list(virt(ind), ind) = orb
                endif
            enddo
            if (elec == nel_loc) exit
        enddo

        ! Give final class count
        store%ClassCountUnocc = OrbClassCount - store%ClassCountOcc
        store%tFilled = .true.
        store%scratch3(1) = -1

        ! Fill in the remainder of the virtuals list
        forall (ind = 1:ScratchSize)
            !if (virt(ind) /= store%ClassCountUnocc(ind)) then
                store%virt_list ( &
                    virt(ind) + 1 : &
                    store%ClassCountUnocc(ind), ind) = &
                SymLabelList2 (&
                    SymLabelCounts2(1, ind) + virt(ind) + &
                        store%ClassCountOcc(ind) : &
                    SymLabelCounts2(1, ind) + OrbClassCount(ind) - 1)
            !endif
        endforall

    end subroutine

    pure subroutine decode_bit_det_chunks (nI, iLut)

        ! This is a routine to take a determinant in bit form and construct
        ! the natural ordered integer form of the det.
        ! If CSFs are enabled, transfer the Yamanouchi symbol as well.

        use FciMCData, only: blank_det

        integer(n_int), intent(in) :: ilut(0:NIftot)
        integer, intent(out) :: nI(:)
        integer :: nopen, i, j, k, val, elec, offset, pos
        logical :: bIsCsf
        integer :: nel_loc

        nel_loc = size(nI)

        ! We need to use the CSF decoding routine if CSFs are enabled, and
        ! we are below a truncation limit if set.

        bIsCsf = .false.
        if (tCSF) then
            if (tTruncateCSF) then
                nopen = count_open_orbs(ilut)
                if (nopen <= csf_trunc_level) then
                    bIsCsf = .true.
                endif
            else
                bIsCsf = .true.
            endif
        endif

        elec = 0
        if (bIsCsf) then
            ! ****************
            ! Currently this just works in the old fashioned way. We aren't
            ! really that worried about CSF efficiency atm.
            ! ****************
            ! Consider the closed shell electrons first
            do i=0,NIfD
                do j=0,bits_n_int-2,2
                    if (btest(iLut(i),j)) then
                        if (btest(iLut(i),j+1)) then
                            ! An electron pair is in this spatial orbital
                            ! (2 matched spin orbitals)
                            elec = elec + 2
                            nI(elec-1) = (bits_n_int*i) + (j+1)
                            nI(elec) = (bits_n_int*i) + (j+2)
                            if (elec == nel_loc) return
                        endif
                    endif
                enddo
            enddo

            ! Now consider the open shell electrons
            ! TODO: can we move in steps of two, to catch unmatched pairs?
            nopen = 0
            do i=0,NIfD
                do j=0,end_n_int
                    if (btest(iLut(i),j)) then
                        if (.not.btest(iLut(i),ieor(j,1))) then
                            elec = elec + 1
                            nI(elec) = (bits_n_int*i) + (j+1)
                            pos = NIfD + 1 + (nopen/bits_n_int)
                            if (btest(iLut(Pos),mod(nopen,bits_n_int))) then
                                nI(elec) = ibset(nI(elec),csf_yama_bit)
                            endif
                            nopen = nopen + 1
                        endif
                    endif
                    if (elec==nel_loc) exit
                enddo
                if (elec==nel_loc) exit
            enddo
            ! If there are any open shell e-, set the csf bit
            nI = ibset(nI, csf_test_bit)
        else
            offset = 0
            do i = 0, NIfD
                do j = 0, bits_n_int - 1, 8
!                    val = iand(ishft(ilut(i), -j), Z'FF')
                    val = int(iand(ishft(ilut(i), -j), int(255,n_int)),sizeof_int)
                    do k = 1, decode_map_arr(0, val)
                        elec = elec + 1
                        nI(elec) = offset + decode_map_arr(k, val)
                        if (elec == nel_loc) return ! exit
                    enddo
                    offset = offset + 8
                enddo
            enddo

        endif

    end subroutine

    pure subroutine decode_bit_det_bitwise (nI, iLut)

        ! This is a routine to take a determinant in bit form and construct
        ! the natural ordered integer forim of the det.
        ! If CSFs are enabled, transfer the yamanouchi symbol as well.

        integer(n_int), intent(in) :: iLut(0:NIfTot)
        integer, intent(out) :: nI(:)
        integer :: i, j, elec, pos, nopen
        integer :: nel_loc
        logical :: bIsCsf

        nel_loc = size(nI)

        ! We need to use the CSF decoding routine if CSFs are enable, and we
        ! are below a truncation limit if set.
        bIsCsf = .false.
        if (tCSF) then
            if (tTruncateCSF) then
                nopen = count_open_orbs(iLut)
                if (nopen <= csf_trunc_level) then
                    bIsCsf = .true.
                endif
            else
                bIsCsf = .true.
            endif
        endif

        elec=0
        if (bIsCsf) then
            ! Consider the closed shell electrons first
            do i=0,NIfD
                do j=0,bits_n_int-2,2
                    if (btest(iLut(i),j)) then
                        if (btest(iLut(i),j+1)) then
                            ! An electron pair is in this spatial orbital
                            ! (2 matched spin orbitals)
                            elec = elec + 2
                            nI(elec-1) = (bits_n_int*i) + (j+1)
                            nI(elec) = (bits_n_int*i) + (j+2)
                            if (elec == nel_loc) return
                        endif
                    endif
                enddo
            enddo

            ! Now consider the open shell electrons
            ! TODO: can we move in steps of two, to catch unmatched pairs?
            nopen = 0
            do i=0,NIfD
                do j=0,end_n_int
                    if (btest(iLut(i),j)) then
                        if (.not.btest(iLut(i),ieor(j,1))) then
                            elec = elec + 1
                            nI(elec) = (bits_n_int*i) + (j+1)
                            pos = NIfD + 1 + (nopen/bits_n_int)
                            if (btest(iLut(Pos),mod(nopen,bits_n_int))) then
                                nI(elec) = ibset(nI(elec),csf_yama_bit)
                            endif
                            nopen = nopen + 1
                        endif
                    endif
                    if (elec==nel_loc) exit
                enddo
                if (elec==nel_loc) exit
            enddo
            ! If there are any open shell e-, set the csf bit
            nI = ibset(nI, csf_test_bit)
        else
            do i=0,NIfD
                do j=0,end_n_int
                    if(btest(iLut(i),j)) then
                        !An electron is at this orbital
                        elec=elec+1
                        nI(elec)=(i*bits_n_int)+(j+1)
                        if (elec == nel_loc) exit
                    endif
                enddo
                if (elec == nel_loc) exit
            enddo
        endif
    end subroutine decode_bit_det_bitwise

    subroutine add_ilut_lists(ndets_1, ndets_2, sorted_lists, list_1, list_2, list_out, ndets_out)

        ! WARNING 1: This routine assumes that both list_1 and list_2 contain no
        ! repeated iluts, even if one of the repeated iluts has zero amplitude.

        ! WARNING 2: If the input lists are not sorted (as defined by ilut_gt)
        ! then sorted_lists should be input as .false., and the lists will then
        ! be sorted. This routine will not work if unsorted lists are passed in
        ! and sorted_list is input as .true.

        use DetBitOps, only: ilut_lt, ilut_gt
        use sort_mod, only: sort
        use util_mod, only: binary_search_custom

        integer, intent(in) :: ndets_1
        integer, intent(in) :: ndets_2
        logical, intent(in) :: sorted_lists
        integer(n_int), intent(inout) :: list_1(0:,1:)
        integer(n_int), intent(inout) :: list_2(0:,1:)
        integer(n_int), intent(inout) :: list_out(0:,1:)
        integer, intent(out) :: ndets_out

        integer :: i, pos, min_ind
        real(dp) :: sign_1(lenof_sign), sign_2(lenof_sign), sign_out(lenof_sign)

        if (.not. sorted_lists) then
            call sort(list_1(:,1:ndets_1), ilut_lt, ilut_gt)
            call sort(list_2(:,1:ndets_2), ilut_lt, ilut_gt)
        end if

        ndets_out = 0
        ! Where to start searching from in list 1:
        min_ind = 1

        do i = 1, ndets_2
            ! If list_2(:,i) is in list 1 then pos will equal the position it
            ! occupies in list 1.
            ! If list_2(:,i) is not in list 1 then -pos will equal the position
            ! that it should go in, to mantain the sorted ordering.
            pos = binary_search_custom(list_1(:,min_ind:ndets_1), list_2(:,i), NIfTot+1, ilut_gt)

            if (pos > 0) then
                ! Move all the states from list 1 before min_ind+pos-1 across
                ! to the combined list.
                list_out(:,ndets_out+1:ndets_out+pos-1) = list_1(:,min_ind:min_ind+pos-2)
                ndets_out = ndets_out + pos - 1

                ndets_out = ndets_out + 1
                call extract_sign(list_1(:, min_ind+pos-1), sign_1)
                call extract_sign(list_2(:,i), sign_2)
                sign_out = sign_1 + sign_2
                list_out(:,ndets_out) = list_2(:,i)
                call encode_sign(list_out(:,ndets_out), sign_out)

                ! Search a smaller section of list_1 next time.
                min_ind = min_ind + pos
            else
                ! We have a state in list 2 which is not in list 1. Its
                ! position, if it were in list 1 would be min_ind-pos-1. Thus,
                ! first copy across all states from min_ind to min_ind-pos-2
                ! from list 1. Then copy across the state from list 2.
                list_out(:,ndets_out+1:ndets_out-pos-1) = list_1(:,min_ind:min_ind-pos-2)
                ndets_out = ndets_out - pos

                list_out(:,ndets_out) = list_2(:,i)

                ! Search a smaller section of list_1 next time.
                min_ind = min_ind - pos - 1
            end if
        end do

    end subroutine add_ilut_lists

!    subroutine init_excitations()
!        ! Allocate and initialise data in excit_mask.
!        use basis, only: bit_lookup, nbasis, basis_length
!        integer :: ibasis, jbasis, pos, el, ierr
!
!        allocate(excit_mask(basis_length, nbasis), stat=ierr)
!        excit_mask = 0
!
!        do ibasis = 1, nbasis
!            ! Set bits corresponding to all orbitals above ibasis.
!            ! Sure, there are quicker ways of doing this, but it's a one-off
!            do jbasis = ibasis+1, nbasis
!                pos = bit_lookup(1, jbasis)
!                el = bit_lookup(2, jbasis)
!                excit_mask(el, ibasis) = ibset(excit_mask(el, ibasis), pos)
!            end do
!        end do
!
!    end subroutine init_excitations
!
!    !This is a JSS routine for calculation a permutation from a double excitation
!    pure subroutine find_excitation_permutation2(f, excitation)
!        ! Find the parity of the permutation required to maximally line up
!        ! a determinant with an excitation of it, as needed for use with the
!        ! Slater--Condon rules.
!        !
!        ! This version is for double excitations of a determinant.
!        !
!        ! In:
!        !    f: bit string representation of the determinant.
!        !    excitation: excit type specifying how the excited determinant is
!        !        connected to the determinant described by f.
!        !        Note that we require the lists of orbitals excited from/into
!        !        to be ordered.
!        ! Out:
!        !    excitation: excit type with the parity of the permutation also
!        !        specified.
!
!        use basis, only: basis_length
!        use bit_utils, only: count_set_bits
!
!        integer(i0), intent(in) :: f(basis_length)
!        type(excit), intent(inout) :: excitation
!
!        integer :: perm
!        integer(i0) :: ia(basis_length), jb(basis_length)
!
!        ! Fast way of getting the parity of the permutation required to align
!        ! two determinants given one determinant and the connecting excitation.
!        ! This is hard to generalise to all cases, but we actually only care
!        ! about single and double excitations.  The idea is quite different from
!        ! that used in get_excitation (where we also need to find the orbitals
!        ! involved in the excitation).
!
!        ! In the following & represents the bitwise and operation; ^ the
!        ! bitwise exclusive or
!        ! operation; xmask is a mask with all bits representing orbitals above
!        ! x set; f is the string representing the determinant from which we
!        ! excite and the excitation is defined by (i,j)->(a,b), where i<j and
!        ! a<b.
!
!        ! imask ^ amask returns a bit string with bits corresponding to all
!        ! orbitals between i and a set, with max(i,a) set and min(i,a) cleared.
!        ! Thus f & (imask ^ amask) returns a bit string with only bits set for
!        ! the occupied orbitals which are between i and a (possibly including i)
!        ! and so the popcount of this gives the number of orbitals between i and
!        ! a (possibly one larger than the actual answer) number of permutations
!        ! needed to
!        ! align i and a in the same 'slot' in the determinant string.  We need
!        ! to subtract one if i>a to correct for the overcounting.
!
!        ! An analagous approach counts the number of permutations required so
!        ! j and b are coincident.
!
!        ! Finally, we need to account for some more overcounting/undercounting.
!        ! If j is between i and a, then it is counted yet j can either be moved
!        ! before i (resulting in the actual number of permutations being one
!        ! less than that counted) or after i (resulting in moving j taking one
!        ! more permutation than counted).  It doesn't matter which we do, as we
!        ! are only interested in whether the number of permutations is odd or
!        ! even.  We similarly need to take into account the case where i is
!        ! between j and b.
!
!        ia =
!ieor(excit_mask(:,excitation%from_orb(1)),excit_mask(:,excitation%to_orb(1)))
!        jb =
!ieor(excit_mask(:,excitation%from_orb(2)),excit_mask(:,excitation%to_orb(2)))
!
!        perm = sum(count_set_bits(iand(f,ia))) +
!sum(count_set_bits(iand(f,jb)))
!
!        if (excitation%from_orb(1) > excitation%to_orb(1)) perm = perm - 1
!        if (excitation%from_orb(1) > excitation%to_orb(2)) perm = perm - 1
!        if (excitation%from_orb(2) > excitation%to_orb(2) .or. &
!            excitation%from_orb(2) < excitation%to_orb(1)) perm = perm - 1
!
!        excitation%perm = mod(perm,2) == 1
!
!    end subroutine find_excitation_permutation2


end module bit_reps
