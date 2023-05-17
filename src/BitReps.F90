#include "macros.h"

module bit_reps
    use FciMCData, only: WalkVecDets, MaxWalkersPart, tLogNumSpawns, blank_det

    use SystemData, only: nel, nbasis, tGUGA

    use CalcData, only: tTruncInitiator, tUseRealCoeffs, tSemiStochastic, &
                        tTrialWavefunction, semistoch_shift_iter, &
                        tStartTrialLater, tPreCond, tReplicaEstimates, tStoredDets

    use constants, only: lenof_sign, end_n_int, bits_n_int, n_int, dp, sizeof_int, stdout, &
        inum_runs_max, inum_runs

    use DetBitOps, only: count_open_orbs, CountBits

    use DetBitOps, only: ilut_lt, ilut_gt

    use bit_rep_data, only: test_flag, flag_initiator, niftot, nIfD, IlutBits, &
        extract_sign, bit_rdm_init, nIfGUGA, IlutBitsParent

    use SymExcitDataMod, only: excit_gen_store_type, tBuildOccVirtList, &
                               tBuildSpinSepLists, &
                               OrbClassCount, ScratchSize, SymLabelList2, &
                               SymLabelCounts2

    use sym_general_mod, only: ClassCountInd

    use util_mod, only: binary_search_custom

    use sort_mod, only: sort

    use global_det_data, only: get_determinant

    use guga_bitRepOps, only: init_guga_bitrep

    use LoggingData, only: tRDMOnfly

    use guga_bitRepOps, only: transfer_stochastic_rdm_info

    use DeterminantData, only: write_det

    use util_mod, only: stop_all, neci_flush

    better_implicit_none
    private
    public :: decode_bit_det, init_bit_rep, nullify_ilut, encode_sign, get_initiator_flag, &
        decode_bit_det_lists, writebitdet, getExcitationType, &
        decode_bit_det_spinsep, set_flag, extract_bit_rep, any_run_is_initiator, &
        all_runs_are_initiator, nullify_ilut_part, encode_part_sign, &
        zero_parent, get_initiator_flag_by_run, encode_parent, clr_flag_multi, &
        clr_flag, extract_flags, encode_bit_rep, extract_part_sign, log_spawn, &
        increase_spawn_counter, encode_spawn_hdiag, extract_spawn_hdiag, &
        get_num_spawns, bit_parent_zero, encode_det, add_ilut_lists, &
        clear_all_flags, extract_run_sign, encode_flags



    ! Structure of a bit representation:

    ! | 0-NIfD: Det | Sign(Re) | Sign(Im) | Flags |
    !
    ! -------
    ! (NIfD + 1) * 64-bits              Orbital rep.
    !  1         * 32-bits              Signs (Re)
    ! (1         * 32-bits if needed)   Signs (Im)
    ! (1         * 32-bits if needed)   Flags

    interface set_flag
        module procedure set_flag_single
        module procedure set_flag_general
    end interface

    ! Which decoding function do we want to use?
    interface decode_bit_det
        module procedure decode_bit_det_chunks
        module procedure decode_bit_det_lists
    end interface

    integer, parameter :: l1(1:33) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 0, 0, 0]
    integer, parameter :: l2(1:33) = [0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 2, 1, 3, 0, 0, 0, 0, 0, 0, 2, 2, 3, 0, 0, 0, 0, 0, 0, 3, 1, 2]
    integer, parameter :: l3(1:33) = [3, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 2, 1, 4, 0, 0, 0, 0, 0, 0, 2, 2, 4, 0, 0, 0, 0, 0, 0]
    integer, parameter :: l4(1:33) = [3, 1, 2, 4, 0, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0, 0, 3, 1, 3, 4, 0, 0, 0, 0, 0, 3, 2, 3, 4, 0, 0]
    integer, parameter :: l5(1:33) = [0, 0, 0, 4, 1, 2, 3, 4, 0, 0, 0, 0, 1, 5, 0, 0, 0, 0, 0, 0, 0, 2, 1, 5, 0, 0, 0, 0, 0, 0, 2, 2, 5]
    integer, parameter :: l6(1:33) = [0, 0, 0, 0, 0, 0, 3, 1, 2, 5, 0, 0, 0, 0, 0, 2, 3, 5, 0, 0, 0, 0, 0, 0, 3, 1, 3, 5, 0, 0, 0, 0, 0]
    integer, parameter :: l7(1:33) = [3, 2, 3, 5, 0, 0, 0, 0, 0, 4, 1, 2, 3, 5, 0, 0, 0, 0, 2, 4, 5, 0, 0, 0, 0, 0, 0, 3, 1, 4, 5, 0, 0]
    integer, parameter :: l8(1:33) = [0, 0, 0, 3, 2, 4, 5, 0, 0, 0, 0, 0, 4, 1, 2, 4, 5, 0, 0, 0, 0, 3, 3, 4, 5, 0, 0, 0, 0, 0, 4, 1, 3]
    integer, parameter :: l9(1:33) = [4, 5, 0, 0, 0, 0, 4, 2, 3, 4, 5, 0, 0, 0, 0, 5, 1, 2, 3, 4, 5, 0, 0, 0, 1, 6, 0, 0, 0, 0, 0, 0, 0]
    integer, parameter :: l10(1:33) = [2, 1, 6, 0, 0, 0, 0, 0, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 3, 1, 2, 6, 0, 0, 0, 0, 0, 2, 3, 6, 0, 0, 0]
    integer, parameter :: l11(1:33) = [0, 0, 0, 3, 1, 3, 6, 0, 0, 0, 0, 0, 3, 2, 3, 6, 0, 0, 0, 0, 0, 4, 1, 2, 3, 6, 0, 0, 0, 0, 2, 4, 6]
    integer, parameter :: l12(1:33) = [0, 0, 0, 0, 0, 0, 3, 1, 4, 6, 0, 0, 0, 0, 0, 3, 2, 4, 6, 0, 0, 0, 0, 0, 4, 1, 2, 4, 6, 0, 0, 0, 0]
    integer, parameter :: l13(1:33) = [3, 3, 4, 6, 0, 0, 0, 0, 0, 4, 1, 3, 4, 6, 0, 0, 0, 0, 4, 2, 3, 4, 6, 0, 0, 0, 0, 5, 1, 2, 3, 4, 6]
    integer, parameter :: l14(1:33) = [0, 0, 0, 2, 5, 6, 0, 0, 0, 0, 0, 0, 3, 1, 5, 6, 0, 0, 0, 0, 0, 3, 2, 5, 6, 0, 0, 0, 0, 0, 4, 1, 2]
    integer, parameter :: l15(1:33) = [5, 6, 0, 0, 0, 0, 3, 3, 5, 6, 0, 0, 0, 0, 0, 4, 1, 3, 5, 6, 0, 0, 0, 0, 4, 2, 3, 5, 6, 0, 0, 0, 0]
    integer, parameter :: l16(1:33) = [5, 1, 2, 3, 5, 6, 0, 0, 0, 3, 4, 5, 6, 0, 0, 0, 0, 0, 4, 1, 4, 5, 6, 0, 0, 0, 0, 4, 2, 4, 5, 6, 0]
    integer, parameter :: l17(1:33) = [0, 0, 0, 5, 1, 2, 4, 5, 6, 0, 0, 0, 4, 3, 4, 5, 6, 0, 0, 0, 0, 5, 1, 3, 4, 5, 6, 0, 0, 0, 5, 2, 3]
    integer, parameter :: l18(1:33) = [4, 5, 6, 0, 0, 0, 6, 1, 2, 3, 4, 5, 6, 0, 0, 1, 7, 0, 0, 0, 0, 0, 0, 0, 2, 1, 7, 0, 0, 0, 0, 0, 0]
    integer, parameter :: l19(1:33) = [2, 2, 7, 0, 0, 0, 0, 0, 0, 3, 1, 2, 7, 0, 0, 0, 0, 0, 2, 3, 7, 0, 0, 0, 0, 0, 0, 3, 1, 3, 7, 0, 0]
    integer, parameter :: l20(1:33) = [0, 0, 0, 3, 2, 3, 7, 0, 0, 0, 0, 0, 4, 1, 2, 3, 7, 0, 0, 0, 0, 2, 4, 7, 0, 0, 0, 0, 0, 0, 3, 1, 4]
    integer, parameter :: l21(1:33) = [7, 0, 0, 0, 0, 0, 3, 2, 4, 7, 0, 0, 0, 0, 0, 4, 1, 2, 4, 7, 0, 0, 0, 0, 3, 3, 4, 7, 0, 0, 0, 0, 0]
    integer, parameter :: l22(1:33) = [4, 1, 3, 4, 7, 0, 0, 0, 0, 4, 2, 3, 4, 7, 0, 0, 0, 0, 5, 1, 2, 3, 4, 7, 0, 0, 0, 2, 5, 7, 0, 0, 0]
    integer, parameter :: l23(1:33) = [0, 0, 0, 3, 1, 5, 7, 0, 0, 0, 0, 0, 3, 2, 5, 7, 0, 0, 0, 0, 0, 4, 1, 2, 5, 7, 0, 0, 0, 0, 3, 3, 5]
    integer, parameter :: l24(1:33) = [7, 0, 0, 0, 0, 0, 4, 1, 3, 5, 7, 0, 0, 0, 0, 4, 2, 3, 5, 7, 0, 0, 0, 0, 5, 1, 2, 3, 5, 7, 0, 0, 0]
    integer, parameter :: l25(1:33) = [3, 4, 5, 7, 0, 0, 0, 0, 0, 4, 1, 4, 5, 7, 0, 0, 0, 0, 4, 2, 4, 5, 7, 0, 0, 0, 0, 5, 1, 2, 4, 5, 7]
    integer, parameter :: l26(1:33) = [0, 0, 0, 4, 3, 4, 5, 7, 0, 0, 0, 0, 5, 1, 3, 4, 5, 7, 0, 0, 0, 5, 2, 3, 4, 5, 7, 0, 0, 0, 6, 1, 2]
    integer, parameter :: l27(1:33) = [3, 4, 5, 7, 0, 0, 2, 6, 7, 0, 0, 0, 0, 0, 0, 3, 1, 6, 7, 0, 0, 0, 0, 0, 3, 2, 6, 7, 0, 0, 0, 0, 0]
    integer, parameter :: l28(1:33) = [4, 1, 2, 6, 7, 0, 0, 0, 0, 3, 3, 6, 7, 0, 0, 0, 0, 0, 4, 1, 3, 6, 7, 0, 0, 0, 0, 4, 2, 3, 6, 7, 0]
    integer, parameter :: l29(1:33) = [0, 0, 0, 5, 1, 2, 3, 6, 7, 0, 0, 0, 3, 4, 6, 7, 0, 0, 0, 0, 0, 4, 1, 4, 6, 7, 0, 0, 0, 0, 4, 2, 4]
    integer, parameter :: l30(1:33) = [6, 7, 0, 0, 0, 0, 5, 1, 2, 4, 6, 7, 0, 0, 0, 4, 3, 4, 6, 7, 0, 0, 0, 0, 5, 1, 3, 4, 6, 7, 0, 0, 0]
    integer, parameter :: l31(1:33) = [5, 2, 3, 4, 6, 7, 0, 0, 0, 6, 1, 2, 3, 4, 6, 7, 0, 0, 3, 5, 6, 7, 0, 0, 0, 0, 0, 4, 1, 5, 6, 7, 0]
    integer, parameter :: l32(1:33) = [0, 0, 0, 4, 2, 5, 6, 7, 0, 0, 0, 0, 5, 1, 2, 5, 6, 7, 0, 0, 0, 4, 3, 5, 6, 7, 0, 0, 0, 0, 5, 1, 3]
    integer, parameter :: l33(1:33) = [5, 6, 7, 0, 0, 0, 5, 2, 3, 5, 6, 7, 0, 0, 0, 6, 1, 2, 3, 5, 6, 7, 0, 0, 4, 4, 5, 6, 7, 0, 0, 0, 0]
    integer, parameter :: l34(1:33) = [5, 1, 4, 5, 6, 7, 0, 0, 0, 5, 2, 4, 5, 6, 7, 0, 0, 0, 6, 1, 2, 4, 5, 6, 7, 0, 0, 5, 3, 4, 5, 6, 7]
    integer, parameter :: l35(1:33) = [0, 0, 0, 6, 1, 3, 4, 5, 6, 7, 0, 0, 6, 2, 3, 4, 5, 6, 7, 0, 0, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 8, 0]
    integer, parameter :: l36(1:33) = [0, 0, 0, 0, 0, 0, 2, 1, 8, 0, 0, 0, 0, 0, 0, 2, 2, 8, 0, 0, 0, 0, 0, 0, 3, 1, 2, 8, 0, 0, 0, 0, 0]
    integer, parameter :: l37(1:33) = [2, 3, 8, 0, 0, 0, 0, 0, 0, 3, 1, 3, 8, 0, 0, 0, 0, 0, 3, 2, 3, 8, 0, 0, 0, 0, 0, 4, 1, 2, 3, 8, 0]
    integer, parameter :: l38(1:33) = [0, 0, 0, 2, 4, 8, 0, 0, 0, 0, 0, 0, 3, 1, 4, 8, 0, 0, 0, 0, 0, 3, 2, 4, 8, 0, 0, 0, 0, 0, 4, 1, 2]
    integer, parameter :: l39(1:33) = [4, 8, 0, 0, 0, 0, 3, 3, 4, 8, 0, 0, 0, 0, 0, 4, 1, 3, 4, 8, 0, 0, 0, 0, 4, 2, 3, 4, 8, 0, 0, 0, 0]
    integer, parameter :: l40(1:33) = [5, 1, 2, 3, 4, 8, 0, 0, 0, 2, 5, 8, 0, 0, 0, 0, 0, 0, 3, 1, 5, 8, 0, 0, 0, 0, 0, 3, 2, 5, 8, 0, 0]
    integer, parameter :: l41(1:33) = [0, 0, 0, 4, 1, 2, 5, 8, 0, 0, 0, 0, 3, 3, 5, 8, 0, 0, 0, 0, 0, 4, 1, 3, 5, 8, 0, 0, 0, 0, 4, 2, 3]
    integer, parameter :: l42(1:33) = [5, 8, 0, 0, 0, 0, 5, 1, 2, 3, 5, 8, 0, 0, 0, 3, 4, 5, 8, 0, 0, 0, 0, 0, 4, 1, 4, 5, 8, 0, 0, 0, 0]
    integer, parameter :: l43(1:33) = [4, 2, 4, 5, 8, 0, 0, 0, 0, 5, 1, 2, 4, 5, 8, 0, 0, 0, 4, 3, 4, 5, 8, 0, 0, 0, 0, 5, 1, 3, 4, 5, 8]
    integer, parameter :: l44(1:33) = [0, 0, 0, 5, 2, 3, 4, 5, 8, 0, 0, 0, 6, 1, 2, 3, 4, 5, 8, 0, 0, 2, 6, 8, 0, 0, 0, 0, 0, 0, 3, 1, 6]
    integer, parameter :: l45(1:33) = [8, 0, 0, 0, 0, 0, 3, 2, 6, 8, 0, 0, 0, 0, 0, 4, 1, 2, 6, 8, 0, 0, 0, 0, 3, 3, 6, 8, 0, 0, 0, 0, 0]
    integer, parameter :: l46(1:33) = [4, 1, 3, 6, 8, 0, 0, 0, 0, 4, 2, 3, 6, 8, 0, 0, 0, 0, 5, 1, 2, 3, 6, 8, 0, 0, 0, 3, 4, 6, 8, 0, 0]
    integer, parameter :: l47(1:33) = [0, 0, 0, 4, 1, 4, 6, 8, 0, 0, 0, 0, 4, 2, 4, 6, 8, 0, 0, 0, 0, 5, 1, 2, 4, 6, 8, 0, 0, 0, 4, 3, 4]
    integer, parameter :: l48(1:33) = [6, 8, 0, 0, 0, 0, 5, 1, 3, 4, 6, 8, 0, 0, 0, 5, 2, 3, 4, 6, 8, 0, 0, 0, 6, 1, 2, 3, 4, 6, 8, 0, 0]
    integer, parameter :: l49(1:33) = [3, 5, 6, 8, 0, 0, 0, 0, 0, 4, 1, 5, 6, 8, 0, 0, 0, 0, 4, 2, 5, 6, 8, 0, 0, 0, 0, 5, 1, 2, 5, 6, 8]
    integer, parameter :: l50(1:33) = [0, 0, 0, 4, 3, 5, 6, 8, 0, 0, 0, 0, 5, 1, 3, 5, 6, 8, 0, 0, 0, 5, 2, 3, 5, 6, 8, 0, 0, 0, 6, 1, 2]
    integer, parameter :: l51(1:33) = [3, 5, 6, 8, 0, 0, 4, 4, 5, 6, 8, 0, 0, 0, 0, 5, 1, 4, 5, 6, 8, 0, 0, 0, 5, 2, 4, 5, 6, 8, 0, 0, 0]
    integer, parameter :: l52(1:33) = [6, 1, 2, 4, 5, 6, 8, 0, 0, 5, 3, 4, 5, 6, 8, 0, 0, 0, 6, 1, 3, 4, 5, 6, 8, 0, 0, 6, 2, 3, 4, 5, 6]
    integer, parameter :: l53(1:33) = [8, 0, 0, 7, 1, 2, 3, 4, 5, 6, 8, 0, 2, 7, 8, 0, 0, 0, 0, 0, 0, 3, 1, 7, 8, 0, 0, 0, 0, 0, 3, 2, 7]
    integer, parameter :: l54(1:33) = [8, 0, 0, 0, 0, 0, 4, 1, 2, 7, 8, 0, 0, 0, 0, 3, 3, 7, 8, 0, 0, 0, 0, 0, 4, 1, 3, 7, 8, 0, 0, 0, 0]
    integer, parameter :: l55(1:33) = [4, 2, 3, 7, 8, 0, 0, 0, 0, 5, 1, 2, 3, 7, 8, 0, 0, 0, 3, 4, 7, 8, 0, 0, 0, 0, 0, 4, 1, 4, 7, 8, 0]
    integer, parameter :: l56(1:33) = [0, 0, 0, 4, 2, 4, 7, 8, 0, 0, 0, 0, 5, 1, 2, 4, 7, 8, 0, 0, 0, 4, 3, 4, 7, 8, 0, 0, 0, 0, 5, 1, 3]
    integer, parameter :: l57(1:33) = [4, 7, 8, 0, 0, 0, 5, 2, 3, 4, 7, 8, 0, 0, 0, 6, 1, 2, 3, 4, 7, 8, 0, 0, 3, 5, 7, 8, 0, 0, 0, 0, 0]
    integer, parameter :: l58(1:33) = [4, 1, 5, 7, 8, 0, 0, 0, 0, 4, 2, 5, 7, 8, 0, 0, 0, 0, 5, 1, 2, 5, 7, 8, 0, 0, 0, 4, 3, 5, 7, 8, 0]
    integer, parameter :: l59(1:33) = [0, 0, 0, 5, 1, 3, 5, 7, 8, 0, 0, 0, 5, 2, 3, 5, 7, 8, 0, 0, 0, 6, 1, 2, 3, 5, 7, 8, 0, 0, 4, 4, 5]
    integer, parameter :: l60(1:33) = [7, 8, 0, 0, 0, 0, 5, 1, 4, 5, 7, 8, 0, 0, 0, 5, 2, 4, 5, 7, 8, 0, 0, 0, 6, 1, 2, 4, 5, 7, 8, 0, 0]
    integer, parameter :: l61(1:33) = [5, 3, 4, 5, 7, 8, 0, 0, 0, 6, 1, 3, 4, 5, 7, 8, 0, 0, 6, 2, 3, 4, 5, 7, 8, 0, 0, 7, 1, 2, 3, 4, 5]
    integer, parameter :: l62(1:33) = [7, 8, 0, 3, 6, 7, 8, 0, 0, 0, 0, 0, 4, 1, 6, 7, 8, 0, 0, 0, 0, 4, 2, 6, 7, 8, 0, 0, 0, 0, 5, 1, 2]
    integer, parameter :: l63(1:33) = [6, 7, 8, 0, 0, 0, 4, 3, 6, 7, 8, 0, 0, 0, 0, 5, 1, 3, 6, 7, 8, 0, 0, 0, 5, 2, 3, 6, 7, 8, 0, 0, 0]
    integer, parameter :: l64(1:33) = [6, 1, 2, 3, 6, 7, 8, 0, 0, 4, 4, 6, 7, 8, 0, 0, 0, 0, 5, 1, 4, 6, 7, 8, 0, 0, 0, 5, 2, 4, 6, 7, 8]
    integer, parameter :: l65(1:33) = [0, 0, 0, 6, 1, 2, 4, 6, 7, 8, 0, 0, 5, 3, 4, 6, 7, 8, 0, 0, 0, 6, 1, 3, 4, 6, 7, 8, 0, 0, 6, 2, 3]
    integer, parameter :: l66(1:33) = [4, 6, 7, 8, 0, 0, 7, 1, 2, 3, 4, 6, 7, 8, 0, 4, 5, 6, 7, 8, 0, 0, 0, 0, 5, 1, 5, 6, 7, 8, 0, 0, 0]
    integer, parameter :: l67(1:33) = [5, 2, 5, 6, 7, 8, 0, 0, 0, 6, 1, 2, 5, 6, 7, 8, 0, 0, 5, 3, 5, 6, 7, 8, 0, 0, 0, 6, 1, 3, 5, 6, 7]
    integer, parameter :: l68(1:33) = [8, 0, 0, 6, 2, 3, 5, 6, 7, 8, 0, 0, 7, 1, 2, 3, 5, 6, 7, 8, 0, 5, 4, 5, 6, 7, 8, 0, 0, 0, 6, 1, 4]
    integer, parameter :: l69(1:33) = [5, 6, 7, 8, 0, 0, 6, 2, 4, 5, 6, 7, 8, 0, 0, 7, 1, 2, 4, 5, 6, 7, 8, 0, 6, 3, 4, 5, 6, 7, 8, 0, 0]
    integer, parameter :: l70(1:27) = [7, 1, 3, 4, 5, 6, 7, 8, 0, 7, 2, 3, 4, 5, 6, 7, 8, 0, 8, 1, 2, 3, 4, 5, 6, 7, 8]

    ! Some (rather nasty) data for the chunkwise decoding
    integer, parameter :: decode_map_arr(0:8, 0:255) = reshape( &
                          [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, &
                            l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31, l32, l33, l34, l35, l36, l37, l38, l39, l40, &
                            l41, l42, l43, l44, l45, l46, l47, l48, l49, l50, l51, l52, l53, l54, l55, l56, l57, l58, l59, l60, &
                            l61, l62, l63, l64, l65, l66, l67, l68, l69, l70], [9, 256])

    private :: l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26
    private :: l27, l28, l29, l30, l31, l32, l33, l34, l35, l36, l37, l38, l39, l40, l41, l42, l43, l44, l45, l46, l47, l48, l49, l50
    private :: l51, l52, l53, l54, l55, l56, l57, l58, l59, l60, l61, l62, l63, l64, l65, l66, l67, l68, l69, l70

contains

    subroutine init_bit_rep()

        ! Set the values of nifd etc.
#ifdef PROG_NUMRUNS_
        character(*), parameter :: this_routine = 'init_bit_rep'
#endif

        ! This indicates the upper-bound for the determinants when expressed
        ! in bit-form. This will equal int(nBasis/32).
        ! The actual total length for a determinant in bit form will be
        nIfD = int(nbasis / bits_n_int)
        IlutBits%len_orb = nifd

        ! The signs array
        IlutBits%ind_pop = IlutBits%len_orb + 1
        IlutBits%len_pop = lenof_sign
#ifdef PROG_NUMRUNS_
        write(stdout, *) 'Calculation supports multiple parallel runs'
#elif defined(DOUBLERUN_)
        write(stdout, *) "Double run in use."
#endif
#if defined(CMPLX_)
        write(stdout, *) "Complex walkers in use."
#endif
        write(stdout, *) 'Number of simultaneous walker distributions: ', inum_runs
        write(stdout, *) 'Number of sign components in bit representation of determinant: ', &
            IlutBits%len_pop

        ! The number of integers used for sorting / other bit manipulations
        ! WD: this is always just nifd.. so remove nifdbo..

#ifdef PROG_NUMRUNS_
        if (inum_runs > inum_runs_max) then
            write(stdout, *) "Maximally", inum_runs_max, "replicas are allowed"
            call stop_all(this_routine, "Requesting more than the maximum number of replicas")
        end if
#endif

! If we are using programattic lenofsign, then we also need to use separate
! integers for the flags, as the number of initiator/parent flags increases
! dramatically!

        IlutBits%ind_flag = IlutBits%ind_pop + IlutBits%len_pop

        ! N.B. Flags MUST be last!!!!!
        !      If we change this bit, then we need to adjust ilut_lt and
        !      ilut_gt.

        ! The total number of bits_n_int-bit integers used - 1
        ! WD: the +1 is for the always used flag entry now
        NIfTot = IlutBits%len_orb + IlutBits%len_pop + 1
        IlutBits%len_tot = IlutBits%len_orb + IlutBits%len_pop + 1

        write(stdout, "(A,I6)") "Setting integer length of determinants as bit-strings to: ", &
            IlutBits%len_tot + 1
        write(stdout, "(A,I6)") "Setting integer bit-length of determinants as bit-strings to: ", &
            bits_n_int

        if (tGUGA) then
            ! set up a nIfGUGA variable to use a similar integer list to
            ! calculate excitations for a given GUGA CSF

            ! Structure of a bit representation:
            ! the parenthesis is for the stochastic GUGA rdm implementation
            ! | 0-NIfD: Det | x0 | x1 | deltaB | (rdm_ind | rdm_x0 | rdm_x1)
            !
            ! -------
            ! (NIfD + 1) * 64-bits              Orbital rep.
            !  1         * 64-bits              x0 matrix element
            !  1         * 64-bits              x1 matrix element
            !  1         * 32-bits              deltaB value

            call init_guga_bitrep(nifd)
            write(stdout, "(A,I6)") "For GUGA calculation set up a integer list of length: ", &
                nIfGUGA + 1

            ! if we use fast guga rdms we also need space in the
            ! 'normal' ilut to store the rdm_ind and x0 and x1..
            ! for communication or? yes.. and also the parent stuff below
            ! need to be adapted..
            ! but this gets changed in rdm_general.. so do the GUGA
            ! changes there!
            ! no.. I also need to change niftot to be able to hold
            ! rdm_ind, x0 and x1 from the excitation generation step..
            ! so I also need to adapt this here!
            ! and I think I just need to store it within niftot!
            ! I do not even need and additional entry in the parent
            ! atleast in the communication within spawnedparts!
            if (tRDMOnfly) then
                IlutBits%ind_rdm_ind = niftot + 1
                IlutBits%ind_rdm_x0 = IlutBits%ind_rdm_ind + 1
                IlutBits%ind_rdm_x1 = IlutBits%ind_rdm_x0 + 1

                niftot = IlutBits%ind_rdm_x1
                IlutBits%len_tot = IlutBits%ind_rdm_x1
            end if
        end if

        ! By default we DO NOT initialise RDM parts of the bit rep now
        bit_rdm_init = .false.

        !
        ! The broadcasted information, used in annihilation, may require more
        ! information to be used.
        ! TODO: We may not always need the flags array. Test that...
        IlutBits%len_bcast = IlutBits%len_tot

        ! sometimes, we also need to store the number of spawn events
        ! in this iteration
        IlutBits%ind_spawn = IlutBits%len_tot + 1

        if (tLogNumSpawns) then
            ! then, there is an extra integer in spawnedparts just behind
            ! the ilut noting the number of spawn events
            IlutBits%len_bcast = IlutBits%len_bcast + 1
        end if

        ! If we need to communicate the diagonal Hamiltonian element
        ! for the spawning
        if (tPreCond .or. tReplicaEstimates) then
            IlutBits%ind_hdiag = IlutBits%len_bcast + 1
            IlutBits%len_bcast = IlutBits%len_bcast + 1
        end if

        ! also store the information for the spawned_parents in the
        ! RDM calculation in this data-stucture! for a nicer overview!
        IlutBitsParent%len_orb = IlutBits%len_orb
        IlutBitsParent%ind_pop = IlutBitsParent%len_orb + 1
        IlutBitsParent%ind_flag = IlutBitsParent%ind_pop + 1
        IlutBitsParent%ind_source = IlutBitsParent%ind_flag + 1

        IlutBitsParent%len_tot = IlutBitsParent%ind_source

        ! and if we use GUGA we have to enlarge this array by 3 entries
        if (tRDMOnfly .and. tGUGA) then
            IlutBitsParent%ind_rdm_ind = IlutBitsParent%ind_source + 1
            IlutBitsParent%ind_rdm_x0 = IlutBitsParent%ind_rdm_ind + 1
            IlutBitsParent%ind_rdm_x1 = IlutBitsParent%ind_rdm_x0 + 1

            IlutBitsParent%len_tot = IlutBitsParent%ind_rdm_x1
        end if

    end subroutine

    subroutine extract_bit_rep(ilut, nI, real_sgn, flags, j, store)

        ! Extract useful terms out of the bit-representation of a walker

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(out) :: nI(nel), flags
        integer, intent(in), optional :: j
        type(excit_gen_store_type), intent(inout), optional :: store
        real(dp), intent(out) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

        if (tStoredDets .and. present(j)) then
            nI = get_determinant(j)
        else if (tBuildOccVirtList .and. present(store)) then
            if (tBuildSpinSepLists) then
                call decode_bit_det_spinsep(nI, ilut, store)
            else
                call decode_bit_det_lists(nI, ilut, store)
            end if
        else
            call decode_bit_det(nI, ilut)
        end if

        sgn = iLut(IlutBits%ind_pop:IlutBits%ind_pop + lenof_sign - 1)
        real_sgn = transfer(sgn, real_sgn)

        flags = int(iLut(IlutBits%ind_flag))

    end subroutine extract_bit_rep

    !Extract all flags as a single integer
    function extract_flags(iLut) result(flags)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer :: flags

        flags = int(ilut(IlutBits%ind_flag))

    end function extract_flags

    !Extract the sign (as a real_dp) for a particular element in the lenof_sign "array"
    pure function extract_part_sign(ilut, part_type) result(real_sgn)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: part_type
        real(dp) :: real_sgn
        real_sgn = transfer(ilut(IlutBits%ind_pop + part_type - 1), real_sgn)
    end function

    pure function extract_run_sign(ilut, run) result(sgn)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: run
        HElement_t(dp) :: sgn

        ! Strange bug in compiler
        unused_var(run)
#ifdef CMPLX_
        sgn = cmplx(extract_part_sign(ilut, min_part_type(run)), extract_part_sign(ilut, max_part_type(run)), dp)
#else
        sgn = extract_part_sign(ilut, min_part_type(run))
#endif
    end function

    !From the determinants, array of signs, and flag integer, create the
    !"packaged walker" representation
    pure subroutine encode_bit_rep(ilut, Det, real_sgn, flag)
        integer(n_int), intent(out) :: ilut(0:nIfTot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        integer(n_int), intent(in) :: Det(0:IlutBits%len_orb)
        integer, intent(in) :: flag
        integer(n_int) :: sgn(lenof_sign)

        iLut(0:IlutBits%len_orb) = Det

        sgn = transfer(real_sgn, sgn)
        iLut(IlutBits%ind_pop:IlutBits%ind_pop + IlutBits%len_pop - 1) = sgn

        ilut(IlutBits%ind_flag) = int(flag, n_int)

    end subroutine encode_bit_rep

    subroutine encode_flags(ilut, flag)

        ! Add new flag information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flag

        iLut(IlutBits%ind_flag) = int(flag, n_int)

    end subroutine encode_flags

    pure function get_initiator_flag(sgn_index) result(flag)
        integer, intent(in) :: sgn_index
        integer :: flag
        ! Strange bug in compiler
        unused_var(sgn_index)
        ! map 1->1, 2->1, 3->3, 4->3, 5->5, 6->5 for complex,
        ! as the initiator flag is stored in the "real" bit
        ! of each run
        flag = flag_initiator(min_part_type(part_type_to_run(sgn_index)))
    end function get_initiator_flag

    pure function get_initiator_flag_by_run(run) result(flag)
        integer, intent(in) :: run
        integer :: flag
        ! Strange bug in compiler
        unused_var(run)
        ! map 1->1, 2->3, 3->5, 4->7 for complex
        flag = flag_initiator(min_part_type(run))
    end function get_initiator_flag_by_run

    pure function any_run_is_initiator(ilut) result(t)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer :: run
        logical :: t
        t = .false.
        do run = 1, inum_runs
            if (test_flag(ilut, get_initiator_flag_by_run(run))) then
                t = .true.
                return
            end if
        end do
    end function any_run_is_initiator

    pure function all_runs_are_initiator(ilut) result(t)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer :: run
        logical :: t
        t = .true.
        do run = 1, inum_runs
            if (.not. test_flag(ilut, get_initiator_flag_by_run(run))) then
                t = .false.
                return
            end if
        end do
    end function all_runs_are_initiator

    subroutine clear_all_flags(ilut)

        ! Clear all of the flags

        integer(n_int), intent(inout) :: ilut(0:niftot)

        ilut(IlutBits%ind_flag) = 0_n_int

    end subroutine clear_all_flags

    subroutine encode_sign(ilut, real_sgn)

        ! Add new sign information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        real(dp), intent(in) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

        sgn = transfer(real_sgn, sgn)
        iLut(IlutBits%ind_pop:IlutBits%ind_pop + IlutBits%len_pop - 1) = sgn

    end subroutine encode_sign

    subroutine encode_part_sign(ilut, real_sgn, part_type)

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
        iLut(IlutBits%ind_pop + part_type - 1) = sgn

    end subroutine encode_part_sign

    subroutine nullify_ilut(ilut)

        ! Sets the sign of a determinant to equal zero.
        integer(n_int), intent(inout) :: ilut(0:NIfTot)

        iLut(IlutBits%ind_pop:IlutBits%ind_pop + IlutBits%len_pop - 1) &
            = transfer(0.0_dp, 0_n_int)

    end subroutine

    subroutine nullify_ilut_part(ilut, part_type)

        ! Sets the sign of the walkers of the specified particle type on
        ! a determinant to equal zero.
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(in) :: part_type

        iLut(IlutBits%ind_pop + part_type - 1) = transfer(0.0_dp, 0_n_int)

    end subroutine

    subroutine set_flag_general(ilut, flg, state)

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
            call set_flag_single(ilut, flg)
        else
            call clr_flag(ilut, flg)
        end if
    end subroutine set_flag_general

    subroutine set_flag_single(ilut, flg)

        ! Set the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to set
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg

        ! This now assumes that we do not have more flags than bits in an
        ! integer.
        ilut(IlutBits%ind_flag) = ibset(ilut(IlutBits%ind_flag), flg)

    end subroutine set_flag_single

    subroutine clr_flag(ilut, flg)

        ! Clear the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to clear
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg

!This now assumes that we do not have more flags than bits in an integer.
        ilut(IlutBits%ind_flag) = ibclr(ilut(IlutBits%ind_flag), flg)

    end subroutine clr_flag

    subroutine clr_flag_multi(ilut, flg)

        ! Clear the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to clear
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg(:)

        integer :: i

        do i = 1, size(flg)
            ilut(IlutBits%ind_flag) = ibclr(ilut(IlutBits%ind_flag), flg(i))
        end do

    end subroutine clr_flag_multi

    function bit_parent_zero(ilut) result(zero)

        ! Used by the RDM functions
        ! Is the communicated parent zero?

        integer(n_int), intent(in) :: ilut(0:IlutBits%len_bcast)
        logical :: zero
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'bit_parent_zero'
#endif

        ASSERT(bit_rdm_init)

        zero = all(ilut(IlutBits%ind_parent:IlutBits%ind_parent + IlutBits%len_orb) == 0)

    end function

    subroutine encode_parent(ilut, ilut_parent, RDMBiasFacCurr)

        integer(n_int), intent(inout) :: ilut(0:IlutBits%len_bcast)
        integer(n_int), intent(in) :: ilut_parent(0:NIfTot)
        real(dp), intent(in) :: RDMBiasFacCurr
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'encode_parent'
#endif

        ASSERT(bit_rdm_init)

        ilut(IlutBits%ind_parent:IlutBits%ind_parent + IlutBits%len_orb) &
            = ilut_parent(0:IlutBits%len_orb)

        ilut(IlutBits%ind_rdm_fac) = &
            transfer(RDMBiasFacCurr, ilut(IlutBits%ind_rdm_fac))
        ! store the flag
        ilut(IlutBits%ind_parent_flag) = ilut_parent(IlutBits%ind_flag)

        ! in GUGA we also need to encode the rdm information
        ! BUT be careful here! the encoding should actually be done with
        ! IlutBits here, as the input ilut is of len_bcast length
        if (tGUGA) then
            call transfer_stochastic_rdm_info(ilut_parent, ilut, &
                                              BitIndex_from=IlutBits, BitIndex_to=IlutBits)
        end if

    end subroutine

    subroutine zero_parent(ilut)

        integer(n_int), intent(inout) :: ilut(0:IlutBits%len_bcast)
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'zero_parent'
#endif

        ASSERT(bit_rdm_init)

        ! is it intentional that the flag does not get zeroed?
        ilut(IlutBits%ind_parent:IlutBits%ind_rdm_fac) = 0_n_int

    end subroutine

    subroutine encode_spawn_hdiag(ilut, hel)

        integer(n_int), intent(inout) :: ilut(0:IlutBits%len_bcast)
        HElement_t(dp), intent(in) :: hel
#ifdef CMPLX_
        routine_name("encode_spawn_hdiag")
#endif

#ifdef CMPLX_
        ! Properly ensure that complex uses two words instead of one
        call stop_all(this_routine, "not implemented for complex")
        unused_var(ilut)
        unused_var(hel)
#else
        ilut(IlutBits%ind_hdiag) = transfer(hel, ilut(IlutBits%ind_hdiag))
#endif

    end subroutine encode_spawn_hdiag

    function extract_spawn_hdiag(ilut) result(hel)

        integer(n_int), intent(in) :: ilut(0:IlutBits%len_bcast)

        HElement_t(dp) :: hel
#ifdef CMPLX_
        routine_name("extract_spawn_hdiag")
#endif

#ifdef CMPLX_
        ! Properly ensure that complex uses two words instead of one
        call stop_all(this_routine, "not implemented for complex")
        unused_var(ilut)
        hel = 0._dp
#else
        hel = transfer(ilut(IlutBits%ind_hdiag), hel)
#endif

    end function extract_spawn_hdiag

    subroutine log_spawn(ilut)

        ! set the spawn counter to 1
        integer(n_int), intent(inout) :: ilut(0:IlutBits%len_bcast)

        ilut(IlutBits%ind_spawn) = 1
    end subroutine log_spawn

    subroutine increase_spawn_counter(ilut)
        ! increase the spawn counter by 1
        integer(n_int), intent(inout) :: ilut(0:IlutBits%len_bcast)

        ilut(IlutBits%ind_spawn) = ilut(IlutBits%ind_spawn) + 1

    end subroutine increase_spawn_counter

    function get_num_spawns(ilut) result(nSpawn)
        ! read the number of spawns to this det so far
        integer(n_int), intent(inout) :: ilut(0:IlutBits%len_bcast)
        integer :: nSpawn

        nSpawn = int(ilut(IlutBits%ind_spawn))

    end function get_num_spawns

    ! function test_flag is in bit_rep_data
    ! This avoids a circular dependence with DetBitOps.

    subroutine encode_det(ilut, Det)

        ! Add new det information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer(n_int), intent(in) :: Det(0:nifd)

        iLut(0:nifd) = Det

    end subroutine encode_det

    subroutine decode_bit_det_lists(nI, iLut, store)

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
                end if
            end do
            if (elec == nel_loc) exit
        end do

        ! Give final class count
        store%ClassCountUnocc = OrbClassCount - store%ClassCountOcc
        store%tFilled = .true.
        store%scratch3(1) = -1

        ! Fill in the remaineder of the virtuals list
        forall (ind=1:ScratchSize)
        !if (virt(ind) /= store%ClassCountUnocc(ind)) then
        !&<
            store%virt_list(virt(ind) + 1 : store%ClassCountUnocc(ind), ind) &
                = SymLabelList2(SymLabelCounts2(1, ind) + virt(ind) + store%ClassCountOcc(ind) : &
                                SymLabelCounts2(1, ind) + OrbClassCount(ind) - 1)
        !&>
        ! end if
        endforall

    end subroutine

    pure function getExcitationType(ExMat, IC) result(exTypeFlag)
        integer, intent(in) :: ExMat(2, ic), IC
        integer :: exTypeFlag

        ! i need to initialize to something..
        exTypeFlag = -1
        if (IC == 1) then
            if (tGUGA) then
                exTypeFlag = 1
            else
                if (is_beta(ExMat(2, 1)) .neqv. is_beta(ExMat(1, 1))) then
                    exTypeFlag = 3
                    return
                else
                    exTypeFlag = 1
                end if
            end if

        else if (IC == 2) then
            if (tGUGA) then
                ! in GUGA avoid the further checks down below
                exTypeFlag = 2
            else
                if (is_beta(ExMat(1, 1)) .and. is_beta(ExMat(1, 2))) then
                    ! elec orbs are both beta
                    if (is_beta(ExMat(2, 1)) .and. is_beta(ExMat(2, 2))) then
                        ! virt orbs are both beta
                        exTypeFlag = 2
                        return
                    else if (is_alpha(ExMat(2, 1)) .and. is_alpha(ExMat(2, 2))) then
                        ! virt orbs are both alpha
                        exTypeFlag = 5
                        return
                    else
                        ! one of the spins changes
                        exTypeFlag = 4
                        return
                    end if
                else if (is_alpha(ExMat(1, 1)) .and. is_alpha(ExMat(1, 2))) then
                    ! elec orbs are both alpha
                    if (is_alpha(ExMat(2, 1)) .and. is_alpha(ExMat(2, 2))) then
                        ! virt orbs are both alpha
                        exTypeFlag = 2
                        return
                    else if (is_beta(ExMat(2, 1)) .and. is_beta(ExMat(2, 2))) then
                        ! virt orbs are both beta
                        exTypeFlag = 5
                        return
                    else
                        ! one of the spins changes
                        exTypeFlag = 4
                        return
                    end if
                else
                    ! elec orb spins are different
                    if (is_beta(ExMat(2, 1)) .neqv. is_beta(ExMat(2, 2))) then
                        ! virt orbs are of opposite spin
                        exTypeFlag = 2
                        return
                    else
                        ! virt orbs are of the same spin
                        exTypeFlag = 4
                        return
                    end if
                end if
            end if
        else if (ic == 3) then
            ! todo! need to consider more maybe!
            exTypeFlag = 6
        end if

    end function

    subroutine decode_bit_det_spinsep(nI, iLut, store)

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
                    if (mod(ind, 2) == 1) then
                        ! alpha
                        store%nel_alpha = store%nel_alpha + 1
                        store%nI_alpha(store%nel_alpha) = orb
                        store%nI_alpha_inds(store%nel_alpha) = elec
                    else
                        store%nI_beta(elec - store%nel_alpha) = orb
                        store%nI_beta_inds(elec - store%nel_alpha) = elec
                    end if

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
                end if
            end do
            if (elec == nel_loc) exit
        end do

        ! Give final class count
        store%ClassCountUnocc = OrbClassCount - store%ClassCountOcc
        store%tFilled = .true.
        store%scratch3(1) = -1

        ! Fill in the remainder of the virtuals list
        forall (ind=1:ScratchSize)
        !if (virt(ind) /= store%ClassCountUnocc(ind)) then
        !&<
            store%virt_list(virt(ind) + 1 : store%ClassCountUnocc(ind), ind) &
                = SymLabelList2(SymLabelCounts2(1, ind) + virt(ind) + store%ClassCountOcc(ind): &
                                SymLabelCounts2(1, ind) + OrbClassCount(ind) - 1)
        !&>
         ! end if
        endforall

    end subroutine

    pure subroutine decode_bit_det_chunks(nI, iLut)

        ! This is a routine to take a determinant in bit form and construct
        ! the natural ordered integer form of the det.
        ! If CSFs are enabled, transfer the Yamanouchi symbol as well.

        integer(n_int), intent(in) :: ilut(0:NIftot)
        integer, intent(out) :: nI(:)
        integer :: i, j, k, val, elec, offset
        integer :: nel_loc

        nel_loc = size(nI)

        elec = 0
        offset = 0
        do i = 0, NIfD
            do j = 0, bits_n_int - 1, 8
                val = int(iand(ishft(ilut(i), -j), int(255, n_int)))
                do k = 1, decode_map_arr(0, val)
                    elec = elec + 1
                    nI(elec) = offset + decode_map_arr(k, val)
                    if (elec == nel_loc) return ! exit
                end do
                offset = offset + 8
            end do
        end do

    end subroutine

    subroutine add_ilut_lists(ndets_1, ndets_2, sorted_lists, list_1, list_2, list_out, &
                              ndets_out, prefactor)

        ! WARNING 1: This routine assumes that both list_1 and list_2 contain no
        ! repeated iluts, even if one of the repeated iluts has zero amplitude.

        ! WARNING 2: If the input lists are not sorted (as defined by ilut_gt)
        ! then sorted_lists should be input as .false., and the lists will then
        ! be sorted. This routine will not work if unsorted lists are passed in
        ! and sorted_list is input as .true.

        integer, intent(in) :: ndets_1
        integer, intent(in) :: ndets_2
        logical, intent(in) :: sorted_lists
        integer(n_int), intent(inout) :: list_1(0:, 1:)
        integer(n_int), intent(inout) :: list_2(0:, 1:)
        integer(n_int), intent(inout) :: list_out(0:, 1:)
        integer, intent(out) :: ndets_out
        real(dp), intent(in), optional :: prefactor

        integer :: i, pos, min_ind
        real(dp) :: sign_1(lenof_sign), sign_2(lenof_sign), sign_out(lenof_sign)
        real(dp) :: prefactor_

        def_default(prefactor_, prefactor, 1.0_dp)

        if (.not. sorted_lists) then
            write(stdout, *) lbound(list_1, 1), ubound(list_1, 1), lbound(list_2, 1), ubound(list_2, 1)
            call neci_flush(stdout)
            call sort(list_1(:, 1:ndets_1), ilut_lt, ilut_gt)
            call sort(list_2(:, 1:ndets_2), ilut_lt, ilut_gt)
        end if
        ndets_out = 0
        ! Where to start searching from in list 1:
        min_ind = 1

        do i = 1, ndets_2
            ! If list_2(:,i) is in list 1 then pos will equal the position it
            ! occupies in list 1.
            ! If list_2(:,i) is not in list 1 then -pos will equal the position
            ! that it should go in, to mantain the sorted ordering.
            pos = binary_search_custom(list_1(:, min_ind:ndets_1), list_2(:, i), niftot, ilut_gt)

            if (pos > 0) then
                ! Move all the states from list 1 before min_ind+pos-1 across
                ! to the combined list.
                list_out(:, ndets_out + 1:ndets_out + pos - 1) = list_1(:, min_ind:min_ind + pos - 2)
                ndets_out = ndets_out + pos - 1

                ndets_out = ndets_out + 1
                call extract_sign(list_1(:, min_ind + pos - 1), sign_1)
                call extract_sign(list_2(:, i), sign_2)
                sign_out = sign_1 + prefactor_ * sign_2
                list_out(:, ndets_out) = list_2(:, i)
                call encode_sign(list_out(:, ndets_out), sign_out)

                ! Search a smaller section of list_1 next time.
                min_ind = min_ind + pos
            else
                ! We have a state in list 2 which is not in list 1. Its
                ! position, if it were in list 1 would be min_ind-pos-1. Thus,
                ! first copy across all states from min_ind to min_ind-pos-2
                ! from list 1. Then copy across the state from list 2.
                list_out(:, ndets_out + 1:ndets_out - pos - 1) = list_1(:, min_ind:min_ind - pos - 2)
                ndets_out = ndets_out - pos

                list_out(:, ndets_out) = list_2(:, i)
                call extract_sign(list_out(:, ndets_out), sign_2)
                call encode_sign(list_out(:, ndets_out), sign_2 * prefactor_)

                ! Search a smaller section of list_1 next time.
                min_ind = min_ind - pos - 1
            end if
        end do

    end subroutine add_ilut_lists

    ! Write bit-determinant NI to unit NUnit.  Set LTerm if to add a newline at end.  Also prints CSFs
    subroutine writebitdet(nunit, ilutni, lterm)
        integer nunit, ni(nel)
        integer(kind=n_int) :: ilutni(0:niftot)
        logical lterm
        call decode_bit_det(ni, ilutni)
        call write_det(nunit, ni, lterm)
    end
end module bit_reps
