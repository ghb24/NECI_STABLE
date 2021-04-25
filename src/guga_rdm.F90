#include "macros.h"

module guga_rdm
    ! RDM module specifically for the GUGA spin-adapted implementation

    use constants, only: n_int, dp, lenof_sign, EPS, sizeof_int, int_rdm, bn2_, &
                         Root2, int64, int_rdm
    use SystemData, only: nel, nSpatOrbs, current_stepvector, currentB_ilut
    use bit_reps, only: extract_bit_rep, decode_bit_det, niftot, nifd, &
                        any_run_is_initiator, all_runs_are_initiator
    use rdm_data, only: one_rdms, two_rdm_spawn, rdmCorrectionFactor
    use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs, rdm_spawn_t, one_rdm_t
    use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, rdm_list_t
    use rdm_data, only: Sing_ExcList, Doub_ExcList, OneEl_Gap, TwoEl_Gap
    use DetBitOps, only: EncodeBitDet, count_open_orbs, DetBitEq, ilut_lt, ilut_gt
    use load_balance_calcnodes, only: DetermineDetNode
    use guga_excitations, only: excitationIdentifier
    use guga_excitations, only: init_singleWeight, calcRemainingSwitches_excitInfo_single
    use guga_excitations, only: createSingleStart, singleUpdate, singleEnd
    use guga_excitations, only: checkCompatibility, print_excitInfo
    use guga_excitations, only: calcDoubleExcitation_withWeight, &
                                calcNonOverlapDouble, calcSingleOverlapLowering, &
                                calcSingleOverlapRaising, calcSingleOverlapMixed, &
                                calcDoubleLowering, calcDoubleRaising, &
                                calcDoubleL2R, calcDoubleR2L, calcFullstopLowering, &
                                calcFullstopRaising, calcFullStopL2R, &
                                calcFullStopR2L, calcFullStartLowering, &
                                calcFulLStartRaising, calcFullStartL2R, &
                                calcFullStartR2L, calcFullStartFullStopAlike, &
                                calcFullStartFullStopMixed, &
                                calcRemainingSwitches_excitInfo_double, &
                                calc_guga_matrix_element
    use guga_data, only: ExcitationInformation_t, tag_tmp_excits, tag_excitations, &
                         excit_type, gen_type, rdm_ind_bitmask, excit_names, &
                         RdmContribList_t
    use guga_data, only: getDoubleMatrixElement, funA_0_2overR2, funA_m1_1_overR2, &
                         funA_3_1_overR2, funA_2_0_overR2, minFunA_2_0_overR2, &
                         minFunA_0_2_overR2, getDoubleContribution, getMixedFullStop
    use guga_types, only: WeightObj_t
    use guga_bitRepOps, only: update_matrix_element, setDeltaB, extract_matrix_element
    use guga_bitRepOps, only: isProperCSF_ilut, isDouble, init_csf_information
    use guga_bitRepOps, only: write_guga_list, write_det_guga, getSpatialOccupation
    use guga_bitRepOps, only: convert_ilut_toGUGA, convert_ilut_toNECI
    use guga_bitRepOps, only: calc_csf_info, add_guga_lists, EncodeBitDet_guga
    use guga_bitRepOps, only: findFirstSwitch, findLastSwitch
    use guga_bitRepOps, only: contract_1_rdm_ind, contract_2_rdm_ind, &
                              extract_1_rdm_ind, extract_2_rdm_ind, &
                              encode_rdm_ind, extract_rdm_ind, &
                              encode_stochastic_rdm_info, &
                              extract_excit_type_rdm, extract_excit_lvl_rdm, &
                              transfer_stochastic_rdm_info, &
                              extract_stochastic_rdm_info
    use MemoryManager, only: LogMemAlloc, LogMemDealloc
    use bit_rep_data, only: GugaBits, IlutBits
    use FciMCData, only: projEDet, CurrentDets, TotWalkers, ilutref, HFDet_True, &
                         iLutHF_True
    use LoggingData, only: ThreshOccRDM, tThreshOccRDMDiag, RDMExcitLevel, &
                           tExplicitAllRDM
    use UMatCache, only: gtID
    use RotateOrbsData, only: SymLabelListInv_rot
    use CalcData, only: tAdaptiveShift
    use Parallel_neci, only: nProcessors, MPIArg, MPIAlltoAll, MPIAlltoAllv, &
                             iProcIndex, MPISumAll, MPIBarrier
    use searching, only: BinSearchParts_rdm
    use rdm_data_utils, only: add_to_rdm_spawn_t, extract_sign_rdm
    use OneEInts, only: GetTMatEl
    use procedure_pointers, only: get_umat_el
    use guga_matrixElements, only: calcDiagExchangeGUGA_nI
    use util_mod, only: operator(.div.), near_zero, operator(.isclose.)
    use sort_mod, only: sort
    use rdm_data, only: rdm_list_t, rdm_definitions_t
    use util_mod, only: get_free_unit
    use dSFMT_interface, only: genrand_real2_dSFMT
    use ParallelHelper, only: root

    implicit none

    private
    public :: calc_rdm_energy_guga, gen_exc_djs_guga, &
              send_proc_ex_djs, &
              calc_all_excits_guga_rdm_singles, calc_explicit_1_rdm_guga, &
              calc_all_excits_guga_rdm_doubles, calc_explicit_2_rdm_guga, &
              calc_explicit_diag_2_rdm_guga, &
              Add_RDM_From_IJ_Pair_GUGA, fill_diag_1rdm_guga, &
              Add_RDM_HFConnections_GUGA, fill_spawn_rdm_diag_guga, &
              combine_x0_x1, pure_rdm_ind, generator_sign, &
              create_all_rdm_contribs, contract_molcas_1_rdm_index, &
              extract_molcas_1_rdm_index, contract_molcas_2_rdm_index, &
              extract_molcas_2_rdm_index, output_molcas_rdms, &
              fill_sings_1rdm_guga, &
              fill_sings_2rdm_guga, add_rdm_from_ij_pair_guga_exact, &
              conjugate_rdm_ind

    interface conjugate_rdm_ind
        module procedure conjugate_rdm_ind_vec
        module procedure conjugate_rdm_ind_scalar
    end interface conjugate_rdm_ind

contains

    pure function conjugate_rdm_ind_vec(rdm_ind, order) result(conj_rdm_ind)
        integer(int_rdm), intent(in) :: rdm_ind(:)
        integer, intent(in) :: order
        integer(int_rdm), allocatable :: conj_rdm_ind(:)

        integer :: i, j, k, l, m

        allocate(conj_rdm_ind(size(rdm_ind)), source=0_int_rdm)

        if (order == 1) then

            do m = 1, size(rdm_ind)
                call extract_1_rdm_ind(rdm_ind(m), i, j)
                conj_rdm_ind(m) = contract_1_rdm_ind(j, i)
            end do
        else if (order == 2) then
            do m = 1, size(rdm_ind)
                call extract_2_rdm_ind(rdm_ind(m), i, j, k, l)
                conj_rdm_ind(m) = contract_2_rdm_ind(k, l, i, j)
            end do
        end if

    end function conjugate_rdm_ind_vec

    pure function conjugate_rdm_ind_scalar(rdm_ind, order) result(conj_rdm_ind)
        integer(int_rdm), intent(in) :: rdm_ind
        integer, intent(in) :: order
        integer(int_rdm) :: conj_rdm_ind

        integer :: i, j, k, l, m

        if (order == 1) then
            call extract_1_rdm_ind(rdm_ind, i, j)
            conj_rdm_ind = contract_1_rdm_ind(j, i)
        else if (order == 2) then
            call extract_2_rdm_ind(rdm_ind, i, j, k, l)
            conj_rdm_ind = contract_2_rdm_ind(k, l, i, j)
        end if

    end function conjugate_rdm_ind_scalar

    pure function contract_molcas_1_rdm_index(i, j) result(ij)
        ! function which uses the molcas RDM index convention
        integer, intent(in) :: i, j
        integer :: ij

        integer :: p, q

        p = max(i, j)
        q = min(i, j)

        ij = q + p * (p - 1) / 2

    end function contract_molcas_1_rdm_index

    pure subroutine extract_molcas_1_rdm_index(pq, p, q)
        ! function which extracts the orbital indices following molcas
        ! convention
        integer, intent(in) :: pq
        integer, intent(out) :: p, q

        p = int(ceiling(-0.5 + sqrt(2.0 * pq)))
        q = pq - p * (p - 1) / 2

    end subroutine extract_molcas_1_rdm_index

    pure function contract_molcas_2_rdm_index(p, q, r, s) result(pqrs)
        ! function using the molcas rdm index convention
        integer, intent(in) :: p, q, r, s
        integer :: pqrs

        integer :: pq, rs, ij, kl

        pq = contract_molcas_1_rdm_index(p, q)
        rs = contract_molcas_1_rdm_index(r, s)

        ij = max(pq, rs)
        kl = min(pq, rs)

        pqrs = kl + ij * (ij - 1) / 2

    end function contract_molcas_2_rdm_index

    pure subroutine extract_molcas_2_rdm_index(pqrs, p, q, r, s, pq_out, rs_out)
        ! function using the molcas 2 rdm index convention
        integer, intent(in) :: pqrs
        integer, intent(out) :: p, q, r, s
        integer, intent(out), optional :: pq_out, rs_out

        integer :: pq, rs

        call extract_molcas_1_rdm_index(pqrs, pq, rs)

        call extract_molcas_1_rdm_index(pq, p, q)
        call extract_molcas_1_rdm_index(rs, r, s)

        if (present(pq_out)) pq_out = pq
        if (present(rs_out)) rs_out = rs

    end subroutine extract_molcas_2_rdm_index

    subroutine output_molcas_rdms(rdm_defs, rdm, rdm_trace)
        ! routine which prints spin-free GUGA RDMs directly in molcas format
        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        character(*), parameter :: this_routine = "output_molcas_rdms"

        real(dp), allocatable :: psmat(:), pamat(:), dmat(:)
        integer :: iunit_psmat, iunit_pamat, iunit_dmat, i
        integer :: p, q, r, s, pq, rs

        ! first fill the molcas rdms
        call fill_molcas_rdms(rdm_defs, rdm, rdm_trace, psmat, pamat, dmat)

        ! then print them
        ! Open the file to be written to.
        iunit_psmat = get_free_unit()
        iunit_pamat = get_free_unit()
        iunit_dmat = get_free_unit()

        ! only print on the root processor
        if (iProcIndex == root) then
            open(iunit_psmat, file='PSMAT', status='replace')
            do i = 1, size(psmat)
                if (abs(psmat(i)) > 1e-12_dp) then
                    !                 call extract_molcas_2_rdm_index(i, p, q, r, s, pq, rs)
                    !                 write(iunit_psmat, '(I6,G25.17,4I6)') i, psmat(i), p, q, r, s
                    write(iunit_psmat, '(I6, G25.17)') i, psmat(i)
                end if
            end do
            close(iunit_psmat)

            open(iunit_pamat, file='PAMAT', status='replace')
            do i = 1, size(pamat)
                if (abs(pamat(i)) > 1e-12_dp) then
                    !                 call extract_molcas_2_rdm_index(i, p, q, r, s, pq, rs)
                    !                 write(iunit_psmat, '(I6,G25.17,4I6)') i, pamat(i), p, q, r, s
                    write(iunit_pamat, '(I6, G25.17)') i, pamat(i)
                end if
            end do
            close(iunit_pamat)

            open(iunit_dmat, file='DMAT', status='replace')
            do i = 1, size(dmat)
                if (abs(dmat(i)) > 1e-12_dp) then
                    !                 call extract_molcas_1_rdm_index(i, p, q)
                    !                 write(iunit_dmat, '(I6, G25.17,2I6)') i, dmat(i), p, q
                    write(iunit_dmat, '(I6, G25.17)') i, dmat(i)
                end if
            end do
            close(iunit_dmat)
        end if

    end subroutine output_molcas_rdms

    subroutine fill_molcas_rdms(rdm_defs, rdm, rdm_trace, &
                                psmat, pamat, dmat)
        type(rdm_definitions_t), intent(in) :: rdm_defs
        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(in) :: rdm_trace(rdm%sign_length)
        real(dp), intent(out), allocatable :: psmat(:), pamat(:), dmat(:)
        character(*), parameter :: this_routine = "fill_molcas_rdms"

        integer :: n_one_rdm, n_two_rdm, iproc, irdm, ielem
        integer(int_rdm) :: pqrs, pq_, rs_
        integer :: pq, rs, pqrs_m, pq_m, rs_m, p, q, r, s, p_m, q_m, r_m, s_m
        real(dp) :: rdm_sign(rdm%sign_length), rdm_sign_, fac
        real(dp), allocatable :: dmat_loc(:), psmat_loc(:), pamat_loc(:)
        integer :: ierr

        n_one_rdm = nSpatorbs * (nSpatorbs + 1) / 2
        n_two_rdm = n_one_rdm * (n_one_rdm + 1) / 2

        allocate(dmat_loc(n_one_rdm), source=0.0_dp)
        allocate(psmat_loc(n_two_rdm), source=0.0_dp)
        allocate(pamat_loc(n_two_rdm), source=0.0_dp)

        ! first
        do iproc = 0, nProcessors - 1
            if (iproc == iProcIndex) then
                ! Loop over all RDMs beings sampled.
                do irdm = 1, rdm_defs%nrdms
                    do ielem = 1, rdm%nelements
                        pqrs = rdm%elements(0, ielem)
                        call extract_2_rdm_ind(pqrs, p, q, r, s, pq_, rs_)
                        pq = int(pq_)
                        rs = int(rs_)
                        pqrs_m = contract_molcas_2_rdm_index(p, q, r, s)
                        call extract_molcas_2_rdm_index(pqrs_m, &
                                                        p_m, q_m, r_m, s_m, pq_m, rs_m)

                        call extract_sign_rdm(rdm%elements(:, ielem), rdm_sign)
                        rdm_sign_ = rdm_sign(1)
                        ! now make the fill logic
                        ! the molcas rdm elements are given by
                        ! if r /= s (and probably here p /= q)
                        ! psmat_loc(pqrs) = (two_rdm(pqrs) + two_rdm(pqsr)) / 2
                        ! pamat_loc(pqrs) = (two_rdm(pqrs) - two_rdm(pqsr)) / 2
                        ! if r == s (and probably .or. p == q)
                        ! psmat_loc(pqrs) = 2 * two_rdm(pqrs)
                        ! es geht eigentlich nur drum wann das element
                        ! negativ zur anti-symmetrischen beiträgt..
                        ! und wenn es nur zur diagonalen beiträgt..

                        if (pq_m == rs_m) then
                            if (p_m == q_m) then
                                psmat_loc(pqrs_m) = psmat_loc(pqrs_m) + rdm_sign_ / 2.0_dp
                            else
                                psmat_loc(pqrs_m) = psmat_loc(pqrs_m) + rdm_sign_ / 4.0_dp
                                pamat_loc(pqrs_m) = pamat_loc(pqrs_m) + &
                                                    molcas_sign(p, q, r, s) * rdm_sign_ / 4.0_dp
                            end if
                        else
                            if (p_m == q_m) then
                                psmat_loc(pqrs_m) = psmat_loc(pqrs_m) + rdm_sign_ / 4.0_dp
                            else
                                psmat_loc(pqrs_m) = psmat_loc(pqrs_m) + rdm_sign_ / 8.0_dp
                                if (r_m /= s_m) then
                                    pamat_loc(pqrs_m) = pamat_loc(pqrs_m) + &
                                                        molcas_sign(p, q, r, s) * rdm_sign_ / 8.0_dp
                                end if
                            end if
                        end if

                        pq_m = contract_molcas_1_rdm_index(p, q)
                        rs_m = contract_molcas_1_rdm_index(r, s)
                        ! convert to 1-RDM
                        if (r == s .and. p == q) then
                            dmat_loc(pq_m) = dmat_loc(pq_m) + rdm_sign_
                            dmat_loc(rs_m) = dmat_loc(rs_m) + rdm_sign_
                        else
                            if (r == s) then
                                dmat_loc(pq_m) = dmat_loc(pq_m) + rdm_sign_ / 2.0_dp
                            else if (p == q) then
                                dmat_loc(rs_m) = dmat_loc(rs_m) + rdm_sign_ / 2.0_dp
                            end if
                        end if
                    end do
                end do
            end if
            ! Wait for the current processor to finish summing its RDM elements
            call MPIBarrier(ierr)

        end do

        psmat_loc = psmat_loc / rdm_trace(1)
        pamat_loc = pamat_loc / rdm_trace(1)
        dmat_loc = dmat_loc / (2.0 * rdm_trace(1) * real(nel - 1, dp))

        allocate(dmat(n_one_rdm), source=0.0_dp)
        allocate(psmat(n_two_rdm), source=0.0_dp)
        allocate(pamat(n_two_rdm), source=0.0_dp)

        call MPISumAll(dmat_loc, dmat)
        call MPISumAll(psmat_loc, psmat)
        call MPISumAll(pamat_loc, pamat)

    end subroutine fill_molcas_rdms

    pure function molcas_sign(p, q, r, s) result(sgn)
        ! gives me the sign to fill the anti-symmetric molcas RDM with
        integer, intent(in) :: p, q, r, s
        real(dp) :: sgn

        sgn = 1.0_dp

        ! i know p /= q and r /= s when coming here

        if (p < q) sgn = -sgn
        if (r < s) sgn = -sgn

    end function molcas_sign

    ! function get_rdm_hash_entry()
    !     ! make a function to easily extract the rdm element stored in the
    !     ! hash table
    !
    ! end function get_rdm_hash_entry

    subroutine print_rdm_ind(rdm_ind, typ, t_newline)
        integer(int_rdm), intent(in) :: rdm_ind
        integer, intent(in) :: typ
        logical, intent(in), optional :: t_newline

        logical :: t_newline_
        integer :: i, j, k, l
        character(3) :: adv

        def_default(t_newline_, t_newline, .true.)

        if (t_newline_) then
            adv = 'yes'
        else
            adv = 'no'
        end if

        if (typ == 1) then
            call extract_1_rdm_ind(rdm_ind, i, j)

            write(6, '(2i3)', advance=adv) i, j

        else if (typ == 2) then
            call extract_2_rdm_ind(rdm_ind, i, j, k, l)

            write(6, '(4i3)', advance=adv) i, j, k, l
        end if

    end subroutine print_rdm_ind

    elemental function pure_rdm_ind(rdm_ind) result(pure_ind)
        ! return 'just' the rdm index without the excit-level and type info
        integer(int_rdm), intent(in) :: rdm_ind
        integer(int_rdm) :: pure_ind

        pure_ind = iand(rdm_ind, rdm_ind_bitmask)

    end function pure_rdm_ind

    subroutine create_all_rdm_contribs(rdm_ind, x0, x1, rdm_ind_out, rdm_mat)
        ! I also need a function which creates me the remaining
        ! index combination of not directly sampled excitation in the
        ! GUGA excitation generation
        integer(int_rdm), intent(in) :: rdm_ind
        real(dp) :: x0, x1
        integer(int_rdm), intent(out), allocatable :: rdm_ind_out(:)
        real(dp), intent(out), allocatable :: rdm_mat(:)
        debug_function_name("create_all_rdm_contribs")
        integer :: ex_lvl, ex_typ, i, j, k, l

        ex_lvl = extract_excit_lvl_rdm(rdm_ind)
        ex_typ = extract_excit_type_rdm(rdm_ind)

        ASSERT(ex_lvl == 1 .or. ex_lvl == 2)
        ASSERT(ex_typ /= excit_type%invalid)

        if (ex_lvl == 1) then
            ! for singles we 'only' have the original index
            allocate(rdm_ind_out(1), source=rdm_ind)
            ! and the element is only the x0 contrib
            ASSERT(near_zero(x1))
            allocate(rdm_mat(1), source=x0)

        else if (ex_lvl == 2) then

            select case (ex_typ)

                ! we only take stuff from stochastic double excitations here..
                ! so this excludes some excitation types.. assert this below!
            case (excit_type%single_overlap_R_to_L, &
                  excit_type%single_overlap_L_to_R, &
                  excit_type%fullstop_lowering, &
                  excit_type%fullstop_raising, &
                  excit_type%fullstart_raising, &
                  excit_type%fullstart_lowering, &
                  excit_type%fullstart_stop_alike)

                ! here we also only have the ij <-> kl conjugation which
                ! can (and should for now) be handled outside of this function!
                allocate(rdm_ind_out(1), source=rdm_ind)
                ! also x1 must be 0 in these cases
                ASSERT(near_zero(x1))
                allocate(rdm_mat(1), source=x0)

                ! group everything together where x0 is 0:
                ! the conjugated elements of those are actually dealt with
                ! in the single excitations (here they would actually be 0
                ! since we enforce a spin-flip in the overlap range!
            case (excit_type%fullstop_L_to_R, &
                  excit_type%fullstop_R_to_L, &
                  excit_type%fullstart_L_to_R, &
                  excit_type%fullstart_R_to_L, &
                  excit_type%fullstart_stop_mixed)

                allocate(rdm_ind_out(1), source=rdm_ind)
                ASSERT(near_zero(x0))
                allocate(rdm_mat(1), source=x1)

            case (excit_type%double_lowering, &
                  excit_type%double_raising)

                ! here we need to do something
                call extract_2_rdm_ind(rdm_ind, i, j, k, l)
                ! a switch of j <-> l will induce a sign change between x0 and x1!
                allocate(rdm_ind_out(2), source=0_int_rdm)
                allocate(rdm_mat(2), source=0.0_dp)

                rdm_ind_out(1) = rdm_ind
                rdm_mat(1) = x0 + generator_sign(i, j, k, l) * x1

                rdm_ind_out(2) = contract_2_rdm_ind(i, l, k, j, &
                                                    excit_lvl=2, excit_typ=ex_typ)

                rdm_mat(2) = x0 + generator_sign(i, l, k, j) * x1

            case (excit_type%double_L_to_R_to_L, &
                  excit_type%double_R_to_L_to_R, &
                  excit_type%double_L_to_R, &
                  excit_type%double_R_to_L)

                ! these cases contain a non-overlap excitaion if no spin-flip
                ! in the overlap range happened (i hope for now that a
                ! zero x0 element indicates such a spin-flip.. although
                ! it should be fine, since when there is a spin-flip it
                ! is 0 definitely and even if there was, x0 is still 0 and this
                ! is all what is encoded for this flipped version!

                allocate(rdm_ind_out(2), source=0_int_rdm)
                allocate(rdm_mat(2), source=0.0_dp)

                ! the first is the original
                rdm_ind_out(1) = rdm_ind
                rdm_mat(1) = x0 + x1

                call extract_2_rdm_ind(rdm_ind, i, j, k, l)
                ! the second one does hopefully not need a type encoded
                ! in the rdm! since this is not available here!!
                rdm_ind_out(2) = contract_2_rdm_ind(i, l, k, j, excit_lvl=2)
                rdm_mat(2) = -2.0_dp * x0

            case default
                ! should not be here for now!
                ASSERT(.false.)

            end select
        else
            ! but should actually not be here now!
            ASSERT(.false.)
            ! other wise no contributions
            allocate(rdm_ind_out(0))
            allocate(rdm_mat(0))
        end if

    end subroutine create_all_rdm_contribs

    function combine_x0_x1(rdm_ind, x0, x1) result(comb)
        ! this function combines the x0 and x1 coupling coefficient
        ! components correctly depending on the provided rdm_index
        integer(int_rdm), intent(in) :: rdm_ind
        real(dp), intent(in) :: x0, x1
        real(dp) :: comb
        debug_function_name("combine_x0_x1")
        integer :: ex_lvl, ex_typ, i, j, k, l

        ex_lvl = extract_excit_lvl_rdm(rdm_ind)
        ex_typ = extract_excit_type_rdm(rdm_ind)

        ASSERT(ex_lvl == 1 .or. ex_lvl == 2)
        ASSERT(ex_typ /= excit_type%invalid)

        if (ex_lvl == 1) then
            ASSERT(near_zero(x1))
            comb = x0
        else if (ex_lvl == 2) then

            select case (ex_typ)

                ! group together everything where x1 is 0
            case (excit_type%single_overlap_L_to_R, &
                  excit_type%single_overlap_R_to_L, &
                  excit_type%fullstop_lowering, &
                  excit_type%fullstop_raising, &
                  excit_type%fullstart_raising, &
                  excit_type%fullstart_lowering, &
                  excit_type%fullstart_stop_alike)

                ASSERT(near_zero(x1))
                comb = x0

                ! group everything together where x0 is 0:
            case (excit_type%fullstop_L_to_R, &
                  excit_type%fullstop_R_to_L, &
                  excit_type%fullstart_L_to_R, &
                  excit_type%fullstart_R_to_L, &
                  excit_type%fullstart_stop_mixed)

                ASSERT(near_zero(x0))
                comb = x1

            case (excit_type%double_lowering, &
                  excit_type%double_raising)
                ! in these two cases i need to determine the sign influence
                ! from the generator ordering

                call extract_2_rdm_ind(rdm_ind, i, j, k, l)
                comb = x0 + generator_sign(i, j, k, l) * x1

            case (excit_type%double_L_to_R_to_L, &
                  excit_type%double_R_to_L_to_R, &
                  excit_type%double_L_to_R, &
                  excit_type%double_R_to_L)
                ! i think for these cases i really need nothing..
                comb = x0 + x1

            case default
                ! should not be here!
                ! maybe in the future i should expand this
                ! by analysing the indices and producing all of them..
                ASSERT(.false.)
            end select

        else
            comb = 0.0_dp
        end if

    end function combine_x0_x1

    pure function generator_sign(i, j, k, l) result(sgn)
        integer, intent(in) :: i, j, k, l
        real(dp) :: sgn

        ! standard value. only for double raising and lowering a -1 can happen
        sgn = 1.0_dp
        if (i < j .and. k < l) then
            ! double raising
            ! make rudimentary now.. think of increasing speed later
            if (i < k) then
                if (l < j) then
                    sgn = -1.0_dp
                end if
            else if (k < i) then
                if (j < l) then
                    sgn = -1.0_dp
                end if
            end if
        else if (j < i .and. l < k) then
            ! double lowering
            if (j < l) then
                if (k < i) then
                    sgn = -1.0_dp
                end if
            else if (l < j) then
                if (i < k) then
                    sgn = -1.0_dp
                end if
            end if
        end if

    end function generator_sign

    subroutine Add_RDM_HFConnections_GUGA(spawn, one_rdms, ilutJ, av_sign_j, &
                                          av_sign_hf, excit_lvl, iter_rdm)
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: ilutJ(0:IlutBits%len_tot)
        integer, intent(in) :: excit_lvl, iter_rdm(:)
        real(dp), intent(in) :: av_sign_j(:), av_sign_hf(:)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "Add_RDM_HFConnections_GUGA"
#endif

        ! damn.. here we need to do the 'slow' implementation i guess..
        ! since NJ does not come from a spawning event but is
        ! done deterministically for the HF connections..
        ! there should be a clever way to do this..
        ! nah.. not for now.. otherwise i have to check everywhere also
        ! if i sampled this already.. for leave it at that and be done with
        ! it!
        if (excit_lvl == 1 .or. excit_lvl == 2) then
            ! for HF -> nJ we do not have csf info intialized so use calc_type = 2
            call add_rdm_from_ij_pair_guga_exact(spawn, one_rdms, iLutHF_True, &
                                                 ilutJ, av_sign_hf(2::2), iter_rdm*av_sign_j(1::2), calc_type=2)

            ! for nJ we have the csf info initialized (or maybe not..)
            call add_rdm_from_ij_pair_guga_exact(spawn, one_rdms, ilutJ, &
                                                 iLutHF_True, av_sign_j(2::2), iter_rdm*av_sign_hf(1::2), calc_type=2)
        end if

    end subroutine Add_RDM_HFConnections_GUGA

    subroutine add_rdm_from_ij_pair_guga_exact(spawn, one_rdms, ilutI, ilutJ, &
                                               sign_i, sign_j, calc_type)
        ! this routine is called for RDM sampling within the semi-stochastic
        ! space or for the connection to the 'true' HF det!
        ! i also need the calc-type input, as in the semi-stochastic space
        ! we initialize the csf info of nI, while for the HF connections we dont
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer(n_int), intent(in) :: ilutI(0:IlutBits%len_tot), &
                                      ilutJ(0:IlutBits%len_tot)
        real(dp), intent(in) :: sign_i(:), sign_j(:)
        integer, intent(in) :: calc_type
        character(*), parameter :: this_routine = "add_rdm_from_ij_pair_guga_exact"
        integer(int_rdm), allocatable :: rdm_ind(:)
        real(dp), allocatable :: rdm_mat(:)
        type(ExcitationInformation_t) :: excitInfo
        HElement_t(dp) :: mat_ele
        integer :: p, q, r, s, n
        real(dp) :: full_sign(spawn%rdm_send%sign_length)

        call calc_guga_matrix_element(ilutI, ilutJ, excitInfo, mat_ele, &
                                      t_hamil=.false., calc_type=calc_type, rdm_ind=rdm_ind, &
                                      rdm_mat=rdm_mat)

        ! i assume sign_i and sign_j are not 0 if we end up here..
        do n = 1, size(rdm_ind)
            if (.not. near_zero(rdm_mat(n))) then
                if (excitInfo%excitLvl == 1) then
                    if (RDMExcitLevel == 1) then
                        call fill_sings_1rdm_guga(one_rdms, sign_i, sign_j, &
                                                  rdm_mat(n), rdm_ind(n))
                    else
                        call fill_sings_2rdm_guga(spawn, ilutI, &
                                                  ilutJ, sign_i, sign_j, rdm_mat(n), rdm_ind(n))
                    end if
                else if (excitInfo%excitLvl == 2 .and. RDMExcitLevel /= 1) then
                    call extract_2_rdm_ind(rdm_ind(n), p, q, r, s)
                    full_sign = sign_i * sign_j * rdm_mat(n)

                    ! here in the 'exact' filling (coming from HF or
                    ! within the semistochastic space I think it makes
                    ! sense to fill symmetrically.. since here no
                    ! stochastic spawning is happening and this does not
                    ! give us information about the hermiticity error!
                    call add_to_rdm_spawn_t(spawn, p, q, r, s, &
                                            full_sign, .true.)

                    if (.not. &
                        (excitInfo%typ == excit_type%fullstart_stop_alike)) then
                        call add_to_rdm_spawn_t(spawn, r, s, p, q, &
                                                full_sign, .true.)
                    end if
                end if
            end if
        end do

    end subroutine add_rdm_from_ij_pair_guga_exact

    subroutine Add_RDM_From_IJ_Pair_GUGA(spawn, one_rdms, nI, nJ, sign_i, &
                                         sign_j, t_bra_to_ket, t_fast, rdm_ind_in, x0, x1)
        ! corresponding GUGA routine from function found in rdm_filling.F90
        type(rdm_spawn_t), intent(inout) :: spawn
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(nel), nJ(nel)
        real(dp), intent(in) :: sign_i(:), sign_j(:)
        logical, intent(in) :: t_bra_to_ket
        logical, intent(in), optional :: t_fast
        integer(int_rdm), intent(in), optional :: rdm_ind_in
        real(dp), intent(in), optional :: x0, x1
#ifdef DEBUG_
        character(*), parameter :: this_routine = "Add_RDM_From_IJ_Pair_GUGA"
#endif
        type(ExcitationInformation_t) :: excitInfo
        HElement_t(dp) :: mat_ele
        integer(int_rdm), allocatable :: rdm_ind(:)
        real(dp), allocatable :: rdm_mat(:)
        integer :: i, j, k, l, n, ex_lvl, ex_typ, pick, m
        real(dp) :: full_sign(spawn%rdm_send%sign_length)
        integer(n_int) :: ilutGi(0:GugaBits%len_tot), ilutGj(0:GugaBits%len_tot)
        integer(int_rdm) :: pure_ind
        logical :: t_fast_

        def_default(t_fast_, t_fast, .true.)

        if (t_fast_) then

            ASSERT(present(rdm_ind_in))
            ASSERT(present(x0))
            ASSERT(present(x1))
            ASSERT(rdm_ind_in /= 0_int_rdm)

            ex_lvl = extract_excit_lvl_rdm(rdm_ind_in)
            ex_typ = extract_excit_type_rdm(rdm_ind_in)

            ASSERT(ex_lvl == 1 .or. ex_lvl == 2)

            pure_ind = pure_rdm_ind(rdm_ind_in)

            if (ex_lvl == 1) then

                if (RDMExcitLevel == 1) then
                    call fill_sings_1rdm_guga(one_rdms, sign_I, sign_J, &
                                              x0, pure_ind)

                else
                    if (t_bra_to_ket) then
                        call EncodeBitDet_guga(nI, ilutGi)
                        call EncodeBitDet_guga(nJ, ilutGj)
                    else
                        call EncodeBitDet_guga(nI, ilutGJ)
                        call EncodeBitDet_guga(nJ, ilutGi)
                    end if

                    call fill_sings_2rdm_guga(spawn, IlutGi, &
                                              ilutGj, sign_i, sign_j, x0, pure_ind)

                end if

            else if (ex_lvl == 2 .and. RDMExcitLevel /= 1) then

                if (t_bra_to_ket) then
                    call EncodeBitDet_guga(nI, ilutGi)
                    call EncodeBitDet_guga(nJ, ilutGj)
                else
                    call EncodeBitDet_guga(nI, ilutGJ)
                    call EncodeBitDet_guga(nJ, ilutGi)
                end if
                call create_all_rdm_contribs(rdm_ind_in, x0, x1, rdm_ind, rdm_mat)

                do n = 1, size(rdm_ind)
                    if (.not. near_zero(rdm_mat(n))) then
                        call extract_2_rdm_ind(rdm_ind(n), i, j, k, l)

                        full_sign = sign_I * sign_J * rdm_mat(n)
                        if (.not. &
                            (ex_typ == excit_type%fullstart_stop_alike)) then
                            full_sign = 2.0_dp * full_sign
                        end if
                        call add_to_rdm_spawn_t(spawn, i, j, k, l, &
                                                full_sign, .true.)
                    end if
                end do
            end if
        else

            if (t_bra_to_ket) then
                call EncodeBitDet_guga(nI, ilutGi)
                call EncodeBitDet_guga(nJ, ilutGj)
            else
                call EncodeBitDet_guga(nI, ilutGJ)
                call EncodeBitDet_guga(nJ, ilutGi)
            end if

            call calc_guga_matrix_element(IlutGi, ilutGj, excitInfo, mat_ele, &
                                          t_hamil=.false., calc_type=2, rdm_ind=rdm_ind, &
                                          rdm_mat=rdm_mat)

            ! i assume sign_i and sign_j are not 0 if we end up here..
            do n = 1, size(rdm_ind)
                if (.not. near_zero(rdm_mat(n))) then
                    if (excitInfo%excitLvl == 1) then
                        if (RDMExcitLevel == 1) then
                            call fill_sings_1rdm_guga(one_rdms, sign_I, sign_J, &
                                                      rdm_mat(n), rdm_ind(n))
                        else
                            call fill_sings_2rdm_guga(spawn, IlutGi, &
                                                      ilutGj, sign_i, sign_j, rdm_mat(n), rdm_ind(n))
                        end if
                    else if (excitInfo%excitLvl == 2 .and. RDMExcitLevel /= 1) then
                        call extract_2_rdm_ind(rdm_ind(n), i, j, k, l)
                        full_sign = sign_I * sign_J * rdm_mat(n)

                        ! here in the 'exact' filling (coming from HF or
                        ! within the semistochastic space I think it makes
                        ! sense to fill symmetrically.. since here no
                        ! stochastic spawning is happening and this does not
                        ! give us information about the hermiticity error!
                        call add_to_rdm_spawn_t(spawn, i, j, k, l, &
                                                full_sign, .true.)

                        if (.not. &
                            (excitInfo%typ == excit_type%fullstart_stop_alike)) then
                            call add_to_rdm_spawn_t(spawn, k, l, i, j, &
                                                    full_sign, .true.)
                        end if
                    end if
                end if
            end do
        end if

    end subroutine Add_RDM_From_IJ_Pair_GUGA

    subroutine gen_exc_djs_guga(ilutI)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        character(*), parameter :: this_routine = "gen_exc_djs_guga"

        integer :: nI(nel), flags_I, n_singles, n_doubles
        real(dp) :: sign_i(lenof_sign), full_sign(1)

        integer(n_int), allocatable :: excits(:, :)
        integer(n_int) :: ilutG(0:GugaBits%len_tot)

        call extract_bit_rep(ilutI, nI, sign_I, flags_I)

        if (RDMExcitLevel == 1) then
            call fill_diag_1rdm_guga(one_rdms, nI, sign_I)
        else
            full_sign = sign_i(1) * sign_I(lenof_sign)
            call fill_spawn_rdm_diag_guga(two_rdm_spawn, nI, full_sign)
        end if

        call convert_ilut_toGUGA(ilutI, ilutG)

        ! one-rdm is always calculated
        ! calculate the excitations here
        ! but with my general two-body excitation routine i do not need
        ! to calculate this in case of more than singles:
        if (RDMExcitLevel == 1) then
            call calc_explicit_1_rdm_guga(ilutG, n_singles, excits)

            ! and then sort them correctly in the communicated array
            call assign_excits_to_proc_guga(n_singles, excits, 1)

            deallocate(excits)
            call LogMemDealloc(this_routine, tag_excitations)
        end if

        ! now to double excitations if requsted:
        if (RDMExcitLevel /= 1) then

            ! if i want to mimic stochastic RDM sampling I also
            ! have to explicitly create single excitations, but
            ! store them in the according 2-RDM entries
            call calc_explicit_1_rdm_guga(ilutG, n_singles, excits)

            ! and then sort them correctly in the communicated array
            call assign_excits_to_proc_guga(n_singles, excits, 1)

            deallocate(excits)
            call LogMemDealloc(this_routine, tag_excitations)

            call calc_explicit_2_rdm_guga(ilutG, n_doubles, excits)

            call assign_excits_to_proc_guga(n_doubles, excits, 2)

            deallocate(excits)
            call LogMemDealloc(this_routine, tag_excitations)
        end if

    end subroutine gen_exc_djs_guga

    subroutine fill_spawn_rdm_diag_guga(spawn, nI, full_sign)
        ! i have to write a routine, which correctly takes into
        ! account all the diagonal contributions from double excitations
        type(rdm_spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: full_sign(spawn%rdm_send%sign_length)
        character(*), parameter :: this_routine = "fill_spawn_rdm_diag_guga"

        integer :: i, s_orbs(nel), s, inc_i, j, inc_j, p
        integer :: ssss, spsp, sp, sp2, a, b, c, d

        real(dp) :: occ_i, occ_j, x0, x1
        logical :: t_test

        ! i have to figure out what exactly contributes to here..
        ! according to the paper about the two-RDM labels the
        ! diagonal elements of the 2-RDM are given by
        ! D_ij,kl = <psi| e_{ki,lj} | psi>
        ! which for the diagonal contribution kl = ij yields
        ! D_{ij,ij} = <psi| e_{ii,jj} |psi>
        ! which is definetly purely diagonal!

        ! i could try to use the add_to_rdm_spawn_t routine in
        ! rdm_data_utils with the according modified matrix elements.
        i = 1
        s_orbs = gtID(nI)

        ! TODO: I could also try to add the diagonal exchange
        ! contributions directly here! so this routine would be
        ! applicable in the stochastic spawning already
        do while (i <= nel)

            s = s_orbs(i)

            if (isDouble(nI, i)) then
                occ_i = 2.0_dp
                inc_i = 2

                ! i think here i have to put the full diagonal
                ! D_{ii,ii} entry, which counts double occupancies
                ! and i will try for now to reuse the following routine:
                ! put with the spatial orbital s
                ! but maybe I have to also add the 'switched' indices
                ! contribution.. which would mean a factor of 2..
                call add_to_rdm_spawn_t(spawn, s, s, s, s, occ_i * full_sign, .true.)

            else
                occ_i = 1.0_dp
                inc_i = 1
            end if

            j = i + inc_i

            do while (j <= nel)

                p = s_orbs(j)

                if (isDouble(nI, j)) then
                    occ_j = 2.0_dp
                    inc_j = 2
                else
                    occ_j = 1.0_dp
                    inc_j = 1
                end if

                ! i could also just multiply by 2 here, since this will
                ! get strored in the same D_{ij,ij} RDM element!
                ! but for now I decided to fill in all disting 2-RDM
                ! elements
                ! since here we do not have spawning this should also
                ! get stored symmetrically! as otherwise this will
                ! junk up the hermiticity error calculation!
                call add_to_rdm_spawn_t(spawn, s, s, p, p, &
                                        occ_i * occ_j * full_sign, .true.)
                call add_to_rdm_spawn_t(spawn, p, p, s, s, &
                                        occ_i * occ_j * full_sign, .true.)

                ! i can also add the fully diagonal exchange contributions
                ! here. this is also necessary to do, if I want to use
                ! this routine in the stochastic sampling

                ! but for open-shell to open-shell exchange excitations
                ! I have to calculate the correct x1 matrix element..

                ! x0 matrix element is easy
                x0 = -occ_i * occ_j / 2.0_dp

                if (inc_i == 1 .and. inc_j == 1) then
                    ! if we have open-shell to open chell
                    x1 = calcDiagExchangeGUGA_nI(i, j, nI) / 2.0_dp

                    call add_to_rdm_spawn_t(spawn, s, p, p, s, &
                                            (x0 - x1) * full_sign, .true.)
                    call add_to_rdm_spawn_t(spawn, p, s, s, p, &
                                            (x0 - x1) * full_sign, .true.)

                else
                    call add_to_rdm_spawn_t(spawn, s, p, p, s, &
                                            x0 * full_sign, .true.)
                    ! and the symmetric version:
                    call add_to_rdm_spawn_t(spawn, p, s, s, p, &
                                            x0 * full_sign, .true.)
                end if

                j = j + inc_j
            end do
            i = i + inc_i
        end do

    end subroutine fill_spawn_rdm_diag_guga

    subroutine fill_diag_1rdm_guga(one_rdms, nI, contrib_sign, &
                                   tCoreSpaceDetIn, RDMItersIn)
        ! I think I have to change these routines to also take into
        ! account the possible coupling coefficient of doubly occupied
        ! orbitals
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        integer, intent(in) :: nI(nel)
        real(dp), intent(in) :: contrib_sign(:)
        logical, optional, intent(in) :: tCoreSpaceDetIn
        integer, optional, intent(in) :: RDMItersIn(:)
        character(*), parameter :: this_routine = "fill_diag_1rdm_guga"

        integer :: i, ind, irdm, s, s_orbs(nel), inc
        integer :: RDMIters(size(one_rdms))
        real(dp) :: ScaleContribFac, final_contrib(size(one_rdms))
        logical :: tCoreSpaceDet

        ScaleContribFac = 1.0_dp

        RDMIters = 1
        if (present(RDMItersIn)) RDMIters = RDMItersIn

        tCoreSpaceDet = .false.
        if (present(tCoreSpaceDetIn)) tCoreSpaceDet = tCoreSpaceDetIn

        ! This is the single-run cutoff being applied (do not use in DR mode):
        if (.not. tCoreSpaceDet) then
            ! Dets in the core space are never removed from main list, so
            ! strictly do not require corrections.
            if (tThreshOccRDMDiag .and. (abs(contrib_sign(1)) <= ThreshOccRDM)) ScaleContribFac = 0.0_dp
        end if

        s_orbs = gtID(nI)

        i = 1

        do while (i <= nel)

            s = s_orbs(i)

            ind = SymLabelListInv_rot(s)

            if (size(contrib_sign) == 1) then
                final_contrib = contrib_sign**2 * RDMIters * ScaleContribFac
            else
                final_contrib = contrib_sign(1::2)*contrib_sign(2::2) * RDMIters * ScaleContribFac
            end if

            if (isDouble(nI, i)) then
                inc = 2

                final_contrib = 2.0_dp * final_contrib

            else
                inc = 1

            end if

            if (tAdaptiveShift .and. all(nI == projEDet(:, 1))) &
                final_contrib = final_contrib + RDMIters * ScaleContribFac * rdmCorrectionFactor
            do irdm = 1, size(one_rdms)
                one_rdms(irdm)%matrix(ind, ind) = one_rdms(irdm)%matrix(ind, ind) + final_contrib(irdm)
            end do

            i = i + inc

        end do

    end subroutine fill_diag_1rdm_guga

    subroutine send_proc_ex_djs()
        ! i also need a specific routint to send the excitations for
        ! the GUGA RDMs, otherwise this clutters up the det-based routines

        ! this routine is very similar to SendProcExcDjs in the rdm_explicit
        ! module. see there for comments
        character(*), parameter :: this_routine = "send_proc_ex_djs"

        integer :: i, j
        integer :: error, MaxSendIndex, MaxIndex
        integer(MPIArg) :: sendcounts(nProcessors), disps(nProcessors)
        integer(MPIArg) :: sing_recvcounts(nProcessors)
        integer(MPIArg) :: sing_recvdisps(nProcessors)
        integer(MPIArg) :: doub_recvcounts(nProcessors), doub_recvdisps(nProcessors)

        if (RDMExcitLevel == 1 .or. RDMExcitLevel == 3) then
            do i = 0, nProcessors - 1
                sendcounts(i + 1) = int(Sing_ExcList(i) - (nint(OneEl_Gap * i) + 1), MPIArg)
                disps(i + 1) = nint(OneEl_Gap * i, MPIArg)
            end do

            MaxSendIndex = Sing_ExcList(nProcessors - 1) - 1

            sing_recvcounts(1:nProcessors) = 0
            call MPIAlltoAll(sendcounts, 1, sing_recvcounts, 1, error)

            sing_recvdisps(1) = 0
            do i = 2, nProcessors
                sing_recvdisps(i) = sing_recvdisps(i - 1) + sing_recvcounts(i - 1)
            end do

            MaxIndex = sing_recvdisps(nProcessors) + sing_recvcounts(nProcessors)

            do i = 1, nProcessors
                sendcounts(i) = sendcounts(i) * (int(GugaBits%len_tot + 1, MPIArg))
                disps(i) = disps(i) * (int(GugaBits%len_tot + 1, MPIArg))
                sing_recvcounts(i) = sing_recvcounts(i) * (int(GugaBits%len_tot + 1, MPIArg))
                sing_recvdisps(i) = sing_recvdisps(i) * (int(GugaBits%len_tot + 1, MPIArg))
            end do

#ifdef USE_MPI
            call MPIAlltoAllv(Sing_ExcDjs(:, 1:MaxSendIndex), sendcounts, disps, &
                              Sing_ExcDjs2, sing_recvcounts, sing_recvdisps, error)
#else
            Sing_ExcDjs2(0:GugaBits%len_tot, 1:MaxIndex) = Sing_ExcDjs(0:GugaBits%len_tot, 1:MaxSendIndex)
#endif

            ! and also write a new routine for the search of occ. dets
            call singles_search_guga(sing_recvcounts, sing_recvdisps)
        end if

        if (RDMExcitLevel /= 1) then

            do i = 0, nProcessors - 1
                ! Sendcounts is the number of singly excited determinants we
                ! want to send for each processor (but goes from 1, not 0).
                sendcounts(i + 1) = int(Doub_ExcList(i) - (nint(TwoEl_Gap * i) + 1), MPIArg)

                ! disps is the first slot for each processor - 1.
                disps(i + 1) = nint(TwoEl_Gap * i, MPIArg)
            end do

            MaxSendIndex = Doub_ExcList(nProcessors - 1) - 1

            ! We now need to calculate the recvcounts and recvdisps -
            ! this is a job for AlltoAll
            doub_recvcounts(1:nProcessors) = 0
            call MPIAlltoAll(sendcounts, 1, doub_recvcounts, 1, error)

            ! We can now get recvdisps from recvcounts, since we want the data
            ! to be contiguous after the move.
            doub_recvdisps(1) = 0
            do i = 2, nProcessors
                doub_recvdisps(i) = doub_recvdisps(i - 1) + doub_recvcounts(i - 1)
            end do

            MaxIndex = doub_recvdisps(nProcessors) + doub_recvcounts(nProcessors)
            ! But the actual number of integers we need to send is the
            ! calculated values * NIfTot+1.
            do i = 1, nProcessors
                sendcounts(i) = sendcounts(i) * (int(GugaBits%len_tot + 1, MPIArg))
                disps(i) = disps(i) * (int(GugaBits%len_tot + 1, MPIArg))
                doub_recvcounts(i) = doub_recvcounts(i) * (int(GugaBits%len_tot + 1, MPIArg))
                doub_recvdisps(i) = doub_recvdisps(i) * (int(GugaBits%len_tot + 1, MPIArg))
            end do

            ! This is the main send of all the single excitations to the
            ! corresponding processors.
#ifdef USE_MPI
            call MPIAlltoAllv(Doub_ExcDjs(:, 1:MaxSendIndex), sendcounts, disps, &
                              Doub_ExcDjs2, doub_recvcounts, doub_recvdisps, error)
#else
            Doub_ExcDjs2(0:GugaBits%len_tot, 1:MaxIndex) = Doub_ExcDjs(0:GugaBits%len_tot, 1:MaxSendIndex)
#endif

            call doubles_search_guga(doub_recvcounts, doub_recvdisps)

        end if

    end subroutine send_proc_ex_djs

    subroutine doubles_search_guga(recvcounts, recvdisps)
        ! this is the adapted explicit doubles search for the
        ! GUGA-RDMs.
        integer(MPIArg), intent(in) :: recvcounts(nProcessors), recvdisps(nProcessors)
        character(*), parameter :: this_routine = "doubles_search_guga"

        integer :: i, NoDets, StartDets, nI(nel), nJ(nel), PartInd, FlagsDj
        integer :: j, a, b, c, d
        integer(int_rdm) :: rdm_ind
        integer(n_int) :: ilutJ(0:GugaBits%len_tot), ilutI(0:GugaBits%len_tot)
        real(dp) :: mat_ele, sign_i(lenof_sign), sign_j(lenof_sign)
        logical :: tDetFound

        do i = 1, nProcessors

            NoDets = recvcounts(i) / (GugaBits%len_tot + 1)
            StartDets = (recvdisps(i) / (GugaBits%len_tot + 1)) + 1

            if (NoDets > 1) then

                ilutI = Doub_ExcDjs2(:, StartDets)

                sign_i = extract_matrix_element(ilutI, 2)

                do j = StartDets + 1, (NoDets + StartDets - 1)

                    ilutJ = Doub_ExcDjs2(:, j)

                    call BinSearchParts_rdm(ilutJ, 1, int(TotWalkers, sizeof_int), &
                                            PartInd, tDetFound)

                    if (tDetFound) then

                        call extract_bit_rep(CurrentDets(:, PartInd), nJ, sign_j, FlagsDj)

                        mat_ele = extract_matrix_element(ilutJ, 1)
                        rdm_ind = extract_rdm_ind(ilutJ)

                        call extract_2_rdm_ind(rdm_ind, a, b, c, d)

                        call add_to_rdm_spawn_t(two_rdm_spawn, a, b, c, d, &
                                                sign_i * sign_j * mat_ele, .true.)
                    end if
                end do
            end if
        end do

    end subroutine doubles_search_guga

    subroutine singles_search_guga(recvcounts, recvdisps)
        ! make a tailored GUGA search for explicit single excitations

        ! this routine is very similar to Sing_SearchOccDets in the rdm_explicit
        ! module. see there for comments
        integer(MPIArg), intent(in) :: recvcounts(nProcessors), recvdisps(nProcessors)
        character(*), parameter :: this_routine = "singles_search_guga"

        integer :: i, NoDets, StartDets, nI(nel), nJ(nel), PartInd, FlagsDj
        integer :: j
        integer(int_rdm) :: rdm_ind
        integer(n_int) :: ilutJ(0:GugaBits%len_tot), ilutI(0:GugaBits%len_tot)
        real(dp) :: mat_ele, sign_i(lenof_sign), sign_j(lenof_sign)
        logical :: tDetFound

        do i = 1, nProcessors

            NoDets = recvcounts(i) / (GugaBits%len_tot + 1)
            StartDets = (recvdisps(i) / (GugaBits%len_tot + 1)) + 1

            if (NoDets > 1) then

                ! extract the determinant, the matrix element and the
                ! RDM index from the single excitation array

                ! the convention is: i encode the coupling coefficient
                ! in the x0 matrix element and the sign of Di in the
                ! x1 element. and the combined RDM index in the
                ! Delta B value position!
                ilutI = Sing_ExcDjs2(:, StartDets)

                sign_i = extract_matrix_element(ilutI, 2)

                do j = StartDets + 1, (NoDets + StartDets - 1)

                    ! apparently D_i is in the first spot and all
                    ! excitations come here afterwards..

                    ilutJ = Sing_ExcDjs2(:, j)

                    call BinSearchParts_rdm(ilutJ, 1, int(TotWalkers, sizeof_int), PartInd, tDetFound)

                    if (tDetFound) then

                        call extract_bit_rep(CurrentDets(:, PartInd), nJ, sign_j, FlagsDj)

                        mat_ele = extract_matrix_element(Sing_ExcDjs2(:, j), 1)
                        rdm_ind = extract_rdm_ind(Sing_ExcDjs2(:, j))

                        if (RDMExcitLevel == 1) then
                            call fill_sings_1rdm_guga(one_rdms, sign_i, sign_j, &
                                                      mat_ele, rdm_ind)

                        else if (RDMExcitLevel == 3) then
                            call fill_sings_2rdm_guga(two_rdm_spawn, ilutI, &
                                                      ilutJ, sign_i, sign_j, mat_ele, rdm_ind)

                        end if
                    end if
                end do
            end if
        end do

    end subroutine singles_search_guga

    subroutine fill_sings_1rdm_guga(one_rdms, sign_i, sign_j, mat_ele, rdm_ind)
        type(one_rdm_t), intent(inout) :: one_rdms(:)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        integer(int_rdm), intent(in) :: rdm_ind
        character(*), parameter :: this_routine = "fill_sings_1rdm_guga"

        integer :: i, a, ind_i, ind_a, irdm

        call extract_1_rdm_ind(rdm_ind, i, a)

        ind_i = SymLabelListInv_rot(i)
        ind_a = SymLabelListInv_rot(a)

        do irdm = 1, size(one_rdms)

            one_rdms(irdm)%matrix(ind_i, ind_a) = one_rdms(irdm)%matrix(ind_i, ind_a) &
                                                  + sign_i(irdm) * sign_j(irdm) * mat_ele

        end do

    end subroutine fill_sings_1rdm_guga

    subroutine fill_mixed_2rdm_guga(spawn, ilutI, ilutJ, sign_i, sign_j, mat_ele, &
                                    rdm_ind, excitInfo_opt)
        ! this is the routine where I test the filling of 2-RDM based on
        ! single excitations to mimic the workflow in the stochastic
        ! RDM sampling
        type(rdm_spawn_t), intent(inout) :: spawn
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        integer(int_rdm), intent(in) :: rdm_ind
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_opt
        character(*), parameter :: this_routine = "fill_mixed_2rdm_guga"

        type(ExcitationInformation_t) :: excitInfo
        integer :: ij, kl, i, j, k, l
        ! mimic the stochastic processes for the explicit excitation
        ! generation

        if (DetBitEQ(ilutI, ilutJ, nifd)) return

        ! for the explicit code, I first need to recalculate the type
        ! of excitation.. I am not sure how the indices will be encoded..
        if (present(excitInfo_opt)) then
            excitInfo = excitInfo_opt
        else
            call extract_2_rdm_ind(rdm_ind, i, j, k, l)

            excitInfo = excitationIdentifier(i, j, k, l)
        end if

        select case (excitInfo%typ)

        case (excit_type%fullstop_L_to_R, &
              excit_type%fullstop_R_to_L)

            call fill_mixed_end(spawn, ilutI, ilutJ, sign_i, sign_j, &
                                mat_ele, excitInfo)

        case (excit_type%fullstart_L_to_R, &
              excit_type%fullstart_R_to_L)

            call fill_mixed_start(spawn, ilutI, ilutJ, sign_i, sign_j, &
                                  mat_ele, excitInfo)

        case (excit_type%fullstart_stop_mixed)
            call fill_mixed_start_end(spawn, ilutI, ilutJ, sign_i, sign_j, &
                                      mat_ele, excitInfo)

        case default
            call print_excitInfo(excitInfo)
            call Stop_All(this_routine, "wrong excitation type here")

        end select

    end subroutine fill_mixed_2rdm_guga

    subroutine fill_mixed_end(spawn, ilutI, ilutJ, sign_i, sign_j, &
                              mat_ele, excitInfo)
        type(rdm_spawn_t), intent(inout) :: spawn
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        type(ExcitationInformation_t), intent(in) :: excitInfo
        character(*), parameter :: this_routine = "fill_mixed_end"

        integer :: step, sw, elecInd, holeInd, i, j
        integer :: step_i(nSpatOrbs), b_i(nSpatOrbs), int_occ(nSpatOrbs)
        real(dp) :: top_cont, tmp_mat, stay_mat, end_mat, real_b(nSpatOrbs)
        real(dp) :: occ_i(nSpatOrbs)
        logical :: above_flag

        ! do as much stuff as possible beforehand
        associate (st => excitInfo%fullStart, se => excitInfo%secondStart, &
                   en => excitInfo%fullEnd, typ => excitInfo%typ)

            if (typ == excit_type%fullstop_L_to_R) then
                elecInd = st
                holeInd = se
            else if (typ == excit_type%fullstop_R_to_L) then
                elecInd = se
                holeInd = st
            else
                call stop_all(this_routine, "should not be here!")
            end if

            sw = findLastSwitch(ilutI, ilutJ, se, en)

            call calc_csf_info(ilutI, step_i, b_i, occ_i)
            int_occ = int(occ_i)

            step = step_i(en)

            real_b = real(b_i, dp)

            if (en < nSpatOrbs) then
                select case (step)
                case (1)
                    if (isOne(ilutJ, en)) then
                        top_cont = -Root2 * sqrt((real_b(en) + 2.0_dp) / &
                                                 real_b(en))

                    else
                        top_cont = -Root2 / sqrt(real_b(en) * (real_b(en) + 2.0_dp))

                    end if
                case (2)
                    if (isOne(ilutJ, en)) then
                        top_cont = -Root2 / sqrt(real_b(en) * (real_b(en) + 2.0_dp))

                    else
                        top_cont = Root2 * sqrt(real_b(en) / &
                                                (real_b(en) + 2.0_dp))
                    end if

                case default
                    call stop_all(this_routine, "wrong stepvalues!")

                end select

                if (.not. near_zero(top_cont)) then

                    above_flag = .false.
                    tmp_mat = 1.0_dp

                    do i = en + 1, nSpatOrbs
                        if (int_occ(i) /= 1) cycle

                        ! then check if thats the last step
                        if (step_i(i) == 2 .and. b_i(i) == 0) then
                            above_flag = .true.
                        end if

                        ! should be able to do that without second loop too!
                        ! figure out!
                        step = step_i(i)

                        call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, real_b(i), &
                                                    1.0_dp, x1_element=stay_mat)

                        call getMixedFullStop(step, step, 0, real_b(i), &
                                              x1_element=end_mat)

                        ! this check should never be true, but just to be sure
                        if (near_zero(stay_mat)) above_flag = .true.

                        if (.not. near_zero(end_mat)) then

                            call add_to_rdm_spawn_t(spawn, holeInd, i, i, elecInd, &
                                                    top_cont * end_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)
                            call add_to_rdm_spawn_t(spawn, i, elecInd, holeInd, i, &
                                                    top_cont * end_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)

                        end if

                        if (above_flag) exit

                        ! otherwise update your running pgen and matrix element vars
                        tmp_mat = tmp_mat * stay_mat

                    end do

                    ! have to figure out what to do with this here:!
                end if
            end if

            if (sw < en) then

                step = step_i(en)

                ! inverse fullstop matrix element
                call getMixedFullStop(step, step, 0, real_b(en), x1_element=tmp_mat)

                if (.not. near_zero(tmp_mat)) then
                    tmp_mat = 1.0_dp / tmp_mat

                    ! have to change the switches before the first cycle:
                    ! but for cycling backwards, thats not so easy.. need todo

                    do i = en - 1, sw + 1, -1

                        if (int_occ(i) /= 1) cycle

                        step = step_i(i)
                        ! update inverse product
                        call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, real_b(i), &
                                                    1.0_dp, x1_element=stay_mat)

                        call getMixedFullStop(step, step, 0, real_b(i), x1_element=end_mat)

                        ! update matrix element
                        tmp_mat = tmp_mat / stay_mat

                        if (.not. near_zero(end_mat)) then

                            call add_to_rdm_spawn_t(spawn, holeInd, i, i, elecInd, &
                                                    end_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)
                            call add_to_rdm_spawn_t(spawn, i, elecInd, holeInd, i, &
                                                    end_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)

                        end if

                    end do

                    ! deal with switch specifically:

                    step = step_i(sw)

                    if (step == 1) then
                        ! then a -2 branch arrived!
                        call getDoubleMatrixElement(2, 1, -2, gen_type%L, gen_type%R, real_b(sw), &
                                                    1.0_dp, x1_element=stay_mat)

                        call getMixedFullStop(2, 1, -2, real_b(sw), x1_element=end_mat)

                    else
                        ! +2 branch arrived!

                        call getDoubleMatrixElement(1, 2, 2, gen_type%L, gen_type%R, real_b(sw), &
                                                    1.0_dp, x1_element=stay_mat)

                        call getMixedFullStop(1, 2, 2, real_b(sw), x1_element=end_mat)
                    end if

                    tmp_mat = tmp_mat * end_mat / stay_mat

                    call add_to_rdm_spawn_t(spawn, holeInd, sw, sw, elecInd, &
                                            tmp_mat * sign_i * sign_j * mat_ele, .true.)
                    call add_to_rdm_spawn_t(spawn, sw, elecInd, holeInd, sw, &
                                            tmp_mat * sign_i * sign_j * mat_ele, .true.)
                end if
            end if

        end associate

    end subroutine fill_mixed_end

    subroutine fill_mixed_start(spawn, ilutI, ilutJ, sign_i, sign_j, &
                                mat_ele, excitInfo)
        type(rdm_spawn_t), intent(inout) :: spawn
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        type(ExcitationInformation_t), intent(in) :: excitInfo
        character(*), parameter :: this_routine = "fill_mixed_start"

        integer :: sw, i, step, elecInd, holeInd
        integer :: step_i(nSpatOrbs), b_i(nSpatOrbs), int_occ(nSpatOrbs)
        real(dp) :: bot_cont, tmp_mat, stay_mat, start_mat, real_b(nSpatOrbs)
        real(dp) :: occ_i(nSpatOrbs)
        logical :: below_flag

        associate (st => excitInfo%fullStart, se => excitInfo%firstEnd, &
                   en => excitInfo%fullEnd, typ => excitInfo%typ)

            ! depending on the type of excitaiton, calculation of orbital pgens
            ! change
            if (typ == excit_type%fullstart_L_to_R) then
                elecInd = en
                holeInd = se
            else if (typ == excit_type%fullstart_R_to_L) then
                elecInd = se
                holeInd = en
            else
                call stop_all(this_routine, "should not be here!")
            end if

            sw = findFirstSwitch(ilutI, ilutJ, st, se)

            call calc_csf_info(ilutI, step_i, b_i, occ_i)
            real_b = real(b_i, dp)
            int_occ = int(occ_i)

            step = step_i(st)

            if (step == 1) then
                if (isOne(ilutJ, st)) then
                    bot_cont = Root2 * sqrt((real_b(st) - 1.0_dp) / &
                                            (real_b(st) + 1.0_dp))
                else
                    bot_cont = -sqrt(2.0_dp / ((real_b(st) - 1.0_dp) * &
                                               (real_b(st) + 1.0_dp)))
                end if
            else
                if (isOne(ilutJ, st)) then
                    bot_cont = -sqrt(2.0_dp / ((real_b(st) + 1.0_dp) * &
                                               (real_b(st) + 3.0_dp)))
                else
                    bot_cont = -Root2 * sqrt((real_b(st) + 3.0_dp) / &
                                             (real_b(st) + 1.0_dp))
                end if
            end if

            if (.not. near_zero(bot_cont)) then

                tmp_mat = 1.0_dp
                below_flag = .false.

                do i = st - 1, 1, -1
                    if (int_occ(i) /= 1) cycle

                    ! then check if thats the last stepvalue to consider
                    if (step_i(i) == 1 .and. b_i(i) == 1) then
                        below_flag = .true.
                    end if

                    ! then deal with the matrix element and branching probabilities
                    step = step_i(i)

                    ! get both start and staying matrix elements -> and update
                    ! matrix element contributions on the fly to avoid second loop!
                    call getDoubleMatrixElement(step, step, -1, gen_type%R, gen_type%L, real_b(i), &
                                                1.0_dp, x1_element=start_mat)

                    call getDoubleMatrixElement(step, step, 0, gen_type%R, gen_type%L, real_b(i), &
                                                1.0_dp, x1_element=stay_mat)

                    if (near_zero(stay_mat)) below_flag = .true.

                    if (.not. near_zero(start_mat)) then

                        call add_to_rdm_spawn_t(spawn, holeInd, i, i, elecInd, &
                                                bot_cont * start_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)
                        call add_to_rdm_spawn_t(spawn, i, elecInd, holeInd, i, &
                                                bot_cont * start_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)

                    end if

                    ! also update matrix element on the fly
                    tmp_mat = stay_mat * tmp_mat

                end do

            end if

            step = step_i(st)

            ! calculate the necarry values needed to formulate everything in terms
            ! of the already calculated quantities:
            call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, real_b(st), &
                                        1.0_dp, x1_element=tmp_mat)

            ! and calc. x1^-1
            ! keep tempWweight as the running matrix element which gets updated
            ! every iteration
            if (.not. near_zero(tmp_mat)) then
                tmp_mat = 1.0_dp / tmp_mat

                do i = st + 1, sw - 1
                    ! the good thing here is, i do not need to loop a second time,
                    ! since i can recalc. the matrix elements and pgens on-the fly
                    ! here the matrix elements should not be 0 or otherwise the
                    ! excitation wouldnt have happended anyways
                    if (int_occ(i) /= 1) cycle

                    step = step_i(i)

                    ! update inverse product
                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, real_b(i), &
                                                1.0_dp, x1_element=stay_mat)

                    tmp_mat = tmp_mat / stay_mat

                    ! and also get starting contribution
                    call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, real_b(i), &
                                                1.0_dp, x1_element=start_mat)

                    ! because the rest of the matrix element is still the same in
                    ! both cases...
                    if (.not. near_zero(start_mat)) then
                        call add_to_rdm_spawn_t(spawn, holeInd, i, i, elecInd, &
                                                start_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)
                        call add_to_rdm_spawn_t(spawn, i, elecInd, holeInd, i, &
                                                start_mat * tmp_mat * sign_i * sign_j * mat_ele, .true.)
                    end if

                end do

                ! handle switch seperately (but only if switch > start)
                if (sw > st) then

                    step = step_i(sw)

                    ! on the switch the original probability is:
                    if (step == 1) then

                        call getDoubleMatrixElement(2, 1, 0, gen_type%L, gen_type%R, real_b(sw), &
                                                    1.0_dp, x1_element=stay_mat)

                        call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, real_b(sw), &
                                                    1.0_dp, x1_element=start_mat)

                    else

                        call getDoubleMatrixElement(1, 2, 0, gen_type%L, gen_type%R, real_b(sw), &
                                                    1.0_dp, x1_element=stay_mat)

                        call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, real_b(sw), &
                                                    1.0_dp, x1_element=start_mat)

                    end if

                    ! update inverse product
                    ! and also get starting contribution
                    tmp_mat = tmp_mat * start_mat / stay_mat

                    ! because the rest of the matrix element is still the same in
                    ! both cases...

                    call add_to_rdm_spawn_t(spawn, holeInd, sw, sw, elecInd, &
                                            tmp_mat * sign_i * sign_j * mat_ele, .true.)
                    call add_to_rdm_spawn_t(spawn, sw, elecInd, holeInd, sw, &
                                            tmp_mat * sign_i * sign_j * mat_ele, .true.)

                end if
            end if

        end associate

    end subroutine fill_mixed_start

    subroutine fill_mixed_start_end(spawn, ilutI, ilutJ, sign_i, sign_j, &
                                    mat_ele, excitInfo)
        type(rdm_spawn_t), intent(inout) :: spawn
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        type(ExcitationInformation_t), intent(in) :: excitInfo
        character(*), parameter :: this_routine = "fill_mixed_start_end"

        integer :: first, last, deltaB(nSpatOrbs), i, j, k, step1, step2
        integer :: step_i(nSpatOrbs), b_i(nSpatOrbs), int_occ(nSpatOrbs)
        integer :: step_j(nSpatOrbs), b_j(nSpatOrbs)
        logical :: above_flag, below_flag
        real(dp) :: inter, tempWeight_1, real_b(nSpatOrbs), tempWeight
        real(dp) :: occ_i(nSpatOrbs), occ_j(nSpatOrbs)

        unused_var(mat_ele)

        first = findFirstSwitch(ilutI, ilutJ, excitInfo%fullStart, excitInfo%fullEnd)
        last = findLastSwitch(ilutI, ilutJ, first, excitInfo%fullEnd)

        call calc_csf_info(ilutI, step_i, b_i, occ_i)
        call calc_csf_info(ilutJ, step_j, b_j, occ_j)
        real_b = real(b_i, dp)
        int_occ = int(occ_i)

        below_flag = .false.
        above_flag = .false.

        deltaB = b_i - b_j

        inter = 1.0_dp

        ! calculate the always involved intermediate matrix element from
        ! first switch to last switch
        do i = first + 1, last - 1
            if (int_occ(i) /= 1) cycle

            step1 = step_i(i)
            step2 = step_j(i)
            call getDoubleMatrixElement(step2, step1, deltaB(i - 1), gen_type%L, gen_type%R, &
                                        real_b(i), 1.0_dp, x1_element=tempWeight)

            inter = inter * tempWeight
        end do

        do j = last, nSpatOrbs
            if (int_occ(j) /= 1) cycle

            ! calculate the remaining switches once for each (j) but do it
            ! for the worst case until i = 1

            ! check if this is the last end needed to consider
            if (step_i(j) == 2 .and. b_i(j) == 0) then
                above_flag = .true.
            end if

            ! i have to reset the below flag each iteration of j..
            below_flag = .false.

            do i = first, 1, -1
                if (int_occ(i) /= 1) cycle

                if (below_flag) exit

                ! if the bottom stepvector d = 1 and b = 1 there is no
                ! additional contribution from below, since the x1 matrix
                ! element is 0
                ! same if d = 2 and b = 0 for fullstop stepvector
                if (step_i(i) == 1 .and. b_i(i) == 1) then
                    below_flag = .true.
                end if

                ! get the starting matrix element
                step1 = step_i(i)
                step2 = step_j(i)
                call getDoubleMatrixElement(step2, step1, -1, gen_type%L, gen_type%R, &
                                            real_b(i), 1.0_dp, x1_element=tempWeight)

                ! loop over excitation range
                ! distinguish between different regimes
                ! if i do it until switch - 1 -> i know that dB = 0 and
                ! the 2 stepvalues are always the same..
                do k = i + 1, first - 1
                    if (int_occ(k) /= 1) cycle

                    step1 = step_i(k)
                    ! only 0 branch here
                    call getDoubleMatrixElement(step1, step1, 0, gen_type%L, gen_type%R, &
                                                real_b(k), 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = tempWeight * tempWeight_1

                end do

                ! then do first switch site seperately, if (i) is not first
                ! and what if (i) is first??
                if (i /= first) then
                    step1 = step_i(first)

                    if (step1 == 1) then
                        ! i know that step2 = 2
                        call getDoubleMatrixElement(2, 1, 0, gen_type%L, gen_type%R, &
                                                    real_b(first), 1.0_dp, x1_element=tempWeight_1)

                    else
                        ! i know that step2 = 1
                        call getDoubleMatrixElement(1, 2, 0, gen_type%L, gen_type%R, real_b(first), &
                                                    1.0_dp, x1_element=tempWeight_1)

                    end if

                    tempWeight = tempWeight * tempWeight_1

                end if

                ! more efficient to do "last" step seperately, since i have to
                ! check deltaB value and also have to consider matrix element
                ! but only of (j) is not last or otherwise already dealt with
                if (j /= last) then

                    if (step_i(last) == 1) then
                        ! then i know step2 = 2 & dB = -2!
                        call getDoubleMatrixElement(2, 1, -2, gen_type%L, gen_type%R, &
                                                    real_b(last), 1.0_dp, x1_element=tempWeight_1)
                    else
                        ! i know step2 == 1 and dB = +2
                        call getDoubleMatrixElement(1, 2, +2, gen_type%L, gen_type%R, &
                                                    real_b(last), 1.0_dp, x1_element=tempWeight_1)

                    end if

                    tempWeight = tempWeight * tempWeight_1
                end if

                ! then do remaining top range, where i know stepvalues are
                ! the same again and dB = 0 always!
                do k = last + 1, j - 1
                    if (int_occ(k) /= 1) cycle

                    step1 = step_i(k)
                    ! only 0 branch here
                    call getDoubleMatrixElement(step1, step1, 0, gen_type%L, gen_type%R, &
                                                real_b(k), 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = tempWeight * tempWeight_1

                end do

                ! and handle fullend
                ! and then do the the end value at j
                step1 = step_i(j)
                step2 = step_j(j)
                call getMixedFullStop(step2, step1, deltaB(j - 1), real_b(j), &
                                      x1_element=tempWeight_1)

                call add_to_rdm_spawn_t(spawn, i, j, j, i, &
                                        tempWeight * tempWeight_1 * inter * sign_i * sign_j, .true.)
                call add_to_rdm_spawn_t(spawn, j, i, i, j, &
                                        tempWeight * tempWeight_1 * inter * sign_i * sign_j, .true.)

                ! check if i deal with that correctly...
                if (below_flag) exit
            end do
            ! todo: i cant use tthat like that.. or else some combinations
            ! of i and j get left out! i have to reinit it somehow..
            ! not yet sure how..
            if (above_flag) exit
        end do

    end subroutine fill_mixed_start_end

    subroutine fill_sings_2rdm_guga(spawn, ilutI, ilutJ, sign_i, sign_j, mat_ele, rdm_ind)
        ! this is the routine where I test the filling of 2-RDM based on
        ! single excitations to mimic the workflow in the stochastic
        ! RDM sampling
        type(rdm_spawn_t), intent(inout) :: spawn
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot), ilutJ(0:GugaBits%len_tot)
        real(dp), intent(in) :: sign_i(:), sign_j(:), mat_ele
        integer(int_rdm), intent(in) :: rdm_ind
        character(*), parameter :: this_routine = "fill_sings_2rdm_guga"

        integer :: i, a, n, st, en, step_i(nSpatOrbs), step_j(nSpatOrbs), &
                   delta_b(nSpatOrbs), gen, d_i, d_j, delta, &
                   b_i(nSpatOrbs), b_j(nSpatOrbs)
        real(dp) :: occ_n, occ_i(nSpatOrbs), occ_j(nSpatOrbs)
        real(dp) :: botCont, topCont, tempWeight, prod, StartCont, EndCont
        integer :: iO, jO, step

        ! here I have to fill W + R/L, single overlap RR/LL
        ! and full-start/stop RL with no change in the double overlap region
        ! i also have to correct the coupling coefficient here effifiently

        ! for this it would be best to have both CSFs I and J involved.

        ! and i have to figure out all correct index combinations here
        ! for all the entries the singles contribute to.

        ! ok this is the final routine i need to consider i guess..
        ! or hope.. I need the matrix elements and indices, to which a
        ! certain type of single excitations contributes to..
        ! I need to effectively recalculate the correct matrix element
        ! and store it in the correct RDM place.

        ! essentially I need to loop over the contracted index and
        ! correctly take the coupling coefficient into account.
        call extract_1_rdm_ind(rdm_ind, i, a)

        st = min(i, a)
        en = max(i, a)

        ! do i have access to current_stepvector here? I think so..
        ! no I dont! or just to be save

        ! this is essentially a mimic of the function
        ! calc_integral_contribution_single in guga_excitations

        ! but wait a minute..
        ! in the stochastic excitation generation, I calculate the
        ! full Hamiltonian matrix element..
        ! I do not have access to the coupling coefficient for
        ! i, a anymore.. damn..
        ! and in general, I always calculate the full matrix element..
        ! i have to take this into account!
        ! or change the stochastic excitation generation, to also yield
        ! this information..
        ! or "just" recalculate everything..
        ! i could use the x1 element storage for this quantity..
        ! this would simplify things!

        ! for simplicity it is best to have these quantities:
        call calc_csf_info(ilutI, step_i, b_i, occ_i)
        call calc_csf_info(ilutJ, step_j, b_j, occ_j)

        delta_b = b_i - b_j

        ! calculate the bottom contribution depending on the excited stepvalue
        select case (step_i(st))
        case (0)
            ! this implicates a raising st:
            if (isOne(ilutJ, st)) then
                call getDoubleMatrixElement(1, 0, 0, gen_type%L, gen_type%R, real(b_i(st), dp), &
                                            1.0_dp, x1_element=botCont)

            else
                call getDoubleMatrixElement(2, 0, 0, gen_type%L, gen_type%R, real(b_i(st), dp), &
                                            1.0_dp, x1_element=botCont)
            end if

            StartCont = 0.0_dp
            gen = gen_type%R

        case (3)
            ! implies lowering st
            if (isOne(ilutJ, st)) then
                ! need tA(0,2)
                botCont = funA_0_2overR2(real(b_i(st), dp))

            else
                ! need -tA(2,0)
                botCont = minFunA_2_0_overR2(real(b_i(st), dp))
            end if

            StartCont = 1.0_dp
            gen = gen_type%L

        case (1)
            botCont = funA_m1_1_overR2(real(b_i(st), dp))
            ! check which generator
            if (isZero(ilutJ, st)) then
                botCont = -botCont

                StartCont = 0.0_dp
                gen = gen_type%L
            else
                StartCont = 1.0_dp
                gen = gen_type%R
            end if

        case (2)
            botCont = funA_3_1_overR2(real(b_i(st), dp))
            if (isThree(ilutJ, st)) then
                botCont = -botCont

                StartCont = 1.0_dp
                gen = gen_type%R
            else
                StartCont = 0.0_dp
                gen = gen_type%L
            end if

        end select

        ! do top contribution also already

        select case (step_i(en))
        case (0)
            if (isOne(ilutJ, en)) then
                topCont = funA_2_0_overR2(real(b_i(en), dp))
            else
                topCont = minFunA_0_2_overR2(real(b_i(en), dp))
            end if

            EndCont = 0.0_dp

        case (3)
            if (isOne(ilutJ, en)) then
                topCont = minFunA_2_0_overR2(real(b_i(en), dp))
            else
                topCont = funA_0_2overR2(real(b_i(en), dp))
            end if

            EndCont = 1.0_dp

        case (1)
            topCont = funA_2_0_overR2(real(b_i(en), dp))
            if (isThree(ilutJ, en)) then
                topCont = -topCont

                EndCont = 1.0_dp
            else
                EndCont = 0.0_dp
            end if

        case (2)
            topCont = funA_0_2overR2(real(b_i(en), dp))
            if (isZero(ilutJ, en)) then
                topCont = -topCont

                EndCont = 0.0_dp
            else
                EndCont = 1.0_dp
            end if

        end select

        ! depending on i and j calulate the corresponding single and double
        ! integral weights and check if they are non-zero...
        ! gets quite involved... :( need to loop over all orbitals
        ! have to reset prod inside the loop each time!

        do iO = 1, st - 1
            ! no contribution if not occupied.
            if (step_i(iO) == 0) cycle
            ! else it gets a contrbution weighted with orbital occupation
            ! first easy part:

            ! i think also singles I should store symmetrically! as
            ! the conjugated element will of course never be stored!
            ! W + R/L contribution
            call add_to_rdm_spawn_t(spawn, iO, iO, i, a, &
                                    occ_i(iO) * sign_i * sign_j * mat_ele, .true.)

            call add_to_rdm_spawn_t(spawn, i, a, iO, iO, &
                                    occ_i(iO) * sign_i * sign_j * mat_ele, .true.)

            ! exchange contribution:
            ! wtf is this b_i == 0 check? that does not make any sense..
            if (step_i(iO) == 3 .or. ((occ_i(iO) .isclose.1.0_dp) .and. b_i(iO) == 0)) then
                ! then it is easy:
                ! just figure out correct indices
                call add_to_rdm_spawn_t(spawn, i, iO, iO, a, &
                                        -occ_i(iO) / 2.0 * sign_i * sign_j * mat_ele, .true.)
                call add_to_rdm_spawn_t(spawn, iO, a, i, iO, &
                                        -occ_i(iO) / 2.0 * sign_i * sign_j * mat_ele, .true.)

            else
                ! otherwise i have to recalc the x1 element

                step = step_i(iO)
                call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, real(b_i(iO), dp), &
                                            1.0_dp, x1_element=prod)

                ! and then do the remaining:
                do jO = iO + 1, st - 1
                    ! need the stepvalue entries to correctly access the mixed
                    ! generator matrix elements
                    step = step_i(jO)
                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, &
                                                real(b_i(jO), dp), 1.0_dp, x1_element=tempWeight)

                    prod = prod * tempWeight
                end do
                prod = prod * botCont

                call add_to_rdm_spawn_t(spawn, i, iO, iO, a, &
                                        (prod - occ_i(iO) / 2.0) * sign_i * sign_j * mat_ele, .true.)
                call add_to_rdm_spawn_t(spawn, iO, a, i, iO, &
                                        (prod - occ_i(iO) / 2.0) * sign_i * sign_j * mat_ele, .true.)
            end if
        end do

        ! start segment: only W + R/L
        ! but this depends on the type of excitation here.. or?
        call add_to_rdm_spawn_t(spawn, st, st, i, a, &
                                StartCont * sign_i * sign_j * mat_ele, .true.)
        call add_to_rdm_spawn_t(spawn, i, a, st, st, &
                                StartCont * sign_i * sign_j * mat_ele, .true.)

        ! loop over excitation range:

        do iO = st + 1, en - 1

            ! W + R/L
            call add_to_rdm_spawn_t(spawn, iO, iO, i, a, &
                                    occ_i(iO) * sign_i * sign_j * mat_ele, .true.)
            call add_to_rdm_spawn_t(spawn, i, a, iO, iO, &
                                    occ_i(iO) * sign_i * sign_j * mat_ele, .true.)

            ! exchange type:
            ! oh damn I need the delta-B value here..
            d_i = step_i(iO)
            d_j = step_j(iO)
            delta = delta_b(iO - 1)

            prod = getDoubleContribution(d_j, d_i, delta, gen, real(b_i(iO), dp))

            call add_to_rdm_spawn_t(spawn, i, iO, iO, a, &
                                    prod * sign_i * sign_j * mat_ele, .true.)
            call add_to_rdm_spawn_t(spawn, iO, a, i, iO, &
                                    prod * sign_i * sign_j * mat_ele, .true.)
        end do

        ! end contribution
        call add_to_rdm_spawn_t(spawn, en, en, i, a, &
                                EndCont * sign_i * sign_j * mat_ele, .true.)
        call add_to_rdm_spawn_t(spawn, i, a, en, en, &
                                EndCont * sign_i * sign_j * mat_ele, .true.)

        ! loop above:
        do iO = en + 1, nSpatOrbs

            if (step_i(iO) == 0) cycle

            ! W + R/L contribution
            call add_to_rdm_spawn_t(spawn, iO, iO, i, a, &
                                    occ_i(iO) * sign_i * sign_j * mat_ele, .true.)
            call add_to_rdm_spawn_t(spawn, i, a, iO, iO, &
                                    occ_i(iO) * sign_i * sign_j * mat_ele, .true.)

            ! exchange
            if (step_i(iO) == 3 .or. (b_i(iO) == 1 .and. step_i(iO) == 1)) then
                ! only x0 contribution
                ! then it is easy:
                ! just figure out correct indices
                call add_to_rdm_spawn_t(spawn, i, iO, iO, a, &
                                        -occ_i(iO) / 2.0 * sign_i * sign_j * mat_ele, .true.)
                call add_to_rdm_spawn_t(spawn, iO, a, i, iO, &
                                        -occ_i(iO) / 2.0 * sign_i * sign_j * mat_ele, .true.)

            else

                prod = 1.0_dp

                do jO = en + 1, iO - 1

                    step = step_i(jO)

                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, real(b_i(jO), dp), &
                                                1.0_dp, x1_element=tempWeight)

                    prod = prod * tempWeight

                end do

                step = step_i(iO)

                call getMixedFullStop(step, step, 0, real(b_i(iO), dp), &
                                      x1_element=tempWeight)

                prod = prod * tempWeight

                prod = prod * topCont

                call add_to_rdm_spawn_t(spawn, i, iO, iO, a, &
                                        (prod - occ_i(iO) / 2.0) * sign_i * sign_j * mat_ele, .true.)
                call add_to_rdm_spawn_t(spawn, iO, a, i, iO, &
                                        (prod - occ_i(iO) / 2.0) * sign_i * sign_j * mat_ele, .true.)
            end if
        end do

    end subroutine fill_sings_2rdm_guga

    subroutine assign_excits_to_proc_guga(n_tot, excits, excit_lvl)
        integer, intent(in) :: n_tot, excit_lvl
        integer(n_int), intent(in), allocatable :: excits(:, :)
        character(*), parameter :: this_routine = "assign_excits_to_proc_guga"

        integer :: i, proc, nJ(nel)
        integer(n_int) :: ilut(0:niftot)

        ASSERT(excit_lvl == 1 .or. excit_lvl == 2)

        if (excit_lvl == 1) then
            do i = 1, n_tot

                call convert_ilut_toNECI(excits(:, i), ilut)

                call decode_bit_det(nJ, ilut)

                proc = DetermineDetNode(nel, nJ, 0)

                Sing_ExcDjs(:, Sing_ExcList(proc)) = excits(:, i)
                Sing_ExcList(proc) = Sing_ExcList(proc) + 1

#ifdef DEBUG_
                if (Sing_ExcList(Proc) > nint(OneEl_Gap * (Proc + 1))) then
                    write(6, *) 'Proc', Proc
                    write(6, *) 'Sing_ExcList', Sing_ExcList
                    write(6, *) 'No. spaces for each proc', nint(OneEl_Gap)
                    call Stop_All('GenExcDjs', 'Too many excitations for space available.')
                end if
#endif
            end do

        else if (excit_lvl == 2) then
            do i = 1, n_tot

                call convert_ilut_toNECI(excits(:, i), ilut)

                call decode_bit_det(nJ, ilut)
                proc = DetermineDetNode(nel, nJ, 0)

                Doub_ExcDjs(:, Doub_ExcList(proc)) = excits(:, i)
                Doub_ExcList(proc) = Doub_ExcList(proc) + 1

#ifdef DEBUG_
                if (Doub_ExcList(Proc) > nint(TwoEl_Gap * (Proc + 1))) then
                    write(6, *) 'Proc', Proc
                    write(6, *) 'Doub_ExcList', Doub_ExcList
                    write(6, *) 'No. spaces for each proc', nint(TwoEl_Gap)
                    call Stop_All('GenExcDjs', 'Too many excitations for space available.')
                end if
#endif
            end do
        end if

    end subroutine assign_excits_to_proc_guga

    subroutine calc_explicit_diag_2_rdm_guga(ilut, n_tot, excitations)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        character(*), parameter :: this_routine = "calc_explicit_diag_2_rdm_guga"

        integer :: i, j, k, l, nMax, ierr, n, n_excits, jl, ik
        integer(n_int), allocatable :: temp_excits(:, :), tmp_all_excits(:, :)
        integer(int_rdm) :: ijkl

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**3 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:GugaBits%len_tot, nMax), stat=ierr)
        call LogMemAlloc('tmp_all_excits', (GugaBits%len_tot + 1) * nMax, 8, this_routine, tag_tmp_excits)

        n_tot = 0

        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs

                if (i == j) cycle

                call calc_all_excits_guga_rdm_doubles(ilut, i, j, j, i, &
                                                      temp_excits, n_excits)
#ifdef DEBUG_
                do n = 1, n_excits
                    ASSERT(isProperCSF_ilut(temp_excits(:, n), .true.))
                end do
#endif

                if (n_excits > 1) then
                    call add_guga_lists_rdm(n_tot, n_excits - 1, &
                                            tmp_all_excits, temp_excits(:, 2:))
                end if

                deallocate(temp_excits)
            end do
        end do

        j = 1
        do i = 1, n_tot
            if (near_zero(extract_matrix_element(tmp_all_excits(:, i), 1))) cycle

            tmp_all_excits(:, j) = tmp_all_excits(:, i)

            j = j + 1

        end do

        n_tot = j - 1

        allocate(excitations(0:GugaBits%len_tot, n_tot), &
                  source=tmp_all_excits(0:GugaBits%len_tot, 1:n_tot), stat=ierr)
        ! hm to log that does not make so much sense.. since it gets called
        ! more than once and is only a temporary array..
        call LogMemAlloc('excitations', n_tot, 8, this_routine, tag_excitations)

        call sort(excitations(:, 1:n_tot), ilut_lt, ilut_gt)

        deallocate(tmp_all_excits)
        call LogMemDealloc(this_routine, tag_tmp_excits)

    end subroutine calc_explicit_diag_2_rdm_guga

    subroutine calc_explicit_2_rdm_guga(ilut, n_tot, excitations)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        character(*), parameter :: this_routine = "calc_explicit_2_rdm_guga"

        integer :: i, j, k, l, nMax, ierr, n, n_excits, jl, ik
        integer(n_int), allocatable :: temp_excits(:, :), tmp_all_excits(:, :)
        integer(int_rdm) :: ijkl

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**3 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:GugaBits%len_tot, nMax), stat=ierr)
        call LogMemAlloc('tmp_all_excits', (GugaBits%len_tot + 1) * nMax, 8, this_routine, tag_tmp_excits)

        n_tot = 0
        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs
                do k = 1, nSpatOrbs
                    do l = 1, nSpatOrbs

                        ! pure diagonal contribution:
                        if (i == j .and. k == l) cycle

                        call calc_all_excits_guga_rdm_doubles(ilut, i, j, k, l, &
                                                              temp_excits, n_excits)

#ifdef DEBUG_
                        do n = 1, n_excits
                            if (.not. isProperCSF_ilut(temp_excits(:, n), .true.)) then
                                print *, "===="
                                call write_det_guga(6, ilut, .true.)
                                call write_det_guga(6, temp_excits(:, n), .true.)
                                print *, i, j, k, l
                            end if
                            ASSERT(isProperCSF_ilut(temp_excits(:, n), .true.))
                        end do
#endif
                        if (i == l .and. j == k) then
                            ! exclude the diagonal exchange here,
                            ! as it is already accounted for in the
                            ! diagonal contribution routine
#ifdef DEBUG_
                            if (n_excits > 0) then
                                if (.not. DetBitEQ(ilut, temp_excits(:, 1), nifd)) then
                                    print *, "not equal!"
                                end if
                            end if
#endif

                            if (n_excits > 1) then
                                call add_guga_lists_rdm(n_tot, n_excits - 1, &
                                                        tmp_all_excits, temp_excits(:, 2:))
                            end if
                        else
                            if (n_excits > 0) then
                                call add_guga_lists_rdm(n_tot, n_excits, tmp_all_excits, temp_excits)
                            end if
                        end if

                        deallocate(temp_excits)

                    end do
                end do
            end do
        end do

        j = 1
        do i = 1, n_tot
            if (near_zero(extract_matrix_element(tmp_all_excits(:, i), 1))) cycle

            tmp_all_excits(:, j) = tmp_all_excits(:, i)

            j = j + 1

        end do

        n_tot = j - 1

        allocate(excitations(0:GugaBits%len_tot, n_tot), &
                  source=tmp_all_excits(0:GugaBits%len_tot, 1:n_tot), stat=ierr)

        ! hm to log that does not make so much sense.. since it gets called
        ! more than once and is only a temporary array..
        call LogMemAlloc('excitations', n_tot, 8, this_routine, tag_excitations)

        deallocate(tmp_all_excits)
        call LogMemDealloc(this_routine, tag_tmp_excits)

    end subroutine calc_explicit_2_rdm_guga

    subroutine calc_explicit_1_rdm_guga(ilut, n_tot, excitations)
        ! routine to calculate the one-RDM explicitly in the GUGA formalism
        ! the one RDM is given by rho_ij = <Psi|E_ij|Psi>
        ! so, similar to the actHamiltonian routine I need to loop over all
        ! the orbital indices (i,j) and calculate the action of a given E_ij
        ! on the sampled wavefunction |Psi> and the resulting overlap with <Psi|
        ! make this routine similar to GenExcDjs() in rdm_explicit.F90
        ! to insert it there to calculate the GUGA RDMs in this case
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(out) :: n_tot
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        character(*), parameter :: this_routine = "calc_explicit_1_rdm_guga"

        integer :: i, j, nMax, ierr, n, n_excits
        integer(n_int), allocatable :: temp_excits(:, :), tmp_all_excits(:, :)

        call init_csf_information(ilut)

        nMax = 6 + 4 * (nSpatOrbs)**2 * (count_open_orbs(ilut) + 1)
        allocate(tmp_all_excits(0:GugaBits%len_tot, nMax), stat=ierr)
        call LogMemAlloc('tmp_all_excits', (GugaBits%len_tot + 1) * nMax, 8, this_routine, tag_tmp_excits)

        n_tot = 0

        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs
                if (i == j) cycle

                call calc_all_excits_guga_rdm_singles(ilut, i, j, temp_excits, &
                                                      n_excits)

#ifdef DEBUG_
                do n = 1, n_excits
                    ASSERT(isProperCSF_ilut(temp_excits(:, n), .true.))
                end do
#endif

                if (n_excits > 0) then
                    call add_guga_lists_rdm(n_tot, n_excits, tmp_all_excits, temp_excits)
                end if

                deallocate(temp_excits)

            end do
        end do

        j = 1
        do i = 1, n_tot
            if (near_zero(extract_matrix_element(tmp_all_excits(:, i), 1))) cycle

            tmp_all_excits(:, j) = tmp_all_excits(:, i)

            j = j + 1

        end do

        n_tot = j - 1

        allocate(excitations(0:GugaBits%len_tot, n_tot), &
                  source=tmp_all_excits(0:GugaBits%len_tot, 1:n_tot), stat=ierr)
        ! hm to log that does not make so much sense.. since it gets called
        ! more than once and is only a temporary array..
        call LogMemAlloc('excitations', n_tot, 8, this_routine, tag_excitations)

        call sort(excitations(:, 1:n_tot), ilut_lt, ilut_gt)

        deallocate(tmp_all_excits)
        call LogMemDealloc(this_routine, tag_tmp_excits)

    end subroutine calc_explicit_1_rdm_guga

    subroutine calc_all_excits_guga_rdm_doubles(ilut, i, j, k, l, excits, n_excits)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: i, j, k, l
        integer(n_int), intent(out), allocatable :: excits(:, :)
        integer, intent(out) :: n_excits
        character(*), parameter :: this_routine = "calc_all_excits_guga_rdm_doubles"

        real(dp) :: umat, posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer :: ierr, n, exlevel
        type(ExcitationInformation_t) :: excitInfo
        logical :: compFlag
        integer(n_int) :: tmp_ilut(0:niftot)
        integer(int_rdm) :: ijkl

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        ! if called with k = l = 0 -> call single version of function
        if (k == 0 .and. l == 0) then
            call calc_all_excits_guga_rdm_singles(ilut, i, j, excits, n_excits)
            return
        end if

        ASSERT(k > 0 .and. k <= nSpatOrbs)
        ASSERT(l > 0 .and. l <= nSpatOrbs)

        ! default:
        n_excits = 0

        ! otherwise get excitation information
        excitInfo = excitationIdentifier(i, j, k, l)

        ! screw it. for now write a function which checks if indices and ilut
        ! are compatible, and not initiate the excitation right away
        ! but check cases again
        ! in checkCompatibility the number of switches is already
        ! calulated to check if the probabilistic weights fit... maybe but
        ! that out and reuse.. to not waste any effort.
        call checkCompatibility(ilut, excitInfo, compFlag, posSwitches, negSwitches)

        ! for mixed type full starts and/or full stops i have to consider
        ! the possible diagonal/single excitations here!
        ! although I think the full-stops are already correctly taken into
        ! account..
        ! and also the the full-starts are maybe correct already..
        ! so it was just the full-start into full-stop mixed!
        if (.not. compFlag) then
            allocate(excits(0, 0), stat=ierr)
            return
        end if

        select case (excitInfo%typ)

        case (excit_type%raising, &
              excit_type%lowering, &
              excit_type%single_overlap_lowering, &
              excit_type%single_overlap_raising)

            ! in the case of mimicking stochasitic
            ! excitation generation, we should abort here for
            ! these type of excitations!
            allocate(excits(0, 0), stat=ierr)
            n_excits = 0
            return

        end select

        select case (excitInfo%typ)
        case (excit_type%single)
            ! shouldnt be here.. onyl single excits and full weight gens
            allocate(excits(0, 0), stat=ierr)
            return

        case (excit_type%raising) ! weight + raising gen.
            ! can be treated almost like a single excitation
            ! essentially the same, except if d(w) == 3 in the excitaton regime

            call calcDoubleExcitation_withWeight(ilut, excitInfo, excits, &
                                                 n_excits, posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%lowering) ! weight + lowering gen
            call calcDoubleExcitation_withWeight(ilut, excitInfo, excits, &
                                                 n_excits, posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%non_overlap) ! non overlap
            call calcNonOverlapDouble(ilut, excitInfo, excits, n_excits, &
                                      posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%single_overlap_lowering) ! single overlap two lowering
            ! how can i efficiently adress that?
            ! can i write that efficiently in one function or do i need more?
            ! probably need more... i already determined
            call calcSingleOverlapLowering(ilut, excitInfo, excits, n_excits, &
                                           posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%single_overlap_raising) ! single overlap raising
            call calcSingleOverlapRaising(ilut, excitInfo, excits, n_excits, &
                                          posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%single_overlap_L_to_R) ! single overlap lowering into raising
            call calcSingleOverlapMixed(ilut, excitInfo, excits, n_excits, &
                                        posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%single_overlap_R_to_L) ! single overlap raising into lowering
            call calcSingleOverlapMixed(ilut, excitInfo, excits, n_excits, &
                                        posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_lowering) ! normal double overlap two lowering
            call calcDoubleLowering(ilut, excitInfo, excits, n_excits, &
                                    posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_raising) ! normal double overlap two raising
            call calcDoubleRaising(ilut, excitInfo, excits, n_excits, &
                                   posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_L_to_R_to_L) ! lowering into raising into lowering
            call calcDoubleRaising(ilut, excitInfo, excits, n_excits, &
                                   posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_R_to_L_to_R) ! raising into lowering into raising
            call calcDoubleLowering(ilut, excitInfo, excits, n_excits, &
                                    posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_L_to_R) ! lowering into raising double
            call calcDoubleL2R(ilut, excitInfo, excits, n_excits, &
                               posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_R_to_L) ! raising into lowering double
            call calcDoubleR2L(ilut, excitInfo, excits, n_excits, &
                               posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstop_lowering) ! full stop 2 lowering
            ! can i write a function for both alike generator combinations
            ! i think i can
            call calcFullstopLowering(ilut, excitInfo, excits, n_excits, &
                                      posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstop_raising) ! full stop 2 raising
            call calcFullstopRaising(ilut, excitInfo, excits, n_excits, &
                                     posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstop_L_to_R) ! full stop lowering into raising
            call calcFullStopL2R(ilut, excitInfo, excits, n_excits, &
                                 posSwitches, negSwitches, t_no_singles_opt=.true.)

            ! in this case there is also the possibility for one single-like
            ! excitation if there is no change in the double overlap region!
            ! todo! how to fix that? is that so important? its only max. 1

            ! depending on the full-stop step-value:
            ! if it is 3, all excitations are single like, and should be
            ! disregarded if we mimic stochastic sampling.

            ! if the end step value is 1 or 2, there is on Delta B = 0
            ! branch with ofc. multiple possible singles
            ! associated with it.
            ! i think i should use a optional flag to indicate that I want
            ! to mimick stochastic exctiations generation
            exlevel = 2

        case (excit_type%fullstop_R_to_L) ! full stop raising into lowering
            call calcFullStopR2L(ilut, excitInfo, excits, n_excits, &
                                 posSwitches, negSwitches, t_no_singles_opt=.true.)

            ! same as for 16
            exlevel = 2

        case (excit_type%fullstart_lowering) ! full start 2 lowering
            call calcFullStartLowering(ilut, excitInfo, excits, n_excits, &
                                       posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstart_raising) ! full start 2 raising
            call calcFulLStartRaising(ilut, excitInfo, excits, n_excits, &
                                      posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstart_L_to_R) ! full start lowering into raising
            call calcFullStartL2R(ilut, excitInfo, excits, n_excits, &
                                  posSwitches, negSwitches, t_no_singles_opt=.true.)

            ! same as for 16
            exlevel = 2

        case (excit_type%fullstart_R_to_L) ! full start raising into lowering
            call calcFullStartR2L(ilut, excitInfo, excits, n_excits, &
                                  posSwitches, negSwitches, t_no_singles_opt=.true.)

            ! same as for 16
            exlevel = 2

        case (excit_type%fullstart_stop_alike) ! full start into full stop alike
            call calcFullStartFullStopAlike(ilut, excitInfo, excits)
            n_excits = 1

            exlevel = 2

        case (excit_type%fullstart_stop_mixed) ! full start into full stop mixed
            call calcFullStartFullStopMixed(ilut, excitInfo, excits, n_excits, &
                                            posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        end select

        ! indicate the level of excitation IC for the remaining NECI code
        if (n_excits > 0) then
            do n = 1, n_excits
                ! for consistency also encode x0 and x1 in the ind_rdm_x0
                ! and ind_rdm_x1 entries
                ! damn.. i already added x0 and x1 in the routines above
                ! i think...
                call encode_rdm_ind(excits(:, n), contract_2_rdm_ind(i, j, k, l))
            end do
        end if

    end subroutine calc_all_excits_guga_rdm_doubles

    subroutine calc_all_excits_guga_rdm_singles(ilut, i, j, excits, n_excits)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: i, j
        integer(n_int), intent(out), allocatable :: excits(:, :)
        integer, intent(out) :: n_excits
        character(*), parameter :: this_routine = "calc_all_excits_guga_rdm_singles"

        type(ExcitationInformation_t) :: excitInfo
        type(WeightObj_t) :: weights
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer :: iEx, iOrb, ierr
        integer(n_int), allocatable :: tempExcits(:, :)
        real(dp) :: minusWeight, plusWeight

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        n_excits = 0

        excitInfo = excitationIdentifier(i, j)

        associate (gen1 => excitInfo%gen1, en => excitInfo%fullEnd, &
                   st => excitInfo%fullStart)

            if (gen1 /= 0) then
                if (current_stepvector(i) == 3 .or. current_stepvector(j) == 0) then
                    allocate(excits(0, 0), stat=ierr)
                    return
                end if
            end if

            call calcRemainingSwitches_excitInfo_single(excitInfo, posSwitches, negSwitches)

            weights = init_singleWeight(ilut, en)
            plusWeight = weights%proc%plus(posSwitches(st), currentB_ilut(st), weights%dat)
            minusWeight = weights%proc%minus(negSwitches(st), currentB_ilut(st), weights%dat)

            ! check compatibility of chosen indices

            if ((current_stepvector(st) == 1 .and. near_zero(plusWeight)) .or. &
                (current_stepvector(st) == 2 .and. near_zero(minusWeight)) .or. &
                near_zero(minusWeight + plusWeight)) then
                allocate(excits(0, 0), stat=ierr)
                return
            end if

            ! have to give probabilistic weight object as input, to deal
            call createSingleStart(ilut, excitInfo, posSwitches, &
                                   negSwitches, weights, tempExcits, n_excits)

            do iOrb = st + 1, en - 1
                call singleUpdate(ilut, iOrb, excitInfo, posSwitches, &
                                  negSwitches, weights, tempExcits, n_excits)
            end do

            call singleEnd(ilut, excitInfo, tempExcits, &
                           n_excits, excits)

            ! encode the combined RDM-ind in the deltaB position for
            ! communication purposes
            do iEx = 1, n_excits
                ! just for consistency all encode the x0 element in the
                ! ind_rdm_x0 entry
                call encode_stochastic_rdm_info(GugaBits, excits(:, iEx), &
                                                rdm_ind=contract_1_rdm_ind(i, j), &
                                                x0=extract_matrix_element(excits(:, iex), 1), x1=0.0_dp)
            end do

        end associate

    end subroutine calc_all_excits_guga_rdm_singles

    subroutine add_guga_lists_rdm(n_dets_tot, n_dets, list_tot, list)
        ! i need a new guga list addition for RDMs as I must not add
        ! identical excitations coming from different RDM index combinations
        ! otherwise I loose the information of those excitations!
        integer, intent(inout) ::        n_dets_tot
        integer, intent(in) ::           n_dets
        integer(n_int), intent(inout) :: list_tot(0:, 1:), list(0:, 1:)
        character(*), parameter :: this_routine = "add_guga_lists_rdm"

        ! here I essentially only need to append the list to the total
        ! list and add the number of new elements to n_dets_tot

        list_tot(:, n_dets_tot + 1:n_dets_tot + n_dets) = list(:, 1:n_dets)
        n_dets_tot = n_dets_tot + n_dets

    end subroutine add_guga_lists_rdm

    subroutine calc_rdm_energy_guga(rdm, rdm_energy_1, rdm_energy_2)
        ! the GUGA RDM version of the energy calculator

        type(rdm_list_t), intent(in) :: rdm
        real(dp), intent(out) :: rdm_energy_1(rdm%sign_length)
        real(dp), intent(out) :: rdm_energy_2(rdm%sign_length)
        character(*), parameter :: this_routine = "calc_rdm_energy_guga"

        integer(int_rdm) :: ijkl
        integer :: ielem, ij, kl, i, j, k, l
        real(dp) :: rdm_sign(rdm%sign_length)

        rdm_energy_1 = 0.0_dp
        rdm_energy_2 = 0.0_dp

        ! Loop over all elements in the 2-RDM.
        do ielem = 1, rdm%nelements
            ijkl = rdm%elements(0, ielem)
            call extract_sign_rdm(rdm%elements(:, ielem), rdm_sign)

            call extract_2_rdm_ind(ijkl, i, j, k, l)

            ! TODO maybe the indices are not correctly accessed here!
            rdm_energy_2 = rdm_energy_2 + rdm_sign * get_umat_el(i, k, j, l) / 2.0_dp

            if (i == j) then
                rdm_energy_1 = rdm_energy_1 + &
                               rdm_sign * GetTMatEl(2 * k, 2 * l) / real(nel - 1, dp) / 2.0_dp
            end if
            if (k == l) then
                rdm_energy_1 = rdm_energy_1 + &
                               rdm_sign * GetTMatEl(2 * i, 2 * j) / real(nel - 1, dp) / 2.0_dp
            end if
        end do

    end subroutine calc_rdm_energy_guga

end module guga_rdm
