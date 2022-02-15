#include "macros.h"
! GUGA module containg as much matrix element calculation functionality as
! possible.
module guga_matrixElements
    use SystemData, only: nEl, nBasis, ECore, t_tJ_model, t_heisenberg_model, t_new_hubbard, &
        t_new_real_space_hubbard, treal, nSpatOrbs, t_mixed_hubbard, ElecPairs, &
        is_init_guga
    use constants, only: dp, n_int, hel_zero, int_rdm, bn2_, Root2
    use bit_reps, only: decode_bit_det
    use OneEInts, only: GetTMatEl
    use procedure_pointers, only: get_umat_el
    use guga_bitRepOps, only: isDouble, calcB_vector_nI, isProperCSF_nI, CSF_Info_t, &
        convert_ilut_toGUGA, identify_excitation, findFirstSwitch, findLastSwitch, &
        calcb_vector_ilut, count_open_orbs_ij
    use guga_bitRepOps, only: contract_1_rdm_ind, contract_2_rdm_ind
    use util_mod, only: binary_search, operator(.isclose.), stop_all, near_zero, operator(.div.)
    use guga_data, only: projE_replica, ExcitationInformation_t, excit_type, gen_type

    use guga_data, only: funA_0_2overR2, minFunA_2_0_overR2, funA_2_0_overR2, &
                         funA_m1_1_overR2, funA_3_1_overR2, minFunA_0_2_overR2, WeightData_t

    use guga_data, only: getdoublematrixelement, getmixedfullstop, getSingleMatrixElement, &
        getdoublecontribution


    use guga_types, only: WeightObj_t

    use bit_rep_data, only: niftot, nifd
    use MPI_wrapper, only: iprocindex
    use CalcData, only: matele_cutoff, t_matele_cutoff, t_trunc_guga_pgen, &
        t_trunc_guga_pgen_noninits, trunc_guga_pgen, t_trunc_guga_matel, trunc_guga_matel
    use bit_rep_data, only: nifguga
    use DetBitOps, only: DetBitEQ


    use guga_procedure_pointers, only: &
        calc_orbital_pgen_contrib_start, calc_orbital_pgen_contrib_end, &
        calc_orbital_pgen_contr


    use FciMCData, only: tFillingStochRDMOnFly

    better_implicit_none

    private
    public :: calc_guga_matrix_element, calc_mixed_contr_integral, calc_integral_contribution_single
    public :: calcDiagMatEleGuga_nI, calcDiagExchangeGUGA_nI, calcDiagMatEleGUGA_ilut
    public :: calcremainingswitches_excitinfo_double, calcremainingswitches_excitinfo_single, &
                calcstartprob, calcstayingprob, endfx, endgx, &
                init_fullstartweight, init_singleweight

    public :: calc_mixed_start_contr_sym, calc_mixed_end_contr_sym, calc_mixed_contr_sym, &
        calc_mixed_start_contr_integral, calc_mixed_start_contr_pgen, &
        calc_mixed_end_contr_pgen, calc_mixed_end_contr_integral, &
        calc_mixed_contr_pgen

    public :: get_forced_zero_double, getminus_double, getminus_semistart, getplus_double, getplus_semistart, init_doubleweight, init_semistartweight

    abstract interface
        ! to make the recalculation of the branch-weights for mixed full-stop
        ! excitations more efficient, set up an array of functions which
        ! give the corresponding branching tree decision to easily recalculate
        ! the branch_pgen for different fullends without having to check the
        ! taken path for every iteration
        function branch_weight_function(weight, bVal, negSwitches, posSwitches) &
            result(prob)
            use constants, only: dp
            use guga_types, only: WeightObj_t
            implicit none
            type(WeightObj_t), intent(in) :: weight
            real(dp), intent(in) :: bval, negSwitches, posSwitches
            real(dp) :: prob
        end function branch_weight_function
    end interface

    type :: BranchWeightArr_t
        procedure(branch_weight_function), pointer, nopass :: ptr => null()
    end type BranchWeightArr_t

contains

    pure subroutine calc_guga_matrix_element(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, t_hamil, &
                                        rdm_ind, rdm_mat)
        ! function which, given the 2 CSFs ilutI/J and the excitation
        ! information, connecting those 2, calculates the Hamiltionian
        ! matrix element between those 2
        ! use a flag to distinguish between only guga-mat_ele calculation
        ! and full hamiltonian matrix element calculation

        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(out) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in) :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_guga_matrix_element"

        integer(n_int) :: tmp_i(0:nifguga)

        mat_ele = h_cast(0.0_dp)

        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))

        ! check diagonal case first
        if (DetBitEQ(ilutI, ilutJ)) then

            call convert_ilut_toGUGA(ilutI, tmp_i)

            ! i think I should do a change here for the Heisenberg or
            ! tJ model..
            mat_ele = calcDiagMatEleGUGA_ilut(tmp_i)
            return
        end if

        excitInfo = identify_excitation(ilutI, ilutJ)

        ! more than a double excitation! leave with 0 matrix element
        if (.not. excitInfo%valid) return

        ! for the hubbard model implementation, depending if it is in the
        ! momentum- or real-space i can get out of here if we identify
        ! excitation types which are definetly 0
        ! i could do that more efficiently if we identify it already
        ! earlier but for now do it here!
        if (t_hamil) then
            ! make this check only if we want the hamiltonian matrix
            ! element. for general coupling coefficients (eg. for RDMs)
            ! i do need this contributions anyway
            if (t_new_hubbard) then
                if (treal .or. t_new_real_space_hubbard) then
                    ! only singles in the real-space hubbard!
                    if (excitInfo%typ /= excit_type%single) return
                else
                    ! only double excitations in the momentum-space hubbard!
                    if (excitInfo%typ == excit_type%single) return
                end if
            end if

            ! make the adjustment for the Heisenberg model
            if (t_heisenberg_model &
                .and. excitInfo%typ /= excit_type%fullstart_stop_mixed) then
                return
            end if

            if (t_tJ_model .and. &
                (.not. (excitInfo%typ == excit_type%single &
                 .or. excitInfo%typ == excit_type%fullstart_stop_mixed))) &
                return
        end if

        ! i think in the excitation identification i can not find out if the
        ! delta B value is abs>2 so i have to do that here.. or specific for
        ! the type of excitations below.. for singles its not allowed
        ! abs > 1 ..
        ! but i think for the double excitations i cannot do that generally
        ! since i have to check for the overlap and non-overlap regions
        ! specifically

        ! then i need a select case to specifically calculate all the
        ! different types of excitations
        ! essentially i can just mimick the stochastic excitation creation
        ! routines but with a fixed chosen excitation
        select case (excitInfo%typ)

        case (excit_type%single)
            ! pure single excitation

            ! but here i have to calculate all the double excitation
            ! influences which can lead to the same excitation(weights etc.)
            call calc_single_excitation_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                           t_hamil, rdm_ind, rdm_mat)

        case (excit_type%single_overlap_L_to_R)
            ! single overlap lowering into raising

            ! maybe i have to check special conditions on the overlap site.
            call calc_single_overlap_mixed_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                              t_hamil, rdm_ind, rdm_mat)

        case (excit_type%single_overlap_R_to_L)
            ! single overlap raising into lowering

            ! maybe i have to check special conditions on the overlap site.
            call calc_single_overlap_mixed_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                              t_hamil, rdm_ind, rdm_mat)

        case (excit_type%double_lowering)
            ! normal double lowering

            ! question is can i combine more functions here since i know
            ! both CSFs.. i think so!

            ! deal with order parameter for switched indices
            call calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)

        case (excit_type%double_raising)
            ! normal double raising

            ! here i have to deal with the order parameter for switched
            ! indices ..
            call calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)

        case (excit_type%double_L_to_R_to_L)
            ! lowering into raising into lowering

            ! can i combine these 4 similar excitations in one routine?
            ! deal with non-overlap if no spin-coupling changes!
            call calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)

        case (excit_type%double_R_to_L_to_R)
            ! raising into lowering into raising

            ! here i have to consider the non-overlap contribution if no
            ! spin-coupling changes in the overlap range
            call calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)

        case (excit_type%double_L_to_R)
            ! lowering into raising double

            ! consider non-overlap if no spin-coupling changes!
            call calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)

        case (excit_type%double_R_to_L)
            ! raising into lowering double

            ! here i also have to consider the non-overlap contribution if no
            ! spin-coupling changes in the overlap range
            call calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstop_lowering)
            ! full-stop 2 lowering

            ! here only x0 matrix element in overlap range!
            ! also combine fullstop-alike
            call calc_fullstop_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                        t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstop_raising)
            ! full-stop 2 raising

            ! here only x0 matrix elment in overlap range!
            call calc_fullstop_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                        t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstop_L_to_R)
            ! full-stop lowering into raising

            ! here i have to consider all the singly occupied orbital
            ! influences ABOVE the last spin-coupling change
            ! this is more of a pain.. do later!
            ! finished the "easy" ones.. now to the annoying..
            call calc_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
                                        t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstop_R_to_L)
            ! full-stop raising into lowering


            ! here i have to consider all the singly occupied orbital
            ! influences ABOVE the last spin-coupling change
            call calc_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
                                        t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstart_lowering)
            ! full-start 2 lowering

            ! here only x0 matrix element in overlap range!
            call calc_fullstart_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                         t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstart_raising)
            ! full-start 2 raising

            ! here only the x0-matrix in the overlap range (this implies no
            ! spin-coupling changes, but i already dealt with that! (hopefully!))
            call calc_fullstart_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                         t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullStart_L_to_R)
            ! full-start lowering into raising

            ! here i have to consider all the other singly occupied orbital
            ! influences BELOW the first spin-coupling change
            call calc_fullstart_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
                                         t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstart_R_to_L)
            ! full-start raising into lowering

            ! here i have to consider all the other singly occupied orbital
            ! influences BELOW the first spin-coupling change
            call calc_fullstart_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
                                         t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstart_stop_alike)
            ! full-start into full-stop alike

            ! here no spin-coupling changes are allowed!
            call calc_fullstart_fullstop_alike_ex(csf_i, excitInfo, &
                                                  mat_ele, t_hamil, rdm_ind, rdm_mat)

        case (excit_type%fullstart_stop_mixed)
            ! full-start into full-stop mixed

            ! here i have to consider all the singly occupied orbitals
            ! below the first spin-change and above the last spin change
            call calc_fullstart_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, &
                                                  mat_ele, t_hamil, rdm_ind, rdm_mat)

        case default
            call stop_all(this_routine, "unexpected excitation type")

        end select

        if (t_matele_cutoff .and. abs(mat_ele) < matele_cutoff) mat_ele = h_cast(0.0_dp)

    end subroutine calc_guga_matrix_element



    pure function calcDiagMatEleGuga_nI(nI) result(hel_ret)
        ! calculates the diagonal Hamiltonian matrix element when a CSF in
        ! nI(nEl) form is provided and returns hElement of type hElement_t
        integer, intent(in) :: nI(nEl)
        HElement_t(dp) :: hel_ret

        ! have to loop over the number of spatial orbitals i , and within
        ! loop again over orbitals j > i, s indicates spatial orbitals
        integer :: iOrb, jOrb, inc1, inc2, sOrb, pOrb
        real(dp) :: nOcc1, nOcc2

        hel_ret = ECore

        iOrb = 1
        ! loop over nI spin orbital entries: good thing is unoccupied orbitals
        ! do not  contribute to the single matrix element part.
        do while (iOrb <= nEl)
            ! spatial orbital index needed for get_umat_el access
            sOrb = (nI(iOrb) + 1) / 2
            ! have to check if orbital is singly or doubly occupied.
            if (isDouble(nI, iOrb)) then ! double has two part. int. contribution
                nOcc1 = 2.0_dp
                hel_ret = hel_ret + nOcc1 * GetTMatEl(nI(iOrb), nI(iOrb)) + &
                          get_umat_el(sOrb, sOrb, sOrb, sOrb)

                ! correctly count through spin orbitals if its a double occ.
                inc1 = 2

            else ! single occupation
                nOcc1 = 1.0_dp
                hel_ret = hel_ret + nOcc1 * GetTMatEl(nI(iOrb), nI(iOrb))
                inc1 = 1

            end if

            ! second loop:
            jOrb = iOrb + inc1
            do while (jOrb <= nEl)
                pOrb = (nI(jOrb) + 1) / 2
                ! again check for double occupancies
                if (isDouble(nI, jOrb)) then
                    nOcc2 = 2.0_dp
                    inc2 = 2

                else
                    nOcc2 = 1.0_dp
                    inc2 = 1

                end if
                ! standard two particle contribution
                if (.not. (t_tJ_model .or. t_heisenberg_model)) then
                    hel_ret = hel_ret + nOcc1 * nOcc2 * ( &
                              get_umat_el(sOrb, pOrb, sOrb, pOrb) - &
                              get_umat_el(sOrb, pOrb, pOrb, sOrb) / 2.0_dp)
                end if

                ! calculate exchange integral part, involving Shavitt
                ! rules for matrix elements, only contributes if both
                ! stepvectors are 1 or 2, still to be implemented..
                if ((nOcc1.isclose.1.0_dp) .and. (nOcc2.isclose.1.0_dp)) then
                    hel_ret = hel_ret - get_umat_el(sOrb, pOrb, pOrb, sOrb) * &
                              calcDiagExchangeGUGA_nI(iOrb, jOrb, nI) / 2.0_dp

                    if (t_tJ_model) then
                        hel_ret = hel_ret + get_umat_el(sOrb, pOrb, pOrb, sOrb) / 2.0
                    end if
                end if

                ! increment the counters
                jOrb = jOrb + inc2

            end do
            iOrb = iOrb + inc1
        end do

    end function calcDiagMatEleGUGA_nI

    pure function calcDiagMatEleGuga_ilut(ilut) result(hElement)
        ! function to calculate the diagonal matrix element if a stepvector
        ! in ilut format is given
        integer(n_int), intent(in) :: ilut(0:niftot)
        HElement_t(dp) :: hElement

        integer :: nI(nEl)

        call decode_bit_det(nI, ilut)

        hElement = calcDiagMatEleGUGA_nI(nI)

    end function calcDiagMatEleGuga_ilut

    pure function calcDiagExchangeGUGA_nI(iOrb, jOrb, nI) result(exchange)
        ! calculates the exchange contribution to the diagonal matrix elements
        ! this is the implementation if only nI is provided
        integer, intent(in) :: iOrb, jOrb, nI(nEl)
        real(dp) :: exchange

        real(dp) :: bVector(nEl)
        integer :: i

        ! the b-vector is also needed for these calculations:
        bVector = calcB_vector_nI(nI)
        ! probably could use current b vector.. or reference b vector even...
        ! yes definitly. no not really since this is also used for general
        ! diagonal matrix elements not only the current determinant in the
        ! fciqmc loop

        exchange = 1.0_dp
        ! then i need the exchange product term between orbitals s and p
        do i = iOrb + 1, jOrb - 1
            if (.not. isDouble(nI, i)) then
                if (is_beta(nI(i))) then
                    exchange = exchange * functionA(bVector(i), 2.0_dp, 0.0_dp) &
                               * functionA(bVector(i), -1.0_dp, 1.0_dp)

                else
                    exchange = exchange * functionA(bVector(i), 0.0_dp, 2.0_dp) * &
                               functionA(bVector(i), 3.0_dp, 1.0_dp)

                end if
            end if
        end do

        if (is_beta(nI(iOrb))) then
            exchange = exchange * functionA(bVector(iOrb), 2.0_dp, 0.0_dp)

            if (is_beta(nI(jOrb))) then
                exchange = exchange * functionA(bVector(jOrb), -1.0_dp, 1.0_dp)

            else
                exchange = -exchange * functionA(bVector(jOrb), 3.0_dp, 1.0_dp)

            end if

        else
            exchange = exchange * functionA(bVector(iOrb), 0.0_dp, 2.0_dp)

            if (is_beta(nI(jOrb))) then
                exchange = -exchange * functionA(bVector(jOrb), -1.0_dp, 1.0_dp)

            else
                exchange = exchange * functionA(bVector(jOrb), 3.0_dp, 1.0_dp)

            end if
        end if

    end function calcDiagExchangeGUGA_nI

    elemental function functionA(bValue, x, y) result(r)
        ! calculated the "A" function used for Shavitts matrix element calc.
        real(dp), intent(in) :: bValue, x, y
        real(dp) :: r

        r = sqrt((bValue + x) / (bValue + y))

    end function functionA

    pure subroutine calc_single_excitation_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                         t_calc_full, rdm_ind, rdm_mat)
        ! routine to exactly calculate the matrix element between so singly
        ! connected CSFs, with the option to output also all the indices and
        ! overlap matrix elements necessary for the rdm calculation
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_calc_full
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_single_excitation_ex"

        integer :: iOrb, db, step1, step2
        real(dp) :: bVal
        HElement_t(dp) :: integral
        logical :: t_calc_full_
        ! just mimick the stochastic single excitation!
        real(dp) :: tmp_mat
        integer :: delta_b(nSpatOrbs)

        def_default(t_calc_full_, t_calc_full, .true.)

        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))

        ! set defaults for the output if we excit early..
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (i => excitInfo%i, j => excitInfo%j, st => excitInfo%fullstart, &
                   en => excitInfo%fullEnd, gen => excitInfo%currentGen)

            if (present(rdm_ind) .and. present(rdm_mat)) then
                allocate(rdm_ind(1), source=0_int_rdm)
                allocate(rdm_mat(1), source=0.0_dp)
                rdm_ind = contract_1_rdm_ind(i, j)
                rdm_mat = 0.0_dp
            end if

            ! this deltaB info can slip through the excitation identifier..
            if (any(abs(delta_b) > 1)) return

            if (t_calc_full_) then
                integral = getTmatEl(2 * i, 2 * j)
            else
                integral = h_cast(1.0_dp)
            end if

            tmp_mat = 1.0_dp

            ! what do i need from the 2 CSFs to calculate all the matrix elements?
            ! the 2 stepvalues: d, d' -> write a routine which only calculates those
            ! the generator type: gen
            ! the currentB value of I, b -> also calc. that in the routine to get d
            ! and the deltaB value: db

            ! since i mostly use this function to calculate the overlap with
            ! the reference determinant it is probably a good idea to
            ! calculate the necessary information for the reference determinant
            ! once and store it.. and maybe use a flag here to check what
            ! kind of usage of this function is.. or outside in the calling
            ! function..

            ! then i can just loop over the whole excitation range without
            ! distinction between start.. (mabye end i have to deal with specifically!)

            ! here i have to calculate the starting element
            !
            step1 = csf_i%stepvector(st)
            step2 = csf_j%stepvector(st)
            bVal = csf_i%b_real(st)
            ! i think i have to take the deltaB for the next orb or?
            ! because i need the outgoing deltaB value..
            ! nah.. i need the incoming!
            ! damn i think i need the deltaB value from the previous orbital..
            ! no.. for the single start i need the outgoing deltab, that the
            ! convention how to access the matrix elements, but then i need the
            ! incoming deltaB value.. or lets check that out how this works
            ! exactly
            db = delta_b(st)

            tmp_mat = tmp_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

            ! i think it can still always happen that the matrix element is 0..
            ! but maybe for the rdm case i have to do something more involved..
            ! to set all the corresponding indices to 0..

            if (near_zero(tmp_mat)) return

            do iOrb = st + 1, en - 1

                step1 = csf_i%stepvector(iOrb)
                step2 = csf_j%stepvector(iOrb)
                bVal = csf_i%b_real(iOrb)
                db = delta_b(iOrb - 1)

                tmp_mat = tmp_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

                if (near_zero(tmp_mat)) return

                if (t_calc_full_) then
                    if (.not. (t_new_real_space_hubbard .or. t_heisenberg_model &
                               .or. t_tJ_model .or. t_mixed_hubbard)) then
                        integral = integral + get_umat_el(i, iOrb, j, iOrb) * csf_i%Occ_real(iOrb)

                        integral = integral + get_umat_el(i, iOrb, iOrb, j) * &
                                   getDoubleContribution(step2, step1, db, gen, bVal)
                    end if
                end if

            end do

            step1 = csf_i%stepvector(en)
            step2 = csf_j%stepvector(en)
            bVal = csf_i%b_real(en)
            db = delta_b(en - 1)

            tmp_mat = tmp_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

            if (near_zero(tmp_mat)) return

            ! i think i could also exclude the treal case here.. try!
            if (t_calc_full_) then
                if (.not. (treal .or. t_new_real_space_hubbard .or. &
                           t_heisenberg_model .or. t_tJ_model .or. t_mixed_hubbard)) then
                    call calc_integral_contribution_single(csf_i, csf_j, i, j, st, en, integral)
                end if
            end if

            mat_ele = tmp_mat * integral

            if (present(rdm_mat)) rdm_mat = tmp_mat

        end associate

    end subroutine calc_single_excitation_ex

    pure subroutine calc_single_overlap_mixed_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                            t_calc_full, rdm_ind, rdm_mat)
        ! routine to exactly calculate the matrix element between 2 CSFs
        ! connected by a single overlap excitation with mixed generators
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_calc_full
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_single_overlap_mixed_ex"

        real(dp) :: bVal, temp_mat, guga_mat
        HElement_t(dp) :: umat
        integer :: i, db, gen2, step1, step2, delta_b(nSpatOrbs)
        logical :: t_calc_full_

        ! in the case of rdm calculation, i know that this type of exitation
        ! only has one (or two, with switches 2-body integrals..??)
        ! rdm-contribution..

        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))
        def_default(t_calc_full_, t_calc_full, .true.)

        ! set some defaults in case of early exit
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, st => excitInfo%fullStart, &
                   ss => excitInfo%secondStart, en => excitInfo%fullEnd, &
                   gen => excitInfo%firstGen, fe => excitInfo%firstEnd, &
                   typ => excitInfo%typ)

            if (present(rdm_ind) .and. present(rdm_mat)) then
                ! i am not sure yet if I will use symmetries in the RDM
                ! calculation (some are also left out in the SD based implo..
                ! so for now sample both combinations
                allocate(rdm_ind(1), source=0_int_rdm)
                allocate(rdm_mat(1), source=0.0_dp)
                rdm_ind(1) = contract_2_rdm_ind(ii, jj, kk, ll)
            end if

            if (any(abs(delta_b) > 1)) return

            if (t_calc_full_) then
                if (typ == excit_type%single_overlap_L_to_R) then

                    umat = (get_umat_el(fe, ss, st, en) + &
                            get_umat_el(ss, fe, en, st)) / 2.0_dp

                else if (typ == excit_type%single_overlap_R_to_L) then
                    umat = (get_umat_el(st, en, fe, ss) + &
                            get_umat_el(en, st, ss, fe)) / 2.0_dp
                else
                    call stop_all(this_routine, "shouldnt be here!")
                end if
            else
                umat = h_cast(1.0_dp)
            end if

            ! for the hamiltonian matrix element i can exit here if umat is 0
            ! but for the rdm-contribution i need to calc. the GUGA element
            if (t_calc_full .and. near_zero(umat)) return

            guga_mat = 1.0_dp

            ! i have to do the start specifically, due to the access of the
            ! single matrix elements
            step1 = csf_i%stepvector(st)
            step2 = csf_j%stepvector(st)
            bVal = csf_i%b_real(st)
            db = delta_b(st)

            guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

            if (near_zero(guga_mat)) return

            do i = st + 1, ss - 1

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                bVal = csf_i%b_real(i)
                db = delta_b(i - 1)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

                if (near_zero(guga_mat)) return

            end do

            step1 = csf_i%stepvector(ss)
            step2 = csf_j%stepvector(ss)
            bVal = csf_i%b_real(ss)
            db = delta_b(ss - 1)
            gen2 = -gen

            call getDoubleMatrixElement(step2, step1, db, gen, gen2, bVal, 1.0_dp, temp_mat)

            guga_mat = guga_mat * temp_mat

            if (near_zero(guga_mat)) return

            do i = ss + 1, en

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                bVal = csf_i%b_real(i)
                db = delta_b(i - 1)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen2, bVal)

                if (near_zero(guga_mat)) return

            end do

        end associate

        mat_ele = guga_mat * umat

        ! both coupling coeffs are the same
        if (present(rdm_mat)) rdm_mat = guga_mat

    end subroutine calc_single_overlap_mixed_ex

    pure subroutine calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                     t_hamil, rdm_ind, rdm_mat)
        ! combined routine to calculate the mixed generator excitations with
        ! 4 different spatial orbitals. here i have to consider if a
        ! spin change happened in the overlap. if no, i also have to calc. the
        ! non-overlap contribution!
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_normal_double_ex"

        integer :: step1, step2, db, i, delta_b(nSpatOrbs)
        real(dp) :: temp_x0, temp_x1, temp_mat0, temp_mat1, bVal, guga_mat
        logical :: t_hamil_

        def_default(t_hamil_, t_hamil, .true.)
        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))

        ! set defaults for early exits
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, start1 => excitInfo%fullStart, &
                   start2 => excitInfo%secondStart, ende1 => excitInfo%firstEnd, &
                   ende2 => excitInfo%fullEnd, gen1 => excitInfo%gen1, &
                   gen2 => excitInfo%gen2, firstgen => excitInfo%firstgen, &
                   lastgen => excitInfo%lastgen, order => excitInfo%order, &
                   order1 => excitInfo%order1, typ => excitInfo%typ)

            if (present(rdm_ind) .and. present(rdm_mat)) then
                allocate(rdm_ind(2), source=0_int_rdm)
                allocate(rdm_mat(2), source=0.0_dp)

                ! this does get tricky now with the rdm inds and mats..
                ! the indices which end up here should always be intertwined..
                ! as we always deal with an overlap range in the excitation
                ! generation
                ASSERT(max(ii, jj) > min(kk, ll))

                ! so the first two entries correspond to the overlap version
                ! of the generators (in the case of mixed R and L)
                rdm_ind(1) = contract_2_rdm_ind(ii, jj, kk, ll)
                ! and the second to the non-overlap (again only in the case of
                ! mixed generator combinations!)
                rdm_ind(2) = contract_2_rdm_ind(ii, ll, kk, jj)

            end if

            ! i have to check if the deltaB value is in its correct bounds for
            ! the specific parts of the excitations
            if (any(abs(delta_b(start1:start2 - 1)) > 1) .or. &
                any(abs(delta_b(start2:ende1 - 1)) > 2) .or. &
                any(abs(delta_b(ende1:ende2)) > 1)) return

            ! depending on which type of excitation the non-overlap has a specific
            ! tbd generator! i should deal with that in the starting if statement
            ! according to my stochastic implementation it is not so hard luckily..

            ! do first single part
            guga_mat = 1.0_dp

            ! have to do the start specifically as i need the outgoing deltaB
            ! value
            step1 = csf_i%stepvector(start1)
            step2 = csf_j%stepvector(start1)
            db = delta_b(start1)
            bVal = csf_i%b_real(start1)

            guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, firstgen, bVal)

            if (near_zero(guga_mat)) return

            do i = start1 + 1, start2 - 1

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                bVal = csf_i%b_real(i)
                db = delta_b(i - 1)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, firstgen, bVal)

                if (near_zero(guga_mat)) return

            end do

            temp_x0 = guga_mat
            temp_x1 = guga_mat

            ! do overlap part
            do i = start2, ende1

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                bVal = csf_i%b_real(i)
                db = delta_b(i - 1)

                call getDoubleMatrixElement(step2, step1, db, gen1, gen2, bVal, 1.0_dp, &
                                            temp_mat0, temp_mat1)

                temp_x0 = temp_x0 * temp_mat0
                temp_x1 = temp_x1 * temp_mat1

                if (near_zero(temp_x0) .and. near_zero(temp_x1)) return
            end do

            ! now do second single overlap part

            do i = ende1 + 1, ende2

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                bVal = csf_i%b_real(i)
                db = delta_b(i - 1)

                temp_mat0 = getSingleMatrixElement(step2, step1, db, lastgen, bVal)

                temp_x0 = temp_x0 * temp_mat0
                temp_x1 = temp_x1 * temp_mat0

                if (near_zero(temp_x0) .and. near_zero(temp_x1)) return
            end do

            ! now the excitation-type and spin question come into play..
            ! i actually could combine all "normal" double excitations in one
            ! routine with a if statement at the end here..

            if (present(rdm_mat)) then
                select case (typ)
                case (excit_type%double_lowering, &
                      excit_type%double_raising)

                    ! in both the orbital picker and also excitation identifier
                    ! the order quantities are setup to be 1.0..
                    ! i hope I did everything right there.. then this would
                    ! mean here
                    ASSERT(order.isclose.1.0_dp)
                    ASSERT(order1.isclose.1.0_dp)
                    ! and now I have to fill in the correct combinations of
                    ! signs here..
                    ! to check that everything is correctly set to the default
                    ! ordering I should assert that here
#ifdef DEBUG_
                    if (typ == excit_type%double_raising) then
                        ASSERT(ii < jj)
                        ASSERT(kk < ll)
                        ASSERT(ii > kk)
                        ASSERT(jj > ll)
                    end if
#endif
                    ! be sure with the rdm_sign function:

                    ! rdm_mat(1) = temp_x0 + generator_sign(ii,jj,kk,ll) * temp_x1
                    ! rdm_mat(2) = temp_x0 + generator_sign(ii,ll,kk,jj) * temp_x1

                    rdm_mat(1) = temp_x0 + temp_x1
                    rdm_mat(2) = temp_x0 - temp_x1

                case (excit_type%double_L_to_R_to_L, &
                      excit_type%double_R_to_L_to_R, &
                      excit_type%double_L_to_R, &
                      excit_type%double_R_to_L)

                    if (excitInfo%spin_change) then
                        ! if we have a spin-change the non-overlap contribution
                        ! must be 0! which is already intiated to 0 above
                        rdm_mat(1) = temp_x1
                    else
                        ! if there is not spin-change the non-overlap is also
                        ! 0, and in this case is -2 the original x0
                        rdm_mat(1) = temp_x0 + temp_x1
                        rdm_mat(2) = -2.0_dp * temp_x0
                    end if
                end select
            end if

            select case (typ)
            case (excit_type%double_lowering)
                ! double lowering
                ! wait a minute.. i have to do that at the end apparently..
                ! since i need to know the x0 and x1 matrix element contributions
                mat_ele = (temp_x0 * (get_umat_el(ende1, ende2, start1, start2) + &
                                      get_umat_el(ende2, ende1, start2, start1) + &
                                      get_umat_el(ende2, ende1, start1, start2) + &
                                      get_umat_el(ende1, ende2, start2, start1)) + &
                           order * order1 *  temp_x1 * ( &
                                      get_umat_el(ende1, ende2, start1, start2) + &
                                      get_umat_el(ende2, ende1, start2, start1) - &
                                      get_umat_el(ende2, ende1, start1, start2) - &
                                      get_umat_el(ende1, ende2, start2, start1))) / 2.0_dp

            case (excit_type%double_raising)
                ! double raising
                mat_ele = (temp_x0 * (get_umat_el(start1, start2, ende1, ende2) + &
                                      get_umat_el(start2, start1, ende2, ende1) + &
                                      get_umat_el(start1, start2, ende2, ende1) + &
                                      get_umat_el(start2, start1, ende1, ende2)) + &
                           order * order1 * temp_x1 * ( &
                                      get_umat_el(start1, start2, ende1, ende2) + &
                                      get_umat_el(start2, start1, ende2, ende1) - &
                                      get_umat_el(start1, start2, ende2, ende1) - &
                                      get_umat_el(start2, start1, ende1, ende2))) / 2.0_dp

            case (excit_type%double_L_to_R_to_L)
                ! L -> R -> L
                if (excitInfo%spin_change) then
                    ! if a spin-change happenend -> no non-overlap!
                    mat_ele = temp_x1 * (get_umat_el(ende2, start2, start1, ende1) + &
                                         get_umat_el(start2, ende2, ende1, start1)) / 2.0_dp

                else
                    mat_ele = (-temp_x0 * (get_umat_el(start2, ende2, start1, ende1) + &
                                           get_umat_el(ende2, start2, ende1, start1)) * 2.0_dp + &
                              (temp_x0 + temp_x1) * (get_umat_el(ende2, start2, start1, ende1) + &
                                                     get_umat_el(start2, ende2, ende1, start1))) / 2.0_dp

                end if

            case (excit_type%double_R_to_L_to_R)
                ! R -> L -> R
                if (excitInfo%spin_change) then
                    mat_ele = temp_x1 * (get_umat_el(start1, ende1, ende2, start2) + &
                                         get_umat_el(ende1, start1, start2, ende2)) / 2.0_dp

                else
                    mat_ele = (-temp_x0 * (get_umat_el(start1, ende1, start2, ende2) + &
                                           get_umat_el(ende1, start1, ende2, start2)) * 2.0_dp + &
                              (temp_x0 + temp_x1) * (get_umat_el(start1, ende1, ende2, start2) + &
                                                     get_umat_el(ende1, start1, start2, ende2))) / 2.0_dp

                end if

            case (excit_type%double_L_to_R)
                ! L -> R
                if (excitInfo%spin_change) then
                    mat_ele = temp_x1 * (get_umat_el(ende1, start2, start1, ende2) + &
                                         get_umat_el(start2, ende1, ende2, start1)) / 2.0_dp

                else
                    mat_ele = (-temp_x0 * (get_umat_el(start2, ende1, start1, ende2) + &
                                           get_umat_el(ende1, start2, ende2, start1)) * 2.0_dp + &
                              (temp_x0 + temp_x1) * (get_umat_el(ende1, start2, start1, ende2) + &
                                                     get_umat_el(start2, ende1, ende2, start1))) / 2.0_dp
                end if

            case (excit_type%double_R_to_L)
                ! R -> L
                if (excitInfo%spin_change) then
                    mat_ele = temp_x1 * (get_umat_el(start1, ende2, ende1, start2) + &
                                         get_umat_el(ende2, start1, start2, ende1)) / 2.0_dp
                else
                    mat_ele = (-temp_x0 * (get_umat_el(start1, ende2, start2, ende1) + &
                                           get_umat_el(ende2, start1, ende1, start2)) * 2.0_dp + &
                              (temp_x0 + temp_x1) * (get_umat_el(start1, ende2, ende1, start2) + &
                                                     get_umat_el(ende2, start1, start2, ende1))) / 2.0_dp
                end if
                ! combine the "normal" double RR/LL also in here, since the rest
                ! of the routine is totally the same!
            end select

        end associate

    end subroutine calc_normal_double_ex

    pure subroutine calc_fullstop_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                      t_hamil, rdm_ind, rdm_mat)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_fullstop_alike_ex"

        real(dp) :: bVal, temp_mat, nOpen, guga_mat
        HElement_t(dp) :: umat
        integer :: i, step1, step2, db, delta_b(nSpatOrbs)
        logical :: t_hamil_

        def_default(t_hamil_, t_hamil, .true.)
        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))

        ! set some defaults in case of early exit
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, typ => excitInfo%typ, st => excitInfo%fullStart, &
                   se => excitInfo%secondStart, gen => excitInfo%gen1, &
                   en => excitInfo%fullEnd)

            if (present(rdm_ind) .and. present(rdm_mat)) then
                allocate(rdm_ind(1), source=0_int_rdm)
                allocate(rdm_mat(1), source=0.0_dp)
                rdm_ind(1) = contract_2_rdm_ind(ii, jj, kk, ll)
            end if

            ! can i exclude every deltaB > 1, since only db = 0 allowed in
            ! double overlap region? i think so..
            if (any(abs(delta_b) > 1)) return

            if (t_hamil_) then
                if (typ == excit_type%fullstop_lowering) then
                    ! LL
                    umat = (get_umat_el(en, en, st, se) + &
                            get_umat_el(en, en, se, st)) / 2.0_dp

                else if (typ == excit_type%fullstop_raising) then
                    ! RR
                    umat = (get_umat_el(st, se, en, en) + &
                            get_umat_el(se, st, en, en)) / 2.0_dp
                end if
            else
                umat = h_cast(1.0_dp)
            end if

            if (t_hamil_ .and. near_zero(umat)) return

            guga_mat = 1.0_dp

            ! have to deal with start specifically
            step1 = csf_i%stepvector(st)
            step2 = csf_j%stepvector(st)
            bVal = csf_i%b_real(st)
            db = delta_b(st)

            guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

            if (near_zero(guga_mat)) return

            do i = st + 1, se - 1

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                db = delta_b(i - 1)
                bVal = csf_i%b_real(i)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen, bval)

                if (near_zero(guga_mat)) return

            end do

            ! do semi-start:
            step1 = csf_i%stepvector(se)
            step2 = csf_j%stepvector(se)
            db = delta_b(se - 1)
            bVal = csf_i%b_real(se)

            call getDoubleMatrixElement(step2, step1, db, gen, gen, bVal, 1.0_dp, temp_mat)

            guga_mat = guga_mat * temp_mat

            if (near_zero(guga_mat)) return

            nOpen = (-1.0_dp)**real(count_open_orbs_ij(csf_i, se + 1, en - 1), dp)

            ! is this the same for both type of gens?
            mat_ele = guga_mat * nOpen * Root2 * umat

            ! since the x1 element is 0 there is no sign influence from the
            ! order of generators!
            if (present(rdm_mat)) rdm_mat = guga_mat * nOpen * Root2

        end associate

    end subroutine calc_fullstop_alike_ex

    pure subroutine calc_fullstart_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_fullstart_alike_ex"

        integer :: i, step1, step2, db, delta_b(nSpatOrbs)
        real(dp) :: bVal, nOpen, temp_mat, guga_mat
        HElement_t(dp) :: umat
        logical :: t_hamil_

        def_default(t_hamil_, t_hamil, .true.)
        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))

        ! set defaults for early exits
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, start => excitInfo%fullstart, &
                   ende => excitInfo%fullEnd, semi => excitInfo%firstEnd, &
                   gen => excitInfo%firstGen, typ => excitInfo%typ)

            if (present(rdm_ind) .and. present(rdm_mat)) then
                allocate(rdm_ind(1), source=0_int_rdm)
                allocate(rdm_mat(1), source=0.0_dp)
                rdm_ind(1) = contract_2_rdm_ind(ii, jj, kk, ll)
            end if

            ! i think i can exclude every deltaB > 1 sinve only dB = 0 branch
            ! allowed for the alike..
            if (any(abs(delta_b) > 1)) return

            if (t_hamil_) then
                if (typ == excit_type%fullstart_lowering) then
                    ! LL
                    umat = (get_umat_el(ende, semi, start, start) + &
                            get_umat_el(semi, ende, start, start)) / 2.0_dp

                else if (typ == excit_type%fullstart_raising) then
                    ! RR
                    umat = (get_umat_el(start, start, semi, ende) + &
                            get_umat_el(start, start, ende, semi)) / 2.0_dp
                end if
            else
                umat = h_cast(1.0_dp)
            end if

            if (t_hamil_ .and. near_zero(umat)) return

            nOpen = real(count_open_orbs_ij(csf_i, start, semi - 1), dp)

            ! do semi-stop
            step1 = csf_i%stepvector(semi)
            step2 = csf_j%stepvector(semi)
            db = delta_b(semi - 1)
            bVal = csf_i%b_real(semi)

            call getDoubleMatrixElement(step2, step1, db, gen, gen, bVal, 1.0_dp, temp_mat)

            if (near_zero(temp_mat)) return

            guga_mat = Root2 * temp_mat * (-1.0_dp)**nOpen

            ! do single range
            do i = semi + 1, ende

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                db = delta_b(i - 1)
                bVal = csf_i%b_real(i)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

                if (near_zero(guga_mat)) return

            end do

            mat_ele = guga_mat * umat

            ! also no influence on coupling coefficient sign from generator order
            if (present(rdm_mat)) rdm_mat = guga_mat

        end associate

    end subroutine calc_fullstart_alike_ex

    pure subroutine calc_fullstart_fullstop_alike_ex(csf_i, excitInfo, &
                                                mat_ele, t_hamil, rdm_ind, rdm_mat)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_fullstart_fullstop_alike_ex"

        real(dp) :: nOpen, guga_mat
        HElement_t(dp) :: umat
        logical :: t_hamil_

        def_default(t_hamil_, t_hamil, .true.)
        ASSERT(present(rdm_ind) .eqv. present(rdm_mat))

        ! set defaults for early exit
        mat_ele = h_cast(0.0_dp)

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, start => excitInfo%fullStart, &
                   ende => excitInfo%fullEnd)

            if (present(rdm_ind) .and. present(rdm_mat)) then
                allocate(rdm_ind(1), source=0_int_rdm)
                allocate(rdm_mat(1), source=0.0_dp)
                ! only one element here, since indices are the same
                rdm_ind(1) = contract_2_rdm_ind(ii, jj, kk, ll)
            end if

            if (t_hamil_) then
                umat = get_umat_el(ii, ii, jj, jj) / 2.0_dp
            else
                umat = h_cast(1.0_dp)
            end if

            if (t_hamil_ .and. near_zero(umat)) return

            nOpen = real(count_open_orbs_ij(csf_i, start, ende), dp)

            guga_mat = 2.0_dp * (-1.0_dp)**nOpen
            mat_ele = guga_mat * umat

            if (present(rdm_mat)) rdm_mat = guga_mat

        end associate

    end subroutine calc_fullstart_fullstop_alike_ex

    pure subroutine calc_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
                                      t_hamil, rdm_ind, rdm_mat)
        ! from the excitInfo i know the first switch position.
        ! this makes things a bit easier for the exact calculation
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)

        integer :: i, step1, step2, db
        real(dp)  :: bVal, temp_mat0, temp_mat1, temp_x0, temp_x1
        HElement_t(dp) :: integral
        integer(n_int) :: tmp_I(0:nifguga), tmp_J(0:nifguga)

        integer :: st, se, en, delta_b(nSpatOrbs)
        real(dp) :: guga_mat
        logical :: t_hamil_

        def_default(t_hamil_, t_hamil, .true.)
        delta_b = csf_i%B_int - csf_j%B_int

        ! i can not associate to all stuff of excitInfo, since it will
        ! get changed later on..
        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, firstGen => excitInfo%firstgen, &
                   typ => excitInfo%typ)

            ! set defaults in case of early exit
            mat_ele = h_cast(0.0_dp)

            if (any(abs(delta_b(st:se - 1)) > 1) .or. &
                any(abs(delta_b(se:en)) > 2)) return

            ! first do single overlap region
            guga_mat = 1.0_dp

            step1 = csf_i%stepvector(st)
            step2 = csf_j%stepvector(st)
            bVal = csf_i%b_real(st)
            db = delta_b(st)

            guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, firstgen, bval)

            if (near_zero(guga_mat)) return

            do i = st + 1, se - 1

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                db = delta_b(i - 1)
                bVal = csf_i%b_real(i)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, firstgen, bVal)

                if (near_zero(guga_mat)) return

            end do

            ! i also do not have to consider if there is a d=3 at the end since
            ! i determine the last spin-coupling change and thus know it is
            ! definetly a singly occupied orbital at the end

            ! and actually the x0 matrix element has to be 0, otherwise it is
            ! not the excitation i thought it was.. do this as a check to
            ! abort if anything else happens
            temp_x0 = guga_mat
            temp_x1 = guga_mat

            do i = se, en - 1

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                db = delta_b(i - 1)
                bVal = csf_i%b_real(i)

                call getDoubleMatrixElement(step2, step1, db, gen_type%R, gen_type%L, bVal, &
                                            1.0_dp, temp_mat0, temp_mat1)

                temp_x0 = temp_x0 * temp_mat0
                temp_x1 = temp_x1 * temp_mat1

                if (near_zero(temp_x0) .and. near_zero(temp_x1)) return
            end do

            ! to the fullstop
            step1 = csf_i%stepvector(en)
            step2 = csf_j%stepvector(en)
            db = delta_b(en - 1)
            bVal = csf_i%b_real(en)

            call getMixedFullStop(step2, step1, db, bVal, temp_mat0, temp_mat1)

            temp_x0 = temp_x0 * temp_mat0
            temp_x1 = temp_x1 * temp_mat1

            ! so if x0 > 0 abort!

            ! and here misuse the stochastic routine to calculate the influence
            ! of all the other singly-occupied orbitals..
            ! i probably should write a function which does that only for the
            ! integral and stores the specific indices and matrix elements for
            ! the rdm calculation, but to that later! todo
            ! but to use this function i have to transform the 2 iluts to
            ! GUGA iluts and store the matrix element in them..

            call convert_ilut_toGUGA(ilutI, tmp_I)
            call convert_ilut_toGUGA(ilutJ, tmp_J)

            ! is the matrix element transformed within these functions?
            ! no apparently not and not event touched.. so no need to encode the
            ! matrix element


            block
            real(dp) :: discard
            discard = 1.0_dp
            if (t_hamil_ .or. (tFillingStochRDMOnFly .and. present(rdm_mat))) then
                if (typ == excit_type%fullstop_L_to_R) then
                    ! L -> R
                    ! what do i have to put in as the branch pgen?? does it have
                    ! an influence on the integral and matrix element calculation?
!                     call calc_mixed_end_contr_sym(tmp_I, csf_i, tmp_J, excitInfo, discard, &
!                                                   also_discard, integral, rdm_ind, rdm_mat)

                    call calc_mixed_end_contr_integral(tmp_I, csf_i, tmp_J, &
                        excitInfo, integral, rdm_ind, rdm_mat)

                    ! need to multiply by x1
                    if (present(rdm_mat)) rdm_mat = rdm_mat * temp_x1

                    mat_ele = temp_x1 * ((get_umat_el(en, se, st, en) + &
                                          get_umat_el(se, en, en, st)) / 2.0_dp + integral)

                else if (typ == excit_type%fullstop_R_to_L) then
                    ! R -> L
!                     call calc_mixed_end_contr_sym(tmp_I, csf_i, tmp_J, excitInfo, discard, &
!                                                   also_discard, integral, rdm_ind, rdm_mat)

                    call calc_mixed_end_contr_integral(tmp_I, csf_i, tmp_J, &
                        excitInfo, integral, rdm_ind, rdm_mat)

                    if (present(rdm_mat)) rdm_mat = rdm_mat * temp_x1

                    mat_ele = temp_x1 * ((get_umat_el(en, st, se, en) + &
                                          get_umat_el(st, en, en, se)) / 2.0_dp + integral)

                end if
            else
                mat_ele = h_cast(temp_x1)
            end if
            end block
        end associate

    end subroutine calc_fullstop_mixed_ex

    pure subroutine calc_fullstart_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
                                       t_hamil, rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)

        integer :: step1, step2, db, i
        real(dp) :: bVal, temp_mat0, temp_mat1, temp_x0, temp_x1
        HElement_t(dp) :: integral
        integer(n_int) :: tmp_I(0:nifguga), tmp_J(0:nifguga)
        logical :: t_hamil_

        integer :: st, en, se, delta_b(nSpatOrbs)
        real(dp) :: guga_mat

        def_default(t_hamil_, t_hamil, .true.)

        ! set defaults for early exit
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        ! i can not associate to all stuff of excitInfo, since it will
        ! get changed later on..
        st = excitInfo%fullStart
        en = excitInfo%fullEnd
        se = excitInfo%firstEnd

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, gen => excitInfo%lastGen, typ => excitInfo%typ)

            if (any(abs(delta_b(st:se - 1)) > 2) .or. &
                any(abs(delta_b(se:en)) > 1)) return

            ! do the full-start, and i know here that it is singly occupied

            step1 = csf_i%stepvector(st)
            step2 = csf_j%stepvector(st)
            ! to indicate the mixed fullstart matrix elements in the routine!
            db = -1
            bVal = csf_i%b_real(st)

            call getDoubleMatrixElement(step2, step1, -1, gen_type%L, gen_type%R, bVal, 1.0_dp, &
                                        temp_x0, temp_x1)

            if (near_zero(temp_x0) .and. near_zero(temp_x1)) return


            ! then do the double overlap range

            do i = st + 1, se

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                db = delta_b(i - 1)
                bVal = csf_i%b_real(i)

                call getDoubleMatrixElement(step2, step1, db, gen_type%L, gen_type%R, bVal, 1.0_dp, &
                                            temp_mat0, temp_mat1)

                temp_x0 = temp_x0 * temp_mat0
                temp_x1 = temp_x1 * temp_mat1

                if (near_zero(temp_x1)) return

            end do

            ! i think here i should also check if the x0 matrix element is non-zero..

            if (.not. near_zero(temp_x0)) return


            guga_mat = temp_x1
            ! do single range

            do i = se + 1, en

                step1 = csf_i%stepvector(i)
                step2 = csf_j%stepvector(i)
                db = delta_b(i - 1)
                bVal = csf_i%b_real(i)

                guga_mat = guga_mat * getSingleMatrixElement(step2, step1, db, gen, bVal)

                if (near_zero(guga_mat)) return

            end do


            ! then also reuse the stochastic routines to get the integral contribs

            call convert_ilut_toGUGA(ilutI, tmp_I)
            call convert_ilut_toGUGA(ilutJ, tmp_J)

            block
            real(dp) :: discard
            discard = 1.0_dp

            if (t_hamil .or. (tFillingStochRDMOnFly .and. present(rdm_mat))) then
!                 call calc_mixed_start_contr_sym(tmp_I, csf_i, tmp_J, excitInfo, discard, &
!                                                 also_discard, integral, rdm_ind, rdm_mat)
                call calc_mixed_start_contr_integral(tmp_I, csf_i, tmp_J, excitInfo, &
                    integral, rdm_ind, rdm_mat)

                if (present(rdm_mat)) rdm_mat = rdm_mat * guga_mat
                if (typ == excit_type%fullstart_L_to_R) then
                    mat_ele = guga_mat * ((get_umat_el(st, se, en, st) + &
                                           get_umat_el(se, st, st, en)) / 2.0_dp + integral)

                else if (typ == excit_type%fullstart_R_to_L) then
                    mat_ele = guga_mat * ((get_umat_el(st, en, se, st) + &
                                           get_umat_el(en, st, st, se)) / 2.0_dp + integral)

                end if
            end if
            end block
        end associate

    end subroutine calc_fullstart_mixed_ex

    pure subroutine calc_fullstart_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, &
                                                mat_ele, t_hamil, rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: mat_ele
        logical, intent(in), optional :: t_hamil
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)

        integer(n_int) :: tmp_I(0:nifguga), tmp_J(0:nifguga)
        logical :: t_hamil_
        integer :: delta_b(nSpatOrbs)

        def_default(t_hamil_, t_hamil, .true.)

        ! set default for early exits
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, start => excitInfo%fullstart, &
                   ende => excitInfo%fullEnd)

            if (any(abs(delta_b) > 2)) return

            call convert_ilut_toGUGA(ilutI, tmp_I)
            call convert_ilut_toGUGA(ilutJ, tmp_J)

            if (t_hamil_ .or. (tFillingStochRDMOnFly .and. present(rdm_mat))) then
                if (present(rdm_mat)) then
                    call calc_mixed_contr_integral(tmp_I, csf_i, tmp_J, start, ende, &
                                                    mat_ele, rdm_ind, rdm_mat)
                else
                    call calc_mixed_contr_integral(tmp_I, csf_i, tmp_J, start, ende, &
                                                mat_ele)
                end if
            end if
        end associate
    end subroutine calc_fullstart_fullstop_mixed_ex


    pure subroutine calc_integral_contribution_single(csf_i, csf_j, i, j, st, en, integral)
        ! calculates the double-excitaiton contribution to a single excitation
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        integer, intent(in) :: i, j, st, en
        HElement_t(dp), intent(inout) :: integral

        real(dp) :: botCont, topCont, tempWeight, prod
        integer :: iO, jO, step

        ! calculate the bottom contribution depending on the excited stepvalue
        select case (csf_i%stepvector(st))
        case (0)
            ! this implicates a raising st:
            if (csf_j%stepvector(st) == 1) then
                call getDoubleMatrixElement(1, 0, 0, gen_type%L, gen_type%R, csf_i%B_real(st), &
                                            1.0_dp, x1_element=botCont)

            else
                call getDoubleMatrixElement(2, 0, 0, gen_type%L, gen_type%R, csf_i%B_real(st), &
                                            1.0_dp, x1_element=botCont)
            end if

        case (3)
            ! implies lowering st
            if (csf_j%stepvector(st) == 1) then
                ! need tA(0,2)
                botCont = funA_0_2overR2(csf_i%B_real(st))

            else
                ! need -tA(2,0)
                botCont = minFunA_2_0_overR2(csf_i%B_real(st))
            end if

        case (1)
            botCont = funA_m1_1_overR2(csf_i%B_real(st))
            ! check which generator
            if (csf_j%stepvector(st) == 0) botCont = -botCont

        case (2)
            botCont = funA_3_1_overR2(csf_i%B_real(st))
            if (csf_j%stepvector(st) == 3) botCont = -botCont
        end select

        ! do top contribution also already

        select case (csf_i%stepvector(en))
        case (0)
            if (csf_j%stepvector(en) == 1) then
                topCont = funA_2_0_overR2(csf_i%B_real(en))
            else
                topCont = minFunA_0_2_overR2(csf_i%B_real(en))
            end if
        case (3)
            if (csf_j%stepvector(en) == 1) then
                topCont = minFunA_2_0_overR2(csf_i%B_real(en))
            else
                topCont = funA_0_2overR2(csf_i%B_real(en))
            end if
        case (1)
            topCont = funA_2_0_overR2(csf_i%B_real(en))
            if (csf_j%stepvector(en) == 3) topCont = -topCont

        case (2)
            topCont = funA_0_2overR2(csf_i%B_real(en))
            if (csf_j%stepvector(en) == 0) topCont = -topCont

        end select

        ! depending on i and j calulate the corresponding single and double
        ! integral weights and check if they are non-zero...
        ! gets quite involved... :( need to loop over all orbitals
        ! have to reset prod inside the loop each time!

        do iO = 1, st - 1
            ! no contribution if not occupied.
            if (csf_i%stepvector(iO) == 0) cycle
            ! else it gets a contrbution weighted with orbital occupation
            ! first easy part:
            integral = integral + get_umat_el(i, iO, j, iO) * csf_i%Occ_real(iO)

            ! also easy is the non-product involving part...
            integral = integral - get_umat_el(i, iO, iO, j) * &
                       csf_i%Occ_real(iO) / 2.0_dp

            ! the product part is annoying actually... but doesnt help... todo
            ! have to do a second loop for the product
            ! for the first loop iteration i have to access the mixed fullstart
            ! elements, with a deltaB = -1 value!!
            ! think about how to implement that !

            ! do i have to do anything for a d = 3 ? since x1-element is 0
            ! in this case anyway.. there should not be an influence.
            ! also if its a d = 2 with b = 0 the matrix element is also 0
            if (csf_i%stepvector(iO) == 3 .or. csf_i%B_int(iO) == 0) cycle

            step = csf_i%stepvector(iO)
            call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, csf_i%B_real(iO), &
                                        1.0_dp, x1_element=prod)

            ! and then do the remaining:
            do jO = iO + 1, st - 1
                ! need the stepvalue entries to correctly access the mixed
                ! generator matrix elements
                step = csf_i%stepvector(jO)
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, &
                                            csf_i%B_real(jO), 1.0_dp, x1_element=tempWeight)

                prod = prod * tempWeight
            end do
            prod = prod * botCont

            integral = integral + get_umat_el(i, iO, iO, j) * prod

        end do

        ! also have to calc the top and bottom contribution to the product terms
        ! BUT they also depend on the excitation -> so i should only calculate
        ! these terms after i started the excitation! todo

        ! at the start of excittaion also certain values.

        ! need for second sum term the occupation of the excited state... ->
        ! need generator type considerations.. todo
        ! at i the occupation of the excited state is necessary ->
        ! which. independently means

        ! did some stupid double counting down there...
        ! still something wrong down there...
        ! i should be able to formulate that in terms of st and en..
        integral = integral + get_umat_el(i, i, j, i) * csf_i%Occ_real(i)
        integral = integral + get_umat_el(i, j, j, j) * (csf_i%Occ_real(j) - 1.0_dp)

        ! have to reset prod here!!!

        ! also do the same loop on the orbitals above the fullEnd. to get
        ! the double contribution
        do iO = en + 1, nSpatOrbs
            ! do stuff
            if (csf_i%stepvector(iO) == 0) cycle

            integral = integral + get_umat_el(i, iO, j, iO) * csf_i%Occ_real(iO)

            integral = integral - get_umat_el(i, iO, iO, j) * csf_i%Occ_real(iO) / 2.0_dp

            ! not necessary to do it for d = 3 or b = 1, d=1 end value! since
            ! top matrix element 0 in this case

            if (csf_i%stepvector(iO) == 3 .or. (csf_i%B_int(iO) == 1 &
                                                  .and. csf_i%stepvector(iO) == 1)) cycle

            ! have to reset prod every loop
            prod = 1.0_dp

            do jO = en + 1, iO - 1
                ! do stuff
                step = csf_i%stepvector(jO)
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, csf_i%B_real(jO), &
                                            1.0_dp, x1_element=tempWeight)

                prod = prod * tempWeight

            end do
            ! have to seperately access the top most mixed full-stop
            step = csf_i%stepvector(iO)
            call getMixedFullStop(step, step, 0, csf_i%B_real(iO), &
                                  x1_element=tempWeight)

            prod = prod * tempWeight

            prod = prod * topCont

            integral = integral + get_umat_el(i, iO, iO, j) * prod
        end do

    end subroutine calc_integral_contribution_single

    pure subroutine calc_mixed_start_contr_integral(ilut, csf_i, t, excitInfo, &
            integral, rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in), value :: excitInfo
        HElement_t(dp), intent(out) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_mixed_start_contr_integral"
        integer :: st, se, en, elecInd, holeInd, sw, step, i, max_num_rdm
        integer :: rdm_count
        real(dp) :: bot_cont
        logical :: below_flag, rdm_flag
        real(dp) :: mat_ele, start_mat, stay_mat
        integer(int_rdm), allocatable :: tmp_rdm_ind(:)
        real(dp), allocatable :: tmp_rdm_mat(:)

        if (present(rdm_ind) .or. present(rdm_mat)) then
            ASSERT(present(rdm_ind) .and. present(rdm_mat))
            rdm_flag = .true.
        else
            rdm_flag = .false.
        end if

        ! whats different here?? what do i have to consider? and how to optimize?
        ! to make it most similar to the full-start into full-stop calc.
        ! i could loop from the first switch downwards and stop at
        ! a d = 1, b = 1 stepvalue and definetly unify pgen and integral
        ! calculation!
        ! to similary reuse the already calculated quantities loop from
        ! switch to start to 1
        st = excitInfo%fullStart
        se = excitInfo%firstEnd
        en = excitInfo%fullEnd
        ! depending on the type of excitaiton, calculation of orbital pgens
        ! change
        if (excitInfo%typ == excit_type%fullStart_L_to_R) then
            elecInd = en
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstart_R_to_L) then
            elecInd = se
            holeInd = en
        else
            call stop_all(this_routine, "should not be here!")
        end if

        sw = findFirstSwitch(ilut, t, st, se)

        if (rdm_flag) then
            max_num_rdm = sw
            allocate(tmp_rdm_ind(max_num_rdm), source=0_int_rdm)
            allocate(tmp_rdm_mat(max_num_rdm), source=0.0_dp)
            rdm_count = 0
        end if

        ! what can i precalculate beforehand?
        step = csf_i%stepvector(st)

        integral = h_cast(0.0_dp)

        if (step == 1) then

            if (isOne(t, st)) then

                bot_cont = Root2 * sqrt((csf_i%B_real(st) - 1.0_dp) / &
                                        (csf_i%B_real(st) + 1.0_dp))

            else

                bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) - 1.0_dp) * &
                                           (csf_i%B_real(st) + 1.0_dp)))

            end if
        else

            if (isOne(t, st)) then
                bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) + 1.0_dp) * &
                                           (csf_i%B_real(st) + 3.0_dp)))

            else
                bot_cont = -Root2 * sqrt((csf_i%B_real(st) + 3.0_dp) / &
                                         (csf_i%B_real(st) + 1.0_dp))
            end if
        end if

        ! loop from start backwards so i can abort at a d=1 & b=1 stepvalue
        ! also consider if bot_cont < EPS to avoid unnecarry calculations
        if (.not. near_zero(bot_cont)) then

            mat_ele = 1.0_dp
            below_flag = .false.

            do i = st - 1, 1, -1
                if (csf_i%Occ_int(i) /= 1) cycle

                ! then check if thats the last stepvalue to consider
                if (csf_i%stepvector(i) == 1 .and. csf_i%B_int(i) == 1) then
                    below_flag = .true.
                end if

                ! then deal with the matrix element and branching probabilities
                step = csf_i%stepvector(i)

                ! get both start and staying matrix elements -> and update
                ! matrix element contributions on the fly to avoid second loop!
                call getDoubleMatrixElement(step, step, -1, gen_type%R, gen_type%L, &
                    csf_i%B_real(i), 1.0_dp, x1_element=start_mat)

                call getDoubleMatrixElement(step, step, 0, gen_type%R, gen_type%L, &
                    csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)

                ! another check.. although this should not happen
                ! except the other d = 1 & b = 1 condition is already met
                ! above, to not continue:
                if (near_zero(stay_mat)) below_flag = .true.

                ! i think i could avoid the second loop over j
                ! if i express everything in terms of already calculated
                ! quantities!
                ! "normally" matrix element shouldnt be 0 anymore... still check
                if (.not. near_zero(start_mat)) then
                    integral = integral + start_mat * mat_ele * &
                        (get_umat_el(i, holeInd, elecInd, i) &
                       + get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(i, elecInd, holeInd, i)
                        tmp_rdm_mat(rdm_count) = start_mat * mat_ele * bot_cont
                    end if
                end if

                if (below_flag) exit

                ! also update matrix element on the fly
                mat_ele = stay_mat * mat_ele

            end do

            ! and update matrix element finally with bottom contribution
            integral = integral * bot_cont

        end if

        ! start to switch loop: here matrix elements are not 0!
        ! and its only db = 0 branch and no stepvalue change!
        ! if the start is the switch nothing happens

        step = csf_i%stepvector(st)

        ! calculate the necarry values needed to formulate everything in terms
        ! of the already calculated quantities:
        call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, &
            csf_i%B_real(st), 1.0_dp, x1_element=mat_ele)

        ! and calc. x1^-1
        ! keep tempWweight as the running matrix element which gets updated
        ! every iteration

        ! for rdms (in this current setup) I need to make a dummy
        ! output if sw == st)
        if (rdm_flag .and. sw == st) then
            rdm_count = rdm_count + 1
            tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(sw, elecInd, holeInd, sw)
            tmp_rdm_mat(rdm_count) = 1.0_dp
        end if

        if (.not. near_zero(abs(mat_ele))) then

            mat_ele = 1.0_dp / mat_ele

            do i = st + 1, sw - 1
                ! the good thing here is, i do not need to loop a second time,
                ! since i can recalc. the matrix elements and pgens on-the fly
                ! here the matrix elements should not be 0 or otherwise the
                ! excitation wouldnt have happended anyways
                if (csf_i%Occ_int(i) /= 1) cycle

                step = csf_i%stepvector(i)

                ! update inverse product
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, &
                    csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)

                ASSERT(.not. near_zero(stay_mat))

                mat_ele = mat_ele / stay_mat

                ! and also get starting contribution
                call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, &
                    csf_i%B_real(i), 1.0_dp, x1_element=start_mat)

                ! because the rest of the matrix element is still the same in
                ! both cases...
                if (.not. near_zero(start_mat)) then
                    integral = integral + mat_ele * start_mat * &
                        (get_umat_el(holeInd, i, i, elecInd) + &
                         get_umat_el(i, holeInd, elecInd, i)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(i, elecInd, holeInd, i)
                        tmp_rdm_mat(rdm_count) = start_mat * mat_ele
                    end if
                end if
            end do

            ! handle switch seperately (but only if switch > start)
            if (sw > st) then

                step = csf_i%stepvector(sw)
                ! on the switch the original probability is:
                if (step == 1) then
                    call getDoubleMatrixElement(2, 1, 0, gen_type%L, gen_type%R, &
                        csf_i%B_real(sw), 1.0_dp, x1_element=stay_mat)

                    call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, &
                        csf_i%B_real(sw), 1.0_dp, x1_element=start_mat)

                else
                    call getDoubleMatrixElement(1, 2, 0, gen_type%L, gen_type%R, &
                        csf_i%B_real(sw), 1.0_dp, x1_element=stay_mat)

                    call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, &
                        csf_i%B_real(sw), 1.0_dp, x1_element=start_mat)

                end if

                ! update inverse product
                ! and also get starting contribution
                ASSERT(.not. near_zero(stay_mat))

                mat_ele = mat_ele * start_mat / stay_mat

                ! because the rest of the matrix element is still the same in
                ! both cases...
                integral = integral + mat_ele * (get_umat_el(holeInd, sw, sw, elecInd) + &
                                                 get_umat_el(sw, holeInd, elecInd, sw)) / 2.0_dp

                if (rdm_flag) then
                    rdm_count = rdm_count + 1
                    tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(sw, elecInd, holeInd, sw)
                    tmp_rdm_mat(rdm_count) = mat_ele
                end if
            end if
        end if

        if (present(rdm_mat)) then
            allocate(rdm_ind(rdm_count), source=tmp_rdm_ind(1:rdm_count))
            allocate(rdm_mat(rdm_count), source=tmp_rdm_mat(1:rdm_count))

            deallocate(tmp_rdm_ind)
            deallocate(tmp_rdm_mat)
        end if

    end subroutine calc_mixed_start_contr_integral

    subroutine calc_mixed_start_contr_pgen(ilut, csf_i, t, excitInfo, branch_pgen, &
            pgen)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), value :: excitInfo
        real(dp), value :: branch_pgen
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "calc_mixed_start_contr_pgen"
        integer :: st, se, en, elecInd, holeInd, sw, step, i
        real(dp) :: orb_pgen, zero_weight, switch_weight, stay_weight
        real(dp) :: bot_cont, posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        real(dp) :: new_pgen, start_mat, stay_mat, start_weight
        type(WeightObj_t) :: weights
        logical :: below_flag

        st = excitInfo%fullStart
        se = excitInfo%firstEnd
        en = excitInfo%fullEnd
        ! depending on the type of excitaiton, calculation of orbital pgens
        ! change
        if (excitInfo%typ == excit_type%fullStart_L_to_R) then
            elecInd = en
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstart_R_to_L) then
            elecInd = se
            holeInd = en
        else
            call stop_all(this_routine, "should not be here!")
        end if

        sw = findFirstSwitch(ilut, t, st, se)

        ! what can i precalculate beforehand?
        step = csf_i%stepvector(st)

        ! do i actually deal with the actual start orbital influence??
        ! fuck i don't think so.. wtf..
        call calc_orbital_pgen_contrib_start(csf_i, [2 * st, 2 * elecInd], &
            holeInd, orb_pgen)

        pgen = orb_pgen * branch_pgen

        ! since weights only depend on the number of switches at the
        ! semistop and semistop and full-end index i can calculate
        ! it beforehand for all?
        excitInfo%fullStart = 1
        excitInfo%secondStart = 1
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitches, negSwitches)

        weights = init_fullStartWeight(csf_i, se, en, negSwitches(se), &
                                       posSwitches(se), csf_i%B_real(se))

        ! determine the original starting weight
        zero_weight = weights%proc%zero(negSwitches(st), posSwitches(st), &
                                        csf_i%B_real(st), weights%dat)

        if (step == 1) then

            switch_weight = weights%proc%plus(posSwitches(st), csf_i%B_real(st), &
                                              weights%dat)

            if (isOne(t, st)) then

                bot_cont = Root2 * sqrt((csf_i%B_real(st) - 1.0_dp) / &
                                        (csf_i%B_real(st) + 1.0_dp))

                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(st))

                start_weight = zero_weight / (zero_weight + switch_weight)

            else

                bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) - 1.0_dp) * &
                                           (csf_i%B_real(st) + 1.0_dp)))

                stay_weight = 1.0_dp - calcStayingProb(zero_weight, switch_weight, &
                                                       csf_i%B_real(st))

                start_weight = switch_weight / (zero_weight + switch_weight)

            end if
        else

            switch_weight = weights%proc%minus(negSwitches(st), &
                                               csf_i%B_real(st), weights%dat)

            if (isOne(t, st)) then
                bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) + 1.0_dp) * &
                                           (csf_i%B_real(st) + 3.0_dp)))

                stay_weight = 1.0_dp - calcStayingProb(zero_weight, switch_weight, &
                                                       csf_i%B_real(st))

                start_weight = switch_weight / (zero_weight + switch_weight)

            else
                bot_cont = -Root2 * sqrt((csf_i%B_real(st) + 3.0_dp) / &
                                         (csf_i%B_real(st) + 1.0_dp))

                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(st))

                start_weight = zero_weight / (zero_weight + switch_weight)
            end if
        end if

        ASSERT(.not. near_zero(start_weight))
        ! update the pgen stumbs here to reuse start_weight variable
        new_pgen = stay_weight * branch_pgen / start_weight

        ! divide out the original starting weight:
        branch_pgen = branch_pgen / start_weight

        ! loop from start backwards so i can abort at a d=1 & b=1 stepvalue
        ! also consider if bot_cont < EPS to avoid unnecarry calculations
        if (.not. near_zero(bot_cont)) then

            below_flag = .false.

            do i = st - 1, 1, -1
                if (csf_i%Occ_int(i) /= 1) cycle

                ! then check if thats the last stepvalue to consider
                if (csf_i%stepvector(i) == 1 .and. csf_i%B_int(i) == 1) then
                    below_flag = .true.
                end if

                ! then i need to calculate the orbital probability
                ! from the fact that this is a lowering into raising fullstart
                ! i know, more about the restrictions...
                ! and that fullend is the electron eg.
                ! depening on the type of excitation (r2l or l2r) the electron
                ! orbitals change here
                call calc_orbital_pgen_contrib_start(csf_i, [2*i, 2*elecInd], &
                    holeInd, orb_pgen)

                ! then deal with the matrix element and branching probabilities
                step = csf_i%stepvector(i)

                ! get both start and staying matrix elements -> and update
                ! matrix element contributions on the fly to avoid second loop!
                call getDoubleMatrixElement(step, step, -1, gen_type%R, gen_type%L, &
                    csf_i%B_real(i), 1.0_dp, x1_element=start_mat)

                call getDoubleMatrixElement(step, step, 0, gen_type%R, gen_type%L, &
                    csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)

                ! another check.. although this should not happen
                ! except the other d = 1 & b = 1 condition is already met
                ! above, to not continue:
                if (near_zero(stay_mat)) below_flag = .true.

                zero_weight = weights%proc%zero(negSwitches(i), posSwitches(i), &
                    csf_i%B_real(i), weights%dat)

                if (step == 1) then
                    switch_weight = weights%proc%plus(posSwitches(i), &
                                                      csf_i%B_real(i), weights%dat)
                else
                    switch_weight = weights%proc%minus(negSwitches(i), &
                                                       csf_i%B_real(i), weights%dat)
                end if

                start_weight = zero_weight / (zero_weight + switch_weight)
                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(i))

                ! i think i could avoid the second loop over j
                ! if i express everything in terms of already calculated
                ! quantities!
                ! "normally" matrix element shouldnt be 0 anymore... still check
                if (.not. near_zero(start_mat)) then
                    pgen = pgen + orb_pgen * start_weight * new_pgen
                end if

                if (below_flag) exit

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                    if (new_pgen < trunc_guga_pgen) then
                        new_pgen = 0.0_dp
                    end if
                end if

                ! update new_pgen for next cycle
                new_pgen = stay_weight * new_pgen

            end do

        end if

        ! start to switch loop: here matrix elements are not 0!
        ! and its only db = 0 branch and no stepvalue change!
        ! if the start is the switch nothing happens

        step = csf_i%stepvector(st)

        ! calculate the necarry values needed to formulate everything in terms
        ! of the already calculated quantities:
        call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, &
            csf_i%B_real(st), 1.0_dp, x1_element=start_mat)

        if (.not. near_zero(abs(start_mat))) then

            do i = st + 1, sw - 1
                ! the good thing here is, i do not need to loop a second time,
                ! since i can recalc. the matrix elements and pgens on-the fly
                ! here the matrix elements should not be 0 or otherwise the
                ! excitation wouldnt have happended anyways
                if (csf_i%Occ_int(i) /= 1) cycle

                ! calculate orbitals pgen first and cycle if 0
                call calc_orbital_pgen_contrib_start(csf_i, [2 * i, 2 * elecInd], &
                    holeInd, orb_pgen)

                step = csf_i%stepvector(i)

                ! update inverse product
#ifdef DEBUG_
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, &
                    csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)
                ASSERT(.not. near_zero(stay_mat))
#endif

                ! and update pgens also
                zero_weight = weights%proc%zero(negSwitches(i), posSwitches(i), &
                    csf_i%B_real(i), weights%dat)

                if (step == 1) then
                    switch_weight = weights%proc%plus(posSwitches(i), &
                                                csf_i%B_real(i), weights%dat)
                else
                    switch_weight = weights%proc%minus(negSwitches(i), &
                                                csf_i%B_real(i), weights%dat)
                end if

                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(i))

                start_weight = zero_weight / (zero_weight + switch_weight)

                branch_pgen = branch_pgen / stay_weight

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then

                    if (branch_pgen * start_weight > trunc_guga_pgen) then
                        pgen = pgen + orb_pgen * branch_pgen * start_weight
                    end if
                else
                    ! and add up correctly
                    pgen = pgen + orb_pgen * branch_pgen * start_weight
                end if

            end do

            ! handle switch seperately (but only if switch > start)
            if (sw > st) then

                ! check orb_pgen otherwise no influencce
                call calc_orbital_pgen_contrib_start(csf_i, [2 * sw, 2 * elecInd], &
                    holeInd, orb_pgen)

                if (.not. near_zero(orb_pgen)) then

                    step = csf_i%stepvector(sw)

                    zero_weight = weights%proc%zero(negSwitches(sw), posSwitches(sw), &
                                                    csf_i%B_real(sw), weights%dat)

                    ! on the switch the original probability is:
                    if (step == 1) then
                        switch_weight = weights%proc%plus(posSwitches(sw), &
                                                csf_i%B_real(sw), weights%dat)

#ifdef DEBUG_
                        call getDoubleMatrixElement(2, 1, 0, gen_type%L, gen_type%R, &
                            csf_i%B_real(sw), 1.0_dp, x1_element=stay_mat)
#endif

                    else
                        switch_weight = weights%proc%minus(negSwitches(sw), &
                                                csf_i%B_real(sw), weights%dat)

#ifdef DEBUG_
                        call getDoubleMatrixElement(1, 2, 0, gen_type%L, gen_type%R, &
                            csf_i%B_real(sw), 1.0_dp, x1_element=stay_mat)
#endif

                    end if

                    ! update inverse product
                    ! and also get starting contribution
                    ASSERT(.not. near_zero(stay_mat))

                    stay_weight = 1.0_dp - calcStayingProb(zero_weight, switch_weight, &
                                                           csf_i%B_real(sw))

                    ! and the new startProb is also the non-b=0 branch
                    start_weight = switch_weight / (zero_weight + switch_weight)

                    if (t_trunc_guga_pgen .or. &
                        (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                        if (branch_pgen * start_weight / start_weight > trunc_guga_pgen) then
                            pgen = pgen + orb_pgen * branch_pgen * start_weight / stay_weight
                        end if

                    else
                        pgen = pgen + orb_pgen * branch_pgen * start_weight / stay_weight
                    end if

                end if
            end if
        end if

        ! i also need to consider the electron pair picking probability..
        pgen = pgen / real(ElecPairs, dp)
        ! and if the second electron is in a double occupied orbital I have
        ! to modify it with 2
        if (csf_i%stepvector(elecInd) == 3) pgen = pgen * 2.0_dp

    end subroutine calc_mixed_start_contr_pgen

    subroutine calc_mixed_start_contr_sym(ilut, csf_i, t, excitInfo, branch_pgen, &
                                          pgen, integral, rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        real(dp), intent(inout) :: branch_pgen
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_mixed_start_contr_sym"

        integer :: sw, i, st, se, step, en, elecInd, holeInd
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), mat_ele, &
                    new_pgen, zero_weight, switch_weight, stay_mat, start_mat, bot_cont, &
                    orb_pgen, start_weight, stay_weight
        type(WeightObj_t) :: weights
        logical :: below_flag
        integer(int_rdm), allocatable :: tmp_rdm_ind(:)
        real(dp), allocatable :: tmp_rdm_mat(:), temp_int
        integer :: rdm_count, max_num_rdm
        logical :: rdm_flag, test_skip

        test_skip = .false.

        if (present(rdm_ind) .or. present(rdm_mat)) then
            ASSERT(present(rdm_ind) .and. present(rdm_mat))
            rdm_flag = .true.
        else
            rdm_flag = .false.
        end if

        ! whats different here?? what do i have to consider? and how to optimize?
        ! to make it most similar to the full-start into full-stop calc.
        ! i could loop from the first switch downwards and stop at
        ! a d = 1, b = 1 stepvalue and definetly unify pgen and integral
        ! calculation!
        ! to similary reuse the already calculated quantities loop from
        ! switch to start to 1
        st = excitInfo%fullStart
        se = excitInfo%firstEnd
        en = excitInfo%fullEnd
        ! depending on the type of excitaiton, calculation of orbital pgens
        ! change
        if (excitInfo%typ == excit_type%fullStart_L_to_R) then
            elecInd = en
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstart_R_to_L) then
            elecInd = se
            holeInd = en
        else
            call stop_all(this_routine, "should not be here!")
        end if

        sw = findFirstSwitch(ilut, t, st, se)

        if (rdm_flag) then
            max_num_rdm = sw
            allocate(tmp_rdm_ind(max_num_rdm), source=0_int_rdm)
            allocate(tmp_rdm_mat(max_num_rdm), source=0.0_dp)
            rdm_count = 0
        end if

        ! what can i precalculate beforehand?
        step = csf_i%stepvector(st)

        integral = h_cast(0.0_dp)

        ! do i actually deal with the actual start orbital influence??
        ! fuck i don't think so.. wtf..
        call calc_orbital_pgen_contrib_start(csf_i, [2 * st, 2 * elecInd], &
            holeInd, orb_pgen)

        pgen = orb_pgen * branch_pgen

        ! since weights only depend on the number of switches at the
        ! semistop and semistop and full-end index i can calculate
        ! it beforehand for all?
        excitInfo%fullStart = 1
        excitInfo%secondStart = 1
        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitches, negSwitches)

        weights = init_fullStartWeight(csf_i, se, en, negSwitches(se), &
                                       posSwitches(se), csf_i%B_real(se))

        ! determine the original starting weight
        zero_weight = weights%proc%zero(negSwitches(st), posSwitches(st), &
                                        csf_i%B_real(st), weights%dat)

        if (step == 1) then

            switch_weight = weights%proc%plus(posSwitches(st), csf_i%B_real(st), &
                                              weights%dat)

            if (isOne(t, st)) then

                bot_cont = Root2 * sqrt((csf_i%B_real(st) - 1.0_dp) / &
                                        (csf_i%B_real(st) + 1.0_dp))

                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(st))

                start_weight = zero_weight / (zero_weight + switch_weight)

            else

                bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) - 1.0_dp) * &
                                           (csf_i%B_real(st) + 1.0_dp)))

                stay_weight = 1.0_dp - calcStayingProb(zero_weight, switch_weight, &
                                                       csf_i%B_real(st))

                start_weight = switch_weight / (zero_weight + switch_weight)

            end if
        else

            switch_weight = weights%proc%minus(negSwitches(st), &
                                               csf_i%B_real(st), weights%dat)

            if (isOne(t, st)) then
                bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) + 1.0_dp) * &
                                           (csf_i%B_real(st) + 3.0_dp)))

                stay_weight = 1.0_dp - calcStayingProb(zero_weight, switch_weight, &
                                                       csf_i%B_real(st))

                start_weight = switch_weight / (zero_weight + switch_weight)

            else
                bot_cont = -Root2 * sqrt((csf_i%B_real(st) + 3.0_dp) / &
                                         (csf_i%B_real(st) + 1.0_dp))

                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(st))

                start_weight = zero_weight / (zero_weight + switch_weight)
            end if
        end if

        ASSERT(.not. near_zero(start_weight))
        ! update the pgen stumbs here to reuse start_weight variable
        new_pgen = stay_weight * branch_pgen / start_weight

        ! divide out the original starting weight:
        branch_pgen = branch_pgen / start_weight

        ! loop from start backwards so i can abort at a d=1 & b=1 stepvalue
        ! also consider if bot_cont < EPS to avoid unnecarry calculations
        if (.not. near_zero(bot_cont)) then

            mat_ele = 1.0_dp
            below_flag = .false.

            do i = st - 1, 1, -1
                if (csf_i%Occ_int(i) /= 1) cycle

                ! then check if thats the last stepvalue to consider
                if (csf_i%stepvector(i) == 1 .and. csf_i%B_int(i) == 1) then
                    below_flag = .true.
                end if

                ! then i need to calculate the orbital probability
                ! from the fact that this is a lowering into raising fullstart
                ! i know, more about the restrictions...
                ! and that fullend is the electron eg.
                ! depening on the type of excitation (r2l or l2r) the electron
                ! orbitals change here
                call calc_orbital_pgen_contrib_start(csf_i, [2*i, 2*elecInd], &
                    holeInd, orb_pgen)

                ! then deal with the matrix element and branching probabilities
                step = csf_i%stepvector(i)

                ! get both start and staying matrix elements -> and update
                ! matrix element contributions on the fly to avoid second loop!
                call getDoubleMatrixElement(step, step, -1, gen_type%R, gen_type%L, &
                    csf_i%B_real(i), 1.0_dp, x1_element=start_mat)

                call getDoubleMatrixElement(step, step, 0, gen_type%R, gen_type%L, &
                    csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)

                ! another check.. although this should not happen
                ! except the other d = 1 & b = 1 condition is already met
                ! above, to not continue:
                if (near_zero(stay_mat)) below_flag = .true.

                ! check if orb_pgen is non-zero
                if (.not. test_skip) then
                    if (near_zero(orb_pgen) .and. (.not. rdm_flag)) then
                        temp_int = (get_umat_el(i, holeInd, elecInd, i) &
                            + get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp
                        if (.not. near_zero(temp_int) .and. near_zero(orb_pgen)) then
                            print *, "start 1 orb_pgen 0, integral: ", temp_int
                        end if
                        ! still have to update matrix element, even if 0 pgen
                        mat_ele = mat_ele * stay_mat

                        cycle
                    end if
                end if

                zero_weight = weights%proc%zero(negSwitches(i), posSwitches(i), &
                    csf_i%B_real(i), weights%dat)

                if (step == 1) then
                    switch_weight = weights%proc%plus(posSwitches(i), &
                                                      csf_i%B_real(i), weights%dat)
                else
                    switch_weight = weights%proc%minus(negSwitches(i), &
                                                       csf_i%B_real(i), weights%dat)
                end if

                start_weight = zero_weight / (zero_weight + switch_weight)
                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(i))

                ! i think i could avoid the second loop over j
                ! if i express everything in terms of already calculated
                ! quantities!
                ! "normally" matrix element shouldnt be 0 anymore... still check
                if (.not. near_zero(start_mat)) then
                    integral = integral + start_mat * mat_ele * &
                        (get_umat_el(i, holeInd, elecInd, i) &
                       + get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(i, elecInd, holeInd, i)
                        tmp_rdm_mat(rdm_count) = start_mat * mat_ele * bot_cont
                    end if

                    if (t_trunc_guga_pgen .or. &
                        (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                        if (new_pgen < trunc_guga_pgen) then
                            new_pgen = 0.0_dp
                        end if
                    end if

                    pgen = pgen + orb_pgen * start_weight * new_pgen
                end if

                if (below_flag) exit

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                    if (new_pgen < trunc_guga_pgen) then
                        new_pgen = 0.0_dp
                    end if
                end if

                if (test_skip) then
                    if (near_zero(orb_pgen)) then
                        print *, "stay_weight(1): ", stay_weight
                        print *, "i, st: ", i, st
                    end if
                end if
                ! update new_pgen for next cycle
                new_pgen = stay_weight * new_pgen

                ! also update matrix element on the fly
                mat_ele = stay_mat * mat_ele

            end do

            ! and update matrix element finally with bottom contribution
            integral = integral * bot_cont

        end if

        ! start to switch loop: here matrix elements are not 0!
        ! and its only db = 0 branch and no stepvalue change!
        ! if the start is the switch nothing happens

        step = csf_i%stepvector(st)

        ! calculate the necarry values needed to formulate everything in terms
        ! of the already calculated quantities:
        call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, &
            csf_i%B_real(st), 1.0_dp, x1_element=mat_ele)

        ! and calc. x1^-1
        ! keep tempWweight as the running matrix element which gets updated
        ! every iteration

        ! for rdms (in this current setup) I need to make a dummy
        ! output if sw == st)
        if (rdm_flag .and. sw == st) then
            rdm_count = rdm_count + 1
            tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(sw, elecInd, holeInd, sw)
            tmp_rdm_mat(rdm_count) = 1.0_dp
        end if

        if (.not. near_zero(abs(mat_ele))) then

            mat_ele = 1.0_dp / mat_ele

            do i = st + 1, sw - 1
                ! the good thing here is, i do not need to loop a second time,
                ! since i can recalc. the matrix elements and pgens on-the fly
                ! here the matrix elements should not be 0 or otherwise the
                ! excitation wouldnt have happended anyways
                if (csf_i%Occ_int(i) /= 1) cycle

                ! calculate orbitals pgen first and cycle if 0
                call calc_orbital_pgen_contrib_start(csf_i, [2 * i, 2 * elecInd], &
                    holeInd, orb_pgen)

                step = csf_i%stepvector(i)

                ! update inverse product
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, &
                    csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)

                ASSERT(.not. near_zero(stay_mat))

                mat_ele = mat_ele / stay_mat

                ! check if orb_pgen is non-zero
                ! still have to update matrix element in this case..
                ! so do the cycle only afterwards..

                if (.not. test_skip) then
                    if (near_zero(orb_pgen) .and. (.not. rdm_flag)) then
                        temp_int = (get_umat_el(i, holeInd, elecInd, i) &
                                  + get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp
                        if (.not. near_zero(temp_int) .and. near_zero(orb_pgen)) then
                            print *, "start 2 orb_pgen 0, integral: ", temp_int
                        end if

                        cycle
                    end if
                end if

                ! and update pgens also
                zero_weight = weights%proc%zero(negSwitches(i), posSwitches(i), &
                    csf_i%B_real(i), weights%dat)

                if (step == 1) then
                    switch_weight = weights%proc%plus(posSwitches(i), &
                                                csf_i%B_real(i), weights%dat)
                else
                    switch_weight = weights%proc%minus(negSwitches(i), &
                                                csf_i%B_real(i), weights%dat)
                end if

                stay_weight = calcStayingProb(zero_weight, switch_weight, &
                                              csf_i%B_real(i))

                start_weight = zero_weight / (zero_weight + switch_weight)


                ! and also get starting contribution
                call getDoubleMatrixElement(step, step, -1, gen_type%L, gen_type%R, &
                    csf_i%B_real(i), 1.0_dp, x1_element=start_mat)

                ! because the rest of the matrix element is still the same in
                ! both cases...
                if (.not. near_zero(start_mat)) then
                    integral = integral + mat_ele * start_mat * &
                        (get_umat_el(holeInd, i, i, elecInd) + &
                         get_umat_el(i, holeInd, elecInd, i)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(i, elecInd, holeInd, i)
                        tmp_rdm_mat(rdm_count) = start_mat * mat_ele
                    end if

                end if

                ! and update probWeight
                ! i think i should not put in the intermediate start_weight
                ! permanently..

                if (test_skip) then
                    if (near_zero(orb_pgen)) then
                        print *, "stay_weight(2): ", stay_weight
                        print *, "i, st, sw:", i, st, sw
                    end if
                end if

                branch_pgen = branch_pgen / stay_weight

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then

                    if (branch_pgen * start_weight > trunc_guga_pgen) then
                        pgen = pgen + orb_pgen * branch_pgen * start_weight
                    end if
                else
                    ! and add up correctly
                    pgen = pgen + orb_pgen * branch_pgen * start_weight
                end if

            end do

            ! handle switch seperately (but only if switch > start)
            if (sw > st) then

                ! check orb_pgen otherwise no influencce
                call calc_orbital_pgen_contrib_start(csf_i, [2 * sw, 2 * elecInd], &
                    holeInd, orb_pgen)

                temp_int = (get_umat_el(sw, holeInd, elecInd, sw) &
                          + get_umat_el(holeInd, sw, sw, elecInd)) / 2.0_dp

                if (near_zero(orb_pgen) .and. .not. near_zero(temp_int)) then
                    print *, "start 3 orb_pgen 0, integral: ", temp_int
                end if

                if (.not. near_zero(orb_pgen) .or. rdm_flag .or. test_skip) then

                    step = csf_i%stepvector(sw)

                    zero_weight = weights%proc%zero(negSwitches(sw), posSwitches(sw), &
                                                    csf_i%B_real(sw), weights%dat)

                    ! on the switch the original probability is:
                    if (step == 1) then
                        switch_weight = weights%proc%plus(posSwitches(sw), &
                                                csf_i%B_real(sw), weights%dat)

                        call getDoubleMatrixElement(2, 1, 0, gen_type%L, gen_type%R, &
                            csf_i%B_real(sw), 1.0_dp, x1_element=stay_mat)

                        call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, &
                            csf_i%B_real(sw), 1.0_dp, x1_element=start_mat)

                    else
                        switch_weight = weights%proc%minus(negSwitches(sw), &
                                                csf_i%B_real(sw), weights%dat)

                        call getDoubleMatrixElement(1, 2, 0, gen_type%L, gen_type%R, &
                            csf_i%B_real(sw), 1.0_dp, x1_element=stay_mat)

                        call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, &
                            csf_i%B_real(sw), 1.0_dp, x1_element=start_mat)

                    end if

                    ! update inverse product
                    ! and also get starting contribution
                    ASSERT(.not. near_zero(stay_mat))

                    mat_ele = mat_ele * start_mat / stay_mat

                    ! because the rest of the matrix element is still the same in
                    ! both cases...
                    integral = integral + mat_ele * (get_umat_el(holeInd, sw, sw, elecInd) + &
                                                     get_umat_el(sw, holeInd, elecInd, sw)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(sw, elecInd, holeInd, sw)
                        tmp_rdm_mat(rdm_count) = mat_ele
                    end if

                    stay_weight = 1.0_dp - calcStayingProb(zero_weight, switch_weight, &
                                                           csf_i%B_real(sw))

                    ! and the new startProb is also the non-b=0 branch
                    start_weight = switch_weight / (zero_weight + switch_weight)

                    if (t_trunc_guga_pgen .or. &
                        (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                        if (branch_pgen * start_weight / start_weight > trunc_guga_pgen) then
                            pgen = pgen + orb_pgen * branch_pgen * start_weight / stay_weight
                        end if

                    else
                        pgen = pgen + orb_pgen * branch_pgen * start_weight / stay_weight
                    end if

                end if
            end if
        end if

        ! i also need to consider the electron pair picking probability..
        pgen = pgen / real(ElecPairs, dp)
        ! and if the second electron is in a double occupied orbital I have
        ! to modify it with 2
        if (csf_i%stepvector(elecInd) == 3) pgen = pgen * 2.0_dp

        if (present(rdm_mat)) then
            allocate(rdm_ind(rdm_count), source=tmp_rdm_ind(1:rdm_count))
            allocate(rdm_mat(rdm_count), source=tmp_rdm_mat(1:rdm_count))

            deallocate(tmp_rdm_ind)
            deallocate(tmp_rdm_mat)
        end if

    end subroutine calc_mixed_start_contr_sym

    subroutine calc_mixed_end_contr_sym(ilut, csf_i, t, excitInfo, branch_pgen, pgen, &
                                        integral, rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        real(dp), intent(inout) :: branch_pgen
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_mixed_end_contr_sym"

        integer :: st, se, en, step, sw, elecInd, holeInd, i, j
        real(dp) :: top_cont, mat_ele, stay_mat, end_mat, orb_pgen, new_pgen, &
                    posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), &
                    tmp_pos(nSpatOrbs), tmp_neg(nSpatOrbs)
        logical :: above_flag
        type(BranchWeightArr_t) :: weight_funcs(nSpatOrbs)
        type(WeightObj_t) :: weights

        integer(int_rdm), allocatable :: tmp_rdm_ind(:)
        real(dp), allocatable :: tmp_rdm_mat(:)
        logical :: rdm_flag, test_skip
        integer :: rdm_count, max_num_rdm

        test_skip = .false.

        if (present(rdm_ind) .or. present(rdm_mat)) then
            ASSERT(present(rdm_ind))
            ASSERT(present(rdm_mat))
            rdm_flag = .true.
        else
            rdm_flag = .false.
        end if

        ! do as much stuff as possible beforehand
        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd
        if (excitInfo%typ == excit_type%fullstop_L_to_R) then
            elecInd = st
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstop_R_to_L) then
            elecInd = se
            holeInd = st
        else
            call stop_all(this_routine, "should not be here!")
        end if

        integral = h_cast(0.0_dp)
        ! also here i didn't consider the actual end contribution or? ...
        call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * en], &
            holeInd, orb_pgen)

        pgen = orb_pgen * branch_pgen

        step = csf_i%stepvector(en)

        sw = findLastSwitch(ilut, t, se, en)

        if (rdm_flag) then
            max_num_rdm = (nSpatOrbs - sw + 1)
            allocate(tmp_rdm_ind(max_num_rdm), source=0_int_rdm)
            allocate(tmp_rdm_mat(max_num_rdm), source=0.0_dp)
            rdm_count = 0
        end if

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitches, negSwitches)

        ! need temporary switch arrays for more efficiently recalcing
        ! weights
        tmp_pos = posSwitches
        tmp_neg = negSwitches
        ! after last switch only dB = 0 branches! consider that
        call setup_weight_funcs(t, csf_i, st, se, weight_funcs)

        if (en < nSpatOrbs) then
            select case (step)
            case (1)
                if (isOne(t, en)) then
                    top_cont = -Root2 * sqrt((csf_i%B_real(en) + 2.0_dp) / &
                                             csf_i%B_real(en))

                else
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                end if
            case (2)
                if (isOne(t, en)) then
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                else
                    top_cont = Root2 * sqrt(csf_i%B_real(en) / &
                                            (csf_i%B_real(en) + 2.0_dp))
                end if

            case default
                call stop_all(this_routine, "wrong stepvalues!")

            end select

            if (.not. near_zero(top_cont)) then

                above_flag = .false.
                mat_ele = 1.0_dp

                ! to avoid to recalc. remaining switches all the time
                ! just increment them correctly
                if (step == 1) then
                    tmp_neg(se:en - 1) = tmp_neg(se:en - 1) + 1.0_dp
                else
                    tmp_pos(se:en - 1) = tmp_pos(se:en - 1) + 1.0_dp
                end if

                do i = en + 1, nSpatOrbs
                    if (csf_i%Occ_int(i) /= 1) cycle

                    ! then check if thats the last step
                    if (csf_i%stepvector(i) == 2 .and. csf_i%B_int(i) == 0) then
                        above_flag = .true.
                    end if

                    ! then calc. orbital probability
                    call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * i], &
                        holeInd, orb_pgen)

                    ! should be able to do that without second loop too!
                    ! figure out!
                    step = csf_i%stepvector(i)

                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, csf_i%B_real(i), &
                                                1.0_dp, x1_element=stay_mat)

                    call getMixedFullStop(step, step, 0, csf_i%B_real(i), &
                                          x1_element=end_mat)

                    if (.not. test_skip) then
                        if (near_zero(orb_pgen) .and. (.not. rdm_flag)) then
                            ! still have to update the switches before cycling
                            ! update the switches
                            if (csf_i%stepvector(i) == 1) then
                                tmp_neg(se:i - 1) = tmp_neg(se:i - 1) + 1.0_dp
                            else
                                tmp_pos(se:i - 1) = tmp_pos(se:i - 1) + 1.0_dp
                            end if

                            ! also have to update the matrix element, even if
                            ! the orb pgen is 0
                            mat_ele = mat_ele * stay_mat

                            cycle
                        end if
                    end if

                    ! this check should never be true, but just to be sure
                    if (near_zero(stay_mat)) above_flag = .true.

                    if (.not. near_zero(end_mat)) then
                        integral = integral + end_mat * mat_ele * &
                                   (get_umat_el(i, holeInd, elecInd, i) + &
                                    get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp

                        if (rdm_flag) then
                            rdm_count = rdm_count + 1
                            tmp_rdm_ind(rdm_count) = &
                                contract_2_rdm_ind(i, elecInd, holeInd, i)
                            tmp_rdm_mat(rdm_count) = top_cont * end_mat * mat_ele
                        end if

                        ! also only recalc. pgen if matrix element is not 0
                        excitInfo%fullEnd = i
                        excitInfo%firstEnd = i

                        weights = init_semiStartWeight(csf_i, se, i, tmp_neg(se), &
                                                       tmp_pos(se), csf_i%B_real(se))

                        new_pgen = 1.0_dp

                        ! deal with the start and semi-start seperately
                        if (csf_i%Occ_int(st) /= 1) then
                            new_pgen = new_pgen * weight_funcs(st)%ptr(weights, &
                                                                       csf_i%B_real(st), tmp_neg(st), tmp_pos(st))
                        end if

                        do j = st + 1, se - 1
                            ! can and do i have to cycle here if its not
                            ! singly occupied??
                            if (csf_i%Occ_int(j) /= 1) cycle

                            new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                                                      csf_i%B_real(j), tmp_neg(j), tmp_pos(j))
                        end do

                        ! then need to reinit double weight
                        weights = weights%ptr

                        ! and also with the semi-start
                        if (csf_i%Occ_int(se) /= 1) then
                            new_pgen = new_pgen * weight_funcs(se)%ptr(weights, &
                                                                       csf_i%B_real(se), tmp_neg(se), tmp_pos(se))
                        end if

                        do j = se + 1, i - 1
                            if (csf_i%Occ_int(j) /= 1) cycle

                            new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                                                      csf_i%B_real(j), tmp_neg(j), tmp_pos(j))
                        end do

                        if (t_trunc_guga_pgen .or. &
                            (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                            if (new_pgen < trunc_guga_pgen) then
                                new_pgen = 0.0_dp
                            end if
                        end if

                        pgen = pgen + new_pgen * orb_pgen

                    end if

                    if (above_flag) exit

                    ! otherwise update your running pgen and matrix element vars
                    mat_ele = mat_ele * stay_mat

                    ! update the switches
                    if (csf_i%stepvector(i) == 1) then
                        tmp_neg(se:i - 1) = tmp_neg(se:i - 1) + 1.0_dp
                    else
                        tmp_pos(se:i - 1) = tmp_pos(se:i - 1) + 1.0_dp
                    end if

                end do

                integral = integral * top_cont
            end if
        end if

        if (rdm_flag .and. sw == en) then
            rdm_count = rdm_count + 1
            tmp_rdm_ind(rdm_count) = &
                contract_2_rdm_ind(sw, elecInd, holeInd, sw)
            tmp_rdm_mat(rdm_count) = 1.0_dp
        end if

        if (sw < en) then

            step = csf_i%stepvector(en)

            ! inverse fullstop matrix element
            call getMixedFullStop(step, step, 0, csf_i%B_real(en), x1_element=mat_ele)

            ASSERT(.not. near_zero(mat_ele))

            mat_ele = 1.0_dp / mat_ele

            ! have to change the switches before the first cycle:
            ! but for cycling backwards, thats not so easy.. need todo

            do i = en - 1, sw + 1, -1

                if (csf_i%Occ_int(i) /= 1) cycle

                ! get orbital pgen
                call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * i], &
                    holeInd, orb_pgen)

                if (csf_i%stepvector(i) == 1) then
                    ! by looping in this direction i have to reduce
                    ! the number of switches at the beginning
                    ! but only to the left or??
                    ! i think i have to rethink that.. thats not so easy..
                    negSwitches(se:i - 1) = negSwitches(se:i - 1) - 1.0_dp

                else
                    posSwitches(se:i - 1) = posSwitches(se:i - 1) - 1.0_dp

                end if

                step = csf_i%stepvector(i)
                ! update inverse product
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, csf_i%B_real(i), &
                                            1.0_dp, x1_element=stay_mat)

                call getMixedFullStop(step, step, 0, csf_i%B_real(i), x1_element=end_mat)

                ! update matrix element
                ASSERT(.not. near_zero(stay_mat))
                mat_ele = mat_ele / stay_mat

                ! dont i still have to atleast update the matrix element
                ! even if the orbital pgen is 0??
                if (.not. test_skip) then
                    if (near_zero(orb_pgen) .and. (.not. rdm_flag)) cycle
                end if

                if (.not. near_zero(end_mat)) then

                    integral = integral + end_mat * mat_ele * &
                               (get_umat_el(i, holeInd, elecInd, i) + &
                                get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = &
                            contract_2_rdm_ind(i, elecInd, holeInd, i)
                        tmp_rdm_mat(rdm_count) = end_mat * mat_ele
                    end if

                    ! only recalc. pgen if matrix element is not 0
                    excitInfo%fullEnd = i
                    excitInfo%firstEnd = i

                    weights = init_semiStartWeight(csf_i, se, i, negSwitches(se), &
                                                   posSwitches(se), csf_i%B_real(se))

                    new_pgen = 1.0_dp

                    ! deal with the start and semi-start seperately
                    if (csf_i%Occ_int(st) /= 1) then
                        new_pgen = new_pgen * weight_funcs(st)%ptr(weights, &
                                                                   csf_i%B_real(st), negSwitches(st), posSwitches(st))
                    end if

                    do j = st + 1, se - 1
                        if (csf_i%Occ_int(j) /= 1) cycle

                        new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                                                  csf_i%B_real(j), negSwitches(j), posSwitches(j))
                    end do

                    ! then need to reinit double weight
                    weights = weights%ptr

                    ! and also with the semi-start
                    if (csf_i%Occ_int(se) /= 1) then
                        new_pgen = new_pgen * weight_funcs(se)%ptr(weights, &
                                                                   csf_i%B_real(se), negSwitches(se), posSwitches(se))
                    end if

                    do j = se + 1, i - 1
                        if (csf_i%Occ_int(j) /= 1) cycle

                        new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                                                  csf_i%B_real(j), negSwitches(j), posSwitches(j))
                    end do

                    if (t_trunc_guga_pgen .or. &
                        (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                        if (new_pgen < trunc_guga_pgen) then
                            new_pgen = 0.0_dp
                        end if
                    end if

                    pgen = pgen + new_pgen * orb_pgen

                end if

            end do

            ! deal with switch specifically:

            ! figure out orbital pgen
            call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * sw], &
                holeInd, orb_pgen)

            if (.not. near_zero(orb_pgen) .or. rdm_flag .or. test_skip) then

                step = csf_i%stepvector(sw)

                if (step == 1) then
                    ! then a -2 branch arrived!
                    call getDoubleMatrixElement(2, 1, -2, gen_type%L, gen_type%R, csf_i%B_real(sw), &
                                                1.0_dp, x1_element=stay_mat)

                    call getMixedFullStop(2, 1, -2, csf_i%B_real(sw), x1_element=end_mat)

                    ! also reduce negative switches then
                    ! only everything to the left or?
                    negSwitches(se:sw - 1) = negSwitches(se:sw - 1) - 1.0_dp

                else
                    ! +2 branch arrived!

                    call getDoubleMatrixElement(1, 2, 2, gen_type%L, gen_type%R, csf_i%B_real(sw), &
                                                1.0_dp, x1_element=stay_mat)

                    call getMixedFullStop(1, 2, 2, csf_i%B_real(sw), x1_element=end_mat)

                    ! reduce positive switchtes otherwise
                    posSwitches(se:sw - 1) = posSwitches(se:sw - 1) - 1.0_dp

                end if

                ASSERT(.not. near_zero(stay_mat))

                mat_ele = mat_ele * end_mat / stay_mat

                integral = integral + mat_ele * (get_umat_el(sw, holeInd, elecInd, sw) + &
                                                 get_umat_el(holeInd, sw, sw, elecInd)) / 2.0_dp

                if (rdm_flag) then
                    rdm_count = rdm_count + 1
                    tmp_rdm_ind(rdm_count) = &
                        contract_2_rdm_ind(sw, elecInd, holeInd, sw)
                    tmp_rdm_mat(rdm_count) = mat_ele
                end if

                ! loop to get correct pgen
                new_pgen = 1.0_dp

                weights = init_semiStartWeight(csf_i, se, sw, negSwitches(se), &
                                               posSwitches(se), csf_i%B_real(se))

                ! deal with the start and semi-start seperately
                if (csf_i%Occ_int(st) /= 1) then
                    new_pgen = new_pgen * weight_funcs(st)%ptr(weights, &
                                                               csf_i%B_real(st), negSwitches(st), posSwitches(st))
                end if

                do j = st + 1, se - 1
                    if (csf_i%Occ_int(j) /= 1) cycle

                    new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                                              csf_i%B_real(j), negSwitches(j), posSwitches(j))
                end do

                weights = weights%ptr

                ! and also with the semi-start
                if (csf_i%Occ_int(se) /= 1) then
                    new_pgen = new_pgen * weight_funcs(se)%ptr(weights, &
                                                               csf_i%B_real(se), negSwitches(se), posSwitches(se))
                end if

                do j = se + 1, sw - 1
                    if (csf_i%Occ_int(j) /= 1) cycle

                    new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                                              csf_i%B_real(j), negSwitches(j), posSwitches(j))
                end do

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                    if (new_pgen < trunc_guga_pgen) then
                        new_pgen = 0.0_dp
                    end if
                end if

                pgen = pgen + new_pgen * orb_pgen

            end if
        end if

        pgen = pgen / real(ElecPairs, dp)

        if (csf_i%stepvector(elecInd) == 3) pgen = pgen * 2.0_dp

        if (rdm_flag) then
            allocate(rdm_ind(rdm_count), source=tmp_rdm_ind(1:rdm_count))
            allocate(rdm_mat(rdm_count), source=tmp_rdm_mat(1:rdm_count))

            deallocate(tmp_rdm_ind)
            deallocate(tmp_rdm_mat)
        end if

    end subroutine calc_mixed_end_contr_sym

    subroutine calc_mixed_end_contr_pgen(ilut, csf_i, t, excitInfo, branch_pgen, &
            pgen)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), value :: excitInfo
        real(dp), value :: branch_pgen
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "calc_mixed_end_contr_pgen"

        integer :: st, se, en, step, sw, elecInd, holeInd, i, j
        real(dp) :: top_cont, stay_mat, end_mat, orb_pgen, new_pgen, &
                    posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), &
                    tmp_pos(nSpatOrbs), tmp_neg(nSpatOrbs)
        logical :: above_flag
        type(BranchWeightArr_t) :: weight_funcs(nSpatOrbs)
        type(WeightObj_t) :: weights

        ! do as much stuff as possible beforehand
        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd
        if (excitInfo%typ == excit_type%fullstop_L_to_R) then
            elecInd = st
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstop_R_to_L) then
            elecInd = se
            holeInd = st
        else
            call stop_all(this_routine, "should not be here!")
        end if

        ! also here i didn't consider the actual end contribution or? ...
        call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * en], &
            holeInd, orb_pgen)

        pgen = orb_pgen * branch_pgen

        step = csf_i%stepvector(en)

        sw = findLastSwitch(ilut, t, se, en)

        call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitches, &
            negSwitches)

        ! need temporary switch arrays for more efficiently recalcing
        ! weights
        tmp_pos = posSwitches
        tmp_neg = negSwitches
        ! after last switch only dB = 0 branches! consider that
        call setup_weight_funcs(t, csf_i, st, se, weight_funcs)

        if (en < nSpatOrbs) then
            select case (step)
            case (1)
                if (isOne(t, en)) then
                    top_cont = -Root2 * sqrt((csf_i%B_real(en) + 2.0_dp) / &
                                             csf_i%B_real(en))

                else
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                end if
            case (2)
                if (isOne(t, en)) then
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                else
                    top_cont = Root2 * sqrt(csf_i%B_real(en) / &
                                            (csf_i%B_real(en) + 2.0_dp))
                end if

            case default
                call stop_all(this_routine, "wrong stepvalues!")

            end select

            if (.not. near_zero(top_cont)) then

                above_flag = .false.

                ! to avoid to recalc. remaining switches all the time
                ! just increment them correctly
                if (step == 1) then
                    tmp_neg(se:en - 1) = tmp_neg(se:en - 1) + 1.0_dp
                else
                    tmp_pos(se:en - 1) = tmp_pos(se:en - 1) + 1.0_dp
                end if

                do i = en + 1, nSpatOrbs
                    if (csf_i%Occ_int(i) /= 1) cycle

                    ! then check if thats the last step
                    if (csf_i%stepvector(i) == 2 .and. csf_i%B_int(i) == 0) then
                        above_flag = .true.
                    end if

                    ! then calc. orbital probability
                    call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * i], &
                        holeInd, orb_pgen)

                    ! should be able to do that without second loop too!
                    ! figure out!
                    step = csf_i%stepvector(i)

                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, csf_i%B_real(i), &
                                                1.0_dp, x1_element=stay_mat)

                    call getMixedFullStop(step, step, 0, csf_i%B_real(i), &
                                          x1_element=end_mat)

                    ! this check should never be true, but just to be sure
                    if (near_zero(stay_mat)) above_flag = .true.

                    if (.not. near_zero(end_mat)) then

                        ! also only recalc. pgen if matrix element is not 0
                        excitInfo%fullEnd = i
                        excitInfo%firstEnd = i

                        weights = init_semiStartWeight(csf_i, se, i, tmp_neg(se), &
                                                       tmp_pos(se), csf_i%B_real(se))

                        new_pgen = 1.0_dp

                        ! deal with the start and semi-start seperately
                        if (csf_i%Occ_int(st) /= 1) then
                            new_pgen = new_pgen * weight_funcs(st)%ptr(weights, &
                                   csf_i%B_real(st), tmp_neg(st), tmp_pos(st))
                        end if

                        do j = st + 1, se - 1
                            ! can and do i have to cycle here if its not
                            ! singly occupied??
                            if (csf_i%Occ_int(j) /= 1) cycle

                            new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                      csf_i%B_real(j), tmp_neg(j), tmp_pos(j))
                        end do

                        ! then need to reinit double weight
                        weights = weights%ptr

                        ! and also with the semi-start
                        if (csf_i%Occ_int(se) /= 1) then
                            new_pgen = new_pgen * weight_funcs(se)%ptr(weights, &
                                   csf_i%B_real(se), tmp_neg(se), tmp_pos(se))
                        end if

                        do j = se + 1, i - 1
                            if (csf_i%Occ_int(j) /= 1) cycle

                            new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                                      csf_i%B_real(j), tmp_neg(j), tmp_pos(j))
                        end do

                        if (t_trunc_guga_pgen .or. &
                            (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                            if (new_pgen < trunc_guga_pgen) then
                                new_pgen = 0.0_dp
                            end if
                        end if

                        pgen = pgen + new_pgen * orb_pgen

                    end if

                    if (above_flag) exit

                    ! update the switches
                    if (csf_i%stepvector(i) == 1) then
                        tmp_neg(se:i - 1) = tmp_neg(se:i - 1) + 1.0_dp
                    else
                        tmp_pos(se:i - 1) = tmp_pos(se:i - 1) + 1.0_dp
                    end if
                end do
            end if
        end if

        if (sw < en) then

            step = csf_i%stepvector(en)

            ! have to change the switches before the first cycle:
            ! but for cycling backwards, thats not so easy.. need todo

            do i = en - 1, sw + 1, -1

                if (csf_i%Occ_int(i) /= 1) cycle

                ! get orbital pgen
                call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * i], &
                    holeInd, orb_pgen)

                if (csf_i%stepvector(i) == 1) then
                    ! by looping in this direction i have to reduce
                    ! the number of switches at the beginning
                    ! but only to the left or??
                    ! i think i have to rethink that.. thats not so easy..
                    negSwitches(se:i - 1) = negSwitches(se:i - 1) - 1.0_dp

                else
                    posSwitches(se:i - 1) = posSwitches(se:i - 1) - 1.0_dp

                end if

                step = csf_i%stepvector(i)

                call getMixedFullStop(step, step, 0, csf_i%B_real(i), x1_element=end_mat)

                if (.not. near_zero(end_mat)) then

                    ! only recalc. pgen if matrix element is not 0
                    excitInfo%fullEnd = i
                    excitInfo%firstEnd = i

                    weights = init_semiStartWeight(csf_i, se, i, negSwitches(se), &
                                                   posSwitches(se), csf_i%B_real(se))

                    new_pgen = 1.0_dp

                    ! deal with the start and semi-start seperately
                    if (csf_i%Occ_int(st) /= 1) then
                        new_pgen = new_pgen * weight_funcs(st)%ptr(weights, &
                               csf_i%B_real(st), negSwitches(st), posSwitches(st))
                    end if

                    do j = st + 1, se - 1
                        if (csf_i%Occ_int(j) /= 1) cycle

                        new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                              csf_i%B_real(j), negSwitches(j), posSwitches(j))
                    end do

                    ! then need to reinit double weight
                    weights = weights%ptr

                    ! and also with the semi-start
                    if (csf_i%Occ_int(se) /= 1) then
                        new_pgen = new_pgen * weight_funcs(se)%ptr(weights, &
                           csf_i%B_real(se), negSwitches(se), posSwitches(se))
                    end if

                    do j = se + 1, i - 1
                        if (csf_i%Occ_int(j) /= 1) cycle

                        new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                              csf_i%B_real(j), negSwitches(j), posSwitches(j))
                    end do

                    if (t_trunc_guga_pgen .or. &
                        (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                        if (new_pgen < trunc_guga_pgen) then
                            new_pgen = 0.0_dp
                        end if
                    end if

                    pgen = pgen + new_pgen * orb_pgen

                end if
            end do

            ! deal with switch specifically:

            ! figure out orbital pgen
            call calc_orbital_pgen_contrib_end(csf_i, [2 * elecInd, 2 * sw], &
                holeInd, orb_pgen)

            if (.not. near_zero(orb_pgen)) then

                step = csf_i%stepvector(sw)

                if (step == 1) then
                    ! then a -2 branch arrived!
                    call getMixedFullStop(2, 1, -2, csf_i%B_real(sw), x1_element=end_mat)

                    ! also reduce negative switches then
                    ! only everything to the left or?
                    negSwitches(se:sw - 1) = negSwitches(se:sw - 1) - 1.0_dp

                else
                    ! +2 branch arrived!
                    call getMixedFullStop(1, 2, 2, csf_i%B_real(sw), x1_element=end_mat)

                    ! reduce positive switchtes otherwise
                    posSwitches(se:sw - 1) = posSwitches(se:sw - 1) - 1.0_dp

                end if

                ! loop to get correct pgen
                new_pgen = 1.0_dp

                weights = init_semiStartWeight(csf_i, se, sw, negSwitches(se), &
                                               posSwitches(se), csf_i%B_real(se))

                ! deal with the start and semi-start seperately
                if (csf_i%Occ_int(st) /= 1) then
                    new_pgen = new_pgen * weight_funcs(st)%ptr(weights, &
                           csf_i%B_real(st), negSwitches(st), posSwitches(st))
                end if

                do j = st + 1, se - 1
                    if (csf_i%Occ_int(j) /= 1) cycle

                    new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                              csf_i%B_real(j), negSwitches(j), posSwitches(j))
                end do

                weights = weights%ptr

                ! and also with the semi-start
                if (csf_i%Occ_int(se) /= 1) then
                    new_pgen = new_pgen * weight_funcs(se)%ptr(weights, &
                           csf_i%B_real(se), negSwitches(se), posSwitches(se))
                end if

                do j = se + 1, sw - 1
                    if (csf_i%Occ_int(j) /= 1) cycle

                    new_pgen = new_pgen * weight_funcs(j)%ptr(weights, &
                              csf_i%B_real(j), negSwitches(j), posSwitches(j))
                end do

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                    if (new_pgen < trunc_guga_pgen) then
                        new_pgen = 0.0_dp
                    end if
                end if

                pgen = pgen + new_pgen * orb_pgen

            end if
        end if

        pgen = pgen / real(ElecPairs, dp)

        if (csf_i%stepvector(elecInd) == 3) pgen = pgen * 2.0_dp

    end subroutine calc_mixed_end_contr_pgen

    pure subroutine calc_mixed_end_contr_integral(ilut, csf_i, t, excitInfo, integral, &
            rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        HElement_t(dp), intent(out) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_mixed_end_contr_integral"

        integer :: st, se, en, step, sw, elecInd, holeInd, i
        real(dp) :: top_cont, mat_ele, stay_mat, end_mat
        logical :: above_flag

        integer(int_rdm), allocatable :: tmp_rdm_ind(:)
        real(dp), allocatable :: tmp_rdm_mat(:)
        logical :: rdm_flag
        integer :: rdm_count, max_num_rdm

        if (present(rdm_ind) .or. present(rdm_mat)) then
            ASSERT(present(rdm_ind))
            ASSERT(present(rdm_mat))
            rdm_flag = .true.
        else
            rdm_flag = .false.
        end if

        ! do as much stuff as possible beforehand
        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd
        if (excitInfo%typ == excit_type%fullstop_L_to_R) then
            elecInd = st
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstop_R_to_L) then
            elecInd = se
            holeInd = st
        else
            call stop_all(this_routine, "should not be here!")
        end if

        integral = h_cast(0.0_dp)

        step = csf_i%stepvector(en)

        sw = findLastSwitch(ilut, t, se, en)

        if (rdm_flag) then
            max_num_rdm = (nSpatOrbs - sw + 1)
            allocate(tmp_rdm_ind(max_num_rdm), source=0_int_rdm)
            allocate(tmp_rdm_mat(max_num_rdm), source=0.0_dp)
            rdm_count = 0
        end if

        if (en < nSpatOrbs) then
            select case (step)
            case (1)
                if (isOne(t, en)) then
                    top_cont = -Root2 * sqrt((csf_i%B_real(en) + 2.0_dp) / &
                                             csf_i%B_real(en))

                else
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                end if
            case (2)
                if (isOne(t, en)) then
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                else
                    top_cont = Root2 * sqrt(csf_i%B_real(en) / &
                                            (csf_i%B_real(en) + 2.0_dp))
                end if

            case default
                call stop_all(this_routine, "wrong stepvalues!")

            end select

            if (.not. near_zero(top_cont)) then

                above_flag = .false.
                mat_ele = 1.0_dp

                do i = en + 1, nSpatOrbs
                    if (csf_i%Occ_int(i) /= 1) cycle

                    ! then check if thats the last step
                    if (csf_i%stepvector(i) == 2 .and. csf_i%B_int(i) == 0) then
                        above_flag = .true.
                    end if

                    ! should be able to do that without second loop too!
                    ! figure out!
                    step = csf_i%stepvector(i)

                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, &
                        csf_i%B_real(i), 1.0_dp, x1_element=stay_mat)

                    call getMixedFullStop(step, step, 0, csf_i%B_real(i), &
                                          x1_element=end_mat)

                    ! this check should never be true, but just to be sure
                    if (near_zero(stay_mat)) above_flag = .true.

                    if (.not. near_zero(end_mat)) then
                        integral = integral + end_mat * mat_ele * &
                                   (get_umat_el(i, holeInd, elecInd, i) + &
                                    get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp

                        if (rdm_flag) then
                            rdm_count = rdm_count + 1
                            tmp_rdm_ind(rdm_count) = &
                                contract_2_rdm_ind(i, elecInd, holeInd, i)
                            tmp_rdm_mat(rdm_count) = top_cont * end_mat * mat_ele
                        end if

                    end if

                    if (above_flag) exit

                    ! otherwise update your running pgen and matrix element vars
                    mat_ele = mat_ele * stay_mat

                end do

                integral = integral * top_cont
            end if
        end if

        if (rdm_flag .and. sw == en) then
            rdm_count = rdm_count + 1
            tmp_rdm_ind(rdm_count) = &
                contract_2_rdm_ind(sw, elecInd, holeInd, sw)
            tmp_rdm_mat(rdm_count) = 1.0_dp
        end if

        if (sw < en) then

            step = csf_i%stepvector(en)

            ! inverse fullstop matrix element
            call getMixedFullStop(step, step, 0, csf_i%B_real(en), x1_element=mat_ele)

            ASSERT(.not. near_zero(mat_ele))

            mat_ele = 1.0_dp / mat_ele

            do i = en - 1, sw + 1, -1

                if (csf_i%Occ_int(i) /= 1) cycle

                step = csf_i%stepvector(i)
                ! update inverse product
                call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, csf_i%B_real(i), &
                                            1.0_dp, x1_element=stay_mat)

                call getMixedFullStop(step, step, 0, csf_i%B_real(i), x1_element=end_mat)

                ! update matrix element
                ASSERT(.not. near_zero(stay_mat))
                mat_ele = mat_ele / stay_mat

                if (.not. near_zero(end_mat)) then

                    integral = integral + end_mat * mat_ele * &
                               (get_umat_el(i, holeInd, elecInd, i) + &
                                get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp

                    if (rdm_flag) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = &
                            contract_2_rdm_ind(i, elecInd, holeInd, i)
                        tmp_rdm_mat(rdm_count) = end_mat * mat_ele
                    end if
                end if
            end do

            step = csf_i%stepvector(sw)

            if (step == 1) then
                ! then a -2 branch arrived!
                call getDoubleMatrixElement(2, 1, -2, gen_type%L, gen_type%R, csf_i%B_real(sw), &
                                            1.0_dp, x1_element=stay_mat)

                call getMixedFullStop(2, 1, -2, csf_i%B_real(sw), x1_element=end_mat)

            else
                ! +2 branch arrived!

                call getDoubleMatrixElement(1, 2, 2, gen_type%L, gen_type%R, csf_i%B_real(sw), &
                                            1.0_dp, x1_element=stay_mat)

                call getMixedFullStop(1, 2, 2, csf_i%B_real(sw), x1_element=end_mat)

            end if

            ASSERT(.not. near_zero(stay_mat))

            mat_ele = mat_ele * end_mat / stay_mat

            integral = integral + mat_ele * (get_umat_el(sw, holeInd, elecInd, sw) + &
                                             get_umat_el(holeInd, sw, sw, elecInd)) / 2.0_dp

            if (rdm_flag) then
                rdm_count = rdm_count + 1
                tmp_rdm_ind(rdm_count) = &
                    contract_2_rdm_ind(sw, elecInd, holeInd, sw)
                tmp_rdm_mat(rdm_count) = mat_ele
            end if
        end if

        if (rdm_flag) then
            allocate(rdm_ind(rdm_count), source=tmp_rdm_ind(1:rdm_count))
            allocate(rdm_mat(rdm_count), source=tmp_rdm_mat(1:rdm_count))

            deallocate(tmp_rdm_ind)
            deallocate(tmp_rdm_mat)
        end if

    end subroutine calc_mixed_end_contr_integral

    subroutine calc_mixed_contr_sym(ilut, csf_i, t, excitInfo, pgen, integral)
        ! new implementation of the pgen contribution calculation for
        ! fullstart into fullstop excitation with mixed generators
        ! this is a specific implementation for the hubbard/ueg model with
        ! full k-point symmetry information, since in this case the condition
        ! ki + kj = ka + kb is always fullfilled since ka = ki, kj = kb or vv.
        ! NEW: combine pgen and matrix element contribution finally!
        ! to optimize!
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: integral

        integer :: first, last, deltaB(nSpatOrbs), i, j, k, step1, step2
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), &
                    zeroWeight, minusWeight, plusWeight, branch_weight, tempWeight, &
                    inter, tempWeight_1, &
                    above_cpt, below_cpt
        type(WeightObj_t) :: weights
        logical :: above_flag, below_flag, test_skip
        HElement_t(dp) :: temp_int
        type(ExcitationInformation_t) :: tmp_excitInfo

        test_skip = .false.

        tmp_excitInfo = excitInfo

        ! also need the pgen contributions from all other index combinations
        ! shich could lead to this excitation
        first = findFirstSwitch(ilut, t, excitInfo%fullStart, excitInfo%fullEnd)
        last = findLastSwitch(ilut, t, first, excitInfo%fullEnd)

        below_flag = .false.
        above_flag = .false.
        pgen = 0.0_dp

        deltaB = int(csf_i%B_real - calcB_vector_ilut(t(0:nifd)))

        inter = 1.0_dp
        integral = h_cast(0.0_dp)

        ! calculate the always involved intermediate matrix element from
        ! first switch to last switch
        do i = first + 1, last - 1
            if (csf_i%Occ_int(i) /= 1) cycle

            step1 = csf_i%stepvector(i)
            step2 = getStepvalue(t, i)
            call getDoubleMatrixElement(step2, step1, deltaB(i - 1), gen_type%L, gen_type%R, &
                                        csf_i%B_real(i), 1.0_dp, x1_element=tempWeight)

            inter = inter * tempWeight
        end do

        ! also have to recalc. the ab-orbital cumulative probability distrib
        ! essentially this should be the only difference to molucular
        ! calculations.. where i additionally have to check if the corresponding
        ! orbitals are symmetry allowed.. todo
        ! cant do that so generally out here.. since this list depends on the
        ! picked electrons too! -> which i have to recalc.

        ! use a global list of open orbitals to reduce computational amount
        ! and also check for (d=1,b=1) / (d=2,b=0) spot below/above there
        ! is no additional contribution, due to 0 matrix element
        ! hm... or is it too much effort to recalc. a list of open orbitals
        ! generally .. since its actually only used here,
        ! or think about storing a persistent list of open orbitals, and
        ! the number of open orbitals for every CSF updated and calculated
        ! during excitation generation...

        ! still think about that and then make a new efficient implementation
        ! todo!
        ! better to first loop over j, since only for each new end, the weights
        ! have to be recalced..

        do j = last, nSpatOrbs
            if (csf_i%Occ_int(j) /= 1) cycle

            ! calculate the remaining switches once for each (j) but do it
            ! for the worst case until i = 1

            ! check if this is the last end needed to consider
            if (csf_i%stepvector(j) == 2 .and. csf_i%B_int(j) == 0) then
                above_flag = .true.
            end if

            excitInfo%fullStart = 1
            excitInfo%secondStart = 1
            excitInfo%fullEnd = j
            excitInfo%firstEnd = j
            ! reinit remainings switches and weights
            ! i think i could also do a on-the fly switch recalculation..
            ! so only the weights have to be reinited
            call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, posSwitches, &
                                                        negSwitches)

            weights = init_doubleWeight(csf_i, j)

            ! i have to reset the below flag each iteration of j..
            below_flag = .false.

            do i = first, 1, -1
                if (csf_i%Occ_int(i) /= 1) cycle

                if (below_flag) exit

                ! this is the only difference for molecular/hubbard/ueg
                ! calculations
                call calc_orbital_pgen_contr(csf_i, [2 * i, 2 * j], above_cpt, &
                                             below_cpt)

                ! yes they can, and then this orbital does not contribute to the
                ! obtained excitaiton -> cycle..
                ! only from the names.. shouldnt here be below_cpt?
                ! ah ok below is the below_flag!

                ! if the bottom stepvector d = 1 and b = 1 there is no
                ! additional contribution from below, since the x1 matrix
                ! element is 0
                ! same if d = 2 and b = 0 for fullstop stepvector
                if (csf_i%stepvector(i) == 1 .and. csf_i%B_int(i) == 1) then
                    below_flag = .true.
                end if

                temp_int = (get_umat_el(i, j, j, i) + get_umat_el(j, i, i, j)) / 2.0_dp

                if (.not. test_skip) then
                    if (near_zero(above_cpt)) then
                        if (.not. near_zero(temp_int)) then
                            print *, "mixed: above is 0, integral: ", temp_int
                        end if
                        cycle
                    end if
                    if (near_zero(below_cpt)) then
                        if (.not. near_zero(temp_int)) then
                            print *, "mixed: below is 0, integral: ", temp_int
                        end if
                        cycle
                    end if

                end if

                ! calculate the branch probability


                if (t_heisenberg_model .or. t_tJ_model) then
                    ! in the heisenberg and t-J i never pick combinations,
                    ! with 0 matrix element..
                    if (near_zero(temp_int)) cycle
                end if

                zeroWeight = weights%proc%zero(negSwitches(i), posSwitches(i), &
                    csf_i%B_real(i), weights%dat)

                ! deal with the start seperately:
                if (csf_i%stepvector(i) == 1) then
                    plusWeight = weights%proc%plus(posSwitches(i), &
                                                   csf_i%B_real(i), weights%dat)
                    if (isOne(t, i)) then
                        branch_weight = zeroWeight / (zeroWeight + plusWeight)
                    else
                        branch_weight = plusWeight / (zeroWeight + plusWeight)
                    end if
                else
                    minusWeight = weights%proc%minus(negSwitches(i), &
                                                     csf_i%B_real(i), weights%dat)
                    if (isTwo(t, i)) then
                        branch_weight = zeroWeight / (zeroWeight + minusWeight)
                    else
                        branch_weight = minusWeight / (zeroWeight + minusWeight)
                    end if
                end if

                ! get the starting matrix element
                step1 = csf_i%stepvector(i)
                step2 = getStepvalue(t, i)
                call getDoubleMatrixElement(step2, step1, -1, gen_type%L, gen_type%R, &
                                            csf_i%B_real(i), 1.0_dp, x1_element=tempWeight)

                ! loop over excitation range
                ! distinguish between different regimes
                ! if i do it until switch - 1 -> i know that dB = 0 and
                ! the 2 stepvalues are always the same..
                do k = i + 1, first - 1
                    if (csf_i%Occ_int(k) /= 1) cycle

                    step1 = csf_i%stepvector(k)
                    ! only 0 branch here
                    call getDoubleMatrixElement(step1, step1, 0, gen_type%L, gen_type%R, &
                                                csf_i%B_real(k), 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = tempWeight * tempWeight_1

                    zeroWeight = weights%proc%zero(negSwitches(k), &
                                                   posSwitches(k), csf_i%B_real(k), weights%dat)

                    if (step1 == 1) then
                        plusWeight = weights%proc%plus(posSwitches(k), &
                                                       csf_i%B_real(k), weights%dat)

                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, plusWeight, csf_i%B_real(k))

                    else
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                         csf_i%B_real(k), weights%dat)

                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, minusWeight, csf_i%B_real(k))
                    end if

                end do

                ! then do first switch site seperately, if (i) is not first
                ! and what if (i) is first??
                if (i /= first) then
                    step1 = csf_i%stepvector(first)

                    zeroWeight = weights%proc%zero(negSwitches(first), &
                           posSwitches(first), csf_i%B_real(first), weights%dat)

                    if (step1 == 1) then
                        ! i know that step2 = 2
                        call getDoubleMatrixElement(2, 1, 0, gen_type%L, gen_type%R, &
                                csf_i%B_real(first), 1.0_dp, x1_element=tempWeight_1)

                        plusWeight = weights%proc%plus(posSwitches(first), &
                                                       csf_i%B_real(first), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                         zeroWeight, plusWeight, csf_i%B_real(first)))

                    else
                        ! i know that step2 = 1
                        call getDoubleMatrixElement(1, 2, 0, gen_type%L, gen_type%R, &
                            csf_i%B_real(first), 1.0_dp, x1_element=tempWeight_1)

                        minusWeight = weights%proc%minus(negSwitches(first), &
                                             csf_i%B_real(first), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                         zeroWeight, minusWeight, csf_i%B_real(first)))

                    end if
                    tempWeight = tempWeight * tempWeight_1

                end if

                ! loop over the range where switch happened
                do k = first + 1, last - 1
                    ! in this region i know, that the matrix element is
                    ! definetly not 0, since otherwise the excitation would
                    ! have been aborted before
                    ! combine stepvalue and deltaB info in select statement

                    if (csf_i%Occ_int(k) /= 1) cycle

                    zeroWeight = weights%proc%zero(negSwitches(k), &
                                       posSwitches(k), csf_i%B_real(k), weights%dat)

                    select case (deltaB(k - 1) + csf_i%stepvector(k))

                    case (1)
                        ! d=1 + b=0 : 1
                        plusWeight = weights%proc%plus(posSwitches(k), &
                                                       csf_i%B_real(k), weights%dat)
                        if (isOne(t, k)) then
                            branch_weight = branch_weight * calcStayingProb( &
                                            zeroWeight, plusWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             zeroWeight, plusWeight, csf_i%B_real(k)))
                        end if

                    case (2)
                        ! d=2 + b=0 : 2
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                         csf_i%B_real(k), weights%dat)

                        if (isTwo(t, k)) then
                            branch_weight = branch_weight * calcStayingProb( &
                                            zeroWeight, minusWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                                 zeroWeight, minusWeight, csf_i%B_real(k)))
                        end if

                    case (-1)
                        ! d=1 + b=-2 : -1
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                 csf_i%B_real(k), weights%dat)

                        if (isOne(t, k)) then
                            branch_weight = branch_weight * calcStayingProb(minusWeight, &
                                                            zeroWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             minusWeight, zeroWeight, csf_i%B_real(k)))
                        end if

                    case (4)
                        ! d=2 + b=2 : 4
                        zeroWeight = weights%proc%zero(negSwitches(k), &
                                           posSwitches(k), csf_i%B_real(k), weights%dat)

                        plusWeight = weights%proc%plus(posSwitches(k), &
                                                   csf_i%B_real(k), weights%dat)

                        if (isTwo(t, k)) then
                            branch_weight = branch_weight * calcStayingProb(plusWeight, &
                                                            zeroWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                                 plusWeight, zeroWeight, csf_i%B_real(k)))
                        end if

                    end select

                end do

                ! more efficient to do "last" step seperately, since i have to
                ! check deltaB value and also have to consider matrix element
                ! but only of (j) is not last or otherwise already dealt with
                if (j /= last) then

                    if (csf_i%stepvector(last) == 1) then
                        ! then i know step2 = 2 & dB = -2!
                        call getDoubleMatrixElement(2, 1, -2, gen_type%L, gen_type%R, &
                                    csf_i%B_real(last), 1.0_dp, x1_element=tempWeight_1)

                        zeroWeight = weights%proc%zero(negSwitches(last), &
                                       posSwitches(last), csf_i%B_real(last), weights%dat)

                        minusWeight = weights%proc%minus(negSwitches(last), &
                                             csf_i%B_real(last), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             minusWeight, zeroWeight, csf_i%B_real(last)))

                    else
                        ! i know step2 == 1 and dB = +2
                        call getDoubleMatrixElement(1, 2, +2, gen_type%L, gen_type%R, &
                                        csf_i%B_real(last), 1.0_dp, x1_element=tempWeight_1)

                        zeroWeight = weights%proc%zero(negSwitches(last), &
                                           posSwitches(last), csf_i%B_real(last), weights%dat)

                        plusWeight = weights%proc%plus(posSwitches(last), &
                                                   csf_i%B_real(last), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             plusWeight, zeroWeight, csf_i%B_real(last)))

                    end if

                    tempWeight = tempWeight * tempWeight_1
                end if

                ! then do remaining top range, where i know stepvalues are
                ! the same again and dB = 0 always!
                do k = last + 1, j - 1
                    if (csf_i%Occ_int(k) /= 1) cycle

                    step1 = csf_i%stepvector(k)
                    ! only 0 branch here
                    call getDoubleMatrixElement(step1, step1, 0, gen_type%L, gen_type%R, &
                                    csf_i%B_real(k), 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = tempWeight * tempWeight_1

                    zeroWeight = weights%proc%zero(negSwitches(k), &
                               posSwitches(k), csf_i%B_real(k), weights%dat)

                    if (step1 == 1) then
                        ! i know step2 = 1 als
                        plusWeight = weights%proc%plus(posSwitches(k), &
                                               csf_i%B_real(k), weights%dat)

                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, plusWeight, csf_i%B_real(k))
                    else
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                         csf_i%B_real(k), weights%dat)
                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, minusWeight, csf_i%B_real(k))
                    end if
                end do

                ! and handle fullend
                ! and then do the the end value at j
                step1 = csf_i%stepvector(j)
                step2 = getStepvalue(t, j)
                call getMixedFullStop(step2, step1, deltaB(j - 1), csf_i%B_real(j), &
                                      x1_element=tempWeight_1)

                temp_int = tempWeight * tempWeight_1 * inter * temp_int

                ! and multiply and add up all contribution elements
                integral = integral + temp_int

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                    if (branch_weight < trunc_guga_pgen) then
                        branch_weight = 0.0_dp
                    end if
                end if

                if (t_trunc_guga_matel) then
                    if (abs(tempWeight * tempWeight_1 * inter) < trunc_guga_matel) then
                        branch_weight = 0.0_dp
                    end if
                end if

                ! add up pgen contributions..
                pgen = pgen + (below_cpt + above_cpt) * branch_weight

                ! check if i deal with that correctly...
                if (below_flag) exit
            end do
            ! todo: i cant use tthat like that.. or else some combinations
            ! of i and j get left out! i have to reinit it somehow..
            ! not yet sure how..
            if (above_flag) exit
        end do

        ! multiply by always same probability to pick the 2 electrons
        if (.not. (t_heisenberg_model .or. t_tJ_model)) then
            pgen = pgen / real(ElecPairs, dp)
        end if

    end subroutine calc_mixed_contr_sym

    pure subroutine calc_mixed_contr_integral(ilut, csf_i, t, start, ende, integral, &
            rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: start, ende
        HElement_t(dp), intent(out) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_mixed_contr_integral"

        real(dp) :: inter, tempWeight, tempWeight_1
        HElement_t(dp) :: temp_int
        integer :: i, j, k, step1, step2, bVector(nSpatOrbs), first, last
        logical :: rdm_flag
        integer :: max_num_rdm, rdm_count
        integer(int_rdm), allocatable :: tmp_rdm_ind(:)
        real(dp), allocatable :: tmp_rdm_mat(:)
        real(dp) :: tmp_mat

        if (present(rdm_ind) .or. present(rdm_mat)) then
            ASSERT(present(rdm_ind))
            ASSERT(present(rdm_mat))
            rdm_flag = .true.
        else
            rdm_flag = .false.
        end if
        ! do it differently... since its always a mixed double overlap region
        ! and only 2 indices are involved.. just recalc all possible
        ! matrix elements in this case..
        ! so first determine the first and last switches and calculate
        ! the overlap matrix element in this region, since atleast thats
        ! always the same
        first = findFirstSwitch(ilut, t, start, ende)
        last = findLastSwitch(ilut, t, first, ende)

        if (rdm_flag) then
            max_num_rdm = first * (nSpatOrbs - last + 1)
            allocate(tmp_rdm_ind(max_num_rdm), source=0_int_rdm)
            allocate(tmp_rdm_mat(max_num_rdm), source=0.0_dp)
            rdm_count = 0
        end if
        ! calc. the intermediate matrix element..
        ! but what is the deltaB value inbetween? calculate it on the fly..
        bVector = int(calcB_vector_ilut(ilut(0:nifd)) - calcB_vector_ilut(t(0:nifd)))

        inter = 1.0_dp
        integral = h_cast(0.0_dp)

        do i = first + 1, last - 1
            if (csf_i%Occ_int(i) /= 1) cycle

            step1 = csf_i%stepvector(i)
            step2 = getStepvalue(t, i)
            call getDoubleMatrixElement(step2, step1, bVector(i - 1), gen_type%L, gen_type%R, &
                                        csf_i%B_real(i), 1.0_dp, x1_element=tempWeight)

            inter = inter * tempWeight
        end do

        ! and then add up all the possible integral contribs

        do i = 1, first
            if (csf_i%Occ_int(i) /= 1) cycle

            do j = last, nSpatOrbs
                if (csf_i%Occ_int(j) /= 1) cycle

                temp_int = (get_umat_el(i, j, j, i) + get_umat_el(j, i, i, j)) / 2.0_dp


                if ((.not. rdm_flag) .and. near_zero(temp_int)) cycle

                ! get the starting matrix element
                step1 = csf_i%stepvector(i)
                step2 = getStepvalue(t, i)
                call getDoubleMatrixElement(step2, step1, -1, gen_type%L, gen_type%R, &
                                            csf_i%B_real(i), 1.0_dp, x1_element=tempWeight)

                ! then calc. the product:
                do k = i + 1, first
                    if (csf_i%Occ_int(k) /= 1) cycle

                    step1 = csf_i%stepvector(k)
                    step2 = getStepvalue(t, k)
                    ! only 0 branch here
                    call getDoubleMatrixElement(step2, step1, 0, gen_type%L, gen_type%R, &
                                                csf_i%B_real(k), 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = tempWeight * tempWeight_1

                end do

                do k = last, j - 1
                    if (csf_i%Occ_int(k) /= 1) cycle

                    step1 = csf_i%stepvector(k)
                    step2 = getStepvalue(t, k)
                    call getDoubleMatrixElement(step2, step1, bVector(k - 1), gen_type%L, gen_type%R, &
                                                csf_i%B_real(k), 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = tempWeight * tempWeight_1
                end do

                ! and then do the the end value at j
                step1 = csf_i%stepvector(j)
                step2 = getStepvalue(t, j)
                call getMixedFullStop(step2, step1, bVector(j - 1), csf_i%B_real(j), &
                                      x1_element=tempWeight_1)

                ! and multiply and add up all contribution elements
                integral = integral + tempWeight * tempWeight_1 * inter * temp_int

                if (rdm_flag) then
                    tmp_mat = tempWeight * tempWeight_1 * inter
                    if (.not. near_zero(tmp_mat)) then
                        rdm_count = rdm_count + 1
                        tmp_rdm_ind(rdm_count) = contract_2_rdm_ind(i, j, j, i)
                        tmp_rdm_mat(rdm_count) = tmp_mat
                    end if
                end if
            end do
        end do

        if (rdm_flag) then
            allocate(rdm_ind(rdm_count), source = tmp_rdm_ind(1:rdm_count))
            allocate(rdm_mat(rdm_count), source = tmp_rdm_mat(1:rdm_count))
        end if
    end subroutine calc_mixed_contr_integral

    subroutine calc_mixed_contr_pgen(ilut, csf_i, t, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), value :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: first, last, deltaB(nSpatOrbs), i, j, k, step1
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), &
                    zeroWeight, minusWeight, plusWeight, branch_weight, &
                    above_cpt, below_cpt
        type(WeightObj_t) :: weights
        logical :: above_flag, below_flag

        ! also need the pgen contributions from all other index combinations
        ! shich could lead to this excitation
        first = findFirstSwitch(ilut, t, excitInfo%fullStart, excitInfo%fullEnd)
        last = findLastSwitch(ilut, t, first, excitInfo%fullEnd)

        below_flag = .false.
        above_flag = .false.
        pgen = 0.0_dp

        deltaB = int(csf_i%B_real - calcB_vector_ilut(t(0:nifd)))

        do j = last, nSpatOrbs
            if (csf_i%Occ_int(j) /= 1) cycle

            ! calculate the remaining switches once for each (j) but do it
            ! for the worst case until i = 1

            ! check if this is the last end needed to consider
            if (csf_i%stepvector(j) == 2 .and. csf_i%B_int(j) == 0) then
                above_flag = .true.
            end if

            excitInfo%fullStart = 1
            excitInfo%secondStart = 1
            excitInfo%fullEnd = j
            excitInfo%firstEnd = j
            ! reinit remainings switches and weights
            ! i think i could also do a on-the fly switch recalculation..
            ! so only the weights have to be reinited
            call calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, &
                posSwitches, negSwitches)

            weights = init_doubleWeight(csf_i, j)

            ! i have to reset the below flag each iteration of j..
            below_flag = .false.

            do i = first, 1, -1
                if (csf_i%Occ_int(i) /= 1) cycle

                if (below_flag) exit

                ! this is the only difference for molecular/hubbard/ueg
                ! calculations
                call calc_orbital_pgen_contr(csf_i, [2 * i, 2 * j], above_cpt, &
                                             below_cpt)

                ! yes they can, and then this orbital does not contribute to the
                ! obtained excitaiton -> cycle..
                ! only from the names.. shouldnt here be below_cpt?
                ! ah ok below is the below_flag!

                ! if the bottom stepvector d = 1 and b = 1 there is no
                ! additional contribution from below, since the x1 matrix
                ! element is 0
                ! same if d = 2 and b = 0 for fullstop stepvector
                if (csf_i%stepvector(i) == 1 .and. csf_i%B_int(i) == 1) then
                    below_flag = .true.
                end if

                zeroWeight = weights%proc%zero(negSwitches(i), posSwitches(i), &
                    csf_i%B_real(i), weights%dat)

                ! deal with the start seperately:
                if (csf_i%stepvector(i) == 1) then
                    plusWeight = weights%proc%plus(posSwitches(i), &
                                                   csf_i%B_real(i), weights%dat)
                    if (isOne(t, i)) then
                        branch_weight = zeroWeight / (zeroWeight + plusWeight)
                    else
                        branch_weight = plusWeight / (zeroWeight + plusWeight)
                    end if
                else
                    minusWeight = weights%proc%minus(negSwitches(i), &
                                                     csf_i%B_real(i), weights%dat)
                    if (isTwo(t, i)) then
                        branch_weight = zeroWeight / (zeroWeight + minusWeight)
                    else
                        branch_weight = minusWeight / (zeroWeight + minusWeight)
                    end if
                end if

                ! loop over excitation range
                ! distinguish between different regimes
                ! if i do it until switch - 1 -> i know that dB = 0 and
                ! the 2 stepvalues are always the same..
                do k = i + 1, first - 1
                    if (csf_i%Occ_int(k) /= 1) cycle

                    step1 = csf_i%stepvector(k)

                    zeroWeight = weights%proc%zero(negSwitches(k), &
                                                   posSwitches(k), csf_i%B_real(k), weights%dat)

                    if (step1 == 1) then
                        plusWeight = weights%proc%plus(posSwitches(k), &
                                                       csf_i%B_real(k), weights%dat)

                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, plusWeight, csf_i%B_real(k))

                    else
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                         csf_i%B_real(k), weights%dat)

                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, minusWeight, csf_i%B_real(k))
                    end if

                end do

                ! then do first switch site seperately, if (i) is not first
                ! and what if (i) is first??
                if (i /= first) then
                    step1 = csf_i%stepvector(first)

                    zeroWeight = weights%proc%zero(negSwitches(first), &
                           posSwitches(first), csf_i%B_real(first), weights%dat)

                    if (step1 == 1) then
                        ! i know that step2 = 2
                        plusWeight = weights%proc%plus(posSwitches(first), &
                                                       csf_i%B_real(first), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                         zeroWeight, plusWeight, csf_i%B_real(first)))

                    else
                        ! i know that step2 = 1
                        minusWeight = weights%proc%minus(negSwitches(first), &
                                             csf_i%B_real(first), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                         zeroWeight, minusWeight, csf_i%B_real(first)))

                    end if

                end if

                ! loop over the range where switch happened
                do k = first + 1, last - 1
                    ! in this region i know, that the matrix element is
                    ! definetly not 0, since otherwise the excitation would
                    ! have been aborted before
                    ! combine stepvalue and deltaB info in select statement

                    if (csf_i%Occ_int(k) /= 1) cycle

                    zeroWeight = weights%proc%zero(negSwitches(k), &
                                       posSwitches(k), csf_i%B_real(k), weights%dat)

                    select case (deltaB(k - 1) + csf_i%stepvector(k))

                    case (1)
                        ! d=1 + b=0 : 1
                        plusWeight = weights%proc%plus(posSwitches(k), &
                                                       csf_i%B_real(k), weights%dat)
                        if (isOne(t, k)) then
                            branch_weight = branch_weight * calcStayingProb( &
                                            zeroWeight, plusWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             zeroWeight, plusWeight, csf_i%B_real(k)))
                        end if

                    case (2)
                        ! d=2 + b=0 : 2
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                         csf_i%B_real(k), weights%dat)

                        if (isTwo(t, k)) then
                            branch_weight = branch_weight * calcStayingProb( &
                                            zeroWeight, minusWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                                 zeroWeight, minusWeight, csf_i%B_real(k)))
                        end if

                    case (-1)
                        ! d=1 + b=-2 : -1
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                 csf_i%B_real(k), weights%dat)

                        if (isOne(t, k)) then
                            branch_weight = branch_weight * calcStayingProb(minusWeight, &
                                                            zeroWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             minusWeight, zeroWeight, csf_i%B_real(k)))
                        end if

                    case (4)
                        ! d=2 + b=2 : 4
                        zeroWeight = weights%proc%zero(negSwitches(k), &
                                           posSwitches(k), csf_i%B_real(k), weights%dat)

                        plusWeight = weights%proc%plus(posSwitches(k), &
                                                   csf_i%B_real(k), weights%dat)

                        if (isTwo(t, k)) then
                            branch_weight = branch_weight * calcStayingProb(plusWeight, &
                                                            zeroWeight, csf_i%B_real(k))
                        else
                            branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                                 plusWeight, zeroWeight, csf_i%B_real(k)))
                        end if

                    end select

                end do

                ! more efficient to do "last" step seperately, since i have to
                ! check deltaB value and also have to consider matrix element
                ! but only of (j) is not last or otherwise already dealt with
                if (j /= last) then

                    if (csf_i%stepvector(last) == 1) then
                        ! then i know step2 = 2 & dB = -2!
                        zeroWeight = weights%proc%zero(negSwitches(last), &
                                       posSwitches(last), csf_i%B_real(last), weights%dat)

                        minusWeight = weights%proc%minus(negSwitches(last), &
                                             csf_i%B_real(last), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             minusWeight, zeroWeight, csf_i%B_real(last)))

                    else
                        ! i know step2 == 1 and dB = +2
                        zeroWeight = weights%proc%zero(negSwitches(last), &
                                           posSwitches(last), csf_i%B_real(last), weights%dat)

                        plusWeight = weights%proc%plus(posSwitches(last), &
                                                   csf_i%B_real(last), weights%dat)

                        branch_weight = branch_weight * (1.0_dp - calcStayingProb( &
                                             plusWeight, zeroWeight, csf_i%B_real(last)))

                    end if
                end if

                ! then do remaining top range, where i know stepvalues are
                ! the same again and dB = 0 always!
                do k = last + 1, j - 1
                    if (csf_i%Occ_int(k) /= 1) cycle

                    step1 = csf_i%stepvector(k)
                    zeroWeight = weights%proc%zero(negSwitches(k), &
                               posSwitches(k), csf_i%B_real(k), weights%dat)

                    if (step1 == 1) then
                        ! i know step2 = 1 als
                        plusWeight = weights%proc%plus(posSwitches(k), &
                                               csf_i%B_real(k), weights%dat)

                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, plusWeight, csf_i%B_real(k))
                    else
                        minusWeight = weights%proc%minus(negSwitches(k), &
                                                         csf_i%B_real(k), weights%dat)
                        branch_weight = branch_weight * calcStayingProb( &
                                        zeroWeight, minusWeight, csf_i%B_real(k))
                    end if
                end do

                if (t_trunc_guga_pgen .or. &
                    (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                    if (branch_weight < trunc_guga_pgen) then
                        branch_weight = 0.0_dp
                    end if
                end if

                ! add up pgen contributions..
                pgen = pgen + (below_cpt + above_cpt) * branch_weight

                ! check if i deal with that correctly...
                if (below_flag) exit
            end do
            ! todo: i cant use tthat like that.. or else some combinations
            ! of i and j get left out! i have to reinit it somehow..
            ! not yet sure how..
            if (above_flag) exit
        end do

        ! multiply by always same probability to pick the 2 electrons
        if (.not. (t_heisenberg_model .or. t_tJ_model)) then
            pgen = pgen / real(ElecPairs, dp)
        end if

    end subroutine calc_mixed_contr_pgen

    function calcStartProb(prob1, prob2) result(ret)
        ! calculate the probability of a starting branch, given two possibilities
        real(dp), intent(in) :: prob1, prob2
        real(dp) :: ret
        character(*), parameter :: this_routine = "calcStartProb"

        ASSERT(prob1 >= 0.0_dp)
        ASSERT(prob2 >= 0.0_dp)

        ret = prob1 / (prob1 + prob2)

    end function calcStartProb

    function calcStayingProb(prob1, prob2, bVal) result(ret)
        ! calculate the probability to stay on a certain excitation branch
        real(dp), intent(in) :: prob1, prob2, bVal
        real(dp) :: ret, tmp
        character(*), parameter :: this_routine = "calcStayingProb"

        ASSERT(prob1 >= 0.0_dp)
        ASSERT(prob2 >= 0.0_dp)
        ASSERT(bVal >= 0.0_dp)

        ! if b == 0 its stupid to use it..
        tmp = max(bVal, 1.0_dp)

        ret = tmp * prob1 / (tmp * prob1 + prob2)

    end function calcStayingProb

    function init_fullStartWeight(csf_i, sOrb, pOrb, negSwitches, &
                                  posSwitches, bVal) result(fullStart)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb
        real(dp), intent(in) :: negSwitches, posSwitches
        real(dp), intent(in) :: bVal
        type(WeightObj_t) :: fullStart
        character(*), parameter :: this_routine = "init_fullStartWeight"

        type(WeightObj_t), target, save :: single
        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)
        ASSERT(pOrb > 0 .and. pOrb <= nSpatOrbs)
        ASSERT(negSwitches >= 0.0_dp)
        ASSERT(posSwitches >= 0.0_dp)

        fullStart%dat%F = endFx(csf_i, sOrb)
        fullStart%dat%G = endGx(csf_i, sOrb)

        ! have to set up a single weight obj.
        single = init_singleWeight(csf_i, pOrb)

        ! try to reuse the already initialized singles weight in the cause of
        ! an excitation, i hope this works with the pointers and stuff.
        fullstart%ptr => single

        fullStart%dat%minus = single%proc%minus(negSwitches, bVal, single%dat)
        fullStart%dat%plus = single%proc%plus(posSwitches, bVal, single%dat)

        fullStart%proc%minus => getMinus_fullStart
        fullStart%proc%plus => getPlus_fullStart
        fullStart%proc%zero => getZero_fullStart

        fullStart%initialized = .true.

    end function init_fullStartWeight

    elemental function endFx(csf_i, sOrb) result(fx)
        ! flag function used in excitation tree generation to check if spatial
        ! orbital sOrb
        ! is 0,1 or 3. Probably possible to implement it on an efficient
        ! bit-rep level, todo
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb
        real(dp) :: fx

        ! always one except d=2 at end
        if (csf_i%stepvector(sOrb) == 2) then
            fx = 0.0_dp
        else
            fx = 1.0_dp
        end if
    end function endFx

    elemental function endGx(csf_i, sOrb) result(gx)
        ! flag function used in excitation tree generation to check if spatial
        ! orbital sOrb is 0,2,3.
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb
        real(dp) :: gx


        if (csf_i%stepvector(sOrb) == 1) then
            gx = 0.0_dp
        else
            ! always one except d=1 at end
            gx = 1.0_dp
        end if
    end function endGx

    ! proabbilistic weight objects:
    elemental function init_singleWeight(csf_i, sOrb) result(singleWeight)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb
        type(WeightObj_t) :: singleWeight
        character(*), parameter :: this_routine = "init_singleWeight"

        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)

        singleWeight%dat%F = endFx(csf_i, sOrb)
        singleWeight%dat%G = endGx(csf_i, sOrb)

        singleWeight%proc%minus => getMinus_single
        singleWeight%proc%plus => getPlus_single

        singleWeight%initialized = .true.

    end function init_singleWeight

    function getMinus_fullStart(nSwitches, bVal, fullStart) result(minusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: fullStart
        real(dp) :: minusWeight
        character(*), parameter :: this_routine = "getMinus_fullStart"

        ASSERT(nSwitches >= 0.0_dp)
        ! here this assert is valid, since the bvalue really never should be
        ! 0.. or else the parent CSF is invalid.
        ! no of course it can still be 0 at a 0/3 start... but also
        ! set it to the technical value only

        ! try an adhoc second order fix..
        minusWeight = fullStart%F * fullStart%minus + nSwitches / max(1.0_dp, bval) &
                      * (fullStart%G * fullStart%minus + fullStart%F * fullStart%plus) &
                      + (max(nSwitches - 1, 0.0_dp) * fullStart%G * fullStart%plus / (max(1.0_dp, bval)**2))

        ASSERT(minusWeight >= 0.0_dp)

    end function getMinus_fullStart

    function getPlus_fullStart(nSwitches, bVal, fullStart) result(plusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: fullStart
        real(dp) :: plusWeight
        character(*), parameter :: this_routine = "getPlus_fullStart"

        ASSERT(nSwitches >= 0.0_dp)

        ! same as above: have to check if < 2
        if (bVal < 2.0_dp) then
            plusWeight = 0.0_dp
        else
            plusWeight = fullStart%G * fullStart%plus + nSwitches / bVal * &
                         (fullStart%G * fullStart%minus + fullStart%F * fullStart%plus) &
                         + (max(nSwitches - 1, 0.0_dp) * fullStart%F * fullStart%minus / (max(1.0_dp, bval)**2))
        end if

        ASSERT(plusWeight >= 0.0_dp)

    end function getPlus_fullStart

    function getMinus_single(nSwitches, bVal, single) result(minusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: single
        real(dp) :: minusWeight
        character(*), parameter :: this_routine = "getMinus_single"
        ASSERT(nSwitches >= 0.0_dp)
        ! change that, to make it independend if b is zero
        if (near_zero(bVal)) then
            ! make it only depend on f and nSwitches
            ! will be normalized to 1 anyway in the calcStayingProb function
            minusWeight = single%F + nSwitches * single%G

        else
            minusWeight = single%F + nSwitches * single%G / bVal

        end if

        ASSERT(minusWeight >= 0.0_dp)

    end function getMinus_single

    function getPlus_single(nSwitches, bVal, single) result(plusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: single
        real(dp) :: plusWeight
        character(*), parameter :: this_routine = "getPlus_single"
        ASSERT(nSwitches >= 0.0_dp)

        if (near_zero(bVal)) then
            plusWeight = 0.0_dp
        else
            plusWeight = single%G + nSwitches * single%F / bVal
        end if

        ASSERT(plusWeight >= 0.0_dp)

    end function getPlus_single

    function plus_start_single(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, minus

        plus = weights%proc%plus(posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = plus / (plus + minus)

    end function plus_start_single

    function minus_start_single(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, minus

        plus = weights%proc%plus(posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = minus / (plus + minus)

    end function minus_start_single

    function minus_staying_single(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, minus

        plus = weights%proc%plus(posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = calcStayingProb(minus, plus, bVal)

    end function minus_staying_single

    function plus_staying_single(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, minus

        plus = weights%proc%plus(posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = calcStayingProb(plus, minus, bVal)

    end function plus_staying_single

    function plus_switching_single(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, minus

        plus = weights%proc%plus(posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = 1.0_dp - calcStayingProb(plus, minus, bVal)

    end function plus_switching_single

    function minus_switching_single(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, minus

        plus = weights%proc%plus(posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = 1.0_dp - calcStayingProb(minus, plus, bVal)

    end function minus_switching_single

    function minus_start_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: zero, minus

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = minus / (zero + minus)

    end function minus_start_double

    function plus_start_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, zero

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        plus = weights%proc%plus(posSwitches, bVal, weights%dat)

        prob = plus / (zero + plus)

    end function plus_start_double

    function zero_plus_start_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, zero

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        plus = weights%proc%plus(posSwitches, bVal, weights%dat)

        prob = zero / (zero + plus)

    end function zero_plus_start_double

    function zero_minus_start_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: zero, minus

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = zero / (zero + minus)

    end function zero_minus_start_double

    function zero_plus_staying_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, zero

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        plus = weights%proc%plus(posSwitches, bVal, weights%dat)

        prob = calcStayingProb(zero, plus, bVal)

    end function zero_plus_staying_double

    function zero_minus_staying_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: zero, minus

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = calcStayingProb(zero, minus, bVal)

    end function zero_minus_staying_double

    function zero_plus_switching_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, zero

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        plus = weights%proc%plus(posSwitches, bVal, weights%dat)

        prob = 1.0_dp - calcStayingProb(zero, plus, bVal)

    end function zero_plus_switching_double

    function zero_minus_switching_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: zero, minus

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = 1.0_dp - calcStayingProb(zero, minus, bVal)

    end function zero_minus_switching_double

    function minus_staying_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: zero, minus

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = calcStayingProb(minus, zero, bVal)

    end function minus_staying_double

    function plus_staying_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, zero

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        plus = weights%proc%plus(posSwitches, bVal, weights%dat)

        prob = calcStayingProb(plus, zero, bVal)

    end function plus_staying_double

    function minus_switching_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: zero, minus

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        minus = weights%proc%minus(negSwitches, bVal, weights%dat)

        prob = 1.0_dp - calcStayingProb(minus, zero, bVal)

    end function minus_switching_double

    function plus_switching_double(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        real(dp) :: plus, zero

        zero = weights%proc%zero(negSwitches, posSwitches, bVal, weights%dat)
        plus = weights%proc%plus(posSwitches, bVal, weights%dat)

        prob = 1.0_dp - calcStayingProb(plus, zero, bVal)

    end function plus_switching_double

    function probability_one(weights, bVal, negSwitches, posSwitches) result(prob)
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: bVal, negSwitches, posSwitches
        real(dp) :: prob

        unused_var(weights)
        unused_var(bVal)
        unused_var(negSwitches)
        unused_var(posSwitches)

        prob = 1.0_dp

    end function probability_one

    function getZero_fullStart(negSwitches, posSwitches, bVal, fullStart) &
        result(zeroWeight)
        real(dp), intent(in) :: negSwitches, posSwitches, bVal
        type(WeightData_t), intent(in) :: fullStart
        real(dp) :: zeroWeight
        character(*), parameter :: this_routine = "getZero_fullStart"

        ASSERT(negSwitches >= 0.0_dp)
        ASSERT(posSwitches >= 0.0_dp)

        zeroWeight = fullStart%F * fullStart%plus + fullStart%G * fullStart%minus + &
                     (negSwitches * fullStart%G * fullStart%plus + &
                      posSwitches * fullStart%F * fullStart%minus) / max(1.0_dp, bval)

        ASSERT(zeroWeight >= 0.0_dp)

    end function getZero_fullStart

    subroutine calcRemainingSwitches_excitInfo_double(csf_i, excitInfo, &
                                                      posSwitches, negSwitches)
        ! subroutine to determine the number of remaining switches for double
        ! excitations between spatial orbitals (i,j,k,l). orbital indices are
        ! given in type(excitationInformation), extra flag is needed to
        ! indicate that this is a double excitaiton then
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp), intent(out) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer :: iOrb, end1
        real(dp) :: oneCount, twoCount

        ! have to calc. the overlap range of the excitations to more
        ! efficiently decide between different kind of double excitations
        ! even better, get all possible information through excitationIdentifier
        ! assume exitInfo already calculated in calling function
        ! update: already given as input
        !excitInfo = excitationIdentifier(i, j, k, l)

        ! intitialize values
        oneCount = 0.0_dp
        twoCount = 0.0_dp
        posSwitches = 0.0_dp
        negSwitches = 0.0_dp

        if (excitInfo%typ == excit_type%raising .or. &
            excitInfo%typ == excit_type%lowering) then

            call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, &
                                                        posSwitches, negSwitches)
        else

            select case (excitInfo%overlap)
            case (0)
                do iOrb = excitInfo%fullEnd - 1, excitInfo%secondStart, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

                ! reset count past second excitations:
                oneCount = 0.0_dp
                twoCount = 0.0_dp

                do iOrb = excitInfo%firstEnd - 1, excitInfo%fullStart, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

            case (1)
                ! not quite sure anymore why, but have to treat single overlap
                ! excitations with alike generators different then mixed
                ! because it is like a single excitation over the whole excitation
                ! range
                if (excitInfo%gen1 /= excitInfo%gen2) then
                    end1 = 0
                else
                    end1 = excitInfo%firstEnd
                end if

                do iOrb = excitInfo%fullEnd - 1, excitInfo%firstEnd, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

                ! reset the switch number if alike generators are present.
                ! now i am confused.. no its fine.. we actually never want to
                ! go into single overlap with alike generators, since it is
                ! actually a single excitation then!
                if (excitInfo%gen1 == excitInfo%gen2) then
                    oneCount = 0.0_dp
                    twoCount = 0.0_dp
                end if

                do iOrb = excitInfo%firstEnd - 1, excitInfo%fullStart, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

            case default
                ! proper overlap ranges:

                ! do all those excitations in the same way, although this means
                ! for some, that too much work is done.. e.g. for full-start
                ! excitations with alike generators. where only the delta b = 0
                ! branch has non-zero matrix elements in the overlap region
                ! but these things can be handled in the excitations calculation

                ! for certain index combinations some loops wont get executed
                do iOrb = excitInfo%fullEnd - 1, excitInfo%firstEnd, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

                oneCount = 0.0_dp
                twoCount = 0.0_dp

                do iOrb = excitInfo%firstEnd - 1, excitInfo%secondStart, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

                oneCount = 0.0_dp
                twoCount = 0.0_dp

                do iOrb = excitInfo%secondStart - 1, excitInfo%fullStart, -1
                    posSwitches(iOrb) = twoCount
                    negSwitches(iOrb) = oneCount

                    select case (csf_i%stepvector(iOrb))
                    case (1)
                        oneCount = oneCount + 1.0_dp
                    case (2)
                        twoCount = twoCount + 1.0_dp
                    end select
                end do

                oneCount = 0.0_dp
                twoCount = 0.0_dp

            end select

        end if

    end subroutine calcRemainingSwitches_excitInfo_double

    subroutine calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, &
                                                      posSwitches, negSwitches)
        ! subroutine to determine the number of remaining switches for single
        ! excitations between orbitals s, p given in type of excitationInformation.
        ! The switches are given
        ! as a list, to access it for each spatial orbital
        ! stepValue = 1 -> positive delta B switch possibility
        ! stepValue = 2 -> negative delta B switch possibility
        ! assume exitInfo is already calculated
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp), intent(out) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer :: iOrb
        real(dp) :: oneCount, twoCount

        ! ignore b > 0 forced switches for now. As they only change the bias
        ! do not make the excitations calculations wrong if ignored.

        oneCount = 0.0_dp
        twoCount = 0.0_dp
        posSwitches = 0.0_dp
        negSwitches = 0.0_dp
        ! have to count from the reversed ilut entries
        do iOrb = excitInfo%fullEnd - 1, excitInfo%fullStart, -1
            posSwitches(iOrb) = twoCount
            negSwitches(iOrb) = oneCount

            select case (csf_i%stepvector(iOrb))
            case (1)
                oneCount = oneCount + 1.0_dp
            case (2)
                twoCount = twoCount + 1.0_dp
            end select
        end do

    end subroutine calcRemainingSwitches_excitInfo_single

    subroutine setup_weight_funcs(t, csf_i, st, se, weight_funcs)
        integer(n_int), intent(in) :: t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: st, se
        type(BranchWeightArr_t), intent(out) :: weight_funcs(nSpatOrbs)

        integer :: i, step, delta_b(nSpatOrbs)

        delta_b = int(csf_i%B_real - calcB_vector_ilut(t(0:nifd)))

        ! i know that a start was possible -> only check what the excitation
        ! stepvalue is
        ! damn.. where are my notes? im not sure about that..
        if (isOne(t, st)) then
            weight_funcs(st)%ptr => minus_start_single
        else if (isTwo(t, st)) then
            weight_funcs(st)%ptr => plus_start_single
            ! i also need to consider an non-choosing start or deal with
            ! that in the routines above..
        end if

        do i = st + 1, se - 1
            if (csf_i%Occ_int(i) /= 1) cycle

            step = csf_i%stepvector(i)

            if (step == 1 .and. delta_b(i - 1) == -1) then
                if (isOne(t, i)) then
                    weight_funcs(i)%ptr => minus_staying_single
                else
                    weight_funcs(i)%ptr => minus_switching_single
                end if
            else if (step == 2 .and. delta_b(i - 1) == 1) then
                if (isTwo(t, i)) then
                    weight_funcs(i)%ptr => plus_staying_single
                else
                    weight_funcs(i)%ptr => plus_switching_single
                end if
                ! here i need a one-prob. if no switch was possible.. damn..
            else
                weight_funcs(i)%ptr => probability_one
            end if

        end do

        ! similar to the start, only need to check  the stepvalue of the
        ! excitaiton, since we know something must have worked
        if (isOne(t, se)) then
            if (delta_b(se - 1) == -1) then
                weight_funcs(se)%ptr => minus_start_double
            else
                weight_funcs(se)%ptr => zero_plus_start_double
            end if
        else if (isTwo(t, se)) then
            if (delta_b(se - 1) == -1) then
                weight_funcs(se)%ptr => zero_minus_start_double
            else
                weight_funcs(se)%ptr => plus_start_double
            end if
        end if

        do i = se + 1, nSpatOrbs
            if (csf_i%Occ_int(i) /= 1) cycle

            step = csf_i%stepvector(i)

            ! also combine step and deltab value in a select case statement
            select case (delta_b(i - 1) + step)
            case (1)
                ! d=1 + b=0 : 1
                if (isOne(t, i)) then
                    weight_funcs(i)%ptr => zero_plus_staying_double
                else
                    weight_funcs(i)%ptr => zero_plus_switching_double
                end if
            case (2)
                ! d=2 + b=0 :2
                if (isTwo(t, i)) then
                    weight_funcs(i)%ptr => zero_minus_staying_double
                else
                    weight_funcs(i)%ptr => zero_minus_switching_double
                end if

            case (-1)
                if (isOne(t, i)) then
                    weight_funcs(i)%ptr => minus_staying_double
                else
                    weight_funcs(i)%ptr => minus_switching_double
                end if

            case (4)
                if (isTwo(t, i)) then
                    weight_funcs(i)%ptr => plus_staying_double
                else
                    weight_funcs(i)%ptr => plus_switching_double
                end if

                ! i also need a case default to prob 1. for a no-choice..
            case default
                weight_funcs(i)%ptr => probability_one

            end select
        end do

    end subroutine setup_weight_funcs

    function init_semiStartWeight(csf_i, sOrb, pOrb, negSwitches, posSwitches, bVal) &
        result(semiStart)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb
        real(dp), intent(in) :: negSwitches, posSwitches, bVal
        type(WeightObj_t) :: semiStart
        character(*), parameter :: this_routine = "init_semiStartWeight"

        type(WeightObj_t), target, save :: double

        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)
        ASSERT(pOrb > 0 .and. pOrb <= nSpatOrbs)
        ASSERT(negSwitches >= 0.0_dp)
        ASSERT(posSwitches >= 0.0_dp)

        semiStart%dat%F = endFx(csf_i, sOrb)
        semiStart%dat%G = endGx(csf_i, sOrb)

        double = init_doubleWeight(csf_i, pOrb)

        semiStart%ptr => double

        semiStart%dat%minus = double%proc%minus(negSwitches, bVal, double%dat)
        semiStart%dat%plus = double%proc%plus(posSwitches, bVal, double%dat)
        semiStart%dat%zero = double%proc%zero(negSwitches, posSwitches, bVal, double%dat)

        semiStart%proc%minus => getMinus_semiStart
        semiStart%proc%plus => getPlus_semiStart

        semiStart%initialized = .true.

    end function init_semiStartWeight


    function init_doubleWeight(csf_i, sOrb) result(doubleWeight)
        ! obj has the same structure as the semi-start weight, reuse them!
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb
        type(WeightObj_t) :: doubleWeight
        character(*), parameter :: this_routine = "init_doubleWeight"
        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)

        doubleWeight%dat%F = endFx(csf_i, sOrb)
        doubleWeight%dat%G = endGx(csf_i, sOrb)

        doubleWeight%proc%minus => getMinus_double
        doubleWeight%proc%plus => getPlus_double
        doubleWeight%proc%zero => getZero_double

        doubleWeight%initialized = .true.

    end function init_doubleWeight


    function getMinus_double(nSwitches, bVal, double) result(minusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: double
        real(dp) :: minusWeight
        character(*), parameter :: this_routine = "getMinus_double"
        ASSERT(nSwitches >= 0.0_dp)

        minusWeight = double%F + nSwitches / max(1.0_dp, bVal)

        ASSERT(minusWeight >= 0.0_dp)
    end function getMinus_double

    function getPlus_double(nSwitches, bVal, double) result(plusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: double
        real(dp) :: plusWeight
        character(*), parameter :: this_routine = "getPlus_double"
        ASSERT(nSwitches >= 0.0_dp)

        ! update: the correct check for the b value in the case of the +2
        ! double excitation branch should be that its not allowed to be
        ! less than 2 on the current orbital with the new bVector
        ! implementation:
        ! or otherwise the excitation would lead to negative be values
        if (bVal < 2.0_dp) then
            plusWeight = 0.0_dp
        else
            plusWeight = double%G + nSwitches / bVal
        end if

        ASSERT(plusWeight >= 0.0_dp)
    end function getPlus_double

    function get_forced_zero_double(negSwitches, posSwitches, bVal, double) &
        result(zeroWeight)
        real(dp), intent(in) :: posSwitches, negSwitches, bVal
        type(WeightData_t), intent(in) :: double
        real(dp) :: zeroWeight

        ! remove the order(1) branch as we want to switch at the end!
        if (near_zero(bVal)) then
            zeroWeight = negSwitches * double%G + posSwitches * double%F
        else
            zeroWeight = 1.0_dp / bVal * (negSwitches * double%G + posSwitches * double%F)
        end if

    end function get_forced_zero_double

    function getZero_double(negSwitches, posSwitches, bVal, double) &
        result(zeroWeight)
        real(dp), intent(in) :: posSwitches, negSwitches, bVal
        type(WeightData_t), intent(in) :: double
        real(dp) :: zeroWeight
        character(*), parameter :: this_routine = "getZero_double"
        ASSERT(negSwitches >= 0.0_dp)
        ASSERT(posSwitches >= 0.0_dp)

        ! UPDATE: cant set it to one, because there are cases, when a 0
        ! branch comes to a 2 in a double excitation, where one has to
        ! compare against the -2 branch, which can also be non-zero
        ! so to have the correct weights, just ignore b
        zeroWeight = 1.0_dp + &
                     (negSwitches * double%G + posSwitches * double%F) / max(1.0_dp, bVal)

        ASSERT(zeroWeight >= 0.0_dp)
    end function getZero_double

    function getMinus_semiStart(nSwitches, bVal, semiStart) result(minusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: semiStart
        real(dp) :: minusWeight
        character(*), parameter :: this_routine = "getMinus_semiStart"

        ASSERT(nSwitches >= 0.0_dp)
        ! change b value treatment, by just checking if excitations is
        ! technically possible if b == 0

        minusWeight = semiStart%F * semiStart%zero + semiStart%G * semiStart%minus + &
                      nSwitches / max(1.0_dp, bval) * (semiStart%F * semiStart%plus + semiStart%G * semiStart%zero)

        ASSERT(minusWeight >= 0.0_dp)

    end function getMinus_semiStart

    function getPlus_semiStart(nSwitches, bVal, semiStart) result(plusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: semiStart
        real(dp) :: plusWeight
        character(*), parameter :: this_routine = "getPlus_semiStart"

        ASSERT(nSwitches >= 0.0_dp)
        ! just set +1 branch probability to 0 if b == 0

        if (near_zero(bVal)) then
            plusWeight = 0.0_dp
        else
            plusWeight = semiStart%G * semiStart%zero + semiStart%F * semiStart%plus + &
                         nSwitches / bVal * (semiStart%G * semiStart%minus + semiStart%F * semiStart%zero)
        end if

        ASSERT(plusWeight >= 0.0_dp)

    end function getPlus_semiStart
end module guga_matrixElements
