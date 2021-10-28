#include "macros.h"
! GUGA module containg as much matrix element calculation functionality as
! possible.
module guga_matrixElements
    use SystemData, only: nEl, nBasis, ECore, t_tJ_model, t_heisenberg_model, t_new_hubbard, &
        t_new_real_space_hubbard, treal, nSpatOrbs, t_mixed_hubbard
    use constants, only: dp, n_int, hel_zero, int_rdm, bn2_, Root2
    use bit_reps, only: decode_bit_det
    use OneEInts, only: GetTMatEl
    use procedure_pointers, only: get_umat_el
    use guga_bitRepOps, only: isDouble, calcB_vector_nI, isProperCSF_nI, CSF_Info_t, &
        convert_ilut_toGUGA, identify_excitation, findFirstSwitch, findLastSwitch, &
        calcb_vector_ilut, count_open_orbs_ij
    use guga_bitRepOps, only: contract_1_rdm_ind, contract_2_rdm_ind
    use util_mod, only: binary_search, operator(.isclose.), stop_all, near_zero
    use guga_data, only: projE_replica, ExcitationInformation_t, excit_type, gen_type

    use guga_data, only: funA_0_2overR2, minFunA_2_0_overR2, funA_2_0_overR2, &
                         funA_m1_1_overR2, funA_3_1_overR2, minFunA_0_2_overR2

    use guga_data, only: getdoublematrixelement, getmixedfullstop, getSingleMatrixElement, &
        getdoublecontribution

    use guga_procedure_pointers, only: &
        calc_mixed_start_l2r_contr, calc_mixed_end_l2r_contr, calc_mixed_end_r2l_contr, &
        calc_mixed_start_r2l_contr

    use bit_rep_data, only: niftot, nifd
    use MPI_wrapper, only: iprocindex
    use CalcData, only: matele_cutoff, t_matele_cutoff
    use bit_rep_data, only: nifguga
    use DetBitOps, only: DetBitEQ


    use FciMCData, only: tFillingStochRDMOnFly

    better_implicit_none

    private
    public :: calc_guga_matrix_element, calcMixedContribution, calc_integral_contribution_single
    public :: calcDiagMatEleGuga_nI, calcDiagExchangeGUGA_nI, calcDiagMatEleGUGA_ilut

contains

    subroutine calc_guga_matrix_element(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, t_hamil, &
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

#ifdef DEBUG_
        if (present(rdm_ind)) then
            ASSERT(present(rdm_mat))
        end if
        if (present(rdm_mat)) then
            ASSERT(present(rdm_ind))
        end if
#endif


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
            call calc_fullstart_fullstop_alike_ex(csf_i, csf_j, excitInfo, &
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



    function calcDiagMatEleGuga_nI(nI) result(hel_ret)
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

    function calcDiagMatEleGuga_ilut(ilut) result(hElement)
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

    subroutine calc_single_excitation_ex(csf_i, csf_j, excitInfo, mat_ele, &
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

        ! set defaults for the output if we excit early..
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (i => excitInfo%i, j => excitInfo%j, st => excitInfo%fullstart, &
                   en => excitInfo%fullEnd, gen => excitInfo%currentGen)

            if (present(rdm_ind) .or. present(rdm_mat)) then
                ASSERT(present(rdm_ind))
                ASSERT(present(rdm_mat))
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

    subroutine calc_single_overlap_mixed_ex(csf_i, csf_j, excitInfo, mat_ele, &
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

        def_default(t_calc_full_, t_calc_full, .true.)

        ! set some defaults in case of early exit
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, st => excitInfo%fullStart, &
                   ss => excitInfo%secondStart, en => excitInfo%fullEnd, &
                   gen => excitInfo%firstGen, fe => excitInfo%firstEnd, &
                   typ => excitInfo%typ)

            if (present(rdm_ind) .or. present(rdm_mat)) then
                ASSERT(present(rdm_ind))
                ASSERT(present(rdm_mat))
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

    subroutine calc_normal_double_ex(csf_i, csf_j, excitInfo, mat_ele, &
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

            if (present(rdm_ind) .or. present(rdm_mat)) then
                ASSERT(present(rdm_ind))
                ASSERT(present(rdm_mat))
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

    subroutine calc_fullstop_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
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

        ! set some defaults in case of early exit
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, typ => excitInfo%typ, st => excitInfo%fullStart, &
                   se => excitInfo%secondStart, gen => excitInfo%gen1, &
                   en => excitInfo%fullEnd)

            if (present(rdm_ind) .or. present(rdm_mat)) then
                ASSERT(present(rdm_ind))
                ASSERT(present(rdm_mat))
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

    subroutine calc_fullstart_alike_ex(csf_i, csf_j, excitInfo, mat_ele, &
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

        ! set defaults for early exits
        mat_ele = h_cast(0.0_dp)
        delta_b = csf_i%B_int - csf_j%B_int

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, start => excitInfo%fullstart, &
                   ende => excitInfo%fullEnd, semi => excitInfo%firstEnd, &
                   gen => excitInfo%firstGen, typ => excitInfo%typ)

            if (present(rdm_ind) .or. present(rdm_mat)) then
                ASSERT(present(rdm_ind))
                ASSERT(present(rdm_mat))
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

    subroutine calc_fullstart_fullstop_alike_ex(csf_i, csf_j, excitInfo, &
                                                mat_ele, t_hamil, rdm_ind, rdm_mat)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
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

        ! set defaults for early exit
        mat_ele = h_cast(0.0_dp)

        associate (ii => excitInfo%i, jj => excitInfo%j, kk => excitInfo%k, &
                   ll => excitInfo%l, start => excitInfo%fullStart, &
                   ende => excitInfo%fullEnd)

            if (present(rdm_ind) .or. present(rdm_mat)) then
                ASSERT(present(rdm_ind))
                ASSERT(present(rdm_mat))
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

    subroutine calc_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
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

            ! need to input dummy variable into the end contr functions
            ! use temp_mat1
            temp_mat1 = 1.0_dp

            if (t_hamil_ .or. (tFillingStochRDMOnFly .and. present(rdm_mat))) then
                if (typ == excit_type%fullstop_L_to_R) then
                    ! L -> R
                    ! what do i have to put in as the branch pgen?? does it have
                    ! an influence on the integral and matrix element calculation?
                    if (present(rdm_mat)) then
                        call calc_mixed_end_l2r_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                      temp_mat0, integral, rdm_ind, rdm_mat)
                        ! need to multiply by x1
                        rdm_mat = rdm_mat * temp_x1
                    else
                        call calc_mixed_end_l2r_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                      temp_mat0, integral)
                    end if

                    mat_ele = temp_x1 * ((get_umat_el(en, se, st, en) + &
                                          get_umat_el(se, en, en, st)) / 2.0_dp + integral)

                else if (typ == excit_type%fullstop_R_to_L) then
                    ! R -> L
                    if (present(rdm_mat)) then
                        call calc_mixed_end_r2l_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                      temp_mat0, integral, rdm_ind, rdm_mat)
                        rdm_mat = rdm_mat * temp_x1
                    else
                        call calc_mixed_end_r2l_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                      temp_mat0, integral)
                    end if

                    mat_ele = temp_x1 * ((get_umat_el(en, st, se, en) + &
                                          get_umat_el(st, en, en, se)) / 2.0_dp + integral)

                end if
            else
                mat_ele = h_cast(temp_x1)
            end if
        end associate

    end subroutine calc_fullstop_mixed_ex

    subroutine calc_fullstart_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, &
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

            ! need to input variable into start_contr routines, misuse temp_mat1
            temp_mat1 = 1.0_dp

            if (t_hamil .or. (tFillingStochRDMOnFly .and. present(rdm_mat))) then
                if (typ == excit_type%fullstart_L_to_R) then
                    ! L -> R
                    if (present(rdm_mat)) then
                        call calc_mixed_start_l2r_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                        temp_mat0, integral, rdm_ind, rdm_mat)
                        ! need to multiply by guga-mat:
                        rdm_mat = rdm_mat * guga_mat
                    else
                        call calc_mixed_start_l2r_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                        temp_mat0, integral, rdm_ind, rdm_mat)
                    end if

                    mat_ele = guga_mat * ((get_umat_el(st, se, en, st) + &
                                           get_umat_el(se, st, st, en)) / 2.0_dp + integral)

                else if (typ == excit_type%fullstart_R_to_L) then

                    if (present(rdm_mat)) then
                        call calc_mixed_start_r2l_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                        temp_mat0, integral, rdm_ind, rdm_mat)
                        rdm_mat = rdm_mat * guga_mat
                    else
                        call calc_mixed_start_r2l_contr(tmp_I, csf_i, tmp_J, excitInfo, temp_mat1, &
                                                        temp_mat0, integral)
                    end if

                    mat_ele = guga_mat * ((get_umat_el(st, en, se, st) + &
                                           get_umat_el(en, st, st, se)) / 2.0_dp + integral)

                end if
            end if
        end associate

    end subroutine calc_fullstart_mixed_ex

    subroutine calc_fullstart_fullstop_mixed_ex(ilutI, csf_i, ilutJ, csf_j, excitInfo, &
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
                    mat_ele = calcMixedContribution(tmp_I, csf_i, tmp_J, start, ende, &
                                                    rdm_ind, rdm_mat)
                else
                    mat_ele = calcMixedContribution(tmp_I, csf_i, tmp_J, start, ende)
                end if
            end if
        end associate
    end subroutine calc_fullstart_fullstop_mixed_ex


    function calcMixedContribution(ilut, csf_i, t, start, ende, rdm_ind, rdm_mat) result(integral)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: start, ende
        HElement_t(dp) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calcMixedContribution"

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
    end function calcMixedContribution


    subroutine calc_integral_contribution_single(csf_i, csf_j, i, j, st, en, integral)
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




end module guga_matrixElements
