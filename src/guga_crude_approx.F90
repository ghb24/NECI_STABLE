#include "macros.h"
! GUGA excitations module:
! contains as much excitation related functionality as possible
module guga_crude_approx_mod

    use SystemData, only: nEl, G1, nSpatOrbs

    use constants, only: dp, n_int, bits_n_int, bn2_

    use bit_rep_data, only: niftot, nifguga, nifd

    use OneEInts, only: GetTMatEl

    use dSFMT_interface, only: genrand_real2_dSFMT

    use util_mod, only: abs_l1, operator(.isclose.), &
                        operator(.div.), near_zero, stop_all

    use SymExcitDataMod, only: SpinOrbSymLabel, OrbClassCount, SymLabelCounts2, &
                               sym_label_list_spat

    use sym_general_mod, only: ClassCountInd

    use umatcache, only: gtID

    use procedure_pointers, only: get_umat_el

    use guga_procedure_pointers, only: pickOrbitals_double

    use guga_types, only: WeightObj_t

    use guga_excitations, only: global_excitInfo, checkcompatibility

    use guga_data, only: ExcitationInformation_t, excit_type, gen_type

    use guga_bitRepOps, only: isProperCSF_ilut, calcB_vector_ilut, &
                              encode_matrix_element, convert_ilut_toNECI, &
                              CSF_Info_t

    use guga_matrixElements, only: calc_guga_matrix_element

    better_implicit_none

    private
    public :: create_crude_guga_single, create_crude_guga_double, &
        create_crude_single, create_crude_double, perform_crude_excitation

contains

    subroutine create_crude_guga_single(ilut, nI, csf_i, exc, pgen, excitInfo_in)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer(n_int), intent(out) :: exc(0:nifguga)
        real(dp), intent(out) :: pgen
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_in
        character(*), parameter :: this_routine = "create_crude_guga_single"

        type(ExcitationInformation_t) :: excitInfo
        integer :: i, j, st, en, start_d, end_d, gen
        real(dp) :: orb_pgen, r, branch_pgen
        HElement_t(dp) :: mat_ele
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)

        ASSERT(isProperCSF_ilut(ilut))

        ! can i use the original orbital picker? no.. since it allows for
        ! switches..

        if (present(excitInfo_in)) then
            excitInfo = excitInfo_in
        else
            call pick_orbitals_single_crude(ilut, nI, csf_i, excitInfo, orb_pgen)
        end if

        if (.not. excitInfo%valid) then
            ! if no valid indices were picked, return 0 excitation and return
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        ! reimplement it from scratch
        i = excitInfo%i
        j = excitInfo%j

        ! try to make a valid, determinant-like single excitation,
        ! respecting the spin-character of the CSF
        st = excitInfo%fullStart
        en = excitInfo%fullEnd
        gen = excitInfo%currentGen
        start_d = csf_i%stepvector(st)
        end_d = csf_i%stepvector(en)

        branch_pgen = 1.0_dp

        exc = ilut
        associate(b => csf_i%B_int)
            if (start_d == 3) then
                ! here we now it is a lowering generator with d_j = 1 or 2
                ! with restrictions however
                ! here a d_j = 2 is theoretically possible
                ASSERT(gen == gen_type%L)
                ASSERT(end_d /= 3)
                if (end_d == 0) then
                    if (all(b(st:en - 1) > 0)) then
                        ! in this case both starts are possible!
                        r = genrand_real2_dSFMT()

                        if (r < 0.5_dp) then
                            ! make a 3 -> 1 and 0 -> 2
                            set_zero(exc, st)
                            set_one(exc, st)

                            set_two(exc, en)

                        else
                            ! make 3 -> 2 and 0 -> 1
                            set_zero(exc, st)
                            set_two(exc, st)

                            set_one(exc, en)

                        end if
                        branch_pgen = 0.5_dp
                    else
                        ! here only 3 > 1 and 0 > 2 is possible
                        ! make
                        set_zero(exc, st)
                        set_one(exc, st)

                        set_two(exc, en)

                    end if

                else if (end_d == 1) then
                    ! then only the 3 > 1 start is possible and b is irrelevant
                    ! make 3 > 1 and 1 > 3
                    set_zero(exc, st)
                    set_one(exc, st)

                    set_three(exc, en)

                else if (end_d == 2) then
                    ! this is only possible if all b are > 0
                    if (all(b(st:en - 1) > 0)) then
                        ! make 3 > 2 and 2  > 3
                        set_zero(exc, st)
                        set_two(exc, st)

                        set_three(exc, en)

                    else
                        pgen = 0.0_dp
                        exc = 0_n_int
                        return
                    end if
                end if

            else if (start_d == 0) then
                ! here we know it is a raising generator
                ASSERT(gen == gen_type%R)
                ASSERT(end_d /= 0)
                if (end_d == 3) then
                    if (all(b(st:en - 1) > 0)) then
                        r = genrand_real2_dSFMT()

                        if (r < 0.5_dp) then
                            ! make 0 > 1 and  3 > 2
                            set_one(exc, st)

                            set_zero(exc, en)
                            set_two(exc, en)

                        else
                            ! make 0 > 2 and 3 > 1
                            set_two(exc, st)

                            set_zero(exc, en)
                            set_one(exc, en)

                        end if
                        branch_pgen = 0.5_dp
                    else
                        ! make  0 > 1 and 3 > 2
                        set_one(exc, st)

                        set_zero(exc, en)
                        set_two(exc, en)

                    end if

                else if (end_d == 1) then
                    ! make 0 > 1 and 1 > 0

                    set_one(exc, st)

                    set_zero(exc, en)

                else if (end_d == 2) then
                    if (all(b(st:en - 1) > 0)) then
                        ! make 0 > 2 and 2 > 0
                        set_two(exc, st)

                        set_zero(exc, en)

                    else
                        pgen = 0.0_dp
                        exc = 0_n_int
                        return
                    end if
                end if
            else if (start_d == 1) then
                if (all(b(st:en - 1) > 0) .and. end_d /= 1) then
                    ! only in this case it is possible
                    if (end_d == 2) then
                        if (gen == gen_type%R) then
                            ! make 1 > 3 and 2 > 0
                            set_three(exc, st)

                            set_zero(exc, en)

                        else if (gen == gen_type%L) then
                            ! make 1 > 0 and 2 > 3
                            set_zero(exc, st)

                            set_three(exc, en)

                        end if

                    else if (end_d == 3) then
                        ASSERT(gen == gen_type%R)
                        ! make 1 > 3 and 3 > 1
                        set_three(exc, st)

                        set_zero(exc, en)
                        set_one(exc, en)

                    else if (end_d == 0) then
                        ASSERT(gen == gen_type%L)
                        ! make 1 > 0 and 0 > 1
                        set_zero(exc, st)

                        set_one(exc, en)

                    end if
                else
                    pgen = 0.0_dp
                    exc = 0_n_int
                    return
                end if

            else if (start_d == 2) then
                ! here I do not have a b restriction or?
                if (end_d == 0) then
                    ASSERT(gen == gen_type%L)
                    ! make 2 > 0 and 0 > 2
                    set_zero(exc, st)

                    set_two(exc, en)

                else if (end_d == 3) then
                    ASSERT(gen == gen_type%R)
                    ! make 2 > 3 and 3 > 2
                    set_three(exc, st)

                    set_zero(exc, en)
                    set_two(exc, en)

                else if (end_d == 1) then
                    if (gen == gen_type%R) then
                        ! make 2 > 3 and 1 > 0
                        set_three(exc, st)

                        set_zero(exc, en)

                    else if (gen == gen_type%L) then
                        ! make 2 > 0 and 1 > 3
                        set_zero(exc, st)

                        set_three(exc, en)

                    end if
                else if (end_d == 2) then
                    pgen = 0.0_dp
                    exc = 0_n_int
                    return
                end if
            end if
        end associate

        ! we also want to check if we produced a valid CSF here
        if (.not. isProperCSF_ilut(exc, .true.)) then
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, .true.)

        if (near_zero(mat_ele)) then
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        call encode_matrix_element(exc, 0.0_dp, 2)
        call encode_matrix_element(exc, mat_ele, 1)

        global_excitInfo = excitInfo

        if (present(excitInfo_in)) then
            ! then the orbitals were already picked before and we only want
            ! to give the branch_pgen here
            pgen = branch_pgen
        else
            ! otherwise we want the full pgen
            pgen = orb_pgen * branch_pgen
        end if

    end subroutine create_crude_guga_single


    subroutine create_crude_guga_double(ilut, nI, csf_i, exc, pgen, excitInfo_in)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer(n_int), intent(out) :: exc(0:nifguga)
        real(dp), intent(out) :: pgen
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_in

        type(ExcitationInformation_t) :: excitInfo
        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        real(dp) :: branch_pgen, orb_pgen
        HElement_t(dp) :: mat_ele
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(WeightObj_t) :: weights
        integer :: elecs(2), orbs(2)

        if (present(excitInfo_in)) then
            excitInfo = excitInfo_in
        else
            call pickOrbitals_double(ilut, nI, csf_i, excitInfo, orb_pgen)

            if (.not. excitInfo%valid) then
                ! if no valid indices were picked, return 0 excitation and return
                exc = 0_n_int
                pgen = 0.0_dp
                return
            end if
            call checkCompatibility(&
                    csf_i, excitInfo, compFlag, posSwitches, negSwitches, weights)

            if (.not. compFlag) then
                exc = 0_n_int
                pgen = 0.0_dp
                return
            end if

        end if

        ! here I have to do the actual crude double excitation..
        ! my idea for now is to create pseudo random spin-orbital from the
        ! picked spatial orbitals
        call create_random_spin_orbs(ilut, csf_i, excitInfo, elecs, orbs, branch_pgen)

        if (any(elecs == 0) .or. any(orbs == 0)) then
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        if (near_zero(branch_pgen)) then
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        exc = ilut

        clr_orb(exc, elecs(1))
        clr_orb(exc, elecs(2))
        set_orb(exc, orbs(1))
        set_orb(exc, orbs(2))

        ! this check is the same at the end of singles:
        ! maybe I also need to do this only if no excitInfo_in is provided..
        ! we also want to check if we produced a valid CSF here
        if (.not. isProperCSF_ilut(exc, .true.)) then
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        ! we also need to calculate the matrix element here!
        call convert_ilut_toNECI(ilut, ilutI)
        call convert_ilut_toNECI(exc, ilutJ)
        call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, .true.)

        if (near_zero(mat_ele)) then
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        call encode_matrix_element(exc, 0.0_dp, 2)
        call encode_matrix_element(exc, mat_ele, 1)

        global_excitInfo = excitInfo

        if (present(excitInfo_in)) then
            ! then the orbitals were already picked before and we only want
            ! to give the branch_pgen here
            pgen = branch_pgen
        else
            ! otherwise we want the full pgen
            pgen = orb_pgen * branch_pgen
        end if

    end subroutine create_crude_guga_double

    subroutine create_crude_single(ilut, csf_i, exc, branch_pgen, excitInfo_in)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer(n_int), intent(inout) :: exc(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_in
        character(*), parameter :: this_routine = "create_crude_single"

        type(ExcitationInformation_t) :: excitInfo
        integer :: elec, orb
        real(dp) :: r

#ifdef DEBUG_
        ! also assert we are not calling it for a weight gen. accidently
        ASSERT(excitInfo%currentGen /= 0)
        ASSERT(isProperCSF_ilut(ilut))
        ! also check if calculated b vector really fits to ilut
        ASSERT(all(csf_i%B_real.isclose.calcB_vector_ilut(ilut(0:nifd))))
        if (excitInfo%currentGen == gen_type%R) then
            ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        else if (excitInfo%currentGen == gen_type%L) then
            ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        end if
#endif

        if (present(excitInfo_in)) then
            excitInfo = excitInfo_in
        else
            call stop_all(this_routine, "not yet implemented without excitInfo_in")
        end if

        elec = excitInfo%j
        orb = excitInfo%i

        exc = ilut

        select case (csf_i%stepvector(elec))

        case (0)
            call stop_all(this_routine, "empty orbital picked as electron!")

        case (1)

            clr_orb(exc, 2 * elec - 1)
            branch_pgen = 1.0_dp

            if (csf_i%stepvector(orb) == 0) then

                set_orb(exc, 2 * orb - 1)

            else if (csf_i%stepvector(orb) == 1) then

                branch_pgen = 0.0_dp

            else if (csf_i%stepvector(orb) == 2) then

                set_orb(exc, 2 * orb - 1)

            end if

        case (2)

            clr_orb(exc, 2 * elec)
            branch_pgen = 1.0_dp

            if (csf_i%stepvector(orb) == 0) then

                set_orb(exc, 2 * orb)

            else if (csf_i%stepvector(orb) == 1) then

                set_orb(exc, 2 * orb)

            else if (csf_i%stepvector(orb) == 2) then

                branch_pgen = 0.0_dp

            end if

        case (3)

            if (csf_i%stepvector(orb) == 0) then

                ! here i have to decide..
                branch_pgen = 0.5_dp

                r = genrand_real2_dSFMT()

                if (r < 0.5_dp) then

                    ! 1 -> 2
                    clr_orb(exc, 2 * elec)
                    set_orb(exc, 2 * orb)

                else
                    ! 2 -> 1
                    clr_orb(exc, 2 * elec - 1)
                    set_orb(exc, 2 * orb - 1)

                end if

            else
                branch_pgen = 1.0_dp

                if (csf_i%stepvector(orb) == 1) then

                    clr_orb(exc, 2 * elec)
                    set_orb(exc, 2 * orb)

                else if (csf_i%stepvector(orb) == 2) then

                    clr_orb(exc, 2 * elec - 1)
                    set_orb(exc, 2 * orb - 1)

                end if
            end if
        end select

        if (.not. isProperCSF_ilut(exc)) then
            ! i have to check if i created a proper CSF
            branch_pgen = 0.0_dp
        end if

    end subroutine create_crude_single


    subroutine create_crude_double(ilut, exc, branch_pgen, excitInfo_in)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer(n_int), intent(inout) :: exc(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        type(ExcitationInformation_t), intent(in), optional :: excitInfo_in
        character(*), parameter :: this_routine = "create_crude_double"

        type(ExcitationInformation_t) :: excitInfo

        call stop_all(this_routine, "not yet implemented ")

        if (present(excitInfo_in)) then
            excitInfo = excitInfo_in
        else
            call stop_all(this_routine, "not yet implemented without excitInfo_in")
        end if

        ! i think i still have to reuse the excitInfo information of the
        ! excitation type.. do I still know where the holes and electrons are
        ! in the double excitations?
        exc = ilut

        select case (excitInfo%typ)

            ! maybe i can combine some of them together
            ! i think i can.. but i have to think about that more clearly!
        case (excit_type%single_overlap_L_to_R, &
              excit_type%single_overlap_R_to_L)
            ! single overlap excitation

        end select

        if (.not. isProperCSF_ilut(exc)) then
            ! i have to check if i created a proper CSF
            branch_pgen = 0.0_dp
        end if

    end subroutine create_crude_double

    subroutine create_random_spin_orbs(ilut, csf_i, excitInfo, elecs, orbs, pgen)
        ! a subroutine to create random spin-orbitals from chosen
        ! spatial orbital of a GUGA excitation.
        ! this is needed in the crude back-spawn approximation to
        ! create determinant-like excitations
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer, intent(out) :: elecs(2)
        integer, intent(out) :: orbs(2)
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "create_random_spin_orbs"

        integer :: elec_1, orb_1, elec_2, orb_2, d_elec, s_elec, d_orb, s_orb
        real(dp) :: r

        ! i essentially need electron and empty spin-orbitals to
        ! create the excitation. the validity of the CSF will be checked
        ! afterwards and then maybe thrown away..

        elecs = 0
        orbs = 0
        pgen = 1.0_dp

        if (excitInfo%typ == excit_type%single) then
            ! the case of a single excitation
            elec_1 = excitInfo%j
            orb_1 = excitInfo%i

            select case (csf_i%stepvector(elec_1))

            case (0)
                call stop_all(this_routine, "empty elec index")

            case (1)
                ! 1 corresponds to a beta orbital
                elecs(1) = 2 * elec_1 - 1
                orbs(1) = 2 * orb_1 - 1

            case (2)
                ! 2 corresponds to alpha orbitals
                elecs(1) = 2 * elec_1
                orbs(1) = 2 * orb_1

            case (3)
                select case (csf_i%stepvector(orb_1))
                case (0)
                    ! now we can have two possibilities
                    r = genrand_real2_dSFMT()

                    if (r < 0.5_dp) then
                        elecs(1) = 2 * elec_1 - 1
                        orbs(1) = 2 * orb_1 - 1
                    else
                        elecs(1) = 2 * elec_1
                        orbs(1) = 2 * orb_1
                    end if

                    pgen = 0.5_dp

                case (1)
                    elecs(1) = 2 * elec_1
                    orbs(1) = 2 * orb_1

                case (2)
                    elecs(1) = 2 * elec_1 - 1
                    orbs(1) = 2 * orb_1 - 1

                end select
            end select

        else
            elec_1 = excitInfo%j
            elec_2 = excitInfo%l

            orb_1 = excitInfo%i
            orb_2 = excitInfo%k

            ASSERT(csf_i%stepvector(elec_1) /= 0)
            ASSERT(csf_i%stepvector(elec_2) /= 0)
            ASSERT(csf_i%stepvector(orb_1) /= 3)
            ASSERT(csf_i%stepvector(orb_2) /= 3)

            select case (excitInfo%typ)
            case (excit_type%single_overlap_L_to_R, &
                  excit_type%fullstop_lowering, &
                  excit_type%fullstart_raising)
                ! here i know the orbital indices are identical

                ASSERT(orb_1 == orb_2)
                orbs(1) = 2 * orb_1 - 1
                orbs(2) = 2 * orb_1

                ! write a general function which gives me valid spin-orbs
                ! for GUGA CSFs (mayb as an input the 'neutral' number 3 here.c.
                ! maybe later..

                if ((csf_i%stepvector(elec_1) == csf_i%stepvector(elec_2)) &
                    .and. csf_i%Occ_int(elec_1) == 1) then
                    ! in this case no crude excitation is possible
                    pgen = 0.0_dp
                    elecs = 0
                    orbs = 0
                    return
                end if

                select case (csf_i%stepvector(elec_1))

                case (0)
                    call stop_all(this_routine, "something wrong happened")

                case (1)

                    elecs(1) = 2 * elec_1 - 1
                    ! then the second must be of 'opposite' spin
                    elecs(2) = 2 * elec_2

                case (2)
                    elecs(1) = 2 * elec_1
                    elecs(2) = 2 * elec_2 - 1

                case (3)
                    if (csf_i%stepvector(elec_2) == 3) then
                        ! here we have a choice
                        r = genrand_real2_dSFMT()

                        if (r < 0.5_dp) then
                            elecs(1) = 2 * elec_1 - 1
                            elecs(2) = 2 * elec_2
                        else
                            elecs(1) = 2 * elec_1
                            elecs(2) = 2 * elec_2 - 1
                        end if

                        pgen = 0.5_dp
                    else if (csf_i%stepvector(elec_2) == 1) then
                        elecs(2) = 2 * elec_2 - 1
                        elecs(1) = 2 * elec_1

                    else if (csf_i%stepvector(elec_2) == 2) then
                        elecs(2) = 2 * elec_2
                        elecs(1) = 2 * elec_1 - 1

                    end if
                end select

            case (excit_type%single_overlap_R_to_L, &
                  excit_type%fullstop_raising, &
                  excit_type%fullstart_lowering)

                ! here i know the electron indices are the same
                ASSERT(elec_1 == elec_2)

                elecs(1) = 2 * elec_1 - 1
                elecs(2) = 2 * elec_1

                if ((csf_i%stepvector(orb_1) == csf_i%stepvector(orb_2)) &
                    .and. csf_i%Occ_int(orb_1) == 1) then
                    pgen = 0
                    elecs = 0
                    orbs = 0
                    return
                end if

                select case (csf_i%stepvector(orb_1))
                case (3)
                    call stop_all(this_routine, "something went wrong")

                case (1)

                    orbs(1) = 2 * orb_1
                    orbs(2) = 2 * orb_2 - 1

                case (2)

                    orbs(1) = 2 * orb_1 - 1
                    orbs(2) = 2 * orb_2

                case (0)
                    if (csf_i%stepvector(orb_2) == 0) then
                        r = genrand_real2_dSFMT()

                        if (r < 0.5_dp) then
                            orbs(1) = 2 * orb_1 - 1
                            orbs(2) = 2 * orb_2
                        else
                            orbs(1) = 2 * orb_1
                            orbs(2) = 2 * orb_2 - 1
                        end if
                        pgen = 0.5_dp

                    else if (csf_i%stepvector(orb_2) == 1) then
                        orbs(2) = 2 * orb_2
                        orbs(1) = 2 * orb_1 - 1

                    else if (csf_i%stepvector(orb_2) == 2) then
                        orbs(2) = 2 * orb_2 - 1
                        orbs(1) = 2 * orb_1

                    end if
                end select

                ! do the same as above, just for two hole indices!

            case (excit_type%fullstop_L_to_R, &
                  excit_type%fullstop_R_to_L, &
                  excit_type%fullStart_L_to_R, &
                  excit_type%fullstart_R_to_L)

                ! here i know one electron and hole index are the same
                ASSERT(elec_1 /= elec_2)
                ASSERT(orb_1 /= orb_2)

                ! this case is very restrictive..
                if (elec_1 == orb_1) then

                    s_elec = elec_1
                    s_orb = orb_1

                    d_elec = elec_2
                    d_orb = orb_2

                else if (elec_1 == orb_2) then

                    s_elec = elec_1
                    s_orb = orb_1

                    d_elec = elec_2
                    d_orb = orb_1

                else if (elec_2 == orb_1) then

                    s_elec = elec_2
                    s_orb = orb_1

                    d_elec = elec_1
                    d_orb = orb_2

                else if (elec_2 == orb_2) then

                    s_elec = elec_2
                    s_orb = orb_2

                    d_elec = elec_1
                    d_orb = orb_1

#ifdef DEBUG_
                else
                    call stop_all(this_routine, "something went wrong")
#endif
                end if

                ASSERT(csf_i%Occ_int(s_elec) == 1)

                if (csf_i%stepvector(s_elec) == 1) then
                    if (IsOcc(ilut, 2 * d_elec) .and. IsNotOcc(ilut, 2 * d_orb - 1)) then
                        ! this is the only case this is possible
                        elecs(1) = 2 * s_elec - 1
                        elecs(2) = 2 * d_elec

                        orbs(1) = 2 * s_orb
                        orbs(2) = 2 * d_orb - 1

                    else
                        pgen = 0.0_dp
                    end if
                else if (csf_i%stepvector(s_elec) == 2) then
                    if (IsOcc(ilut, 2 * d_elec - 1) .and. IsNotOcc(ilut, 2 * d_orb)) then
                        elecs(1) = 2 * s_elec
                        elecs(2) = 2 * d_elec - 1

                        orbs(1) = 2 * s_orb - 1
                        orbs(2) = 2 * d_orb
                    else
                        pgen = 0.0_dp
                    end if
                end if

            case (excit_type%fullstart_stop_mixed)
                ! full start full stop mixed
                ASSERT(elec_1 /= elec_2)
                ASSERT(orb_1 /= orb_2)
                ASSERT(csf_i%Occ_int(elec_1) == 1)
                ASSERT(csf_i%Occ_int(elec_2) == 1)

                if (csf_i%stepvector(elec_1) == csf_i%stepvector(elec_2)) then
                    pgen = 0.0_dp
                else if (csf_i%stepvector(elec_1) == 1) then
                    elecs(1) = 2 * elec_1 - 1
                    elecs(2) = 2 * elec_2

                    orbs(1) = 2 * elec_2 - 1
                    orbs(2) = 2 * elec_1

                else if (csf_i%stepvector(elec_1) == 2) then
                    elecs(1) = 2 * elec_1
                    elecs(2) = 2 * elec_2 - 1

                    orbs(1) = 2 * elec_2
                    orbs(2) = 2 * elec_1 - 1
                end if

            case (excit_type%fullstart_stop_alike)
                ! full-start full-stop alike
                ASSERT(elec_1 == elec_2)
                ASSERT(orb_1 == orb_2)

                elecs(1) = 2 * elec_1 - 1
                elecs(2) = 2 * elec_1

                orbs(1) = 2 * orb_1 - 1
                orbs(2) = 2 * orb_1

            case default
                ! now the general 4-index double excitations..
                ! this can be nasty again..

                call pick_random_4ind(csf_i, elec_1, elec_2, orb_1, orb_2, elecs, orbs, pgen)

            end select
        end if

    end subroutine create_random_spin_orbs



    subroutine pick_orbitals_single_crude(ilut, nI, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pick_orbitals_single_crude"

        integer :: elec, cc_i, ierr, nOrb, orb_i
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: elec_factor

        unused_var(ilut)

        ! first pick completely random from electrons only!
        elec = 1 + floor(genrand_real2_dSFMT() * nEl)
        ! have to adjust pgen if it is a doubly occupied orbital afterwards
        ! -> since twice the chance to pick that orbital then!

        ! pick associated "spin orbital"
        orb_i = nI(elec)

        ! get the symmetry index:
        ! since there is no spin restriction here have to consider both
        ! again
        cc_i = ClassCountInd(1, SpinOrbSymLabel(orb_i), G1(orb_i)%Ml)

        ! get the number of orbitals in this symmetry sector
        nOrb = OrbClassCount(cc_i)
        allocate(cum_arr(nOrb), stat=ierr)
        select case (csf_i%stepvector(gtID(orb_i)))
            ! der stepvalue sagt mir auch, ob es ein alpha oder beta
            ! elektron war..
            ! i have to change this for the crude implementation
        case (1)
            elec_factor = 1.0_dp
            call gen_crude_guga_single_1(nI, csf_i, orb_i, cc_i, cum_arr)

        case (2)
            ! to do
            elec_factor = 1.0_dp
            call gen_crude_guga_single_2(nI, csf_i, orb_i, cc_i, cum_arr)

        case (3)
            ! adjust pgen, the chance to pick a doubly occupied with
            ! spinorbitals is twice as high..
            elec_factor = 2.0_dp
            call gen_crude_guga_single_3(nI, csf_i, orb_i, cc_i, cum_arr)

        case default
            call stop_all(this_routine, "should not have picked empty orbital")

        end select

        call stop_all(this_routine, "TODO")
        pgen = 0.0_dp

    end subroutine pick_orbitals_single_crude

    subroutine perform_crude_excitation(ilut, csf_i, excitInfo, excitation, compFlag)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: excitation(0:nifguga)
        logical, intent(out) :: compFlag

        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        HElement_t(dp) :: mat_ele
        type(ExcitationInformation_t) :: dummy

        excitation = ilut

        select case (excitInfo%typ)
        case (excit_type%fullstart_stop_mixed)
            ! fully exchange is easy, just switch involved step-vectors
            if (csf_i%stepvector(excitInfo%fullstart) == 1) then
                if (csf_i%stepvector(excitInfo%fullEnd) == 1) then
                    ! not valid if the same step-vectors
                    compFlag = .false.
                    return
                end if

                ! 1 -> 2 at start
                set_two(excitation, excitInfo%fullstart)
                clr_one(excitation, excitInfo%fullstart)

                ! 2 -> 1 at end
                set_one(excitation, excitInfo%fullEnd)
                clr_two(excitation, excitInfo%fullEnd)

            else
                if (csf_i%stepvector(excitInfo%fullEnd) == 2) then
                    compFlag = .false.
                    return
                end if

                ! 2 > 1 at start
                set_one(excitation, excitInfo%fullStart)
                clr_two(excitation, excitInfo%fullstart)

                ! 1 -> 2 at end
                clr_one(excitation, excitInfo%fullEnd)
                set_two(excitation, excitInfo%fullEnd)

            end if

        case (excit_type%fullstop_R_to_L)
            ! full stop raising into lowering
            ! the start and semi-start step-values have to be
            ! different than the full-stop, where a switch is enforced.
            if (csf_i%stepvector(excitInfo%fullEnd) == 1) then
                ! the full-start and semi-start are not allowed to have
                ! the same step-vector as the full-stop
                if (csf_i%stepvector(excitInfo%secondStart) == 1 &
                    .or. csf_i%stepvector(excitInfo%fullStart) == 1) then
                    compFlag = .false.
                    return
                end if

                ! in the case that the end is d = 1:
                set_one(excitation, excitInfo%fullStart)

                clr_two(excitation, excitInfo%secondStart)

                clr_one(excitation, excitInfo%fullEnd)
                set_two(excitation, excitInfo%fullEnd)

            else
                if (csf_i%stepvector(excitInfo%secondStart) == 2 &
                    .or. csf_i%stepvector(excitInfo%fullStart) == 2) then
                    compFlag = .false.
                    return
                end if

                ! in the case of d = 2 at end

                ! in the case that the end is d = 1:
                set_two(excitation, excitInfo%fullStart)

                clr_one(excitation, excitInfo%secondStart)

                clr_two(excitation, excitInfo%fullEnd)
                set_one(excitation, excitInfo%fullEnd)

            end if

        case (excit_type%fullstop_L_to_R)
            ! full-stop lowering into raising

            if (csf_i%stepvector(excitInfo%fullEnd) == 1) then
                ! the full-start and semi-start are not allowed to have
                ! the same step-vector as the full-stop
                if (csf_i%stepvector(excitInfo%secondStart) == 1 &
                    .or. csf_i%stepvector(excitInfo%fullStart) == 1) then
                    compFlag = .false.
                    return
                end if

                ! in the case that the end is d = 1:
                set_one(excitation, excitInfo%secondStart)

                clr_two(excitation, excitInfo%fullStart)

                clr_one(excitation, excitInfo%fullEnd)
                set_two(excitation, excitInfo%fullEnd)

            else
                if (csf_i%stepvector(excitInfo%secondStart) == 2 &
                    .or. csf_i%stepvector(excitInfo%fullStart) == 2) then
                    compFlag = .false.
                    return
                end if

                ! in the case of d = 2 at end

                ! in the case that the end is d = 1:
                set_two(excitation, excitInfo%secondStart)

                clr_one(excitation, excitInfo%fullStart)

                clr_two(excitation, excitInfo%fullEnd)
                set_one(excitation, excitInfo%fullEnd)

            end if

        case (excit_type%fullStart_L_to_R)
            ! full-start lowering into raising

            if (csf_i%stepvector(excitInfo%fullStart) == 1) then
                ! the full-start and semi-start are not allowed to have
                ! the same step-vector as the full-stop
                if (csf_i%stepvector(excitInfo%firstEnd) == 1 &
                    .or. csf_i%stepvector(excitInfo%fullEnd) == 1) then
                    compFlag = .false.
                    return
                end if

                ! in the case that the end is d = 1:
                set_one(excitation, excitInfo%firstEnd)

                clr_two(excitation, excitInfo%fullEnd)

                clr_one(excitation, excitInfo%fullStart)
                set_two(excitation, excitInfo%fullStart)

            else
                if (csf_i%stepvector(excitInfo%firstEnd) == 2 &
                    .or. csf_i%stepvector(excitInfo%fullEnd) == 2) then
                    compFlag = .false.
                    return
                end if

                ! in the case of d = 2 at end

                ! in the case that the end is d = 1:
                set_two(excitation, excitInfo%firstEnd)

                clr_one(excitation, excitInfo%fullEnd)

                clr_two(excitation, excitInfo%fullStart)
                set_one(excitation, excitInfo%fullStart)

            end if

        case (excit_type%fullstart_R_to_L)
            ! full-start raising into lowering

            if (csf_i%stepvector(excitInfo%fullStart) == 1) then
                ! the full-start and semi-start are not allowed to have
                ! the same step-vector as the full-stop
                if (csf_i%stepvector(excitInfo%firstEnd) == 1 &
                    .or. csf_i%stepvector(excitInfo%fullEnd) == 1) then
                    compFlag = .false.
                    return
                end if

                ! in the case that the end is d = 1:
                set_one(excitation, excitInfo%fullEnd)

                clr_two(excitation, excitInfo%firstEnd)

                clr_one(excitation, excitInfo%fullStart)
                set_two(excitation, excitInfo%fullStart)

            else
                if (csf_i%stepvector(excitInfo%firstEnd) == 2 &
                    .or. csf_i%stepvector(excitInfo%fullEnd) == 2) then
                    compFlag = .false.
                    return
                end if

                ! in the case of d = 2 at end

                ! in the case that the end is d = 1:
                set_two(excitation, excitInfo%fullEnd)

                clr_one(excitation, excitInfo%firstEnd)

                clr_two(excitation, excitInfo%fullStart)
                set_one(excitation, excitInfo%fullStart)

            end if

        end select

        compFlag = isProperCSF_ilut(excitation, .true.)

        if (.not. compFlag) then
            return
        end if

        ! and then recalculate the matrix element

        call convert_ilut_toNECI(ilut, ilutI)
        call convert_ilut_toNECI(excitation, ilutJ)

        call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), dummy, mat_ele, .true.)

        if (near_zero(mat_ele)) then
            compFlag = .false.
            excitation = 0_n_int

            return
        end if

        call encode_matrix_element(excitation, 0.0_dp, 2)
        call encode_matrix_element(excitation, mat_ele, 1)

    end subroutine perform_crude_excitation


    subroutine pick_random_4ind(csf_i, elec_1, elec_2, orb_1, orb_2, elecs, orbs, pgen)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: elec_1, elec_2, orb_1, orb_2
        integer, intent(out) :: elecs(2), orbs(2)
        real(dp), intent(out) :: pgen

        real(dp) :: r

        elecs = 0
        orbs = 0

        pgen = 1.0_dp

        select case (csf_i%stepvector(elec_1))
        case (1)
            elecs(1) = 2 * elec_1 - 1

        case (2)
            elecs(1) = 2 * elec_1

        case (3)
            r = genrand_real2_dSFMT()

            if (r < 0.5_dp) then
                elecs(1) = 2 * elec_1 - 1
            else
                elecs(1) = 2 * elec_1
            end if

            pgen = 0.5_dp

        end select

        select case (csf_i%stepvector(elec_2))
        case (1)
            elecs(2) = 2 * elec_2 - 1

        case (2)
            elecs(2) = 2 * elec_2

        case (3)
            r = genrand_real2_dSFMT()

            if (r < 0.5_dp) then
                elecs(2) = 2 * elec_2 - 1
            else
                elecs(2) = 2 * elec_2
            end if

            pgen = pgen * 0.5_dp

        end select

        select case (csf_i%stepvector(orb_1))
        case (0)
            if (same_spin(elecs(1), elecs(2))) then
                if (is_beta(elecs(1))) then
                    orbs(1) = 2 * orb_1 - 1
                else
                    orbs(1) = 2 * orb_1
                end if
            else
                r = genrand_real2_dSFMT()

                if (r < 0.5_dp) then
                    orbs(1) = 2 * orb_1 - 1
                else
                    orbs(1) = 2 * orb_1
                end if

                pgen = 0.5_dp * pgen
            end if

        case (1)
            orbs(1) = 2 * orb_1

        case (2)
            orbs(1) = 2 * orb_1 - 1

        end select

        ! i think i can restrict the last one or?..

        if (same_spin(elecs(1), elecs(2))) then
            if (is_beta(elecs(1))) then
                if (csf_i%stepvector(orb_2) == 1) then
                    pgen = 0.0_dp
                    elecs = 0
                    orbs = 0
                    return
                else
                    orbs(2) = 2 * orb_2 - 1
                end if
            else
                if (csf_i%stepvector(orb_2) == 2) then
                    pgen = 0.0_dp
                    elecs = 0
                    orbs = 0
                    return
                else
                    orbs(2) = 2 * orb_2
                end if
            end if
        else
            if (is_beta(orbs(1))) then
                if (csf_i%stepvector(orb_2) == 2) then
                    pgen = 0.0_dp
                    elecs = 0
                    orbs = 0
                    return
                else
                    orbs(2) = 2 * orb_2
                end if
            else
                if (csf_i%stepvector(orb_2) == 1) then
                    pgen = 0.0_dp
                    elecs = 0
                    orbs = 0
                    return
                else
                    orbs(2) = 2 * orb_2 - 1
                end if
            end if
        end if

    end subroutine pick_random_4ind




    subroutine gen_crude_guga_single_1(nI, csf_i, orb_i, cc_i, cum_arr)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i, cc_i
        real(dp), intent(out) :: cum_arr(OrbClassCount(cc_i))

        integer :: nOrb, i, label_index, j, n_id(nEl), id_i, s_orb
        real(dp) :: hel, cum_sum

        nOrb = OrbClassCount(cc_i)
        label_index = SymLabelCounts2(1, cc_i)
        n_id = gtID(nI)
        id_i = gtID(orb_i)

        cum_sum = 0.0_dp

        do i = 1, nOrb
            s_orb = sym_label_list_spat(label_index + i - 1)

            if (s_orb == id_i) then
                cum_arr(i) = cum_sum
                cycle
            end if

            hel = 0.0_dp

            select case (csf_i%stepvector(s_orb))
            case (0)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))

                do j = 1, nEl

                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    ! exclude initial orbital, since this case gets
                    ! contributed already outside of loop over electrons!
                    ! but only spin-orbital or spatial??
                    if (n_id(j) == id_i) cycle
                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            case (2)
                ! no restrictions for 2 -> 1 excitations
                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))
                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)
                hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                do j = 1, nEl
                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            end select

            cum_sum = cum_sum + abs_l1(hel)
            cum_arr(i) = cum_sum

        end do

    end subroutine gen_crude_guga_single_1

    subroutine gen_crude_guga_single_2(nI, csf_i, orb_i, cc_i, cum_arr)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i, cc_i
        real(dp), intent(out) :: cum_arr(OrbClassCount(cc_i))

        integer :: nOrb, i, label_index, j, n_id(nEl), id_i, s_orb
        real(dp) :: cum_sum, hel

        nOrb = OrbClassCount(cc_i)
        label_index = SymLabelCounts2(1, cc_i)
        n_id = gtID(nI)
        id_i = gtID(orb_i)

        cum_sum = 0.0_dp

        do i = 1, nOrb
            s_orb = sym_label_list_spat(label_index + i - 1)

            if (s_orb == id_i) then
                cum_arr(i) = cum_sum
                cycle
            end if

            hel = 0.0_dp

            select case (csf_i%stepvector(s_orb))
            case (0)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))

                do j = 1, nEl

                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    ! exclude initial orbital, since this case gets
                    ! contributed already outside of loop over electrons!
                    ! but only spin-orbital or spatial??
                    if (n_id(j) == id_i) cycle
                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))

                    ! now depending on generator and relation of j to
                    ! st and en -> i know sign or don't

                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            case (1)
                ! no restrictions for 2 -> 1 excitations
                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))
                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)
                hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                do j = 1, nEl
                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            end select
            cum_sum = cum_sum + abs_l1(hel)
            cum_arr(i) = cum_sum
        end do

    end subroutine gen_crude_guga_single_2

    subroutine gen_crude_guga_single_3(nI, csf_i, orb_i, cc_i, cum_arr)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i, cc_i
        real(dp), intent(out) :: cum_arr(OrbClassCount(cc_i))

        integer :: nOrb, i, label_index, j, n_id(nEl), id_i, s_orb
        real(dp) :: cum_sum, hel

        nOrb = OrbClassCount(cc_i)
        label_index = SymLabelCounts2(1, cc_i)
        n_id = gtID(nI)
        id_i = gtID(orb_i)

        cum_sum = 0.0_dp

        do i = 1, nOrb
            s_orb = sym_label_list_spat(label_index + i - 1)

            if (s_orb == id_i) then
                cum_arr(i) = cum_sum
                cycle
            end if

            hel = 0.0_dp

            select case (csf_i%stepvector(s_orb))
            case (0)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))

                hel = hel + abs(get_umat_el(id_i, id_i, s_orb, id_i))

                do j = 1, nEl

                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    ! exclude initial orbital, since this case gets
                    ! contributed already outside of loop over electrons!
                    ! but only spin-orbital or spatial??
                    if (n_id(j) == id_i) cycle
                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            case (1)
                ! no restrictions for 2 -> 1 excitations
                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))
                ! do the loop over all the other electrons

                hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))
                hel = hel + abs(get_umat_el(id_i, id_i, s_orb, id_i))

                do j = 1, nEl
                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            case (2)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb))
                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)

                hel = hel + abs(get_umat_el(id_i, id_i, s_orb, id_i))
                hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                do j = 1, nEl

                    ! todo: finish all contributions later for now only do
                    ! those which are the same for all
                    if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle

                    hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                    hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                end do

            end select

            cum_sum = cum_sum + abs_l1(hel)
            cum_arr(i) = cum_sum

        end do

    end subroutine gen_crude_guga_single_3






end module guga_crude_approx_mod
