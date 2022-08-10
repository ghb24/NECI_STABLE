#include "macros.h"
submodule(tau_search_conventional) tau_search_conventional_impls
    use util_mod, only: stop_all

    use SystemData, only: nEl, tHPHF, nBasis, max_ex_level, G1, tKPntSym, t_uniform_excits

    use HPHFRandExcitMod, only: ReturnAlphaOpenDet, CalcPGenHPHF, &
                                CalcNonUniPGen
    use HPHF_integrals, only: hphf_off_diag_helement_norm

    use k_space_hubbard, only: calc_pgen_k_space_hubbard_uniform_transcorr, &
                               calc_pgen_k_space_hubbard_transcorr, calc_pgen_k_space_hubbard

    use GenRandSymExcitNUMod, only: &
        construct_class_counts, init_excit_gen_store, clean_excit_gen_store

    use SymExcitDataMod, only: excit_gen_store_type

    use SymExcit3, only: GenExcitations3

    use sym_general_mod, only: SymAllowedExcit

    use Determinants, only: get_helement

    use bit_reps, only: decode_bit_det

    use bit_rep_data, only: NIfTot

    use fortran_strings, only: str

    use DetBitOps, only: FindBitExcitLevel, TestClosedShellDet, &
                         EncodeBitDet, GetBitExcitation

    implicit none

contains

    module subroutine find_tau_from_refdet_conn()

        ! Routine to find an upper bound to tau, by consideration of the
        ! singles and doubles connected to the reference determinant
        !
        ! Obviously, this make assumptions about the possible range of pgen,
        ! so may actually give a tau that is too SMALL for the latest
        ! excitation generators, which is exciting!

        use neci_intfce
        use SymExcit4, only: GenExcitations4, ExcitGenSessionType
        type(excit_gen_store_type) :: store, store2
        logical :: tAllExcitFound, tParity, tSameFunc, tSwapped, tSign
        character(len=*), parameter :: t_r = "find_tau_from_refdet_conn"
        character(len=*), parameter :: this_routine = "find_tau_from_refdet_conn"
        integer :: ex(2, maxExcit), ex2(2, maxExcit), exflag, iMaxExcit, nStore(6), nExcitMemLen(1)
        integer, allocatable :: Excitgen(:)
        real(dp) :: nAddFac, MagHel, pGen, pGenFac
        HElement_t(dp) :: hel
        integer :: ic, nJ(nel), nJ2(nel), ierr, iExcit, ex_saved(2, maxExcit)
        integer(kind=n_int) :: iLutnJ(0:niftot), iLutnJ2(0:niftot)

        type(ExcitGenSessionType) :: session

        integer(n_int), allocatable :: det_list(:, :)
        integer :: n_excits, i, ex_3(2, 3)

        if (tGUGA) then
            call stop_all(this_routine, "Not implemented for GUGA")
        end if

        if (MaxWalkerBloom.isclose.-1._dp) then
            !No MaxWalkerBloom specified
            !Therefore, assume that we do not want blooms larger than n_add if initiator,
            !or 5 if non-initiator calculation.
            if (tTruncInitiator) then
                nAddFac = InitiatorWalkNo
            else
                nAddFac = 5.0_dp    !Won't allow more than 5 particles at a time
            end if
        else
            nAddFac = real(MaxWalkerBloom, dp) !Won't allow more than MaxWalkerBloom particles to spawn in one event.
        end if

        call assign_value_to_tau(clamp(1000.0_dp, min_tau, max_tau), 'Initial assignment')

        ! NOTE: test if the new real-space implementation works with this
        ! function! maybe i also have to use a specific routine for this !
        ! since it might be necessary in the transcorrelated approach to
        ! the real-space hubbard

        ! bypass everything below for the new k-space hubbard implementation
        if (t_k_space_hubbard) then
            if (tHPHF) then
                call Stop_All(this_routine, &
                              "not yet implemented with HPHF, since gen_all_excits not atapted to it!")
            end if

            call gen_all_excits_k_space_hubbard(ProjEDet(:, 1), n_excits, det_list)

            ! now loop over all of them and determine the worst case H_ij/pgen ratio
            do i = 1, n_excits
                call decode_bit_det(nJ, det_list(:, i))
                ! i have to take the right direction in the case of the
                ! transcorrelated, due to non-hermiticity..
                ic = FindBitExcitlevel(det_list(:, i), ilutRef(:, 1))
                ASSERT(ic == 2 .or. ic == 3)
                if (ic == 2) then
                    call GetBitExcitation(ilutRef(:, 1), det_list(:, i), ex, tParity)
                else if (ic == 3) then
                    call GetBitExcitation(ilutRef(:, 1), det_list(:, i), ex_3, tParity)
                end if

                MagHel = abs(get_helement_lattice(nJ, ProjEDet(:, 1)))
                ! and also get the generation probability
                if (t_trans_corr_2body) then
                    if (t_uniform_excits) then
                        ! i have to setup pDoubles and the other quantities
                        ! before i call this functionality!
                        pgen = calc_pgen_k_space_hubbard_uniform_transcorr(ex_3, ic)
                    else
                        pgen = calc_pgen_k_space_hubbard_transcorr( &
                               ProjEDet(:, 1), ilutRef(:, 1), ex_3, ic)
                    end if
                else
                    pgen = calc_pgen_k_space_hubbard( &
                           ProjEDet(:, 1), ilutRef(:, 1), ex, ic)
                end if

                if (MagHel > EPS) then
                    pGenFac = pgen * nAddFac / MagHel

                    if (tau > pGenFac .and. pGenFac > EPS) then
                        call assign_value_to_tau(pGenFac, this_routine)
                    end if
                end if
            end do

            if (tau > 0.075_dp) then
                call assign_value_to_tau(0.075_dp, this_routine)
                write(stdout, "(A,F8.5,A)") "Small system. Setting initial timestep to be ", Tau, " although this &
                                                &may be inappropriate. Care needed"
            else
                write(stdout, "(A,F18.10)") "From analysis of reference determinant and connections, &
                                         &an upper bound for the timestep is: ", Tau
            end if

            return
        end if

        tAllExcitFound = .false.
        Ex_saved(:, :) = 0
        exflag = 3
        tSameFunc = .false.
        call init_excit_gen_store(store)
        call init_excit_gen_store(store2)
        store%tFilled = .false.
        store2%tFilled = .false.
        CALL construct_class_counts(ProjEDet(:, 1), store%ClassCountOcc, &
                                    store%ClassCountUnocc)
        store%tFilled = .true.
        if (tKPntSym) then
            !TODO: It REALLY needs to be fixed so that we don't need to do this!!
            !Setting up excitation generators that will work with kpoint sampling
            iMaxExcit = 0
            nStore(:) = 0
            CALL GenSymExcitIt2(ProjEDet(:, 1), NEl, G1, nBasis, .TRUE., nExcitMemLen, nJ, iMaxExcit, nStore, exFlag)
            allocate(EXCITGEN(nExcitMemLen(1)), stat=ierr)
            IF (ierr /= 0) CALL Stop_All(t_r, "Problem allocating excitation generator")
            EXCITGEN(:) = 0
            CALL GenSymExcitIt2(ProjEDet(:, 1), NEl, G1, nBasis, .TRUE., EXCITGEN, nJ, iMaxExcit, nStore, exFlag)
        end if

        do while (.not. tAllExcitFound)
            if (tKPntSym) then
                call GenSymExcitIt2(ProjEDet(:, 1), nel, G1, nBasis, .false., EXCITGEN, nJ, iExcit, nStore, exFlag)
                if (nJ(1) == 0) exit
                !Calculate ic, tParity and Ex
                call EncodeBitDet(nJ, iLutnJ)
                Ex(:, :) = 0
                ic = FindBitExcitlevel(iLutnJ, iLutRef(:, 1), 2)
                ex(1, 1) = ic
                call GetExcitation(ProjEDet(:, 1), nJ, Nel, ex, tParity)
            else
                if (tReltvy) then
                    call GenExcitations4(session, ProjEDet(:, 1), nJ, exflag, ex_saved, tParity, tAllExcitFound, .false.)
                else
                    CALL GenExcitations3(ProjEDet(:, 1), iLutRef(:, 1), nJ, exflag, Ex_saved, tParity, tAllExcitFound, .false.)
                end if

                IF (tAllExcitFound) EXIT
                Ex(:, :) = Ex_saved(:, :)
                if (Ex(2, 2) == 0) then
                    ic = 1
                else
                    ic = 2
                end if
                call EncodeBitDet(nJ, iLutnJ)
            end if

            ! Exclude an excitation if it isn't symmetry allowed.
            ! Note that GenExcitations3 is not perfect, especially if there
            ! additional restrictions, such as LzSymmetry.
            if (.not. SymAllowedExcit(ProjEDet(:, 1), nJ, ic, ex)) &
                cycle

            if (tHPHF) then
                if (.not. TestClosedShellDet(iLutnJ)) then
                    CALL ReturnAlphaOpenDet(nJ, nJ2, iLutnJ, iLutnJ2, .true., .true., tSwapped)
                    if (tSwapped) then
                        !Have to recalculate the excitation matrix.
                        ic = FindBitExcitLevel(iLutnJ, iLutRef(:, 1), 2)
                        ex(:, :) = 0
                        if (ic <= max_ex_level) then
                            ex(1, 1) = ic
                            call GetBitExcitation(iLutRef(:, 1), iLutnJ, Ex, tParity)
                        end if
                    end if
                end if
                hel = hphf_off_diag_helement_norm(ProjEDet(:, 1), nJ, iLutRef(:, 1), iLutnJ)
            else
                hel = get_helement(ProjEDet(:, 1), nJ, ic, ex, tParity)
            end if

            MagHel = abs(hel)

            !Find pGen (nI -> nJ)
            if (tHPHF) then
                call CalcPGenHPHF(ProjEDet(:, 1), iLutRef(:, 1), nJ, iLutnJ, ex, store%ClassCountOcc, &
                                  store%ClassCountUnocc, pDoubles, pGen, tSameFunc)
            else
                call CalcNonUnipGen(ProjEDet(:, 1), ilutRef(:, 1), ex, ic, store%ClassCountOcc, store%ClassCountUnocc, pDoubles, pGen)
            end if
            if (tSameFunc) cycle
            if (MagHel > 0.0_dp) then
                pGenFac = pGen * nAddFac / MagHel
                if (Tau > pGenFac .and. pGenFac > EPS) then
                    call assign_value_to_tau(pGenFac, this_routine)
                end if
            end if

            !Find pGen(nJ -> nI)
            CALL construct_class_counts(nJ, store2%ClassCountOcc, &
                                        store2%ClassCountUnocc)
            store2%tFilled = .true.
            if (tHPHF) then
                ic = FindBitExcitLevel(iLutnJ, iLutRef(:, 1), 2)
                ex2(:, :) = 0
                if (ic <= max_ex_level) then
                    ex2(1, 1) = ic

                    call GetBitExcitation(iLutnJ, iLutRef(:, 1), Ex2, tSign)
                end if
                call CalcPGenHPHF(nJ, iLutnJ, ProjEDet(:, 1), iLutRef(:, 1), ex2, store2%ClassCountOcc, &
                                  store2%ClassCountUnocc, pDoubles, pGen, tSameFunc)
            else
                ex2(1, :) = ex(2, :)
                ex2(2, :) = ex(1, :)
                call CalcNonUnipGen(nJ, ilutnJ, ex2, ic, store2%ClassCountOcc, store2%ClassCountUnocc, pDoubles, pGen)
            end if
            if (tSameFunc) cycle
            if (MagHel > 0.0_dp) then
                pGenFac = pGen * nAddFac / MagHel
                if (Tau > pGenFac .and. pGenFac > EPS) then
                    call assign_value_to_tau(pGenFac, this_routine)
                end if
            end if

        end do

        call clean_excit_gen_store(store)
        call clean_excit_gen_store(store2)
        if (tKPntSym) deallocate(EXCITGEN)

        if (tau > 0.075_dp) then
            call assign_value_to_tau(0.075_dp, this_routine)
            write(stdout, "(A,F8.5,A)") "Small system. Setting initial timestep to be ", Tau, " although this &
                                            &may be inappropriate. Care needed"
        else
            write(stdout, "(A,F18.10)") "From analysis of reference determinant and connections, &
                                     &an upper bound for the timestep is: ", Tau
        end if

        ! if (tau < min_tau .or. tau > max_tau) then
        !     call stop_all(this_routine, "The determined tau "//str(tau, 4)//" is smaller than min_tau or larger than max_tau")
        ! end if

    end subroutine find_tau_from_refdet_conn

end submodule
