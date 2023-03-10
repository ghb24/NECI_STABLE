#include "macros.h"
submodule (tau_main) tau_main_impls
    use constants, only: n_int, maxexcit

    use FciMCData, only: ProjEDet, ilutRef, pDoubles, pParallel

    use CalcData, only: InitiatorWalkNo, tTruncInitiator, tReadPops, tWalkContGrow

    use SystemData, only: t_k_space_hubbard, t_trans_corr_2body, tReltvy, tGUGA, &
        nOccAlpha, nOccBeta

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard

    use util_mod, only: operator(.isclose.)

    use lattice_mod, only: get_helement_lattice

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

    use fortran_strings, only: str, to_lower

    use DetBitOps, only: FindBitExcitLevel, TestClosedShellDet, &
                         EncodeBitDet, GetBitExcitation

    use neci_intfce, only: GenSymExcitIt2

    use excit_mod, only: GetExcitation

    use SymExcit4, only: GenExcitations4, ExcitGenSessionType

    use tau_search_conventional, only: init_tau_search_conventional, finalize_tau_search_conventional

    use tau_search_hist, only: init_hist_tau_search, finalize_hist_tau_search, t_fill_frequency_hists

    better_implicit_none

contains

    module subroutine init_tau()
        ! And what is the maximum death-component found
        max_death_cpt = 0

        if (tau_start_val == possible_tau_start%refdet_connections) then
            call find_tau_from_refdet_conn()
        end if

        if (tau_search_method == possible_tau_search_methods%CONVENTIONAL) then
            call init_tau_search_conventional()
        else if (tau_search_method == possible_tau_search_methods%HISTOGRAMMING) then
            call init_hist_tau_search()
        else
            ! Add a couple of checks for sanity
            if (nOccAlpha == 0 .or. nOccBeta == 0) then
                pParallel = 1.0_dp
            end if
            if (nOccAlpha == 1 .and. nOccBeta == 1) then
                pParallel = 0.0_dp
            end if
        end if


        write(stdout, *) ">>> Initial tau from source: ", to_lower(tau_start_val%str), " is ", tau, "."
        if (tau_search_method /= possible_tau_search_methods%OFF) then
            write(stdout, *) ">>> Tau-search activated. Using ", to_lower(tau_search_method%str), " algorithm. ", &
                "Stop if ", to_lower(tau_stop_method%str), '.'
        else
            write(stdout, *) ">>> Tau-search off."
        end if

    end subroutine

    module subroutine stop_tau_search(stop_method)
        type(StopMethod_t), intent(in) :: stop_method
        if (tau_search_method == possible_tau_search_methods%HISTOGRAMMING) then
            t_fill_frequency_hists = .false.
        end if
        write(stdout, *)
        write(stdout, *) 'The current tau search method is: ', trim(tau_search_method%str)
        write(stdout, *) 'It is switched off now because of: ', trim(stop_method%str)
        write(stdout, *)
        tau_search_method = possible_tau_search_methods%OFF
    end subroutine

    module subroutine finalize_tau()
        call finalize_tau_main()
        call finalize_tau_search_conventional()
        call finalize_hist_tau_search()
    end subroutine

    module subroutine find_tau_from_refdet_conn()

        !! Routine to find an upper bound to tau, by consideration of the
        !! singles and doubles connected to the reference determinant
        !!
        !! Obviously, this make assumptions about the possible range of pgen,
        !! so may actually give a tau that is too SMALL for the latest
        !! excitation generators, which is exciting!
        character(len=*), parameter :: this_routine = "find_tau_from_refdet_conn"

        if (tGUGA) then
            call stop_all(this_routine, "Not implemented for GUGA")
        else if (t_k_space_hubbard) then
            call hubbard_find_tau_from_refdet_conn()
        else
            call ab_initio_find_tau_from_refdet_conn()
        end if
    end subroutine find_tau_from_refdet_conn

    subroutine ab_initio_find_tau_from_refdet_conn()

        type(excit_gen_store_type) :: store, store2
        logical :: tAllExcitFound, tParity, tSameFunc, tSwapped, tSign
        character(len=*), parameter :: this_routine = "find_tau_from_refdet_conn"
        integer :: ex(2, maxExcit), ex2(2, maxExcit), exflag, iMaxExcit, nStore(6), nExcitMemLen(1)
        integer, allocatable :: Excitgen(:)
        real(dp) :: nAddFac, MagHel, pGen, pGenFac
        HElement_t(dp) :: hel
        integer :: ic, nJ(nel), nJ2(nel), iExcit, ex_saved(2, maxExcit)
        integer(kind=n_int) :: iLutnJ(0:niftot), iLutnJ2(0:niftot)
        real(dp) :: new_tau

        type(ExcitGenSessionType) :: session

        ASSERT(.not. tGUGA)
        ASSERT(.not. t_k_space_hubbard)

        new_tau = huge(new_tau)

        nAddFac = MaxWalkerBloom

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
            allocate(EXCITGEN(nExcitMemLen(1)))
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
                new_tau = clamp(&
                    merge(new_tau, min(pGenFac, new_tau), near_zero(pGenFac)), &
                    min_tau, max_tau)
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
                new_tau = clamp(&
                    merge(new_tau, min(pGenFac, new_tau), near_zero(pGenFac)), &
                    min_tau, max_tau)
            end if
        end do


        call clean_excit_gen_store(store)
        call clean_excit_gen_store(store2)
        if (tKPntSym) deallocate(EXCITGEN)

        if (new_tau > 0.075_dp) then
            new_tau = clamp(0.075_dp, min_tau, max_tau)
            write(stdout, "(A,F8.5,A)") "Small system. Setting initial timestep to be ", Tau, " although this &
                                            &may be inappropriate. Care needed"
        end if
        call assign_value_to_tau(new_tau, this_routine)

    end subroutine ab_initio_find_tau_from_refdet_conn


    subroutine hubbard_find_tau_from_refdet_conn()

        ! Routine to find an upper bound to tau, by consideration of the
        ! singles and doubles connected to the reference determinant
        !
        ! Obviously, this make assumptions about the possible range of pgen,
        ! so may actually give a tau that is too SMALL for the latest
        ! excitation generators, which is exciting!

        character(len=*), parameter :: this_routine = "find_tau_from_refdet_conn"
        integer :: ex(2, maxExcit), ic, nJ(nel), n_excits, i, ex_3(2, 3)
        real(dp) :: nAddFac, MagHel, pGen, pGenFac
        logical :: tParity
        integer(n_int), allocatable :: det_list(:, :)
        real(dp) :: new_tau

        ASSERT(t_k_space_hubbard)
        ASSERT(.not. tGUGA)

        ! NOTE: test if the new real-space implementation works with this
        ! function! maybe i also have to use a specific routine for this !
        ! since it might be necessary in the transcorrelated approach to
        ! the real-space hubbard

        new_tau = huge(new_tau)

        nAddFac = MaxWalkerBloom

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
                new_tau = clamp(&
                    merge(new_tau, min(pGenFac, new_tau), near_zero(pGenFac)), &
                    min_tau, max_tau)
            end if
        end do

        if (new_tau > 0.075_dp) then
            new_tau = clamp(0.075_dp, min_tau, max_tau)
            write(stdout, "(A,F8.5,A)") "Small system. Setting initial timestep to be ", Tau, " although this &
                                            &may be inappropriate. Care needed"
        end if
        call assign_value_to_tau(new_tau, this_routine)
    end subroutine hubbard_find_tau_from_refdet_conn

    subroutine finalize_tau_main()
        !! Reset the values
        tau = 0._dp

        tau_search_method = possible_tau_search_methods%OFF
        if (allocated(input_tau_search_method)) deallocate(input_tau_search_method)
        tau_stop_method = possible_tau_stop_methods%var_shift
        if (allocated(tau_start_val)) deallocate(tau_start_val)

        search_data = TauSearchData_t()
        stop_options = StopOptions_t()

        min_tau = 0._dp
        max_tau = huge(max_tau)
        taufactor = 0._dp
        scale_tau_to_death_triggered = .false.
        t_scale_tau_to_death = .false.
        max_death_cpt = 0._dp
        readpops_but_tau_not_from_popsfile = .false.
    end subroutine



end submodule
