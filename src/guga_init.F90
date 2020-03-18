#include "macros.h"
! guga module containing all necessary functionality needed to initialize
! a guga simulation
module guga_init
    ! module use statements
    use SystemData, only: tSPN, tHPHF, lNoSymmetry, STOT, nEl, &
        nBasis, tGUGA, tNoBrillouin, tExactSizeSpace, tUHF, tUEGNewGenerator, &
                          tPickVirtUniform, tGenHelWeighted, tGen_4ind_2, tGen_4ind_weighted, &
                          tGen_4ind_reverse, tGen_sym_guga_ueg, tGen_sym_guga_mol, &
                          tGen_nosym_guga, nSpatOrbs, t_consider_diff_bias, &
                          current_stepvector, currentOcc_ilut, currentOcc_int, &
                          currentB_ilut, currentB_int, current_cum_list, &
                          ref_stepvector, ref_b_vector_int, ref_occ_vector, &
                          ref_b_vector_real, treal, tHUB, t_guga_noreorder, tgen_guga_crude, &
                          t_new_real_space_hubbard, t_heisenberg_model, &
                          t_tJ_model

    use CalcData, only: tUseRealCoeffs, tRealCoeffByExcitLevel, RealCoeffExcitThresh, &
                        t_guga_mat_eles, t_hist_tau_search, tSpinProject, &
                        tReplicaEstimates, tPreCond

    use hist_data, only: tHistSpawn

    use LoggingData, only: tCalcFCIMCPsi, tPrintOrbOcc, tRDMonfly, tExplicitAllRDM

    use bit_rep_data, only: tUseFlags, nifd

    use guga_data, only: init_guga_data_procPtrs, orbitalIndex, t_slow_guga_rdms, &
                         t_fast_guga_rdms, n_excit_info_bits

    use guga_procedure_pointers, only: pickOrbitals_single, pickOrbitals_double, &
                        calc_orbital_pgen_contr, calc_mixed_contr, calc_mixed_start_l2r_contr, &
                        calc_mixed_start_r2l_contr, calc_mixed_end_r2l_contr, calc_mixed_end_l2r_contr, &
                        pick_first_orbital, orb_pgen_contrib_type_2, orb_pgen_contrib_type_3, &
                        calc_off_diag_guga_ref

    use guga_excitations, only: pickOrbs_sym_uniform_ueg_single, pickOrbs_sym_uniform_ueg_double, &
                        pickOrbs_sym_uniform_mol_single, pickOrbs_sym_uniform_mol_double, &
                        pickOrbitals_nosym_single, pickOrbitals_nosym_double, &
                        calc_orbital_pgen_contr_ueg, calc_orbital_pgen_contr_mol, &
                        calc_mixed_contr_sym, calc_mixed_contr_nosym, calc_mixed_start_l2r_contr_nosym, &
                        calc_mixed_start_r2l_contr_nosym, calc_mixed_start_contr_sym,&
                        calc_mixed_x2x_ueg, calc_mixed_end_l2r_contr_nosym, calc_mixed_end_r2l_contr_nosym, &
                        calc_mixed_end_contr_sym, pick_first_orbital_nosym_guga_diff, &
                        pick_first_orbital_nosym_guga_uniform, orb_pgen_contrib_type_2_diff, &
                        orb_pgen_contrib_type_3_diff, orb_pgen_contrib_type_2_uniform, &
                        orb_pgen_contrib_type_3_uniform, temp_step_i, temp_step_j, &
                        temp_delta_b, temp_occ_i, temp_b_real_i, calc_off_diag_guga_ref_direct, &
                        pickOrbs_real_hubbard_single, pickOrbs_real_hubbard_double

    use guga_matrixElements, only: calc_off_diag_guga_ref_list

    use FciMCData, only: pExcit2, pExcit4, pExcit2_same, pExcit3_same, tSearchTau

    use constants, only: dp, int_rdm, n_int

    use ParallelHelper, only: iProcIndex

    use tj_model, only: pick_orbitals_guga_heisenberg, calc_orbital_pgen_contr_heisenberg, &
                        init_get_helement_tj_guga, init_get_helement_heisenberg_guga, &
                        init_guga_tJ_model, init_guga_heisenberg_model, &
                        pick_orbitals_guga_tJ

    use back_spawn, only: setup_virtual_mask

    use guga_rdm, only: init_guga_rdm

    use guga_bitRepOps, only: init_guga_bitrep

    ! variable declaration
    implicit none

    private

    public :: checkInputGUGA, init_guga

contains

    subroutine init_guga_orbital_pickers()
        ! this routine, depending on the input set the orbital pickers
        ! to differentiate between the different excitation generators

        ! now i have to differentiate between the real- and momentum space
        ! hubbard models..
        if (tGen_sym_guga_ueg) then
            if (.not. (treal .or. t_new_real_space_hubbard)) then
                pickOrbitals_single => pickOrbs_sym_uniform_ueg_single
                pickOrbitals_double => pickOrbs_sym_uniform_ueg_double
                calc_mixed_contr => calc_mixed_contr_sym
                calc_orbital_pgen_contr => calc_orbital_pgen_contr_ueg
                calc_mixed_start_r2l_contr => calc_mixed_x2x_ueg
                calc_mixed_start_l2r_contr => calc_mixed_x2x_ueg
                calc_mixed_end_r2l_contr => calc_mixed_x2x_ueg
                calc_mixed_end_l2r_contr => calc_mixed_x2x_ueg
            else
                pickOrbitals_single => pickOrbs_real_hubbard_single
                pickOrbitals_double => pickOrbs_real_hubbard_double
                ! what about the contributions? do i need dummy functions?
            end if

        else if (tGen_sym_guga_mol) then

            pickOrbitals_single => pickOrbs_sym_uniform_mol_single
            pickOrbitals_double => pickOrbs_sym_uniform_mol_double
            calc_orbital_pgen_contr => calc_orbital_pgen_contr_mol
            calc_mixed_start_l2r_contr => calc_mixed_start_contr_sym
            calc_mixed_start_r2l_contr => calc_mixed_start_contr_sym
            calc_mixed_end_l2r_contr => calc_mixed_end_contr_sym
            calc_mixed_end_r2l_contr => calc_mixed_end_contr_sym
            calc_mixed_contr => calc_mixed_contr_sym

        else if (tGen_nosym_guga) then
            pickOrbitals_single => pickOrbitals_nosym_single
            pickOrbitals_double => pickOrbitals_nosym_double
            calc_mixed_contr => calc_mixed_contr_nosym
            calc_mixed_start_l2r_contr => calc_mixed_start_l2r_contr_nosym
            calc_mixed_start_r2l_contr => calc_mixed_start_r2l_contr_nosym
            calc_mixed_end_l2r_contr => calc_mixed_end_l2r_contr_nosym
            calc_mixed_end_r2l_contr => calc_mixed_end_r2l_contr_nosym

            ! for specific probability updates:
            if (t_consider_diff_bias) then
                pick_first_orbital => pick_first_orbital_nosym_guga_diff
                orb_pgen_contrib_type_2 => orb_pgen_contrib_type_2_diff
                orb_pgen_contrib_type_3 => orb_pgen_contrib_type_3_diff
            else
                pick_first_orbital => pick_first_orbital_nosym_guga_uniform
                orb_pgen_contrib_type_2 => orb_pgen_contrib_type_2_uniform
                orb_pgen_contrib_type_3 => orb_pgen_contrib_type_3_uniform
            end if

        else if (t_heisenberg_model) then
!             pickOrbitals_single => pickOrbs_sym_uniform_ueg_single
            pickOrbitals_double => pick_orbitals_guga_heisenberg
            calc_orbital_pgen_contr => calc_orbital_pgen_contr_heisenberg
            calc_mixed_contr => calc_mixed_contr_sym

        else if (t_tJ_model) then
            pickOrbitals_single => pick_orbitals_guga_tJ
            pickOrbitals_double => pick_orbitals_guga_heisenberg
            calc_mixed_contr => calc_mixed_contr_sym
            calc_orbital_pgen_contr => calc_orbital_pgen_contr_heisenberg

        else ! standardly also use nosymmetry version
            pickOrbitals_single => pickOrbitals_nosym_single
            pickOrbitals_double => pickOrbitals_nosym_double
            calc_mixed_contr => calc_mixed_contr_nosym
            calc_mixed_start_l2r_contr => calc_mixed_start_l2r_contr_nosym
            calc_mixed_start_r2l_contr => calc_mixed_start_r2l_contr_nosym
            calc_mixed_end_l2r_contr => calc_mixed_end_l2r_contr_nosym
            calc_mixed_end_r2l_contr => calc_mixed_end_r2l_contr_nosym

            ! for specific probability updates:
            if (t_consider_diff_bias) then
                pick_first_orbital => pick_first_orbital_nosym_guga_diff
                orb_pgen_contrib_type_2 => orb_pgen_contrib_type_2_diff
                orb_pgen_contrib_type_3 => orb_pgen_contrib_type_3_diff
            else
                pick_first_orbital => pick_first_orbital_nosym_guga_uniform
                orb_pgen_contrib_type_2 => orb_pgen_contrib_type_2_uniform
                orb_pgen_contrib_type_3 => orb_pgen_contrib_type_3_uniform
            end if

        end if

    end subroutine init_guga_orbital_pickers

    ! in Fortran no executable code is allowed to be in the module header part
    ! so a initialization subroutine is needed, which has to be called in the
    ! other modules using the guga_data module
    subroutine init_guga()
        integer :: i, ierr
        ! main initialization routine

   ! this routine is called in SysInit() of System_neci.F90
        write(6,*) ' ************ Using the GUGA-CSF implementation **********'
        write(6,*) ' Restricting the total spin of the system, tGUGA : ', tGUGA
        write(6,'(A,I5)') '  Restricting total spin S in units of h/2 to ', STOT
        write(6,*) ' So eg. S = 1 corresponds to one unpaired electron '
        write(6,*) ' not quite sure yet how to deal with extensively used m_s quantum number..'
        write(6,*) ' NOTE: for now, although SPIN-RESTRICT is set off, internally m_s(LMS) '
        write(6,*) ' is set to STOT, to make use of reference determinant creations already implemented'
        write(6,*) ' Since NECI always seems to take the beta orbitals first for open shell or '
        write(6,*) ' spin restricted systems, associate those to positively coupled +h/2 orbitals '
        write(6,*) ' to always ensure a S >= 0 value!'
        write(6,*) ' *********************************************************'

        ! initialize the procedure pointer arrays, needed in the matrix
        ! element calculation
        call init_guga_data_procPtrs()

        ! initialize and point the excitation generator functions to the
        ! correct ones
        call init_guga_orbital_pickers()


        ! also have to set tRealCoeffs true to use it in excitation creation
        ! dont actually need realCoeffs anymore since i changed the accessing
        ! of the ilut lists used for excitation creation
        ! but probably have to set a few other things so it works with
        ! reals
        tUseRealCoeffs = .true.

        tUseFlags = .true.

        ! define global variable of spatial orbitals
        ! do that in a more general setup routine! where nBasis is defined
        ! eg
        ! i have to all this routine again from a point after freezing
        ! where the new number of NBasis is determined already..
        nSpatOrbs = nBasis / 2

        if (allocated(orbitalIndex)) deallocate(orbitalIndex)
        ! but also have to set up the global orbitalIndex list
        allocate(orbitalIndex(nSpatOrbs), stat = ierr)
        orbitalIndex = [ (i, i = 1, nSpatOrbs)]

        ! maybe more to come...
        ! also allocate the currrent_ quantities
        if (allocated(current_stepvector)) deallocate(current_stepvector)
        if (allocated(currentB_ilut))      deallocate(currentB_ilut)
        if (allocated(currentOcc_ilut))    deallocate(currentOcc_ilut)
        if (allocated(currentB_int))       deallocate(currentB_int)
        if (allocated(currentOcc_int))     deallocate(currentOcc_int)

        allocate(current_stepvector(nSpatOrbs), stat = ierr)
        allocate(currentB_ilut(nSpatOrbs), stat = ierr)
        allocate(currentOcc_ilut(nSpatOrbs), stat = ierr)
        allocate(currentB_int(nSpatOrbs), stat = ierr)
        allocate(currentOcc_int(nSpatOrbs), stat = ierr)

        if (allocated(current_cum_list)) deallocate(current_cum_list)

        allocate(current_cum_list(nSpatOrbs), stat = ierr)


        ! also allocate the temporary variables used in the matrix element
        ! calculation and also the similar variables for the reference
        ! determinant!
        if (allocated(temp_step_i))   deallocate(temp_step_i)
        if (allocated(temp_step_j))   deallocate(temp_step_j)
        if (allocated(temp_delta_b))  deallocate(temp_delta_b)
        if (allocated(temp_occ_i))    deallocate(temp_occ_i)
        if (allocated(temp_b_real_i)) deallocate(temp_b_real_i)

        allocate(temp_step_i(nSpatOrbs))
        allocate(temp_step_j(nSpatOrbs))
        allocate(temp_delta_b(nSpatOrbs))
        allocate(temp_occ_i(nSpatOrbs))
        allocate(temp_b_real_i(nSpatOrbs))

        if (allocated(ref_stepvector))    deallocate(ref_stepvector)
        if (allocated(ref_b_vector_int))  deallocate(ref_b_vector_int)
        if (allocated(ref_b_vector_real)) deallocate(ref_b_vector_real)
        if (allocated(ref_occ_vector))    deallocate(ref_occ_vector)

        allocate(ref_stepvector(nSpatOrbs))
        allocate(ref_b_vector_int(nSpatOrbs))
        allocate(ref_b_vector_real(nSpatOrbs))
        allocate(ref_occ_vector(nSpatOrbs))

        ! also initiate the pExcit values here.. or otherwise they dont
        ! get initialized without tau-search option on
        if (tGen_nosym_guga) then
            pExcit4 = (1.0_dp - 1.0_dp / real(nSpatOrbs, dp))
            pExcit2 = 1.0_dp / real(nSpatOrbs - 1, dp)

            root_print "initial pExcit4 set to: ", pExcit4
            root_print "initial pExcit2 set to: ", pExcit2

            if (t_consider_diff_bias) then
                pExcit2_same = 0.9_dp
                pExcit3_same = 0.9_dp
            else
                pExcit2_same = 1.0_dp
                pExcit3_same = 1.0_dp
            end if

            root_print "initial pExcit2_same set to: ", pExcit2_same
            root_print "initial pExcit3_same set to: ", pExcit3_same
        end if

        ! for now (time/iteration comparison) reasons, decide which
        ! reference energy calculation method we use
        if (t_guga_mat_eles) then
            ! use the new "direct" calculation method
            calc_off_diag_guga_ref => calc_off_diag_guga_ref_direct

        else
            ! use the "old" with the projected energy list
            calc_off_diag_guga_ref => calc_off_diag_guga_ref_list

        end if

        ! make checks for the RDM calculation
        if (tRDMonfly) then
            call check_rdm_guga_setup()
            call init_guga_rdm()
        end if

        ! make a unified bit rep initializer:
        call init_guga_bitrep(nifd)

    end subroutine init_guga

    subroutine check_rdm_guga_setup
        character(*), parameter :: this_routine = "check_rdm_guga_setup"


        ! check if the integer types fit for out setup
        if (bit_size(0_n_int) /= bit_size(0_int_rdm)) then
            call stop_all(this_routine, "n_int and int_rdm have different size!")
        end if

        ! we use some bits in the rdm_ind for other information..
        ! check if we still have enough space for all the indices..
        if (nSpatOrbs ** 4 > 2 ** (bit_size(int_rdm) - n_excit_info_bits - 1) - 1) then
            call stop_all(this_routine, "cannot store enough indices in rdm_ind!")
        end if

    end subroutine check_rdm_guga_setup


    subroutine checkInputGUGA()
        ! routine to check if all the input parameters given are consistent
        ! and otherwise stops the excecution
        ! is called inf checkinput() in file readinput.F90
        character(*), parameter :: this_routine = 'checkInputGUGA'


        ! test flags: to be removed after tests pass
        if (tRDMonfly .and. .not. tExplicitAllRDM) then
            if (t_slow_guga_rdms .and. t_fast_guga_rdms) then
                call stop_all(this_routine, "not both slow and fast GUGA RDMs!")
            end if
            if (.not. (t_slow_guga_rdms .or. t_fast_guga_rdms)) then
                call stop_all(this_routine, "neither slow nor fast GUGA rdm implo chosen!")
            end if
        end if

        if (tSPN) then
            call stop_all(this_routine, &
                "GUGA not yet implemented with spin restriction SPIN-RESTRICT!")
        end if

        if (tHPHF) then
            call stop_all(this_routine, &
                "GUGA not compatible with HPHF option!")
        end if

        if (tSpinProject) then
            call stop_all(this_routine, &
                "GUGA not compatible with tSpinProject!")
        end if

        ! with the new UEG/Hubbard implementation of the excitation generator
        ! i need symmetry actually!! or otherwise its wrong
        ! have to somehow find out how to check if k-point symmetry is
        ! provided
        if (tGen_sym_guga_ueg .and. lNoSymmetry .and. .not. treal) then
            call stop_all(this_routine, &
                "UEG/Hubbard implementation of GUGA excitation generator needs symmetry but NOSYMMETRY set! abort!")
        end if

        ! in the real-space do not reorder the orbitals!
        if (treal) t_guga_noreorder = .true.

        if (tExactSizeSpace) then
            call stop_all(this_routine, &
                "calculation of exact Hilbert space size not yet implemented with GUGA!")
        end if

        if (tUEGNewGenerator) then
            call stop_all(this_routine, &
                "wrong input: ueg excitation generator chosen! abort!")
        end if

        if (tPickVirtUniform) then
            call stop_all(this_routine, &
                "wrong input: tPickVirtUniform excitation generator chosen! abort!")
        end if

        if (tGenHelWeighted) then
            call stop_all(this_routine, &
                "wrong input: tGenHelWeighted excitation generator chosen with GUGA! abort!")
        end if

        if (tGen_4ind_2) then
            call stop_all(this_routine, &
                "wrong input: tGen_4ind_2 excitation generator chosen with GUGA! abort!")
        end if

        if (tGen_4ind_weighted) then
            call stop_all(this_routine, &
                "wrong input: tGen_4ind_weighted excitation generator chosen with GUGA! abort!")
        end if

        if (tGen_4ind_reverse) then
            call stop_all(this_routine, &
                "wrong input: tGen_4ind_reverse excitation generator chosen with GUGA! abort!")
        end if

        if (.not. tNoBrillouin) then
            call stop_all(this_routine, &
                "Brillouin theorem not valid for GUGA approach!(I think atleast...)")
        end if

        ! also check if provided input values match:
        ! CONVENTION: STOT in units of h/2!
        if (STOT .gt. nEl) then
            call stop_all(this_routine, &
                "total spin S in units of h/2 cannot be higher than number of electrons!")
        end if

        if (mod(STOT,2) /= mod(nEl,2)) then
            call stop_all(this_routine, &
                "number of electrons nEl and total spin S in units of h/2 must have same parity!")
        end if

        ! maybe more to come...
        ! UHF basis is also not compatible with guga? or not... or atleast
        ! i am not yet implementing it in such a way so it can work
        if (tUHF) then
            call stop_all(this_routine, &
                "GUGA approach and UHF basis not yet (or never?) compatible!")
        end if

        if (tRealCoeffByExcitLevel) then
            if (RealCoeffExcitThresh > 2) then
                call stop_all(this_routine, &
                    "can only determine up to excit level 2 in GUGA for now!")
            end if
        end if

        ! avoid using both the old and the new tau search functionality
        if (tSearchTau .and. t_hist_tau_search) then
            call stop_all(this_routine, &
                "can't use both the old and the new tau search option at the same time!")
        end if


        if (tReplicaEstimates) then
            call stop_all(this_routine, &
                "'replica-estimates' not yet implemented with GUGA")
        end if

        if (tPreCond) then
            call stop_all(this_routine, &
                "'precond' not yet implemented with GUGA. mostly because of communication")
        end if
        ! assert that tUseFlags is set, to be able to encode deltaB values
        ! in the ilut representation for excitation generation
        !if (.not.tUseFlags) then
        !    call stop_all(this_routine, &
        !        "tUseFlags has to be .true. to encode deltaB values in ilut!")
        !end if
        ! cannot do assertion here, since flag is set afterwards, but changed
        ! corresponding code, so flag is set.

    end subroutine checkInputGUGA


end module guga_init
