
# List the files that ought to be compiled

set( SRC_FILES

    adi_references.F90
    adi_initiators.F90
    adi_data.F90
    aliasSampling.F90
    Annihilation.F90
    bit_rep_data.F90
    BitReps.F90
    CalcData.F90
    Calc.F90
    calcrho.F
    CMakeLists.txt
    compiler_issues.txt
    cont_time_fcimc.F90
    cont_time_rates.F90
	core_space_util.F90
    CPMDData.F90
    cpmdinit.F
    cpmdstub.F
    davidson_neci.F90
    davidson_semistoch.F90
    default_sets.F90
    DetBitOps.F90
    DetCalcData.F90
    DetCalc.F90
    DeterminantData.F90
    Determinants.F90
    determ_proj.F90
    detham.F
    direct_ci.F90
    enumerate_excitations.F90
    ErrorAnalysis.F90
    exact_diag.F90
    exact_spectrum.F90
    excit.F
    excit_parity.F90
    fast_determ_hamil.F90
    FciMCData.F90
    fcimc_helper.F90
    fcimc_iter_utilities.F90
    FciMCLogging.F90
    fcimc_output.F90
    FciMCPar.F90
    fcimc_pointed_fns.F90
    fcoul.F
    frsblk.F
    ftlm_neci.F90
	gdata_io.F90
    gen_coul.F
    gen_coul_ueg.F90
    global_utilities.F90
    gndts_blk.F
    gndts.F90
    gndwork.F90
    hash.F90
    hamiltonian_linalg.F90
    hdf5_popsfile.F90
    hdiag.F
    hdiag_from_excit.F90
    HElem.F90
    hfbasis.F
    HFCalc.F90
    hilbert_space_size.F90
    hist_data.F90
    HPHFRandExcit.F90
    hubbard.F
	impurity_models.F90
	index_rhash.F90
    init_coul2D.F
    init_coul.F
    initial_trial_states.F90
    initiator_space_procs.F90
    input.F90
    IntegralsData.F90
    Integrals_neci.F90
    ISKRandExcit.F90
    kp_fciqmc_data.F90
    kp_fciqmc.F90
    kp_fciqmc_init.F90
    kp_fciqmc_procs.F90
    kp_fciqmc_proj_est.F90
    lanczos_general.F90
    lanczos_wrapper.F90
    legacy_data.F90
    lineup.F
    LMat_calc.F90
	LMat_class.F90
	LMat_freeze.F90
    LMat_indexing.F90
    LMat_mod.F90
    load_balance_calcnodes.F90
    LoggingData.F90
    Logging.F90
    Makefile
    matmul.F
    modifiedHub.F90
    molcas_interface.F90
    MolproPlugin.F90
    PluginGuest.cpp
    NatOrbs.F90
    NECICore.F90
    neci_intfce.F90
    OneEInts.F90
    orthogonalise.F90
    Orthonorm.F90
    MPI_wrapper.F90
    Parallel_Calc.F90
    parallel_workarounds.F90
    pchb_excitgen.F90
    perturbations.F90
    Popsfile.F90
    proc_ptrs.F90
    ras_data.F90
    ras.F90
    rdm_data.F90
    rdm_data_utils.F90
    rdm_estimators.F90
    rdm_explicit.F90
    rdm_filling.F90
    rdm_finalising.F90
    rdm_integral_fns.F90
    rdm_nat_orbs.F90
    rdm_reading.F90
    readdf.F90
    readinput.F90
    readint.F90
    read_psi.F
    readri.F90
    replica_data.F90
    replica_estimates.F90
    rhodiag.F
    RhoElement.F90
    RotateOrbsData.F90
    RotateOrbs.F90
    RPA.F90
    scrtransf.F
    searching.F90
    semi_stoch_gen.F90
    semi_stoch_procs.F90
	shared_rhash.F90
    sparse_arrays.F90
    spectral_data.F90
    spectral_lanczos.F90
    SymData.F90
    symexcit2.F90
    symexcit3.F90
    symexcit4.F90
    symexcitData.F90
    symexcit.F
    symexcit_old.F
    sym.F90
    sym_general.F90
    symrandexcit1.F
    symrandexcit2.F90
    symrandexcit3.F90
    symrandexcit_Ex_Mag.F90
    symrandexcit4.F90
    symrandexcit5.F90
    symrandexcit.F
    symrandexcit_ueg.F90
    SystemData.F90
    System_neci.F90
    tau_search.F90
    tc_three_body_data.F90
    tc_three_body_excitgen.F90
    trial_wf_gen.F90
    trial_ht_procs.F90
    UMatCache.F90
    tc_precomputed_excitgen.F90
    unit_test_helper_excitgen.F90
    vasp_init.F90
    vasp_interface.F90
    vasp_neci_interface_stub.F90
    double_occ.F90
    tau_search_hist.F90
    real_space_hubbard.F90
    back_spawn.F90
    back_spawn_excit_gen.F90
    guga_bitRepOps.F90
    guga_data.F90
    guga_excitations.F90
    guga_init.F90
    guga_write_H_matrix.F90
    guga_matrixElements.F90
    guga_procedure_pointers.F90
    guga_tausearch.F90
    guga_testsuite.F90
    guga_types.F90
    lattice_mod.F90
    tJ_model.F90
    k_space_hubbard.F90
    lattice_models_utils.F90
    analyse_wf_symmetry.F90
    guga_rdm.F90
    unit_test_helpers.F90
    real_time_data.F90
    cc_amplitudes.F90
    cepa_shifts.F90
    real_time.F90
    real_time_aux.F90
    real_time_init.F90
    real_time_procs.F90
    real_time_read_input.F90
    verlet_aux.F90
    CAS_distribution_init.F90
    excitation_generators.F90
    exc_gen_class_wrappers.F90
    exc_gen_classes.F90

    lib/backtrace.cpp
    lib/cons_neci.F90
    lib/dSFMT.cpp
    lib/dSFMT_interface.F90
    lib/environment_report.F90
    lib/error_handling_neci.F90
    lib/hdf5_util.F90
    lib/iso_c_hack.F90
    lib/MemoryManagement.F90
    lib/ranlux.F
    lib/signals.F90
    lib/signals_helper.cpp
    lib/soft_exit.F90
    lib/timing.F90
    lib/util_mod_fmt.F90
    lib/fortran_strings.F90
    lib/matrix_util.F90
    lib/display_matrices.F90
    lib/binomial_lookup.F90
)

set( RELAX_WARNINGS_SRCS

    Annihilation.F90
    BitReps.F90
    CalcData.F90
    Calc.F90
    calcrho.F
    CMakeLists.txt
    compiler_issues.txt
    cont_time_fcimc.F90
    cont_time_rates.F90
	core_space_util.F90
    CPMDData.F90
    cpmdinit.F
    cpmdstub.F
    davidson_neci.F90
    davidson_semistoch.F90
    default_sets.F90
    DetBitOps.F90
    DetCalcData.F90
    DetCalc.F90
    DeterminantData.F90
    Determinants.F90
    determ_proj.F90
    detham.F
    direct_ci.F90
    enumerate_excitations.F90
    ErrorAnalysis.F90
    exact_diag.F90
    exact_spectrum.F90
    excit.F
    excit_parity.F90
    fast_determ_hamil.F90
    FciMCData.F90
    fcimc_helper.F90
    fcimc_iter_utilities.F90
    FciMCLogging.F90
    fcimc_output.F90
    FciMCPar.F90
    fcimc_pointed_fns.F90
    fcoul.F
    frsblk.F
    ftlm_neci.F90
	gdata_io.F90
    gen_coul.F
    gen_coul_ueg.F90
    global_utilities.F90
    gndts_blk.F
    gndts.F90
    gndwork.F90
    hash.F90
    hamiltonian_linalg.F90
    hdf5_popsfile.F90
    hdiag.F
    hdiag_from_excit.F90
    HElem.F90
    hfbasis.F
    HFCalc.F90
    hilbert_space_size.F90
    hist_data.F90
    HPHFRandExcit.F90
    hubbard.F
	impurity_models.F90
	index_rhash.F90
    init_coul2D.F
    init_coul.F
    initial_trial_states.F90
    initiator_space_procs.F90
    input.F90
    IntegralsData.F90
    Integrals_neci.F90
    ISKRandExcit.F90
    kp_fciqmc_data.F90
    kp_fciqmc.F90
    kp_fciqmc_init.F90
    kp_fciqmc_procs.F90
    kp_fciqmc_proj_est.F90
    lanczos_general.F90
    lanczos_wrapper.F90
    legacy_data.F90
    lineup.F
    LMat_calc.F90
	LMat_class.F90
	LMat_freeze.F90
    LMat_indexing.F90
    LMat_mod.F90
    load_balance_calcnodes.F90
    LoggingData.F90
    Logging.F90
    Makefile
    matmul.F
    modifiedHub.F90
    molcas_interface.F90
    MolproPlugin.F90
    NatOrbs.F90
    NECICore.F90
    neci_intfce.F90
    OneEInts.F90
    orthogonalise.F90
    Orthonorm.F90
    Parallel_Calc.F90
    parallel_workarounds.F90
    pchb_excitgen.F90
    perturbations.F90
    Popsfile.F90
    proc_ptrs.F90
    ras_data.F90
    ras.F90
    rdm_data.F90
    rdm_data_utils.F90
    rdm_estimators.F90
    rdm_explicit.F90
    rdm_filling.F90
    rdm_finalising.F90
    rdm_integral_fns.F90
    rdm_nat_orbs.F90
    rdm_reading.F90
    readdf.F90
    readinput.F90
    readint.F90
    read_psi.F
    readri.F90
    replica_data.F90
    replica_estimates.F90
    rhodiag.F
    RhoElement.F90
    RotateOrbsData.F90
    RotateOrbs.F90
    RPA.F90
    scrtransf.F
    searching.F90
    semi_stoch_gen.F90
    semi_stoch_procs.F90
	shared_rhash.F90
    sparse_arrays.F90
    spectral_data.F90
    spectral_lanczos.F90
    SymData.F90
    symexcit2.F90
    symexcit3.F90
    symexcit4.F90
    symexcitData.F90
    symexcit.F
    symexcit_old.F
    sym.F90
    sym_general.F90
    symrandexcit1.F
    symrandexcit2.F90
    symrandexcit3.F90
    symrandexcit_Ex_Mag.F90
    symrandexcit4.F90
    symrandexcit5.F90
    symrandexcit.F
    symrandexcit_ueg.F90
    SystemData.F90
    System_neci.F90
    tau_search.F90
    tc_three_body_data.F90
    tc_three_body_excitgen.F90
    trial_wf_gen.F90
    trial_ht_procs.F90
    UMatCache.F90
    tc_precomputed_excitgen.F90
    unit_test_helper_excitgen.F90
    vasp_init.F90
    vasp_interface.F90
    vasp_neci_interface_stub.F90
    double_occ.F90
    tau_search_hist.F90
    real_space_hubbard.F90
    back_spawn.F90
    back_spawn_excit_gen.F90
    lattice_mod.F90
    guga_testsuite.F90
    tJ_model.F90
    k_space_hubbard.F90
    lattice_models_utils.F90
    analyse_wf_symmetry.F90
    guga_rdm.F90
    unit_test_helpers.F90
    real_time_data.F90
    cc_amplitudes.F90
    cepa_shifts.F90
    real_time.F90
    real_time_aux.F90
    real_time_init.F90
    real_time_procs.F90
    real_time_read_input.F90
    verlet_aux.F90
    CAS_distribution_init.F90
)

set( FYPP_SRC_FILES
    growing_buffers.fpp
    fcimc_initialisation.fpp
    HPHFIntegrals.fpp
    sltcnd.fpp
    gasci.fpp
    gasci_disconnected.fpp
    gasci_discarding.fpp
    gasci_general.fpp
    gasci_pchb.fpp
    gasci_supergroup_index.fpp
    gasci_util.fpp
    orb_idx.fpp
    excitation_types.fpp
	shared_array.fpp
    load_balancer.fpp
    global_det_data.fpp
    rdm_general.fpp
    hist.fpp
    lib/util_mod.fpp
    lib/util_mod_epsilon_close.fpp
    sets_mod.fpp
	shared_ragged_array.fpp
)

set( TEMPLATED_SRC_FILES

    BinSearch.F90.template
    Parallel.F90.template
    shared_memory_mpi.F90.template
    scalar_shared_memory_mpi.F90.template
    lib/quicksort.F90.template
    lib/util_mod_byte_size.F90.template
    lib/util_mod_comparisons.F90.template
    lib/util_mod_cpts.F90.template
    lib/util_mod_numerical.F90.template
)


# ---------------------------------------------------------------------------------

# Don't use generator expressions for including libraries, as this causes problems when there are list
# expansions.
set( _conditional_libraries "" )
if (HAVE_FFTW)
    list(APPEND _conditional_libraries ${FFTW_LIBRARIES} )
endif()
if (HAVE_HDF5)
    list(APPEND _conditional_libraries ${HDF5_NECI_LIBRARIES} )
endif()

# ---------------------------------------------------------------------------------

set( _library_link_libs
    "${LAPACK_LIBRARIES}"
    "${LIBRT_LIBRARIES}"
    "${_conditional_libraries}"
    "${MPI_LIBRARIES}"
    "${MPI_Fortran_LIBRARIES}")

if( HAVE_SHARED )
    set( _link_type SHARED )
    set( _exe_link_libs "" )
else()
    set( _link_type STATIC )
    set( _exe_link_libs "${_library_link_libs}" )
endif()

# ---------------------------------------------------------------------------------
# Build the main parts as static libraries
# --> They don't need to be recompiled to be included in unit tests/python/...

foreach(_config ${${PROJECT_NAME}_CONFIGURATIONS})

    message(DEBUG "Adding library for configuration ${_config}")
    message(DEBUG "Library definitions: ${${PROJECT_NAME}_${_config}_DEFINITIONS}")
    neci_add_library(

        TARGET lib${_config}
        TYPE ${_link_type}
        OUTPUT_NAME ${_config}
        SOURCES ${SRC_FILES}
        RELAX_WARNINGS ${RELAX_WARNINGS_SRCS}
        TEMPLATED_SOURCES ${TEMPLATED_SRC_FILES}
        FYPP_SOURCES ${FYPP_SRC_FILES}
        LINKER_LANGUAGE Fortran

        PRIVATE_INCLUDES
            "${MPI_Fortran_INCLUDE_PATH}"
            "${HDF5_NECI_INCLUDE_PATH}"

        DEFINITIONS
            "${NECI_GLOBAL_DEFINES}"
            "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
            "${HDF5_NECI_DEFINITIONS}"

        LIBS "${_library_link_libs}"
	    WARNERR ${${PROJECT_NAME}_${_config}_WARNERR}
    )

    message(DEBUG "Adding executable for configuration ${_config}")
    neci_add_executable(
        TARGET ${_config}
        SOURCES necimain.F90
        LIBS
            lib${_config}
            "${_exe_link_libs}"

        LINKER_LANGUAGE Fortran
    )

endforeach()
