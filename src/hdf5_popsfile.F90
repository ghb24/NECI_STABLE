#include "macros.h"

module hdf5_popsfile

    ! Read and write POPSFILES using hdf5 as the data format
    !
    ! Data structure:
    !
    !   --> Don't assume that only walker-data is going to be relevant
    !   --> Create groups appropriately.
    !   --> Note that we use collective writing, so we write all data from all
    !       nodes (except where explicitly managed), and the HDF library
    !       ensures the writes happen in a sensible way.

        !namelist /POPSHEAD/ Pop64Bit,PopHPHF,PopLz,PopNEl, &
        !            PopCyc,PopNIfFlag,PopNIfTot, &
        !            PopTau,PopiBlockingIter,PopRandomHash&
        !            PopNNodes, PopWalkersOnNodes, &
        !            PopMultiSumNoatHF, PopMultiSumENum, PopBalanceBlocks
    !
    ! A: vcs_ver             - The SHA ID of the git commit
    ! A: compiled_at         - The time of code compilation
    ! A: compiled_on         - The date of code compilation
    ! A: date                - Date/time of calculation stored in popsfile
    ! A: seq_no              - The sequence number of the calculation (i.e.
    !                          how many restarts have there been).
    ! A: config              - The build configuration
    ! A: status              - (opt) indicates if local changes beyond SHA ID
    !
    ! /system/               - Details of the system that is being restarted
    !
    ! /calculation/          - Details of the calculation
    !     /random_hash/      - Random values used in the orbital mapping
    !     /tau_search/       - Values used in timestep optimisation
    !         /gamma_sing/
    !         /gamma_doub/
    !         /gamma_opp/
    !         /gamma_par/
    !         /enough_sing/
    !         /enough_doub/
    !         /enough_opp/
    !         /enough_par/
    !         /cnt_sing/
    !         /cnt_doub/
    !         /cnt_opp/
    !         /cnt_par/
    !         /max_death_cpt/
    !
    !         /psingles/     - And the values which have been optimised
    !         /pdoubles/
    !         /pparallel/
    !         /tau/
    !     /accumulators/     - Accumulated (output) data
    !         /sum_no_ref/
    !         /sum_enum/
    !     /completed_iters/  - How many iterations have already been completed
    !     /tot_imag_time/    - Total amount of imaginary time completed
    !     /shift/            - The diagshift value (absolute, invarient to a
    !                          change of reference)
    ! 
    ! /wavefunction/         - Details of a determinental Hilbert space
    !     A: width           - Width of the bit-rep in 64-bit integers
    !     A: num_dets        - Number of determinants in the file
    !     A: lenof_sign      - Number of elements in the sign data
    !     A: norm_sqr        - sum(sgn**2) for each sign element
    !     A: num_parts       - sum(abs(sgn)) for each sign element
    !
    !     /ilut/             - The bit representations of the determinants
    !     /sgns/             - The occupation of the determinants

    use ParallelHelper
    use Parallel_neci
    use constants
    use hdf5_util
    use util_mod
#ifdef __USE_HDF
    use hdf5
#endif
    implicit none
    private

    ! Constants for naming various sections
    character(*), parameter :: &

            nm_vcs_ver = 'vcs_ver', &
            nm_comp_time = 'compiled_at', &
            nm_comp_date = 'compiled_on', &
            nm_date = 'date', &
            nm_seq_no = 'seq_no', &
            nm_comp_config = 'config', &
            nm_comp_status = 'status', &

            nm_calc_grp = 'calculation', &
            nm_random_hash = 'random_hash', &
            nm_iters = 'completed_iters', &
            nm_tot_imag = 'tot_imag_time', &
            nm_shift = 'shift', &

            nm_tau_grp = 'tau_search', &
            nm_gam_sing = 'gamma_sing', &
            nm_gam_doub = 'gamma_doub', &
            nm_gam_opp = 'gamma_opp', &
            nm_gam_par = 'gamma_par', &
            nm_en_sing = 'enough_sing', &
            nm_en_doub = 'enough_doub', &
            nm_en_opp = 'enough_opp', &
            nm_en_par = 'enough_par', &
            nm_cnt_sing = 'cnt_sing', &
            nm_cnt_doub = 'cnt_doub', &
            nm_cnt_opp = 'cnt_opp', &
            nm_cnt_par = 'cnt_par', &
            nm_max_death = 'max_death', &
            nm_psingles = 'psingles', &
            nm_pdoubles = 'pdoubles', &
            nm_pparallel = 'pparallel', &
            nm_tau = 'tau', &
            ! [W.D.]: 
            ! can i just add another entry without breaking anything?
            nm_hist_tau = 'hist_tau_search', &

            nm_acc_grp = 'accumulators', &
            nm_sum_no_ref = 'sum_no_ref', &
            nm_sum_enum = 'sum_enum', &
            
            nm_wfn_grp = 'wavefunction', &
            nm_rep_width = 'width', &
            nm_sgn_len = 'lenof_sign', &
            nm_num_dets = 'num_dets', &
            nm_ilut = 'ilut', &
            nm_sgns = 'sgns', &
            nm_norm_sqr = 'norm_sqr', &
            nm_num_parts = 'num_parts'

#ifdef __USE_HDF
    ! hsize_t is only defined in the hdf5 library
    integer(hsize_t), dimension(:,:), allocatable :: receivebuff
#else
    integer, dimension(:,:), allocatable :: receivebuff
#endif
    integer:: receivebuff_tag

    public :: write_popsfile_hdf5, read_popsfile_hdf5
    public :: add_pops_norm_contrib

contains

    subroutine write_popsfile_hdf5()

        use CalcData, only: iPopsFileNoWrite
        use LoggingData, only: tIncrementPops

        ! TODO:
        ! 1) Deal with multiple filenames
        ! 2) Deal with build configurations without HDF5
        ! 3) Deal with HDF5 build configurations without MPIO
        ! 4) Should we in some way make incrementpops default?

        character(*), parameter :: t_r = 'write_popsfile_hdf5'
#ifdef __USE_HDF
        integer(hid_t) :: plist_id, file_id, err
        character(255) :: filename

        ! Get a unique filename for this popsfile. This needs to be done on
        ! the head node to avoid collisions.
        if (iProcIndex == 0) &
            call get_unique_filename('popsfile', tIncrementPops, .true., &
                                     iPopsFileNoWrite, filename, ext='.h5')
        call MPIBCast(filename)

        write(6,*)
        write(6,*) "============== Writing HDF5 popsfile =============="
        write(6,*) "File name: ", trim(filename)

        ! Initialise the hdf5 fortran interface
        call h5open_f(err)

        ! Set up a property list to ensure file handling across all nodes.
        ! TODO: Check if we should be using a more specific communicator
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
        call h5pset_fapl_mpio_f(plist_id, CommGlobal, mpiInfoNull, err)

        ! TODO: Do sensible file handling here...
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, &
                         access_prp=plist_id)
        call h5pclose_f(plist_id, err)
        write(6,*) "writing metadata"
        call write_metadata(file_id)
        write(6,*) "writing calc_data"
        call write_calc_data(file_id)

        call MPIBarrier(err)
        write(6,*) "writing walkers"
        call write_walkers(file_id)

        call MPIBarrier(err)
        write(6,*) "closing popsfile"
        ! And we are done!
        call h5fclose_f(file_id, err)
        call h5close_f(err)

        call MPIBarrier(err)
        write(6,*) "popsfile write successful"
#else
        call stop_all(t_r, 'HDF5 support not enabled at compile time')
#endif

    end subroutine write_popsfile_hdf5


    function read_popsfile_hdf5(dets) result(CurrWalkers)

        use CalcData, only: iPopsFileNoRead
        use LoggingData, only: tIncrementPops

        ! Read a popsfile in, prior to running a new calculation
        ! TODO: Integrate with CheckPopsParams

        ! n.b. This reads into the specified array, to allow use of popsfiles
        !      for initialising perturbations, etc.
        integer(n_int), intent(out) :: dets(:, :)
        integer(int64) :: CurrWalkers
        character(*), parameter :: t_r = 'write_popsfile_hdf5'
#ifdef __USE_HDF
        integer(hid_t) :: err, file_id, plist_id
        integer :: tmp
        character(255) :: filename

        ! Get the name for the popsfile to read in
        if (iProcIndex == 0) &
            call get_unique_filename('popsfile', tIncrementPops, .false., &
                                     iPopsFileNoRead, filename, ext='.h5')
        call MPIBcast(filename)

        write(6,*)
        write(6,*) "========== Reading in from HDF5 popsfile =========="
        write(6,*) 'File name: ', trim(filename)

        ! Initialise the hdf5 fortran interface
        call h5open_f(err)

        ! Set up a property list to ensure file handling across all nodes.
        ! TODO: Check if we should be using a more specific communicator
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
        call h5pset_fapl_mpio_f(plist_id, CommGlobal, mpiInfoNull, err)

        ! Open the popsfile
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err, &
                       access_prp=plist_id)

        call read_metadata(file_id)
        call read_walkers(file_id, dets, CurrWalkers)
        call read_calc_data(file_id)

        ! And we are done
        call h5pclose_f(plist_id, err)
        call h5fclose_f(file_id, err)
        call h5close_f(err)

        call neci_flush(6)
        call MPIBarrier(tmp)
#else
        CurrWalkers = 0
        call stop_all(t_r, 'HDF5 support not enabled at compile time')
#endif

    end function


#ifdef __USE_HDF
    subroutine write_metadata(parent)

        use CalcData, only: calc_seq_no

        ! Output macroscopic metadata applicable to a restart file, which may
        ! be used for establishing providence of calculations, etc.

        integer(hid_t), intent(in) :: parent
        character(19) :: date_str
        integer :: date_values(8)

        ! TODO: Run by
        call write_string_attribute(parent, nm_vcs_ver, _VCS_VER)
        call write_string_attribute(parent, nm_comp_date, __DATE__)
        call write_string_attribute(parent, nm_comp_time, __TIME__)
        call write_string_attribute(parent, nm_comp_config, _CONFIG)
#ifdef _WORKING_DIR_CHANGES
        call write_string_attribute(parent, nm_comp_status, &
                                    "Working directory contains local changes")
#endif

        ! How many calculations have been run to get to this popsfile pt?
        call write_int32_attribute(parent, nm_seq_no, int(calc_seq_no, int32))

        ! When are we running this?
        call date_and_time(values=date_values)

        write(date_str, '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') &
            date_values(1:3), date_values(5:7)
        call write_string_attribute(parent, nm_date, date_str)

    end subroutine

    subroutine read_metadata(parent)

        use CalcData, only: calc_seq_no

        ! Read in the macroscopic metadata applicable to the restart file.

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: err, attribute

        logical :: exists
        character(100) :: str_buf

        write(6,*) 'Previous calculation'

        call read_string_attribute(parent, nm_date, str_buf, exists)
        if (exists) write(6,*) 'Date: ', trim(str_buf)
        call read_int32_attribute(parent, nm_seq_no,  calc_seq_no, &
                                  default=1_int32)
        write(6,*) 'Sequence no.:', calc_seq_no

        ! Output nice details for usability
        call read_string_attribute(parent, nm_vcs_ver, str_buf, exists)
        if (exists) write(6,*) 'VCS ver: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_config, str_buf, exists)
        if (exists) write(6,*) 'Build config: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_status, str_buf, exists)
        if (exists) write(6,*) 'Build status: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_date, str_buf, exists)
        if (exists) write(6,*) 'Build date: ', trim(str_buf)
        call read_string_attribute(parent, nm_comp_time, str_buf, exists)
        if (exists) write(6,*) 'Build time: ', trim(str_buf)

        ! Update values for the new calculation
        calc_seq_no = calc_seq_no + 1

    end subroutine

    subroutine write_calc_data(parent)

        use load_balance_calcnodes, only: RandomOrbIndex
        use FciMCData, only: Iter, PreviousCycles, TotImagTime, Hii
        use CalcData, only: DiagSft

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: calc_grp, err

        ! Firstly create the group for storing calculation-related data
        call h5gcreate_f(parent, nm_calc_grp, calc_grp, err)

        ! Write out the random orbital mapping index
        call write_int64_1d_dataset(calc_grp, nm_random_hash, RandomOrbIndex)

        call MPIBcast(PreviousCycles)
        call write_int64_scalar(calc_grp, nm_iters, iter + PreviousCycles)
        call write_dp_scalar(calc_grp, nm_tot_imag, TotImagTime)
        call write_dp_1d_dataset(calc_grp, nm_shift, DiagSft + Hii)

        ! Output the values used for tau optimisation. Only output non-zero
        ! (i.e. used) values.
        call write_tau_opt(calc_grp)

        ! Output accumulator data
        call write_accumulator_data(calc_grp)

        ! Clear stuff up
        call h5gclose_f(calc_grp, err)

    end subroutine

    subroutine write_tau_opt(parent)
    
        use tau_search, only: gamma_sing, gamma_doub, gamma_opp, gamma_par, &
                              enough_sing, enough_doub, enough_opp, &
                              enough_par, cnt_sing, cnt_doub, cnt_opp, &
                              cnt_par, max_death_cpt
        use FciMCData, only: pSingles, pDoubles, pParallel
        use CalcData, only: tau, t_hist_tau_search_option, t_previous_hist_tau

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: tau_grp, err

        real(dp) :: max_gam_sing, max_gam_doub, max_gam_opp, max_gam_par
        real(dp) :: max_max_death_cpt
        logical :: all_en_sing, all_en_doub, all_en_opp, all_en_par
        integer :: max_cnt_sing, max_cnt_doub, max_cnt_opp, max_cnt_par

        real(dp) :: all_pdoub, all_psing, all_ppar, all_tau

        ! Create the group
        call h5gcreate_f(parent, nm_tau_grp, tau_grp, err)

        ! We want to use the maximised values across all the processors
        ! (there is nothing ensuring that all the processors are adjusted to
        ! the same values...)
        call MPIAllReduce(gamma_sing, MPI_MAX, max_gam_sing)
        call MPIAllReduce(gamma_doub, MPI_MAX, max_gam_doub)
        call MPIAllReduce(gamma_opp, MPI_MAX, max_gam_opp)
        call MPIAllReduce(gamma_par, MPI_MAX, max_gam_par)
        call MPIAllReduce(max_death_cpt, MPI_MAX, max_max_death_cpt)
        call MPIAllLORLogical(enough_sing, all_en_sing)
        call MPIAllLORLogical(enough_doub, all_en_doub)
        call MPIAllLORLogical(enough_opp, all_en_opp)
        call MPIAllLORLogical(enough_par, all_en_par)
        call MPIAllReduce(cnt_sing, MPI_MAX, max_cnt_sing)
        call MPIAllReduce(cnt_doub, MPI_MAX, max_cnt_doub)
        call MPIAllReduce(cnt_opp, MPI_MAX, max_cnt_opp)
        call MPIAllReduce(cnt_par, MPI_MAX, max_cnt_par)

        if (max_gam_sing /= 0) &
            call write_dp_scalar(tau_grp, nm_gam_sing, max_gam_sing)
        if (max_gam_doub /= 0) &
            call write_dp_scalar(tau_grp, nm_gam_doub, max_gam_doub)
        if (max_gam_opp /= 0) &
            call write_dp_scalar(tau_grp, nm_gam_opp, max_gam_opp)
        if (max_gam_par /= 0) &
            call write_dp_scalar(tau_grp, nm_gam_par, max_gam_par)
        if (max_max_death_cpt /= 0) &
            call write_dp_scalar(tau_grp, nm_max_death, max_max_death_cpt)
        if (all_en_sing) &
            call write_log_scalar(tau_grp, nm_en_sing, all_en_sing)
        if (all_en_doub) &
            call write_log_scalar(tau_grp, nm_en_doub, all_en_doub)
        if (all_en_opp) &
            call write_log_scalar(tau_grp, nm_en_opp, all_en_opp)
        if (all_en_par) &
            call write_log_scalar(tau_grp, nm_en_par, all_en_par)
        if (max_cnt_sing /= 0) &
            call write_int64_scalar(tau_grp, nm_cnt_sing, max_cnt_sing)
        if (max_cnt_doub /= 0) &
            call write_int64_scalar(tau_grp, nm_cnt_doub, max_cnt_doub)
        if (max_cnt_opp /= 0) &
            call write_int64_scalar(tau_grp, nm_cnt_opp, max_cnt_opp)
        if (max_cnt_par /= 0) &
            call write_int64_scalar(tau_grp, nm_cnt_par, max_cnt_par)

        ! Use the probability values from the head node
        all_psing = pSingles; all_pdoub = pDoubles; all_ppar = pParallel
        all_tau = tau
        call MPIBcast(all_psing)
        call MPIBcast(all_pdoub)
        call MPIBcast(all_ppar)
        call MPIBcast(all_tau)

        call write_dp_scalar(tau_grp, nm_psingles, all_psing)
        call write_dp_scalar(tau_grp, nm_pdoubles, all_pdoub)
        call write_dp_scalar(tau_grp, nm_pparallel, all_ppar)
        call write_dp_scalar(tau_grp, nm_tau, all_tau)

        ! [W.D.]:
        ! for the new hist-tau search i essentially only need to indicat 
        ! that a histogramming tau-search was used: 
        if (t_hist_tau_search_option .or. t_previous_hist_tau) then 
            call write_log_scalar(tau_grp, nm_hist_tau, .true.)
        end if

        ! Clear up
        call h5gclose_f(tau_grp, err)

    end subroutine write_tau_opt


    subroutine write_accumulator_data(parent)

        use FciMCData, only: AllSumNoatHF, AllSumENum

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: acc_grp, err

        ! Create group
        call h5gcreate_f(parent, nm_acc_grp, acc_grp, err)

        ! Write the energy accumulator values
        ! (n.b. ensure values on all procs)
        call MPIBcast(AllSumENum)
        call MPIBcast(AllSumNoatHF)
#ifdef __CMPLX
        call write_cplx_1d_dataset(acc_grp, nm_sum_enum, AllSumENum)
#else
        call write_dp_1d_dataset(acc_grp, nm_sum_enum, AllSumENum)
#endif
        call write_dp_1d_dataset(acc_grp, nm_sum_no_ref, AllSumNoatHF)

        ! Clear up
        call h5gclose_f(acc_grp, err)

    end subroutine

    subroutine read_calc_data(parent)

        use load_balance_calcnodes, only: RandomOrbIndex
        use FciMCData, only: PreviousCycles, Hii, TotImagTime, tSearchTauOption, &
                             tSearchTau, pSingles, pDoubles, pParallel
        use CalcData, only: DiagSft, tWalkContGrow, tau, t_hist_tau_search, &
                            hdf5_diagsft
        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: grp_id, err
        logical :: exists

        call h5gopen_f(parent, nm_calc_grp, grp_id, err)

        ! Read out the random orbital mapping index
        call read_int64_1d_dataset(grp_id, nm_random_hash, RandomOrbIndex, &
                                   required=.true.)

        ! Previous iteration data.
        call read_int64_scalar(grp_id, nm_iters, PreviousCycles, &
                               default=0_int64, exists=exists)
        if (exists) &
            write(6,*) 'Completed iterations: ', PreviousCycles

        call read_dp_scalar(grp_id, nm_tot_imag, TotImagTime, default=0.0_dp, &
                            exists=exists)
        if (exists) &
            write(6,*) 'Resuming calculation after ', TotImagTime, ' a.u.'

        ! Read in the diagsft. Note that it uses the absolute value (to be
        ! independent of choice of reference), so we must subtract out the
        ! reference
        ! TODO: Do scale up from 1 --> 2 runs for RDMs
        if (.not. tWalkContGrow) then
            call read_dp_1d_dataset(grp_id, nm_shift, DiagSft, required=.true.)
            DiagSft = DiagSft - Hii
            tSinglePartPhase = (abs(DiagSft(1)) < 1.0e-6)
            write(6,*) 'Initial shift: ', DiagSft
        else
            tSinglePartPhase = .true.
            
            ! i still want to capture the diagshift in a temporary file 
            ! atleast 
            call read_dp_1d_dataset(grp_id, nm_shift, hdf5_diagsft, required=.true.)
            hdf5_diagsft = hdf5_diagsft - Hii

        end if

        ! [W.D.]:
        ! i think i also want to read in pSingles etc. even if we do not 
        ! tau-search anymore in a restarted run..
        ! but i guess i have to be careful to set the appropriate 
        ! default, if no tau-search was used and then is restarted..
        ! well, even if the tau-search is not turned on, the 
        ! values are written anyway.. so i can also read them in, but 
        ! not use the read-in tau, but the one specified in the input! 
        ! except the hist_tau was used the we want to use the 
        ! one in the popsfile all the time
!         if (tSearchTauOption) then
        call read_tau_opt(grp_id)
!         else
!             write(6,*) 'Skipping tau optimisation data as tau optimisation is &
!                        &disabled'
!         end if
        ! and also output some info: 
        write(6,*) "read-in tau optimization data: "
        write(6,*) "time-step: ", tau 
        write(6,*) "pSingles: ", pSingles
        write(6,*) "pDoubles: ", pDoubles
        write(6,*) "pParallel: ", pParallel
        if (tSearchTau .or. t_hist_tau_search) then
            write(6,*) "continuing tau-search!"
        else
            write(6,*) "Do not continue tau-search!"
        end if

        call read_accumulator_data(grp_id)

        ! TODO: Read nbasis, nel, ms2, etc.
        !       --> Check that these values haven't changed from the
        !           previous run

        call h5gclose_f(grp_id, err)

    end subroutine

    subroutine read_tau_opt(parent)

        use tau_search, only: gamma_sing, gamma_doub, gamma_opp, gamma_par, &
                              enough_sing, enough_doub, enough_opp, &
                              enough_par, cnt_sing, cnt_doub, cnt_opp, &
                              cnt_par, max_death_cpt, update_tau
        use FciMCData, only: pSingles, pDoubles, pParallel, tSearchTau, &
                             tSearchTauOption 
        use CalcData, only: tau, t_previous_hist_tau, t_restart_hist_tau, &
                            t_hist_tau_search, t_hist_tau_search_option, &
                            t_fill_frequency_hists
        use LoggingData, only: t_print_frq_histograms
        use tau_search_hist, only: deallocate_histograms

        ! Read accumulator values for the timestep optimisation
        ! TODO: Add an option to reset these values...

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: grp_id, err
        logical :: ppar_set, tau_set, hist_tau, temp_previous

        real(dp) :: temp_tau

        call h5gopen_f(parent, nm_tau_grp, grp_id, err)

        ! These are all optional things to have in the popsfile. If they don't
        ! exist, then they will be left unchanged.
        call read_dp_scalar(grp_id, nm_gam_sing, gamma_sing)
        call read_dp_scalar(grp_id, nm_gam_doub, gamma_doub)
        call read_dp_scalar(grp_id, nm_gam_opp, gamma_opp)
        call read_dp_scalar(grp_id, nm_gam_par, gamma_par)
        call read_dp_scalar(grp_id, nm_max_death, max_death_cpt)
        call read_log_scalar(grp_id, nm_en_sing, enough_sing)
        call read_log_scalar(grp_id, nm_en_doub, enough_doub)
        call read_log_scalar(grp_id, nm_en_opp, enough_opp)
        call read_log_scalar(grp_id, nm_en_par, enough_par)
        call read_int64_scalar(grp_id, nm_cnt_sing, cnt_sing)
        call read_int64_scalar(grp_id, nm_cnt_doub, cnt_doub)
        call read_int64_scalar(grp_id, nm_cnt_opp, cnt_opp)
        call read_int64_scalar(grp_id, nm_cnt_par, cnt_par)

        call read_dp_scalar(grp_id, nm_psingles, psingles)
        call read_dp_scalar(grp_id, nm_pdoubles, pdoubles)
        call read_dp_scalar(grp_id, nm_pparallel, pparallel, exists=ppar_set)
        ! here i want to make the distinction if we want to tau-search 
        ! or not
        call read_dp_scalar(grp_id, nm_tau, temp_tau, exists=tau_set)

        call read_log_scalar(grp_id, nm_hist_tau, temp_previous, &
            exists = hist_tau)

        call h5gclose_f(grp_id, err)

        if (tSearchTauOption .and. tau_set) then
           tau = temp_tau 
        end if

        ! also set if previous hist-tau
        if (hist_tau) then
            tau = temp_tau
            ! and turn off if i dont want to force restart! 
            if (.not. t_restart_hist_tau) then
                t_previous_hist_tau = temp_previous

                if (t_previous_hist_tau) then
                    tSearchTau = .false.
                    tSearchTauOption = .false.

                    if (t_hist_tau_search) then
                        call deallocate_histograms()
                        t_hist_tau_search = .false.
                        t_fill_frequency_hists = .false.

                        t_hist_tau_search_option = .true.
                        t_print_frq_histograms = .false.
                    end if
                end if
            end if
        end if

        ! if tau is 0, because no input provided, use the one here too
        if (tau < EPS .and. (.not. temp_tau < EPS)) then 
            tau = temp_tau
        end if

        ! Deal with a previous bug, that leads to popsfiles existing with all
        ! the optimising parameters excluding tau, such that the first
        ! iteration results in chaos
        ! [W.D]: this if should suffice or?
        if (.not. hist_tau) then
            t_previous_hist_tau = .false.
            if (ppar_set .and. .not. tau_set) &
                call update_tau()
        end if

    end subroutine

    subroutine read_accumulator_data(parent)

        use FciMCData, only: AllSumNoatHF, AllSumENum

        integer(hid_t), intent(in) :: parent
        integer(hid_t) :: grp_id, err

        call h5gopen_f(parent, nm_acc_grp, grp_id, err)
        call read_dp_1d_dataset(grp_id, nm_sum_no_ref, AllSumNoatHF, &
                                required=.true.)
#ifdef __CMPLX
        call read_cplx_1d_dataset(grp_id, nm_sum_enum, AllSumENum, &
                                  required=.true.)
#else
        call read_dp_1d_dataset(grp_id, nm_sum_enum, AllSumENum, &
                                required=.true.)
#endif

        call h5gclose_f(grp_id, err)
        
    end subroutine

    subroutine write_walkers(parent)

        use iso_c_hack
        use bit_rep_data, only: NIfD, NIfTot, NOffSgn
        use FciMCData, only: AllTotWalkers, CurrentDets, MaxWalkersPart, &
                             TotWalkers
        use CalcData, only: tUseRealCoeffs

        ! Output the wavefunction information to the relevant groups in the
        ! wavefunction.

        integer(hid_t), intent(in) :: parent
        type(c_ptr) :: cptr
        integer(int32), pointer :: ptr(:)
        integer(int32) :: boop

        character(*), parameter :: t_r = 'write_walkers'

        integer(hid_t) :: wfn_grp_id, dataspace, dataset, err, memspace
        integer(hid_t) :: plist_id

        integer(hsize_t) :: counts(0:nProcessors-1)
        integer(hsize_t) :: all_count

        integer(int32) :: bit_rep_width
        integer(hsize_t) :: mem_offset(2), write_offset(2)
        integer(hsize_t) :: dims(2), hyperdims(2)
        real(dp) :: all_parts(lenof_sign), all_norm_sqr(lenof_sign)
        integer(hsize_t) :: block_size, block_start, block_end
        integer(hsize_t), dimension(:,:), allocatable :: temp_dets
        integer :: ierr

        ! TODO: Add a (slower) fallback routine for weird cases, odd HDF libs

        ! Firstly create the group for storing wavefunction info
        call h5gcreate_f(parent, nm_wfn_grp, wfn_grp_id, err)

        ! TODO: Refactor these chunks into their own little subroutines.
        ! We fix the format of the binary file. Thus if we are on a 32-bit
        ! build, we need to convert the data into 64-bit compatibile chunks.
        if (build_64bit) then
            bit_rep_width = NIfD + 1
        else
            bit_rep_width = 2 * (NIfD + 1)
            call stop_all(t_r, "Needs manual, careful, testing")
        end if

        ! How many occuiped determinants are there on each of the processors
        call MPIAllGather(TotWalkers, counts, ierr)
        all_count = sum(counts)
        write_offset = [0_hsize_t, sum(counts(0:iProcIndex-1))]

        ! Output the bit-representation data
        call write_int32_attribute(wfn_grp_id, nm_rep_width, bit_rep_width)
        call write_int32_attribute(wfn_grp_id, nm_sgn_len, &
                                   int(lenof_sign, int32))
        call write_int64_attribute(wfn_grp_id, nm_num_dets, all_count)
        ! TODO: Check these values. May need to sum them explicitly

        ! Accumulated values only valid on head node. collate_iter_data
        ! has not yet run.
        all_parts = AllTotParts
        call MPISumAll(norm_psi_squared, all_norm_sqr)
        call MPIBcast(all_parts)
        call write_dp_1d_attribute(wfn_grp_id, nm_norm_sqr, all_norm_sqr)
        call write_dp_1d_attribute(wfn_grp_id, nm_num_parts, all_parts)

        !we do an explicitly buffered write to avoid performance problems with
        !complicated hyperslabs + collective buffering
        ! Write out the determinant bit-representations
        call write_2d_multi_arr_chunk_buff( &
                wfn_grp_id, nm_ilut, H5T_NATIVE_INTEGER_8, &
                CurrentDets, arr_2d_dims(CurrentDets), &
                [int(nifd+1, hsize_t), int(TotWalkers, hsize_t)], & ! dims
                [0_hsize_t, 0_hsize_t], & ! offset
                [int(nifd+1, hsize_t), all_count], & ! all dims
                [0_hsize_t, sum(counts(0:iProcIndex-1))] & ! output offset
        )

        ! Write out the sign values on each of the processors
!        if (.not. tUseRealCoeffs) &
!            call stop_all(t_r, "This could go badly...")

        call write_2d_multi_arr_chunk_buff( &
                wfn_grp_id, nm_sgns, H5T_NATIVE_REAL_8, &
                CurrentDets, arr_2d_dims(CurrentDets), &
                [int(lenof_sign, hsize_t), int(TotWalkers, hsize_t)], & ! dims
                [int(nOffSgn, hsize_t), 0_hsize_t], & ! offset
                [int(lenof_sign, hsize_t), all_count], & ! all dims
                [0_hsize_t, sum(counts(0:iProcIndex-1))] & ! output offset
        )

        ! And we are done
        call h5gclose_f(wfn_grp_id, err)

    end subroutine write_walkers

    subroutine read_walkers(parent, dets, CurrWalkers)

        use bit_rep_data, only: NIfD, NIfTot
        use CalcData, only: pops_norm
        use FciMCData, only: InitialSpawnedSlots

        ! This is the routine that has complexity!!!
        !
        ! We read chunks of the walkers datasets on _all_ of the processors,
        ! and communicate the walkers to the correct node after each of the
        ! blocks. This is essentially a giant annihilation step
        !
        ! We used the spawnedparts arrays as a buffer. These are split up per
        ! processor in the same way as normally used for annihilation. We
        ! keep looping over however much we can fit in this array unil
        ! we are done.
        !
        ! --> This is quite complicated, but should equally be quite efficient

        integer(hid_t), intent(in) :: parent
        integer(n_int), intent(out) :: dets(:, :)
        integer(int64), intent(out) :: CurrWalkers
        character(*), parameter :: t_r = 'read_walkers'

        integer :: proc, nreceived
        integer(hid_t) :: grp_id, err
        integer(hid_t) :: ds_sgns, ds_ilut
        integer(int64) :: nread_walkers

        integer(int32) :: bit_rep_width, tmp_lenof_sign
        integer(hsize_t) :: all_count, block_size, counts(0:nProcessors-1)
        integer(hsize_t) :: offsets(0:nProcessors-1)
        integer(hsize_t) :: block_start, block_end, last_part
        integer(hsize_t) :: this_block_size
        real(dp) :: pops_num_parts(lenof_sign), pops_norm_sqr(lenof_sign)
        real(dp) :: norm(lenof_sign), parts(lenof_sign)
        logical :: running, any_running
        integer(hsize_t), dimension(:,:), allocatable :: temp_ilut, temp_sgns
        integer :: temp_ilut_tag, temp_sgns_tag, rest

        ! TODO:
        ! - Read into a relatively small buffer. Make this such that all the
        !   particles could be broadcast to a single processor and everything
        !   would be fine (this is the only way to _ensure_ no overflows).
        ! - Loop over these read walkers, decode and determine which proc
        !   they are for
        ! - MPIAllToAllV to the correct processor
        ! - Ensure that these end up in the correct place in the main list,
        !   and that we don't overflow anywhere!
        ! - Writing and reading of flags, etc.
        !    -- For flags, include attribute info with which bit is which
        !       flag for supportability.

        call h5gopen_f(parent, nm_wfn_grp, grp_id, err)

        ! TODO: Make sure that this plays nicely with 32bit!!!
        if (.not. build_64bit) &
            call stop_all(t_r, "Needs manual, careful testing and adjustment")

        ! Get attributes about the sizes of stuff in the popsfile.
        ! ========

        ! We can only read in bit representations that are the right size!
        call read_int32_attribute(grp_id, nm_rep_width, bit_rep_width)
        if (bit_rep_width /= NIfD + 1) &
            call stop_all(t_r, "Mismatched bit representations")

        ! TODO: Deal with increasing the number of runs (e.g. for seeding RDMs)
        call read_int32_attribute(grp_id, nm_sgn_len, tmp_lenof_sign)
        if (lenof_sign /= tmp_lenof_sign) &
            call stop_all(t_r, "Mismatched sign length")

        call read_int64_attribute(grp_id, nm_num_dets, all_count, &
                                  required=.true.)
        call read_dp_1d_attribute(grp_id, nm_norm_sqr, pops_norm_sqr, &
                                  required=.true.)
        call read_dp_1d_attribute(grp_id, nm_num_parts, pops_num_parts, &
                                  required=.true.)
        write(6,*)
        write(6,*) 'Reading in ', all_count, ' determinants'

        ! How many particles should each processor read in (try and distribute
        ! this as uniformly as possible. Also calculate the associated data
        ! offsets
        counts = all_count / nProcessors
        rest=mod(all_count,nProcessors)
        if(rest.gt.0) counts(0:rest-1)=counts(0:rest-1)+1

        if (sum(counts) /= all_count .or. any(counts < 0)) &
            call stop_all(t_r, "Invalid particles counts")
        do proc = 0, nProcessors - 1
            offsets(proc) = sum(counts(0:proc-1))
        end do

        write(6,*) 'Reading in ', counts(iProcIndex), &
                   ' determinants on this process ...'

        ! TODO: Split this up into more manageable chunks!

        ! Open the relevant datasets
        call h5dopen_f(grp_id, nm_ilut, ds_ilut, err)
        call h5dopen_f(grp_id, nm_sgns, ds_sgns, err)

        ! Check that these datasets look like we expect them to.
        call check_dataset_params(ds_ilut, nm_ilut, 8_hsize_t, H5T_INTEGER_F, &
                                  [int(bit_rep_width, hsize_t), all_count])
        call check_dataset_params(ds_sgns, nm_sgns, 8_hsize_t, H5T_FLOAT_F, &
                                  [int(lenof_sign, hsize_t), all_count])

        !limit the buffer size per MPI task to 50MB or MaxSpawned entries
        block_size=50000000/(bit_rep_width*lenof_sign)/sizeof(SpawnedParts(1,1))
        block_size = min(block_size,MaxSpawned)
#if 0
        do proc = 0, nProcessors - 2
            block_size = min(block_size, &
                InitialSpawnedSlots(proc+1) - InitialSpawnedSlots(proc))
        end do
        block_size = min(block_size, &
                MaxSpawned - InitialSpawnedSlots(nProcessors-1))
#endif

        ! Initialise relevant counters
        CurrWalkers = 0
        nread_walkers = 0
        pops_norm = 0
        norm = 0
        parts = 0

        ! Note that reading in the HDF5 library is zero based, not one based
        ! TODO: We need to keep looping until all the processes are done, not
        !       just the one!
        block_start = offsets(iProcIndex)
        if (iProcIndex == nProcessors - 1) then
            last_part = all_count - 1
        else
            last_part = offsets(iProcIndex + 1) - 1
        end if
        block_end = min(block_start + block_size - 1, last_part)
        this_block_size = block_end - block_start + 1

        running = .true.
        any_running = .true.

        allocate(temp_ilut(int(bit_rep_width),int(this_block_size)))
        call LogMemAlloc('temp_ilut',size(temp_ilut),sizeof(temp_ilut(1,1)),'read_walkers',temp_ilut_tag,err)

        allocate(temp_sgns(int(lenof_sign),int(this_block_size)))
        call LogMemAlloc('temp_sgns',size(temp_sgns),sizeof(temp_sgns(1,1)),'read_walkers',temp_sgns_tag,err)

        do while (any_running)

            ! If this is the last block, its size will differ from the biggest
            ! one allowed.
            if (running) then
                this_block_size = block_end - block_start + 1
            else
                this_block_size = 0
            end if

            call read_walker_block_buff(ds_ilut, ds_sgns, block_start, &
                                   this_block_size, bit_rep_width, temp_ilut, temp_sgns)

            call distribute_and_add_walkers(this_block_size, temp_ilut, temp_sgns, dets, &
                 nreceived, CurrWalkers, norm, parts)

            nread_walkers = nread_walkers + nreceived
            
            ! And update for the next block
            if (running) then
                block_start = block_end + 1
                block_end = min(block_start + block_size - 1, last_part)

                if (block_start > last_part) running = .false.
            end if

            ! Test if _all_ of the processes have finished. If they have
            ! not then
            call MPIAllLORLogical(running, any_running)

        end do

        deallocate(temp_ilut, temp_sgns)
        call LogMemDeAlloc('read_walkers',temp_ilut_tag)
        call LogMemDeAlloc('read_walkers',temp_sgns_tag)


        call h5dclose_f(ds_sgns, err)
        call h5dclose_f(ds_ilut, err)
        call h5gclose_f(grp_id, err)

        ! Do some checking
        call check_read_particles(nread_walkers, norm, parts, all_count, &
                                  pops_num_parts, pops_norm_sqr)

        write(6,*) "... done"
        write(6,*)

    end subroutine read_walkers

    subroutine read_walker_block_buff(ds_ilut, ds_sgns, block_start, block_size, &
                                 bit_rep_width, temp_ilut, temp_sgns)

        use bit_rep_data, only: NIfD
        use FciMCData, only: SpawnedParts2

        ! Read the walkers into the array spawnedparts2
        !
        ! N.B. This routine is quite sensitive to the particular structure
        !      of the bit representations determined in BitReps.F90
        !
        ! --> It would also be possible to read into scratch arrays, and then
        !     do some transferring.

        integer(hid_t), intent(in) :: ds_ilut, ds_sgns
        integer(hsize_t), intent(in) :: block_start, block_size
        integer(int32), intent(in) :: bit_rep_width
        integer(hsize_t), dimension(:,:) :: temp_ilut, temp_sgns
        integer(hid_t) :: plist_id

#ifdef __INT64

        call read_2d_multi_chunk( &
                ds_ilut, temp_ilut, H5T_NATIVE_INTEGER_8, &
                [int(bit_rep_width, hsize_t), block_size], &
                [0_hsize_t, block_start], &
                [0_hsize_t, 0_hsize_t])

        call read_2d_multi_chunk( &
             ds_sgns, temp_sgns, H5T_NATIVE_REAL_8, &
             [int(lenof_sign, hsize_t), block_size], &
             [0_hsize_t, block_start], &
             [0_hsize_t, 0_hsize_t])

#else
        call stop_all("read_walker_block", "32-64bit conversion not yet implemented")
#endif

        ! TODO: Flags here!!!

    end subroutine read_walker_block_buff


    subroutine distribute_and_add_walkers(block_size, temp_ilut, temp_sgns, dets, &
         nreceived, CurrWalkers, norm, parts)
      integer(n_int), intent(out) :: dets(:, :)
      integer(int64), intent(inout) :: CurrWalkers
      integer(hsize_t) :: block_size
      integer:: nreceived
      real(dp), intent(inout) :: norm(lenof_sign), parts(lenof_sign)
      integer(hsize_t):: temp_ilut(:,:), temp_sgns(:,:)
      integer:: sendcount(0:nProcessors-1), nlocal=0

      call assign_dets_to_procs_buff(block_size, temp_ilut, temp_sgns, sendcount)

      !add elements that are on the right processor already
#define localfirst
!#undef localfirst
#ifdef localfirst
      nlocal=sendcount(iProcIndex)
      call add_new_parts(dets, nlocal, CurrWalkers, norm, parts)      
      sendcount(iProcIndex)=0
#endif
      !communicate the remaining elements
      nreceived = communicate_read_walkers_buff(sendcount)
      call add_new_parts(dets, nreceived, CurrWalkers, norm, parts)

      if (allocated(receivebuff)) then
         deallocate(receivebuff)
         call LogMemDeAlloc('distribute_and_add',receivebuff_tag)
      end if

      nreceived=nreceived+nlocal

    end subroutine distribute_and_add_walkers

    subroutine assign_dets_to_procs_buff(block_size, temp_ilut, temp_sgns, sendcount)

        use load_balance_calcnodes, only: DetermineDetNode
        use bit_reps, only: decode_bit_det, extract_sign
        use FciMCData, only: SpawnedParts2, SpawnedParts, ValidSpawnedList

        use Determinants, only: write_det
        use bit_rep_data, only: NIfD, NIFBCast
        use SystemData, only: nel

        integer(hsize_t), intent(in) :: block_size
        character(*), parameter :: t_r = 'distribute_walkers_from_block'
        integer(hsize_t), dimension(:,:) :: temp_ilut, temp_sgns
        integer(hsize_t) :: onepart(0:NIfBCast)
        integer :: det(nel), p, j, proc, ierr, sizeilut, targetproc(block_size)
        integer:: sendcount(0:nProcessors-1), index, index2
        logical :: list_full

        sizeilut=size(temp_ilut,1)

        ! Iterate through walkers in temp_ilut+temp_sgns and determine the target processor. 
        onepart=0
        sendcount=0
        do j = 1, block_size
            onepart(0:sizeilut-1)=temp_ilut(:,j)
            onepart(sizeilut:sizeilut+int(lenof_sign)-1)=temp_sgns(:,j)
            ! Which processor does this determinant live on?
            call decode_bit_det(det, onepart)
            proc = DetermineDetNode(nel, det, 0)
            targetproc(j)=proc
            sendcount(proc)=sendcount(proc)+1
        end do
        
        ! Write the elements to SpawnedParts in the correct order for sending
        index=1
        index2=1
        do p = 0, nProcessors-1
#ifdef localfirst
           if (p.eq.iProcIndex) then
#else
           if (.false.) then
#endif
              !elements that don't have to be communicated are written to SpawnedParts2 
              do j = 1, block_size
                 if(targetproc(j).eq.p) then
                    onepart(0:sizeilut-1)=temp_ilut(:,j)
                    onepart(sizeilut:sizeilut+int(lenof_sign)-1)=temp_sgns(:,j)
                    SpawnedParts2(:,index2)=onepart
                    index2=index2+1
                 end if
              end do
           else
              !elements that have to be sent to other procs are written to SpawnedParts 
              do j = 1, block_size
                 if(targetproc(j).eq.p) then
                    onepart(0:sizeilut-1)=temp_ilut(:,j)
                    onepart(sizeilut:sizeilut+int(lenof_sign)-1)=temp_sgns(:,j)
                    SpawnedParts(:,index)=onepart
                    index=index+1
                 end if
              end do
           end if
        end do

    end subroutine assign_dets_to_procs_buff


    function communicate_read_walkers_buff(sendcounts) result(num_received)
        integer:: sendcounts(0:nProcessors-1)
        integer :: num_received
        integer(int64) :: lnum_received

        integer(MPIArg) :: recvcounts(0:nProcessors-1)
        integer(MPIArg) :: disps(0:nProcessors-1), recvdisps(0:nProcessors-1)

        integer :: j, ierr


        !offsets for data to the different procs
        disps(0) = 0
        do j = 1, nProcessors-1
            disps(j) = disps(j-1)+sendcounts(j-1)
        end do

        ! Communicate the number of particles that need to go to each proc
        call MPIAllToAll(sendcounts, 1, recvcounts, 1, ierr)


        ! We want the data to be contiguous after the move. So calculate the
        ! offsets
        recvdisps(0) = 0
        do j = 1, nProcessors-1
            recvdisps(j) = recvdisps(j-1) + recvcounts(j-1)
        end do
        num_received = recvdisps(nProcessors-1) + recvcounts(nProcessors-1)
        lnum_received = recvdisps(nProcessors-1) + recvcounts(nProcessors-1)

        ! Adjust offsets so that they match the size of the array
        recvcounts = recvcounts * size(SpawnedParts, 1)
        recvdisps = recvdisps * size(SpawnedParts, 1)
        sendcounts = sendcounts * size(SpawnedParts, 1)
        disps = disps * size(SpawnedParts, 1)

        if (num_received.gt.size(SpawnedParts2,2)) then
           ! there could in principle be a memory problem because we are not limiting the
           ! size of receivebuff.
           write(6,*) 'Allocating additional buffer for communication on Processor ', iProcIndex, 'with ', &
                num_received*size(SpawnedParts,1)*sizeof(SpawnedParts(1,1))/1000000, 'MB'
           allocate(receivebuff(size(SpawnedParts,1),num_received))
           call LogMemAlloc('receivebuff',size(receivebuff),sizeof(receivebuff(1,1)),&
                'communicate_read_walkers',receivebuff_tag,ierr)
           call MPIAllToAllV(SpawnedParts, sendcounts, disps, receivebuff, &
                recvcounts, recvdisps, ierr)
        else
           call MPIAllToAllV(SpawnedParts, sendcounts, disps, SpawnedParts2, &
                recvcounts, recvdisps, ierr)
        end if
      end function communicate_read_walkers_buff

    subroutine add_new_parts(dets, nreceived, CurrWalkers, norm, parts)

        use CalcData, only: iWeightPopRead
        use bit_reps, only: extract_sign

        ! Integrate the just-read block of walkers into the main list.

        integer(n_int), intent(inout) :: dets(:, :)
        integer, intent(in) :: nreceived
        integer(int64), intent(inout) :: CurrWalkers
        real(dp), intent(inout) :: norm(lenof_sign), parts(lenof_sign)

        integer(int64) :: j
        real(dp) :: sgn(lenof_sign)

            ! TODO: inum_runs == 2, PopNIfSgn == 1
        if (allocated(receivebuff)) then
           
           do j = 1, nreceived
              
              ! Check that the site is occupied, and passes the relevant
              ! thresholds before adding it to the system.
              call extract_sign(receivebuff(: ,j), sgn)
              if (any(abs(sgn) >= iWeightPopRead) .and. .not. IsUnoccDet(sgn)) then
                 
                 ! Add this site to the main list
                 CurrWalkers = CurrWalkers + 1
                 dets(:, CurrWalkers) = receivebuff(:, j)
                 call add_pops_norm_contrib(dets(:, CurrWalkers))
                 call extract_sign(receivebuff(:,j), sgn)
                 norm = norm + sgn**2
                 parts = parts + abs(sgn)
              end if
           end do
        else
           do j = 1, nreceived
              
              ! Check that the site is occupied, and passes the relevant
              ! thresholds before adding it to the system.
              call extract_sign(SpawnedParts2(: ,j), sgn)
              if (any(abs(sgn) >= iWeightPopRead) .and. .not. IsUnoccDet(sgn)) then
                 
                 ! Add this site to the main list
                 CurrWalkers = CurrWalkers + 1
                 dets(:, CurrWalkers) = SpawnedParts2(:, j)
                 call add_pops_norm_contrib(dets(:, CurrWalkers))
                 call extract_sign(SpawnedParts2(:,j), sgn)
                 norm = norm + sgn**2
                 parts = parts + abs(sgn)
              end if
           end do
           
        endif
        ! TODO: Add check that we have read in the correct number of parts

    end subroutine

    subroutine check_read_particles(nread_walkers, norm, parts, &
                                    pops_det_count, pops_num_parts, &
                                    pops_norm_sqr)

        use CalcData, only: pops_norm

        ! Check that the values received in these routines are valid
        !
        ! nread_walkers  - Number of determinant/walker lines read (this proc)
        ! norm           - Norm (this proc)
        ! parts          - Coefficient weight (this proc)
        ! pops_det_count - Total dets in popsfile (from header)
        ! pops_num_parts - Total particle weight (from header)
        ! pops_norm_sqr  - Total norm of wavefunction (from header)

        integer(int64), intent(in) :: nread_walkers, pops_det_count
        real(dp), intent(in) :: pops_num_parts(lenof_sign)
        real(dp), intent(in) :: pops_norm_sqr(lenof_sign)
        real(dp), intent(in) :: norm(lenof_sign), parts(lenof_sign)
        character(*), parameter :: t_r = 'check_read_particles'

        integer(int64) :: all_read_walkers, tot_walkers
        real(dp) :: all_norm(lenof_sign), all_parts(lenof_sign)

        ! Have all the sites been correctly read in from the file
        ! n.b. CurrWalkers may not equal pops_det_count, as any unoccupied
        !      sites, or sites below the specified threshold will have been
        !      dropped.
        call MPISum(nread_walkers, all_read_walkers)
        if (iProcIndex == 0 .and. all_read_walkers /= pops_det_count) then
            write(6,*) 'Determinants in popsfile header: ', pops_det_count
            write(6,*) 'Determinants read in: ', all_read_walkers
            call stop_all(t_r, 'Particle number mismatch')
        end if

        ! Is the total number of walkers (sum of the sign values) correct?
        !
        ! This is relative to the total number to account for relative errors
        call MPISumAll(parts, all_parts)
        if (any(abs(all_parts - pops_num_parts) > (pops_num_parts * 1.0e-10_dp))) then
            write(6,*) 'popsfile particles: ', pops_num_parts
            write(6,*) 'read particles: ', all_parts
            call stop_all(t_r, 'Incorrect particle weight read from popsfile')
        end if

        ! Is the total norm of the wavefunction correct
        ! The test condition is rather lax, as the hugely differing magnitudes
        ! of the total numbers on each processor combined with the varying
        ! summation orders can lead to larger errors here
        !
        ! Same relative error test as before
        call MPISumAll(norm, all_norm)
        if (any(abs(all_norm - pops_norm_sqr) > (pops_norm_sqr * 1.0e-10_dp))) then
            write(6,*) 'popsfile norm**2: ', pops_norm_sqr
            write(6,*) 'read norm**2: ', all_norm
            call stop_all(t_r, 'Wavefunction norm incorrect')
        end if

        ! If the absolute sum, and the sum of the squares is correct, we can
        ! be fairly confident that they have all been read in!...
        
        ! [W.D.]
        ! on behalf of sasha bring back the feature that turns off the walker 
        ! grow even if walkcontgrow was set unintentionally but the number of 
        ! read-in walkers already exceeds or is close to the target number 
        ! so i guess it is enough to set the global AllTotParts here 
        ! of walkers 
        AllTotParts = all_parts

    end subroutine
#endif

    !
    ! This is only here for dependency circuit breaking
    subroutine add_pops_norm_contrib(ilut)

        use bit_rep_data, only: NIfTot
        use bit_reps, only: extract_sign
        use CalcData, only: pops_norm

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp) :: real_sign(lenof_sign)

        call extract_sign(ilut, real_sign)

#ifdef __DOUBLERUN
        pops_norm = pops_norm + real_sign(1)*real_sign(2)
#elif __CMPLX
        pops_norm = pops_norm + real_sign(1)**2 + real_sign(2)**2
#else
        pops_norm = pops_norm + real_sign(1)*real_sign(1)
#endif

    end subroutine add_pops_norm_contrib

end module
