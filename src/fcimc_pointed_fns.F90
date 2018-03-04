#include "macros.h"

module fcimc_pointed_fns

    use SystemData, only: nel, tGUGA, tGen_4ind_weighted, tGen_4ind_2, tGen_nosym_guga, &
                          tGen_sym_guga_mol, t_consider_diff_bias, nSpatOrbs, thub, & 
                          tUEG, tGen_4ind_reverse, nBasis
    use LoggingData, only: tHistExcitToFrom, FciMCDebug
    use CalcData, only: RealSpawnCutoff, tRealSpawnCutoff, tAllRealCoeff, &
                        RealCoeffExcitThresh, AVMcExcits, tau, DiagSft, &
                        tRealCoeffByExcitLevel, InitiatorWalkNo, &
                        t_fill_frequency_hists, t_truncate_spawns, n_truncate_spawns, & 
                        t_matele_cutoff, matele_cutoff 
    use DetCalcData, only: FciDetIndex, det
    use procedure_pointers, only: get_spawn_helement, log_spawn_magnitude
    use fcimc_helper, only: CheckAllowedTruncSpawn
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, count_open_orbs
    use bit_rep_data, only: NIfTot, test_flag
    use bit_reps, only: get_initiator_flag
    use tau_search, only: log_death_magnitude, log_spawn_magnitude
    use rdm_general, only: calc_rdmbiasfac
    use hist, only: add_hist_excit_tofrom
    use searching, only: BinSearchParts2
    use util_mod
    use FciMCData
    use constants
    use bit_reps, only: nifguga

#ifndef __CMPLX
#ifdef __DEBUG
    use guga_bitRepOps, only: convert_ilut_toGUGA, write_det_guga
    use guga_excitations, only: print_excitInfo, global_excitInfo
#endif
#endif
    
    use tau_search_hist, only: fill_frequency_histogram_4ind, &
                               fill_frequency_histogram_sd, &
                               fill_frequency_histogram

    use excit_gen_5, only: pgen_select_a_orb

    implicit none

    contains

    function attempt_create_trunc_spawn (DetCurr,&
                                         iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type, &
                                         AvSignCurr, RDMBiasFacCurr) result(child)
        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t(dp), intent(in) :: HElGen

        if (CheckAllowedTruncSpawn (walkExcitLevel, nJ, iLutnJ, IC)) then
            child = attempt_create_normal (DetCurr, &
                               iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, ic, ex, &
                               tParity, walkExcitLevel, part_type, AvSignCurr, RDMBiasFacCurr)
        else
            child = 0
        endif
    end function

!Decide whether to spawn a particle at nJ from DetCurr. (bit strings iLutnJ and iLutCurr respectively).  
!  ic and ex specify the excitation of nJ from DetCurr, along with the sign change tParity.
!  part_type:           Is the parent real (1) or imaginary (2)
!  wSign:               wSign gives the sign of the particle we are trying to spawn from
!                          if part_type is 1, then it will only use wsign(1)
!                                          2,                       wsign(2)
!                       Only the sign, not magnitude is used.
!  prob:                prob is the generation probability of the excitation in order to unbias.
!                       The probability of spawning is divided by prob to do this.
!  HElGen:              If the matrix element has already been calculated, it is sent in here.
!  get_spawn_helement:  A function pointer for looking up or calculating the relevant matrix element.
!  walkExcitLevel:      Is Unused
! 
!  child:      A lenof_sign array containing the particles spawned.
    function att_create_trunc_spawn_enc (DetCurr,&
                                         iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, &
                                         ic, ex, tparity, walkExcitLevel, part_type, &
                                         AvSignCurr,RDMBiasFacCurr) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type 
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t(dp) , intent(in) :: HElGen

        call EncodeBitDet (nJ, iLutnJ)
        if (CheckAllowedTruncSpawn (walkExcitLevel, nJ, iLutnJ, IC)) then
            child = attempt_create_normal (DetCurr, &
                               iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, ic, ex, &
                               tParity, walkExcitLevel, part_type, AvSignCurr, RDMBiasFacCurr)
        else
            child = 0
        endif
    end function

    function attempt_create_normal (DetCurr, iLutCurr, &
                                    RealwSign, nJ, iLutnJ, prob, HElGen, ic, ex, tParity,&
                                    walkExcitLevel, part_type, AvSignCurr, RDMBiasFacCurr) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! odd = Real parent particle, even = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp) , dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        HElement_t(dp) , intent(in) :: HElGen
        character(*), parameter :: this_routine = 'attempt_create_normal'

        real(dp) :: rat, r, walkerweight, pSpawn, nSpawn, MatEl, p_spawn_rdmfac
        integer :: extracreate, tgt_cpt, component, i, iUnused
        integer :: TargetExcitLevel
        logical :: tRealSpawning
        HElement_t(dp) :: rh, rh_used
#ifdef __DEBUG
        integer :: nOpen
        integer(n_int) :: ilutTmpI(0:nifguga), ilutTmpJ(0:nifguga)
#endif
        logical :: t_par
        real(dp) :: temp_prob, pgen_a, dummy_arr(nBasis), cum_sum
        integer :: ispn

        ! Just in case
        child = 0.0_dp

        ! If each walker does not have exactly one spawning attempt
        ! (if AvMCExcits /= 1.0_dp) then the probability of an excitation
        ! having been chosen, prob, must be altered accordingly.
        prob = prob * AvMCExcits

        ! In the case of using HPHF, and when tGenMatHEl is on, the matrix
        ! element is calculated at the time of the excitation generation, 
        ! and returned in HElGen. In this case, get_spawn_helement simply
        ! returns HElGen, rather than recomputing the matrix element.
        rh = get_spawn_helement (DetCurr, nJ, iLutCurr, iLutnJ, ic, ex, &
                                 tParity, HElGen)

!         if (abs(rh) > EPS) then
!             print *, "HElGen: ", rh
!             print *, "prob: ", prob
!         end if
        ! We actually want to calculate Hji - take the complex conjugate, 
        ! rather than swap around DetCurr and nJ.
#ifdef __CMPLX
        rh_used = conjg(rh)
#else
        rh_used = rh
#endif
        
        ! [W.D.]
        ! if the matrix element happens to be zero, i guess i should 
        ! abort as early as possible? so check that here already, or even 
        ! earlier.. 
!         if (abs(rh_used) < EPS) then 
!             child = 0.0_dp
!             return
!         end if
! 
!         if (t_matele_cutoff) then
!             if (abs(rh_used) < EPS) then
!                 child = 0.0_dp
!             end if
!         end if

        !write(6,*) 'p,rh', prob, rh

        ! essentially here i have all the information for my frequency 
        ! analysis of the H_ij/pgen ratio so call the routine here
        ! but i have to remember to keep it parallel! so dont forget to 
        ! sum up all the contributions from different cores! 
        ! and divide prob by AvMCExcits again to get correct pgen! 
        if (t_frequency_analysis .and. t_fill_frequency_hists) then
            ! use specific ones for different types of excitation gens
            if (tGen_4ind_weighted .or. tGen_4ind_2) then
                ! determine if excitation was parallel or anti-parallel
                ! ex(1,1) and ex(1,2) are the electrons 
                t_par = (is_beta(ex(1,1)) .eqv. is_beta(ex(1,2)))

                call fill_frequency_histogram_4ind(abs(rh), prob / AvMCExcits, &
                    ic, t_par)

            else if (tGen_nosym_guga) then 
                ! have to also check if diff bias is considered 
                if (t_consider_diff_bias) then 
                    call fill_frequency_histogram_nosym_diff(abs(rh), prob / AvMCExcits, & 
                        ic, ex(1,1), ex(1,2))
                else 
                    call fill_frequency_histogram_nosym_nodiff(abs(rh), &
                        prob / AvMCExcits, ic, ex(1,1))
                end if

            else if (tGen_sym_guga_mol) then
                call fill_frequency_histogram_sd(abs(rh), prob / AvMCExcits, ic)

            else 
                ! for any other excitation generator just use one histogram 
                ! for all the excitations.. 
                call fill_frequency_histogram(abs(rh), prob / AvMCExcits)

            end if
        end if

        ! The following is useful for debugging the contributions of single
        ! excitations, and double excitations of spin-paired/opposite
        ! electron pairs to the value of tau.
!        if (ic == 2) then
!            if (G1(ex(1,1))%Ms /= G1(ex(1,2))%Ms) then
!                write(6,*) 'OPP', rh, prob
!            else
! !                write(6,*) 'SAM', rh, prob
!            end if
!        else
!            write(6,*) 'IC1', rh, prob
!        end if

        ! fill in the frequency histograms here! 
        ! [Werner Dobrautz 4.4.2017:]
        if (t_fill_frequency_hists) then 
            if (tHUB .or. tUEG) then 
                call fill_frequency_histogram(abs(rh_used), prob / AvMCExcits)
            else 
                if (tGen_4ind_2 .or. tGen_4ind_weighted .or. tGen_4ind_reverse) then 
                    t_par = (is_beta(ex(1,1)) .eqv. is_beta(ex(1,2)))

                    ! not sure about the AvMCExcits!! TODO
                    call fill_frequency_histogram_4ind(abs(rh_used), prob / AvMCExcits, &
                        ic, t_par, ex)

                else

                    call fill_frequency_histogram_sd(abs(rh_used), prob / AvMCExcits, ic)
                    
                end if
            end if
        end if
        ! Are we doing real spawning?
        
        tRealSpawning = .false.
        if (tAllRealCoeff) then
            tRealSpawning = .true.
        elseif (tRealCoeffByExcitLevel) then
            if (tGUGA) call stop_all(this_routine,&
                "excit level does not work with GUGA here...")

            TargetExcitLevel = FindBitExcitLevel (iLutRef(:,1), iLutnJ)

            if (TargetExcitLevel <= RealCoeffExcitThresh) &
                tRealSpawning = .true.
        endif

        ! Spawn to real and imaginary particles. Note that spawning from
        ! imaginary parent particles has slightly different rules:
        !       - Attempt to spawn REAL walkers with prob +AIMAG(Hij)/P
        !       - Attempt to spawn IMAG walkers with prob -REAL(Hij)/P



#if !defined(__CMPLX) && (defined(__PROG_NUMRUNS) || defined(__DOUBLERUN))
        child = 0.0_dp
        tgt_cpt = part_type
        walkerweight = sign(1.0_dp, RealwSign(part_type))
        matEl = real(rh_used, dp)
        if (t_matele_cutoff) then
            if (abs(matEl) < matele_cutoff) matel = 0.0_dp
        end if
#else
        do tgt_cpt = 1, (lenof_sign/inum_runs)

            ! Real, single run:    inum_runs=1, lenof_sign=1 --> 1 loop
            ! Real, double run:    inum_runs=2, lenof_sign=1 --> 1 loop
            ! Complex, single run: inum_runs=1, lenof_sign=2 --> 2 loops
            ! Complex, double run: inum_runs=2, lenof_sign=4 --> 2 loops
            ! Complex, multiple run: inum_runs=m, lenof_sign=2*m --> 2 loops

            ! For spawning from imaginary particles, we cross-match the 
            ! real/imaginary matrix-elements/target-particles.


#if defined(__CMPLX) && (defined(__PROG_NUMRUNS) || defined(__DOUBLERUN))
            component = part_type+tgt_cpt-1
            if (.not. btest(part_type,0)) then
                ! even part_type => imag replica =>  map 4->3,4 ; 6->5,6 etc.
                component = part_type - tgt_cpt + 1
            endif
#else
            component = tgt_cpt
            if ((part_type.eq.2).and.(inum_runs.eq.1)) component = 3 - tgt_cpt
#endif

            ! Get the correct part of the matrix element
            walkerweight = sign(1.0_dp, RealwSign(part_type))
            if (btest(component,0)) then
                ! real component
                MatEl = real(rh_used, dp)
                if (t_matele_cutoff) then
                    if (abs(MatEl) < matele_cutoff) MatEl = 0.0_dp
                end if
            else
#ifdef __CMPLX
                MatEl = real(aimag(rh_used), dp)
                if (t_matele_cutoff) then
                    if (abs(MatEl) < matele_cutoff) MatEl = 0.0_dp
                end if
                ! n.b. In this case, spawning is of opposite sign.
                if (.not. btest(part_type,0)) then
                    ! imaginary parent -> imaginary child
                    walkerweight = -walkerweight
                endif
#endif
            end if
#endif
            nSpawn = - tau * MatEl * walkerweight / prob

            ! so.. does a |nSpawn| > 1 in my new tau-search implitly mean 
            ! that the tau is to small for this kind of exciation? 
            ! i guess so yeah.. so do it really brute force to start with 
            if (t_truncate_spawns .and. abs(nSpawn) > n_truncate_spawns) then
                

                ! in debug mode i should output some additional information 
                ! to analyze the type of excitations and how many open-orbitals
                ! etc. are used 
#ifndef __CMPLX
#ifdef __DEBUG 
            if (abs(nSpawn) > 10.0_dp) then
            if (tGUGA) then
                write(iout,*) "=================================================="
                call convert_ilut_toGUGA(iLutCurr, ilutTmpI)
                call convert_ilut_toGUGA(ilutnj, ilutTmpJ)
                write(iout,*) "nSpawn > n_truncate_spawns!", nSpawn
                write(iout,*) "limit the number of spawned walkers to: ", n_truncate_spawns
                write(iout,*) "for spawn from determinant: "
                call write_det_guga(6,ilutTmpI) 
                write(iout,*) "to: " 
                call write_det_guga(iout, ilutTmpJ)
                nOpen = count_open_orbs(iLutCurr)
                write(iout,*) "# of openshell orbitals: ", nOpen, count_open_orbs(ilutnj)
                write(iout,*) "open/spatial: ", nOpen/real(nSpatOrbs,dp)
                write(iout,*) "(t |H_ij|/pgen) / #open ratio: ", abs(nSpawn) / real(nOpen,dp)
                write(iout,*) " H_ij, pgen: ", MatEl, prob
                write(iout,*) "=================================================="
                ! excitation type would be cool too.. but how do i get it 
                ! to here?? do i still have global_excitInfo??
                call print_excitInfo(global_excitInfo)
                call neci_flush(iout)
            end if
            end if
#endif
#endif

                ! [Werner Dobrautz 4.4.2017:]
                ! apply the spawn truncation, when using histogramming tau-search
                nSpawn = sign(n_truncate_spawns, nSpawn)

            end if
            
            ! n.b. if we ever end up with |walkerweight| /= 1, then this
            !      will need to ffed further through.
            if (tSearchTau .and. (.not. tFillingStochRDMonFly)) then
                ! in the back-spawning i have to adapt the probabilites 
                ! back, to be sure the time-step covers the changed 
                ! non-initiators spawns! 

                call log_spawn_magnitude (ic, ex, matel, prob)

            end if

            ! Keep track of the biggest spawn this cycle
            max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)
            
            if (tRealSpawning) then
                ! Continuous spawning. Add in acceptance probabilities.
                
                if (tRealSpawnCutoff .and. &
                    abs(nSpawn) < RealSpawnCutoff) then
                    p_spawn_rdmfac=abs(nSpawn)/RealSpawnCutoff
                    nSpawn = RealSpawnCutoff &
                           * stochastic_round (nSpawn / RealSpawnCutoff)
               else
                    p_spawn_rdmfac=1.0_dp !The acceptance probability of some kind of child was equal to 1
               endif
            else
                if(abs(nSpawn).ge.1) then
                    p_spawn_rdmfac=1.0_dp !We were certain to create a child here.
                    ! This is the special case whereby if P_spawn(j | i) > 1, 
                    ! then we will definitely spawn from i->j.
                    ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                    ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                    ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                else
                    p_spawn_rdmfac=abs(nSpawn)
                endif
                
                ! How many children should we spawn?

                ! And round this to an integer in the usual way
                ! HACK: To use the same number of random numbers for the tests.
                nSpawn = real(stochastic_round (nSpawn), dp)
                
            endif
            ! And create the parcticles
#ifdef __CMPLX
            child((part_type_to_run(part_type)-1)*2+tgt_cpt) = nSpawn
#else
            child(tgt_cpt) = nSpawn
#endif

#if defined(__CMPLX) || !defined(__PROG_NUMRUNS) && !defined(__DOUBLERUN)
        enddo
#endif

       
        if(tFillingStochRDMonFly) then
            if (child(part_type).ne.0) then
                !Only add in contributions for spawning events within population 1
                !(Otherwise it becomes tricky in annihilation as spawnedparents doesn't tell you which population
                !the event came from at present)
                call calc_rdmbiasfac(p_spawn_rdmfac, prob, realwSign(part_type), RDMBiasFacCurr) 
            else
                RDMBiasFacCurr = 0.0_dp
            endif
        else
            ! Not filling the RDM stochastically, bias is zero.
            RDMBiasFacCurr = 0.0_dp
        endif

        ! Avoid compiler warnings
        iUnused = walkExcitLevel

    end function

    ! 
    ! This is a null routine for encoding spawned sites
    ! --> DOES NOTHING!!!
    subroutine null_encode_child (ilutI, ilutJ, ic, ex)
        implicit none
        integer(kind=n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: ic, ex(2,2)
        integer(kind=n_int), intent(inout) :: ilutj(0:niftot)

        ! Avoid compiler warnings
        integer :: iUnused
        integer(n_int) :: iUnused2
        iLutJ(0) = iLutJ(0); iUnused = IC; iUnused = ex(2,2)
        iUnused2 = iLutI(0)
    end subroutine

    subroutine new_child_stats_hist_hamil (iter_data, iLutI, nJ, iLutJ, ic, &
                                           walkExLevel, child, parent_flags, &
                                           part_type)
        ! Based on old AddHistHamilEl. Histograms the hamiltonian matrix, and 
        ! then calls the normal statistics routine.

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        real(dp), dimension(lenof_sign) , intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'new_child_stats_hist_hamil'
        integer :: partInd, partIndChild, childExLevel
        logical :: tSuccess

        if (walkExLevel == nel) then
            call BinSearchParts2 (iLutI, FCIDetIndex(walkExLevel), Det, &
                                  PartInd, tSuccess)
        else
            call BinSearchParts2 (iLutI, FCIDetIndex(walkExLevel), &
                                  FciDetIndex(walkExLevel+1)-1, partInd, &
                                  tSuccess)
        endif

        if (.not. tSuccess) &
            call stop_all (this_routine, 'Cannot find determinant nI in list')

        childExLevel = FindBitExcitLevel (iLutHF, iLutJ, nel)
        if (tGUGA) call stop_all(this_routine, &
            "excit level does not work with GUGA here...")

        if (childExLevel == nel) then
            call BinSearchParts2 (iLutJ, FCIDetIndex(childExLevel), Det, &
                                  partIndChild, tSuccess)
        elseif (childExLevel == 0) then
            partIndChild = 1
            tSuccess = .true.
        else
            call BinSearchParts2 (iLutJ, FCIDetIndex(childExLevel), &
                                  FciDetIndex(childExLevel+1)-1, &
                                  partIndChild, tSuccess)
        endif

        histHamil (partIndChild, partInd) = &
                histHamil (partIndChild, partInd) + (1.0_dp * child(1))
        histHamil (partInd, partIndChild) = &
                histHamil (partInd, partIndChild) + (1.0_dp * child(1))
        avHistHamil (partIndChild, partInd) = &
                avHistHamil (partIndChild, partInd) + (1.0_dp * child(1))
        avHistHamil (partInd, partIndChild) = &
                avHistHamil (partInd, partIndChild) + (1.0_dp * child(1))

        ! Call the normal stats routine
        call new_child_stats_normal (iter_data, iLutI, nJ, iLutJ, ic, &
                                     walkExLevel, child, parent_flags, &
                                     part_type)

    end subroutine

    subroutine new_child_stats_normal (iter_data, iLutI, nJ, iLutJ, ic, &
                                       walkExLevel, child, parent_flags, &
                                       part_type)

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        real(dp), dimension(lenof_sign), intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(n_int) :: iUnused
        integer :: run
        integer :: i

        ! Write out some debugging information if asked
        IFDEBUG(FCIMCDebug,3) then
            write(iout,"(A)",advance='no') "Creating "
            do i = 1,lenof_sign
                write(iout,"(f10.5)",advance='no') child(i)
            enddo
            write(iout,"(A)",advance='no') " particles: "
            write(iout,"(A,2I4,A)",advance='no') &
                                      "Parent flag: ", parent_flags, part_type
            call writebitdet (iout, ilutJ, .true.)
            call neci_flush(iout)
        endif

        ! Count the number of children born
#ifdef __CMPLX
        do run = 1, inum_runs
            NoBorn(run) = NoBorn(run) + sum(abs(child(min_part_type(run):max_part_type(run))))
            if (ic == 1) SpawnFromSing(run) = SpawnFromSing(run) + sum(abs(child(min_part_type(run):max_part_type(run))))

        
           ! Count particle blooms, and their sources
            if (sum(abs(child(min_part_type(run):max_part_type(run)))) > InitiatorWalkNo) then
                bloom_count(ic) = bloom_count(ic) + 1
                bloom_sizes(ic) = max(real( sum(abs(child(min_part_type(run):max_part_type(run)))),dp), bloom_sizes(ic))
            end if
        enddo
#else
        NoBorn = NoBorn + abs(child)
        if (ic == 1) SpawnFromSing = SpawnFromSing + abs(child)

        ! Count particle blooms, and their sources
        if (abs(child(part_type)) > InitiatorWalkNo) then
            bloom_count(ic) = bloom_count(ic) + 1
            bloom_sizes(ic) = max(real((abs(child(part_type))), dp), bloom_sizes(ic))
        end if
#endif
        iter_data%nborn = iter_data%nborn + abs(child)

        ! Histogram the excitation levels as required
        if (tHistExcitToFrom) &
            call add_hist_excit_tofrom(ilutI, ilutJ, child)

        ! Avoid compiler warnings
        iUnused = iLutI(0); iUnused = iLutJ(0)

    end subroutine

    function attempt_die_normal (DetCurr, Kii, realwSign, WalkExcitLevel) result(ndie)
        
        ! Should we kill the particle at determinant DetCurr. 
        ! The function allows multiple births (if +ve shift), or deaths from
        ! the same particle. The returned number is the number of deaths if
        ! positive, and the
        !
        ! In:  DetCurr - The determinant to consider
        !      Kii     - The diagonal matrix element of DetCurr (-Ecore)
        !      wSign   - The sign of the determinant being considered. If
        !                |wSign| > 1, attempt to die multiple particles at
        !                once (multiply probability of death by |wSign|)
        ! Ret: ndie    - The number of deaths (if +ve), or births (If -ve).

        integer, intent(in) :: DetCurr(nel)
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        real(dp), intent(in) :: Kii
        real(dp), dimension(lenof_sign) :: ndie
        integer, intent(in) :: WalkExcitLevel
        character(*), parameter :: t_r = 'attempt_die_normal'

        real(dp) :: probsign, r
        real(dp), dimension(inum_runs) :: fac
        integer :: i, run, iUnused
#ifdef __CMPLX
        real(dp) :: rat(2)
#else
        real(dp) :: rat(1)
#endif        

        do i=1, inum_runs
            fac(i)=tau*(Kii-DiagSft(i))

            ! And for tau searching purposes
            call log_death_magnitude (Kii - DiagSft(i))
        enddo

        if(any(fac > 1.0_dp)) then
            if (any(fac > 2.0_dp)) then
                if (tSearchTau) then
                    ! If we are early in the calculation, and are using tau
                    ! searching, then this is not a big deal. Just let the
                    ! searching deal with it
                    write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                    write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                    do i = 1, inum_runs
                        fac(i) = min(2.0_dp, fac(i))
                    end do
                else
                    call stop_all(t_r, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                end if
            else
                write(iout,'("** WARNING ** Death probability > 1: Creating Antiparticles. "&
                    & //"Timestep errors possible: ")',advance='no')
                do i = 1, inum_runs
                    write(iout,'(1X,f13.7)',advance='no') fac(i)
                end do
                write(iout,'()')
            endif
        endif


        if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel .le. RealCoeffExcitThresh)) &
            .or. tAllRealCoeff ) then
            do run=1, inum_runs
                ndie(min_part_type(run))=fac(run)*abs(realwSign(min_part_type(run)))
#ifdef __CMPLX
                ndie(max_part_type(run))=fac(run)*abs(realwSign(max_part_type(run)))
#endif
            enddo
        else
            do run=1,inum_runs
                
                ! Subtract the current value of the shift, and multiply by tau.
                ! If there are multiple particles, scale the probability.
                
                rat(:) = fac(run) * abs(realwSign(min_part_type(run):max_part_type(run)))

                ndie(min_part_type(run):max_part_type(run)) = real(int(rat), dp)
                rat(:) = rat(:) - ndie(min_part_type(run):max_part_type(run))

                ! Choose to die or not stochastically
                r = genrand_real2_dSFMT() 
                if (abs(rat(1)) > r) ndie(min_part_type(run)) = &
                    ndie(min_part_type(run)) + real(nint(sign(1.0_dp, rat(1))), dp)
#ifdef __CMPLX
                r = genrand_real2_dSFMT() 
                if (abs(rat(2)) > r) ndie(max_part_type(run)) = &
                    ndie(max_part_type(run)) + real(nint(sign(1.0_dp, rat(2))), dp)
#endif               
            enddo
        endif

        ! Avoid compiler warnings
        iUnused = DetCurr(1)

    end function

end module
