#include "macros.h"

module fcimc_pointed_fns

    use SystemData, only: nel, tGen_4ind_2, tGen_4ind_weighted, tHub, tUEG, &
                          tGen_4ind_reverse,  nBasis
    use LoggingData, only: tHistExcitToFrom, FciMCDebug
    use CalcData, only: RealSpawnCutoff, tRealSpawnCutoff, tAllRealCoeff, &
                        RealCoeffExcitThresh, AVMcExcits, tau, DiagSft, &
                        tRealCoeffByExcitLevel, InitiatorWalkNo, &
                        t_fill_frequency_hists, t_truncate_spawns, n_truncate_spawns, & 
                        t_matele_cutoff, matele_cutoff 
    use DetCalcData, only: FciDetIndex, det
    use procedure_pointers, only: get_spawn_helement
    use fcimc_helper, only: CheckAllowedTruncSpawn
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use bit_rep_data, only: NIfTot, test_flag
    use bit_reps, only: get_initiator_flag
    use tau_search, only: log_death_magnitude, log_spawn_magnitude
    use rdm_general, only: calc_rdmbiasfac
    use hist, only: add_hist_excit_tofrom
    use searching, only: BinSearchParts2
    use util_mod
    use FciMCData
    use constants
#ifdef __REALTIME
    use real_time_data, only: NoBorn_1, SpawnFromSing_1, bloom_count_1, &
                              bloom_sizes_1, runge_kutta_step
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

        ! The following is useful for debugging the contributions of single
        ! excitations, and double excitations of spin-paired/opposite
        ! electron pairs to the value of tau.
!        if (ic == 2) then
!            if (G1(ex(1,1))%Ms /= G1(ex(1,2))%Ms) then
!                write(6,*) 'OPP', rh, prob
!            else
!                write(6,*) 'SAM', rh, prob
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
            TargetExcitLevel = FindBitExcitLevel (iLutRef(:,1,1), iLutnJ)
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
            if (mod(part_type,2) == tgt_cpt) then
                ! even part_type => imag replica =>  map 4->3,4 ; 6->5,6 etc.
                component = 1
             else
                component = 2
            endif
#else
            component = tgt_cpt
            if ((part_type.eq.2).and.(inum_runs.eq.1)) component = 3 - tgt_cpt
#endif

            ! Get the correct part of the matrix element
            walkerweight = sign(1.0_dp, RealwSign(part_type))
            if (mod(component,2) == 1) then
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
                if(mod(part_type,2) == 0) then
                   walkerweight = - walkerweight
                endif
#endif
            end if
#endif
            nSpawn = - tau * MatEl * walkerweight / prob
!            write(66,*) part_type, nspawn, RealSpawnCutoff, RealSpawnCutoff, stochastic_round (nSpawn / RealSpawnCutoff)
!            write(66,*) part_type, nspawn, RealSpawnCutoff, RealSpawnCutoff, stochastic_round (nSpawn / RealSpawnCutoff)
            ! [Werner Dobrautz 4.4.2017:]
            ! apply the spawn truncation, when using histogramming tau-search
            if (t_truncate_spawns .and. abs(nspawn) > n_truncate_spawns) then
                ! TODO: add some additional output if this event happens
!                 write(iout,*) "Truncating spawn magnitude from: ", abs(nspawn), " to ", n_truncate_spawns
                nSpawn = sign(n_truncate_spawns, nspawn)


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
        ! in the real-time fciqmc, it is probably good to keep track of the 
        ! stats of the 2 distinct RK loops .. use a global variable for the step
        ! RT_M_Merge: Merge conflict with master due to introduction of kmneci resolved. 
        ! For rmneci, the __REALTIME case will have to be adapted to inum_runs>1

        ! rmneci_setup: Added multirun support for real-time case
#if defined(__REALTIME) 
        do run = 1, inum_runs
           if (runge_kutta_step == 1) then
              NoBorn_1(run) = NoBorn_1(run) + sum(abs(child))
              if (ic == 1) SpawnFromSing_1(run) = SpawnFromSing_1(run) + sum(abs(child))

              if (sum(abs(child)) > InitiatorWalkNo) then
                 bloom_count_1(ic) = bloom_count_1(ic) + 1
                 bloom_sizes_1(ic) = max(real(sum(abs(child)),dp), bloom_sizes_1(ic))
              end if

           else if (runge_kutta_step == 2) then
              NoBorn(run) = NoBorn(run) + sum(abs(child))
              if (ic == 1) SpawnFromSing(run) = SpawnFromSing(run) + sum(abs(child))

              ! Count particle blooms, and their sources
              if (sum(abs(child)) > InitiatorWalkNo) then
                 bloom_count(ic) = bloom_count(ic) + 1
                 bloom_sizes(ic) = max(real(sum(abs(child)), dp), bloom_sizes(ic))
              end if
           end if
        enddo
#elif defined( __CMPLX) && !defined(__REALTIME)
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
