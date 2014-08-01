#include "macros.h"

module tau_search

    use SystemData, only: AB_elec_pairs, par_elec_pairs, tGen_4ind_weighted, &
                          tHPHF, tCSF, tKpntSym, tMomInv, nel, G1, nbasis, &
                          AB_hole_pairs, par_hole_pairs, tGen_4ind_reverse
    use CalcData, only: tTruncInitiator, tReadPops, MaxWalkerBloom, tau, &
                        InitiatorWalkNo, tWalkContGrow
    use FciMCData, only: tRestart, pSingles, pDoubles, pParallel, &
                         ProjEDet, ilutRef
    use GenRandSymExcitNUMod, only: construct_class_counts, &
                                    init_excit_gen_store, clean_excit_gen_store
    use SymExcit3, only: GenExcitations3
    use Determinants, only: get_helement
    use HPHFRandExcitMod, only: ReturnAlphaOpenDet, CalcPGenHPHF, &
                                CalcNonUniPGen
    use HPHF_integrals, only: hphf_off_diag_helement_norm
    use SymExcitDataMod, only: excit_gen_store_type
    use bit_rep_data, only: NIfTot
    use DetBitOps, only: FindBitExcitLevel, TestClosedShellDet, &
                         EncodeBitDet
    use Parallel_neci
    use constants
    implicit none

    real(dp) :: gamma_sing, gamma_doub, gamma_opp, gamma_par, max_death_cpt
    real(dp) :: max_permitted_spawn
    integer :: cnt_sing, cnt_doub, cnt_opp, cnt_par
    integer :: n_opp, n_par
    logical :: enough_sing, enough_doub, enough_opp, enough_par
    logical :: consider_par_bias

contains

    subroutine init_tau_search ()

        ! N.B. This must be called BEFORE a popsfile is read in, otherwise
        !      we screw up the gamma values that have been carefully read in.

        ! We want to start off with zero-values
        gamma_sing = 0
        gamma_doub = 0
        gamma_opp = 0
        gamma_par = 0

        ! And what is the maximum death-component found
        max_death_cpt = 0

        ! And the counts are used to make sure we don't update anything too
        ! early
        cnt_sing = 0
        cnt_doub = 0
        cnt_opp = 0
        cnt_par = 0
        enough_sing = .false.
        enough_doub = .false.
        enough_opp = .false.
        enough_par = .false.

        ! Unless it is already specified, set an initial value for tau
        if (.not. tRestart .and. .not. tReadPops .and. tau == 0) &
            call FindMaxTauDoubs()
        write(6,*) 'Using initial time-step: ', tau

        
        ! Set the maximum spawn size
        if (MaxWalkerBloom == -1) then
            ! No maximum manually specified, so we set the limit of spawn
            ! size to either the initiator criterion, or to 5 otherwise
            if (tTruncInitiator) then
                max_permitted_spawn = InitiatorWalkNo
            else
                max_permitted_spawn = 5.0_dp
            end if
        else
            ! This is specified manually
            max_permitted_spawn = real(MaxWalkerBloom, dp)
        end if

        if (.not. (tReadPops .and. .not. tWalkContGrow)) then
            write(iout, "(a,f10.5)") "Will dynamically update timestep to &
                         &limit spawning probability to", max_permitted_spawn
        end if

        ! Are we considering parallel-spin bias in the doubles?
        ! Do this logic here, so that if we add opposite spin bias to more
        ! excitation generators, then there is only one place that this logic
        ! needs to be updated!
        if (tGen_4ind_weighted) then
            consider_par_bias = .true.
            n_opp = AB_elec_pairs
            n_par = par_elec_pairs
        else if (tGen_4ind_reverse) then
            consider_par_bias = .true.
            n_opp = AB_hole_pairs
            n_par = par_hole_pairs
        else
            consider_par_bias = .false.
        end if

    end subroutine

    subroutine log_spawn_magnitude (ic, ex, matel, prob)

        integer, intent(in) :: ic, ex(2,2)
        real(dp), intent(in) :: prob, matel
        real(dp) :: tmp_gamma, tmp_prob
        integer, parameter :: cnt_threshold = 50

        if (ic == 1) then

            ! Log the details if necessary!
            tmp_prob = prob / pSingles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_sing) &
                gamma_sing = tmp_gamma
            
            ! And keep count!
            if (.not. enough_sing) then
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
            end if

        else
            
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoubles

            ! We need to deal with the doubles
            if (consider_par_bias) then

                ! In this case, distinguish between parallal and oppisite spins
                if (is_beta(ex(1,1)) .eqv. is_beta(ex(1,2))) then
                    tmp_prob = tmp_prob / pParallel
                    tmp_gamma = abs(matel) / tmp_prob
                    if (tmp_gamma > gamma_par) &
                        gamma_par = tmp_gamma

                    ! And keep count
                    if (.not. enough_par) then
                        cnt_par = cnt_par + 1
                        if (cnt_par > cnt_threshold) enough_par = .true.
                        if (enough_opp .and. enough_par) enough_doub = .true.
                    end if
                else
                    tmp_prob = tmp_prob / (1.0_dp - pParallel)
                    tmp_gamma = abs(matel) / tmp_prob
                    if (tmp_gamma > gamma_opp) &
                        gamma_opp = tmp_gamma

                    ! And keep count
                    if (.not. enough_opp) then
                        cnt_opp = cnt_opp + 1
                        if (cnt_opp > cnt_threshold) enough_opp = .true.
                        if (enough_opp .and. enough_par) enough_doub = .true.
                    end if
                end if
            else
                ! We are not playing around with the same/opposite spin bias
                ! then we should just treat doubles like the singles
                tmp_gamma = abs(matel) / tmp_prob
                if (tmp_gamma > gamma_doub) &
                    gamma_doub = tmp_gamma

                ! And keep count
                if (.not. enough_doub) then
                    cnt_doub = cnt_doub + 1
                    if (cnt_doub > cnt_threshold) enough_doub = .true.
                end if
            end if

        end if

    end subroutine

    subroutine log_death_magnitude (mult)

        ! The same as above, but for particle death

        real(dp) :: mult

        if (mult > max_death_cpt) &
            max_death_cpt = mult

    end subroutine

    subroutine update_tau ()

        real(dp) :: psingles_new, tau_new, mpi_tmp, tau_death, pParallel_new
        logical :: mpi_ltmp

        ! What needs doing depends on the number of parametrs that are being
        ! updated.
        if (consider_par_bias) then

            ! Considering two types of double exctitaion...
            call MPIAllReduce (gamma_sing, MPI_MAX, mpi_tmp)
            gamma_sing = mpi_tmp
            call MPIAllReduce (gamma_opp, MPI_MAX, mpi_tmp)
            gamma_opp = mpi_tmp
            call MPIAllReduce (gamma_par, MPI_MAX, mpi_tmp)
            gamma_par = mpi_tmp
            call MPIAllReduce (enough_opp, MPI_LOR, mpi_ltmp)
            enough_opp = mpi_ltmp
            call MPIAllReduce (enough_par, MPI_LOR, mpi_ltmp)
            enough_par = mpi_ltmp

            if (enough_sing .and. enough_doub) then
                pparallel_new = gamma_par / (gamma_opp + gamma_par)
                psingles_new = gamma_sing * pparallel_new &
                             / (gamma_par + gamma_sing * pparallel_new)
                tau_new = psingles_new * max_permitted_spawn &
                              / gamma_sing
            else
                pparallel_new = pParallel
                psingles_new = pSingles
                tau_new = max_permitted_spawn * &
                        min(pSingles / gamma_sing, &
                        min(pDoubles * pParallel / gamma_par, &
                            pDoubles * pParallel / gamma_opp))
            end if

            ! We only want to update the opposite spins bias here, as we only
            ! consider it here!
            if (enough_opp .and. enough_par) then
                if (abs(pParallel_new-pParallel) / pParallel > 0.0001_dp) then
                    root_print "Updating parallel-spin bias; new pParallel = ", &
                        pParallel_new
                end if
                pParallel = pParallel_new
            end if

        else

            ! Only considering a direct singles/doubles bias
            call MPIAllReduce (gamma_sing, MPI_MAX, mpi_tmp)
            gamma_sing = mpi_tmp
            call MPIAllReduce (gamma_doub, MPI_MAX, mpi_tmp)
            gamma_doub = mpi_tmp

            ! Get the values of pSingles and tau that correspond to the stored
            ! values
            if (enough_sing .and. enough_doub) then
                psingles_new = gamma_sing / (gamma_doub + gamma_sing)
                tau_new = max_permitted_spawn / (gamma_doub + gamma_sing)
            else
                psingles_new = pSingles
                tau_new = max_permitted_spawn * &
                            min(pSingles / gamma_sing, pDoubles / gamma_doub)
            end if

        end if
        call MPIAllReduce (enough_sing, MPI_LOR, mpi_ltmp)
        enough_sing = mpi_ltmp
        call MPIAllReduce (enough_doub, MPI_LOR, mpi_ltmp)
        enough_doub = mpi_ltmp

        ! The range of tau is restricted by particle death. It MUST be <=
        ! the value obtained to restrict the maximum death-factor to 1.0.
        call MPIAllReduce (max_death_cpt, MPI_MAX, mpi_tmp)
        max_death_cpt = mpi_tmp
        tau_death = 1.0_dp / max_death_cpt
        if (tau_death < tau_new) &
            tau_new = tau_death

        ! If the calculated tau is less than the current tau, we should ALWAYS
        ! update it. Once we have a reasonable sample of excitations, then we
        ! can permit tau to increase if we have started too low.
        if (tau_new < tau .or. (enough_sing .and. enough_doub)) then
            if (abs(tau - tau_new) / tau > 0.001_dp) then
                root_print "Updating time-step. New time-step = ", tau_new
            end if
            tau = tau_new
        end if

        ! Make sure that we have at least some of both singles and doubles
        ! before we allow ourselves to change the probabilities too much...
        if (enough_sing .and. enough_doub .and. psingles_new > 1e-5_dp &
            .and. psingles_new < (1.0_dp - 1e-5_dp)) then

            if (abs(psingles - psingles_new) / psingles > 0.0001_dp) then
                root_print "Updating singles/doubles bias. pSingles = ", &
                    psingles_new, ", pDoubles = ", 1.0_dp - psingles_new
            end if
            pSingles = psingles_new
            pDoubles = 1.0_dp - pSingles
        end if

    end subroutine


    subroutine FindMaxTauDoubs()

        ! Routine to find an upper bound to tau, by consideration of the
        ! singles and doubles connected to the reference determinant
        ! 
        ! Obviously, this make assumptions about the possible range of pgen,
        ! so may actually give a tau that is too SMALL for the latest
        ! excitation generators, which is exciting!

        use neci_intfce
        type(excit_gen_store_type) :: store, store2
        logical :: tAllExcitFound,tParity,tSameFunc,tSwapped,tSign
        character(len=*), parameter :: t_r="FindMaxTauDoubs"
        integer :: ex(2,2),ex2(2,2),exflag,iMaxExcit,nStore(6),nExcitMemLen(1)
        integer, allocatable :: Excitgen(:)
        real(dp) :: nAddFac,MagHel,pGen,pGenFac
        HElement_t :: hel
        integer :: ic,nJ(nel),nJ2(nel),ierr,iExcit,ex_saved(2,2)
        integer(kind=n_int) :: iLutnJ(0:niftot),iLutnJ2(0:niftot)

        if(tCSF.or.tMomInv) call stop_all(t_r,"TauSearching needs fixing to work with CSFs or MI funcs")

        if(MaxWalkerBloom.eq.-1) then
            !No MaxWalkerBloom specified
            !Therefore, assume that we do not want blooms larger than n_add if initiator,
            !or 5 if non-initiator calculation.
            if(tTruncInitiator) then
                nAddFac = InitiatorWalkNo
            else
                nAddFac = 5.0_dp    !Won't allow more than 5 particles at a time
            endif
        else
            nAddFac = real(MaxWalkerBloom,dp) !Won't allow more than MaxWalkerBloom particles to spawn in one event. 
        endif

        Tau = 1000.0_dp
        tAllExcitFound=.false.
        Ex_saved(:,:)=0
        exflag=3
        tSameFunc=.false.
        call init_excit_gen_store(store)
        call init_excit_gen_store(store2)
        store%tFilled = .false.
        store2%tFilled = .false.
        CALL construct_class_counts(ProjEDet, store%ClassCountOcc, &
                                    store%ClassCountUnocc)
        store%tFilled = .true.
        if(tKPntSym) then
            !TODO: It REALLY needs to be fixed so that we don't need to do this!!
            !Setting up excitation generators that will work with kpoint sampling
            iMaxExcit=0
            nStore(:)=0
            CALL GenSymExcitIt2(ProjEDet,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
            ALLOCATE(EXCITGEN(nExcitMemLen(1)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(t_r,"Problem allocating excitation generator")
            EXCITGEN(:)=0
            CALL GenSymExcitIt2(ProjEDet,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
        !    CALL GetSymExcitCount(EXCITGEN,DetConn)
        endif

        do while (.not.tAllExcitFound)
            if(tKPntSym) then
                call GenSymExcitIt2(ProjEDet,nel,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
                if(nJ(1).eq.0) exit
                !Calculate ic, tParity and Ex
                call EncodeBitDet (nJ, iLutnJ)
                Ex(:,:)=0
                Ex(1,1)=FindBitExcitlevel(iLutnJ,iLutRef,2)
                call GetExcitation(ProjEDet,nJ,Nel,ex,tParity)
            else
                CALL GenExcitations3(ProjEDet,iLutRef,nJ,exflag,Ex_saved,tParity,tAllExcitFound,.false.)
                IF(tAllExcitFound) EXIT
                Ex(:,:) = Ex_saved(:,:)
                if(Ex(2,2).eq.0) then
                    ic=1
                else
                    ic=2
                endif
                call EncodeBitDet (nJ, iLutnJ)
            endif

            !Find Hij
            if(tHPHF) then
                if(.not.TestClosedShellDet(iLutnJ)) then
                    CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutnJ2,.true.,.true.,tSwapped)
                    if(tSwapped) then
                        !Have to recalculate the excitation matrix.
                        ic = FindBitExcitLevel(iLutnJ, iLutRef, 2)
                        ex(:,:) = 0
                        if(ic.le.2) then
                            ex(1,1) = ic
                            call GetBitExcitation(iLutRef,iLutnJ,Ex,tParity)
                        endif
                    endif
                endif
                hel = hphf_off_diag_helement_norm(ProjEDet,nJ,iLutRef,iLutnJ)
            else
                hel = get_helement(ProjEDet,nJ,ic,ex,tParity)
            endif

            MagHel = abs(hel)

            !Find pGen (nI -> nJ)
            if(tHPHF) then
                call CalcPGenHPHF(ProjEDet,iLutRef,nJ,iLutnJ,ex,store%ClassCountOcc,    &
                            store%ClassCountUnocc,pDoubles,pGen,tSameFunc)
            else
                call CalcNonUnipGen(ProjEDet,ilutRef,ex,ic,store%ClassCountOcc,store%ClassCountUnocc,pDoubles,pGen)
            endif
            if(tSameFunc) cycle
            if(MagHel.gt.0.0_dp) then
                pGenFac = pGen*nAddFac/MagHel
                if(Tau.gt.pGenFac) then
                    Tau = pGenFac
                endif
            endif

            !Find pGen(nJ -> nI)
            CALL construct_class_counts(nJ, store2%ClassCountOcc, &
                                        store2%ClassCountUnocc)
            store2%tFilled = .true.
            if(tHPHF) then
                ic = FindBitExcitLevel(iLutnJ, iLutRef, 2)
                ex2(:,:) = 0
                if(ic.le.2) then
                    ex2(1,1) = ic
                    call GetBitExcitation(iLutnJ,iLutRef,Ex2,tSign)
                endif
                call CalcPGenHPHF(nJ,iLutnJ,ProjEDet,iLutRef,ex2,store2%ClassCountOcc,    &
                            store2%ClassCountUnocc,pDoubles,pGen,tSameFunc)
            else
                ex2(1,:) = ex(2,:)
                ex2(2,:) = ex(1,:)
                call CalcNonUnipGen(nJ,ilutnJ,ex2,ic,store2%ClassCountOcc,store2%ClassCountUnocc,pDoubles,pGen)
            endif
            if(tSameFunc) cycle
            if(MagHel.gt.0.0_dp) then
                pGenFac = pGen*nAddFac/MagHel
                if(Tau.gt.pGenFac) then
                    Tau = pGenFac
                endif
            endif

        enddo
                
        call clean_excit_gen_store (store)
        call clean_excit_gen_store (store2)
        if(tKPntSym) deallocate(EXCITGEN)

        if(tau.gt.0.075_dp) then
            tau=0.075_dp
            write(iout,"(A,F8.5,A)") "Small system. Setting initial timestep to be ",Tau," although this &
                                            &may be inappropriate. Care needed"
        else
            write(iout,"(A,F18.10)") "From analysis of reference determinant and connections, &
                                     &an upper bound for the timestep is: ",Tau
        endif

    end subroutine FindMaxTauDoubs


end module
