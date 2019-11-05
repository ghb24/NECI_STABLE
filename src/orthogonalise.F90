#include "macros.h"
module orthogonalise

    use FciMCData, only: TotWalkers, CurrentDets, all_norm_psi_squared, &
                         NoBorn, NoDied, fcimc_iter_data,  replica_overlaps_real, &
#ifdef __CMPLX
                         replica_overlaps_imag, &
#endif
                         HolesInList, iter
    use CalcData, only: OccupiedThresh, tOrthogonaliseSymmetric, tSemiStochastic, &
                        tPairedReplicas, t_test_overlap, overlap_eps, n_stop_ortho
    use dSFMT_interface, only: genrand_real2_dSFMT
    use load_balance, only: CalcHashTableStats
    use bit_reps, only: extract_sign, encode_sign
    use semi_stoch_procs, only: check_determ_flag
    use Parallel_neci
    use constants
    use util_mod
    implicit none

contains

    subroutine orthogonalise_replicas (iter_data)

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.
        !
        ! |psi_2'> = |psi_2> - (|psi_1><psi_1|psi_2>)/(<psi_1|psi_1>)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: tgt_run, src_run, run, j, TotWalkersNew
        integer :: tgt_sgn, src_sgn
        real(dp) :: all_overlaps_real(inum_runs, inum_runs)
        integer :: low_loop_bound, high_loop_bound, loop_step
        real(dp) :: norms(inum_runs), overlaps_real(inum_runs, inum_runs)
        real(dp) :: sgn(lenof_sign), sgn_orig, delta_real, r
#ifdef __CMPLX
        ! the min_part_type and max_part_type macros return the index of the sgns
        ! corresponding to the real and imag parts of a replica respectively
        real(dp) :: overlaps_imag(inum_runs, inum_runs), all_overlaps_imag(inum_runs, inum_runs)
        real(dp) :: delta_imag
#endif
        logical :: tCoreDet
        character(*), parameter :: this_routine = 'orthogonalise_replicas'

#ifndef __PROG_NUMRUNS
        call stop_all(this_routine, "orthogonalise replicas requires mneci or kmneci")
#endif

        ! If we are using the symmetric orthogonaliser, then bypass this
        ! routine...
        if (tOrthogonaliseSymmetric) then
            call orthogonalise_replicas_lowdin(iter_data)
            return
        end if

        norms = 0.0_dp
        overlaps_real = 0.0_dp
#ifdef __CMPLX
        overlaps_imag = 0.0_dp
#endif

        do tgt_run = 1, inum_runs

            HolesInList = 0
            do j = 1, int(TotWalkers, sizeof_int)

                ! n.b. We are using a non-contiguous list (Hash algorithm)
                call extract_sign(CurrentDets(:,j), sgn)
                tCoreDet = check_determ_flag(CurrentDets(:,j))
                if (IsUnoccDet(sgn) .and. (.not. tCoreDet)) then
                    HolesInList = HolesInList + 1
                    cycle
                end if

                ! Loop over source runs, and subtract out components
                sgn_orig = sgn(tgt_run)

                ! set loop parameters
                if (tPairedReplicas) then
                    low_loop_bound = 2-mod(tgt_run, 2)
                    high_loop_bound = tgt_run-2
                    loop_step = 2
                else
                    low_loop_bound = 1
                    high_loop_bound = tgt_run-1
                    loop_step = 1
                endif

                do src_run = low_loop_bound, high_loop_bound, loop_step
#ifdef __CMPLX
                    ! (a + ib)*(c + id) = ac - bd + i(ad +bc)
                    delta_real = - ( sgn(min_part_type(src_run)) &
                                    * all_overlaps_real(src_run, tgt_run) &
                                 - sgn(max_part_type(src_run)) &
                                    * all_overlaps_imag(src_run, tgt_run)) &
                                           / all_norms(src_run)
                    sgn(min_part_type(tgt_run)) = sgn(min_part_type(tgt_run)) + delta_real

                    delta_imag = - ( sgn(min_part_type(src_run)) &
                                    * all_overlaps_imag(src_run, tgt_run) &
                                 + sgn(max_part_type(src_run)) &
                                    * all_overlaps_real(src_run, tgt_run)) &
                                           / all_norms(src_run)
                    sgn(max_part_type(tgt_run)) = sgn(max_part_type(tgt_run)) + delta_imag

#else
                    delta_real = - sgn(src_run) * all_overlaps_real(src_run, tgt_run) &
                                           / all_norms(src_run)
                    ! test if a small remaining overlap still gives the 
                    ! correct shift.. 
                    if (t_test_overlap .and. tgt_run == inum_runs) then
                        ! draw a random sign and change by the 
                        ! specified epsilon
                        ! or choose just a random number out of [-eps,+eps]
!                         delta_real = delta_real + rand_sign * overlap_eps
!                         delta_real = delta_real + (rand()-0.5_dp)*overlap_eps
                        ! or try a prefactor..
                        delta_real = overlap_eps*delta_real
                        if (n_stop_ortho > 0 .and. iter > n_stop_ortho) then
                            delta_real  = 0.0_dp
                        end if
                    end if

                    sgn(tgt_run) = sgn(tgt_run) + delta_real
#endif
                end do

                call encode_sign(CurrentDets(:,j), sgn)

                ! Note that we shouldn't be able to kill all particles on a
                ! site, as we can only change a run if there are particles in
                ! a lower indexed run to change...
                ! The exception is when using semi-stochastic, where
                ! unoccupied determinants can be stored.
                if (.not. tSemiStochastic) then
                    ASSERT(.not. IsUnoccDet(sgn))
                end if

                ! Now we need to our accounting to make sure that NoBorn/
                ! Died/etc. counters remain reasonable.
                !
                ! n.b. we don't worry about the delta=0 case, as adding 0 to
                !      the accumulators doesn't cause errors...
                if (sgn(tgt_run) >= 0 .neqv. sgn_orig >= 0) then
                    NoDied(tgt_run) = NoDied(tgt_run) + abs(sgn_orig)
                    iter_data%ndied(tgt_run) = iter_data%ndied(tgt_run) &
                                             + abs(sgn_orig)
                    NoBorn(tgt_run) = NoBorn(tgt_run) + abs(sgn(tgt_run))
                    iter_data%nborn(tgt_run) = iter_data%nborn(tgt_run) &
                                             + abs(sgn(tgt_run))
                else if (abs(sgn(tgt_run)) >= abs(sgn_orig)) then
                    NoBorn(tgt_run) = NoBorn(tgt_run) &
                                    + abs(sgn(tgt_run) - sgn_orig)
                    iter_data%nborn(tgt_run) = iter_data%nborn(tgt_run) &
                                             + abs(sgn(tgt_run) - sgn_orig)
                else
                    NoDied(tgt_run) = NoDied(tgt_run) &
                                    + abs(sgn(tgt_run) - sgn_orig)
                    iter_data%ndied(tgt_run) = iter_data%ndied(tgt_run) &
                                             + abs(sgn(tgt_run) - sgn_orig)
                end if

                ! We should be accumulating the norm for this run, so it can
                ! be used in processing later runs. We should also be
                ! accumulating the overlap terms for doing the same
                !
                ! n.b. These are the values _after_ orthogonalisation with
                !      the previous runs
#ifdef __CMPLX
                ! (a+ib)*(a-ib) = a^2 + b^2
                norms(tgt_run) = norms(tgt_run) &
                    + sgn(min_part_type(tgt_run))**2 + sgn(max_part_type(tgt_run))**2
#else
                norms(tgt_run) = norms(tgt_run) + sgn(tgt_run)**2
#endif
                do run = tgt_run + 1, inum_runs
#ifdef __CMPLX
                    ! we want the inner product of the vectors:
                    !        ___
                    ! <n| =  \  ` [ a(n,r) - ib(n,r) ] <r|
                    !        /__,
                    !          r
                    !        ___
                    ! |m> =  \  ` [ a(n,r) + ib(n,r) ] |r>
                    !        /__,
                    !          r
                    !          ___
                    ! <n|m> =  \  ` a(n,r)*a(m,r) + b(n,r)*b(m,r)
                    !          /__,   + i [ a(n,r)*b(m,r) - b(n,r)*a(m,r) ]
                    !            r
                    !
                    ! Note that symmetry may no longer be assumed for the overlap matrix
                    ! i.e. overlaps_real(n,m) =  overlaps_real(m,n), but
                    !      overlaps_imag(n,m) = -overlaps_imag(m,n)
                    ! only storing overlaps for run(m) > tgt_run(n):
                    ! <n|m> =
                    !  _             _
                    ! |               |
                    ! | x             |
                    ! | x x           |
                    ! | x x x         |
                    ! | x x x x       |
                    ! | x x x x x     |
                    ! |_x x x x x x  _|

                    overlaps_real(tgt_run, run) = overlaps_real(tgt_run, run) &
                                           + sgn(min_part_type(run))*sgn(min_part_type(tgt_run)) &
                                           + sgn(max_part_type(run))*sgn(max_part_type(tgt_run))
                    overlaps_imag(tgt_run, run) = overlaps_imag(tgt_run, run) &
                                           + sgn(min_part_type(run))*sgn(max_part_type(tgt_run)) &
                                           + sgn(max_part_type(run))*sgn(min_part_type(tgt_run))
                    overlaps_real(run, tgt_run) = overlaps_real(tgt_run, run)
                    overlaps_imag(run, tgt_run) = -overlaps_imag(tgt_run, run)
#else
                    overlaps_real(tgt_run, run) = overlaps_real(tgt_run, run) &
                                           + sgn(run)*sgn(tgt_run)
                    overlaps_real(run, tgt_run) = overlaps_real(tgt_run, run)
#endif
                end do

            end do

            ! And ensure that the norm/overlap data is accumulated onto all
            ! of the processors.
            call MPISumAll(norms, all_norms)
            call MPISumAll(overlaps_real, all_overlaps_real)
#ifdef __CMPLX
            call MPISumAll(overlaps_imag, all_overlaps_imag)
#endif

        end do

        ! Store a normalised overlap matrix for each of the replicas.
        do src_run = 1, inum_runs - 1
            do tgt_run = src_run + 1, inum_runs
               if(all_norms(src_run)*all_norms(tgt_run) > EPS) then
                replica_overlaps_real(src_run, tgt_run) = &
                    all_overlaps_real(src_run, tgt_run) / &
                    sqrt(all_norms(src_run) * all_norms(tgt_run))
#ifdef __CMPLX
                replica_overlaps_imag(src_run, tgt_run) = &
                    all_overlaps_imag(src_run, tgt_run) / &
                    sqrt(all_norms(src_run) * all_norms(tgt_run))
#endif
             endif
            end do
        end do

        ! Store in global data for perturbation accumulation.
        all_overlaps = all_overlaps_real

        ! We now need to aggregate statistics here, rather than at the end
        ! of annihilation, as they have been modified by this routine...
        !
        ! n.b. TotWalkersNew is integer, TotWalkers is int64 to ensure that it
        !      can be collected using AllTotWalkers (which may end up with far
        !      too many walkers).
        TotWalkersNew = int(TotWalkers)
        call CalcHashTableStats(TotWalkersNew, iter_data)
        TotWalkers = TotWalkersNew

    end subroutine orthogonalise_replicas


    subroutine orthogonalise_replica_pairs (iter_data)

        ! complex walkers not supported here. This routine is never called: deprecated code?

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.

        ! This routine works in cases where there are two replicas representing
        ! each excited state. Whenever an overlap has to be calculated between
        ! different excited states, both pairs are used and averaged. When a
        ! norm needs to be calculated, the overlap is taken between the two
        ! replicas for the same excited state.

        ! Underscores denote the excited state label and uperscores represent
        ! the replica label (i.e. 1 or 2):

        ! |psi_2^1'> = |psi_2^1> - (|psi_1^1><psi_1^2|psi_2^1> +
        !                           |psi_1^2><psi_1^1|psi_2^1>)/(2<psi_1^1|psi_1^2>)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: tgt_state, src_state, run, j, irep, imod1, imod2, TotWalkersNew
        real(dp) :: norms(inum_runs .div. 2), overlaps(inum_runs, inum_runs)
        real(dp) :: sgn(lenof_sign), sgn_orig(2), delta, r
        logical :: tCoreDet
        character(len=*), parameter :: this_routine = "orthogonalise_replica_pairs"

        ASSERT(inum_runs == lenof_sign)
#ifndef __PROG_NUMRUNS
        call unused(iter_data % update_iters)
        call stop_all(this_routine, "orthogonalise replicas requires mneci.x")
#else

        norms = 0.0_dp
        overlaps = 0.0_dp

        ! Loop over all excited states (thus over all replica *pairs*).
        do tgt_state = 1, inum_runs/2

            HolesInList = 0
            do j = 1, int(TotWalkers, sizeof_int)

                ! n.b. We are using a non-contiguous list (Hash algorithm).
                call extract_sign(CurrentDets(:,j), sgn)
                tCoreDet = check_determ_flag(CurrentDets(:,j))
                if (IsUnoccDet(sgn) .and. (.not. tCoreDet)) then
                    HolesInList = HolesInList + 1
                    cycle
                end if

                if (tgt_state > 1) then

                ! The two replica signs for this excited state.
                sgn_orig(1) = sgn(tgt_state*2-1)
                sgn_orig(2) = sgn(tgt_state*2)

                ! Loop over both replicas for this excited state.
                do irep = 1, 2

                    ! These variables equal (1,2) and (2,1) when irep = 1 and
                    ! 2, repsectively. They are used to access the correct
                    ! indices of the sign arrays.
                    imod1 = mod(irep, 2)
                    imod2 = 1 - imod1

                    do src_state = 1, tgt_state - 1
                        delta = - sgn(src_state*2-imod1) * &
                                  all_overlaps(src_state*2-imod2, tgt_state*2-imod1) / &
                                  (2 * all_norms(src_state))
                        sgn(tgt_state*2-imod1) = sgn(tgt_state*2-imod1) + delta

                        delta = - sgn(src_state*2-imod2) * &
                                  all_overlaps(src_state*2-imod1, tgt_state*2-imod1) / &
                                  (2 * all_norms(src_state))
                        sgn(tgt_state*2-imod1) = sgn(tgt_state*2-imod1) + delta
                    end do

                end do

                call encode_sign(CurrentDets(:,j), sgn)

                ! Note that we shouldn't be able to kill all particles on a
                ! site, as we can only change a run if there are particles in
                ! a lower indexed run to change...
                ! The exception is when using semi-stochastic, where
                ! unoccupied determinants can be stored.
                if (.not. tSemiStochastic) then
                    ASSERT(.not. IsUnoccDet(sgn))
                end if

                ! Now we need to our accounting to make sure that NoBorn/
                ! Died/etc. counters remain reasonable.
                !
                ! n.b. we don't worry about the delta=0 case, as adding 0 to
                !      the accumulators doesn't cause errors...
                do irep = 1, 2
                    run = tgt_state*2 - mod(irep,2)

                    if (sgn(run) >= 0 .neqv. sgn_orig(irep) >= 0) then
                        NoDied(run) = NoDied(run) + abs(sgn_orig(irep))
                        iter_data%ndied(run) = iter_data%ndied(run) &
                                                 + abs(sgn_orig(irep))
                        NoBorn(run) = NoBorn(run) + abs(sgn(run))
                        iter_data%nborn(run) = iter_data%nborn(run) &
                                                 + abs(sgn(run))
                    else if (abs(sgn(run)) >= abs(sgn_orig(irep))) then
                        NoBorn(run) = NoBorn(run) &
                                        + abs(sgn(run) - sgn_orig(irep))
                        iter_data%nborn(run) = iter_data%nborn(run) &
                                                 + abs(sgn(run) - sgn_orig(irep))
                    else
                        NoDied(run) = NoDied(run) &
                                        + abs(sgn(run) - sgn_orig(irep))
                        iter_data%ndied(run) = iter_data%ndied(run) &
                                                 + abs(sgn(run) - sgn_orig(irep))
                    end if
                end do

                end if

                ! We should be accumulating the norm for this run, so it can
                ! be used in processing later runs. We should also be
                ! accumulating the overlap terms for doing the same
                !
                ! n.b. These are the values _after_ orthogonalisation with
                !      the previous runs
                norms(tgt_state) = norms(tgt_state) + sgn(2*tgt_state-1)*sgn(2*tgt_state)

                ! Calculate overlaps for all 4 combinations of the 2 replicas
                ! of the 2 excited states.
                do run = tgt_state + 1, inum_runs/2
                    overlaps(tgt_state*2-1, run*2-1) = overlaps(tgt_state*2-1, run*2-1) &
                                                     + sgn(tgt_state*2-1) * sgn(run*2-1)

                    overlaps(tgt_state*2-1, run*2) = overlaps(tgt_state*2-1, run*2) &
                                                     + sgn(tgt_state*2-1) * sgn(run*2)

                    overlaps(tgt_state*2, run*2-1) = overlaps(tgt_state*2, run*2-1) &
                                                     + sgn(tgt_state*2) * sgn(run*2-1)

                    overlaps(tgt_state*2, run*2) = overlaps(tgt_state*2, run*2) &
                                                     + sgn(tgt_state*2) * sgn(run*2)
                end do

            end do

            ! And ensure that the norm/overlap data is accumulated onto all
            ! of the processors.
            call MPISumAll(norms, all_norms)
            call MPISumAll(overlaps, all_overlaps)

        end do

        ! We now need to aggregate statistics here, rather than at the end
        ! of annihilation, as they have been modified by this routine...
        !
        ! n.b. TotWalkersNew is integer, TotWalkers is int64 to ensure that it
        !      can be collected using AllTotWalkers (which may end up with far
        !      too many walkers).
        TotWalkersNew = int(TotWalkers)
        call CalcHashTableStats(TotWalkersNew, iter_data)
        TotWalkers = TotWalkersNew

#endif

    end subroutine orthogonalise_replica_pairs

    subroutine orthogonalise_replicas_2runs (iter_data)

        ! complex walkers not supported here. This routine is never called: deprecated code?

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.
        !
        ! |psi_2'> = |psi_2> - (|psi_1><psi_1|psi_2>)/(<psi_1|psi_1>)

        integer :: j, TotWalkersNew
        real(dp) :: sgn(lenof_sign), delta, scal_prod, all_scal_prod, r
        real(dp) :: sgn_orig
        real(dp) :: psi_squared(lenof_sign), all_psi_squared(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'orthogonalise_replicas_2runs'

        ASSERT(inum_runs == 2)
        ASSERT(lenof_sign == 2)

#ifndef __PROG_NUMRUNS
        call stop_all(this_routine, "orthogonalise replicas requires mneci.x")
#else

        ! We need the norm of the wavefunction to do anything. Don't trust
        ! the global values here, as they aren't valid until after the
        ! communication --> not set yet.
        !
        ! We want them to be valid for the new psi...
        psi_squared = 0
        scal_prod = 0
        do j = 1, int(TotWalkers, sizeof_int)
            call extract_sign(CurrentDets(:, j), sgn)
            if (IsUnoccDet(sgn)) cycle
            psi_squared = psi_squared + sgn**2
            scal_prod = scal_prod + sgn(1) * sgn(2)
        end do
        call MPISumAll(psi_squared, all_psi_squared)
        call MPISumAll(scal_prod, all_scal_prod)

        ! Calculate the change
        do j = 1, int(TotWalkers, sizeof_int)

            ! Adjust the wavefunctions
            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

            delta = - sgn(1) * all_scal_prod / all_psi_squared(1)
            sgn_orig = sgn(2)
            sgn(2) = sgn(2) + delta

            ! And stochastically round, so that the minimum particle sign
            ! is maintained in an unbiased way.
            if (abs(sgn(2)) < 1.0_dp) then
                r = genrand_real2_dSFMT()
                if (r > abs(sgn(2))) then
                    sgn(2) = sign(1.0_dp, sgn(2))
                else
                    sgn(2) = 0.0_dp
                endif
            end if

            call encode_sign(CurrentDets(:,j), sgn)

            ! Note that we shouldn't be able to kill all particles on a site,
            ! as we can only change run 2 if run 1 is occupied, and run 1
            ! doesn't change...
            ! The exception is when using semi-stochastic, where
            ! unoccupied determinants can be stored.
            if (.not. tSemiStochastic) then
                ASSERT(.not. IsUnoccDet(sgn))
            end if

            ! Now we need to do our accounting to make sure that the NoBorn/
            ! Died/etc. counters remain reasonable.
            !
            ! n.b. we don't worry about the delta=0 case, as adding 0 to the
            !      accumulators doesn't cause errors.
            if (sgn(2) >= 0 .neqv. sgn_orig >= 0) then
                NoDied(2) = NoDied(2) + abs(sgn_orig)
                iter_data%ndied(2) = iter_data%ndied(2) + abs(sgn_orig)
                NoBorn(2) = NoBorn(2) + abs(sgn(2))
                iter_data%nborn(2) = iter_data%nborn(2) + abs(sgn(2))
            else if (abs(sgn(2)) >= abs(sgn_orig)) then
                NoBorn(2) = NoBorn(2) + abs(sgn(2) - sgn_orig)
                iter_data%nborn(2) = iter_data%nborn(2) + abs(sgn(2) -sgn_orig)
            else
                NoDied(2) = NoDied(2) + abs(sgn(2) - sgn_orig)
                iter_data%ndied(2) = iter_data%ndied(2) + abs(sgn(2) -sgn_orig)
            end if


        end do

#endif

        ! We now need to aggregate statistics here, rather than at the end
        ! of annihilation, as they have been modified by this routine...
        !
        ! n.b. TotWalkersNew is integer, TotWalkers is int64 to ensure that it
        !      can be collected using AllTotWalkers (which may end up with far
        !      to many walkers).
        TotWalkersNew = int(TotWalkers)
        call CalcHashTableStats(TotWalkersNew, iter_data)
        TotWalkers = TotWalkersNew

    end subroutine

    subroutine write_mat(mat)

        real(dp) :: mat(:,:)
        integer :: i, j

        do i = lbound(mat, 1), ubound(mat, 1)
            do j = lbound(mat, 2), ubound(mat, 2)
                write(6, '(f17.9, " ")', advance='no') mat(i, j)
            end do
            write(6,*)
        end do

    end subroutine


    subroutine orthogonalise_replicas_lowdin(iter_data)

        ! Perform a symmetric (Lowdin) orthogonalisation

        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine='orthogonalise_replicas_lowdin'

        real(dp) :: S(inum_runs, inum_runs), S_all(inum_runs, inum_runs)
        real(dp) :: evecs(inum_runs, inum_runs), evecs_t(inum_runs, inum_runs)
        real(dp) :: S_half(inum_runs, inum_runs)
        real(dp) :: elem, sgn(lenof_sign), sgn_orig(lenof_sign)
        real(dp) :: work(3*inum_runs-1), evals(inum_runs)
        real(dp) :: sgn_orig_norm(lenof_sign), sgn_norm(lenof_sign)
        real(dp) :: norm(inum_runs)
        integer :: j, run, runa, runb, info, TotWalkersNew
        logical :: tCoreDet

        ! Not implemented for complex (yet)
        ASSERT(inum_runs == lenof_sign)
#ifndef __PROG_NUMRUNS
        call stop_all(this_routine, "orthogonalise replicas requires mneci.x")
#else

        ! Generate the overlap matrix (unnormalised)
        S = 0
        do j = 1, int(TotWalkers, sizeof_int)

            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

            do runa = 1, inum_runs
                do runb = runa, inum_runs
                    elem = sgn(runa) * sgn(runb)
                    S(runa, runb) = S(runa, runb) + elem
                    if (runa /= runb) &
                        S(runb, runa) = S(runb, runa) + elem
                end do
            end do
        end do
        call MPISumAll(S, S_all)

        ! Normalise everything (the diagonal terms give the normalisation
        ! constants)
        S = S_all
        do run = 1, inum_runs
            norm(run) = sqrt(S_all(run, run))
            S(run, :) = S(run, :) / norm(run)
            S(:, run) = S(:, run) / norm(run)
        end do
        evecs = S

        ! Diagonalise the S matrix
        call dsyev('V', 'U', inum_runs, evecs, inum_runs, evals, work, &
                   size(work), info)

        if (any(evals < 0)) then
            write(6,*) '*** WARNING ***'
            write(6,*) "Not orthogonalising this iteration."
            write(6,*) 'Negative eigenvalue of overlap matrix found'
            return
        end if

        ! Take the square roots of the eigenvalues
        evals = 1.0_dp / sqrt(evals)
        evecs_t = transpose(evecs)

        ! Multiply through by the square root, and obtain S^{-0.5}
        do j = 1, inum_runs
            evecs_t(j, :) = evecs_t(j, :) * evals(j)
        end do

        S_half = matmul(evecs, evecs_t)

        if (any(isnan_neci(s_half))) &
            call stop_all(this_routine, "NaNs found")

        ! Go through and update the values!
        HolesInList = 0
        do j = 1, int(TotWalkers, sizeof_int)

            ! n.b. We are using a non-contiguous list (Hash algorith)
            call extract_sign(CurrentDets(:,j), sgn_orig)
            tCoreDet = check_determ_flag(CurrentDets(:,j))
            if (IsUnoccDet(sgn_orig) .and. (.not. tCoreDet)) then
                HolesInList = HolesInList + 1
                cycle
            end if

            ! Ultimately, we want to use normalised signs for the Lowdin
            ! expressions
            sgn_orig_norm = sgn_orig / norm

            ! Obtain the new sign values
            sgn_norm = matmul(S_half, sgn_orig_norm)
            sgn = sgn_norm * norm
            call encode_sign(CurrentDets(:,j), sgn)

            ! We should not be able to kill all particles on a site. This is
            ! a rotation. The exception is when using semi-stochastic, where
            ! unoccupied determinants can be stored.
            if (.not. tSemiStochastic) then
                ASSERT(.not. IsUnoccDet(sgn))
            end if

            ! Do some particle accounting
            do run = 1, inum_runs
                if (sgn(run) >= 0 .neqv. sgn_orig(run) >= 0) then
                    NoDied(run) = NoDied(run) + abs(sgn_orig(run))
                    iter_data%ndied(run) = iter_data%ndied(run) &
                                         + abs(sgn_orig(run))
                    NoBorn(run) = NoBorn(run) + abs(sgn(run))
                    iter_data%nborn(run) = iter_data%nborn(run) + abs(sgn(run))
                else if (abs(sgn(run)) >= abs(sgn_orig(run))) then
                    NoBorn(run) = NoBorn(run) + abs(sgn(run) - sgn_orig(run))
                    iter_data%nborn(run) = iter_data%nborn(run) &
                                         + abs(sgn(run) - sgn_orig(run))
                else
                    NoDied(run) = NoDied(run) + abs(sgn(run) - sgn_orig(run))
                    iter_data%ndied(run) = iter_data%ndied(run) &
                                         + abs(sgn(run) - sgn_orig(run))
                end if
            end do

        end do

#endif

        ! We now need to aggregate statistics here, rather than at the end
        ! of annihilation, as they have been modified by this routine...
        !
        ! n.b. TotWalkersNew is integer, TotWalkers is int64 to ensure that it
        !      can be collected using AllTotWalkers (which may end up with far
        !      to many walkers).
        TotWalkersNew = int(TotWalkers)
        call CalcHashTableStats(TotWalkersNew, iter_data)
        TotWalkers = TotWalkersNew

    end subroutine

    subroutine calc_replica_overlaps()

        ! A routine for just calculating the overlap, in cases where
        ! orthogonalisation is not being performed.

        integer :: j, run, tgt_run, src_run
        real(dp) :: sgn(lenof_sign)
        real(dp) :: norms(inum_runs), overlaps(inum_runs, inum_runs)
        character(*), parameter :: this_routine = 'calc_replica_overlaps'

        norms = 0.0_dp
        overlaps = 0.0_dp
        do j = 1, int(TotWalkers, sizeof_int)

            ! n.b. We are using a non-contiguous list (Hash algorithm)
            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

#ifndef __CMPLX
            norms = norms + sgn*sgn
#endif

            do tgt_run = 1, inum_runs
                do run = tgt_run + 1, inum_runs
                    overlaps(tgt_run, run) = overlaps(tgt_run, run) &
                                           + sgn(tgt_run) * sgn(run)
                    overlaps(run, tgt_run) = 99999999.0_dp ! invalid
                end do
            end do

        end do

        ! And ensure that the norm/overlap data is accumulated onto all
        ! of the processors.
        call MPISumAll(norms, all_norms)
        call MPISumAll(overlaps, all_overlaps)

        ! Store a normalised overlap matrix for each of the replicas.
        do src_run = 1, inum_runs - 1
            do tgt_run = src_run + 1, inum_runs
               if(all_norms(src_run)*all_norms(tgt_run) > EPS) then
                  replica_overlaps_real(src_run, tgt_run) = &
                       all_overlaps(src_run, tgt_run) / &
                       sqrt(all_norms(src_run) * all_norms(tgt_run))
               endif
                replica_overlaps_real(src_run, tgt_run) = &
                    replica_overlaps_real(src_run, tgt_run)
            end do
        end do

    end subroutine calc_replica_overlaps

end module
