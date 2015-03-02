#include "macros.h"
module orthogonalise

    use FciMCData, only: TotWalkers, CurrentDets, all_norm_psi_squared, &
                         NoBorn, NoDied, fcimc_iter_data
    use dSFMT_interface, only: genrand_real2_dSFMT
    use bit_reps, only: extract_sign, encode_sign
    use CalcData, only: OccupiedThresh
    use Parallel_neci
    use constants
    implicit none

contains

    subroutine orthogonalise_replicas (iter_data)

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.
        !
        ! |psi_2'> = |psi_2> - (|psi_1><psi_1|psi_2>)/(<psi_1|psi_1>)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: tgt_run, src_run, run, j
        real(dp) :: norms(inum_runs), overlaps(inum_runs, inum_runs)
        real(dp) :: all_norms(inum_runs), all_overlaps(inum_runs, inum_runs)
        real(dp) :: sgn(lenof_sign), sgn_orig, delta, r

        ASSERT(inum_runs == lenof_sign)
#ifndef __PROG_NUMRUNS
        call stop_all(this_routine, "orthogonalise replicas requires mneci.x")
#else

        norms = 0
        overlaps = 0
        do tgt_run = 1, inum_runs

            do j = 1, int(TotWalkers, sizeof_int)

                ! n.b. We are using a non-contiguous list (Hash algorith)
                call extract_sign(CurrentDets(:,j), sgn)
                if (IsUnoccDet(sgn)) cycle

                ! Loop over source runs, and subtract out components
                sgn_orig = sgn(tgt_run)
                do src_run = 1, tgt_run - 1
                    delta = - sgn(src_run) * all_overlaps(src_run, tgt_run) &
                                           / all_norms(src_run)
                    sgn(tgt_run) = sgn(tgt_run) + delta
                end do

                ! TODO: Consider InitiatorOccupiedThresh?
                if (abs(sgn(tgt_run)) < OccupiedThresh) then
                    r = genrand_real2_dSFMT()
                    if (r > abs(sgn(tgt_run)) / OccupiedThresh) then
                        sgn(tgt_run) = sign(OccupiedThresh, sgn(tgt_run))
                    else
                        sgn(tgt_run) = 0.0_dp
                    end if
                end if

                call encode_sign(CurrentDets(:,j), sgn)

                ! Note that we shouldn't be able to kill all particles on a
                ! site, as we can only change a run if there are particles in
                ! a lower indexed run to change...
                ASSERT(.not. IsUnoccDet(sgn))

                ! Now we need to to our accounting to make sure that NoBorn/
                ! Died/etc. counters remain reasonable.
                !
                ! n.b. we don't worry about the delta=0 case, as adding 0 to
                !      the accumulators doesn't cause errors...
                if (sgn(tgt_run) >= 0 .eqv. sgn_orig >= 0) then
                    NoDied(tgt_run) = NoDied(tgt_run) + abs(sgn_orig)
                    iter_data%ndied(tgt_run) = iter_data%ndied(tgt_run) &
                                             + abs(sgn_orig)
                    NoBorn(tgt_run) = NoBorn(tgt_run) + abs(sgn(tgt_run))
                    iter_data%nborn(tgt_run) = iter_data%nborn(tgt_run) &
                                             + abs(sgn(tgt_run))
                else if (abs(sgn(tgt_run)) >= abs(sgn_orig)) then
                    NoDied(tgt_run) = NoBorn(tgt_run) &
                                    + abs(sgn(tgt_run) - sgn_orig)
                    iter_data%ndied(tgt_run) = iter_data%ndied(tgt_run) &
                                             + abs(sgn(tgt_run) - sgn_orig)
                else
                    NoBorn(tgt_run) = NoBorn(tgt_run) &
                                    + abs(sgn(tgt_run) - sgn_orig)
                    iter_data%nborn(tgt_run) = iter_data%nborn(tgt_run) &
                                             + abs(sgn(tgt_run) - sgn_orig)
                end if

                ! We should be accumulating the norm for this run, so it can
                ! be used in processing later runs. We should also be
                ! accumulating the overlap terms for doing the same
                !
                ! n.b. These are the values _after_ orthogonalisation with
                !      the previous runs
                norms(tgt_run) = norms(tgt_run) + sgn(tgt_run)**2
                do run = tgt_run + 1, inum_runs
                    overlaps(tgt_run, run) = overlaps(tgt_run, run) &
                                           * sgn(tgt_run) * sgn(run)
                    overlaps(run, tgt_run) = 99999999.0_dp ! invalid
                end do

            end do

            ! And ensure that the norm/overlap data is accumulated onto all
            ! of the processors
            call MPISumAll(norms, all_norms)
            call MPISumAll(overlaps, all_overlaps)

        end do
#endif

    end subroutine

    subroutine orthogonalise_replicas_2runs (iter_data)

        ! Apply a Gram Schmidt orthogonalisation to the different system
        ! replicas.
        !
        ! |psi_2'> = |psi_2> - (|psi_1><psi_1|psi_2>)/(<psi_1|psi_1>)

        integer :: j
        real(dp) :: sgn(lenof_sign), delta, scal_prod, all_scal_prod, r
        real(dp) :: sgn_orig
        real(dp) :: psi_squared(lenof_sign), all_psi_squared(lenof_sign)
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'orthogonalise_replicas'

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
            ASSERT(.not. IsUnoccDet(sgn))

            ! Now we need to do our accounting to make sure that the NoBorn/
            ! Died/etc. counters remain reasonable.
            !
            ! n.b. we don't worry about the delta=0 case, as adding 0 to the
            !      accumulators doesn't cause errors.
            if (sgn(2) >= 0 .eqv. sgn_orig >= 0) then
                NoDied(2) = NoDied(2) + abs(sgn_orig)
                iter_data%ndied(2) = iter_data%ndied(2) + abs(sgn_orig)
                NoBorn(2) = NoBorn(2) + abs(sgn(2))
                iter_data%nborn(2) = iter_data%nborn(2) + abs(sgn(2))
            else if (abs(sgn(2)) >= abs(sgn_orig)) then
                NoDied(2) = NoBorn(2) + abs(sgn(2) - sgn_orig)
                iter_data%ndied(2) = iter_data%ndied(2) + abs(sgn(2) -sgn_orig)
            else
                NoBorn(2) = NoBorn(2) + abs(sgn(2) - sgn_orig)
                iter_data%nborn(2) = iter_data%nborn(2) + abs(sgn(2) -sgn_orig)
            end if


        end do

#endif

    end subroutine

end module
