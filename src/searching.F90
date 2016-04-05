! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include  "macros.h"

module searching

    use bit_rep_data, only: nIfDBO, NIfTot, flag_trial, flag_connected, test_flag
    use bit_reps, only: decode_bit_det, set_flag
    use constants
    use DetBitOps, only: DetBitLt, ilut_gt, DetBitEq
    use FciMCData, only: CurrentDets, trial_space, min_trial_ind, trial_space_size, trial_wfs, &
                         con_space, min_conn_ind, con_space_size, con_space_vecs, &
                         trial_numerator, trial_denom, Trial_Search_Time, ntrial_excits
    use hash, only: FindWalkerHash
    use sparse_arrays, only: trial_ht, con_ht
    use SystemData, only: nel
    use timing_neci, only: set_timer, halt_timer
    use util_mod, only: binary_search_custom

    implicit none

contains

    SUBROUTINE LinSearchParts(DetArray,iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER(KIND=n_int) :: iLut(0:NIfTot),DetArray(0:NIfTot,1:MaxInd)
        INTEGER :: N,Comp
        LOGICAL :: tSuccess

        N=MinInd
        do while(N.le.MaxInd)
            Comp=DetBitLT(DetArray(:,N),iLut(:),NIfDBO)
            IF(Comp.eq.1) THEN
                N=N+1
            ELSEIF(Comp.eq.-1) THEN
                PartInd=N-1
                tSuccess=.false.
                RETURN
            ELSE
                tSuccess=.true.
                PartInd=N
                RETURN
            ENDIF
        enddo
        tSuccess=.false.
        PartInd=MaxInd-1

    END SUBROUTINE LinSearchParts

!Do a binary search in CurrentDets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchParts(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF

        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower
!bound of the search up to N.  However, if the lower bound is already equal to N then the two bounds are 
!consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry 
!of the list and leave.  We need to check the last index.
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfDBO)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN

                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts

    subroutine bin_search_trial(ilut, amp, tTrial, tCon)

        integer(n_int), intent(in) :: ilut(0:)
        HElement_t(dp), intent(out) :: amp(:)
        logical, intent(out) :: tTrial, tCon

        integer :: pos

        amp = 0.0_dp
        tTrial = .false.
        tCon = .false.

        ! Search both the trial space and connected space to see if this state exists in either list.
        ! First the trial space:
        pos = binary_search_custom(trial_space(:, min_trial_ind:trial_space_size), ilut, NIfTot+1, ilut_gt)

        if (pos > 0) then
            amp = trial_wfs(:,pos+min_trial_ind-1)
            min_trial_ind = min_trial_ind + pos
            tTrial = .true.
        else
            ! The state is not in the trial space. Just update min_trial_ind accordingly.
            min_trial_ind = min_trial_ind - pos - 1
        end if

        ! If pos > 0 then the state is in the trial space. A state cannot be in both the trial and
        ! connected space, so, unless pos < 0, don't bother doing the following binary search.
        if (pos < 0) then
            pos = binary_search_custom(con_space(:, min_conn_ind:con_space_size), ilut, NIfTot+1, ilut_gt)

            if (pos > 0) then
                amp = con_space_vecs(:,pos+min_conn_ind-1)
                min_conn_ind = min_conn_ind + pos
                tCon = .true.
            else
                min_conn_ind = min_conn_ind - pos - 1
            end if
        end if

    end subroutine bin_search_trial

    subroutine hash_search_trial(ilut, nI, amp, tTrial, tCon)

        integer(n_int), intent(in) :: ilut(0:)
        integer, intent(in) :: nI(nel)
        HElement_t(dp), intent(out) :: amp(:)
        logical, intent(out) :: tTrial, tCon

        integer :: i, hash_val

        amp = 0.0_dp
        tTrial = .false.
        tCon = .false.

        ! Note we search the trial space first, and don't add a contribution
        ! from the connected space if we are also in the trial space.

        if (trial_space_size > 0) then
            ! Find the hash value of this state.
            hash_val = FindWalkerHash(nI, trial_space_size)
            do i = 1, trial_ht(hash_val)%nclash
                if (DetBitEq(ilut, trial_ht(hash_val)%states(0:NIfDBO,i))) then
                    tTrial = .true.
                    amp = transfer(trial_ht(hash_val)%states(NIfDBO+1:,i), amp)
                    return
                end if
            end do
        end if
        
        ! If it wasn't in the trial space, check to see if it is in the connected space.
        if (con_space_size > 0) then
            hash_val = FindWalkerHash(nI, con_space_size)
            ! Loop over all hash clashes for this hash value.
            do i = 1, con_ht(hash_val)%nclash
                if (DetBitEq(ilut, con_ht(hash_val)%states(0:NIfDBO,i))) then
                    tCon = .true.
                    amp = transfer(con_ht(hash_val)%states(NIfDBO+1:,i), amp)
                    return
                end if
            end do
        end if

    end subroutine hash_search_trial

    subroutine get_con_amp_trial_space(ilut, amps)

        ! WARNING: This routines expects that the state passed in, ilut, is
        ! definitely in the trial space, and performs a stop_all if not.

        integer(n_int), intent(in) :: ilut(0:)
        HElement_t(dp), intent(out) :: amps(:)

        integer :: i, hash_val
        integer :: nI(nel)

        amps = 0.0_dp
        call decode_bit_det(nI, ilut)

        hash_val = FindWalkerHash(nI, con_space_size)
        ! Loop over all hash clashes for this hash value.
        do i = 1, con_ht(hash_val)%nclash
            if (all(ilut(0:NIfDBO) == con_ht(hash_val)%states(0:NIfDBO,i))) then
                amps = transfer(con_ht(hash_val)%states(NIfDBO+1:,i), amps)
                return
            end if
        end do

        call stop_all("get_trial_wf_amp","The input state is not in the trial space.")

    end subroutine get_con_amp_trial_space

    subroutine add_trial_energy_contrib(ilut, RealwSign, ireplica)
    
        integer(n_int), intent(in) :: ilut(0:)
        real(dp), intent(in) :: RealwSign
        integer, intent(in) :: ireplica

        integer :: i, hash_val
        integer :: nI(nel)
        real(dp) :: amp(ntrial_excits)

        call decode_bit_det(nI, ilut)

        ! First search the trial space.
        if (trial_space_size > 0) then
            hash_val = FindWalkerHash(nI, trial_space_size)
            do i = 1, trial_ht(hash_val)%nclash
                if (DetBitEq(trial_ht(hash_val)%states(0:NIfDBO,i), ilut)) then
                    amp = transfer(trial_ht(hash_val)%states(NIfDBO+1:,i), amp)
                    if (ntrial_excits == 1) then
                        trial_denom(ireplica) = trial_denom(ireplica) + amp(1)*RealwSign
                    else if (ntrial_excits == lenof_sign) then
                        trial_denom(ireplica) = trial_denom(ireplica) + amp(ireplica)*RealwSign
                    end if
                    return
                end if
            end do
        end if

        ! If not in the trial space, search the connected space.
        if (con_space_size > 0) then
            hash_val = FindWalkerHash(nI, con_space_size)
            do i = 1, con_ht(hash_val)%nclash
                if (DetBitEq(con_ht(hash_val)%states(0:NIfDBO,i), ilut)) then
                    amp = transfer(con_ht(hash_val)%states(NIfDBO+1:,i), amp)
                    if (ntrial_excits == 1) then
                        trial_numerator(ireplica) = trial_numerator(ireplica) + amp(1)*RealwSign
                    else if (ntrial_excits == lenof_sign) then
                        trial_numerator(ireplica) = trial_numerator(ireplica) + amp(ireplica)*RealwSign
                    end if
                    return
                end if
            end do
        end if

    end subroutine add_trial_energy_contrib

    ! This is the same as BinSearchParts1, but this time, it searches though the 
    ! full list of determinants created by the full diagonalizer when the 
    ! histogramming option is on. 
    !
    ! This is outside the module so it is accessible to AnnihilateMod
    subroutine BinSearchParts2(iLut, MinInd, MaxInd, PartInd, tSuccess)

        use DetCalcData , only : FCIDets
        use DetBitOps, only: DetBitLT
        use constants, only: n_int
        use bit_reps, only: NIfTot,NIfDBO

        IMPLICIT NONE
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

    !    WRITE(iout,*) "Binary searching between ",MinInd, " and ",MaxInd
    !    CALL neci_flush(iout)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(FCIDets(:,MaxInd),iLut(:),NIfDBO)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
    !        WRITE(iout,*) i,j,n

            ! Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is 
            ! more or 0 if they are the same
            Comp=DetBitLT(FCIDets(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
    !Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
                ! The value of the determinant at N is LESS than the determinant 
                ! we're looking for. Therefore, move the lower bound of the 
                ! search up to N. However, if the lower bound is already equal to
                ! N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN

                IF(i.eq.MaxInd-1) THEN
                    ! This deals with the case where we are interested in the 
                    ! final/first entry in the list. Check the final entry of the
                    ! list and leave. We need to check the last index.
                    Comp=DetBitLT(FCIDets(:,i+1),iLut(:),NIfDBO)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
    !final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
    !The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
    !We have failed - exit loop
                i=j
            ENDIF

        enddo

    !If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    end subroutine BinSearchParts2

    subroutine BinSearchParts_rdm(iLut,MinInd,MaxInd,PartInd,tSuccess)

        ! Do a binary search in CurrentDets, between the indices of MinInd and
        ! MaxInd. If successful, tSuccess will be true and  PartInd will be a
        ! coincident determinant. If there are multiple values, the chosen one
        ! may be any of them... If failure, then the index will be one less than
        ! the index that the particle would be in if it was present in the list.
        ! (or close enough!)

        integer(kind=n_int) :: iLut(0:NIfTot)
        integer :: MinInd, MaxInd, PartInd
        integer :: i, j, N, Comp
        logical :: tSuccess

        i = MinInd
        j = MaxInd
        if (i-j .eq. 0) then
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO)
            if (Comp .eq. 0) then
                tSuccess = .true.
                PartInd = MaxInd
                return
            else
                tSuccess = .false.
                PartInd = MinInd
            end if
        end if

        do while(j-i .gt. 0)  ! End when the upper and lower bound are the same.
            N = (i+j)/2       ! Find the midpoint of the two indices.

            ! Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is
            ! more or 0 if they are the same
            Comp = DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO)

            if (Comp .eq. 0) then
                ! Praise the lord, we've found it!
                tSuccess = .true.
                PartInd = N
                return
            else if ((Comp .eq. 1) .and. (i .ne. N)) then
                ! The value of the determinant at N is LESS than the determinant
                ! we're looking for. Therefore, move the lower bound of the search
                ! up to N. However, if the lower bound is already equal to N then
                ! the two bounds are consecutive and we have failed...
                i = N
            else if (i .eq. N) then


                if (i .eq. MaxInd-1) then
                    ! This deals with the case where we are interested in the
                    ! final/first entry in the list. Check the final entry of
                    ! the list and leave. We need to check the last index.
                    Comp = DetBitLT(CurrentDets(:,i+1), iLut(:), NIfDBO)
                    if (Comp .eq. 0) then
                        tSuccess = .true.
                        PartInd = i + 1
                        return
                    else if (Comp .eq. 1) then
                        ! final entry is less than the one we want.
                        tSuccess = .false.
                        PartInd = i + 1
                        return
                    else
                        tSuccess = .false.
                        PartInd = i
                        return
                    end if

                else if (i .eq. MinInd) then
                    tSuccess = .false.
                    PartInd = i
                    return
                else
                    i = j
                end if

            else if (Comp .eq. -1) then
                ! The value of the determinant at N is MORE than the determinant
                ! we're looking for. Move the upper bound of the search down to N.
                j = N
            else
                ! We have failed - exit loop.
                i = j
            end if

        end do

        ! If we have failed, then we want to find the index that is one less
        ! than where the particle would have been.
        tSuccess = .false.
        PartInd = max(MinInd,i-1)

    end subroutine BinSearchParts_rdm

    subroutine remove_repeated_states(list, list_size)

        use DetBitOps, only: ilut_lt, ilut_gt
        use sort_mod, only: sort

        integer, intent(inout) :: list_size
        integer(n_int), intent(inout) :: list(0:,:)
        integer :: i, counter

        if (list_size > 0) then
            ! Annihilation-like steps to remove repeated states.
            call sort(list(:, 1:list_size), ilut_lt, ilut_gt)
            counter = 1
            do i = 2, list_size
                ! If this state and the previous one were identical, don't add this state to the
                ! list so that repeats aren't included.
                if (.not. all(list(0:NIfDBO, i-1) == list(0:NIfDBO, i)) ) then
                    counter = counter + 1
                    list(:, counter) = list(:, i)
                end if
            end do

            list_size = counter
        end if

    end subroutine remove_repeated_states

end module searching
