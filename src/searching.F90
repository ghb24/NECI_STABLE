#include  "macros.h"

module searching

    use bit_rep_data, only: nIfDBO, NIfTot, flag_trial, flag_connected, test_flag
    use bit_reps, only: decode_bit_det, set_flag
    use constants
    use DetBitOps, only: DetBitLt, ilut_gt, DetBitEq
    use FciMCData, only: CurrentDets, trial_space, min_trial_ind, trial_space_size, trial_wf, &
                         con_space, min_conn_ind, con_space_size, con_space_vector, &
                         trial_numerator, trial_denom, Trial_Search_Time, trial_temp, con_temp
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
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO,.false.)
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
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO,.false.)

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
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfDBO,.false.)
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

    subroutine bin_search_trial(ilut, amp)

        integer(n_int), intent(inout) :: ilut(0:)
        real(dp), intent(out) :: amp
        integer :: i, pos

        amp = 0.0_dp

        ! Search both the trial space and connected space to see if this state exists in either list.
        ! First the trial space:
        pos = binary_search_custom(trial_space(:, min_trial_ind:trial_space_size), ilut, NIfTot+1, ilut_gt)

        if (pos > 0) then
            amp = trial_wf(pos+min_trial_ind-1)
            min_trial_ind = min_trial_ind + pos
            call set_flag(ilut, flag_trial)
        else
            ! The state is not in the trial space. Just update min_trial_ind accordingly.
            min_trial_ind = min_trial_ind - pos - 1
        end if

        ! If pos > 0 then the state is in the trial space. A state cannot be in both the trial and
        ! connected space, so, unless pos < 0, don't bother doing the following binary search.
        if (pos < 0) then

            pos = binary_search_custom(con_space(:, min_conn_ind:con_space_size), ilut, NIfTot+1, ilut_gt)

            if (pos > 0) then
                amp = con_space_vector(pos+min_conn_ind-1)
                min_conn_ind = min_conn_ind + pos
                call set_flag(ilut, flag_connected)
            else
                min_conn_ind = min_conn_ind - pos - 1
            end if
        end if

    end subroutine bin_search_trial

    subroutine hash_search_trial(ilut, nI, amp)

        integer(n_int), intent(inout) :: ilut(0:)
        integer, intent(in) :: nI(nel)
        real(dp), intent(out) :: amp
        integer :: i, hash_val

        amp = 0.0_dp
        
        if (con_space_size > 0) then
            ! Find the hash value of this state.
            hash_val = FindWalkerHash(nI, con_space_size)
            ! Loop over all hash clashes for this hash value.
            do i = 1, con_ht(hash_val)%nclash
                if (DetBitEq(ilut, con_ht(hash_val)%states(0:NIfDBO,i))) then
                    call set_flag(ilut, flag_connected)
                    amp = transfer(con_ht(hash_val)%states(NIfDBO+1,i), amp)
                    return
                end if
            end do
        end if

        if (trial_space_size > 0) then
            ! If it wasn't in the connected space, check to see if it is in the trial space.
            hash_val = FindWalkerHash(nI, trial_space_size)
            do i = 1, trial_ht(hash_val)%nclash
                if (DetBitEq(ilut, trial_ht(hash_val)%states(0:NIfDBO,i))) then
                    call set_flag(ilut, flag_trial)
                    amp = transfer(trial_ht(hash_val)%states(NIfDBO+1,i), amp)
                    return
                end if
            end do
        end if

    end subroutine hash_search_trial

    subroutine add_trial_energy_contrib(ilut, RealwSign)
    
        integer(n_int), intent(in) :: ilut(0:)
        real(dp), intent(in) :: RealwSign
        integer :: i, hash_val
        integer :: nI(nel)
        real(dp) :: amp

        call decode_bit_det(nI, ilut)

        ! First search the connected space.
        if (con_space_size > 0) then
            hash_val = FindWalkerHash(nI, con_space_size)
            do i = 1, con_ht(hash_val)%nclash
                if (DetBitEq(con_ht(hash_val)%states(0:NIfDBO,i), ilut)) then
                    amp = transfer(con_ht(hash_val)%states(NIfDBO+1,i), amp)
                    trial_numerator = trial_numerator + amp*RealwSign
                    return
                end if
            end do
        end if

        ! If not in the connected space, search the trial space.
        if (trial_space_size > 0) then
            hash_val = FindWalkerHash(nI, trial_space_size)
            do i = 1, trial_ht(hash_val)%nclash
                if (DetBitEq(trial_ht(hash_val)%states(0:NIfDBO,i), ilut)) then
                    amp = transfer(trial_ht(hash_val)%states(NIfDBO+1,i), amp)
                    trial_denom = trial_denom + amp*RealwSign
                    return
                end if
            end do
        end if

    end subroutine add_trial_energy_contrib

    subroutine find_trial_and_con_states_bin(num_states, ilut_list, ntrial, ncon)

        integer(int64), intent(in) :: num_states
        integer(n_int), intent(inout) :: ilut_list(:,:)
        integer, intent(out) :: ntrial, ncon
        integer :: pos
        integer(int64) :: i

        min_trial_ind = 1
        min_conn_ind = 1
        ntrial = 0
        ncon = 0

        call set_timer(Trial_Search_Time)

        do i = 1, num_states
            ! Search both the trial space and connected space to see if this state exists in either list.
            if (min_trial_ind <= trial_space_size) then
                ! First the trial space:
                pos = binary_search_custom(trial_space(:, min_trial_ind:trial_space_size), &
                                           ilut_list(:,i), NIfTot+1, ilut_gt)

                if (pos > 0) then
                    ntrial = ntrial + 1
                    trial_temp(ntrial) = trial_wf(pos+min_trial_ind-1)
                    min_trial_ind = min_trial_ind + pos
                    call set_flag(ilut_list(:,i), flag_trial)
                else
                    ! The state is not in the trial space. Just update min_trial_ind accordingly.
                    min_trial_ind = min_trial_ind - pos - 1
                end if
            else
                ! To make sure that the connected space can be searched next.
                pos = -1
            end if

            ! If pos > 0 then the state is in the trial space. A state cannot be in both the trial and
            ! connected space, so, unless pos < 0, don't bother doing the following binary search.
            if (pos < 0 .and. min_conn_ind <= con_space_size) then

                pos = binary_search_custom(con_space(:, min_conn_ind:con_space_size), &
                                           ilut_list(:,i), NIfTot+1, ilut_gt)

                if (pos > 0) then
                    ncon = ncon + 1
                    con_temp(ncon) = con_space_vector(pos+min_conn_ind-1)
                    min_conn_ind = min_conn_ind + pos
                    call set_flag(ilut_list(:,i), flag_connected)
                else
                    min_conn_ind = min_conn_ind - pos - 1
                end if
            end if
        end do

        call halt_timer(Trial_Search_Time)

    end subroutine find_trial_and_con_states_bin

    subroutine find_trial_and_con_states_hash(num_states, ilut_list, ntrial, ncon)

        integer(int64), intent(in) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer, intent(out) :: ntrial, ncon
        integer :: pos
        integer(int64) :: i
        integer :: nI(nel)
        real(dp) :: amp

        ntrial = 0
        ncon = 0

        call set_timer(Trial_Search_Time)

        do i = 1, num_states
            call decode_bit_det(nI, ilut_list(:,i))
            ! Search the trial and connected list to see if this state exists in
            ! either. If it is, this routine sets the corresponding flag and returns the
            ! corresponding amplitude.
            call hash_search_trial(ilut_list(:,i), nI, amp)

            if (test_flag(ilut_list(:,i),flag_trial)) then
                ! If this state is in the trial space.
                ntrial = ntrial + 1
                trial_temp(ntrial) = amp
            else if(test_flag(ilut_list(:,i),flag_connected)) then
                ! If this state is in the connected space.
                ncon = ncon + 1
                con_temp(ncon) = amp
            end if
        end do

        call halt_timer(Trial_Search_Time)

    end subroutine find_trial_and_con_states_hash

    ! This is the same as BinSearchParts1, but this time, it searches though the 
    ! full list of determinants created by the full diagonalizer when the 
    ! histogramming option is on. 
    !
    ! This is outside the module so it is accessible to AnnihilateMod
    SUBROUTINE BinSearchParts2(iLut,MinInd,MaxInd,PartInd,tSuccess)
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

    END SUBROUTINE BinSearchParts2

end module searching
