! This code will merge two sorted lists of determinants in bit-format.
! On input, nlist1 is the length of list1. List1 must be ordered in increasing
! order (this is not checked here).
! The lists are of integers of length 0:NIfTot for each element.
! nlist2 is the length of list2. List2 must be ordered in non-decreasing
! order (this is not checked here)
! list1 must be at least of dimension nlist1+nlist2
! list2 must be at least of dimension nlist2
! on output, list1 contains the merged list
! nlist1 is the length of the merged list
! list2 is not overwritten
! Along with each array, is a seperate array, signlist(1/2) and the elements
! in this array are moved with the elements in the main arrays.
! The HList real array is an optional argument. If specified, it will find the diagonal
! matrix elements for the elements it is merging into the main list.
! The list1 will be binary searched to find insertion points. Generally, if list2 > list1/2,
! a linear search would be quicker.
!nlist1 = TotWalkersNew
!nlist2 = ValidSpawned
!list2 = SpawnedParts(:,1:nlist2)
    SUBROUTINE MergeListswH(nlist1,nlist2,list2)
        USE FciMCParMOD , only : Hii,CurrentDets
        use FciMCData , only : tFillingStochRDMonFly, InstNoatHF, ntrial_occ, &
                               ncon_occ, occ_trial_amps, occ_con_amps, &
                               trial_temp, con_temp, tTrialHash, iLutHF_True, &
                               iter
        use SystemData, only: nel, tHPHF
        use bit_rep_data, only: extract_sign, flag_trial, flag_connected
        use bit_reps, only: NIfTot, NIfDBO, decode_bit_det, test_flag, &
                            encode_first_iter
        USE Determinants , only : get_helement
        use DetBitOps, only: DetBitEQ
        use hphf_integrals, only: hphf_diag_helement
        USE CalcData , only : tTruncInitiator, tTrialWavefunction
        USE HElem
        use constants, only: dp,n_int,lenof_sign
        use util_mod, only: binary_search_custom
        use searching, only: find_trial_and_con_states_bin, find_trial_and_con_states_hash
        use global_det_data, only: global_determinant_data, set_det_diagH, &
                                   det_diagH, set_av_sgn, set_iter_occ
        IMPLICIT NONE
        INTEGER :: nlisto,nlist1,nlist2,i
        INTEGER(KIND=n_int) :: list2(0:NIfTot,1:nlist2),DetCurr(0:NIfTot) 
        INTEGER :: ips
        HElement_t :: HDiagTemp
        real(dp) :: HDiag
        INTEGER :: nJ(NEl),j,ntrial,ncon,ncon_old,ntrial_old
        logical :: tTrialState, tConState

        if (tTrialWavefunction .and. nlist2 > 0) then
            ! If using a trial wavefunction, count the number of states to be merged
            ! in which are in the trial and connected space. This routine also stores
            ! the corresponding amplitudes from trial_wf and con_space_vector in the
            ! arrays trial_temp and con_temp, and sets the flags of the new states.
            if (tTrialHash) then
                call find_trial_and_con_states_hash(int(nlist2,8), list2, ntrial, ncon)
            else
                call find_trial_and_con_states_bin(int(nlist2,8), list2, ntrial, ncon)
            end if
            ! The number of trial and connected states currently in CurrentDets.
            ntrial_old = ntrial_occ
            ncon_old = ncon_occ
            ! Update the new number of trial and connected states, for after the merge.
            ntrial_occ = ntrial_old + ntrial
            ncon_occ = ncon_old + ncon
        end if

!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
       nlisto=nlist1
       do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
           DetCurr(:)=list2(:,i)

           ! Does this state belong to the trial or connected space?
           if (tTrialWavefunction) then
               tTrialState = test_flag(list2(:,i), flag_trial)
               tConState = test_flag(list2(:,i), flag_connected)
           end if

!..ips1 is the position in list1 which num is to be inserted
! nlisto is the last element in list1 which might have to be moved to
! accommodate an element from list2.
           if (nlisto == 0) then
! No more elements in list1 to be moved.  All remaining elements in list2 need
! to be inserted with the same positions into list1.
               ips = 1
           else
               call search(nlisto,DetCurr,ips)
           end if

! Move all elements that between DetCurr and nlisto to their position in the
! completely merged list.  We know that there must be i elements still to be
! inserted...
           do j=nlisto,ips,-1
               if (tTrialWavefunction) then
                   if (test_flag(CurrentDets(:,j), flag_trial)) then
                       ! If a trial state, shift the corresponding trial vector
                       ! element. There are still ntrial elements to be inserted.
                       occ_trial_amps(ntrial_old+ntrial) = occ_trial_amps(ntrial_old)
                       ntrial_old = ntrial_old - 1
                   else if (test_flag(CurrentDets(:,j), flag_connected)) then
                       ! If a connected state, shift the corresponding connected vector
                       ! element. There are still ncon elements to be inserted.
                       occ_con_amps(ncon_old+ncon) = occ_con_amps(ncon_old)
                       ncon_old = ncon_old - 1
                   end if
               end if

               CurrentDets(:,j+i)=CurrentDets(:,j)
               global_determinant_data(:,j+i) = global_determinant_data(:,j)
           enddo

! Insert DetCurr into its position in the completely merged list (i-1 elements
! below it still to be inserted).
           CurrentDets(:,ips+i-1)=list2(:,i)
           
           ! Merge in the new trial and connected vector amplitudes.
           if (tTrialWavefunction) then
               if (tTrialState) then
                   occ_trial_amps(ntrial_old+ntrial) = trial_temp(ntrial)
                   ntrial = ntrial - 1
               else if (tConState) then
                   occ_con_amps(ncon_old+ncon) = con_temp(ncon)
                   ncon = ncon - 1
               end if
           end if
           
           ! We want to calculate the diagonal hamiltonian matrix element for
           ! the new particle to be merged.
           call decode_bit_det (nJ, list2(:,i))
           if (tHPHF) then
               HDiagTemp = hphf_diag_helement (nJ, list2(:,i))
           else
               HDiagTemp = get_helement (nJ, nJ, 0)
           endif
           HDiag=(REAL(HDiagTemp,dp))-Hii
           call set_det_diagH(ips+i-1, HDiag)
           if(DetBitEQ(CurrentDets(:,ips+i-1),iLutHF_True,NIfDBO)) call extract_sign(CurrentDets(:,ips+i-1),InstNoatHF)

           ! Store the initial iteration for this created particle
           call encode_first_iter(CurrentDets(:,ips+i-1), iter)
           
           if (tFillingStochRDMonFly) then
               if (lenof_sign == 2) then
                   call set_av_sgn(ips+i-1, (/0.0_dp, 0.0_dp/))
                   call set_iter_occ(ips+i-1, (/0.0_dp, 0.0_dp/))
               else
                   call set_av_sgn(ips+i-1, (/0.0_dp/))
                   call set_iter_occ(ips+i-1, (/0.0_dp/))
               end if
           end if
               
! Next element to be inserted must be smaller than DetCurr, so must be inserted
! at (at most) at ips-1.
! If nlisto=0 then all remaining elements in list2 must be inserted directly
! into list1---there are no more elements in list1 to be moved.
           nlisto=ips-1
        enddo
        nlist1=nlist1+nlist2

        if (tTrialWavefunction .and. nlist2 > 0) then
            if (ntrial /= 0 .or. ncon /= 0) then
                write(6,'(a7,i7,a5,i7)') "ntrial:", ntrial, "ncon:", ncon
                call neci_flush(6)
                call stop_all("MergeListswH","Incorrect number of trial or connected states merged.")
            end if
        end if

        return
    END SUBROUTINE MergeListswH

!This routine is the same as MergeListswH, but will not generate the diagonal 
!hamiltonian matrix elements to go with the inserted determinants
    SUBROUTINE MergeLists(nlist1,nlist2,list2)
        USE FciMCParMOD , only : Hii,CurrentDets
        use SystemData, only: nel
        use bit_reps, only: NIfTot
        USE CalcData , only : tTruncInitiator
        USE HElem
        use constants, only : n_int, dp
        IMPLICIT NONE
        INTEGER :: nlisto,nlist1,nlist2,i
        INTEGER(KIND=n_int) :: list2(0:NIfTot,1:nlist2),DetCurr(0:NIfTot) 
        INTEGER :: ips
        INTEGER :: j
!        LOGICAL :: tbin
!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
       nlisto=nlist1
       do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
           DetCurr(0:NIfTot)=list2(0:NIfTot,i)
           if(nlisto.eq.0) then
               ips=1
           else
               call search(nlisto,DetCurr,ips)
           endif
!          write(6,*) 'position in list1 to be inserted:',ips1
!..ips1 is the position in list1 which num is to be inserted
!           write(6,*) ' Going to insert into position:',ips
!           write(6,*) ' Copy elements from:',ips,' to',nlisto
!..if ips is less than nlisto, then no elements will be copied.
           do j=nlisto,ips,-1
             CurrentDets(0:NIfTot,j+i)=CurrentDets(0:NIfTot,j)
           enddo
!.elements of list1 which were copied over started
!.. from nlisto and went up to ips
!           write(6,*) ' position labels of newly enlarged list:'
!           write(6,'(20i15)') (j,j=1,nlist1+nlist2)
!           write(6,*) 'new enlarged list on step:',i
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)

           CurrentDets(0:NIfTot,ips+i-1)=list2(0:NIfTot,i)
               
!           write(6,*) ' newly inserted member on position:'                             &
!     & ,ips+i-1,' of value:',list2(:,i)
!           write(6,*) 'new list with insertion in it:'
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)
!..new end of list position is the next insertion point
!           nlisto=min(nlisto,ips+i-2)
           nlisto=ips-1
!           write(6,*) ' new end of list position:',nlisto
        enddo
        nlist1=nlist1+nlist2
        return
    END SUBROUTINE MergeLists
!..............................................................................
!..find the position in list such that 
!.. list(0:NIfTot,ipos-1) < DetCurr(0:NIfTot)
!.. list(0:NIfTot,ipos) ge DetCurr(0:NIfTot)
!..list is assumed to be in increasing order
    SUBROUTINE search(n,DetCurr,ipos)
        use bit_reps, only: NIfTot, NIfDBO
        use DetBitOps, only: DetBitLT
        USE FciMCParMOD , only : CurrentDets
        use constants, only: n_int
        implicit none
        integer(n_int), intent(in) :: DetCurr(0:NifTot)
        integer, intent(in) :: n
        integer, intent(out) :: ipos
        integer :: nlo, nup, i, ncurr, CompPart
!        logical :: tbin
!        if(.not.tbin) goto 200
!.......................................................................
!..binary seach
        nlo=1
        nup=n

        IF(n.eq.0) THEN
            ipos=1
            return
        ENDIF
 100    continue
!..if num is larger than the last element of list,
!.. return ipos as nup+1
        if(DetBitLT(CurrentDets(0:NIfTot,nup),DetCurr(:),NIfDBO).eq.1) then 
           ipos=nup+1
           return
        endif
!..if num is le the first element of the list
!.. return ipos as nlo
        if(DetBitLT(DetCurr(:),CurrentDets(0:NIfTot,nlo),NIfDBO).eq.1) then
           ipos=nlo
           return
        endif
!..at this point num is within the range of list
!.. take the mid-point and see how num compares
        ncurr=(nlo+nup)/2
        IF(nlo.eq.ncurr) THEN
!Insertion point is between nlo and nup...
            ipos=nlo+1
            return
        ENDIF

        CompPart=DetBitLT(CurrentDets(0:NIfTot,ncurr),DetCurr(:),NIfDBO)

!.. if list(ncurr) gt num then the upper bound to the 
!.. list can be shifted to nup
!        if(list(ncurr).gt.num) nup=ncurr
        if(CompPart.eq.-1) THEN
            nup=ncurr
        else
!..if list(ncurr) le num then the lower bound can be 
!.. can be shifted to ncurr
            nlo=ncurr
        endif
!..
!.. list(ncurr).eq.num
        if(CompPart.eq.0) then 
           write(6,"(I8)",advance='no') ncurr 
           call WriteBitDet(6,CurrentDets(:,ncurr),.true.)
           call stop_all("Search","When merging lists, no entries should exist on both lists")
!..check to see if the previous member is less than num.
!.. if so, return ipose=ncurr
           if(DetBitLT(CurrentDets(0:NIfTot,ncurr-1),DetCurr(:),NIfDBO).eq.1) then
!           if(list(ncurr-1).lt.num) then
              ipos=ncurr 
              return
           endif
!..check to see if the next member of list is ge num
!.. if so, return ipos=ncurr
           if(DetBitLT(CurrentDets(0:NIfTot,ncurr+1),DetCurr(:),NIfDBO).eq.-1) then
!           if(list(ncurr+1).gt.num) then
              ipos=ncurr 
              return
           endif
        endif
        goto 100
!...........................................................................
!        continue
!..simple linear search. At the moment, you cannot get here.
!        do i=1,n
!           if(DetBitLT(CurrentDets(0:NIfTot,i),DetCurr(:),NIfDBO).ne.1) then 
!             ipos=i
!             return
!           endif
!        enddo
!        ipos=n+1
    END SUBROUTINE Search
!..............................................................................
!..find the position in list such that 
!.. list(0:NIfTot,ipos-1) < DetCurr(0:NIfTot)
!.. list(0:NIfTot,ipos) ge DetCurr(0:NIfTot)
!..list is assumed to be in increasing order
    SUBROUTINE searchgen(n,list,DetCurr,ipos)
        use bit_reps, only: NIfTot, NIfDBO
        use DetBitOps, only: DetBitLT
        use constants, only: n_int
        IMPLICIT NONE
        INTEGER :: nlo,nup,i,ipos,ncurr,CompPart,n
        INTEGER(KIND=n_int) :: DetCurr(0:NIfTot),list(0:NIfTot,n)
!        logical :: tbin
!        if(.not.tbin) goto 200
!.......................................................................
!..binary seach
        nlo=1
        nup=n
        IF(n.eq.0) THEN
            ipos=1
            return
        ENDIF
 100    continue
!..if num is larger than the last element of list,
!.. return ipos as nup+1
        if(DetBitLT(list(:,nup),DetCurr(:),NIfDBO).eq.1) then 
           ipos=nup+1
           return
        endif
!..if num is le the first element of the list
!.. return ipos as nlo
        if(DetBitLT(DetCurr(:),list(:,nlo),NIfDBO).eq.1) then
           ipos=nlo
           return
        endif
!..at this point num is within the range of list
!.. take the mid-point and see how num compares
        ncurr=(nlo+nup)/2
        IF(nlo.eq.ncurr) THEN
!Insertion point is between nlo and nup...
            ipos=nlo+1
            return
        ENDIF

        CompPart=DetBitLT(list(:,ncurr),DetCurr(:),NIfDBO)

!.. if list(ncurr) gt num then the upper bound to the 
!.. list can be shifted to nup
!        if(list(ncurr).gt.num) nup=ncurr
        if(CompPart.eq.-1) THEN
            nup=ncurr
        else
!..if list(ncurr) le num then the lower bound can be 
!.. can be shifted to ncurr
            nlo=ncurr
        endif
!..
!.. list(ncurr).eq.num
        if(CompPart.eq.0) then 
!..check to see if the previous member is less than num.
!.. if so, return ipose=ncurr
           if(DetBitLT(list(:,ncurr-1),DetCurr(:),NIfDBO).eq.1) then
!           if(list(ncurr-1).lt.num) then
              ipos=ncurr 
              return
           endif
!..check to see if the next member of list is ge num
!.. if so, return ipos=ncurr
           if(DetBitLT(list(:,ncurr+1),DetCurr(:),NIfDBO).eq.-1) then
!           if(list(ncurr+1).gt.num) then
              ipos=ncurr 
              return
           endif
        endif
        goto 100
!...........................................................................
!        continue
!!..simple linear search. At the moment, you cannot get here.
!        do i=1,n
!           if(DetBitLT(list(:,i),DetCurr(:),NIfDBO).ne.1) then 
!             ipos=i
!             return
!           endif
!        enddo
!        ipos=n+1
    END SUBROUTINE Searchgen



    ! Insert the (already sorted) items of n2 into the sorted array
    ! n1 which only contains items up to position n1, but is large
    ! enough to contain list2
    subroutine int_list_merge (list1, list2, n1, n2)
        integer, intent(in) :: n1, n2, list2(1:n2)
        integer, intent(inout) :: list1(1:n1+n2)
        integer i, j, pos

        ! Work backwards through the lists, therefore never have to move
        ! items out of the way.
        ! j contains the reading position from list1
        ! pos contains the writing position into list1
        pos = n1+n2+1
        j = n1
        do i=n2,1,-1
            do j=j,1,-1
                pos = pos - 1
                if (list2(i) > list1(j)) then
                    list1(pos) = list2(i)
                    exit
                else
                    list1(pos) = list1(j)
                endif
            enddo
            if (j == 0) exit
        enddo

        list1(1:i) = list2(1:i)
    end subroutine

