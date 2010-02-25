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
    SUBROUTINE MergeListswH(nlist1,nlist1max,nlist2,list2,SignList2)
        USE FciMCParMOD , only : iLutHF,Hii,CurrentDets,CurrentSign,CurrentH
        USE SystemData , only : NEl,tHPHF,NIfTot,NIfDBO
        USE Determinants , only : get_helement
        use DetBitOps, only: DecodeBitDet, DetBitEQ
        USE HElem
        IMPLICIT NONE
!        INTEGER :: list1(0:NIfTot,nlist1max),list2(0:NIfTot,1:nlist2)
        INTEGER :: list2(0:NIfTot,1:nlist2)
        INTEGER :: nlisto,nlist1,nlist2,nlo,i,DetCurr(0:NIfTot) 
        INTEGER :: ips,ips1,SignList2(nlist2)!,SignList1(nlist1max),
!        REAL*8 :: HList(nlist1max)
        TYPE(HElement) :: HDiagTemp
        REAL*8 :: HDiag
        INTEGER :: nJ(NEl),j,nlist1max
!        LOGICAL :: tbin
!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
       nlisto=nlist1
       nlo=nlist1
       do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
           DetCurr(:)=list2(:,i)
           call search(nlisto,DetCurr,ips1)
!          write(6,*) 'position in list1 to be inserted:',ips1
!..ips1 is the position in list1 which num is to be inserted
           ips=ips1      
!           write(6,*) ' Going to insert into position:',ips
!           write(6,*) ' Copy elements from:',ips,' to',nlisto
!..if ips is less than nlisto, then no elements will be copied.
           do j=nlisto,ips,-1
              if(j.le.nlo) then 
!                 write(6,*) j,'->',j+i
                 CurrentDets(:,j+i)=CurrentDets(:,j)
                 CurrentSign(j+i)=CurrentSign(j)
                 CurrentH(j+i)=CurrentH(j)
              endif
           enddo
!.elements of list1 which were copied over started
!.. from nlisto and went up to ips
           nlo=ips
!           write(6,*) ' position labels of newly enlarged list:'
!           write(6,'(20i15)') (j,j=1,nlist1+nlist2)
!           write(6,*) 'new enlarged list on step:',i
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)
           CurrentDets(:,ips+i-1)=list2(:,i)
           CurrentSign(ips+i-1)=SignList2(i)
!We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
           IF(DetBitEQ(list2(:,i),iLutHF,NIfDBO)) THEN
!We know we are at HF - HDiag=0
               HDiag=0.D0
!               IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                   CALL DecodeBitDet(nJ,list2(0:NIfTot,i))
!                   HDiagTemp=GetHElement3(nJ,nJ,0)
!                   HDiag=(REAL(HDiagTemp%v,8))
!               ENDIF
           ELSE
               CALL DecodeBitDet(nJ,list2(:,i))
               IF(tHPHF) THEN
                   CALL HPHFGetDiagHElement(nJ,list2(:,i),HDiagTemp)
               ELSE
                   HDiagTemp = get_helement (nJ, nJ, 0)
               ENDIF
               HDiag=(REAL(HDiagTemp%v,8))-Hii
           ENDIF
           CurrentH(ips+i-1)=HDiag
               
!           write(6,*) ' newly inserted member on position:'                             &
!     & ,ips+i-1,' of value:',list2(:,i)
!           write(6,*) 'new list with insertion in it:'
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)
!..new end of list position is the next insertion point
!           nlisto=min(nlisto,ips+i-2)
           nlisto=min(nlisto,ips-1)
!           write(6,*) ' new end of list position:',nlisto
        enddo
        nlist1=nlist1+nlist2
        return
    END SUBROUTINE MergeListswH

                              
! This is pretty much the same as MergeListswH, however in this case, as well as sorting the determinants
! and signs, the Hii and Hij and parents are taken with the determinants too.
   SUBROUTINE MergeListswH2(nlist1,nlist1max,nlist2,list2,list3,SignList2)
        USE FciMCParMOD , only : iLutHF,Hii,MinorStarDets,MinorStarSign,MinorStarParent,MinorStarHii,MinorStarHij
        USE SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,nMsh,tHPHF,NIfTot
        USE Determinants , only : get_helement
        USE IntegralsData , only : fck,NMax,UMat
        USE HElem
        use DetBitOps, only: DecodeBitDet
        IMPLICIT NONE
        INTEGER :: list2(0:NIfTot,1:nlist2),list3(0:NIfTot,1:nlist2)
        INTEGER :: nlisto,nlist1,nlist2,nlo,i,DetCurr(0:NIfTot),DetCurr2(0:NIfTot) 
        INTEGER :: ips,ips1,SignList2(nlist2)
        TYPE(HElement) :: HDiagTemp,HOffDiagTemp
        REAL*8 :: HDiag,HOffDiag
        INTEGER :: nJ(NEl),j,nlist1max,ExcitLevel,nK(NEl)
        LOGICAL :: Det2BitEQ
!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
        IF(nlist1.gt.0) THEN
           nlisto=nlist1
           nlo=nlist1
           do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
               DetCurr(0:NIfTot)=list2(0:NIfTot,i)
               DetCurr2(0:NIfTot)=list3(0:NIfTot,i)
               call searchminor(nlisto,DetCurr,DetCurr2,ips1)
!              write(6,*) 'position in list1 to be inserted:',ips1
!..ips1 is the position in list1 which num is to be inserted
               ips=ips1      
!           write(6,*) ' Going to insert into position:',ips
!           write(6,*) ' Copy elements from:',ips,' to',nlisto
!..if ips is less than nlisto, then no elements will be copied.
               do j=nlisto,ips,-1
                  if(j.le.nlo) then 
!                 write(6,*) j,'->',j+i
                     MinorStarDets(0:NIfTot,j+i)=MinorStarDets(0:NIfTot,j)
                     MinorStarParent(0:NIfTot,j+i)=MinorStarParent(0:NIfTot,j)
                     MinorStarSign(j+i)=MinorStarSign(j)
                     MinorStarHii(j+i)=MinorStarHii(j)
                     MinorStarHij(j+i)=MinorStarHij(j+1)
                  endif
               enddo
!.elements of list1 which were copied over started
!.. from nlisto and went up to ips
               nlo=ips
!           write(6,*) ' position labels of newly enlarged list:'
!           write(6,'(20i15)') (j,j=1,nlist1+nlist2)
!           write(6,*) 'new enlarged list on step:',i
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)
               MinorStarDets(0:NIfTot,ips+i-1)=list2(0:NIfTot,i)
               MinorStarSign(ips+i-1)=SignList2(i)
               MinorStarParent(0:NIfTot,ips+i-1)=list3(0:NIfTot,i)

! Want to calculate the diagonal and off diagonal H elements of the particle to be merged.           
               CALL DecodeBitDet(nJ,list2(:,i))
               IF(tHPHF) THEN
                   CALL HPHFGetDiagHElement(nJ,list2(:,i),HDiagTemp)
               ELSE
                   HDiagTemp = get_helement (nJ, nJ, 0)
               ENDIF
               HDiag=(REAL(HDiagTemp%v,8))-Hii
               MinorStarHii(ips+i-1)=HDiag

               CALL DecodeBitDet(nK,list3(:,i))
               IF(tHPHF) THEN
                   CALL HPHFGetOffDiagHElement(nJ,nK,list2(:,i),list3(:,i),HOffDiagTemp)
               ELSE
                   HOffDiagTemp = get_helement(nJ, nK, list2(:,i), list3(:,i))
               ENDIF
               HOffDiag=(REAL(HOffDiagTemp%v,8))
               MinorStarHij(ips+i-1)=HOffDiag


               nlisto=min(nlisto,ips-1)
!           write(6,*) ' new end of list position:',nlisto
            enddo
            nlist1=nlist1+nlist2
            return
        ELSE
! If there are no entries in the star arrays to merge with, just copy the spawned walkers straight over to star array            
            do j=1,nlist2
                MinorStarDets(:,j)=list2(:,j)
                MinorStarSign(j)=SignList2(j)
                MinorStarParent(:,j)=list3(:,j)

                CALL DecodeBitDet(nJ,list2(:,j))
                IF(tHPHF) THEN
                    CALL HPHFGetDiagHElement(nJ,list2(:,j),HDiagTemp)
                ELSE
                    HDiagTemp = get_helement (nJ, nJ, 0)
                ENDIF
                HDiag=(REAL(HDiagTemp%v,8))-Hii
                MinorStarHii(j)=HDiag

                CALL DecodeBitDet(nK,list3(:,j))
                IF(tHPHF) THEN
                    CALL HPHFGetOffDiagHElement(nJ,nK,list2(:,j),list3(:,j),HOffDiagTemp)
                ELSE
                    HOffDiagTemp = get_helement(nJ, nK, list2(:,j),list3(:,j))
                ENDIF
                HOffDiag=(REAL(HOffDiagTemp%v,8))
                MinorStarHij(j)=HOffDiag
            enddo
            nlist1=nlist2
            return
        ENDIF


    END SUBROUTINE MergeListswH2



!This routine is the same as MergeListswH, but will not generate the diagonal hamiltonian matrix elements to go with the inserted determinants
    SUBROUTINE MergeLists(nlist1,nlist1max,nlist2,list2,SignList2)
        USE FciMCParMOD , only : iLutHF,Hii,CurrentDets,CurrentSign
        USE SystemData , only : NEl, NIfTot
        USE HElem
        IMPLICIT NONE
        INTEGER :: list2(0:NIfTot,1:nlist2)
        INTEGER :: nlisto,nlist1,nlist2,nlo,i,DetCurr(0:NIfTot) 
        INTEGER :: ips,ips1,SignList2(nlist2)!,SignList1(nlist1max)
        REAL*8 :: HDiag
        INTEGER :: nJ(NEl),j,nlist1max
!        LOGICAL :: tbin
!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
       nlisto=nlist1
       nlo=nlist1
       do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
           DetCurr(0:NIfTot)=list2(0:NIfTot,i)
           call search(nlisto,DetCurr,ips1)
!          write(6,*) 'position in list1 to be inserted:',ips1
!..ips1 is the position in list1 which num is to be inserted
           ips=ips1      
!           write(6,*) ' Going to insert into position:',ips
!           write(6,*) ' Copy elements from:',ips,' to',nlisto
!..if ips is less than nlisto, then no elements will be copied.
           do j=nlisto,ips,-1
              if(j.le.nlo) then 
!                 write(6,*) j,'->',j+i
                 CurrentDets(0:NIfTot,j+i)=CurrentDets(0:NIfTot,j)
                 CurrentSign(j+i)=CurrentSign(j)
              endif
           enddo
!.elements of list1 which were copied over started
!.. from nlisto and went up to ips
           nlo=ips
!           write(6,*) ' position labels of newly enlarged list:'
!           write(6,'(20i15)') (j,j=1,nlist1+nlist2)
!           write(6,*) 'new enlarged list on step:',i
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)
           CurrentDets(0:NIfTot,ips+i-1)=list2(0:NIfTot,i)
           CurrentSign(ips+i-1)=SignList2(i)
               
!           write(6,*) ' newly inserted member on position:'                             &
!     & ,ips+i-1,' of value:',list2(:,i)
!           write(6,*) 'new list with insertion in it:'
!           do j=1,nlist1+nlist2
!               write(6,'(20i15)') j,list1(:,j)
!           enddo
!           write(6,'(20i15)') (list1(:,j),j=1,nlist1+nlist2)
!..new end of list position is the next insertion point
!           nlisto=min(nlisto,ips+i-2)
           nlisto=min(nlisto,ips-1)
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
        use SystemData, only: NIfTot,NIfDBO
        use DetBitOps, only: DetBitLT
        USE FciMCParMOD , only : CurrentDets
        IMPLICIT NONE
        INTEGER :: n,DetCurr(0:NIfTot)!,list(0:NIFTot,n)
        INTEGER :: nlo,nup,i,ipos,ncurr,CompPart
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
 200    continue
!..simple linear search. At the moment, you cannot get here.
        do i=1,n
           if(DetBitLT(CurrentDets(0:NIfTot,i),DetCurr(:),NIfDBO).ne.1) then 
             ipos=i
             return
           endif
        enddo
        ipos=n+1
    END SUBROUTINE Search
!..............................................................................
!..find the position in list such that 
!.. list(0:NIfTot,ipos-1) < DetCurr(0:NIfTot)
!.. list(0:NIfTot,ipos) ge DetCurr(0:NIfTot)
!..list is assumed to be in increasing order
    SUBROUTINE searchgen(n,list,DetCurr,ipos)
        use SystemData, only: NIfTot,NIfDBO
        use DetBitOps, only: DetBitLT
        IMPLICIT NONE
        INTEGER :: n,DetCurr(0:NIfTot)
        INTEGER :: nlo,nup,i,ipos,ncurr,CompPart
        INTEGER :: list(0:NIfTot,n)
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
 200    continue
!..simple linear search. At the moment, you cannot get here.
        do i=1,n
           if(DetBitLT(list(:,i),DetCurr(:),NIfDBO).ne.1) then 
             ipos=i
             return
           endif
        enddo
        ipos=n+1
    END SUBROUTINE Searchgen


!..............................................................................
!..find the position in list such that 
!.. list(0:NIfTot,ipos-1) le DetCurr(0:NIfTot)
!.. list(0:NIfTot,ipos) ge DetCurr(0:NIfTot)
!.. AND
!.. list2(0:NIfTot,ipos-1) < DetCurr2(0:NIfTot)
!.. list2(0:NIfTot,ipos) ge DetCurr2(0:NIfTot)
!..list is assumed to be in increasing order
!..i.e inserting an entry in two lists in the correct position relative to both lists.
    SUBROUTINE searchminor(n,DetCurr,DetCurr2,ipos)
        use SystemData, only: NIfTot,NIfDBO
        use DetBitOps, only: Det2BitLT
        USE FciMCParMOD , only : MinorStarDets,MinorStarParent
        IMPLICIT NONE
        INTEGER :: n,DetCurr(0:NIfTot),DetCurr2(0:NIfTot)!,list(0:NIFd,n)
        INTEGER :: nlo,nup,i,ipos,ncurr,CompPart
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
!        WRITE(6,*) 'in 100 loop'
        if(Det2BitLT(MinorStarDets(:,nup),DetCurr(:),MinorStarParent(:,nup),DetCurr2(:),NIfDBO).eq.1) then 
           ipos=nup+1
           return
        endif
!..if num is le the first element of the list
!.. return ipos as nlo
        if(Det2BitLT(DetCurr(:),MinorStarDets(:,nlo),DetCurr2(:),MinorStarParent(:,nlo),NIfDBO).eq.1) then
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

        CompPart=Det2BitLT(MinorStarDets(:,ncurr),DetCurr(:),MinorStarParent(:,ncurr),DetCurr2(:),NIfDBO)
        ! Compares determinants w regards to both the determinants and the parents.


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
           if(Det2BitLT(MinorStarDets(:,ncurr-1),DetCurr(:),MinorStarParent(:,ncurr-1),DetCurr2(:),NIfDBO).eq.1) then
!           if(list(ncurr-1).lt.num) then
              ipos=ncurr 
              return
           endif
!..check to see if the next member of list is ge num
!.. if so, return ipos=ncurr
           if(Det2BitLT(MinorStarDets(:,ncurr+1),DetCurr(:),MinorStarParent(:,ncurr+1),DetCurr2(:),NIfDBO).eq.-1) then
!           if(list(ncurr+1).gt.num) then
              ipos=ncurr 
              return
           endif
        endif
        goto 100
!...........................................................................
 200    continue
!..simple linear search. At the moment, you cannot get here.
        do i=1,n
           if(Det2BitLT(MinorStarDets(:,i),DetCurr(:),MinorStarParent(:,i),DetCurr2(:),NIfDBO).ne.1) then 
             ipos=i
             return
           endif
        enddo
        ipos=n+1
    END SUBROUTINE searchminor

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
