! This code will merge two sorted lists of determinants in bit-format.
! On input, nlist1 is the length of list1. List1 must be ordered in increasing
! order (this is not checked here).
! The lists are of integers of length 0:NIfD for each element.
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
    SUBROUTINE MergeListswH(nlist1,nlist1max,nlist2,list2,SignList2,NIfD)
        USE FciMCParMOD , only : iLutHF,Hii,CurrentDets,CurrentSign,CurrentH
        USE SystemData , only : NEl
        USE Determinants , only : GetHElement3
        USE HElem
        IMPLICIT NONE
!        INTEGER :: list1(0:NIfD,nlist1max),list2(0:NIfD,1:nlist2)
        INTEGER :: list2(0:NIfD,1:nlist2)
        INTEGER :: nlisto,nlist1,nlist2,NIfD,nlo,i,DetCurr(0:NIfD) 
        INTEGER :: ips,ips1,SignList2(nlist2)!,SignList1(nlist1max),
!        REAL*8 :: HList(nlist1max)
        TYPE(HElement) :: HDiagTemp,GetHElement2
        REAL*8 :: HDiag
        INTEGER :: nJ(NEl),j,nlist1max
        LOGICAL :: DetBitEQ
!        LOGICAL :: tbin
!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
       nlisto=nlist1
       nlo=nlist1
       do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
           DetCurr(0:NIfD)=list2(0:NIfD,i)
           call search(nlisto,DetCurr,ips1,NIfD)
!          write(6,*) 'position in list1 to be inserted:',ips1
!..ips1 is the position in list1 which num is to be inserted
           ips=ips1      
!           write(6,*) ' Going to insert into position:',ips
!           write(6,*) ' Copy elements from:',ips,' to',nlisto
!..if ips is less than nlisto, then no elements will be copied.
           do j=nlisto,ips,-1
              if(j.le.nlo) then 
!                 write(6,*) j,'->',j+i
                 CurrentDets(0:NIfD,j+i)=CurrentDets(0:NIfD,j)
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
           CurrentDets(0:NIfD,ips+i-1)=list2(0:NIfD,i)
           CurrentSign(ips+i-1)=SignList2(i)
!We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
           IF(DetBitEQ(list2(0:NIfD,i),iLutHF,NIfD)) THEN
!We know we are at HF - HDiag=0
               HDiag=0.D0
!               IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                   CALL DecodeBitDet(nJ,list2(0:NIfD,i),NEl,NIfD)
!                   HDiagTemp=GetHElement3(nJ,nJ,0)
!                   HDiag=(REAL(HDiagTemp%v,8))
!               ENDIF
           ELSE
               CALL DecodeBitDet(nJ,list2(0:NIfD,i),NEl,NIfD)
               HDiagTemp=GetHElement3(nJ,nJ,0)
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

!This routine is the same as MergeListswH, but will not generate the diagonal hamiltonian matrix elements to go with the inserted determinants
    SUBROUTINE MergeLists(nlist1,nlist1max,nlist2,list2,SignList2,NIfD)
        USE FciMCParMOD , only : iLutHF,Hii,CurrentDets,CurrentSign
        USE SystemData , only : NEl
        USE Determinants , only : GetHElement3
        USE HElem
        IMPLICIT NONE
        INTEGER :: list1(0:NIfD,nlist1max),list2(0:NIfD,1:nlist2)
        INTEGER :: nlisto,nlist1,nlist2,NIfD,nlo,i,DetCurr(0:NIfD) 
        INTEGER :: ips,ips1,SignList2(nlist2)!,SignList1(nlist1max)
        TYPE(HElement) :: HDiagTemp,GetHElement2
        REAL*8 :: HDiag
        INTEGER :: nJ(NEl),j,nlist1max
        LOGICAL :: DetBitEQ
!        LOGICAL :: tbin
!.................................................................
!..starting from the end of the list, expand list1 to accomodate
!.. elements of list2
       nlisto=nlist1
       nlo=nlist1
       do i=nlist2,1,-1
!.. find the positions in list1 which the list2 would be inserted
           DetCurr(0:NIfD)=list2(0:NIfD,i)
           call search(nlisto,DetCurr,ips1,NIfD)
!          write(6,*) 'position in list1 to be inserted:',ips1
!..ips1 is the position in list1 which num is to be inserted
           ips=ips1      
!           write(6,*) ' Going to insert into position:',ips
!           write(6,*) ' Copy elements from:',ips,' to',nlisto
!..if ips is less than nlisto, then no elements will be copied.
           do j=nlisto,ips,-1
              if(j.le.nlo) then 
!                 write(6,*) j,'->',j+i
                 CurrentDets(0:NIfD,j+i)=CurrentDets(0:NIfD,j)
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
           CurrentDets(0:NIfD,ips+i-1)=list2(0:NIfD,i)
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
!.. list(0:NIfD,ipos-1) < DetCurr(0:NIfD)
!.. list(0:NIfD,ipos) ge DetCurr(0:NIfD)
!..list is assumed to be in increasing order
    SUBROUTINE search(n,DetCurr,ipos,NIfD)
        USE FciMCParMOD , only : CurrentDets
        IMPLICIT NONE
        INTEGER :: n,NIfD,DetCurr(0:NIfD)!,list(0:NIFd,n)
        INTEGER :: nlo,nup,DetBitLT,i,ipos,ncurr,CompPart
!        logical :: tbin
!        if(.not.tbin) goto 200
!.......................................................................
!..binary seach
        nlo=1
        nup=n
 100    continue
!..if num is larger than the last element of list,
!.. return ipos as nup+1
        if(DetBitLT(CurrentDets(0:NIfD,nup),DetCurr(0:NIfD),NIfD).eq.1) then 
           ipos=nup+1
           return
        endif
!..if num is le the first element of the list
!.. return ipos as nlo
        if(DetBitLT(DetCurr(0:NIfD),CurrentDets(0:NIfD,nlo),NIfD).eq.1) then
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

        CompPart=DetBitLT(CurrentDets(0:NIfD,ncurr),DetCurr(0:NIfD),NIfD)

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
           if(DetBitLT(CurrentDets(0:NIfD,ncurr-1),DetCurr(0:NIfD),NIfD).eq.1) then
!           if(list(ncurr-1).lt.num) then
              ipos=ncurr 
              return
           endif
!..check to see if the next member of list is ge num
!.. if so, return ipos=ncurr
           if(DetBitLT(CurrentDets(0:NIfD,ncurr+1),DetCurr(0:NIfD),NIfD).eq.-1) then
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
           if(DetBitLT(CurrentDets(0:NIfD,i),DetCurr(0:NIfD),NIfD).ne.1) then 
             ipos=i
             return
           endif
        enddo
        ipos=n+1
    END SUBROUTINE Search
