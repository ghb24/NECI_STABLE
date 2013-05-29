! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module CCMCData
   use constants, only: dp,n_int,end_n_int,lenof_sign
   implicit none
   save
   real(dp)   dT1SqCuml
   logical  tExactCluster  ! Go through all combinations of excitors to make all clusters
   logical  tExactSpawn    ! For each cluster, go through all connected dets, and spawn there
   integer  nSpawnings     ! The number of spawning events per cluster if not tExactSpawn
   integer  nClustSelections !The number of cluster selections per cycle
   real(dp)   dClustSelectionRatio !The number of cluster selections per excitor (if nClustSelections==-1)
   logical  tCCMCFCI       ! Run CCMC code without excitation clusters, recovering the FCIMC result
   logical  tAmplitudes    ! Use real numbers to indicate the amplitudes rather than stochastically sampling
   real(dp)   dInitAmplitude ! Specify the initial amplitude for use in CCMC amplitude calculations.
   real(dp)   dProbSelNewExcitor !The probability that the cluster selection algorithm terminates after each addition of an excitor.
   LOGICAL  tSpawnProp     ! Set if we use spawning proportional to the cluster amplitude rather than equally

   LOGICAL  tCCBuffer      ! Buffer the CC Amplitudes - this is useful when there are many cluster selections which 
                !lead to the same collapsed det. It creates a combined amplitude of the det first and spawns from that.
   LOGICAL  tExactEnergy   ! Do we calculate projected energy exactly rather than through sampling?
   LOGICAL  tSharedExcitors !Do we use shared memory for the excitor list?
   LOGICAL  tCCNoCuml        ! Do we use a new algorithm which doesn't create a cumulative amplitude list
   

!This contains information as to a chosen Cluster
TYPE Cluster 
   INTEGER(KIND=n_int), allocatable :: SelectedExcitors(:,:)      !(0:NIfTot,nEl)  !The excitors which make up this cluster
   INTEGER, allocatable :: SelectedExcitorIndices(:)  !(nEl)     !The indices in a list of excitors, of the excitors 
                                                                !which make up this cluster
   INTEGER(KIND=n_int), allocatable :: iLutDetCurr(:)             !(0:NIfTot) The determinant made from 
                                                                !collapsing this cluster in bit representation
   INTEGER, allocatable :: DetCurr(:)                 !(nEl) The determinant made from collapsing this cluster.
   INTEGER  iSize
   INTEGER,dimension(lenof_sign) :: iSgn                                      !The sign of the 
                                                        !determinant after collapsing the cluster
   INTEGER iExcitLevel                                !The excitation level of the resultant det

   INTEGER initFlag                                   !Zero if this cluster is an initiator, or 1 if it isn't
   real(dp)   dAbsAmplitude
! dAbsAmplitude is the product of the coefficients of the excitors with the relevant normalizations.
! i.e. abs ( N0  (tI/N0) (tJ/N0) ... )
   real(dp)   dSelectionProb
! dSelectionProb is the probability that the cluster was selected
   
   real(dp) dSelectionNorm



!The following are old
   real(dp)   dProbNorm
! dProbNorm is the prob that a cluster in this level would've been chosen had they been equally weighted
   real(dp)   dClusterProb
!dClusterProb is the Probability of having chosen this cluster excitor (normalized such that <1/dClusterProb> = 1)
   real(dp)   dClusterNorm                           
!  dClusterNorm is the probability that this cluster was chosen, given the level had already 
!been selected.  This includes multiple selections of the same excitor as well as combinations of excitors which produce a 0 sign.
END TYPE Cluster

TYPE ClustSelector 
   INTEGER iIndex
   LOGICAL tFull     !Set if we are to generate all possible clusters
   INTEGER iMaxSize  !The maximum size of a cluster
   INTEGER nSelects  !If we're stochastically sampling the cluster space, this is the number of samples we take
   real(dp) dRatio     !If tDynamic then dRatio is the ratio of selections to excitors
   real(dp) dProbSelNewExcitor  !The probability that we quit at every stage of selecting a new excitor for a cluster  
   INTEGER iRefPos   !The Location in teh amplitude list of the reference det
   LOGICAL tDynamic  !If set, we choose as many clusters as there are excitors.
   LOGICAL tInitiators !Set if we are using initiators
   real(dp) dInitiatorThresh !Threshold for creating intiator cluster
   TYPE(Cluster) C

END TYPE ClustSelector

TYPE Spawner 
   LOGICAL tFull        !Set if we go through all possible spawnees sequentially
   INTEGER nSpawnings   !The number of spawning events to attempt if we are randomly spawning
   INTEGER iIndex       !The index of the current spawning event
   INTEGER iMaxExcitLevel !The Max level away from the reference that an excit canbe
!   INTEGER, allocatable :: nI(:)           !The det from which to spawn
!   INTEGER, allocatable :: iLutnI(:)       !The det from which to spawn
   LOGICAL bValid       !Set if a valid det was found
   Type(Cluster), pointer:: C
   INTEGER, allocatable :: Scratch1(:)
   INTEGER, allocatable :: Scratch2(:)
   INTEGER, allocatable :: nJ(:)           !The det which is spawned to.
   INTEGER(KIND=n_int), allocatable :: iLutnJ(:)       !The det from which to spawn
   INTEGER iExcitLevel                     !The excitation level of the resultant det from the composite cluster
   HElement_t       :: HIJ
   real(dp)               :: dProbSpawn      !Prob that we spawned here (including the number of spawning events)
   INTEGER              :: ExcitMat(2,2)   !Internal data corresponding to the excitation matrix of the last generated det.
   INTEGER              :: ExFlag
END TYPE Spawner 

TYPE CCTransitionLog
   real(dp), allocatable :: dProbTransition(:,:,:,:) !(2,2,nClust,nClust)
   real(dp), allocatable :: dProbClust(:,:)  !(2,nClust)
   real(dp), allocatable :: dProbUniqClust(:,:)  !(2,-1:nClust)
   INTEGER nExcitors    !Number of excitors
   INTEGER nMaxSize     !Largest cluster Size
   INTEGER MaxIndex     !Last possible index of a cluster (nClust)
   LOGICAL tNonUniq     !Set if we're doing non-uniq logging
END TYPE CCTransitionLog

contains

!Excitors and determinants are different although both can be specified by a LookUpTable.
!  Excitors are of the form
!  e.g. t_ij^ab  means   a^+_a a^+_b a_j a_i  (where i<j and a<b).
!  our determinants are specified by an ordered list of occupied orbitals:
!    i,j,k  (i<j<k)
!  This correpsonds to a^+_i a^+_j a^+_k |0>

!  Applying the excitor to the reference det may lead to a change in sign.  That is calculated here.

FUNCTION ExcitToDetSign(iLutRef,iLutDet,iLevel)
   use SystemData, only: nEl
   use bit_rep_data, only: NIfDBO, NIfD, NIfTot
   IMPLICIT NONE
   INTEGER ExcitToDetSign
   INTEGER iLevel
   INTEGER(KIND=n_int) iLutRef(0:nIfTot),iLutDet(0:nIfTot)
   INTEGER iSgn,i,j
   INTEGER(KIND=n_int) mask
   INTEGER iAnnihil, iCreation
   iSgn=1
   iAnnihil=iLevel
   iCreation=iLevel
!   write(6,*) "Excitation level ",iLevel
!   write(6,*) "Ref",iLutRef
!   write(6,*) "Det",iLutDet
   DO i=0,nIfD
      mask=ieor(iLutRef(i),iLutDet(i))
      Do j=0,end_n_int
         if(btest(iLutRef(i),j)) then
! electron is in ref det
!            WRITE(6,*) "Bit ",i*31+j," set in ref"
            if(btest(mask,j)) then
!               WRITE(6,*) "Bit ",i*31+j," not set in det"
               ! electron not in this det, so we annihilate
               iAnnihil=iAnnihil-1
            else
!propagate the rest of the operators through
               if(iand(iAnnihil+iCreation,1).ne.0) iSgn=-iSgn  
            endif
         else
!virtual orb
            if(btest(mask,j)) then
               ! electron is in this det, so we create
               iCreation=iCreation-1
            endif
         endif
      enddo
   enddo
   if(iAnnihil+iCreation.ne.0) then
      call WriteBitEx(6,iLutRef,iLutDet,.false.)
      write(6,"(A)",advance='no') "Ref, Det:"
      write(6,"(8Z17)") iLutRef, iLutDet
      write(6,*) "bits/byte", end_n_int+1
      write(6,*) "Level:", iLevel
      write(6,"(A,Z17,A,Z17)") " left over annihil: ", iAnnihil," left over creation: ", iCreation
      call Stop_All("Failed to normal order excitor.","ExcitToDetSign")
   endif
!   call WriteBitEx(6,iLutRef,iLutDet,.false.)
!   write(6,*) " Ex2Det sign: ",iSgn
   ExcitToDetSign=iSgn
   return
end function ExcitToDetSign            
!Add the excitation in iLutnJ to iLutnI and return it in iLutnI.  iSgn is
!updated with the relevant permutation or set to zero if the excitation is
!disallowed.
SUBROUTINE AddBitExcitor(iLutnI,iLutnJ,iLutRef,iSgn)
   use SystemData, only : nEl
   use DetBitOps, only: FindBitExcitLevel
   use bit_rep_data, only: NIfDBO,NIfD,NIfTot
   IMPLICIT NONE
   INTEGER(KIND=n_int) iLutnI(0:nIfTot), iLutnJ(0:nIfTot),iLutRef(0:nIfTot)
   INTEGER(KIND=n_int) iLutTmp(0:nIfTot)
   INTEGER(KIND=n_int) T1,T2,T3
   INTEGER iSgn
! We need to run through the bits of J and I concurrently, setting bits of I
   INTEGER i,j

   integer iIlevel,iJlevel,iTmpLevel,iI,iJ
   iIlevel = FindBitExcitLevel(iLutRef, iLutnI, nEl)
   iJlevel = FindBitExcitLevel(iLutRef, iLutnJ, nEl)
   iI=iIlevel
   iJ=iJlevel
   iTmpLevel=0
!   WRITE(6,'(A,I)',advance='no') 'x1,x2,s: ',iSgn
!   call WriteBitEx(6,iLutRef,iLutnI,.false.)
!   call WriteBitEx(6,iLutRef,iLutnJ,.true.)

!We specify the definitions of the excitors in terms of creation and annihilation operators:

!  e.g. t_ij^ab  means   a^+_a a^+_b a_j a_i  (where i<j and a<b).
!  our determinants are specified by an ordered list of occupied orbitals:
!    i,j,k  (i<j<k)
!  This correpsonds to a^+_i a^+_j a^+_k |0>

!To compose two excitors, we must move the annihilation operators to the right,
! and the creation operators to the left.  Each switch of operators incurs a sign flip.

!NB excitors and dets are different (there may be a sign change).  See ExcitToDetSign

!First the occ
   do i=0,nIfD
! First mask so we only have 'occupied' orbitals.  The occupieds which are zero
! mean they have been excited from.
      T1=IAND(iLutnI(i),iLutRef(i))
      T2=IAND(iLutnJ(i),iLutRef(i))
      IF(IAND(IEOR(T1,iLutRef(i)),IEOR(T2,iLutRef(i))).NE.0) THEN
!The two excitors are exciting from at least one of the same orbitals.  This
!gives a sign of zero.
!         WRITE(6,'(A,B,B)') 'T1,T2: ',T1,T2
!         WRITE(6,'(A,B)') "Overlapping from",IAND(IEOR(T1,iLutRef(i)),IEOR(T2,iLutRef(i)))
         iSgn=0
         return
      ENDIF
      iLutTmp(i)=IAND(T1,T2)  !Combine the excitors
!We now (very naively, but I can't think of a better way to do it) work out the sign change for the occupied reorder.
!  Tmp lies  between I and J, and will receive the eventual complete excitor.
!  Tmp contains no creation operators as yet, so the right of I can move freely to the left of Tmp
!  J contains a string of iILevel creation operators which must also be taken into account.
! We step through both I and J simultaneously, and move any annihilation operators from either to Tmp,
! with the appropriate sign change.
! iI,iJ and iTmpLevel indicate the number of annihilators left in I J or Tmp.
      T1=IAND(NOT(iLutnI(i)),iLutRef(i))
      T2=IAND(NOT(iLutnJ(i)),iLutRef(i))
      T3=0
      do j=0,end_n_int
!      WRITE(6,'(A,B33.32,B33.32,B33.32)') 'T1,T3,T2: ',T1,T3,T2
         if(BTEST(T1,j)) THEN
            iTmpLevel=iTmpLevel+1
            iI=iI-1
            T1=IBCLR(T1,j)
            T3=IBSET(T3,j)
         else if(BTest(T2,j)) then
!Shift the relevant bit from J to tmp.
            if(btest(iJLevel+iJ+iTmpLevel,0)) iSgn=-iSgn
            iJ=iJ-1
            T2=IBCLR(T2,j)
            T3=IBSET(T3,j)
         endif
!         Write(6,*) j,iTmpLevel,iI,iJ,iSgn
      enddo
   enddo
   iI=iILevel
   iJ=iJLevel 
   iTmpLevel=0
!   WRITE(6,*) "Creation"
   do i=nIfD,0,-1  !Go backwards for the sign
! Now mask so we only have 'virtual' orbitals.  The virtuals which are set mean they have been excited to.
      T3=NOT(iLutRef(i))
      T1=IAND(iLutnI(i),T3)
      T2=IAND(iLutnJ(i),T3)
      IF(IAND(T1,T2).NE.0) THEN
!The two excitors are exciting to at least one of the same orbitals.  This
!gives a sign of zero.
!         WRITE(6,'(A,B)') "Overlapping to",IAND(T1,T2)
         iSgn=0
      ENDIF
      iLutnI(i)=IOR(iLutTmp(i),IOR(T1,T2)) 
!Combine the excitors and the excitor from and put them back into where they should go
!We now (very naively, but I can't think of a better way to do it) work out the sign change for the occupied reorder.
!  Tmp lies  between I and J, and will receive the eventual complete excitor.
!  The right of I can still move freely to the left of Tmp
!  J contains a string of iJ creation operators which must also be taken into account,
!   as well as the iJLevel+iILevel annihilation operators and iTmpLevel creation operators in Tmp
!  totalling  (iILevel+iJLevel)+iI+iTmpLeel
! We step through both I and J simultaneously, and move any annihilation operators from either to Tmp,
! with the appropriate sign change.
! iI,iJ and iTmpLevel indicate the number of annihilators left in I J or Tmp.
      T3=0
      do j=end_n_int,0,-1
!      WRITE(6,'(A,B33.32,B33.32,B33.32)') 'T1,T3,T2: ',T1,IOR(T3,ieor(iLutTmp(i),iLutRef(i))),T2
         if(BTEST(T1,j)) THEN
            iTmpLevel=iTmpLevel+1
            iI=iI-1
            T1=IBCLR(T1,j)
            T3=IBSET(T3,j)
         else if(BTest(T2,j)) then
!Shift the relevant bit from J to tmp. 
            if(btest((iILevel+iJLevel+iJ)+iTmpLevel,0)) iSgn=-iSgn
            iJ=iJ-1
            T2=IBCLR(T2,j)
            T3=IBSET(T3,j)
         endif
!         Write(6,*) j,iTmpLevel,iI,iJ,iSgn
      enddo
   enddo
!   write(6,*) 'x final sign:',iSgn
   return
END SUBROUTINE AddBitExcitor

!Write a Cluster.
   SUBROUTINE WriteCluster(iUnit,C,lTerm)
      use FCIMCData, only: iLutHF
      IMPLICIT NONE
      TYPE(Cluster) C
      INTEGER i
      LOGICAL lTerm
      INTEGER iUnit
      WRITE(iUnit,'(A)',advance='no') '['
      do i=1,C%iSize
         call WriteBitEx(iUnit,iLutHF,C%SelectedExcitors(:,i),.false.)
         write(iUnit,'(A)',advance='no') ','
      enddo
      WRITE(iUnit,'(A)',advance='no') ']'
      if(lTerm) WRITE(iUnit,*)
   END SUBROUTINE WriteCluster
end module CCMCData
