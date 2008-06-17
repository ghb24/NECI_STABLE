!  Do Monte Carlo in pure determinant space.
#include "macros.h"


MODULE MCDets
   use HElem
   IMPLICIT NONE
   logical,parameter    :: tLogExGens=.false.           !Turn on memory logging for excitation generators?  Best not to as slow

!This holds an excitation generator and reference counting associated with it
Type ExcitGen
   integer, pointer     :: exData(:)       !Its excitation generator (allocation and deallocation managed by the calling routine)
   integer                 tagExData    !A memory tag for its excitation generator
   integer                 iRefCount   !A reference count for this excitation generator
End Type ExcitGen
! Setup Particle and Particle list formalism

! Particle holds information about a given collection of subparticles all at the same det

! ParticleData holds the non-array info for each particle
Type ParticleData
   type(HDElement)         hHii        !Its diagonal Hamil element
   integer                 iWeight     !Its weight (the number of subparticles in it)
   integer                 iSgn        !Its sign
   integer                 iParent     !The index of the parent node
   integer                 iBefore     !The child node from this which is less than it 
   integer                 iAfter      !The child node from this which is greater than it
   Type(ExcitGen),pointer::exGen       !A pointer to the excitation generator which manages its own existence with reference counting
   integer                 iProgeny    !  Set to the number of children we have spawned at this particle
End Type
integer, parameter :: ParticleDataSize=HDElementSize+6
integer, parameter :: ParticleDataSizeB=ParticleDataSize*8

! Particle contains all info about each particle
Type Particle
   integer, allocatable :: nI(:)       !The determinant which this corresponds to
   Type(ParticleData)      d
end Type


!ParticleList is a binary tree which holds information about particles.  Adding a particle to the tree will either add a new node, or increment or decrement a current node
Type ParticleList
   integer nParticles                           !Number of nodes in this list
   integer nMaxParticles                        !Max number of nodes in this list
   integer iNextFreeParticle                    !The next position to be allocated
   integer nSubParticles                        !Sum of the number of subparticles at each node
   integer nEl                                  !Number of electrons in a det - used to allocate ParticleDets
   integer, allocatable :: ParticleDets(:,:)    !The Determinant corresponding to each node ! (nel, nMaxParticles)
   integer tagParticleDets
   Type(ParticleData), allocatable :: ParticleData(:)      !All other info for each particle
   integer tagParticleData
   integer iHeadParticle                        !The first node
End Type
contains

! Particle management functions

!Allocate the determinant within a particle
subroutine AllocParticle(P,nEl)
   implicit none
   integer nEl
   Type(Particle) P
   allocate(P%nI(nEl))
end subroutine

!DeAllocate the determinant within a particle
subroutine DeallocParticle(P)
   implicit none
   integer nEl
   Type(Particle) P
   deallocate(P%nI)
end subroutine

!Particle List management functions

!Allocate arrays within a particle list, and clear then
subroutine AllocParticleList(PL,nMaxParticles, nEl)
   Use MemoryManager, only: LogMemAlloc, LogMemDealloc
   implicit none
   Type(ParticleList) PL
   integer nMaxParticles
   integer nEl
   integer ierr
   character(25), parameter :: this_routine='AllocParticleList'
   PL%nEl=nEl
   PL%nMaxParticles=nMaxParticles
   allocate(PL%ParticleDets(nEl,nMaxParticles), stat=ierr)
   LogAlloc(ierr,'ParticleDets',nMaxParticles*nEl,4,PL%tagParticleDets)
   allocate(PL%ParticleData(nMaxParticles), stat=ierr)
   LogAlloc(ierr,'ParticleData',nMaxParticles,ParticleDataSizeB,PL%tagParticleData)
   call ClearParticleList(PL)
end subroutine AllocParticleList

subroutine DeallocParticleList(PL)
   Use MemoryManager, only: LogMemAlloc, LogMemDealloc
   implicit none
   Type(ParticleList) PL
   integer ierr
   character(25), parameter :: this_routine='DeallocParticleList'
   LogDealloc(PL%tagParticleDets)
   Deallocate(PL%ParticleDets)
   LogDealloc(PL%tagParticleData)
   Deallocate(PL%ParticleData)
end subroutine DeallocParticleList

subroutine AddParticle(PL,nI,iWeight,iSgn,hHii,exGen,iPreExisted)
!Add a particle to this list
! First find a node to put it at
   implicit none
   Type(ParticleList) PL
   integer nI(PL%nEl)
   integer iWeight,iSgn
   Type(HDElement) hHii
   Type(ExcitGen), pointer :: ExGen
   integer iNode,iAction
   integer iPreExisted
   iPreExisted=0
   call FindParticleNode(PL,nI,iNode,iAction)
   if(iAction.eq.0) iPreExisted=1
   call AddParticleInt(PL,nI,iWeight,iSgn,hHii,ExGen,iNode,iAction)
end subroutine AddParticle
subroutine AddParticleInt(PL,nI,iWeight,iSgn,hHii,ExGen,iNode,iAction)
!Add a particle to this list
! First find a node to put it at
   implicit none
   Type(ParticleList) PL
   integer nI(PL%nEl)
   integer iWeight,iSgn
   Type(HDElement) hHii
   Type(ExcitGen), pointer :: ExGen
   integer iNode,iAction
   if(iAction.eq.0) then
! This particle matches this node.  Just change its Weight
      PL%nSubParticles=PL%nSubParticles-PL%ParticleData(iNode)%iWeight
      PL%ParticleData(iNode)%iWeight=PL%ParticleData(iNode)%iWeight*PL%ParticleData(iNode)%iSgn+iSgn*iWeight
      if(PL%ParticleData(iNode)%iWeight.lt.0) then
         PL%ParticleData(iNode)%iWeight=-PL%ParticleData(iNode)%iWeight
         PL%ParticleData(iNode)%iSgn=-1
      else
         PL%ParticleData(iNode)%iSgn=1
      endif
      PL%nSubParticles=PL%nSubParticles+PL%ParticleData(iNode)%iWeight
      return
   elseif(iAction.eq.1) then
! This particle is greater than this node, and we're at the end of the line
      if(iNode.ne.0) PL%ParticleData(iNode)%iAfter=PL%iNextFreeParticle
      if(iNode.eq.0) PL%iHeadParticle=PL%iNextFreeParticle
      PL%ParticleData(PL%iNextFreeParticle)%iParent=iNode
      iNode=PL%iNextFreeParticle
      PL%iNextFreeParticle=iNode+1
      PL%nParticles=PL%nParticles+1
      if(PL%iNextFreeParticle.gt.PL%nMaxParticles) call stopgm("Particle List Full")
   elseif(iAction.eq.-1) then
! This particle is less than this node, and we're at the end of the line
      if(iNode.ne.0) PL%ParticleData(iNode)%iBefore=PL%iNextFreeParticle
      if(iNode.eq.0) PL%iHeadParticle=PL%iNextFreeParticle
      PL%ParticleData(PL%iNextFreeParticle)%iParent=iNode
      iNode=PL%iNextFreeParticle
      PL%iNextFreeParticle=iNode+1
      PL%nParticles=PL%nParticles+1
      if(PL%iNextFreeParticle.gt.PL%nMaxParticles) call stopgm("Particle List Full")
   endif   
   PL%ParticleData(iNode)%iWeight=iWeight
   PL%ParticleData(iNode)%iSgn=iSgn
   PL%ParticleData(iNode)%hHii=hHii
   PL%ParticleData(iNode)%ExGen=>ExGen
   call AddReference(ExGen)
   PL%ParticleData(iNode)%iAfter=0
   PL%ParticleData(iNode)%iBefore=0
   PL%ParticleDets(:,iNode)=nI(:)
   PL%nSubParticles=PL%nSubParticles+iWeight
end subroutine AddParticleInt

!Zero a particle List
subroutine ClearParticleList(PL)
   implicit none
   Type(ParticleList) PL
   PL%iNextFreeParticle=1
   PL%nSubParticles=0
   PL%nParticles=0
   PL%iHeadParticle=0
   PL%ParticleDets(:,:)=0
end subroutine ClearParticleList

!Get the next particle with non-zero weight.  We traverse the tree as follows.  
!  Start at HEAD node
!    The next node after CURRENT is
!     iBefore (if its weight is non-zero)
!     iAfter (if its weight is non-zero)
!     Recursing up the tree, the first non-zero iAfter, otherwise 0

!The Int version does not do a copy to P
!  If there are no more particles, it returns iParticle=0
subroutine GetNextParticleInt(PL,iParticle)
   implicit none
   Type(ParticleList) PL
   integer iParticle
   integer iOldPart
   if(iParticle.eq.0) then
! Just get the head if there is one
      if(PL%nParticles.gt.0) then
         iParticle=PL%iHeadParticle
      endif
   else
      if(PL%ParticleData(iParticle)%iBefore.ne.0) then
         iParticle=PL%ParticleData(iParticle)%iBefore
      else
! We need to look upwards in the tree
         iOldPart=0
         do while(iParticle.gt.0.and.PL%ParticleData(iParticle)%iAfter.eq.iOldPart)
!            WRITE(6,*) iOldPart, iParticle
            iOldPart=iParticle
            iParticle=PL%ParticleData(iParticle)%iParent
            if(PL%ParticleData(iParticle)%iBefore.eq.iOldPart) iOldPart=0
!            WRITE(6,*) iOldPart, iParticle
         enddo
         if(iParticle.gt.0) iParticle=PL%ParticleData(iParticle)%iAfter
      endif
   endif
end subroutine GetNextParticleInt
! GetnextParticle Gets the particle after iParticle.  If iParticle=0, it returns the head node.
!  It puts the Particle info into P.  If there are no more particles, it returns iParticle=0
subroutine GetNextParticle(PL,iParticle,P)
   implicit none
   Type(ParticleList) PL
   Type(Particle) P
   integer iParticle
   call GetNextParticleInt(PL,iParticle)
   call GetParticle(PL,iParticle,P)
end subroutine GetNextParticle

!Extract the data from a particle
subroutine GetParticle(PL,iParticle,P)
   implicit none
   Type(ParticleList) PL
   Type(Particle) P
   integer iParticle
   if(iParticle.ne.0) then
      P%d=PL%ParticleData(iParticle)
      P%nI(:)=PL%ParticleDets(:,iParticle)
   endif
end subroutine GetParticle

!run GetNextParticle until we find one with zero weight
subroutine GetNextZeroParticle(PL,iParticle,P)
   implicit none
   Type(ParticleList) PL
   Type(Particle) P
   integer iParticle
   call GetNextParticleInt(PL,iParticle)
   do while(iParticle.gt.0.and.PL%ParticleData(iParticle)%iWeight.ne.0)
      call GetNextParticleInt(PL,iParticle)
   enddo
   call GetParticle(PL,iParticle,P)
end Subroutine GetNextZeroParticle

! Remove a particle from the node list, and fix up the list
subroutine DeleteParticle(PL,iParticle)
   implicit none
   Type(ParticleList) PL
   integer iParticle
   integer iParent
   integer iHanging

!  We remove this node first
!   If the particle has a parent
   iParent=PL%ParticleData(iParticle)%iParent
   if(iParent.ne.0) then
      if(PL%ParticleData(iParent)%iBefore.eq.iParticle) then
!  We're a before, so we link in our own before, and leave the iAfter hanging
         PL%ParticleData(iParent)%iBefore=PL%ParticleData(iParticle)%iBefore
         PL%ParticleData(PL%ParticleData(iParticle)%iBefore)%iParent=iParent
         iHanging=PL%ParticleData(iParticle)%iAfter
      else
!  We're an after, so we link in our own after, and leave the iBefore hanging
         PL%ParticleData(iParent)%iAfter=PL%ParticleData(iParticle)%iAfter
         PL%ParticleData(PL%ParticleData(iParticle)%iAfter)%iParent=iParent
         iHanging=PL%ParticleData(iParticle)%iBefore
      endif
   else
!  We're the parent, so we make our before the headif we can, and leave the after hanging
      PL%iHeadParticle=PL%ParticleData(iParticle)%iBefore
      iHanging=PL%ParticleData(iParticle)%iAfter
      if(PL%iHeadParticle.eq.0) then
!  We don't have a before, so try our after
         PL%iHeadParticle=PL%ParticleData(iParticle)%iAfter
         PL%ParticleData(iParticle)%iParent=iParent
         iHanging=0
      else
         PL%ParticleData(PL%iHeadParticle)%iParent=0
      endif
      if(PL%iHeadParticle.eq.0) then
! no before or after, so we must be top of the list
         PL%iHeadParticle=0
      else
         PL%ParticleData(PL%iHeadParticle)%iParent=0
      endif
   endif
!  Deal with the hanging
   if(iHanging) then
      call AddParticleNode(PL,iHanging)
   endif
!  Now note that we are zero
   PL%ParticleDets(1,iParticle)=0
! consolidate the particle list a little
   do while(PL%ParticleDets(1,PL%iNextFreeParticle).eq.0)
      PL%iNextFreeParticle=PL%iNextFreeParticle-1
   enddo
end subroutine DeleteParticle

subroutine AddParticleNode(PL,iParticle)
   implicit none
   Type(ParticleList) PL
   integer iParticle
   integer iNode,iAction
   integer iPreExisted
   call FindParticleNode(PL,PL%ParticleDets(:,iParticle),iNode,iAction)
   if(iAction.eq.0) stop 'Found duplicate particle in list'
   if(iAction.eq.1) then
! This particle is greater than this node, and we're at the end of the line
      if(iNode.ne.0) PL%ParticleData(iNode)%iAfter=iParticle
      if(iNode.eq.0) PL%iHeadParticle=iParticle
      PL%ParticleData(iParticle)%iParent=iNode
      iNode=iParticle
   elseif(iAction.eq.-1) then
! This particle is less than this node, and we're at the end of the line
      if(iNode.ne.0) PL%ParticleData(iNode)%iBefore=iParticle
      if(iNode.eq.0) PL%iHeadParticle=iParticle
      PL%ParticleData(iParticle)%iParent=iNode
      iNode=iParticle
   endif   
end subroutine AddParticleNode

!returns iNode as the node nI is in if already there (and Action=0).  Otherwise Action=-1 for needing a new node before the returned iNode, or +1 for after.
subroutine FindParticleNode(PL,nI,iNode,iAction)
   implicit none
   Type(ParticleList) PL
   integer nI(*)
   integer iNode
   integer iAction,iCmp
   integer iCmpDets  !external function

   iNode=PL%iHeadParticle
   if(iNode.eq.0) then
      iAction=-1
      return 
   endif
   iCmp=iCmpDets(nI,PL%ParticleDets(:,iNode),PL%nEl)
   iAction=0
   do while(iCmp.ne.0)
      if(iCmp.lt.0) then
!our det is less than the current, so look at before
         if(PL%ParticleData(iNode)%iBefore.eq.0) then
            iAction=-1
            iCmp=0
         else
            iNode= PL%ParticleData(iNode)%iBefore
         endif
      else
!our det is more than the current, so look at after
         if(PL%ParticleData(iNode)%iAfter.eq.0) then
            iAction=1
            iCmp=0
         else
            iNode= PL%ParticleData(iNode)%iAfter
         endif
      endif
      if(iCmp.ne.0)    iCmp=iCmpDets(nI,PL%ParticleDets(:,iNode),PL%nEl)
   enddo
   return
end subroutine FindParticleNode

subroutine DumpParticleList(iunit,PL)
   implicit none
   Type(ParticleList) PL
   integer iUnit
   Type(Particle) P
   integer iParticle
   iParticle=0
   Allocate(P%nI(PL%nEl))
   call GetNextParticle(PL,iParticle,P)
   do while (iParticle.ne.0)
      write(iUnit,"(I,$)") iParticle
      call WriteDet(iUnit,P%nI,PL%nEl,.false.)
      write(iUnit,"(7I)") P%d%iWeight,P%d%iSgn,P%d%iParent, P%d%iBefore,P%d%iAfter,P%d%ExGen%tagExData,P%d%ExGen%iRefCount
      call GetNextParticle(PL,iParticle,P)
   enddo
   Deallocate(P%nI)
end subroutine DumpParticleList
   
subroutine AllocateExGen(exGen,nI,nEl,this_routine)
   Use MemoryManager, only: LogMemAlloc, LogMemDealloc
   Type(ExcitGen), pointer :: exGen
   integer nEl 
   integer nI(nEl),nJ(nEl)
   integer Store(6)  !Used for Excit gen
   integer iLen,iC     
   integer ierr
   character(*) this_routine


   Allocate(exGen)
   exGen%iRefCount=0

!.. Setup the spin excit generator
   Store(1)=0
   ! .true. for
   ! 3 for both singles and doubles. (Brill theorem is taken into account if needed however)
   call GenSymExcitIt3(nI,.true.,iLen,nJ,iC,0,Store,3)
   ExGen%tagExData=0
   Allocate(exGen%exData(iLen),STAT=ierr)
   if(tLogExGens) LogAlloc(ierr,'exGen',iLen,8,exGen%tagExData)
   ExGen%exData(1)=0
!Now setup the generator (Store contains some necessary info)
   call GenSymExcitIt3(nI,.true.,exGen%exData,nJ,iC,0,Store,3)
   return
end subroutine AllocateExGen
   
subroutine DeallocateExGen(exGen,this_routine)
   Use MemoryManager, only: LogMemAlloc, LogMemDealloc
   Type(ExcitGen) exGen
   integer ierr
   character(*) this_routine
   if(exGen%iRefCount.ne.0) then
      stop 'Attempting to deallocate referenced excitation generator'
   endif
!   write(6,*) "Dealloc",ExGen%tagExData
   if(tLogExGens) LogDealloc(ExGen%tagExData)
   Deallocate(ExGen%ExData)
end subroutine DeallocateExGen

subroutine AddReference(exGen)
   Type(ExcitGen) exGen
   exGen%iRefCount=exGen%iRefCount+1
end subroutine

subroutine DelReference(exGen)
   Type(ExcitGen) exGen
   exGen%iRefCount=exGen%iRefCount-1
end subroutine
   


subroutine MCDetsCalc(nI,iSeed,nCycles,dTau,dMu,nMaxParticles,nInitParticles,iStep,dInitShift)
   Use HElem
   Use MemoryManager, only: LogMemAlloc, LogMemDealloc
   Use Determinants, only: GetHElement3
   Use System, only: BasisFN,nEl
   IMPLICIT NONE
   character(25), parameter :: this_routine='MCDets'
   integer nI(nEl)   !The root determinant
   integer iSeed     !Random number seed
   integer nCycles   !Number of cycles we should do
   Type(ExcitGen), pointer :: exGen
   integer nJ(nEl)   !Store generated determinants
   real*8  pGen      !Store generation prob of dets
   integer iC        ! The level of excitation returned
   integer ierr      !errors in memory allocation
   integer iCycle    !Cycle counter
   Type(Helement) hHij
   Type(Helement) hRhIJ  
   Type(Helement) hRhJJ  
 
   Type(HDElement) EShift, EHF,hHjj
   integer iParticle

   
   Type(ParticleList) PL(2)      ! 2 particle lists - one for past cycle and one created this cycle
   Type(Particle) P
   
   integer iPL,nPL
   integer iSubPart,iNewCount
   integer nMaxParticles,isgn,iNewWeight
   real*8 r,rnum,dtau,dMu
   real*8 RAN2 ! external
   integer iPreExisted
   integer nInitParticles
   integer iOldCount
   integer iStep
   integer iCount
   integer iNode,iAction         !Internal loop variables documented later
   integer ICMin,ICMax           ! Min and Max excitation level we encounter in a step
   real*8 dInitShift
   real*8 ICMean                 !Mean excit level 
   Type(HDElement) ECumlTot,ECumlTotSq ! Used to give stats on the energy fluctuation
   integer nECumls               ! "
   integer iTot
!   dTau=1e-2
!   dMu=1e-1
!   nMaxParticles=1000
!   nInitParticles=1000
!   iStep=10
   
   ECumlTot=0.d0
   nECumls=0
   call writeDet(6,nI,nEl,.true.)
! We need to allocate space for our current particle
   call AllocParticle(P,nEl)   
   
   EShift=dInitShift
   EHF=GetHElement3(nI,nI,0)
   EShift=EHF+EShift

   call AllocateExGen(exGen,nI,nEl,this_routine)

! We switch the current particle list between 1 and 2  iPL is current.  nPL is the one we're creating
   iPL=1
   nPL=2
!  Allocate particle lists in PL
   call AllocParticleList(PL(1),nMaxParticles,nEl)
   call AllocParticleList(PL(2),nMaxParticles,nEl)
!Setup single particle with weight nInitParticles at the HF det
   call AddParticle(PL(iPL),nI,nInitParticles,+1,EHF,exGen,iPreExisted)

   write(6,*) "nCycles:",nCycles
   iOldCount=nInitParticles

   do iCycle=1,nCycles
      call ClearParticleList(PL(nPL))
      call SpawnParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P)
      call PropagateParticles(PL,iPL,nPL,dTau,EShift,ICmean,ICmin,ICMax,nI,nEl,iSeed,ExGen,P)
      call CleanupParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P)
!      write(6,*) "Cycle",iCycle
!      call DumpParticleList(6,PL(nPL))
         

      if(PL(nPL)%nSubParticles.eq.0) then
         write(6,*) "All particles have died."
         return
      endif
      if(mod(iCycle,iStep).eq.0) then
         EShift=EShift+HDElement(dMu*log((iOldCount+0.d0)/PL(nPL)%nSubParticles)/(dTau*iStep))
         iOldCount=PL(nPL)%nSubParticles
         WRITE(6,'(2I9,G20.12,F10.5,2I7)') iCycle, PL(nPL)%nSubParticles,EShift,ICMean,ICMin,ICMax
         WRITE(76,"(2I,2G,2I)") iCycle, PL(nPL)%nSubParticles,DReal(EShift),ICMean,ICMin,ICMax
         call flush(76)
         if(mod(iCycle,iStep*10).eq.0) then
            ECumlTot=ECumlTot+EShift
            ECumlTotSq=ECumlTotSq+EShift*EShift
            nECumls=nECumls+1
         endif
      endif
! Switch to the other particle list
      nPL=iPL
      iPL=3-iPL
   end do !iCycle

   write(6,*) "Average EShift:", DReal(ECumlTot)/nECumls
   write(6,*) "Std of EShift:", sqrt((DReal(ECumlTotSq)/nECumls-(DReal(ECumlTot)/nECumls)**2)/(nECumls-1))
! We need to deallocate all the excitation generators in the last particle list we created 
   iParticle=0
   call GetNextParticle(PL(iPL),iParticle,P)
   do while (iParticle.gt.0)
      call DeallocateExGen(P%d%ExGen,this_routine)
      call GetNextParticle(PL(iPL),iParticle,P)
   enddo
   call DeallocParticleList(PL(1))
   call DeallocParticleList(PL(2))
   call DeallocParticle(P)
end subroutine MCDetsCalc

subroutine PropagateParticles(PL,iPL,nPL,dTau,EShift,ICmean,ICmin,ICMax,nI,nEl,iSeed,ExGen,P)
   Type(ParticleList) PL(2)
   integer iPL,nPL
   integer icmin,icmax
   real*8 icmean
   real*8 dTau
   Type(HDElement) EShift
   integer iSeed
   Type(ExcitGen),pointer :: ExGen
   Type(Particle) P

   integer iParticle
   integer iC 
   integer iNewCount,iNewWeight
   integer iSubPart
   integer nEl,nI(nEl)
   real*8 r,rnum
   integer iPreExisted
   
   integer iGetExcitLevel
   real*8 RAN2 ! external
   character(25), parameter :: this_routine='PropagateParticles'

   iParticle=0 !Just get the head node first and its info in P
   call GetNextParticle(PL(iPL),iParticle,P)
   ICMean=0.D0
   ICMin=10
   ICMax=0
   do while (iParticle.gt.0)  !i.e. while we're at a valid particle
! Each Particle has iWeight of subparticles at the same det
!         write(6,*) iParticle,P%d%exGen%tagExData
!         write(56,'(I,$)')  P%d%iWeight

      IC=iGetExcitLevel(Ni,P%Ni,NEl)
      ICmean=ICmean+IC*(P%d%iWeight)
      IF(IC.lt.ICMin) ICMin=IC
      IF(IC.gt.ICMax) ICMax=IC
      
      iNewCount=0
      do iSubPart=1,P%d%iWeight
! With this current subparticle we do two processes.  Firstly we deplete the current particle according to 
! Hii and tau
         r=DREAL(P%d%hHii-EShift)*dTau !Prob of death
         r=1-r ! Prob of survival
         iNewWeight=int(r)  ! Number of newly created definitely
         r=r-iNewWeight
!            if(r.lt.0) THEN write(6,*) "Warning: Self-destruction  probabililty <0"
! This particle doesn't self-destruct if a random num is > r
         rnum=ran2(iSeed)
         if(r.gt.rnum) then
!  Add it to the new list
            inewWeight=iNewWeight+1
         endif
         if(iNewWeight.gt.0) then
!Copies accross to the nPL list (if actually creating more weight at particle)
            call AddParticle(PL(nPL),P%nI,iNewWeight,+P%d%iSgn,P%d%hHii,P%d%exGen,iPreExisted)
!Number of subparticles of this particle which have are present at end of creation/death process
            iNewCount=iNewCount+iNewWeight
         endif
      enddo !iSubPart
      call DelReference(P%d%exGen)
!Save the number of children we've spawned at this det
      PL(iPL)%ParticleData(iParticle)%iProgeny=iNewCount
      call GetNextParticle(PL(iPL),iParticle,P)
   enddo !iParticle
!      call DumpParticleList(6,PL(nPL))
   ICMean=ICMean/(Pl(NPl)%NSubparticles)

end Subroutine PropagateParticles
subroutine SpawnParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P)
   Use Determinants, only: GetHElement3
   Type(ParticleList) PL(2)
   integer iPL,nPL
   integer nEl
   integer iSeed
   real*8 dTau
   Type(ExcitGen),pointer :: exGen
   Type(Particle) P

   integer iParticle
   integer iSubPart
   integer iAction
   Type(HElement) hHij
   Type(HDelement) HHjj
   integer nJ(nEl)
   integer iC,iCount
   integer iNode
   integer iSgn
   real*8 pGen
   real*8 r,rnum

   real*8 ran2
   character(25), parameter :: this_routine='SpawnParticles'

!      write(56,*)
! Secondly we might create positive or negative particle at an excitation
   iParticle=0
   call GetNextParticle(PL(iPL),iParticle,P)
   do while (iParticle.gt.0)
!         write(6,*) iParticle,P%d%exGen%tagExData
!         call WriteDet(6,P%nI,nEl,.false.)
!         write(6,*) P%d%iProgeny
      
      do iSubPart=1,P%d%iWeight
! generate a random excitation of the root
         call GenRandSymExcitIt3(P%nI,P%d%ExGen%exData,nJ,iSeed,iC,0,pGen,iCount)
! The excitation is in nJ, and its gen prob in pGen.
!  Just print out for now
!            write(6,*)  pGen
         hHij=GetHElement3(P%nI,nJ,iC)

! see if we create a particle in it
         
         r=DREAL(abs(hHij))*dTau/pGen
!            WRITE(6,*) r,hHij
         rnum=ran2(iSeed)
! This a new particle is createdif r > a random num
         if(r.gt.1) write(6,*) "Warning: Child creation probabililty>1"
         if(r.gt.rnum) then
!  See if this particle's already in the list.  iNode is the node it's in in PL(nPL) if it exists, or a new clean node ready to receive it.  iAction is set if 
            call FindParticleNode(PL(nPL),nJ,iNode,iAction)
            if(iAction.ne.0) then
!  If it doesn't already exist,
!  Now we have to get info about this particle
               hHjj=GetHElement3(nJ,nJ,0)
               call AllocateExGen(exGen,nJ,nEl,this_routine)
            else
! Already exists, so we don't recalculate info
            endif
!  Add it to the new list
!   If the Hij element is -ve then keep the same sign, otherwise flip the sign
            iSgn=-P%d%iSgn
            if(DREAL(hHij).lt.0) iSgn=-iSgn
!  Adds on the weight if iAction=0 or creates a new node w.r.t iNode and iAction if iAction +-1
            call AddParticleInt(PL(nPL),nJ,1,iSgn,hHjj,ExGen,iNode,iAction)
         endif
      enddo !iSubPart
      call GetNextParticle(PL(iPL),iParticle,P)
   enddo !iParticle
end subroutine SpawnParticles
subroutine CleanupParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P)
   Use Determinants, only: GetHElement3
   Type(ParticleList) PL(2)
   integer iPL,nPL
   integer nEl
   integer iSeed
   real*8 dTau
   Type(ExcitGen),pointer :: exGen
   Type(Particle) P

   integer iParticle
   integer iSubPart
   character(25), parameter :: this_routine='SpawnParticles'
   iParticle=0
   call GetNextParticle(PL(iPL),iParticle,P)
   do while (iParticle.gt.0)
      if(P%d%exGen%iRefCount.eq.0) then
!            WRITE(6,*) "Progeny 0",iParticle
         call DeallocateExGen(P%d%ExGen,this_routine)
      endif
! Just see if we need to keep our old excitation generator around, otherwise kill it
      call GetNextParticle(PL(iPL),iParticle,P)
   enddo !iParticle
end subroutine CleanupParticles
end Module MCDets
