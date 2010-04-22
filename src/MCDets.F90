! TO DO
! Separate routines for death/spawning to allow multiple growth/death
! Sort output for consistency
!** Add optional energy calculations (both types) - store IC for each particle (print mean & max/min)
! Impliment culling
!** Impliment fixed-sign approx (extra memory slots for H0 and Hij)
! Updatesft in new routine
!** Net positive particles
!** Read & Dump pops

!  Do Monte Carlo in pure determinant space.
#include "macros.h"


MODULE MCDets
!   use constants, only: dp
!   IMPLICIT NONE
!   logical,parameter    :: tLogExGens=.false.           !Turn on memory logging for excitation generators?  Best not to as slow
!
!!This holds an excitation generator and reference counting associated with it
!Type ExcitGen
!   integer, pointer     :: exData(:)       !Its excitation generator (allocation and deallocation managed by the calling routine)
!   integer                 tagExData    !A memory tag for its excitation generator
!   integer                 iRefCount   !A reference count for this excitation generator. We only want to deallocate the excitgen if
!!all particles have finished using it. If a node is killed, it can still be present in the new particle list, since other excitations
!!may have generated it in the current cycle. Therefore we want to keep the excitation generator.
!
!End Type ExcitGen
!! Setup Particle and Particle list formalism
!
!! Particle holds information about a given collection of subparticles all at the same det
!
!! ParticleData holds the non-array info for each particle
!Type ParticleData
!   real(dp)         hHii        !Its diagonal Hamil element
!   integer                 iWeight     !Its weight (the number of subparticles in it)
!   integer                 iSgn        !Its sign
!   integer                 IC0         !Its excitation level away from HF determinant
!   HElement_t          hHi0        !Its Hamiltonian matrix element between itself and HF determinant
!   integer                 iParent     !The index of the parent node
!   integer                 iBefore     !The child node from this which is less than it 
!   integer                 iAfter      !The child node from this which is greater than it
!   Type(ExcitGen),pointer::exGen       !A pointer to the excitation generator which manages its own existence with reference counting
!   integer                 iProgeny    !  Set to the number of children we have spawned at this particle
!End Type
!integer, parameter :: ParticleDataSize=Size+HElement_t_size+7
!integer, parameter :: ParticleDataSizeB=ParticleDataSize*8
!
!! Particle contains all info about each particle
!Type Particle
!   integer, allocatable :: nI(:)       !The determinant which this corresponds to
!   Type(ParticleData)      d
!end Type
!
!
!!ParticleList is a binary tree which holds information about particles.  Adding a particle to the tree will either add a new node, or increment or decrement a current node
!Type ParticleList
!   integer nParticles                           !Number of nodes in this list
!   integer nMaxParticles                        !Max number of nodes in this list
!   integer iNextFreeParticle                    !The next position to be allocated
!   integer nSubParticles                        !Sum of the number of subparticles at each node
!   integer nEl                                  !Number of electrons in a det - used to allocate ParticleDets
!   integer, allocatable :: ParticleDets(:,:)    !The Determinant corresponding to each node ! (nel, nMaxParticles)
!   integer tagParticleDets
!   Type(ParticleData), allocatable :: ParticleData(:)      !All other info for each particle
!   integer tagParticleData
!   integer iHeadParticle                        !The first node
!End Type

!Summary of four types:
!ExcitGen - holds excitgen and reference counting associated with it.
!ParticleData - non-array info for each particle
!ParticleList - Contains the whole binary tree - all particles (including particleData type and determinant for each one), and details about the overall tree
!Particle - holds data for a single particle - a single determinant and its data
contains

!subroutine MCDetsCalc(nI,iSeed,nCycles,dTau,dMu,nMaxParticles,nInitParticles,iStep,dInitShift,GrowMaxFactor,CullFactor)
!   use constants, only: dp
!   Use global_utilities
!   Use Determinants, only: GetHElement3
!   use SystemData, only: BasisFN,nEl
!   IMPLICIT NONE
!   character(25), parameter :: this_routine='MCDets'
!   integer nI(nEl)   !The root determinant
!   integer iSeed     !Random number seed
!   integer nCycles   !Number of cycles we should do
!   Type(ExcitGen), pointer :: exGen
!   integer nJ(nEl)   !Store generated determinants
!   real*8  pGen      !Store generation prob of dets
!   integer iC        ! The level of excitation returned
!   integer ierr      !errors in memory allocation
!   integer iCycle    !Cycle counter
!   HElement_t hHij
!   HElement_t hRhIJ  
!   HElement_t hRhJJ  
! 
!   real(dp) EShift, EHF,hHjj
!   integer iParticle
!
!   
!   Type(ParticleList) PL(2)      ! 2 particle lists - one for past cycle and one created this cycle
!   Type(Particle) P
!   
!   integer iPL,nPL
!   integer iSubPart,iNewCount
!   integer nMaxParticles,isgn,iNewWeight
!   real*8 r,rnum,dtau,dMu,GrowRate
!   real*8 RAN2 ! external
!   integer iPreExisted
!   integer nInitParticles
!   integer iOldCount
!   integer iStep
!   integer iCount
!   integer NoCulls               !The number of culls/growths in a given shift cycle
!!CullInfo is the number of walkers before and after the cull (columns 1&2), and the third element is the previous number of steps before this cull...
!!Only 15 culls/growth increases are allowed in a given shift cycle
!   integer CullInfo(15,3)   
!   real*8 GrowMaxFactor,CullFactor      !Info for the culling routine
!   integer iNode,iAction         !Internal loop variables documented later
!   integer ICMin,ICMax           ! Min and Max excitation level we encounter in a step
!   real*8 dInitShift
!   real*8 ICMean                 !Mean excit level 
!   integer iTot
!!   dTau=1e-2
!!   dMu=1e-1
!!   nMaxParticles=1000
!!   nInitParticles=1000
!!   iStep=10

!    CALL Stop_All("MCDetsCalc","This code has been commented out.")

!    NoCulls=0
!    CullInfo=0
!
!    OPEN(15,file='FCIMCStats',status='unknown')     !Open same file as FCIMC.F90 for consistency
!
!    WRITE(6,*) ""
!    WRITE(6,*) "Performing MC Dets...."
!    write(6,*) "nCycles:",nCycles
!    WRITE(6,*) "Initial number of walkers chosen to be: ", nInitParticles
!    WRITE(6,*) "Damping parameter for Diag Shift set to: ", dMu
!    WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", dInitShift
!    WRITE(6,*) "       Step  Shift  ParticleChange  GrowRate  TotParticles        Proj.E      Net+veWalk     Proj.E-Inst"
!    WRITE(15,*) "#       Step  Shift  ParticleChange  GrowRate  TotParticles        Proj.E      Net+veWalk     Proj.E-Inst"
!
!   
!!   call writeDet(6,nI,nEl,.true.)
!! We need to allocate space for our current particle
!   call AllocParticle(P,nEl)   
!   
!   EShift=dInitShift
!   EHF=GetHElement3(nI,nI,0)
!   EShift=EHF+EShift
!
!!   WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
!!   WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7,I9,G16.7)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE,TotWalkers,ProjectionEInst
!
!   call AllocateExGen(exGen,nI,nEl,this_routine)
!
!! We switch the current particle list between 1 and 2  iPL is current.  nPL is the one we're creating
!   iPL=1
!   nPL=2
!!  Allocate particle lists in PL
!   call AllocParticleList(PL(1),nMaxParticles,nEl)
!   call AllocParticleList(PL(2),nMaxParticles,nEl)
!!Setup single particle with weight nInitParticles at the HF det - Current PL, FDet, Initial Particles, Sign, HF Energy (diag elem), Exitgen for nI, Returned value whether node was in list.
!   call AddParticle(PL(iPL),nI,nInitParticles,+1,EHF,(EHF),0,exGen,iPreExisted)
!
!   iOldCount=nInitParticles
!
!   do iCycle=1,nCycles
!      call ClearParticleList(PL(nPL))   !Zero list contained in nPL
!      call SpawnParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P,nI,EHF)   !Run through all subparticles and determine if they want to spawn weight to neighbouring particles
!      call PropagateParticles(PL,iPL,nPL,dTau,EShift,ICmean,ICmin,ICMax,nI,nEl,iSeed,ExGen,P)   !Attempt to kill all subparticles
!      call CleanupParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P)      !Search through all particles and deallocate exgens if no longer used
!!      write(6,*) "Cycle",iCycle
!!      call DumpParticleList(6,PL(nPL))
!
!      if(PL(nPL)%nSubParticles.eq.0) then
!         write(6,*) "All particles have died."
!         return
!      endif
!
!      IF(PL(nPL)%nSubParticles.gt.(nInitParticles*GrowMaxFactor)) THEN
!!Particle number has grown too large - cull particles randomly
!         WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to over ",GrowMaxFactor," times initial number..."
!         call CullParticles(PL(nPL),.true.,iSeed,NoCulls,CullInfo,iCycle,iStep,CullFactor)
!
!      ELSEIF(PL(nPL)%nSubParticles.lt.(nInitParticles/2)) THEN
!!Particle number has dropped too low - double all particles
!         WRITE(6,*) "Total number of particles has dropped to half initial number - doubling all particles..."
!         call CullParticles(PL(nPL),.false.,iSeed,NoCulls,CullInfo,iCycle,iStep,CullFactor)
!      ENDIF
!
!      if(mod(iCycle,iStep).eq.0) then
!!Update shift & print out info
!
!         call UpdateShift(PL(nPL),EShift,dMu,dTau,iStep,CullInfo,NoCulls,iOldCount,GrowRate)
!
!         WRITE(6,'(I9,G20.12,I10,F10.5,I12)') iCycle, DREAL(EShift), PL(nPL)%nSubParticles-iOldCount,GrowRate,PL(nPL)%nSubParticles
!         WRITE(15,'(I9,G20.12,I10,F10.5,I12)') iCycle, DREAL(EShift), PL(nPL)%nSubParticles-iOldCount,GrowRate,PL(nPL)%nSubParticles
!         call flush(15)
!         call flush(6)
!         iOldCount=PL(nPL)%nSubParticles
!
!      endif
!
!! Switch the current and new particle lists around
!      nPL=iPL
!      iPL=3-iPL
!   end do !iCycle
!
!!   write(6,*) "Average EShift:", DReal(ECumlTot)/nECumls
!!   write(6,*) "Std of EShift:", sqrt((DReal(ECumlTotSq)/nECumls-(DReal(ECumlTot)/nECumls)**2)/(nECumls-1))
!
!!Final tidy...
!! We need to deallocate all the excitation generators in the last particle list we created 
!   iParticle=0
!   call GetNextParticle(PL(iPL),iParticle,P)
!   do while (iParticle.gt.0)
!      P%d%ExGen%iRefCount=0     !Need to set all references to zero before it will let you deallocate
!      call DeallocateExGen(P%d%ExGen,this_routine)
!      call GetNextParticle(PL(iPL),iParticle,P)
!   enddo
!   call DeallocParticleList(PL(1))
!   call DeallocParticleList(PL(2))
!   call DeallocParticle(P)
!
!   CLOSE(15)

!end subroutine MCDetsCalc


!!Routine to update the shift
!subroutine UpdateShift(PL,EShift,dMu,dTau,iStep,CullInfo,NoCulls,iOldCount,GrowRate)
!   Type(ParticleList) PL
!   integer iStep,NoCulls,iOldCount,GrowthSteps,j,k
!   real*8 dMu,dTau,GrowRate
!   real(dp) EShift
!!CullInfo is the number of walkers before and after the cull (columns 1&2), and the third element is the previous number of steps before this cull...
!!Only 15 culls/growth increases are allowed in a given shift cycle
!   integer CullInfo(15,3)   
!
!   IF(NoCulls.eq.0) THEN
!       GrowRate=(PL%nSubParticles+0.D0)/(iOldCount+0.D0)
!   ELSE
!!We have had culling in this update of the shift which needs to be taken into account
!       GrowRate=((CullInfo(1,3)+0.D0)/(iStep+0.D0))*((CullInfo(1,1)+0.D0)/(iOldCount+0.D0))
!       do j=2,NoCulls
!
!!This is needed since the steps between culls are stored cumulativly
!           GrowthSteps=CullInfo(j,3)-CullInfo(j-1,3)
!           GrowRate=GrowRate+((GrowthSteps+0.D0)/(iStep+0.D0))*((CullInfo(j,1)+0.D0)/(CullInfo(j-1,2)+0.D0))
!
!       enddo
!
!       GrowthSteps=iStep-CullInfo(NoCulls,3)
!       GrowRate=GrowRate+((GrowthSteps+0.D0)/(iStep+0.D0))*((PL%nSubParticles+0.D0)/(CullInfo(NoCulls,2)+0.D0))
!
!       NoCulls=0
!       CullInfo=0
!
!   ENDIF
!
!   EShift=EShift-(dMu*log(GrowRate)/(dTau*iStep))
!
!   RETURN
!
!end subroutine UpdateShift
!   
!
!!The number of particles has grown/shrunk too much - reduce/double the number of subparticles
!!This is done by running through all particles, and stochastically reducing/doubling their weights
!subroutine CullParticles(PL,tHighLow,iSeed,NoCulls,CullInfo,iCycle,iStep,CullFactor)
!   Type(ParticleList) PL
!   Logical tHighLow     !Indicates whether we are culling, or doubling weights
!   integer iCycle,iStep    !MC Cycle we are on, and the frequency with which the shift is updated
!   integer iSeed,iParticle,Culled,ToCull,NoCulls
!!CullInfo is the number of walkers before and after the cull (columns 1&2), and the third element is the previous number of steps before this cull...
!!Only 15 culls/growth increases are allowed in a given shift cycle
!   integer CullInfo(15,3)   
!   real*8 r,rnum,ran2,CullFactor
!   type(Particle) P
!
!   NoCulls=NoCulls+1      !Log the fact we have made a cull
!   IF(NoCulls.gt.15) CALL Stop_All("CullParticles","Too many culls/growths in a given shift cycle")
!     
!!CullInfo(:,1) is total subparticles before cull
!   CullInfo(NoCulls,1)=PL%nSubParticles
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!   CullInfo(NoCulls,3)=mod(iCycle,iStep)
!
!   Culled=0
!   iParticle=0 !Just get the head node first and its info in P
!   call GetNextParticle(PL,iParticle,P)
!
!   do while (iParticle.gt.0)    !i.e. while we're at a valid particle
!
!       IF(tHighLow) THEN
!!There are too many particles - cull all weights with CullFactor
!           r=(P%d%iWeight+0.D0)/CullFactor
!           ToCull=INT(r)
!           r=r-(ToCull+0.D0)
!           rnum=ran2(iSeed)
!           IF(r.gt.rnum) THEN
!!Cull an additional subparticle
!               ToCull=ToCull+1
!           ENDIF
!           
!           PL%ParticleData(iParticle)%iWeight=PL%ParticleData(iParticle)%iWeight-ToCull
!
!           Culled=Culled+ToCull
!!           WRITE(6,*) iParticle,P%d%iWeight,PL%ParticleData(iParticle)%iWeight
!
!       ELSE
!!There are too few particles - non-stochastically double the weight on every particle
!
!           PL%ParticleData(iParticle)%iWeight=PL%ParticleData(iParticle)%iWeight*2
!
!       ENDIF
!
!       call GetNextParticle(PL,iParticle,P)
!
!   enddo
!
!   IF(tHighLow) THEN
!       WRITE(6,*) "Total particle number culled from ", PL%nSubParticles," by ", Culled," particles."
!       PL%nSubParticles=PL%nSubParticles-Culled
!   ELSE
!       WRITE(6,*) "Total particle number increased from ", PL%nSubParticles," to ",PL%nSubParticles*2
!       PL%nSubParticles=PL%nSubParticles*2
!   ENDIF
!   
!   CullInfo(NoCulls,2)=PL%nSubParticles     !The second column of CullInfo contains info about populations after culling
!   
!   RETURN
!
!end subroutine CullParticles
!
!
!!Run through all particles, and attempt to kill them in turn
!subroutine PropagateParticles(PL,iPL,nPL,dTau,EShift,ICmean,ICmin,ICMax,nI,nEl,iSeed,ExGen,P)
!   Type(ParticleList) PL(2)
!   integer iPL,nPL
!   integer icmin,icmax
!   real*8 icmean
!   real*8 dTau
!   real(dp) EShift
!   integer iSeed
!   Type(ExcitGen),pointer :: ExGen
!   Type(Particle) P
!
!   integer iParticle
!   integer iC 
!   integer iNewCount,iNewWeight
!   integer iSubPart
!   integer nEl,nI(nEl)
!   real*8 r,rnum
!   integer iPreExisted
!   
!   integer iGetExcitLevel
!   real*8 RAN2 ! external
!   character(25), parameter :: this_routine='PropagateParticles'
!
!   iParticle=0 !Just get the head node first and its info in P
!   call GetNextParticle(PL(iPL),iParticle,P)
!   ICMean=0.D0
!   ICMin=10
!   ICMax=0
!   do while (iParticle.gt.0)  !i.e. while we're at a valid particle
!! Each Particle has iWeight of subparticles at the same det
!!         write(6,*) iParticle,P%d%exGen%tagExData
!!         write(56,'(I)',advance='no')  P%d%iWeight
!
!      IC=iGetExcitLevel(Ni,P%Ni,NEl)
!      ICmean=ICmean+IC*(P%d%iWeight)
!      IF(IC.lt.ICMin) ICMin=IC
!      IF(IC.gt.ICMax) ICMax=IC
!      
!      iNewCount=0
!      do iSubPart=1,P%d%iWeight
!! With this current subparticle we do two processes.  Firstly we deplete the current particle according to 
!! Hii and tau
!         r=DREAL(P%d%hHii-EShift)*dTau !Prob of death
!         r=(1.D0)-r ! Prob of survival
!         iNewWeight=int(r)  ! Number of newly created definitely
!         r=r-iNewWeight
!!            if(r.lt.0) THEN write(6,*) "Warning: Self-destruction  probabililty <0"
!! This particle doesn't self-destruct if a random num is > r
!         rnum=ran2(iSeed)
!         if(r.gt.rnum) then
!!  Add it to the new list
!            iNewWeight=iNewWeight+1
!         endif
!         if(iNewWeight.gt.0) then
!!Copies accross to the nPL list (if actually creating more weight at particle) - sign is same sign as before
!!Will need to search though the nPL copy of PL to see if particle already exists. If it does then iPreExisted=1 
!            call AddParticle(PL(nPL),P%nI,iNewWeight,+P%d%iSgn,P%d%hHii,P%d%hHi0,P%d%IC0,P%d%exGen,iPreExisted)
!!Number of subparticles of this particle which have are present at end of creation/death process
!            iNewCount=iNewCount+iNewWeight
!         endif
!      enddo !iSubPart
!
!      call DelReference(P%d%exGen)   !Reduce the reference for this exgen by one. 
!!Even though this node is finished with in iPL, it may be present in nPL, so keep reference
!
!!Save the number of children we've spawned at this det
!      PL(iPL)%ParticleData(iParticle)%iProgeny=iNewCount
!      call GetNextParticle(PL(iPL),iParticle,P)
!   enddo !iParticle
!!      call DumpParticleList(6,PL(nPL))
!   ICMean=ICMean/(Pl(nPL)%NSubparticles)
!
!end Subroutine PropagateParticles
!
!!Arguments: Particle Lists(2),Current index, new index, Tau, NEl, Seed, ExGen for ..., Empty particle slot,HFDet,HF Energy
!subroutine SpawnParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P,nI,EHF)
!   Use Determinants, only: GetHElement3
!   Type(ParticleList) PL(2)
!   integer iPL,nPL
!   integer nEl
!   integer iSeed
!   real*8 dTau
!   Type(ExcitGen),pointer :: exGen
!   Type(Particle) P
!
!   integer iParticle
!   integer iSubPart
!   integer iAction
!   HElement_t hHij,hHi0
!   real(dp) HHjj,EHF
!   integer nJ(nEl),nI(nEl),IC0
!   integer iC,iCount,ExtraCreate
!   integer iNode,iGetExcitLevel
!   integer iSgn
!   real*8 pGen
!   real*8 r,rnum
!
!   real*8 ran2
!   character(25), parameter :: this_routine='SpawnParticles'
!
!!      write(56,*)
!! Secondly we might create positive or negative particle at an excitation
!   iParticle=0      ! setting iParticle to zero initialises run through nodes in list
!   call GetNextParticle(PL(iPL),iParticle,P)    !Finds next particle and puts it into P
!   do while (iParticle.gt.0)
!!         write(6,*) iParticle,P%d%exGen%tagExData
!!         call WriteDet(6,P%nI,nEl,.false.)
!!         write(6,*) P%d%iProgeny
!      
!      do iSubPart=1,P%d%iWeight     !Run through the integer weight on each particle
!! generate a random excitation of iParticle
!         call GenRandSymExcitIt3(P%nI,P%d%ExGen%exData,nJ,iSeed,iC,0,pGen,iCount)
!! The excitation is in nJ, and its gen prob in pGen.
!!            write(6,*)  pGen
!         hHij=GetHElement3(P%nI,nJ,iC)  !This is the connection strength to nJ from P
!
!! see if we create a particle in nJ, which we do with probability given by r
!         r=DREAL(abs(hHij))*dTau/pGen
!!            WRITE(6,*) r,hHij
!
!!If probability > 1, then we can just create change the weights of particle created at nJ by more than +-1
!         ExtraCreate=INT(r)
!         r=r-real(ExtraCreate)
!
!         rnum=ran2(iSeed)
!! There is new weight created if r > a random num
!!         if(r.gt.1) write(6,*) "Warning: Child creation probabililty>1"
!         if(r.gt.rnum) ExtraCreate=ExtraCreate+1
!         
!         if(ExtraCreate.gt.0) then
!! Particle is created, or if particle is already in the list, its weight is increased by ExtraCreate
!!  See if this particle's already in the list.  iNode is the node it's in in PL(nPL) if it exists, or a new clean node ready to receive it.  iAction is set 
!! to 1 if it doesn't already exist.
!            call FindParticleNode(PL(nPL),nJ,iNode,iAction)
!            if(iAction.ne.0) then
!
!!  Now we have to get info about this particle to add it to the particle list
!               hHjj=GetHElement3(nJ,nJ,0)   !Find Diagonal element
!               call AllocateExGen(exGen,nJ,nEl,this_routine)    !Find excitation generator for nJ
!
!               IC0=iGetExcitLevel(nI,nJ,nEl)
!               IF(IC0.eq.0) THEN
!!Particle created is on the HF det - hHi0 is HF Energy
!                  hHi0=(EHF)
!               ELSEIF(IC0.eq.2) THEN
!!Created particle is a double excit - calculate connection back to HF
!                  hHi0=GetHElement3(nI,nJ,IC0)
!               ELSE
!!Created particle is a single/ > double, therefore connection back to HF = 0
!                  hHi0=(0.D0)
!               ENDIF
!
!            else
!! Already exists, so we don't recalculate info
!            endif
!!  Add it to the new list
!!   If the Hij element is -ve then keep the same sign, otherwise flip the sign
!            iSgn=-P%d%iSgn
!            if(DREAL(hHij).lt.0) iSgn=-iSgn
!!  Adds on the weight if iAction=0 or creates a new node w.r.t iNode and iAction if iAction +-1
!!  PL to add the weight to, determinant to add (if iAction+-1), weight to add, sign of weight to add, Diag Elem (if iAction+-1),
!!  Exgen for nJ (if iAction+-1), Node in PL(nPL) to add weight, +- = New particle
!            call AddParticleInt(PL(nPL),nJ,ExtraCreate,iSgn,hHjj,hHi0,IC0,ExGen,iNode,iAction)
!         endif
!
!      enddo !iSubPart
!      call GetNextParticle(PL(iPL),iParticle,P)
!   enddo !iParticle
!end subroutine SpawnParticles
!
!
!! Particle management functions
!
!!Allocate the determinant within a particle
!subroutine AllocParticle(P,nEl)
!   implicit none
!   integer nEl
!   Type(Particle) P
!   allocate(P%nI(nEl))
!end subroutine
!
!!DeAllocate the determinant within a particle
!subroutine DeallocParticle(P)
!   implicit none
!   integer nEl
!   Type(Particle) P
!   deallocate(P%nI)
!end subroutine
!
!!Particle List management functions
!
!!Allocate arrays within a particle list, and clear then
!subroutine AllocParticleList(PL,nMaxParticles, nEl)
!   Use global_utilities
!   implicit none
!   Type(ParticleList) PL
!   integer nMaxParticles
!   integer nEl
!   integer ierr
!   character(25), parameter :: this_routine='AllocParticleList'
!   PL%nEl=nEl
!   PL%nMaxParticles=nMaxParticles
!   allocate(PL%ParticleDets(nEl,nMaxParticles), stat=ierr)
!   LogAlloc(ierr,'ParticleDets',nMaxParticles*nEl,4,PL%tagParticleDets)
!   allocate(PL%ParticleData(nMaxParticles), stat=ierr)
!   LogAlloc(ierr,'ParticleData',nMaxParticles,ParticleDataSizeB,PL%tagParticleData)
!   call ClearParticleList(PL)
!end subroutine AllocParticleList
!
!subroutine DeallocParticleList(PL)
!   Use global_utilities
!   implicit none
!   Type(ParticleList) PL
!   integer ierr
!   character(25), parameter :: this_routine='DeallocParticleList'
!   LogDealloc(PL%tagParticleDets)
!   Deallocate(PL%ParticleDets)
!   LogDealloc(PL%tagParticleData)
!   Deallocate(PL%ParticleData)
!end subroutine DeallocParticleList
!
!subroutine AddParticle(PL,nI,iWeight,iSgn,hHii,hHi0,IC0,exGen,iPreExisted)
!!Add a particle to this list
!! First find a node to put it at
!   implicit none
!   Type(ParticleList) PL
!   integer nI(PL%nEl)
!   integer iWeight,iSgn,IC0
!   real(dp) hHii
!   HElement_t hHi0
!   Type(ExcitGen), pointer :: ExGen
!   integer iNode,iAction
!   integer iPreExisted
!   iPreExisted=0
!   call FindParticleNode(PL,nI,iNode,iAction)
!   if(iAction.eq.0) iPreExisted=1
!   call AddParticleInt(PL,nI,iWeight,iSgn,hHii,hHi0,IC0,ExGen,iNode,iAction)
!end subroutine AddParticle
!
!
!!  Adds on the weight if iAction=0 or creates a new node w.r.t iNode and iAction if iAction +-1
!!  PL to add the weight to, determinant to add (if iAction+-1), weight to add, sign of weight to add, Diag Elem (if iAction+-1),
!!  Connection to HF(if iAction+-1), ExcitLevel from HF, Exgen for nJ (if iAction+-1), Node in PL(nPL) to add weight, +- = New particle
!subroutine AddParticleInt(PL,nI,iWeight,iSgn,hHii,hHi0,IC0,ExGen,iNode,iAction)
!!Add a particle to this list
!! First find a node to put it at
!   implicit none
!   Type(ParticleList) PL
!   integer nI(PL%nEl)
!   integer iWeight,iSgn,IC0
!   real(dp) hHii
!   HElement_t hHi0
!   Type(ExcitGen), pointer :: ExGen
!   integer iNode,iAction
!   if(iAction.eq.0) then
!! This particle matches this node.  Just change its Weight
!      PL%nSubParticles=PL%nSubParticles-PL%ParticleData(iNode)%iWeight  !Remove original weight on node from total weight
!      PL%ParticleData(iNode)%iWeight=PL%ParticleData(iNode)%iWeight*PL%ParticleData(iNode)%iSgn+iSgn*iWeight  !Calculate new weight on the particle
!      if(PL%ParticleData(iNode)%iWeight.lt.0) then
!!Find the sign of the weight
!         PL%ParticleData(iNode)%iWeight=-PL%ParticleData(iNode)%iWeight     !Weight is always positive value
!         PL%ParticleData(iNode)%iSgn=-1     !Update the sign of the particle
!      else
!         PL%ParticleData(iNode)%iSgn=1
!      endif
!      PL%nSubParticles=PL%nSubParticles+PL%ParticleData(iNode)%iWeight    !Update the new total weight of the system
!      return
!
!   elseif(iAction.eq.1) then
!! This particle is greater than this node, and we're at the end of the line
!      if(iNode.ne.0) PL%ParticleData(iNode)%iAfter=PL%iNextFreeParticle
!      if(iNode.eq.0) PL%iHeadParticle=PL%iNextFreeParticle
!      PL%ParticleData(PL%iNextFreeParticle)%iParent=iNode
!      iNode=PL%iNextFreeParticle
!      PL%iNextFreeParticle=iNode+1
!      PL%nParticles=PL%nParticles+1
!      if(PL%iNextFreeParticle.gt.PL%nMaxParticles) call Stop_All("Particle List Full")
!   elseif(iAction.eq.-1) then
!! This particle is less than this node, and we're at the end of the line
!      if(iNode.ne.0) PL%ParticleData(iNode)%iBefore=PL%iNextFreeParticle
!      if(iNode.eq.0) PL%iHeadParticle=PL%iNextFreeParticle
!      PL%ParticleData(PL%iNextFreeParticle)%iParent=iNode
!      iNode=PL%iNextFreeParticle
!      PL%iNextFreeParticle=iNode+1
!      PL%nParticles=PL%nParticles+1
!      if(PL%iNextFreeParticle.gt.PL%nMaxParticles) call Stop_All("Particle List Full")
!   endif   
!   PL%ParticleData(iNode)%iWeight=iWeight
!   PL%ParticleData(iNode)%iSgn=iSgn
!   PL%ParticleData(iNode)%hHii=hHii
!   PL%ParticleData(iNode)%hHi0=hHi0
!   PL%ParticleData(iNode)%IC0=IC0
!   PL%ParticleData(iNode)%ExGen=>ExGen
!   call AddReference(ExGen)
!   PL%ParticleData(iNode)%iAfter=0
!   PL%ParticleData(iNode)%iBefore=0
!   PL%ParticleDets(:,iNode)=nI(:)
!   PL%nSubParticles=PL%nSubParticles+iWeight
!end subroutine AddParticleInt
!
!!Zero a particle List
!subroutine ClearParticleList(PL)
!   implicit none
!   Type(ParticleList) PL
!   PL%iNextFreeParticle=1
!   PL%nSubParticles=0
!   PL%nParticles=0
!   PL%iHeadParticle=0
!   PL%ParticleDets(:,:)=0
!end subroutine ClearParticleList
!
!!Get the next particle with non-zero weight.  We traverse the tree as follows.  
!!  Start at HEAD node
!!    The next node after CURRENT is
!!     iBefore (if its weight is non-zero)
!!     iAfter (if its weight is non-zero)
!!     Recursing up the tree, the first non-zero iAfter, otherwise 0
!
!!The Int version does not do a copy to P
!!  If there are no more particles, it returns iParticle=0
!subroutine GetNextParticleInt(PL,iParticle)
!   implicit none
!   Type(ParticleList) PL
!   integer iParticle
!   integer iOldPart
!   if(iParticle.eq.0) then
!! Just get the head if there is one
!      if(PL%nParticles.gt.0) then
!         iParticle=PL%iHeadParticle
!      endif
!   else
!      if(PL%ParticleData(iParticle)%iBefore.ne.0) then
!         iParticle=PL%ParticleData(iParticle)%iBefore
!      else
!! We need to look upwards in the tree
!         iOldPart=0
!         do while(iParticle.gt.0.and.PL%ParticleData(iParticle)%iAfter.eq.iOldPart)
!!            WRITE(6,*) iOldPart, iParticle
!            iOldPart=iParticle
!            iParticle=PL%ParticleData(iParticle)%iParent
!            if(PL%ParticleData(iParticle)%iBefore.eq.iOldPart) iOldPart=0
!!            WRITE(6,*) iOldPart, iParticle
!         enddo
!         if(iParticle.gt.0) iParticle=PL%ParticleData(iParticle)%iAfter
!      endif
!   endif
!end subroutine GetNextParticleInt
!
!
!! GetnextParticle Gets the particle after iParticle.  If iParticle=0, it returns the head node.
!!  It puts the Particle info into P.  If there are no more particles, it returns iParticle=0
!subroutine GetNextParticle(PL,iParticle,P)
!   implicit none
!   Type(ParticleList) PL
!   Type(Particle) P
!   integer iParticle
!   call GetNextParticleInt(PL,iParticle)
!   call GetParticle(PL,iParticle,P)
!end subroutine GetNextParticle
!
!!Extract the data from a particle
!subroutine GetParticle(PL,iParticle,P)
!   implicit none
!   Type(ParticleList) PL
!   Type(Particle) P
!   integer iParticle
!   if(iParticle.ne.0) then
!      P%d=PL%ParticleData(iParticle)
!      P%nI(:)=PL%ParticleDets(:,iParticle)
!   endif
!end subroutine GetParticle
!
!!run GetNextParticle until we find one with zero weight
!subroutine GetNextZeroParticle(PL,iParticle,P)
!   implicit none
!   Type(ParticleList) PL
!   Type(Particle) P
!   integer iParticle
!   call GetNextParticleInt(PL,iParticle)
!   do while(iParticle.gt.0.and.PL%ParticleData(iParticle)%iWeight.ne.0)
!      call GetNextParticleInt(PL,iParticle)
!   enddo
!   call GetParticle(PL,iParticle,P)
!end Subroutine GetNextZeroParticle
!
!! Remove a particle from the node list, and fix up the list
!subroutine DeleteParticle(PL,iParticle)
!   implicit none
!   Type(ParticleList) PL
!   integer iParticle
!   integer iParent
!   integer iHanging
!
!!  We remove this node first
!!   If the particle has a parent
!   iParent=PL%ParticleData(iParticle)%iParent
!   if(iParent.ne.0) then
!      if(PL%ParticleData(iParent)%iBefore.eq.iParticle) then
!!  We're a before, so we link in our own before, and leave the iAfter hanging
!         PL%ParticleData(iParent)%iBefore=PL%ParticleData(iParticle)%iBefore
!         PL%ParticleData(PL%ParticleData(iParticle)%iBefore)%iParent=iParent
!         iHanging=PL%ParticleData(iParticle)%iAfter
!      else
!!  We're an after, so we link in our own after, and leave the iBefore hanging
!         PL%ParticleData(iParent)%iAfter=PL%ParticleData(iParticle)%iAfter
!         PL%ParticleData(PL%ParticleData(iParticle)%iAfter)%iParent=iParent
!         iHanging=PL%ParticleData(iParticle)%iBefore
!      endif
!   else
!!  We're the parent, so we make our before the headif we can, and leave the after hanging
!      PL%iHeadParticle=PL%ParticleData(iParticle)%iBefore
!      iHanging=PL%ParticleData(iParticle)%iAfter
!      if(PL%iHeadParticle.eq.0) then
!!  We don't have a before, so try our after
!         PL%iHeadParticle=PL%ParticleData(iParticle)%iAfter
!         PL%ParticleData(iParticle)%iParent=iParent
!         iHanging=0
!      else
!         PL%ParticleData(PL%iHeadParticle)%iParent=0
!      endif
!      if(PL%iHeadParticle.eq.0) then
!! no before or after, so we must be top of the list
!         PL%iHeadParticle=0
!      else
!         PL%ParticleData(PL%iHeadParticle)%iParent=0
!      endif
!   endif
!!  Deal with the hanging
!   if(iHanging.ne.0) then
!      call AddParticleNode(PL,iHanging)
!   endif
!!  Now note that we are zero
!   PL%ParticleDets(1,iParticle)=0
!! consolidate the particle list a little
!   do while(PL%ParticleDets(1,PL%iNextFreeParticle).eq.0)
!      PL%iNextFreeParticle=PL%iNextFreeParticle-1
!   enddo
!end subroutine DeleteParticle
!
!subroutine AddParticleNode(PL,iParticle)
!   implicit none
!   Type(ParticleList) PL
!   integer iParticle
!   integer iNode,iAction
!   integer iPreExisted
!   call FindParticleNode(PL,PL%ParticleDets(:,iParticle),iNode,iAction)
!   if(iAction.eq.0) stop 'Found duplicate particle in list'
!   if(iAction.eq.1) then
!! This particle is greater than this node, and we're at the end of the line
!      if(iNode.ne.0) PL%ParticleData(iNode)%iAfter=iParticle
!      if(iNode.eq.0) PL%iHeadParticle=iParticle
!      PL%ParticleData(iParticle)%iParent=iNode
!      iNode=iParticle
!   elseif(iAction.eq.-1) then
!! This particle is less than this node, and we're at the end of the line
!      if(iNode.ne.0) PL%ParticleData(iNode)%iBefore=iParticle
!      if(iNode.eq.0) PL%iHeadParticle=iParticle
!      PL%ParticleData(iParticle)%iParent=iNode
!      iNode=iParticle
!   endif   
!end subroutine AddParticleNode
!
!!returns iNode as the node nI is in if already there (and Action=0).  Otherwise Action=-1 for needing a new node before the returned iNode, or +1 for after.
!subroutine FindParticleNode(PL,nI,iNode,iAction)
!   implicit none
!   Type(ParticleList) PL
!   integer nI(*)
!   integer iNode
!   integer iAction,iCmp
!   integer iCmpDets  !external function
!
!   iNode=PL%iHeadParticle
!   if(iNode.eq.0) then
!      iAction=-1
!      return 
!   endif
!   iCmp=iCmpDets(nI,PL%ParticleDets(:,iNode),PL%nEl)
!   iAction=0
!   do while(iCmp.ne.0)
!      if(iCmp.lt.0) then
!!our det is less than the current, so look at before
!         if(PL%ParticleData(iNode)%iBefore.eq.0) then
!            iAction=-1
!            iCmp=0
!         else
!            iNode= PL%ParticleData(iNode)%iBefore
!         endif
!      else
!!our det is more than the current, so look at after
!         if(PL%ParticleData(iNode)%iAfter.eq.0) then
!            iAction=1
!            iCmp=0
!         else
!            iNode= PL%ParticleData(iNode)%iAfter
!         endif
!      endif
!      if(iCmp.ne.0)    iCmp=iCmpDets(nI,PL%ParticleDets(:,iNode),PL%nEl)
!   enddo
!   return
!end subroutine FindParticleNode
!
!subroutine DumpParticleList(iunit,PL)
!   implicit none
!   Type(ParticleList) PL
!   integer iUnit
!   Type(Particle) P
!   integer iParticle
!   iParticle=0
!   Allocate(P%nI(PL%nEl))
!   call GetNextParticle(PL,iParticle,P)
!   do while (iParticle.ne.0)
!      write(iUnit,"(I12)",advance='no') iParticle
!      call WriteDet(iUnit,P%nI,PL%nEl,.false.)
!      write(iUnit,"(7I12)") P%d%iWeight,P%d%iSgn,P%d%iParent, P%d%iBefore,P%d%iAfter,P%d%ExGen%tagExData,P%d%ExGen%iRefCount
!      call GetNextParticle(PL,iParticle,P)
!   enddo
!   Deallocate(P%nI)
!end subroutine DumpParticleList
!   
!subroutine AllocateExGen(exGen,nI,nEl,this_routine)
!   Use global_utilities
!   Type(ExcitGen), pointer :: exGen
!   integer nEl 
!   integer nI(nEl),nJ(nEl)
!   integer Store(6)  !Used for Excit gen
!   integer iLen,iC     
!   integer ierr
!   character(*) this_routine
!
!
!   Allocate(exGen)
!   exGen%iRefCount=0
!
!!.. Setup the spin excit generator
!   Store(1)=0
!   ! .true. for
!   ! 3 for both singles and doubles. (Brill theorem is taken into account if needed however)
!   call GenSymExcitIt3(nI,.true.,iLen,nJ,iC,0,Store,3)
!   ExGen%tagExData=0
!   Allocate(exGen%exData(iLen),STAT=ierr)
!   if(tLogExGens) LogAlloc(ierr,'exGen',iLen,8,exGen%tagExData)
!   ExGen%exData(1)=0
!!Now setup the generator (Store contains some necessary info)
!   call GenSymExcitIt3(nI,.true.,exGen%exData,nJ,iC,0,Store,3)
!   return
!end subroutine AllocateExGen
!   
!subroutine DeallocateExGen(exGen,this_routine)
!   Use global_utilities
!   Type(ExcitGen) exGen
!   integer ierr
!   character(*) this_routine
!   if(exGen%iRefCount.ne.0) then
!      stop 'Attempting to deallocate referenced excitation generator'
!   endif
!!   write(6,*) "Dealloc",ExGen%tagExData
!   if(tLogExGens) LogDealloc(ExGen%tagExData)
!   Deallocate(ExGen%ExData)
!end subroutine DeallocateExGen
!
!subroutine AddReference(exGen)
!   Type(ExcitGen) exGen
!   exGen%iRefCount=exGen%iRefCount+1
!end subroutine
!
!subroutine DelReference(exGen)
!   Type(ExcitGen) exGen
!   exGen%iRefCount=exGen%iRefCount-1
!end subroutine
!   
!
!subroutine CleanupParticles(PL,iPL,nPL,dTau,nEl,iSeed,ExGen,P)
!   Use Determinants, only: GetHElement3
!   Type(ParticleList) PL(2)
!   integer iPL,nPL
!   integer nEl
!   integer iSeed
!   real*8 dTau
!   Type(ExcitGen),pointer :: exGen
!   Type(Particle) P
!
!   integer iParticle
!   integer iSubPart
!   character(25), parameter :: this_routine='CleanupParticles'
!   iParticle=0
!   call GetNextParticle(PL(iPL),iParticle,P)
!   do while (iParticle.gt.0)
!      if(P%d%exGen%iRefCount.eq.0) then
!!            WRITE(6,*) "Progeny 0",iParticle
!         call DeallocateExGen(P%d%ExGen,this_routine)
!      endif
!! Just see if we need to keep our old excitation generator around, otherwise kill it
!      call GetNextParticle(PL(iPL),iParticle,P)
!   enddo !iParticle
!end subroutine CleanupParticles

end Module MCDets
