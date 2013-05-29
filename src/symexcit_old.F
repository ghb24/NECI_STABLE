! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
!! NB SYMSETUPEXCITS uses a store which has integer(int64) to cope with 64-bit machines, and passes this directly to SYMSETUPEXCITS2
!!  Iterators which use GENSYMEXCITIT2 use an INTEGER STORE, and create an integer(int64) STORE2 internally, and pass this to SYMSETUPEXCITS2

!! SYMSETUPEXCITS takes an integer(int64) STORE(6) which stores the pointers used by SYMSETUPEXCITS2
!!   This STORE can then be passed to SYMGENEXCITS.
C.. Setup for symmetry routine below.
      SUBROUTINE SYMSETUPEXCITS(NI,NEL,NBASIS,STORE,TCOUNT,
     &      ICOUNT,ILEVEL,iMinElec1,iMaxElec1)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB,g1
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         use legacy_data, only: irat
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         integer, POINTER :: DSTORE(:)
         INTEGER STORE(6)
         INTEGER ICOUNT
         LOGICAL TCOUNT
         INTEGER ILEVEL
         INTEGER iMinElec1,iMaxElec1

         STORE(1:6)=0
         allocate(DSTORE(SymClassSize*NEL+(nBasis/32)+1
     &      +SymmetrySize*(NEL*NEL+1)))
         STORE(1)=1 !Meaning the first element of DSTORE
         ! STORE(1) --> 
         CALL SYMSETUPEXCITS3(NI,NEL,G1,NBASIS,STORE,STORE(1),
     &      STORE(1),STORE(1),STORE(1),TCOUNT,
     &      ICOUNT,DSTORE(1),
     &      DSTORE(SymClassSize*NEL+1),
     &      DSTORE(SymClassSize*NEL+1+(nBasis/32)+1),ILEVEL,iMinElec1,
     &      iMaxElec1)
C.. If we're just counting, we don't need to keep DSTORE
         IF(TCOUNT) deallocate(dstore)
      END

!We wish to calculate what excitation class the excitation NI->NJ falls into, with the appropriate 
! IFROM class and index within that, IFROMINDEX, and ITO class, and index within that. ITOINDEX
! After that, we generate the probability that nJ would be an excitation from nI.
! We require the excitation generator of I, which is used to pass info to GenExcitProbInternal (in symexcit2.F90)
      SUBROUTINE GenExcitProb(nI,nJ,nEl,nIExcitor,G1,nBasisMax,Arr,
     &      nBasis,pGen)
         use constants, only: dp
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         IMPLICIT NONE
         INTEGER nEl,nI(nEl),nJ(nEl),nBasis,nBasisMax(*)
         INTEGER, target :: nIExcitor(*)
         TYPE(BasisFN) G1(nBasis) 
         INTEGER, pointer :: DStore(:)
         real(dp) pGen
         real(dp) Arr(nBasis,2)
         DSTORE=>nIExcitor(nIExcitor(1):niExcitor(2)-1)
         CALL GenExcitProb2(nI,nJ,nEl,G1,nBasisMax,Arr,nBasis,
     &      nIExcitor(nIExcitor(4)),nIExcitor(nIExcitor(5)),
     &      DStore(SymClassSize*nEl+1),
     &      DStore(SymClassSize*nEl+1+(nBasis/32)+1),
     &      nIExcitor(nIExcitor(2)),nIExcitor(23),pGen)
      END

!  Must be called after SYMSETUPEXCITS.  Needs integer(int64) STORE and frees memory therein
      SUBROUTINE SYMGENEXCITS(NI,NEL,STORE,LSTE,DSTORE,ExcitTypes,
     &   nAllowPPS,OrbPairs,SymProdInd,
     &   ICLIST,ICOUNT)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use SymData, only: SymClassSize
         use global_utilities
         IMPLICIT NONE
         INTEGER NEL,NI(NEL)
         INTEGER STORE(6)
         INTEGER  ORBPAIRS(2,*)
         INTEGER SYMPRODIND(2,3,1:*)
         INTEGER NALLOWPPS(3,*)

         INTEGER EXCITTYPES(5,*)
         INTEGER DSTORE(*)
         type(timer), save :: proc_timer
         INTEGER LSTE(NEL,*),ICLIST(*),ICOUNT
         proc_timer%timer_name='SYMGNEXCIT'
         call set_timer(proc_timer,65)

         CALL SYMGENALLEXCITS(NI,NEL,EXCITTYPES,STORE(6),DSTORE(1),
     &        SymProdInd,DSTORE(SymClassSize*NEL+1),ORBPAIRS,
     &            LSTE,ICLIST,ICOUNT)
   
!         CALL N_FREEM(IP_EXCITTYPES)
!         CALL N_FREEM(IP_NALLOWPPS)
!         CALL N_FREEM(IP_ORBPAIRS)
!         CALL N_FREEM(IP_SYMPRODIND)
!         CALL N_FREEM(IP_DSTORE)

         call halt_timer(proc_timer)
 
      END

C.. IF(TSETUP) Generate an iterator which allows up to double excitations to be generated
C.. one at a time (in an unordered fashion) from a given det.  THis needs to be called twice,
C.. first with STORE(1)=0, and then again.  Finally it can be called with TSETUP=.FALSE. 
C.. to actually generate the excitations
!! First run.  TSETUP=.TRUE.  STORE(1)=0.  NMEM can just be an INTEGER
!!  requires INTEGER STORE(6) to hold initialization data.  NMEM is returned containing the memory required to hold the excitation generator (call this nLength)
!! Second run,  TSETUP=.TRUE.  STORE(:) contains data from first run.  NMEM now is an array of nLenght INTEGERs.  The excitation generator data is stored in this.
!! Third run,   TSETUP=.FALSE.  will generate an excitation of NI, and put it in NEL.  

      SUBROUTINE GENSYMEXCITIT2(NI,NEL,G1,NBASIS,
     &         TSETUP,NMEM,NJ,IC,STORE,ILEVEL)
         use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
         use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB
         use legacy_data, only: irat
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NBASIS
         INTEGER G1(nBasis)
!  STORE contains lengths of various components of the excitation generator
         INTEGER STORE(6)
!  STORE will contain the addesses of various components of the excitation generator, and is passed to SYMSETUPEXCITS2
         INTEGER, target :: NMEM(*)
         INTEGER NJ(NEL),IC
         LOGICAL TSETUP
         INTEGER ILEVEL
         CALL GenSymExcitIt2Par(NI,NEL,G1,NBASIS,
     &         TSETUP,NMEM,NJ,IC,STORE,ILEVEL,1,NEL)
      END
