subroutine VaspSystemInit(ArrLEN)
   use System, only: Symmetry,SymmetrySize,SymmetrySizeB
   use System, only: BasisFN,BasisFNSize,BasisFNSizeB
   use vasp_interface
   implicit none
   include 'sym.inc'
   integer :: ArrLEN
   integer :: i,ik

   ArrLEN=nStates*2
   nRot=nKP
   PropBitLen=15
   nProp(:)=KPMsh(:)
   tAbelian=.true.
   write (6,*) tAbelian,nProp,PropBitLen,nRot,ArrLEN

   call Memory(IP_KPntSym,3*nKP,'KPntSym')
   do ik=1,nKP
      do i=1,3
         KPntSym(i,ik)=kpnts(i,ik)*nProp(i)
      end do
   end do

   write (6,*) tAbelian,nProp,PropBitLen,nRot,ArrLEN
   return
end subroutine VaspSystemInit

subroutine VASPInitIntegrals(nOrbUsed,ECore,tOrder)
   use HElem
   use System, only: BasisFN,nEl
   use OneEInts, only: TMatSym, TMatInd
   use vasp_interface
   use UMatCache, only: SetupUMatCache,UMat2D
   use MemoryManager, only: LogMemAlloc
   implicit none
   integer :: nOrbUsed
   real(q) ::  ECore
   logical :: tOrder
   integer :: iSub,I,J,II,A,B,nStatesUsed,ierr
   type(HElement) :: HarXC,HarXCSum
   character(*), parameter :: thisroutine='VASPInitIntegrals'
   
   call TISET('VASPInitInts',ISUB)
   ! ECore=EIonIon???
   ECore=0.d0
   write (6,*) 'Core Energy: ',ECORE
   nStatesUsed=nOrbUsed/2

   call SetupUMatCache(nStatesUsed,NSTATESUSED.NE.NSTATES)

   HarXCSum=dcmplx(0.d0,0.d0)
   write (6,*) "Calculating TMAT"
   open(10,file='TMAT',status='unknown')
   do I=1,nStatesUsed
      ! Subtract out the double counting. Assume closed-shell.
      HarXC=dcmplx(0.d0,0.d0)
      do J=1,nEl/2
         if (I.ne.J) then
            A=min(I,J)
            B=max(I,J)
            HarXC=HarXC-HElement(2)*UMat2D(A,B)+UMat2D(B,A)
         end if
         if (i.le.j) HarXCSum=HarXCSum+HarXC
      end do
      II=I*2-1
      TMATSYM(TMatInd(II+1,II+1))=HElement(eigv(I))-HarXC
      write (10,*) I,J,TMATSYM(TMatInd(II+1,II+1))
   end do
   

   if (tOrder) then
      call Stop_All('VASPInitIntegrals','tOrder not implemented in VASP interface yet.')
   end if
   write (6,*) "Finished TMAT"
   close (10)

   call TIHALT('VASPInitInts',ISUB)
   
   return
end subroutine VASPInitIntegrals

subroutine VASPBasisInit(ARR,BRR,G1,LEN)
   ! Largely lifted from the CPMD analogue.  Should be doing (roughly) the same
   ! thing to get going!
   use System, only: Symmetry,SymmetrySize,SymmetrySizeB
   use System, only: BasisFN,BasisFNSize,BasisFNSizeB,nBASISMax
   use vasp_interface, only: q,nStates,nKP,KPntInd,eigv
   implicit none
   include 'sym.inc'
   real(q) :: ARR(LEN,2)
   integer :: BRR(LEN),LEN
   type(BasisFN) :: G1(LEN)
   integer :: i
   type(Symmetry) :: iDecomp
   integer(8) :: ComposeAbelianSym

   NBASISMAX(1:3,1:2)=0
   NBASISMAX(1,3)=2
   NBASISMAX(4,1)=-1
   NBASISMAX(4,2)=1

!  set ISPINSKIP=0, to tell the SCRs that there's no UMAT
   NBASISMAX(2,3)=0

   call GenKPtIrreps(NKP,NKP,KPNTIND,NSTATES)

   call IAZZERO(G1,LEN*BasisFNSize)

   do i=1,nStates
      IDECOMP%s=ComposeAbelianSym(KpntSym(:,KPntInd(I)))
      G1(I*2-1)%SYM=IDECOMP
      G1(I*2)%SYM=IDECOMP
      G1(I*2-1)%MS=-1
      G1(I*2)%MS=1
      G1(I*2-1)%K(1:3)=0
      G1(I*2)%K(1:3)=0
      ARR(2*I-1,1)=eigv(I)
      ARR(2*I,1)=eigv(I)
      ARR(2*I-1,2)=eigv(I)
      ARR(2*I,2)=eigv(I)
      BRR(2*I-1)=2*I-1
      BRR(2*I)=2*I
   end do

   WRITE (6,*) 'Using Abelian symmetry formulation.'
   NBASISMAX(5,1)=0
   NBASISMAX(5,2)=NSYM-1

! Show it's a generic spatial basis
   NBASISMAX(3,3)=1

   return
end subroutine VASPBasisInit
