subroutine VaspSystemInit(ArrLEN)
   use global_utilities
   use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
   use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,NullBasisFn
   use vasp_interface
   use SymData, only: nRot,PropBitLen,tAbelian,nProp,KPntSym,tagKPntSym
   use constants, only: dp,sizeof_int
   implicit none
   integer :: ArrLEN
   integer :: i,ik
   character(*),parameter :: this_routine='VaspSystemInit'

   ArrLEN=nStates*2
   nRot=nKP
   PropBitLen=15
   nProp(:)=KPMsh(:)
   tAbelian=.true.
   write (6,*) tAbelian,nProp,PropBitLen,nRot,ArrLEN

   allocate(KPntSym(3,nKP))
   call LogMemAlloc('KPntSym',3*nKP,4,this_routine,tagKPntSym)
   do ik=1,nKP
      do i=1,3
         KPntSym(i,ik)=int(kpnts(i,ik)*real(nProp(i),dp),sizeof_int)
      end do
   end do

   write (6,*) tAbelian,nProp,PropBitLen,nRot,ArrLEN
   return
end subroutine VaspSystemInit

subroutine VASPInitIntegrals(nOrbUsed,ECore,tOrder)
   use constants, only: dp
   use SystemData, only: BasisFN,nEl
   use OneEInts, only: TMatSym, TMatInd
   use vasp_interface
   use UMatCache, only: SetupUMatCache,UMat2D
   use global_utilities
   use constants, only: dp
   implicit none
   integer :: nOrbUsed
   real(dp) ::  ECore
   logical :: tOrder
   type(timer), save :: proc_timer
   integer :: I,J,II,A,B,nStatesUsed
   HElement_t(dp) :: HarXC,HarXCSum
   
   proc_timer%timer_name='VASPInitInts'
   call set_timer(proc_timer)
   ! ECore=EIonIon???
   ECore=0.0_dp
   write (6,*) 'Core Energy: ',ECORE
   nStatesUsed=nOrbUsed/2

   call SetupUMatCache(nStatesUsed,NSTATESUSED.NE.NSTATES)

#ifdef __CMPLX
   HarXCSum=cmplx(0.0_dp,0.0_dp,dp)
#else
   HarXCSum=0.0_dp
#endif
   write (6,*) "Calculating TMAT"
   open(10,file='TMAT',status='unknown')
   do I=1,nStatesUsed
      ! Subtract out the double counting. Assume closed-shell.
#ifdef __CMPLX
      HarXC=cmplx(0.0_dp,0.0_dp,dp)
#else
      HarXC=0.0_dp
#endif
      do J=1,nEl/2
         if (I.ne.J) then
            A=min(I,J)
            B=max(I,J)
            HarXC=HarXC-(2)*UMat2D(A,B)+UMat2D(B,A)
         end if
         if (i.le.j) HarXCSum=HarXCSum+HarXC
      end do
      II=I*2-1
      TMATSYM(TMatInd(II+1,II+1))=(eigv(I))-HarXC
      write (10,*) I,J,TMATSYM(TMatInd(II+1,II+1))
   end do
   

   if (tOrder) then
      call Stop_All('VASPInitIntegrals','tOrder not implemented in VASP interface yet.')
   end if
   write (6,*) "Finished TMAT"
   close (10)

   call halt_timer(proc_timer)
   
   return
end subroutine VASPInitIntegrals

subroutine VASPBasisInit(ARR,BRR,G1,LEN)
   ! Largely lifted from the CPMD analogue.  Should be doing (roughly) the same
   ! thing to get going!
   use SystemData, only: Symmetry,SymmetrySize,SymmetrySizeB
   use SystemData, only: BasisFN,BasisFNSize,BasisFNSizeB,nBASISMax,NullBasisFn
   use vasp_interface, only: nStates,nKP,KPntInd,eigv
   use SymData, only: KPntSym,nSym
   use constants, only: dp
   use sym_mod
   implicit none
   integer :: LEN
   real(dp) :: ARR(LEN,2)
   integer :: BRR(LEN)
   type(BasisFN) :: G1(LEN)
   integer :: i
   type(Symmetry) :: iDecomp

   NBASISMAX(1:3,1:2)=0
   NBASISMAX(1,3)=2
   NBASISMAX(4,1)=-1
   NBASISMAX(4,2)=1

!  set ISPINSKIP=0, to tell the SCRs that there's no UMAT
   NBASISMAX(2,3)=0

   call GenKPtIrreps(NKP,NKP,KPNTIND,NSTATES)

   G1(1:LEN)=NullBasisFn

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
