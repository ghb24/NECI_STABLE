#ifdef PARALLEL
!Module with Parallel information etc.
!
!  NECI will run as a standalone parallel app or in conjunction with CPMD.
!    Only standalone is currently implemented.

!  Parallelization is done over occupied electrons.  
!    These are split among each processor such that there are an approximately equal number of pairs of electrons on each processor.
!    Each processor is given a set of 'Electron 1'.  From this it can generate a set of single excitations as well
!      as a set of double excitations.  Double excitations can have any electron past Electron 1 as Electron 2.
!      This means that the lower electron numbers will have more possible pairs.
!    Parallelization is supported by the symmetry excitation generators, interfaced through GenSymExcitIt3Par
!    There is no clean way to automatically parallelize high-vertex graph code, so each parallel routine must
!      be specifically written.
!
!  ROUTINES
!     MPIInit     Setup MPI and init Nodefiles if we're standalone
!     MPIEnd      Shutdown MPI
!     MPIStopAll  Abort all processors
!     MPIISum     Sum an array of integers among all processors, and distribute the results array to each.
!     MPIDSum     As MPIISum but for real*8
!     MPIHelSum   As MPIDSum but for Type(HElement)


module Parallel
! This is mixed case to avoid it being picked up by the the Configure script, as it doesn't require a module file.
   uSE mpi
   IMPLICIT NONE
   save
   integer iProcIndex
   integer nProcessors
Contains

! MPIInit
!     tExternal    (in) Set if using VASP/CPMD's MPI interface, so we don't have to init our own.
!
!Determine the number of processors, and fork each off to its own NodeFile output file

Subroutine MPIInit(tExternal)
   Use Determinants, only: FDet
   implicit none
   logical tExternal
   integer numtasks, rank, ierr, rc
   integer a,b,g
   character*20 NodeFile
   if(tExternal) then
     write(6,*) 'Using CPMD MPI configuration'
   else 
     write(6,*) 'Initing MPI'

     call MPI_INIT(ierr)
     if (ierr .ne. MPI_SUCCESS) then
        print *,'Error starting MPI program. Terminating.'
        call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
     end if
   endif
   call MPI_COMM_RANK(MPI_COMM_WORLD, iProcIndex, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcessors, ierr)
   if(tExternal) then
      write(6,*) "NECI Processor ",iProcIndex+1,'/',nProcessors
   else
      if(iProcIndex.eq.0) then
         write(6,*) "Processor ",iProcIndex+1,'/',nProcessors, ' as head node.'
      else

         write(6,*) "Processor ",iProcIndex+1,'/',nProcessors, ' moving to local output.'
         if (iProcIndex.le.9) then
             write (NodeFile,'(a,i1)') 'NodeFile',iProcIndex+1
         else
             write (NodeFile,'(a,i2)') 'NodeFile',iProcIndex+1
         end if
         write(6,*) "outfile=",NodeFile
         close(6,status="keep")
         open(6,file=NodeFile)
         write(6,*) "Processor ",iProcIndex+1,'/',nProcessors, ' on local output.'
      endif
!  Just synchronize everything briefly
      a=iProcIndex+1
      call MPIISum(a,1,g)
      WRITE(6,*) "Sum: ",g
   endif
   
end Subroutine

! MPIEnd
!
!  tExternal    (in)  Set if using an external program's MPI interface
!  (currently CPMD or VASP).
!
!  Shutdown our MPI Interface if we're not using CPMD/VASP's
Subroutine MPIEnd(tExternal)
   uSe mpi
   implicit none
   logical tExternal
   integer ierr
   if(.not.tExternal) then
      call MPI_FINALIZE(ierr)
   endif
end subroutine

! MPIStopAll
!  
!  n     (in)  MPI Abort parameter (???)
!  
!  Abort all processors
Subroutine MPIStopAll(n)
   integer n,rc
   call MPI_ABORT(MPI_COMM_WORLD, rc, n)
   stop 999
end

! MPIISum
!
!  iValues(iLen)  (in)  Array of integers.  The corresponding elements for each processor are summed and returnd into iReturn(iLen)
!  iLen           (in)  Length of the arrays
!  iReturn(iLen)  (out) Array of integers to get the results.

Subroutine MPIISum(iValues, iLen, iReturn)
   integer iValues, iReturn, iLen
   integer g, ierr,rc
   g=MPI_COMM_WORLD
   call MPI_ALLREDUCE(iValues,iReturn,iLen*4,MPI_INTEGER,MPI_SUM,g,ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if
end Subroutine

! MPIDSum
!
!  dValues(iLen)  (in)  Array of real*8s.  The corresponding elements for each processor are summed and returnd into dReturn(iLen)
!  iLen           (in)  Length of the arrays
!  dReturn(iLen)  (out) Array of real*8s to get the results.
Subroutine MPIDSum(dValues, iLen, dReturn)
   real*8 dValues, dReturn
   integer iLen
   integer g, ierr,rc
   g=MPI_COMM_WORLD
   call MPI_ALLREDUCE(dValues,dReturn,iLen*8,MPI_DOUBLE_PRECISION,MPI_SUM,g,ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if
end Subroutine

!MPIDSumRoot
!Same as above, but only updates the value on the root processor
Subroutine MPIDSumRoot(dValues,iLen,dReturn,Root)
    REAL*8 :: dValues,dReturn
    INTEGER :: iLen
    INTEGER :: g,ierr,rc,Root
    g=MPI_COMM_WORLD
    call MPI_REDUCE(dValues,dReturn,iLen,MPI_DOUBLE_PRECISION,MPI_SUM,Root,g,ierr)
    if(ierr.ne.MPI_SUCCESS) then
        print *,'Error summing values in MPIDSumRoot. Terminating.'
        call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
    endif
end subroutine MPIDSumRoot

! MPIHelSum
!
!  dValues(iLen)  (in)  Array of Type(HElement).  The corresponding elements for each processor are summed and returnd into dReturn(iLen)
!  iLen           (in)  Length of the arrays
!  dReturn(iLen)  (out) Array of Type(HElement) to get the results.
Subroutine MPIHelSum(dValues, iLen, dReturn)
   Use HElem
   Type(HElement) dValues(:), dReturn(:)
   integer iLen
   integer g, ierr,rc
   g=MPI_COMM_WORLD
   call MPI_ALLREDUCE(dValues,dReturn,iLen*8*HElementSize,MPI_DOUBLE_PRECISION,MPI_SUM,g,ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if
end Subroutine
End Module Parallel


! GetProcElectrons
!
!  iProcIndex  (in)     Index of this processor (starting at 1)
!  nProcessors (in)     Total number of processors
!  iMinElec    (out)    First electron to allocate to this processor
!  iMaxElec    (out)    Last electron to allocate to this processor
!
!  Choose min and max electrons such that ordered pairs are distributed evenly across processors
Subroutine GetProcElectrons(iProcIndex,nProcessors, iMinElec, iMaxElec)
   Use System, only: nEl
   implicit none
   integer iProcIndex,nProcessors, iMinElec,iMaxElec
   real*8 nCur
!Invert X=n(n-1)/2
   nCur=((nProcessors+1-iProcIndex)*nEl*(nEl-1.d0)/nProcessors)
   nCur=nEl+1-(1+sqrt(1.d0+4*nCur))/2
   iMinElec=ceiling(nCur)
   if(floor(nCur).eq.ceiling(nCur)) then
! We hit smack bang on an integer for nCur, so to avoid duplicating this electron here and in the previous processor, we choose the next one.
      iMinElec=iMinElec+1
   endif
   nCur=((nProcessors-iProcIndex)*nEl*(nEl-1.d0)/nProcessors)
   nCur=nEl+1-(1+sqrt(1.d0+4*nCur))/2
   iMaxElec=floor(nCur)
End

! Par2vSum
!
!  nI(nEl)        The root determinant of the 2v sum.
!
! A parallel version of the 2-vertex sum
!
! 
Subroutine Par2vSum(nI)
   USE HElem
   uSE MPI
   Use Parallel, only : iProcIndex, nProcessors,MPIHElSum
   Use System, only: nEl,Beta
   Use Determinants, only: HElement, GetHElement3
   IMPLICIT NONE
   Integer nI(nEl)
   integer iMinElec, iMaxElec
   real*8 nCur
   integer i
   integer nPr
   integer store(6)
   integer ic,exlen,iC0
   integer, pointer :: Ex(:)
   integer nJ(nEl)
   TYPE(HElement) dU
   Type(HDElement) dE1,dE2
   Type(HElement) dEw,dw,dEwtot,dwtot,dTots(2),dTots2(2)
   iC0=0
   nPr=nProcessors
   i=iProcIndex+1
   write(6,*) "Proc ",i,"/",nProcessors
   call GetProcElectrons(iProcIndex+1,nProcessors,iMinElec,iMaxElec)
   Write(6,*) "Electrons ", iMinElec, " TO ", iMaxElec

!  The root's energy
   dE1=GetHElement3(nI,nI,iC0)

!  Initialize.  If we're the first processor then we add in the 1-vertex graph.
   if(iProcIndex.EQ.0) THEN
      dEwTot=dE1
      dwTot=HElement(1.d0)
   ELSE
      dEwTot=0.D0
      dwTot=HElement(0.d0)
   ENDIF

! Now enumerate all 2v graphs
!.. Setup the spin excit generator
   STORE(1)=0
   CALL GENSYMEXCITIT3Par(NI,.TRUE.,EXLEN,nJ,IC,0,STORE,3,iMinElec,iMaxElec)
   Allocate(Ex(exLen))
   EX(1)=0
   CALL GENSYMEXCITIT3Par(NI, .TRUE.,EX,nJ,IC,0,STORE,3,iMinElec,iMaxElec)


!  Generate the first excitation
   CALL GENSYMEXCITIT3Par(NI, .False.,EX,nJ,IC,0,STORE,3,iMinElec,iMaxElec)
   i=0
!NJ(1) is zero when there are no more excitations.
   DO WHILE(NJ(1).NE.0)
      i=i+1
      dU=GetHElement3(nI,nJ,iC)
      dE2=GetHElement3(nJ,nJ,iC0)
      call Get2vWeightEnergy(dE1,dE2,dU,Beta,dw,dEw)
      dEwTot=dEwTot+dEw
      dwTot=dwTot+dw
      CALL GENSYMEXCITIT3Par(NI, .False.,EX,nJ,IC,0,STORE,3,iMinElec,iMaxElec)
   ENDDO 
   write(6,*) I
   write(6,*) dEwTot,dwTot, dEwTot/dwTot
   dTots(1)=dwTot
   dTots(2)=dEwTot
   Call MPIHelSum(dTots,2,dTots2)
   write(6,*) dTots2(2),dTots2(2),dTots2(2)/dTots2(1)
   deallocate(Ex)
End Subroutine


!Get the MP2 energy contribution for an excitation 1->2
subroutine getMP2E(dE1,dE2,dU,dE)
   implicit none
   real*8 dE1,dE2,dU,dE
   dE=dU*dU/(dE1-dE2)
end subroutine


!Get the two-vertex contribution from a graph given its matrix elements.  Weights are divided by exp(-beta E1)

!Get the two-vertex contribution from a graph given its matrix elements.  Weights are divided by exp(-beta E1)
subroutine Get2vWeightEnergy(dE1,dE2,dU,dBeta,dw,dEt)
   implicit none
   real*8 dE1,dE2,dU,dEt,dw,dBeta
   real*8 dEp,dEm,dD,dEx,dD2
   real*8 dTmp
   if(dU.eq.0) then
      dw=0.D0
      dEt=0.D0
      return 
   endif
   dD=(dE1+dE2)**2-4*dE1*dE2+4*dU*dU
   dD=sqrt(dD)/2
   dEp=(dE1+dE2)/2
   dEm=dEp-dD
   dEp=dEp+dD
!  The normalized first coefficient of an eigenvector is U/sqrt((dE1-dEpm)**2+U**2)
   dD=1/sqrt((dE1-dEp)**2+dU**2)
   dD2=dD*(dEp-dE1)
   dD=dD*dU
!dD is the eigenvector component
   dEx=exp(-dBeta*(dEp-dE1))
   dw=dD*dD*dEx
   dEt=dE1*dD*dD*dEx
   dEt=dEt+dU*dD*dD2*dEx
!   write(6,*) dEp,dD,dD2,dw,dEx,dBeta
!  This can be numerically unstable when dE1 is v close to dEm
!   Instead we just swap dD2 and dD around
!   dD=1/sqrt((dE1-dEm)**2+dU**2)
!   dD2=dD*(dEm-dE1)
!   dD=dD*dU
   dTmp=dD
   dD=dD2
   dD2=-dTmp
   dEx=exp(-dBeta*(dEm-dE1))
   dw=dw+dD*dD*dEx
!   write(6,*) dEm,dD,dD2,dw,dEx
   dEt=dEt+dE1*dD**2*dEx
   dEt=dEt+dU*dD*dD2*dEx
   dEt=dEt-dE1
   dw=dw-1
end subroutine
#else
module Parallel
!Dummy so we actually create a module
   IMPLICIT NONE
   save
   integer iProcIndex
   integer nProcessors
end module
#endif

