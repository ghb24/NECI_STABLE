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
!  In:
!     error_str: parameter string containing error used as argument to STOP.
!  
!  Abort all processors
Subroutine MPIStopAll(error_str)
   character(3) :: error_str
   integer error_code,ierror
   ! errorcode: Error returned to invoking environment.
   ! ierror: error status (of abort: was abort successful?)
   ! Currently neither are analysed.
   call MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
   stop error_str
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
! JSS: This causes problems when the nEl/nProc is an integer.
! More thinking needed Alex! :-)
!   if(floor(nCur).eq.ceiling(nCur)) then
! We hit smack bang on an integer for nCur, so to avoid duplicating this electron here and in the previous processor, we choose the next one.
!      iMinElec=iMinElec+1
!   endif
   nCur=((nProcessors-iProcIndex)*nEl*(nEl-1.d0)/nProcessors)
   nCur=nEl+1-(1+sqrt(1.d0+4*nCur))/2
   iMaxElec=floor(nCur)
End subroutine GetProcElectrons



subroutine ParMP2(nI)
! Calculate the MP2 energy in parallel.
! Designed for CPMD calculations.
! Parallel used only for distribution of the wavefunctions and to
! compute the integrals.
! In:
!    nI(nEl) : list of occupied spin orbitals in the reference 
!              determinant.
! Prints out the <D_0|H|D_0>+MP2 energy.

! To do: 
! 1. Write a custom routine to calculate the <D_0|H|D_1> elements.
! This would then enable the excitations ij->ab,iJ->aB,Ij-Ab and IJ->AB
! to be handled at once, removing the need to any integrals to be stored.
! 2. Implement for standalone NECI calculations (molecular, Hubbard etc).

   USE HElem
   uSE MPI
   Use Parallel, only : iProcIndex, nProcessors,MPIHElSum
   Use System, only: nBasisMax,nEl,Beta,ARR,nBasis,ECore,G1,AreSameSpatialOrb,tCPMD,Symmetry
   use Calc, only: NWHTAY
   use Integrals, only: GetUMatEl2
   use UMatCache, only: GTID
   use OneEInts, only: GetTMatEl
   Use Determinants, only: HElement, GetHElement3,GetH0Element3
   use MemoryManager, only: LogMemAlloc,LogMemDealloc
   use SymData, only: SymLabels
   IMPLICIT NONE
   integer :: nI(nEl)
   integer :: iMinElec, iMaxElec
   real*8 :: nCur
   integer :: i,j
   integer :: IA,JA,AA,BA
   integer :: store(6),Excit(2,2)
   integer :: ic,exlen,iC0,ExcitOrbs(2,2),ExLevel
   integer, pointer :: Ex(:)
   integer :: nJ(nEl),weight
   TYPE(HElement) dU(2)
   Type(HDElement) :: dE1,dE2
   Type(HElement) :: dE,dEtot(2),dE0
   ! MPIHelSum requires arrays as input/output.
   Type(HElement) :: dEarr(2),dEres(2)
   type(Symmetry) :: iSym1,iSym2
   type(Symmetry) :: iSym1Conj,iSym2Conj
   type(Symmetry) :: SymConj
   logical :: tSign
   integer :: isub,ierr,tag_Ex
   character(*), parameter :: this_routine='ParMP2'

   call TiSet('ParMP2    ',isub)
   
   select case(IAND(nWHTay(1,1),24))
   case(0)
       ! Generate both single and double excitations.  CPMD works in a Kohn--Sham
       ! basis, and so Brillouin's theorem does not apply.
       write (6,*) 'ParMP2: using single and double excitation.'
       ExLevel=3
   case(8)
       write (6,*) 'ParMP2: using only single excitations.'
       ExLevel=1
   case(16)
       write (6,*) 'ParMP2: using only double excitations.'
       ExLevel=2
   case(24)
       call stop_all('ParMP2','Invalid combination of flags in nWHTay.  Invalid EXCITAIONS specification?')
   end select

   iC0=0

   write(6,*) "Proc ",iProcIndex+1,"/",nProcessors
   if (tCPMD) then
       ! For CPMD jobs, we actually want each processor to do the full sum, as each
       ! integral is divided across the processors.
       iMinElec=1
       iMaxElec=nEl
   else
       ! For other calculations, the sum is split over processors.
       call GetProcElectrons(iProcIndex+1,nProcessors,iMinElec,iMaxElec)
   end if
   write(6,*) "Electrons ", iMinElec, " TO ", iMaxElec

!  The root's "energy"---sum over the eigenvalues of occ. spin orbitals and the
!  core energy. Only used it the H0Element formulation is used (see below). 
   dE1=GetH0Element3(nI)

!  Initialize: get the contribution from the reference determinant.
   dE0=GetHElement3(nI,nI,iC0)

!  Now enumerate all 2v graphs
!  Setup the spin excit generator
   STORE(1)=0
!  IC is the excitation level (relative to the reverence det).
   CALL GENSYMEXCITIT3Par(NI,.TRUE.,EXLEN,nJ,IC,0,STORE,ExLevel,iMinElec,iMaxElec)
   Allocate(Ex(exLen),stat=ierr)
   call LogMemAlloc('Ex',Exlen,4,this_routine,tag_Ex,ierr)
   EX(1)=0
   CALL GENSYMEXCITIT3Par(NI, .TRUE.,EX,nJ,IC,0,STORE,ExLevel,iMinElec,iMaxElec)

!  Generate the first excitation
   CALL GENSYMEXCITIT3Par(NI, .False.,EX,nJ,IC,0,STORE,ExLevel,iMinElec,iMaxElec)
   i=0
   j=0

!NJ(1) is zero when there are no more excitations.
   DO WHILE(NJ(1).NE.0)

      i=i+1

! Quickest attempt (and useful for debugging).
! Not as efficient as the code below, but clearer and useful for debugging.
! Also, this should be used for unrestricted calculations.
      ! MP2 theory refers to the unperturbed excited determinant
      ! => use GetH0Element3 rather than GetHElement3.
!      dE2=GetH0Element3(nJ)
!      dU=GetHElement3(nI,nJ,iC)
!      call getMP2E(dE1,dE2,dU,dE)
!      dETot=dETot+dE

! Alternatively, calculate the energy of the excited determinant
! in reference to that of the reference determinant (i.e. setting dE1=0).
      Excit(1,1)=2
      call GetExcitation(nI,nJ,nEl,Excit,tSign)

      ! Assuming a restricted calculation.
      weight=0
      if (Excit(1,2).eq.0) then
          ! Single excitation
          if (G1(Excit(1,1))%Ms.eq.-1) then
              ! alpha -> alpha single excitation.
              ! count for beta->beta as well.
              weight=2
          end if
      else
          ! Double excitation.
          if (G1(Excit(1,1))%Ms.eq.-1.and.G1(Excit(1,2))%Ms.eq.-1) then
              ! alpha,alpha -> alpha,alpha double excitation.
              ! count for beta,beta  -> beta,beta as well.
              weight=2
          else if (G1(Excit(1,1))%Ms.eq.-1.and.G1(Excit(1,2))%Ms.eq.1) then
              ! alpha,beta -> alpha,beta double excitation.
              ! We consider just these, and bring in the beta,alpha -> beta,alpha
              ! excitations via spin symmetry.
              if (AreSameSpatialOrb(Excit(1,1),Excit(1,2))) then
                  ! Excitations from the spin orbitals of the same spatial
                  ! orbital are unique.
                  ! e.g (1a,1b) -> (2a,2b)
                  ! e.g (1a,1b) -> (2a,3b)
                  weight=1
              else
                  ! Excitations from different spatial orbitals.  Use spin
                  ! symmetry:
                  ! (1a,2b) -> (3a,3b) has an identical contribution to the MP2
                  ! as (1b,2a) -> (3a,3b).
                  ! Similarly (1a,2b) -> (3a,4b) and (1b,2a) -> (3b,4a).
                  ! Count twice for excitations from beta,alpha spin-orbitals.
                  if (AreSameSpatialOrb(Excit(2,1),Excit(2,2))) then
                      weight=2
                      !write (6,'(2i4,a2,2i4,2f17.8)') Excit(1,:),'->',Excit(2,:),GetHElement3(nI,nJ,iC)
                  end if

              end if
          end if
      end if

      write (6,'(2i4,a2,2i4)') Excit(1,:),'->',Excit(2,:)
      call flush(6)
      if (weight.eq.0) then
          ! Have counted current excitation elsewhere.
          ! Generate next excitation.
          write (6,*) 'skipping excitation'
          CALL GENSYMEXCITIT3Par(NI,.false.,EX,nJ,IC,0,STORE,ExLevel,iMinElec,iMaxElec)
          cycle
      end if

      j=j+1
      dE2=HDElement(Arr(Excit(2,1),2)-Arr(Excit(1,1),2))
      dU=HElement(0.d0)
      if (Excit(2,2).ne.0) then

          ! Double excitation
          dE2=dE2+HDElement(Arr(Excit(2,2),2)-Arr(Excit(1,2),2))

          call GTID(nBasisMax,Excit(1,1),IA)
          call GTID(nBasisMax,Excit(1,2),JA)
          call GTID(nBasisMax,Excit(2,1),AA)
          call GTID(nBasisMax,Excit(2,2),BA)
          if (G1(Excit(1,1))%Ms.eq.G1(Excit(1,2))%Ms) then
              dU(2)=GetUMatEl2(IA,JA,AA,BA)
              dU(1)=-GetUMatEl2(IA,JA,BA,AA)
          else if (G1(Excit(1,1))%Ms.eq.G1(Excit(2,1))%Ms) then
              dU(1)=GetUMatEl2(IA,JA,AA,BA)
          else
              dU(1)=-GetUMatEl2(IA,JA,BA,AA)
          end if
          !write (6,'(2i4,a2,2i4,2f17.8)') Excit(1,:),'->',Excit(2,:),dU(2)
          if (tSign) dU(1)=-dU(1)
      else

          ! Single excitation.
          ! dU=\sum_J 2<IJ|AJ>-<IJ|JA> (in terms of spatial orbitals).
          call GTID(nBasisMax,Excit(1,1),IA)
          call GTID(nBasisMax,Excit(2,1),AA)
          do JA=1,nEl/2
              if (JA.eq.IA) then
                  dU(1)=dU(1)+GetUMatEl2(IA,JA,AA,JA)
              else
                  dU(1)=dU(1)+HElement(2)*GetUMatEl2(IA,JA,AA,JA)-GetUMatEl2(IA,JA,JA,AA)
              end if
          end do
          dU(1)=dU(1)+GetTMATEl(Excit(1,1),Excit(2,1))
          if (tSign) dU(1)=-dU(1)

      end if

!      write (6,'(2i3,a2,2i3,6f12.7)') Excit(1,:),'->',Excit(2,:),dU,dE2,dE

      if (Excit(1,2).eq.0) then
          call getMP2E(HDElement(0.d0),dE2,dU(1),dE)
          ! Singles contribution.
          dETot(1)=dETot(1)+HElement(weight)*dE
      else
          if (dU(2).agt.0.d0) then
              write (6,'(2i4,a2,2i4,2f17.8)') Excit(1,:),'->',Excit(2,:),dU(1)
              call getMP2E(HDElement(0.d0),dE2,dU(1),dE)
              dETot(2)=dETot(2)+HElement(weight)*dE
              write (6,'(2i4,a2,2i4,2f17.8)') Excit(1,:),'->',Excit(2,:),dU(2)
              call getMP2E(HDElement(0.d0),dE2,dU(2),dE)
              dETot(2)=dETot(2)+HElement(weight)*dE
!              write (6,*) 'dE',dETot(2)
          end if
          dU(1)=dU(1)-dU(2)
          ! Doubles contribution.
          call getMP2E(HDElement(0.d0),dE2,dU(1),dE)
          dETot(2)=dETot(2)+HElement(weight)*dE
          !write (6,'(2i4,a2,2i4)') Excit(1,:),'->',Excit(2,:)
          !write (6,*) 'dE',dETot(2)
      end if
! END  of more efficient approach.

      ! Get next excitation.
      CALL GENSYMEXCITIT3Par(NI,.false.,EX,nJ,IC,0,STORE,ExLevel,iMinElec,iMaxElec)

   ENDDO 

   write(6,*) 'No. of excitations=',I
   write(6,*) 'No. of spin-unique excitations=',J
   if (.not.tCPMD) then
       write (6,'(a28,i3,a1,2f15.8)') 'Contribution from processor',iProcIndex+1,':',dEtot
       dEarr=dETot
       call MPIHElSum(dEArr,2,dEres)
   end if
   if (iand(ExLevel,1).eq.1) write(6,*) 'MP2 SINGLES=',dETot(1)+dE0
   if (iand(ExLevel,2).eq.2) write(6,*) 'MP2 DOUBLES=',dETot(2)+dE0
   write(6,*) 'MP2 ENERGY =',dETot(1)+dETot(2)+dE0

   deallocate(Ex)
   call LogMemDealloc(this_routine,tag_Ex)

   call TiHalt('ParMP2    ',isub)

end subroutine ParMP2



! Par2vSum
!
!  nI(nEl)        The root determinant of the 2v sum.
!
! A parallel version of the 2-vertex sum
!
! This is not quite stable yet.
! Issues:
!  * Some problems remain with how electrons are distributed over processors.
!  * Doesn't work for CPMD calculations.
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
   integer store(6)
   integer ic,exlen,iC0
   integer, pointer :: Ex(:)
   integer nJ(nEl)
   TYPE(HElement) dU
   Type(HDElement) dE1,dE2
   Type(HElement) dEw,dw,dEwtot,dwtot,dTots(2),dTots2(2)
   iC0=0
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
!  IC is the excitation level (relative to the reverence det).
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
   write(6,*) dEwTot,dwTot,dEwTot/dwTot
   dTots(1)=dwTot
   dTots(2)=dEwTot
   Call MPIHelSum(dTots,2,dTots2)
   write(6,*) dTots2(2),dTots2(1),dTots2(2)/dTots2(1)
   deallocate(Ex)
End Subroutine


subroutine getMP2E(dE1,dE2,dU,dE)
   ! Get the MP2 energy contribution for an excitation 1->2
   ! In:
   !    dE1 energy of reference determinant.
   !    dE2 energy of excited determinant.
   !    dU  Cross term of Hamiltonian matrix, < 1 | H | 2 >.
   ! Out:
   !    dE = |< 1 | H | 2 >|^2 / (dE2 - dE1)
   !         contribution to the MP2 energy.
   use HElem
   implicit none
   type(HDElement) dE1,dE2
   type(HElement) dU,dE
   dE=Sq(dU)/(dE1%v-dE2%v)
end subroutine


! Get the two-vertex contribution from a graph containing the reference
! determinant, 0, and a connected determinant, i, given its matrix elements.
! Weights are divided by exp(-beta E1)

! Each two vertex graph is represented by a 2x2 (Hermitian) matrix
!  | dE1   dU  |
!  | dU*   dE2 |
! The eigenvalues are [ (dE1+dE2) +- \Sqrt( (dE1+dE2) - 4dE1*dE2 + 4|dU|^2 ) ]/2.
! where:
!   dE1=<D_0|H|D_0>
!   dE2=<D_i|H|D_i>
!   dU =<D_0|H|D_i>
! Denoting the eigenvalues as dEp and dEm, the normalised eigenvectors are:
!   \frac{1}{\Sqrt{ dU^2 + (dE1-dEp)^2 } ( U, dEp-dE1 )
!   \frac{1}{\Sqrt{ dU^2 + (dE1-dEm)^2 } ( U, dEm-dE1 )

subroutine Get2vWeightEnergy(dE1,dE2,dU,dBeta,dw,dEt)
   ! In:
   !    dE1    <D_0|H|D_0>
   !    dE2    <D_i|H|D_i>
   !    dU     <D_0|H|D_i>
   !    dBeta  beta
   ! Out:
   !    dw     weight of the graph, w_i[G]
   !    dEt    weighted contribution of the graph, w_i[G] \tilde{E}_i[G]
   use HElem
   implicit none
   type(HDElement) dE1,dE2
   type(HElement) dU,dEt,dw
   real*8 dBeta
   type(HElement) dEp,dEm,dD,dEx,dD2,dTmp
   if(abs(dU%v).eq.0.d0) then
      ! Determinants are not connected.
      ! => zero contribution.
      dw=0.D0
      dEt=0.D0
      return 
   endif

!  Calculate eigenvalues.
!   write (6,*) dE1,dE2,dU

   dD=HElement((dE1%v+dE2%v)**2-4*(dE1%v*dE2%v)+4*abs(dU%v)**2)

!   write (6,*) dD

   dD=sqrt(dD%v)/2
   dEp=(dE1+dE2)/2
   dEm=dEp-dD
   dEp=dEp+dD

!   write (6,*) dD,dEp,dEm

!  The normalized first coefficient of an eigenvector is U/sqrt((dE1-dEpm)**2+U**2)
   dD=1/sqrt((dE1%v-dEp%v)**2+abs(dU%v)**2) ! normalisation factor
   dD2=dD%v*(dEp%v-dE1%v)               ! second coefficient
   dD=dD*dU                           ! first coefficient
   
!   write (6,*) dD,dD2

!dD is the eigenvector component
   dEx=exp(-dBeta*(dEp%v-dE1%v))
   dw=sq(abs(dD))*dEx%v
   dEt=dE1%v*sq(abs(dD))*dEx%v
   dEt=dEt+dU*dD*dconjg(dD2)*dEx
   
!   write (6,*) dEx,dw,dEt

!   write(6,*) dEp,dD,dD2,dw,dEx,dBeta
!  This can be numerically unstable when dE1 is v close to dEm:
!      dD=1/sqrt((dE1-dEm)**2+dU**2)
!      dD2=dD*(dEm-dE1)
!      dD=dD*dU
!  Instead we just swap dD2 and dD around
   dTmp=dD
   dD=dconjg(dD2)
   dD2=-dconjg(dTmp)
   dEx=exp(-dBeta*(dEm%v-dE1%v))
   dw=dw+HElement(sq(abs(dD))*dEx%v)
!   write(6,*) dEm,dD,dD2,dw,dEx
   dEt=dEt+HElement(dE1%v*sq(abs(dD))*dEx%v)
   dEt=dEt+dU*dD*dconjg(dD2)*dEx

!  Now, we already have calculated the effect of the one-vertex graph, so
!  we need to remove the double-counting of this from the 2-vertex graph.
!  For the one-vertex graph containing just i:
!     w_i[i] = exp(-\beta E_i)
!     w_i[i] \tilde{E}_i[i] = exp(-\beta E_i) E_i
!  As we factor out exp(-\beta E_i), this becomes just:
   dEt=dEt-HElement(dE1%v)
   dw=dw-HElement(1)
   write (6,*) 'wE,E',dEt,dw
   
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

