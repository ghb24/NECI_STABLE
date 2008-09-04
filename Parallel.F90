module Parallel

!=  MPI interface and helper routines for parallel work.
!=  If compiled without the PARALLEL c pre-processing statement, then
!=  contains dummy routines to enable the parallel calculation algorithms to
!=  be used in serial.  This is useful for testing, development and as (in some 
!=  cases) the parallel algorithms are much more efficient that the serial
!=  analogues, as the latter are much more general.

!=  NECI will run as a standalone parallel app or in conjunction with CPMD.
!=  Only standalone is currently implemented.

!=  Parallelization is done over occupied electrons.  These are split among each
!=  processor such that there are an approximately equal number of pairs of
!=  electrons on each processor.
!=
!=  Each processor is given a set of 'Electron 1'.  From this it can generate
!=  a set of single excitations as well as a set of double excitations.  Double
!=  excitations can have any electron past Electron 1 as Electron 2.  This means
!=  that the lower electron numbers will have more possible pairs.
!=
!=  Parallelization is supported by the symmetry excitation generators,
!=  interfaced through GenSymExcitIt3Par There is no clean way to automatically
!=  parallelize high-vertex graph code, so each parallel routine must be
!=  specifically written.
!=
!=  ROUTINES
!=     MPIInit     Setup MPI and init Nodefiles if we're standalone
!=     MPIEnd      Shutdown MPI
!=     MPIStopAll  Abort all processors
!=     MPIISum     Sum an array of integers among all processors, and distribute the results array to each.
!=     MPIDSum     As MPIISum but for real*8
!=     MPIHelSum   As MPIDSum but for Type(HElement)

   ! mpi USE is in mixed case to avoid it being picked up by the the Configure
   ! script, as it doesn't require a module file.
#ifdef PARALLEL
   uSE mpi
#endif
   IMPLICIT NONE
   save
   integer iProcIndex
   integer nProcessors

Contains



Subroutine MPIInit(tExternal)
   != Determine the number of processors, and fork each off to its own NodeFile output file
   !=
   != In:
   !=   tExternal True if using VASP/CPMD's MPI interface, so we don't have to initialise our own.
   Use Determinants, only: FDet
   implicit none
   logical tExternal
   integer numtasks, rank, ierr, rc
   integer a,b,g
   character*20 NodeFile
#if PARALLEL
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
#endif
   
end Subroutine MPIInit



Subroutine MPIEnd(tExternal)
   !=  Shutdown our MPI Interface if we're not using CPMD/VASP's
   !=
   != In:
   !=   tExternal Set if using an external program's MPI interface
   !=             (currently CPMD or VASP), in which case the external
   !=             program handles MPI termination.
   implicit none
   logical tExternal
   integer ierr
#if PARALLEL
   if(.not.tExternal) then
      call MPI_FINALIZE(ierr)
   endif
#endif
end subroutine MPIEnd



Subroutine MPIStopAll(error_str)
   !=  Abort all processors.
   !=  
   !=  In:
   !=     error_str: parameter string containing error used as argument to STOP.
   character(3) :: error_str
   integer error_code,ierror
#if PARALLEL
   ! errorcode: Error returned to invoking environment.
   ! ierror: error status (of abort: was abort successful?)
   ! Currently neither are analysed.
   call MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
   stop error_str
#endif
end subroutine MPIStopAll



Subroutine MPIISum(iValues, iLen, iReturn)
   !=  In:
   !=     iValues(iLen)  Array of integers.  The corresponding elements for each
   !=                    processor are summed over the processors and returned in 
   !=                    iReturn(iLen).
   !=     iLen           Length of the arrays.
   !=  Out:
   !=     iReturn(iLen)  Array of integers to get the results.
   integer iValues, iReturn, iLen
   integer g, ierr,rc
#if PARALLEL
   g=MPI_COMM_WORLD
   call MPI_ALLREDUCE(iValues,iReturn,iLen*4,MPI_INTEGER,MPI_SUM,g,ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if
#else
   iReturn=iValues
#endif
end Subroutine MPIISum



Subroutine MPIDSum(dValues, iLen, dReturn)
   !=  In:
   !=     dValues(iLen)  Array of real*8s.  The corresponding elements for each
   !=                    processor are summed and returnd into dReturn(iLen).
   !=     iLen           Length of the arrays.
   !=  Out:
   !=     dReturn(iLen)  Array of real*8s to get the results.
   real*8 dValues, dReturn
   integer iLen
   integer g, ierr,rc
#ifdef PARALLEL
   g=MPI_COMM_WORLD
   call MPI_ALLREDUCE(dValues,dReturn,iLen*8,MPI_DOUBLE_PRECISION,MPI_SUM,g,ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if
#else
    dReturn=dValues
#endif
end Subroutine MPIDSum



Subroutine MPIDSumRoot(dValues,iLen,dReturn,Root)
    !=  Same as above, but only updates the value on the root processor.
    REAL*8 :: dValues,dReturn
    INTEGER :: iLen
    INTEGER :: g,ierr,rc,Root
#ifdef PARALLEL
    g=MPI_COMM_WORLD
    call MPI_REDUCE(dValues,dReturn,iLen,MPI_DOUBLE_PRECISION,MPI_SUM,Root,g,ierr)
    if(ierr.ne.MPI_SUCCESS) then
        print *,'Error summing values in MPIDSumRoot. Terminating.'
        call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
    endif
#else
    dReturn=dValues
    Root=dValues
#endif
end subroutine MPIDSumRoot



Subroutine MPIHElSum(dValues, iLen, dReturn)
   !=  In:
   !=     dValues(iLen)  Array of Type(HElement).  The corresponding elements
   !=                    for each processor are summed and returnd into dReturn(iLen).
   !=     iLen           Length of the arrays.
   !=  Out:
   !=     dReturn(iLen)  (out) Array of Type(HElement) to get the results.
   Use HElem
   Type(HElement) dValues(:), dReturn(:)
   integer iLen
   integer g, ierr,rc
#ifdef PARALLEL
   g=MPI_COMM_WORLD
   call MPI_ALLREDUCE(dValues,dReturn,iLen*8*HElementSize,MPI_DOUBLE_PRECISION,MPI_SUM,g,ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if
#else
   dReturn=dValues
#endif
end Subroutine MPIHElSum



Subroutine GetProcElectrons(iProcIndex,nProcessors, iMinElec, iMaxElec)
   !=  Choose min and max electrons such that ordered pairs are distributed evenly across processors
   !=
   !=  In:
   !=     iProcIndex  Index of this processor (starting at 1).
   !=     nProcessors Total number of processors.
   !=  Out:
   !=     iMinElec    First electron to allocate to this processor.
   !=     iMaxElec    Last electron to allocate to this processor.
   use SystemData, only: nEl
   implicit none
   integer iProcIndex,nProcessors, iMinElec,iMaxElec
   real*8 nCur
#ifdef PARALLEL
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
#else
   ! Serial calculation: all electrons on one processor.
   iMinElec=1
   iMaxElec=1
#endif
End subroutine GetProcElectrons



End Module Parallel
