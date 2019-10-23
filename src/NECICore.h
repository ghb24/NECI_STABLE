! Can be included into free- and fixed-format code


      interface
      subroutine necicore(icacheflag, tcpmd, tvasp, tmolpro_local,      &
     &      call_as_lib, int_name, filename_in, MemSize)
         use, intrinsic :: iso_fortran_env, only: int64
!= NECICore is the main outline of the NECI Program.
!= It provides a route for calling NECI when accessed as a library,
!= rather than as a standalone program.
!= In:
!=   iCacheFlag: controls the behaviour of the 4-index integral cache.
!=               Currently relevant only for CPMD and VASP calculations.
!=               iCacheFlag=0: initialise and destroy the cache.
!=               iCacheFlag=1: initialise but don't destroy the cache.
!=               iCacheFlag=2: reuse and destroy the cache.
!=               iCacheFlag=3: reuse and keep the cache.
!=   tCPMD: True if doing a CPMD-based calculation.
!=   tVASP: True if doing a VASP-based calculation.
!=   call_as_lib: True if called as subroutine from external code.
!=   int_name is the name of the integral file to read in if necessary
!=   filename is the name of the input file to read in if necessary
        implicit none
        integer, intent(in), optional :: iCacheFlag
        logical,intent(in), optional ::                                 &
     &      tCPMD,tVASP,tMolpro_local,call_as_lib
        character(*), intent(in), optional :: filename_in, int_name
        integer(int64), intent(in), optional :: MemSize
      end subroutine
      end interface
