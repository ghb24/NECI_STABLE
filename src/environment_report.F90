subroutine environment_report(tCPMD)

!= Print out a summary of the environment:
!=   * When the code was compiled
!=   * The version control system (VCS) repository id of the codebase. 
!=   * Whether the codebase contains local changes.
!=   * The working directory.
!=   * The host computer.
!=
!= The VCS information is added via preproccessing and requires various settings
!= in the makefile.
!=
!= This has been tested with gfortran, g95, ifort, pgf90.
!= Note that hostnm and getcwd are intrinsic functions (at least according to gnu
!= documentation).  There should also exist the intrinsic subroutines, but pgf90
!= and ifort gave segmentation faults when they were used.

!= In:
!=   tCPMD: true if doing a CPMD-based calculation.  CPMD already prints out the 
!=          directory and host names, so we avoid printing duplicate information.

implicit none
logical :: tCPMD
integer :: stat,getcwd,hostnm
character(255) :: dirname,host

write (6,'(/,1X,64("="))')
write (6,'(a13,a,a4,a)') 'Compiled on ',__DATE__,'at ',__TIME__
write (6,'(a29,/,5X,a)') 'VCS BASE repository version:',_VCS_VER
#ifdef _WORKING_DIR_CHANGES
write (6,'(a42)') 'Working directory contains local changes.'
#endif
stat=getcwd(dirname)
if (stat.eq.0.and..not.tCPMD) then
    write (6,'(a20)') 'Working directory: '
    write (6,'(5X,a)') trim(dirname)
end if
stat=hostnm(host)
if (stat.eq.0.and..not.tCPMD) then
    write (6,'(a13,a)') 'Running on: ',trim(host)
end if
write (6,'(1X,64("="),/)')

return 
end subroutine environment_report
