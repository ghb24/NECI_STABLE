subroutine environment_report(tCPMD)

!= Print out a summary of the environment:
!=   * When the code was compiled
!=   * The version control system (VCS) repository id of the codebase. 
!=   * Whether the codebase contains local changes.
!=   * The working directory.
!=   * The host computer.
!=   * The time the calculation started.
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

#ifdef NAGF95
USe f90_unix_dir
use constants, only: int32
#endif
implicit none
logical :: tCPMD
integer :: stat,hostnm
#ifndef NAGF95
integer :: getcwd
#else
integer(kind=int32) :: stat_dum
#endif
character(255) :: dirname,host
integer :: date_values(8)

#ifndef MOLPRO
write (6,'(/,1X,64("="))')
write (6,'(a13,a,a4,a)') 'Compiled on ',__DATE__,'at ',__TIME__
write (6,'(a30,/,5X,a)') 'Compiled using configuration:',_CONFIG
write (6,'(a29,/,5X,a)') 'VCS BASE repository version:',_VCS_VER
#endif
#ifdef _WORKING_DIR_CHANGES
write (6,'(a42)') 'Working directory contains local changes.'
#endif
#ifdef NAGF95
call getcwd(path=dirname,errno=stat_dum)
stat=stat_dum
#else
stat=getcwd(dirname)
#endif
if (stat.eq.0.and..not.tCPMD) then
    write (6,'(a20)') 'Working directory: '
    write (6,'(5X,a)') trim(dirname)
end if
#if defined(NAGF95) || defined(_WIN32_)
!Can't find a hostnm intrinsic equivalent in the nag system modules
    stat=1
#else
    stat=hostnm(host)
#endif
if (stat.eq.0.and..not.tCPMD) then
    write (6,'(a13,a)') 'Running on: ',trim(host)
end if

call date_and_time(VALUES=date_values)

write (6,'(1X,"Started running on",1X,i2.2,"/",i2.2,"/",i4.4,1X,"at",1X,i2.2,2(":",i2.2))') date_values(3:1:-1), date_values(5:7)

write (6,'(1X,64("="),/)')

return 
end subroutine environment_report
