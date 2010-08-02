!Miscellaneous MODULEs and SUBROUTINES common to the code go here.
!-----------------------------------------------------------------

module common_routines
use precision
!use global_data     !AJWT
private
public :: check_info, check_a_file, check_allocate, check_deallocate
public :: getunit, my_timer, get_available_file
public :: name_tmp_file, name_mol_file
public :: internal_error, general_error, close_file_units
public :: warning
contains

subroutine check_info(info,from_routine,called_by_routine,line_no,do_i_stop)
implicit none
integer, intent(in) :: info
character(len=*), intent(in) :: from_routine
character(len=*), intent(in) :: called_by_routine
integer, intent(in) :: line_no
logical, intent(in), optional :: do_i_stop
!------------
logical :: i_stop
!------------
if (present(do_i_stop)) then
  i_stop = do_i_stop
else
  i_stop = .true.
endif
!
if (info.lt.0) then
  write(6,*)'ERROR code ',info,' from routine ',trim(from_routine)
  write(6,*)'which was called by routine ',trim(called_by_routine)
  write(6,*)'at line number ',line_no
  if (i_stop) then
    call internal_error('check_info',__LINE__,'INFO flagged error')
  endif
endif
!if ((info.gt.0).and.g_debug) then     AJWT
!  print *,'WARNING code ',info,' passed by ',trim(from_routine)
!  write(6,*)'which was called by routine ',trim(called_by_routine)
!  write(6,*)'at line number ',line_no
!endif
return
end subroutine check_info

!--------------------------------------------------------------------------

subroutine check_a_file(fname,info,flag)
!flag = 0 : normal working
!     = 1 : supress warnings
implicit none
character(len=*), intent(in) :: fname
integer, intent(out) :: info
integer, intent(in), optional :: flag
!------------
logical :: iexist
integer :: lflag
!------------
if (present(flag)) then
  lflag = flag
else
  lflag = 0
endif
!
inquire(file=fname,exist=iexist)
if (iexist) then
  info = 0
else
  info = -1
  select case(lflag)
  case(0)
    write(6,*)'WARNING: File ',trim(fname),' is needed but is not present'
  case(1)
  case default
  end select
endif
return
end subroutine check_a_file

!----------------------------------------------------------------------

subroutine get_available_file(file_prefix,file_suffix,filename,fileunit,info,&
     & overwrite)
!INPUT:
!  file_prefix,file_suffix: prefix and suffix of file name
!OUTPUT:
!  filename: = file_prefix//file_suffix unless this file already exists in which
!           case it sill be
!              file_prefixNNNfile_suffix 
!           where NNN is a 3 digit integer
!  fileunit: unit number available for file to be opened
!  info: 0 = all's well
!       -1 = no available file for NNN between 0 and 999
!            In this case, NNN=999
!  overwrite: LOGICAL (optional) : If TRUE, if the default file name exists,
!            overwrite it.
!use global_data, only : g_overwrite   !AJWT
implicit none
character(80), intent(in) :: file_prefix, file_suffix
character(80), intent(out) :: filename
integer, intent(out) :: fileunit
integer, intent(out) :: info
logical, optional, intent(in) :: overwrite
!------------
character(80) :: myfilename
character(3) :: string_i
integer :: i, myinfo, iunit
logical :: got_file, overwrite_file
!------------
if (present(overwrite)) then
  overwrite_file = overwrite
!else       !AJWT
!  overwrite_file = g_overwrite
endif
!
filename = ''
fileunit = -1
got_file = .false.
!
myfilename = trim(file_prefix)//trim(file_suffix)
call check_a_file(myfilename,myinfo,flag=1)
if ((myinfo.eq.-1).or.overwrite_file) then
  !got a file name not in use or has been requested to be overwritten so
  !we have our filename
  filename = myfilename
  !If it is to be overwritten, first close it if it is open:
  inquire(file=filename,number=iunit)
  if (iunit.gt.0) close(iunit)
  got_file = .true.
else
  !Look for an alternative name
  file_loop: do i = 1, 999
    call number2ascii(i,string_i)
    myfilename = trim(file_prefix)//'_'//trim(string_i)//trim(file_suffix)
    call check_a_file(myfilename,myinfo,flag=1)
    if (myinfo.eq.-1) then
      !got a file name not in use
      got_file = .true.
      exit file_loop
    endif
  enddo file_loop
  filename = myfilename
  !
endif
!
call getunit(fileunit)
!
if (got_file) then
  info = 0
else
  info = -1
endif
!
return
end subroutine get_available_file

!----------------------------------------------------------------------

subroutine name_mol_file(mol_name,filename,info)
!INPUT:
!  mol_name = name of the molecule 
!OUTPUT:
!  filename: = <mol_name>_NNN.mol   where NNN = 3 digit integer
!  info: 0 = all's well
!       -1 = no available file for NNN between 0 and 999
!            In this case, NNN=999
implicit none
character(80), intent(in) :: mol_name
character(80), intent(out) :: filename
integer, intent(out) :: info
!------------
character(4), parameter :: file_suffix = '.mol'
character(80) :: myfilename
character(3) :: string_i
integer :: i, myinfo
!------------
filename = 'NO_FILE_SET_FOR_MOLECULE'
info = -1
!
file_loop: do i = 0, 999
  call number2ascii(i,string_i)
  myfilename = trim(mol_name)//'_'//trim(string_i)//trim(file_suffix)
  call check_a_file(myfilename,myinfo,flag=1)
  if (myinfo.eq.-1) then
    !got a file name not in use
    filename = myfilename
    info = 0
    exit file_loop
  endif
enddo file_loop
!
return
end subroutine name_mol_file

!----------------------------------------------------------------------

subroutine name_tmp_file(tmpfile,operation,info,filetag,inlist)
!USAGE:
!  call name_tmp_file(tmpfile,'get') = get a name. See note below.
!  call name_tmp_file(tmpfile,'rel') = release a name
!  call name_tmp_file(tmpfile,'add') = add a name to the list. Enforces the
!                                      chosen name and flags an error if there
!                                      is a conflict.
!  call name_tmp_file(tmpfile,'inq',inuse) = inquire if the file is in use.
!                                      inlist = .true. if file is in list 
!                                             = .false. if not.
!  call name_tmp_file(tmpfile,'out') = write out file information to standard
!                                      output.
!The optional third argument governs stopping:
!If present, no stops are executed here, rather, if an error occurs, then 
!info < 0 is returned, else, info = 0 is returned.
!
!INPUT/OUTPUT: 
!  tmpfile    character(10)
!INPUT:
!  operation  character(3)       Operation to be performed
!  info       integer            Returns error status (optional)
!  filetag    character(80)      Adds tag to file list to make identification
!                                possible. (optional)
!
!Given an input filename <tmpfile>, this routine attempts to find a suitable
!file name that can be used for a temporary file. The criteria used are:
! tmpfile = 
! 'TMP_______' :: find an available file of the form TMP_<number>
!                 where <number> is a 3 digit number starting from 0.
! ''           :: same as above
!
! Default      :: Retain the name of the TMP file unless there is a
!                 conflict with another file. If so, replace the last 3
!                 characters of the name by numbers starting from 000 till an
!                 unused name is found.
!                 This is similar to using 'add' as the operation, only, 'add'
!                 enforces the same name and flags an error if there is a
!                 conflict.
use parameters, only : max_tmp_files, tmp_file_name
implicit none
character(10), intent(inout) :: tmpfile
character(3), intent(in) :: operation
integer, intent(out), optional :: info
character(80), intent(in), optional :: filetag
logical, intent(out), optional :: inlist
!------------
!Data types to be saved:
character(10), dimension(max_tmp_files), save :: tmp_file_list = ''
character(80), dimension(max_tmp_files), save :: tmp_file_tags = ''
integer, save :: num_files = 0
!End save.
character(len=*), parameter :: this_routine = 'name_tmp_file'
character(7), parameter :: tmp_name7 = 'TMP____'
character(10) :: filename
character(3) :: name3
character(80) :: my_tag
integer, parameter :: max_count = 999 !maximum count possible due to subroutine
                                      !number2ascii
integer :: indx, at_indx, my_info
logical :: got_file, it_exists, no_stops
!------------
if (present(info)) then
  no_stops = .true.
else
  no_stops = .false.
endif
if (present(filetag)) then
  my_tag = filetag
else
  my_tag = ''
endif
!
my_info = 0
got_file = .false.
!
select case(operation)
case('get','GET')
  select case(tmpfile)
  case('',tmp_file_name)
    !Find a temporary file 
    loop1: do indx = 0, max_count
      call number2ascii(indx,name3)
      filename = tmp_name7//name3
      call check_if_exists(filename,it_exists,at_indx)
      if (.not.it_exists) then
        tmpfile = filename
        got_file = .true.
        exit loop1
      endif
    enddo loop1
  case default
    !Find a file with name similar to the one supplied
    call check_if_exists(tmpfile,it_exists,at_indx)
    if (it_exists) then
      !Find a alternative file name
      loop2: do indx = 0, max_count
        call number2ascii(indx,name3)
        filename = tmp_name7//name3
        call check_if_exists(filename,it_exists,at_indx)
        if (.not.it_exists) then
          tmpfile = filename
          got_file = .true.
          exit loop2
        endif
      enddo loop2
    else
      !This one (tmpfile) is OK
      got_file = .true.
    endif
  end select
  !
  if (got_file) then
    call add_file_to_list(tmpfile,my_info,my_tag)
    num_files = num_files + 1
  endif
  !
case('add','ADD')
  select case(tmpfile)
  case('',tmp_file_name)
    !These names are not allowed. Flag an error.
    print *,'INTERNAL ERROR in ',trim(this_routine)
    print *,'Illegal choice of name chosen is :',trim(tmpfile)
    if (no_stops) then
      my_info = -1
    else
      call internal_error(this_routine,__LINE__,'Illegal file name chosen')
    endif
  case default
    !Check to see if exactly the supplied name is available
    call check_if_exists(tmpfile,it_exists,at_indx)
    if (it_exists) then
      print *,'INTERNAL ERROR: File requested is in use. : ',trim(tmpfile)
      if (no_stops) then
        my_info = -1
      else
        call internal_error(this_routine,__LINE__,'File in use')
      endif
    else
      !This one (tmpfile) is OK
      got_file = .true.
    endif
  end select
  !
  if (got_file) then
    call add_file_to_list(tmpfile,my_info,my_tag)
    num_files = num_files + 1
  endif
  !
case('rel','REL')
  call remove_file_from_list(tmpfile,my_info)
  num_files = num_files - 1
case('inq','INQ')
  call check_if_exists(tmpfile,it_exists,at_indx)
  if (present(inlist)) then
    inlist = it_exists
  else
    print *,'INTERNAL ERROR: Variable INLIST not passed in call'
    if (no_stops) then
      my_info = -1
    else
      call internal_error(this_routine,__LINE__,'Variable missing in call')
    endif
  endif
case('out','OUT')
  call write_file_list
case default
  if (no_stops) then
    my_info = -1
  else
    print *,'INTERNAL ERROR: Illegal value of operation =',operation
    print *,'Temp file name passed is ',tmpfile
    call internal_error(this_routine,__LINE__,'Illegal value of operation')
  endif
end select
!
if (present(info)) info = my_info
!
return
contains
  subroutine check_if_exists(filename,it_exists,at_indx)
  implicit none
  character(*), intent(in) :: filename
  logical, intent(out) :: it_exists
  integer, intent(out) :: at_indx
  !------------
  integer :: indx
  !------------
  it_exists = .false.
  at_indx = -1
  loop1: do indx = 1, max_tmp_files
    if (filename.eq.tmp_file_list(indx)) then
      it_exists = .true.
      at_indx = indx
      exit loop1
    endif
  enddo loop1
  !
  return
  end subroutine check_if_exists

  subroutine add_file_to_list(filename,info,filetag)
  implicit none
  character(*), intent(in) :: filename
  integer, intent(inout) :: info
  character(80), intent(in) :: filetag
  !------------
  integer :: i
  logical :: added_file
  !------------
  added_file = .false.
  loop: do i = 1, max_tmp_files
    if (tmp_file_list(i).eq.'') then
      tmp_file_list(i) = filename
      tmp_file_tags(i) = filetag
      added_file = .true.
      exit loop
    endif
  enddo loop
  !
  if (.not.added_file) then
    print *,'ERROR: Number of temporary files exceeds limit set by'
    print *,'parameter MAX_TMP_FILES = ',max_tmp_files
    print *,'Re-compile code with larger value'
    info = -1
  endif
  !
  return
  end subroutine add_file_to_list

  subroutine remove_file_from_list(filename,info)
  implicit none
  character(*), intent(in) :: filename
  integer, intent(inout) :: info
  !------------
  integer :: at_indx
  logical :: it_exists
  !character(40) :: sys_command
  !integer :: sys_result
  !------------
  call check_if_exists(filename,it_exists,at_indx)
  if (it_exists) then
    tmp_file_list(at_indx) = ''
    tmp_file_tags(at_indx) = ''
    !Do not use these till record_handler is re-written. It doesn't like to find
    !file missing at present.
    ! sys_command = 'rm -f '//trim(filename)
    ! call system(sys_command,sys_result)
  else
    print *,'WARNING: Temporary file to be removed does not exist in the list'
    print *,'File name is ',trim(filename)
  endif
  info = -1
  !
  return
  end subroutine remove_file_from_list

  subroutine write_file_list
  implicit none
  integer :: i, indx
  !------------
  print *,''
  print *,'List of temporary file names and tags: '
  print '(1x,''Number of files currently in the list = '',i4)',num_files
  print *,'======================================'
  indx = 0
  do i = 1, max_tmp_files
    if (tmp_file_list(i).ne.'') then
      indx = indx + 1
      print '(1x,i4,1x,a10,4x,a60)',indx,trim(tmp_file_list(i)),&
        & trim(tmp_file_tags(i))
    endif
  enddo
  print *,'--------------------------------------'
  print *,''
  !
  return
  end subroutine write_file_list

end subroutine name_tmp_file

!----------------------------------------------------------------------

subroutine check_allocate(info,object,from_routine,line_no,do_i_stop)
implicit none
integer, intent(in) :: info
character(len=*), intent(in) :: from_routine
character(len=*), intent(in) :: object
integer, intent(in), optional :: line_no
logical, intent(in), optional :: do_i_stop
!------------
logical :: istop
!------------
if (present(do_i_stop)) then
  istop = do_i_stop
else
  istop = .true.
endif
!
if (info.ne.0) then
  write(6,*)'Allocation error for object ',trim(object)
  write(6,*)'Allocation attempted by routine ',trim(from_routine)
  if (present(line_no)) then
    write(6,*)'at line number ',line_no
  endif
  if (istop) then
    call internal_error('check_allocate',__LINE__,'ALLOCATE error')
  endif
endif
return
end subroutine check_allocate

subroutine check_deallocate(info,from_routine,object,line_no,do_i_stop)
implicit none
integer, intent(in) :: info
character(len=*), intent(in) :: from_routine
character(len=*), intent(in) :: object
integer, intent(in), optional :: line_no
logical, intent(in), optional :: do_i_stop
!------------
logical :: istop
!------------
if (present(do_i_stop)) then
  istop = do_i_stop
else
  istop = .true.
endif
!
if (info.ne.0) then
  write(6,*)'De-allocation error for object ',trim(object)
  write(6,*)'De-allocation attempted by routine ',trim(from_routine)
  if (present(line_no)) then
    write(6,*)'at line number ',line_no
  endif
  if (istop) then
    call internal_error('check_allocate',__LINE__,'DE-ALLOCATE error')
  endif
endif
return
end subroutine check_deallocate

!--------------------------------------------------------------------------

subroutine close_file_units
!------------------------
!Closes open file units (except special units)
implicit none
integer, parameter :: maxunit=100
integer, parameter :: isp=10 !Max number of special restricted unit nos.
integer :: ispecial(isp)
integer :: i,j
logical :: isitopen
data ispecial/5,6,8*0/
!------------
unit: do i = 1, maxunit
 do j = 1, isp
    if (i.eq.ispecial(j)) cycle unit
 enddo
 inquire(unit=i,opened=isitopen)
 if (isitopen) close(i)
enddo unit
!
return
end subroutine close_file_units

subroutine getunit(iunit)
!------------------------
!Returns the first available unused unit number.
!A very handy routine for getting unit numbers for temporary
!usage. Avoids conflicts with main unit numbers used in 
!the calling code.
!NOTE: Successive calls to this routine must be alternated with
!      OPEN commands. Else the same unit number will be returned 
!      each time!!!
!
!NOTE: Put special unit numbers in ispecial() below. These will not
!      be used.
!
!OUTPUT: 
!-------
!IUNIT = INTEGER : first available free unit number.
!
use parameters, only : maximum_unit_numbers
implicit none
integer, intent(out) :: iunit
!------------
character(len=*), parameter :: this_routine='getunit'
integer, parameter :: isp=10 !Max number of special restricted unit nos.
integer :: ispecial(isp)
integer :: i,j
logical :: isitopen
data ispecial/5,6,8*0/
!------------
!
iunit = -1
!
unit: do i = 1, maximum_unit_numbers
 special: do j = 1, isp
  !check to see if this is a special unit number
  if (i.eq.ispecial(j)) cycle unit
 enddo special
 inquire(unit=i,opened=isitopen)
 if (.not.isitopen) then
  !this is the unit number we need
  iunit=i 
  exit unit
 endif
enddo unit
!
select case(iunit)
 case (-1)
   !if you come here then no unit is available (amazing!)
   print '(1x,a,i5,a)','Subroutine GETUNIT: No unit number <=',&
      & maximum_unit_numbers,' available'
   print *,'Why the heck do you have so many files open?'
   print *,'Stopping: change value of MAXIMUM_UNIT_NUMBERS in '
   print *,'module parameters and recompile'
   call internal_error(this_routine,__LINE__,'Too many open files')
 case default
end select
!
return
end subroutine getunit

!----------------------------------------------------------------------

  subroutine my_timer(operation,title) !,debug)
  !Based on Wojtek's timing routine.
  !operation : One of 
  !           'enter' : entering routine/function
  !           'exit'  : exiting
  !           'report': write it all out nicely
  !title     : 20 character title of the routine/function you are timing.
  !debug     : (OPTIONAL) LOGICAL
  !            If present and TRUE, a single line stating the name of the
  !            subroutine entered is printed.
  !
!  use global_data, only : g_timer_debug  !AJWT
  use precision
  implicit none
  save
  character(len=*), intent(in) :: operation
  character(len=*), intent(in) :: title
!  logical, intent(in), optional :: debug
  !------------
  character(20) :: my_operation, my_title
  logical :: iam_init = .FALSE.  !Initial value. Lost after updating
  integer :: current_len = 0     !---do---
  integer :: posn
  integer, parameter :: maxlen = 100
  character(20), dimension(maxlen) :: titles
  integer, dimension(maxlen) :: num_calls
  real(dp), dimension(maxlen) :: entry_time, accumulated_time
  real(dp) :: current_time
  character(len=*), parameter :: fmt1='(1x,''Subroutine'',10x,&
  & ''Number of Calls'',3x,''Time (seconds)'')'
  character(len=*), parameter :: fmt2='(1x,A20,2x,i8,10x,f10.2)'
  !------------
  my_operation = trim(operation)
  my_title = trim(title)
!  if (present(debug)) then      !AJWT
!    if (debug.and.g_timer_debug) print *,'Subroutine ',trim(my_title),&
!      & '  ',trim(my_operation)
!  endif
  !
  if (.not.iam_init) then
    titles = ''
    posn = 0  !Current position in arrays
    !Now update the value of IAM_INIT. Contrary to expectation, the initial
    !value defined in the declaration statement above, will be over-written.
    iam_init = .true.
  endif
  !
  select case(my_operation)
    case('enter')
      posn = find_title(my_title) 
      if (posn.eq.-1) then
        !this is a new entry in the table
        current_len = current_len + 1
        if (current_len.GT.maxlen) then
          call error(1)
          !Set current length back one place
          current_len = current_len - 1
          return
        endif
        posn = current_len
        titles(posn) = my_title
        num_calls(posn) = 0
        accumulated_time(posn) = 0.0_dp
      endif
      num_calls(posn) = num_calls(posn) + 1
      entry_time(posn) = my_cpu_time()
      !
    case('exit')
      posn = find_title(my_title)
      if (posn.eq.-1) then
        call error(2)
        return
      endif
      current_time = my_cpu_time()
      accumulated_time(posn) = accumulated_time(posn) + &
                              & ( current_time - entry_time(posn) )
      !
    case('report')
      write(6,*)
      write(6,*)' Timing Report '
      write(6,*)'==============='
      write(6,fmt1)
      do posn = 1, current_len
         write(6,fmt2)titles(posn),num_calls(posn),accumulated_time(posn)
      enddo
      write(6,*)'===================================='
      write(6,*)
      !
    case default
      write(6,*)'MY_TIMER: Wrong operation command given'
      write(6,*)'Received: operation=',operation,' and title =',title
  end select
  !
  return
  contains
    function find_title(my_title)
    !Looks for TITLE in the list of TITLES(). If present, the position of TITLE
    !in the list is returned, else -1 is returned.
    implicit none
    character(len=*), intent(in) :: my_title
    !------------
    integer :: find_title
    integer :: found, pos
    !------------
    found = -1
    loop: do pos = 1, current_len
      if (my_title.eq.titles(pos)) then
        found = pos
        exit loop
      endif
    enddo loop
    !
    find_title = found
    !
    return
    end function find_title

    subroutine error(indx)
    implicit none
    integer indx
    !------------
    if (indx.eq.1) then
      write(6,*)'WARNING: In Timing routine MY_TIMER: MAXLEN exceeded.'
      write(6,*)'Skipping timing routine ',my_title
    elseif (indx.eq.2) THEN
      write(6,*)'WARNING: In Timing routine MY_TIMER: Cannot find TITLE'
      write(6,*)'Missing title is =',my_title
    else
      write(6,*)'WARNING: Wrong INDX value passed to SUBROUTINE ERROR'
    endif
    return
    end subroutine error

    function my_cpu_time()
    implicit none
    real(dp) :: my_cpu_time
    !------------
#if defined(G77)
    my_cpu_time = second()
#endif
#if defined(DECALPHA) || defined(SUN) || defined(SUNF90)
    real(4), dimension(2) :: tarray
    real(4) :: etime
    my_cpu_time = etime(tarray)
#endif
#if (SGI) || (RS6K) || (IBM32) || (IBM64)
    my_cpu_time = float(mclock())*1.0E-2
#endif
#if (IFC) || (IFORT) || (LAHEY) || (PGF90) || (G95) || (GFORTRAN)
    real(dp) :: cputime
    call cpu_time(cputime)
    my_cpu_time = cputime
#endif
    !
    return
    end function my_cpu_time

  end subroutine my_timer

  subroutine warning(routine_name,line_no,warn_str)
  implicit none
  character(len=*), intent(in) :: routine_name
  integer, intent(in) :: line_no
  character(len=*), intent(in) :: warn_str
  !------------
  print *,'WARNING: in subroutine ',trim(routine_name),' line ',line_no
  print *,trim(warn_str)
  call flush(6)
  stop
  end subroutine warning

  subroutine internal_error(routine_name,line_no,error_str)
  implicit none
  character(len=*), intent(in) :: routine_name
  integer, intent(in) :: line_no
  character(len=*), intent(in) :: error_str
  !------------
  character(10) :: tmpfile = ''
  !------------
  print *,'INTERNAL ERROR: in subroutine ',trim(routine_name)
  print *,'at line number ',line_no
  print *,trim(error_str)
  print *,'Contact program developers with this problem.'
  print *,''
  print *,'Closing open files'
  call close_file_units
  call name_tmp_file(tmpfile,'out')
  print *,'Stopping'
  call flush(6)
  stop
  end subroutine internal_error

  subroutine general_error(routine_name,line_no)
  implicit none
  character(len=*), intent(in) :: routine_name
  integer, intent(in) :: line_no
  !------------
  character(10) :: tmpfile = ''
  !------------
  print *,'ERROR: in subroutine ',trim(routine_name)
  print *,'at line number ',line_no
  print *,'This error is likely due to problems with your input file'
  print *,'or data files (hessians, MO files, etc.) or, if you have'
  print *,'restarted from an earlier calculation, there could be a possible'
  print *,'inconsistency. '
  print *,''
  print *,'Closing open files'
  call close_file_units
  call name_tmp_file(tmpfile,'out')
  print *,'Stopping'
  call flush(6)
  stop
  end subroutine general_error

end module common_routines

subroutine number2ascii(nr,str)
!Written by Wojtek.
!Takes in a three digit (max) integer and returns a three
!character string of that integer. Other languages would have a
!command like CHAR(). Apparently the Lahey F95 compiler (from
!Fijitsu) has this command (called ACHAR), but G77 doesn't hence
!this bit of code...
use common_routines
use parameters, only : max_tmp_files
implicit none
integer, intent(in) :: nr
character(3), intent(out) :: str
!------------
character(len=*), parameter :: this_routine = 'number2ascii'
integer :: i, j, k
character(1), dimension(0:9) :: digits
data digits/'0','1','2','3','4','5','6','7','8','9'/
!------------
if (nr.ge.max_tmp_files) then
   write (6,*) 'ERROR: Number too big for number2ascii'
   call internal_error(this_routine,__LINE__,'Probably too many files open')
endif
i=nr/100
j=(nr-100*i)/10
k=nr-100*i-10*j
!
str(1:1)=digits(i)
str(2:2)=digits(j)
str(3:3)=digits(k)
! 
return
end subroutine number2ascii
