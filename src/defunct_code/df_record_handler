!RECORD HANDLER
!==============
!A. J. Misquitta
!3rd June 2004
!University of Cambridge
!-----------------------
!Routines to handle writing and reading of direct access un-formatted file 
!that are arbritrarily large
!USAGE:
!------
!Inlcude the module with the statement:
!
! USE record_handler, &
! & ONLY : init_record_handler, write_record, read_record, &
! &        leave_record_handler, query_record_handler
!
!Then call routines as described below:
!--------------------------------------
!NOTE: All calls take the optional last argument 
!      printinfo=.true. or .false.
!This controls printing of comments.
!
!To query files (no printinfo option):
!  call query_record_handler(filename,filestatus,len_first_rec,num_records)
!
!performs some queries on the file filename.
!INPUT:
!     filename = character(10)
!OUTPUT: 
!     (All optional)
!     filestatus = 'NEW' or 'OLD'
!     len_first_rec = length of the first record
!     num_records = number of records written
!On exit the file is closed.
!
!To write a file
!===============
!  (1) call init_record_handler(filename,filestatus,info,num_records)
!      where filename = character(10) string
!            filestatus = NEW/new or OVR/ovr
!                        NEW: specifies if the file is to be created (i.e., it
!                        is new)
!                        OVR: specifies that the file exists and is to be
!                        overwritten.
!            info     = (OUTPUT) 0=normal exit, <0=error
!            num_records = (OUTPUT, OPTIONAL) set to the number of records.
!                     (use query_record_handler for this instead)
!      This step is optional as, if found necessary, this call can also be made
!      by the reading and writing routines (where the filestatus is assumed to
!      be 'NEW' in the write routine, and 'OLD' in the read routine).
!
!  (2) call write_record(filename,X,indicesX,labelX,info)
!      where 
!         filename = character(10) string
!         indexX = integer: the index of the data set to be written. This index
!                  will be used to identify the data set while reading.
!         X = real(dp) array : The data set values
!         indicesX = integer array : The indices of the data set values
!         labelX = character(8) label to facilitate identification of the data
!                  set. Most likely this won't be used very much. But I put it
!                  here for compatibility with the SAPT integrals.
!                  NOTE: If set to '        ', then no checks are made. Useful
!                  when the label is unknown.
!         info     = (OUTPUT) 0=normal exit, <0=error
!      Notice how there is no mention of packing/unpacking indices. This is left
!      to the user. All this routine does is keep track of records.
!      Call this routine as many times as is needed.
!
!  (3) call leave_record_handler(filename,info)
!      where
!           filename = character(10) filename. 
!           info     = (OUTPUT) 0=normal exit, <0=error
!      This MUST be called or else the record_handler doesn't know that you've
!      completed writing into the file. 
!
!Once a file has been written, another can be written using the same set of
!routines. However, multiple files cannot be written simultaneously.
!           
!Reading data from a file
!========================
!  (1) call init_record_handler(filename,filestatus,info,num_records)
!      where filename = character(10) string
!            filestatus = OLD/old
!                        This specified that the file exists and it is to 
!                        be read from disk (i.e., it is old)
!                        ADD/add: add a file created by record_handler in an 
!                        earlier run. 
!            info     = (OUTPUT) 0=normal exit, <0=error
!            num_records = number of records in file if file is OLD.
!                        (optional argument)
!      This step is optional as, if found necessary, this call can also be made
!      by the reading and writing routines (where the filestatus is assumed to
!      be 'NEW' in the write routine, and 'OLD' in the read routine).
!
!  (2) read_record(filename,indexX,labelX,X,info)
!      where
!           filename = character(10) string
!           indexX = integer: The index of the data set that you wish to be read
!                    in (it is used to figure out the record on which the data
!                    set resides).
!           labelX = character(8) label of the data set required. 
!           X = real(dp) array into which the data set will be put.
!           info     = (OUTPUT) 0=normal exit, <0=error
!     Once again, no packing/unpacking routines are employed. The indices saved
!     with the data set are used to fill the array X() provided.
!     This routine can be called as many times as is needed.
!
!  (3) call leave_record_handler(filename,info)
!      where
!           filename = character(10) filename. 
!           info     = (OUTPUT) 0=normal exit, <0=error
!      This MUST be called or else the record_handler doesn't know that you've
!      completed reading the file and that it should release accumulated
!      information about the file.
!
!***Multiple files can now be read/written simultaneously. 
!The number of files that can be handled by this module is determined by the
!parameter maxfiles (see MODULE record_handler_arrays). 
!
!One drawback: At the moment there is no way of removing a file from the list of
!files. So, during the course of the calling program, as the number of files
!handled by record_handler increases, the list keeps increasing (upto maxfiles).
!There is no provision for removing a file from the list.
!
!Caveat (1): Don't try to write and read the same file without calling
!leave_record_handler between the writes and reads. This is because the
!Table-of-Contents are written to the file only intermittently. So reading is
!not possible untill the file is closed and the last TOC is written.
!
!Caveat (2):While you can read two files like so:
!           call read_record(file1,...)
!           call read_record(file2,...)
!           call read_record(file1,...)
!           call read_record(file2,...)
!           etc.
!this is wasteful as the TOC is erased and re-read at each read. Also, the
!output of the code is quite a bit. 
!
!To Be Done: Some means of adding records to a file that already exists.
!
!A brief programming detail: The capacity of the code to handle multiple files
!was carried out quite simply: 
!   (1) Identify the variables that determine the status of a file. 
!   (2) Make arrays of these variables - each entry in these arrays determines
!   the variables for a given file.
!   (3) On entry (at any of the four entry points), find out if the file is on
!   list. If so, identify its position in the list. If not, add it to the list
!   if there is space.
!   (4) The using pointers, set the values of the variables identified in (1)
!   using the arrays in (2). 
!   (5) Now use the code as before (i.e., as it was in the single-file
!   implementation).
!   (6) The full TOC is not replicated for each file, so this needs to be
!   re-created for each file read.
!
module record_handler_arrays
!The arrays here store a snapshot of the relevant variables keeping track of the
!status of a file. The variable names have the prefix `ar_' added.
!
use precision
use common_routines
SAVE
integer, parameter :: maxfiles = 1000
character(10), dimension(maxfiles), target :: ar_filename = ''
integer, dimension(maxfiles), target :: ar_fileunit = 0
integer, dimension(maxfiles), target :: ar_curr_rec = 0
integer, dimension(maxfiles), target :: ar_next_rec = 0
!Logical flags for initialization etc. (also set their initial values)
logical, dimension(maxfiles), target :: ar_iam_init_toc = .false.
logical, dimension(maxfiles), target :: ar_iam_init_record_handler = .false.
logical, dimension(maxfiles), target :: ar_iam_init_bookkeeper = .false.
logical, dimension(maxfiles), target :: ar_file_is_new = .false.
!Parameters governing all things related to record size.
!reclen determines maxlen and len_one_toc.
! 16 + 12*maxlen <= reclen (solve for largest integer)
! 4*len_one_toc = reclen
integer, parameter :: reclen = 32760     !record length
integer, parameter :: maxlen = 2728      !max len of arrays x() n() in a rec
integer, parameter :: len_one_toc = 8190 !length of each toc (see below)
!These arrays are used by bookkeeper
!-----------------------------------
!curr_toc : contains the current table of contents for a file
!rec_curr_toc: is the record at which the curr_toc is to be written
!num_in_curr_toc : is the number of entries in the curr_toc
integer, dimension(len_one_toc,maxfiles), target :: ar_curr_toc
integer, dimension(maxfiles), target :: ar_rec_curr_toc
integer, dimension(maxfiles), target :: ar_num_in_curr_toc
end module record_handler_arrays

module record_handler_pointers
use precision
SAVE
integer :: file_indx         !index of the current file
integer :: numfiles = 0      !number of files (should be  <= maxfiles) 
integer, pointer :: fileunit
!Logical flags for initialization etc.
logical, pointer :: iam_init_toc
logical, pointer :: iam_init_record_handler
logical, pointer :: iam_init_bookkeeper
logical, pointer :: file_is_new
end module record_handler_pointers

module bookkeeper_pointers
!These arrays are used by bookkeeper
SAVE
integer, pointer :: curr_rec, next_rec
integer, dimension(:), pointer :: curr_toc
integer, pointer :: rec_curr_toc
integer, pointer :: num_in_curr_toc
end module bookkeeper_pointers

!--------------------
module record_handler
!--------------------
!A set of routines to read and write large direct-access files.
!Records are assumed to be in a similar format as those used in the SAPT 
!suite of codes. A linked-list like format is used. 
!Each record is of the form:
!
! next_rec len_rec label x(len_rec) n(len_rec)
!    4        4      8     len_rec*(8 + 4)         + 8 (overhead)
! 
! This means 2728 elements of x() and n() can be written into one record. 
!
!where
!   next_rec = next record to be read in. If zero, then no more records left.
!   len_rec  = number of elements to be read from present record.
!   label    = character(8) label of current record
!   x(len_rec) = real(dp) array 
!   n(len_rec) = integer(4) indices of elements x(). This is needed as x(i) 
!                could be the j^th element in an array A(). n(i) tells you the
!                position of x(i) in array A(). This is useful to have as often
!                not all elements of A() are stored. Only non-zero values. So we
!                need some way of knowing which elements have been skipped and
!                what the next non-zero element is.
!
!Record length: 32KB - 8 bytes= 32760 bytes 
!               (NOTE: each record has 8 bytes of overhead storage required)
!
!The Table of Contents will be put starting at record 1. As many records as are
!needed will be reserved for the table of contents. A linked-list is used to
!expand the storage for the table of contents. See below for details.
!                
use precision
use record_handler_arrays
use record_handler_pointers
use common_routines
IMPLICIT NONE
PRIVATE
SAVE
public :: init_record_handler, write_record, read_record, leave_record_handler
public :: query_record_handler
!Pertaining to the Full Table of Contents (toc)
integer(4), dimension(:), allocatable :: toc !The full Table of Contents 
integer :: len_full_toc                   !length of the full TOC
!Pertaining to each record (maxlen is set in record_handler_arrays)
real(dp), dimension(maxlen) :: values    !values stored in record
integer(4), dimension(maxlen) :: indices !indices of those values
!Misc
logical :: l_print  !controls printing of info (not errors)
character(len=*), parameter :: this_routine='record_handler'
contains
  subroutine init_record_handler(filename,filestatus,info,&
                                & num_records,printinfo)
  !Read the Table of Contents of a direct access file
  !INPUT: 
  !    filename = character(10) file name
  !    filestatus = character(3): OLD/old or NEW/new or OVR/ovr or ADD/add
  !                 NEW: This is a new file to be opened
  !                 OLD: This is an old file already opened and closed in the
  !                      current run.
  !                 OVR: This file has been written and closed, overwrite it.
  !                 ADD: This file has been created in an earlier run. Add it to
  !                      record_handler.
  !    info      =  0 : normal exit
  !               < 0 : Error
  !    num_records = integer (optional) is set equal to the number of records
  !    printinfo = LOGICAL (optional) controls printing.
  !
  !filestatus tells the record_handler whether the file has already been created
  !(old or OLD) or whether it is to be written into for the first time (new or
  !NEW). In the latter case, any existing file by the same name is over-written.
  !
  !Table of Contents (TOC):
  !------------------------
  !Unlike the case with SAPT where only a small (less than 90) number of types
  !of integrals needed to be written to disk, the routines in this module are
  !designed to handle a large (really large!) number of types of integrals. Each
  !type needs its own entry in the TOC. This means that the TOC can be quite
  !large and may well span many records itself. This, of course, means that we
  !must have some mechanism for expanding the TOC to arbritrary length. This
  !will be done by using the linked-list idea which is described next.
  !  Record no.     Contents
  !     1           TOC(1) 
  !     2           integrals
  !     ...         ...
  !    8189         integrals
  !    8190         TOC(2)
  !    8191         integrals
  !     ...         ...
  !    etc
  !
  !There are reclen = 32760 bytes in a record. Thus, each TOC can store 8190
  !integer(4) numbers. 8189 of these will point to records containing integrals.
  !The last, the 8190^th will contain the record of the next TOC. In this way
  !the TOC can be expanded arbritrarliy. The user need not worry about this. 
  !
  implicit none
  character(10), intent(in) :: filename
  character(3), intent(in) :: filestatus
  integer, intent(out) :: info
  integer, optional, intent(out) :: num_records
  logical, optional, intent(in) :: printinfo
  !------------
  character(3) :: actual_file_status
  logical :: openit, readit, iamopen, iexist, error
  integer :: stat
  logical :: overwrite_file, warn, initialize, add_file
  integer(4) :: reclenrat !used to fix differences in record lengths between compilers
  inquire(iolength=reclenrat) reclenrat
!  AJWT 
!  some compilers give reclenrat as 4, and others 1
!  We want to count in bytes not words, so we change reclenrat to 4/reclenrat
!  later we shall use reclen/reclenrat as our length scale
  reclenrat=4/reclenrat
  !------------
  !
  info = 0
  !
  if (present(printinfo)) then
    l_print = printinfo
  else
    l_print = .false.
  endif
  !Set up the record_handler variables for this file
  call find_file(filename,file_indx,info,actual_file_status)
  if (info.lt.0) return
  call set_pointers(file_indx)
  !
  openit = .false.
  readit = .false.
  error = .false.
  warn = .false.
  add_file = .false.
  initialize=.false.
  select case(filestatus)
    case('new','NEW')
      openit = .true.
      file_is_new = .true.
      initialize=.true.
      if (actual_file_status.ne.'NEW') error = .true.
    case('old','OLD')
      readit = .true.
      file_is_new = .false.
      if (actual_file_status.ne.'OLD') error = .true.
    case('ovr','OVR')
      openit = .true.
      overwrite_file = .true.
      file_is_new = .true.
      initialize=.true.
      if (actual_file_status.ne.'OLD') warn = .true.
    case('add','ADD')
      openit = .true.
      add_file = .true.
    case default
      write(6,*)'ERROR: record_handler: Incorrect filestatus provided'
      write(6,*)'filestatus = ',filestatus
      print *,'Filename = ',trim(filename)
      info = -1
      return
  end select
  !
  if (error) then
    write(6,*)'init_record_handler: ERROR: Inconsistent file status'
    write(6,*)'File ',filename,' was expected to be ',filestatus
    write(6,*)'but was found to be ',actual_file_status
    info = -1
    return
  endif
  !
  if (warn.and.l_print) then
    write(6,*)'init_record_handler: WARNING: File to be overwritten seems'
    write(6,*)'to be new to record_handler. File is ',filename
  endif
  !
  if (add_file) then
    if (actual_file_status.eq.'OLD') then
      print *,'init_record_handler:ERROR: File to be added is already present'
      print *,'in record_handler. Possible name conflict?'
      print *,'File name is ',filename
      info = -1
      return
    else
      readit = .true.
      openit = .false.
    endif
  endif
  !
  if (initialize) then
    !Some initializations
    fileunit = -1 !nonsense values
    iam_init_toc = .false.
  endif
  !
  if (openit) then
    inquire(file=filename,opened=iamopen,exist=iexist,number=fileunit)
    if (overwrite_file) then
      if (iamopen) close(fileunit)
      call initialize_file_variables
    else
      if (iamopen.or.iexist) then
        write(6,*)'ERROR: init_record_handler: Expecting to open a new file'
        write(6,*)'but the file ',filename,' seems to exist. '
        info = -1
        return
      endif
    endif
    call getunit(fileunit)
    
    open (fileunit,file=filename,status='unknown',access='direct',&
         & recl=reclen/reclenrat,form='unformatted',action='readwrite',iostat=stat)
    call io_check(filename,stat,info,openingfile=.true.)
    if (info.lt.0) return
    !
    if (present(num_records)) then
      !The file is new so no records written. Set a nonsense value.
      num_records = -1
    endif
  endif
  !
  if (readit) then
    inquire(file=filename,opened=iamopen,exist=iexist)
    if (.not.iexist) then
      write(6,*)'ERROR: record_handler: Bad news. Your file ',filename
      write(6,*)'does not exist.'
      info = -1
      return
    endif
    if (.not.iamopen) then
      call getunit(fileunit)
      open (fileunit,file=filename,status='old',access='direct',&
           & recl=reclen/reclenrat,form='unformatted',action='read',iostat=stat)
      call io_check(filename,stat,info,openingfile=.true.)
      if (info.lt.0) return
    else
      !if the file is opened, determine the unit number it is opened on
      inquire(file=filename,number=fileunit)
    endif  
    !Now read the TOC
    call read_toc(filename,info)
    !
    if (present(num_records)) then
      num_records = len_full_toc
    endif
  endif
  !
  iam_init_record_handler = .true.
  if (l_print) write(6,*)'...initialization complete for file ',filename
  return
  end subroutine init_record_handler

  subroutine query_record_handler(filename,info,filestatus,&
                                 & len_first_rec,num_records,printinfo)
  !performs some queries on the file filename.
  !INPUT:
  !     filename = character(10)
  !     printinfo = T or F (optional)
  !OUTPUT: 
  !     info = 0  Normal exit
  !          < 0  Errors
  !     (All optional)
  !     filestatus = 'NEW' or 'OLD'
  !     len_first_rec = length of the first record after the TOC
  !     num_records = number of records written
  !On exit the file is closed.
  implicit none
  character(10), intent(in) :: filename
  integer, intent(out) :: info
  logical, optional, intent(in) :: printinfo
  character(3), optional, intent(out) :: filestatus
  integer, optional, intent(out) :: len_first_rec, num_records
  !------------
  character(3) :: file_status
  logical :: l_check_status, l_check_len_rec1, l_check_num_recs
  logical :: l_file_open, l_file_exists
  integer(4) :: stat, curr, next, length
  integer :: i
  integer(4) :: reclenrat !used to fix differences in record lengths between compilers
  inquire(iolength=reclenrat) reclenrat
!  AJWT 
!  some compilers give reclenrat as 4, and others 1
!  We want to count in bytes not words, so we change reclenrat to 4/reclenrat
!  later we shall use reclen/reclenrat as our length scale
  reclenrat=4/reclenrat
  !------------
  info = 0
  !
  if (present(printinfo)) then
    l_print = printinfo
  else
    l_print = .false.
  endif
  if (present(filestatus)) then
    l_check_status = .true.
  else
    l_check_status = .false.
  endif
  if (present(len_first_rec)) then
    l_check_len_rec1 = .true.
    l_check_status = .true.
  else
    l_check_len_rec1 = .false.
  endif
  if (present(num_records)) then
    l_check_num_recs = .true.
    l_check_status = .true.
  else
    l_check_num_recs = .false.
  endif
  !
  if (l_check_status) then
    call find_file(filename,file_indx,info,actual_file_status=file_status,&
                  & query=.true.)
    if (info.lt.0) return
  endif
  if (present(filestatus)) then
    filestatus = file_status
  endif
  !
  l_file_open = .false.
  l_file_exists = .false.
  !
  if (l_check_len_rec1.or.l_check_num_recs) then
    if (file_status.eq.'NEW') then
      write(6,*)'query_record_handler: Cannot check a new file :',filename
      if (l_check_len_rec1) len_first_rec = -1
      if (l_check_num_recs) num_records = -1
    else
      call set_pointers(file_indx)
      inquire(file=filename,opened=l_file_open,exist=l_file_exists)
      if (.not.l_file_exists) then
        write(6,*)'query_record_handler: ERROR nonexistent file: ',filename
        if (l_check_len_rec1) len_first_rec = -1
        if (l_check_num_recs) num_records = -1
        return
      endif
      if (.not.l_file_open) then
        call getunit(fileunit)
        open (fileunit,file=filename,status='old',access='direct',&
             &recl=reclen/reclenrat,form='unformatted',action='read',iostat=stat)
        call io_check(filename,stat,info,openingfile=.true.)
        if (info.lt.0) return
      else
        !if the file is opened, determine the unit number it is opened on
        inquire(file=filename,number=fileunit)
      endif  
      !
      call read_toc(filename,info)
      if (info.lt.0) return
      !
      if (l_check_num_recs) then
        !The number of records stored is the length of the TOC
        num_records = len_full_toc
      endif
      if (l_check_len_rec1) then
        !Read in the first bits of record 1. Note: The conceptual record 1
        !is stored in potentially many records of maximum length maxlen starting
        !from physical record 2. We have to read in all physical records that go
        !to form the conceptual record 1 before we have the actual length of
        !this record.
        if (len_full_toc.lt.1) then
          write(6,*)'query_record_handler: No records in file ',filename
          len_first_rec = -1
        else
          len_first_rec = 0
          curr = 2 !start from record 2
          rec_loop: do
          if (curr.eq.0) exit rec_loop
            read(fileunit,rec=curr,iostat=stat)next,length
            call io_check(filename,stat,info,record=curr)
            if (info.lt.0) return
            curr = next
            len_first_rec = len_first_rec + length
          enddo rec_loop
        endif
      endif
      !
    endif
  endif
  !
  if (.not.l_file_open) then
    !The file wasn't open, so close it before exiting
    call leave_record_handler(filename,info,printinfo=l_print)
    if (info.lt.0) then
      print *,'query_record_handler: Internal error'
      print *,'Error while closing file : ',trim(filename)
      print *,'Could be due to an earlier error.'
      info = -1
    endif
  endif
  !
  if (info.lt.0) then
    write(6,*)' Index       filename'
    write(6,*)'---------------------'
    do i = 1, numfiles
      write(6,*)i,ar_filename(i) 
    enddo
    write(6,*)'====================='
  endif
  !
  return
  end subroutine query_record_handler

  subroutine read_toc(filename,info)
  !read in the entire linked-list TOC into memory. Since the total length of the
  !TOC is not known (unless you read it), we first perform a dummy read to
  !determine the number of elements in the TOC, and only then allocate space for
  !the TOC and read it all in. 
  implicit none
  character(10), intent(in) :: filename
  integer, intent(out) :: info
  !------------
  integer(4), dimension(len_one_toc) :: one_toc
  integer :: i
  integer :: stat, record
  integer :: entries_in_toc
  logical :: last
  integer :: ierr
  !------------
  integer, save :: toc_of_file_indx
  !-----------
  if (allocated(toc)) then
    if (toc_of_file_indx.eq.file_indx) then
      !then the toc in memory belongs to the file we want to read
      return
    else
      !the toc in memory belongs to another file...
      deallocate(toc,stat=ierr)
      call check_deallocate(ierr,'toc',this_routine,__LINE__,.false.)
      if (ierr.ne.0) then; info = -1; return; endif
      !...and set the toc initilization flag to false
      ar_iam_init_toc(toc_of_file_indx) = .false.
    endif
  endif
  !
  !set the index that states which file this is the toc of...
  toc_of_file_indx = file_indx
  !  
  last = .false.
  entries_in_toc = 0
  record = 1
  !
  toc_loop1 : do
  if (record.eq.0) exit toc_loop1
    !
    read(fileunit,rec=record,iostat=stat) one_toc
    call io_check(filename,stat,info,record=record)
    if (info.lt.0) return
    !read all but the last element (which is reserved for the record number of
    !the next toc - if present) and check to see where the first zero is
    !present. This signifies the end of the TOC.
    entries1: do i = 1, (len_one_toc - 1)
      if (one_toc(i).eq.0) then
        last = .true.
        exit entries1
      else
        entries_in_toc = entries_in_toc + 1
      endif
    enddo entries1
    !
    if (last) then
      !we have found the last element in the TOC
      len_full_toc = entries_in_toc
      exit toc_loop1
    else
      !there is more to read. The next toc is located on record number
      !one_toc(len_one_toc) (i.e., the last entry in one_toc).
      record = one_toc(len_one_toc)
    endif
    !
  enddo toc_loop1
  !
  !Now that we have the number of elements in the TOC, allocate it and read in
  !all the entries...
  allocate(toc(len_full_toc),stat=ierr)
  call check_allocate(ierr,'toc',this_routine,__LINE__,.false.)
  if (ierr.ne.0) then
    print *,'record_handler: Filename = ',trim(filename)
    info = -1
    return
  endif
  !
  last = .false.
  record = 1
  entries_in_toc = 0
  toc_loop2 : do
  if (record.eq.0) exit toc_loop2
    !
    read(fileunit,rec=record,iostat=stat) one_toc
    call io_check(filename,stat,info,record=record)
    if (info.lt.0) return
    entries2: do i = 1, (len_one_toc - 1)
      if (one_toc(i).eq.0) then
        last = .true.
        exit entries2
      else
        entries_in_toc = entries_in_toc + 1
        toc(entries_in_toc) = one_toc(i)
      endif
    enddo entries2
    !
    if (last) then 
      exit toc_loop2
    else
      !there is more to read. The next toc is located on record number
      !one_toc(len_one_toc) (i.e., the last entry in one_toc).
      record = one_toc(len_one_toc)
    endif
    !
  enddo toc_loop2
  !
  iam_init_toc = .true.
  !
  call check_toc
  !
  return
  contains
    subroutine check_toc
    !perform some simple checks on the toc
    implicit none
    if (l_print) write(6,*)'Length of the Table of Contents is ',len_full_toc
    if (len_full_toc.eq.0) then
      if (l_print) then
        write(6,*)'WARNING: Odd! There are no records to be read in this file.'
      endif
    else
      if (toc(1).ne.2) then
        write(6,*)'ERROR: This Table of Contents appears to be incorrect'
        write(6,*)'Was the file ',filename,' written using record_handler'
        write(6,*)'routine write_record ?'
        info = -1
        return
      endif
    endif
    return
    end subroutine check_toc

  end subroutine read_toc

  subroutine read_record(filename,indexX,labelX,X,info,printinfo)
  implicit none
  character(10), intent(in) :: filename       !file to be read from
  integer, intent(in) :: indexX               !position of record to be read in
  character(8), intent(in) :: labelX          !label of the record 
  real(dp), dimension(:), intent(out) :: X    !array to be filled up
  integer, intent(out) :: info
  logical, optional, intent(in) :: printinfo
  !------------
  character(3) :: filestatus
  integer :: stat
  integer(4) :: num_entries
  integer :: rec_num, i
  character(8) :: readlabel
  integer(4) :: next_rec
  !------------
  info = 0
  !
  if (present(printinfo)) then
    l_print = printinfo
  else
    l_print = .false.
  endif
  !set up the record-handler for this file
  call find_file(filename,file_indx,info)
  if (info.lt.0) return
  call set_pointers(file_indx)
  !
  if (.not.allocated(toc)) then
    iam_init_toc = .false.
  endif
  !
  if (.not.iam_init_record_handler) then
    if (l_print) then
      write(6,*)'WARNING: Record_handler not initialized for file ',filename
    endif
    filestatus = 'OLD'
    call init_record_handler(filename,filestatus,info)
    if (info.lt.0) then
      print *,'read_record: Internal error while trying to initialize'
      print *,'file ',trim(filename)
      info = -1
      return
    endif
  else
    call read_toc(filename,info)
    if (info.lt.0) then
     print *,'read_record: Internal error while reading TOC for ',trim(filename)
     info = -1
     return
    endif
  endif
  !
  X = 0.0_dp
  !first a check on indexX 
  if ((indexX.le.0).or.(indexX.gt.len_full_toc)) then
    write(6,*)'ERROR: record_handler: Illegal value of record position',indexX
    write(6,*)'Number of records in file ',filename,' is ',len_full_toc
    info = -1
    return
  endif
  !Determine the record number
  rec_num = toc(indexX)
  !One more check...just in case
  if (rec_num.eq.0) then
    write(6,*)'MODULE record_handler: Internal error. Got zero record number'
    write(6,*)'for record position ',indexX
    write(6,*)'This should never have happened if the file being read was'
    write(6,*)'written by record_handler. Was it?'
    print *,'Filename is ',trim(filename)
    info = -1
    return
  endif
  !
  rec_loop: do
  if (rec_num.eq.0) exit rec_loop
    !
    indices = 0
    values = 0.0_dp
    read(fileunit,rec=rec_num,iostat=stat)next_rec,num_entries,readlabel,&
    & values,indices
    !
    call io_check(filename,stat,info,record=rec_num)
    if (info.lt.0) return
    if (readlabel.ne.labelX) then
      if (labelX.eq.'        ') then
        !Blank label. Do nothing.
      else
        write(6,*)'ERROR record_handler:The record read in had the wrong label'
        write(6,*)'Label wanted is ',labelX
        write(6,*)'Label of record read in is ',readlabel
        print *,'Filename is ',trim(filename)
        info = -1
        return
      endif
    endif
    !
    !fill up the X() array
    do i = 1, num_entries
      X(indices(i)) = values(i)
    enddo
    !
    !update next record to be read...
    rec_num = next_rec
    !
  enddo rec_loop
  !
  return
  end subroutine read_record

  subroutine write_record(filename,X,indicesX,labelX,info,printinfo)
  ! filename = character(10) name of file to write to
  ! indexX   = index of array X in the TOC
  ! X        = read(dp) array to be written
  ! indicesX = integer array of indices of X
  ! info     = 0  Normal exit
  !          < 0  Errors
  ! 
  implicit none
  character(10), intent(in) :: filename
  real(dp), dimension(:), intent(in) :: X
  integer(4), dimension(:), intent(in) :: indicesX
  character(8), intent(in) :: labelX
  integer, intent(out) :: info
  logical, optional, intent(in) :: printinfo
  !------------
  character(3) :: filestatus, newold
  integer :: lenX, len_indxX
  integer :: elements_written
  logical :: last
  integer :: start_indx, end_indx, num_to_be_written
  integer :: curr_rec, next_rec
  !------------
  !
  info = 0
  !
  if (present(printinfo)) then
    l_print = printinfo
  else
    l_print = .false.
  endif
  !set up the record-handler for this file
  call find_file(filename,file_indx,info)
  if (info.lt.0) return
  call set_pointers(file_indx)
  !
  if (.not.iam_init_record_handler) then
    if (l_print) then
      write(6,*)'WARNING: Record_handler not initialized for file ',filename
    endif
    filestatus = 'NEW'
    call init_record_handler(filename,filestatus,info)
    if (info.lt.0) then
      print *,'write_record: Internal error while trying to initialize'
      print *,'file ',trim(filename)
      info = -1
      return
    endif
  endif
  !
  lenX = SIZE(X)
  len_indxX = SIZE(indicesX)
  if (lenX.ne.len_indxX) then
    write(6,*)'ERROR: record_handler: arrays X and indicesX have different'
    write(6,*)'lengths. Lengths are ',lenX, len_indxX
    print *,'Filename is ',trim(filename)
    info = -1
    return
  endif
  !
  newold = 'NEW'
  elements_written = 0
  start_indx = 1    !starting index of X and indicesX to be written
  end_indx = maxlen !ending indx of X and indicesX to be written
  last = .false.    !flag signaling last record of the data set
  Xloop : do 
  if (last) exit Xloop
    !Initialize arrays that will be used to write info into records
    values = 0.0_dp
    indices = 0
    !
    num_to_be_written = maxlen
    if (end_indx.gt.lenX) then
      !This is the last bit to be written
      end_indx = lenX
      num_to_be_written = lenX - start_indx + 1
      last = .true.
    endif
    elements_written = elements_written + num_to_be_written
    !
    values(1:num_to_be_written) = X(start_indx:end_indx)
    indices(1:num_to_be_written) = indicesX(start_indx:end_indx)
    !
    call bookkeeper(curr_rec,next_rec,newold,info)
    if (info.lt.0) then
      print *,'write_record: Internal error in bookkeeper'
      print *,'File is ',trim(filename)
      info = -1
      return
    endif
    newold = 'OLD'
    !
    if (last) then
      !If this is the last record, then there are no more records to be read.
      !Signal this by setting next_rec = 0
      next_rec = 0
    endif
    !
    write(fileunit,rec=curr_rec,ERR=1)&
    & next_rec, num_to_be_written, labelX, values, indices
    !
    !determine the new start and end index of X and indicesX to be written
    start_indx = end_indx + 1
    end_indx = start_indx + maxlen - 1
    !
  enddo Xloop

  call flush(fileunit)

  return
  !
  !Errors come here
1 write(6,*)'MODULE record_handler: A write-error has occured (not EOF or EOR)'
  write(6,*)'file name =',filename,' record ',curr_rec
  info = -1
  return
  !
  end subroutine write_record

  subroutine leave_record_handler(filename,info,printinfo)
  implicit none
  character(10), intent(in) :: filename
  integer, intent(out) :: info
  logical, optional, intent(in) :: printinfo
  !------------
  character(3) :: expected_status, actual_status
  integer :: i, j
  integer :: ierr
  !------------
  info = 0
  !
  if (present(printinfo)) then
    l_print = printinfo
  else
    l_print = .false.
  endif
  !set up the record_handler variables for this file
  expected_status = 'OLD'
  call find_file(filename,file_indx,info,actual_file_status=actual_status)
  if (info.lt.0) return
  if (expected_status.ne.actual_status) then
    write(6,*)'ERROR: record_handler: In leave_record_handler'
    write(6,*)'Unexpected file status. Was expecting an OLD file'
    write(6,*)'but found a ',actual_status,' file.'
    print *,'Filename is ',trim(filename)
    info = -1
    return
  endif
  call set_pointers(file_indx)
  !
  !Close the toc
  if (file_is_new) then
    call bookkeeper(i,j,'END',info)
    if (info.lt.0) then
      print *,'leave_record_handler: Internal error in bookkeeper'
      print *,'File is ',trim(filename)
      info = -1
      return
    endif
    file_is_new = .false.
  endif
  call flush(fileunit)
  close(fileunit)
  !
  if (allocated(toc)) then
    deallocate(toc,stat=ierr)
    call check_deallocate(ierr,'toc',this_routine,__LINE__,.false.)
    if (ierr.ne.0) then; info = -1; return; endif
  endif
  !
  iam_init_toc = .false.
  iam_init_record_handler = .false.
  iam_init_bookkeeper = .false.
  !
  call nullify_pointers
  !
  if (l_print) write(6,*)'Leaving record_handler. File ',filename,' is closed.'
  !
  return
  end subroutine leave_record_handler

  subroutine io_check(filename,stat,info,record,openingfile)
  implicit none
  character(10), intent(in) :: filename
  integer, intent(in) :: stat
  integer, intent(out) :: info
  integer, optional, intent(in) :: record
  logical, optional, intent(in) :: openingfile
  logical :: lrec, lfile
  !-------
  lrec = .false.
  lfile = .false.
  if (present(record)) lrec = .true.
  if (present(openingfile)) then
    if (openingfile) lfile = .true.
  endif
  select case(stat)
  case(:-1)
    write(6,*)'ERROR: record_handler: End-of-file or End-of-record'
    if (lrec) then
      write(6,*)'encountered while reading file ',filename
      write(6,*)'Record being read was ',record
    elseif (lfile) then
      write(6,*)'encountered while opening file ',filename
      write(6,*)'Strange! How can this happen? '
    else
      write(6,*)'Not sure what the operation was when this error occured.'
    endif
    info = -1
    return
  case(1:)
    write(6,*)'ERROR: record_handler: File error other than EOR or EOF'
    if (lrec) then
      write(6,*)'encountered while reading file ',filename
      write(6,*)'Record being read was ',record
    elseif (lfile) then
      write(6,*)'encountered while opening file ',filename
      write(6,*)'Strange! How can this happen? '
    else
      write(6,*)'Not sure what the operation was when this error occured.'
    endif
    info = -1
    return
  case default
    !nothing. 
  end select
  !
  return
  end subroutine io_check

  subroutine bookkeeper(lcurr_rec,lnext_rec,newold,info)
  !update the current TOC (curr_toc), provide an updated current record that can
  !be written into...
  use bookkeeper_pointers
  implicit none
  save
  integer, intent(out) :: lcurr_rec, lnext_rec
  character(3), intent(in) :: newold        !beginning a 'NEW' record, or
                                            !continuing with an 'OLD' record
                                            !At end of file call with 'END'
  integer, intent(out) :: info
  integer, parameter :: max_in_curr_toc = len_one_toc - 1
                                            !max entries in current toc
  !------------
  if (newold.eq.'NEW') then
    if (.not.iam_init_bookkeeper) then
      !This is the beginning of a new type of data set to be written. 
      !perform the initialization
      !Zero out the current TOC matrix
      curr_toc = 0
      !the record that can be written into is the second one. 
      !The first is for the current TOC
      curr_rec = 2        !data can be written starting from record 2
      next_rec = 3        !next record available for writing
      num_in_curr_toc = 1 !number of entries in current toc so far
      curr_toc(num_in_curr_toc) = curr_rec
      rec_curr_toc = 1    !this toc is to be written into record 1
      iam_init_bookkeeper = .true.
      if (l_print) write(6,*)'Bookkeeper initialized'
    else
      !Perform the following tasks
      !(0) Update the curr_toc to specify the starting record for this data set
      !(1) update the current record 
      !(2) decide if there is more space in the current TOC or if a new one 
      !    has to be created. 
      !(3) If there is more space, update the next record position (next_rec) 
      !    that can be written to.
      !(4) If a new current TOC needs to be created, write curr_toc to disk, 
      !    and update curr_rec and next_rec.
      !
      curr_rec = next_rec
      if (num_in_curr_toc.lt.max_in_curr_toc) then
        !we have space remaining in the current toc...
        next_rec = curr_rec + 1
        num_in_curr_toc = num_in_curr_toc + 1 
        !Put entry into the current TOC
        curr_toc(num_in_curr_toc) = curr_rec
      else
        !There is no more space left.
        !Set the last entry in curr_toc to point to the record of the next toc
        !which is the record pointed at by curr_rec
        curr_toc(len_one_toc) = curr_rec
        !write out the current toc (curr_toc) to file
        call write_curr_toc(rec_curr_toc,info)
        if (info.lt.0) return
        !set a new value for rec_curr_toc
        rec_curr_toc = curr_rec
        !...and, of course, set new values for...
        curr_rec = curr_rec + 1
        next_rec = curr_rec + 1
        !initialize the new curr_toc
        curr_toc = 0
        num_in_curr_toc = 1
        curr_toc(num_in_curr_toc) = curr_rec
      endif
      !
    endif
    !
  elseif (newold.eq.'OLD') then
    !This means we are just writing more of an old data set into an additional
    !record (most data sets will utilize more than one record)
    curr_rec = next_rec
    next_rec = curr_rec + 1
  elseif (newold.eq.'END') then
    !A check
    if (.not.iam_init_bookkeeper) then
      curr_toc = 0
      rec_curr_toc = 1
    endif
    !The writing is complete. Write out the existing curr_toc
    call write_curr_toc(rec_curr_toc,info)
    if (info.lt.0) return
  else
    write(6,*)'Wrong value of newold in MODULE record_handler. ',newold
    info = -1
    return
  endif
  !
  lcurr_rec = curr_rec
  lnext_rec = next_rec
  !
  return
  contains
    subroutine write_curr_toc(rec_curr_toc,info)
    implicit none
    integer, intent(in) :: rec_curr_toc
    integer, intent(out) :: info
    !------------
    character(10) :: filename
    !------------
    ! 
    filename = ar_filename(file_indx)
    !
    write(fileunit,rec=rec_curr_toc,ERR=1) curr_toc
    !
    return
    !
    !Errors come here
1   write(6,*)'MODULE record_handler: TOC: Write-error (not EOF or EOR)'
    write(6,*)'file name =',filename,' record ',rec_curr_toc
    info = -1
    return
    !
    end subroutine write_curr_toc
    
  end subroutine bookkeeper

  subroutine find_file(filename,file_indx,info,actual_file_status,query)
  !determine the index of the file in the list. If it is not present, add it to
  !the list and return the file_indx. actual_file_status is the (perceived)
  !actual status of the file (OLD if it is already on the list, NEW otherwise).
  use common_routines
  implicit none
  character(10), intent(in) :: filename
  integer, intent(out) :: file_indx
  integer, intent(out) :: info
  character(3), intent(out), optional :: actual_file_status
  logical, intent(in), optional :: query
  !------------
  integer :: i
  logical :: l_query, l_status
  character(10) :: tmpfile
  !------------
  if (present(query)) then
    l_query = query
  else
    l_query = .false.
  endif
  if (present(actual_file_status)) then
    l_status = .true.
  else
    l_status = .false.
  endif
  !
  file_indx = 0
  if (l_status) actual_file_status = 'NEW' !default
  !
  if (numfiles.gt.0) then
    files: do i = 1, numfiles
      if (filename.eq.ar_filename(i)) then
        file_indx = i
        if (l_status) actual_file_status = 'OLD'
        exit files
      endif
    enddo files
  endif
  !
  if ((file_indx.eq.0).and.(.not.l_query)) then
    !the file was not found. Add it to the list if there is space.
    if (numfiles.lt.maxfiles) then
      numfiles = numfiles + 1
      file_indx = numfiles
      ar_filename(file_indx) = filename
      if (l_status) actual_file_status = 'NEW'
    else
      write(6,*)'MODULE record_handler: No space left in arrays to handle'
      write(6,*)'the file ',filename
      write(6,*)'Consider increasing the value of maxfiles in this module'
      call name_tmp_file(tmpfile,'out')
      !
      print *,'List of files in record_handler:'
      write(6,*)' Index       filename'
      write(6,*)'---------------------'
      do i = 1, maxfiles
        write(6,*)i,ar_filename(i) 
      enddo
      write(6,*)'====================='
      info = -1
      return
    endif
  endif
  !
  return
  end subroutine find_file

  subroutine initialize_file_variables
  use bookkeeper_pointers, only : curr_rec, next_rec
  implicit none
  iam_init_toc = .false.
  iam_init_record_handler = .false.
  iam_init_bookkeeper = .false.
  curr_rec = -1
  next_rec = -1
  return
  end subroutine initialize_file_variables

  subroutine set_pointers(file_indx)
  use bookkeeper_pointers
  implicit none
  integer, intent(in) :: file_indx
  !------------
  call nullify_pointers
  !
  fileunit => ar_fileunit(file_indx)
  curr_rec => ar_curr_rec(file_indx)
  next_rec => ar_next_rec(file_indx)
  iam_init_toc => ar_iam_init_toc(file_indx)
  iam_init_record_handler => ar_iam_init_record_handler(file_indx)
  iam_init_bookkeeper => ar_iam_init_bookkeeper(file_indx)
  file_is_new => ar_file_is_new(file_indx)
  curr_toc => ar_curr_toc(:,file_indx)
  rec_curr_toc => ar_rec_curr_toc(file_indx)
  num_in_curr_toc => ar_num_in_curr_toc(file_indx)
  !
  return
  end subroutine set_pointers

  subroutine nullify_pointers
  use bookkeeper_pointers
  implicit none
  nullify(fileunit,curr_rec,next_rec)
  nullify(iam_init_toc,iam_init_record_handler,iam_init_bookkeeper,file_is_new)
  nullify(curr_toc,rec_curr_toc,num_in_curr_toc)
  return
  end subroutine nullify_pointers

!------------------------
end module record_handler
!------------------------
