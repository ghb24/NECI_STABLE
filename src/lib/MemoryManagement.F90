#include "macros.h"
module MemoryManager
use constants , only : sizeof_int, dp, int64, int32

! JSS.  Memory book-keeping routines.  Contains a few elements of the initialisation,
! output and structure of the memory_manager module from CamCASP (formerly SITUS),
! written by Alston Misquitta, with permission.

!=====================================================================================

! USAGE
!
! See individual routines for more detail and optional arguments for
! error checking the allocation/deallocation.
!
! 1. Initialise the logger:
!      call InitMemoryManager(MemSize,print_err)
!    where the optional arguments are:
!      MemSize: max amount of memory available in MB. Default 1GB.
!      print_err: print error messages. Default true.
!
! 2. For each array to be logged, define an integer "tag" with an initial value of 0.
!      real(8), allocatable :: testarr(:)
!      integer :: tag_testarr=0
!    Make sure you save the tag (and the array!) if used for an array stored in a
!    module which can go out of scope between allocation and deallocation.
!
! 3. At allocation, call the allocation logger:
!      allocate(testarr(N))
!      call LogMemAlloc('testarr',N,8,'allocating_routine',tag_testarr)
!    where:
!      'testarr' is the name of the array being logged;
!       N is the size of the array;
!       8 is the number of bytes per element in the array;
!       'allocating_routine' is the routine in which the allocation is performed;
!       tag_testarr is the "tag" of the array being logged.
!
! 4. At deallocation, call the deallocation logger:
!      deallocate(testarr)
!      call LogMemDealloc('deallocating_routine',tag_testarr)
!    where:
!      'deallocating_routine' is the routine in which the deallocation is performed.
!
! 5. If, at any point you want to print (to STOUT) the allocated elements of the memory log:
!      call PrintMemory()
!    Optional arguments exist for printing out the deallocated elements and
!    printing to other units.
!
! 6. At the end of the calculation, print out memory usage:
!      call LeaveMemoryManager()

!=====================================================================================

! Log memory usage in one of two ways:
!   1. Store everything.  The size of the log (MaxLen) had best be suitably large.
!   2. Store active allocations.  When the top most slot in use in the log is
!   deallocated, free up all the slots at the top of the log that have been
!   deallocated for logging later actions (i.e. psuedo-LIFO [last in, first
!   out]).  This does allow a fractured log to occur, but memory
!   allocation-deallocation is often LIFO and it allows us to be efficient in
!   storing the actions without wasting much effort searching the log.
! This is controlled by the CachingMemLog (true==approach 2) flag.

! Using a "tag" for each allocated routine makes searching the log for an array
! trivial (and fast).

implicit none

!private
private
public :: MemoryLeft, MemoryUsed, MaxMemory, LookupPointer, PrintMemory
! Allow users to do the potentially dangerous thing of changing how the log is run.
! We'll hope they'll only use this for good...
public :: CachingMemLog
! Routines that need to be accessible.
public :: InitMemoryManager,LogMemAlloc,LogMemDealloc,LeaveMemoryManager
public :: TagIntType



integer, parameter :: TagIntType = sizeof_int   !This is for CPMD which needs to know what type of integer to pass as a tag
! Configuration.
integer, parameter :: MaxLen = 500000   ! size of memory log (max number of arrays
                                        ! that can be logged at any one time if
                                        ! CachingMemLog=.true., else max total number of
                                        ! arrays that can be logged in a calculation).
integer, parameter :: MaxWarn = 10     ! maximum number of low memory warning messages to be printed.
integer, parameter :: nLargeObjects = 10 ! maximum number of the largest memory allocations remember.
logical, save :: CachingMemLog = .true. ! See above for how MemLog is used.
logical, save :: MemUnitsBytes = .true. ! If true, then output object size in bytes/KB/MB.
                                        ! If false, then output in words (an ugly unit!).

logical, save :: err_output = .true. ! Print error messages.

type :: MemLogEl
    character(len=25) :: ObjectName=''
    character(len=25) :: AllocRoutine=''
    character(len=25) :: DeallocRoutine='not deallocated'
    integer(int64) :: ObjectSize=0
end type MemLogEl

! All memory variables stored in bytes.
integer(int64), save :: MaxMemory
integer(int64), save :: MemoryUsed
integer(int64), save :: MemoryLeft
integer(int64), save :: MaxMemoryUsed

! Warnings, debug flags, output parameters.
logical, save :: initialised = .false.
logical, save :: warned = .false.
integer, save :: nWarn = 0
logical, save :: debug = .false.

! Log of memory allocations.
type(MemLogEl), allocatable, save :: MemLog(:)
integer, save :: ipos=1  ! Next available empty slot in the log.

! Capture the state of the MemLog at peak usage.  Currently not outputted, but
! useful for diagnostics.
type(MemLogEl), allocatable, save :: PeakMemLog(:)

type(MemLogEl), save :: LargeObjLog(nLargeObjects) ! Store the largest allocations.
integer, save :: ismall=1 ! The smallest large object (remember to avoid repeating minloc again and again...)

! For backwards compatibility with the existing CPMD scheme, where the IP address is
! stored as the tag. Use long integer (int64) so can handle POINTER8.
integer(int64), allocatable, save :: LookupPointer(:)

! Log a memory allocation.
! INPUT:
!       ObjectName - Name of object.
!       ObjectSize - Number of elements in object.
!       ElementSize - Number of bytes per element.
!       AllocRoutine - routine in which object is allocated.
!       err (optional) - error output from allocate statement (checked if present).
! OUTPUT:
!       tag - position in memory log the object is stored at.
!             If -1, then the log is full and it's not been stored.
! IN/OUT:
!       nCalls (optional) -  increments nCalls: counts the number of times
!       a routine has called the LogMemAlloc routine (useful for tracking
!       repeated allocations in debugging).
!
! Details:
!       ObjectSize and ElementSize can either be both int32 or int64
interface LogMemAlloc
    module procedure LogMemAlloc_int32, LogMemAlloc_int64
end interface

contains

    subroutine InitMemoryManager(MemSize, print_err)
    ! Initialise memory manager.

    ! In:
    !    MemSize (optional) : max amount of memory available in MB.  Default: MaxMemLimit
    !    print_err (optional): print all error messages from memory manager. Default: true.

    ! MAXMEM must be set via c pre-processing or set to be an integer.

    implicit none

    integer(int64), intent(in), optional :: MemSize
    logical, intent(in), optional :: print_err
    integer(int64) :: MaxMemBytes
    ! Obtained via CPP in the makefile. MAXMEM in MB.
    integer(int64), parameter :: MaxMemLimit = MAXMEM

    def_default(err_output, print_err, .true.)
    if (present(MemSize)) then
        MaxMemBytes = MemSize * 1024**2
    else
        MaxMemBytes = MaxMemLimit * 1024**2
    end if

    if (initialised) then
        if (err_output) write (6,*) 'Already initialised memory manager.  Not re-initialsing.'
    else
        if (MaxMemBytes.le.0) then
            if (err_output) then
                write (6,*) 'Illegal maximum memory value passed to memorymanager.'
                write (6,*) 'MaxMemgbytes = ',real(MaxMemBytes,dp)/(1024**2)
                write (6,*) 'Setting maximum memory available to 1GB.'
            end if
            MaxMemBytes=1024**3
        endif


        if(.not.allocated(MemLog)) allocate(MemLog(MaxLen))
        if(.not.allocated(PeakMemLog)) allocate(PeakMemLog(MaxLen))
        if(.not.allocated(LookupPointer)) allocate(LookupPointer(MaxLen))
        lookuppointer = 0
        MaxMemory = MaxMemBytes
        MemoryUsed = 0
        MemoryLeft = MaxMemory
        MaxMemoryUsed = 0
        initialised = .true.
        nWarn = 0
!       Deal with debug options at a later date.
!       debug = gmemdebug

        write (6,'(a33,f8.1,a3)') ' Memory Manager initialised with ',real(MaxMemBytes,dp)/(1024**2),' MB'
    end if

    return
    end subroutine InitMemoryManager


    subroutine LogMemAlloc_int32(ObjectName,ObjectSize,ElementSize,AllocRoutine,tag,err,nCalls)
        character(len=*), intent(in) :: ObjectName, AllocRoutine
        integer(int32), intent(in) :: ObjectSize
        integer(int32), intent(in) :: ElementSize
        integer(TagIntType), intent(out) :: tag
        integer, intent(in), optional :: err
        integer, intent(inout), optional :: nCalls
        call LogMemAlloc(ObjectName, int(ObjectSize, int64), int(ElementSize, int64), &
                         AllocRoutine, tag, err, nCalls)
    end subroutine

    subroutine LogMemAlloc_int64(ObjectName,ObjectSize,ElementSize,AllocRoutine,tag,err,nCalls)
    implicit none

    character(len=*), intent(in) :: ObjectName, AllocRoutine
    integer(int64), intent(in) :: ObjectSize
    integer(int64), intent(in) :: ElementSize
    integer(TagIntType), intent(out) :: tag
    integer, intent(in), optional :: err
    integer, intent(inout), optional :: nCalls

    integer(int64) :: ObjectSizeBytes
    integer :: ismallloc(1)
    character(*), parameter :: this_routine = 'LogMemAlloc'
    external :: warning_neci
    if (present(nCalls)) nCalls=nCalls+1

    if (.not. initialised) then
        write (6,*) 'Memory manager not initialised. Doing so now with 1GB limit.'
        call InitMemoryManager()
    end if

    ObjectSizeBytes=ObjectSize*ElementSize

    MemoryUsed=MemoryUsed+ObjectSizeBytes
    MaxMemoryUsed=max(MemoryUsed,MaxMemoryUsed)
    MemoryLeft=MaxMemory-MemoryUsed

    if (MemoryLeft.lt.0.and.nWarn.lt.MaxWarn) then
        if (err_output) write (6,*) 'WARNING: Memory used exceeds maximum memory set',MemoryLeft
        nWarn=nWarn+1
    end if

    if (present(err)) then
        if (err /= 0) then
            call stop_all(this_routine, 'Failure to allocate array '//ObjectName//' in '//AllocRoutine)
        end if
    end if

    if (ipos.gt.MaxLen) then
        if (.not.warned) then
            warned=.true.
            if (err_output) then
                write (6,*) 'Warning: Array capacity of memory manager exceeded.'
                write (6,*) 'Required array length is ',ipos
                write (6,*) 'Max memory used is likely to be incorrect.'
            end if
        end if
        tag=-1
        ! If we're not putting it in the log, test if it's a huge array:
        ! it's always the biggest fishes that get away!
        if (ObjectSizeBytes > LargeObjLog(ismall)%ObjectSize) then
            LargeObjLog(ismall)%ObjectName=ObjectName
            LargeObjLog(ismall)%AllocRoutine=AllocRoutine
            LargeObjLog(ismall)%ObjectSize=ObjectSizeBytes
            ismallloc=minloc(LargeObjLog(:)%ObjectSize)
            ismall=ismallloc(1)
        end if
    else
        MemLog(ipos)%ObjectName=ObjectName
        MemLog(ipos)%AllocRoutine=AllocRoutine
        MemLog(ipos)%DeallocRoutine='not deallocated' ! In case this slot has already been used and abandoned.
        MemLog(ipos)%ObjectSize=ObjectSizeBytes
        tag=ipos
        ipos=ipos+1
    end if

    if (debug) then
        write (6,"(A,I6,I12,' ',A,' ',A,' ',I12)") 'Allocating memory: ',tag,ObjectSizeBytes,ObjectName,AllocRoutine,MemoryUsed
    end if

    return
    end subroutine LogMemAlloc_int64



    subroutine LogMemDealloc(DeallocRoutine,tag,err)
    ! Log a memory deallocation.
    ! INPUT:
    !       DeallocRoutine - routine in which object is deallocated.
    !       tag - position in memory log the object is stored at.
    !       err (optional) - error output from deallocate statement (checked if present).
    ! OUTPUT:
    !       tag - 0 if successfully logged (or noted that it wasn't stored in
    !             the log in the first place).

    implicit none

    character(len=*), intent(in) :: DeallocRoutine
    integer(TagIntType), intent(inout) :: tag
    integer, intent(in), optional :: err
    integer :: i,ismallloc(1)
    character(len=25) :: ObjectName
    external :: stop_all
    if (.not.initialised) then
        if (err_output) write (6,*) 'Memory manager not initialised. Cannot log deallocation.'
        return
    end if

    ObjectName='Unknown'

    if (tag.eq.0) then
        if (err_output) write (6,*) 'Warning: attempting to log deallocation but never logged allocation.'
        tag=-1
    else if(tag.gt.MaxLen.or.tag.lt.-1) then
        if (err_output) write (6,*) 'Warning: attempting to log deallocation but tag does not exist: ',tag
        tag=-1
    else

        if (MemoryUsed.eq.MaxMemoryUsed) then
            ! Are at peak memory usage.  Copy the memory log.
            ! Useful to see what's around when memory usage is at a maximum.
            PeakMemLog(:)=MemLog(:)
        end if

        if (tag.eq.-1) then
            ! No record of it in the log: can only print out a debug message.
            if (debug) then
                write (6,"(2A,I5)") 'Deallocating memory in: ',DeallocRoutine,tag
            end if
        else
            MemoryUsed=MemoryUsed-MemLog(tag)%ObjectSize
            MemoryLeft=MaxMemory-MemoryUsed

            ! Object was stored in the cache.
            MemLog(tag)%DeallocRoutine=DeallocRoutine
            ObjectName=MemLog(tag)%ObjectName

            ! Check to see if object is larger than the smallest of the large
            ! objects: if so, keep a record of it.
            if (MemLog(tag)%ObjectSize.gt.LargeObjLog(ismall)%ObjectSize) then
                LargeObjLog(ismall)=MemLog(tag)
                ismallloc=minloc(LargeObjLog(:)%ObjectSize)
                ismall=ismallloc(1)
            end if

            if (CachingMemLog) then
                ! Then we free up this slot and slots of all objects directly below it in
                ! the log that have also been deallocated.  This is not the most
                ! efficient storage (we still can have a fractured log) but
                ! works well for LIFO approaches, which are most common for us.
                MemLog(tag)%ObjectSize=0 ! Nothing to see here now.
                if (tag.eq.ipos-1) then
                    do i=tag,1,-1
                        if (MemLog(i)%ObjectSize.eq.0) then
                            ipos=ipos-1
                        else
                            exit
                        end if
                    end do
                end if
            end if
            if (debug) then
                write (6,"(A,I5,' ',A,' ',A,' ',A,' ',I12)") 'Deallocating memory: ',tag,MemLog(tag)
            end if
        end if

        ! Set tag to zero: there was no problem with the logging deallocation
        ! (apart from maybe a too small cache).
        tag=0

    end if

    if (present(err)) then
        if (err.ne.0) then
            call Stop_All('LogMemAlloc','Failed to deallocate array '//ObjectName//' in '//DeallocRoutine)
        end if
    end if


    return
    end subroutine LogMemDealloc



    subroutine LeaveMemoryManager()
    ! Call this to print out the largest memory allocations.
    ! If debug flag is on, then the full memory log is dumped to file.

    implicit none

    integer :: iunit,iobjloc(1),iobj,i
    integer(int64), allocatable :: ObjectSizes(:)
    type(MemLogEl), allocatable :: AllMemEl(:)
    character(len=*), parameter :: memoryfile = 'TMPMemoryusage.dat'
    character(len=*), parameter :: fmt1='(3a19)'

    if (.not.initialised) then
        if (err_output) write (6,*) 'Memory manager not initialised. Cannot leave memory manager.'
        return
    end if

    allocate(ObjectSizes(nLargeObjects+MaxLen))
    allocate(AllMemEl(nLargeObjects+MaxLen))

    if (MemoryUsed.eq.MaxMemoryUsed) then
        ! Peak memory usage is now.
        PeakMemLog(:)=MemLog(:)
    end if

    call WriteMemLogHeader(6)

    if (CachingMemLog) then
        ! Large objects might be residing in the MemLog, but not deallocated
        ! (and so haven't been moved to the large object store).
        AllMemEl(1:MaxLen)=MemLog
        AllMemEl(MaxLen+1:MaxLen+nLargeObjects)=LargeObjLog
    else
        ! Everything really ought to be held in just the MemLog: if not, then
        ! this is a "feature".
        AllMemEl(1:MaxLen)=MemLog(:)
    end if

    ! Copy the sizes to an integer array: we use maxloc on the copy.  This
    ! allows us to check for arrays of the same size without writing over
    ! information in our log.
    ObjectSizes(:)=AllMemEl(:)%ObjectSize
    iobjloc(:)=maxloc(ObjectSizes)
    iobj=iobjloc(1)
    ObjectSizes(iobj)=ObjectSizes(iobj)+1
    do i=2,nLargeObjects+1
        ! Print out i-1 large object.
        write (6,fmt1,advance='no') ' '//AllMemEl(iobj)%ObjectName,AllMemEl(iobj)%AllocRoutine,AllMemEl(iobj)%DeallocRoutine
        call WriteMemSize(6,AllMemEl(iobj)%ObjectSize)
        ! Find the next large object.
        iobjloc=maxloc(ObjectSizes,mask=ObjectSizes.lt.ObjectSizes(iobj))
        iobj=iobjloc(1)
        if (AllMemEl(iobj)%ObjectName.eq.''                                &
             .and.AllMemEl(iobj)%AllocRoutine.eq.''                        &
             .and.AllMemEl(iobj)%DeallocRoutine.eq.'not deallocated'       &
             .and.AllMemEl(iobj)%ObjectSize.eq.0) then
             ! Have logged less than nLargeObjects allocations.
             exit
         end if
         ObjectSizes(iobj)=ObjectSizes(iobj)+1 ! So we don't find this object next time round.
    end do
    if (warned) then
        write (6,*) '== NOTE: Length of logging arrays exceeded. Length needed is ',ipos
    endif
    write (6,*) '================================================================'

    if (debug) then
        ! Dump entire memory log to file.
        iunit=93
!        call get_free_unit(iunit)  !Avoid circular dependancies - hack.
        open(unit=iunit,file=memoryfile,form='formatted',status='unknown')
        call PrintMemory(.true.,iunit)
        close(iunit)
    end if

    initialised=.false.
    deallocate(MemLog)
    deallocate(PeakMemLog)
    deallocate(LookupPointer)

    deallocate(ObjectSizes)
    deallocate(AllMemEl)

    end subroutine LeaveMemoryManager



    subroutine PrintMemory(PrintDeallocated,iunit)
    ! Print out the memory log.  If using the cache memory log, then it will
    ! only print out the elements stored, *not* any which have been
    ! over-written.  By default, only the active allocations are printed out to
    ! STOUT.
    !
    ! INPUT:
    !      PrintDeallocated (optional, default=.false.) - print out objects which
    !                       have been deallocated (but are still present in the cache).
    !      iunit (optional, default=6) - unit to output to.

    implicit none
    logical, intent(in), optional :: PrintDeallocated
    integer, intent(in), optional :: iunit

    logical :: pd
    integer :: io,iobj
    character(len=*), parameter :: fmt1='(3a19)'

    if (present(PrintDeallocated)) then
        pd=PrintDeallocated
    else
        pd=.false.
    end if

    if (present(iunit)) then
        io=iunit
    else
        io=6
    end if

    call WriteMemLogHeader(io)
    do iobj = 1, min(ipos,MaxLen)
        if (pd.or.MemLog(iobj)%DeallocRoutine.eq.'not deallocated') then
            write (io,fmt1,advance='no') ' '//MemLog(iobj)%ObjectName,MemLog(iobj)%AllocRoutine,MemLog(iobj)%DeallocRoutine
            call WriteMemSize(io,MemLog(iobj)%ObjectSize)
        end if
    enddo
    if (ipos.gt.MaxLen) then
        write (io,*) '== NOTE: Length of logging arrays exceeded. Length needed is ',ipos
    endif
    write (io,*) '================================================================'

    return
    end subroutine PrintMemory



    subroutine WriteMemLogHeader(iunit)
    implicit none
    integer, intent(in) :: iunit
    write (iunit,*)
    write (iunit,*) '================================================================'
    write (iunit,*) 'Memory usage'
    write (iunit,'(a34,f9.1)') ' Maximum memory defined is (MB) : ',real(MaxMemory,dp)/1024.0_dp**2.0_dp
    write (iunit,'(a34,f9.1)') ' Maximum memory used is    (MB) : ',real(MaxMemoryUsed,dp)/1024.0_dp**2.0_dp
    if (nWarn.gt.0) then
        write (iunit,*)'Maximum memory exceeded ',nWarn,' times.'
    endif
    write (iunit,*) ''
    if (iunit.eq.6) then
        write (iunit,*) 'Large memory allocations:'
        write (iunit,*) ''
    end if
    write (iunit,*) 'Name              Allocated in       Deallocated in         Size'
    write (iunit,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    return
    end subroutine WriteMemLogHeader



    subroutine WriteMemSize(iunit, MemSize)
        ! Write out a human-readable amount of memory.  MemSize is in bytes.
        implicit none
        integer, intent(in) :: iunit
        integer(int64), intent(in) :: MemSize
        character(len=*), parameter :: fmt1='(f6.1, a2)'
        character(len=*), parameter :: fmt2='(i7, a1)'
        if (MemUnitsBytes) then
            if (MemSize < 1024) then
                ! output in KB.
                write (iunit, fmt1) real(MemSize, dp) / 1024, 'KB'
            else if (MemSize < 1024**2) then
                ! output in MB.
                write (iunit, fmt1) real(MemSize, dp) / 1024**2, 'MB'
            else if (MemSize < 1024**3) then
                ! output in GB.
                write (iunit, fmt1) real(MemSize, dp) / 1024**3, 'GB'
            else if (MemSize < 1024_int64**4) then
                ! output in GB.
                write (iunit, fmt1) real(MemSize, dp) / 1024_int64**4, 'TB'
            else
                write (iunit, '(A)') '> 1 PB'
            end if
        else
            write (iunit, fmt2) MemSize / 8, 'W'
        end if
    end subroutine WriteMemSize



end module MemoryManager
