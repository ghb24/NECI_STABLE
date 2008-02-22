! JSS.  Memory book-keeping routines.  Contains some elements of the memory_manager
! module from CamCASP (formerly SITUS), written by Alston Misquitta, with
! permission.

module MemoryManager

implicit none

private
public :: MemoryLeft, MemoryUsed,MaxMemory,li,LookupPointer
public :: InitMemoryManager,LogMemAlloc,LogMemDealloc,LeaveMemoryManager

integer, parameter :: li = selected_int_kind(18) !ints between +-10^18

type MemLogEl
    character(len=25) :: ObjectName=''
    character(len=25) :: AllocRoutine=''
    character(len=25) :: DeallocRoutine='not deallocated'
    integer(li) :: ObjectSize=0
end type MemLogEl

! All memory variables stored in bytes.
integer(li), save :: MaxMemory
integer(li), save :: MemoryUsed
integer(li), save :: MemoryLeft
integer(li), save :: MaxMemoryUsed

! Warnings, debug flags, output parameters.
logical, save :: initialised = .false.
logical, save :: warned = .false.
integer, parameter :: MaxWarn = 10 ! maximum number of low memory warning messages to be printed.
integer, parameter :: nLargeObjects = 10 ! maximum number of large memory allocations to print out.
integer, save :: nWarn
logical, save :: debug = .false.

! Log of memory allocations.
integer, parameter :: MaxLen = 5000
type(MemLogEl), save :: MemLog(MaxLen)
integer, save :: ipos=1  ! Next available empty slot in the log.
! For backwards compatibility with the existing scheme, where the IP address is
! stored.
integer(li), save :: LookupPointer(MaxLen)

contains

    subroutine InitMemoryManager(MaxMemBytes)
    ! Initialise memory manager.

    use common_routines, only: internal_error
    implicit none

    integer(li), intent(in) :: MaxMemBytes
    character(len=*), parameter :: ThisRoutine = 'InitMemoryManager'

    if (MaxMemBytes.le.0) then
        write (6,*) 'Illegal maximum memory value passed to memorymanager.'
        write (6,*) 'MaxMemgbytes = ',dfloat(MaxMemBytes)/(1024**2)
        call internal_error(ThisRoutine,__LINE__,&
           & 'Illegal maximum memory. Check MEMORY in your input file.')
    endif

    MaxMemory = MaxMemBytes
    MemoryUsed = 0
    MemoryLeft = MaxMemory
    MaxMemoryUsed = 0
    initialised = .true.
    nWarn = 0
!   Deal with debug options at a later date.
!   debug = gmemdebug

    write (6,'(a33,f8.1,a3)') ' Memory Manager initialised with ',dfloat(MaxMemBytes)/(1024**2),' MB'

    return
    end subroutine InitMemoryManager



    subroutine LogMemAlloc(ObjectName,ObjectSize,ElementSize,AllocRoutine,tag)
    ! Log a memory allocation.
    ! INPUT:
    !       ObjectName - Name of object.
    !       ObjectSize - Number of elements in object.
    !       ElementSize - Number of bytes per element.
    !       AllocRoutine - routine in which object is allocated.
    ! OUTPUT:
    !       tag - position in memory log the object is stored at.
    !             If -1, then the log is full and it's not been stored.

    implicit none

    character(len=*),intent(in) :: ObjectName,AllocRoutine
    integer, intent(in) :: ObjectSize
    integer, intent(in) :: ElementSize
    integer, intent(out) :: tag
 
    integer(li), parameter :: DefaultMem=1024**3
    integer :: ObjectSizeBytes

    if (.not.initialised) then
        write (6,*) 'Memory manager not initialised. Doing so now with 1GB limit.'
        call InitMemoryManager(DefaultMem)
    end if

    ObjectSizeBytes=ObjectSize*ElementSize

    MemoryUsed=MemoryUsed+ObjectSizeBytes
    MaxMemoryUsed=max(MemoryUsed,MaxMemoryUsed)
    MemoryLeft=MaxMemory-MemoryUsed

    if (MemoryLeft.lt.0.and.nWarn.lt.MaxWarn) then
        write (6,*) 'WARNING: Memory used exceeds maximum memory set',MemoryLeft
        nWarn=nWarn+1
    end if

    if (debug) then
        write (6,*) 'Allocating memory: ',tag,ObjectSizeBytes,ObjectName,AllocRoutine,MemoryUsed
    end if

    if (ipos.gt.MaxLen) then
        if (.not.warned) then
            warned=.true.
            write (6,*) 'Warning: Array capacity of memory manager exceeded.'
            write (6,*) 'Required array length is ',ipos
            write (6,*) 'Max memory used is likely to be incorrect.'
        end if
        tag=-1
    else
        MemLog(ipos)%ObjectName=ObjectName
        MemLog(ipos)%AllocRoutine=AllocRoutine
        MemLog(ipos)%ObjectSize=ObjectSizeBytes
        tag=ipos
        ipos=ipos+1
    end if

    return
    end subroutine LogMemAlloc



    subroutine LogMemDealloc(DeallocRoutine,tag)
    ! Log a memory deallocation.
    ! INPUT:
    !       DeallocRoutine - routine in which object is deallocated.
    !       tag - position in memory log the object is stored at.
    ! OUTPUT:
    !       tag - 0 if successfully logged (or noted that it wasn't stored in
    !             the log in the first place).

    implicit none

    character(len=*), intent(in) :: DeallocRoutine
    integer, intent(inout) :: tag

    if (tag.eq.0) then
        write (6,*) 'Warning: attempting to log deallocation but never logged allocation.'
        tag=-1
    else if(tag.gt.MaxLen.or.tag.lt.-1) then
        write (6,*) 'Warning: attempting to log deallocation but tag does not exist: ',tag
        tag=-1
    else if (tag.eq.-1) then
        tag=0
    else
        MemoryUsed=MemoryUsed-MemLog(tag)%ObjectSize
        MemoryLeft=MaxMemory-MemoryUsed
        MemLog(tag)%DeallocRoutine=DeallocRoutine
        if (debug) then
            write (6,*) 'Deallocating memory: ',tag,MemLog(tag)
        end if
        tag=0
    end if

    return
    end subroutine LogMemDealloc



    subroutine LeaveMemoryManager
    ! Call this to print out the largest memory allocations.
    ! If debug flag is on, then the full memory log is dumped to file.

    use common_routines, only: getunit
    implicit none

    integer :: iunit,iobjloc(1),iobj,i
    integer(li) :: ObjectSizes(MaxLen)
    character(len=*), parameter :: memoryfile = 'TMPMemoryusage.dat'
    character(len=*), parameter :: fmt1='(3a19,f8.1)'

    call WriteMemLogHeader(6)
    ObjectSizes(:)=MemLog(:)%ObjectSize
    iobjloc(:)=maxloc(ObjectSizes)
    iobj=iobjloc(1)
    ObjectSizes(iobj)=ObjectSizes(iobj)+1
    do i=2,min(nLargeObjects,ipos)
        write (6,fmt1) ' '//MemLog(iobj)%ObjectName,MemLog(iobj)%AllocRoutine,MemLog(iobj)%DeallocRoutine,dfloat(MemLog(iobj)%ObjectSize)/1024**2
        iobjloc=maxloc(ObjectSizes,mask=ObjectSizes.lt.ObjectSizes(iobj))
        iobj=iobjloc(1)
        ObjectSizes(iobj)=ObjectSizes(iobj)+1
    end do
    if (warned) then
        write (6,*) '== NOTE: Length of logging arrays exceeded. Length needed is ',ipos
    endif
    write (6,*) '================================================================'

    if (debug) then
        ! Dump entire memory log to file.
        call getunit(iunit)
        open(unit=iunit,file=memoryfile,form='formatted',status='unknown')
        call WriteMemLogHeader(iunit)
        do iobj = 1, min(ipos,MaxLen)
            write (6,fmt1) ' '//MemLog(iobj)%ObjectName,MemLog(iobj)%AllocRoutine,MemLog(iobj)%DeallocRoutine,dfloat(MemLog(iobj)%ObjectSize)/1024**2
        enddo
        if (warned) then
            write (iunit,*) '== NOTE: Length of logging arrays exceeded. Length needed is ',ipos
        endif
        write (iunit,*) '================================================================'
        close(iunit)
    end if

    return
    end subroutine LeaveMemoryManager

    

    subroutine WriteMemLogHeader(iunit)
    implicit none
    integer :: iunit
    write(iunit,*)
    write (iunit,*) '================================================================'
    write (iunit,*) 'Memory usage'
    write (iunit,'(a34,f9.1)') ' Maximum memory defined is (MB) : ',dfloat(MaxMemory)/1024**2
    write (iunit,'(a34,f9.1)') ' Maximum memory used is    (MB) : ',dfloat(MaxMemoryUsed)/1024**2
    if (nWarn.gt.0) then
        write (iunit,*)'Maximum memory exceeded ',nWarn,' times.'
    endif
    write (iunit,*) ''
    if (iunit.eq.6) then
        write (iunit,*) 'Large memory allocations:'
        write (iunit,*) ''
    end if
    write (iunit,*) 'Name              Allocated in       Deallocated in     Size(MB)'
    write (iunit,*) '----------------------------------------------------------------'
    return
    end subroutine WriteMemLogHeader



end module MemoryManager
