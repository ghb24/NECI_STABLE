MODULE MolproPlugin
 USE iso_c_binding
 IMPLICIT NONE
 PRIVATE
 PUBLIC :: MolproPluginInit, MolproPluginterm, MolproPluginResult
 PUBLIC :: molpro_plugin, molpro_plugin_fcidumpname, molpro_plugin_datafilename
 LOGICAL :: molpro_plugin
 INTEGER, PARAMETER :: PATHLEN = 1024
 CHARACTER(PATHLEN) :: molpro_plugin_fcidumpname, molpro_plugin_datafilename

! we have to drive the C implementation of PluginGuest because neci's Fortran wrappers for MPI are not compliant with Fortran MPI code
 INTERFACE
  SUBROUTINE PluginGuestOpen(host) BIND(C, name='PluginGuestOpen')
   USE iso_c_binding
   CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: host
  END SUBROUTINE PluginGuestOpen
  FUNCTION PluginGuestActive() BIND(C, name='PluginGuestActive')
   USE iso_c_binding
   INTEGER(c_int) :: PluginGuestActive
  END FUNCTION PluginGuestActive
  FUNCTION PluginGuestSend(value) BIND(C, name='PluginGuestSend')
   USE iso_c_binding
   CHARACTER(c_char), DIMENSION(*), INTENT(in) :: value
   INTEGER(c_int) :: PluginGuestSend
  END FUNCTION PluginGuestSend
  SUBROUTINE PluginGuestReceive(str,len) BIND(C, name='PluginGuestReceive')
   USE iso_c_binding
   character(kind=c_char), dimension(*) :: str
   integer(kind=c_size_t), value :: len
  END SUBROUTINE PluginGuestReceive
  SUBROUTINE PluginGuestClose() BIND(C, name='PluginGuestClose')
  END SUBROUTINE PluginGuestClose
 END INTERFACE
CONTAINS
 FUNCTION PluginGuestReceiveF()
  CHARACTER(:), ALLOCATABLE :: PluginGuestReceiveF
  integer :: length
  character(kind=c_char,len=1), allocatable, dimension(:) :: cstring
  integer(kind=c_size_t) :: clen
   
  clen=int(PATHLEN,kind=c_size_t)
  allocate(cstring(PATHLEN+1)) 
  call PluginGuestReceive(cstring,clen)
  DO length=1,PATHLEN
   IF (cstring(length).EQ.c_null_char) EXIT
  END DO
  length = length-1
  ALLOCATE (CHARACTER(len=length) :: PluginGuestReceiveF)
  DO length=1,LEN(PluginGuestReceiveF)
   PluginGuestReceiveF(length:length) = cstring(length)
  END DO
 END FUNCTION PluginGuestReceiveF
 FUNCTION PluginGuestSendF(value)
  LOGICAL :: PluginGuestSendF
  CHARACTER(*), INTENT(in) :: value
  CHARACTER(kind=c_char), DIMENSION(:), ALLOCATABLE :: cstring
  INTEGER :: i
  ALLOCATE(cstring(LEN_TRIM(value)+1))
  DO i=1,LEN_TRIM(value)
   cstring(i) = value(i:i)
  END DO
  cstring(LEN_TRIM(value)+1) = c_null_char
  PluginGuestSendF = PluginGuestSend(cstring).NE.0
  DEALLOCATE(cstring)
 END FUNCTION PluginGuestSendF
SUBROUTINE MolproPluginInit (tMolpro)
 IMPLICIT NONE
 LOGICAL, INTENT(inout) :: tMolpro
 CHARACTER(1024) :: id
 CHARACTER(c_char), DIMENSION(7) :: host = ['M','O','L','P','R','O',c_null_char]
 ! is this a Molpro plugin?
 CALL PluginGuestOpen(host)
 tMolpro = PluginGuestActive().NE.0
 molpro_plugin=.FALSE.
 IF (.NOT. tMolpro) RETURN
 molpro_plugin=.TRUE.
 
! ask for an FCIDUMP
 IF (.NOT. PluginGuestSendF('GIVE OPERATOR HAMILTONIAN FCIDUMP')) STOP 'plugin request has failed'
 molpro_plugin_fcidumpname = PluginGuestReceiveF()

! ask for a data file
 IF (.NOT. PluginGuestSendF('GIVE INPUT NECI')) STOP 'plugin request has failed'
 molpro_plugin_datafilename = PluginGuestReceiveF()

 IF (.FALSE.) THEN ! debugging
  WRITE (6, '(''Dump file: '',A)') TRIM(molpro_plugin_fcidumpname)
  WRITE (6, '(''Input file: '',A)') TRIM(molpro_plugin_datafilename)
  OPEN(1,file=molpro_plugin_datafilename,status='OLD')
  DO WHILE (.TRUE.)
   READ(1,'(A)',END=99) id
   WRITE (6, '(A)') TRIM(id)
  END DO
99 CLOSE(1)
 END IF
 
 END SUBROUTINE MolproPluginInit

SUBROUTINE MolproPluginTerm(signal)
 USE iso_c_binding, ONLY : c_int
 INTEGER, INTENT(in) :: signal
 INTERFACE
  SUBROUTINE fsleep(seconds) BIND(C,name="sleep")
   IMPORT
   INTEGER(c_int), value :: seconds
  END SUBROUTINE fsleep
 END INTERFACE
 IF (PluginGuestActive().NE.0) THEN
! Graceful exit if Molpro server
  CALL PluginGuestClose
  ! without this print, then MPI gets lost ???
  WRITE (6,*) 'Stopping Molpro plugin, signal =',signal; FLUSH(6)
! doesn't look like slave threads ever make it here, so do not have a barrier
!  IF (signal.EQ.0) THEN
!   CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
!  ELSE
  ! for abnormal termination, signal.lt.0, and then don't risk a barrier
   CALL fsleep(1_c_int) ! give the message time to arrive
!  END IF
 END IF
END SUBROUTINE MolproPluginTerm

SUBROUTINE MolproPluginResult(property,values)
 IMPLICIT NONE
 CHARACTER(*), INTENT(in) :: property
 DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: values
 CHARACTER(:), ALLOCATABLE :: buffer
 IF (PluginGuestActive().EQ.0) RETURN
 IF (PluginGuestSendF('TAKE PROPERTY '//TRIM(property)//c_null_char)) THEN
  ALLOCATE(CHARACTER(24*size(values)) :: buffer)
  WRITE (buffer,'(1000(G23.16,1X))') values
  IF (.NOT. PluginGuestSendF(buffer//c_null_char)) STOP 'Failure to send results to plugin host'
  DEALLOCATE(buffer)
 END IF
END SUBROUTINE MolproPluginResult

END MODULE MolproPlugin
