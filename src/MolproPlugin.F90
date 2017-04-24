MODULE MolproPlugin
 USE iso_c_binding
 IMPLICIT NONE
 PRIVATE
 PUBLIC :: MolproPluginInit, MolproPluginterm, MolproPluginResult
 PUBLIC :: molpro_plugin, molpro_plugin_fcidumpname, molpro_plugin_datafilename
 LOGICAL :: molpro_plugin
 CHARACTER(1024) :: molpro_plugin_fcidumpname, molpro_plugin_datafilename

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
  FUNCTION PluginGuestMaster() BIND(C, name='PluginGuestMaster')
   USE iso_c_binding
   INTEGER(c_int) :: PluginGuestMaster
  END FUNCTION PluginGuestMaster
  FUNCTION PluginGuestSend(value) BIND(C, name='PluginGuestSend')
   USE iso_c_binding
   CHARACTER(c_char), DIMENSION(:), INTENT(in) :: value
   INTEGER(c_int) :: PluginGuestSend
  END FUNCTION PluginGuestSend
  FUNCTION PluginGuestReceive() BIND(C, name='PluginGuestReceive')
   USE iso_c_binding
   !CHARACTER(c_char), DIMENSION(:), ALLOCATABLE :: PluginGuestReceive
   TYPE(c_ptr) :: PluginGuestReceive
  END FUNCTION PluginGuestReceive
  SUBROUTINE PluginGuestClose() BIND(C, name='PluginGuestClose')
  END SUBROUTINE PluginGuestClose
 END INTERFACE
CONTAINS
 FUNCTION PluginGuestReceiveF()
  CHARACTER(:), ALLOCATABLE :: PluginGuestReceiveF
  TYPE(c_ptr) :: p
  CHARACTER, POINTER, DIMENSION(:) :: fp
  integer :: length
  p = PluginGuestReceive()
  CALL c_f_pointer(p,fp, [1])
  DO length=1,1000000
   IF (fp(length).EQ.c_null_char) EXIT
  END DO
  length = length-1
  ALLOCATE (CHARACTER(len=length) :: PluginGuestReceiveF)
  DO length=1,LEN(PluginGuestReceiveF)
   PluginGuestReceiveF(length:length) = fp(length)
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
 USE Parallel_neci, ONLY : MPIBcast
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
 
 IF (PluginGuestMaster().NE.0) THEN
! ask for an FCIDUMP
  IF (.NOT. PluginGuestSendF('GIVE OPERATOR HAMILTONIAN FCIDUMP')) STOP 'plugin request has failed'
  molpro_plugin_fcidumpname = PluginGuestReceiveF()
 END IF
 CALL MPIBcast(molpro_plugin_fcidumpname)

 IF (PluginGuestMaster().NE.0) THEN
! ask for a data file
  IF (.NOT. PluginGuestSendF('GIVE INPUT NECI')) STOP 'plugin request has failed'
  molpro_plugin_datafilename = PluginGuestReceiveF()
   IF (.FALSE.) THEN
    WRITE (6, '(''Input file: '',A)') TRIM(molpro_plugin_datafilename)
    OPEN(1,file=molpro_plugin_datafilename,status='OLD')
    DO WHILE (.TRUE.)
     READ(1,'(A)',END=99) id
     WRITE (6, '(A)') TRIM(id)
    END DO
99  CLOSE(1)
   END IF
  END IF
  CALL MPIBcast(molpro_plugin_datafilename)
  
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
  IF (PluginGuestMaster().NE.0) CALL PluginGuestClose
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
 IF (PluginGuestMaster().EQ.0) RETURN
 IF (PluginGuestSendF('TAKE PROPERTY '//TRIM(property)//c_null_char)) THEN
  ALLOCATE(CHARACTER(24*size(values)) :: buffer)
  WRITE (buffer,'(1000(G23.16,1X))') values
  IF (.NOT. PluginGuestSendF(buffer//c_null_char)) STOP 'Failure to send results to plugin host'
  DEALLOCATE(buffer)
 END IF
END SUBROUTINE MolproPluginResult

END MODULE MolproPlugin
