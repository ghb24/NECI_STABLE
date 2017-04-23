MODULE MolproPlugin
 USE PluginGuestF
 IMPLICIT NONE
#include "mpif.h"
 INTEGER :: molpro_plugin
 CHARACTER(1024) :: molpro_plugin_fcidumpname, molpro_plugin_datafilename
CONTAINS
SUBROUTINE MolproPluginInit (tMolpro)
 IMPLICIT NONE
 LOGICAL, INTENT(inout) :: tMolpro
 INTEGER :: rank
 INTEGER :: ierr
 CHARACTER(1024) :: id
 CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
 ! is this a Molpro plugin?
 CALL PluginGuestF_open('MOLPRO')
 tMolpro = PluginGuestF_active()
 molpro_plugin=0
 IF (.NOT. tMolpro) RETURN
 molpro_plugin=1
 
 IF (rank.EQ.0) THEN
! ask for an FCIDUMP
  IF (.NOT. PluginGuestF_send('GIVE OPERATOR HAMILTONIAN FCIDUMP')) STOP 'plugin request has failed'
  molpro_plugin_fcidumpname = PluginGuestF_receive()
 END IF
 CALL MPI_Bcast(molpro_plugin_fcidumpname,LEN(molpro_plugin_fcidumpname),MPI_Char,0,MPI_COMM_WORLD,ierr)

  IF (rank.EQ.0) THEN
! ask for a data file
   IF (.NOT. PluginGuestF_send('GIVE INPUT NECI')) STOP 'plugin request has failed'
   molpro_plugin_datafilename = PluginGuestF_receive()
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
  CALL MPI_Bcast(molpro_plugin_datafilename,LEN(molpro_plugin_datafilename),MPI_Char,0,MPI_COMM_WORLD,ierr)
  
 END SUBROUTINE MolproPluginInit

SUBROUTINE MolproPluginTerm(signal)
 USE iso_c_binding, ONLY : c_int
 INTEGER, INTENT(in) :: signal
 INTEGER :: ierr
 INTERFACE
  SUBROUTINE fsleep(seconds) BIND(C,name="sleep")
   IMPORT
   INTEGER(c_int), VALUE :: seconds
  END SUBROUTINE fsleep
 END INTERFACE
 INTEGER :: rank
 IF (PluginGuestF_active()) THEN
! Graceful exit if Molpro server
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  IF (rank.EQ.0) CALL PluginGuestF_close
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
 INTEGER :: rank, ierr
 IF (.NOT. PluginGuestF_active()) RETURN
 CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
 IF (rank.EQ.0) RETURN
 IF (PluginGuestF_send('TAKE PROPERTY '//TRIM(property))) THEN
  ALLOCATE(CHARACTER(24*size(values)) :: buffer)
  WRITE (buffer,'(1000(G23.16,1X))') values
  IF (.NOT. PluginGuestF_send(buffer)) stop 'Failure to send results to plugin host'
 END IF
END SUBROUTINE MolproPluginResult
END MODULE MolproPlugin
