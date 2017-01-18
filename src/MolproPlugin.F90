MODULE MolproPlugin
#include "mpif.h"
 INTEGER :: molpro_plugin_intercomm
 INTEGER :: molpro_plugin
 CHARACTER(1024) :: molpro_plugin_fcidumpname, molpro_plugin_datafilename
 LOGICAL :: molpro_plugin_redirect_output = .TRUE. ! normally .true., maybe .false. for debugging crashes
CONTAINS
SUBROUTINE MolproPluginInit (tMolpro)
 USE iso_c_binding, ONLY : c_null_char
 IMPLICIT NONE
 LOGICAL, INTENT(inout) :: tMolpro
 INTEGER :: size
 INTEGER :: rank
 INTEGER :: ierr
 CHARACTER(4) :: molpro_version
 CHARACTER(1024) :: id
 INTEGER length
 INTEGER, DIMENSION(MPI_Status_size) :: status
 CALL MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
 CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
 ! is this a Molpro plugin?
  molpro_plugin=0
  CALL MPI_Comm_get_parent(molpro_plugin_intercomm, ierr)
  IF (rank.EQ.0 .AND. molpro_plugin_intercomm .NE. MPI_COMM_NULL) THEN
! expect plugin server to identify itself
   CALL MPI_Recv(length,1,MPI_INT,0,0,molpro_plugin_intercomm,status,ierr)
   CALL MPI_Recv(id,length,MPI_CHAR,0,1,molpro_plugin_intercomm,status,ierr)
   id(length:)=' '
   !PRINT '(''Plugin server: '',A)',trim(id)
   IF (id(1:6).EQ.'MOLPRO') molpro_plugin=1
   IF (molpro_plugin.GT.0) THEN
    molpro_version=id(8:11)
   END IF
  END IF
  CALL MPI_Bcast(molpro_plugin,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

  IF (molpro_plugin.EQ.0) THEN
   CONTINUE ! IF (rank.EQ.0) PRINT '(''Not running in plugin mode'')' ! or stay silent
  ELSE
   if (molpro_plugin_redirect_output) CLOSE(6)
   tMolpro = .TRUE.
   CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   IF (rank.EQ.0) THEN ! communication should be handled by just the root process,
    if (molpro_plugin_redirect_output) OPEN(6,file=TRIM(id(13:)))
! exchanging messages with the root process on the server
    WRITE (6, '(''Plugin for Molpro version '',A,'' running on '',I3,'' processors'')') molpro_version,size
   END IF
  END IF

  IF (molpro_plugin.GT.0 .AND. rank.EQ.0) THEN
! ask for an FCIDUMP
    id="GIVE OPERATOR HAMILTONIAN FCIDUMP"//c_null_char
    length=len_trim(id)
    CALL MPI_Send(length, 1,MPI_INT,0,0,molpro_plugin_intercomm,ierr)
    CALL MPI_Send(id,length, MPI_CHAR,0,1,molpro_plugin_intercomm,ierr)
    CALL MPI_Recv(length,1,MPI_INT,0,0,molpro_plugin_intercomm,status,ierr)
    IF (length.EQ.0) STOP 'plugin request has failed'
    CALL MPI_Recv(molpro_plugin_fcidumpname,length,MPI_CHAR,0,1,molpro_plugin_intercomm,status,ierr)
    molpro_plugin_fcidumpname(length+1:)=' '
   END IF

  IF (molpro_plugin.GT.0 .AND. rank.EQ.0) THEN
! ask for a data file
    id="GIVE INPUT NECI"//c_null_char
    length=len_trim(id)
    CALL MPI_Send(length, 1,MPI_INT,0,0,molpro_plugin_intercomm,ierr)
    CALL MPI_Send(id,length, MPI_CHAR,0,1,molpro_plugin_intercomm,ierr)
    CALL MPI_Recv(length,1,MPI_INT,0,0,molpro_plugin_intercomm,status,ierr)
    IF (length.EQ.0) STOP 'plugin request has failed'
    CALL MPI_Recv(molpro_plugin_datafilename,length,MPI_CHAR,0,1,molpro_plugin_intercomm,status,ierr)
    molpro_plugin_datafilename(length+1:)=' '
    IF (.FALSE.) THEN
     WRITE (6, '(''Input file: '',A)') TRIM(molpro_plugin_datafilename)
     OPEN(1,file=molpro_plugin_datafilename,status='OLD')
     DO WHILE (.TRUE.)
      READ(1,'(A)',END=99) id
      WRITE (6, '(A)') TRIM(id)
     END DO
99   CLOSE(1)
    END IF
   END IF

END SUBROUTINE MolproPluginInit
END MODULE MolproPlugin
