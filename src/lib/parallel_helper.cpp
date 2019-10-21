#ifdef _MOLCAS_
#include "molcas_wrapper.h"
#endif

#ifdef CBINDMPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif
#include <stdio.h>
#include <string.h>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#else
#include <unistd.h>
#endif

#ifdef CBINDMPI

int g_argc;
char** g_argv;
//
// We need an entry point which can be found by a c++ based linker
extern "C" void neci_main_c ();

int main (int argc, char ** argv)
{
	g_argc = argc;
	g_argv = argv;
	neci_main_c ();
}

// We need C-linkage
extern "C" {

//
// The next three functions return information about the command line
// arguments to fortran, so that we can run properly when compiled
// with CBINDMPI.
//
int c_argc ()
{
	return g_argc;
}

int c_getarg_len (int i)
{
	return strlen(g_argv[i]);
}

void c_getarg (int i, char * str)
{
	strcpy(str, g_argv[i]);
}

MPI_Fint mpicommworld_c2f ()
{
	return MPI_Comm_c2f(MPI_COMM_WORLD);
}


//
// The MPI specification has several constants and specific datatypes
// which are not available to fortran
//
// --> We need to create a mapping between constants in fortran and C.
//     See ParallelHelper.F90 for the fortran equivalents.
//
// --> We need to store generated Comms and Groups within the C code and
//     pass an index to those lists between fortran and C.
//

//
// Match the datatype labels between fortran and C
const MPI_Datatype dtype_map[] = {MPI_INTEGER4,
                                  MPI_INTEGER8,
                                  MPI_DOUBLE_PRECISION,
                                  MPI_DOUBLE_COMPLEX,
                                  MPI_2INTEGER,
                                  MPI_INTEGER,
                                  MPI_CHARACTER,
                                  MPI_2DOUBLE_PRECISION};

//
// Match the MPI operations between fortran and C
const MPI_Op op_map[] = {MPI_SUM,
                         MPI_MAX,
                         MPI_MIN,
                         MPI_LOR,
                         MPI_MINLOC,
                         MPI_MAXLOC};


//
// Wrapper for MPI_INIT
void mpi_init_wrap (int * ierr)
{
    *ierr = MPI_Init (NULL, NULL);
}

void mpi_initialized_wrap(int * flag, int * ierr)
{
  *ierr = MPI_Initialized(flag);
}


//
// Wrapper for MPI_Finalize
void mpi_finalize_wrap (int *ierr)
{
    *ierr = MPI_Finalize();
}


//
// Wrapper for MPI_Abort
void mpi_abort_wrap (MPI_Fint comm, int err, int * ierr)
{
    *ierr = MPI_Abort (MPI_Comm_f2c(comm), err);
}


//
// Wrapper for MPI_Barrier
void mpi_barrier_wrap (MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Barrier (MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Comm_rank
void mpi_comm_rank_wrap (MPI_Fint comm, int * rank, int * ierr)
{
    *ierr = MPI_Comm_rank (MPI_Comm_f2c(comm), rank);
}


//
// Wrapper for MPI_Comm_size
void mpi_comm_size_wrap (MPI_Fint comm, int * size, int * ierr)
{
    *ierr = MPI_Comm_size (MPI_Comm_f2c(comm), size);
}


//
// Wrapper for MPI_Comm_group
void mpi_comm_group_wrap (MPI_Fint comm, MPI_Fint * group, int * ierr)
{
	// Get the comm handle
	MPI_Comm comm_handle = MPI_Comm_f2c(comm);

	// Call the MPI routine
	MPI_Group grp_handle;
	*ierr = MPI_Comm_group (comm_handle, &grp_handle);

	// And convert into an integer we can push back to fortran
	*group = MPI_Group_c2f(grp_handle);
}


//
// Wrapper for MPI_Comm_create
void mpi_comm_create_wrap (MPI_Fint comm, MPI_Fint group, MPI_Fint * ncomm, int * ierr)
{
	MPI_Comm comm_handle = MPI_Comm_f2c(comm);
	MPI_Group grp_handle = MPI_Group_f2c(group);
	MPI_Comm new_comm;

	*ierr = MPI_Comm_create (comm_handle, grp_handle, &new_comm);

	// And convert to an integer to push back to fortran
	*ncomm = MPI_Comm_c2f(new_comm);
}


//
// Wrapper for MPI_Group_incl
void mpi_group_incl_wrap (MPI_Fint group, int n, int * ranks, 
		                  MPI_Fint * ogroup, int * ierr)
{
	MPI_Group grp_handle = MPI_Group_f2c(group);
	MPI_Group new_group;

	*ierr = MPI_Group_incl (grp_handle, n, ranks, &new_group);

	// And convert to an integer fortran will understand
	*ogroup = MPI_Group_c2f(new_group);
}


//
// Wrapper for MPI_Error_string
void mpi_error_string_wrap (int err, char * str, int * len, int * ierr)
{
    *ierr = MPI_Error_string (err, str, len);
}


//
// Wrapper for MPI_Reduce
void mpi_reduce_wrap (void * sbuf, void * rbuf, int count, int dtype,
                      int op, int root, MPI_Fint comm, int * ierr)
{

    *ierr = MPI_Reduce (sbuf ? sbuf : MPI_IN_PLACE,
                        rbuf ? rbuf : MPI_IN_PLACE, count, dtype_map[dtype],
                        op_map[op], root, MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Allreduce
void mpi_allreduce_wrap (double * sbuf, double * rbuf, int count, int dtype,
                         int op, MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Allreduce (sbuf ? sbuf : MPI_IN_PLACE,
                           rbuf ? rbuf : MPI_IN_PLACE, count,
                           dtype_map[dtype], op_map[op], MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Bcast
void mpi_bcast_wrap (void * buf, int count, int dtype, int root, 
                     MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Bcast (buf, count, dtype_map[dtype], root, 
	                   MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Alltoall
void mpi_alltoall_wrap (void * sbuf, int scount, int stype, void * rbuf,
                        int rcount, int rtype, MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Alltoall (sbuf, scount, dtype_map[stype], rbuf, rcount,
                          dtype_map[rtype], MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_AlltoallV
void mpi_alltoallv_wrap (void * sbuf, int * scount, int * sdispl, int stype,
                         void * rbuf, int * rcount, int * rdispl, int rtype,
                         MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Alltoallv (sbuf, scount, sdispl, dtype_map[stype], rbuf,
                           rcount, rdispl, dtype_map[rtype], 
						   MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Gather
void mpi_gather_wrap (void * sbuf, int scount, int stype, void * rbuf,
                      int rcount, int rtype, int root, MPI_Fint comm, 
					  int * ierr)
{
    *ierr = MPI_Gather (sbuf, scount, dtype_map[stype], rbuf, rcount,
                        dtype_map[rtype], root, MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_GatherV
void mpi_gatherv_wrap (void * sbuf, int scount, int stype, void * rbuf,
                      int * rcount, int * displs, int rtype, int root,
                      MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Gatherv (sbuf, scount, dtype_map[stype], rbuf, rcount,
                         displs, dtype_map[rtype], root, MPI_Comm_f2c(comm));
}


// Wrapper for MPI_Allgather
void mpi_allgather_wrap (void * sbuf, int scount, int stype, void * rbuf,
                         int rcount, int rtype, MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Allgather (sbuf, scount, dtype_map[stype], rbuf, rcount,
                           dtype_map[rtype], MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_ScatterV
void mpi_scatterv_wrap (void * sbuf, int * scount, int * displs, int stype,
                        void * rbuf, int rcount, int rtype, int root,
                        MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Scatterv (sbuf, scount, displs, dtype_map[stype], rbuf,
                          rcount, dtype_map[rtype], root, MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_AllGatherV
void mpi_allgatherv_wrap (void * sbuf, int scount, int stype, void * rbuf,
		                  int * rcount, int * displs, int rtype, MPI_Fint comm,
						  int * ierr)
{
	*ierr = MPI_Allgatherv (sbuf, scount, dtype_map[stype], rbuf, rcount,
                            displs, dtype_map[rtype], MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Scatter
void mpi_scatter_wrap (void * sbuf, int scount, int stype, void * rbuf,
		               int rcount, int rtype, int root, MPI_Fint comm,
					   int * ierr)
{
    *ierr = MPI_Scatter (sbuf, scount, dtype_map[stype], rbuf,
                         rcount, dtype_map[rtype], root, MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Send
void mpi_send_wrap (void * buf, int count, int dtype, int dest, int tag,
                    MPI_Fint comm, int * ierr)
{
    *ierr = MPI_Send (buf, count, dtype_map[dtype], dest, tag,
                      MPI_Comm_f2c(comm));
}


//
// Wrapper for MPI_Recv
void mpi_recv_wrap (void * buf, int count, int dtype, int src, int tag,
                    MPI_Fint comm, int * stat_ignore, int * ierr)
{
    MPI_Status stat;
    *ierr = MPI_Recv (buf, count, dtype_map[dtype], src, tag,
                      MPI_Comm_f2c(comm), &stat);
}




// End C-linkage
}

#endif // CBINDMPI
