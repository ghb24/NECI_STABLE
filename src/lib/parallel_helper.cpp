#ifdef CBINDMPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unistd.h>

#ifdef CBINDMPI

int g_argc;
char** g_argv;
#ifndef MOLPRO
//
// We need an entry point which can be found by a c++ based linker
extern "C" void neci_main_c ();

int main (int argc, char ** argv)
{
	g_argc = argc;
	g_argv = argv;
	neci_main_c ();
}
#endif

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
// A comm list.
// This is initialised so that the first value is MPI_COMM_WORLD (i.e. the
// default communicator).
typedef std::vector<MPI_Comm> comm_vec_t;
comm_vec_t new_comm_vec () {
	comm_vec_t v;
	v.push_back(MPI_COMM_WORLD);
	return v;
}
comm_vec_t comm_vec = new_comm_vec();

//
// A group list
typedef std::vector<MPI_Group> group_vec_t;
group_vec_t group_vec;


//
// Wrapper for MPI_INIT
void mpi_init_wrap (int * ierr)
{
    *ierr = MPI_Init (NULL, NULL);
}


//
// Wrapper for MPI_Finalize
void mpi_finalize_wrap (int * ierr)
{
    *ierr = MPI_Finalize ();
}


//
// Wrapper for MPI_Abort
void mpi_abort_wrap (int comm, int err, int * ierr)
{
    *ierr = MPI_Abort (comm_vec[comm], err);
}


//
// Wrapper for MPI_Barrier
void mpi_barrier_wrap (int comm, int * ierr)
{
    *ierr = MPI_Barrier (comm_vec[comm]);
}


//
// Wrapper for MPI_Comm_rank
void mpi_comm_rank_wrap (int comm, int * rank, int * ierr)
{
    *ierr = MPI_Comm_rank (comm_vec[comm], rank);
}


//
// Wrapper for MPI_Comm_size
void mpi_comm_size_wrap (int comm, int * size, int * ierr)
{
    *ierr = MPI_Comm_size (comm_vec[comm], size);
}


//
// Wrapper for MPI_Comm_group
void mpi_comm_group_wrap (int comm, int * group, int * ierr)
{
	// Get the comm handle
	MPI_Comm comm_handle = comm_vec[comm];

	// Call the MPI routine
	MPI_Group grp_handle;
	*ierr = MPI_Comm_group (comm_handle, &grp_handle);

	// We need to return an integer, so store group in a vector and return
	// its index. See if it is already in there first...
	// n.b. Don't use an iterator, as we actually want the index...
	for (int i = 0; i < group_vec.size(); ++i) {
		if (group_vec[i] == grp_handle) {
			*group = i;
			return;
		}
	}

	// Not in the list, so we need to append it
	group_vec.push_back(grp_handle);
	*group = group_vec.size() - 1;
}


//
// Wrapper for MPI_Comm_create
void mpi_comm_create_wrap (int comm, int group, int * ncomm, int * ierr)
{
	MPI_Comm comm_handle = comm_vec[comm];
	MPI_Group grp_handle = group_vec[group];
	MPI_Comm new_comm;

	*ierr = MPI_Comm_create (comm_handle, grp_handle, &new_comm);

	// Add this comm to the list, and return its index
	comm_vec.push_back(new_comm);
	*ncomm = comm_vec.size() - 1;
}


//
// Wrapper for MPI_Group_incl
void mpi_group_incl_wrap (int group, int n, int * ranks, int * ogroup,
                          int * ierr)
{
	MPI_Group grp_handle = group_vec[group];
	MPI_Group new_group;

	*ierr = MPI_Group_incl (grp_handle, n, ranks, &new_group);

	// Add this group to the list, and return its index
	group_vec.push_back (new_group);
	*ogroup = group_vec.size() - 1;
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
                      int op, int root, int comm, int * ierr)
{

    *ierr = MPI_Reduce (sbuf ? sbuf : MPI_IN_PLACE,
                        rbuf ? rbuf : MPI_IN_PLACE, count, dtype_map[dtype],
                        op_map[op], root, comm_vec[comm]);
}


//
// Wrapper for MPI_Allreduce
void mpi_allreduce_wrap (double * sbuf, double * rbuf, int count, int dtype,
                         int op, int comm, int * ierr)
{
    *ierr = MPI_Allreduce (sbuf ? sbuf : MPI_IN_PLACE,
                           rbuf ? rbuf : MPI_IN_PLACE, count,
                           dtype_map[dtype], op_map[op], comm_vec[comm]);
}


//
// Wrapper for MPI_Bcast
void mpi_bcast_wrap (void * buf, int count, int dtype, int root, int comm,
                     int * ierr)
{
    *ierr = MPI_Bcast (buf, count, dtype_map[dtype], root, comm_vec[comm]);
}


//
// Wrapper for MPI_Alltoall
void mpi_alltoall_wrap (void * sbuf, int scount, int stype, void * rbuf,
                        int rcount, int rtype, int comm, int * ierr)
{
    *ierr = MPI_Alltoall (sbuf, scount, dtype_map[stype], rbuf, rcount,
                          dtype_map[rtype], comm_vec[comm]);
}


//
// Wrapper for MPI_AlltoallV
void mpi_alltoallv_wrap (void * sbuf, int * scount, int * sdispl, int stype,
                         void * rbuf, int * rcount, int * rdispl, int rtype,
                         int comm, int * ierr)
{
    *ierr = MPI_Alltoallv (sbuf, scount, sdispl, dtype_map[stype], rbuf,
                           rcount, rdispl, dtype_map[rtype], comm_vec[comm]);
}


//
// Wrapper for MPI_Gather
void mpi_gather_wrap (void * sbuf, int scount, int stype, void * rbuf,
                      int rcount, int rtype, int root, int comm, int * ierr)
{
    *ierr = MPI_Gather (sbuf, scount, dtype_map[stype], rbuf, rcount,
                        dtype_map[rtype], root, comm_vec[comm]);
}


//
// Wrapper for MPI_GatherV
void mpi_gatherv_wrap (void * sbuf, int scount, int stype, void * rbuf,
                      int * rcount, int * displs, int rtype, int root,
                      int comm, int * ierr)
{
    *ierr = MPI_Gatherv (sbuf, scount, dtype_map[stype], rbuf, rcount,
                         displs, dtype_map[rtype], root, comm_vec[comm]);
}


// Wrapper for MPI_Allgather
void mpi_allgather_wrap (void * sbuf, int scount, int stype, void * rbuf,
                         int rcount, int rtype, int comm, int * ierr)
{
    *ierr = MPI_Allgather (sbuf, scount, dtype_map[stype], rbuf, rcount,
                           dtype_map[rtype], comm_vec[comm]);
}


//
// Wrapper for MPI_ScatterV
void mpi_scatterv_wrap (void * sbuf, int * scount, int * displs, int stype,
                        void * rbuf, int rcount, int rtype, int root,
                        int comm, int * ierr)
{
    *ierr = MPI_Scatterv (sbuf, scount, displs, dtype_map[stype], rbuf,
                          rcount, dtype_map[rtype], root, comm_vec[comm]);
}


//
// Wrapper for MPI_Scatter
void mpi_scatter_wrap (void * sbuf, int scount, int stype, void * rbuf,
		               int rcount, int rtype, int root, int comm, int * ierr)
{
    *ierr = MPI_Scatter (sbuf, scount, dtype_map[stype], rbuf,
                          rcount, dtype_map[rtype], root, comm_vec[comm]);
}


//
// Wrapper for MPI_Send
void mpi_send_wrap (void * buf, int count, int dtype, int dest, int tag,
                    int comm, int * ierr)
{
    *ierr = MPI_Send (buf, count, dtype_map[dtype], dest, tag,
                      comm_vec[comm]);
}


//
// Wrapper for MPI_Recv
void mpi_recv_wrap (void * buf, int count, int dtype, int src, int tag,
                    int comm, int * stat_ignore, int * ierr)
{
    MPI_Status stat;
    *ierr = MPI_Recv (buf, count, dtype_map[dtype], src, tag,
                      comm_vec[comm], &stat); // comm
}




// End C-linkage
}

#endif // CBINDMPI
