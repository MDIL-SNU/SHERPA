#ifndef __MY_MPI_H__
#define __MY_MPI_H__

#ifdef LMP
#include <mpi.h>
#endif
#ifdef VASP
#include <stdio.h>
#include <stdlib.h>

#define MPI_COMM_WORLD 0
#define MPI_SUCCESS 0
#define MPI_FAILURE 1

#define MPI_Comm int
#define MPI_Datatype int
#define MPI_Op int

#define MPI_IN_PLACE NULL

#define MPI_DATATYPE_NULL 0
#define MPI_INT 1 
#define MPI_DOUBLE 2

#define MPI_SUM 1

int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm);
int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcount, int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Barrier(MPI_Comm comm);
int MPI_Bcast(void *data, int n, MPI_Datatype datatype, int node, MPI_Comm comm);
int MPI_Comm_free(MPI_Comm *comm);
int MPI_Comm_rank(MPI_Comm comm, int *rank);
int MPI_Comm_size(MPI_Comm comm, int *size);
int MPI_Comm_split(MPI_Comm comm, int icolor, int ikey, MPI_Comm *new_comm);
int MPI_Finalize(void);
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount,
               MPI_Datatype datatype, int root, MPI_Comm comm);
int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int *recvcount, int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Init(int *argc, char **argv[]);


#define MPI_Win int
#define MPI_Aint int
#define MPI_Info int
#define MPI_INFO_NULL 0
#define MPI_LOCK_EXCLUSIVE 1

int MPI_Fetch_and_op(void *origin_addr, void *result_addr,
                     MPI_Datatype datatype, int target_rank,
                     MPI_Aint target_disp, MPI_Op op, MPI_Win win);
int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                     MPI_Comm comm, void *baseptr, MPI_Win *win);
int MPI_Win_free(MPI_Win *win);
int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
int MPI_Win_unlock(int target_rank, MPI_Win win);

#endif
#endif
