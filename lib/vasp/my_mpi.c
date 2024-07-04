#include <string.h>
#include <time.h>
#include "my_mpi.h"


int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
{
    int ierror;
    if (sendbuf == MPI_IN_PLACE) {
        ierror = MPI_SUCCESS;
        return ierror;
    } else {
        if (sendtype == MPI_INT) {
            memcpy(recvbuf, sendbuf, sizeof(int) * sendcount);
            ierror = MPI_SUCCESS;
        } else if (sendtype == MPI_DOUBLE) {
            memcpy(recvbuf, sendbuf, sizeof(double) * sendcount);
            ierror = MPI_SUCCESS;
        } else {
            ierror = MPI_FAILURE;
        }
    }
    return ierror;
}


int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcount, int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm)
{
    int ierror;
    if (sendbuf == MPI_IN_PLACE) {
        ierror = MPI_SUCCESS;
        return ierror;
    } else {
        if (sendtype == MPI_INT) {
            memcpy(recvbuf, sendbuf, sizeof(int) * sendcount);
            ierror = MPI_SUCCESS;
        } else if (sendtype == MPI_DOUBLE) {
            memcpy(recvbuf, sendbuf, sizeof(double) * sendcount);
            ierror = MPI_SUCCESS;
        } else {
            ierror = MPI_FAILURE;
        }
    }
    return ierror;
}


int MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ierror;
    if (sendbuf == MPI_IN_PLACE) {
        ierror = MPI_SUCCESS;
        return ierror;
    } else {
        if (datatype == MPI_INT) {
            memcpy(recvbuf, sendbuf, sizeof(int) * count);
            ierror = MPI_SUCCESS;
        } else if (datatype == MPI_DOUBLE) {
            memcpy(recvbuf, sendbuf, sizeof(double) * count);
            ierror = MPI_SUCCESS;
        } else {
            ierror = MPI_FAILURE;
        }
    }
    return ierror;
}


int MPI_Barrier(MPI_Comm comm)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Bcast(void *data, int n, MPI_Datatype datatype, int node, MPI_Comm comm)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Comm_free(MPI_Comm *comm)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
    int ierror;
    ierror = MPI_SUCCESS;
    *rank = 0;
    return ierror;
}


int MPI_Comm_size(MPI_Comm comm, int *size)
{
    int ierror;
    ierror = MPI_SUCCESS;
    *size = 1;
    return ierror;
}


int MPI_Comm_split(MPI_Comm comm, int icolor, int ikey, MPI_Comm *new_comm)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Finalize(void)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ierror;
    if (sendbuf == MPI_IN_PLACE) {
        ierror = MPI_SUCCESS;
        return ierror;
    } else {
        if (sendtype == MPI_INT) {
            memcpy(recvbuf, sendbuf, sizeof(int) * sendcount);
            ierror = MPI_SUCCESS;
        } else if (sendtype == MPI_DOUBLE) {
            memcpy(recvbuf, sendbuf, sizeof(double) * sendcount);
            ierror = MPI_SUCCESS;
        } else {
            ierror = MPI_FAILURE;
        }
    }
    return ierror;
}


int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int *recvcount, int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ierror;
    if (sendbuf == MPI_IN_PLACE) {
        ierror = MPI_SUCCESS;
        return ierror;
    } else {
        if (sendtype == MPI_INT) {
            memcpy(recvbuf, sendbuf, sizeof(int) * sendcount);
            ierror = MPI_SUCCESS;
        } else if (sendtype == MPI_DOUBLE) {
            memcpy(recvbuf, sendbuf, sizeof(double) * sendcount);
            ierror = MPI_SUCCESS;
        } else {
            ierror = MPI_FAILURE;
        }
    }
    return ierror;
}


int MPI_Init(int *argc, char **argv[])
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Fetch_and_op(void *origin_addr, void *result_addr,
                     MPI_Datatype datatype, int target_rank,
                     MPI_Aint target_disp, MPI_Op op, MPI_Win win)
{
    int ierror;
    if (op == MPI_SUM) {
        ierror = MPI_SUCCESS;
        if (datatype == MPI_INT) {
            *(int *)result_addr += *(int *)origin_addr;
        } else if (datatype == MPI_DOUBLE) {
            *(double *)result_addr += *(double *)origin_addr;
        } else {
            ierror = MPI_FAILURE;
        }
    } else if (op == MPI_REPLACE) {
        ierror = MPI_SUCCESS;
        if (datatype == MPI_INT) {
            *(int *)result_addr = *(int *)origin_addr;
        } else if (datatype == MPI_DOUBLE) {
            *(double *)result_addr = *(double *)origin_addr;
        } else {
            ierror = MPI_FAILURE;
        }
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Win_create(void *baseptr, MPI_Aint size, int disp_unit,
                   MPI_Info info, MPI_Comm comm, MPI_Win *win)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Win_free(MPI_Win *win)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Win_unlock(int target_rank, MPI_Win win)
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


double MPI_Wtime()
{
    return (double)clock();
}
