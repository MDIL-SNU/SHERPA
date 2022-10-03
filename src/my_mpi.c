#ifdef LMP
#include <mpi.h>
#endif
#ifdef VASP
#include "my_mpi.h"


int MPI_Copy_int(int *data1, int *data2, int n)
{
    int i, ierror;
    for (i = 0; i < n; ++i) {
        data2[i] = data1[i];
    }
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Copy_double(double *data1, double *data2, int n)
{
    int i, ierror;
    for (i = 0; i < n; ++i) {
        data2[i] = data1[i];
    }
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Reduce_int(int *data1, int *data2, int n, MPI_Op op)
{
    int i, ierror;
    if (op == MPI_SUM) {
        for (i = 0; i < n; ++i) {
            data2[i] = data1[i];
        }
        ierror = MPI_SUCCESS;
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Reduce_double(double *data1, double *data2, int n, MPI_Op op)
{
    int i, ierror;
    if (op == MPI_SUM) {
        for (i = 0; i < n; ++i) {
            data2[i] = data1[i];
        }
        ierror = MPI_SUCCESS;
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
{
    int ierror;
    if (sendtype == MPI_INT) {
        if (sendbuf != NULL) {
            ierror = MPI_Copy_int((int *)sendbuf, (int *)recvbuf, sendcount);
        }
    } else if (sendtype == MPI_DOUBLE) {
        if (sendbuf != NULL) {
            ierror = MPI_Copy_double((double *)sendbuf, (double *)recvbuf, sendcount);
        }
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcount, int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm)
{
    int ierror;
    if (sendtype == MPI_INT) {
        if (sendbuf != NULL) {
            ierror = MPI_Copy_int((int *)sendbuf, (int *)recvbuf, sendcount);
        }
    } else if (sendtype == MPI_DOUBLE) {
        if (sendbuf != NULL) {
            ierror = MPI_Copy_double((double *)sendbuf, (double *)recvbuf, sendcount);
        }
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ierror;
    if (datatype == MPI_INT) {
        ierror = MPI_Reduce_int((int *)sendbuf, (int *)recvbuf, count, op);
    } else if (datatype == MPI_DOUBLE) {
        ierror = MPI_Reduce_double((double *)sendbuf, (double *)recvbuf, count, op);
    } else {
        ierror = MPI_FAILURE;
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


int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, const int recvcount,
               MPI_Datatype datatype, int root, MPI_Comm comm)
{
    int ierror;
    if (datatype == MPI_INT) {
        ierror = MPI_Copy_int((int *)sendbuf, (int *)recvbuf, sendcount);
    } else if (datatype == MPI_DOUBLE) {
        ierror = MPI_Copy_double((double *)sendbuf, (double *)recvbuf, sendcount);
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int *recvcount, const int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ierror;
    if (sendtype == MPI_INT) {
        ierror = MPI_Copy_int((int *)sendbuf, (int *)recvbuf, sendcount);
    } else if (sendtype == MPI_DOUBLE) {
        ierror = MPI_Copy_double((double *)sendbuf, (double *)recvbuf, sendcount);
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Init(int *argc, char **argv[])
{
    int ierror;
    ierror = MPI_SUCCESS;
    return ierror;
}


int MPI_Fetch_and_op(const void *origin_addr, void *result_addr,
                     MPI_Datatype datatype, int target_rank,
                     MPI_Aint target_disp, MPI_Op op, MPI_Win win)
{
    int ierror;
    if (op == MPI_SUM) {
        ierror = MPI_SUCCESS;
        if (datatype == MPI_INT) {
            int sum = *((int *)result_addr) + *((int *)origin_addr);
            result_addr = &sum;
        } else if (datatype == MPI_DOUBLE) {
            double sum = *((double *)result_addr) + *((double *)origin_addr);
            result_addr = &sum;
        } else {
            ierror = MPI_FAILURE;
        }
    } else {
        ierror = MPI_FAILURE;
    }
    return ierror;
}


int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                     MPI_Comm comm, void *baseptr, MPI_Win *win)
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
#endif
