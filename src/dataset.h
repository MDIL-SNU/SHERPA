#ifndef __DATASET_H__
#define __DATASET_H__
#include <mpi.h>
#include "config.h"
#include "input.h"


typedef struct _Data
{
    struct _Data *next;

    int index;
    int *type;
    double *saddle;
    double *eigenmode;
} Data;

typedef struct _Dataset
{
    int numdata;
    struct _Data *head;
} Dataset;

void insert_data(Dataset *, int, int, int *, double *, double *);
void recycle_data(Config *, Config *, Input *, Data *, Config *, double *, MPI_Comm);
void build_dataset(Dataset *, Input *, int, int *, int);
void free_dataset(Dataset *);
#endif
