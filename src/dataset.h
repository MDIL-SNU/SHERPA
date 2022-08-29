#ifndef __DATASET_H__
#define __DATASET_H__
#include "config.h"


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
void build_dataset(Dataset *, Config *, Input *, long long);
void free_dataset(Dataset *);
#endif
