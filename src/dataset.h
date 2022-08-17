#ifndef __DATASET_H__
#define __DATASET_H__
#include "config.h"


typedef struct _Data
{
    struct _Data *next;

    int natom;
    int *type;
    int index;

    double *initial;
    double *saddle;
    double *eigenmode;
} Data;

typedef struct _Dataset
{
    int numdata;
    struct _Data *head;
} Dataset;

void insert_data(Dataset *, int, int *, int, int *, int *,
                 double *, double *, double *);
Data *search_data(Dataset *, Config *, long long, int *, int, int);
void delete_data(Dataset *);
void free_dataset(Dataset *);
#endif
