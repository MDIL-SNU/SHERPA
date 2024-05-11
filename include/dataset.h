#ifndef __DATASET_H__
#define __DATASET_H__


typedef struct _Data
{
    struct _Data *next;
    int index;
    double *eigenmode;
} Data;

typedef struct _Dataset
{
    struct _Data *head;
    struct _Data *tail;
} Dataset;

void insert_data(Dataset *dataset, int n, int index, double *eigenmode);
void build_dataset(Dataset *dataset, int n);
void free_dataset(Dataset *dataset);
#endif
