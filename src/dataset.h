#ifndef __DATASET_H__
#define __DATASET_H__
#include "config.h"
#include "input.h"


typedef struct _Data
{
    struct _Data *next;
    int index;
    double *eigenmode;
} Data;

typedef struct _Dataset
{
    struct _Data *head;
} Dataset;

void insert_data(Dataset *, int, int, double *);
void build_dataset(Dataset *, Input *, int);
void free_dataset(Dataset *);
#endif
