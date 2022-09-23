#ifndef __ART_NOUVEAU__
#define __ART_NOUVEAU__
#include <mpi.h>
#include "config.h"
#include "dataset.h"
#include "input.h"


int lanczos(Config *, Input *, int, int *, double *, double *, MPI_Comm);
void uphill_push(Config *, Input *, int, int *, double *, double, double *, MPI_Comm);
void normal_relax(Config *, Input *, int, int *, double, double *, MPI_Comm);
int art_nouveau(Config *, Config *, Input *, Data *, int, int, double *, MPI_Comm);

#endif
