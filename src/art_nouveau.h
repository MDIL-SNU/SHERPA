#ifndef __ART_NOUVEAU__
#define __ART_NOUVEAU__
#include "config.h"
#include "dataset.h"
#include "input.h"
#include "my_mpi.h"


void lanczos(Config *, Input *, int, int *, double *, double *, MPI_Comm);
double *uphill_push(Config *, Input *, int, int *, double, double *, double *, MPI_Comm);
void perp_relax(Config *, Input *, int, int *, double, double *, int, int, MPI_Comm);
int art_nouveau(Config *, Config *, Input *, Data *, int, int, double *, MPI_Comm);

#endif
