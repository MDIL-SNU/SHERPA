#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include <mpi.h>
#include "config.h"
#include "input.h"

void *lmp_init(Config *, Input *, int, char **, MPI_Comm);
double global_oneshot(Config *, Input *, MPI_Comm);
double local_oneshot(Config *, Input *, int, MPI_Comm);
double atom_relax(Config *, Input *, MPI_Comm);
#endif
