#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include <mpi.h>
#include "config.h"
#include "input.h"

void *lmp_init(Config *, Input *, int, char **, MPI_Comm);
double oneshot(Config *, Input *, MPI_Comm);
void oneshot_disp(Config *, Input *, double *, double *, int, int *, MPI_Comm);
double atom_relax(Config *, Input *, MPI_Comm);
#endif
