#ifndef __LMP_CALCULATOR_H__
#define __LMP_CALCULATOR_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"

void *lmp_init(Config *, Input *, int, char **, MPI_Comm);
double oneshot(Config *, Input *, MPI_Comm);
void oneshot_disp(Config *, Input *, double *, double *, int, int *, MPI_Comm);
double atom_relax(Config *, Input *, MPI_Comm);
#endif
