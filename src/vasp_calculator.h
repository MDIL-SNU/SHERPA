#ifndef __VASP_CALCULATOR_H__
#define __VASP_CALCULATOR_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"

void oneshot(Config *, Input *, double *, MPI_Comm);
void oneshot_disp(Config *, Input *, double *, double *, int, int *, MPI_Comm);
void atom_relax(Config *, Input *, double *, MPI_Comm);
#endif
