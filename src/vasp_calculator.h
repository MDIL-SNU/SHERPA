#ifndef __VASP_CALCULATOR_H__
#define __VASP_CALCULATOR_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


void oneshot(Config *config, Input *input, double *energy, double *force,
             MPI_Comm comm);
void atom_relax(Config *config, Input *input, double *energy, MPI_Comm);
#endif
