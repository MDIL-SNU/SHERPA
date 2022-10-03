#ifndef __DIMER_H__
#define __DIMER_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


int dimer(Config *initial, Config *final, Input *input, double *full_eigenmode,
          int count, int index, double *Ea, MPI_Comm comm);
#endif
