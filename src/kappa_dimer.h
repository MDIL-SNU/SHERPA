#ifndef __KAPPA_DIMER_H__
#define __KAPPA_DIMER_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


int kappa_dimer(Config *initial, Config *saddle, Config *final, Input *input,
                double *full_eigenmode, int count, int index, double *Ea,
                MPI_Comm comm);
#endif
