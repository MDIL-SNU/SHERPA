#ifndef __DIMER_H__
#define __DIMER_H__
#include "calculator.h"
#include "config.h"
#include "input.h"
#include "my_mpi.h"


int dimer(Calc *calc, Config *initial, Config *saddle, Config *final,
          Input *input, double *full_eigenmode, int count, int index,
          double *Ea, MPI_Comm comm);
#endif
