#ifndef __SNC_DIMER_H__
#define __SNC_DIMER_H__
#include "config.h"
#include "dimer.h"
#include "input.h"
#include "my_mpi.h"


int snc_dimer(Config *, Config *, Input *, Data *, int, int, double *, MPI_Comm);
#endif
