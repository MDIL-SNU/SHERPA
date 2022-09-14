#ifndef __KAPPA_DIMER_H__
#define __KAPPA_DIMER_H__
#include <mpi.h>
#include "config.h"
#include "dimer.h"
#include "input.h"


int kappa_dimer(Config *, Config *, Config *, Input *, double *, int, int, double *, MPI_Comm);
#endif
