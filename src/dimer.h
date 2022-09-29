#ifndef __DIMER_H__
#define __DIMER_H__
#include <mpi.h>
#include "config.h"
#include "dataset.h"
#include "input.h"


void rotate_vector(double *, double *, double **, double **, int, double);
double *get_rot_force(Input *, double *, double *, double *, int);
int dimer(Config *, Config *, Input *, Data *, int, int, double *, MPI_Comm);
#endif
