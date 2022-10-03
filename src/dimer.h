#ifndef __DIMER_H__
#define __DIMER_H__
#include "config.h"
#include "dataset.h"
#include "input.h"
#include "my_mpi.h"


void rotate_vector(double *, double *, double **, double **, int, double);
double *get_rot_force(Input *, double *, double *, double *, int);
int dimer(Config *, Config *, Input *, Data *, int, int, double *, MPI_Comm);
#endif
