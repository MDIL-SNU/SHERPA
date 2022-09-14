#ifndef __DIMER_H__
#define __DIMER_H__
#include <mpi.h>
#include "config.h"
#include "input.h"


double normal_random(double, double);
double norm(double *, int);
double *normalize(double *, int);
double dot(double *, double *, int);
double *parallel_vector(double *, double *, int);
double *perpendicular_vector(double *, double *, int);
void rotate_vector(double *, double *, double **, double **, int, double);
void cut_sphere(Config *, Input *, int, int *);
double *gen_eigenmode(Input *, int, MPI_Comm);
void gen_list(Config *, Input *, double *, int *, int **, int *, int **, MPI_Comm);
double *get_rot_force(Input *, double *, double *, double *, int);
void get_cg_direction(double *, double *, double *, int);
int dimer(Config *, Config *, Config *, Input *, double *, int, int, double *, MPI_Comm);
#endif
