#ifndef __DIMER_H__
#define __DIMER_H__
#include "config.h"
#include "input.h"


double normal_random(double, double);
double **parallel_vector(double **, double **, int);
double **perpendicular_vector(double **, double **, int);
double **gen_eigenmode(Input *, int);
double **get_rot_force(Input *, double **, double **, int, int *, double **);
void cut_sphere(Config *, Input *, double *);
void rotate(Config *, Config *, Input *, int, int *, double **);
void translate();
void dimer(Config *, Input *, int);
#endif
