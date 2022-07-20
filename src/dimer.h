#ifndef __DIMER_H__
#define __DIMER_H__
#include "config.h"
#include "input.h"


double normal_random(double, double);
double *normalize(double *, int);
double *parallel_vector(double *, double *, int);
double *perpendicular_vector(double *, double *, int);
double *gen_eigenmode(Input *, int);
double *get_rot_force(Input *, double *, double *, double *, int);
double *projected_force(double *, double *, double, int);
double *displace(Input *, int, int, int *, int);
void rotate_vector(double *, double *, double **, double **, int, double);
void cut_sphere(Config *, Input *, double *);
void rotate(Config *, Input *, int, int *, double *, int, int);
void translate(Config *, Input *, int, int *, double *, double *, double *);
double dimer(Config *, Input *, int, int);
#endif
