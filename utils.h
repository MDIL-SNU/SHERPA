#ifndef __UTILS_H__
#define __UTILS_H__
#include <mpi.h>
#include "config.h"
#include "input.h"

extern double norm(double *vec);
extern double dot(double *vec1, double *vec2);
extern void cross(double *vec1, double *vec2, double *vec3);
extern double det(double (*mat)[3]);
extern void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz);
extern double normal_random(double, double);
int get_atom_num(char *);
double get_mass(int);
char *get_symbol(int);
int get_mask(Config *, Input *, char *, char *, int *, int, MPI_Comm);
#endif
