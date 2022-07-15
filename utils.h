#ifndef __UTILS_H__
#define __UTILS_H__
#include <mpi.h>

extern double norm(double *vec);
extern double dot(double *vec1, double *vec2);
extern void cross(double *vec1, double *vec2, double *vec3);
extern double det(double (*mat)[3]);
extern void minimum_image(double *del, double *boxlo, double *boxhi,
                          double xy, double yz, double xz);
int get_atom_num(char *);
double get_mass(int);
char *get_symbol(int);
#endif
