#ifndef __UTILS_H__
#define __UTILS_H__
#include <dirent.h>
#include <mpi.h>
#include "config.h"
#include "input.h"


void int_sort_decrease(int *, int);
void int_sort_increase(int *, int);
extern void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz);
int get_atom_num(char *);
double get_mass(int);
char *get_symbol(int);
int name_filter(const struct dirent *);
int check_unique(Config *config, Input *input, char *);
double normal_random(double, double);
double norm(double *, int);
double *normalize(double *, int);
double dot(double *, double *, int);
double *parallel_vector(double *, double *, int);
double *perpendicular_vector(double *, double *, int);
void get_cg_direction(double *, double *, double *, int);
void cut_sphere(Config *, Input *, int, int *);
double *get_eigenmode(Input *, int, MPI_Comm);
void gen_list(Config *, Input *, double *, int *, int **, int *, int **, MPI_Comm);
int split_configs(Config *, Config *, Config *, Input *,
                  double *, int, int, int, int *, int, int *, MPI_Comm);
#endif
