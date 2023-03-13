#ifndef __SPS_UTILS_H__
#define __SPS_UTILS_H__
#include <dirent.h>
#include "config.h"
#include "input.h"
#include "my_mpi.h"


extern void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz);
int get_atom_num(char *symbol);
double get_mass(int atom_num);
char *get_symbol(int atom_num);
int name_filter(const struct dirent *info);
int check_unique(Config *config, Input *input, char *self);
double *get_rot_force(Input *input, double *force1, double *force2,
                      double *eigenmode, int n);
void get_cg_direction(double *direction, double *direction_old,
                      double *cg_direction, int n);
double *get_eigenmode(Input *input, int n, MPI_Comm comm);
void get_sphere_list(Config *config, Input *input, double *center, double cutoff,
                     int *atom_num, int **atom_list, MPI_Comm comm);
void expand_active_volume(Config *initial, Config *saddle, Input *input,
                          int *active_num, int *active_list, int *max_index,
                          MPI_Comm comm);
int postprocess(Config *initial, Config *final, Input *input,
                double *Ea, double *eigenmode, int count, int index,
                int global_num, int *global_list, double time, MPI_Comm comm);
#endif
