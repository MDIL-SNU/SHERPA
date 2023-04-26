#ifndef __UTILS_H__
#define __UTILS_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


extern void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz);
int check_unique(Config *config, Input *input);
void concat_files(char *filename1, char *filename2);
double *get_rot_force(Input *input, double *force1, double *force2,
                      double *eigenmode, int n);
void get_cg_direction(double *direction, double *direction_old,
                      double *cg_direction, int n);
double *get_eigenmode(Input *input, int n, MPI_Comm comm);
void get_sphere_list(Config *config, Input *input, double *center, double cutoff,
                     int *atom_num, int **atom_list, MPI_Comm comm);
void expand_active_volume(Config *initial, Config *saddle, Input *input,
                          int *active_num, int *active_list, MPI_Comm comm);
int diff_config(Config *config1, Config *config2, double tol);
int split_config(Config *initial, Config *saddle, Config *final, Input *input,
                double *Ea, double *dE, double eigenvalue, double *eigenmode,
                int active_num, int *active_list, int count, int index, MPI_Comm comm);
#endif
