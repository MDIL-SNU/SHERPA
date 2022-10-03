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
void trim_atoms(Config *config, int update_num, int *update_list);
double *get_eigenmode(Input *input, int n, MPI_Comm comm);
void set_active_volume(Config *config, Input *input, double *center,
                       int *update_num, int **update_list,
                       int *extract_num, int **extract_list, MPI_Comm comm);
int split_configs(Config *initial, Config *final, Config *config0, Input *input,
                  double *eigenmode, int count, int index,
                  int update_num, int *update_list,
                  int disp_num, int *disp_list,
                  MPI_Comm comm);
#endif
