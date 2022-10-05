#ifndef __VASP_CALCULATOR_H__
#define __VASP_CALCULATOR_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


void modify_incar(Input *input, char *filename, int ibrion);
void read_outcar(char *filename, double *pos, double *energy, double *force);
void append_outcar(char *filename1, char *filename2);
void oneshot(Config *config, Input *input, double *energy, double *force,
             MPI_Comm comm);
void oneshot_disp(Config *config, Input *input, double *energy, double *force,
                  int disp_num, int *disp_list, MPI_Comm comm);
void atom_relax(Config *config, Input *input, double *energy, MPI_Comm);
#endif
