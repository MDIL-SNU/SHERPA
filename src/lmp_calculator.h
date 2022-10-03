#ifndef __LMP_CALCULATOR_H__
#define __LMP_CALCULATOR_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"

void *lmp_init(Config *config, Input *input, int lmpargc, char **lmpargv,
               MPI_Comm comm);
void oneshot(Config *config, Input *input, double *energy, double *force,
             MPI_Comm comm);
void oneshot_disp(Config *config, Input *input, double *energy, double *force,
                  int disp_num, int *disp_list, MPI_Comm comm);
void atom_relax(Config *config, Input *input, double *energy, MPI_Comm comm);
#endif
