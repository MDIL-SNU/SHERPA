#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vasp_calculator.h"


void oneshot(Config *config, Input *input, double *energy, MPI_Comm comm)
{
    *energy = 0.0;
}


void oneshot_disp(Config *config, Input *input, double *energy, double *force,
                  int disp_num, int *disp_list, MPI_Comm comm)
{
    *energy = 0.0;
}


void atom_relax(Config *config, Input *input, double *energy, MPI_Comm comm)
{
    *energy = 0.0;
}
