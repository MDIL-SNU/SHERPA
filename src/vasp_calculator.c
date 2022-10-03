#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vasp_calculator.h"
#include "utils.h"


double oneshot(Config *config, Input *input, MPI_Comm comm)
{
    double energy = 0.0;
    return energy;
}


void oneshot_disp(Config *config, Input *input, double *energy, double *force,
                  int disp_num, int *disp_list, MPI_Comm comm)
{
    return;
}


double atom_relax(Config *config, Input *input, MPI_Comm comm)
{
    double pe = 0.0;
    return pe;
}
