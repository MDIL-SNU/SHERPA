#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


typedef struct _Calc
{
    void *vasp;
    int initialized;
} Calc;

void oneshot(Calc *calc, Config *config, Input *input,
             double *energy, double *force, MPI_Comm comm);
void atom_relax(Calc *calc, Config *config, Input *input,
                double *energy, MPI_Comm comm);
void free_calc(Calc *calc);
#endif
