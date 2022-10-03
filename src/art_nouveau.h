#ifndef __ART_NOUVEAU__
#define __ART_NOUVEAU__
#include "config.h"
#include "input.h"
#include "my_mpi.h"


int art_nouveau(Config *initial, Config *final, Input *input,
                double *full_eigenmode, int count, int index, double *Ea,
                MPI_Comm comm);
#endif
