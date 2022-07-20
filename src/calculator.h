#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include "config.h"
#include "input.h"

void *lmp_init(Config *, Input *, int, char **);
void oneshot(Config *, Input *, double *, double *, int, int *);
double atom_relax(Config *, Input *);
#endif
