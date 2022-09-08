#ifndef __TARGET_H__
#define __TARGET_H__
#include "config.h"
#include "input.h"

int read_target(Config *, Input *, int **, int *, int *);
void write_target(Input *, int *, int);
int get_index(Config *, Config *);
#endif
