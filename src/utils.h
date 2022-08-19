#ifndef __UTILS_H__
#define __UTILS_H__
#include <dirent.h>
#include "config.h"
#include "input.h"


void int_sort(int *, int);
extern void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz);
int get_atom_num(char *);
double get_mass(int);
char *get_symbol(int);
int name_filter(const struct dirent *);
int check_unique(Config *config, Input *input);
#endif
