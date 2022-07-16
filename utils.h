#ifndef __UTILS_H__
#define __UTILS_H__

extern void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz);
int get_atom_num(char *);
double get_mass(int);
char *get_symbol(int);
#endif
