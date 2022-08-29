#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"


/* large to small */
void int_sort(int *int_list, int list_len)
{
    int i;
    while (1) {
        int done = 1;
        for (i = 1; i < list_len; ++i) {
            if (int_list[i - 1] < int_list[i]) {
                int tmp_int = int_list[i - 1];
                int_list[i - 1] = int_list[i];
                int_list[i] = tmp_int;
                done = 0;
            }
        }
        if (done) {
            break;
        }
    }
}


inline void get_minimum_image(double *del, double *boxlo, double *boxhi,
                              double xy, double yz, double xz)
{
    double xprd, yprd, zprd;
    double xprd_half, yprd_half, zprd_half;
    xprd = boxhi[0] - boxlo[0];
    yprd = boxhi[1] - boxlo[1];
    zprd = boxhi[2] - boxlo[2];
    xprd_half = xprd * 0.5;
    yprd_half = yprd * 0.5;
    zprd_half = zprd * 0.5;
    while (fabs(del[2]) > zprd_half) {
        if (del[2] < 0.0) {
            del[2] += zprd;
            del[1] += yz;
            del[0] += xz;
        } else {
            del[2] -= zprd;
            del[1] -= yz;
            del[0] -= xz;
        }
    }
    while (fabs(del[1]) > yprd_half) {
        if (del[1] < 0.0) {
            del[1] += yprd;
            del[0] += xy;
        } else {
            del[1] -= yprd;
            del[0] -= xy;
        }
    }
    while (fabs(del[0]) > xprd_half) {
        if (del[0] < 0.0) {
            del[0] += xprd;
        } else {
            del[0] -= xprd;
        }
    }
    return;
}


int get_atom_num(char *symbol)
{
    int i;
    const char *name[92] = {"H", "He", "Li", "Be", "B",
    "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
    "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};
    for (i = 0; i < 92; ++i) {
        if (strncmp(symbol, name[i], strlen(symbol)) == 0) {
            return i + 1;
        }
    }
    return 0;
}


double get_mass(int atom_num)
{
  const double mass[92] = {1.0079, 4.0026, 6.941, 9.0122, 10.811,
  12.011, 14.007, 15.999, 18.998, 20.180, 22.990, 24.305, 26.982,
  28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956,
  47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546,
  65.380, 69.723, 72.640, 74.922, 78.960, 79.904, 83.798, 85.468,
  87.620, 88.906, 91.224, 92.906, 95.960, 98.000, 101.07, 102.91,
  106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90,
  131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 145.00,
  150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93,
  173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22,
  195.08, 196.97, 200.59, 204.38, 207.20, 208.98, 209.00, 210.00,
  222.00, 223.00, 226.00, 227.00, 232.04, 231.04, 238.03};
  return mass[atom_num - 1];
}


char *get_symbol(int atom_num)
{
    char *name[92] = {"H", "He", "Li", "Be", "B",
    "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
    "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};
    return name[atom_num - 1];
}


int name_filter(const struct dirent *info)
{
    if (strncmp(info->d_name, "Final", 5) == 0) {
        return 1;
    } else {
        return 0;
    }
}


/* nonpositive: not unique, positive: unique */
int check_unique(Config *config, Input *input, char *self, long long step)
{
    int i, j, errno, unique;
    struct dirent **namelist;

    char path[128];
    sprintf(path, "%s/%lld", input->output_dir, step);
    int count = scandir(path, &namelist, name_filter, NULL);
    if (count > 0) {
        for (i = 0; i < count; ++i) {
            if (strcmp(self, namelist[i]->d_name) == 0) {
                continue;
            }
            char filename[128];
            sprintf(filename, "%s/%lld/%s",
                    input->output_dir, step, namelist[i]->d_name);
            Config *tmp_config = (Config *)malloc(sizeof(Config));
            errno = read_config(tmp_config, input, filename);
            /* already deleted */
            if (errno > 0) {
                free(tmp_config);
                continue;
            }
            /* 0: identical, 1: different */
            unique = diff_config(tmp_config, config, 2 * input->max_step);
            free_config(tmp_config);
            if (unique == 0) {
                strtok(namelist[i]->d_name, "_");
                int index = atoi(strtok(NULL, "."));
                for (j = 0; j < count; ++j) {
                    free(namelist[j]);
                }
                free(namelist);
                return -index;
            }
        }
        for (j = 0; j < count; ++j) {
            free(namelist[j]);
        }
        free(namelist);
        return 1;
    } else {
        return 1;
    }
}


int count_unique(Input *input, long long step)
{
    int i;
    struct dirent **namelist;

    char path[128];
    sprintf(path, "%s/%lld", input->output_dir, step);
    int count = scandir(path, &namelist, name_filter, NULL);
    if (count > 0) {
        for (i = 0; i < count; ++i) {
            free(namelist[i]);
        }
        free(namelist);
    }
    return count;
}
