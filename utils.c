#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

inline double norm(double *vec) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}


inline double dot(double *vec1, double *vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}


inline void cross(double *vec1, double *vec2, double *vec3) {
    vec3[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vec3[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vec3[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}


inline double det(double (*mat)[3]) {
    int i;
    double cofactor[3];
    cofactor[0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    cofactor[1] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    cofactor[2] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];

    double det = 0.0;
    for (i = 0; i < 3; i++) {
        det += cofactor[i] * mat[0][i];
    }
    return det;
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


inline double normal_random(double mean, double std)
{
    double u, v, s;
    do {
        u = ((double)rand() / RAND_MAX) * 2 - 1;
        v = ((double)rand() / RAND_MAX) * 2 - 1;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    s = sqrt(-2 * log(s) / s);
    return mean + std * u * s;
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


int get_mask(Config *config, Input *input, char *rlx_cmd, char *fix_cmd,
             int *mask, int index, MPI_Comm comm)
{
    int i;
    int rank, size;
    double e_cutoff, f_cutoff, dist, ref[3], del[3];
    char tmp_cmd[64];
    int nummask = 0;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int q = config->tot_num / size;
    int r = config->tot_num % size;
    int begin = rank * q + ((rank > r) ? r : rank);
    int end = (r > rank) ? begin + q + 1 : begin + q;

    int *rlx = (int *)malloc(sizeof(int) * (end - begin));
    int *fix = (int *)malloc(sizeof(int) * (end - begin));
    int num_rlx = 0;
    int num_fix = 0;

    e_cutoff = input->cutoff;
    f_cutoff = e_cutoff + input->cutoff;

    ref[0] = config->pos[index * 3 + 0]; 
    ref[1] = config->pos[index * 3 + 1]; 
    ref[2] = config->pos[index * 3 + 2]; 
    for (i = begin; i < end; ++i) {
        del[0] = config->pos[i * 3 + 0] - ref[0]; 
        del[1] = config->pos[i * 3 + 1] - ref[1]; 
        del[2] = config->pos[i * 3 + 2] - ref[2]; 
        get_minimum_image(del, config->boxlo, config->boxhi,
                          config->xy, config->yz, config->xz);
        dist = norm(del);
        if (dist > f_cutoff) {
            continue;
        } else if (dist > e_cutoff) {
            fix[num_fix] = i + 1;
            num_fix++;
        } else {
            rlx[num_rlx] = i + 1;
            num_rlx++;
        }
    }

    /* rlx */
    int *counts = (int *)malloc(sizeof(int) * size);
    MPI_Allgather(&num_rlx, 1, MPI_INT, counts, 1, MPI_INT, comm);
    int *displ = (int *)malloc(sizeof(int) * size);
    displ[0] = 0;
    if (size > 1) {
        for (i = 1; i < size; ++i) {
            displ[i] = displ[i - 1] + counts[i - 1];
        } 
    }
    int tot_rlx;
    MPI_Allreduce(&num_rlx, &tot_rlx, 1, MPI_INT, MPI_SUM, comm);
    int *rlx_id = (int *)malloc(sizeof(int) * tot_rlx);
    MPI_Allgatherv(rlx, num_rlx, MPI_INT, rlx_id, counts, displ, MPI_INT, comm);
    
    sprintf(rlx_cmd, "group rlx id");
    for (i = 0; i < tot_rlx; ++i) {
        sprintf(tmp_cmd, " %d", rlx_id[i]);
        strcat(rlx_cmd, tmp_cmd);
        mask[nummask] = rlx_id[i];
        nummask++;
    }

    /* fix */
    MPI_Allgather(&num_fix, 1, MPI_INT, counts, 1, MPI_INT, comm);
    if (size > 1) {
        for (i = 1; i < size; ++i) {
            displ[i] = displ[i - 1] + counts[i - 1];
        }
    }
    int tot_fix;
    MPI_Allreduce(&num_fix, &tot_fix, 1, MPI_INT, MPI_SUM, comm);
    int *fix_id = (int *)malloc(sizeof(int) * tot_fix);
    MPI_Allgatherv(fix, num_fix, MPI_INT, fix_id, counts, displ, MPI_INT, comm);

    sprintf(fix_cmd, "group fix id");
    for (i = 0; i < tot_fix; ++i) {
        sprintf(tmp_cmd, " %d", fix_id[i]);
        strcat(fix_cmd, tmp_cmd);
    }

    free(rlx);
    free(fix);
    free(rlx_id);
    free(fix_id);
    free(counts);
    free(displ);

    return nummask;
}

