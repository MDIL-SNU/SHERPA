#include "config.h"
#include "linalg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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


/* convert into lmp basis */
static void convert_basis(Config *config)
{
    int i, j;
    double A_norm, AxB_norm, vol, tmp_value, tmp_pos[3];
    double A[3], B[3], C[3];
    double Ahat[3], AxB[3], AxBhat[3], AhatxB[3], AxBhatxAhat[3];
    double frac[3][3], trans[3][3];

    /* rotate bases into triangular matrix */
    for (i = 0; i < 3; i++) {
        A[i] = config->cell[0][i];
        B[i] = config->cell[1][i];
        C[i] = config->cell[2][i];
    }

    A_norm = norm(A, 1);
    Ahat[0] = A[0] / A_norm;
    Ahat[1] = A[1] / A_norm;
    Ahat[2] = A[2] / A_norm;

    cross(A, B, AxB);
    AxB_norm = norm(AxB, 1);
    AxBhat[0] = AxB[0] / AxB_norm;
    AxBhat[1] = AxB[1] / AxB_norm;
    AxBhat[2] = AxB[2] / AxB_norm;

    cross(Ahat, B, AhatxB);
    cross(AxBhat, Ahat, AxBhatxAhat);

    /* column vector (a b c) */
    config->cell[0][0] = A_norm;
    config->cell[0][1] = dot(B, Ahat, 1);
    config->cell[0][2] = dot(C, Ahat, 1);
    config->cell[1][0] = 0.0;
    config->cell[1][1] = norm(AhatxB, 1);
    config->cell[1][2] = dot(C, AxBhatxAhat, 1);
    config->cell[2][0] = 0.0;
    config->cell[2][1] = 0.0;
    config->cell[2][2] = dot(C, AxBhat, 1);

    /* edge & tilting */
    config->boxlo[0] = 0.0;
    config->boxlo[1] = 0.0;
    config->boxlo[2] = 0.0;
    config->boxhi[0] = config->cell[0][0];
    config->boxhi[1] = config->cell[1][1];
    config->boxhi[2] = config->cell[2][2];
    config->xy = config->cell[0][1];
    config->xz = config->cell[0][2];
    config->yz = config->cell[1][2];

    /* transformation matrix */
    vol = det(config->cell);
    cross(B, C, frac[0]);
    cross(C, A, frac[1]);
    cross(A, B, frac[2]);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            frac[i][j] /= vol;
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            trans[i][j] = config->cell[i][0] * frac[0][j] +
                          config->cell[i][1] * frac[1][j] +
                          config->cell[i][2] * frac[2][j];
        }
    }

    /* convert position */
    for (i = 0; i < config->tot_num; i++) {
        for (j = 0; j < 3; j++) {
            tmp_pos[j] = trans[j][0] * config->pos[i * 3 + 0] +
                         trans[j][1] * config->pos[i * 3 + 1] +
                         trans[j][2] * config->pos[i * 3 + 2];
        }
        config->pos[i * 3 + 0] = tmp_pos[0];
        config->pos[i * 3 + 1] = tmp_pos[1];
        config->pos[i * 3 + 2] = tmp_pos[2];
    }

    /* transpose lattice vector */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 3; j++) {
            tmp_value = config->cell[i][j];
            config->cell[i][j] = config->cell[j][i];
            config->cell[j][i] = tmp_value;
        }
    }
}


/* remove one atom from config */
void remove_atom(Config *config, int index)
{
    int i, j, ntype, tot_num;
    int acc_num = 0;
    for (i = 0; i < config->ntype; ++i) {
        acc_num += config->each_num[i];
        if (acc_num > index) {
            config->each_num[i]--;
            if (config->each_num[i] == 0) {
                for (j = i; j < config->ntype - 1; ++j) {
                    config->atom_num[i] = config->atom_num[i + 1];
                    config->each_num[i] = config->each_num[i + 1];
                }
                config->ntype--;
                ntype = config->ntype;
                config->atom_num = (int *)realloc(config->atom_num, sizeof(int) * ntype);
                config->each_num = (int *)realloc(config->each_num, sizeof(int) * ntype);
            }
            break;
        }
    }

    for (i = index; i < config->tot_num - 1; ++i) {
        config->fix[i] = config->fix[i + 1];
        config->pos[i * 3 + 0] = config->pos[(i + 1) * 3 + 0];
        config->pos[i * 3 + 1] = config->pos[(i + 1) * 3 + 1];
        config->pos[i * 3 + 2] = config->pos[(i + 1) * 3 + 2];
    }
    config->tot_num--;
    tot_num = config->tot_num;
    config->fix = (int *)realloc(config->fix, sizeof(int) * tot_num);
    config->pos = (double *)realloc(config->pos, sizeof(double) * tot_num * 3);
}


int read_config(Config *config, char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    int i, j, k;
    double tmp_pos[3];
    char line[1024], tmp_line[1024], *ptr;

    /* system name */
    ptr = fgets(line, 1024, fp);

    /* scale */
    ptr = fgets(line, 1024, fp);
    double scale = atof(line);

    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        ptr = fgets(line, 1024, fp);
        config->cell[i][0] = atof(strtok(line, " \n")) * scale;
        config->cell[i][1] = atof(strtok(NULL, " \n")) * scale;
        config->cell[i][2] = atof(strtok(NULL, " \n")) * scale;
    }

    /* the number of type */
    ptr = fgets(line, 1024, fp);
    strncpy(tmp_line, line, 1024);
    config->ntype = 0;
    ptr = strtok(line, " \r\n");
    while (ptr != NULL) {
        if (strlen(ptr) > 0) {
            config->ntype++;
        }
        ptr = strtok(NULL, " \r\n");
    }

    /* atomic number type */
    ptr = strtok(tmp_line, " \r\n");
    config->atom_num = (int *)malloc(sizeof(int) * config->ntype);
    for (i = 0; i < config->ntype; ++i) {
        config->atom_num[i] = get_atom_num(ptr);
        ptr = strtok(NULL, " \r\n");
    }

    /* each number of type */
    config->each_num = (int *)malloc(sizeof(int) * config->ntype);
    config->tot_num = 0;
    ptr = fgets(line, 1024, fp);
    config->each_num[0] = atoi(strtok(line, " \n"));
    config->tot_num += config->each_num[0];
    if (config->ntype > 1) {
        for (i = 1; i < config->ntype; ++i) {
            config->each_num[i] = atoi(strtok(NULL, " \n"));
            config->tot_num += config->each_num[i];
        }
    }

    /* positions and constraint */
    ptr = fgets(line, 1024, fp);
    config->fix = (int *)calloc(config->tot_num, sizeof(int));
    config->pos = (double *)malloc(sizeof(double) * config->tot_num * 3);
    int constraint = 0;
    if (strncasecmp(line, "S", 1) == 0) {
        constraint = 1;
        ptr = fgets(line, 1024, fp);
    }
    if (strncasecmp(line, "D", 1) == 0) {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, 1024, fp);
            tmp_pos[0] = atof(strtok(line, " \n"));
            tmp_pos[1] = atof(strtok(NULL, " \n"));
            tmp_pos[2] = atof(strtok(NULL, " \n"));
            for (j = 0; j < 3; ++j) {
                config->pos[i * 3 + j] = tmp_pos[0] * config->cell[0][j]
                                       + tmp_pos[1] * config->cell[1][j]
                                       + tmp_pos[2] * config->cell[2][j];
            }
            if (constraint > 0) {
                ptr = strtok(NULL, " \n");
                if (strcmp(ptr, "F") == 0) {
                    config->fix[i] = 1;
                }
            }
        }
    } else {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, 1024, fp);
            config->pos[i * 3 + 0] = atof(strtok(line, " \n"));
            config->pos[i * 3 + 1] = atof(strtok(NULL, " \n"));
            config->pos[i * 3 + 2] = atof(strtok(NULL, " \n"));
            if (constraint > 0) {
                ptr = strtok(NULL, " \n");
                if (strcmp(ptr, "F") == 0) {
                    config->fix[i] = 1;
                }
            }
        }
    }
    convert_basis(config);
    fclose(fp);
    return 0;
}


void write_config(Config *config, char *filename, char *header, char *mode)
{
    int i;
    FILE *fp;
    fp = fopen(filename, mode);

    /* title */
    fprintf(fp, "%s\n", header);

    /* scale */
    fputs("1.0\n", fp);
    
    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        fprintf(fp, " %.15f %.15f %.15f\n",
                config->cell[i][0], config->cell[i][1], config->cell[i][2]);
    }

    /* symbols */
    for (i = 0; i < config->ntype; ++i) {
        fprintf(fp, " %s", get_symbol(config->atom_num[i]));
    }
    fputs("\n", fp);

    /* the number of each type */
    for (i = 0; i < config->ntype; ++i) {
        fprintf(fp, " %d", config->each_num[i]);
    }
    fputs("\n", fp);

    /* positions and constraint */
    fputs("Selective dynamics\n", fp);
    fputs("Cartesian\n", fp);
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] > 0) {
            fprintf(fp, "  %.15f  %.15f  %.15f  F F F\n",
                    config->pos[i * 3 + 0],
                    config->pos[i * 3 + 1],
                    config->pos[i * 3 + 2]);
        } else {
            fprintf(fp, "  %.15f  %.15f  %.15f  T T T\n",
                    config->pos[i * 3 + 0],
                    config->pos[i * 3 + 1],
                    config->pos[i * 3 + 2]);
        }
    }
    fclose(fp);
}


void copy_config(Config *config2, Config *config1)
{
    int i;

    config2->ntype = config1->ntype;
    config2->tot_num = config1->tot_num;
    config2->cell[0][0] = config1->cell[0][0];
    config2->cell[0][1] = config1->cell[0][1];
    config2->cell[0][2] = config1->cell[0][2];
    config2->cell[1][0] = config1->cell[1][0];
    config2->cell[1][1] = config1->cell[1][1];
    config2->cell[1][2] = config1->cell[1][2];
    config2->cell[2][0] = config1->cell[2][0];
    config2->cell[2][1] = config1->cell[2][1];
    config2->cell[2][2] = config1->cell[2][2];

    config2->boxlo[0] = config1->boxlo[0];
    config2->boxlo[1] = config1->boxlo[1];
    config2->boxlo[2] = config1->boxlo[2];
    config2->boxhi[0] = config1->boxhi[0];
    config2->boxhi[1] = config1->boxhi[1];
    config2->boxhi[2] = config1->boxhi[2];
    config2->xy = config1->xy;
    config2->xz = config1->xz;
    config2->yz = config1->yz;

    config2->atom_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->each_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->fix = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->pos = (double *)malloc(sizeof(double) * config1->tot_num * 3);

    for (i = 0; i < config1->ntype; ++i) {
        config2->atom_num[i] = config1->atom_num[i];
        config2->each_num[i] = config1->each_num[i];
    }
    for (i = 0; i < config1->tot_num; ++i) {
        config2->fix[i] = config1->fix[i];
        config2->pos[i * 3 + 0] = config1->pos[i * 3 + 0];
        config2->pos[i * 3 + 1] = config1->pos[i * 3 + 1];
        config2->pos[i * 3 + 2] = config1->pos[i * 3 + 2];
    }
}


void free_config(Config *config)
{
    free(config->atom_num);
    free(config->each_num);
    free(config->fix);
    free(config->pos);
    free(config);
}
