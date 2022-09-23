#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "utils.h"


void int_sort_decrease(int *int_list, int list_len)
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


void int_sort_increase(int *int_list, int list_len)
{
    int i;
    while (1) {
        int done = 1;
        for (i = 1; i < list_len; ++i) {
            if (int_list[i - 1] > int_list[i]) {
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
int check_unique(Config *config, Input *input, char *self)
{
    int i, j, errno, unique;
    struct dirent **namelist;

    int ncount = scandir(input->output_dir, &namelist, name_filter, NULL);
    if (ncount > 0) {
        for (i = 0; i < ncount; ++i) {
            if (strcmp(self, namelist[i]->d_name) == 0) {
                continue;
            }
            char filename[512];
            sprintf(filename, "%s/%s",
                    input->output_dir, namelist[i]->d_name);
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
                int count = atoi(strtok(NULL, "_"));
                for (j = 0; j < ncount; ++j) {
                    free(namelist[j]);
                }
                free(namelist);
                return -count;
            }
        }
        for (j = 0; j < ncount; ++j) {
            free(namelist[j]);
        }
        free(namelist);
        return 1;
    } else {
        return 1;
    }
}


double normal_random(double mean, double std)
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


double norm(double *vec, int n)
{
    int i;
    double output = 0.0;
    for (i = 0; i < n; ++i) {
        output += vec[i * 3 + 0] * vec[i * 3 + 0];
        output += vec[i * 3 + 1] * vec[i * 3 + 1];
        output += vec[i * 3 + 2] * vec[i * 3 + 2];
    }
    return sqrt(output);
}


double *normalize(double *vec, int n)
{
    int i;
    double *output = (double *)malloc(sizeof(double) * n * 3); 
    double magnitude = norm(vec, n);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = vec[i * 3 + 0] / magnitude;
        output[i * 3 + 1] = vec[i * 3 + 1] / magnitude;
        output[i * 3 + 2] = vec[i * 3 + 2] / magnitude;
    }
    return output;
}


double dot(double *vec1, double *vec2, int n)
{
    int i;
    double output = 0.0;
    for (i = 0; i < n; ++i) {
        output += vec1[i * 3 + 0] * vec2[i * 3 + 0];
        output += vec1[i * 3 + 1] * vec2[i * 3 + 1];
        output += vec1[i * 3 + 2] * vec2[i * 3 + 2];
    }
    return output;
}


double *parallel_vector(double *vector, double *unit, int n)
{
    int i; 
    double magnitude = dot(vector, unit, n);
    double *output = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = magnitude * unit[i * 3 + 0];
        output[i * 3 + 1] = magnitude * unit[i * 3 + 1];
        output[i * 3 + 2] = magnitude * unit[i * 3 + 2];
    } 
    return output;
}


double *perpendicular_vector(double *vector, double *unit, int n)
{
    int i;
    double *tmp_vector = parallel_vector(vector, unit, n);
    double *output = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = vector[i * 3 + 0] - tmp_vector[i * 3 + 0];
        output[i * 3 + 1] = vector[i * 3 + 1] - tmp_vector[i * 3 + 1];
        output[i * 3 + 2] = vector[i * 3 + 2] - tmp_vector[i * 3 + 2];
    } 
    free(tmp_vector);
    return output;
}


void get_cg_direction(double *direction, double *direction_old,
                      double *cg_direction, int n)
{
    int i;
    double old_norm = norm(direction_old, n);
    double betaPR;
    if (fabs(old_norm) > 1e-8) {
        double *ddirection = (double *)malloc(sizeof(double) * n * 3);
        for (i = 0; i < n; ++i) {
            ddirection[i * 3 + 0] = direction[i * 3 + 0]
                                  - direction_old[i * 3 + 0];
            ddirection[i * 3 + 1] = direction[i * 3 + 1] 
                                  - direction_old[i * 3 + 1];
            ddirection[i * 3 + 2] = direction[i * 3 + 2] 
                                  - direction_old[i * 3 + 2];
        }
        betaPR = dot(direction, ddirection, n) / old_norm; 
        free(ddirection);
    } else {
        betaPR = 0.0;
    }
    if (betaPR < 0.0) {
        betaPR = 0.0;
    }
    for (i = 0; i < n; ++i) {
        cg_direction[i * 3 + 0] = direction[i * 3 + 0]
                                + cg_direction[i * 3 + 0] * betaPR;
        cg_direction[i * 3 + 1] = direction[i * 3 + 1] 
                                + cg_direction[i * 3 + 1] * betaPR;
        cg_direction[i * 3 + 2] = direction[i * 3 + 2] 
                                + cg_direction[i * 3 + 2] * betaPR;
        direction_old[i * 3 + 0] = direction[i * 3 + 0];
        direction_old[i * 3 + 1] = direction[i * 3 + 1];
        direction_old[i * 3 + 2] = direction[i * 3 + 2];
    }
}


/* not normalized */
double *get_eigenmode(Input *input, int n, MPI_Comm comm)
{
    int i, rank, size;
    double *eigenmode = (double *)malloc(sizeof(double) * n * 3);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    int q = n / input->ncore;
    int r = n % input->ncore;
    int begin = local_rank * q + ((local_rank > r) ? r : local_rank);
    int end = begin + q;
    if (r > local_rank) {
        end++;
    }
    for (i = begin; i < end; ++i) {
        eigenmode[i * 3 + 0] = normal_random(0, input->stddev);
        eigenmode[i * 3 + 1] = normal_random(0, input->stddev);
        eigenmode[i * 3 + 2] = normal_random(0, input->stddev);
    }
    int count = (end - begin) * 3;
    int *counts = (int *)malloc(sizeof(int) * input->ncore);
    MPI_Allgather(&count, 1, MPI_INT, counts, 1, MPI_INT, comm);
    int *disp = (int *)malloc(sizeof(int) * input->ncore);
    disp[0] = 0;
    if (input->ncore > 1) {
        for (i = 1; i < input->ncore; ++i) {
            disp[i] = disp[i - 1] + counts[i - 1];
        }
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   eigenmode, counts, disp, MPI_DOUBLE, comm); 
    free(disp);
    free(counts);
    return eigenmode;
}


/* update_list < 2 * cutoff
   extract_list < acti_cutoff */
void gen_list(Config *config, Input *input, double *center,
              int *update_num, int **update_list,
              int *extract_num, int **extract_list, MPI_Comm comm)
{
    int i, rank, size;
    double del[3];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    int q = config->tot_num / input->ncore;
    int r = config->tot_num % input->ncore;
    int begin = local_rank * q + ((local_rank > r) ? r : local_rank);
    int end = begin + q;
    if (r > local_rank) {
        end++;
    }
    int tmp_update_num = 0;
    int *tmp_update_list = (int *)malloc(sizeof(int) * config->tot_num);
    int tmp_extract_num = 0;
    int *tmp_extract_list = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = begin; i < end; ++i) {
        del[0] = config->pos[i * 3 + 0] - center[0];
        del[1] = config->pos[i * 3 + 1] - center[1];
        del[2] = config->pos[i * 3 + 2] - center[2];
        get_minimum_image(del, config->boxlo, config->boxhi,
                          config->xy, config->yz, config->xz);
        double dist = sqrt(del[0] * del[0]
                         + del[1] * del[1]
                         + del[2] * del[2]);
        if (dist < 2 * input->pair_cutoff) {
            tmp_update_list[tmp_update_num] = i;
            tmp_update_num++; 
            if (dist < input->acti_cutoff) {
                tmp_extract_list[tmp_extract_num] = i;
                tmp_extract_num++;
            }
        }
    }
    MPI_Allreduce(&tmp_update_num, update_num, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&tmp_extract_num, extract_num, 1, MPI_INT, MPI_SUM, comm);
    *update_list = (int *)malloc(sizeof(int) * (*update_num));
    *extract_list = (int *)malloc(sizeof(int) * (*extract_num));

    int *counts = (int *)malloc(sizeof(int) * input->ncore);
    int *disp = (int *)malloc(sizeof(int) * input->ncore);

    /* update list */
    MPI_Allgather(&tmp_update_num, 1, MPI_INT, counts, 1, MPI_INT, comm);
    disp[0] = 0;
    if (input->ncore > 1) {
        for (i = 1; i < input->ncore; ++i) {
            disp[i] = disp[i - 1] + counts[i - 1];
        }
    }
    MPI_Allgatherv(tmp_update_list, tmp_update_num, MPI_INT,
                   *update_list, counts, disp, MPI_INT, comm);

    /* extract list */
    MPI_Allgather(&tmp_extract_num, 1, MPI_INT, counts, 1, MPI_INT, comm);
    disp[0] = 0;
    if (input->ncore > 1) {
        for (i = 1; i < input->ncore; ++i) {
            disp[i] = disp[i - 1] + counts[i - 1];
        }
    }
    MPI_Allgatherv(tmp_extract_list, tmp_extract_num, MPI_INT,
                   *extract_list, counts, disp, MPI_INT, comm);

    free(tmp_update_list);
    free(tmp_extract_list);
    free(counts);
    free(disp);

    /* sort */
    int_sort_increase(*update_list, *update_num);
    int_sort_increase(*extract_list, *extract_num);
}


int split_configs(Config *initial, Config *final, Config *config0, Input *input,
                  double *eigenmode, int count, int index,  
                  int update_num, int *update_list,
                  int disp_num, int *disp_list,
                  MPI_Comm comm)
{
    int i, rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    Config *config1 = (Config *)malloc(sizeof(Config));
    Config *config2 = (Config *)malloc(sizeof(Config));
    int trial = 1;
    while (1) {
        copy_config(config1, config0);
        copy_config(config2, config0);
        for (i = 0; i < disp_num; ++i) {
            config1->pos[disp_list[i] * 3 + 0] = config0->pos[disp_list[i] * 3 + 0]
                                               + 0.1 * trial * eigenmode[i * 3 + 0];
            config1->pos[disp_list[i] * 3 + 1] = config0->pos[disp_list[i] * 3 + 1]
                                               + 0.1 * trial * eigenmode[i * 3 + 1];
            config1->pos[disp_list[i] * 3 + 2] = config0->pos[disp_list[i] * 3 + 2]
                                               + 0.1 * trial * eigenmode[i * 3 + 2];
            config2->pos[disp_list[i] * 3 + 0] = config0->pos[disp_list[i] * 3 + 0]
                                               - 0.1 * trial * eigenmode[i * 3 + 0];
            config2->pos[disp_list[i] * 3 + 1] = config0->pos[disp_list[i] * 3 + 1]
                                               - 0.1 * trial * eigenmode[i * 3 + 1];
            config2->pos[disp_list[i] * 3 + 2] = config0->pos[disp_list[i] * 3 + 2]
                                               - 0.1 * trial * eigenmode[i * 3 + 2];
        }
        atom_relax(config1, input, comm); 
        atom_relax(config2, input, comm); 
        if (diff_config(config1, config2, 2 * input->max_step) == 1) {
            break;
        } else {
            free_config(config1);
            free_config(config2);
            config1 = (Config *)malloc(sizeof(Config));
            config2 = (Config *)malloc(sizeof(Config));
        }
        trial++;
    }

    int diff1 = diff_config(initial, config1, 2 * input->max_step);
    int diff2 = diff_config(initial, config2, 2 * input->max_step);
    /* log */
    if (diff1 * diff2 > 0) {
        if (local_rank == 0) {
            char filename[128];
            sprintf(filename, "%s/Dimer_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            fputs("----------------------------------------------------------------------------\n", fp);
            fputs(" Saddle state: disconnected\n", fp);
            fclose(fp);
        }
        free_config(config1);
        free_config(config2);
        return 1;
    } else {
        if (local_rank == 0) {
            char  filename[128];
            sprintf(filename, "%s/Dimer_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            fputs("----------------------------------------------------------------------------\n", fp);
            fputs(" Saddle state: connected\n", fp);
            fclose(fp);
        }
        if (diff1 == 0) {
            for (i = 0; i < config0->tot_num; ++i) {
                final->pos[update_list[i] * 3 + 0] = config2->pos[i * 3 + 0];
                final->pos[update_list[i] * 3 + 1] = config2->pos[i * 3 + 1];
                final->pos[update_list[i] * 3 + 2] = config2->pos[i * 3 + 2];
            }
        } else {
            for (i = 0; i < config0->tot_num; ++i) {
                final->pos[update_list[i] * 3 + 0] = config1->pos[i * 3 + 0];
                final->pos[update_list[i] * 3 + 1] = config1->pos[i * 3 + 1];
                final->pos[update_list[i] * 3 + 2] = config1->pos[i * 3 + 2];
            }
        }
        if (local_rank == 0) {
            char filename[128];
            sprintf(filename, "%s/Final_%d_%d.POSCAR",
                    input->output_dir, count, index);
            write_config(final, filename, "w");
        }
        free_config(config1);
        free_config(config2);
        return 0;
    }
}



