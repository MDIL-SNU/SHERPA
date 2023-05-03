#include "utils.h"
#include "calculator.h"
#include "linalg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


void get_minimum_image(double *del, double *boxlo, double *boxhi,
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


/* -1: unique, nonnegative: degenerate */
int check_unique(Config *config, Input *input)
{
    int i, j, unique, count, errno;
    char filename[128];
    FILE *rp = fopen("./Saddle.POSCAR", "r");
    /* first open */
    if (rp == NULL) {
        return -1;
    } else {
        char line[1024], header[128], *ptr;
        while (1) {
            if (fgets(header, 128, rp) == NULL) {
                break;
            }
            FILE *wp = fopen("./tmp_Saddle.POSCAR", "w");
            fputs(header, wp);
            for (i = 0; i < 8; ++i) {
                fgets(line, 1024, rp);
                fputs(line, wp);
            }
            for (i = 0; i < config->tot_num; ++i) {
                fgets(line, 1024, rp);
                fputs(line, wp);
            }
            fclose(wp);
            Config *tmp_config = (Config *)malloc(sizeof(Config));
            read_config(tmp_config, "./tmp_Saddle.POSCAR");
            int diff = diff_config(config, tmp_config, input->diff_tol);
            /* 0: identical */
            if (diff == 0) {
                fclose(rp);
                remove("./tmp_Saddle.POSCAR");
                free_config(tmp_config);
                return atoi(strtok(header, "_"));
            }
            free_config(tmp_config);
        }
        fclose(rp);
        remove("./tmp_Saddle.POSCAR");
        return -1;
    }
}


/* filename1 <- filename2 */
void concat_files(char *filename1, char *filename2)
{
    FILE *wp = fopen(filename1, "a");
    FILE *rp = fopen(filename2, "r");
    char line[1024];
    while (fgets(line, 1024, rp) != NULL) {
        fputs(line, wp);
    }
    fclose(wp);
    fclose(rp);
}


double *get_rot_force(Input *input, double *force1, double *force2,
                      double *eigenmode, int n)
{
    int i;
    double *dforce = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        dforce[i * 3 + 0] = force1[i * 3 + 0] - force2[i * 3 + 0];
        dforce[i * 3 + 1] = force1[i * 3 + 1] - force2[i * 3 + 1];
        dforce[i * 3 + 2] = force1[i * 3 + 2] - force2[i * 3 + 2];
    }
    double *rot_force = perpendicular_vector(dforce, eigenmode, n);
    for (i = 0; i < n; ++i) {
        rot_force[i * 3 + 0] /= 2 * input->finite_diff;
        rot_force[i * 3 + 1] /= 2 * input->finite_diff;
        rot_force[i * 3 + 2] /= 2 * input->finite_diff;
    }
    free(dforce);
    return rot_force;
}


void get_cg_direction(double *direction, double *direction_old,
                      double *cg_direction, int n)
{
    int i;
    double a1 = dot(direction, direction_old, n);
    double a2 = norm(direction_old, n);
    double gamma;
    if ((a1 < 0.5 * a2) && (fabs(a2) > 1e-8)) {
        double *ddirection = (double *)malloc(sizeof(double) * n * 3);
        for (i = 0; i < n; ++i) {
            ddirection[i * 3 + 0] = direction[i * 3 + 0]
                                  - direction_old[i * 3 + 0];
            ddirection[i * 3 + 1] = direction[i * 3 + 1] 
                                  - direction_old[i * 3 + 1];
            ddirection[i * 3 + 2] = direction[i * 3 + 2] 
                                  - direction_old[i * 3 + 2];
        }
        gamma = dot(direction, ddirection, n) / a2; 
        free(ddirection);
    } else {
        gamma = 0.0;
    }
    for (i = 0; i < n; ++i) {
        cg_direction[i * 3 + 0] = direction[i * 3 + 0]
                                + cg_direction[i * 3 + 0] * gamma;
        cg_direction[i * 3 + 1] = direction[i * 3 + 1] 
                                + cg_direction[i * 3 + 1] * gamma;
        cg_direction[i * 3 + 2] = direction[i * 3 + 2] 
                                + cg_direction[i * 3 + 2] * gamma;
        direction_old[i * 3 + 0] = direction[i * 3 + 0];
        direction_old[i * 3 + 1] = direction[i * 3 + 1];
        direction_old[i * 3 + 2] = direction[i * 3 + 2];
    }
}


/* not normalized */
double *get_random_vector(Input *input, int n, MPI_Comm comm)
{
    int i, rank;
    double *vector = (double *)malloc(sizeof(double) * n * 3);

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
        vector[i * 3 + 0] = normal_random(0, 1);
        vector[i * 3 + 1] = normal_random(0, 1);
        vector[i * 3 + 2] = normal_random(0, 1);
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
                   vector, counts, disp, MPI_DOUBLE, comm);
    free(disp);
    free(counts);
    return vector;
}


void get_sphere_list(Config *config, Input *input, double *center,
                     double cutoff, int *atom_num, int **atom_list,
                     MPI_Comm comm)
{
    int i, rank;
    double del[3];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    int q = config->tot_num / input->ncore;
    int r = config->tot_num % input->ncore;
    int begin = local_rank * q + ((local_rank > r) ? r : local_rank);
    int end = begin + q;
    if (r > local_rank) {
        end++;
    }
    int tmp_num = 0;
    int *tmp_list = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = begin; i < end; ++i) {
        del[0] = config->pos[i * 3 + 0] - center[0];
        del[1] = config->pos[i * 3 + 1] - center[1];
        del[2] = config->pos[i * 3 + 2] - center[2];
        get_minimum_image(del, config->boxlo, config->boxhi,
                          config->xy, config->yz, config->xz);
        if (norm(del, 1) < cutoff) {
            tmp_list[tmp_num] = i;
            tmp_num++;
        }
    }
    MPI_Allreduce(&tmp_num, atom_num, 1, MPI_INT, MPI_SUM, comm);
    *atom_list = (int *)malloc(sizeof(int) * (*atom_num));

    int *counts = (int *)malloc(sizeof(int) * input->ncore);
    int *disp = (int *)malloc(sizeof(int) * input->ncore);

    /* all gather */
    MPI_Allgather(&tmp_num, 1, MPI_INT, counts, 1, MPI_INT, comm);
    disp[0] = 0;
    if (input->ncore > 1) {
        for (i = 1; i < input->ncore; ++i) {
            disp[i] = disp[i - 1] + counts[i - 1];
        }
    }
    MPI_Allgatherv(tmp_list, tmp_num, MPI_INT,
                   *atom_list, counts, disp, MPI_INT, comm);
    free(tmp_list);
    free(counts);
    free(disp);
}


void expand_active_volume(Config *initial, Config *saddle, Input *input,
                          double cutoff, int *active_num, int *active_list,
                          MPI_Comm comm)
{
    int i, j, max_index;
    double del[3];
    double dmax = 0.0;
    /* maximally displaced atom */
    for (i = 0; i < (*active_num); ++i) {
        del[0] = saddle->pos[active_list[i] * 3 + 0]
               - initial->pos[active_list[i] * 3 + 0];
        del[1] = saddle->pos[active_list[i] * 3 + 1]
               - initial->pos[active_list[i] * 3 + 1];
        del[2] = saddle->pos[active_list[i] * 3 + 2]
               - initial->pos[active_list[i] * 3 + 2];
        if (norm(del, 1) > dmax) {
            dmax = norm(del, 1);
            max_index = active_list[i];
        }
    }

    /* generate lists */
    double center[3];
    center[0] = saddle->pos[max_index * 3 + 0];
    center[1] = saddle->pos[max_index * 3 + 1];
    center[2] = saddle->pos[max_index * 3 + 2];
    int tmp_num;
    int *tmp_list;
    get_sphere_list(saddle, input, center, cutoff,
                    &tmp_num, &tmp_list, comm);
    /* append active_list */
    int tmp_active_num = *active_num;
    for (i = 0; i < tmp_num; ++i) {
        int new = 1;
        for (j = 0; j < tmp_active_num; ++j) {
            if (active_list[j] == tmp_list[i]) {
                new = 0;
                break;
            }
        }
        if ((new == 1) && (saddle->fix[tmp_list[i]] == 0)) {
            active_list[*active_num] = tmp_list[i];
            (*active_num)++;
        }
    }
    free(tmp_list);
}


/* 0: identical, 1: different */
int diff_config(Config *config1, Config *config2, double tol)
{
    int i;
    double del[3];
    if (config1->ntype != config2->ntype) {
        return 1;
    }
    for (i = 0; i < config1->ntype; ++i) {
        if (config1->each_num[i] != config2->each_num[i]) {
            return 1;
        }
    }
    for (i = 0; i < config1->tot_num; ++i) {
        del[0] = config2->pos[i * 3 + 0] - config1->pos[i * 3 + 0];
        del[1] = config2->pos[i * 3 + 1] - config1->pos[i * 3 + 1];
        del[2] = config2->pos[i * 3 + 2] - config1->pos[i * 3 + 2];
        get_minimum_image(del, config1->boxlo, config1->boxhi,
                          config1->xy, config1->yz, config1->xz);
        if (norm(del, 1) > tol) {
            return 1;
        }
    }
    return 0;
}


int split_config(Config *initial, Config *saddle, Config *final, Input *input,
                double *Ea, double *dE, double eigenvalue, double *eigenmode,
                int active_num, int *active_list, int count, int index,
                MPI_Comm comm)
{
    int i, j, rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * saddle->tot_num * 3);
    /* initial oneshot */
    oneshot(initial, input, &energy0, force0, comm);
    double initial_energy = energy0;
    /* saddle oneshot */
    oneshot(saddle, input, &energy0, force0, comm);
    double saddle_energy = energy0;
    *Ea = saddle_energy - initial_energy;
    free(force0);

    /* step size */
    double dr = -2 * input->f_tol / eigenvalue;

    /* forward image */
    Config *config1 = (Config *)malloc(sizeof(Config));
    copy_config(config1, saddle);
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * config1->tot_num * 3);
    for (i = 0; i < 10; ++i) {
        for (j = 0; j < active_num; ++j) {
            config1->pos[active_list[j] * 3 + 0] += dr * eigenmode[j * 3 + 0];
            config1->pos[active_list[j] * 3 + 1] += dr * eigenmode[j * 3 + 1];
            config1->pos[active_list[j] * 3 + 2] += dr * eigenmode[j * 3 + 2];
        }
        oneshot(config1, input, &energy1, force1, comm);
        if (energy1 < energy0) {
            break;
        }
    }
    free(force1);
    atom_relax(config1, input, &energy1, comm);
    int diff1 = diff_config(initial, config1, input->diff_tol);

    /* backward image */
    Config *config2 = (Config *)malloc(sizeof(Config));
    copy_config(config2, saddle);
    double energy2;
    double *force2 = (double *)malloc(sizeof(double) * config2->tot_num * 3);
    for (i = 0; i < 10; ++i) {
        for (j = 0; j < active_num; ++j) {
            config2->pos[active_list[j] * 3 + 0] -= dr * eigenmode[j * 3 + 0];
            config2->pos[active_list[j] * 3 + 1] -= dr * eigenmode[j * 3 + 1];
            config2->pos[active_list[j] * 3 + 2] -= dr * eigenmode[j * 3 + 2];
        }
        oneshot(config2, input, &energy2, force2, comm);
        if (energy2 < energy0) {
            break;
        }
    }
    free(force2);
    atom_relax(config2, input, &energy2, comm);
    int diff2 = diff_config(initial, config2, input->diff_tol);

    /* disconnected */
    if ((diff1 == 1) && (diff2 == 1)) {
        free_config(config1);
        free_config(config2);
        return 1;
    /* not splited */
    } else if ((diff1  == 0) && (diff2 == 0)) {
        free_config(config1);
        free_config(config2);
        return 2;
    /* connected */
    } else {
        if (diff1 == 0) {
            *dE = energy2 - initial_energy;
            for (i = 0; i < final->tot_num; ++i) {
                final->pos[i * 3 + 0] = config2->pos[i * 3 + 0];
                final->pos[i * 3 + 1] = config2->pos[i * 3 + 1];
                final->pos[i * 3 + 2] = config2->pos[i * 3 + 2];
            }
        } else {
            *dE = energy1 - initial_energy;
            for (i = 0; i < final->tot_num; ++i) {
                final->pos[i * 3 + 0] = config1->pos[i * 3 + 0];
                final->pos[i * 3 + 1] = config1->pos[i * 3 + 1];
                final->pos[i * 3 + 2] = config1->pos[i * 3 + 2];
            }
        }
        free_config(config1);
        free_config(config2);
        return 0;
    }
}
