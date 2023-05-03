#include "dimer.h"
#include "calculator.h"
#include "config.h"
#include "linalg.h"
#include "utils.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


static double *projected_force(double *force0, double *eigenmode,
                               double curvature, int n)
{
    int i;
    double *parallel_force = parallel_vector(force0, eigenmode, n);
    double *output = (double *)malloc(sizeof(double) * n * 3);
    if (curvature > 0) {
        for (i = 0; i < n; ++i) {
            output[i * 3 + 0] = -parallel_force[i * 3 + 0];
            output[i * 3 + 1] = -parallel_force[i * 3 + 1];
            output[i * 3 + 2] = -parallel_force[i * 3 + 2];
        }
    } else {
        for (i = 0; i < n; ++i) {
            output[i * 3 + 0] = force0[i * 3 + 0] - 2 * parallel_force[i * 3 + 0];
            output[i * 3 + 1] = force0[i * 3 + 1] - 2 * parallel_force[i * 3 + 1];
            output[i * 3 + 2] = force0[i * 3 + 2] - 2 * parallel_force[i * 3 + 2];
        }
    }
    free(parallel_force);
    return output;
}


static void rotate(Config *config0, Input *input,
                   int active_num, int *active_list,
                   double *curvature, double *eigenmode,
                   double energy0, double *force0,
                   int count, int sps_step, MPI_Comm comm)
{
    int i, rank;
    double magnitude, cmin;
    char filename[128];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * active_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * active_num * 3);
    double *full_force = (double *)malloc(sizeof(double) * config0->tot_num * 3);

    int rot_step;
    for (rot_step = 0; rot_step < input->max_num_rot; ++rot_step) {
        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (i = 0; i < active_num; ++i) {
            config1->pos[active_list[i] * 3 + 0] += input->finite_diff
                                                  * eigenmode[i * 3 + 0];
            config1->pos[active_list[i] * 3 + 1] += input->finite_diff
                                                  * eigenmode[i * 3 + 1];
            config1->pos[active_list[i] * 3 + 2] += input->finite_diff
                                                  * eigenmode[i * 3 + 2];
        }
        oneshot(config1, input, &energy1, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force1[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force1[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force1[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        free_config(config1);
        for (i = 0; i < active_num; ++i) {
            force2[i * 3 + 0] = 2 * force0[i * 3 + 0] - force1[i * 3 + 0];
            force2[i * 3 + 1] = 2 * force0[i * 3 + 1] - force1[i * 3 + 1];
            force2[i * 3 + 2] = 2 * force0[i * 3 + 2] - force1[i * 3 + 2];
        }
        double *f_rot_A = get_rot_force(input, force1, force2, eigenmode, active_num);
        /* no rotation */
        if (norm(f_rot_A, active_num) < input->f_rot_min) {
            if (local_rank == 0) {
                char filename[128];
                sprintf(filename, "./%d.log", count);
                FILE *fp = fopen(filename, "a");
                fprintf(fp, " %10d   %8d   %16f   ---------   %14f\n",
                        sps_step, rot_step, energy0, norm(f_rot_A, active_num));
                fclose(fp);
            }
            free(f_rot_A);
            break;
        }
        double *rot_unit_A = normalize(f_rot_A, active_num);
        /* curvature */
        double *dforce = (double *)malloc(sizeof(double) * active_num * 3);
        double *n_A = (double *)malloc(sizeof(double) * active_num * 3);
        for (i = 0; i < active_num; ++i) {
            n_A[i * 3 + 0] = eigenmode[i * 3 + 0];
            n_A[i * 3 + 1] = eigenmode[i * 3 + 1];
            n_A[i * 3 + 2] = eigenmode[i * 3 + 2];
            dforce[i * 3 + 0] = force2[i * 3 + 0] - force1[i * 3 + 0];
            dforce[i * 3 + 1] = force2[i * 3 + 1] - force1[i * 3 + 1];
            dforce[i * 3 + 2] = force2[i * 3 + 2] - force1[i * 3 + 2];
        }
        magnitude = dot(dforce, eigenmode, active_num);
        double c0 = magnitude / (2.0 * input->finite_diff);
        *curvature = c0;
        magnitude = dot(dforce, rot_unit_A, active_num);
        double c0d = magnitude / input->finite_diff;
        /* trial rotation */
        double *n_B, *rot_unit_B;
        rotate_vector(n_A, rot_unit_A, &n_B, &rot_unit_B,
                      active_num, 3.1415926535897932384626 * 0.25);
        Config *trial_config1 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config1, config0);
        for (i = 0; i < active_num; ++i) {
            trial_config1->pos[active_list[i] * 3 + 0] += n_B[i * 3 + 0]
                                                        * input->finite_diff;
            trial_config1->pos[active_list[i] * 3 + 1] += n_B[i * 3 + 1]
                                                        * input->finite_diff;
            trial_config1->pos[active_list[i] * 3 + 2] += n_B[i * 3 + 2]
                                                        * input->finite_diff;
        } 
        /* derivative of curvature */
        oneshot(trial_config1, input, &energy1, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force1[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force1[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force1[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        free_config(trial_config1);
        for (i = 0; i < active_num; ++i) {
            force2[i * 3 + 0] = 2 * force0[i * 3 + 0] - force1[i * 3 + 0];
            force2[i * 3 + 1] = 2 * force0[i * 3 + 1] - force1[i * 3 + 1];
            force2[i * 3 + 2] = 2 * force0[i * 3 + 2] - force1[i * 3 + 2];
            dforce[i * 3 + 0] = force2[i * 3 + 0] - force1[i * 3 + 0];
            dforce[i * 3 + 1] = force2[i * 3 + 1] - force1[i * 3 + 1];
            dforce[i * 3 + 2] = force2[i * 3 + 2] - force1[i * 3 + 2];
        }
        magnitude = dot(dforce, rot_unit_B, active_num);
        double c1d = magnitude / input->finite_diff;
        /* fourier coefficients */
        double a1 = (c0d * cos(2 * 3.1415926535897932384626 * 0.25) - c1d)
                  / (2 * sin(2 * 3.1415926535897932384626 * 0.25));
        double b1 = 0.5 * c0d;
        double a0 = 2 * (c0 - a1);
        /* rotational angle */
        double rotangle = 0.5 * atan(b1 / a1);
        cmin = 0.5 * a0 + a1 * cos(2 * rotangle) + b1 * sin(2 * rotangle);
        if (c0 < cmin) {
            rotangle += 3.1415926535897932384626 * 0.5;
        }
        double *new_eigenmode, *tmp_force;
        rotate_vector(n_A, rot_unit_A, &new_eigenmode, &tmp_force,
                      active_num, rotangle);
        for (i = 0; i < active_num; ++i) {
            eigenmode[i * 3 + 0] = new_eigenmode[i * 3 + 0];
            eigenmode[i * 3 + 1] = new_eigenmode[i * 3 + 1];
            eigenmode[i * 3 + 2] = new_eigenmode[i * 3 + 2];
        }
        free(n_A);
        free(n_B);
        free(dforce);
        free(rot_unit_A);
        free(rot_unit_B);
        free(new_eigenmode);
        free(tmp_force);
        if (local_rank == 0) {
            char filename[128];
            sprintf(filename, "./%d.log", count);
            FILE *fp = fopen(filename, "a");
            fprintf(fp, " %10d   %8d   %16f   %9f   %14f\n",
                    sps_step, rot_step + 1, energy0, cmin,
                    norm(f_rot_A, active_num));
            fclose(fp);
        }
        if (norm(f_rot_A, active_num) < input->f_rot_max) {
            free(f_rot_A);
            break;
        }
        free(f_rot_A);
    }
    free(force1);
    free(force2);
    free(full_force);
}


static void translate(Config *config0, Input *input, int active_num, int *active_list,
                      double *eigenmode, double *force0, double *direction_old,
                      double *cg_direction, int count, int index, int sps_step,
                      MPI_Comm comm)
{
    int i, rank;
    char filename[128];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double magnitude;
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * active_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * active_num * 3);
    double *full_force = (double *)malloc(sizeof(double) * config0->tot_num * 3);
    Config *config1 = (Config *)malloc(sizeof(Config));
    copy_config(config1, config0);
    for (i = 0; i < active_num; ++i) {
        config1->pos[active_list[i] * 3 + 0] += input->finite_diff
                                              * eigenmode[i * 3 + 0];
        config1->pos[active_list[i] * 3 + 1] += input->finite_diff
                                              * eigenmode[i * 3 + 1];
        config1->pos[active_list[i] * 3 + 2] += input->finite_diff
                                              * eigenmode[i * 3 + 2];
    }
    /* curvature */
    oneshot(config1, input, &energy1, full_force, comm);
    for (i = 0; i < active_num; ++i) {
        force1[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
        force1[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
        force1[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
    }
    free_config(config1);
    for (i = 0; i < active_num; ++i) {
        force2[i * 3 + 0] = 2 * force0[i * 3 + 0] - force1[i * 3 + 0];
        force2[i * 3 + 1] = 2 * force0[i * 3 + 1] - force1[i * 3 + 1];
        force2[i * 3 + 2] = 2 * force0[i * 3 + 2] - force1[i * 3 + 2];
    }
    double *dforce = (double *)malloc(sizeof(double) * active_num * 3);
    for (i = 0; i < active_num; ++i) {
        dforce[i * 3 + 0] = force2[i * 3 + 0] - force1[i * 3 + 0];
        dforce[i * 3 + 1] = force2[i * 3 + 1] - force1[i * 3 + 1];
        dforce[i * 3 + 2] = force2[i * 3 + 2] - force1[i * 3 + 2];
    }
    magnitude = dot(dforce, eigenmode, active_num);
    double curvature = magnitude / (2.0 * input->finite_diff);
    /* projected force */
    double *f0p = projected_force(force0, eigenmode, curvature, active_num);
    /* cg_direction */
    get_cg_direction(f0p, direction_old, cg_direction, active_num);
    double *direction = normalize(cg_direction, active_num);
    /* step */
    double *step = (double *)malloc(sizeof(double) * active_num * 3);
    if (curvature > 0) {
        for (i = 0; i < active_num; ++i) {
            step[i * 3 + 0] = direction[i * 3 + 0] * input->max_move;
            step[i * 3 + 1] = direction[i * 3 + 1] * input->max_move;
            step[i * 3 + 2] = direction[i * 3 + 2] * input->max_move;
        }
    } else {
        Config *trial_config0 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config0, config0);
        for (i = 0; i < active_num; ++i) {
            trial_config0->pos[active_list[i] * 3 + 0] += direction[i * 3 + 0]
                                                        * input->trial_move;
            trial_config0->pos[active_list[i] * 3 + 1] += direction[i * 3 + 1]
                                                        * input->trial_move;
            trial_config0->pos[active_list[i] * 3 + 2] += direction[i * 3 + 2]
                                                        * input->trial_move;
        }
        double trial_energy0;
        double *trial_force0 = (double *)malloc(sizeof(double) * active_num * 3);
        double *tmp_force = (double *)malloc(sizeof(double) * active_num * 3);
        oneshot(trial_config0, input, &trial_energy0, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            trial_force0[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            trial_force0[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            trial_force0[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        double *f0tp = projected_force(trial_force0, eigenmode,
                                       curvature, active_num);
        for (i = 0; i < active_num; ++i) {
            tmp_force[i * 3 + 0] = f0tp[i * 3 + 0] + f0p[i * 3 + 0];
            tmp_force[i * 3 + 1] = f0tp[i * 3 + 1] + f0p[i * 3 + 1];
            tmp_force[i * 3 + 2] = f0tp[i * 3 + 2] + f0p[i * 3 + 2];
        }
        double F = dot(tmp_force, direction, active_num) / 2.0;
        for (i = 0; i < active_num; ++i) {
            tmp_force[i * 3 + 0] = f0tp[i * 3 + 0] - f0p[i * 3 + 0];
            tmp_force[i * 3 + 1] = f0tp[i * 3 + 1] - f0p[i * 3 + 1];
            tmp_force[i * 3 + 2] = f0tp[i * 3 + 2] - f0p[i * 3 + 2];
        }
        double C = dot(tmp_force, direction, active_num) / input->trial_move;
        double coeff = -F / C + input->trial_move * 0.5;
        for (i = 0; i < active_num; ++i) {
            step[i * 3 + 0] = coeff * direction[i * 3 + 0];
            step[i * 3 + 1] = coeff * direction[i * 3 + 1];
            step[i * 3 + 2] = coeff * direction[i * 3 + 2];
        }
        free(trial_force0);
        free(tmp_force);
        free(f0tp);
        free_config(trial_config0);
        if (norm(step, active_num) > input->max_move) {
            for (i = 0; i < active_num; ++i) {
                step[i * 3 + 0] = direction[i * 3 + 0] * input->max_move;
                step[i * 3 + 1] = direction[i * 3 + 1] * input->max_move;
                step[i * 3 + 2] = direction[i * 3 + 2] * input->max_move;
            }
        }
    }
    /* check nan */
    int nan = 0;
    for (i = 0; i < active_num; ++i) {
        if (isnan(step[i * 3 + 0]) != 0) {
            nan = 1;
        }
        if (isnan(step[i * 3 + 1]) != 0) {
            nan = 1;
        }
        if (isnan(step[i * 3 + 2]) != 0) {
            nan = 1;
        }
    }
    if (nan > 0) {
        free(dforce);
        free(force1);
        free(force2);
        free(full_force);
        free(f0p);
        free(direction);
        free(step);
        return;
    }

    for (i = 0; i < active_num; ++i) {
        config0->pos[active_list[i] * 3 + 0] += step[i * 3 + 0];
        config0->pos[active_list[i] * 3 + 1] += step[i * 3 + 1];
        config0->pos[active_list[i] * 3 + 2] += step[i * 3 + 2];
    } 
    /* trajectory */
    if (local_rank == 0) {
        char header[128];
        sprintf(header, "%d_%d %d", count, index, sps_step);
        sprintf(filename, "./%d.XDATCAR", count);
        write_config(config0, filename, header, "a");
    }

    free(dforce);
    free(force1);
    free(force2);
    free(full_force);
    free(f0p);
    free(direction);
    free(step);
}


int dimer(Config *initial, Config *saddle, Config *final, Input *input,
          double *full_eigenmode, int count, int index, double *Ea,
          MPI_Comm comm)
{
    int i, j, rank;
    int conv = -1;
    char filename[128];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    /* active region */
    Config *config0 = (Config *)malloc(sizeof(Config));
    copy_config(config0, initial);
    int active_num = 0;
    int *active_list = (int *)malloc(sizeof(int) * config0->tot_num);
    for (i = 0; i < config0->tot_num; ++i) {
        if (config0->fix[i] == 0) {
            active_list[active_num] = i;
            active_num++;
        }
    }

    /* eigenmode */
    if (full_eigenmode == NULL) {
        full_eigenmode = get_eigenmode(input, config0->tot_num, comm);
    }
    double *eigenmode = (double *)calloc(active_num * 3, sizeof(double));
    for (i = 0; i < active_num; ++i) {
        eigenmode[i * 3 + 0] = full_eigenmode[active_list[i] * 3 + 0];
        eigenmode[i * 3 + 1] = full_eigenmode[active_list[i] * 3 + 1];
        eigenmode[i * 3 + 2] = full_eigenmode[active_list[i] * 3 + 2];
    }
    memset(full_eigenmode, 0, sizeof(double) * config0->tot_num * 3);
    double *tmp_eigenmode = normalize(eigenmode, active_num);
    for (i = 0; i < active_num; ++i) {
        eigenmode[i * 3 + 0] = tmp_eigenmode[i * 3 + 0];
        eigenmode[i * 3 + 1] = tmp_eigenmode[i * 3 + 1];
        eigenmode[i * 3 + 2] = tmp_eigenmode[i * 3 + 2];
    }
    free(tmp_eigenmode);

    /* initial */
    int tmp_num;
    int *tmp_list;
    double center[3] = {config0->pos[index * 3 + 0],
                        config0->pos[index * 3 + 1],
                        config0->pos[index * 3 + 2]};
    get_sphere_list(config0, input, center, input->init_cutoff,
                    &tmp_num, &tmp_list, comm);
    if (input->init_disp > 0) {
        double *tmp_init_disp = (double *)calloc(active_num * 3, sizeof(double));
        for (i = 0; i < tmp_num; ++i) {
            if (config0->fix[tmp_list[i]] == 0) {
                tmp_init_disp[i * 3 + 0] = normal_random(0, 1);
                tmp_init_disp[i * 3 + 0] = normal_random(0, 1);
                tmp_init_disp[i * 3 + 0] = normal_random(0, 1);
            }
        }
        double *init_disp = normalize(tmp_init_disp, tmp_num);
        free(tmp_init_disp);
        for (i = 0; i < tmp_num; ++i) {
            config0->pos[tmp_list[i] * 3 + 0] += input->disp_move
                                               * init_disp[i * 3 + 0];
            config0->pos[tmp_list[i] * 3 + 1] += input->disp_move
                                               * init_disp[i * 3 + 1];
            config0->pos[tmp_list[i] * 3 + 2] += input->disp_move
                                               * init_disp[i * 3 + 2];
        }
    }
    free(tmp_list);

    /* cg optimization */
    double *direction_old = (double *)calloc(active_num * 3, sizeof(double));
    double *cg_direction = (double *)calloc(active_num * 3, sizeof(double));

    /* run */
    if (local_rank == 0) {
        sprintf(filename, "./%d.log", count);
        FILE *fp = fopen(filename, "w");
        fprintf(fp, " %d_%d\n", count, index);
        fputs("-----------------------------------------------------------------------\n", fp);
        fputs(" Trans step   Rot step   Potential energy   Curvature   Rotation force\n", fp);
        fputs("-----------------------------------------------------------------------\n", fp);
        fclose(fp);
        char header[128];
        sprintf(header, "%d_%d %d", count, index, 0);
        sprintf(filename, "./%d.XDATCAR", count);
        write_config(config0, filename, header, "w");
    }

    int sps_step;
    double curvature = 1.0;
    double energy0;
    double *force0 = (double *)calloc(config0->tot_num * 3, sizeof(double));
    double *full_force = (double *)malloc(sizeof(double) * config0->tot_num * 3);
    clock_t start = clock();
    for (sps_step = 1; sps_step <= input->max_num_tls; ++sps_step) {
        oneshot(config0, input, &energy0, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force0[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force0[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force0[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        /* rotation */
        rotate(config0, input, active_num, active_list, &curvature,
               eigenmode, energy0, force0, count, sps_step, comm);
        if (curvature < 0) {
            /* force criteria */
            double fmax = 0.0;
            for (i = 0; i < active_num; ++i) {
                double tmpf = force0[i * 3 + 0] * force0[i * 3 + 0]
                            + force0[i * 3 + 1] * force0[i * 3 + 1]
                            + force0[i * 3 + 2] * force0[i * 3 + 2];
                tmpf = sqrt(tmpf);
                if (tmpf > fmax) {
                    fmax = tmpf;
                }
            }
            if (fmax < input->f_tol) {
                conv = 0;
            }
        }
        /* translation */
        translate(config0, input, active_num, active_list, eigenmode, force0,
                  direction_old, cg_direction, count, index, sps_step, comm);
    }
    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    free(force0);
    free(direction_old);
    free(cg_direction);
    if (local_rank == 0) {
        sprintf(filename, "./%d.log", count);
        FILE *fp = fopen(filename, "a");
        fputs("-----------------------------------------------------------------------\n", fp);
        fputs(" State of the saddle   Barrier energy   Reaction energy   Elapsed time\n", fp);
        fputs("-----------------------------------------------------------------------\n", fp);
        if (conv < 0) {
            fprintf(fp, "         Unconverged   --------------   ---------------   %12f\n", time);
            fputs("-----------------------------------------------------------------------\n", fp);
        }
        fclose(fp);
    }
    if (conv < 0) {
        free_config(config0);
        free(active_list);
        free(full_eigenmode);
        free(full_force);
        free(eigenmode);
        return conv;
    }

    /* saddle update */
    for (i = 0; i < active_num; ++i) {
        saddle->pos[active_list[i] * 3 + 0] = config0->pos[active_list[i] * 3 + 0];
        saddle->pos[active_list[i] * 3 + 1] = config0->pos[active_list[i] * 3 + 1];
        saddle->pos[active_list[i] * 3 + 2] = config0->pos[active_list[i] * 3 + 2];
        full_eigenmode[active_list[i] * 3 + 0] = eigenmode[i * 3 + 0];
        full_eigenmode[active_list[i] * 3 + 1] = eigenmode[i * 3 + 1];
        full_eigenmode[active_list[i] * 3 + 2] = eigenmode[i * 3 + 2];
    }

    /* postprocess */
    double dE;
    conv = split_config(initial, saddle, final, input, Ea, &dE,
                        curvature, eigenmode, active_num, active_list,
                        count, index, comm);
    if (local_rank == 0) {
        char filename[128];
        /* modecar */
        sprintf(filename, "./%d.MODECAR", count);
        FILE *fp = fopen(filename, "w");
        fprintf(fp, "%d_%d\n", count, index);
        for (i = 0; i < config0->tot_num; ++i) {
            fprintf(fp, "%f %f %f\n",
                    full_eigenmode[i * 3 + 0],
                    full_eigenmode[i * 3 + 1],
                    full_eigenmode[i * 3 + 2]);
        }
        fclose(fp);
        /* log */
        sprintf(filename, "./%d.log", count);
        fp = fopen(filename, "a");
        if (conv == 0) {
            fprintf(fp, "           Connected   %14f   %15f   %12f\n", *Ea, dE, time);
        } else if (conv == 1) {
            fprintf(fp, "        Disconnected   %14f   ---------------   %12f\n", *Ea, time);
        } else if (conv == 2) {
            fprintf(fp, "         Not splited   %14f   ---------------   %12f\n", *Ea, time);
        } else {
        }
        fputs("-----------------------------------------------------------------------\n", fp);
        fclose(fp);
    }

    free_config(config0);
    free(active_list);
    free(full_eigenmode);
    free(full_force);
    free(eigenmode);
    return conv;
}
