#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef LMP
#include "lmp_calculator.h"
#endif
#ifdef VASP
#include "vasp_calculator.h"
#endif
#include "alg_utils.h"
#include "config.h"
#include "kappa_dimer.h"
#include "sps_utils.h"


static double *projected_force(double *force0, double *eigenmode,
                               double kappa, int n)
{
    int i;
    double beta = 5.0;
    double gamma1 = 2.0 / (1.0 + exp(beta * kappa)) - 1.0;
    double gamma2 = 1.0 - 1.0 / (1.0 + exp(beta * kappa));
    double *parallel_force = parallel_vector(force0, eigenmode, n);
    double *perpendicular_force = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        perpendicular_force[i * 3 + 0] = force0[i * 3 + 0] - parallel_force[i * 3 + 0];
        perpendicular_force[i * 3 + 1] = force0[i * 3 + 1] - parallel_force[i * 3 + 1];
        perpendicular_force[i * 3 + 2] = force0[i * 3 + 2] - parallel_force[i * 3 + 2];
    }
    double *output = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = gamma2 * perpendicular_force[i * 3 + 0]
                          - gamma1 * parallel_force[i * 3 + 0];
        output[i * 3 + 1] = gamma2 * perpendicular_force[i * 3 + 1]
                          - gamma1 * parallel_force[i * 3 + 1];
        output[i * 3 + 2] = gamma2 * perpendicular_force[i * 3 + 2]
                          - gamma1 * parallel_force[i * 3 + 2];
    }
    free(parallel_force);
    free(perpendicular_force);
    return output;
}


static void rotate(Config *config0, double energy0, double *force0,
                   Input *input, int disp_num, int *disp_list,
                   double *eigenmode, int count, int dimer_step, MPI_Comm comm)
{
    int i, j, rank, size;
    double magnitude, cmin;
    char filename[128];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < input->max_num_rot; ++i) {
        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (j = 0; j < disp_num; ++j) {
            config1->pos[disp_list[j] * 3 + 0] += input->disp_dist
                                                * eigenmode[j * 3 + 0];
            config1->pos[disp_list[j] * 3 + 1] += input->disp_dist 
                                                * eigenmode[j * 3 + 1];
            config1->pos[disp_list[j] * 3 + 2] += input->disp_dist 
                                                * eigenmode[j * 3 + 2];
        }
        oneshot_disp(config1, input, &energy1, force1, disp_num, disp_list, comm);
        free_config(config1);
        for (j = 0; j < disp_num; ++j) {
            force2[j * 3 + 0] = 2 * force0[j * 3 + 0] - force1[j * 3 + 0];
            force2[j * 3 + 1] = 2 * force0[j * 3 + 1] - force1[j * 3 + 1];
            force2[j * 3 + 2] = 2 * force0[j * 3 + 2] - force1[j * 3 + 2];
        }
        double *f_rot_A = get_rot_force(input, force1, force2, eigenmode, disp_num);
        /* no rotation */
        if (norm(f_rot_A, disp_num) < input->f_rot_min) {
            if (local_rank == 0) {
                char filename[128];
                sprintf(filename, "%s/SPS_%d.log",
                        input->output_dir, count);
                FILE *fp = fopen(filename, "a");
                fprintf(fp, " %8d   %8d   %16f   ---------   ---------   %9f\n",
                        dimer_step, i, energy0, norm(f_rot_A, disp_num));
                fclose(fp);
            }
            free(f_rot_A);
            break;
        }
        double *rot_unit_A = normalize(f_rot_A, disp_num);
        /* curvature */
        double *dforce = (double *)malloc(sizeof(double) * disp_num * 3);
        double *n_A = (double *)malloc(sizeof(double) * disp_num * 3);
        for (j = 0; j < disp_num; ++j) {
            n_A[j * 3 + 0] = eigenmode[j * 3 + 0];
            n_A[j * 3 + 1] = eigenmode[j * 3 + 1];
            n_A[j * 3 + 2] = eigenmode[j * 3 + 2];
            dforce[j * 3 + 0] = force2[j * 3 + 0] - force1[j * 3 + 0];
            dforce[j * 3 + 1] = force2[j * 3 + 1] - force1[j * 3 + 1];
            dforce[j * 3 + 2] = force2[j * 3 + 2] - force1[j * 3 + 2];
        }
        magnitude = dot(dforce, eigenmode, disp_num);
        double c0 = magnitude / (2.0 * input->disp_dist);
        magnitude = dot(dforce, rot_unit_A, disp_num);
        double c0d = magnitude / input->disp_dist;
        /* trial rotation */
        double *n_B, *rot_unit_B;
        rotate_vector(n_A, rot_unit_A, &n_B, &rot_unit_B,
                      disp_num, 3.1415926535897932384626 * 0.25);
        Config *trial_config1 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config1, config0);
        for (j = 0; j < disp_num; ++j) {
            trial_config1->pos[disp_list[j] * 3 + 0] += n_B[j * 3 + 0]
                                                      * input->disp_dist;
            trial_config1->pos[disp_list[j] * 3 + 1] += n_B[j * 3 + 1]
                                                      * input->disp_dist;
            trial_config1->pos[disp_list[j] * 3 + 2] += n_B[j * 3 + 2]
                                                      * input->disp_dist;
        } 
        /* derivative of curvature */
        oneshot_disp(trial_config1, input, &energy1, force1,
                     disp_num, disp_list, comm);
        free_config(trial_config1);
        for (j = 0; j < disp_num; ++j) {
            force2[j * 3 + 0] = 2 * force0[j * 3 + 0] - force1[j * 3 + 0];
            force2[j * 3 + 1] = 2 * force0[j * 3 + 1] - force1[j * 3 + 1];
            force2[j * 3 + 2] = 2 * force0[j * 3 + 2] - force1[j * 3 + 2];
            dforce[j * 3 + 0] = force2[j * 3 + 0] - force1[j * 3 + 0];
            dforce[j * 3 + 1] = force2[j * 3 + 1] - force1[j * 3 + 1];
            dforce[j * 3 + 2] = force2[j * 3 + 2] - force1[j * 3 + 2];
        }
        magnitude = dot(dforce, rot_unit_B, disp_num);
        double c1d = magnitude / input->disp_dist;
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
                      disp_num, rotangle);
        for (j = 0; j < disp_num; ++j) {
            eigenmode[j * 3 + 0] = new_eigenmode[j * 3 + 0];
            eigenmode[j * 3 + 1] = new_eigenmode[j * 3 + 1];
            eigenmode[j * 3 + 2] = new_eigenmode[j * 3 + 2];
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
            sprintf(filename, "%s/SPS_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            fprintf(fp, " %8d   %8d   %16f   %9f   %9f   %9f\n",
                    dimer_step, i + 1, energy0, cmin,
                    rotangle * 180 / 3.1415926535897932384626,
                    norm(f_rot_A, disp_num));
            fclose(fp);
        }
        if (norm(f_rot_A, disp_num) < input->f_rot_max) {
            free(f_rot_A);
            break;
        }
        free(f_rot_A);
    }
    free(force1);
    free(force2);
}


static double constrained_rotate(Config *config0, double *force0,
                                 Input *input, int disp_num, int *disp_list,
                                 double *eigenmode, MPI_Comm comm)
{
    int i, j;
    double magnitude, cmin;
    double kappa = 0.0;
    double *tmp_eigenmode, *new_eigenmode;

    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *unit_force0 = normalize(force0, disp_num);
    for (i = 0; i < input->max_num_rot; ++i) {
        /* let eigenmode normal to force */
        tmp_eigenmode = perpendicular_vector(eigenmode, unit_force0, disp_num);
        new_eigenmode = normalize(tmp_eigenmode, disp_num);
        for (j = 0; j < disp_num; ++j) {
            eigenmode[j * 3 + 0] = new_eigenmode[j * 3 + 0];
            eigenmode[j * 3 + 1] = new_eigenmode[j * 3 + 1];
            eigenmode[j * 3 + 2] = new_eigenmode[j * 3 + 2];
        }
        free(tmp_eigenmode);
        free(new_eigenmode);

        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (j = 0; j < disp_num; ++j) {
            config1->pos[disp_list[j] * 3 + 0] += input->disp_dist
                                                * eigenmode[j * 3 + 0];
            config1->pos[disp_list[j] * 3 + 1] += input->disp_dist 
                                                * eigenmode[j * 3 + 1];
            config1->pos[disp_list[j] * 3 + 2] += input->disp_dist 
                                                * eigenmode[j * 3 + 2];
        }
        oneshot_disp(config1, input, &energy1, force1, disp_num, disp_list, comm);
        free_config(config1);
        for (j = 0; j < disp_num; ++j) {
            force2[j * 3 + 0] = 2 * force0[j * 3 + 0] - force1[j * 3 + 0];
            force2[j * 3 + 1] = 2 * force0[j * 3 + 1] - force1[j * 3 + 1];
            force2[j * 3 + 2] = 2 * force0[j * 3 + 2] - force1[j * 3 + 2];
        }
        double *f_rot_A = get_rot_force(input, force1, force2, eigenmode, disp_num);
        /* let rotational force normal to force */
        double *new_f_rot_A = perpendicular_vector(f_rot_A, unit_force0, disp_num);
        for (j = 0; j < disp_num; ++j) {
            f_rot_A[j * 3 + 0] = new_f_rot_A[j * 3 + 0];
            f_rot_A[j * 3 + 1] = new_f_rot_A[j * 3 + 1];
            f_rot_A[j * 3 + 2] = new_f_rot_A[j * 3 + 2];
        }
        free(new_f_rot_A);

        /* no rotation */
        if (norm(f_rot_A, disp_num) < input->f_rot_min) {
            free(f_rot_A);
            break;
        }
        double *rot_unit_A = normalize(f_rot_A, disp_num);
        /* curvature */
        double *dforce = (double *)malloc(sizeof(double) * disp_num * 3);
        double *n_A = (double *)malloc(sizeof(double) * disp_num * 3);
        for (j = 0; j < disp_num; ++j) {
            n_A[j * 3 + 0] = eigenmode[j * 3 + 0];
            n_A[j * 3 + 1] = eigenmode[j * 3 + 1];
            n_A[j * 3 + 2] = eigenmode[j * 3 + 2];
            dforce[j * 3 + 0] = force2[j * 3 + 0] - force1[j * 3 + 0];
            dforce[j * 3 + 1] = force2[j * 3 + 1] - force1[j * 3 + 1];
            dforce[j * 3 + 2] = force2[j * 3 + 2] - force1[j * 3 + 2];
        }
        magnitude = dot(dforce, eigenmode, disp_num);
        double c0 = magnitude / (2.0 * input->disp_dist);
        magnitude = dot(dforce, rot_unit_A, disp_num);
        double c0d = magnitude / input->disp_dist;
        /* trial rotation */
        double *n_B, *rot_unit_B;
        rotate_vector(n_A, rot_unit_A, &n_B, &rot_unit_B,
                      disp_num, 3.1415926535897932384626 * 0.25);
        Config *trial_config1 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config1, config0);
        for (j = 0; j < disp_num; ++j) {
            trial_config1->pos[disp_list[j] * 3 + 0] += n_B[j * 3 + 0]
                                                      * input->disp_dist;
            trial_config1->pos[disp_list[j] * 3 + 1] += n_B[j * 3 + 1]
                                                      * input->disp_dist;
            trial_config1->pos[disp_list[j] * 3 + 2] += n_B[j * 3 + 2]
                                                      * input->disp_dist;
        } 
        /* derivative of curvature */
        oneshot_disp(trial_config1, input, &energy1, force1, disp_num, disp_list, comm);
        free_config(trial_config1);
        for (j = 0; j < disp_num; ++j) {
            force2[j * 3 + 0] = 2 * force0[j * 3 + 0] - force1[j * 3 + 0];
            force2[j * 3 + 1] = 2 * force0[j * 3 + 1] - force1[j * 3 + 1];
            force2[j * 3 + 2] = 2 * force0[j * 3 + 2] - force1[j * 3 + 2];
            dforce[j * 3 + 0] = force2[j * 3 + 0] - force1[j * 3 + 0];
            dforce[j * 3 + 1] = force2[j * 3 + 1] - force1[j * 3 + 1];
            dforce[j * 3 + 2] = force2[j * 3 + 2] - force1[j * 3 + 2];
        }
        magnitude = dot(dforce, rot_unit_B, disp_num);
        double c1d = magnitude / input->disp_dist;
        /* fourier coefficients */
        double a1 = (c0d * cos(2 * 3.1415926535897932384626 * 0.25) - c1d)
                  / (2 * sin(2 * 3.1415926535897932384626 * 0.25));
        double b1 = 0.5 * c0d;
        double a0 = 2 * (c0 - a1);
        /* rotational angle */
        double rotangle = 0.5 * atan(b1 / a1);
        cmin = 0.5 * a0 + a1 * cos(2 * rotangle) + b1 * sin(2 * rotangle);
        kappa = -cmin / norm(force0, disp_num);
        if (c0 < cmin) {
            rotangle += 3.1415926535897932384626 * 0.5;
        }
        double *tmp_force;
        rotate_vector(n_A, rot_unit_A, &new_eigenmode, &tmp_force,
                      disp_num, rotangle);
        for (j = 0; j < disp_num; ++j) {
            eigenmode[j * 3 + 0] = new_eigenmode[j * 3 + 0];
            eigenmode[j * 3 + 1] = new_eigenmode[j * 3 + 1];
            eigenmode[j * 3 + 2] = new_eigenmode[j * 3 + 2];
        }
        free(n_A);
        free(n_B);
        free(dforce);
        free(rot_unit_A);
        free(rot_unit_B);
        free(tmp_force);
        free(new_eigenmode);
        if (norm(f_rot_A, disp_num) < input->f_rot_max) {
            free(f_rot_A);
            break;
        }
        free(f_rot_A);
    }
    free(force1);
    free(force2);
    free(unit_force0);
    return kappa;
}


static void translate(Config *config0, double *force0, Input *input,
                      int disp_num, int *disp_list, double *eigenmode,
                      double *direction_old, double *cg_direction,
                      int dimer_step, double kappa, MPI_Comm comm)
{
    int i;
    double magnitude;
    char filename[128];
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    Config *config1 = (Config *)malloc(sizeof(Config));
    copy_config(config1, config0);
    for (i = 0; i < disp_num; ++i) {
        config1->pos[disp_list[i] * 3 + 0] += input->disp_dist
                                            * eigenmode[i * 3 + 0];
        config1->pos[disp_list[i] * 3 + 1] += input->disp_dist 
                                            * eigenmode[i * 3 + 1];
        config1->pos[disp_list[i] * 3 + 2] += input->disp_dist 
                                            * eigenmode[i * 3 + 2];
    }
    /* curvature */
    oneshot_disp(config1, input, &energy1, force1, disp_num, disp_list, comm);
    free_config(config1);
    for (i = 0; i < disp_num; ++i) {
        force2[i * 3 + 0] = 2 * force0[i * 3 + 0] - force1[i * 3 + 0];
        force2[i * 3 + 1] = 2 * force0[i * 3 + 1] - force1[i * 3 + 1];
        force2[i * 3 + 2] = 2 * force0[i * 3 + 2] - force1[i * 3 + 2];
    }
    double *dforce = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        dforce[i * 3 + 0] = force2[i * 3 + 0] - force1[i * 3 + 0];
        dforce[i * 3 + 1] = force2[i * 3 + 1] - force1[i * 3 + 1];
        dforce[i * 3 + 2] = force2[i * 3 + 2] - force1[i * 3 + 2];
    }
    magnitude = dot(dforce, eigenmode, disp_num);
    double curvature = magnitude / (2.0 * input->disp_dist);
    /* projected force */
    double *f0p = projected_force(force0, eigenmode, kappa, disp_num);
    /* cg_direction */
    get_cg_direction(f0p, direction_old, cg_direction, disp_num);
    double *direction = normalize(cg_direction, disp_num);
    /* step */
    double *step = (double *)malloc(sizeof(double) * disp_num * 3);
    if (curvature > 0) {
        for (i = 0; i < disp_num; ++i) {
            step[i * 3 + 0] = direction[i * 3 + 0] * input->max_move;
            step[i * 3 + 1] = direction[i * 3 + 1] * input->max_move;
            step[i * 3 + 2] = direction[i * 3 + 2] * input->max_move;
        }
    } else {
        Config *trial_config0 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config0, config0);
        for (i = 0; i < disp_num; ++i) {
            trial_config0->pos[disp_list[i] * 3 + 0] += direction[i * 3 + 0]
                                                      * input->trial_move;
            trial_config0->pos[disp_list[i] * 3 + 1] += direction[i * 3 + 1]
                                                      * input->trial_move;
            trial_config0->pos[disp_list[i] * 3 + 2] += direction[i * 3 + 2]
                                                      * input->trial_move;
        }
        double trial_energy0;
        double *trial_force0 = (double *)malloc(sizeof(double) * disp_num * 3);
        double *tmp_force = (double *)malloc(sizeof(double) * disp_num * 3);
        oneshot_disp(trial_config0, input, &trial_energy0, trial_force0,
                     disp_num, disp_list, comm); 
        double *f0tp = projected_force(trial_force0, eigenmode,
                                       kappa, disp_num);
        for (i = 0; i < disp_num; ++i) {
            tmp_force[i * 3 + 0] = f0tp[i * 3 + 0] + f0p[i * 3 + 0];
            tmp_force[i * 3 + 1] = f0tp[i * 3 + 1] + f0p[i * 3 + 1];
            tmp_force[i * 3 + 2] = f0tp[i * 3 + 2] + f0p[i * 3 + 2];
        }
        double F = dot(tmp_force, direction, disp_num) / 2.0;
        for (i = 0; i < disp_num; ++i) {
            tmp_force[i * 3 + 0] = f0tp[i * 3 + 0] - f0p[i * 3 + 0];
            tmp_force[i * 3 + 1] = f0tp[i * 3 + 1] - f0p[i * 3 + 1];
            tmp_force[i * 3 + 2] = f0tp[i * 3 + 2] - f0p[i * 3 + 2];
        }
        double C = dot(tmp_force, direction, disp_num) / input->trial_move;
        double coeff = -F / C + input->trial_move * 0.5;
        for (i = 0; i < disp_num; ++i) {
            step[i * 3 + 0] = coeff * direction[i * 3 + 0];
            step[i * 3 + 1] = coeff * direction[i * 3 + 1];
            step[i * 3 + 2] = coeff * direction[i * 3 + 2];
        }
        free(trial_force0);
        free(tmp_force);
        free(f0tp);
        free_config(trial_config0);
        if (norm(step, disp_num) > input->max_move) {
            for (i = 0; i < disp_num; ++i) {
                step[i * 3 + 0] = direction[i * 3 + 0] * input->max_move;
                step[i * 3 + 1] = direction[i * 3 + 1] * input->max_move;
                step[i * 3 + 2] = direction[i * 3 + 2] * input->max_move;
            }
        }
    }
    /* check nan */
    int nan = 0;
    for (i = 0; i < disp_num; ++i) {
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
        free(f0p);
        free(direction);
        free(step);
        return;
    }

    for (i = 0; i < disp_num; ++i) {
        config0->pos[disp_list[i] * 3 + 0] += step[i * 3 + 0];
        config0->pos[disp_list[i] * 3 + 1] += step[i * 3 + 1];
        config0->pos[disp_list[i] * 3 + 2] += step[i * 3 + 2];
    } 
    free(dforce);
    free(force1);
    free(force2);
    free(f0p);
    free(direction);
    free(step);
}


int kappa_dimer(Config *initial, Config *final, Input *input,
                double *full_eigenmode, int count, int index, double *Ea,
                MPI_Comm comm)
{
    int i, j, rank, size;
    int conv = 1;
    char filename[128];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    /* generate lists */
    int tmp_num;
    int *tmp_list;
    double center[3] = {initial->pos[index * 3 + 0],
                        initial->pos[index * 3 + 1],
                        initial->pos[index * 3 + 2]};
    /* extract list */
    get_sphere_list(initial, input, center, input->acti_cutoff,
                    &tmp_num, &tmp_list, comm);
    int extract_num = 0;
    int *extract_list = (int *)malloc(sizeof(int) * initial->tot_num);
    for (i = 0; i < tmp_num; ++i) {
        if (initial->fix[tmp_list[i]] == 0) {
            extract_list[extract_num] = tmp_list[i];
            extract_num++;
        }
    }
    free(tmp_list);
    /* update list */
    get_sphere_list(initial, input, center, 2 * input->pair_cutoff,
                    &tmp_num, &tmp_list, comm);
    int update_num = 0;
    int *update_list = (int *)malloc(sizeof(int) * initial->tot_num);
    int inner_num = 0;
    int *inner_mask = (int *)calloc(initial->tot_num, sizeof(int));
    for (i = 0; i < tmp_num; ++i) {
        update_list[update_num] = tmp_list[i];
        update_num++;
        inner_mask[tmp_list[i]] = 1;
    }
    free(tmp_list);
    /* sphere cut */
    int outer_num = 0;
    int *outer_list = (int *)malloc(sizeof(int) * initial->tot_num);
    for (i = initial->tot_num - 1; i >= 0; --i) {
        if (inner_mask[i] == 0) {
            outer_list[outer_num] = i;
            outer_num++;
        }
    }
    for (i = 0; i < outer_num; ++i) {
        extract_atom(initial, outer_list[i]);
    }
    free(inner_mask);
    free(outer_list);
    /* starting configuration */
    Config *config0 = (Config *)malloc(sizeof(Config));
    copy_config(config0, initial);
    /* displacement list */
    get_sphere_list(config0, input, center, input->acti_cutoff,
                    &tmp_num, &tmp_list, comm);
    int disp_num = 0;
    int *disp_list = (int *)malloc(sizeof(int) * config0->tot_num);
    for (i = 0; i < tmp_num; ++i) {
        if (config0->fix[tmp_list[i]] == 0) {
            disp_list[disp_num] = tmp_list[i];
            disp_num++;
        }
    }
    free(tmp_list);

    /* eigenmode */
    if (full_eigenmode == NULL) {
        full_eigenmode = get_eigenmode(input, final->tot_num, comm); 
    }
    double *tmp_eigenmode = (double *)malloc(sizeof(double) * extract_num * 3);
    for (i = 0; i < extract_num; ++i) {
        tmp_eigenmode[i * 3 + 0] = full_eigenmode[extract_list[i] * 3 + 0];
        tmp_eigenmode[i * 3 + 1] = full_eigenmode[extract_list[i] * 3 + 1];
        tmp_eigenmode[i * 3 + 2] = full_eigenmode[extract_list[i] * 3 + 2];
    }
    memset(full_eigenmode, 0, sizeof(double) * final->tot_num * 3);

    /* initial perturbation */
    if (input->init_disp > 0) {
        get_sphere_list(config0, input, center, input->disp_cutoff,
                        &tmp_num, &tmp_list, comm);
        for (i = 0; i < disp_num; ++i) {
            for (j = 0; j < tmp_num; ++j) {
                if (disp_list[i] == tmp_list[j]) {
                    config0->pos[tmp_list[j] * 3 + 0] += tmp_eigenmode[i * 3 + 0];
                    config0->pos[tmp_list[j] * 3 + 1] += tmp_eigenmode[i * 3 + 1];
                    config0->pos[tmp_list[j] * 3 + 2] += tmp_eigenmode[i * 3 + 2];
                    break;
                }
            }
        }
        free(tmp_list);
    }
    /* normalize */
    double *eigenmode = normalize(tmp_eigenmode, extract_num);
    free(tmp_eigenmode);
    double *init_direction = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        init_direction[i * 3 + 0] = eigenmode[i * 3 + 0];
        init_direction[i * 3 + 1] = eigenmode[i * 3 + 1];
        init_direction[i * 3 + 2] = eigenmode[i * 3 + 2];
    }

    /* cg optimization */
    double *direction_old = (double *)calloc(disp_num * 3, sizeof(double));
    double *cg_direction = (double *)calloc(disp_num * 3, sizeof(double));

    /* run */
    double kappa;
    int dimer_step;
    if (local_rank == 0) {
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "w");
        fputs("----------------------------------------------------------------------------\n", fp);
        fputs(" Opt step   Rot step   Potential energy   Curvature   Rot angle   Rot force\n", fp);
        fputs("----------------------------------------------------------------------------\n", fp);
        fclose(fp);
        sprintf(filename, "%s/SPS_%d.XDATCAR",
                input->output_dir, count);
        write_config(config0, filename, "w");
    }

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
    for (dimer_step = 1; dimer_step <= input->max_num_tls; ++dimer_step) {
        rotate(config0, energy0, force0, input, disp_num, disp_list,
               eigenmode, count, dimer_step, comm);
        /* test */
        if ((local_rank == 0) && (input->write_mode)) {
            for (i = 0; i < disp_num; ++i) {
                full_eigenmode[extract_list[i] * 3 + 0] = eigenmode[i * 3 + 0];
                full_eigenmode[extract_list[i] * 3 + 1] = eigenmode[i * 3 + 1];
                full_eigenmode[extract_list[i] * 3 + 2] = eigenmode[i * 3 + 2];
            }
            sprintf(filename, "%s/%d_%d.MODECAR",
                    input->output_dir, count, dimer_step);
            FILE *fp = fopen(filename, "w");
            for (i = 0; i < final->tot_num; ++i) {
                fprintf(fp, "%f %f %f\n",
                        full_eigenmode[i * 3 + 0],
                        full_eigenmode[i * 3 + 1],
                        full_eigenmode[i * 3 + 2]);
            }
            fclose(fp);
        }
        /* kappa-dimer */
        double *tmp_eigenmode = (double *)malloc(sizeof(double) * disp_num * 3);
        for (i = 0; i < disp_num; ++i) {
            tmp_eigenmode[i * 3 + 0] = eigenmode[i * 3 + 0];
            tmp_eigenmode[i * 3 + 1] = eigenmode[i * 3 + 1];
            tmp_eigenmode[i * 3 + 2] = eigenmode[i * 3 + 2];
        }
        kappa = constrained_rotate(config0, force0, input, disp_num, disp_list,
                                   tmp_eigenmode, comm);
        free(tmp_eigenmode);
        translate(config0, force0, input, disp_num, disp_list, eigenmode,
                  direction_old, cg_direction, dimer_step, kappa, comm);
        oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
        double fmax = 0.0;
        for (i = 0; i < disp_num; ++i) {
            double tmpf = force0[i * 3 + 0] * force0[i * 3 + 0]
                        + force0[i * 3 + 1] * force0[i * 3 + 1]
                        + force0[i * 3 + 2] * force0[i * 3 + 2];
            tmpf = sqrt(tmpf);
            if (tmpf > fmax) {
                fmax = tmpf;
            } 
        }
        /* trajectory */
        if (local_rank == 0) {
            sprintf(filename, "%s/SPS_%d.XDATCAR",
                    input->output_dir, count);
            write_config(config0, filename, "a");
        }
        if (fmax < input->f_tol) {
            conv = 0;
            break;
        }
    }
    free(force0);
    free(direction_old);
    free(cg_direction);
    if (local_rank == 0) {
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "a");
        fputs("----------------------------------------------------------------------------\n", fp);
        if (conv > 0) {
            fputs(" Saddle state: not converged\n", fp);
        }
        fclose(fp);
    }
    if (conv > 0) {
        free_config(config0);
        free(extract_list);
        free(update_list);
        free(disp_list);
        free(full_eigenmode);
        free(eigenmode);
        return conv;
    }

    /* saddle update */
    for (i = 0; i < disp_num; ++i) {
        /* final stores saddle */
        final->pos[extract_list[i] * 3 + 0] = config0->pos[disp_list[i] * 3 + 0];
        final->pos[extract_list[i] * 3 + 1] = config0->pos[disp_list[i] * 3 + 1];
        final->pos[extract_list[i] * 3 + 2] = config0->pos[disp_list[i] * 3 + 2];
        full_eigenmode[extract_list[i] * 3 + 0] = eigenmode[i * 3 + 0];
        full_eigenmode[extract_list[i] * 3 + 1] = eigenmode[i * 3 + 1];
        full_eigenmode[extract_list[i] * 3 + 2] = eigenmode[i * 3 + 2];
    }
    if (local_rank == 0) {
        sprintf(filename, "%s/Saddle_%d_%d.POSCAR",
                input->output_dir, count, index);
        write_config(final, filename, "w");
        sprintf(filename, "%s/%d.MODECAR",
                input->output_dir, count);
        FILE *fp = fopen(filename, "w");
        for (i = 0; i < final->tot_num; ++i) {
            fprintf(fp, "%f %f %f\n",
                    full_eigenmode[i * 3 + 0],
                    full_eigenmode[i * 3 + 1],
                    full_eigenmode[i * 3 + 2]);
        }
        fclose(fp);
    }

    /* split */
    conv = split_configs(initial, final, config0, input,
                         Ea, eigenmode, count, index,
                         update_num, update_list, disp_num, disp_list, comm);

    free_config(config0);
    free(extract_list);
    free(update_list);
    free(disp_list);
    free(full_eigenmode);
    free(eigenmode);
    return conv;
}
