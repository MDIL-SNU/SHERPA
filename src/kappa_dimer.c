#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "config.h"
#include "kappa_dimer.h"
#include "utils.h"


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


static void rotate(Config *config0, Input *input, int disp_num, int *disp_list,
                   double *eigenmode, int count, int dimer_step, MPI_Comm comm)
{
    int i, j, rank, size;
    double magnitude, cmin;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy0, energy1;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
    for (i = 0; i < input->max_num_rot; ++i) {
        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (j = 0; j < disp_num; ++j) {
            config1->pos[disp_list[j] * 3 + 0] += input->dimer_dist
                                                * eigenmode[j * 3 + 0];
            config1->pos[disp_list[j] * 3 + 1] += input->dimer_dist 
                                                * eigenmode[j * 3 + 1];
            config1->pos[disp_list[j] * 3 + 2] += input->dimer_dist 
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
                sprintf(filename, "%s/Dimer_%d.log",
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
        double c0 = magnitude / (2.0 * input->dimer_dist);
        magnitude = dot(dforce, rot_unit_A, disp_num);
        double c0d = magnitude / input->dimer_dist;
        /* trial rotation */
        double *n_B, *rot_unit_B;
        rotate_vector(n_A, rot_unit_A, &n_B, &rot_unit_B,
                      disp_num, input->trial_angle); 
        Config *trial_config1 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config1, config0);
        for (j = 0; j < disp_num; ++j) {
            trial_config1->pos[disp_list[j] * 3 + 0] += n_B[j * 3 + 0]
                                                      * input->dimer_dist;
            trial_config1->pos[disp_list[j] * 3 + 1] += n_B[j * 3 + 1]
                                                      * input->dimer_dist;
            trial_config1->pos[disp_list[j] * 3 + 2] += n_B[j * 3 + 2]
                                                      * input->dimer_dist;
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
        double c1d = magnitude / input->dimer_dist;
        /* fourier coefficients */
        double a1 = (c0d * cos(2 * input->trial_angle) - c1d) 
                  / (2 * sin(2 * input->trial_angle));
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
            sprintf(filename, "%s/Dimer_%d.log",
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
    free(force0);
    free(force1);
    free(force2);
}


static double constrained_rotate(Config *config0, Input *input,
                                 int disp_num, int *disp_list,
                                 double *eigenmode, MPI_Comm comm)
{
    int i, j;
    double magnitude, cmin;
    double kappa = 0.0;
    double *tmp_eigenmode, *new_eigenmode;

    double energy0, energy1;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
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
            config1->pos[disp_list[j] * 3 + 0] += input->dimer_dist
                                                * eigenmode[j * 3 + 0];
            config1->pos[disp_list[j] * 3 + 1] += input->dimer_dist 
                                                * eigenmode[j * 3 + 1];
            config1->pos[disp_list[j] * 3 + 2] += input->dimer_dist 
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
        double c0 = magnitude / (2.0 * input->dimer_dist);
        magnitude = dot(dforce, rot_unit_A, disp_num);
        double c0d = magnitude / input->dimer_dist;
        /* trial rotation */
        double *n_B, *rot_unit_B;
        rotate_vector(n_A, rot_unit_A, &n_B, &rot_unit_B,
                      disp_num, input->trial_angle); 
        Config *trial_config1 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config1, config0);
        for (j = 0; j < disp_num; ++j) {
            trial_config1->pos[disp_list[j] * 3 + 0] += n_B[j * 3 + 0]
                                                      * input->dimer_dist;
            trial_config1->pos[disp_list[j] * 3 + 1] += n_B[j * 3 + 1]
                                                      * input->dimer_dist;
            trial_config1->pos[disp_list[j] * 3 + 2] += n_B[j * 3 + 2]
                                                      * input->dimer_dist;
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
        double c1d = magnitude / input->dimer_dist;
        /* fourier coefficients */
        double a1 = (c0d * cos(2 * input->trial_angle) - c1d) 
                  / (2 * sin(2 * input->trial_angle));
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
    free(force0);
    free(force1);
    free(force2);
    free(unit_force0);
    return kappa;
}


static void translate(Config *config0, Input *input,
                      int disp_num, int *disp_list, double *eigenmode,
                      double *direction_old, double *cg_direction,
                      int dimer_step, double kappa, MPI_Comm comm)
{
    int i;
    double magnitude;
    double energy0, energy1;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
    Config *config1 = (Config *)malloc(sizeof(Config));
    copy_config(config1, config0);
    for (i = 0; i < disp_num; ++i) {
        config1->pos[disp_list[i] * 3 + 0] += input->dimer_dist
                                            * eigenmode[i * 3 + 0];
        config1->pos[disp_list[i] * 3 + 1] += input->dimer_dist 
                                            * eigenmode[i * 3 + 1];
        config1->pos[disp_list[i] * 3 + 2] += input->dimer_dist 
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
    double curvature = magnitude / (2.0 * input->dimer_dist);
    /* projected force */
    double *f0p = projected_force(force0, eigenmode, kappa, disp_num);
    /* cg_direction */
    if (dimer_step == 1) {
        for (i = 0; i < disp_num; ++i) {
            direction_old[i * 3 + 0] = f0p[i * 3 + 0];
            direction_old[i * 3 + 1] = f0p[i * 3 + 1];
            direction_old[i * 3 + 2] = f0p[i * 3 + 2];
            cg_direction[i * 3 + 0] = f0p[i * 3 + 0];
            cg_direction[i * 3 + 1] = f0p[i * 3 + 1];
            cg_direction[i * 3 + 2] = f0p[i * 3 + 2];
        }
    }
    get_cg_direction(f0p, direction_old, cg_direction, disp_num);
    double *direction = normalize(cg_direction, disp_num);
    /* step */
    double *step = (double *)malloc(sizeof(double) * disp_num * 3);
    if (curvature > 0) {
        for (i = 0; i < disp_num; ++i) {
            step[i * 3 + 0] = direction[i * 3 + 0] * input->max_step;
            step[i * 3 + 1] = direction[i * 3 + 1] * input->max_step;
            step[i * 3 + 2] = direction[i * 3 + 2] * input->max_step;
        }
    } else {
        Config *trial_config0 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config0, config0);
        for (i = 0; i < disp_num; ++i) {
            trial_config0->pos[disp_list[i] * 3 + 0] += direction[i * 3 + 0]
                                                      * input->trial_step;
            trial_config0->pos[disp_list[i] * 3 + 1] += direction[i * 3 + 1]
                                                      * input->trial_step;
            trial_config0->pos[disp_list[i] * 3 + 2] += direction[i * 3 + 2]
                                                      * input->trial_step;
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
        double C = dot(tmp_force, direction, disp_num) / input->trial_step;
        double coeff = -F / C + input->trial_step * 0.5;
        for (i = 0; i < disp_num; ++i) {
            step[i * 3 + 0] = coeff * direction[i * 3 + 0];
            step[i * 3 + 1] = coeff * direction[i * 3 + 1];
            step[i * 3 + 2] = coeff * direction[i * 3 + 2];
        }
        free(trial_force0);
        free(tmp_force);
        free(f0tp);
        free_config(trial_config0);
        if (norm(step, disp_num) > input->max_step) {
            for (i = 0; i < disp_num; ++i) {
                step[i * 3 + 0] = direction[i * 3 + 0] * input->max_step;
                step[i * 3 + 1] = direction[i * 3 + 1] * input->max_step;
                step[i * 3 + 2] = direction[i * 3 + 2] * input->max_step;
            }
        }
    }
    for (i = 0; i < disp_num; ++i) {
        config0->pos[disp_list[i] * 3 + 0] += step[i * 3 + 0];
        config0->pos[disp_list[i] * 3 + 1] += step[i * 3 + 1];
        config0->pos[disp_list[i] * 3 + 2] += step[i * 3 + 2];
    } 
    free(dforce);
    free(force0);
    free(force1);
    free(force2);
    free(f0p);
    free(direction);
    free(step);
}


// TODO: orthogonalization
int kappa_dimer(Config *initial, Config *saddle, Config *final, Input *input,
                double *full_eigenmode, int count, int index,
                double *Ea, MPI_Comm comm)
{
    int i, j, rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    /* generate lists */
    double center[3] = {initial->pos[index * 3 + 0],
                        initial->pos[index * 3 + 1],
                        initial->pos[index * 3 + 2]};
    int update_num;
    int extract_num;
    int *update_list;
    int *extract_list;
    gen_list(initial, input, center, &update_num, &update_list,
             &extract_num, &extract_list, comm);

    /* cut far atoms */ 
    cut_sphere(initial, input, update_num, update_list);

    /* set dimer space */
    double del[3];
    int disp_num = 0;
    int *disp_list = (int *)malloc(sizeof(int) * initial->tot_num);
    for (i = 0; i < initial->tot_num; ++i) {
        del[0] = initial->pos[i * 3 + 0] - center[0];
        del[1] = initial->pos[i * 3 + 1] - center[1];
        del[2] = initial->pos[i * 3 + 2] - center[2];
        get_minimum_image(del, initial->boxlo, initial->boxhi,
                          initial->xy, initial->yz, initial->xz);
        double dist = sqrt(del[0] * del[0]
                         + del[1] * del[1] 
                         + del[2] * del[2]);
        if (dist < input->disp_cutoff) {
            disp_list[disp_num] = i;
            disp_num++;
        }
    }

    /* cut far atoms for starting dimer */ 
    Config *config0 = (Config *)malloc(sizeof(Config));
    copy_config(config0, saddle);
    cut_sphere(config0, input, update_num, update_list);

    /* oneshot for first force */
    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);

    /* eigenmode */
    double *tmp_eigenmode = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        tmp_eigenmode[i * 3 + 0] = full_eigenmode[extract_list[i] * 3 + 0];
        tmp_eigenmode[i * 3 + 1] = full_eigenmode[extract_list[i] * 3 + 1];
        tmp_eigenmode[i * 3 + 2] = full_eigenmode[extract_list[i] * 3 + 2];
    }
    double *eigenmode = normalize(tmp_eigenmode, disp_num);

    /* cg optimization */
    double *direction_old = (double *)malloc(sizeof(double) * disp_num * 3);
    double *cg_direction = (double *)malloc(sizeof(double) * disp_num * 3);

    /* run */
    if (local_rank == 0) {
        char filename[128];
        sprintf(filename, "%s/Dimer_%d.XDATCAR",
                input->output_dir, count);
        write_config(config0, filename, "w");
    }
    double fmax, kappa;
    int converge = 0;
    int dimer_step;
    if (local_rank == 0) {
        char filename[128];
        sprintf(filename, "%s/Dimer_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "w");
        fputs(" Opt step   Rot step   Potential energy   Curvature   Rot angle   Rot force\n", fp);
        fclose(fp);
    }

    for (dimer_step = 1; dimer_step <= 1000; ++dimer_step) {
        rotate(config0, input, disp_num, disp_list,
               eigenmode, count, dimer_step, comm);
        /* kappa-dimer */
        for (i = 0; i < disp_num; ++i) {
            tmp_eigenmode[i * 3 + 0] = eigenmode[i * 3 + 0];
            tmp_eigenmode[i * 3 + 1] = eigenmode[i * 3 + 1];
            tmp_eigenmode[i * 3 + 2] = eigenmode[i * 3 + 2];
        }
        kappa = constrained_rotate(config0, input, disp_num, disp_list,
                                   tmp_eigenmode, comm);
        translate(config0, input, disp_num, disp_list, eigenmode,
                  direction_old, cg_direction, dimer_step, kappa, comm);
        oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);     
        fmax = 0.0;
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
            char filename[128];
            sprintf(filename, "%s/Dimer_%d.XDATCAR",
                    input->output_dir, count);
            write_config(config0, filename, "a");
        }
        if (fmax < input->f_tol) {
            converge = 1;
            break;
        }
    }
    free(tmp_eigenmode);
    free(direction_old);
    free(cg_direction);
    if (converge == 0) {
        if (local_rank == 0) {
            char filename[128];
            sprintf(filename, "%s/Dimer_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            fputs("----------------------------------------------------------------------------\n", fp);
            fputs(" Saddle state: not converged\n", fp);
            fclose(fp);
        }
        free_config(config0);
        free(force0);
        free(eigenmode);
        free(update_list);
        free(extract_list);
        return 1;
    }
    /* relax initial structure and barrier energy */
    atom_relax(initial, input, comm);
    oneshot_disp(initial, input, &energy0, force0, disp_num, disp_list, comm);
    double i_energy = energy0;
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
    double ts_energy = energy0;
    free(force0);
    *Ea = ts_energy - i_energy;

    /* saddle update */
    for (i = 0; i < config0->tot_num; ++i) {
        saddle->pos[update_list[i] * 3 + 0] = config0->pos[i * 3 + 0];
        saddle->pos[update_list[i] * 3 + 1] = config0->pos[i * 3 + 1];
        saddle->pos[update_list[i] * 3 + 2] = config0->pos[i * 3 + 2];
    }
    for (i = 0; i < disp_num; ++i) {
        full_eigenmode[extract_list[i] * 3 + 0] = eigenmode[i * 3 + 0];
        full_eigenmode[extract_list[i] * 3 + 0] = eigenmode[i * 3 + 0];
        full_eigenmode[extract_list[i] * 3 + 0] = eigenmode[i * 3 + 0];
    }
    if (local_rank == 0) {
        char filename[128];
        sprintf(filename, "%s/Saddle_%d_%d.POSCAR",
                input->output_dir, count, index);
        write_config(saddle, filename, "w");
        sprintf(filename, "%s/%d.MODECAR",
                input->output_dir, count);
        FILE *fp = fopen(filename, "wb");     
        fwrite(full_eigenmode, sizeof(double), saddle->tot_num * 3, fp);
        fclose(fp);
    }

    /* split */
    Config *config1 = (Config *)malloc(sizeof(Config));
    Config *config2 = (Config *)malloc(sizeof(Config));
    int trial = 1;
    while (1) {
        copy_config(config1, config0);
        copy_config(config2, config0);
        for (j = 0; j < disp_num; ++j) {
            config1->pos[disp_list[j] * 3 + 0] = config0->pos[disp_list[j] * 3 + 0]
                                               + 0.1 * trial * eigenmode[j * 3 + 0];
            config1->pos[disp_list[j] * 3 + 1] = config0->pos[disp_list[j] * 3 + 1]
                                               + 0.1 * trial * eigenmode[j * 3 + 1];
            config1->pos[disp_list[j] * 3 + 2] = config0->pos[disp_list[j] * 3 + 2]
                                               + 0.1 * trial * eigenmode[j * 3 + 2];
            config2->pos[disp_list[j] * 3 + 0] = config0->pos[disp_list[j] * 3 + 0]
                                               - 0.1 * trial * eigenmode[j * 3 + 0];
            config2->pos[disp_list[j] * 3 + 1] = config0->pos[disp_list[j] * 3 + 1]
                                               - 0.1 * trial * eigenmode[j * 3 + 1];
            config2->pos[disp_list[j] * 3 + 2] = config0->pos[disp_list[j] * 3 + 2]
                                               - 0.1 * trial * eigenmode[j * 3 + 2];
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
    free(disp_list);
    free(extract_list);
    free(eigenmode);

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
        free_config(config0);
        free_config(config1);
        free_config(config2);
        free(update_list);
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
        free_config(config0);
        free_config(config1);
        free_config(config2);
        free(update_list);
        return 0;
    }
}
