#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "config.h"
#include "dimer.h"
#include "utils.h"


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


void rotate_vector(double *vec1i, double *vec2i,
                   double **vec1o, double **vec2o,
                   int n, double angle)
{
    int i;
    double cAng = cos(angle);
    double sAng = sin(angle);
    double *tmp_vec1o = (double *)malloc(sizeof(double) * n * 3);
    double *tmp_vec2o = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        tmp_vec1o[i * 3 + 0] = vec1i[i * 3 + 0] * cAng + vec2i[i * 3 + 0] * sAng;
        tmp_vec1o[i * 3 + 1] = vec1i[i * 3 + 1] * cAng + vec2i[i * 3 + 1] * sAng;
        tmp_vec1o[i * 3 + 2] = vec1i[i * 3 + 2] * cAng + vec2i[i * 3 + 2] * sAng;
        tmp_vec2o[i * 3 + 0] = vec2i[i * 3 + 0] * cAng - vec1i[i * 3 + 0] * sAng;
        tmp_vec2o[i * 3 + 1] = vec2i[i * 3 + 1] * cAng - vec1i[i * 3 + 1] * sAng;
        tmp_vec2o[i * 3 + 2] = vec2i[i * 3 + 2] * cAng - vec1i[i * 3 + 2] * sAng;
    }
    double magnitude1 = norm(vec1i, n);
    double magnitude2 = norm(vec2i, n);
    *vec1o = normalize(tmp_vec1o, n);
    *vec2o = normalize(tmp_vec2o, n);
    for (i = 0; i < n; ++i) {
        (*vec1o)[i * 3 + 0] *= magnitude1;
        (*vec1o)[i * 3 + 1] *= magnitude1;
        (*vec1o)[i * 3 + 2] *= magnitude1;
        (*vec2o)[i * 3 + 0] *= magnitude2;
        (*vec2o)[i * 3 + 1] *= magnitude2;
        (*vec2o)[i * 3 + 2] *= magnitude2;
    }
    free(tmp_vec1o);
    free(tmp_vec2o);
}


double *projected_force(double *force0, double *eigenmode,
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


void cut_sphere(Config *config, Input *input, double *center, int *update_list)
{
    int i;
    double del[3];
    int cut_num = 0;
    int update_num = 0;
    int *cut_list = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = 0; i < config->tot_num; ++i) {
        del[0] = config->pos[i * 3 + 0] - center[0];
        del[1] = config->pos[i * 3 + 1] - center[1];
        del[2] = config->pos[i * 3 + 2] - center[2];
        get_minimum_image(del, config->boxlo, config->boxhi,
                          config->xy, config->yz, config->xz);
        double dist = sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);
        if (dist > input->cutoff * 2) {
            cut_list[cut_num] = i;
            cut_num++;
        } else {
            update_list[update_num] = i;
            update_num++;
        }
    }
    /* sort */
    while (1) {
        int done = 1;
        for (i = 1; i < cut_num; ++i) {
            if (cut_list[i - 1] < cut_list[i]) {
                int tmp_cut = cut_list[i - 1];
                cut_list[i - 1] = cut_list[i];
                cut_list[i] = tmp_cut;
                done = 0;
            }
        }
        if (done) {
            break;
        }
    }
    for (i = 0; i < cut_num; ++i) {
        extract_atom(config, cut_list[i]);
    }
    free(cut_list);
}


double *displace(Input *input, int tot_num,
                 int disp_num, int *disp_list, int count)
{
    int i, rank, size;
    double *disp = (double *)malloc(sizeof(double) * disp_num * 3);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        double *tmp_disp = (double *)calloc(tot_num * 3, sizeof(double));
        char line[128];
        FILE *fp;
        if (input->init_mode > 0) {
            fp = fopen("./MODECAR", "r");
            for (i = 0; i < tot_num; ++i) {
                fgets(line, 128, fp);
                tmp_disp[i * 3 + 0] = atof(strtok(line, " \n"));
                tmp_disp[i * 3 + 1] = atof(strtok(NULL, " \n"));
                tmp_disp[i * 3 + 2] = atof(strtok(NULL, " \n"));
            }
        } else {
            for (i = 0; i < disp_num; ++i) {
                tmp_disp[disp_list[i] * 3 + 0] = normal_random(0, input->stddev); 
                tmp_disp[disp_list[i] * 3 + 1] = normal_random(0, input->stddev); 
                tmp_disp[disp_list[i] * 3 + 2] = normal_random(0, input->stddev); 
            }
            char filename[128];
            sprintf(filename, "%s/Eigenmode_%d.dat", input->output_dir, count);
            fp = fopen(filename, "w");     
            for (i = 0; i < tot_num; ++i) {
                sprintf(line, "%e %e %e\n",
                        tmp_disp[i * 3 + 0],
                        tmp_disp[i * 3 + 1],
                        tmp_disp[i * 3 + 2]);
                fputs(line, fp);
            }
        }
        fclose(fp);
        for (i = 0; i < disp_num; ++i) {
            disp[i * 3 + 0] = tmp_disp[disp_list[i] * 3 + 0];
            disp[i * 3 + 1] = tmp_disp[disp_list[i] * 3 + 1];
            disp[i * 3 + 2] = tmp_disp[disp_list[i] * 3 + 2];
        }
        free(tmp_disp);
    }
    MPI_Bcast(disp, disp_num * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return disp;
}


double *get_rot_force(Input *input, double *force1, double *force2,
                      double *eigenmode, int disp_num)
{
    int i;
    double *dforce = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        dforce[i * 3 + 0] = force1[i * 3 + 0] - force2[i * 3 + 0];
        dforce[i * 3 + 1] = force1[i * 3 + 1] - force2[i * 3 + 1];
        dforce[i * 3 + 2] = force1[i * 3 + 2] - force2[i * 3 + 2];
    }
    double *rot_force = perpendicular_vector(dforce, eigenmode, disp_num);
    for (i = 0; i < disp_num; ++i) {
        rot_force[i * 3 + 0] /= 2 * input->dimer_dist;
        rot_force[i * 3 + 1] /= 2 * input->dimer_dist;
        rot_force[i * 3 + 2] /= 2 * input->dimer_dist;
    }
    free(dforce);
    return rot_force;
}


void rotate(Config *config0, Input *input, int disp_num, int *disp_list,
            double *eigenmode, int count, int dimer_step)
{
    int i, j, rank, size;
    double magnitude;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double energy0, energy1, energy2;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot(config0, input, &energy0, force0, disp_num, disp_list); 
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
        oneshot(config1, input, &energy1, force1, disp_num, disp_list);
        free_config(config1);
        for (j = 0; j < disp_num; ++j) {
            force2[j * 3 + 0] = 2 * force0[j * 3 + 0] - force1[j * 3 + 0];
            force2[j * 3 + 1] = 2 * force0[j * 3 + 1] - force1[j * 3 + 1];
            force2[j * 3 + 2] = 2 * force0[j * 3 + 2] - force1[j * 3 + 2];
        }
        double *f_rot_A = get_rot_force(input, force1, force2,
                                        eigenmode, disp_num);
        /* no rotation */
        if (norm(f_rot_A, disp_num) < input->f_rot_min) {
            if (rank == 0) {
                char line[128], filename[128];
                sprintf(filename, "%s/Dimer_%d.log", input->output_dir, count);
                FILE *fp = fopen(filename, "a");
                sprintf(line, " %8d   %8d   %16f   ---------   ---------   %9f\n",
                        dimer_step, i, energy0, norm(f_rot_A, disp_num));
                fputs(line, fp);
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
        oneshot(trial_config1, input, &energy1, force1, disp_num, disp_list);
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
        double cmin = 0.5 * a0 + a1 * cos(2 * rotangle) + b1 * sin(2 * rotangle);
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
        if (rank == 0) {
            char line[128], filename[128];
            sprintf(filename, "%s/Dimer_%d.log", input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            sprintf(line, " %8d   %8d   %16f   %9f   %9f   %9f\n",
                    dimer_step, i + 1, energy0, cmin,
                    rotangle * 180 / 3.1415926535897932384624,
                    norm(f_rot_A, disp_num));
            fputs(line, fp);
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


void translate(Config *config0, Input *input, int disp_num, int *disp_list,
               double *eigenmode, double *direction_old, double *cg_direction,
               int dimer_step)
{
    int i;
    double magnitude;
    double energy0, energy1, energy2;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot(config0, input, &energy0, force0, disp_num, disp_list); 
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
    oneshot(config1, input, &energy1, force1, disp_num, disp_list);
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
    double *f0p = projected_force(force0, eigenmode, curvature, disp_num);
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
        oneshot(trial_config0, input, &trial_energy0, trial_force0,
                disp_num, disp_list); 
        double *f0tp = projected_force(trial_force0, eigenmode,
                                       curvature, disp_num);
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
double dimer(Config *config0, Input *input, int count, int ii)
{
    int i, j, rank, size;
    int tot_num = config0->tot_num;
    double del[3];
    double center[3] = {config0->pos[ii * 3 + 0],
                        config0->pos[ii * 3 + 1],
                        config0->pos[ii * 3 + 2]};

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* First, cut far atoms */
    Config *template = (Config *)malloc(sizeof(Config));
    copy_config(template, config0);
    int *update_list = (int *)malloc(sizeof(int) * config0->tot_num);
    cut_sphere(config0, input, center, update_list);
    atom_relax(config0, input); 

    /* Save initial config */
    Config *origin = (Config *)malloc(sizeof(Config));
    copy_config(origin, config0);

    /* Second, set dimer space */
    int disp_num = 0;
    int *disp_list = (int *)malloc(sizeof(int) * config0->tot_num);
    for (i = 0; i < config0->tot_num; ++i) {
        del[0] = config0->pos[i * 3 + 0] - center[0];
        del[1] = config0->pos[i * 3 + 1] - center[1];
        del[2] = config0->pos[i * 3 + 2] - center[2];
        get_minimum_image(del, config0->boxlo, config0->boxhi,
                          config0->xy, config0->yz, config0->xz);
        double dist = sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);
        if (dist < input->disp_cutoff) {
            disp_list[disp_num] = i;
            disp_num++;
        }
    }
    
    /* Third, initial energy */
    double initial_energy;
    double *initial_force = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot(config0, input, &initial_energy, initial_force, disp_num, disp_list);     
    free(initial_force);

    /* Fourth, displacement & eigenmode */
    double *disp = displace(input, tot_num, disp_num, disp_list, count);
    for (i = 0; i < disp_num; ++i) {
        config0->pos[disp_list[i] * 3 + 0] += disp[i * 3 + 0]; 
        config0->pos[disp_list[i] * 3 + 1] += disp[i * 3 + 1]; 
        config0->pos[disp_list[i] * 3 + 2] += disp[i * 3 + 2]; 
    }
    double *eigenmode = normalize(disp, disp_num);

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot(config0, input, &energy0, force0, disp_num, disp_list);     
    double *direction_old = (double *)malloc(sizeof(double) * disp_num * 3);
    double *cg_direction = (double *)malloc(sizeof(double) * disp_num * 3);

    /* Fifth, run */
    /* trajectory */
    if (rank == 0) {
        char filename[128];
        sprintf(filename, "%s/Dimer_%d.XDATCAR", input->output_dir, count);
        write_config(config0, filename);
    }
    double fmax;
    int dimer_step = 1;
    do {
        rotate(config0, input, disp_num, disp_list, eigenmode,
               count, dimer_step);
        translate(config0, input, disp_num, disp_list, eigenmode,
                  direction_old, cg_direction, dimer_step);
        oneshot(config0, input, &energy0, force0, disp_num, disp_list);     
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
        if (rank == 0) {
            char filename[128];
            sprintf(filename, "%s/Dimer_%d.XDATCAR", input->output_dir, count);
            write_config(config0, filename);
        }
        dimer_step++;
    } while (fmax > input->f_tol);
    oneshot(config0, input, &energy0, force0, disp_num, disp_list);     
    double saddle_energy = energy0;

    free(disp);
    free(direction_old);
    free(cg_direction);
    free(force0);

    /* saddle update */
    for (i = 0; i < config0->tot_num; ++i) {
        template->pos[update_list[i] * 3 + 0] = config0->pos[i * 3 + 0];
        template->pos[update_list[i] * 3 + 1] = config0->pos[i * 3 + 1];
        template->pos[update_list[i] * 3 + 2] = config0->pos[i * 3 + 2];
    }
    if (rank == 0) {
        char filename[128];
        sprintf(filename, "%s/Saddle_%d.POSCAR", input->output_dir, count);
        write_config(template, filename);
    }

    /* split */
    Config *config1 = (Config *)malloc(sizeof(Config));
    copy_config(config1, config0);
    Config *config2 = (Config *)malloc(sizeof(Config));
    copy_config(config2, config0);
    for (i = 1; i < 11; ++i) {
        for (j = 0; j < disp_num; ++j) {
            config1->pos[disp_list[j] * 3 + 0] += 0.2 * i * eigenmode[j * 3 + 0];
            config1->pos[disp_list[j] * 3 + 1] += 0.2 * i * eigenmode[j * 3 + 1];
            config1->pos[disp_list[j] * 3 + 2] += 0.2 * i * eigenmode[j * 3 + 2];
            config2->pos[disp_list[j] * 3 + 0] -= 0.2 * i * eigenmode[j * 3 + 0];
            config2->pos[disp_list[j] * 3 + 1] -= 0.2 * i * eigenmode[j * 3 + 1];
            config2->pos[disp_list[j] * 3 + 2] -= 0.2 * i * eigenmode[j * 3 + 2];
        }
        atom_relax(config1, input); 
        atom_relax(config2, input); 
        if (diff_config(config1, config2, input->max_step) > 0) {
            if (rank == 0) {
                char line[128], filename[128];
                sprintf(filename, "%s/Dimer_%d.log", input->output_dir, count);
                FILE *fp = fopen(filename, "a");
                fputs("---------------------------------------------------------------------------\n", fp);
                sprintf(line, " Split success: %f\n", 0.2 * i);
                fputs(line, fp);
                fclose(fp);
            }
            break;
        };
    }
    free(disp_list);
    free(eigenmode);

    /* log */
    if (rank == 0) {
        int diff1 = diff_config(origin, config1, input->max_step);
        int diff2 = diff_config(origin, config2, input->max_step);
        char line[128], filename[128];
        sprintf(filename, "%s/Dimer_%d.log", input->output_dir, count);
        FILE *fp = fopen(filename, "a");
        if (diff1 * diff2 > 0) {
            fputs(" Saddle state: disconnected\n", fp);
            sprintf(filename, "%s/Initial_%d.POSCAR", input->output_dir, count);
            for (i = 0; i < config0->tot_num; ++i) {
                template->pos[update_list[i] * 3 + 0] = config1->pos[i * 3 + 0];
                template->pos[update_list[i] * 3 + 1] = config1->pos[i * 3 + 1];
                template->pos[update_list[i] * 3 + 2] = config1->pos[i * 3 + 2];
            }
            write_config(template, filename);
            sprintf(filename, "%s/Final_%d.POSCAR", input->output_dir, count);
            for (i = 0; i < config0->tot_num; ++i) {
                template->pos[update_list[i] * 3 + 0] = config2->pos[i * 3 + 0];
                template->pos[update_list[i] * 3 + 1] = config2->pos[i * 3 + 1];
                template->pos[update_list[i] * 3 + 2] = config2->pos[i * 3 + 2];
            }
            write_config(template, filename);
        } else {
            fputs(" Saddle state: connected\n", fp);
            if (diff1 == 0) {
                sprintf(filename, "%s/Initial_%d.POSCAR", input->output_dir, count);
                for (i = 0; i < config0->tot_num; ++i) {
                    template->pos[update_list[i] * 3 + 0] = config1->pos[i * 3 + 0];
                    template->pos[update_list[i] * 3 + 1] = config1->pos[i * 3 + 1];
                    template->pos[update_list[i] * 3 + 2] = config1->pos[i * 3 + 2];
                }
                write_config(template, filename);
                sprintf(filename, "%s/Final_%d.POSCAR", input->output_dir, count);
                for (i = 0; i < config0->tot_num; ++i) {
                    template->pos[update_list[i] * 3 + 0] = config2->pos[i * 3 + 0];
                    template->pos[update_list[i] * 3 + 1] = config2->pos[i * 3 + 1];
                    template->pos[update_list[i] * 3 + 2] = config2->pos[i * 3 + 2];
                }
                write_config(template, filename);
            } else {
                sprintf(filename, "%s/Initial_%d.POSCAR", input->output_dir, count);
                for (i = 0; i < config0->tot_num; ++i) {
                    template->pos[update_list[i] * 3 + 0] = config2->pos[i * 3 + 0];
                    template->pos[update_list[i] * 3 + 1] = config2->pos[i * 3 + 1];
                    template->pos[update_list[i] * 3 + 2] = config2->pos[i * 3 + 2];
                }
                write_config(template, filename);
                sprintf(filename, "%s/Final_%d.POSCAR", input->output_dir, count);
                for (i = 0; i < config0->tot_num; ++i) {
                    template->pos[update_list[i] * 3 + 0] = config1->pos[i * 3 + 0];
                    template->pos[update_list[i] * 3 + 1] = config1->pos[i * 3 + 1];
                    template->pos[update_list[i] * 3 + 2] = config1->pos[i * 3 + 2];
                }
                write_config(template, filename);
            }
        }
        fclose(fp);
    }
    free_config(template);
    free_config(config1);
    free_config(config2);
    free_config(origin);
    free(update_list);

    return saddle_energy - initial_energy;
}
