#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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


double norm(double **vec, int n)
{
    int i;
    double output = 0.0;
    for (i = 0; i < n; ++i) {
        output += vec[i][0] * vec[i][0];
        output += vec[i][1] * vec[i][1];
        output += vec[i][2] * vec[i][2];
    }
    return sqrt(output);
}


double **normalize(double **vec, int n)
{
    int i;
    double **output = (double **)malloc(sizeof(double *) * n); 
    double magnitude = norm(vec, n);
    for (i = 0; i < n; ++i) {
        output[i] = (double *)malloc(sizeof(double) * 3);
        output[i][0] = vec[i][0] / magnitude;
        output[i][1] = vec[i][1] / magnitude;
        output[i][2] = vec[i][2] / magnitude;
    }
    return output;
}


double dot(double **vec1, double **vec2, int n)
{
    int i;
    double output = 0.0;
    for (i = 0; i < n; ++i) {
        output += vec1[i][0] * vec2[i][0];
        output += vec1[i][1] * vec2[i][1];
        output += vec1[i][2] * vec2[i][2];
    }
    return output;
}


double **parallel_vector(double **vector, double **unit, int n)
{
    int i; 
    double magnitude = dot(vector, unit, n);
    double **output = (double **)malloc(sizeof(double *) * n);
    for (i = 0; i < n; ++i) {
        output[i] = (double *)malloc(sizeof(double) * 3);
        output[i][0] = magnitude * unit[i][0];
        output[i][1] = magnitude * unit[i][1];
        output[i][2] = magnitude * unit[i][2];
    } 
    return output;
}


double **perpendicular_vector(double **vector, double **unit, int n)
{
    int i;
    double **tmp_vector = parallel_vector(vector, unit, n);
    double **output = (double **)malloc(sizeof(double *) * n);
    for (i = 0; i < n; ++i) {
        output[i] = (double *)malloc(sizeof(double) * 3); 
        output[i][0] = vector[i][0] - tmp_vector[i][0];
        output[i][1] = vector[i][1] - tmp_vector[i][1];
        output[i][2] = vector[i][2] - tmp_vector[i][2];
    } 
    for (i = 0; i < n; ++i) {
        free(tmp_vector[i]);
    }
    free(tmp_vector);
    return output;
}


void rotate_vector(double **vec0i, double **vec0j,
                   double ***vec1i, double ***vec1j,
                   int n, double angle)
{
    int i;
    double cAng = cos(angle);
    double sAng = sin(angle);
    double **tmp_vec1i = (double **)malloc(sizeof(double *) * n);
    double **tmp_vec1j = (double **)malloc(sizeof(double *) * n);
    for (i = 0; i < n; ++i) {
        tmp_vec1i[i] = (double *)malloc(sizeof(double) * 3);
        tmp_vec1i[i][0] = vec0i[i][0] * cAng + vec0j[i][0] * sAng;
        tmp_vec1i[i][1] = vec0i[i][1] * cAng + vec0j[i][1] * sAng;
        tmp_vec1i[i][2] = vec0i[i][2] * cAng + vec0j[i][2] * sAng;
        tmp_vec1j[i] = (double *)malloc(sizeof(double) * 3);
        tmp_vec1j[i][0] = vec0j[i][0] * cAng - vec0i[i][0] * sAng;
        tmp_vec1j[i][1] = vec0j[i][1] * cAng - vec0i[i][1] * sAng;
        tmp_vec1j[i][2] = vec0j[i][2] * cAng - vec0i[i][2] * sAng;
    }
    double magnitude1 = norm(vec0i, n);
    double magnitude2 = norm(vec0j, n);
    *vec1i = normalize(tmp_vec1i, n);
    *vec1j = normalize(tmp_vec1j, n);
    for (i = 0; i < n; ++i) {
        (*vec1i)[i][0] *= magnitude1;
        (*vec1i)[i][1] *= magnitude1;
        (*vec1i)[i][2] *= magnitude1;
        free(tmp_vec1i[i]);
        (*vec1j)[i][0] *= magnitude2;
        (*vec1j)[i][1] *= magnitude2;
        (*vec1j)[i][2] *= magnitude2;
        free(tmp_vec1j[i]);
    }
    free(tmp_vec1i);
    free(tmp_vec1j);
}


void cut_sphere(Config *config, Input *input, double *center)
{
    int i;
    double del[3];

    int cut_num = 0;
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


double **displace(Input *input, int n)
{
    int i, j;
    double **disp = (double **)malloc(sizeof(double *) * n);
    if (input->init_mode > 0) {
        // TODO: MODECAR
    } else {
        // TODO: orthogonalization
        for (i = 0; i < n; ++i) {
            disp[i] = (double *)malloc(sizeof(double) * 3);
            disp[i][0] = normal_random(0, input->stddev); 
            disp[i][1] = normal_random(0, input->stddev); 
            disp[i][2] = normal_random(0, input->stddev); 
        }
    }
    return disp;
}


double **get_rot_force(Input *input, double **force1, double **force2,
                       double **eigenmode, int disp_num)
{
    int i;
    double **dforce = (double **)malloc(sizeof(double *) * disp_num);
    for (i = 0; i < disp_num; ++i) {
        dforce[i] = (double *)malloc(sizeof(double) * 3);
        dforce[i][0] = force1[i][0] - force2[i][0];
        dforce[i][1] = force1[i][1] - force2[i][1];
        dforce[i][2] = force1[i][2] - force2[i][2];
    }
    double **rot_force = perpendicular_vector(dforce, eigenmode, disp_num);
    for (i = 0; i < disp_num; ++i) {
        rot_force[i][0] /= 2 * input->dimer_dist;
        rot_force[i][1] /= 2 * input->dimer_dist;
        rot_force[i][2] /= 2 * input->dimer_dist;
        free(dforce[i]);
    }
    free(dforce);
    return rot_force;
}


void rotate(Config *config0, Input *input, int disp_num, int *disp_list,
            double **eigenmode)
{
    int i;
    double magnitude;
    double energy0, energy1, energy2;
    double **force0 = (double **)malloc(sizeof(double *) * disp_num);
    double **force1 = (double **)malloc(sizeof(double *) * disp_num);
    double **force2 = (double **)malloc(sizeof(double *) * disp_num);
    for (i = 0; i < disp_num; ++i) {
        force0[i] = (double *)malloc(sizeof(double) * 3);
        force1[i] = (double *)malloc(sizeof(double) * 3);
        force2[i] = (double *)malloc(sizeof(double) * 3);
    }
    oneshot(config0, input, &energy0, &force0, disp_num, disp_list); 
    for (i = 0; i < input->max_num_rot; ++i) {
        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (i = 0; i < disp_num; ++i) {
            config1->pos[disp_list[i] * 3 + 0] += input->dimer_dist
                                                * eigenmode[i][0];
            config1->pos[disp_list[i] * 3 + 1] += input->dimer_dist 
                                                * eigenmode[i][1];
            config1->pos[disp_list[i] * 3 + 2] += input->dimer_dist 
                                                * eigenmode[i][2];
        }
        oneshot(config1, input, &energy1, &force1, disp_num, disp_list);
        free_config(config1);
        for (i = 0; i < disp_num; ++i) {
            force2[i][0] = 2 * force0[i][0] - force1[i][0];
            force2[i][1] = 2 * force0[i][1] - force1[i][1];
            force2[i][2] = 2 * force0[i][2] - force1[i][2];
        }
        double **f_rot_A = get_rot_force(input, force1, force2,
                                         eigenmode, disp_num);
        if (norm(f_rot_A, disp_num) < input->f_rot_min) {
            for (i = 0; i < disp_num; ++i) {
                free(f_rot_A[i]);
            }
            free(f_rot_A);
            break;
        }
        double **n_A = (double **)malloc(sizeof(double *) * disp_num);
        double **rot_unit_A = normalize(f_rot_A, disp_num);
        /* curvature */
        double **dforce = (double **)malloc(sizeof(double *) * disp_num);
        for (i = 0; i < disp_num; ++i) {
            n_A[i] = (double *)malloc(sizeof(double) * 3);
            n_A[i][0] = eigenmode[i][0];
            n_A[i][1] = eigenmode[i][1];
            n_A[i][2] = eigenmode[i][2];
            dforce[i] = (double *)malloc(sizeof(double) * 3);
            dforce[i][0] = force2[i][0] - force1[i][0];
            dforce[i][1] = force2[i][1] - force1[i][1];
            dforce[i][2] = force2[i][2] - force1[i][2];
        }
        magnitude = dot(dforce, eigenmode, disp_num);
        double c0 = magnitude / (2.0 * input->dimer_dist);
        magnitude = dot(dforce, rot_unit_A, disp_num);
        double c0d = magnitude / input->dimer_dist;
        /* trial rotation */
        double **n_B, **rot_unit_B;
        rotate_vector(n_A, rot_unit_A, &n_B, &rot_unit_B,
                      disp_num, input->trial_angle); 
        Config *trial_config1 = (Config *)malloc(sizeof(Config));
        copy_config(trial_config1, config0);
        for (i = 0; i < disp_num; ++i) {
            trial_config1->pos[disp_list[i] * 3 + 0] += n_B[i][0]
                                                      * input->dimer_dist;
            trial_config1->pos[disp_list[i] * 3 + 1] += n_B[i][1]
                                                      * input->dimer_dist;
            trial_config1->pos[disp_list[i] * 3 + 2] += n_B[i][2]
                                                      * input->dimer_dist;
        } 
        /* derivative of curvature */
        oneshot(trial_config1, input, &energy1, &force1, disp_num, disp_list);
        free_config(trial_config1);
        for (i = 0; i < disp_num; ++i) {
            force2[i][0] = 2 * force0[i][0] - force1[i][0];
            force2[i][1] = 2 * force0[i][1] - force1[i][1];
            force2[i][2] = 2 * force0[i][2] - force1[i][2];
            dforce[i][0] = force2[i][0] - force1[i][0];
            dforce[i][1] = force2[i][1] - force1[i][1];
            dforce[i][2] = force2[i][2] - force1[i][2];
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
        double **new_eigenmode, **tmp_force;
        rotate_vector(n_A, rot_unit_A, &new_eigenmode, &tmp_force,
                      disp_num, rotangle);
        for (i = 0; i < disp_num; ++i) {
            eigenmode[i][0] = new_eigenmode[i][0];
            eigenmode[i][1] = new_eigenmode[i][1];
            eigenmode[i][2] = new_eigenmode[i][2];
            free(n_A[i]);
            free(n_B[i]);
            free(dforce[i]);
            free(rot_unit_A[i]);
            free(rot_unit_B[i]);
            free(new_eigenmode[i]);
            free(tmp_force[i]);
        }
        free(n_A);
        free(n_B);
        free(dforce);
        free(rot_unit_A);
        free(rot_unit_B);
        free(new_eigenmode);
        free(tmp_force);
        if (norm(f_rot_A, disp_num) < input->f_rot_max) {
            for (i = 0; i < disp_num; ++i) {
                free(f_rot_A[i]);
            }
            free(f_rot_A);
            break;
        }
        for (i = 0; i < disp_num; ++i) {
            free(f_rot_A[i]);
        }
        free(f_rot_A);
    }
    for (i = 0; i < disp_num; ++i) {
        free(force0[i]);
        free(force1[i]);
        free(force2[i]);
    }
    free(force0);
    free(force1);
    free(force2);
}


void translate(Config *config0, Input *input, int disp_num, int *disp_list,
               double **eigenmode)
{
    int i;
    double magnitude;
    double energy0, energy1, energy2;
    double **force0 = (double **)malloc(sizeof(double *) * disp_num);
    double **force1 = (double **)malloc(sizeof(double *) * disp_num);
    double **force2 = (double **)malloc(sizeof(double *) * disp_num);
    for (i = 0; i < disp_num; ++i) {
        force0[i] = (double *)malloc(sizeof(double) * 3);
        force1[i] = (double *)malloc(sizeof(double) * 3);
        force2[i] = (double *)malloc(sizeof(double) * 3);
    }
    oneshot(config0, input, &energy0, &force0, disp_num, disp_list); 
    Config *config1 = (Config *)malloc(sizeof(Config));
    copy_config(config1, config0);
    for (i = 0; i < disp_num; ++i) {
        config1->pos[disp_list[i] * 3 + 0] += input->dimer_dist
                                            * eigenmode[i][0];
        config1->pos[disp_list[i] * 3 + 1] += input->dimer_dist 
                                            * eigenmode[i][1];
        config1->pos[disp_list[i] * 3 + 2] += input->dimer_dist 
                                            * eigenmode[i][2];
    }
    oneshot(config1, input, &energy1, &force1, disp_num, disp_list);
    free_config(config1);
    for (i = 0; i < disp_num; ++i) {
        force2[i][0] = 2 * force0[i][0] - force1[i][0];
        force2[i][1] = 2 * force0[i][1] - force1[i][1];
        force2[i][2] = 2 * force0[i][2] - force1[i][2];
    }
    double **dforce = (double **)malloc(sizeof(double *) * disp_num);
    for (i = 0; i < disp_num; ++i) {
        dforce[i] = (double *)malloc(sizeof(double) * 3);
        dforce[i][0] = force2[i][0] - force1[i][0];
        dforce[i][1] = force2[i][1] - force1[i][1];
        dforce[i][2] = force2[i][2] - force1[i][2];
    }
    magnitude = dot(dforce, eigenmode, disp_num);
    double curvature = magnitude / (2.0 * input->dimer_dist);
    for (i = 0; i < disp_num; ++i) {
        free(dforce[i]);
        free(force0[i]);
        free(force1[i]);
        free(force2[i]);
    } 
    free(dforce);
    free(force0);
    free(force1);
    free(force2);
}


void dimer(Config *config0, Input *input, int ii)
{
    int i;
    double del[3];
    double center[3] = {config0->pos[ii * 3 + 0],
                        config0->pos[ii * 3 + 1],
                        config0->pos[ii * 3 + 2]};

    /* First, cut far atoms */
    cut_sphere(config0, input, center);

    /* Second, set dimer */
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

    double **disp = displace(input, disp_num);
    for (i = 0; i < disp_num; ++i) {
        config0->pos[disp_list[i] * 3 + 0] += disp[i][0]; 
        config0->pos[disp_list[i] * 3 + 1] += disp[i][1]; 
        config0->pos[disp_list[i] * 3 + 2] += disp[i][2]; 
    }
    double **eigenmode = normalize(disp, disp_num);

    double fmax = 0.0001;
    /* TODO: direction old */
    do {
        rotate(config0, input, disp_num, disp_list, eigenmode);
        translate(config0, input, disp_num, disp_list, eigenmode);
        /* test */
    } while (fmax > input->f_rot_min);

    for (i = 0; i < disp_num; ++i) {
        free(disp[i]);
        free(eigenmode[i]);
    }
    free(disp);
    free(eigenmode);
    free(disp_list);
}
