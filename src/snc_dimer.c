#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef LMP
#include "lmp_calculator.h"
#endif
#ifdef VASP
#include "vasp_calculator.h"
#endif
#include "alg_utils.h"
#include "config.h"
#include "snc_dimer.h"
#include "sps_utils.h"


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


static double *get_hessian(Config *config, Input *input,
                           int disp_num, int *disp_list, MPI_Comm comm)
{
    int i, j, k, l;
    double *H = (double *)calloc(3 * disp_num * 3 * disp_num, sizeof(double));
    double energy;
    double *force_1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force_2 = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        for (j = 0; j < 3; ++j) {
            config->pos[disp_list[i] * 3 + j] += 0.001;
            oneshot_disp(config, input, &energy, force_1, disp_num, disp_list, comm);     
            config->pos[disp_list[i] * 3 + j] -= 2 * 0.001;
            oneshot_disp(config, input, &energy, force_2, disp_num, disp_list, comm);     
            config->pos[disp_list[i] * 3 + j] += 0.001;
            /* column-major lower triangle matrix */
            for (k = i * 3 + j; k < 3 * disp_num; ++k) {
                /* f = -dE/dx */
                double dforce = force_2[k] - force_1[k];
                H[i * 3 * disp_num * 3 + j * disp_num * 3 + k] = dforce / (2 * 0.001);
            }
        }
    }
    free(force_1);
    free(force_2);
    return H;
}


static double *get_eigenvalue(double *H, int disp_num)
{
    MKL_INT n = disp_num * 3;
    MKL_INT info, lwork;
    double wkopt;
    double *w = (double *)malloc(sizeof(double) * disp_num * 3);
    lwork = -1;
    /* query and allocate the optimal workspace */
    dsyev("V", "L", &n, H, &n, w, &wkopt, &lwork, &info);
    lwork = (MKL_INT)wkopt;
    double *work = (double *)malloc(lwork * sizeof(double));
    /* solve */
    /* w: eigenvalues in ascending order
       H: the columns of H contain the orthonormal eigenvectors */
    dsyev("V", "L", &n, H, &n, w, work, &lwork, &info);
    free(work);
    return w;
}


/* cartesian -> snc */
static void transform_disp(double *vector, double *H, double *w,
                           int disp_num, int *disp_list)
{
    int i;
    double *x = (double *)malloc(sizeof(double) * 3 * disp_num);
    for (i = 0; i < disp_num; ++i) {
        x[i * 3 + 0] = vector[disp_list[i] * 3 + 0];
        x[i * 3 + 1] = vector[disp_list[i] * 3 + 1];
        x[i * 3 + 2] = vector[disp_list[i] * 3 + 2];
    }
    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasTrans;
    MKL_INT n = 3 * disp_num;
    double alpha = 1.0;
    double beta = 0.0;
    double *y = (double *)malloc(sizeof(double) * 3 * disp_num);
    cblas_dgemv(layout, trans, n, n, alpha, H, n, x, 1, beta, y, 1);
    for (i = 0; i < disp_num; ++i) {
        vector[disp_list[i] * 3 + 0] = y[i * 3 + 0] * sqrt(w[i * 3 + 0]);
        vector[disp_list[i] * 3 + 1] = y[i * 3 + 1] * sqrt(w[i * 3 + 1]);
        vector[disp_list[i] * 3 + 2] = y[i * 3 + 2] * sqrt(w[i * 3 + 2]);
    }
    free(x);
    free(y);
}
 

/* snc -> cartesian */
static void inv_transform_disp(double *vector, double *H, double *w,
                               int disp_num, int *disp_list)
{
    int i;
    double *x = (double *)malloc(sizeof(double) * 3 * disp_num);
    for (i = 0; i < disp_num; ++i) {
        x[i * 3 + 0] = vector[disp_list[i] * 3 + 0] / sqrt(w[i * 3 + 0]);
        x[i * 3 + 1] = vector[disp_list[i] * 3 + 1] / sqrt(w[i * 3 + 1]);
        x[i * 3 + 2] = vector[disp_list[i] * 3 + 2] / sqrt(w[i * 3 + 2]);
    }
    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    MKL_INT n = 3 * disp_num;
    double alpha = 1.0;
    double beta = 0.0;
    double *y = (double *)malloc(sizeof(double) * 3 * disp_num);
    cblas_dgemv(layout, trans, n, n, alpha, H, n, x, 1, beta, y, 1);
    for (i = 0; i < disp_num; ++i) {
        vector[disp_list[i] * 3 + 0] = y[i * 3 + 0];
        vector[disp_list[i] * 3 + 1] = y[i * 3 + 1];
        vector[disp_list[i] * 3 + 2] = y[i * 3 + 2];
    }
    free(x);
    free(y);
}
 

/* snc -> cartesian */
static void inv_transform(double *x, double *H, double *w, int disp_num)
{
    int i;
    for (i = 0; i < disp_num; ++i) {
        x[i * 3 + 0] /= sqrt(w[i * 3 + 0]);
        x[i * 3 + 1] /= sqrt(w[i * 3 + 1]);
        x[i * 3 + 2] /= sqrt(w[i * 3 + 2]);
    }
    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    MKL_INT n = 3 * disp_num;
    double alpha = 1.0;
    double beta = 0.0;
    double *y = (double *)malloc(sizeof(double) * 3 * disp_num);
    cblas_dgemv(layout, trans, n, n, alpha, H, n, x, 1, beta, y, 1);
    for (i = 0; i < disp_num; ++i) {
        x[i * 3 + 0] = y[i * 3 + 0];
        x[i * 3 + 1] = y[i * 3 + 1];
        x[i * 3 + 2] = y[i * 3 + 2];
    }
    free(y);
}
 

/* cartesian -> snc */
static void transform(double *x, double *H, double *w, int disp_num)
{
    int i;
    CBLAS_LAYOUT layout = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasTrans;
    MKL_INT n = 3 * disp_num;
    double alpha = 1.0;
    double beta = 0.0;
    double *y = (double *)malloc(sizeof(double) * 3 * disp_num);
    cblas_dgemv(layout, trans, n, n, alpha, H, n, x, 1, beta, y, 1);
    for (i = 0; i < disp_num; ++i) {
        x[i * 3 + 0] = y[i * 3 + 0] / sqrt(w[i * 3 + 0]);
        x[i * 3 + 1] = y[i * 3 + 1] / sqrt(w[i * 3 + 1]);
        x[i * 3 + 2] = y[i * 3 + 2] / sqrt(w[i * 3 + 2]);
    }
    free(y);
}


static void rotate(Config *config0, Input *input, int disp_num, int *disp_list,
                   double *eigenmode, double *eigenvector, double *eigenvalue,
                   int count, int dimer_step, MPI_Comm comm)
{
    int i, j, rank, size;
    double magnitude, cmin;
    char filename[128];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy0, energy1;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    /* convert coordinate */
    inv_transform_disp(config0->pos, eigenvector, eigenvalue, disp_num, disp_list);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
    transform_disp(config0->pos, eigenvector, eigenvalue, disp_num, disp_list);
    transform(force0, eigenvector, eigenvalue, disp_num);
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
        /* convert coordinate */
        inv_transform_disp(config1->pos, eigenvector, eigenvalue, disp_num, disp_list);
        oneshot_disp(config1, input, &energy1, force1, disp_num, disp_list, comm);
        transform_disp(config1->pos, eigenvector, eigenvalue, disp_num, disp_list);
        transform(force1, eigenvector, eigenvalue, disp_num);
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
                      disp_num, input->trial_angle); 
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
        /* convert coordinate */
        inv_transform_disp(trial_config1->pos, eigenvector, eigenvalue, disp_num, disp_list);
        oneshot_disp(trial_config1, input, &energy1, force1,
                     disp_num, disp_list, comm);
        transform_disp(trial_config1->pos, eigenvector, eigenvalue, disp_num, disp_list);
        transform(force1, eigenvector, eigenvalue, disp_num);
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
    free(force0);
    free(force1);
    free(force2);
}


static void translate(Config *config0, Input *input, int disp_num, int *disp_list,
                      double *eigenmode, double *eigenvector, double *eigenvalue,
                      double *direction_old, double *cg_direction,
                      int dimer_step, MPI_Comm comm)
{
    int i;
    double magnitude;
    char filename[128];
    double energy0, energy1;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    double *force2 = (double *)malloc(sizeof(double) * disp_num * 3);
    /* convert coordinate */
    inv_transform_disp(config0->pos, eigenvector, eigenvalue, disp_num, disp_list);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
    transform_disp(config0->pos, eigenvector, eigenvalue, disp_num, disp_list);
    transform(force0, eigenvector, eigenvalue, disp_num);
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
    inv_transform_disp(config1->pos, eigenvector, eigenvalue, disp_num, disp_list);
    oneshot_disp(config1, input, &energy1, force1, disp_num, disp_list, comm);
    transform_disp(config1->pos, eigenvector, eigenvalue, disp_num, disp_list);
    transform(force1, eigenvector, eigenvalue, disp_num);
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
        /* convert coordinate */
        inv_transform_disp(trial_config0->pos, eigenvector, eigenvalue, disp_num, disp_list);
        oneshot_disp(trial_config0, input, &trial_energy0, trial_force0,
                     disp_num, disp_list, comm); 
        transform_disp(trial_config0->pos, eigenvector, eigenvalue, disp_num, disp_list);
        transform(trial_force0, eigenvector, eigenvalue, disp_num);
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


int snc_dimer(Config *initial, Config *final, Input *input,
              double *full_eigenmode, int count, int index, double *Ea,
              MPI_Comm comm)
{
    int i, j, rank, size;
    char filename[128];

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
    set_active_volume(initial, input, center, &update_num, &update_list,
                      &extract_num, &extract_list, comm);
    trim_atoms(initial, update_num, update_list);

    /* starting dimer */ 
    Config *config0 = (Config *)malloc(sizeof(Config));
    copy_config(config0, initial);
    int tmp_num;
    int disp_num;
    int *tmp_list;
    int *disp_list;
    set_active_volume(config0, input, center, &tmp_num, &tmp_list,
                      &disp_num, &disp_list, comm);
    free(tmp_list);

    /* eigenmode */
    if (full_eigenmode == NULL) {
        full_eigenmode = get_eigenmode(input, final->tot_num, comm); 
    }

    /* normalize */
    double *tmp_eigenmode = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        tmp_eigenmode[i * 3 + 0] = full_eigenmode[extract_list[i] * 3 + 0];
        tmp_eigenmode[i * 3 + 1] = full_eigenmode[extract_list[i] * 3 + 1];
        tmp_eigenmode[i * 3 + 2] = full_eigenmode[extract_list[i] * 3 + 2];
    }
    double *eigenmode = normalize(tmp_eigenmode, disp_num);

    /* perturbate starting config */
    if (input->init_disp > 0) {
        for (i = 0; i < disp_num; ++i) {
            config0->pos[disp_list[i] * 3 + 0] += input->stddev * eigenmode[i * 3 + 0];
            config0->pos[disp_list[i] * 3 + 1] += input->stddev * eigenmode[i * 3 + 1];
            config0->pos[disp_list[i] * 3 + 2] += input->stddev * eigenmode[i * 3 + 2];
        }
    }

    /* cg optimization */
    double *direction_old = (double *)calloc(disp_num * 3, sizeof(double));
    double *cg_direction = (double *)calloc(disp_num * 3, sizeof(double));

    /* lower triangle matrix in col-major */
    double *H = get_hessian(initial, input, disp_num, disp_list, comm); 
    double *w = get_eigenvalue(H, disp_num);

    /* run */
    double fmax;
    int converge = 0;
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

    /* convert coordinate */
    transform_disp(config0->pos, H, w, disp_num, disp_list);
    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    for (dimer_step = 1; dimer_step <= 1000; ++dimer_step) {
        rotate(config0, input, disp_num, disp_list,
               eigenmode, H, w, count, dimer_step, comm);
        /* test */
        if ((local_rank == 0) && (input->write_mode)) {
            inv_transform(eigenmode, H, w, disp_num);
            for (i = 0; i < extract_num; ++i) {
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
            transform(eigenmode, H, w, disp_num);
        }
        translate(config0, input, disp_num, disp_list, eigenmode, H, w,
                  direction_old, cg_direction, dimer_step, comm);
        /* convert coordinate */
        inv_transform_disp(config0->pos, H, w, disp_num, disp_list);
        oneshot_disp(config0, input, &energy0, force0,
                     disp_num, disp_list, comm);     
        /* trajectory */
        if (local_rank == 0) {
            char filename[128];
            sprintf(filename, "%s/SPS_%d.XDATCAR",
                    input->output_dir, count);
            write_config(config0, filename, "a");
        }
        transform_disp(config0->pos, H, w, disp_num, disp_list);
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
        if (fmax < input->f_tol) {
            converge = 1;
            break;
        }
    }
    free(tmp_eigenmode);
    free(direction_old);
    free(cg_direction);
    if (local_rank == 0) {
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "a");
        fputs("----------------------------------------------------------------------------\n", fp);
        fclose(fp);
    }
    if (converge == 0) {
        if (local_rank == 0) {
            sprintf(filename, "%s/SPS_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            fputs(" Saddle state: not converged\n", fp);
            fclose(fp);
        }
        free_config(config0);
        free(force0);
        free(eigenmode);
        free(update_list);
        free(extract_list);
        free(H);
        free(w);
        return 1;
    }
    /* relax initial structure and barrier energy */
    atom_relax(initial, input, &energy0, comm);
    oneshot_disp(initial, input, &energy0, force0, disp_num, disp_list, comm);
    double i_energy = energy0;
    /* convert coordinate */
    inv_transform_disp(config0->pos, H, w, disp_num, disp_list);
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
    double ts_energy = energy0;
    free(force0);
    *Ea = ts_energy - i_energy;

    /* saddle update */
    for (i = 0; i < update_num; ++i) {
        final->pos[update_list[i] * 3 + 0] = config0->pos[i * 3 + 0];
        final->pos[update_list[i] * 3 + 1] = config0->pos[i * 3 + 1];
        final->pos[update_list[i] * 3 + 2] = config0->pos[i * 3 + 2];
    }
    for (i = 0; i < extract_num; ++i) {
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
    inv_transform(eigenmode, H, w, disp_num);
    tmp_eigenmode = normalize(eigenmode, disp_num);
    for (i = 0; i < disp_num; ++i) {
        eigenmode[i * 3 + 0] = tmp_eigenmode[i * 3 + 0];
        eigenmode[i * 3 + 1] = tmp_eigenmode[i * 3 + 1];
        eigenmode[i * 3 + 2] = tmp_eigenmode[i * 3 + 2];
    }
    free(tmp_eigenmode);
    int conv = split_configs(initial, final, config0, input,
                             eigenmode, count, index,
                             update_num, update_list,
                             disp_num, disp_list, comm);

    if ((local_rank == 0) && (conv == 0)) {
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "a");
        fprintf(fp, " Barrier energy: %f eV\n", *Ea);
        fclose(fp);
    }

    free_config(config0);
    free(extract_list);
    free(update_list);
    free(disp_list);
    free(full_eigenmode);
    free(eigenmode);
    return conv;
}
