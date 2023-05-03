#include "art_nouveau.h"
#include "calculator.h"
#include "config.h"
#include "linalg.h"
#include "utils.h"
#include <float.h>
#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


static void lanczos(Config *config, Input *input, int active_num, int *active_list,
                    double *eigenvalue, double *eigenmode, int *lanczos_step,
                    double *force0, MPI_Comm comm)
{
    int i, j;
    int size = 32;
    double *alpha = (double *)malloc(sizeof(double) * size);
    double *beta = (double *)malloc(sizeof(double) * size);
    double *eigenvector;

    /* V is transpose of column matrix */
    double *V = (double *)malloc(sizeof(double) * size * active_num * 3);
    double *HL = (double *)malloc(sizeof(double) * active_num * 3);
    /* stap 1 */
    alpha[0] = 0.0;
    for (i = 0; i < active_num; ++i) {
        V[i * 3 + 0] = eigenmode[i * 3 + 0];
        V[i * 3 + 1] = eigenmode[i * 3 + 1];
        V[i * 3 + 2] = eigenmode[i * 3 + 2];
        config->pos[active_list[i] * 3 + 0] += input->finite_diff * V[i * 3 + 0];
        config->pos[active_list[i] * 3 + 1] += input->finite_diff * V[i * 3 + 1];
        config->pos[active_list[i] * 3 + 2] += input->finite_diff * V[i * 3 + 2];
    }
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * active_num * 3);
    double *full_force = (double *)malloc(sizeof(double) * config->tot_num * 3);
    oneshot(config, input, &energy1, full_force, comm);
    for (i = 0; i < active_num; ++i) {
        force1[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
        force1[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
        force1[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
    }
    for (i = 0; i < active_num; ++i) {
        config->pos[active_list[i] * 3 + 0] -= input->finite_diff * V[i * 3 + 0];
        config->pos[active_list[i] * 3 + 1] -= input->finite_diff * V[i * 3 + 1];
        config->pos[active_list[i] * 3 + 2] -= input->finite_diff * V[i * 3 + 2];
        HL[i * 3 + 0] = (force0[i * 3 + 0] - force1[i * 3 + 0]) / input->finite_diff;
        HL[i * 3 + 1] = (force0[i * 3 + 1] - force1[i * 3 + 1]) / input->finite_diff;
        HL[i * 3 + 2] = (force0[i * 3 + 2] - force1[i * 3 + 2]) / input->finite_diff;
        alpha[0] += V[i * 3 + 0] * HL[i * 3 + 0];
        alpha[0] += V[i * 3 + 1] * HL[i * 3 + 1];
        alpha[0] += V[i * 3 + 2] * HL[i * 3 + 2];
    }
    for (i = 0; i < active_num; ++i) {
        V[active_num * 3 + i * 3 + 0] = HL[i * 3 + 0] - alpha[0] * V[i * 3 + 0];
        V[active_num * 3 + i * 3 + 1] = HL[i * 3 + 1] - alpha[0] * V[i * 3 + 1];
        V[active_num * 3 + i * 3 + 2] = HL[i * 3 + 2] - alpha[0] * V[i * 3 + 2];
    }
    beta[0] = norm(&V[active_num * 3], active_num);
    for (i = 0; i < active_num; ++i) {
        V[active_num * 3 + i * 3 + 0] /= beta[0];
        V[active_num * 3 + i * 3 + 1] /= beta[0];
        V[active_num * 3 + i * 3 + 2] /= beta[0];
    }

    /* step 2 */
    double lambda_old = 0;
    double lambda_new;
    double criteria;
    int k = 1;
    while (1) {
        if (k + 1 >= size) {
            size = size << 1;
            alpha = (double *)realloc(alpha, sizeof(double) * size);
            beta = (double *)realloc(beta, sizeof(double) * size);
            V = (double *)realloc(V, sizeof(double) * size * active_num * 3);
        }
        double *L0 = &V[(k - 1) * active_num * 3];
        double *L1 = &V[k * active_num * 3];
        double *L2 = &V[(k + 1) * active_num * 3];
        for (i = 0; i < active_num; ++i) {
            config->pos[active_list[i] * 3 + 0] += input->finite_diff * L1[i * 3 + 0];
            config->pos[active_list[i] * 3 + 1] += input->finite_diff * L1[i * 3 + 1];
            config->pos[active_list[i] * 3 + 2] += input->finite_diff * L1[i * 3 + 2];
        } 
        oneshot(config, input, &energy1, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force1[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force1[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force1[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        alpha[k] = 0.0;
        for (i = 0; i < active_num; ++i) {
            config->pos[active_list[i] * 3 + 0] -= input->finite_diff * L1[i * 3 + 0];
            config->pos[active_list[i] * 3 + 1] -= input->finite_diff * L1[i * 3 + 1];
            config->pos[active_list[i] * 3 + 2] -= input->finite_diff * L1[i * 3 + 2];
            HL[i * 3 + 0] = (force0[i * 3 + 0] - force1[i * 3 + 0]) / input->finite_diff;
            HL[i * 3 + 1] = (force0[i * 3 + 1] - force1[i * 3 + 1]) / input->finite_diff;
            HL[i * 3 + 2] = (force0[i * 3 + 2] - force1[i * 3 + 2]) / input->finite_diff;
            alpha[k] += L1[i * 3 + 0] * HL[i * 3 + 0];
            alpha[k] += L1[i * 3 + 1] * HL[i * 3 + 1];
            alpha[k] += L1[i * 3 + 2] * HL[i * 3 + 2];
        }
        for (i = 0; i < active_num; ++i) {
            L2[i * 3 + 0] = HL[i * 3 + 0] - alpha[k] * L1[i * 3 + 0]
                          - beta[k - 1] * L0[i * 3 + 0];
            L2[i * 3 + 1] = HL[i * 3 + 1] - alpha[k] * L1[i * 3 + 1]
                          - beta[k - 1] * L0[i * 3 + 1];
            L2[i * 3 + 2] = HL[i * 3 + 2] - alpha[k] * L1[i * 3 + 2]
                          - beta[k - 1] * L0[i * 3 + 2];
        }
        /* gram-schmidt orthogonalization */
        for (i = 0; i < k + 1; ++i) {
            double dot_product = dot(L2, &V[i * active_num * 3], active_num);
            for (j = 0; j < active_num; ++j) {
                L2[j * 3 + 0] -= dot_product * V[i * active_num * 3 + j * 3 + 0];
                L2[j * 3 + 1] -= dot_product * V[i * active_num * 3 + j * 3 + 1];
                L2[j * 3 + 2] -= dot_product * V[i * active_num * 3 + j * 3 + 2];
            }
            double L2_norm = norm(L2, active_num);
            for (j = 0; j < active_num; ++j) {
                L2[j * 3 + 0] /= L2_norm;
                L2[j * 3 + 1] /= L2_norm;
                L2[j * 3 + 2] /= L2_norm;
            }
        }
        beta[k] = dot(L2, HL, active_num);
        /* caution index */
        k++;
        eigenvector = (double *)calloc(k * k, sizeof(double));
        eigenvector[0] = alpha[0];
        for (i = 1; i < k; ++i) {
            eigenvector[i * k + i - 1] = beta[i - 1];
            eigenvector[i * k + i] = alpha[i];
        }

        /* eigenvector consists of orthonormal columns */
        MKL_INT n = k;
        MKL_INT info, lwork;
        double wkopt;
        double *w = (double *)malloc(sizeof(double) * k);
        lwork = -1;
        dsyev("V", "U", &n, eigenvector, &n, w, &wkopt, &lwork, &info);
        lwork = (MKL_INT)wkopt;
        double *work = (double *)malloc(lwork * sizeof(double));
        dsyev("V", "U", &n, eigenvector, &n, w, work, &lwork, &info);
        lambda_new = w[0]; 
        criteria = fabs((lambda_new - lambda_old) / lambda_new);
        lambda_old = lambda_new;
        free(work);
        free(w);
        if (criteria < input->lambda_conv) {
            break;
        } else {
            free(eigenvector);
        }
    }

    *eigenvalue = lambda_new;
    *lanczos_step = k;
    double *y = (double *)calloc(active_num * 3, sizeof(double));
    for (i = 0; i < active_num; ++i) {
        for (j = 0; j < k; ++j) {
            y[i * 3 + 0] += V[j * active_num * 3 + i * 3 + 0] * eigenvector[j];
            y[i * 3 + 1] += V[j * active_num * 3 + i * 3 + 1] * eigenvector[j];
            y[i * 3 + 2] += V[j * active_num * 3 + i * 3 + 2] * eigenvector[j];
        }
        eigenmode[i * 3 + 0] = y[i * 3 + 0];
        eigenmode[i * 3 + 1] = y[i * 3 + 1];
        eigenmode[i * 3 + 2] = y[i * 3 + 2];
    }

    free(eigenvector);
    free(full_force);
    free(force1);
    free(alpha);
    free(beta);
    free(V);
    free(HL);
    free(y);
}


static void uphill_push(Config *config, Input *input, int active_num, int *active_list,
                        double eigenvalue, double *eigenmode, double **push_direction,
                        double *init_direction, double *force, int negative, MPI_Comm comm)
{
    int i;
    if (eigenvalue < 0) {
        double *parallel_force = parallel_vector(force, eigenmode, active_num);
        double f_norm = norm(parallel_force, active_num);
        double alpha = f_norm / fabs(eigenvalue);
        double dr = input->max_move < alpha ? input->max_move : alpha;
        double ratio = negative > input->mixing_step ?
                       1.0 : (double)negative / input->mixing_step;
        free(parallel_force);
        double *push_vector = (double *)malloc(sizeof(double) * active_num * 3);
        if (dot(force, eigenmode, active_num) > 0) {
            for (i = 0; i < active_num; ++i) {
                push_vector[i * 3 + 0] = init_direction[i * 3 + 0] * (1 - ratio)
                                       - eigenmode[i * 3 + 0] * ratio;
                push_vector[i * 3 + 1] = init_direction[i * 3 + 1] * (1 - ratio)
                                       - eigenmode[i * 3 + 1] * ratio;
                push_vector[i * 3 + 2] = init_direction[i * 3 + 2] * (1 - ratio)
                                       - eigenmode[i * 3 + 2] * ratio;
            }
        } else {
            for (i = 0; i < active_num; ++i) {
                push_vector[i * 3 + 0] = init_direction[i * 3 + 0] * (1 - ratio)
                                       + eigenmode[i * 3 + 0] * ratio;
                push_vector[i * 3 + 1] = init_direction[i * 3 + 1] * (1 - ratio)
                                       + eigenmode[i * 3 + 1] * ratio;
                push_vector[i * 3 + 2] = init_direction[i * 3 + 2] * (1 - ratio)
                                       + eigenmode[i * 3 + 2] * ratio;
            }
        }
        *push_direction = normalize(push_vector, active_num);
        free(push_vector);
        for (i = 0; i < active_num; ++i) {
            config->pos[active_list[i] * 3 + 0] += dr * (*push_direction)[i * 3 + 0];
            config->pos[active_list[i] * 3 + 1] += dr * (*push_direction)[i * 3 + 1];
            config->pos[active_list[i] * 3 + 2] += dr * (*push_direction)[i * 3 + 2];
        }
    } else {
        *push_direction = (double *)malloc(sizeof(double) * active_num * 3);
        for (i = 0; i < active_num; ++i) {
            (*push_direction)[i * 3 + 0] = init_direction[i * 3 + 0];
            (*push_direction)[i * 3 + 1] = init_direction[i * 3 + 1];
            (*push_direction)[i * 3 + 2] = init_direction[i * 3 + 2];
        }
        for (i = 0; i < active_num; ++i) {
            config->pos[active_list[i] * 3 + 0] += input->max_move
                                                 * (*push_direction)[i * 3 + 0];
            config->pos[active_list[i] * 3 + 1] += input->max_move
                                                 * (*push_direction)[i * 3 + 1];
            config->pos[active_list[i] * 3 + 2] += input->max_move
                                                 * (*push_direction)[i * 3 + 2];
        }
    }
}


static void perp_relax(Config *initial, Config *config0, Input *input,
                       int active_num, int *active_list, double eigenvalue,
                       double *push_direction, int count, int index,
                       int sps_step, int lanczos_step, int negative, MPI_Comm comm)
{
    int i, rank;
    char filename[128];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * active_num * 3);
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * active_num * 3);
    double *full_force = (double *)malloc(sizeof(double) * config0->tot_num * 3);

    double *disp_vector = (double *)malloc(sizeof(double) * active_num * 3);
    for (i = 0; i < active_num; ++i) {
        disp_vector[i * 3 + 0] = config0->pos[active_list[i] * 3 + 0]
                               - initial->pos[active_list[i] * 3 + 0];
        disp_vector[i * 3 + 1] = config0->pos[active_list[i] * 3 + 1]
                               - initial->pos[active_list[i] * 3 + 1];
        disp_vector[i * 3 + 2] = config0->pos[active_list[i] * 3 + 2]
                               - initial->pos[active_list[i] * 3 + 2];
    }
    double *perp_disp = perpendicular_vector(disp_vector, push_direction, active_num);
    free(disp_vector);
    double *disp_direction = normalize(perp_disp, active_num);
    free(perp_disp);

    /* cg optimization */
    double *direction_old = (double *)calloc(active_num * 3, sizeof(double));
    double *cg_direction = (double *)calloc(active_num * 3, sizeof(double));

    int relax_step;
    int end_step = eigenvalue < 0 ? 500 : input->max_num_rlx;
    for (relax_step = 0; relax_step < end_step; ++relax_step) {
        oneshot(config0, input, &energy0, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force0[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force0[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force0[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        /* trajectory */
        if (local_rank == 0) {
            sprintf(filename, "./%d.log", count);
            FILE *fp = fopen(filename, "a");
            if (sps_step > input->delay_step) {
                fprintf(fp, " %9d   %10d   %12d   %16f   %10f\n",
                        sps_step, relax_step, lanczos_step, energy0, eigenvalue);
            } else {
                fprintf(fp, " %9d   %10d   %12d   %16f   ----------\n",
                        sps_step, relax_step, lanczos_step, energy0);
            }
            fclose(fp);
            char header[128];
            sprintf(header, "%d_%d %d", count, index, sps_step);
            sprintf(filename, "./%d.XDATCAR", count);
            write_config(config0, filename, header, "a");
        }
    
        double *perp_force0;
        if (negative > input->hyper_step) {
            perp_force0 = perpendicular_vector(force0, push_direction, active_num);
        } else {
            double *push_parallel_force0 = parallel_vector(force0, push_direction, active_num);
            double *disp_parallel_force0 = parallel_vector(force0, disp_direction, active_num);
            perp_force0 = (double *)malloc(sizeof(double) * active_num * 3);
            for (i = 0; i < active_num; ++i) {
                perp_force0[i * 3 + 0] = force0[i * 3 + 0]
                                       - push_parallel_force0[i * 3 + 0]
                                       - disp_parallel_force0[i * 3 + 0];
                perp_force0[i * 3 + 1] = force0[i * 3 + 1]
                                       - push_parallel_force0[i * 3 + 1]
                                       - disp_parallel_force0[i * 3 + 1];
                perp_force0[i * 3 + 2] = force0[i * 3 + 2]
                                       - push_parallel_force0[i * 3 + 2]
                                       - disp_parallel_force0[i * 3 + 2];
            }
            free(push_parallel_force0);
            free(disp_parallel_force0);
        }
        if (eigenvalue < 0) {
            double *parallel_force0 = (double *)malloc(sizeof(double) * active_num * 3);
            for (i = 0; i < active_num; ++i) {
                parallel_force0[i * 3 + 0] = force0[i * 3 + 0] - perp_force0[i * 3 + 0];
                parallel_force0[i * 3 + 1] = force0[i * 3 + 1] - perp_force0[i * 3 + 1];
                parallel_force0[i * 3 + 2] = force0[i * 3 + 2] - perp_force0[i * 3 + 2];
            }
            if (norm(perp_force0, active_num) < norm(parallel_force0, active_num)) {
                free(parallel_force0);
                break;
            }
            free(parallel_force0);
        }

        /* cg direction */
        get_cg_direction(perp_force0, direction_old, cg_direction, active_num);
        double *direction = normalize(cg_direction, active_num);

        /* trial step */
        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (i = 0; i < active_num; ++i) {
            config1->pos[active_list[i] * 3 + 0] += direction[i * 3 + 0]
                                                  * input->trial_move;
            config1->pos[active_list[i] * 3 + 1] += direction[i * 3 + 1]
                                                  * input->trial_move;
            config1->pos[active_list[i] * 3 + 2] += direction[i * 3 + 2]
                                                  * input->trial_move;
        }
        oneshot(config1, input, &energy1, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force1[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force1[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force1[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }

        double *perp_force1;
        if (negative > input->hyper_step) {
            perp_force1 = perpendicular_vector(force1, push_direction, active_num);
        } else {
            double *push_parallel_force1 = parallel_vector(force1, push_direction, active_num);
            double *disp_parallel_force1 = parallel_vector(force1, disp_direction, active_num);
            perp_force1 = (double *)malloc(sizeof(double) * active_num * 3);
            for (i = 0; i < active_num; ++i) {
                perp_force1[i * 3 + 0] = force1[i * 3 + 0]
                                       - push_parallel_force1[i * 3 + 0]
                                       - disp_parallel_force1[i * 3 + 0];
                perp_force1[i * 3 + 1] = force1[i * 3 + 1]
                                       - push_parallel_force1[i * 3 + 1]
                                       - disp_parallel_force1[i * 3 + 1];
                perp_force1[i * 3 + 2] = force1[i * 3 + 2]
                                       - push_parallel_force1[i * 3 + 2]
                                       - disp_parallel_force1[i * 3 + 2];
            }
            free(push_parallel_force1);
            free(disp_parallel_force1);
        }

        double *tmp_force = (double *)malloc(sizeof(double) * active_num * 3);
        for (i = 0; i < active_num; ++i) {
            tmp_force[i * 3 + 0] = perp_force1[i * 3 + 0] + perp_force0[i * 3 + 0];
            tmp_force[i * 3 + 1] = perp_force1[i * 3 + 1] + perp_force0[i * 3 + 1];
            tmp_force[i * 3 + 2] = perp_force1[i * 3 + 2] + perp_force0[i * 3 + 2];
        }
        double F = dot(tmp_force, direction, active_num) / 2.0;
        for (i = 0; i < active_num; ++i) {
            tmp_force[i * 3 + 0] = perp_force1[i * 3 + 0] - perp_force0[i * 3 + 0];
            tmp_force[i * 3 + 1] = perp_force1[i * 3 + 1] - perp_force0[i * 3 + 1];
            tmp_force[i * 3 + 2] = perp_force1[i * 3 + 2] - perp_force0[i * 3 + 2];
        }
        double C = dot(tmp_force, direction, active_num) / input->trial_move;
        double coeff = -F / C + input->trial_move * 0.5;
        double *step = (double *)malloc(sizeof(double) * active_num * 3);
        for (i = 0; i < active_num; ++i) {
            step[i * 3 + 0] = coeff * direction[i * 3 + 0];
            step[i * 3 + 1] = coeff * direction[i * 3 + 1];
            step[i * 3 + 2] = coeff * direction[i * 3 + 2];
        }
        if (norm(step, active_num) > input->max_move) {
            for (i = 0; i < active_num; ++i) {
                step[i * 3 + 0] = direction[i * 3 + 0] * input->max_move;
                step[i * 3 + 1] = direction[i * 3 + 1] * input->max_move;
                step[i * 3 + 2] = direction[i * 3 + 2] * input->max_move;
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
            free(perp_force0);
            free(perp_force1);
            free_config(config1);
            free(direction);
            free(tmp_force);
            free(step);
            break;
        }

        for (i = 0; i < active_num; ++i) {
            config0->pos[active_list[i] * 3 + 0] += step[i * 3 + 0];
            config0->pos[active_list[i] * 3 + 1] += step[i * 3 + 1];
            config0->pos[active_list[i] * 3 + 2] += step[i * 3 + 2];
        }
        free(perp_force0);
        free(perp_force1);
        free_config(config1);
        free(direction);
        free(tmp_force);
        free(step);
    }
    free(force0);
    free(force1);
    free(full_force);
    free(direction_old);
    free(cg_direction);
    free(disp_direction);
}


int art_nouveau(Config *initial, Config *saddle, Config *final, Input *input,
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

    /* initial */
    int tmp_num;
    int *tmp_list;
    double center[3] = {config0->pos[index * 3 + 0],
                        config0->pos[index * 3 + 1],
                        config0->pos[index * 3 + 2]};
    get_sphere_list(config0, input, center, input->disp_cutoff,
                    &tmp_num, &tmp_list, comm);
    double *init_direction = (double *)calloc(active_num * 3, sizeof(double));
    for (i = 0; i < active_num; ++i) {
        for (j = 0; j < tmp_num; ++j) {
            if (active_list[i] == tmp_list[j]) {
                init_direction[i * 3 + 0] = eigenmode[i * 3 + 0];
                init_direction[i * 3 + 1] = eigenmode[i * 3 + 1];
                init_direction[i * 3 + 2] = eigenmode[i * 3 + 2];
            }
        }
    }
    double *tmp_init_direction = normalize(init_direction, active_num);
    for (i = 0; i < active_num; ++i) {
        init_direction[i * 3 + 0] = tmp_init_direction[i * 3 + 0];
        init_direction[i * 3 + 1] = tmp_init_direction[i * 3 + 1];
        init_direction[i * 3 + 2] = tmp_init_direction[i * 3 + 2];
    }
    free(tmp_init_direction);
    if (input->init_disp > 0) {
        double *tmp_init_disp = (double *)calloc(active_num * 3, sizeof(double));
        for (i = 0; i < tmp_num; ++i) {
            if (config0->fix[tmp_list[i]] > 0) {
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

    /* normalize */
    double *tmp_eigenmode = normalize(eigenmode, active_num);
    for (i = 0; i < active_num; ++i) {
        eigenmode[i * 3 + 0] = tmp_eigenmode[i * 3 + 0];
        eigenmode[i * 3 + 1] = tmp_eigenmode[i * 3 + 1];
        eigenmode[i * 3 + 2] = tmp_eigenmode[i * 3 + 2];
    }
    free(tmp_eigenmode);

    if (local_rank == 0) {
        sprintf(filename, "./%d.log", count);
        FILE *fp = fopen(filename, "w");
        fprintf(fp, " %d_%d\n", count, index);
        fputs("-----------------------------------------------------------------------\n", fp);
        fputs(" Push step   Relax step   Lanczos step   Potential energy   Eigenvalue\n", fp);
        fputs("-----------------------------------------------------------------------\n", fp);
        fclose(fp);
        char header[128];
        sprintf(header, "%d_%d %d", count, index, 0);
        sprintf(filename, "./%d.XDATCAR", count);
        write_config(config0, filename, header, "w");
    }

    int sps_step;
    double eigenvalue = 1.0;
    int negative = 0;
    int lanczos_step = 0;
    double energy0;
    double *force0 = (double *)calloc(active_num * 3, sizeof(double));
    double *full_force = (double *)malloc(sizeof(double) * config0->tot_num * 3);
    clock_t start = clock();
    for (sps_step = 1; sps_step <= 500; ++sps_step) {
        oneshot(config0, input, &energy0, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force0[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force0[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force0[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        /* lanczos */
        if (sps_step > input->delay_step) {
            lanczos(config0, input, active_num, active_list,
                    &eigenvalue, eigenmode, &lanczos_step, force0, comm);
        }
        if (eigenvalue < 0) {
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
                break;
            }
            /* mixing & hyper */
            negative++;
        } else {
            negative = 0;
        }
        /* uphill push */
        double *push_direction;
        uphill_push(config0, input, active_num, active_list, eigenvalue, eigenmode,
                    &push_direction, init_direction, force0, negative, comm);
        /* normal relax */
        perp_relax(initial, config0, input, active_num, active_list, eigenvalue,
                   push_direction, count, index, sps_step, lanczos_step, negative, comm);
        free(push_direction);
    }
    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    free(force0);
    free(init_direction);
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
    double dE;
    conv = split_config(initial, saddle, final, input, Ea, &dE,
                        eigenvalue, eigenmode, active_num, active_list,
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
