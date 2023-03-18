#include "art_nouveau.h"
#include "calculator.h"
#include "config.h"
#include "linalg.h"
#include "utils.h"
#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


static int lanczos(Config *config, Input *input, int active_num, int *active_list,
                   double *eigenvalue, double *eigenmode, double *force0, MPI_Comm comm)
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

    return k;
}


static void uphill_push(Config *config, Input *input, int active_num, int *active_list,
                        double eigenvalue, double *eigenmode, double *push_vector,
                        double *init_direction, double *force, int negative, MPI_Comm comm)
{
    int i;
    if (eigenvalue < input->lambda_crit) {
        double *parallel_force = parallel_vector(force, eigenmode, active_num);
        double f_norm = norm(parallel_force, active_num);
        double divisor = fabs(eigenvalue) > fabs(input->lambda_crit) ?
                         fabs(eigenvalue) : fabs(input->lambda_crit);
        double alpha = f_norm / divisor;
        double dr = input->max_move < alpha ? input->max_move : alpha;
        double ratio = negative > input->art_mixing ?
                       1.0 : (double)negative / input->art_mixing;
        free(parallel_force);
        if (dot(force, eigenmode, active_num) > 0) {
            for (i = 0; i < active_num; ++i) {
                push_vector[i * 3 + 0] = dr * (init_direction[i * 3 + 0] * (1 - ratio)
                                                  - eigenmode[i * 3 + 0] * ratio);
                push_vector[i * 3 + 1] = dr * (init_direction[i * 3 + 1] * (1 - ratio)
                                                  - eigenmode[i * 3 + 1] * ratio);
                push_vector[i * 3 + 2] = dr * (init_direction[i * 3 + 2] * (1 - ratio)
                                                  - eigenmode[i * 3 + 2] * ratio);
            }
        } else {
            for (i = 0; i < active_num; ++i) {
                push_vector[i * 3 + 0] = dr * (init_direction[i * 3 + 0] * (1 - ratio)
                                                  + eigenmode[i * 3 + 0] * ratio);
                push_vector[i * 3 + 1] = dr * (init_direction[i * 3 + 1] * (1 - ratio)
                                                  + eigenmode[i * 3 + 1] * ratio);
                push_vector[i * 3 + 2] = dr * (init_direction[i * 3 + 2] * (1 - ratio)
                                                  + eigenmode[i * 3 + 2] * ratio);
            }
        }
    } else {
        for (i = 0; i < active_num; ++i) {
            push_vector[i * 3 + 0] = input->max_move * init_direction[i * 3 + 0];
            push_vector[i * 3 + 1] = input->max_move * init_direction[i * 3 + 1];
            push_vector[i * 3 + 2] = input->max_move * init_direction[i * 3 + 2];
        }
    }
    for (i = 0; i < active_num; ++i) {
        config->pos[active_list[i] * 3 + 0] += push_vector[i * 3 + 0];
        config->pos[active_list[i] * 3 + 1] += push_vector[i * 3 + 1];
        config->pos[active_list[i] * 3 + 2] += push_vector[i * 3 + 2];
    }
}


static void perp_relax(Config *config0, Input *input, int active_num, int *active_list,
                       double eigenvalue, double *push_direction, int count, int index,
                       int sps_step, int lanczos_step, MPI_Comm comm)
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

    /* cg optimization */
    double *direction_old = (double *)calloc(active_num * 3, sizeof(double));
    double *cg_direction = (double *)calloc(active_num * 3, sizeof(double));

    int relax_step;
    int end_step = eigenvalue < input->lambda_crit ? 500 : input->max_num_rlx;
    for (relax_step = 0; relax_step < end_step; ++relax_step) {
        oneshot(config0, input, &energy0, full_force, comm);
        for (i = 0; i < active_num; ++i) {
            force0[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
            force0[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
            force0[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
        }
        /* trajectory */
        if (local_rank == 0) {
            sprintf(filename, "%s/%d.log", input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            if (sps_step > input->art_delay) {
                fprintf(fp, " %9d   %10d   %12d   %16f   %10f\n",
                        sps_step, relax_step, lanczos_step, energy0, eigenvalue);
            } else {
                fprintf(fp, " %9d   %10d   %12d   %16f   ----------\n",
                        sps_step, relax_step, lanczos_step, energy0);
            }
            fclose(fp);
            char header[128];
            sprintf(header, "%d_%d %d", count, index, sps_step);
            sprintf(filename, "%s/%d.XDATCAR", input->output_dir, count);
            write_config(config0, filename, header, "a");
        }
    
        double *perp_force0 = perpendicular_vector(force0, push_direction, active_num);
        if (eigenvalue < input->lambda_crit) {
            double *parallel_force0 = parallel_vector(force0, push_direction, active_num);
            if (norm(perp_force0, active_num) < norm(parallel_force0, active_num)) {
                free(perp_force0);
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
        double *perp_force1 = perpendicular_vector(force1, push_direction, active_num);

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

    /* starting configuration */
    Config *config0 = (Config *)malloc(sizeof(Config));
    copy_config(config0, initial);
    /* generate lists */
    int tmp_num;
    int *tmp_list;
    double center[3] = {config0->pos[index * 3 + 0],
                        config0->pos[index * 3 + 1],
                        config0->pos[index * 3 + 2]};
    get_sphere_list(config0, input, center, input->acti_cutoff,
                    &tmp_num, &tmp_list, comm);
    int active_num = 0;
    int *active_list = (int *)malloc(sizeof(int) * config0->tot_num);
    for (i = 0; i < tmp_num; ++i) {
        if (config0->fix[tmp_list[i]] == 0) {
            active_list[active_num] = tmp_list[i];
            active_num++;
        }
    }
    free(tmp_list);

    /* eigenmode */
    if (full_eigenmode == NULL) {
        full_eigenmode = get_eigenmode(input, config0->tot_num, comm);
    }
    double *eigenmode = (double *)calloc(config0->tot_num * 3, sizeof(double));
    for (i = 0; i < active_num; ++i) {
        eigenmode[i * 3 + 0] = full_eigenmode[active_list[i] * 3 + 0];
        eigenmode[i * 3 + 1] = full_eigenmode[active_list[i] * 3 + 1];
        eigenmode[i * 3 + 2] = full_eigenmode[active_list[i] * 3 + 2];
    }
    memset(full_eigenmode, 0, sizeof(double) * config0->tot_num * 3);

    /* initial perturbation */
    if (input->init_disp > 0) {
        get_sphere_list(config0, input, center, input->disp_cutoff,
                        &tmp_num, &tmp_list, comm);
        for (i = 0; i < active_num; ++i) {
            for (j = 0; j < tmp_num; ++j) {
                if (active_list[i] == tmp_list[j]) {
                    config0->pos[tmp_list[j] * 3 + 0] += eigenmode[i * 3 + 0];
                    config0->pos[tmp_list[j] * 3 + 1] += eigenmode[i * 3 + 1];
                    config0->pos[tmp_list[j] * 3 + 2] += eigenmode[i * 3 + 2];
                    break;
                }
            }
        }
        free(tmp_list);
    }
    /* normalize */
    double *tmp_eigenmode = normalize(eigenmode, active_num);
    for (i = 0; i < active_num; ++i) {
        eigenmode[i * 3 + 0] = tmp_eigenmode[i * 3 + 0];
        eigenmode[i * 3 + 1] = tmp_eigenmode[i * 3 + 1];
        eigenmode[i * 3 + 2] = tmp_eigenmode[i * 3 + 2];
    }
    free(tmp_eigenmode);
    double *init_direction = (double *)calloc(config0->tot_num * 3, sizeof(double));
    for (i = 0; i < active_num; ++i) {
        init_direction[i * 3 + 0] = eigenmode[i * 3 + 0];
        init_direction[i * 3 + 1] = eigenmode[i * 3 + 1];
        init_direction[i * 3 + 2] = eigenmode[i * 3 + 2];
    }

    if (local_rank == 0) {
        sprintf(filename, "%s/%d.log", input->output_dir, count);
        FILE *fp = fopen(filename, "w");
        fprintf(fp, " %d_%d\n", count, index);
        fputs("-----------------------------------------------------------------------\n", fp);
        fputs(" Push step   Relax step   Lanczos step   Potential energy   Eigenvalue\n", fp);
        fputs("-----------------------------------------------------------------------\n", fp);
        fclose(fp);
        char header[128];
        sprintf(header, "%d_%d %d", count, index, 0);
        sprintf(filename, "%s/%d.XDATCAR", input->output_dir, count);
        write_config(config0, filename, header, "w");
    }

    int sps_step;
    int max_index = index;
    double eigenvalue = 1.0;
    int negative = 0;
    int lanczos_step = 0;
    double energy0;
    double *force0 = (double *)calloc(config0->tot_num * 3, sizeof(double));
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
        if (sps_step > input->art_delay) {
            lanczos_step = lanczos(config0, input, active_num, active_list,
                                   &eigenvalue, eigenmode, force0, comm);
        }
        if (eigenvalue < input->lambda_crit) {
            /* force criteria */
            double fmax = 0.0;
            for (i = 0; i < active_num; ++i) {
                force0[i * 3 + 0] = full_force[active_list[i] * 3 + 0];
                force0[i * 3 + 1] = full_force[active_list[i] * 3 + 1];
                force0[i * 3 + 2] = full_force[active_list[i] * 3 + 2];
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
            } else {
                /* mixing */
                negative++;
            }
        } else {
            negative = 0;
        }
        /* uphill push */
        double *push_vector = (double *)malloc(sizeof(double) * active_num * 3);
        uphill_push(config0, input, active_num, active_list, eigenvalue,
                    eigenmode, push_vector, init_direction, force0, negative, comm);
        double *push_direction = normalize(push_vector, active_num);
        free(push_vector);
        /* normal relax */
        perp_relax(config0, input, active_num, active_list, eigenvalue,
                   push_direction, count, index, sps_step, lanczos_step, comm);
        free(push_direction);

        /* change active volume */
        if ((sps_step > input->acti_nevery) &&
            ((sps_step - 1) % input->acti_nevery == 0)) {
            expand_active_volume(initial, config0, input,
                                 &active_num, active_list, &max_index, comm);
        }
    }
    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    free(force0);
    free(init_direction);
    if (local_rank == 0) {
        sprintf(filename, "%s/%d.log", input->output_dir, count);
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
    conv = split_config(initial, saddle, final, input, Ea, &dE, eigenmode,
                        count, index, active_num, active_list, comm);
    if (local_rank == 0) {
        char filename[128];
        /* modecar */
        sprintf(filename, "%s/%d.MODECAR", input->output_dir, count);
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
        sprintf(filename, "%s/%d.log", input->output_dir, count);
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
