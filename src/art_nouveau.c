#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include "art_nouveau.h"
#ifdef LMP
#include "lmp_calculator.h"
#endif
#ifdef VASP
#include "vasp_calculator.h"
#endif
#include "alg_utils.h"
#include "config.h"
#include "sps_utils.h"


static int lanczos(Config *config, Input *input, int disp_num, int *disp_list,
                   double *eigenvalue, double *eigenmode, MPI_Comm comm)
{
    int i, j;
    int size = 32;
    double *alpha = (double *)malloc(sizeof(double) * size);
    double *beta = (double *)malloc(sizeof(double) * size);
    double *eigenvector;

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config, input, &energy0, force0, disp_num, disp_list, comm);

    /* V is transpose of column matrix */
    double *V = (double *)malloc(sizeof(double) * size * disp_num * 3);
    double *HL = (double *)malloc(sizeof(double) * disp_num * 3);
    /* stap 1 */
    alpha[0] = 0.0;
    for (i = 0; i < disp_num; ++i) {
        V[i * 3 + 0] = eigenmode[i * 3 + 0];
        V[i * 3 + 1] = eigenmode[i * 3 + 1];
        V[i * 3 + 2] = eigenmode[i * 3 + 2];
        config->pos[disp_list[i] * 3 + 0] += input->disp_dist * V[i * 3 + 0];
        config->pos[disp_list[i] * 3 + 1] += input->disp_dist * V[i * 3 + 1];
        config->pos[disp_list[i] * 3 + 2] += input->disp_dist * V[i * 3 + 2];
    }
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config, input, &energy1, force1, disp_num, disp_list, comm);
    for (i = 0; i < disp_num; ++i) {
        config->pos[disp_list[i] * 3 + 0] -= input->disp_dist * V[i * 3 + 0];
        config->pos[disp_list[i] * 3 + 1] -= input->disp_dist * V[i * 3 + 1];
        config->pos[disp_list[i] * 3 + 2] -= input->disp_dist * V[i * 3 + 2];
        HL[i * 3 + 0] = (force0[i * 3 + 0] - force1[i * 3 + 0]) / input->disp_dist;
        HL[i * 3 + 1] = (force0[i * 3 + 1] - force1[i * 3 + 1]) / input->disp_dist;
        HL[i * 3 + 2] = (force0[i * 3 + 2] - force1[i * 3 + 2]) / input->disp_dist;
        alpha[0] += V[i * 3 + 0] * HL[i * 3 + 0];
        alpha[0] += V[i * 3 + 1] * HL[i * 3 + 1];
        alpha[0] += V[i * 3 + 2] * HL[i * 3 + 2];
    }
    for (i = 0; i < disp_num; ++i) {
        V[disp_num * 3 + i * 3 + 0] = HL[i * 3 + 0] - alpha[0] * V[i * 3 + 0];
        V[disp_num * 3 + i * 3 + 1] = HL[i * 3 + 1] - alpha[0] * V[i * 3 + 1];
        V[disp_num * 3 + i * 3 + 2] = HL[i * 3 + 2] - alpha[0] * V[i * 3 + 2];
    }
    beta[0] = norm(&V[disp_num * 3], disp_num);
    for (i = 0; i < disp_num; ++i) {
        V[disp_num * 3 + i * 3 + 0] /= beta[0];
        V[disp_num * 3 + i * 3 + 1] /= beta[0];
        V[disp_num * 3 + i * 3 + 2] /= beta[0];
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
            V = (double *)realloc(V, sizeof(double) * size * disp_num * 3);
        }
        double *L0 = &V[(k - 1) * disp_num * 3];
        double *L1 = &V[k * disp_num * 3];
        double *L2 = &V[(k + 1) * disp_num * 3];
        for (i = 0; i < disp_num; ++i) {
            config->pos[disp_list[i] * 3 + 0] += input->disp_dist * L1[i * 3 + 0];
            config->pos[disp_list[i] * 3 + 1] += input->disp_dist * L1[i * 3 + 1];
            config->pos[disp_list[i] * 3 + 2] += input->disp_dist * L1[i * 3 + 2];
        } 
        oneshot_disp(config, input, &energy1, force1, disp_num, disp_list, comm);
        alpha[k] = 0.0;
        for (i = 0; i < disp_num; ++i) {
            config->pos[disp_list[i] * 3 + 0] -= input->disp_dist * L1[i * 3 + 0];
            config->pos[disp_list[i] * 3 + 1] -= input->disp_dist * L1[i * 3 + 1];
            config->pos[disp_list[i] * 3 + 2] -= input->disp_dist * L1[i * 3 + 2];
            HL[i * 3 + 0] = (force0[i * 3 + 0] - force1[i * 3 + 0]) / input->disp_dist;
            HL[i * 3 + 1] = (force0[i * 3 + 1] - force1[i * 3 + 1]) / input->disp_dist;
            HL[i * 3 + 2] = (force0[i * 3 + 2] - force1[i * 3 + 2]) / input->disp_dist;
            alpha[k] += L1[i * 3 + 0] * HL[i * 3 + 0];
            alpha[k] += L1[i * 3 + 1] * HL[i * 3 + 1];
            alpha[k] += L1[i * 3 + 2] * HL[i * 3 + 2];
        }
        for (i = 0; i < disp_num; ++i) {
            L2[i * 3 + 0] = HL[i * 3 + 0]
                          - alpha[k] * L1[i * 3 + 0]
                          - beta[k - 1] * L0[i * 3 + 0];
            L2[i * 3 + 1] = HL[i * 3 + 1]
                          - alpha[k] * L1[i * 3 + 1]
                          - beta[k - 1] * L0[i * 3 + 1];
            L2[i * 3 + 2] = HL[i * 3 + 2]
                          - alpha[k] * L1[i * 3 + 2]
                          - beta[k - 1] * L0[i * 3 + 2];
        }
        /* gram-schmidt orthogonalization */
        for (i = 0; i < k + 1; ++i) {
            double dot_product = dot(L2, &V[i * disp_num * 3], disp_num); 
            for (j = 0; j < disp_num; ++j) {
                L2[j * 3 + 0] -= dot_product * V[i * disp_num * 3 + j * 3 + 0];
                L2[j * 3 + 1] -= dot_product * V[i * disp_num * 3 + j * 3 + 1];
                L2[j * 3 + 2] -= dot_product * V[i * disp_num * 3 + j * 3 + 2];
            }
        }
        double L2_norm = norm(L2, disp_num);
        for (j = 0; j < disp_num; ++j) {
            L2[j * 3 + 0] /= L2_norm;
            L2[j * 3 + 1] /= L2_norm;
            L2[j * 3 + 2] /= L2_norm;
        }
        beta[k] = dot(L2, HL, disp_num);
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
        /* test */
        for (i = 0; i < k; ++i) {
            printf("eigenvalue i %d: %f\n", i, w[i]);
        }
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
    double *y = (double *)calloc(disp_num * 3, sizeof(double));
    for (i = 0; i < disp_num; ++i) {
        for (j = 0; j < k; ++j) {
            y[i * 3 + 0] += V[j * disp_num * 3 + i * 3 + 0] * eigenvector[j];
            y[i * 3 + 1] += V[j * disp_num * 3 + i * 3 + 1] * eigenvector[j];
            y[i * 3 + 2] += V[j * disp_num * 3 + i * 3 + 2] * eigenvector[j];
        }
        eigenmode[i * 3 + 0] = y[i * 3 + 0];
        eigenmode[i * 3 + 1] = y[i * 3 + 1];
        eigenmode[i * 3 + 2] = y[i * 3 + 2];
    }

    free(eigenvector);
    free(force0);
    free(force1);
    free(alpha);
    free(beta);
    free(V);
    free(HL);
    free(y);

    return k;
}


static double *uphill_push(Config *config, Input *input,
                           int disp_num, int *disp_list,
                           double eigenvalue, double *eigenmode,
                           double *init_direction, int negative, MPI_Comm comm)
{
    int i;
    double energy;
    double *force = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config, input, &energy, force, disp_num, disp_list, comm);
    double *push_vector = (double *)malloc(sizeof(double) * disp_num * 3);
    if (eigenvalue < input->lambda_crit) {
        double *parallel_force = parallel_vector(force, eigenmode, disp_num);
        double f_norm = norm(parallel_force, disp_num);
        double divisor = fabs(eigenvalue) > fabs(input->lambda_crit) ?
                         fabs(eigenvalue) : fabs(input->lambda_crit);
        double alpha = f_norm / divisor;
        double dr = input->max_step < alpha ? input->max_step : alpha;
        double ratio = negative > input->art_mixing ?
                       1.0 : (double)negative / input->art_mixing;
        free(parallel_force);
        if (dot(force, eigenmode, disp_num) > 0) {
            for (i = 0; i < disp_num; ++i) {
                push_vector[i * 3 + 0] = dr * (init_direction[i * 3 + 0] * (1 - ratio)
                                                  - eigenmode[i * 3 + 0] * ratio);
                push_vector[i * 3 + 1] = dr * (init_direction[i * 3 + 1] * (1 - ratio)
                                                  - eigenmode[i * 3 + 1] * ratio);
                push_vector[i * 3 + 2] = dr * (init_direction[i * 3 + 2] * (1 - ratio)
                                                  - eigenmode[i * 3 + 2] * ratio);
            }
        } else {
            for (i = 0; i < disp_num; ++i) {
                push_vector[i * 3 + 0] = dr * (init_direction[i * 3 + 0] * (1 - ratio)
                                                  + eigenmode[i * 3 + 0] * ratio);
                push_vector[i * 3 + 1] = dr * (init_direction[i * 3 + 1] * (1 - ratio)
                                                  + eigenmode[i * 3 + 1] * ratio);
                push_vector[i * 3 + 2] = dr * (init_direction[i * 3 + 2] * (1 - ratio)
                                                  + eigenmode[i * 3 + 2] * ratio);
            }
        }
    } else {
        for (i = 0; i < disp_num; ++i) {
            push_vector[i * 3 + 0] = input->max_step * init_direction[i * 3 + 0];
            push_vector[i * 3 + 1] = input->max_step * init_direction[i * 3 + 1];
            push_vector[i * 3 + 2] = input->max_step * init_direction[i * 3 + 2];
        }
    }
    for (i = 0; i < disp_num; ++i) {
        config->pos[disp_list[i] * 3 + 0] += push_vector[i * 3 + 0];
        config->pos[disp_list[i] * 3 + 1] += push_vector[i * 3 + 1];
        config->pos[disp_list[i] * 3 + 2] += push_vector[i * 3 + 2];
    }
    double *push_direction = normalize(push_vector, disp_num);
    free(push_vector);
    free(force);
    return push_direction;
}


static void perp_relax(Config *config0, Input *input,
                       int disp_num, int *disp_list,
                       double eigenvalue, double *push_direction,
                       int count, int art_step, int lanczos_step, MPI_Comm comm)
{
    int i, rank, size;
    char filename[128];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);

    /* cg optimization */
    double *direction_old = (double *)calloc(disp_num * 3, sizeof(double));
    double *cg_direction = (double *)calloc(disp_num * 3, sizeof(double));

    int relax_step;
    int end_step = eigenvalue < input->lambda_crit ? 1000 : input->max_num_rlx;
    for (relax_step = 0; relax_step < end_step; ++relax_step) {
        oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
        double *perp_force0 = perpendicular_vector(force0, push_direction, disp_num);

        if (local_rank == 0) {
            sprintf(filename, "%s/SPS_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            if (art_step > input->art_delay) {
                fprintf(fp, " %8d   %8d   %12d   %16f   %10f\n",
                        art_step, relax_step, lanczos_step, energy0, eigenvalue);
            } else {
                fprintf(fp, " %8d   %8d   %12d   %16f   ----------\n",
                        art_step, relax_step, lanczos_step, energy0);
            }
            fclose(fp);
            sprintf(filename, "%s/SPS_%d.XDATCAR",
                    input->output_dir, count);
            write_config(config0, filename, "a");
        }
    
        if (eigenvalue < input->lambda_crit) {
            double *parallel_force0 = parallel_vector(force0, push_direction, disp_num);
            if (norm(perp_force0, disp_num) < norm(parallel_force0, disp_num)) {
                free(perp_force0);
                free(parallel_force0);
                break;
            }
        }

        /* cg direction */
        get_cg_direction(perp_force0, direction_old, cg_direction, disp_num);
        double *direction = normalize(cg_direction, disp_num);

        /* trial step */
        Config *config1 = (Config *)malloc(sizeof(Config));
        copy_config(config1, config0);
        for (i = 0; i < disp_num; ++i) {
            config1->pos[disp_list[i] * 3 + 0] += direction[i * 3 + 0]
                                                * input->trial_step;
            config1->pos[disp_list[i] * 3 + 1] += direction[i * 3 + 1]
                                                * input->trial_step;
            config1->pos[disp_list[i] * 3 + 2] += direction[i * 3 + 2]
                                                * input->trial_step;
        }
        oneshot_disp(config1, input, &energy1, force1, disp_num, disp_list, comm); 
        double *perp_force1 = perpendicular_vector(force1, push_direction, disp_num);

        double *tmp_force = (double *)malloc(sizeof(double) * disp_num * 3);
        for (i = 0; i < disp_num; ++i) {
            tmp_force[i * 3 + 0] = perp_force1[i * 3 + 0] + perp_force0[i * 3 + 0];
            tmp_force[i * 3 + 1] = perp_force1[i * 3 + 1] + perp_force0[i * 3 + 1];
            tmp_force[i * 3 + 2] = perp_force1[i * 3 + 2] + perp_force0[i * 3 + 2];
        }
        double F = dot(tmp_force, direction, disp_num) / 2.0;
        for (i = 0; i < disp_num; ++i) {
            tmp_force[i * 3 + 0] = perp_force1[i * 3 + 0] - perp_force0[i * 3 + 0];
            tmp_force[i * 3 + 1] = perp_force1[i * 3 + 1] - perp_force0[i * 3 + 1];
            tmp_force[i * 3 + 2] = perp_force1[i * 3 + 2] - perp_force0[i * 3 + 2];
        }
        double C = dot(tmp_force, direction, disp_num) / input->trial_step;
        double coeff = -F / C + input->trial_step * 0.5;
        double *step = (double *)malloc(sizeof(double) * disp_num * 3);
        for (i = 0; i < disp_num; ++i) {
            step[i * 3 + 0] = coeff * direction[i * 3 + 0];
            step[i * 3 + 1] = coeff * direction[i * 3 + 1];
            step[i * 3 + 2] = coeff * direction[i * 3 + 2];
        }
        if (norm(step, disp_num) > input->max_step) {
            for (i = 0; i < disp_num; ++i) {
                step[i * 3 + 0] = direction[i * 3 + 0] * input->max_step;
                step[i * 3 + 1] = direction[i * 3 + 1] * input->max_step;
                step[i * 3 + 2] = direction[i * 3 + 2] * input->max_step;
            }
        }
        for (i = 0; i < disp_num; ++i) {
            config0->pos[disp_list[i] * 3 + 0] += step[i * 3 + 0];
            config0->pos[disp_list[i] * 3 + 1] += step[i * 3 + 1];
            config0->pos[disp_list[i] * 3 + 2] += step[i * 3 + 2];
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
    free(direction_old);
    free(cg_direction);
}


int art_nouveau(Config *initial, Config *final, Input *input,
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
    free(tmp_eigenmode);

    double *init_direction = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        init_direction[i * 3 + 0] = eigenmode[i * 3 + 0];
        init_direction[i * 3 + 1] = eigenmode[i * 3 + 1];
        init_direction[i * 3 + 2] = eigenmode[i * 3 + 2];
    }
    /* perturbate starting config */
    if (input->init_disp > 0) {
        for (i = 0; i < disp_num; ++i) {
            config0->pos[disp_list[i] * 3 + 0] += input->stddev * init_direction[i * 3 + 0];
            config0->pos[disp_list[i] * 3 + 1] += input->stddev * init_direction[i * 3 + 1];
            config0->pos[disp_list[i] * 3 + 2] += input->stddev * init_direction[i * 3 + 2];
        }
    }

    double fmax;
    int art_step;
    int lanczos_step = 0;
    int converge = 0;
    double eigenvalue = 1;
    if (local_rank == 0) {
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "w");
        fputs("--------------------------------------------------------------------\n", fp);
        fputs(" Art step   Opt step   Lanczos step   Potential energy   Eigenvalue\n", fp);
        fputs("--------------------------------------------------------------------\n", fp);
        fclose(fp);
        sprintf(filename, "%s/SPS_%d.XDATCAR",
                input->output_dir, count);
        write_config(config0, filename, "w");
    }

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    int negative = 0;
    for (art_step = 1; art_step <= 1000; ++art_step) {
        /* lanczos */
        if (art_step > input->art_delay) {
            lanczos_step = lanczos(config0, input, disp_num, disp_list,
                                   &eigenvalue, eigenmode, comm);
        }
        /* test */
        if ((local_rank == 0) && (input->write_mode)) {
            sprintf(filename, "%s/%d_%d.MODECAR",
                    input->output_dir, count, art_step);
            FILE *fp = fopen(filename, "w");
            for (i = 0; i < disp_num; ++i) {
                fprintf(fp, "%f %f %f\n",
                        eigenmode[i * 3 + 0],
                        eigenmode[i * 3 + 1],
                        eigenmode[i * 3 + 2]);
            }
            fclose(fp);
        }
        /* mixing */
        if (eigenvalue < 0) {
            negative++;
        } else {
            negative = 0;
        }
        /* uphill push */
        double *push_direction = uphill_push(config0, input, disp_num, disp_list,
                                             eigenvalue, eigenmode,
                                             init_direction, negative, comm);
        /* normal relax */
        perp_relax(config0, input, disp_num, disp_list,
                   eigenvalue, push_direction, count, art_step, lanczos_step, comm);
        free(push_direction);

        /* force criteria */
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
        if (fmax < input->f_tol) {
            converge = 1;
            break;
        }
    }
    free(init_direction);
    if (local_rank == 0) {
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "a");
        fputs("--------------------------------------------------------------------\n", fp);
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
        free(force0);
        free_config(config0);
        free(disp_list);
        free(update_list);
        free(extract_list);
        free(full_eigenmode);
        free(eigenmode);
        return 1;
    }

    /* relax initial structure and barrier energy */
    atom_relax(initial, input, &energy0, comm);
    oneshot_disp(initial, input, &energy0, force0, disp_num, disp_list, comm);
    double i_energy = energy0;
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
    double ts_energy = energy0;
    free(force0);
    *Ea = ts_energy - i_energy;

    /* saddle update */
    for (i = 0; i < config0->tot_num; ++i) {
        final->pos[update_list[i] * 3 + 0] = config0->pos[i * 3 + 0];
        final->pos[update_list[i] * 3 + 1] = config0->pos[i * 3 + 1];
        final->pos[update_list[i] * 3 + 2] = config0->pos[i * 3 + 2];
    }
    for (i = 0; i < disp_num; ++i) {
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
    int conv = split_configs(initial, final, config0, input,
                             eigenmode, count, index,
                             update_num, update_list,
                             disp_num, disp_list, comm);

    if ((local_rank == 0) && (conv == 0)) {
        char filename[128];
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
