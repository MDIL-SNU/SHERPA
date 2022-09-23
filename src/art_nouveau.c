#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include "art_nouveau.h"
#include "calculator.h"
#include "config.h"
#include "utils.h"


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


void lanczos(Config *config, Input *input, int disp_num, int *disp_list,
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

    /* distortion */
    for (i = 0; i < disp_num; ++i) {
        config->pos[disp_list[i] * 3 + 0] += input->disp_dist * eigenmode[i * 3 + 0];
        config->pos[disp_list[i] * 3 + 1] += input->disp_dist * eigenmode[i * 3 + 1];
        config->pos[disp_list[i] * 3 + 2] += input->disp_dist * eigenmode[i * 3 + 2];
    } 
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);
    oneshot_disp(config, input, &energy1, force1, disp_num, disp_list, comm);
    for (i = 0; i < disp_num; ++i) {
        config->pos[disp_list[i] * 3 + 0] -= input->disp_dist * eigenmode[i * 3 + 0];
        config->pos[disp_list[i] * 3 + 1] -= input->disp_dist * eigenmode[i * 3 + 1];
        config->pos[disp_list[i] * 3 + 2] -= input->disp_dist * eigenmode[i * 3 + 2];
    } 

    /* V is transpose of column matrix */
    double *V = (double *)malloc(sizeof(double) * size * disp_num * 3);
    double *Av = (double *)malloc(sizeof(double) * disp_num * 3);
    double *w = (double *)malloc(sizeof(double) * disp_num * 3);
    /* stap 1 */
    alpha[0] = 0.0;
    for (i = 0; i < disp_num; ++i) {
        V[i * 3 + 0] = eigenmode[i * 3 + 0];
        V[i * 3 + 1] = eigenmode[i * 3 + 1];
        V[i * 3 + 2] = eigenmode[i * 3 + 2];
        Av[i * 3 + 0] = (force0[i * 3 + 0] - force1[i * 3 + 0])
                      / input->disp_dist;
        Av[i * 3 + 1] = (force0[i * 3 + 1] - force1[i * 3 + 1])
                      / input->disp_dist;
        Av[i * 3 + 2] = (force0[i * 3 + 2] - force1[i * 3 + 2])
                      / input->disp_dist;
        alpha[0] += V[i * 3 + 0] * Av[i * 3 + 0];
        alpha[0] += V[i * 3 + 1] * Av[i * 3 + 1];
        alpha[0] += V[i * 3 + 2] * Av[i * 3 + 2];
    }
    for (i = 0; i < disp_num; ++i) {
        w[i * 3 + 0] = Av[i * 3 + 0] - alpha[0] * V[i * 3 + 0];
        w[i * 3 + 1] = Av[i * 3 + 1] - alpha[0] * V[i * 3 + 1];
        w[i * 3 + 2] = Av[i * 3 + 2] - alpha[0] * V[i * 3 + 2];
    }

    /* step 2 */
    double lambda_old = 0;
    double lambda_new;
    double criteria;
    int k = 0;
    while (1) {
        k++;
        if (k >= size) {
            size = size << 1;
            alpha = (double *)realloc(alpha, sizeof(double) * size);
            beta = (double *)realloc(beta, sizeof(double) * size);
            V = (double *)realloc(V, sizeof(double) * size * disp_num * 3);
        }
        beta[k] = norm(w, disp_num); 
        for (i = 0; i < disp_num; ++i) {
            V[k * disp_num * 3 + i * 3 + 0] = w[i * 3 + 0] / beta[k]; 
            V[k * disp_num * 3 + i * 3 + 1] = w[i * 3 + 1] / beta[k]; 
            V[k * disp_num * 3 + i * 3 + 2] = w[i * 3 + 2] / beta[k]; 
            /* distortion */
            config->pos[disp_list[i] * 3 + 0] += input->disp_dist
                                               * V[k * disp_num * 3 + i * 3 + 0];
            config->pos[disp_list[i] * 3 + 1] += input->disp_dist
                                               * V[k * disp_num * 3 + i * 3 + 1];
            config->pos[disp_list[i] * 3 + 2] += input->disp_dist
                                               * V[k * disp_num * 3 + i * 3 + 2];
        } 
        oneshot_disp(config, input, &energy1, force1, disp_num, disp_list, comm);
        alpha[k] = 0.0;
        for (i = 0; i < disp_num; ++i) {
            config->pos[disp_list[i] * 3 + 0] -= input->disp_dist
                                               * V[k * disp_num * 3 + i * 3 + 0];
            config->pos[disp_list[i] * 3 + 1] -= input->disp_dist
                                               * V[k * disp_num * 3 + i * 3 + 1];
            config->pos[disp_list[i] * 3 + 2] -= input->disp_dist
                                               * V[k * disp_num * 3 + i * 3 + 2];
            Av[i * 3 + 0] = (force0[i * 3 + 0] - force1[i * 3 + 0])
                          / input->disp_dist;
            Av[i * 3 + 1] = (force0[i * 3 + 1] - force1[i * 3 + 1])
                          / input->disp_dist;
            Av[i * 3 + 2] = (force0[i * 3 + 2] - force1[i * 3 + 2])
                          / input->disp_dist;
            alpha[k] += V[k * disp_num * 3 + i * 3 + 0] * Av[i * 3 + 0];
            alpha[k] += V[k * disp_num * 3 + i * 3 + 1] * Av[i * 3 + 1];
            alpha[k] += V[k * disp_num * 3 + i * 3 + 2] * Av[i * 3 + 2];
        }
        for (i = 0; i < disp_num; ++i) {
            w[i * 3 + 0] = Av[i * 3 + 0]
                         - alpha[k] * V[k * disp_num * 3 + i * 3 + 0]
                         - beta[k] * V[(k - 1) * disp_num * 3 + i * 3 + 0];
            w[i * 3 + 1] = Av[i * 3 + 1]
                         - alpha[k] * V[k * disp_num * 3 + i * 3 + 1]
                         - beta[k] * V[(k - 1) * disp_num * 3 + i * 3 + 1];
            w[i * 3 + 2] = Av[i * 3 + 2]
                         - alpha[k] * V[k * disp_num * 3 + i * 3 + 2]
                         - beta[k] * V[(k - 1) * disp_num * 3 + i * 3 + 2];
        }
        double *tmp_alpha = (double *)malloc(sizeof(double) * k);
        double *tmp_beta = (double *)malloc(sizeof(double) * k);
        tmp_alpha[0] = alpha[0];
        for (i = 1; i < k; ++i) {
            tmp_alpha[i] = alpha[i];
            tmp_beta[i - 1] = beta[i];
        }
        eigenvector = (double *)malloc(sizeof(double) * k * k);
        /* eigenvector consists of orthonormal columns */
        LAPACKE_dstev(LAPACK_ROW_MAJOR, 'V', k,
                      tmp_alpha, tmp_beta, eigenvector, k); 
        lambda_new = tmp_alpha[0];
        criteria = fabs((lambda_new - lambda_old) / lambda_new);
        lambda_old = lambda_new;
        free(tmp_alpha);
        free(tmp_beta);

        if (criteria < input->lambda_conv) {
            break;
        } else {
            free(eigenvector);
        }
    }

    *eigenvalue = lambda_new;
    double *x = (double *)malloc(sizeof(double) * k);
    for (i = 0; i < k; ++i) {
        x[i] = eigenvector[k * i];
    }
    double *y = (double *)malloc(sizeof(double) * disp_num * 3);
    double *basis = (double *)malloc(sizeof(double) * disp_num * 3 * k);
    for (i = 0; i < k; ++i) {
        for (j = 0; j < disp_num; ++j) {
            basis[(j * 3 + 0) * k + i] = V[i * disp_num * 3 + j * 3 + 0];
            basis[(j * 3 + 1) * k + i] = V[i * disp_num * 3 + j * 3 + 1];
            basis[(j * 3 + 2) * k + i] = V[i * disp_num * 3 + j * 3 + 2];
        }
    }
    cblas_dgemv(CblasRowMajor, CblasNoTrans, disp_num * 3, k,
                1.0, basis, k, x, 1, 0.0, y, 1);

    /* y is already almost norm, but... */
    double *norm_y = normalize(y, disp_num);
    for (i = 0; i < disp_num; ++i) {
        eigenmode[i * 3 + 0] = norm_y[i * 3 + 0];
        eigenmode[i * 3 + 1] = norm_y[i * 3 + 1];
        eigenmode[i * 3 + 2] = norm_y[i * 3 + 2];
    }

    free(eigenvector);
    free(basis);
    free(force0);
    free(force1);
    free(alpha);
    free(beta);
    free(V);
    free(Av);
    free(w);
    free(x);
    free(y);
    free(norm_y);
}


void uphill_push(Config *config, Input *input, int disp_num, int *disp_list,
                 double *init_direction, double eigenvalue, double *eigenmode,
                 MPI_Comm comm)
{
    int i;
    if (eigenvalue < input->lambda_crit) {
        double energy;
        double *force = (double *)malloc(sizeof(double) * disp_num * 3);
        oneshot_disp(config, input, &energy, force, disp_num, disp_list, comm);
        double *parallel_force = parallel_vector(force, eigenmode, disp_num);
        double f_norm = norm(parallel_force, disp_num);
        double divisor = abs(eigenvalue) > 0.5 ? abs(eigenvalue) : 0.5;
        double alpha = f_norm / divisor;
        double dr = input->max_step < alpha ? input->max_step : alpha;
        for (i = 0; i < disp_num; ++i) {
            config->pos[disp_list[i] * 3 + 0] += dr * eigenmode[i * 3 + 0];
            config->pos[disp_list[i] * 3 + 1] += dr * eigenmode[i * 3 + 1];
            config->pos[disp_list[i] * 3 + 2] += dr * eigenmode[i * 3 + 2];
        } 
        free(force);
        free(parallel_force);
    } else {
        for (i = 0; i < disp_num; ++i) {
            config->pos[disp_list[i] * 3 + 0] += init_direction[i * 3 + 0];
            config->pos[disp_list[i] * 3 + 1] += init_direction[i * 3 + 1];
            config->pos[disp_list[i] * 3 + 2] += init_direction[i * 3 + 2];
        }
    }
}


void normal_relax(Config *config0, Input *input, int disp_num, int *disp_list,
                  double eigenvalue, double *eigenmode, int count, int art_step,
                  MPI_Comm comm)
{
    int i, rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int local_rank = rank % input->ncore;


    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    double energy1;
    double *force1 = (double *)malloc(sizeof(double) * disp_num * 3);

    /* cg optimization */
    double *direction_old = (double *)malloc(sizeof(double) * disp_num * 3);
    double *cg_direction = (double *)malloc(sizeof(double) * disp_num * 3);

    int relax_step = 0;
    while (1) {
        oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm);
        double *perp_force0 = perpendicular_vector(force0, eigenmode, disp_num);

        if (local_rank == 0) {
            char filename[128];
            sprintf(filename, "%s/SPS_%d.log",
                    input->output_dir, count);
            FILE *fp = fopen(filename, "a");
            fprintf(fp, " %8d   %8d   %16f   %10f\n",
                    art_step, relax_step, energy0, eigenvalue);
            fclose(fp);
            sprintf(filename, "%s/SPS_%d.XDATCAR",
                    input->output_dir, count);
            write_config(config0, filename, "a");
        }
    
        if (eigenvalue < input->lambda_crit) {
            double *parallel_force0 = parallel_vector(force0, eigenmode, disp_num);
            if (norm(perp_force0, disp_num) < norm(parallel_force0, disp_num)) {
                free(perp_force0);
                free(parallel_force0);
                break;
            }
        } else {
            if (relax_step >= input->max_num_rlx) {
                free(perp_force0);
                break;
            }
        }

        if (relax_step == 1) {
            for (i = 0; i < disp_num; ++i) {
                direction_old[i * 3 + 0] = perp_force0[i * 3 + 0];
                direction_old[i * 3 + 1] = perp_force0[i * 3 + 1];
                direction_old[i * 3 + 2] = perp_force0[i * 3 + 2];
                cg_direction[i * 3 + 0] = perp_force0[i * 3 + 0];
                cg_direction[i * 3 + 1] = perp_force0[i * 3 + 1];
                cg_direction[i * 3 + 2] = perp_force0[i * 3 + 2];
            }
        }
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
        double *perp_force1 = perpendicular_vector(force1, eigenmode, disp_num);

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
        relax_step++;
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


int art_nouveau(Config *initial, Config *final, Input *input, Data *data,
                int count, int index, double *Ea, MPI_Comm comm)
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
    cut_sphere(initial, input, update_num, update_list);

    /* starting dimer */ 
    Config *config0 = (Config *)malloc(sizeof(Config));
    copy_config(config0, initial);

    /* set dimer space */
    double del[3];
    int disp_num = 0;
    int *disp_list = (int *)malloc(sizeof(int) * config0->tot_num);
    for (i = 0; i < config0->tot_num; ++i) {
        del[0] = config0->pos[i * 3 + 0] - center[0];
        del[1] = config0->pos[i * 3 + 1] - center[1];
        del[2] = config0->pos[i * 3 + 2] - center[2];
        get_minimum_image(del, config0->boxlo, config0->boxhi,
                          config0->xy, config0->yz, config0->xz);
        double dist = sqrt(del[0] * del[0]
                         + del[1] * del[1] 
                         + del[2] * del[2]);
        if (dist < input->acti_cutoff) {
            disp_list[disp_num] = i;
            disp_num++;
        }
    }

    /* eigenmode */
    double *full_eigenmode;
    if (data == NULL) {
        full_eigenmode = get_eigenmode(input, final->tot_num, comm); 
    } else {
        full_eigenmode = (double *)malloc(sizeof(double) * final->tot_num * 3);
        for (i = 0; i < final->tot_num; ++i) {
            full_eigenmode[i * 3 + 0] = data->eigenmode[i * 3 + 0];
            full_eigenmode[i * 3 + 1] = data->eigenmode[i * 3 + 1];
            full_eigenmode[i * 3 + 2] = data->eigenmode[i * 3 + 2];
        }
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

    /* perturbate starting config */
    double *init_direction = (double *)malloc(sizeof(double) * disp_num * 3);
    for (i = 0; i < disp_num; ++i) {
        init_direction[i * 3 + 0] = input->stddev * eigenmode[i * 3 + 0];
        init_direction[i * 3 + 1] = input->stddev * eigenmode[i * 3 + 1];
        init_direction[i * 3 + 2] = input->stddev * eigenmode[i * 3 + 2];
        config0->pos[disp_list[i] * 3 + 0] += init_direction[i * 3 + 0];
        config0->pos[disp_list[i] * 3 + 1] += init_direction[i * 3 + 1];
        config0->pos[disp_list[i] * 3 + 2] += init_direction[i * 3 + 2];
    }

    if (local_rank == 0) {
        char filename[128];
        sprintf(filename, "%s/SPS_%d.log",
                input->output_dir, count);
        FILE *fp = fopen(filename, "w");
        fputs("-----------------------------------------------------\n", fp);
        fputs(" Art step   Opt step   Potential energy   Eigenvalue\n", fp);
        fputs("-----------------------------------------------------\n", fp);
        fclose(fp);
        sprintf(filename, "%s/SPS_%d.XDATCAR",
                input->output_dir, count);
        write_config(config0, filename, "w");
    }

    int art_step;
    double eigenvalue;
    /* initial perturb */
    for (art_step = 1; art_step <= 1000; ++art_step) {
        /* lanczos */
        lanczos(config0, input, disp_num, disp_list,
                &eigenvalue, eigenmode, comm);

        /* test */
        double *H = get_hessian(config0, input, disp_num, disp_list, comm);
        double *w = get_eigenvalue(H, disp_num);

        printf("rank %d lanczos %f hessian %f\n", rank, eigenvalue, w[0]);
        fflush(stdout);
        free(H);
        free(w);

        /* uphill push */
        uphill_push(config0, input, disp_num, disp_list, init_direction,
                    eigenvalue, eigenmode, comm);
        /* normal relax */
        normal_relax(config0, input, disp_num, disp_list,
                     eigenvalue, eigenmode, count, art_step, comm);

        double energy;
        double *force = (double *)malloc(sizeof(double) * disp_num * 3);
        oneshot_disp(config0, input, &energy, force, disp_num, disp_list, comm);
        if (norm(force, disp_num) < input->f_tol) {
            free(force);
            break;
        } else {
            free(force);
        }
    }
    free(init_direction);

    /* split */
    int conv = split_configs(initial, final, config0, input,
                             eigenmode, count, index,
                             update_num, update_list,
                             disp_num, disp_list, comm);

    if ((local_rank == 0) && (conv == 0)) {
        char filename[128];
        sprintf(filename, "%s/ARTn_%d.log",
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
