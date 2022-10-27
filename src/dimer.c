#include <math.h>
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
#include "dimer.h"
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


static void rotate(Config *config0, Input *input, int disp_num, int *disp_list,
                   double *eigenmode, int count, int dimer_step, MPI_Comm comm)
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
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
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
        /* test */
        if ((local_rank == 0) && (input->write_mode)) {
            sprintf(filename, "%s/%d_%d_%d.MODECAR",
                    input->output_dir, count, dimer_step, i);
            FILE *fp = fopen(filename, "w");
            for (j = 0; j < disp_num; ++j) {
                fprintf(fp, "%f %f %f\n",
                        eigenmode[j * 3 + 0],
                        eigenmode[j * 3 + 1],
                        eigenmode[j * 3 + 2]);
            }
            fclose(fp);
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
    free(force0);
    free(force1);
    free(force2);
}


static void translate(Config *config0, Input *input,
                      int disp_num, int *disp_list, double *eigenmode,
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
    oneshot_disp(config0, input, &energy0, force0, disp_num, disp_list, comm); 
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
    double *f0p = projected_force(force0, eigenmode, curvature, disp_num);
    /* cg_direction */
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


int dimer(Config *initial, Config *final, Input *input, double *full_eigenmode,
          int count, int index, double *Ea, MPI_Comm comm)
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

    double energy0;
    double *force0 = (double *)malloc(sizeof(double) * disp_num * 3);
    for (dimer_step = 1; dimer_step <= 1000; ++dimer_step) {
        rotate(config0, input, disp_num, disp_list,
               eigenmode, count, dimer_step, comm);
        translate(config0, input, disp_num, disp_list, eigenmode,
                  direction_old, cg_direction, dimer_step, comm);
        oneshot_disp(config0, input, &energy0, force0,
                     disp_num, disp_list, comm);     
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
            sprintf(filename, "%s/SPS_%d.XDATCAR",
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
