#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "calculator.h"
#include "config.h"
#include "dataset.h"
#include "dimer.h"
#include "input.h"
#include "target.h"
#include "utils.h"


void recycle(Config *saddle, Config *config_new, Config *config_old,
             Input *input, Data *data, double *eigenmode, MPI_Comm comm)
{
    int i, rank, size;
    double del[3];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int group_size = size / input->ncore;
    int group_rank = rank / input->ncore;
    int local_rank = rank % input->ncore;

    int q = saddle->tot_num / input->ncore;
    int r = saddle->tot_num % input->ncore;
    int begin = local_rank * q + ((local_rank > r) ? r : local_rank);
    int end = begin + q;
    if (r > local_rank) {
        end++;
    }

    for (i = begin; i < end; ++i) {
        del[0] = config_new->pos[i * 3 + 0] - config_old->pos[i * 3 + 0];
        del[1] = config_new->pos[i * 3 + 1] - config_old->pos[i * 3 + 1];
        del[2] = config_new->pos[i * 3 + 2] - config_old->pos[i * 3 + 2];
        get_minimum_image(del, saddle->boxlo, saddle->boxhi,
                          saddle->xy, saddle->yz, saddle->xz);
        double dist = sqrt(del[0] * del[0]
                         + del[1] * del[1] 
                         + del[2] * del[2]);
        if (dist < 2 * input->max_step) {
            saddle->pos[i * 3 + 0] = data->saddle[i * 3 + 0];
            saddle->pos[i * 3 + 1] = data->saddle[i * 3 + 1];
            saddle->pos[i * 3 + 2] = data->saddle[i * 3 + 2];
            eigenmode[i * 3 + 0] = data->eigenmode[i * 3 + 0];
            eigenmode[i * 3 + 1] = data->eigenmode[i * 3 + 1];
            eigenmode[i * 3 + 2] = data->eigenmode[i * 3 + 2];
        } else {
            saddle->pos[i * 3 + 0] = config_new->pos[i * 3 + 0];
            saddle->pos[i * 3 + 1] = config_new->pos[i * 3 + 1];
            saddle->pos[i * 3 + 2] = config_new->pos[i * 3 + 2];
            eigenmode[i * 3 + 0] = 0.0;
            eigenmode[i * 3 + 1] = 0.0;
            eigenmode[i * 3 + 2] = 0.0;
        }
    }
    int count = (end - begin) * 3;
    int *counts = (int *)malloc(sizeof(int) * input->ncore);
    MPI_Allgather(&count, 1, MPI_INT, counts, 1, MPI_INT, comm);
    int *disp = (int *)malloc(sizeof(int) * input->ncore);
    disp[0] = 0;
    if (input->ncore > 1) {
        for (i = 1; i < input->ncore; ++i) {
            disp[i] = disp[i - 1] + counts[i - 1];
        }
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   saddle->pos, counts, disp, MPI_DOUBLE, comm);
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   eigenmode, counts, disp, MPI_DOUBLE, comm);
    free(counts);
    free(disp);
}


int main(int argc, char *argv[])
{
    int i, j, errno, rank, size;
    int buffer_size = 4096;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* read input */
    Input *input = (Input *)malloc(sizeof(Input));
    errno = read_input(input, "./INPUT");
    if (errno > 0) {
        printf("ERROR in INPUT FILE!\n");
        free(input);
        MPI_Finalize();
        return 1;
    }

    /* make directory */
    if (rank == 0) {
        errno = mkdir(input->output_dir, 0775);
    }
    MPI_Bcast(&errno, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (errno != 0) {
        if (rank == 0) {
            printf("OUTPUT_DIR exists!\n");
            printf("We will not overwrite your data.\n");
        }
        free_input(input);
        MPI_Finalize();
        return 1;
    }

    /* random number */
    srand(input->random_seed);
    if (rank == 0) {
        char line[64], filename[64];
        sprintf(filename, "%s/REPORT", input->output_dir);
        FILE *fp = fopen(filename, "w"); 
        sprintf(line, "Random seed: %d\n", input->random_seed);
        fputs(line, fp);
        fclose(fp);
    }

    /* read config */
    Config *config = (Config *)malloc(sizeof(Config));
    errno = read_config(config, input, input->init_config);
    if (errno > 0) {
        printf("ERROR in INIT_CONFIG FILE!\n");
        free_input(input);
        free(config);
        MPI_Finalize();
        return 1;
    }

    /* initial relax */
    if (input->init_relax > 0) {
        atom_relax(config, input, MPI_COMM_WORLD);
    }

    long long step;
    double time = 0.0;
    /* log */
    if (rank == 0) {
        char line[128], filename[128];
        sprintf(filename, "%s/kMC.log", input->output_dir);
        FILE *fp = fopen(filename, "a");
        fputs("--------------------------------------------------------------\n", fp);
        fputs(" kMC step   Potential energy   Reaction barrier      kMC time\n", fp);
        fputs("--------------------------------------------------------------\n", fp);
        fclose(fp);
        printf("--------------------------------------------------------------\n");
        printf(" kMC step   Reaction index   Reaction barrier   Reaction rate\n");
    }

    /* one-sided communication */
    int group_size = size / input->ncore;
    int group_rank = rank / input->ncore;
    int local_rank = rank % input->ncore;
    MPI_Comm local_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group_rank, rank, &local_comm);
    int head;
    if (local_rank == 0) {
        head = 1;
    } else {
        head = 0;
    }
    MPI_Comm group_comm;
    MPI_Comm_split(MPI_COMM_WORLD, head, rank, &group_comm);

    MPI_Win count_win;
    int *global_count;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_count, &count_win);

    MPI_Win redundant_win;
    int *global_redundant;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_redundant, &redundant_win);

    /* for recycle */
    Config *config_old = (Config *)malloc(sizeof(Config));
    copy_config(config_old, config);

    /* dataset */
    Dataset *dataset = (Dataset *)malloc(sizeof(Dataset));
    dataset->numdata = 0;
    dataset->head = NULL;

    /* kmc loop */
    for (step = 1; step <= input->end_step; ++step) {
        int unique;
        double Ea;
        double *eigenmode;

        /* read target */
        int target_num = 0;
        int list_size = 64;
        int *target_list = (int *)malloc(sizeof(int) * list_size);
        errno = get_target(config, input, &target_list, &target_num, &list_size);
        if (errno > 0) {
            printf("ERROR in TARGET FILE!\n");
            free_input(input);
            free_config(config);
            MPI_Finalize();
            return 1;
        }
        /* dimer loop */
        int unique_num = (dataset->numdata > 0) ? dataset->numdata : 1;
        while (1) {
            if (unique_num < target_num) {
                unique_num = unique_num << 1;
            } else {
                break;
            }
        }
        int local_count;
        *global_count = 0;
        int count = 1;
        int zero = 0;
        int local_reac_num = 0;
        double local_rate_sum = 0.0;
        int *local_reac_list = (int *)malloc(sizeof(int) * unique_num);
        double *local_acti_list = (double *)malloc(sizeof(double) * unique_num);
        double *local_rate_list = (double *)malloc(sizeof(double) * unique_num);
        /* all data in dataset are unique */
        /* do not have to count redundant search */
        if (dataset->numdata > 0) {
            while (1) {
                /* one-side communication */
                if (local_rank == 0) {
                    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, count_win);
                    MPI_Fetch_and_op(&count, &local_count, MPI_INT,
                                     0, (MPI_Aint)0, MPI_SUM, count_win);
                    MPI_Win_unlock(0, count_win);
                }
                MPI_Bcast(&local_count, 1, MPI_INT, 0, local_comm);
                Data *data = dataset->head;
                if (local_count > 0) {
                    for (i = 0; i < local_count; ++i) {
                        data = data->next;
                        if (data == NULL) {
                            break;
                        }
                    }
                    if (data == NULL) {
                        break;
                    }
                }
                /* initial/saddle/final configuration */
                Config *initial = (Config *)malloc(sizeof(Config));
                copy_config(initial, config);
                Config *saddle = (Config *)malloc(sizeof(Config));
                copy_config(saddle, config);
                Config *final = (Config *)malloc(sizeof(Config));
                copy_config(final, config);
                eigenmode = (double *)malloc(sizeof(double) * config->tot_num * 3);
                /* recycle */
                recycle(saddle, config, config_old,
                        input, data, eigenmode, local_comm);
                int conv = dimer(initial, saddle, final, input, local_count,
                                 data->index, eigenmode, &Ea, local_comm);
                /* conv == 0 -> success */
                if (conv == 0) {
                    if (local_rank == 0) {
                        char line[128], filename[128];
                        sprintf(filename, "%s/Dimer_%d.log",
                                input->output_dir, local_count);
                        FILE *fp = fopen(filename, "a");
                        sprintf(line, " Barrier energy: %f eV\n", Ea);
                        fputs(line, fp);
                        fclose(fp);
                    }
                    unique = check_unique(final, input);
                    if (unique > 0) {
                        local_reac_list[local_reac_num] = local_count;
                        local_acti_list[local_reac_num] = Ea;
                        double kT = 8.61733034e-5 * input->temperature;
                        local_rate_list[local_reac_num] = input->att_freq * exp(- Ea / kT);
                        local_rate_sum += local_rate_list[local_reac_num];
                        local_reac_num++;
                        if (local_reac_num > unique_num) {
                            unique_num = unique_num << 1;
                            local_reac_list = (int *)realloc(local_reac_list,
                                                             sizeof(int) * unique_num);
                            local_acti_list = (double *)realloc(local_acti_list,
                                                             sizeof(double) * unique_num);
                            local_rate_list = (double *)realloc(local_rate_list,
                                                             sizeof(double) * unique_num);
                        }
                        if (local_rank == 0) {
                            char filename[128];
                            sprintf(filename, "%s/Final_%d.POSCAR",
                                    input->output_dir, local_count);
                            write_config(final, filename, "w");
                        }
                    }
                }
                free(initial);
                free(saddle);
                free(final);
                free(eigenmode);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int local_redundant;
        *global_redundant = 0;
        while (1) {
            /* confidence check */
            if (local_rank == 0) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
                MPI_Fetch_and_op(&zero, &local_redundant, MPI_INT,
                                 0, (MPI_Aint)0, MPI_SUM, redundant_win);
                MPI_Win_unlock(0, redundant_win);
            }
            MPI_Bcast(&local_redundant, 1, MPI_INT, 0, local_comm);
            if (local_redundant >= input->nredundant) {
                break;
            }
            /* one-sided communication */
            if (local_rank == 0) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, count_win);
                MPI_Fetch_and_op(&count, &local_count, MPI_INT,
                                 0, (MPI_Aint)0, MPI_SUM, count_win);
                MPI_Win_unlock(0, count_win);
            }
            MPI_Bcast(&local_count, 1, MPI_INT, 0, local_comm);
            /* atom index */
            int atom_index = target_list[rand() % target_num];
            MPI_Bcast(&atom_index, 1, MPI_INT, 0, local_comm);
            /* initial/saddle/final configuration */
            Config *initial = (Config *)malloc(sizeof(Config));
            copy_config(initial, config);
            Config *saddle = (Config *)malloc(sizeof(Config));
            copy_config(saddle, config);
            Config *final = (Config *)malloc(sizeof(Config));
            copy_config(final, config);
            /* generate not normalized eigenmode */
            eigenmode = gen_eigenmode(input, config->tot_num, local_comm);
            int conv = dimer(initial, saddle, final, input, local_count,
                             atom_index, eigenmode, &Ea, local_comm);
            /* conv == 0 -> success */
            if (conv == 0) {
                if (local_rank == 0) {
                    char line[128], filename[128];
                    sprintf(filename, "%s/Dimer_%d.log",
                            input->output_dir, local_count);
                    FILE *fp = fopen(filename, "a");
                    sprintf(line, " Barrier energy: %f eV\n", Ea);
                    fputs(line, fp);
                    fclose(fp);
                }
                unique = check_unique(final, input);
                if (unique > 0) {
                    if (local_rank == 0) {
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
                        MPI_Fetch_and_op(&zero, &local_redundant, MPI_INT,
                                         0, (MPI_Aint)0, MPI_REPLACE, redundant_win);
                        MPI_Win_unlock(0, redundant_win);
                    }
                    local_reac_list[local_reac_num] = local_count;
                    local_acti_list[local_reac_num] = Ea;
                    double kT = 8.61733034e-5 * input->temperature;
                    local_rate_list[local_reac_num] = input->att_freq * exp(-Ea / kT);
                    local_rate_sum += local_rate_list[local_reac_num];
                    local_reac_num++;
                    if (local_reac_num > unique_num) {
                        unique_num = unique_num << 1;
                        local_reac_list = (int *)realloc(local_reac_list,
                                                         sizeof(int) * unique_num);
                        local_acti_list = (double *)realloc(local_acti_list,
                                                         sizeof(double) * unique_num);
                        local_rate_list = (double *)realloc(local_rate_list,
                                                         sizeof(double) * unique_num);
                    }
                    if (local_rank == 0) {
                        char filename[128];
                        sprintf(filename, "%s/Final_%d.POSCAR",
                                input->output_dir, local_count);
                        write_config(final, filename, "w");
                    }
                } else {
                    if (local_rank == 0) {
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
                        MPI_Fetch_and_op(&count, &local_redundant, MPI_INT,
                                         0, (MPI_Aint)0, MPI_SUM, redundant_win);
                        MPI_Win_unlock(0, redundant_win);
                    }
                }
            }
            free_config(initial);
            free_config(saddle);
            free_config(final);
            free(eigenmode);
        }

        int total_reac_num;
        double total_rate_sum;
        if (local_rank == 0) {
            MPI_Allreduce(&local_rate_sum, &total_rate_sum, 1, MPI_DOUBLE,
                          MPI_SUM, group_comm);
            MPI_Allreduce(&local_reac_num, &total_reac_num, 1, MPI_INT,
                          MPI_SUM, group_comm);
        }
        MPI_Bcast(&total_reac_num, 1, MPI_INT, 0, local_comm);

        int *global_reac_num = (int *)malloc(sizeof(int) * group_size);
        int *global_reac_list = (int *)malloc(sizeof(int) * total_reac_num);
        double *global_acti_list = (double *)malloc(sizeof(double) * total_reac_num);
        double *global_rate_list = (double *)malloc(sizeof(double) * total_reac_num);
        if (local_rank == 0) {
            MPI_Allgather(&local_reac_num, 1, MPI_INT,
                          global_reac_num, 1, MPI_INT, group_comm);
            int *disp = (int *)malloc(sizeof(int) * group_size);
            disp[0] = 0;
            if (group_size > 1) {
                for (int i = 1; i < group_size; ++i) {
                    disp[i] = disp[i - 1] + global_reac_num[i - 1];
                }
            }
            MPI_Gatherv(local_acti_list, local_reac_num, MPI_DOUBLE,
                        global_acti_list, global_reac_num, disp, MPI_DOUBLE,
                        0, group_comm);
            MPI_Gatherv(local_rate_list, local_reac_num, MPI_DOUBLE,
                        global_rate_list, global_reac_num, disp, MPI_DOUBLE,
                        0, group_comm);
            MPI_Gatherv(local_reac_list, local_reac_num, MPI_INT,
                        global_reac_list, global_reac_num, disp, MPI_INT,
                        0, group_comm);
            free(disp);
        }

        /* selection */
        int reac_index;
        if (rank == 0) {
            double random1 = (double)rand() / RAND_MAX;
            double random2 = (double)rand() / RAND_MAX;
            double rate_acc = 0.0;
            for (i = 0; i < total_reac_num; ++i) {
                rate_acc += global_rate_list[i] / total_rate_sum;
                if (rate_acc > random1) {
                    reac_index = i;
                    break;
                } 
            }
            time -= log(random2) / total_rate_sum;
        }

        MPI_Bcast(&time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double energy = oneshot(config, input, MPI_COMM_WORLD);
        /* log */
        if (rank == 0) {
            char line[128], filename[128];
            sprintf(filename, "%s/kMC.log", input->output_dir);
            FILE *fp = fopen(filename, "a");
            sprintf(line, " %8lld   %16f   %16f   %8e\n",
                    step, energy, global_acti_list[reac_index], time);
            fputs(line, fp);
            fclose(fp);
            sprintf(filename, "%s/kMC.XDATCAR", input->output_dir);
            write_config(config, filename, "a");
            printf("--------------------------------------------------------------\n");
            for (i = 0; i < total_reac_num; ++i) {
                printf(" %8lld   %14d   %16f   %13e\n",
                       step, global_reac_list[i],
                       global_acti_list[i], global_rate_list[i]);
            }
        }
        /* free and build dataset */
        free_dataset(dataset);
        dataset = (Dataset *)malloc(sizeof(Dataset));
        dataset->numdata = 0;
        dataset->head = NULL;
        build_dataset(dataset, config, input);

        free_config(config_old);
        config_old = (Config *)malloc(sizeof(Config));
        copy_config(config_old, config);
        free_config(config);
        free(target_list);
        if (time > input->end_time) {
            free(local_reac_list);
            free(local_acti_list);
            free(local_rate_list);
            free(global_reac_num);
            free(global_reac_list);
            free(global_acti_list);
            free(global_rate_list);
            break;
        }
        char filename[128];
        if (rank == 0) {
            sprintf(filename, "%s/Final_%d.POSCAR",
                    input->output_dir, global_reac_list[reac_index]);
        }
        MPI_Bcast(filename, 128, MPI_CHAR, 0, MPI_COMM_WORLD);

        free(local_reac_list);
        free(local_acti_list);
        free(local_rate_list);
        free(global_reac_num);
        free(global_reac_list);
        free(global_acti_list);
        free(global_rate_list);

        config = (Config *)malloc(sizeof(Config));
        errno = read_config(config, input, filename);
        if (errno > 0) {
            printf("ERROR in INIT_CONFIG FILE!\n");
            free_input(input);
            free(config);
            MPI_Win_free(&count_win);
            MPI_Win_free(&redundant_win);
            MPI_Comm_free(&local_comm);
            MPI_Comm_free(&group_comm);
            MPI_Finalize();
            return 1;
        }
    }
    free_input(input);
    free_config(config_old);
    MPI_Win_free(&count_win);
    MPI_Win_free(&redundant_win);
    MPI_Comm_free(&local_comm);
    MPI_Comm_free(&group_comm);
    MPI_Finalize();
    return 0;
}
