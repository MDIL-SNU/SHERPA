#include "art_nouveau.h"
#include "config.h"
#include "dataset.h"
#include "dimer.h"
#include "input.h"
#include "kappa_dimer.h"
#include "linalg.h"
#include "my_mpi.h"
#include "target.h"
#include "utils.h"
#include "calculator.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>


int main(int argc, char *argv[])
{
    int i, j, atom_index, errno, rank, size;
    FILE *fp;

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
    /* write_input */
    if (rank == 0) {
        write_input(input);
    }

    /* random number */
    srand(input->random_seed + rank);

    /* read config */
    Config *config = (Config *)malloc(sizeof(Config));
    errno = read_config(config, "./POSCAR");
    if (errno > 0) {
        printf("ERROR in POSCAR FILE!\n");
        free_input(input);
        free(config);
        MPI_Finalize();
        return 1;
    }

    /* read target */
    int target_num = 0;
    int list_size = 64;
    int *target_list = (int *)malloc(sizeof(int) * list_size);
    errno = read_target(config, &target_num, &target_list, &list_size);
    MPI_Bcast(target_list, target_num, MPI_INT, 0, MPI_COMM_WORLD);
    if (errno > 0) {
        printf("ERROR in TARGET FILE!\n");
        free_input(input);
        free_config(config);
        MPI_Finalize();
        return 1;
    }
    if (rank == 0) {
        write_target(target_num, target_list);
    }

    /* continue */
    if (input->cont > 0) {
        fp = fopen("./Statistics.log", "r");
        if (fp == NULL) {
            printf("Cannot find Statistics.log.");
            free_input(input);
            free_config(config);
            free(target_list);
            MPI_Finalize();
            return 1;
        } else {
            fclose(fp);
        }
        fp = fopen("./Event.log", "r");
        if (fp == NULL) {
            printf("Cannot find Event.log");
            free_input(input);
            free_config(config);
            free(target_list);
            MPI_Finalize();
            return 1;
        } else {
            fclose(fp);
        }
    }

    /* initial relax */
    if (input->init_relax > 0) {
        double energy;
        atom_relax(config, input, &energy, MPI_COMM_WORLD);
        if (rank == 0) {
            write_config(config, "./POSCAR_read", "INIT_CONFIG", "w");
        }
    }

    int group_size = size / input->ncore;
    int group_rank = rank / input->ncore;
    int local_rank = rank % input->ncore;
    /* intra-group communication */
    MPI_Comm local_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group_rank, rank, &local_comm);
    /* inter-group communication */
    MPI_Comm group_comm;
    MPI_Comm_split(MPI_COMM_WORLD, local_rank, rank, &group_comm);
    /* reset ncore for the last group */
    if (rank >= group_size * input->ncore) {
        input->ncore = size - group_size * input->ncore;
    }
    /* one-sided communication */
    MPI_Win count_win;
    int global_count = 0;
    MPI_Win_create(&global_count, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &count_win);

    MPI_Win redundant_win;
    int global_redundant = 0;
    MPI_Win_create(&global_redundant, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &redundant_win);

    MPI_Win done_win;
    int global_done = 0;
    MPI_Win_create(&global_done, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &done_win);

    MPI_Win conv_win;
    int global_conv = 0;
    MPI_Win_create(&global_conv, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &conv_win);

    MPI_Win unique_win;
    int global_unique = 0;
    MPI_Win_create(&global_unique, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &unique_win);

    MPI_Win write_win;
    int global_write = 0;
    MPI_Win_create(&global_write, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &write_win);

    MPI_Win exit_win;
    int global_exit = 0;
    MPI_Win_create(&global_exit, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, group_comm, &exit_win);

    int zero = 0;
    int one = 1;
    int local_count;
    int local_done;
    int local_conv;
    int local_unique;
    int local_write;
    int local_redundant;
    int local_exit;
    int local_reac_num = 0;
    int local_freq_num = 0;
    int *local_reac_list = (int *)malloc(sizeof(int) * list_size);
    int *local_freq_list = (int *)malloc(sizeof(int) * list_size);
    double *local_acti_list = (double *)malloc(sizeof(double) * list_size);

    /* log */
    if (rank == 0) {
        if (input->cont == 0) {
            fp = fopen("./Redundancy.log", "w");
            fputs("-----------------\n", fp);
            fputs(" Earlier   Later\n", fp);
            fputs("-----------------\n", fp);
            fclose(fp);
            fp = fopen("./Statistics.log", "w");
            fputs("------------------------------------------\n", fp);
            fputs(" Unique events   Relevant events   Trials\n", fp);
            fputs("------------------------------------------\n", fp);
            fclose(fp);
        } else {
            char tmp_line[1024], line[1024], *ptr;
            fp = fopen("./Statistics.log", "r");
            fgets(tmp_line, 1024, fp);
            fgets(tmp_line, 1024, fp);
            fgets(tmp_line, 1024, fp);
            ptr = fgets(tmp_line, 1024, fp);
            while (ptr != NULL) {
                strcpy(line, tmp_line);
                ptr = fgets(tmp_line, 1024, fp);
            }
            global_unique = atoi(strtok(line, " \n"));
            global_conv = atoi(strtok(NULL, " \n"));
            global_count = atoi(strtok(NULL, " \n"));
            global_done = global_count;
            fclose(fp);
            fp = fopen("./Event.log", "r");
            fgets(line, 1024, fp);
            fgets(line, 1024, fp);
            fgets(line, 1024, fp);
            ptr = fgets(line, 1024, fp);
            while (ptr != NULL) {
                local_reac_list[local_reac_num] = atoi(strtok(line, " \n"));
                local_acti_list[local_reac_num] = atof(strtok(NULL, " \n"));
                int frequency = atoi(strtok(NULL, " \n"));
                for (i = 0; i < frequency; ++i) {
                    local_freq_list[local_freq_num] = local_reac_list[local_reac_num];
                    local_freq_num++;
                    if (local_freq_num >= list_size) {
                        list_size = list_size << 1;
                        local_freq_list = (int *)realloc(local_freq_list,
                                                 sizeof(int) * list_size);
                    }
                }
                local_reac_num++;
                if (local_reac_num >= list_size) {
                    list_size = list_size << 1;
                    local_reac_list = (int *)realloc(local_reac_list,
                                             sizeof(int) * list_size);
                    local_acti_list = (double *)realloc(local_acti_list,
                                                sizeof(double) * list_size);
                }
                ptr = fgets(line, 1024, fp);
            }
            fclose(fp);
        }
        fp = fopen("./Event.log", "w");
        fputs("---------------------------------------------\n", fp);
        fputs(" Reaction index   Barrier energy   Frequency\n", fp);
        fputs("---------------------------------------------\n", fp);
        fclose(fp);
    }
    MPI_Bcast(&global_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
    input->max_search += global_done;

    /* dataset */
    Dataset *dataset = (Dataset *)malloc(sizeof(Dataset));
    dataset->head = NULL;
    if (input->init_mode > 0) {
        build_dataset(dataset, config->tot_num);
    }

    int conv, unique;
    double Ea, *eigenmode;
    while (1) {
        if (local_rank == 0) {
            /* increase count */
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, count_win);
            MPI_Fetch_and_op(&one, &local_count, MPI_INT,
                             0, (MPI_Aint)0, MPI_SUM, count_win);
            MPI_Win_unlock(0, count_win);
            /* check exit condition */
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
            MPI_Fetch_and_op(&zero, &local_redundant, MPI_INT,
                             0, (MPI_Aint)0, MPI_SUM, redundant_win);
            MPI_Win_unlock(0, redundant_win);
            if ((local_redundant >= input->nredundant) ||
                (local_count >= input->max_search)) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, exit_win);
                MPI_Fetch_and_op(&one, &local_exit, MPI_INT,
                                 0, (MPI_Aint)0, MPI_SUM, exit_win);
                MPI_Win_unlock(0, exit_win);
            }
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, exit_win);
            MPI_Fetch_and_op(&zero, &local_exit, MPI_INT,
                             0, (MPI_Aint)0, MPI_SUM, exit_win);
            MPI_Win_unlock(0, exit_win);
        }
        MPI_Bcast(&local_exit, 1, MPI_INT, 0, local_comm);
        if (local_exit > 0) {
            break;
        }
        MPI_Bcast(&local_count, 1, MPI_INT, 0, local_comm);

        /* atom_index */
        Data *data = dataset->head;
        if ((data != NULL) && (local_count > 0)) {
            for (i = 0; i < local_count; ++i) {
                data = data->next;
                if (data == NULL) {
                    break;
                }
            }
        }
        if (data == NULL) {
            atom_index = target_list[local_count % target_num];
            eigenmode = NULL;
        } else {
            atom_index = data->index;
            eigenmode = data->eigenmode;
        }

        /* initial/saddle configuration */
        Config *initial = (Config *)malloc(sizeof(Config));
        copy_config(initial, config);
        Config *saddle = (Config *)malloc(sizeof(Config));
        copy_config(saddle, config);
        Config *final = (Config *)malloc(sizeof(Config));
        copy_config(final, config);
        if (input->art_nouveau > 0) {
            conv = art_nouveau(initial, saddle, final, input, eigenmode,
                               local_count, atom_index, &Ea, local_comm);
        } else if (input->kappa_dimer > 0) {
            conv = kappa_dimer(initial, saddle, final, input, eigenmode,
                               local_count, atom_index, &Ea, local_comm);
        } else {
            conv = dimer(initial, saddle, final, input, eigenmode,
                         local_count, atom_index, &Ea, local_comm);
        }
        if (local_rank == 0) {
            /* conv == -1 -> Unconverged  */
            /* conv ==  0 -> Connected    */
            /* conv ==  1 -> Disconnected */
            /* conv ==  2 -> Not splited  */
            while (1) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, write_win);
                MPI_Fetch_and_op(&one, &local_write, MPI_INT,
                                 0, (MPI_Aint)0, MPI_SUM, write_win);
                MPI_Win_unlock(0, write_win);
                if (local_write == 0) { 
                    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, done_win);
                    MPI_Fetch_and_op(&one, &local_done, MPI_INT,
                                     0, (MPI_Aint)0, MPI_SUM, done_win);
                    MPI_Win_unlock(0, done_win);
                    char filename[128], header[128];
                    sprintf(header, "%d_%d", local_count, atom_index);
                    if (conv >= 0) {
                        if (conv == 0) {
                            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, conv_win);
                            MPI_Fetch_and_op(&one, &local_conv, MPI_INT,
                                             0, (MPI_Aint)0, MPI_SUM, conv_win);
                            MPI_Win_unlock(0, conv_win);
                            /* unique < 0 -> unique */
                            unique = check_unique(saddle, input);
                            if (unique < 0) {
                                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, unique_win);
                                MPI_Fetch_and_op(&one, &local_unique, MPI_INT,
                                                 0, (MPI_Aint)0, MPI_SUM, unique_win);
                                MPI_Win_unlock(0, unique_win);
                                local_reac_list[local_reac_num] = local_count;
                                local_acti_list[local_reac_num] = Ea;
                                local_reac_num++;
                                if (local_reac_num >= list_size) {
                                    list_size = list_size << 1;
                                    local_reac_list = (int *)realloc(local_reac_list,
                                                             sizeof(int) * list_size);
                                    local_acti_list = (double *)realloc(local_acti_list,
                                                                sizeof(double) * list_size);
                                }
                                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
                                MPI_Fetch_and_op(&zero, &local_redundant, MPI_INT,
                                                 0, (MPI_Aint)0, MPI_REPLACE, redundant_win);
                                MPI_Win_unlock(0, redundant_win);
                            } else {
                                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
                                MPI_Fetch_and_op(&one, &local_redundant, MPI_INT,
                                                 0, (MPI_Aint)0, MPI_SUM, redundant_win);
                                MPI_Win_unlock(0, redundant_win);
                                local_freq_list[local_freq_num] = unique;
                                local_freq_num++;
                                if (local_freq_num >= list_size) {
                                    list_size = list_size << 1;
                                    local_freq_list = (int *)realloc(local_freq_list,
                                                             sizeof(int) * list_size);
                                }
                                fp = fopen("./Redundancy.log", "a");
                                fprintf(fp, " %7d   %5d\n", unique, local_count);
                                fclose(fp);
                            }
                            /* Final */
                            write_config(final, "./Final.POSCAR", header, "a");
                        }
                        /* Saddle */
                        write_config(saddle, "./Saddle.POSCAR", header, "a");
                        /* MODECAR */
                        sprintf(filename, "./%d.MODECAR", local_count);
                        concat_files("./SPS.MODECAR", filename);
                        remove(filename);
                    }
                    /* XDATCAR */
                    sprintf(filename, "./%d.XDATCAR", local_count);
                    concat_files("./SPS.XDATCAR", filename);
                    remove(filename);
                    /* log */
                    sprintf(filename, "./%d.log", local_count);
                    concat_files("./SPS.log", filename);
                    remove(filename);
                    /* update statistics */
                    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, done_win);
                    MPI_Fetch_and_op(&zero, &local_done, MPI_INT,
                                     0, (MPI_Aint)0, MPI_SUM, done_win);
                    MPI_Win_unlock(0, done_win);
                    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, conv_win);
                    MPI_Fetch_and_op(&zero, &local_conv, MPI_INT,
                                     0, (MPI_Aint)0, MPI_SUM, conv_win);
                    MPI_Win_unlock(0, conv_win);
                    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, unique_win);
                    MPI_Fetch_and_op(&zero, &local_unique, MPI_INT,
                                     0, (MPI_Aint)0, MPI_SUM, unique_win);
                    MPI_Win_unlock(0, unique_win);
                    fp = fopen("./Statistics.log", "a");
                    fprintf(fp, " %13d   %15d   %6d\n",
                            local_unique, local_conv, local_done);
                    fclose(fp);
                    /* write done */
                    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, write_win);
                    MPI_Fetch_and_op(&zero, &local_write, MPI_INT,
                                     0, (MPI_Aint)0, MPI_REPLACE, write_win);
                    MPI_Win_unlock(0, write_win);
                    break;
                } else {
                    sleep(1);
                }
            }
        }
        free_config(initial);
        free_config(saddle);
        free_config(final);
    }
    if (local_rank == 0) {
        int total_reac_num;
        int total_freq_num;
        MPI_Allreduce(&local_reac_num, &total_reac_num, 1, MPI_INT,
                      MPI_SUM, group_comm);
        MPI_Allreduce(&local_freq_num, &total_freq_num, 1, MPI_INT,
                      MPI_SUM, group_comm);

        int *global_reac_num = (int *)malloc(sizeof(int) * group_size);
        int *global_freq_num = (int *)malloc(sizeof(int) * group_size);
        int *global_reac_list = (int *)malloc(sizeof(int) * total_reac_num);
        int *global_freq_list = (int *)malloc(sizeof(int) * total_freq_num);
        double *global_acti_list = (double *)malloc(sizeof(double) * total_reac_num);
        /* reaction */
        MPI_Allgather(&local_reac_num, 1, MPI_INT,
                      global_reac_num, 1, MPI_INT, group_comm);
        int *disp = (int *)malloc(sizeof(int) * group_size);
        disp[0] = 0;
        if (group_size > 1) {
            for (i = 1; i < group_size; ++i) {
                disp[i] = disp[i - 1] + global_reac_num[i - 1];
            }
        }
        MPI_Gatherv(local_acti_list, local_reac_num, MPI_DOUBLE,
                    global_acti_list, global_reac_num, disp, MPI_DOUBLE,
                    0, group_comm);
        MPI_Gatherv(local_reac_list, local_reac_num, MPI_INT,
                    global_reac_list, global_reac_num, disp, MPI_INT,
                    0, group_comm);

        /* frequency */
        MPI_Allgather(&local_freq_num, 1, MPI_INT,
                      global_freq_num, 1, MPI_INT, group_comm);
        disp[0] = 0;
        if (group_size > 1) {
            for (i = 1; i < group_size; ++i) {
                disp[i] = disp[i - 1] + global_freq_num[i - 1];
            }
        }
        MPI_Gatherv(local_freq_list, local_freq_num, MPI_INT,
                    global_freq_list, global_freq_num, disp, MPI_INT,
                    0, group_comm);
        free(disp);

        /* log */
        if (rank == 0) {
            int *freq_num = (int *)calloc(total_reac_num, sizeof(int));
            for (i = 0; i < total_reac_num; ++i) {
                for (j = 0; j < total_freq_num; ++j) {
                    if (global_freq_list[j] == global_reac_list[i]) {
                        freq_num[i]++;
                    }
                }
            }
            fp = fopen("./Redundancy.log", "a");
            fclose(fp);
            fp = fopen("./Statistics.log", "a");
            fclose(fp);
            fp = fopen("./Event.log", "a");
            for (i = 0; i < total_reac_num; ++i) {
                fprintf(fp, " %14d   %14f   %9d\n",
                        global_reac_list[i], global_acti_list[i], freq_num[i] + 1);
            }
            fclose(fp);
            free(freq_num);
        }
        free(global_reac_num);
        free(global_freq_num);
        free(global_reac_list);
        free(global_acti_list);
        free(global_freq_list);
    }

    free(target_list);
    free(local_reac_list);
    free(local_acti_list);
    free(local_freq_list);

    free_dataset(dataset);
    free_config(config);
    free_input(input);

    MPI_Win_free(&count_win);
    MPI_Win_free(&redundant_win);
    MPI_Win_free(&done_win);
    MPI_Win_free(&conv_win);
    MPI_Win_free(&unique_win);
    MPI_Win_free(&write_win);
    MPI_Win_free(&exit_win);

    MPI_Comm_free(&local_comm);
    MPI_Comm_free(&group_comm);
    MPI_Finalize();
    return 0;
}
