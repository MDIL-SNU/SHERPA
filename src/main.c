#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "art_nouveau.h"
#include "alg_utils.h"
#include "config.h"
#include "dataset.h"
#include "dimer.h"
#include "input.h"
#include "kappa_dimer.h"
#include "my_mpi.h"
#include "sps_utils.h"
#include "target.h"
#ifdef LMP
#include "lmp_calculator.h"
#endif
#ifdef VASP
#include "vasp_calculator.h"
#endif


int main(int argc, char *argv[])
{
    int i, j, atom_index, errno, rank, size;

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

    /* write_input */
    if (rank == 0) {
        write_input(input);
    }

    /* random number */
    srand(input->random_seed + rank);
    if (rank == 0) {
        char filename[64];
        sprintf(filename, "%s/REPORT", input->output_dir);
        FILE *fp = fopen(filename, "w"); 
        fprintf(fp, "Random seed: %d\n", input->random_seed);
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

    /* read target */
    int target_num = 0;
    int list_size = 64;
    int *target_list = (int *)malloc(sizeof(int) * list_size);
    errno = read_target(config, input, &target_num, &target_list, &list_size);
    if (errno > 0) {
        printf("ERROR in TARGET FILE!\n");
        free_input(input);
        free_config(config);
        MPI_Finalize();
        return 1;
    }
    if (rank == 0) {
        write_target(input, target_num, target_list);
    }

    /* initial relax */
    if (input->init_relax > 0) {
        double energy;
        atom_relax(config, input, &energy, MPI_COMM_WORLD);
    }

    /* log */
    if (rank == 0) {
        FILE *fp;
        char filename[128];
        sprintf(filename, "%s/Statistics.log", input->output_dir);
        fp = fopen(filename, "w");
        fputs("----------------------------------------------------------------\n", fp);
        fputs(" Recycled events   New unique events   Relevant events   Trials\n", fp);    
        fputs("----------------------------------------------------------------\n", fp);
        fclose(fp);
        sprintf(filename, "%s/Event.log", input->output_dir);
        fp = fopen(filename, "w");
        fputs("----------------------------------------------------------------\n", fp);
        fputs(" Reaction index   Reaction barrier   Reaction rate   Degeneracy\n", fp);
        fputs("----------------------------------------------------------------\n", fp);
        fclose(fp);
        sprintf(filename, "%s/POSCAR", input->output_dir);
        write_config(config, filename, "w");
    }
    MPI_Barrier(MPI_COMM_WORLD);

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
    /* reset ncore for the last group */
    if (rank >= group_size * input->ncore) {
        input->ncore = size - group_size * input->ncore;
    }

    /* one-sided communication */
    MPI_Win count_win;
    int *global_count;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_count, &count_win);

    MPI_Win redundant_win;
    int *global_redundant;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_redundant, &redundant_win);

    MPI_Win recycle_win;
    int *global_recycle;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_recycle, &recycle_win);

    MPI_Win conv_win;
    int *global_conv;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_conv, &conv_win);

    MPI_Win exit_win;
    int *global_exit;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_exit, &exit_win);

    #ifdef LMP
    *global_count = 0;
    *global_recycle = 0;
    *global_conv = 0;
    *global_redundant = 0;
    *global_exit = 0;
    #endif

    int zero = 0;
    int one = 1;
    int mone = -1;
    int local_count = -1;
    int local_recycle = 0;
    int local_conv = 0;
    int local_redundant = 0;
    int local_exit = 0;
    int local_reac_num = 0;
    int local_dege_num = 0;
    double local_rate_sum = 0.0;
    int *local_reac_list = (int *)malloc(sizeof(int) * list_size);
    int *local_dege_list = (int *)malloc(sizeof(int) * list_size);
    double *local_acti_list = (double *)malloc(sizeof(double) * list_size);
    double *local_rate_list = (double *)malloc(sizeof(double) * list_size);
    /* all data in dataset are unique */
    /* do not have to count redundant search */
    Dataset *dataset = (Dataset *)malloc(sizeof(Dataset));
    dataset->head = NULL;
    /* dataset */
    Data *data;
    if (input->restart > 0) {
        build_dataset(dataset, input, config->tot_num);
    }

    int conv, unique;
    double Ea;
    double *eigenmode;
    while (1) {
        /* check exit condition */
        if (local_rank == 0) {
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, count_win);
            MPI_Fetch_and_op(&one, &local_count, MPI_INT,
                             0, (MPI_Aint)0, MPI_SUM, count_win);
            MPI_Win_unlock(0, count_win);
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
        MPI_Bcast(&local_redundant, 1, MPI_INT, 0, local_comm);

        /* atom_index */
        data = dataset->head;
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

        /* initial/saddle/final configuration */
        Config *initial = (Config *)malloc(sizeof(Config));
        copy_config(initial, config);
        Config *final = (Config *)malloc(sizeof(Config));
        copy_config(final, config);
        if (input->art_nouveau > 0) {
            conv = art_nouveau(initial, final, input, eigenmode,
                               local_count, atom_index, &Ea, local_comm);
        } else if (input->kappa_dimer > 0) {
            conv = kappa_dimer(initial, final, input, eigenmode,
                               local_count, atom_index, &Ea, local_comm);
        } else {
            conv = dimer(initial, final, input, eigenmode,
                         local_count, atom_index, &Ea, local_comm);
        }
        if (local_rank == 0) {
            /* conv == 0 -> success */
            if (conv == 0) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, conv_win);
                MPI_Fetch_and_op(&one, &local_conv, MPI_INT,
                                 0, (MPI_Aint)0, MPI_SUM, conv_win);
                MPI_Win_unlock(0, conv_win);
                char filename[128];
                sprintf(filename, "Final_%d_%d.POSCAR", local_count, atom_index);
                unique = check_unique(final, input, filename);
                /* unique == 1 -> unique */
                if (unique > 0) {
                    local_reac_list[local_reac_num] = local_count;
                    local_acti_list[local_reac_num] = Ea;
                    double kT = 8.61733034e-5 * input->temperature;
                    local_rate_list[local_reac_num] = input->frequency
                                                    * exp(- Ea / kT);
                    local_rate_sum += local_rate_list[local_reac_num];
                    local_reac_num++;
                    if (local_reac_num >= list_size) {
                        list_size = list_size << 1;
                        local_reac_list = (int *)realloc(local_reac_list,
                                                 sizeof(int) * list_size);
                        local_acti_list = (double *)realloc(local_acti_list,
                                                    sizeof(double) * list_size);
                        local_rate_list = (double *)realloc(local_rate_list,
                                                    sizeof(double) * list_size);
                    }
                    if (data != NULL) {
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, recycle_win);
                        MPI_Fetch_and_op(&one, &local_recycle, MPI_INT,
                                         0, (MPI_Aint)0, MPI_SUM, recycle_win);
                        MPI_Win_unlock(0, recycle_win);
                    }
                } else {
                    char old_filename[128];
                    sprintf(old_filename, "%s/Final_%d_%d.POSCAR",
                            input->output_dir, local_count, atom_index);
                    char new_filename[128];
                    sprintf(new_filename, "%s/%d_Final_%d_%d.POSCAR",
                            input->output_dir, -unique, local_count, atom_index);
                    rename(old_filename, new_filename);
                    if (data == NULL) {
                        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, redundant_win);
                        MPI_Fetch_and_op(&one, &local_redundant, MPI_INT,
                                         0, (MPI_Aint)0, MPI_SUM, redundant_win);
                        MPI_Win_unlock(0, redundant_win);
                    }
                    local_dege_list[local_dege_num] = -unique;
                    local_dege_num++;
                    if (local_dege_num >= list_size) {
                        list_size = list_size << 1;
                        local_dege_list = (int *)realloc(local_dege_list,
                                                 sizeof(int) * list_size);
                    }
                }
            }
            #ifdef VASP
            char old_filename[128];
            sprintf(old_filename, "%s/OUTCAR", input->output_dir);
            char new_filename[128];
            sprintf(new_filename, "%s/OUTCAR_%d", input->output_dir, local_count);
            rename(old_filename, new_filename);
            # endif
        }
        free_config(initial);
        free_config(final);
    }
    int recycle_num = 0;
    if (rank == 0) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, recycle_win);
        MPI_Fetch_and_op(&zero, &local_recycle, MPI_INT,
                         0, (MPI_Aint)0, MPI_SUM, recycle_win);
        MPI_Win_unlock(0, recycle_win);
        recycle_num = local_recycle;
    }
    /* overestimated local_count */
    if (local_rank == 0) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, count_win);
        MPI_Fetch_and_op(&mone, &local_count, MPI_INT,
                         0, (MPI_Aint)0, MPI_SUM, count_win);
        MPI_Win_unlock(0, count_win);
    }
    int total_reac_num;
    int total_dege_num;
    double total_rate_sum;
    if (local_rank == 0) {
        MPI_Allreduce(&local_rate_sum, &total_rate_sum, 1, MPI_DOUBLE,
                      MPI_SUM, group_comm);
        MPI_Allreduce(&local_reac_num, &total_reac_num, 1, MPI_INT,
                      MPI_SUM, group_comm);
        MPI_Allreduce(&local_dege_num, &total_dege_num, 1, MPI_INT,
                      MPI_SUM, group_comm);
    }
    MPI_Bcast(&total_reac_num, 1, MPI_INT, 0, local_comm);

    if (rank == 0) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, conv_win);
        MPI_Fetch_and_op(&zero, &local_conv, MPI_INT,
                         0, (MPI_Aint)0, MPI_SUM, conv_win);
        MPI_Win_unlock(0, conv_win);
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, count_win);
        MPI_Fetch_and_op(&zero, &local_count, MPI_INT,
                         0, (MPI_Aint)0, MPI_SUM, count_win);
        MPI_Win_unlock(0, count_win);
    }
    if (rank == 0) {
        char filename[128];
        sprintf(filename, "%s/Statistics.log", input->output_dir);
        FILE *fp = fopen(filename, "a");
        fprintf(fp, " %15d   %17d   %15d   %6d\n",
                recycle_num, total_reac_num - recycle_num,
                local_conv, local_count);
        fclose(fp);
    }

    int *global_reac_num = (int *)malloc(sizeof(int) * group_size);
    int *global_dege_num = (int *)malloc(sizeof(int) * group_size);
    int *global_reac_list = (int *)malloc(sizeof(int) * total_reac_num);
    int *global_dege_list = (int *)malloc(sizeof(int) * total_dege_num);
    double *global_acti_list = (double *)malloc(sizeof(double) * total_reac_num);
    double *global_rate_list = (double *)malloc(sizeof(double) * total_reac_num);
    if (local_rank == 0) {
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
        MPI_Gatherv(local_rate_list, local_reac_num, MPI_DOUBLE,
                    global_rate_list, global_reac_num, disp, MPI_DOUBLE,
                    0, group_comm);
        MPI_Gatherv(local_reac_list, local_reac_num, MPI_INT,
                    global_reac_list, global_reac_num, disp, MPI_INT,
                    0, group_comm);

        /* degeneracy */
        MPI_Allgather(&local_dege_num, 1, MPI_INT,
                      global_dege_num, 1, MPI_INT, group_comm);
        disp[0] = 0;
        if (group_size > 1) {
            for (i = 1; i < group_size; ++i) {
                disp[i] = disp[i - 1] + global_dege_num[i - 1];
            }
        }
        MPI_Gatherv(local_dege_list, local_dege_num, MPI_INT,
                    global_dege_list, global_dege_num, disp, MPI_INT,
                    0, group_comm);
        free(disp);
    }

    /* log */
    if (rank == 0) {
        int *dege_num = (int *)calloc(total_reac_num, sizeof(int));
        for (i = 0; i < total_reac_num; ++i) {
            for (j = 0; j < total_dege_num; ++j) {
                if (global_dege_list[j] == global_reac_list[i]) {
                    dege_num[i]++;
                }
            }
        }
        char filename[128];
        sprintf(filename, "%s/Event.log", input->output_dir);
        FILE *fp;
        fp = fopen(filename, "a");
        for (i = 0; i < total_reac_num; ++i) {
            fprintf(fp, " %14d   %16f   %13e   %10d\n",
                    global_reac_list[i], global_acti_list[i],
                    global_rate_list[i], dege_num[i]);
        }
        fclose(fp);
        free(dege_num);
    }

    free(target_list);
    free(local_reac_list);
    free(local_acti_list);
    free(local_rate_list);
    free(local_dege_list);
    free(global_reac_num);
    free(global_reac_list);
    free(global_acti_list);
    free(global_rate_list);
    free(global_dege_list);

    free_dataset(dataset);
    free_config(config);
    free_input(input);

    MPI_Win_free(&count_win);
    MPI_Win_free(&redundant_win);
    MPI_Win_free(&recycle_win);
    MPI_Win_free(&conv_win);
    MPI_Win_free(&exit_win);

    MPI_Comm_free(&local_comm);
    MPI_Comm_free(&group_comm);
    MPI_Finalize();
    return 0;
}
