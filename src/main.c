#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "calculator.h"
#include "config.h"
#include "dimer.h"
#include "input.h"
#include "target.h"
#include "utils.h"


int main(int argc, char *argv[])
{
    int i, ii, j, errno, rank, size;
    int buffer_size = 4096;
    double del[3];

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

    MPI_Win win;
    int *global_index;
    MPI_Win_allocate((MPI_Aint)sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &global_index, &win);

    /* kmc loop */
    for (step = 1; step <= input->end_step; ++step) {
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
        int local_index;
        *global_index = 0;
        int count = 1;
        int local_reac_num = 0;
        int *global_reac_num = (int *)malloc(sizeof(int) * group_size);
        double local_rate_sum = 0.0;
        double global_rate_sum;
        int *local_reac_list = (int *)malloc(sizeof(int) * target_num);
        double *local_acti_list = (double *)malloc(sizeof(double) * target_num);
        double *local_rate_list = (double *)malloc(sizeof(double) * target_num);
        int *global_reac_list = (int *)malloc(sizeof(int) * target_num);
        double *global_acti_list = (double *)malloc(sizeof(double) * target_num);
        double *global_rate_list = (double *)malloc(sizeof(double) * target_num);
        //TODO: recycle
        while (1) {
            if (local_rank == 0) {
                MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
                MPI_Fetch_and_op(&count, &local_index, MPI_INT, 0, (MPI_Aint)0,
                                 MPI_SUM, win);
                MPI_Win_unlock(0, win);
            }
            MPI_Bcast(&local_index, 1, MPI_INT, 0, local_comm);
            if (local_index >= target_num) {
                break;
            }
            Config *initial_config = (Config *)malloc(sizeof(Config));
            copy_config(initial_config, config);
            Config *final_config = NULL;
            double Ea;
            int conv = dimer(initial_config, &final_config, input,
                             local_index, target_list[local_index], &Ea, local_comm);
            if (conv == 0) {
                if (local_rank == 0) {
                    char line[128], filename[128];
                    sprintf(filename, "%s/Dimer_%d.log", input->output_dir, local_index);
                    FILE *fp = fopen(filename, "a");
                    sprintf(line, " Barrier energy: %f eV\n", Ea);
                    fputs(line, fp);
                    fclose(fp);
                }
                local_reac_list[local_reac_num] = local_index;
                local_acti_list[local_reac_num] = Ea;
                double kT = 8.61733034e-5 * input->temperature;
                local_rate_list[local_reac_num] = input->att_freq * exp(- Ea / kT);
                local_rate_sum += local_rate_list[local_reac_num];
                local_reac_num++;
            }
            free_config(initial_config);
            if (final_config) { 
                free_config(final_config);
            }
        }

        int total_reac_num;
        if (local_rank == 0) {
            MPI_Allgather(&local_reac_num, 1, MPI_INT,
                          global_reac_num, 1, MPI_INT, group_comm);
            int *disp = (int *)calloc(group_size, sizeof(int));
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
            MPI_Reduce(&local_rate_sum, &global_rate_sum, 1, MPI_DOUBLE,
                       MPI_SUM, 0, group_comm);
            MPI_Reduce(&local_reac_num, &total_reac_num, 1, MPI_INT,
                       MPI_SUM, 0, group_comm);
            free(disp);
        }

        /* selection */
        //TODO: unique saddle
        int reac_index;
        if (rank == 0) {
            double random1 = (double)rand() / RAND_MAX;
            double random2 = (double)rand() / RAND_MAX;
            double rate_acc = 0.0;
            for (i = 0; i < total_reac_num; ++i) {
                rate_acc += global_rate_list[i] / global_rate_sum;
                if (rate_acc > random1) {
                    reac_index = i;
                    break;
                } 
            }
            time -= log(random2) / global_rate_sum;
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
            write_config(config, filename, 1);
            printf("--------------------------------------------------------------\n");
            for (i = 0; i < total_reac_num; ++i) {
                printf(" %8lld   %14d   %16f   %13e\n",
                       step, global_reac_list[i],
                       global_acti_list[i], global_rate_list[i]);
            }
        }
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
        config = (Config *)malloc(sizeof(Config));
        errno = read_config(config, input, filename);
        MPI_Barrier(MPI_COMM_WORLD);
        free(local_reac_list);
        free(local_acti_list);
        free(local_rate_list);
        free(global_reac_num);
        free(global_reac_list);
        free(global_acti_list);
        free(global_rate_list);
        if (errno > 0) {
            printf("ERROR in INIT_CONFIG FILE!\n");
            free_input(input);
            free(config);
            MPI_Finalize();
            return 1;
        }
    }
    free_input(input);
    MPI_Win_free(&win);
    MPI_Comm_free(&local_comm);
    MPI_Comm_free(&group_comm);
    MPI_Finalize();
    return 0;
}
