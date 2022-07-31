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
        sprintf(filename, "%s/seed.dat", input->output_dir);
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
        atom_relax(config, input);
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
        sprintf(filename, "%s/kMC.XDATCAR", input->output_dir);
        write_config(config, filename);
    }
    /* kmc loop */
    for (step = 0; step < input->end_step; ++step) {
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
        int reac_num = 0;
        double rate_sum = 0.0;
        int *reac_list = (int *)malloc(sizeof(int) * target_num);
        double *acti_list = (double *)malloc(sizeof(double) * target_num);
        double *rate_list = (double *)malloc(sizeof(double) * target_num);
        for (i = 0; i < target_num; ++i) {
            ii = target_list[i];
            Config *initial_config = (Config *)malloc(sizeof(Config));
            copy_config(initial_config, config);
            if (rank == 0) {
                char line[128], filename[128];
                sprintf(filename, "%s/Dimer_%d.log", input->output_dir, i);
                FILE *fp = fopen(filename, "a");
                fputs(" Opt step   Rot step   Potential energy   Curvature   Rot angle   Rot force\n", fp);
                fclose(fp);
            }
            Config *final_config = NULL;
            double Ea;
            int conv = dimer(initial_config, &final_config, input, i, ii, &Ea);
            if (conv == 0) {
                if (rank == 0) {
                    char line[128], filename[128];
                    sprintf(filename, "%s/Dimer_%d.log", input->output_dir, i);
                    FILE *fp = fopen(filename, "a");
                    sprintf(line, " Barrier energy: %f eV\n", Ea);
                    fputs(line, fp);
                    fclose(fp);
                }
                reac_list[reac_num] = ii;
                acti_list[reac_num] = Ea;
                double kT = 8.61733034e-5 * input->temperature;
                rate_list[reac_num] = input->att_freq * exp(- Ea / kT);
                rate_sum += rate_list[reac_num];
                reac_num++;
            }
            free_config(initial_config);
            if (final_config) { 
                free_config(final_config);
            }
        }
        /* selection */
        //TODO: unique saddle
        double random1 = (double)rand() / RAND_MAX;
        double random2 = (double)rand() / RAND_MAX;
        MPI_Bcast(&random1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&random2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double rate_acc = 0.0;
        int reac_index = 0;
        for (i = 0; i < reac_num; ++i) {
            rate_list[i] /= rate_sum;
            rate_acc += rate_list[i];
            if (rate_acc > random1) {
                reac_index = i;
                break;
            } 
        }
        time -= log(random2) / rate_sum;
        double energy = oneshot(config, input);
        /* log */
        if (rank == 0) {
            char line[128], filename[128];
            sprintf(filename, "%s/kMC.log", input->output_dir);
            FILE *fp = fopen(filename, "a");
            sprintf(line, " %8lld   %16f   %16f   %8e\n",
                    step, energy, acti_list[reac_index], time);
            fputs(line, fp);
            fclose(fp);
            sprintf(filename, "%s/kMC.XDATCAR", input->output_dir);
            write_config(config, filename);
        }
        free_config(config);
        free(target_list);
        if (time > input->end_time) {
            free(reac_list);
            free(acti_list);
            free(rate_list);
            break;
        }
        config = (Config *)malloc(sizeof(Config));
        char filename[128];
        sprintf(filename, "%s/Final_%d.POSCAR",
                input->output_dir, reac_list[reac_index]);
        errno = read_config(config, input, filename);
        free(reac_list);
        free(acti_list);
        free(rate_list);
        if (errno > 0) {
            printf("ERROR in INIT_CONFIG FILE!\n");
            free_input(input);
            free(config);
            MPI_Finalize();
            return 1;
        }
    }
    free_input(input);
    MPI_Finalize();
    return 0;
}
