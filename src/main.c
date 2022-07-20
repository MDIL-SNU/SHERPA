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

    /* make directory */
    if (rank == 0) {
        errno = mkdir("./output", 0775);
    }
    MPI_Bcast(&errno, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (errno != 0) {
        if (rank == 0) {
            printf("'output' directory exists!\n");
            printf("We will not overwrite your data.\n");
        }
        MPI_Finalize();
        return 1;
    }

    /* read input */
    Input *input = (Input *)malloc(sizeof(Input));
    errno = read_input(input, "./INPUT");
    if (errno > 0) {
        printf("ERROR in INPUT FILE!\n");
        free(input);
        exit(1);
    }
    srand(input->random_seed);
    if (rank == 0) {
        char line[64], filename[64];
        sprintf(filename, "output/seed.dat");
        FILE *fp = fopen(filename, "w"); 
        sprintf(line, "Random seed: %d\n", input->random_seed);
        fputs(line, fp);
        fclose(fp);
    }

    /* read config */
    Config *config = (Config *)malloc(sizeof(Config));
    errno = read_config(config, input, input->init_config);

    if (errno > 0) {
        printf("ERROR in CONFIG FILE!\n");
        free_input(input);
        exit(1);
    }
    
    /* read target */
    int target_num = 0;
    int *target_list;
    errno = gen_target(config, input, &target_list, &target_num);

    if (input->init_relax > 0) {
        atom_relax(config, input);
    }

    /* main loop */
    //for (i = 0; i < target_num; ++i) {
    for (i = 0; i < 1; ++i) {
        ii = target_list[i];
        Config *tmp_config = (Config *)malloc(sizeof(Config));
        copy_config(tmp_config, config);
        if (rank == 0) {
            char line[128], filename[128];
            sprintf(filename, "output/Dimer_%d.log", i);
            FILE *fp = fopen(filename, "a");
            fputs(" Opt step   Rot step   Potential energy   Curvature   Rot angle   Rot force\n", fp);
            fclose(fp);
        }
        double barrier_E = dimer(tmp_config, input, i, ii);
        if (rank == 0) {
            char line[128], filename[128];
            sprintf(filename, "output/Dimer_%d.log", i);
            FILE *fp = fopen(filename, "a");
            fputs("---------------------------------------------------------------------------\n", fp);
            sprintf(line, " Barrier energy: %f eV\n", barrier_E);
            fputs(line, fp);
            fclose(fp);
        }
        free_config(tmp_config);
    }
    free(target_list);
    free_config(config);
    free_input(input);
    
    MPI_Finalize();
    return 0;
}
