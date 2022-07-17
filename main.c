#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
        exit(1);
    }
    srand(input->random_seed);

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

    // TODO: fix atom
    atom_relax(config, input);

    /* main loop */
    //for (i = 0; i < target_num; ++i) {
    for (i = 0; i < 1; ++i) {
        ii = target_list[i];
        Config *tmp_config = (Config *)malloc(sizeof(Config));
        copy_config(tmp_config, config);
        dimer(tmp_config, input, ii);
        free_config(tmp_config);
    }
    free(target_list);
    free_config(config);
    free_input(input);
    
    MPI_Finalize();
    return 0;
}
