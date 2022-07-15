#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "config.h"
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
    int *target_index;
    errno = gen_target(config, input, &target_index, &target_num);

    // TODO: fix atom
    atom_relax(config, input, MPI_COMM_WORLD);

    /* main loop */
    //for (i = 0; i < target_num; ++i) {
    for (i = 0; i < 1; ++i) {
        ii = target_index[i];
        Config *tmp_config = (Config *)malloc(sizeof(Config));
        copy_config(tmp_config, config);
        for (j = 0; j < tmp_config->tot_num; ++j) {
            del[0] = tmp_config->pos[j * 3 + 0] - tmp_config->pos[ii * 3 + 0];
            del[1] = tmp_config->pos[j * 3 + 1] - tmp_config->pos[ii * 3 + 1];
            del[2] = tmp_config->pos[j * 3 + 2] - tmp_config->pos[ii * 3 + 2];
            get_minimum_image(del, tmp_config->boxlo, tmp_config->boxhi,
                              tmp_config->xy, tmp_config->yz, tmp_config->xz);
            double dist = norm(del);
            if (dist < input->cutoff) {
                tmp_config->pos[j * 3 + 0] += normal_random(0, input->stddev);
                tmp_config->pos[j * 3 + 1] += normal_random(0, input->stddev);
                tmp_config->pos[j * 3 + 2] += normal_random(0, input->stddev);
            }
        }
        free_config(tmp_config);
    }
    free(target_index);
    free_config(config);
    free_input(input);
    
    return 0;
}
