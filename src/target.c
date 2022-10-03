#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "my_mpi.h"
#include "target.h"
#include "utils.h"


int read_target(Config *config, Input *input,
               int **target_list, int *target_num, int *list_size)
{
    int i;
    FILE *fp;
    fp = fopen(input->target_list, "r");
    if (fp == NULL) {
        return 1;
    }
    char line[4096], tmp_line[4096], *ptr;
    while (1) {
        ptr = fgets(line, 4096, fp);
        if (ptr == NULL) {
            break;
        } else if (strcmp(ptr, "\n") == 0 || strncmp(ptr, "#", 1) == 0) {
            continue;
        } else if (strncmp(ptr, "I", 1) == 0) {
            strcpy(tmp_line, line);
            strtok(line, " \n\t");
            ptr = strtok(NULL, " \n\t");
            int tmp_target_num = 0;
            while (ptr != NULL) {
                tmp_target_num++;
                ptr = strtok(NULL, " \n\t");
            }
            if (tmp_target_num + (*target_num) > (*list_size)) {
                do {
                    (*list_size) = (*list_size) << 1;
                } while (tmp_target_num + (*target_num) > (*list_size));
                *target_list = (int *)realloc(*target_list, sizeof(int) * (*list_size));
            }
            strtok(tmp_line, " \n\t");
            for (i = 0; i < tmp_target_num; ++i) {
                (*target_list)[*target_num] = atoi(strtok(NULL, " \n\t"));
                (*target_num)++;
            }
        } else if (strncmp(ptr, "A", 1) == 0) {
            //gen_all();
        } else if (strncmp(ptr, "R", 1) == 0) {
            //gen_random();
        } else {
            return 1;
        }
    }
    return 0;
}


//TODO: A, R
void write_target(Input *input, int *target_list, int target_num)
{
    int i;
    char filename[128];
    sprintf(filename, "%s/TARGET", input->output_dir);
    FILE *fp = fopen(filename, "w");
    fputs("I", fp);
    for (i = 0; i < target_num; ++i) {
        fprintf(fp, " %d", target_list[i]);
    }
    fclose(fp);
}


int get_index(Config *config_new, Config *config_old)
{
    int i, rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = config_new->tot_num;
    int q = n / size;
    int r = n % size;
    int begin = rank * q + ((rank > r) ? r : rank);
    int end = begin + q;
    if (r > rank) {
        end++;
    }

    double del[3];
    double max_dist = 0.0;
    int max_index = 0;
    for (i = begin; i < end; ++i) {
        del[0] = config_new->pos[i * 3 + 0] - config_old->pos[i * 3 + 0]; 
        del[1] = config_new->pos[i * 3 + 1] - config_old->pos[i * 3 + 1]; 
        del[2] = config_new->pos[i * 3 + 2] - config_old->pos[i * 3 + 2]; 
        get_minimum_image(del, config_new->boxlo, config_new->boxhi,
                          config_new->xy, config_new->yz, config_new->xz);
        double dist = sqrt(del[0] * del[0]
                         + del[1] * del[1]
                         + del[2] * del[2]);
        if (dist > max_dist) {
            max_dist = dist;
            max_index = i;
        }
    }
    double *dist_list = (double *)malloc(sizeof(double) * size);
    MPI_Gather(&max_dist, 1, MPI_DOUBLE,
               dist_list, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int *index_list = (int *)malloc(sizeof(int) * size);
    MPI_Gather(&max_index, 1, MPI_INT,
               index_list, 1, MPI_INT, 0, MPI_COMM_WORLD);

    max_dist = 0.0;
    max_index = 0;
    for (i = 0; i < size; ++i) {
        if (dist_list[i] > max_dist) {
            max_dist = dist_list[i];
            max_index = index_list[i];
        }
    }
    free(dist_list);
    free(index_list);

    MPI_Bcast(&max_index, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return max_index;
}
