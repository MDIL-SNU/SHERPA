#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "my_mpi.h"
#include "target.h"


int read_target(Config *config, Input *input,
               int *target_num, int **target_list, int *list_size)
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
            if (tmp_target_num >= (*list_size)) {
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
        } else if (strncmp(ptr, "T", 1) == 0) {
            strtok(line, " \n\t");
            ptr = strtok(NULL, " \n\t");
            while (ptr != NULL) {
                int type = atoi(ptr);
                for (i = 0; i < config->tot_num; ++i) {
                    if (config->type[i] == type) {
                        (*target_list)[*target_num] = i;
                        (*target_num)++;
                        if ((*target_num) >= (*list_size)) {
                            (*list_size) = (*list_size) << 1;
                            *target_list = (int *)realloc(*target_list, sizeof(int) * (*list_size));
                        }
                    }
                }
                ptr = strtok(NULL, " \n\t");
            }
        } else if (strncmp(ptr, "A", 1) == 0) {
            while (config->tot_num >= (*list_size)) {
                (*list_size) = (*list_size) << 1;
            }
            *target_list = (int *)malloc(sizeof(int) * (*list_size));
            for (i = 0; i < config->tot_num; ++i) {
                (*target_list)[*target_num] = i;
                (*target_num)++;
            }
        } else {
            return 1;
        }
    }
    return 0;
}


void write_target(Input *input, int target_num, int *target_list)
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
