#include "target.h"
#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int read_target(Config *config, int *target_num, int **target_list,
                int *list_size)
{
    int i;
    FILE *fp;
    fp = fopen("./TARGET", "r");
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
        } else {
            int random = 0;
            if (strchr(ptr, 'R') != NULL ) {
                random = 1;
            }
            if (strncmp(ptr, "I", 1) == 0) {
                strcpy(tmp_line, line);
                strtok(line, " \n\t");
                ptr = strtok(NULL, " \n\t");
                int tmp_target_num = 0;
                while (ptr != NULL) {
                    tmp_target_num++;
                    ptr = strtok(NULL, " \n\t");
                }
                while (tmp_target_num + (*target_num) > (*list_size)) {
                    (*list_size) = (*list_size) << 1;
                }
                *target_list = (int *)realloc(*target_list, sizeof(int) * (*list_size));
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
                    while (config->each_num[type - 1] + (*target_num) >= (*list_size)) {
                        (*list_size) = (*list_size) << 1;
                    }
                    int begin = 0;
                    if (type > 1) {
                        for (i = 0; i < type - 1; ++i) {
                            begin += config->each_num[i];
                        }
                    }
                    for (i = begin; i < begin + config->each_num[type - 1]; ++i) {
                        (*target_list)[*target_num] = i;
                        (*target_num)++;
                    }
                    ptr = strtok(NULL, " \n\t");
                }
            } else if (strncmp(ptr, "A", 1) == 0) {
                while (config->tot_num + (*target_num) >= (*list_size)) {
                    (*list_size) = (*list_size) << 1;
                }
                *target_list = (int *)realloc(*target_list, sizeof(int) * (*list_size));
                for (i = 0; i < config->tot_num; ++i) {
                    (*target_list)[*target_num] = i;
                    (*target_num)++;
                }
            } else {
                fclose(fp);
                return 1;
            }
            if (random > 0) {
                int_shuffle(*target_num, *target_list);
            }
        }
    }
    fclose(fp);
    return 0;
}


void write_target(int target_num, int *target_list)
{
    int i;
    FILE *fp = fopen("./TARGET_read", "w");
    fputs("I", fp);
    for (i = 0; i < target_num; ++i) {
        fprintf(fp, " %d", target_list[i]);
    }
    fclose(fp);
}
