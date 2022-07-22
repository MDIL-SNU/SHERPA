#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "target.h"


int get_target(Config *config, Input *input,
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
