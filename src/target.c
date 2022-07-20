#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "target.h"


int gen_target(Config *config, Input *input, int **target_list, int *target_num)
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
            while (ptr != NULL) {
                (*target_num)++;
                ptr = strtok(NULL, " \n\t");
            }
            *target_list = (int *)malloc(sizeof(int) * (*target_num));
            strtok(tmp_line, " \n\t");
            for (i = 0; i < *target_num; ++i) {
                (*target_list)[i] = atoi(strtok(NULL, " \n\t"));
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


