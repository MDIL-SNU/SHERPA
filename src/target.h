#ifndef __TARGET_H__
#define __TARGET_H__
#include "config.h"

int read_target(Config *config, char *filename,
                int *target_num, int **target_list, int *list_size);
void write_target(int target_num, int *target_list);
#endif
