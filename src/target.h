#ifndef __TARGET_H__
#define __TARGET_H__
#include "config.h"
#include "input.h"

int read_target(Config *config, Input *input,
                int *target_num, int **target_list, int *list_size);
void write_target(Input *input, int target_num, int *target_list);
#endif
