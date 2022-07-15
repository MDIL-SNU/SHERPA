#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "input.h"


int main(int argc, char *argv[])
{
    int errno;

    /* read input */
    Input *input = (Input *)malloc(sizeof(Input));
    errno = read_input(input, "./INPUT");
    if (errno) {
        printf("ERROR in INPUT FILE!\n");
        free(input);
        exit(1);
    }

    /* read config */
    FILE *fp;
    fp = fopen(input->init_config, "r");
    if (fp == NULL) {
        printf("NO CONFIG FILE!\n");
        free_input(input);
        exit(1);
    }
    free_input(input);
    
    return 0;
}
