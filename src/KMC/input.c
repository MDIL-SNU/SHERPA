#include "input.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


int input_int(int *var, char *tag, char *filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    char line[128], *ptr;
    while (1) {
        ptr = fgets(line, 128, fp);
        if (ptr == NULL) {
            break;
        } else if (strcmp(ptr, "\n") == 0 || strncmp(ptr, "#", 1) == 0) {
            continue;
        } else if (strncmp(ptr, tag, strlen(tag)) == 0) {
            strtok(line, " \n\t");
            strtok(NULL, " \n\t");
            *var = atoi(strtok(NULL, "\n"));
            fclose(fp);
            return 0;
        } else {
            continue;
        }
    }
    fclose(fp);
    return 1;
}


int input_double(double *var, char *tag, char *filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    char line[128], *ptr;
    while (1) {
        ptr = fgets(line, 128, fp);
        if (ptr == NULL) {
            break;
        } else if (strcmp(ptr, "\n") == 0 || strncmp(ptr, "#", 1) == 0) {
            continue;
        } else if (strncmp(ptr, tag, strlen(tag)) == 0) {
            strtok(line, " \n\t");
            strtok(NULL, " \n\t");
            *var = atof(strtok(NULL, "\n"));
            fclose(fp);
            return 0;
        } else {
            continue;
        }
    }
    fclose(fp);
    return 1;
}


int input_char(char **var, char *tag, char *filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    char line[128], *ptr;
    while (1) {
        ptr = fgets(line, 128, fp);
        if (ptr == NULL) {
            break;
        } else if (strcmp(ptr, "\n") == 0 || strncmp(ptr, "#", 1) == 0) {
            continue;
        } else if (strncmp(ptr, tag, strlen(tag)) == 0) {
            strtok(line, " \n\t");
            strtok(NULL, " \n\t");
            *var = (char *)malloc(sizeof(char) * 65536);
            strcpy(*var, strtok(NULL, "\n"));
            fclose(fp);
            return 0;
        } else {
            continue;
        }
    }
    fclose(fp);
    return 1;
}


int input_char_arr(char ***var, char *tag, int n, char *filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    char line[128], *ptr;
    while (1) {
        ptr = fgets(line, 128, fp);
        if (ptr == NULL) {
            break;
        } else if (strcmp(ptr, "\n") == 0 || strncmp(ptr, "#", 1) == 0) {
            continue;
        } else if (strncmp(ptr, tag, strlen(tag)) == 0) {
            strtok(line, " \n\t");
            strtok(NULL, " \n\t");
            *var = (char **)malloc(sizeof(char *) * n);
            for (int i = 0; i < n; ++i) {
                ptr = strtok(NULL, " \n\t");
                (*var)[i] = (char *)calloc(16, sizeof(char));
                strcpy((*var)[i], ptr);
            }
            fclose(fp);
            return 0;
        } else {
            continue;
        }
    }
    fclose(fp);
    return 1;
}


int read_input(Input *input, char *filename)
{
    int errno;
    errno = input_int(&(input->kmc_step), "KMC_STEP", filename);
    if (errno) {
        input->kmc_step = 1000;
    }
    errno = input_double(&(input->temperature), "TEMPERATURE", filename);
    if (errno) {
        input->temperature = 298;
    }
    errno = input_double(&(input->att_freq), "ATT_FREQ", filename);
    if (errno) {
        input->att_freq = 1e13;
    }
    errno = input_char(&(input->extractor), "EXTRACTOR", filename);
    if (errno) {
        printf("EXTRACTOR is necessary.\n");
        return 1;
    }
    errno = input_char(&(input->sherpa_path), "SHERPA_CMD", filename);
    if (errno) {
        printf("SHERPA_CMD is necessary.\n");
        return 1;
    }
    errno = input_int(&(input->num_inputs), "NUM_INPUTS", filename);
    if (errno) {
        printf("NUM_INPUTS is necessary.\n");
        return 1;
    }
    errno = input_char_arr(&(input->inputs), "INPUTS", input->num_inputs, filename);
    if (errno) {
        printf("INPUTS is necessary.\n");
        return 1;
    }
    errno = input_int(&(input->random_seed), "RANDOM_SEED", filename);
    if (errno) {
        input->random_seed = (unsigned int)time(NULL);
    }
    return 0;
}


void write_input(Input *input)
{
    int i;
    FILE *fp = fopen("./INPUT_KMC_read", "w");

    fputs("# general parameter #\n", fp);
    fprintf(fp, "KMC_STEP\t= %d\n", input->kmc_step);
    fprintf(fp, "TEMPERATURE\t= %f\n", input->temperature);
    fprintf(fp, "ATT_FREQ\t= %f\n", input->att_freq);
    fputs("\n", fp);

    fputs("# path parameter #\n", fp);
    fprintf(fp, "POT_PATH\t= %s\n", input->pot_path);
    fprintf(fp, "EXT_PATH\t= %s\n", input->ext_path);
    fprintf(fp, "SHERPA_CMD\t= %s\n", input->sherpa_cmd);
    fputs("\n", fp);

    fputs("# random parameter #\n", fp);
    fprintf(fp, "RANDOM_SEED\t= %d\n", input->random_seed);
    fputs("\n", fp);

    fclose(fp);
}


void free_input(Input *input)
{
    int i;
    for (i = 0; i < input->num_inputs; ++i) {
        free(input->inputs[i]);
    }
    free(input->inputs);
    free(input->extractor);
    free(input->sherpa_cmd);
    free(input);
}
