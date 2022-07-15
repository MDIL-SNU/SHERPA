#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"


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
            *var = (char *)malloc(sizeof(char) * 32);
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
                (*var)[i] = (char *)malloc(sizeof(char) * 16);
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
    errno = input_int(&(input->nelem), "NELEMENT", filename);
    if (errno) {
        return 1;
    }
    errno = input_int(&(input->init_mode), "INIT_MODE", filename);
    if (errno) {
        return 1;
    }
    errno = input_int(&(input->random_seed), "RANDOM_SEED", filename);
    if (errno) {
        return 1;
    }
    errno = input_char_arr(&(input->atom_type), "ATOM_TYPE", input->nelem, filename);
    if (errno) {
        return 1;
    }
    errno = input_char(&(input->pair_style), "PAIR_STYLE", filename);
    if (errno) {
        return 1;
    }
    errno = input_char(&(input->pair_coeff), "PAIR_COEFF", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->cutoff), "CUTOFF", filename);
    if (errno) {
        return 1;
    }
    errno = input_char(&(input->init_config), "INIT_CONFIG", filename);
    if (errno) {
        return 1;
    }
    errno = input_char(&(input->target), "TARGET", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->dimer_dist), "DIMER_DIST", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->fmax), "FMAX", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->ftol), "FTOL", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->trial_angle), "TRIAL_ANGLE", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->disp_cutoff), "DISP_CUTOFF", filename);
    if (errno) {
        return 1;
    }
    errno = input_double(&(input->stddev), "STDDEV", filename);
    if (errno) {
        return 1;
    }
    return 0;
}


void free_input(Input *input)
{
    for (int i = 0; i < input->nelem; ++i) {
        free(input->atom_type[i]);
    }
    free(input->atom_type);
    free(input->pair_style);
    free(input->pair_coeff);
    free(input->init_config);
    free(input);
}
