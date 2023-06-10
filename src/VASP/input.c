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
    errno = input_double(&(input->acti_cutoff), "ACTI_CUTOFF", filename);
    if (errno) {
        input->acti_cutoff = 6.0;
    }
    errno = input_int(&(input->acti_nevery), "ACTI_NEVERY", filename);
    if (errno) {
        input->acti_nevery = 3;
    }
    errno = input_double(&(input->finite_diff), "FINITE_DIFF", filename);
    if (errno) {
        input->finite_diff = 0.01;
    }
    errno = input_double(&(input->f_tol), "F_TOL", filename);
    if (errno) {
        input->f_tol = 0.01;
    }
    errno = input_double(&(input->diff_tol), "DIFF_TOL", filename);
    if (errno) {
        input->diff_tol = 0.4;
    }
    errno = input_double(&(input->max_move), "MAX_MOVE", filename);
    if (errno) {
        input->max_move = 0.1;
    }
    errno = input_double(&(input->trial_move), "TRIAL_MOVE", filename);
    if (errno) {
        input->trial_move = 0.01;
    }
    errno = input_double(&(input->confidence), "CONFIDENCE", filename);
    if (errno) {
        input->confidence = 0.99;
    }
    errno = input_int(&(input->max_search), "MAX_SEARCH", filename);
    if (errno) {
        input->max_search = 100;
    }
    errno = input_int(&(input->cont), "CONTINUE", filename);
    if (errno) {
        input->cont = 0;
    }
    errno = input_int(&(input->nelem), "NELEMENT", filename);
    if (errno) {
        printf("NELEMENT is necessary.\n");
        return 1;
    }
    errno = input_char_arr(&(input->atom_type), "ATOM_TYPE", input->nelem, filename);
    if (errno) {
        printf("ATOM_TYPE is necessary.\n");
        return 1;
    }
    errno = input_int(&(input->init_relax), "INIT_RELAX", filename);
    if (errno) {
        input->init_relax = 1;
    }
    errno = input_int(&(input->init_disp), "INIT_DISP", filename);
    if (errno) {
        input->init_disp = 0;
    }
    errno = input_double(&(input->disp_cutoff), "DISP_CUTOFF", filename);
    if (errno) {
        input->disp_cutoff = 3.0;
    }
    errno = input_double(&(input->disp_move), "DISP_MOVE", filename);
    if (errno) {
        input->disp_move = 0.1;
    }
    errno = input_int(&(input->init_mode), "INIT_MODE", filename);
    if (errno) {
        input->init_mode = 0;
    }
    errno = input_char(&(input->vasp_cmd), "VASP_CMD", filename);
    if (errno) {
        printf("VASP_CMD is necessary.\n");
        return 1;
    }
    errno = input_int(&(input->kappa_dimer), "KAPPA_DIMER", filename);
    if (errno) {
        input->kappa_dimer = 0;
    }
    errno = input_double(&(input->f_rot_min), "F_ROT_MIN", filename);
    if (errno) {
        input->f_rot_min = 0.1;
    }
    errno = input_double(&(input->f_rot_max), "F_ROT_MAX", filename);
    if (errno) {
        input->f_rot_max = 1.0;
    }
    errno = input_int(&(input->max_num_rot), "MAX_NUM_ROT", filename);
    if (errno) {
        input->max_num_rot = 4;
    }
    errno = input_int(&(input->max_num_tls), "MAX_NUM_TLS", filename);
    if (errno) {
        input->max_num_tls = 500;
    }
    errno = input_int(&(input->art_nouveau), "ART_NOUVEAU", filename);
    if (errno) {
        input->art_nouveau = 1;
    }
    errno = input_double(&(input->lambda_conv), "LAMBDA_CONV", filename);
    if (errno) {
        input->lambda_conv = 0.01;
    }
    errno = input_int(&(input->max_num_rlx), "MAX_NUM_RLX", filename);
    if (errno) {
        input->max_num_rlx = 4;
    }
    errno = input_int(&(input->delay_step), "DELAY_STEP", filename);
    if (errno) {
        input->delay_step = 0;
    }
    errno = input_int(&(input->mixing_step), "MIXING_STEP", filename);
    if (errno) {
        input->mixing_step = 0;
    }
    errno = input_int(&(input->hyper_step), "HYPER_STEP", filename);
    if (errno) {
        input->hyper_step = 3;
    }
    errno = input_int(&(input->random_seed), "RANDOM_SEED", filename);
    if (errno) {
        input->random_seed = (unsigned int)time(NULL);
    }
    input->nredundant = (int)round(1 / (1 - input->confidence));
    input->ncore = 1;
    if (input->kappa_dimer + input->art_nouveau > 1) {
        printf("Choose only one algirhtm.\n");
        return 1;
    }
    return 0;
}


void write_input(Input *input)
{
    int i;
    FILE *fp = fopen("./INPUT_read", "w");

    fputs("# general parameter #\n", fp);
    fprintf(fp, "ACTI_CUTOFF\t= %f\n", input->acti_cutoff);
    fprintf(fp, "ACTI_NEVERY\t= %d\n", input->acti_nevery);
    fprintf(fp, "FINITE_DIFF\t= %f\n", input->finite_diff);
    fprintf(fp, "F_TOL\t\t= %f\n", input->f_tol);
    fprintf(fp, "DIFF_TOL\t= %f\n", input->diff_tol);
    fprintf(fp, "MAX_MOVE\t= %f\n", input->max_move);
    fprintf(fp, "TRIAL_MOVE\t= %f\n", input->trial_move);
    fprintf(fp, "CONFIDENCE\t= %f\n", input->confidence);
    fprintf(fp, "MAX_SEARCH\t= %d\n", input->max_search);
    fprintf(fp, "CONTINUE\t= %d\n", input->cont);
    fputs("\n", fp);

    fputs("# initial structure parameter #\n", fp);
    fprintf(fp, "NELEMENT\t= %d\n", input->nelem);
    fputs("ATOM_TYPE\t=", fp);
    for (i = 0; i < input->nelem; ++i) {
        fprintf(fp, " %s", input->atom_type[i]);
    }
    fputs("\n", fp);
    fprintf(fp, "INIT_RELAX\t= %d\n", input->init_relax);
    fprintf(fp, "INIT_DISP\t= %d\n", input->init_disp);
    fprintf(fp, "DISP_CUTOFF\t= %f\n", input->disp_cutoff);
    fprintf(fp, "DISP_MOVE\t= %f\n", input->disp_move);
    fprintf(fp, "INIT_MODE\t= %d\n", input->init_mode);
    fputs("\n", fp);

    fputs("# VASP parameter #\n", fp);
    fprintf(fp, "VASP_CMD\t= %s\n", input->vasp_cmd);
    fputs("\n", fp);

    fputs("# dimer parameter #\n", fp);
    fprintf(fp, "KAPPA_DIMER\t= %d\n", input->kappa_dimer);
    fprintf(fp, "F_ROT_MIN\t= %f\n", input->f_rot_min);
    fprintf(fp, "F_ROT_MAX\t= %f\n", input->f_rot_max);
    fprintf(fp, "MAX_NUM_ROT\t= %d\n", input->max_num_rot);
    fprintf(fp, "MAX_NUM_TLS\t= %d\n", input->max_num_tls);
    fputs("\n", fp);

    fputs("# art_nouveau parameter #\n", fp);
    fprintf(fp, "ART_NOUVEAU\t= %d\n", input->art_nouveau);
    fprintf(fp, "LAMBDA_CONV\t= %f\n", input->lambda_conv);
    fprintf(fp, "MAX_NUM_RLX\t= %d\n", input->max_num_rlx);
    fprintf(fp, "DELAY_STEP\t= %d\n", input->delay_step);
    fprintf(fp, "MIXING_STEP\t= %d\n", input->mixing_step);
    fprintf(fp, "HYPER_STEP\t= %d\n", input->hyper_step);
    fputs("\n", fp);

    fputs("# random parameter #\n", fp);
    fprintf(fp, "RANDOM_SEED\t= %d\n", input->random_seed);
    fputs("\n", fp);

    fclose(fp);
}


void free_input(Input *input)
{
    for (int i = 0; i < input->nelem; ++i) {
        free(input->atom_type[i]);
    }
    free(input->atom_type);
    free(input->vasp_cmd);
    free(input);
}
