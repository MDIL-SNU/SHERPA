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
    errno = input_char(&(input->target_list), "TARGET_LIST", filename);
    if (errno) {
        printf("TARGET_LIST is missing.\n");
        return 1;
    }
    errno = input_double(&(input->finite_diff), "FINITE_DIFF", filename);
    if (errno) {
        printf("FINITE_DIFF is missing.\n");
        return 1;
    }
    errno = input_double(&(input->acti_cutoff), "ACTI_CUTOFF", filename);
    if (errno) {
        printf("ACTI_CUTOFF is missing.\n");
        return 1;
    }
    errno = input_int(&(input->acti_nevery), "ACTI_NEVERY", filename);
    if (errno) {
        printf("ACTI_NEVERY is missing.\n");
        return 1;
    }
    errno = input_double(&(input->f_tol), "F_TOL", filename);
    if (errno) {
        printf("F_TOL is missing.\n");
        return 1;
    }
    errno = input_double(&(input->diff_tol), "DIFF_TOL", filename);
    if (errno) {
        printf("DIFF_TOL is missing.\n");
        return 1;
    }
    errno = input_double(&(input->max_move), "MAX_MOVE", filename);
    if (errno) {
        printf("MAX_MOVE is missing.\n");
        return 1;
    }
    errno = input_double(&(input->trial_move), "TRIAL_MOVE", filename);
    if (errno) {
        printf("TRIAL_MOVE is missing.\n");
        return 1;
    }
    errno = input_double(&(input->confidence), "CONFIDENCE", filename);
    if (errno) {
        printf("CONFIDENCE is missing.\n");
        return 1;
    }
    errno = input_int(&(input->max_search), "MAX_SEARCH", filename);
    if (errno) {
        printf("MAX_SEARCH is missing.\n");
        return 1;
    }
    errno = input_int(&(input->nelem), "NELEMENT", filename);
    if (errno) {
        printf("NELEMENT is missing.\n");
        return 1;
    }
    errno = input_char_arr(&(input->atom_type), "ATOM_TYPE", input->nelem, filename);
    if (errno) {
        printf("ATOM_TYPE is missing.\n");
        return 1;
    }
    errno = input_char(&(input->init_config), "INIT_CONFIG", filename);
    if (errno) {
        printf("INIT_CONFIG is missing.\n");
        return 1;
    }
    errno = input_int(&(input->init_relax), "INIT_RELAX", filename);
    if (errno) {
        printf("INIT_RELAX is missing.\n");
        return 1;
    }
    errno = input_int(&(input->init_disp), "INIT_DISP", filename);
    if (errno) {
        printf("INIT_DISP is missing.\n");
        return 1;
    }
    errno = input_double(&(input->disp_cutoff), "DISP_CUTOFF", filename);
    if (errno) {
        printf("DISP_CUTOFF is missing.\n");
        return 1;
    }
    errno = input_double(&(input->disp_stddev), "DISP_STDDEV", filename);
    if (errno) {
        printf("DISP_STDDEV is missing.\n");
        return 1;
    }
    errno = input_int(&(input->init_mode), "INIT_MODE", filename);
    if (errno) {
        printf("INIT_MODE is missing.\n");
        return 1;
    }
    errno = input_char(&(input->mode_list), "MODE_LIST", filename);
    if (errno) {
        printf("MODE_LIST is missing.\n");
        return 1;
    }
    errno = input_char(&(input->pair_style), "PAIR_STYLE", filename);
    if (errno) {
        printf("PAIR_STYLE is missing.\n");
        return 1;
    }
    errno = input_char(&(input->pair_coeff), "PAIR_COEFF", filename);
    if (errno) {
        printf("PAIR_COEFF is missing.\n");
        return 1;
    }
    errno = input_int(&(input->ncore), "NCORE", filename);
    if (errno) {
        printf("NCORE is missing.\n");
        return 1;
    }
    errno = input_int(&(input->kappa_dimer), "KAPPA_DIMER", filename);
    if (errno) {
        printf("KAPPA_DIMER is missing.\n");
        return 1;
    }
    errno = input_double(&(input->f_rot_min), "F_ROT_MIN", filename);
    if (errno) {
        printf("F_ROT_MIN is missing.\n");
        return 1;
    }
    errno = input_double(&(input->f_rot_max), "F_ROT_MAX", filename);
    if (errno) {
        printf("F_ROT_MAX is missing.\n");
        return 1;
    }
    errno = input_int(&(input->max_num_rot), "MAX_NUM_ROT", filename);
    if (errno) {
        printf("MAX_NUM_ROT is missing.\n");
        return 1;
    }
    errno = input_int(&(input->max_num_tls), "MAX_NUM_TLS", filename);
    if (errno) {
        printf("MAX_NUM_TLS is missing.\n");
        return 1;
    }
    errno = input_int(&(input->art_nouveau), "ART_NOUVEAU", filename);
    if (errno) {
        printf("ART_NOUVEAU is missing.\n");
        return 1;
    }
    errno = input_double(&(input->lambda_conv), "LAMBDA_CONV", filename);
    if (errno) {
        printf("LAMBDA_CONV is missing.\n");
        return 1;
    }
    errno = input_int(&(input->max_num_rlx), "MAX_NUM_RLX", filename);
    if (errno) {
        printf("MAX_NUM_RLX is missing.\n");
        return 1;
    }
    errno = input_int(&(input->delay_step), "DELAY_STEP", filename);
    if (errno) {
        printf("DELAY_STEP is missing.\n");
        return 1;
    }
    errno = input_int(&(input->mixing_step), "MIXING_STEP", filename);
    if (errno) {
        printf("MIXING_STEP is missing.\n");
        return 1;
    }
    errno = input_int(&(input->hyper_step), "HYPER_STEP", filename);
    if (errno) {
        printf("HYPER_STEP is missing.\n");
        return 1;
    }
    errno = input_int(&(input->random_seed), "RANDOM_SEED", filename);
    if (errno) {
        printf("RANDOM_SEED is missing.\n");
        return 1;
    }
    if (input->random_seed == -1) {
        input->random_seed = (unsigned int)time(NULL);
    }
    input->nredundant = (int)round(1 / (1 - input->confidence));
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
    fprintf(fp, "TARGET_LIST\t= %s\n", input->target_list);
    fprintf(fp, "FINITE_DIFF\t= %f\n", input->finite_diff);
    fprintf(fp, "ACTI_CUTOFF\t= %f\n", input->acti_cutoff);
    fprintf(fp, "ACTI_NEVERY\t= %d\n", input->acti_nevery);
    fprintf(fp, "F_TOL\t\t= %f\n", input->f_tol);
    fprintf(fp, "DIFF_TOL\t= %f\n", input->diff_tol);
    fprintf(fp, "MAX_MOVE\t= %f\n", input->max_move);
    fprintf(fp, "TRIAL_MOVE\t= %f\n", input->trial_move);
    fprintf(fp, "CONFIDENCE\t= %f\n", input->confidence);
    fprintf(fp, "MAX_SEARCH\t= %d\n", input->max_search);
    fputs("\n", fp);

    fputs("# initial structure parameter #\n", fp);
    fprintf(fp, "NELEMENT\t= %d\n", input->nelem);
    fputs("ATOM_TYPE\t=", fp);
    for (i = 0; i < input->nelem; ++i) {
        fprintf(fp, " %s", input->atom_type[i]);
    }
    fputs("\n", fp);
    fprintf(fp, "INIT_CONFIG\t= %s\n", input->init_config);
    fprintf(fp, "INIT_RELAX\t= %d\n", input->init_relax);
    fprintf(fp, "INIT_DISP\t= %d\n", input->init_disp);
    fprintf(fp, "DISP_CUTOFF\t= %f\n", input->disp_cutoff);
    fprintf(fp, "DISP_STDDEV\t= %f\n", input->disp_stddev);
    fprintf(fp, "INIT_MODE\t= %d\n", input->init_mode);
    fprintf(fp, "MODE_LIST\t= %s\n", input->mode_list);
    fputs("\n", fp);

    fputs("# LAMMPS parameter #\n", fp);
    fprintf(fp, "PAIR_STYLE\t= %s\n", input->pair_style);
    fprintf(fp, "PAIR_COEFF\t= %s\n", input->pair_coeff);
    fprintf(fp, "NCORE\t\t= %d\n", input->ncore);
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
    free(input->init_config);
    free(input->target_list);
    free(input->pair_style);
    free(input->pair_coeff);
    free(input->mode_list);
    free(input);
}
