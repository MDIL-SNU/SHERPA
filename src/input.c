#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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
            *var = (char *)malloc(sizeof(char) * 128);
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
    errno = input_char(&(input->target_list), "TARGET_LIST", filename);
    if (errno) {
        printf("TARGET_LIST is missing.\n");
        return 1;
    }
    errno = input_double(&(input->disp_dist), "DISP_DIST", filename);
    if (errno) {
        printf("DISP_DIST is missing.\n");
        return 1;
    }
    errno = input_double(&(input->acti_cutoff), "ACTI_CUTOFF", filename);
    if (errno) {
        printf("ACTI_CUTOFF is missing.\n");
        return 1;
    }
    errno = input_double(&(input->f_tol), "F_TOL", filename);
    if (errno) {
        printf("F_TOL is missing.\n");
        return 1;
    }
    errno = input_double(&(input->stddev), "STDDEV", filename);
    if (errno) {
        printf("STDDEV is missing.\n");
        return 1;
    }
    errno = input_double(&(input->max_step), "MAX_STEP", filename);
    if (errno) {
        printf("MAX_STEP is missing.\n");
        return 1;
    }
    errno = input_double(&(input->trial_step), "TRIAL_STEP", filename);
    if (errno) {
        printf("TRIAL_STEP is missing.\n");
        return 1;
    }
    errno = input_int(&(input->init_relax), "INIT_RELAX", filename);
    if (errno) {
        printf("INIT_RELAX is missing.\n");
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
    errno = input_int(&(input->write_mode), "WRITE_MODE", filename);
    if (errno) {
        printf("WRITE_MODE is missing.\n");
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
    errno = input_double(&(input->pair_cutoff), "PAIR_CUTOFF", filename);
    if (errno) {
        printf("PAIR_CUTOFF is missing.\n");
        return 1;
    }
    errno = input_int(&(input->ncore), "NCORE", filename);
    if (errno) {
        printf("NCORE is missing.\n");
        return 1;
    }
    errno = input_char(&(input->vasp_cmd), "VASP_CMD", filename);
    if (errno) {
        printf("VASP_CMD is missing.\n");
        return 1;
    }
    errno = input_int(&(input->istart), "ISTART", filename);
    if (errno) {
        printf("ISTART is missing.\n");
        return 1;
    }
    errno = input_int(&(input->kappa_dimer), "KAPPA_DIMER", filename);
    if (errno) {
        printf("KAPPA_DIMER is missing.\n");
        return 1;
    }
    errno = input_int(&(input->snc_dimer), "SNC_DIMER", filename);
    if (errno) {
        printf("SNC_DIMER is missing.\n");
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
    errno = input_double(&(input->trial_angle), "TRIAL_ANGLE", filename);
    if (errno) {
        printf("TRIAL_ANGLE is missing.\n");
        return 1;
    }
    errno = input_int(&(input->art_nouveau), "ART_NOUVEAU", filename);
    if (errno) {
        printf("ART_NOUVEAU is missing.\n");
        return 1;
    }
    errno = input_double(&(input->lambda_crit), "LAMBDA_CRIT", filename);
    if (errno) {
        printf("LAMBDA_CRIT is missing.\n");
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
    errno = input_int(&(input->art_delay), "ART_DELAY", filename);
    if (errno) {
        printf("ART_DELAY is missing.\n");
        return 1;
    }
    errno = input_int(&(input->art_mixing), "ART_MIXING", filename);
    if (errno) {
        printf("ART_MIXING is missing.\n");
        return 1;
    }
    errno = input_double(&(input->frequency), "FREQUENCY", filename);
    if (errno) {
        printf("FREQUENCY is missing.\n");
        return 1;
    }
    errno = input_double(&(input->temperature), "TEMPERATURE", filename);
    if (errno) {
        printf("TEMPERATURE is missing.\n");
        return 1;
    }
    errno = input_int(&(input->random_seed), "RANDOM_SEED", filename);
    if (errno) {
        printf("RANDOM_SEED is missing.\n");
        return 1;
    }
    errno = input_char(&(input->output_dir), "OUTPUT_DIR", filename);
    if (errno) {
        printf("OUTPUT_DIR is missing.\n");
        return 1;
    }
    errno = input_int(&(input->restart), "RESTART", filename);
    if (errno) {
        printf("RESTART is missing.\n");
        return 1;
    }
    errno = input_char(&(input->restart_dir), "RESTART_DIR", filename);
    if (errno) {
        printf("RESTART_DIR is missing.\n");
        return 1;
    }
    if (input->random_seed == -1) {
        input->random_seed = (unsigned int)time(NULL);
    }
    #ifdef VASP
    input->ncore = 1;
    input->pair_cutoff = 100.0;
    #endif
    input->trial_angle *= 3.1415926535897932384626 / 180;
    input->nredundant = (int)round(1 / (1 - input->confidence));
    if (input->kappa_dimer + input->snc_dimer + input->art_nouveau > 1) {
        return 1;
    }
    return 0;
}


void write_input(Input *input)
{
    int i;
    char filename[128];
    sprintf(filename, "%s/INPUT", input->output_dir);
    FILE *fp = fopen(filename, "w");

    fputs("# general parameter #\n", fp);
    fprintf(fp, "NELEMENT\t= %d\n", input->nelem);
    fputs("ATOM_TYPE\t=", fp);
    for (i = 0; i < input->nelem; ++i) {
        fprintf(fp, " %s", input->atom_type[i]);
    }
    fputs("\n", fp);
    fprintf(fp, "INIT_CONFIG\t= %s\n", input->init_config);
    fprintf(fp, "TARGET_LIST\t= %s\n", input->target_list);
    fprintf(fp, "DISP_DIST\t= %f\n", input->disp_dist);
    fprintf(fp, "ACTI_CUTOFF\t= %f\n", input->acti_cutoff);
    fprintf(fp, "F_TOL\t\t= %f\n", input->f_tol);
    fprintf(fp, "STDDEV\t\t= %f\n", input->stddev);
    fprintf(fp, "MAX_STEP\t= %f\n", input->max_step);
    fprintf(fp, "TRIAL_STEP\t= %f\n", input->trial_step);
    fprintf(fp, "INIT_RELAX\t= %d\n", input->init_relax);
    fprintf(fp, "CONFIDENCE\t= %f\n", input->confidence);
    fprintf(fp, "MAX_SEARCH\t= %d\n", input->max_search);
    fprintf(fp, "WRITE_MODE\t= %d\n", input->write_mode);
    fputs("\n", fp);

    fputs("# LAMMPS parameter #\n", fp);
    fprintf(fp, "PAIR_STYLE\t= %s\n", input->pair_style);
    fprintf(fp, "PAIR_COEFF\t= %s\n", input->pair_coeff);
    fprintf(fp, "PAIR_CUTOFF\t= %f\n", input->pair_cutoff);
    fprintf(fp, "NCORE\t\t= %d\n", input->ncore);
    fputs("\n", fp);

    fputs("# VASP parameter #\n", fp);
    fprintf(fp, "VASP_CMD\t= %s\n", input->vasp_cmd);
    fprintf(fp, "ISTART\t\t= %d\n", input->istart);
    fputs("\n", fp);

    fputs("# dimer parameter #\n", fp);
    fprintf(fp, "KAPPA_DIMER\t= %d\n", input->kappa_dimer);
    fprintf(fp, "SNC_DIMER\t= %d\n", input->snc_dimer);
    fprintf(fp, "F_ROT_MIN\t= %f\n", input->f_rot_min);
    fprintf(fp, "F_ROT_MAX\t= %f\n", input->f_rot_max);
    fprintf(fp, "MAX_NUM_ROT\t= %d\n", input->max_num_rot);
    fprintf(fp, "TRIAL_ANGLE\t= %f\n", input->trial_angle * 180 / 3.141592);
    fputs("\n", fp);

    fputs("# art_nouveau parameter #\n", fp);
    fprintf(fp, "ART_NOUVEAU\t= %d\n", input->art_nouveau);
    fprintf(fp, "LAMBDA_CRIT\t= %f\n", input->lambda_crit);
    fprintf(fp, "LAMBDA_CONV\t= %f\n", input->lambda_conv);
    fprintf(fp, "MAX_NUM_RLX\t= %d\n", input->max_num_rlx);
    fprintf(fp, "ART_DELAY\t= %d\n", input->art_delay);
    fprintf(fp, "ART_MIXING\t= %d\n", input->art_mixing);
    fputs("\n", fp);

    fputs("# system parameter #\n", fp);
    fprintf(fp, "FREQUENCY\t= %f\n", input->frequency);
    fprintf(fp, "TEMPERATURE\t= %f\n", input->temperature);
    fputs("\n", fp);

    fputs("# random parameter #\n", fp);
    fprintf(fp, "RANDOM_SEED\t= %d\n", input->random_seed);
    fputs("\n", fp);

    fputs("# directory parameter #\n", fp);
    fprintf(fp, "OUTPUT_DIR\t= %s\n", input->output_dir);
    fputs("\n", fp);

    fputs("# restart parameter #\n", fp);
    fprintf(fp, "RESTART\t\t= %d\n", input->restart);
    fprintf(fp, "RESTART_DIR\t= %s\n", input->restart_dir);
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
    free(input->vasp_cmd);
    free(input->output_dir);
    free(input->restart_dir);
    free(input);
}
