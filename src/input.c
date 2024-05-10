#include "input.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


int read_input(Input *input, char *filename)
{
    input->acti_cutoff = 6.0;
    input->acti_nevery = 3;
    input->finite_diff = 0.01;
    input->f_tol = 0.01;
    input->diff_tol = 0.4;
    input->max_move = 0.1;
    input->trial_move = 0.01;
    input->max_search = 100;
    input->write_traj = 1;
    input->cont = 0;

    input->nelem = 0;
    input->atom_type = NULL;
    input->init_relax = 1;
    input->init_disp = 0;
    input->disp_cutoff = 3.0;
    input->disp_move = 0.1;
    input->init_mode = 0;

    input->pair_style = NULL;
    input->pair_coeff = NULL;
    input->ncore = 1;

    input->vasp_cmd = NULL;

    input->ase_calc = NULL;
    input->model_path = NULL;

    input->kappa_dimer = 0;
    input->f_rot_min = 0.1;
    input->f_rot_max = 1.0;
    input->max_num_rot = 4;
    input->max_num_tls = 500;

    input->art_nouveau = 1;
    input->lambda_conv = 0.01;
    input->max_num_itr = 500;
    input->max_num_rlx = 1;
    input->delay_step = 0;
    input->mixing_step = 0;
    input->hyper_step = 3;

    input->random_seed = (unsigned int)time(NULL);

    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Check INPUT directory.\n");
        return 1;
    }
    char line[128], *ptr;
    while (1) {
        ptr = fgets(line, 128, fp);
        if (ptr == NULL) {
            break;
        } else if (strcmp(ptr, "\n") == 0 || strncmp(ptr, "#", 1) == 0) {
            continue;
        } else {
            ptr = strtok(line, " \n\t");
            if (strcmp(ptr, "ACTI_CUTOFF") == 0) {
                strtok(NULL, " \n\t");
                input->acti_cutoff = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "ACTI_NEVERY") == 0) {
                strtok(NULL, " \n\t");
                input->acti_nevery = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "FINITE_DIFF") == 0) {
                strtok(NULL, " \n\t");
                input->finite_diff = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "F_TOL") == 0) {
                strtok(NULL, " \n\t");
                input->f_tol = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "DIFF_TOL") == 0) {
                strtok(NULL, " \n\t");
                input->diff_tol = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MAX_MOVE") == 0) {
                strtok(NULL, " \n\t");
                input->max_move = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "TRIAL_MOVE") == 0) {
                strtok(NULL, " \n\t");
                input->trial_move = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MAX_SEARCH") == 0) {
                strtok(NULL, " \n\t");
                input->max_search = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "WRITE_TRAJ") == 0) {
                strtok(NULL, " \n\t");
                input->write_traj = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "CONTINUE") == 0) {
                strtok(NULL, " \n\t");
                input->cont = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "NELEMENT") == 0) {
                strtok(NULL, " \n\t");
                input->nelem = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "ATOM_TYPE") == 0) {
                strtok(NULL, " \n\t");
                input->atom_type = (char **)malloc(sizeof(char *) * input->nelem);
                for (int i = 0; i < input->nelem; ++i) {
                    input->atom_type[i] = (char *)calloc(16, sizeof(char));
                    strcpy(input->atom_type[i], strtok(NULL, " \n\t"));
                }
            } else if (strcmp(ptr, "INIT_RELAX") == 0) {
                strtok(NULL, " \n\t");
                input->init_relax = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "INIT_DISP") == 0) {
                strtok(NULL, " \n\t");
                input->init_disp = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "DISP_CUTOFF") == 0) {
                strtok(NULL, " \n\t");
                input->disp_cutoff = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "DISP_MOVE") == 0) {
                strtok(NULL, " \n\t");
                input->disp_move = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "INIT_MODE") == 0) {
                strtok(NULL, " \n\t");
                input->init_mode = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "PAIR_STYLE") == 0) {
                strtok(NULL, " \n\t");
                input->pair_style = (char *)malloc(sizeof(char) * 65536);
                strcpy(input->pair_style, strtok(NULL, "\n"));
            } else if (strcmp(ptr, "PAIR_COEFF") == 0) {
                strtok(NULL, " \n\t");
                input->pair_coeff = (char *)malloc(sizeof(char) * 65536);
                strcpy(input->pair_coeff, strtok(NULL, "\n"));
            } else if (strcmp(ptr, "NCORE") == 0) {
                strtok(NULL, " \n\t");
                input->ncore = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "VASP_CMD") == 0) {
                strtok(NULL, " \n\t");
                input->vasp_cmd = (char *)malloc(sizeof(char) * 65536);
                strcpy(input->vasp_cmd, strtok(NULL, "\n"));
            } else if (strcmp(ptr, "ASE_CALC") == 0) {
                strtok(NULL, " \n\t");
                input->ase_calc = (char *)malloc(sizeof(char) * 65536);
                strcpy(input->ase_calc, strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MODEL_PATH") == 0) {
                strtok(NULL, " \n\t");
                input->model_path = (char *)malloc(sizeof(char) * 65536);
                strcpy(input->model_path, strtok(NULL, "\n"));
            } else if (strcmp(ptr, "KAPPA_DIMER") == 0) {
                strtok(NULL, " \n\t");
                input->kappa_dimer = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "F_ROT_MIN") == 0) {
                strtok(NULL, " \n\t");
                input->f_rot_min = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "F_ROT_MAX") == 0) {
                strtok(NULL, " \n\t");
                input->f_rot_max = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MAX_NUM_ROT") == 0) {
                strtok(NULL, " \n\t");
                input->max_num_rot = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MAX_NUM_TLS") == 0) {
                strtok(NULL, " \n\t");
                input->max_num_tls = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "ART_NOUVEAU") == 0) {
                strtok(NULL, " \n\t");
                input->art_nouveau = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "LAMBDA_CONV") == 0) {
                strtok(NULL, " \n\t");
                input->lambda_conv = atof(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MAX_NUM_ITR") == 0) {
                strtok(NULL, " \n\t");
                input->max_num_itr = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MAX_NUM_RLX") == 0) {
                strtok(NULL, " \n\t");
                input->max_num_rlx = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "DELAY_STEP") == 0) {
                strtok(NULL, " \n\t");
                input->delay_step = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "MIXING_STEP") == 0) {
                strtok(NULL, " \n\t");
                input->mixing_step = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "HYPER_STEP") == 0) {
                strtok(NULL, " \n\t");
                input->hyper_step = atoi(strtok(NULL, "\n"));
            } else if (strcmp(ptr, "RANDOM_SEED") == 0) {
                strtok(NULL, " \n\t");
                input->random_seed = atoi(strtok(NULL, "\n"));
            } else {
                printf("Unidentified tag. %s\n", ptr);
            }
        }
    }
    fclose(fp);
    if (input->pair_style == NULL && input->pair_coeff == NULL
       && input->vasp_cmd == NULL && input->ase_calc == NULL) {
        printf("Choose one method. LAMMPS, VASP, or ASE.\n");
        return 1;
    }
    /* TODO: dimer is not a default option */
    if (input->kappa_dimer + input->art_nouveau > 1) {
        printf("Choose only one algirhtm. Kappa-dimer or ARTn\n");
        return 1;
    }
    if ((input->vasp_cmd != NULL || input->ase_calc != NULL) && input->ncore != 1) {
        input->ncore = 1;
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
    fprintf(fp, "FINITE_DIFF\t= %f\n", input->finite_diff);
    fprintf(fp, "F_TOL\t\t= %f\n", input->f_tol);
    fprintf(fp, "DIFF_TOL\t= %f\n", input->diff_tol);
    fprintf(fp, "MAX_MOVE\t= %f\n", input->max_move);
    fprintf(fp, "TRIAL_MOVE\t= %f\n", input->trial_move);
    fprintf(fp, "MAX_SEARCH\t= %d\n", input->max_search);
    fprintf(fp, "WRITE_TRAJ\t= %d\n", input->write_traj);
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

    if (input->pair_style != NULL && input->pair_coeff != NULL) {
        fputs("# LAMMPS parameter #\n", fp);
        fprintf(fp, "PAIR_STYLE\t= %s\n", input->pair_style);
        fprintf(fp, "PAIR_COEFF\t= %s\n", input->pair_coeff);
        fprintf(fp, "NCORE\t\t= %d\n", input->ncore);
        fputs("\n", fp);
    }
    if (input->vasp_cmd != NULL) {
        fputs("# VASP parameter #\n", fp);
        fprintf(fp, "VASP_CMD\t= %s\n", input->vasp_cmd);
        fputs("\n", fp);
    }
    if (input->ase_calc != NULL) {
        fputs("# ASE parameter #\n", fp);
        fprintf(fp, "ASE_CALC\t= %s\n", input->ase_calc);
        fputs("\n", fp);
    }

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
    fprintf(fp, "MAX_NUM_ITR\t= %d\n", input->max_num_itr);
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
    if (input->pair_style != NULL) {
        free(input->pair_style);
    }
    if (input->pair_coeff != NULL) {
        free(input->pair_coeff);
    }
    if (input->vasp_cmd != NULL) {
        free(input->vasp_cmd);
    }
    free(input);
}
