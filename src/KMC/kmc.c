#include "input.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void copy_files(char *filename2, char *filename1)
{
    FILE *rp = fopen(filename1, "r");
    FILE *wp = fopen(filename2, "w");
    char line[1024];
    while (fgets(line, 1024, rp) != NULL) {
        fputs(line, wp);
    }
    fclose(rp);
    fclose(wp);
}


int main(int argc, char *argv[])
{
    int i, kmc_step, errno;
    char directory[1024], filename1[1024], filename2[1024];

    /* read input */
    Input *input = (Input *)malloc(sizeof(Input));
    errno = read_input(input, "./INPUT_KMC");
    if (errno > 0) {
        printf("ERROR in INPUT FILE!\n");
        free(input);
        return 1;
    }
    /* write_input */
    write_input(input);

    for (kmc_step = 1; kmc_step <= input->kmc_step; ++kmc_step) {
        /* mkdir */
        sprintf(directory, "%d", kmc_step);
        mkdir(directory, 0775);
        /* cp INPUTS */
        for (i = 0; i < input->num_inputs; ++i) {
            sprintf(filename1, "%s", input->inputs[i]);
            sprintf(filename2, "%d/%s", kmc_step, input->inputs[i]);
            copy(filename2, filename1);
        }
        /* chdir */
        chdir(directory);
        /* run sherpa */
        FILE *pp = popen(input->sherpa_cmd, "r");
        FILE *fp;
        if (pp != NULL) {
            char filename[1024];
            sprintf(filename, "%s/Event.log");
            while (1) {
                if (fopen(filename, "r") != NULL) {
                    break;
                }
            }
            pclose(pp);
        } else {
            printf("Check SHERPA path\n");
            free_input(input);
            return 1;
        }
        /* select */
        /* chdir */
        /* log */
    }

    char filename[1024];
    sprintf(filename, "%s/Event.log", argv[1]);
    FILE *rp = fopen(filename, "r");
    if (rp == NULL) {
        printf("Check your directory.");
        return 1;
    }

    char line[1024], tmp_line[1024], *ptr;
    double Ea, rate;
    /* skip header */
    fgets(line, 1024, rp);
    fgets(line, 1024, rp);
    /* read Event.log */
    fgets(line, 1024, rp);
    int list_size = 1024;
    int reac_num = 0;
    int *reac_list = (int *)malloc(sizeof(int) * list_size);
    double *acti_list = (double *)malloc(sizeof(double) * list_size);
    double *rate_list = (double *)malloc(sizeof(double) * list_size);
    double kT = 8.61733326E-5 * atof(argv[2]);
    while (fgets(line, 1024, rp) != NULL) {
        reac_list[reac_num] = atoi(strtok(line, " \n"));
        Ea = atof(strtok(NULL, " \n"));
        acti_list[reac_num] = Ea;
        rate_list[reac_num] = atof(argv[3]) * exp(-Ea / kT);
        reac_num++;
        if (reac_num >= list_size) {
            list_size = list_size << 1;
            reac_list = (int *)realloc(reac_list, sizeof(int) * list_size);
            acti_list = (double *)realloc(acti_list, sizeof(double) * list_size);
            rate_list = (double *)realloc(rate_list, sizeof(double) * list_size);
        }
    }
    fclose(rp);
    /* reaction index */
    int i;
    double rate_sum = 0.0;
    for (i = 0; i < reac_num; ++i) {
        rate_sum += rate_list[i];
    }
    int index, tmp_index;
    double acc_rate_sum = 0.0;
    double rn = (double)random() / RAND_MAX;
    for (i = 0; i < reac_num; ++i) {
        rate_list[i] /= rate_sum;
        acc_rate_sum += rate_list[i];
        if (acc_rate_sum > rn) {
            index = reac_list[i];    
            Ea = acti_list[i];
            rate = rate_list[i] * rate_sum;
            break;
        }
    }
    printf("%d %f %f\n", index, Ea, rate);

    free(reac_list);
    free(acti_list);
    free(rate_list);
    return 0;
}
