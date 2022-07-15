#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "utils.h"


#define MAXLINE 128
int read_config(Config *config, Input *input, char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    int i, j, k;
    double tmp_pos[3];
    char line[MAXLINE], tmp_line[MAXLINE], *ptr;

    /* system name */
    ptr = fgets(line, MAXLINE, fp);

    /* scale */
    ptr = fgets(line, MAXLINE, fp);
    double scale = atof(line);

    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        ptr = fgets(line, MAXLINE, fp);
        config->cell[i][0] = atof(strtok(line, " \n")) * scale;
        config->cell[i][1] = atof(strtok(NULL, " \n")) * scale;
        config->cell[i][2] = atof(strtok(NULL, " \n")) * scale;
    }

    /* the number of type */
    ptr = fgets(line, MAXLINE, fp);
    strncpy(tmp_line, line, MAXLINE);
    config->ntype = 0;
    ptr = strtok(line, " \r\n");
    while (ptr != NULL) {
        if (strlen(ptr) > 0) {
            config->ntype++;
        }
        ptr = strtok(NULL, " \r\n");
    }

    /* atomic number type */
    ptr = strtok(tmp_line, " \n");
    config->atom_num = (int *)malloc(sizeof(int) * config->ntype);
    for (i = 0; i < config->ntype; ++i) {
        config->atom_num[i] = get_atom_num(ptr);
        ptr = strtok(NULL, " \n");
    }

    /* each number of type */
    config->each_num = (int *)malloc(sizeof(int) * config->ntype);
    int *start_idx = (int *)calloc(config->ntype, sizeof(int));
    config->tot_num = 0;
    ptr = fgets(line, MAXLINE, fp);
    config->each_num[0] = atoi(strtok(line, " \n"));
    config->tot_num += config->each_num[0];
    for (i = 1; i < config->ntype; ++i) {
        config->each_num[i] = atoi(strtok(NULL, " \n"));
        config->tot_num += config->each_num[i];
        start_idx[i] = config->each_num[i - 1] + start_idx[i - 1];
    }

    /* type index of atom */
    int count = 0;
    config->type = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = 0; i < config->ntype; ++i) {
        for (j = 0; j < input->nelem; ++j) {
            if (config->atom_num[i] == get_atom_num(input->atom_type[j])) {
                for (k = 0; k < config->each_num[i]; ++k) {
                    config->type[count] = j + 1;
                    count++;
                }
            }
        }
    }
    free(start_idx);

    /* positions and constraint */
    ptr = fgets(line, MAXLINE, fp);
    config->pos = (double *)malloc(sizeof(double) * config->tot_num * 3);
    if (strncasecmp(line, "S", 1) == 0) {
        ptr = fgets(line, MAXLINE, fp);
    }
    if (strncasecmp(line, "D", 1) == 0) {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, MAXLINE, fp);
            tmp_pos[0] = atof(strtok(line, " \n"));
            tmp_pos[1] = atof(strtok(NULL, " \n"));
            tmp_pos[2] = atof(strtok(NULL, " \n"));
            for (j = 0; j < 3; ++j) {
                config->pos[i * 3 + j] = tmp_pos[0] * config->cell[0][j]
                                       + tmp_pos[1] * config->cell[1][j]
                                       + tmp_pos[2] * config->cell[2][j];
            }
        }
    } else {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, MAXLINE, fp);
            config->pos[i * 3 + 0] = atof(strtok(line, " \n"));
            config->pos[i * 3 + 1] = atof(strtok(NULL, " \n"));
            config->pos[i * 3 + 2] = atof(strtok(NULL, " \n"));
        }
    }
    fclose(fp);
    return 0;
}


void write_config(Config *config, char *filename)
{
    int i;
    char line[MAXLINE];
    FILE *fp = fopen(filename, "w");

    /* title */
    fputs("POSCAR from MINK\n", fp);

    /* scale */
    fputs("1.0\n", fp);
    
    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        sprintf(line, " %19.16f", config->cell[i][0]);
        fputs(line, fp);
        sprintf(line, " %19.16f", config->cell[i][1]);
        fputs(line, fp);
        sprintf(line, " %19.16f\n", config->cell[i][2]);
        fputs(line, fp);
    }

    /* symbols */
    for (i = 0; i < config->ntype; ++i) {
        sprintf(line, " %s", get_symbol(config->atom_num[i]));
        fputs(line, fp);
    }
    fputs("\n", fp);

    /* the number of each type */
    for (i = 0; i < config->ntype; ++i) {
        sprintf(line, " %d", config->each_num[i]);
        fputs(line, fp);
    }
    fputs("\n", fp);

    /* positions and constraint */
    fputs("Cartesian\n", fp);
    for (i = 0; i < config->tot_num; ++i) {
        sprintf(line, "  %19.16f", config->pos[i * 3 + 0]);
        fputs(line, fp);
        sprintf(line, "  %19.16f", config->pos[i * 3 + 1]);
        fputs(line, fp);
        sprintf(line, "  %19.16f\n", config->pos[i * 3 + 2]);
        fputs(line, fp);
    }
    fclose(fp);
}


void copy_config(Config *config2, Config *config1)
{
    int i;

    config2->ntype = config1->ntype;
    config2->tot_num = config1->tot_num;
    config2->cell[0][0] = config1->cell[0][0];
    config2->cell[0][1] = config1->cell[0][1];
    config2->cell[0][2] = config1->cell[0][2];
    config2->cell[1][0] = config1->cell[1][0];
    config2->cell[1][1] = config1->cell[1][1];
    config2->cell[1][2] = config1->cell[1][2];
    config2->cell[2][0] = config1->cell[2][0];
    config2->cell[2][1] = config1->cell[2][1];
    config2->cell[2][2] = config1->cell[2][2];

    config2->atom_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->each_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->type = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->pos = (double *)malloc(sizeof(double) * config1->tot_num * 3);

    for (i = 0; i < config1->ntype; ++i) {
        config2->atom_num[i] = config1->atom_num[i];
        config2->each_num[i] = config1->each_num[i];
    }
    for (i = 0; i < config1->tot_num; ++i) {
        config2->type[i] = config1->type[i];
        config2->pos[i * 3 + 0] = config1->pos[i * 3 + 0];
        config2->pos[i * 3 + 1] = config1->pos[i * 3 + 1];
        config2->pos[i * 3 + 2] = config1->pos[i * 3 + 2];
    }
}


void free_config(Config *config)
{
    free(config->atom_num);
    free(config->each_num);
    free(config->type);
    free(config->pos);
    free(config);
}