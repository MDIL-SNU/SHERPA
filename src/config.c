#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "utils.h"


double norm3(double *vec)
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}


double dot3(double *vec1, double *vec2)
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}


void cross3(double *vec1, double *vec2, double *vec3)
{
    vec3[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vec3[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vec3[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}


double det(double (*mat)[3])
{
    int i;
    double cofactor[3];
    cofactor[0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    cofactor[1] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    cofactor[2] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];

    double det = 0.0;
    for (i = 0; i < 3; i++) {
        det += cofactor[i] * mat[0][i];
    }
    return det;
}


/* convert into lmp basis */
void convert_basis(Config *config)
{
    int i, j;
    double A_norm, AxB_norm, vol, tmp_value, tmp_pos[3];
    double A[3], B[3], C[3];
    double Ahat[3], AxB[3], AxBhat[3], AhatxB[3], AxBhatxAhat[3];
    double frac[3][3], trans[3][3];

    /* rotate bases into triangular matrix */
    for (i = 0; i < 3; i++) {
        A[i] = config->cell[0][i];
        B[i] = config->cell[1][i];
        C[i] = config->cell[2][i];
    }

    A_norm = norm3(A);
    Ahat[0] = A[0] / A_norm;
    Ahat[1] = A[1] / A_norm;
    Ahat[2] = A[2] / A_norm;

    cross3(A, B, AxB);
    AxB_norm = norm3(AxB);
    AxBhat[0] = AxB[0] / AxB_norm;
    AxBhat[1] = AxB[1] / AxB_norm;
    AxBhat[2] = AxB[2] / AxB_norm;

    cross3(Ahat, B, AhatxB);
    cross3(AxBhat, Ahat, AxBhatxAhat);

    /* column vector (a b c) */
    config->cell[0][0] = A_norm;
    config->cell[0][1] = dot3(B, Ahat);
    config->cell[0][2] = dot3(C, Ahat);
    config->cell[1][0] = 0.0;
    config->cell[1][1] = norm3(AhatxB);
    config->cell[1][2] = dot3(C, AxBhatxAhat);
    config->cell[2][0] = 0.0;
    config->cell[2][1] = 0.0;
    config->cell[2][2] = dot3(C, AxBhat);

    /* edge & tilting */
    config->boxlo[0] = 0.0;
    config->boxlo[1] = 0.0;
    config->boxlo[2] = 0.0;
    config->boxhi[0] = config->cell[0][0];
    config->boxhi[1] = config->cell[1][1];
    config->boxhi[2] = config->cell[2][2];
    config->xy = config->cell[0][1];
    config->xz = config->cell[0][2];
    config->yz = config->cell[1][2];

    /* transformation matrix */
    vol = det(config->cell);
    cross3(B, C, frac[0]);
    cross3(C, A, frac[1]);
    cross3(A, B, frac[2]);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            frac[i][j] /= vol;
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            trans[i][j] = config->cell[i][0] * frac[0][j] +
                          config->cell[i][1] * frac[1][j] +
                          config->cell[i][2] * frac[2][j];
        }
    }

    /* convert position */
    for (i = 0; i < config->tot_num; i++) {
        for (j = 0; j < 3; j++) {
            tmp_pos[j] = trans[j][0] * config->pos[i * 3 + 0] +
                         trans[j][1] * config->pos[i * 3 + 1] +
                         trans[j][2] * config->pos[i * 3 + 2];
        }
        config->pos[i * 3 + 0] = tmp_pos[0];
        config->pos[i * 3 + 1] = tmp_pos[1];
        config->pos[i * 3 + 2] = tmp_pos[2];
    }

    /* transpose lattice vector */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 3; j++) {
            tmp_value = config->cell[i][j];
            config->cell[i][j] = config->cell[j][i];
            config->cell[j][i] = tmp_value;
        }
    }
}


/* remove one atom from config */
void extract_atom(Config *config, int idx)
{
    int i, j, ntype, tot_num;
    int acc_num = 0;
    for (i = 0; i < config->ntype; ++i) {
        acc_num += config->each_num[i];
        if (acc_num > idx) {
            config->each_num[i]--;
            if (config->each_num[i] == 0) {
                for (j = i; j < config->ntype - 1; ++j) {
                    config->atom_num[i] = config->atom_num[i + 1];
                    config->each_num[i] = config->each_num[i + 1];
                }
                config->ntype--;
                ntype = config->ntype;
                config->atom_num = (int *)realloc(config->atom_num, sizeof(int) * ntype);
                config->each_num = (int *)realloc(config->each_num, sizeof(int) * ntype);
            }
            break;
        }
    }

    for (i = idx; i < config->tot_num - 1; ++i) {
        config->fix[i] = config->fix[i + 1];
        config->type[i] = config->type[i + 1];
        config->pos[i * 3 + 0] = config->pos[(i + 1) * 3 + 0];
        config->pos[i * 3 + 1] = config->pos[(i + 1) * 3 + 1];
        config->pos[i * 3 + 2] = config->pos[(i + 1) * 3 + 2];
    }
    config->tot_num--;
    tot_num = config->tot_num;
    config->id = (int *)realloc(config->id, sizeof(int) * tot_num);
    config->fix = (int *)realloc(config->fix, sizeof(int) * tot_num);
    config->type = (int *)realloc(config->type, sizeof(int) * tot_num);
    config->pos = (double *)realloc(config->pos, sizeof(double) * tot_num * 3);
}


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
    config->id = (int *)malloc(sizeof(int) * config->tot_num);
    config->type = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = 0; i < config->ntype; ++i) {
        for (j = 0; j < input->nelem; ++j) {
            if (config->atom_num[i] == get_atom_num(input->atom_type[j])) {
                for (k = 0; k < config->each_num[i]; ++k) {
                    config->type[count] = j + 1;
                    config->id[count] = count + 1;
                    count++;
                }
            }
        }
    }
    free(start_idx);

    /* positions and constraint */
    ptr = fgets(line, MAXLINE, fp);
    config->fix = (int *)calloc(config->tot_num, sizeof(int));
    config->pos = (double *)malloc(sizeof(double) * config->tot_num * 3);
    int constraint = 0;
    if (strncasecmp(line, "S", 1) == 0) {
        constraint = 1;
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
            if (constraint > 0) {
                ptr = strtok(NULL, " \n");
                if (strcmp(ptr, "F") == 0) {
                    config->fix[i] = 1;
                }
            }
        }
    } else {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, MAXLINE, fp);
            config->pos[i * 3 + 0] = atof(strtok(line, " \n"));
            config->pos[i * 3 + 1] = atof(strtok(NULL, " \n"));
            config->pos[i * 3 + 2] = atof(strtok(NULL, " \n"));
            if (constraint > 0) {
                ptr = strtok(NULL, " \n");
                if (strcmp(ptr, "F") == 0) {
                    config->fix[i] = 1;
                }
            }
        }
    }
    convert_basis(config);
    fclose(fp);
    return 0;
}


void write_config(Config *config, char *filename, char *mode)
{
    int i;
    char line[MAXLINE];
    FILE *fp;
    fp = fopen(filename, mode);

    /* title */
    fputs("POSCAR\n", fp);

    /* scale */
    fputs("1.0\n", fp);
    
    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        sprintf(line, " %.15f", config->cell[i][0]);
        fputs(line, fp);
        sprintf(line, " %.15f", config->cell[i][1]);
        fputs(line, fp);
        sprintf(line, " %.15f\n", config->cell[i][2]);
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
    fputs("Selective dynamics\n", fp);
    fputs("Cartesian\n", fp);
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] > 0) {
            sprintf(line, "  %.15f", config->pos[i * 3 + 0]);
            fputs(line, fp);
            sprintf(line, "  %.15f", config->pos[i * 3 + 1]);
            fputs(line, fp);
            sprintf(line, "  %.15f", config->pos[i * 3 + 2]);
            fputs(line, fp);
            fputs(" F F F\n", fp);
        } else {
            sprintf(line, "  %.15f", config->pos[i * 3 + 0]);
            fputs(line, fp);
            sprintf(line, "  %.15f", config->pos[i * 3 + 1]);
            fputs(line, fp);
            sprintf(line, "  %.15f", config->pos[i * 3 + 2]);
            fputs(line, fp);
            fputs(" T T T\n", fp);
        }
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

    config2->boxlo[0] = config1->boxlo[0];
    config2->boxlo[1] = config1->boxlo[1];
    config2->boxlo[2] = config1->boxlo[2];
    config2->boxhi[0] = config1->boxhi[0];
    config2->boxhi[1] = config1->boxhi[1];
    config2->boxhi[2] = config1->boxhi[2];
    config2->xy = config1->xy;
    config2->xz = config1->xz;
    config2->yz = config1->yz;

    config2->atom_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->each_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->id = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->fix = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->type = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->pos = (double *)malloc(sizeof(double) * config1->tot_num * 3);

    for (i = 0; i < config1->ntype; ++i) {
        config2->atom_num[i] = config1->atom_num[i];
        config2->each_num[i] = config1->each_num[i];
    }
    for (i = 0; i < config1->tot_num; ++i) {
        config2->id[i] = config1->id[i];
        config2->fix[i] = config1->fix[i];
        config2->type[i] = config1->type[i];
        config2->pos[i * 3 + 0] = config1->pos[i * 3 + 0];
        config2->pos[i * 3 + 1] = config1->pos[i * 3 + 1];
        config2->pos[i * 3 + 2] = config1->pos[i * 3 + 2];
    }
}


/* 0: identical, 1: different */
int diff_config(Config *config1, Config *config2, double tol)
{
    int i;
    for (i = 0; i < config1->tot_num; ++i) {
        if (config1->type[i] != config2->type[i]) {
            return 1;
        };
        double disp[3] = {config2->pos[i * 3 + 0] - config1->pos[i * 3 + 0],
                          config2->pos[i * 3 + 1] - config1->pos[i * 3 + 1],
                          config2->pos[i * 3 + 2] - config1->pos[i * 3 + 2]};
        get_minimum_image(disp, config1->boxlo, config1->boxhi,
                          config1->xy, config1->yz, config1->xz);
        if (norm3(disp) > tol) {
            return 1;
        }
    }
    return 0;
}


void free_config(Config *config)
{
    free(config->atom_num);
    free(config->each_num);
    free(config->id);
    free(config->fix);
    free(config->type);
    free(config->pos);
    free(config);
}
