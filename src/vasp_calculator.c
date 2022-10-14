#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "config.h"
#include "my_mpi.h"
#include "vasp_calculator.h"


void vasp_run(Input *input)
{
    FILE *fp;
    char line[1024];
    chdir("tmp");
    fp = popen(input->vasp_cmd, "r");
    if (fp != NULL) {
        char line[1024];
        while (fgets(line, 1024, fp)) {
            if (strstr(line, "General timing") != NULL) {
                break;
            }
        }
        pclose(fp);
    }
    chdir("../");
}


void modify_incar(Input *input, char *filename, int ibrion)
{
    FILE *wp = fopen(filename, "w");
    char line[1024];
    fputs("# Overwritten by SPS #\n", wp);
    sprintf(line, "ISTART    =    %d\n", input->istart);
    fputs(line, wp);
    if (input->istart > 0) {
        fputs("ICHARG    =    0\n", wp);
        fputs("LWAVE     =    TRUE\n", wp);
    }
    sprintf(line, "IBRION    =    %d\n", ibrion);
    fputs(line, wp);
    if (ibrion == -1) {
        fputs("NSW       =    0\n", wp);
    } else {
        fputs("NSW       =    1000\n", wp);
        fputs("POTIM     =    0.5\n", wp);
    }
    sprintf(line, "EDIFFG    =    -%f\n", input->f_tol);
    fputs(line, wp);
    fputs("\n", wp);
    FILE *rp = fopen("INCAR", "r");
    while (fgets(line, 1024, rp)) {
        fputs(line, wp);
    }
    fclose(rp);
    fclose(wp);
}


void read_outcar(char *filename, double *pos, double *energy, double *force)
{
    int i, n;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("open fail?\n");
    }
    char line[1024], *ptr;
    n = 0;
    while (1) {
        ptr = fgets(line, 1024, fp);
        if (ptr == NULL) {
            break;
        }
        /* # of ions */
        ptr = strstr(line, "NIONS");
        if (ptr != NULL) {
            strtok(ptr, " \n\t");
            strtok(NULL, " \n\t");
            n = atoi(strtok(NULL, "\n"));
        }
        /* position and force */
        ptr = strstr(line, "POSITION");
        if (ptr != NULL) {
            fgets(line, 1024, fp);
            for (i = 0; i < n; ++i) {
                fgets(line, 1024, fp);
                pos[i * 3 + 0] = atof(strtok(line, " \n\t"));
                pos[i * 3 + 1] = atof(strtok(NULL, " \n\t"));
                pos[i * 3 + 2] = atof(strtok(NULL, " \n\t"));
                force[i * 3 + 0] = atof(strtok(NULL, " \n\t"));
                force[i * 3 + 1] = atof(strtok(NULL, " \n\t"));
                force[i * 3 + 2] = atof(strtok(NULL, " \n\t"));
            }
        }
        /* energy */
        ptr = strstr(line, "free  ");
        if (ptr != NULL) {
            strtok(ptr, " \n\t");
            strtok(NULL, " \n\t");
            strtok(NULL, " \n\t");
            strtok(NULL, " \n\t");
            *energy = atof(strtok(NULL, " \n\t"));
        }
    }
    fclose(fp);
}


void append_outcar(char *filename1, char *filename2)
{
    FILE *rp = fopen(filename1, "r");
    FILE *ap = fopen(filename2, "a");
    char line[1024];
    while (fgets(line, 1024, rp)) {
        fputs(line, ap);
    }
    fclose(rp);
    fclose(ap);
}


void oneshot(Config *config, Input *input, double *energy, double *force,
             MPI_Comm comm)
{
    FILE *fp;
    char cmd[1024], filename[128];
    mkdir("tmp", 0775);
    modify_incar(input, "tmp/INCAR", -1);
    write_config(config, "tmp/POSCAR", "w");
    fp = popen("cp POTCAR tmp", "r");
    pclose(fp);
    fp = popen("cp KPOINTS tmp", "r");
    pclose(fp);
    vasp_run(input);
    read_outcar("tmp/OUTCAR", config->pos, energy, force);
    sprintf(filename, "%s/OUTCAR", input->output_dir);
    append_outcar("tmp/OUTCAR", filename);
}


void oneshot_disp(Config *config, Input *input, double *energy, double *force,
                  int disp_num, int *disp_list, MPI_Comm comm)
{
    int i;
    FILE *fp;
    char cmd[1024], filename[128];
    mkdir("tmp", 0775);
    modify_incar(input, "tmp/INCAR", -1);
    write_config(config, "tmp/POSCAR", "w");
    fp = popen("cp POTCAR tmp", "r");
    pclose(fp);
    fp = popen("cp KPOINTS tmp", "r");
    pclose(fp);
    vasp_run(input);
    double *tmp_force = (double *)malloc(sizeof(double) * config->tot_num * 3);
    read_outcar("tmp/OUTCAR", config->pos, energy, tmp_force);
    for (i = 0; i < disp_num; ++i) {
        force[i * 3 + 0] = tmp_force[disp_list[i] * 3 + 0];
        force[i * 3 + 1] = tmp_force[disp_list[i] * 3 + 1];
        force[i * 3 + 2] = tmp_force[disp_list[i] * 3 + 2];
    }
    free(tmp_force);
    sprintf(filename, "%s/OUTCAR", input->output_dir);
    append_outcar("tmp/OUTCAR", filename);
}


void atom_relax(Config *config, Input *input, double *energy, MPI_Comm comm)
{
    FILE *fp;
    char cmd[1024], filename[128];
    mkdir("tmp", 0775);
    modify_incar(input, "tmp/INCAR", 2);
    write_config(config, "tmp/POSCAR", "w");
    fp = popen("cp POTCAR tmp", "r");
    pclose(fp);
    fp = popen("cp KPOINTS tmp", "r");
    pclose(fp);
    vasp_run(input);
    double *tmp_force = (double *)malloc(sizeof(double) * config->tot_num * 3);
    read_outcar("tmp/OUTCAR", config->pos, energy, tmp_force);
    free(tmp_force);
    sprintf(filename, "%s/OUTCAR", input->output_dir);
    append_outcar("tmp/OUTCAR", filename);
}
