#include "calculator.h"
#include "my_mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>


void vasp_run(Input *input)
{
    FILE *fp;
    char line[1024];
    chdir("./VASP_tmp");
    fp = popen(input->vasp_cmd, "r");
    if (fp != NULL) {
        while (fgets(line, 1024, fp)) {
            if (strstr(line, "General timing") != NULL) {
                break;
            }
        }
        pclose(fp);
    }
    chdir("../");
}


static void write_incar(Input *input, char *filename, int ibrion)
{
    FILE *wp = fopen(filename, "w");
    char line[1024];
    fputs("# Written by SHERPA #\n", wp);
    fputs("ISTART    =    1\n", wp);
    fputs("ICHARG    =    1\n", wp);
    fputs("LWAVE     =    TRUE\n", wp);
    fputs("LCHARG    =    TRUE\n", wp);
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


void oneshot(Calc *calc, Config *config, Input *input,
             double *energy, double *force, MPI_Comm comm)
{
    FILE *fp;
    mkdir("./VASP_tmp", 0775);
    write_incar(input, "./VASP_tmp/INCAR", -1);
    write_config(config, "./VASP_tmp/POSCAR", "poscar", "w");
    fp = popen("cp POTCAR VASP_tmp", "r");
    pclose(fp);
    fp = popen("cp KPOINTS VASP_tmp", "r");
    pclose(fp);
    vasp_run(input);
    read_outcar("./VASP_tmp/OUTCAR", config->pos, energy, force);
    fp = popen("cat VASP_tmp/OUTCAR >> OUTCAR_concat", "r");
    pclose(fp);
}


void atom_relax(Calc *calc, Config *config, Input *input,
                double *energy, MPI_Comm comm)
{
    FILE *fp;
    mkdir("./VASP_tmp", 0775);
    write_incar(input, "VASP_tmp/INCAR", 2);
    write_config(config, "VASP_tmp/POSCAR", "poscar", "w");
    fp = popen("cp POTCAR VASP_tmp", "r");
    pclose(fp);
    fp = popen("cp KPOINTS VASP_tmp", "r");
    pclose(fp);
    vasp_run(input);
    double *tmp_force = (double *)malloc(sizeof(double) * config->tot_num * 3);
    read_outcar("./VASP_tmp/OUTCAR", config->pos, energy, tmp_force);
    free(tmp_force);
    fp = popen("cat VASP_tmp/OUTCAR >> OUTCAR_concat", "r");
    pclose(fp);
}


void free_calc(Calc *calc)
{
    free(calc);
}
