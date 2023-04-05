#include "calculator.h"
#define LAMMPS_LIB_MPI
#include "library.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void *lmp_init(Config *config, Input *input, MPI_Comm comm)
{
    /* create LAMMPS instance */
    int i, j, k;
    void *lmp;
    char cmd[65536], tmp_cmd[65536];
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
//    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lammps_open(lmpargc, lmpargv, comm, NULL);
    if (lmp == NULL) {
        printf("LAMMPS initialization failed");
    }
    /* basic */
    const char *cmds[] = {"units metal",
                          "neigh_modify every 1 delay 0 check yes",
                          "atom_modify map array sort 0 0.0",
                          "box tilt large"};
    lammps_commands_list(lmp, sizeof(cmds) / sizeof(const char *), cmds);
    /* box */
    sprintf(cmd, "region cell prism 0 %f 0 %f 0 %f %f %f %f units box",
            config->cell[0][0], config->cell[1][1], config->cell[2][2],
            config->cell[1][0], config->cell[2][0], config->cell[2][1]);
    lammps_command(lmp, cmd);
    sprintf(cmd, "create_box %d cell", input->nelem);
    lammps_command(lmp, cmd);
    /* atoms */
    int count = 0;
    int *id = (int *)malloc(sizeof(int) * config->tot_num);
    int *type = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = 0; i < config->ntype; ++i) {
        for (j = 0; j < input->nelem; ++j) {
            if (config->atom_num[i] == get_atom_num(input->atom_type[j])) {
                for (k = 0; k < config->each_num[i]; ++k) {
                    type[count] = j + 1;
                    id[count] = count + 1;
                    count++;
                }
            }
        }
    }
    lammps_create_atoms(lmp, config->tot_num, id,
                        type, config->pos, NULL, NULL, 0);
    free(id);
    free(type);
    /* mass */
    for (i = 0; i < input->nelem; ++i) {
        sprintf(cmd, "mass %d %f", i + 1,
                get_mass(get_atom_num(input->atom_type[i])));
        lammps_command(lmp, cmd);
    }
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    char *ptr;
    char tmp_pair_coeff[65536];
    memcpy(tmp_pair_coeff, input->pair_coeff, sizeof(char) * 65536);
    ptr = strtok(tmp_pair_coeff, "|\n");
    while (ptr != NULL) {
        sprintf(cmd, "pair_coeff %s", ptr);
        lammps_command(lmp, cmd);
        ptr = strtok(NULL, "|\n");
    }
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 20 1.0");
    /* fix */
    int fix = 0;
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] > 0) {
            fix++;
            break;
        }
    }
    if (fix > 0) {
        sprintf(cmd, "group freeze id");
        for (i = 0; i < config->tot_num; ++i) {
            if (config->fix[i] > 0) {
                sprintf(tmp_cmd, " %d", i + 1);
                strcat(cmd, tmp_cmd);
            }
        }
        lammps_command(lmp, cmd);
        lammps_command(lmp, "fix int freeze setforce 0.0 0.0 0.0");
    }
    return lmp;
}


void oneshot(Config *config, Input *input, double *energy, double *force,
             MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    lmp = lmp_init(config, input, comm);
    /* oneshot */
    lammps_command(lmp, "run 0");
    *energy = lammps_get_thermo(lmp, "pe");
    lammps_gather_atoms(lmp, "f", 1, 3, force);
    /* delete LAMMPS instance */
    lammps_close(lmp);
}


void atom_relax(Config *config, Input *input, double *energy, MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    lmp = lmp_init(config, input, comm);
    /* minimize */
    sprintf(cmd, "minimize 0 %f 10000 100000", input->f_tol);
    lammps_command(lmp, cmd);
    *energy = lammps_get_thermo(lmp, "pe");
    /* update positions */
    lammps_gather_atoms(lmp, "x", 1, 3, config->pos);
    /* delete LAMMPS instance */
    lammps_close(lmp);
}
