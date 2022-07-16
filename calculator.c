#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#define LAMMPS_LIB_MPI
#include "library.h"
#include "utils.h"


void *lmp_init(Config *config, Input *input, int lmpargc, char **lmpargv)
{
    /* create LAMMPS instance */
    int i;
    void *lmp;
    char cmd[1024];
    lmp = lammps_open(lmpargc, lmpargv, MPI_COMM_WORLD, NULL);
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
    lammps_create_atoms(lmp, config->tot_num, config->id,
                        config->type, config->pos, NULL, NULL, 0);
    for (i = 0; i < input->nelem; ++i) {
        sprintf(cmd, "mass %d %f", i + 1, get_mass(get_atom_num(input->atom_type[i])));
        lammps_command(lmp, cmd);
    }
    return lmp;
}


void oneshot(Config *config, Input *input, double *energy, double ***force,
             int disp_num, int *disp_list)
{
    int i;
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    //char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* oneshot */
    lammps_command(lmp, "run 0");
    *energy = lammps_get_thermo(lmp, "pe");
    double **tmp_force = (double **)lammps_extract_atom(lmp, "f");
    for (i = 0;  i < disp_num; ++i) {
        (*force)[i][0] = tmp_force[disp_list[i]][0];
        (*force)[i][1] = tmp_force[disp_list[i]][1];
        (*force)[i][2] = tmp_force[disp_list[i]][2];
    }
    /* delete LAMMPS instance */
    lammps_close(lmp);
}


double atom_relax(Config *config, Input *input)
{
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    //char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* minimize */
    sprintf(cmd, "minimize 0 %f 10000 100000", input->ftol);
    lammps_command(lmp, cmd);
    double pe = lammps_get_thermo(lmp, "pe");
    /* update positions */
    lammps_gather_atoms(lmp, "x", 2, 3, config->pos);
    /* delete LAMMPS instance */
    lammps_close(lmp);

    return pe;
}
