#ifndef __CONFIG_H__
#define __CONFIG_H__
#include "input.h"

typedef struct _Config
{
    int ntype;            /* the number of types */
    int tot_num;          /* the total number of atoms */
    double cell[3][3];    /* lattice vector */

    int *atom_num;        /* atomic number for symbols */
    int *each_num;        /* the number of each types */

    int *type;            /* type index starting from 1 */
    double *pos;          /* [N * 3] dimension of positions */
} Config;

int read_config(Config *, Input *, char *);
void write_config(Config *, char *);
void copy_config(Config *, Config *);
void free_config(Config *);
#endif