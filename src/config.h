#ifndef __CONFIG_H__
#define __CONFIG_H__

typedef struct _Config
{
    int ntype;            /* the number of types */
    int tot_num;          /* the total number of atoms */
    double cell[3][3];    /* lattice vector */
    double boxlo[3];      /* edge */
    double boxhi[3];
    double xy;            /* tilting */
    double xz;
    double yz;

    int *atom_num;        /* atomic number for symbols */
    int *each_num;        /* the number of each types */

    int *fix;
    double *pos;          /* [N * 3] dimension of positions */
} Config;

int get_atom_num(char *symbol);
double get_mass(int atom_num);
char *get_symbol(int atom_num);
void remove_atom(Config *config, int index);
int read_config(Config *config, char *filename);
void write_config(Config *config, char *filename, char *header, char *mode);
void copy_config(Config *config2, Config *config1);
void free_config(Config *config);
#endif
