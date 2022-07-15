#ifndef __INITIAL_H__
#define __INITIAL_H__

typedef struct _Input
{
    int nelem;
    int init_mode;
    int random_seed;

    char *pair_style;
    char *pair_coeff;
    char *init_config;
    char *target;
    char **atom_type;

    double cutoff;
    double dimer_dist;
    double fmax;
    double ftol;
    double trial_angle;
    double disp_cutoff;
    double stddev;
} Input;

int input_int(int *, char *, char *);
int input_double(double *, char *, char *);
int input_char(char **, char *, char *);
int input_char_arr(char ***, char *, int, char *);
int read_input(Input *, char *);
void free_input(Input *);
#endif
