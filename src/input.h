#ifndef __INITIAL_H__
#define __INITIAL_H__

typedef struct _Input
{
    int nelem;
    int init_relax;
    int max_num_rot;
    int random_seed;
    int ncore;
    int nredundant;
    int restart;

    char *pair_style;
    char *pair_coeff;
    char *init_config;
    char *target_list;
    char **atom_type;
    char *output_dir;
    char *dataset_dir;

    double pair_cutoff;
    double dimer_dist;
    double f_tol;
    double f_rot_min;
    double f_rot_max;
    double trial_angle;
    double disp_cutoff;
    double stddev;
    double max_step;
    double trial_step;
    double confidence;
    double temperature;
    double frequency;
} Input;

int input_int(int *, char *, char *);
int input_double(double *, char *, char *);
int input_char(char **, char *, char *);
int input_char_arr(char ***, char *, int, char *);
int read_input(Input *, char *);
void write_input(Input *);
void free_input(Input *);
#endif
