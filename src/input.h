#ifndef __INITIAL_H__
#define __INITIAL_H__

typedef struct _Input
{
    int nelem;
    char **atom_type;
    char *pair_style;
    char *pair_coeff;
    double pair_cutoff;

    char *init_config;
    char *target_list;
    double disp_dist;
    double acti_cutoff;
    double f_tol;
    double stddev;
    double max_step;
    double trial_step;
    int init_relax;
    double confidence;

    double f_rot_min;
    double f_rot_max;
    int max_num_rot;
    double trial_angle;
    int kappa_dimer;
    int snc_dimer;

    int art_nouveau;
    double lambda_crit;
    double lambda_conv;
    int max_num_rlx;
    int art_delay;

    double frequency;
    double temperature;

    int random_seed;

    char *output_dir;

    int restart;
    char *restart_dir;

    int ncore;

    int nredundant;
} Input;

int input_int(int *var, char *tag, char *filename);
int input_double(double *var, char *tag, char *filename);
int input_char(char **var, char *tag, char *filename);
int input_char_arr(char ***var, char *tag, int n, char *filename);
int read_input(Input *input, char *filename);
void write_input(Input *input);
void free_input(Input *input);
#endif
