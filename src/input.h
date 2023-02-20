#ifndef __INITIAL_H__
#define __INITIAL_H__

typedef struct _Input
{
    int nelem;
    char **atom_type;
    char *init_config;
    char *target_list;
    double disp_dist;
    double acti_cutoff;
    double calc_cutoff;
    double f_tol;
    double max_move;
    double trial_move;
    int init_relax;
    int init_disp;
    double disp_cutoff;
    double disp_stddev;
    double confidence;
    int max_search;
    int write_mode;

    char *pair_style;
    char *pair_coeff;
    int ncore;

    char *vasp_cmd;

    int kappa_dimer;
    double f_rot_min;
    double f_rot_max;
    int max_num_rot;
    int max_num_tls;

    int art_nouveau;
    double lambda_crit;
    double lambda_conv;
    int max_num_rlx;
    int art_delay;
    int art_mixing;

    double frequency;
    double temperature;

    int random_seed;

    char *output_dir;

    int restart;
    char *restart_dir;

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
