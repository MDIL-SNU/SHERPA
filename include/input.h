#ifndef __INPUT_H__
#define __INPUT_H__

typedef struct _Input
{
    char algorithm;
    double acti_cutoff;
    int acti_nevery;
    double finite_diff;
    double f_tol;
    double diff_tol;
    double max_move;
    double trial_move;
    int max_search;
    int write_traj;
    int cont;

    int nelem;
    char **atom_type;
    int init_relax;
    int init_disp;
    double disp_cutoff;
    double disp_move;
    int init_mode;

    char *pair_style;
    char *pair_coeff;
    int ncore;
    char *vasp_cmd;
    char *ase_calc;
    char *model_path;
    char *device;

    double f_rot_min;
    double f_rot_max;
    int max_num_rot;
    int max_num_tls;

    double lambda_conv;
    int max_num_itr;
    int max_num_rlx;
    int delay_step;
    int mixing_step;
    int hyper_step;

    int random_seed;
} Input;

int read_input(Input *input, char *filename);
void write_input(Input *input);
void free_input(Input *input);
#endif
