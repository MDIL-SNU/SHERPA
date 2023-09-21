#ifndef __INPUT_H__
#define __INPUT_H__

typedef struct _Input
{
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

    int ncore;
    char *vasp_cmd;

    int kappa_dimer;
    double f_rot_min;
    double f_rot_max;
    int max_num_rot;
    int max_num_tls;

    int art_nouveau;
    double lambda_conv;
    int max_num_itr;
    int max_num_rlx;
    int delay_step;
    int mixing_step;
    int hyper_step;

    int random_seed;
} Input;

int input_int(int *var, char *tag, char *filename);
int input_double(double *var, char *tag, char *filename);
int input_char(char **var, char *tag, char *filename);
int input_char_arr(char ***var, char *tag, int n, char *filename);
int read_input(Input *input, char *filename);
void write_input(Input *input);
void free_input(Input *input);
#endif
