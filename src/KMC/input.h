#ifndef __INPUT_H__
#define __INPUT_H__

typedef struct _Input
{
    int kmc_step;
    double temperature;
    double att_freq;

    char *extractor;

    char *sherpa_cmd;
    int num_inputs;
    char **inputs;

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
