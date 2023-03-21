#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{
    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        printf("Usage: KMC {OLD_DIRECTORY} {NEW_DIRECTORY} {TEMPERATURE} {ATTEMPT_FREQUENCY}\n");
        printf("Output: new INPUT file and time.log\n");
        return 0;
    }
    char filename[1024];
    sprintf(filename, "%s/Event.log", argv[1]);
    FILE *rp = fopen(filename, "r");
    if (rp == NULL) {
        printf("Check your directory.");
        return 1;
    }

    char line[1024], tmp_line[1024], *ptr;
    double Ea;
    /* skip header */
    fgets(line, 1024, rp);
    fgets(line, 1024, rp);
    fgets(line, 1024, rp);
    /* read Event.log */
    fgets(line, 1024, rp);
    int list_size = 64;
    int reac_num = 0;
    int *reac_list = (int *)malloc(sizeof(int) * list_size);
    double *acti_list = (double *)malloc(sizeof(double) * list_size);
    double *rate_list = (double *)malloc(sizeof(double) * list_size);
    double kT = 8.61733326E-5 * atof(argv[3]);
    while (strchr(line, '-') == NULL) {
        reac_list[reac_num] = atoi(strtok(line, " \n"));
        Ea = atof(strtok(NULL, " \n"));
        acti_list[reac_num] = Ea;
        rate_list[reac_num] = atof(argv[4]) * exp(-Ea / kT);
        reac_num++;
        if (reac_num >= list_size) {
            list_size = list_size << 1;
            reac_list = (int *)realloc(reac_list, sizeof(int) * list_size);
            acti_list = (double *)realloc(acti_list, sizeof(double) * list_size);
            rate_list = (double *)realloc(rate_list, sizeof(double) * list_size);
        }
        fgets(line, 1024, rp);
    }
    fclose(rp);
    /* reaction index */
    int i;
    double rate_sum = 0.0;
    for (i = 0; i < reac_num; ++i) {
        rate_sum += rate_list[i];
    }
    int index, tmp_index;
    double acc_rate_sum = 0.0;
    double rn = (double)random() / RAND_MAX;
    for (i = 0; i < reac_num; ++i) {
        rate_list[i] /= rate_sum;
        acc_rate_sum += rate_list[i];
        if (acc_rate_sum > rn) {
            index = reac_list[i];    
            Ea = acti_list[i];
            break;
        }
    }
    /* time log */
    double dt = -log((double)random() / RAND_MAX) / rate_sum;
    rp = fopen("./time.log", "r");
    if (rp == NULL) {
        rp = fopen("./time.log", "w");
        fputs("--------------------------------------------------\n", rp);
        fputs(" kMC step   Index   Barrier energy   kMC time (s)\n", rp);
        fputs("--------------------------------------------------\n", rp);
        fputs("        0   -----   --------------   0.0000000000\n", rp);
    }
    fclose(rp);
    rp = fopen("./time.log", "a+");
    while (fgets(line, 1024, rp) != NULL) {
        strcpy(tmp_line, line);
    }
    int step = atoi(strtok(tmp_line, " \n"));
    strtok(NULL, " \n");
    ptr = strtok(NULL, " \n");
    fprintf(rp, " %8d   %5d   %14f   %12e\n",
            step + 1, index, Ea, atof(strtok(NULL, " \n")) + dt);
    fclose(rp);
    /* next configuration */
    sprintf(filename, "%s/Final.POSCAR", argv[1]);
    rp = fopen(filename, "r");
    if (rp == NULL) {
        printf("Check your directory.");
        return 1;
    }
    FILE *wp = fopen("./tmp_POSCAR", "w");
    int flag = 0;
    while (fgets(tmp_line, 1024, rp) != NULL) {
        strcpy(line, tmp_line);
        if ((strchr(line, '_')) != NULL) {
            if (atoi(strtok(tmp_line, "_")) == index) {
               flag = 1; 
            } else {
                if (flag > 0) {
                    break;
                }
            }
        }
        if (flag > 0) {
            fputs(line, wp);
        }
    }
    fclose(rp);
    fclose(wp);

    free(reac_list);
    free(acti_list);
    free(rate_list);
    return 0;
}
