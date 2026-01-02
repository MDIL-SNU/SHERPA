#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "files.h"


int main(int argc, char *argv[])
{
    int i;
    long long step = 0;
    long long initial_step = 0;
    double Ea = 0.0;
    double kmc_time = 0.0;
    char directory[1024], filename1[1024], filename2[1024];
    FILE *fp, *pp;

    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        printf("Usage: KMC --att_freq {att_freq} --temperature {temperature} --kmc_step {kmc_step} --low_cut {low_cut} --restart {restart} --inputs_path {inputs_path} --sherpa_cmd {sherpa_cmd} --random_seed {random_seed}\n");
        return 0;
    }
    double att_freq = 1e13;
    double temperature = 298;
    long long kmc_step = 10000;
    double low_cut = 0.0;
    long long restart = 0;
    char *sherpa_cmd = NULL;
    char inputs_path[1024] = "./INPUTS";
    sprintf(inputs_path, "./INPUTS");
    unsigned int random_seed = (unsigned int)time(NULL);
    /* read input */
    for (i = 1; i < argc; i = i + 2) {
        if (strcmp(argv[i], "--att_freq") == 0) {
            att_freq = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "--temperature") == 0) {
            temperature = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "--kmc_step") == 0) {
            kmc_step = atoll(argv[i + 1]);
        }
        if (strcmp(argv[i], "--low_cut") == 0) {
            low_cut = atof(argv[i + 1]);
        }
        if (strcmp(argv[i], "--restart") == 0) {
            restart = atoll(argv[i + 1]);
        }
        if (strcmp(argv[i], "--inputs_path") == 0) {
            strncpy(inputs_path, argv[i + 1], sizeof(inputs_path) - 1);
        }
        if (strcmp(argv[i], "--sherpa_cmd") == 0) {
            sherpa_cmd = argv[i + 1];
        }
        if (strcmp(argv[i], "--random_seed") == 0) {
            random_seed = atoi(argv[i + 1]);
        }
    }
    if (sherpa_cmd == NULL) {
        printf("--sherpa_cmd is missing.\n");
        return 1;
    }
    srand(random_seed);

    if (restart == 0) {
        /* write log */
        fp = fopen("./KMC.log", "w");
        fputs("-----------------------------------------------------\n", fp);
        fputs("   KMC steps      Barrier energy      Elapsed time   \n", fp);
        fputs("-----------------------------------------------------\n", fp);
        fputs("           0      --------------               0.0   \n", fp);
        fclose(fp);
        sprintf(filename1, "%s/POSCAR", inputs_path);
        copy_files("./KMC.XDATCAR", filename1);
        kmc_time = 0.0;
        initial_step = 1;
    } else {
        /* read log */
        fp = fopen("./KMC.log", "r");
        FILE *tmp_fp = fopen("./KMC_tmp.log", "w");
        char *ptr;
        char line[1024], tmp_line[1024];
        ptr = fgets(tmp_line, 1024, fp);
        while (ptr != NULL) {
            memcpy(line, ptr, sizeof(char) * (strlen(ptr) + 1));
            step = atoll(strtok(tmp_line, " \n"));
            fputs(line, tmp_fp);
            if (step == restart) {
                strtok(NULL, " \n");
                kmc_time = atof(strtok(NULL, " \n"));
                break;
            }
            ptr = fgets(tmp_line, 1024, fp);
        }
        fclose(fp);
        fclose(tmp_fp);
        if (step != restart) {
            printf("Check restart step.\n");
            return 1;
        }
        rename("./KMC.log", "./KMC_old.log");
        rename("./KMC_tmp.log", "./KMC.log");
        initial_step = step + 1;
        /* copy POSCAR */
        DIR *dp = NULL;
        struct dirent *entry = NULL;
        sprintf(directory, "%lld", step);
        if ((dp = opendir(directory)) == NULL) {
            printf("Check restart step.\n");
            return 1;
        }
        while ((entry = readdir(dp)) != NULL) {
            if (strncmp(entry->d_name, "Final_", 6) == 0) {
                sprintf(filename1, "%lld/%s", step, entry->d_name);
                break;
            }
        }
        sprintf(filename2, "%s/POSCAR", inputs_path);
        copy_files(filename2, filename1);
    }

    for (step = initial_step; step < initial_step + kmc_step; ++step) {
        /* mkdir */
        sprintf(directory, "%lld", step);
        mkdir(directory, 0775);
        /* cp INPUTS */
        DIR *direct;
        struct dirent *entry;
        direct = opendir(inputs_path);
        if (direct == NULL) {
            printf("Check inputs_path.\n");
            return 1;
        }
        while ((entry = readdir(direct)) != NULL) {
            /* only files */
            if (entry->d_type == DT_REG) {
                sprintf(filename1, "%s/%s", inputs_path, entry->d_name);
                sprintf(filename2, "%lld/%s", step, entry->d_name);
                copy_files(filename2, filename1);
            }
        }
        /* chdir */
        chdir(directory);
        int list_size = 1024;
        int reac_num = 0;
        int *reac_list = (int *)malloc(sizeof(int) * list_size);
        double *acti_list = (double *)malloc(sizeof(double) * list_size);
        double *rate_list = (double *)malloc(sizeof(double) * list_size);
        do {
            /* run sherpa */
            pp = popen(sherpa_cmd, "r");
            if (pp != NULL) {
                while (1) {
                    fp = fopen("./Time.log", "r");
                    if (fp != NULL) {
                        printf("Time found\n");
                        fclose(fp);
                        break;
                    }
                }
                pclose(pp);
            } else {
                printf("Please provide absolute SHERPA path.\n");
                return 1;
            }
            /* read log */
            fp = fopen("./Event.log", "r");
            if (fp == NULL) {
                printf("Check your SHERPA.");
                return 1;
            }
            char line[1024];
            fgets(line, 1024, fp);
            fgets(line, 1024, fp);
            fgets(line, 1024, fp);
            double kT = 8.61733326E-5 * temperature;
            while (fgets(line, 1024, fp) != NULL) {
                int tmp_index = atoi(strtok(line, " \n"));
                Ea = atof(strtok(NULL, " \n"));
                if (Ea < low_cut) {
                    continue;
                }
                reac_list[reac_num] = tmp_index;
                acti_list[reac_num] = Ea;
                rate_list[reac_num] = att_freq * exp(-Ea / kT);
                reac_num++;
                if (reac_num >= list_size) {
                    list_size = list_size << 1;
                    reac_list = (int *)realloc(reac_list, sizeof(int) * list_size);
                    acti_list = (double *)realloc(acti_list, sizeof(double) * list_size);
                    rate_list = (double *)realloc(rate_list, sizeof(double) * list_size);
                }
            }
            fclose(fp);
        } while (reac_num == 0);
        /* select event */
        double rate_sum = 0.0;
        for (i = 0; i < reac_num; ++i) {
            rate_sum += rate_list[i];
        }
        int count = 0;
        double acc_rate_sum = 0.0;
        double rn = (double)random() / RAND_MAX;
        for (i = 0; i < reac_num; ++i) {
            rate_list[i] /= rate_sum;
            acc_rate_sum += rate_list[i];
            if (acc_rate_sum > rn) {
                count = reac_list[i];
                Ea = acti_list[i];
                kmc_time += -log((double)random() / RAND_MAX) / rate_sum;
                break;
            }
        }
        free(reac_list);
        free(acti_list);
        free(rate_list);
        /* extract */
        sprintf(filename1, "Final.POSCAR");
        extract_files(filename1, count);
        sprintf(filename1, "./Final_%d.POSCAR", count); 
        sprintf(filename2, "../%s/POSCAR", inputs_path);
        copy_files(filename2, filename1);
        /* log */
        fp = fopen("../KMC.log", "a");
        fprintf(fp, "   %9lld      %14.3f      %12e\n", step, Ea, kmc_time);
        fclose(fp);
        append_files("../KMC.XDATCAR", filename2);
        /* chdir */
        chdir("../");
    }
    return 0;
}
