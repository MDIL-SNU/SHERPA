#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dataset.h"
#include "utils.h"


void insert_data(Dataset *dataset, int n, int index,
                 int *type, double *saddle, double *eigenmode)
{
    int i;
    Data *data = (Data *)malloc(sizeof(Data));
    data->index = index;
    data->type = (int *)malloc(sizeof(int) * n);
    data->saddle = (double *)malloc(sizeof(double) * n * 3);
    data->eigenmode = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        data->type[i] = type[i];
        data->saddle[i * 3 + 0] = saddle[i * 3 + 0];
        data->saddle[i * 3 + 1] = saddle[i * 3 + 1];
        data->saddle[i * 3 + 2] = saddle[i * 3 + 2];
        data->eigenmode[i * 3 + 0] = eigenmode[i * 3 + 0];
        data->eigenmode[i * 3 + 1] = eigenmode[i * 3 + 1];
        data->eigenmode[i * 3 + 2] = eigenmode[i * 3 + 2];
    }
    if (dataset->head == NULL) {
        dataset->head = data;
        dataset->head->next = NULL;
    } else {
        data->next = dataset->head;
        dataset->head = data;
    }
    dataset->numdata++;
}


void build_dataset(Dataset *dataset, Config *config, Input *input)
{
    int i, j, errno;
    double del[3];
    char *ptr;
    FILE *fp;
    struct dirent **namelist;

    int count = scandir(input->output_dir, &namelist, name_filter, NULL); 
    if (count > 0) {
        for (i = 0; i < count; ++i) {
            /* load saddle */
            char filename[128];
            strtok(namelist[i]->d_name, "_");
            ptr = strtok(NULL, ".");
            sprintf(filename, "%s/Saddle_%s.POSCAR", input->output_dir, ptr);
            Config *tmp_config = (Config *)malloc(sizeof(Config));
            errno = read_config(tmp_config, input, filename);
            if (errno > 0) {
                printf("Cannot find %s\n", filename);
            }
            /* load eigenmode */
            int n = config->tot_num;
            double max_dist = 0.0;
            int index = 0;
            for (j = 0; j < n; ++j) {
                del[0] = config->pos[j * 3 + 0] - tmp_config->pos[j * 3 + 0];
                del[1] = config->pos[j * 3 + 1] - tmp_config->pos[j * 3 + 1];
                del[2] = config->pos[j * 3 + 2] - tmp_config->pos[j * 3 + 2];
                if ((abs(del[0]) < max_dist) || (abs(del[1]) < max_dist) || (abs(del[2]) < max_dist)) {
                    continue;
                }
                get_minimum_image(del, config->boxlo, config->boxhi,
                                  config->xy, config->yz, config->xz);
                double dist = sqrt(del[0] * del[0]
                                 + del[1] * del[1] 
                                 + del[2] * del[2]);
                if (dist > max_dist) {
                    index = j;
                    max_dist = dist;
                }                 
            }
            double *eigenmode = (double *)malloc(sizeof(double) * n * 3);
            sprintf(filename, "%s/%s.MODECAR", input->output_dir, ptr);
            fp = fopen(filename, "rb");
            int why = fread(eigenmode, sizeof(double), n * 3, fp);
            fclose(fp);
            /* insert data */
            insert_data(dataset, n, index, tmp_config->type, tmp_config->pos, eigenmode);
            free(tmp_config);
            free(eigenmode);
        }
    }
    for (i = 0; i < count; ++i) {
        free(namelist[i]);
    }
    free(namelist);
}


void free_dataset(Dataset *dataset)
{
    Data *ptr = dataset->head;
    if (ptr != NULL) {
        Data *next = ptr->next;
        while (1) {
            free(ptr->type);
            free(ptr->saddle);
            free(ptr->eigenmode);
            free(ptr);
            if (next == NULL) {
                break;
            }
            ptr = next;
            next = ptr->next;
        }
    }
    free(dataset);
}
