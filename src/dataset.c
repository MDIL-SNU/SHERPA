#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dataset.h"
#include "my_mpi.h"
#include "sps_utils.h"


void insert_data(Dataset *dataset, int n, int index, double *eigenmode)
{
    int i;
    Data *data = (Data *)malloc(sizeof(Data));
    data->index = index;
    data->eigenmode = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
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
}


void build_dataset(Dataset *dataset, Input *input, int n)
{
    int i, j, count, index, errno;
    FILE *fp;
    struct dirent **namelist;

    int ncount = scandir(input->restart_dir, &namelist, name_filter, NULL); 
    if (ncount > 0) {
        for (i = 0; i < ncount; ++i) {
            /* load eigenmode */
            strtok(namelist[i]->d_name, "_");
            count = atoi(strtok(NULL, "_"));
            index = atoi(strtok(NULL, "."));
            MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&index, 1, MPI_INT, 0, MPI_COMM_WORLD);
            char filename[128];
            double *eigenmode = (double *)malloc(sizeof(double) * n * 3);
            sprintf(filename, "%s/%d.MODECAR",
                    input->restart_dir, count);
            fp = fopen(filename, "r");
            char line[1024];
            for (j = 0; j < n; ++j) {
                fgets(line, 1024, fp);
                eigenmode[j * 3 + 0] = atof(strtok(line, " \n"));
                eigenmode[j * 3 + 1] = atof(strtok(NULL, " \n"));
                eigenmode[j * 3 + 2] = atof(strtok(NULL, " \n"));
            }
            fclose(fp);
            /* insert data */
            insert_data(dataset, n, index, eigenmode);
            free(eigenmode);
        }
    }
    for (i = 0; i < ncount; ++i) {
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
            /* already freed in sps */
            //free(ptr->eigenmode);
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
