#include "dataset.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void insert_data(Dataset *dataset, int n, int index, double *eigenmode)
{
    int i;
    Data *data = (Data *)malloc(sizeof(Data));
    data->next = NULL;
    data->index = index;
    data->eigenmode = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        data->eigenmode[i * 3 + 0] = eigenmode[i * 3 + 0];
        data->eigenmode[i * 3 + 1] = eigenmode[i * 3 + 1];
        data->eigenmode[i * 3 + 2] = eigenmode[i * 3 + 2];
    }
    if ((dataset->head == NULL) && (dataset->tail == NULL)) {
        dataset->head = dataset->tail = data;
    } else {
        dataset->tail->next = data;
        dataset->tail = data;
    }
}


void build_dataset(Dataset *dataset, int n)
{
    int i, index;
    FILE *fp = fopen("./Initial.MODECAR", "r");
    char line[1024];
    while (fgets(line, 1024, fp) != NULL) {
        strtok(line, "_\n");
        index = atoi(strtok(NULL, "_\n"));
        double *eigenmode = (double *)malloc(sizeof(double) * n * 3);
        for (i = 0; i < n; ++i) {
            fgets(line, 1024, fp);
            eigenmode[i * 3 + 0] = atof(strtok(line, " \n"));
            eigenmode[i * 3 + 1] = atof(strtok(NULL, " \n"));
            eigenmode[i * 3 + 2] = atof(strtok(NULL, " \n"));
        }
        insert_data(dataset, n, index, eigenmode);
        free(eigenmode);
    }
    fclose(fp);
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
