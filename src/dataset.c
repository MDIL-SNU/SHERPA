#include <stdlib.h>
#include "dataset.h"
#include "utils.h"


void insert_data(Dataset *dataset, int key_len, int *key,
                 int recycle_num, int *recycle_list,
                 int *type, double *initial, double *saddle, double *eigenmode)
{
    int i;
    Data *data = (Data *)malloc(sizeof(Data));
    data->recycle_num = recycle_num;
    data->key = (int *)malloc(sizeof(int) * key_len);
    data->type = (int *)malloc(sizeof(int) * recycle_num);
    data->initial = (double *)malloc(sizeof(double) * recycle_num * 3);
    data->saddle = (double *)malloc(sizeof(double) * recycle_num * 3);
    for (i = 0; i < natom; ++i) {
        data->type[i] = type[i];
        data->initial[i * 3 + 0] = initial[index[i] * 3 + 0];
        data->initial[i * 3 + 1] = initial[index[i] * 3 + 1];
        data->initial[i * 3 + 2] = initial[index[i] * 3 + 2];
        data->saddle[i * 3 + 0] = saddle[index[i] * 3 + 0];
        data->saddle[i * 3 + 1] = saddle[index[i] * 3 + 1];
        data->saddle[i * 3 + 2] = saddle[index[i] * 3 + 2];
    }
    if (dataset->head == NULL) {
        dataset->head = data;
    } else {
        data->next = dataset->head;
        dataset->head = data;
    }
    dataset->numdata++;
    if (dataset->numdata > 5000) {
        delete_data(dataset);
    }
}


void delete_data(Dataset *dataset) {
    Data *data = dataset->head;
    while (data->next->next != NULL) {
        data = data->next;
    }
    free(data->next->key)
    free(data->next->type);
    free(data->next->initial);
    free(data->next->saddle);
    free(data->next);
    data->next = NULL;
    dataset->numdata--; 
}


Data *search_data(Dataset *dataset, Config *config, long long key,
                  int *index, int numindex, int index)
{
    Data *ptr = dataset->head;
    while (1) {
        if (ptr == NULL) {
            return ptr;
        } else if ((ptr->key == key) && (numindex == ptr->numatom)) {
            int i;
            int numatom = 0;
            double vec[3];
            for (i = 0; i < numindex; ++i) {
                if (ptr->type[i] != config->type[index[i] - 1]) {
                    break;
                }
                vec[0] = ptr->initial[i * 3 + 0] - config->pos[(index[i] - 1) * 3 + 0];
                vec[1] = ptr->initial[i * 3 + 1] - config->pos[(index[i] - 1) * 3 + 1];
                vec[2] = ptr->initial[i * 3 + 2] - config->pos[(index[i] - 1) * 3 + 2];
                if ((norm(vec) > 1.0) && (index[i] - 1 != index)) {
                    break;
                }
                numatom++;
            }
            if (numatom == numindex) {
                return ptr;
            } else {
                ptr = ptr->next;
            }
        } else {
            ptr = ptr->next;
        }
    }
}


void free_dataset(Dataset *dataset)
{
    Data *ptr = dataset->head;
    if (ptr != NULL) {
        Data *next = ptr->next;
        while (1) {
            free(ptr->type);
            free(ptr->initial);
            free(ptr->saddle);
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
