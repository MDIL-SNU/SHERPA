#ifndef __ALG_UTILS_H__
#define __ALG_UTILS_H__
#include "config.h"
#include "input.h"


void int_sort_decrease(int list_len, int *int_list);
void int_sort_increase(int list_len, int *int_list);
double normal_random(double mean, double std);
double norm(double *vec, int n);
double *normalize(double *vec, int n);
double dot(double *vec1, double *vec2, int n);
void cross(double *vec1, double *vec2, double *vec3);
double det(double (*mat)[3]);
double *parallel_vector(double *vec, double *unit, int n);
double *perpendicular_vector(double *vec, double *unit, int n);
void rotate_vector(double *vec1i, double *vec2i,
                   double **vec10, double **vec20,
                   int n, double angle);
#endif
