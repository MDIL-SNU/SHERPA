#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "alg_utils.h"


void int_shuffle(int list_len, int *int_list)
{
    int i, tmp;
    for (i = 0; i < list_len - 1; ++i) {
        int random = i + rand() % (list_len - i);
        tmp = int_list[i];
        int_list[i] = int_list[random];
        int_list[random] = tmp;
    }
}


void int_sort_decrease(int list_len, int *int_list)
{
    int i;
    while (1) {
        int done = 1;
        for (i = 1; i < list_len; ++i) {
            if (int_list[i - 1] < int_list[i]) {
                int tmp_int = int_list[i - 1];
                int_list[i - 1] = int_list[i];
                int_list[i] = tmp_int;
                done = 0;
            }
        }
        if (done) {
            break;
        }
    }
}


void int_sort_increase(int list_len, int *int_list)
{
    int i;
    while (1) {
        int done = 1;
        for (i = 1; i < list_len; ++i) {
            if (int_list[i - 1] > int_list[i]) {
                int tmp_int = int_list[i - 1];
                int_list[i - 1] = int_list[i];
                int_list[i] = tmp_int;
                done = 0;
            }
        }
        if (done) {
            break;
        }
    }
}


double normal_random(double mean, double std)
{
    double u, v, s;
    do {
        u = ((double)rand() / RAND_MAX) * 2 - 1;
        v = ((double)rand() / RAND_MAX) * 2 - 1;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    s = sqrt(-2 * log(s) / s);
    return mean + std * u * s;
}


double norm(double *vec, int n)
{
    int i;
    double output = 0.0;
    for (i = 0; i < n; ++i) {
        output += vec[i * 3 + 0] * vec[i * 3 + 0];
        output += vec[i * 3 + 1] * vec[i * 3 + 1];
        output += vec[i * 3 + 2] * vec[i * 3 + 2];
    }
    return sqrt(output);
}


double *normalize(double *vec, int n)
{
    int i;
    double *output = (double *)malloc(sizeof(double) * n * 3); 
    double magnitude = norm(vec, n);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = vec[i * 3 + 0] / magnitude;
        output[i * 3 + 1] = vec[i * 3 + 1] / magnitude;
        output[i * 3 + 2] = vec[i * 3 + 2] / magnitude;
    }
    return output;
}


double dot(double *vec1, double *vec2, int n)
{
    int i;
    double output = 0.0;
    for (i = 0; i < n; ++i) {
        output += vec1[i * 3 + 0] * vec2[i * 3 + 0];
        output += vec1[i * 3 + 1] * vec2[i * 3 + 1];
        output += vec1[i * 3 + 2] * vec2[i * 3 + 2];
    }
    return output;
}


void cross(double *vec1, double *vec2, double *vec3)
{
    vec3[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vec3[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vec3[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}


double det(double (*mat)[3])
{
    int i;
    double cofactor[3];
    cofactor[0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    cofactor[1] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    cofactor[2] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];

    double det = 0.0;
    for (i = 0; i < 3; i++) {
        det += cofactor[i] * mat[0][i];
    }
    return det;
}


double *parallel_vector(double *vec, double *unit, int n)
{
    int i; 
    double magnitude = dot(vec, unit, n);
    double *output = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = magnitude * unit[i * 3 + 0];
        output[i * 3 + 1] = magnitude * unit[i * 3 + 1];
        output[i * 3 + 2] = magnitude * unit[i * 3 + 2];
    } 
    return output;
}


double *perpendicular_vector(double *vec, double *unit, int n)
{
    int i;
    double *tmp_vec = parallel_vector(vec, unit, n);
    double *output = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        output[i * 3 + 0] = vec[i * 3 + 0] - tmp_vec[i * 3 + 0];
        output[i * 3 + 1] = vec[i * 3 + 1] - tmp_vec[i * 3 + 1];
        output[i * 3 + 2] = vec[i * 3 + 2] - tmp_vec[i * 3 + 2];
    } 
    free(tmp_vec);
    return output;
}


void rotate_vector(double *vec1i, double *vec2i,
                   double **vec1o, double **vec2o,
                   int n, double angle)
{
    int i;
    double cAng = cos(angle);
    double sAng = sin(angle);
    double *tmp_vec1o = (double *)malloc(sizeof(double) * n * 3);
    double *tmp_vec2o = (double *)malloc(sizeof(double) * n * 3);
    for (i = 0; i < n; ++i) {
        tmp_vec1o[i * 3 + 0] = vec1i[i * 3 + 0] * cAng + vec2i[i * 3 + 0] * sAng;
        tmp_vec1o[i * 3 + 1] = vec1i[i * 3 + 1] * cAng + vec2i[i * 3 + 1] * sAng;
        tmp_vec1o[i * 3 + 2] = vec1i[i * 3 + 2] * cAng + vec2i[i * 3 + 2] * sAng;
        tmp_vec2o[i * 3 + 0] = vec2i[i * 3 + 0] * cAng - vec1i[i * 3 + 0] * sAng;
        tmp_vec2o[i * 3 + 1] = vec2i[i * 3 + 1] * cAng - vec1i[i * 3 + 1] * sAng;
        tmp_vec2o[i * 3 + 2] = vec2i[i * 3 + 2] * cAng - vec1i[i * 3 + 2] * sAng;
    }
    double magnitude1 = norm(vec1i, n);
    double magnitude2 = norm(vec2i, n);
    *vec1o = normalize(tmp_vec1o, n);
    *vec2o = normalize(tmp_vec2o, n);
    for (i = 0; i < n; ++i) {
        (*vec1o)[i * 3 + 0] *= magnitude1;
        (*vec1o)[i * 3 + 1] *= magnitude1;
        (*vec1o)[i * 3 + 2] *= magnitude1;
        (*vec2o)[i * 3 + 0] *= magnitude2;
        (*vec2o)[i * 3 + 1] *= magnitude2;
        (*vec2o)[i * 3 + 2] *= magnitude2;
    }
    free(tmp_vec1o);
    free(tmp_vec2o);
}
