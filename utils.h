#ifndef UTILS_H
#define UTILS_H

#include "spmat.h"
#include <math.h>

void print_vector_int(int *vector, int size);

void print_vector(double *vector, int size);

char handle_errors(Error error, char *name);

double dot_product(double *row1, double *row2, int size);

/*input file -----> adjacency matrix
 *save the degree of each vertex in degree (by increasing order)*/
Error read_input(FILE *input, spmat *A, int *degree, int nof_vertex);

/**/
Error compute_modularity_matrix_row(spmat *A, int row, group *g, int *degrees, double M, double *B_g_row);

double compute_modularity_value(spmat *B_g, double *s);

Error power_iteration(spmat *mat, double *vector);

double calculate_eigen_value(spmat *mat, double *eigen_vector, group *g, int *degrees, double M, double *B_g_row, double B_1norm);

/*determine the partition into 2 groups, and calc the s vector*/
void eigen2s(double *eigen, group *g1, group *g2, double *s, int size);

#endif