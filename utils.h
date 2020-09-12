#ifndef UTILS_H
#define UTILS_H

#include "spmat.h"
#include <math.h>

void print_vector_int(int *vector, int size);

void print_vector(double *vector, int size);

void print_group(group *g);

void print_stack(group_set *s);

Error compute_1norm(spmat *A, group *g, int *degrees, double M, double *res);

void print_errors(Error error, char *name, char *func);

double dot_product(double *row1, double *row2, int size);

/*input file -----> adjacency matrix
 *save the degree of each vertex in degree (by increasing order)*/
Error read_input(FILE *input, spmat *A, int *degree, int nof_vertex);

/**/
Error compute_modularity_matrix_row(spmat *A, int row, group *g, int *degrees, double M, double *B_g_row, int g_count);

Error compute_modularity_value(spmat *A, group *g, int *degrees, double *s, double M, double *B_g_row, double *mult_vector, double *res);

Error compute_for_improved_score(spmat *A, int A_index, int g_index, group *g, double *s, double M, int *degrees, double *score, double *B_g_row);

Error power_iteration(spmat *mat, double *vector);

Error calculate_eigen_value(spmat *mat, double *eigen_vector, group *g, int *degrees, double M, double *B_g_row, double B_1norm, double *res);

/*determine the partition into 2 groups, and calc the s vector*/
int eigen2s(double *eigen, group *g, double *s);

Error construct_g1g2(group *g, double *s, group *g1, group *g2, int size);

Error write2_output_file(FILE *output, group_set *O);

Error print_output(FILE *output, int nof_vertex);

#endif