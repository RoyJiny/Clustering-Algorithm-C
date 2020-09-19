#ifndef UTILS_H
#define UTILS_H

#include "spmat.h"
#include <math.h>

void print_vector_int(int *vector, int size);

void print_vector(double *vector, int size);

void print_group(group *g);

void print_stack(group_set *s);

/*computes the 1norm of the entire 'A' matrix, the results stored in res.*/
void compute_1norm(spmat *A, group *g, int *degrees, double M, double *res);

/*return the dot product of row1 with row2.
 *assuming both of the vectors are at least of size 'size.'*/
double dot_product(double *row1, double *row2, int size);

/*convert the information in 'input' file to an adjacency matrix stored in'A'.
 *save the degree of each vertex in 'degree' (by increasing order).*/
void read_input(FILE *input, spmat *A, int *degree, int nof_vertex);

/*computes row 'row' of the modularity matrix corresponding to group g.
 *the result is stored in 'B_g_row'.*/
void compute_modularity_matrix_row(spmat *A, int row, group *g, int *degrees, double M, double *B_g_row, int g_count);

double compute_D_row(int A_row, group *g, int *degrees, double M, double *B_g_row);

/*computes the modularity value according to the partition represented in 's'.
 *the result is stored in 'res'*/
void compute_modularity_value(spmat *A, group *g, int *degrees, double *s, double M, double *B_g_row, double *mult_vector, double *res);

/*computes the difference between the initial modularity and the modularity. 
 *we would get if we move the vertex 'A_index' to the other group.
 *the result is stored in 'score'.*/
void compute_score(spmat *A, int A_index, int g_index, group *g, double *s, double M, int *degrees, double *score, double *B_g_row);

/*computes the eigen value ,corresponding to 'eigen_vector', of matrix 'mat'.
 *the result is stored in 'res'.*/
void calculate_eigen_value(spmat *mat, double *eigen_vector, group *g, int *degrees, double M, double *B_g_row, double B_1norm, double *res);

/*computes the partition 's' according to 'eigen'.
 *the result is stored in 's'.
 *the function returns the number of vertex that will be in the first group (or g1) after the partition.*/
int eigen2s(double *eigen, group *g, double *s);

/*update the 2 groups('g1','g2') according to the partition of 'g' represented in 's'.*/
void construct_g1g2(group *g, double *s, group *g1, group *g2, int size);

/*write to 'output' file the final partition into groups.*/
void write_output_file(FILE *output, group_set *O);

void print_output(FILE *output, int nof_vertex);

#endif