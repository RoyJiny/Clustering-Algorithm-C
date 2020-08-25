#ifndef UTILS_H
#define UTILS_H


#include "spmat.h"

char handle_errors(Error error, char *name);

/*input file -----> adjacency matrix
 *save the degree of each vertex in degree (by increasing order)*/
Error read_input(FILE *input, spmat *A, int *degree, int nof_vertex);

/*adjacency matrix -----> modularity matrix of group g*/
Error compute_modularity_matrix(spmat *A, int *g, int *degree, double M, spmat *B_g);

double compute_modularity_value(spmat *B_g, double *s);

Error power_iteration(spmat *mat, double *vector );

int calculate_eigen_value(spmat *mat, double *eigen_vector);/*TODO: substrat C_1norm*/

/*f_g is a vector of all f_i_g*/
Error compute_f_g(spmat* B_g , double* f_g);

/*modifying B_g to B_hat*/
Error convert2_B_hat (spmat* B_g , double* f_g); /*delta is not needed*/

Error matrix_shifting (spmat* B_hat ,double C_1norm);

Error matrix_deshifting (spmat* B_hat ,double C_1norm);

/*determine the partition into 2 groups, and calc the s vector*/
void eigen2s(double *eigen, int *g1 ,int *g2 , double* s, int g_size);

#endif