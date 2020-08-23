#include "spmat.h"
typedef enum
{
    NONE,
    ALLOCATION_FAILED,
    READ_FAILED
} Error;

char handle_errors(Error error, char *name);

/*input file -----> adjacency matrix
 *save the degree of each vertex in degree (by increasing order)*/
Error read_input(FILE *input, spmat *A, int *degree, int nof_vertex);

/*adjacency matrix -----> modularity matrix of group g*/
Error compute_modularity_matrix(spmat *A, double *g, int *degree, double M, spmat *B_g);

double compute_modularity_value(spmat *B_g, int *g);