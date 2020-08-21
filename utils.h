#include "spmat.h"

void read_input(char *path, spmat *A);

double compute_C_prime(spmat *C, spmat *C_prime);

void compute_modularity_matrix(spmat *A, int *g, spmat *B_g); /*TODO: replace name*/

void power_iteration(spmat *C_prime, double *vector, double epsilon, double c_1norm); /*vector will be the eigen vector*/

double dot_product(double *row1, double *row2, int size);

void eigen2s(double *eigen, int *s); /*TODO: replace name*/

double get_eigen_value(spmat *A, double *eigen);

void compute_f_g(spmat *B_g, double *f_g);

void compute_B_hat(spmat *B_g, double *f_g, spmat *B_hat);

int compute_two_groups_modularity(spmat *B, int *s); /* 0.5*(s^T)*B*s */

void modularity_maximization(spmat *B_g, int *g);
/*TODO: maybe add enum of errors*/