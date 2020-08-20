#include "spmat.h"


/*input file -----> adjacency matrix*/
errors read_input(FILE* input , spmat* A , int* degree ,int nof_vertex); /*TDOD: give a file. to support error in main*/
/*degree is a pointer of size nof_vertex*sizeof(int). (to save K_i)*/

/*computes C_prime and C_1norm*/
errors compute_C_prime (spmat* C , spmat* C_prime, double* C_1norm);

/*adjacency matrix -----> modularity matrix of group g*/
errors compute_modularity_matrix(spmat* A , int* g ,int* degree, int M, spmat* B_g);


void power_iteration(spmat* C_prime, double* vector, double epsilon ,double c_1norm);/*vector will be the eigen vector*/
/*TODO: 1.use IS_POSITIVE
*       2.save the previous to the final eigen vector (b_k-1)*/

double dot_product(double* row1, double* row2, int size);

/*s[i]=1 if eigen[i]>0 else s[i]=-1*/
void eigen2s (double* eigen , int* s ,int n);

/*computes eigen value*/
double get_eigen_value (double* eigen , double* prev_eigen, int size);

/*f_g is a vector of all f_i_g*/
errors compute_f_g(spmat* B_g , double* f_g);

/*B_g ------> B_hat_g*/
errors compute_B_hat (spmat* B_g , double* f_g , spmat* B_hat); /*delta is not needed*/


void modularity_maximization (spmat* B_g , int* g1 , int* g2);



typedef enum {NONE , ALLOCATION_FAILED ,READ_FAILED} errors;