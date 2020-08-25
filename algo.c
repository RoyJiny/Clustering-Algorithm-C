#include <stdio.h>
#include <stdlib.h>

#include "algo.h"
/*TODO: need to do computations each row at a time ):*/
Error algo_2(spmat *A, int* degrees, double* init_vector, int g_size, int* g, int* g1, int* g2)
{
    int *temp_i;
    double M , C_1norm ,modularity_value ,eigen_value ,*s;
    Error error;
    spmat *B_g;
    double *f_g ,*eigen_vector = init_vector;
    
    /*------------------------ALLOCATIONS------------------------*/
    B_g = spmat_allocate_list(g_size);
    if(!B_g){
        return ALLOCATION_FAILED;
    }

    f_g = (double*) malloc(g_size*sizeof(double));
    if(!f_g){
        return ALLOCATION_FAILED;
    }

    s = (double*) malloc(g_size*sizeof(double));
    if(!s){
        return ALLOCATION_FAILED;
    }

    /*---------------------computing M----------------------------*/
	for (temp_i = degrees; temp_i < degrees + A->n; temp_i++)
	{
		M += *temp_i;
	}

    /*---------------------compute B_g from A----------------------*/
	error = compute_modularity_matrix(A, g, degrees, M, B_g);
    if(error != NONE){
        printf("failed in compute_modularity_matrix\n");
        return error;
    }

    /*------------compute f_g (needed to compute B_hat)-------------*/
	error = compute_f_g(B_g ,f_g);
	if(error != NONE){
        printf("failed in compute_f_g\n");
        return error;
    }

    /*---------------------convert B_g to B_hat---------------------*/
    error = convert2_B_hat(B_g, f_g);/*includ matrix shifting*/
	if(error != NONE){
        printf("failed in convert2_B_hat\n");
        return error;
    }

    /*-------------------------matrix shifting----------------------*/
    C_1norm = B_g->compute_1norm(B_g);
    error = matrix_shifting(B_g , C_1norm);
	if(error != NONE){
        printf("failed in matrix_shifting\n");
        return error;
    }

    /*-----------------computing leading eigen pair---------*/
    error = power_iteration(B_g, eigen_vector);
    if(error != NONE){
        printf("failed in power_iteration\n");
        return error;
    }
    eigen_value = calculate_eigen_value(B_g, eigen_vector);
    /*TODO: something in the nirmool (see forum), and need to sub C_1norm*/

    /*------------------------matrix deshifting----------------------*/
    error = matrix_deshifting(B_g , C_1norm); /*check if it working, it should be OK*/
	if(error != NONE){
        printf("failed in matrix_deshifting\n");
        return error;
    }

    /*--------------------decide the right partition-----------------*/
    if (!IS_POSITIVE(eigen_value)){ /*g is indivisible*/
        /*TODO : update g1 & g2*/
    }
    else{
        eigen2s(eigen_vector, g1, g2, s, g_size);
        modularity_value = compute_modularity_value(B_g, s);
        if(!IS_POSITIVE(modularity_value)){ /*g is indivisible*/
             /*TODO : update g1 & g2*/
        }
        /*else: g is divisible and the correct g1 & g2 are already computed in eigen2s*/
    }


free(f_g);
free(s);
B_g->free(B_g);
return NONE;
}


