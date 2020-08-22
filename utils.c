#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "utils.h"

double dot_product(double *row1, double *row2, int size)
{
	double sum = 0;
	double *end = row1 + size;
	for (; row1 < end;)
	{
		sum += (*row1) * (*row2);
		row1++;
		row2++;
	}

	return sum;
}

void power_iteration(spmat *A, double *vector, double epsilon, double c_1norm)
{
	int stop = 0;
	double magnitude;
	double *mul_result;
	double *ptr_v, *ptr_r; /*pointers to vector and mul_result*/
	mul_result = (double *)malloc((A->n) * sizeof(double));
	if (!mul_result)
	{
		return;
	}
	while (stop == 0)
	{

		A->mult(A, vector, mul_result); /*result=A*vector*/
		magnitude = sqrt(dot_product(mul_result, mul_result, A->n));

		stop = 1;
		ptr_r = mul_result;
		for (ptr_v = vector; ptr_v < vector + A->n; ptr_v++)
		{
			if (fabs(*ptr_v - ((*(ptr_r)) / magnitude)) > epsilon)
			{
				stop = 0;
			}
			*ptr_v = (*(ptr_r)) / magnitude; /*updating vector to the next vector*/
			ptr_r++;
		}
	}
	free(mul_result);
}

void modularity_maximization(spmat *B_g, int *g)
{
	/*g is a single vector that can represnt g1,g2 - devision to 2 groups*/
	int size = B_g->n;
	/*save a linked list of all verticies for the beginning of every iteration*/

	/*while modularity of the best state keeps increasing:*/

	/*get a copy of the verticies list*/
	/*save the initial modularity*/
	/*create the new g1,g2 to represent the current state*/
	/*(at first the original g1,g2, after that - the best state from the previous process)*/
	/*create the new g1,g2 to represent the new best state (same as above)*/

	/*while verticies list is not empty: (to make sure that each vertex moves once)*/

	/*create variables for the best move: vertex, modularity value, division(g1,g2)*/

	/*while we havent reached the end of the list*/

	/*move the current vertex to the opposite group*/
	/*modify g1,g2*/
	/*calculate the new modularity*/

	/*if the new modularity is better than the highest so far:*/
	/*save the new modularity value as the highest so far this iteration*/
	/*save the new g1,g2 as the highest g1,g2*/
	/*save the respective vertex as part of the best move*/

	/*compare with best state overall, and update if better*/

	/*move the head of verticies list forward*/
}
/*remove the vertex of the best move from the list*/
/*assume that the number of vertexes as already been read from input file*/
errors read_input(FILE *input, spmat *A, int *degree, int nof_vertex)
{
	int *curr_row;
	int *start_row; /*save the start of the row to add_row function in spmat*/
	int *temp;
	int *start_temp; /*save the start of temp for reuse*/
	int *curr_vertex;
	int i, j;

	curr_vertex = degree; /*fill degree*/
	start_row = (int *)malloc(nof_vertex * sizeof(int));
	if (!start_row)
	{
		return ALLOCATION_FAILED;
	}
	start_temp = (int *)malloc(nof_vertex * sizeof(int));
	if (!start_temp)
	{
		return ALLOCATION_FAILED;
	}

	for (i = 0; i < nof_vertex; i++)
	{
		curr_row = start_row;
		temp = start_temp;
		/*try read k_i*/
		if (fread(curr_vertex, sizeof(int), 1, input) != 1)
		{
			return READ_FAILED;
		}
		/*try read the k_i neighbors*/
		if (fread(temp, sizeof(int), *curr_vertex, input) != *curr_vertex)
		{
			return READ_FAILED;
		}
		for (j = 0; j < nof_vertex; j++)
		{
			if (j == *temp)
			{
				*curr_row = 1;
				/*TODO: maybe here count nnz for array imp*/
				temp++;
			}
			else
			{
				*curr_row = 0;
			}
			curr_row++;
		}
		A->add_row(A, start_row, i);
		curr_vertex++;
	}
	free(curr_row);
	free(temp);
	return NONE;
}

/*TODO:check g is sorted (assume it in this function)*/
errors compute_modularity_matrix(spmat *A, int *g, int *degree, int M, spmat *B_g)
{
	int current_row; /*from group g*/
	int i, j, index = 0;
	double *row, *start_row;
	double *expected_nof_edges_row, *start_expected_nof_edges_row; /*(k_i*k_j)/M*/
	int *temp_i, *temp_j;

	start_row = (double *)malloc((B_g->n) * sizeof(double));
	if (!start_row)
	{
		return ALLOCATION_FAILED;
	}
	start_expected_nof_edges_row = (double *)malloc((B_g->n) * sizeof(double));
	if (!start_expected_nof_edges_row)
	{
		return ALLOCATION_FAILED;
	}

	expected_nof_edges_row = start_expected_nof_edges_row;
	temp_i = g;
	temp_j = g;
	for (i = 0; i < (A->n); i++)
	{
		if (*temp_i == 1)
		{ /* i is in g*/
			for (j = 0; j < (A->n); j++)
			{
				if (*temp_j == 1)
				{ /* j is in g*/
					*expected_nof_edges_row = -((degree[i] * degree[j]) / M);
					expected_nof_edges_row++;
				}
				temp_j++;
			}
			temp_j = g; /*reset temp_j*/
			A->sum_rows(A, i, start_expected_nof_edges_row, start_row);
			B_g->add_row(B_g, start_row, index);
			index++;
		}
		temp_i++;
	}

	free(start_row);
	free(start_expected_nof_edges_row);
	return NONE;
}

void eigen2s(double *eigen, int *s, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		*s = IS_POSITIVE(*eigen) ? 1 : -1;
		s++;
		eigen++;
	}
}

errors compute_f_g(spmat *B_g, double *f_g)
{
	int i;
	double *unit_vector, *temp;
	unit_vector = (double *)malloc((B_g->n) * sizeof(double));
	if (!unit_vector)
	{
		return ALLOCATION_FAILED;
	}
	temp = unit_vector;
	while (temp < unit_vector + B_g->n)
	{ /*setting every value of unit_vector to 1*/
		*temp = 1;
		temp++;
	}
	B_g->mult(B_g, unit_vector, f_g);

	free(unit_vector);
	return NONE;
}
/*TODO: check if delta_f is OK*/
errors compute_B_hat(spmat *B_g, double *f_g, spmat *B_hat)
{
	int i;
	double *start_row;
	double *temp_delta_f, *delta_f; /*represent delta*f */
	double *temp_f_g;
	start_row = (double *)malloc((B_g->n) * sizeof(double));
	if (!start_row)
	{
		return ALLOCATION_FAILED;
	}
	delta_f = (double *)malloc((B_g->n) * sizeof(double));
	if (!delta_f)
	{
		return ALLOCATION_FAILED;
	}
	temp_delta_f = delta_f;
	while (temp_delta_f < delta_f + B_g->n)
	{
		*temp_delta_f = 0;
		temp_delta_f++;
	}
	temp_f_g = f_g;
	temp_delta_f = delta_f;
	for (int i = 0; i < B_g->n; i++)
	{
		*temp_delta_f = *temp_f_g;
		B_g->sum_rows(B_g, i, delta_f, start_row);
		B_hat->add_row(B_hat, start_row, i);

		*temp_delta_f = 0;
		temp_delta_f++;
		temp_f_g++;
	}

	free(start_row);
	free(delta_f);
	return NONE;
}

double get_eigen_value(double *eigen, double *prev_eigen, int size)
{
	double mone = dot_product(eigen, prev_eigen, size);
	double mehane = dot_product(prev_eigen, prev_eigen, size);
	return mone / mehane;
}

errors compute_C_prime(spmat *C, spmat *C_prime, double *C_1norm)
{
	double *sum_of_col;
	double *row;
	double *temp, *norm_vector;
	int i;
	sum_of_col = (double *)malloc((C->n) * sizeof(double));
	if (!sum_of_col)
	{
		return ALLOCATION_FAILED;
	}
	row = (double *)malloc((C->n) * sizeof(double));
	if (!row)
	{
		return ALLOCATION_FAILED;
	}
	norm_vector = (double *)malloc((C->n) * sizeof(double));
	if (!norm_vector)
	{
		return ALLOCATION_FAILED;
	}
	for (temp = norm_vector; temp < norm_vector + C->n; temp++)
	{
		*temp = 0;
	}

	*C_1norm = C->sum_of_largest_column(C, sum_of_col);
	temp = norm_vector;
	for (i = 0; i < C->n; i++)
	{
		*temp = *C_1norm;
		C->sum_rows(C, i, norm_vector, row);
		C_prime->add_row(C_prime, row, i);
		*temp = 0;
		temp++;
	}

	free(sum_of_col);
	free(row);
	free(temp);
	return NONE;
}