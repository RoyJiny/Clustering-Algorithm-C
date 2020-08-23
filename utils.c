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

/*the result is stored in vector - at first it should be initialized with random numbers*/
/*calculates the dominant eigen vector*/
/*PROBLEM: currently the loop doesnt stop*/
void power_iteration(spmat *mat, double *vector, double epsilon)
{
	int stop = 0;
	double magnitude;
	double *mul_result;
	double *ptr_v, *ptr_r; /*pointers to vector and mul_result*/
	mul_result = (double *)malloc((mat->n) * sizeof(double));
	if (!mul_result)
	{
		return;
	}
	while (stop == 0)
	{
		mat->mult(mat, vector, mul_result); /*result=A*vector*/
		magnitude = sqrt(dot_product(mul_result, mul_result, mat->n));
		stop = 1;
		ptr_r = mul_result;
		for (ptr_v = vector; ptr_v < vector + mat->n; ptr_v++)
		{
			if (fabs(*ptr_v - ((*ptr_r) / magnitude)) > epsilon)
			{
				stop = 0;
			}
			*ptr_v = (*(ptr_r)) / magnitude; /*updating vector to the next vector*/
			ptr_r++;
		}
	}
	free(mul_result);
}

/*calculate the eigen vector with the matrix and the vector from the last iteration of power iteration*/
int calculate_eigen_value(spmat *mat, double *eigen_vector)
{
	int size = mat->n;
	double *Ab;
	double numerator;
	double denominator;
	Ab = malloc(size * sizeof(double));
	if (!Ab)
	{
		printf("malloc failed on pointer Ab\n");
		return 5;
	}

	mat->mult(mat, eigen_vector, Ab);
	numerator = dot_product(eigen_vector, Ab, size);
	denominator = dot_product(eigen_vector, eigen_vector, size);
	return numerator / denominator;

	free(Ab);
}

/*return 1 if there was an error and 0 otherwise*/
char handle_errors(Error error, char *name)
{
	switch (error)
	{
	case ALLOCATION_FAILED:
		printf("allocation failed in: %s\n", name);
		return 1;
	case READ_FAILED:
		printf("read failed in: %s\n", name);
		return 1;
	default:
		return 0;
	}
}

Error read_input(FILE *input, spmat *A, int *degree, int nof_vertex)
{
	double *curr_row;
	double *start_row; /*save the start of the row to add_row function in spmat*/
	int *temp;
	int *start_temp; /*save the start of temp for reuse*/
	int *curr_vertex;
	int i, j;

	curr_vertex = degree; /*fill degree*/
	start_row = (double *)malloc(nof_vertex * sizeof(double));
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
		if ((signed int)fread(temp, sizeof(int), *curr_vertex, input) != *curr_vertex)
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
	free(start_row);
	free(start_temp);
	return NONE;
}

Error compute_modularity_matrix(spmat *A, double *g, int *degree, double M, spmat *B_g)
{
	int i, j, index = 0;
	double *start_row;
	double *expected_nof_edges_row, *start_expected_nof_edges_row; /*(k_i*k_j)/M*/
	double *temp_i, *temp_j;

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

	printf("size of A:%d \n", A->n);
	printf("size of B_g:%d \n", B_g->n);

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
					*expected_nof_edges_row = -(((double)degree[i] * (double)degree[j]) / M);
					expected_nof_edges_row++;
				}
				temp_j++;
			}
			temp_j = g; /*reset temp_j*/
			A->sum_rows(A, i, start_expected_nof_edges_row, start_row);
			B_g->add_row(B_g, start_row, index);
			index++;
		}
		expected_nof_edges_row = start_expected_nof_edges_row;
		temp_i++;
	}

	free(start_row);
	free(start_expected_nof_edges_row);
	return NONE;
}

double compute_modularity_value(spmat *B_g, double *g)
{
	/*s and g are the same - represent the group division*/
	int size = B_g->n;
	double *Bs;

	Bs = (double *)malloc(size * sizeof(double));
	if (!Bs)
	{
		printf("malloc failed on pointer Bs\n");
		return -1;
	}

	B_g->mult(B_g, g, Bs);

	return 0.5 * dot_product(g, Bs, size);
}
