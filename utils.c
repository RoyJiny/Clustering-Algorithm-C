#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "utils.h"

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

Error compute_modularity_matrix(spmat *A, int *g, int *degree, int M, spmat *B_g)
{
	int i, j, index = 0;
	double *start_row;
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

	printf("size of A:%d ", A->n);
	printf("size of B_g:%d ", B_g->n);

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
			printf("start of func\n");
			A->sum_rows(A, i, start_expected_nof_edges_row, start_row);
			printf("end of func\n");
			B_g->add_row(B_g, start_row, index);
			printf("end of func2\n");
			index++;
			expected_nof_edges_row = start_expected_nof_edges_row;
		}
		temp_i++;
	}

	free(start_row);
	free(start_expected_nof_edges_row);
	return NONE;
}