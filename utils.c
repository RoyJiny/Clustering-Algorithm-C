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
		printf("allocation failed in: %s", name);
		return 1;
	case READ_FAILED:
		printf("read failed in: %s", name);
		return 1;
	case NONE:
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