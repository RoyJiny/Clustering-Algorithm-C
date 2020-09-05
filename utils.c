#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
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

void print_vector(double *vector, int size)
{
	int i;
	printf("(");
	for (i = 0; i < size - 1; i++)
	{
		printf("%f ,", *vector);
		vector++;
	}
	printf("%f)\n", *vector);
}

void print_vector_int(int *vector, int size)
{
	int i;
	printf("(");
	for (i = 0; i < size - 1; i++)
	{
		printf("%d ,", *vector);
		vector++;
	}
	printf("%d)\n", *vector);
}

void print_group(group *g, int size)
{
	int i;
	char *g_mem = g->members;
	printf("(");
	for (i = 0; i < size - 1; i++)
	{
		printf("%d ,", *g_mem);
		g_mem++;
	}
	printf("%d)\n", *g_mem);
}

void print_stack(group_set *s, int size)
{
	int counter = 0;
	group_node *curr = s->first;
	while (curr != NULL)
	{
		printf("group %d is\n:", counter);
		print_group(curr->value, size);
		printf("\n");
		curr = curr->next;
		counter++;
	}
}

double compute_1norm(spmat *A, int *degrees, double M)
{
	int i, j;
	Error error;
	group *g;
	double norm = 0, *col_sums, *runner, *B_row;

	col_sums = (double *)malloc((A->n) * sizeof(double));
	if (!col_sums)
	{
		printf("allocation failed");
	}

	B_row = (double *)malloc((A->n) * sizeof(double));
	if (!B_row)
	{
		printf("allocation failed");
	}

	g = (group *)malloc(sizeof(group));
	if (!g)
	{
		printf("allocation failed");
	}
	g->members = (char *)malloc((A->n) * sizeof(char));
	if (!g->members)
	{
		printf("allocation failed");
	}
	g->size = A->n;
	memset(g->members, 1, A->n);

	memset(col_sums, 0, A->n);
	for (i = 0; i < A->n; i++)
	{
		error = compute_modularity_matrix_row(A, i, g, degrees, M, B_row, i);
		if (error != NONE)
		{
			printf("failed in compute_modularity_matrix_row - for B\n");
			return error;
		}

		j = 0;
		for (runner = col_sums; runner < col_sums + A->n; runner++)
		{
			*runner += fabs(*(B_row + j));
			j++;
		}
	}
	for (i = 0; i < A->n; i++)
	{
		if (*(col_sums + i) > norm)
		{
			norm = *(col_sums + i);
		}
	}

	free(col_sums);
	free(g->members);
	free(g);
	free(B_row);

	return norm;
}

/*calculate the eigen vector with the matrix and the vector from the last iteration of power iteration*/
double calculate_eigen_value(spmat *mat, double *eigen_vector, group *g, int *degrees, double M, double *B_g_row, double B_1norm)
{
	Error error;
	int i, g_count = 0;
	double numerator, denominator;
	char *g_members = g->members;
	double *mult_vector, *runner;

	mult_vector = (double *)malloc((g->size) * sizeof(double));
	if (!mult_vector)
	{
		printf("allocation failed in calculate eigen value");
	}
	runner = mult_vector;
	for (i = 0; i < mat->n; i++)
	{
		if (*g_members)
		{ /*curr vertex in g*/
			error = compute_modularity_matrix_row(mat, i, g, degrees, M, B_g_row, g_count);
			if (error != NONE)
			{
				printf("failed in compute_modularity_matrix_row\n");
				return error;
			}
			B_g_row[g_count] += B_1norm;
			*runner = dot_product(B_g_row, eigen_vector, g->size);
			runner++;
			g_count++;
		}
		g_members++;
	}
	numerator = dot_product(eigen_vector, mult_vector, g->size);
	denominator = dot_product(eigen_vector, eigen_vector, g->size);
	free(mult_vector);
	return numerator / denominator;
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
	case DIVISION_BY_ZERO:
		printf("division by zero in: %s\n", name);
		return 1;
	case WRITE_FAILED:
		printf("write failed in: %s\n", name);
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
	int i, j, counter;

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
		counter = 0;
		for (j = 0; j < nof_vertex; j++)
		{
			if (counter < *curr_vertex && j == *temp)
			{
				*curr_row = 1;
				/*TODO: maybe here count nnz for array imp*/
				temp++;
				counter++;
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

/*calculate a specific row of the modularity matrix, for group g (B_hat)*/
/*assuming g->members is a nof_vertex size vector with 1 indicates in g and 0 indicates not in g*/
Error compute_modularity_matrix_row(spmat *A, int row, group *g, int *degrees, double M, double *B_g_row, int g_count)
{
	int i;
	double row_sum, *temp = B_g_row;
	char *g_members = g->members;
	double row_degree = (double)degrees[row];

	if (!M)
	{
		return DIVISION_BY_ZERO;
	}

	for (i = 0; i < A->n; i++)
	{
		if (*g_members)
		{ /*the current vertex in g*/
			*temp = -(row_degree * (double)(*degrees)) / M;
			temp++;
		}
		degrees++;
		g_members++;
	}
	row_sum = A->add_to_row(A, row, B_g_row, g);
	B_g_row[g_count] -= row_sum; /*for B_hat*/
	return NONE;
}

double compute_modularity_value(spmat *A, group *g, int *degrees, double *s, double M, double *B_g_row, double *mult_vector)
{
	double *runner;
	char *g_members;
	Error error;
	int i, g_count;
	g_members = g->members;
	runner = mult_vector;
	g_count = 0;
	for (i = 0; i < A->n; i++)
	{
		/*do only if the vertex is in g*/
		if (*g_members)
		{
			error = compute_modularity_matrix_row(A, i, g, degrees, M, B_g_row, g_count);
			if (error != NONE)
			{
				printf("failed in compute_modularity_matrix_row\n");
				return error;
			}
			*runner = dot_product(B_g_row, s, g->size);
			runner++;
			g_count++;
		}
		g_members++;
	}
	return 0.5 * dot_product(mult_vector, s, g->size);
}

void eigen2s(double *eigen, group *g, double *s, int size)
{
	int i;
	char *g_members = g->members;
	for (i = 0; i < size; i++)
	{
		if (*g_members)
		{
			if (IS_POSITIVE(*eigen))
			{
				*s = 1;
			}
			else
			{
				*s = -1;
			}
			eigen++;
			s++;
		}
		g_members++;
	}
}

void construct_g1g2(group *g, double *s, group *g1, group *g2, int size)
{
	int i;
	char *g1_members = g1->members, *g2_members = g2->members, *g_members = g->members;
	g1->size = 0;
	g2->size = 0;
	for (i = 0; i < size; i++)
	{
		if(*g_members){
			if (*s == 1)
			{
				*g1_members = 1;
				*g2_members = 0;
				g1->size++;
			}
			else
			{
				*g2_members = 1;
				*g1_members = 0;
				g2->size++;
			}
			s++;
		}
		else{
			*g1_members = 0;
			*g2_members = 0;
		}
		g1_members++;
		g2_members++;
		g_members++;
	}
}

Error write2_output_file(FILE *output, group_set *O, int nof_vertex)
{
	group *g;
	int i, nof_vertex_in_group, nof_groups = O->size;
	int *curr, *runner;
	char *g_members;
	/*----------------allocations----------------*/
	curr = (int *)malloc(nof_vertex * sizeof(int));
	if (!curr)
	{
		return ALLOCATION_FAILED;
	}
	/*-------------------------------------------*/
	/*first number = number of groups*/
	if (fwrite(&nof_groups, sizeof(int), 1, output) != 1)
	{
		return WRITE_FAILED;
	}
	while (!(O->is_empty(O)))
	{
		g = O->pop(O);
		g_members = g->members;
		nof_vertex_in_group = g->size;
		/*number of vertex in the current group*/
		if (fwrite(&nof_vertex_in_group, sizeof(int), 1, output) != 1)
		{
			return WRITE_FAILED;
		}
		runner = curr;
		for (i = 0; i < nof_vertex; i++)
		{
			if (*g_members)
			{
				*runner = i;
				runner++;
			}
			g_members++;
		}
		/*the actual vertex in the curr group*/
		if ((signed int)fwrite(curr, sizeof(int), nof_vertex_in_group, output) != nof_vertex_in_group)
		{
			return WRITE_FAILED;
		}
		free(g->members);
		free(g);
	}
	free(curr);
	return NONE;
}

Error print_output(FILE *output, int nof_vertex)
{
	int i, nof_groups, nof_vertex_in_group, *curr;

	/*----------------allocations----------------*/
	curr = (int *)malloc(nof_vertex * sizeof(int));
	if (!curr)
	{
		return ALLOCATION_FAILED;
	}

	if (fread(&nof_groups, sizeof(int), 1, output) != 1)
	{
		return READ_FAILED;
	}
	printf("result:\ndivision into %d groups\n", nof_groups);
	for (i = 0; i < nof_groups; i++)
	{
		if (fread(&nof_vertex_in_group, sizeof(int), 1, output) != 1)
		{
			return READ_FAILED;
		}
		if ((signed int)fread(curr, sizeof(int), nof_vertex_in_group, output) != nof_vertex_in_group)
		{
			return READ_FAILED;
		}
		printf("group number %d is:\n", i);
		print_vector_int(curr, nof_vertex_in_group);
		printf("\n");
	}

	free(curr);
	return NONE;
}