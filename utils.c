#include <stdio.h>
#include <stdlib.h>
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

void print_vector(double *vector, int size)
{
	int i;
	printf("(");
	for (i = 0; i < size - 1; i++)
	{
		printf("%f ,", *vector);
		vector++;
	}
	printf("%f)\n\n\n", *vector);
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
	printf("%d)\n\n\n", *vector);
}

/*the result is stored in vector - at first it should be initialized with random numbers*/
/*calculates the dominant eigen vector*/
/*PROBLEM: currently the loop doesnt stop*/
Error power_iteration(spmat *mat, double *vector)
{
	int nof_iterations = 0;
	int stop = 0;
	double magnitude;
	double *mul_result;
	double *ptr_v, *ptr_r; /*pointers to vector and mul_result*/
	mul_result = (double *)malloc((mat->n) * sizeof(double));
	if (!mul_result)
	{
		return ALLOCATION_FAILED;
	}
	while (stop == 0 && nof_iterations <= 200)
	{
		mat->mult(mat, vector, mul_result); /*result=A*vector*/
		magnitude = sqrt(dot_product(mul_result, mul_result, mat->n));
		printf("--------------magnitute is:%f, mul result is:----------------\n", magnitude);
		print_vector(mul_result, mat->n);
		stop = 1;
		ptr_r = mul_result;
		for (ptr_v = vector; ptr_v < vector + mat->n; ptr_v++)
		{
			if (IS_POSITIVE(fabs(*ptr_v - ((*ptr_r) / magnitude))))
			{
				stop = 0;
			}
			*ptr_v = (*(ptr_r)) / magnitude; /*updating vector to the next vector*/
			ptr_r++;
		}
		printf("-----------------------vector is:-----------------\n");
		print_vector(vector, mat->n);
		nof_iterations++;
	}

	free(mul_result);
	return NONE;
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
			error = compute_modularity_matrix_row(mat, i, g, degrees, M, B_g_row);
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
	/*
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
	*/
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
Error compute_modularity_matrix_row(spmat *A, int row, group *g, int *degrees, double M, double *B_g_row)
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
		}
		degrees++;
		g_members++;
		temp++;
	}
	printf("Bg row: \n");
	print_vector(B_g_row, g->size);
	row_sum = A->add_to_row(A, row, B_g_row, g);
	B_g_row[row] -= row_sum; /*for B_hat*/
	return NONE;
}

double compute_modularity_value(spmat *B_g, double *s)
{
	/*s and g are the same - represent the group division*/
	double *Bs;
	int size = B_g->n;

	Bs = (double *)malloc(size * sizeof(double));
	if (!Bs)
	{
		printf("malloc failed on pointer Bs\n");
		return -1;
	}

	B_g->mult(B_g, s, Bs);

	return 0.5 * dot_product(s, Bs, size);
}

void eigen2s(double *eigen, group *g1, group *g2, double *s, int size)
{
	int i;
	char *g1_members = g1->members, *g2_members = g2->members;
	g1->size = 0;
	g2->size = 0;
	for (i = 0; i < size; i++)
	{
		if (IS_POSITIVE(*eigen))
		{
			*g1_members = 1;
			(g1->size)++;
			*g2_members = 0;
			*s = 1;
		}
		else
		{
			*g2_members = 1;
			(g2->size)++;
			*g1_members = 0;
			*s = -1;
		}
		eigen++;
		g1_members++;
		g2_members++;
		s++;
	}
}
