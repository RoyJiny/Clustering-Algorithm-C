#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "algo.h"

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
	printf("%d)\n\n\n", *g_mem);
}

void test_input_read(spmat *mat, FILE *compare)
{
	unsigned int n, i;
	double *row;
	spmat *A = spmat_allocate_list(n);
	if (!A)
	{
		printf("error allocating sparse matrix\n");
		return;
	}
	if (fread(&n, sizeof(int), 1, compare) != 1)
	{
		printf("reading error - 1\n");
		return;
	}
	if (fread(&n, sizeof(int), 1, compare) != 1)
	{
		printf("reading error - 2\n");
		return;
	}
	row = (double *)malloc(n * sizeof(double));
	if (!row)
	{
		printf("error allocating row\n");
		return;
	}

	for (i = 0; i < n; i++)
	{
		if (fread(row, sizeof(int), n, compare) != n)
		{ /* read row*/
			printf("reading error - 3\n");
			return;
		}
		A->add_row(A, row, i);
	}

	if (A->equal2(mat, A))
	{
		printf("success\n");
	}
	else
	{
		printf("failed\n");
	}
}

int main(int argc, char *argv[])
{
	/*TODO: move it to algo_2*/
	/*int debug = 1;*/		/*for short debug prints*/
	/*int deep_debug = 1;*/ /*for the long debug prints*/
	FILE *input_file;
	spmat *A;
	int nof_vertex;
	int *degrees;
	double *eigen_vector, *temp_e;
	group *g, *g1, *g2;
	char *temp_g;
	Error error;

	srand(time(NULL));
	printf("argc: %d\n", argc);

	/*--------------------try to open the input file---------------------*/
	input_file = fopen(argv[1], "r");
	if (!input_file)
	{
		printf("input file is invalid\n");
		return 5;
	}

	/*-----------------try to read the number of vertexes-----------------*/
	if (fread(&nof_vertex, sizeof(int), 1, input_file) != 1)
	{
		printf("read from input file failed\n");
		return 5;
	}

	/*---------allocate vector to save the degree of each vector (k_i)-----*/
	degrees = (int *)malloc(nof_vertex * sizeof(int));
	if (!degrees)
	{
		printf("malloc failed on pointer degree\n");
		return 5;
	}

	/*---------------try allocate sparse matrix using list imp--------------*/
	A = spmat_allocate_list(nof_vertex);
	if (!A)
	{
		printf("sparse matrix allocation failed on A\n");
		return 5;
	}

	/*------------------compute the adj matrix (as spars matrix)-------------*/
	error = read_input(input_file, A, degrees, nof_vertex);
	if (handle_errors(error, "read_input\n"))
	{
		return 5;
	}
	/*this g is just for now*/
	g = (group *)malloc(sizeof(group)); /*remember that g will have a different size for each iteration*/
	if (!g)
	{
		printf("malloc failed on struct g\n");
		return 5;
	}
	g->members = (char *)malloc(nof_vertex * sizeof(char));
	if (!(g->members))
	{
		printf("malloc failed on pointer g->members\n");
		return 5;
	}
	g->size = nof_vertex;
	for (temp_g = g->members; temp_g < g->members + nof_vertex; temp_g++)
	{
		*temp_g = 1;
	}

	/*----calculate initial random vector to the power iteration algorithem----*/
	eigen_vector = malloc(nof_vertex * sizeof(double));
	if (!eigen_vector)
	{
		printf("malloc failed on pointer eigen_vector\n");
		return 5;
	}
	for (temp_e = eigen_vector; temp_e < eigen_vector + nof_vertex; temp_e++)
	{
		*temp_e = rand();
	}

	/*-----------allocate 2 groups to later hold the division into 2------------*/
	g1 = (group *)malloc(sizeof(group)); /*remember that g will have a different size for each iteration*/
	if (!g1)
	{
		printf("malloc failed on struct g\n");
		return 5;
	}
	g1->members = (char *)malloc(nof_vertex * sizeof(char));
	if (!(g1->members))
	{
		printf("malloc failed on pointer g->members\n");
		return 5;
	}

	g2 = (group *)malloc(sizeof(group)); /*remember that g will have a different size for each iteration*/
	if (!g2)
	{
		printf("malloc failed on struct g\n");
		return 5;
	}
	g2->members = (char *)malloc(nof_vertex * sizeof(char));
	if (!(g2->members))
	{
		printf("malloc failed on pointer g->members\n");
		return 5;
	}

	error = algo_2(A, degrees, eigen_vector, g, g1, g2);
	if (handle_errors(error, "algo_2"))
	{
		return 5;
	}

	printf("g1:\n");
	print_group(g1, nof_vertex);
	printf("g2:\n");
	print_group(g2, nof_vertex);

	A->free(A);
	free(degrees);
	free(g->members);
	free(g);
	free(g1->members);
	free(g1);
	free(g2->members);
	free(g2);
	free(eigen_vector);
	fclose(input_file);

	printf("FINISHED\n");

	return 0;
}