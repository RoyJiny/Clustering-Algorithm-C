#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "algo.h"

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
	group_set *P, *O;
	group *g;
	spmat *A;
	int nof_vertex;
	int *degrees;
	Error error;

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
	/*---------------------------alocate group set P & O----------------------*/
	g = (group*) malloc(sizeof(group));
	if(!g){
		printf("allocation failed on g");
		return 5;
	}
	g->members = (char*) malloc(nof_vertex*sizeof(char));
	if(!(g->members)){
		printf("allocation failed on g->members");
		return 5;
	}
	memset(g->members, 1, nof_vertex); /*initial group of all the vertex*/
	P = allocate_group_set();
	O = allocate_group_set();

	error = algo_3(A, degrees, P, O, nof_vertex);
	if (handle_errors(error, "algo_3"))
	{
		return 5;
	}

	print_stack(O, nof_vertex);
	A->free(A);
	free(degrees);
	fclose(input_file);
	free(g->members);
	free(g);
	P->free(P);
	O->free(O);

	printf("FINISHED\n");

	return 0;
}