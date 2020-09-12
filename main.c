#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "algo.h"

int run_num = 0;

void test_input_read(spmat *mat, FILE *compare)
{
	unsigned int n, i;
	char *row;
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
	row = (char *)malloc(n * sizeof(char));
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

Error create_graph(FILE *file, int numOfVertex, char empty, char full)
{
	int i, *vector, *runner, j;
	int random, temp;
	char **mat, **runner3, *runner2;
	/*--------------------allocations---------------------------*/
	vector = (int *)malloc(numOfVertex * sizeof(int));
	if (!vector)
	{
		return ALLOCATION_FAILED;
	}

	if (fwrite(&numOfVertex, sizeof(int), 1, file) != 1)
	{
		return WRITE_FAILED;
	}
	if (empty)
	{ /*creates a graph with isolated vertexes*/
		printf("generating empty graph of size %d\n", numOfVertex);
		for (runner = vector; runner < vector + numOfVertex; runner++)
		{
			*runner = 0;
		}
		if ((signed int)fwrite(vector, sizeof(int), numOfVertex, file) != numOfVertex)
		{
			return WRITE_FAILED;
		}
	}
	else if (full)
	{ /*creates a clique*/
		printf("generating full graph of size %d\n", numOfVertex);
		numOfVertex--;
		for (i = 0; i < numOfVertex + 1; i++)
		{
			runner = vector;
			for (j = 0; j < numOfVertex + 1; j++)
			{
				if (j != i)
				{
					*runner = j;
					runner++;
				}
			}
			if (fwrite(&numOfVertex, sizeof(int), 1, file) != 1)
			{
				return WRITE_FAILED;
			}
			if ((signed int)fwrite(vector, sizeof(int), numOfVertex, file) != numOfVertex)
			{
				return WRITE_FAILED;
			}
		}
	}
	else
	{ /*random connections*/
		printf("generating random graph of size %d\n", numOfVertex);
		mat = (char **)malloc(numOfVertex * sizeof(char *));
		if (!mat)
		{
			return ALLOCATION_FAILED;
		}
		runner3 = mat;
		for (i = 0; i < numOfVertex; i++)
		{
			*runner3 = (char *)malloc(numOfVertex * sizeof(char));
			if (!(*runner3))
			{
				return ALLOCATION_FAILED;
			}
			runner3++;
		}
		for (i = 0; i < numOfVertex; i++)
		{
			for (j = 0; j < i; j++)
			{
				random = rand() % 25;
				if (random == 0)
				{
					mat[i][j] = 1;
					mat[j][i] = 1;
				}
				else
				{
					mat[i][j] = 0;
					mat[j][i] = 0;
				}
			}
		}
		runner3 = mat;
		for (i = 0; i < numOfVertex; i++)
		{
			runner2 = *runner3;
			runner = vector;
			temp = 0;
			for (j = 0; j < numOfVertex; j++)
			{
				if (*runner2)
				{
					temp++;
					*runner = j;
					runner++;
				}
				runner2++;
			}
			if (fwrite(&temp, sizeof(int), 1, file) != 1)
			{
				return WRITE_FAILED;
			}
			if ((signed int)fwrite(vector, sizeof(int), temp, file) != temp)
			{
				return WRITE_FAILED;
			}
			runner3++;
		}
	}
	for (runner3 = mat; runner3 < mat + numOfVertex; runner3++)
	{
		free(*runner3);
	}
	free(mat);
	free(vector);
	printf("done generating random graph of size %d\n", numOfVertex);
	return NONE;
}

Error test_create_graph(char *name, int numOfVertex, char empty, char full)
{
	FILE *file;
	/*Error error;
	spmat *A;
	int *degrees, i;
	A = spmat_allocate_list(numOfVertex);
	if (!A)
	{
		printf("sparse matrix allocation failed on A\n");
		return ALLOCATION_FAILED;
	}
	degrees = (int *)malloc(numOfVertex * sizeof(int));
	if (!degrees)
	{
		printf("malloc failed on pointer degree\n");
		return ALLOCATION_FAILED;
	}*/

	file = fopen(name, "w");
	if (!file)
	{
		printf("file is invalid\n");
		return WRITE_FAILED;
	}
	create_graph(file, numOfVertex, empty, full);
	fclose(file);
	/*file = fopen(name, "r");
	if (!file)
	{
		printf("file is invalid2\n");
		return WRITE_FAILED;
	}
	if (fread(&i, sizeof(int), 1, file) != 1)
	{
		return READ_FAILED;
	}
	error = read_input(file, A, degrees, numOfVertex);
	if (handle_errors(error, "read_input"))
	{
		return error;
	}
	A->print_matrix(A);
	fclose(file);
	A->free(A);
	free(degrees);*/
	return NONE;
}

int main(int argc, char *argv[])
{
	/*TODO: move it to algo_2*/
	/*int debug = 1;*/		/*for short debug prints*/
	/*int deep_debug = 1;*/ /*for the long debug prints*/
	FILE *input_file, *output_file;
	group_set *P, *O;
	group *g;
	spmat *A;
	int nof_vertex;
	int *degrees;
	Error error;
	clock_t start;

	if (!argc)
	{
		return 5;
	}
	srand(time(0));
	if (argv[1] == "--test")
	{
		test_create_graph(argv[2], atoi(argv[4]), 0, 0);
	}

	start = clock();
	/*--------------------try to open the input file---------------------*/
	input_file = fopen(argv[2], "r");
	if (!input_file)
	{
		printf("[main]: input file is invalid\n");
		return 5;
	}
	/*-----------------try to read the number of vertexes-----------------*/
	if (fread(&nof_vertex, sizeof(int), 1, input_file) != 1)
	{
		printf("[main]: read from input file failed\n");
		return 5;
	}

	/*---------allocate vector to save the degree of each vector (k_i)-----*/
	degrees = (int *)malloc(nof_vertex * sizeof(int));
	if (!degrees)
	{
		printf("[main]: malloc failed on pointer degree\n");
		return 5;
	}

	/*---------------try allocate sparse matrix using list imp--------------*/
	A = spmat_allocate_list(nof_vertex);
	if (!A)
	{
		printf("[main]: sparse matrix allocation failed on A\n");
		return 5;
	}

	/*------------------compute the adj matrix (as spars matrix)-------------*/
	error = read_input(input_file, A, degrees, nof_vertex);
	if (error != NONE)
	{
		return 5;
	}
	/*---------------------------alocate group set P & O----------------------*/
	g = (group *)malloc(sizeof(group));
	if (!g)
	{
		printf("[main]: allocation failed on g");
		return 5;
	}
	g->members = (int *)malloc(nof_vertex * sizeof(int));
	if (!(g->members))
	{
		printf("[main]: allocation failed on g->members");
		return 5;
	}

	/*-------------------------------initial calculations----------------------------*/
	g->size = nof_vertex;
	memset(g->members, 1, nof_vertex); /*trivial division*/
	P = allocate_group_set();
	O = allocate_group_set();
	P->push(P, g);
	error = algo_3(A, degrees, P, O, nof_vertex);
	if (error != NONE)
	{
		return 5;
	}
	/*if (handle_errors(error, "algo_3"))
	{
		return 5;
	}*/
	/*printf("stack_is :\n");
	print_stack(O, nof_vertex);*/
	/*----------------------------write the division in the output_file--------------*/
	output_file = fopen(argv[3], "w");
	if (!output_file)
	{
		printf("[main]: output file is invalid\n");
		return 5;
	}
	error = write2_output_file(output_file, O, nof_vertex);
	if (error != NONE)
	{
		return 5;
	}
	/*if (handle_errors(error, "write2_output_file"))
	{
		return 5;
	}*/
	fclose(output_file);
	/*-----------printing the output file-----------------*/
	output_file = fopen(argv[3], "r");
	if (!output_file)
	{
		printf("[main]: output file is invalid2\n");
		return 5;
	}
	printf("[main]: actual: \n");
	print_output(output_file, nof_vertex);

	/*result = fopen(argv[3], "r");
	if (!result)
	{
		printf("output file is invalid2\n");
		return 5;
	}
	printf("expected: \n");
	print_output(result, nof_vertex);
	fclose(result);*/

	A->free(A);
	free(degrees);
	fclose(input_file);
	fclose(output_file);
	/*free(g->members);*/ /*freed in algo_3 or in O->free_set*/
	/*free(g);*/
	/*maybe its not neccessary because the stack getting empty in write2_output file*/
	P->free_set(P);
	O->free_set(O);

	printf("FINISHED- %ld\n", (clock() - start) / CLOCKS_PER_SEC);

	return 0;
}