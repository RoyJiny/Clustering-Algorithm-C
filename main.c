#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "utils.h"

void print_vector(double *v, int n)
{
	double *ptr = v;
	for (; ptr < v + n; ptr++)
	{
		printf("%f	", *ptr);
	}
}

void test_input_read(spmat *mat, FILE *compare)
{
	int n, i, *row;
	spmat *A = spmat_allocate_list(n);
	if (!A)
	{
		printf("error allocating sparse matrix");
		return;
	}
	if (fread(&n, sizeof(int), 1, compare) != 1)
	{
		printf("reading error");
		return;
	}
	if (fread(&n, sizeof(int), 1, compare) != 1)
	{
		printf("reading error");
		return;
	}
	row = (int *)malloc(n * sizeof(int));
	if (!row)
	{
		printf("error allocating row");
		return;
	}

	for (i = 0; i < n; i++)
	{
		if (fread(row, sizeof(int), n, compare) != n)
		{ /* read row*/
			printf("reading error");
			return;
		}
		A->add_row(A, row, i);
	}

	if (A->equal2(mat, A))
	{
		printf("success");
	}
	else
	{
		printf("failed");
	}
}

int main(int argc, char *argv[])
{
	FILE *input;
	spmat *A;
	int nof_vertex;
	int *degrees;
	int *temp;
	int M = 0; /*sum of degrees*/
	Error error;

	/*for tests:*/
	FILE *compare;
	compare = fopen(argv[2], "r");
	if (!compare)
	{
		printf("compare file is invalid");
		return 5;
	}

	/*try to open the input file*/
	input = fopen(argv[1], "r");
	if (!input)
	{
		printf("input file is invalid");
		return 5;
	}

	/*try to read the number of vertexes*/
	if (fread(&nof_vertex, sizeof(int), 1, input) != 1)
	{
		printf("read from input file failed");
		return 5;
	}

	/*allocate vector to save the degree of each vector (k_i)*/
	degrees = (int *)malloc(nof_vertex * sizeof(int));
	if (!degrees)
	{
		printf("malloc failed on pointer degree");
		return 5;
	}

	/*try allocate sparse matrix using list imp*/
	A = spmat_allocate_list(nof_vertex);
	if (!A)
	{
		printf("sparse matrix allocation failed on A");
		return 5;
	}

	error = read_input(input, A, degrees, nof_vertex);

	if (handle_errors(error, "read_input"))
	{
		return 5;
	}
	/*computing M*/
	for (temp = degrees; temp < degrees + nof_vertex; temp++)
	{
		M += *temp;
	}

	/*------------------test start (read input)---------------------*/
	test_input_read(A, compare);
	/*------------------test end (read input)---------------------*/
	return 0;
}