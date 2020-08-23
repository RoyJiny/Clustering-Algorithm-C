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
		printf("%f	\n", *ptr);
	}
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
	int debug = 0;		/*for short debug prints*/
	int deep_debug = 0; /*for the long debug prints*/
	FILE *input_file;
	spmat *A;
	spmat *B_g;
	int nof_vertex;
	int *degrees;
	int *temp_i;
	double *temp_d;
	double *g;
	double *eigen_vector;
	double eigen_value;
	int M = 0; /*sum of degrees*/
	Error error;

	/*for tests:*/
	double Q;
	FILE *result_file;
	result_file = fopen(argv[2], "r");
	if (!result_file)
	{
		printf("result file is invalid\n");
		return 5;
	}

	printf("argc: %d\n", argc);

	/*try to open the input file*/
	input_file = fopen(argv[1], "r");
	if (!input_file)
	{
		printf("input file is invalid\n");
		return 5;
	}

	/*try to read the number of vertexes*/
	if (fread(&nof_vertex, sizeof(int), 1, input_file) != 1)
	{
		printf("read from input file failed\n");
		return 5;
	}

	/*allocate vector to save the degree of each vector (k_i)*/
	degrees = (int *)malloc(nof_vertex * sizeof(int));
	if (!degrees)
	{
		printf("malloc failed on pointer degree\n");
		return 5;
	}

	/*try allocate sparse matrix using list imp*/
	A = spmat_allocate_list(nof_vertex);
	if (!A)
	{
		printf("sparse matrix allocation failed on A\n");
		return 5;
	}
	B_g = spmat_allocate_list(nof_vertex); /*remember that B_g will have a different size for each iteration*/
	if (!B_g)
	{
		printf("sparse matrix allocation failed on B_g\n");
		return 5;
	}

	error = read_input(input_file, A, degrees, nof_vertex);

	if (handle_errors(error, "read_input\n"))
	{
		return 5;
	}
	/*computing M*/
	for (temp_i = degrees; temp_i < degrees + nof_vertex; temp_i++)
	{
		M += *temp_i;
	}

	/*------------------test start (read input)---------------------*/
	/*test_input_read(A, result_file);*/
	/*------------------test end (read input)---------------------*/

	g = (double *)malloc(nof_vertex * sizeof(double));
	if (!g)
	{
		printf("malloc failed on pointer g\n");
		return 5;
	}
	for (temp_d = g; temp_d < g + nof_vertex; temp_d++)
	{
		*temp_d = 1;
	}

	error = compute_modularity_matrix(A, g, degrees, (double)M, B_g);
	if (handle_errors(error, "compute_modularity_matrix\n"))
	{
		return 5;
	}

	if (debug == 1 && deep_debug == 1)
	{
		printf("A:\n");
		A->print_matrix(A);
		printf("B_g:\n");
		B_g->print_matrix(B_g);
	}

	Q = compute_modularity_value(B_g, g);
	if (debug == 1)
	{
		printf("Modularity value for B_g is: %f\n", Q);
	}

	/*calculate the eigen vector with power iteration*/
	eigen_vector = malloc(nof_vertex * sizeof(double));
	if (!eigen_vector)
	{
		printf("malloc failed on pointer eigen_vector\n");
		return 5;
	}
	for (temp_d = eigen_vector; temp_d < eigen_vector + nof_vertex; temp_d++)
	{
		*temp_d = rand();
	}
	power_iteration(B_g, eigen_vector, 0.00001); /*TODO: change epsilon to a variable*/

	if (debug == 1)
	{
		printf("dominant eigen vector: [");
		for (temp_d = eigen_vector; temp_d < eigen_vector + nof_vertex; temp_d++)
		{
			printf("%f,", *temp_d);
		}
		printf("]\n");
	}

	/*calculate the corresponding eigen value*/
	eigen_value = calculate_eigen_value(B_g, eigen_vector);

	if (debug == 1)
	{
		printf("and the corresponding eigen value is: %f\n", eigen_value);
	}

	A->free(A);
	B_g->free(B_g);
	free(degrees);
	free(g);
	fclose(input_file);
	fclose(result_file);

	printf("FINISHED\n");

	return 0;
}