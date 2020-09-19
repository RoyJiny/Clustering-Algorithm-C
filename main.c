#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "algo.h"


int main(int argc, char *argv[])
{
	FILE *input_file, *output_file;
	group_set *P, *O;
	group *g;
	group *test;
	double *test_vec, *test_vec2, *test_vec3;
	spmat *A;
	int nof_vertex, i = 0;
	int *degrees, *g_members;
	clock_t start;

	if (!argc)
	{
		exit(4);
	}
	srand(time(0));
	start = clock();

	/*read input*/
	input_file = fopen(argv[1], "r");
	if (!input_file)
	{
		handle_errors(READ_FAILED,"main",argv[1]);
	}	
	if (fread(&nof_vertex, sizeof(int), 1, input_file) != 1)
	{
		handle_errors(READ_FAILED,"main","input_file");
	}

	alloc(degrees,int,nof_vertex,"main","degrees");
	A = spmat_allocate_list(nof_vertex);
	read_input(input_file, A, degrees, nof_vertex);

	printf("A is:\n");
	A->print_matrix(A);
	printf("degrees:\n");
	print_vector_int(degrees, A->n);
	alloc(test,group,1,"main","test");
	alloc(test->members,int,4,"main","test->members");
	test->size = 4;
	alloc(test_vec, double, 4, "main", "test_vec");
	alloc(test_vec2, double, 4, "main", "test_vec2");
	alloc(test_vec3, double, 4, "main", "test_vec3");
	g_members = test->members;
	*g_members = 0;
	*(g_members+1) = 2;
	*(g_members+2) = 8;
	*(g_members+3) = 19;
	*test_vec = 2;
	*(test_vec+1) = 2;
	*(test_vec+2) = 2;
	*(test_vec+3) = 2;
	A->mult(A, test_vec, test_vec2, test_vec3, test);
	printf("mult result is:\n");
	print_vector(test_vec2, 4);
	printf("nof elements is:\n");
	print_vector(test_vec3, 4);


	alloc(g,group,1,"main","g");
	alloc(g->members,int,nof_vertex,"main","g->members");

	/*trivial division*/
	g->size = nof_vertex;
	g_members = g->members;
	for (; i < g->size; i++)
	{
		*g_members = i;
		g_members++;
	}
	P = allocate_group_set();
	O = allocate_group_set();
	P->push(P, g);

	/*run*/
	algo_3(A, degrees, P, O, nof_vertex);

	/*write the division in the output file*/
	output_file = fopen(argv[2], "w");
	if (!output_file)
	{
		handle_errors(WRITE_FAILED,"main",argv[2]);
	}
	write_output_file(output_file, O);
	fclose(output_file);

	/*-----------printing the output file-----------------TODO-reomve at submission*/
	output_file = fopen(argv[2], "r");
	if (!output_file)
	{
		printf("[main]: output file is invalid2\n");
		return 5;
	}
	printf("[main]: actual: \n");
	print_output(output_file, nof_vertex);
	/*---------------------------------------------------------------------------*/

	A->free(A);
	free(degrees);
	fclose(input_file);
	fclose(output_file);
	P->free_set(P);
	O->free_set(O);

	printf("FINISHED- %ld\n", (clock() - start) / CLOCKS_PER_SEC); /*TODO: remove*/

	return 0;
}
