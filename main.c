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
	/*--------------------try to open the input file---------------------*/
	input_file = fopen(argv[1], "r");
	if (!input_file)
	{
		handle_errors(READ_FAILED,"main",argv[1]);
	}
	/*-----------------try to read the number of vertexes-----------------*/
	if (fread(&nof_vertex, sizeof(int), 1, input_file) != 1)
	{
		handle_errors(READ_FAILED,"main","input_file");
	}

	/*---------allocate vector to save the degree of each vector (k_i)-----*/
	alloc(degrees,int,nof_vertex,"main","degrees",5);
	/*---------------try allocate sparse matrix using list imp--------------*/
	A = spmat_allocate_list(nof_vertex);

	/*------------------compute the adj matrix (as spars matrix)-------------*/
	read_input(input_file, A, degrees, nof_vertex);

	/*---------------------------alocate group set P & O----------------------*/
	alloc(g,group,1,"main","g",5);
	alloc(g->members,int,nof_vertex,"main","g->members",5);

	/*-------------------------------initial calculations----------------------------*/
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
	algo_3(A, degrees, P, O, nof_vertex);

	/*----------------------------write the division in the output_file--------------*/
	output_file = fopen(argv[2], "w");
	if (!output_file)
	{
		handle_errors(WRITE_FAILED,"main",argv[2]);
	}
	write2_output_file(output_file, O);
	fclose(output_file);
	/*-----------printing the output file-----------------*/
	output_file = fopen(argv[2], "r");
	if (!output_file)
	{
		printf("[main]: output file is invalid2\n");
		return 5;
	}
	printf("[main]: actual: \n");
	print_output(output_file, nof_vertex);

	A->free(A);
	free(degrees);
	fclose(input_file);
	fclose(output_file);
	P->free_set(P);
	O->free_set(O);

	printf("FINISHED- %ld\n", (clock() - start) / CLOCKS_PER_SEC);

	return 0;
}
