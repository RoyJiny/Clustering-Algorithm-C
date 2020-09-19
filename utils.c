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

void print_group(group *g)
{
	int i;
	int *g_mem = g->members;
	printf("(");
	for (i = 0; i < g->size - 1; i++)
	{
		printf("%d ,", *g_mem);
		g_mem++;
	}
	printf("%d)\n", *g_mem);
}

void print_stack(group_set *s)
{
	int counter = 0;
	group_node *curr = s->first;
	while (curr != NULL)
	{
		printf("group %d is\n:", counter);
		print_group(curr->value);
		printf("\n");
		curr = curr->next;
		counter++;
	}
}

void compute_1norm(spmat *A, group *g, int *degrees, double M, double *res)
{
	int i, j;
	double norm = 0, *col_sums, *runner, *B_row;

	alloc(col_sums,double,A->n,"compute_1norm","col_sums");
	alloc(B_row,double,A->n,"compute_1norm","B_row");

	for (runner = col_sums; runner < col_sums + (A->n); runner++)
	{
		*runner = 0;
	}

	for (i = 0; i < A->n; i++)
	{
		compute_modularity_matrix_row(A, i, g, degrees, M, B_row, i);

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
	*res = norm;

	free(col_sums);
	free(B_row);
}

void calculate_eigen_value(spmat *mat, double *eigen_vector, group *g, int *degrees, double M, double *B_g_row, double B_1norm, double *res)
{
	int i;
	double numerator;
	int *g_members = g->members;
	double *mult_vector, *runner, *B_g_runner;

	alloc(mult_vector,double,g->size,"calculate_eigen_value","mult_vector");

	runner = mult_vector;
	B_g_runner = B_g_row;
	for (i = 0; i < g->size; i++)
	{
		compute_modularity_matrix_row(mat, *g_members, g, degrees, M, B_g_row, i);
		*B_g_runner += B_1norm;
		*runner = dot_product(B_g_row, eigen_vector, g->size);
		runner++;
		B_g_runner++;
		g_members++;
	}

	numerator = dot_product(eigen_vector, mult_vector, g->size);
	*res = numerator;

	free(mult_vector);
}

void read_input(FILE *input, spmat *A, int *degree, int nof_vertex)
{
	char *curr_row;
	char *start_row; /*save the start of the row to add_row function in spmat*/
	int *temp;
	int *start_row_for_read; /*save the start of temp for reuse*/
	int *vertex_degree;
	int i, j, counter;

	vertex_degree = degree;
	alloc(start_row,char,nof_vertex,"read_input","start_row");
	alloc(start_row_for_read,int,nof_vertex,"read_input","start_row_for_read");

	for (i = 0; i < nof_vertex; i++)
	{
		curr_row = start_row;
		temp = start_row_for_read;
		if (fread(vertex_degree, sizeof(int), 1, input) != 1) /*try read k_i*/
		{
			handle_errors(READ_FAILED,"read_input", "input_file");
		}		
		if ((signed int)fread(temp, sizeof(int), *vertex_degree, input) != *vertex_degree) /*try read the k_i neighbors*/
		{
			handle_errors(READ_FAILED,"read_input", "input_file");
		}
		counter = 0;
		for (j = 0; j < nof_vertex; j++)
		{
			if (counter < *vertex_degree && j == *temp)
			{
				*curr_row = 1;
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
		vertex_degree++;
	}
	free(start_row);
	free(start_row_for_read);
}

void compute_modularity_matrix_row(spmat *A, int A_row, group *g, int *degrees, double M, double *B_g_row, int g_row)
{
	int i;
	double row_sum, *runner, *temp = B_g_row;
	int g_vertex, g_prev_vertex = 0;
	int *g_members = g->members;
	double row_degree = (double)degrees[A_row];

	if (!M)
	{
		handle_errors(DIVISION_BY_ZERO, "compute_modularity_matrix_row", "M");
	}
	runner = B_g_row;
	for (i = 0; i < g->size; i++)
	{
		g_vertex = *g_members;
		degrees = degrees + (g_vertex - g_prev_vertex);
		*temp = -(row_degree * (double)(*degrees)) / M;
		temp++;
		g_members++;
		g_prev_vertex = g_vertex;
		if(i<g_row) runner++;
	}

	row_sum = A->add_to_row(A, A_row, B_g_row, g);
	*runner -= row_sum;
}

void compute_score(spmat *A, int A_index, int g_index, group *g, double *s, double M, int *degrees, double *score, double *B_g_row)
{
	int i;
	double row_sum, *temp;
	int g_vertex, g_prev_vertex = 0;
	int *g_members = g->members;
	double row_degree = (double)degrees[A_index];
	int deg = *(degrees + A_index);
	temp = B_g_row;

	if (!M)
	{
		handle_errors(DIVISION_BY_ZERO, "compute_score", "M");
	}

	for (i = 0; i < g->size; i++)
	{
		g_vertex = *g_members;
		degrees = degrees + (g_vertex - g_prev_vertex);
		*temp = -(row_degree * (double)(*degrees)) / M;
		temp++;
		g_members++;
		g_prev_vertex = g_vertex;
	}

	row_sum = A->add_to_row(A, A_index, B_g_row, g);
	if (row_sum != row_sum) { printf("using 'unused return value'\n"); }

	*score = dot_product(B_g_row, s, g->size) * 4 * (*(s + g_index));
	*score += 4 * ((deg) * (deg)) / M;
	*score = 0.5 * (*score);
}

void compute_modularity_value(spmat *A, group *g, int *degrees, double *s, double M, double *B_g_row, double *mult_vector, double *res)
{
	double *runner;
	int *g_members;
	int i;
	g_members = g->members;
	runner = mult_vector;

	for (i = 0; i < g->size; i++)
	{
		compute_modularity_matrix_row(A, *g_members, g, degrees, M, B_g_row, i);
		*runner = dot_product(B_g_row, s, g->size);
		runner++;
		g_members++;
	}

	*res = 0.5 * dot_product(mult_vector, s, g->size);
}

int eigen2s(double *eigen, group *g, double *s)
{
	int i, g1_counter = 0;
	for (i = 0; i < g->size; i++)
	{
		if (IS_POSITIVE(*eigen))
		{
			*s = 1;
			g1_counter++;
		}
		else
		{
			*s = -1;
		}
		eigen++;
		s++;
	}
	return g1_counter;
}

void construct_g1g2(group *g, double *s, group *g1, group *g2, int g1_count)
{
	int i;
	int *g1_members, *g2_members, *g_members = g->members;

	alloc(g1->members,int,g1_count,"construct_g1g2","g1->members");
	alloc(g2->members,int,g->size-g1_count,"construct_g1g2","g2->members");
	g1_members = g1->members;
	g2_members = g2->members;

	g1->size = g1_count;
	g2->size = g->size - g1_count;
	for (i = 0; i < g->size; i++)
	{
		if (*s == 1)
		{
			*g1_members = *g_members;
			g1_members++;
		}
		else
		{
			*g2_members = *g_members;
			g2_members++;
		}
		s++;
		g_members++;
	}
}

void write_output_file(FILE *output, group_set *O)
{
	group *g;
	int nof_vertex_in_group, nof_groups = O->size;
	int *g_members;
	
	/*first number = number of groups*/
	if (fwrite(&nof_groups, sizeof(int), 1, output) != 1)
	{
		handle_errors(WRITE_FAILED, "write2_output_file", "output_file");
	}
	while (!(O->is_empty(O)))
	{
		g = O->pop(O);
		g_members = g->members;
		nof_vertex_in_group = g->size;
		/*number of vertex in the current group*/
		if (fwrite(&nof_vertex_in_group, sizeof(int), 1, output) != 1)
		{
			handle_errors(WRITE_FAILED, "write2_output_file", "output_file");
		}
		/*the actual vertex in the curr group*/
		if ((signed int)fwrite(g_members, sizeof(int), nof_vertex_in_group, output) != nof_vertex_in_group)
		{
			handle_errors(WRITE_FAILED, "write2_output_file", "output_file");
		}
		free(g->members);
		free(g);
	}
}

void print_errors(Error error, char *name, char *func)
{
	switch (error)
	{
	case ALLOCATION_FAILED:
		printf("[%s]: allocation failed for: %s\n", func, name);
		return;
	case READ_FAILED:
		printf("[%s]: read failed on %s\n", func, name);
		return;
	case DIVISION_BY_ZERO:
		printf("[%s]: division by zero, %s is zero\n", func, name);
		return;
	case WRITE_FAILED:
		printf("[%s]: write failed on: %s\n", func, name);
		return;
	case ENDLESS_LOOP:
		printf("[%s]: reached infinite loop threshold at: %s\n",func,name);
		return;
	default:
		return;
	}
}

void print_output(FILE *output, int nof_vertex)
{
	int i, nof_groups, nof_vertex_in_group, *curr;

	alloc(curr,int,nof_vertex,"print_output","curr");

	if (fread(&nof_groups, sizeof(int), 1, output) != 1)
	{
		exit(5);
	}
	printf("result:\ndivision into %d groups\n", nof_groups);
	for (i = 0; i < nof_groups; i++)
	{
		if (fread(&nof_vertex_in_group, sizeof(int), 1, output) != 1)
		{
			exit(5);
		}
		if ((signed int)fread(curr, sizeof(int), nof_vertex_in_group, output) != nof_vertex_in_group)
		{
			exit(5);
		}
		printf("group number %d is:\n", i);
		print_vector_int(curr, nof_vertex_in_group);
		printf("\n");
	}

	free(curr);
}
