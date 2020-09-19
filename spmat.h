#ifndef _SPMAT_H
#define _SPMAT_H

#include "group.h"

typedef struct _spmat
{
	/* Matrix size (n*n) */
	int n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, const char *row, int i);

	/* Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/* compute the sum between A[row] and row2add and puts it in result.*/
	void (*sum_rows)(const struct _spmat *A, int row, double *row2add, double *result);

	/*print the matrix for testing*/
	void (*print_matrix)(const struct _spmat *mat);

	/*add n to mat[row][index] , return 1 if successful and 0 otherwise*/
	char (*add_by_index)(const struct _spmat *mat, int row, int index, double n);

	/*add A[row] to row*/
	double (*add_to_row)(const struct _spmat *A, int row_index, double *row, group *g);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void *handle;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat *spmat_allocate_list(int n);

#endif
