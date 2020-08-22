#ifndef _SPMAT_H
#define _SPMAT_H

#define IS_POSITIVE(x) ((x) > 0.00001)

typedef struct _spmat
{
	/* Matrix size (n*n) */
	int n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, const double *row, int i);

	/* Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void (*mult)(const struct _spmat *A, const double *v, double *result);

	/* compute the sum between A[row] and row2add and puts it in result.*/
	void (*sum_rows)(const spmat *A, int row, double *row2add, double *result);

	/* compute the sum of the largest colomn.
	* sum_col is a vector of size A->n which assumed to be empty*/
	double (*sum_of_largest_column)(const spmat *A, double *sum_col);

	/* return 1 if A==B and 0 otherwise.*/
	char (*equal2)(const spmat *A, const spmat *B);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void *handle;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat *spmat_allocate_list(int n);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat *spmat_allocate_array(int n, int nnz);

typedef struct
{
	double val;
	int index;
	struct node *next;
} node;

typedef struct
{
	node **rows;
} list;

typedef struct
{
	double *values;
	int *cols;
	int *rows;
	int curr;
} arrays;

#endif
