#include "spmat.h"
#include <stdlib.h>
#include <stdio.h>
#define print 0

typedef struct node
{
	double val;
	int index;
	struct node *next;
} node;

typedef struct arrays
{
	double *values;
	int *cols;
	int *rows;
	int curr;
} arrays;

typedef struct list
{
	node **rows;
} list;

void print_list(node *head)
{
	if (head == NULL)
	{
		printf("NULL\n");
	}
	else
	{
		printf("%f->", head->val);
		print_list(head->next);
	}
}
/*-------------------------------------------------------ARRAY_IMP----------------------------------------------------------*/

void addRow_arrays(spmat *A, const double *row, int i)
{
	const double *p = row;
	int j = 0;
	char is_first = 1;
	arrays *helper = (arrays *)(A->handle);
	for (; p < row + A->n; p++)
	{
		if (*p != 0)
		{
			(helper->values)[helper->curr] = *p;
			(helper->cols)[helper->curr] = j;
			if (is_first)
			{
				(helper->rows)[i] = helper->curr;
				is_first = 0;
			}
			(helper->curr)++;
		}
		j++;
	}
	if (is_first)
	{
		(helper->rows)[i] = helper->curr;
		is_first = 0;
	}
	if (i == (A->n) - 1)
	{
		(helper->rows)[i + 1] = helper->curr;
	}
	if (print)
	{
		printf("%d\n", (helper->rows)[i]);
	}
}

void free_arrays(spmat *A)
{
	arrays *helper = (arrays *)(A->handle);
	free(helper->values);
	free(helper->cols);
	free(helper->rows);
	free(helper);
	free(A);
}

void mult_arrays(const spmat *A, const double *v, double *result)
{
	int *row, *col;
	int count;
	double *vals;
	arrays *helper;
	double sum;

	helper = (arrays *)(A->handle);
	vals = helper->values;
	col = helper->cols;
	row = helper->rows;
	sum = 0;

	for (; row < (helper->rows) + (A->n); row++)
	{
		count = *(row + 1) - (*row);
		while (count > 0)
		{
			sum += (*vals) * (*(v + *col));
			vals++;
			col++;
			count--;
		}
		*result = sum;
		sum = 0;
		result++;
	}
}

spmat *spmat_allocate_array(int n, int nnz)
{
	spmat *spm;
	arrays *arrs;
	spm = (spmat *)malloc(sizeof(spmat));
	if (!spm)
	{
		return NULL;
	}
	arrs = (arrays *)malloc(sizeof(arrays));
	if (!arrs)
	{
		free(spm);
		return NULL;
	}
	(arrs->values) = (double *)malloc(nnz * sizeof(double));
	if (!(arrs->values))
	{
		free(spm);
		free(arrs);
		return NULL;
	}
	(arrs->cols) = (int *)malloc(nnz * sizeof(int));
	if (!(arrs->cols))
	{
		free(spm);
		free((arrs->values));
		free(arrs);
		return NULL;
	}
	(arrs->rows) = (int *)malloc((n + 1) * sizeof(int));
	if (!(arrs->rows))
	{
		free(spm);
		free((arrs->values));
		free((arrs->cols));
		free(arrs);
		return NULL;
	}
	(arrs->curr) = 0;
	(spm->handle) = arrs;
	(spm->n) = n;
	(spm->add_row) = addRow_arrays;
	(spm->free) = free_arrays;
	(spm->mult) = mult_arrays;
	return spm;
}

/*--------------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------LIST_IMP----------------------------------------------------------*/
node *create_list(const double *row, int n)
{
	node *head = NULL, *tail = NULL;
	const double *p = row;
	char first = 1;
	int counter = 0;
	if (print)
	{
		printf("%d create_list\n", n);
	}
	for (; p < row + n; p++)
	{
		if (*p != 0)
		{

			if (first)
			{
				head = (node *)malloc(sizeof(node));
				if (!head)
				{
					return NULL;
				}
				head->val = *p;
				head->index = counter;
				tail = head;
				first = 0;
				if (print)
				{
					printf("first succeeded");
				}
			}
			else
			{
				tail->next = (node *)malloc(sizeof(node));
				tail = tail->next;
				if (!tail)
				{
					return NULL;
				}
				tail->val = *p;
				tail->index = counter;
			}
		}
		counter++;
	}
	if (tail != NULL)
	{
		tail->next = NULL;
	}
	if (print)
	{
		printf("\n");
		print_list(head);
		printf("\n");
	}
	return head;
}

void addRow_list(spmat *A, const double *row, int i)
{
	if (print)
	{
		printf("addRow_list\n");
	}
	*((((list *)(A->handle))->rows) + i) = create_list(row, A->n);
}
/*TODO: maybe replace with iterative version*/
void delete_list(node *l)
{
	if (l != NULL)
	{
		if (print)
		{
			printf("delete node %f\n", l->val);
		}
		delete_list(l->next);
		free(l);
	}
}

void free_list(spmat *A)
{
	node **l = ((list *)(A->handle))->rows;
	if (print)
	{
		printf("free_list\n");
	}
	for (; l < (((list *)(A->handle))->rows) + A->n; l++)
	{
		if (print)
		{
			print_list(*l);
		}
		delete_list(*l);
	}
	if (print)
	{
		printf("pass delete\n");
	}
	free(((list *)(A->handle))->rows);
	free(A->handle);
	free(A);
}

void mult_list(const spmat *A, const double *v, double *result)
{
	double sum;
	node **currRow = ((list *)(A->handle))->rows; /* current row*/
	node *currElem = *currRow;					  /* current element*/
	double *currRes = result;					  /* current result element*/
	int i;

	for (i = 0; i < A->n; i++)
	{
		sum = 0;
		while (currElem != NULL)
		{
			sum += (currElem->val) * (*(v + (currElem->index)));
			currElem = currElem->next;
		}
		*(currRes) = sum;
		currRes++;
		currRow++;
		currElem = *currRow;
	}
}

char equal2_list(const spmat *A, const spmat *B)
{
	node *temp_A, *temp_B;
	int i;
	node **rows_A = ((list *)(A->handle))->rows;
	node **rows_B = ((list *)(B->handle))->rows;

	if (A->n != B->n)
	{
		return 0;
	}
	for (i = 0; i < A->n; i++)
	{
		temp_A = *rows_A;
		temp_B = *rows_B;
		while (temp_A != NULL && temp_B != NULL)
		{
			if (temp_A->index != temp_B->index || temp_A->val != temp_B->val)
			{
				return 0;
			}
			temp_A++;
			temp_B++;
		}
		if (temp_A != NULL || temp_B != NULL)
		{
			return 0;
		}
		rows_A++;
		rows_B++;
	}
	return 1;
}

void sum_rows_list(const spmat *A, int row, double *row2add, double *result)
{
	node **rows = ((list *)(A->handle))->rows;
	node *req_row = *(rows + row);
	int i;
	double *end = result + A->n;
	if (req_row == NULL)
	{ /*req_row is all zeroes*/
		while (result < end)
		{
			*result = *row2add;
			result++;
		}
	}
	else
	{ /*req_row is not all zeroes*/
		for (i = 0; i < A->n; i++)
		{
			if (i == req_row->index)
			{
				*result = req_row->val + *row2add;
				req_row++;
			}
			else
			{
				*result = *row2add;
			}
			result++;
		}
	}
	return;
}

void print_matrix_list(const spmat *mat)
{
	print_list(*(((list *)mat->handle)->rows));
	print_list(*(((list *)mat->handle)->rows + 1));
}

spmat *spmat_allocate_list(int n)
{
	list *l;
	spmat *spm;
	l = (list *)malloc(sizeof(list));
	if (!l)
	{
		return NULL;
	}
	l->rows = (node **)malloc(n * sizeof(node *));
	if (!(l->rows))
	{
		free(l);
		return NULL;
	}
	spm = (spmat *)malloc(sizeof(spmat));
	if (!spm)
	{
		free(l->rows);
		free(l);
		return NULL;
	}
	spm->n = n;
	spm->handle = l;

	spm->add_row = addRow_list;
	spm->free = free_list;
	spm->mult = mult_list;
	spm->equal2 = equal2_list;
	spm->sum_rows = sum_rows_list;
	spm->print_matrix = print_matrix_list;
	return spm;
}
/*--------------------------------------------------------------------------------------------------------------------*/