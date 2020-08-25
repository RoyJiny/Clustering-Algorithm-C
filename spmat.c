#include "spmat.h"
#include <stdlib.h>
#include <math.h>
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
	double *cols_sum;
} list;

void print_list(node *head)
{
	if (head == NULL)
	{
		printf("NULL\n");
	}
	else
	{
		printf(" (%d)%f->", head->index, head->val);
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
node *create_list(const double *row, int n ,double *sum_of_cols)
{
	node *head = NULL, *tail = NULL;
	const double *p = row;
	double* temp = sum_of_cols;
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
			*temp += fabs(tail->val);
		}
		counter++;
		temp++;
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
	node **rows = ((list *)(A->handle))->rows;
	double *sum_of_cols = ((list *)(A->handle))->cols_sum;
	if (print)
	{
		printf("addRow_list\n");
	}

	*(rows + i) = create_list(row, A->n , sum_of_cols);
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
	free(((list *)(A->handle))->cols_sum);
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
	printf("start mult_list\n");
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
	printf("end mult_list\n");
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
			temp_A = temp_A->next;
			temp_B = temp_B->next;
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
	double *temp = result;
	if (req_row == NULL)
	{ /*req_row is all zeroes*/
		while (temp < result + A->n)
		{
			*temp = *row2add;
			temp++;
			row2add++;
		}
	}
	else
	{ /*req_row is not all zeroes*/
		for (i = 0; i < A->n; i++)
		{
			if(req_row != NULL){
				if (i == req_row->index)
				{
					*result = req_row->val + *row2add;
					req_row = req_row->next;
				}
				else
				{
					*result = *row2add;
				}
				result++;
				row2add++;
			}
			else{
				*result = *row2add;
				result++;
				row2add++;
			}
		}
	}
	return;
}

void print_matrix_list(const spmat *mat)
{
	int i;
	node** rows = ((list *)(mat->handle))->rows;
	for (i = 0; i < mat->n; i++)
	{
		if (i < 10)
		{
			printf(" ");
		}
		printf("%d| ", i);
		print_list(*rows);
		printf("\n");
		rows++;
	}
}
/*return 1 if insert succeed and 0 otherwise*/
char insert2list(node* n1, double n ,int index){
	node* next = n1->next;
	node* new_node = (node*) malloc(sizeof(node));
	if(!new_node){
		return 0;
	}
	printf("in insert2list\n");
	new_node->val = n;
	new_node->index = index;

	n1->next = new_node;
	new_node->next = next;
	return 1;
}
/*return 1 if insert succeed and 0 otherwise*/
char add_by_index_list(const spmat *mat , int row, int index, double n)
{
	char res;
	node** rows = ((list *)(mat->handle))->rows;
	double* sum_of_cols = ((list *)(mat->handle))->cols_sum;
	node *req_row = *(rows+row), *temp ,*prev;
	temp = req_row;
	prev = temp;

	printf("\n");
	printf("index recieved is: %d\n",index);
	print_list(req_row);
	printf("\n");
	while(temp!=NULL){
		if(temp->index >= index){
			break;
		}
		prev = temp;
		printf("%d->",prev->index);
		temp = temp->next;
	}
	printf("\nprev index is: %d, prev val is: %f\n", prev->index, prev->val);
	if(temp == NULL || (temp->index) > index){
		res = insert2list(prev,n,index);
	}
	else{ /*temp->index == index*/
		temp->val = (temp->val) + n;

		if(temp->val == 0){
			printf("zero");
			prev->next = temp->next;
			sum_of_cols[index] -= fabs(temp->val);
			free(temp);
			return 1;
		}
		res = 1;
	}
	sum_of_cols[index] += fabs(n);
	return res;
}

double compute_1norm_list(const spmat *mat){
	double* sum_of_cols = ((list *)(mat->handle))->cols_sum;
	double max = *sum_of_cols;
	int i;
	sum_of_cols++;
	for(i=1; i<mat->n; i++){
		if(*sum_of_cols > max){
			max = *sum_of_cols;
		}
		sum_of_cols++;
	}
	return max;
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
	l->cols_sum = (double*)calloc(n, sizeof(double));/*all set to 0*/
	if(!(l->cols_sum))
	{
		free(l->rows);
		free(l);
		return NULL;
	}
	spm = (spmat *)malloc(sizeof(spmat));
	if (!spm)
	{
		free(l->rows);
		free(l->cols_sum);
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
	spm->add_by_index = add_by_index_list;
	spm->compute_1norm = compute_1norm_list;
	return spm;
}
/*--------------------------------------------------------------------------------------------------------------------*/