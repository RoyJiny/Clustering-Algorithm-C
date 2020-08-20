#include "spmat.h"
#include <stdlib.h>
#include <stdio.h>
#define print 0


void print_list(node* head) {
	if (head == NULL) {
		printf("NULL\n");
	}
	else {
		printf("%f->", head->val);
		print_list(head->next);
	}
}
/*-------------------------------------------------------ARRAY_IMP----------------------------------------------------------*/

void addRow_arrays(spmat* A, const double* row, int i) {
	const double* p = row;
	int j = 0;
	char is_first = 1;
	arrays* helper = (arrays*)(A->handle);
	for (; p < row + A->n; p++) {
		if (*p != 0) {
			(helper->values)[helper->curr] = *p;
			(helper->cols)[helper->curr] = j;
			if (is_first) {
				(helper->rows)[i] = helper->curr;
				is_first = 0;
			}
			(helper->curr)++;
		}
		j++;
	}
	if (is_first) {
		(helper->rows)[i] = helper->curr;
		is_first = 0;
	}
	if (i == (A->n) - 1) {
		(helper->rows)[i+1] = helper->curr;
	}
	if (print) {
		printf("%d\n", (helper->rows)[i]);
	}
}

void free_arrays(spmat* A) {
	arrays* helper = (arrays*)(A->handle);
	free(helper->values);
	free(helper->cols);
	free(helper->rows);
	free(helper);
	free(A);
}

void mult_arrays(const spmat* A, const double* v, double* result) {
	int* row, * col;
	int count;
	double* vals;
	arrays* helper;
	double sum;

	helper = (arrays*)(A->handle);
	vals = helper->values;
	col = helper->cols;
	row = helper->rows;
	sum = 0;

	for (; row < (helper->rows) + (A->n); row++) {
		count = *(row + 1) - (*row);
		while (count > 0) {
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
/* return  1 if all the elements in row are 0, and 0 otherwise*/
char is_row_empty_arrays(const spmat* A , int row){
	arrays* helper = (arrays*)(A->handle);
	int* rows = ((helper->rows)+row);
	return *(rows+1) - *rows == 0 ? 1:0;
}

/*TODO: implement sum_rows_arrays*/

/*TODO: implement sum_of_largest_colomn_arrays*/

spmat* spmat_allocate_array(int n, int nnz) {
	spmat* spm;
	arrays* arrs;
	spm = (spmat*)malloc(sizeof(spmat));
	if (!spm) {
		return 0;
	}
	arrs = (arrays*)malloc(sizeof(arrays));
	if (!arrs) {
		free(spm);
		return 0;
	}
	(arrs->values) = (double*)malloc(nnz * sizeof(double));
	if (!(arrs->values)) {
		free(spm);
		free(arrs);
		return 0;
	}
	(arrs->cols) = (int*)malloc(nnz * sizeof(int));
	if (!(arrs->cols)) {
		free(spm);
		free((arrs->values));
		free(arrs);
		return 0;
	}
	(arrs->rows) = (int*)malloc((n + 1) * sizeof(int));
	if (!(arrs->rows)) {
		free(spm);
		free((arrs->values));
		free((arrs->cols));
		free(arrs);
		return 0;
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
node* create_list(const double* row, int n) {
	node* head=NULL, * tail=NULL;
	const double* p = row;
	char first = 1;
	int counter = 0;
	if (print) {
		printf("%d create_list\n", n);
	}
	for (; p < row + n; p++) {
		if (*p != 0) {

			if (first) {
				head = (node*)malloc(sizeof(node));
				if (!head) {
					return NULL;
				}
				head->val = *p;
				head->index = counter;
				tail = head;
				first = 0;
				if (print) {
					printf("first succeeded");
				}
			}
			else {
				tail->next = (node*)malloc(sizeof(node));
				tail = tail->next;
				if (!tail) {
					return NULL;
				}
				tail->val = *p;
				tail->index = counter;
			}
		}
		counter++;
	}
	if (tail != NULL) {
		tail->next = NULL;
	}
	if (print) {
		printf("\n");
		print_list(head);
		printf("\n");
	}
	return head;
}
/*TODO: maybe its not neccessary*/
void addRow_list(spmat* A, const double* row, int i) {
	if (print) {
		printf("addRow_list\n");
	}
	*((((list*)(A->handle))->rows)+i) = create_list(row, A->n);
}
/*TODO: maybe replace with iterative version*/
void delete_list(node* l) {
	if (l != NULL) {
		if (print) {
			printf("delete node %f\n", l->val);
		}
		delete_list(l->next);
		free(l);
	}
}

void free_list(spmat* A) {
	node** l=((list*)(A->handle))->rows;
	if (print) {
		printf("free_list\n");
	}
	for (; l < (((list*)(A->handle))->rows) + A->n; l++) {
		if (print) {
			print_list(*l);
		}
		delete_list(*l);
	}
	if (print) {
		printf("pass delete\n");
	}
	free(((list*)(A->handle))->rows);
	free(A->handle);
	free(A);
}


void mult_list(const spmat* A, const double* v, double* result) {
	double sum;
	node** currRow=((list*)(A->handle))->rows;  /* current row*/
	node* currElem = *currRow;  /* current element*/
	double* currRes = result;  /* current result element*/
	int i;

	for (i = 0; i < A->n; i++) {
		sum = 0;
		while (currElem != NULL) {
			sum += (currElem->val) * (*(v + (currElem->index)));
			currElem = currElem->next;
		}
		*(currRes) = sum;
		currRes++;
		currRow++;
		currElem = *currRow;
	}
}
/* return  1 if all the elements in row are 0, and 0 otherwise*/
char is_row_empty_list(const spmat* A , int row){
	node** rows=((list*)(A->handle))->rows;
	return *(rows+row) == NULL ? 1:0;
}

void sum_rows_list(const spmat* A ,int row , double* row2add , double* result){
	node** rows=((list*)(A->handle))->rows;
	node* req_row = *(rows+row);
	int i, index;
	int end = result + A->n;
	if (req_row == NULL){ /*req_row is all zeroes*/
		while(result < end){
			*result = *row2add;
			result++;
		}
	}
	else{/*req_row is not all zeroes*/
		for(i=0; i<A->n; i++){
			if(i == req_row->index){
				*result = req_row->val + *row2add;
				req_row++;
			}
			else{
				*result = *row2add;
			}
			result++;
		}
	}
	return;
}

double sum_of_largest_column_list(const spmat* A, double* sum_col){
	double curr_max;
	int i,j;
	double *temp_sum_cols;
	node *curr_row ,*temp_row;
	node** rows=((list*)(A->handle))->rows;
	curr_row = *rows;
	/*computes the sum of each colomn*/
	for(i=0; i<A->n; i++){
		if(curr_row != NULL){
			temp_row = curr_row;
			temp_sum_cols = sum_col;
			j=0;
			while(temp_row != NULL){
				if(j==temp_row->index){
					*temp_sum_cols = *temp_sum_cols + temp_row->val;
					temp_row++;
				}
				j++;
				temp_sum_cols++;
			}
		}
		curr_row++;
	}
	/*computes the max sum*/
	temp_sum_cols = sum_col;
	for(i=0; i<A->n; i++){
		curr_max = (*temp_sum_cols > curr_max) ? *temp_sum_cols:curr_max;
	}
	return curr_max;
}

spmat* spmat_allocate_list(int n) {
	list* l;
	spmat* spm;
	l =(list*) malloc(sizeof(list));
	if (!l) {
		return 0;
	}
	l->rows=(node**) malloc(n*sizeof(node*));
	if (!(l->rows)) {
		free(l);
		return 0;
	}
	spm = (spmat*)malloc(sizeof(spmat));
	if (!spm) {
		free(l->rows);
		free(l);
		return 0;
	}
	spm->n = n;
	spm->handle = l;
	spm->add_row = addRow_list;
	spm->free = free_list;
	spm->mult = mult_list;
	spm->sum_rows = sum_rows_list;
	spm->sum_of_largest_column = sum_of_largest_column_list;
	return spm;
}
/*--------------------------------------------------------------------------------------------------------------------*/