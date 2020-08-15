#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "utils.h"


double dot_product(double* row1, double* row2, int size)
{
	double sum = 0;
	double* end = row1 + size;
	for (; row1 < end;) {
		sum += (*row1) * (*row2);
		row1++;
		row2++;
	}

	return sum;
}

void power_iteration(spmat* A, double* vector, double epsilon) {
	int stop = 0;
	double magnitude;
	double* mul_result;
	double* ptr_v, * ptr_r;  /*pointers to vector and mul_result*/
	mul_result= (double*)malloc((A->n) * sizeof(double));
	if (!mul_result) {
		return;
	}
	while (stop == 0) {
		
		A->mult(A, vector, mul_result); /*result=A*vector*/
		magnitude = sqrt(dot_product(mul_result, mul_result, A->n));

		stop = 1;
		ptr_r = mul_result;
		for (ptr_v = vector; ptr_v < vector + A->n; ptr_v++) {
			if (fabs(*ptr_v - ((*(ptr_r)) / magnitude)) > epsilon) {
				stop = 0;
			}
			*ptr_v = (*(ptr_r)) / magnitude;  /*updating vector to the next vector*/
			ptr_r++;
		}
	}
	free(mul_result);
}