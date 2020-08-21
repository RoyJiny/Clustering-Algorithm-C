#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

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

void power_iteration(spmat *A, double *vector, double epsilon, double c_1norm)
{
	int stop = 0;
	double magnitude;
	double *mul_result;
	double *ptr_v, *ptr_r; /*pointers to vector and mul_result*/
	mul_result = (double *)malloc((A->n) * sizeof(double));
	if (!mul_result)
	{
		return;
	}
	while (stop == 0)
	{

		A->mult(A, vector, mul_result); /*result=A*vector*/
		magnitude = sqrt(dot_product(mul_result, mul_result, A->n));

		stop = 1;
		ptr_r = mul_result;
		for (ptr_v = vector; ptr_v < vector + A->n; ptr_v++)
		{
			if (fabs(*ptr_v - ((*(ptr_r)) / magnitude)) > epsilon)
			{
				stop = 0;
			}
			*ptr_v = (*(ptr_r)) / magnitude; /*updating vector to the next vector*/
			ptr_r++;
		}
	}
	free(mul_result);
}

void modularity_maximization(spmat *B_g, int *g)
{
	/*g is a single vector that can represnt g1,g2 - devision to 2 groups*/
	int size = B_g->n;
	/*save a linked list of all verticies for the beginning of every iteration*/

	/*while modularity of the best state keeps increasing:*/

		/*get a copy of the verticies list*/
		/*save the initial modularity*/
		/*create the new g1,g2 to represent the current state*/
			/*(at first the original g1,g2, after that - the best state from the previous process)*/
		/*create the new g1,g2 to represent the new best state (same as above)*/

		/*while verticies list is not empty: (to make sure that each vertex moves once)*/

			/*create variables for the best move: vertex, modularity value, division(g1,g2)*/

			/*while we havent reached the end of the list*/

				/*move the current vertex to the opposite group*/
				/*modify g1,g2*/
				/*calculate the new modularity*/

				/*if the new modularity is better than the highest so far:*/
					/*save the new modularity value as the highest so far this iteration*/
					/*save the new g1,g2 as the highest g1,g2*/
					/*save the respective vertex as part of the best move*/

				/*compare with best state overall, and update if better*/

				/*move the head of verticies list forward*/

			/*remove the vertex of the best move from the list*/
}