#ifndef _ALGO_H
#define _ALGO_H
#include "utils.h"

/*i*/
DivisionResult algo_2(spmat *A, int *degrees, double *init_vector, group *g, group *g1, group *g2, double B_1norm, double M);

void algo_3(spmat *A,int *degrees, group_set *P, group_set *O ,int nof_vertex);
#endif
