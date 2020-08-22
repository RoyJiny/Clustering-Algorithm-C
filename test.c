#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "utils.h"

void test(spmat *mat, FILE *compare)
{
    int n, i, *row;
    spmat *A = spmat_allocate_list(n);
    if (!A)
    {
        printf("error allocating sparse matrix");
        return;
    }
    if (fread(&n, sizeof(int), 1, compare) != 1)
    {
        printf("reading error");
        return;
    }
    if (fread(&n, sizeof(int), 1, compare) != 1)
    {
        printf("reading error");
        return;
    }
    row = (int *)malloc(n * sizeof(int));
    if (!row)
    {
        printf("error allocating row");
        return;
    }

    for (i = 0; i < n; i++)
    {
        if (fread(row, sizeof(int), n, compare) != n)
        { /* read row*/
            printf("reading error");
            return;
        }
        A->add_row(A, row, i);
    }

    if (A->equal2(mat, A))
    {
        printf("success");
    }
    else
    {
        printf("failed");
    }
}