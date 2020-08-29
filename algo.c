#include <stdio.h>
#include <stdlib.h>
#include "algo.h"

Error algo_2(spmat *A, int *degrees, double *eigen_vector, group *g, group *g1, group *g2)
{
    int i, *temp_i, g_count;
    char stop = 0, *g_members;
    double *B_row, B_1norm;
    double *B_g_row, *runner1, *runner2, *runner3, *mult_vector;
    double M, modularity_value, eigen_value, *s, magnitude;
    double *unnormalized_eigen_vector;
    Error error;

    /*------------------------ALLOCATIONS------------------------*/

    B_g_row = (double *)malloc((g->size) * sizeof(double));
    if (!B_g_row)
    {
        return ALLOCATION_FAILED;
    }

    B_row = (double *)malloc((A->n) * sizeof(double));
    if (!B_row)
    {
        return ALLOCATION_FAILED;
    }

    unnormalized_eigen_vector = (double *)malloc((g->size) * sizeof(double));
    if (!unnormalized_eigen_vector)
    {
        return ALLOCATION_FAILED;
    }

    mult_vector = (double *)malloc((g->size) * sizeof(double));
    if (!mult_vector)
    {
        return ALLOCATION_FAILED;
    }

    s = (double *)malloc((g->size) * sizeof(double));
    if (!s)
    {
        return ALLOCATION_FAILED;
    }

    /*---------------------computing M----------------------------*/
    for (temp_i = degrees; temp_i < degrees + A->n; temp_i++)
    {
        M += *temp_i;
    }

    /*---------------compute the 1norm for initial B----------------*/
    B_1norm = compute_1norm(A, degrees, M);
    printf("the 1norm for B is: %f\n", B_1norm);

    /*---------------------power iteration-------------------------*/
    while (!stop)
    {
        runner1 = mult_vector;
        g_members = g->members;
        g_count = 0;

        for (i = 0; i < A->n; i++)
        {
            if (*g_members)
            { /*curr vertex in g*/
                error = compute_modularity_matrix_row(A, i, g, degrees, M, B_g_row);
                if (error != NONE)
                {
                    printf("failed in compute_modularity_matrix_row\n");
                    return error;
                }

                B_g_row[g_count] += B_1norm;
                *runner1 = dot_product(B_g_row, eigen_vector, g->size);
                runner1++;
                g_count++;
            }
            g_members++;
        }
        magnitude = sqrt(dot_product(mult_vector, mult_vector, g->size));
        stop = 1;
        runner2 = eigen_vector;
        runner3 = unnormalized_eigen_vector;
        for (runner1 = mult_vector; runner1 < mult_vector + g->size; runner1++)
        {
            if (IS_POSITIVE(fabs(*runner2 - (*runner1 / magnitude))))
            {
                stop = 0;
            }
            *runner2 = *runner1 / magnitude;
            *runner3 = *runner1;
            runner2++;
            runner3++;
        }
    }

    printf("the eigen vector is:\n");
    print_vector(eigen_vector, g->size);

    /*---------------------computing leading eigen value-------------*/
    eigen_value = calculate_eigen_value(A, unnormalized_eigen_vector, g, degrees, M, B_g_row, B_1norm);
    printf("the eigen value is: %f\n", eigen_value);
    /*TODO: something in the nirmool (see forum), and need to sub B_1norm*/

    /*--------------------decide the right partition-----------------*/
    if (!IS_POSITIVE(eigen_value))
    {
        printf("g is indivisible\n");
        printf("eigen value is not visible, value: %f\n", eigen_value);
        /*TODO : update g1 & g2*/
        return NONE;
    }

    eigen2s(eigen_vector, g1, g2, s, A->n);

    /*computing the modularity value*/
    g_members = g->members;
    runner1 = mult_vector;
    for (i = 0; i < A->n; i++)
    {
        if (*g_members)
        { /*curr vertex in g*/
            error = compute_modularity_matrix_row(A, i, g, degrees, M, B_g_row);
            if (error != NONE)
            {
                printf("failed in compute_modularity_matrix_row\n");
                return error;
            }
            *runner1 = dot_product(B_g_row, s, g->size);
            runner1++;
        }
        g_members++;
    }
    modularity_value = 0.5 * dot_product(mult_vector, s, g->size);
    printf("modularity value is: %f\n", modularity_value);
    if (!IS_POSITIVE(modularity_value))
    { /*g is indivisible*/
        printf("g is indivisible\n");
        /*TODO : update g1 & g2*/
        return NONE;
    }

    /*remember: if there is a division of g then "eigen2s already computed the division"*/

    free(s);
    free(mult_vector);
    free(B_g_row);
    return NONE;
}
