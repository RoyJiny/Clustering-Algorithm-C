#include <stdio.h>
#include <stdlib.h>
#include "algo.h"

Error algo_2(spmat *A, int *degrees, double *eigen_vector, group *g, group *g1, group *g2, double B_1norm, double M)
{
    int i, g_count;
    char stop = 0, *g_members;
    double *B_row;
    double *B_g_row, *runner1, *runner2, *runner3, *mult_vector;
    double modularity_value, eigen_value, *s, magnitude;
    double *unnormalized_eigen_vector;
    Error error;

    /*------------------------ALLOCATIONS------------------------*/
    B_g_row = (double *)malloc((g->size) * sizeof(double));
    if (!B_g_row)
    {
        printf("alocation failed in gb row");
        return ALLOCATION_FAILED;
    }

    B_row = (double *)malloc((A->n) * sizeof(double));
    if (!B_row)
    {
        printf("alocation failed in b row");
        return ALLOCATION_FAILED;
    }

    unnormalized_eigen_vector = (double *)malloc((g->size) * sizeof(double));
    if (!unnormalized_eigen_vector)
    {
        printf("alocation failed in unmrlzd ev");
        return ALLOCATION_FAILED;
    }

    mult_vector = (double *)malloc((g->size) * sizeof(double));
    if (!mult_vector)
    {
        printf("alocation failed in mult");
        return ALLOCATION_FAILED;
    }

    s = (double *)malloc((g->size) * sizeof(double));
    if (!s)
    {
        printf("alocation failed in s");
        return ALLOCATION_FAILED;
    }

    /*---------------------power iteration-------------------------*/
    while (!stop)
    {
        runner1 = mult_vector;
        g_members = g->members;
        g_count = 0;

        for (i = 0; i < A->n; i++)
        {
            /*do only if the vertex is in g*/
            if (*g_members)
            {
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

    /*--------------------decide the right partition-----------------*/
    if (!IS_POSITIVE(eigen_value))
    {
        printf("g is indivisible\n");
        printf("eigen value is not non-positive, value: %f\n", eigen_value);
        /*TODO : update g1 & g2*/
        return NONE;
    }

    eigen2s(eigen_vector, g1, g2, s, A->n);

    /*computing the modularity value*/
    g_members = g->members;
    runner1 = mult_vector;
    for (i = 0; i < A->n; i++)
    {
        /*do only if the vertex is in g*/
        if (*g_members)
        {
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
    {
        printf("g is indivisible\n");
        /*TODO : update g1 & g2*/
        return NONE;
    }

    /*remember: if there is a division of g then "eigen2s already computed the division"*/

    free(s);
    free(mult_vector);
    free(B_g_row);
    free(B_row);
    free(unnormalized_eigen_vector);
    return NONE;
}

/*P and O are the initial division*/
Error algo_3(spmat *A, int *degrees, group_set *P, group_set *O, int nof_vertex)
{
    group *g, *g1, *g2;
    double B_1norm, *init_vector, M, *runner_d;
    int *runner_i;
    Error error;

    /*init vector is set to the maximum size*/
    init_vector = malloc(nof_vertex * sizeof(double));
    if (!init_vector)
    {
        return ALLOCATION_FAILED;
    }
    /*------------------------randomize initial vector--------------------*/
    for (runner_d = init_vector; runner_d < init_vector + nof_vertex; runner_d++)
    {
        *runner_d = rand();
    }

    /*---------------------------computing M-----------------------------*/
    for (runner_i = degrees; runner_i < degrees + A->n; runner_i++)
    {
        M += *runner_i;
    }

    /*-------------------compute the 1norm for initial B------------------*/
    B_1norm = compute_1norm(A, degrees, M);
    printf("the 1norm for B is: %f\n", B_1norm);

    /*-----------------------------run-------------------------------------*/
    while (!(P->is_empty(P)))
    {
        /*-------------------------Allocations-------------------------*/
        /*we allocate g1,g2 the maximum possible size, although in future runs they won't*/
        /*acutally use the full allocated size*/
        g1 = (group *)malloc(sizeof(group));
        if (!g1)
        {
            return ALLOCATION_FAILED;
        }
        g1->members = (char *)malloc(nof_vertex * sizeof(char));
        if (!(g1->members))
        {
            return ALLOCATION_FAILED;
        }

        g2 = (group *)malloc(sizeof(group));
        if (!g2)
        {
            return ALLOCATION_FAILED;
        }
        g2->members = (char *)malloc(nof_vertex * sizeof(char));
        if (!(g2->members))
        {
            return ALLOCATION_FAILED;
        }

        g = P->pop(P);
        printf("size of g: %d, the value is:\n", g->size);
        print_group(g, nof_vertex);
        error = algo_2(A, degrees, init_vector, g, g1, g2, B_1norm, M);
        printf("finished algo 2 run\n\n");

        if (handle_errors(error, "algo_2"))
        {
            return error;
        }

        /*modularity maximization*/

        if (g1->size == 0 || g2->size == 0)
        {
            O->push(O, g);
            printf("one group is empty, pushing g to O, freeing g1,g2\n");
            /*g1,g2 won't be used anymore*/
            free(g1->members);
            free(g2->members);
            free(g1);
            free(g2);
        }
        else
        {
            printf("groups are not empty, pushing to P and O\n");
            g1->size == 1 ? O->push(O, g1) : P->push(P, g1);
            g2->size == 1 ? O->push(O, g2) : P->push(P, g2);
        }
    }

    /*free(g1->members);
    free(g1);
    free(g2->members);
    free(g2);*/
    free(init_vector);
    return NONE;
}
