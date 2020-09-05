#include <stdio.h>
#include <stdlib.h>
#include "algo.h"
#include "list.h"

extern int run_num;

Error modularity_maximization(spmat *A, int *degrees, double *s, double M, group *g)
{
    dynamic_list *unmoved;
    dynamic_node *node_runner;
    double max_score, Q0, new_score, max_improve, delta_Q;
    double *mult_vector, *B_g_row;
    int i, j, max_score_index, max_improve_index;
    double *improve, *improve_runner;
    int *indices, *indices_runner;
    char stop = 0;

    /*----------------------Alocations--------------------*/
    
    B_g_row = (double *)malloc(g->size * sizeof(double));
    if (!B_g_row)
    {
        return ALLOCATION_FAILED;
    }
    mult_vector = (double *)malloc(g->size * sizeof(double));
    if (!mult_vector)
    {
        return ALLOCATION_FAILED;
    }
    improve = (double *)malloc(g->size * sizeof(double));
    if (!improve)
    {
        return ALLOCATION_FAILED;
    }
    indices = (int *)malloc(g->size * sizeof(int));
    if (!indices)
    {
        return ALLOCATION_FAILED;
    }
    /*------------------------------------------------------*/

    while (!stop)
    {
        indices_runner = indices;
        improve_runner = improve;
        unmoved = allocate_dynamic_list(g->size);
        for (i = 0; i < g->size; i++)
        {
            Q0 = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector);
            node_runner = unmoved->head;

            /*first run*/
            
            *(s + node_runner->vertex) = -*(s + node_runner->vertex);
            max_score = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector) - Q0;
            max_score_index = node_runner->vertex;
            *(s + node_runner->vertex) = -*(s + node_runner->vertex);
            node_runner = node_runner->next;

            while (node_runner != NULL)
            {
                *(s + node_runner->vertex) = -*(s + node_runner->vertex);
                new_score = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector) - Q0;
                if (new_score > max_score)
                {
                    max_score = new_score;
                    max_score_index = node_runner->vertex;
                }
                *(s + node_runner->vertex) = -*(s + node_runner->vertex);
                node_runner = node_runner->next;
            }

            *(s + max_score_index) = -*(s + max_score_index);
            *indices_runner = max_score_index;
            if (i == 0)
            {
                *improve_runner = max_score;
                max_improve_index = 0;
                max_improve = max_score;
            }
            else
            {
                *improve_runner = *(improve_runner - 1) + max_score;
                if (*improve_runner > max_improve)
                {
                    max_improve = *improve_runner;
                    max_improve_index = i;
                }
            }
            delete_node_by_index(unmoved, max_score_index);
            indices_runner++;
            improve_runner++;
        }
        
        indices_runner = indices + (g->size - 1);
        for (i = g->size - 1; i > max_improve_index; i--)
        {
            j = *indices_runner;
            *(s + j) = -*(s + j);
            indices_runner--;
        }
        if (max_improve_index == g->size - 1)
        {
            delta_Q = 0;
        }
        else
        {
            delta_Q = max_improve;
        }

        if (!IS_POSITIVE(delta_Q))
        {
            stop = 1;
        }
        free(unmoved);
    }
    free(mult_vector);
    free(improve);
    free(indices);

    return NONE;
}

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
                /*g_count is the relative row number for the sub matrix of g*/
                error = compute_modularity_matrix_row(A, i, g, degrees, M, B_g_row, g_count);
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

    /*printf("the eigen vector is:\n");
    print_vector(eigen_vector, g->size);*/
    /*---------------------computing leading eigen value-------------*/
    eigen_value = calculate_eigen_value(A, unnormalized_eigen_vector, g, degrees, M, B_g_row, B_1norm);
    /*printf("the eigen value is: %f\n", eigen_value);*/

    /*--------------------decide the right partition-----------------*/
    if (!IS_POSITIVE(eigen_value))
    {
       /* printf("g is indivisible\n");
        printf("eigen value is not non-positive, value: %f\n", eigen_value);*/
        /*TODO : update g1 & g2*/
        free(mult_vector);
        free(B_row);
        free(unnormalized_eigen_vector);
        free(B_g_row);
        free(s);
        return INDIVISIBLE;
    }

    eigen2s(eigen_vector, g, s, A->n);

    /*printf("s before the max:\n");
    print_vector(s,g->size);*/

    modularity_maximization(A, degrees, s, M, g);

    /*printf("s after the max:\n");
    print_vector(s,g->size);*/

    construct_g1g2(g, s, g1, g2, A->n);

    /*computing the modularity value*/
    modularity_value = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector);
    /*printf("modularity value is: %f\n", modularity_value);*/

    if (!IS_POSITIVE(modularity_value))
    {
        /*printf("g is indivisible\n");*/
        free(mult_vector);
        free(B_row);
        free(unnormalized_eigen_vector);
        free(B_g_row);
        free(s);
        return INDIVISIBLE;
    }

    /*remember: if there is a division of g then "eigen2s already computed the division"*/

    free(mult_vector);
    free(B_row);
    free(unnormalized_eigen_vector);
    free(B_g_row);
    free(s);
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

    /*-----------------------------run-------------------------------------*/
    while (!(P->is_empty(P)))
    {
        /*printf("-------------------------------------------------------------------\n");
        printf("starting a loop run of algorithm 3 (run %d)\n", run_num);*/
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
        /*printf("popped g. size of g: %d, the value is:\n", g->size);
        print_group(g, nof_vertex);*/

        error = algo_2(A, degrees, init_vector, g, g1, g2, B_1norm, M);

        if (handle_errors(error, "algo_2"))
        {
            return error;
        }

        if (g1->size == 0 || g2->size == 0 || error == INDIVISIBLE)
        {
            O->push(O, g);
            /*printf("one group is empty, pushing g to O, freeing g1,g2\n");*/
            /*g1,g2 won't be used anymore*/
            free(g1->members);
            free(g2->members);
            free(g1);
            free(g2);
        }
        else
        {
            /*printf("after algo_2\ng1 is:\n");
            print_group(g1, nof_vertex);
            printf("\ng2 is:\n");
            print_group(g2, nof_vertex);
            printf("groups are not empty, pushing to P and O\n");*/
            g2->size == 1 ? O->push(O, g2) : P->push(P, g2);
            g1->size == 1 ? O->push(O, g1) : P->push(P, g1);
            free(g->members);
            free(g);
        }
        run_num++;
    }

    /*free(g1->members);
    free(g1);
    free(g2->members);
    free(g2);*/
    free(init_vector);
    return NONE;
}
