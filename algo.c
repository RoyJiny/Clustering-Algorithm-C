#include <stdio.h>
#include <stdlib.h>
#include "algo.h"
#include "list.h"

extern int run_num;

Error modularity_maximization(spmat *A, int *degrees, double *s, double M, group *g, int *g1_counter)
{
    dynamic_list *unmoved;
    dynamic_node *node_runner;
    double max_score, new_score, max_improve, delta_Q;
    double *mult_vector, *B_g_row;
    int i, j, max_score_index, max_improve_index;
    double *improve, *improve_runner;
    double *d_pointer;
    int *indices, *indices_runner;
    char stop = 0;
    int current_vertex_index;
    Error error;
    int iteration_counter = 0, counter = 0;

    /*----------------------Alocations--------------------*/

    B_g_row = (double *)malloc(g->size * sizeof(double));
    if (!B_g_row)
    {
        print_errors(ALLOCATION_FAILED, "B_g_row", "modularity_maximization");
        return ALLOCATION_FAILED;
    }
    mult_vector = (double *)malloc(g->size * sizeof(double));
    if (!mult_vector)
    {
        print_errors(ALLOCATION_FAILED, "mult_vector", "modularity_maximization");
        return ALLOCATION_FAILED;
    }
    improve = (double *)malloc(g->size * sizeof(double));
    if (!improve)
    {
        print_errors(ALLOCATION_FAILED, "improve", "modularity_maximization");
        return ALLOCATION_FAILED;
    }
    indices = (int *)malloc(g->size * sizeof(int));
    if (!indices)
    {
        print_errors(ALLOCATION_FAILED, "indices", "modularity_maximization");
        return ALLOCATION_FAILED;
    }
    unmoved = (dynamic_list *)malloc(sizeof(dynamic_list));
    if (!unmoved)
    {
        print_errors(ALLOCATION_FAILED, "unmoved", "modularity_maximization");
        return ALLOCATION_FAILED;
    }
    /*------------------------------------------------------*/

    while (!stop)
    {
        iteration_counter++;
        indices_runner = indices;
        improve_runner = improve;
        if (!create_dynamic_list(unmoved, g->size))
        {
            print_errors(ALLOCATION_FAILED, "create_dynamic_list", "modularity_maximization");
            return ALLOCATION_FAILED;
        }

        for (i = 0; i < g->size; i++)
        {

            /*error = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector, &Q0);
            if (error != NONE)
            {
                return error;
            }*/
            node_runner = unmoved->head;

            /*first run*/
            current_vertex_index = node_runner->vertex;
            *(s + current_vertex_index) = -*(s + current_vertex_index);

            /*error = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector, &max_score);*/
            error = compute_for_improved_score(A, *(g->members + current_vertex_index), current_vertex_index, g, s, M, degrees, &max_score);
            if (error != NONE)
            {
                return error;
            }

            max_score_index = current_vertex_index;
            *(s + current_vertex_index) = -*(s + current_vertex_index);
            node_runner = node_runner->next;

            while (node_runner != NULL)
            {
                counter++;
                current_vertex_index = node_runner->vertex;
                *(s + current_vertex_index) = -*(s + current_vertex_index);
                /*error = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector, &new_score);*/
                error = compute_for_improved_score(A, *(g->members + current_vertex_index), current_vertex_index, g, s, M, degrees, &new_score);
                if (error != NONE)
                {
                    return error;
                }

                if (new_score > max_score)
                {
                    max_score = new_score;
                    max_score_index = current_vertex_index;
                }
                *(s + current_vertex_index) = -*(s + current_vertex_index);
                node_runner = node_runner->next;
            }

            counter = 0;
            d_pointer = s + max_score_index;
            *(d_pointer) = -*(d_pointer);
            (*(d_pointer) == 1) ? (*g1_counter)++ : (*g1_counter)--;
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
            d_pointer = s + j;
            *(d_pointer) = -*(d_pointer);
            (*(d_pointer) == 1) ? (*g1_counter)++ : (*g1_counter)--;
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
    }
    free(mult_vector);
    free(improve);
    free(indices);
    free(unmoved);

    return NONE;
}

Error algo_2(spmat *A, int *degrees, double *eigen_vector, group *g, group *g1, group *g2, double B_1norm, double M)
{
    int i, g1_count;
    char stop = 0;
    int *g_members;
    double *B_row; /*TODO: check if it used*/
    double *B_g_row, *runner1, *runner2, *mult_vector;
    double modularity_value, eigen_value, *s, magnitude;
    /*double *unnormalized_eigen_vector;*/
    Error error;
    int iteration_counter = 0;
    /*------------------------ALLOCATIONS------------------------*/
    B_g_row = (double *)malloc((g->size) * sizeof(double));
    if (!B_g_row)
    {
        print_errors(ALLOCATION_FAILED, "B_g_row", "algo_2");
        return ALLOCATION_FAILED;
    }

    B_row = (double *)malloc((A->n) * sizeof(double));
    if (!B_row)
    {
        print_errors(ALLOCATION_FAILED, "B_row", "algo_2");
        return ALLOCATION_FAILED;
    }

    /*unnormalized_eigen_vector = (double *)malloc((g->size) * sizeof(double));
    if (!unnormalized_eigen_vector)
    {
        print_errors(ALLOCATION_FAILED, "unnormalized_eigen_vector", "algo_2");
        return ALLOCATION_FAILED;
    }*/

    mult_vector = (double *)malloc((g->size) * sizeof(double));
    if (!mult_vector)
    {
        print_errors(ALLOCATION_FAILED, "mult_vector", "algo_2");
        return ALLOCATION_FAILED;
    }

    s = (double *)malloc((g->size) * sizeof(double));
    if (!s)
    {
        print_errors(ALLOCATION_FAILED, "s", "algo_2");
        return ALLOCATION_FAILED;
    }

    /*---------------------power iteration-------------------------*/
    while (!stop)
    {
        iteration_counter++;
        runner1 = mult_vector;
        g_members = g->members;

        for (i = 0; i < g->size; i++)
        {
            error = compute_modularity_matrix_row(A, *g_members, g, degrees, M, B_g_row, i);
            if (error != NONE)
            {
                return error;
            }
            B_g_row[i] += B_1norm;
            *runner1 = dot_product(B_g_row, eigen_vector, g->size);
            runner1++;
            g_members++;
        }

        magnitude = sqrt(dot_product(mult_vector, mult_vector, g->size));
        stop = 1;
        runner2 = eigen_vector;
        /*runner3 = unnormalized_eigen_vector;*/
        for (runner1 = mult_vector; runner1 < mult_vector + g->size; runner1++)
        {
            if (IS_POSITIVE(fabs(*runner2 - (*runner1 / magnitude))))
            {
                stop = 0;
            }
            *runner2 = *runner1 / magnitude;
            /**runner3 = *runner1;*/
            runner2++;
            /*runner3++;*/
        }
    }
    /*printf("the eigen vector is:\n");
    print_vector(eigen_vector, g->size);*/
    /*---------------------computing leading eigen value-------------*/
    error = calculate_eigen_value(A, eigen_vector, g, degrees, M, B_g_row, B_1norm, &eigen_value);
    if (error != NONE)
    {
        return error;
    }
    /*error = calculate_eigen_value(A, unnormalized_eigen_vector, g, degrees, M, B_g_row, B_1norm, &c);
    if(error != NONE){
        return error;
    }
    if(IS_POSITIVE(c - eigen_value)){
        printf("fail: unnrz = %f, norm = %f\n", c, eigen_value);
    }*/
    /*printf("the eigen value is: %f\n", eigen_value);*/

    /*--------------------decide the right partition-----------------*/
    if (!IS_POSITIVE(eigen_value))
    {
        printf("g is indivisible\n");
        printf("eigen value is non-positive, value: %f\n", eigen_value);
        free(mult_vector);
        free(B_row);
        free(B_g_row);
        free(s);
        free(g1);
        free(g2);
        return INDIVISIBLE;
    }

    g1_count = eigen2s(eigen_vector, g, s);

    /*printf("s before the max:\n");
    print_vector(s,g->size);*/
    error = modularity_maximization(A, degrees, s, M, g, &g1_count);
    if (error != NONE)
    {
        return error;
    }
    /*printf("s after the max:\n");
    print_vector(s,g->size);*/
    /*computing the modularity value*/
    /*printf("A is:\n");
    A->print_matrix(A);
    printf("\ndegrees is:\n");
    print_vector_int(degrees, A->n);*/
    error = compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector, &modularity_value);
    if (error != NONE)
    {
        return error;
    }
    /*printf("modularity value is: %f\n", modularity_value);*/

    if (!IS_POSITIVE(modularity_value))
    {
        printf("g is indivisible\n");
        printf("modularity value is non-positive, value: %f\n", modularity_value);
        free(mult_vector);
        free(B_row);
        free(B_g_row);
        free(s);
        free(g1);
        free(g2);
        return INDIVISIBLE;
    }

    error = construct_g1g2(g, s, g1, g2, g1_count);
    if (error != NONE)
    {
        return error;
    }

    /*remember: if there is a division of g then "eigen2s already computed the division"*/

    free(mult_vector);
    free(B_row);
    /*free(unnormalized_eigen_vector);*/
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
        print_errors(ALLOCATION_FAILED, "init_vector", "algo_3");
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
    g = P->top(P); /*TODO: maybe check if P isn't empty*/
    error = compute_1norm(A, g, degrees, M, &B_1norm);
    if (error != NONE)
    {
        return error;
    }
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
            print_errors(ALLOCATION_FAILED, "g1", "algo_3");
            return ALLOCATION_FAILED;
        }
        /*g1->members = (int *)malloc(nof_vertex * sizeof(int));
        if (!(g1->members))
        {
            print_errors(ALLOCATION_FAILED, "g1->members", "algo_3");
            return ALLOCATION_FAILED;
        }*/

        g2 = (group *)malloc(sizeof(group));
        if (!g2)
        {
            print_errors(ALLOCATION_FAILED, "g1", "algo_3");
            return ALLOCATION_FAILED;
        }
        /*g2->members = (int *)malloc(nof_vertex * sizeof(int));
        if (!(g2->members))
        {
            print_errors(ALLOCATION_FAILED, "g2->members", "algo_3");
            return ALLOCATION_FAILED;
        }*/

        g = P->pop(P);
        /*printf("popped g. size of g: %d, the value is:\n", g->size);
        print_group(g, nof_vertex);*/

        error = algo_2(A, degrees, init_vector, g, g1, g2, B_1norm, M);
        if (error != NONE && error != INDIVISIBLE)
        {
            return error;
        }
        if (g1->size == 0 || g2->size == 0 || error == INDIVISIBLE)
        {
            O->push(O, g);
            /*printf("one group is empty, pushing g to O, freeing g1,g2\n");*/
            /*g1,g2 won't be used anymore*/
            /*free(g1->members);
            free(g2->members);
            free(g1);
            free(g2);*/
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
