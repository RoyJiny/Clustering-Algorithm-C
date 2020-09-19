#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "algo.h"
#include "list.h"

void modularity_maximization(spmat *A, int *degrees, double *s, double M, group *g, int *g1_counter)
{
    double max_score, new_score, max_improve, delta_Q;
    double *mult_vector, *B_g_row;
    int i, j, max_score_index, max_improve_index;
    double *improve, *improve_runner;
    double *d_pointer;
    int *indices, *indices_runner;
    char stop = 0;
    int current_vertex_index;
    int iteration_counter = 0, counter = 0;
    int *unmoved ,*unmoved_runner, *unmoved_last, *vertex2delete;
    int unmoved_size;

    alloc(unmoved,int,g->size,"modularity_maximization","unmoved");
    alloc(B_g_row,double,g->size,"modularity_maximization","B_g_row");
    alloc(mult_vector,double,g->size,"modularity_maximization","mult_vector");
    alloc(improve,double,g->size,"modularity_maximization","improve");
    alloc(indices,int,g->size,"modularity_maximization","indices");

    for(unmoved_runner = unmoved; unmoved_runner<unmoved+g->size; unmoved_runner++){
        *unmoved_runner = i;
        i++;
    }

    while (!stop)
    {
        iteration_counter++;
        indices_runner = indices;
        improve_runner = improve;
        i=0;
        unmoved_size = g->size;
        unmoved_last = unmoved + unmoved_size -1;

        for (i = 0; i < g->size; i++)
        {
            unmoved_runner = unmoved;
            
            /*first run*/
            current_vertex_index = *unmoved_runner;
            d_pointer = s + current_vertex_index;
            *d_pointer = -(*d_pointer);

            compute_score(A, *(g->members + current_vertex_index), current_vertex_index, g, s, M, degrees, &max_score, B_g_row);

            max_score_index = current_vertex_index;
            *d_pointer = -(*d_pointer);
            vertex2delete = unmoved_runner;
            unmoved_runner++;

            for(j=1; j<unmoved_size; j++){
                counter++;
                current_vertex_index = *unmoved_runner;
                d_pointer = s + current_vertex_index;
                *d_pointer = -(*d_pointer);

                compute_score(A, *(g->members + current_vertex_index), current_vertex_index, g, s, M, degrees, &new_score, B_g_row);
                if (new_score > max_score)
                {
                    max_score = new_score;
                    max_score_index = current_vertex_index;
                    vertex2delete = unmoved_runner;
                   
                }
                *d_pointer = -(*d_pointer);
                unmoved_runner++;
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
            if (unmoved_size > 0){
                *vertex2delete = *unmoved_last;
                unmoved_size--;
                *unmoved_last = unmoved_size; /*update for next while iteration*/
                if (unmoved_size) unmoved_last--;
            }
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
        printf("delta q = %f\n",delta_Q); /*TODO:remove*/
    }

    free(mult_vector);
    free(improve);
    free(indices);
    free(unmoved);
    free(B_g_row);
}

DivisionResult algo_2(spmat *A, int *degrees, double *eigen_vector, group *g, group *g1, group *g2, double B_1norm, double M)
{
    int i, g1_count;
    char stop = 0;
    int *g_members;
    double *B_g_row, *runner1, *runner2, *mult_vector;
    double modularity_value, eigen_value, *s, magnitude;
    double row_sum = 0, *row_result, *A_sizes, *runner3, *runner4;
    int iteration_counter = MAX_NOF_ITERATIONS(g->size);
    time_t start ;

    alloc(B_g_row,double,g->size,"algo_2","B_g_row");
    alloc(mult_vector,double,g->size,"algo_2","mult_vector");
    alloc(row_result,double,g->size,"algo_2","row_result");
    alloc(A_sizes,double,g->size,"algo_2","A_sizes");
    alloc(s,double,g->size,"algo_2","s");

    /*---------------------power iteration-------------------------*/
    start = clock();
    while (!stop)
    {
        iteration_counter--;
        runner1 = mult_vector;
        g_members = g->members;
        runner3 = A_sizes;
        magnitude = 0;
        A->mult(A, eigen_vector, B_g_row, A_sizes, g);
        for (i = 0; i < g->size; i++)
        {
            row_sum = compute_D_row(*g_members, g, degrees, M, row_result);
            *runner3 += row_sum;
            *runner1 = dot_product(row_result, eigen_vector, g->size);
            runner1++;
            g_members++;
            runner3++;
        }
        runner1 = mult_vector;
        runner2 = B_g_row;
        runner3 = eigen_vector;
        runner4 = A_sizes;
        for (i = 0; i < g->size; i++)
        {
            *runner1 += *runner2 + (B_1norm + *(runner4))*(*runner3);
            magnitude += (*runner1) * (*runner1);
            runner1 ++;
            runner2 ++;
            runner3 ++;
            runner4 ++;
        }
        magnitude = sqrt(magnitude);
        stop = 1;
        runner2 = eigen_vector;
        for (runner1 = mult_vector; runner1 < mult_vector + g->size; runner1++)
        {
            if (IS_POSITIVE(fabs(*runner2 - (*runner1 / magnitude))))
            {
                stop = 0;
            }
            *runner2 = *runner1 / magnitude;
            runner2++;
        }
        if(!iteration_counter){
            handle_errors(ENDLESS_LOOP,"algo_2","power iteration");
        }
    }
    printf("the size of g -%d\n", g->size); /*TODO-remove*/
    printf("done power interation run -%d\n", MAX_NOF_ITERATIONS(g->size) - iteration_counter);

    calculate_eigen_value(A, eigen_vector, g, degrees, M, B_g_row, B_1norm, &eigen_value);

    if (!IS_POSITIVE(eigen_value))
    {
        free(mult_vector);
        free(B_g_row);
        free(s);
        free(g1);
        free(g2);
        return INDIVISIBLE;
    }

    g1_count = eigen2s(eigen_vector, g, s);

    start = clock();
    modularity_maximization(A, degrees, s, M, g, &g1_count);
    printf("done max- %ld\n", (clock() - start) / CLOCKS_PER_SEC); /*TODO-remove*/

    compute_modularity_value(A, g, degrees, s, M, B_g_row, mult_vector, &modularity_value);

    if (!IS_POSITIVE(modularity_value))
    {
        free(mult_vector);
        free(B_g_row);
        free(s);
        free(g1);
        free(g2);
        return INDIVISIBLE;
    }

    construct_g1g2(g, s, g1, g2, g1_count);

    free(mult_vector);
    free(B_g_row);
    free(s);
    return DIVISIBLE;
}

/*P and O are the initial division*/
void algo_3(spmat *A, int *degrees, group_set *P, group_set *O, int nof_vertex)
{
    group *g, *g1, *g2;
    double B_1norm, *init_vector, M, *runner_d;
    int *runner_i;
    DivisionResult is_divisible;

    /*initial vector*/
    alloc(init_vector,double,nof_vertex,"algo_3","init_vector");
    for (runner_d = init_vector; runner_d < init_vector + nof_vertex; runner_d++)
    {
        *runner_d = rand();
    }
    
    for (runner_i = degrees; runner_i < degrees + A->n; runner_i++)
    {
        M += *runner_i;
    }

    g = P->top(P);
    compute_1norm(A, g, degrees, M, &B_1norm);

    /*run*/
    while (!(P->is_empty(P)))
    {
        alloc(g1,group,1,"algo_3","g1");
        alloc(g2,group,1,"algo_3","g2");

        g = P->pop(P);
       
        is_divisible = algo_2(A, degrees, init_vector, g, g1, g2, B_1norm, M);
        if (g1->size == 0 || g2->size == 0 || is_divisible == INDIVISIBLE)
        {
            O->push(O, g);
        }
        else
        {
            g2->size == 1 ? O->push(O, g2) : P->push(P, g2);
            g1->size == 1 ? O->push(O, g1) : P->push(P, g1);
            free(g->members);
            free(g);
        }
    }

    free(init_vector);
}