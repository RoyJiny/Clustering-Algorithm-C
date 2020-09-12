#ifndef _LIST_H
#define _LIST_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef struct _dynamic_node
{
    int vertex;
    struct _dynamic_node *next;
} dynamic_node;

typedef struct _dynamic_list
{
    dynamic_node *head;
    int size;
} dynamic_list;

void print_dynamic_list(dynamic_list *list)
{
    dynamic_node *runner = list->head;
    while(runner){
        printf("%d->",runner->vertex);
        runner  = runner->next;
    }
    printf("\n\n");
}

char create_dynamic_list(dynamic_list *list, int size)
{
    dynamic_node *head = NULL, *tail = NULL;
    int i = 1;
    if (size == 0)
    {
        return 0;
    }
    head = (dynamic_node *)malloc(sizeof(dynamic_node));
    if (!head)
    {
        printf("alocation failed in gb row");
        return 0;
    }
    head->vertex = 0;
    head->next = NULL;
    tail = head;
    for (; i < size; i++)
    {
        tail->next = (dynamic_node *)malloc(sizeof(dynamic_node));
        tail = tail->next;
        if (!tail)
        {
            printf("alocation failed in gb row");
            return 0;
        }
        tail->vertex = i;
        tail->next = NULL;
    }
    if(tail != NULL){
        tail->next = NULL;
    }
    list->head = head;
    list->size = size;
    print_dynamic_list(list);
    return 1;
}

void delete_node_by_index(dynamic_list *list, int index)
{
    dynamic_node *next, *runner = list->head;
    if (runner->vertex == index)
    {
        next = runner->next;
        free(runner);
        /*runner = next;*/
        list->head = next;
        list->size--;
        return;
    }
    while (runner->next != NULL)
    {
        if (runner->next->vertex == index)
        {
            next = runner->next;
            runner->next = runner->next->next;
            free(next);
            list->size--;
            return;
        }
        runner = runner->next;
    }
    printf("didn't find the node of index %d\n", index);
}

#endif