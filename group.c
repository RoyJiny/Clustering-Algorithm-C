#include "group.h"
#include <stdlib.h>
#include <stdio.h>

void push(group_set *s, group *g)
{
    group_node *new_node;
    new_node = (group_node *)malloc(sizeof(group_node));
    if (!new_node)
    {
        printf("malloc failed on new_node");
    }
    new_node->value = g;
    new_node->next = s->first;
    s->first = new_node;
    s->size++;
}

group *pop(group_set *s)
{
    group *g;
    group_node *node;
    g = s->first->value;
    node = s->first;
    s->first = s->first->next;
    s->size--;
    free(node);
    return g;
}

group *top(const group_set *s)
{
    return s->first->value;
}

char is_empty(const group_set *s)
{
    return (s->size) <= 0;
}

void free_set(group_set *s)
{
    group_node *node;
    while (s->first != NULL)
    {
        node = s->first;
        s->first = s->first->next;
        free(node->value->members);
        free(node->value);
        free(node);
        s->size--;
    }
    free(s);
}

group_set *allocate_group_set()
{
    group_set *s;
    s = (group_set *)malloc(sizeof(group_set));
    if (!s)
    {
        printf("allocation faild on s");
    }
    s->size = 0;
    s->first = NULL;
    s->top = top;
    s->free_set = free_set;
    s->is_empty = is_empty;
    s->pop = pop;
    s->push = push;

    return s;
}
