#ifndef _GROUP_H
#define _GROUP_H

#include "param.h"

typedef struct _group
{
	int *members; /*i in g iff members[i] == 1*, size is always nof_vertex*/
	int size;
} group;

typedef struct _group_node
{
	group *value;
	struct _group_node *next;
} group_node;

typedef struct _group_set
{
	group_node *first;
	int size;

	void (*push)(struct _group_set *s, group *g);

	group *(*pop)(struct _group_set *s);

	group *(*top)(const struct _group_set *s);

	char (*is_empty)(const struct _group_set *s);

	void (*free_set)(struct _group_set *s);

} group_set;

group_set *allocate_group_set();

#endif