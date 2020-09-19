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

	/*push a group to the set*/
	void (*push)(struct _group_set *s, group *g);

	/*pop a group from the set*/
	group *(*pop)(struct _group_set *s);

	/*get top group of the set*/
	group *(*top)(const struct _group_set *s);

	/*check if there are any groups in the set*/
	char (*is_empty)(const struct _group_set *s);

	/*free the set itself and all of its groups*/
	void (*free_set)(struct _group_set *s);

} group_set;

/*allocate a group set with no groups*/
group_set *allocate_group_set();

void print_group(group *g);

#endif