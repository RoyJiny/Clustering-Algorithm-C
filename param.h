#ifndef PARAM
#define PARAM

#define IS_POSITIVE(x) ((x) > 0.00001)

typedef enum
{
    DIVISIBLE,
    INDIVISIBLE
} DivisionResult;

typedef enum
{
	NONE,
	ALLOCATION_FAILED,
	READ_FAILED,
	WRITE_FAILED,
	DIVISION_BY_ZERO,
    ENDLESS_LOOP
} Error;

void print_errors(Error error, char *name, char *func)
{
	switch (error)
	{
	case ALLOCATION_FAILED:
		printf("[%s]: allocation failed on: %s\n", func, name);
		return;
	case READ_FAILED:
		printf("[%s]: read failed on %s\n", func, name);
		return;
	case DIVISION_BY_ZERO:
		printf("[%s]: division by zero, %s is zero\n", func, name);
		return;
	case WRITE_FAILED:
		printf("[%s]: write failed on: %s\n", func, name);
		return;
	case ENDLESS_LOOP:
		printf("[%s]: suspision of an endless loop in: %s\n",func,name);
		return;
	default:
		return;
	}
}

#define alloc(var,type,size,func_name,var_name)\
do{\
    (var = (type*)malloc((size)*(sizeof(type))));\
    if(!var)\
    {\
        handle_errors(ALLOCATION_FAILED,func_name,var_name);\
    }\
}while(0)

#define handle_errors(error,func_name,var_name)\
do{\
    if(error != NONE)\
    {\
        print_errors(error,var_name,func_name);\
        exit(4);\
    }\
}while(0)

#define MAX_NOF_ITERATIONS(size) ((size)*(size)*(size)*(size) + (size)*100)

#endif