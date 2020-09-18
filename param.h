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

void print_errors(Error error, char *name, char *func);

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