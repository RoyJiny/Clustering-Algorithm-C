
#ifndef PARAM
#define PARAM
/*#include <stdarg.h>*/

#define IS_POSITIVE(x) ((x) > 0.00001)

#define alloc(var,type,size,func_name,var_name,return_if_fail)\
do{\
    (var = (type*)malloc((size)*(sizeof(type))));\
    if(!var)\
    {\
        print_errors(ALLOCATION_FAILED,var_name,func_name);\
        return return_if_fail;\
    }\
}while(0)

#define handle_errors(error)\
do{\
    if(error != NONE && error != INDIVISIBLE)\
    {\
        return error;\
    }\
}while(0)

#define exit_if_error(error) if(error != NONE)\
{\
exit(5);\
}

#define MAX_NOF_ITERATIONS(size) ((size)*(size)*(size)*(size) + (size)*100)

/*#define FREE_ALL(...)\
do {\
    int i=0;\
    void *pta[] = {__VA_ARGS__};\
    for(i=0; i < sizeof(pta)/sizeof(void*); i++)\
    {\
        free(pta[i]);\
    }\
} while(0)*/

typedef enum
{
	NONE,
	ALLOCATION_FAILED,
	READ_FAILED,
	WRITE_FAILED,
	DIVISION_BY_ZERO,
	INDIVISIBLE,
    ENDLESS_LOOP
} Error;



#endif