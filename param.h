
#ifndef PARAM
#define PARAM
#include <varargs.h>

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
    exit(error);\
}

#define MAX_NOF_ITERATIONS(size) ((size)*(size)*(size)*(size) + (size)*100)

/*#define FREE_ALL(va_alist)\
do {\
    va_list list;
    int i=0;\
    void *p;
    va_start(va_list);
    for(i=0; i < sizeof(pta)/sizeof(void*); i++)\
    {\
        free(pta[i]);\
    }\
} while(0)*/

/*#define FREE_ALL( ... ) Free( &free_stop , __VA_ARGS__ , &free_stop )

int free_stop;

void Free( void* point , ... )  
{
    va_list list;
    void* p;
    if( !point )
        return;
    va_start( list , point );

   p = va_arg( list ,void*);
    while( p != point ) 
    {
        free( p ) ;
        p = va_arg( list , void* );
    }
    va_end( list ) ;
}*/

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