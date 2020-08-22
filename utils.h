#include "spmat.h"
typedef enum
{
    NONE,
    ALLOCATION_FAILED,
    READ_FAILED
} Error;

char handle_errors(Error error, char *name);

Error read_input(FILE *input, spmat *A, int *degree, int nof_vertex);