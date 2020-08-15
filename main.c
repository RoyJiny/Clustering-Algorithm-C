#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "utils.h"


void print_vector(double* v,int n) {
	double* ptr = v;
	for (; ptr < v + n; ptr++) {
		printf("%f	", *ptr);
	}
}

int main(int argc, char* argv[]) {
	return 0;
}