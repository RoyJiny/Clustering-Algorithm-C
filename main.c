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
/*return 1 if there was an error and 0 otherwise*/
char handle_errors(errors error , char* name){
	switch (error)
	{
	case ALLOCATION_FAILED:
		printf("allocation failed in: %s",name);
		return 1;
	case READ_FAILED:
		printf("read failed in: %s",name);
		return 1;
	case NONE:
		return 0;
	}
}

int main(int argc, char* argv[]) {
	/*variables*/
	FILE* input;
	spmat *A , *B_g1, *B_g2;
	int nof_vertex;
	int* degree;
	int* temp;
	int M=0;
	int *g1, *g2; /*the 2 groups after the partition*/
	int sof_g1, sof_g2;
	errors error;

	/*try to open the input file*/
	input = fopen(argv[1] ,"r");
	if(!input){
		printf("input file is invalid");
		return 5;
	}
	/*try to read the number of vertexes*/
	if(fread(&nof_vertex, sizeof(int), 1, input) != 1){ /*TODO: check if file pass to a function is rewined??*/
		printf("read from input file failed");
		return 5;
	}
	/*allocate vector to save the degree of each vector (k_i)*/
	degree =(int*) malloc(nof_vertex * sizeof(int));
	if(!degree){
		printf("malloc failed on pointer degree");
		return 5;
	}
	/*try allocate sparse matrix using list imp*/
	A = spmat_allocate_list(nof_vertex);
	if(!A){
		printf("sparse matrix allocation failed on A");
		return 5;
	}

	error = read_input(input, A , degree, nof_vertex);
	if(handle_errors(error ,"read_input")){
		return 5;
	}
	/*computing M*/
	for(temp=degree; temp<degree+nof_vertex; temp++){
		M += *temp;
	}
	g1 = (int*) malloc(nof_vertex*sizeof(int));
	if(!g1){
		printf("malloc failed on g1");
		return 5;
	}
	g2 = (int*) malloc(nof_vertex*sizeof(int));
	if(!g2){
		printf("malloc failed on g2");
		return 5;
	}

	/*-------partition algorithem-------
	needs to compute g1 ,g2 , sof_g1 , sof_g2
	group g is a vector of size nof_vertex, if vertex i is in g then g[i]=1 else g[i]=0
	& allocate a sparse matrix : B_g1 and B_g2*/
	B_g1 = spmat_allocate_list(sof_g1);
	if(!B_g1){
		printf("sparse matrix allocation failed on B_g1");
		return 5;
	}
	error=compute_modularity_matrix(A, g1, degree, M, B_g1);
	if(handle_errors(error ,"compute_modularity_matrix: g1")){
		return 5;
	}

	B_g2 = spmat_allocate_list(sof_g2);
	if(!B_g2){
		printf("sparse matrix allocation failed on B_g2");
		return 5;
	}
	error=compute_modularity_matrix(A, g2, degree, M, B_g2);
	if(handle_errors(error ,"compute_modularity_matrix: g2")){
		return 5;
	}


	free(degree);
	free(g1);
	free(g2);
	fclose(input);
	A->free;
	return 0;
}