#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0; i<v->size; i++){
		printf("%10g",gsl_vector_get(v,i));
		printf("\n");
	}
}

void matrix_print(char s[], gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0; i<A->size1; i++){
		for(int j=0; j<A->size2; j++){
		printf("%10g",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

int main(){
	int n = 3;
	double a[] = {6.13, -2.90, 5.86, 8.08, -6.31, -3.89, -4.39, 1.00, 0.19};
	double bval[] = {6.23, 5.37, 2.29};
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_calloc(n);
	for(int i = 0; i<(*A).size1; i++){
		for(int j = 0; j<(*A).size2; j++){
			gsl_matrix_set(A,i,j,a[3*i+j]);
		}
	}
	for(int i=0; i< b->size; i++){
		gsl_vector_set(b,i,bval[i]);
	}
	matrix_print("Matrix A:", A);
	gsl_matrix_memcpy(Acopy,A);
	gsl_linalg_HH_solve(Acopy,b,x);
	vector_print("vector b:", b);
	vector_print("vector x:", x);
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
	vector_print("vector Ax:", y);
return 0;
}
