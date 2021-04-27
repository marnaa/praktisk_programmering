#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <gsl/gsl_vector.h>



void matrix_print(gsl_matrix* A,FILE* fil){
    for(int i=0; i<(A->size1);i++){
        for(int j=0; j<(A->size2); j++){
            double Aij = gsl_matrix_get(A,i,j);
            fprintf(fil, "%0.3g ",Aij);
        }
        fprintf(fil,"\n");
    }
}
//Gram-Schmidt orthogonalisation function
void gs_decomp(gsl_matrix* A, gsl_matrix* R) {
    gsl_vector *ai = gsl_vector_alloc((A->size1));
    gsl_vector *qi = gsl_vector_alloc((A->size1));
    //Putting normalised ai's into Ri
    for (int i = 0; i < (A->size2); i++) {
        //finding norm of a_i
        gsl_matrix_get_col(ai, A, i);
        gsl_vector_memcpy(qi,ai);
        double ai_norm = gsl_blas_dnrm2(ai);
        assert(ai_norm>0);
        //Setting R_ii = norm(a_i, a_i)
        gsl_matrix_set(R, i, i, ai_norm);
        gsl_vector_scale(qi,1/ai_norm);
        gsl_matrix_set_col(A,i,qi);
        for(int j=i+1; j<(A->size2); j++){
            //a_i==a_j
            gsl_matrix_get_col(ai,A,j);
            double Rij=0;
            //Rij=q^T a_i
            gsl_blas_ddot(qi,ai,&Rij);
            //a_j=-q_i+a_i
            gsl_blas_daxpy(-Rij,qi,ai);
            //set a_j coloumn
            gsl_matrix_set(R,i,j,Rij);
            gsl_matrix_set(R,j,i,0.);
            gsl_matrix_set_col(A,j,ai);
        }

    }
    gsl_vector_free(ai);
    gsl_vector_free(qi);
}

void backsub(gsl_matrix* R, gsl_vector* x){
    for(int i=(x->size)-1; i>=0;i--){
        double QTbi = gsl_vector_get(x,i);
        double Rii = gsl_matrix_get(R,i,i);
        double xi = QTbi;
        for(int j=i+1; j<x->size;j++){
            double Rij = gsl_matrix_get(R,i,j);
            double xj = gsl_vector_get(x,j);
            xi-=Rij*xj;
        }
        xi /= Rii;
        gsl_vector_set(x,i,xi);
    }
}

void gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
    gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
    backsub(R,x);
}

void gs_invers(gsl_matrix *Q, gsl_matrix* R, gsl_matrix* X){
    //Need to solve QRB=1, in principle also BQR=1, but not nice form
    gsl_matrix_transpose_memcpy(X,Q);
    for(int i=0; i<(Q->size2); i++){
        gsl_vector_view qi = gsl_matrix_column(X,i);
        backsub(R, &qi.vector);
    }
}