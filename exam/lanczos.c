//
// Created by Martin on 15-06-2021.
//

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>


void matrix_print(gsl_matrix* A,FILE* fil){
    for(int i=0; i<(A->size1);i++){
        for(int j=0; j<(A->size2); j++){
            double Aij = gsl_matrix_get(A,i,j);
            fprintf(fil, "%0.3f  ",Aij);
        }
        fprintf(fil,"\n");
    }
}

void GS_ortho_arnoldi(gsl_vector* q, gsl_matrix* Q, gsl_matrix* H, int k){
    for(int i=0; i <= k; i++){
        gsl_vector_view qi = gsl_matrix_column(Q,i);
        double h_ik;
        gsl_blas_ddot(&qi.vector,q,&h_ik);
        //printf("(%i,%i), h = %g \n",i,k,h_ik);
        gsl_matrix_set(H,i,k,h_ik);
        gsl_blas_daxpy(-h_ik,&qi.vector,q);
    }
}

void GS_ortho(gsl_vector* q, gsl_matrix* Q, gsl_matrix* H, int k) {
    if (k == 0) {
        gsl_vector_view qi = gsl_matrix_column(Q, k);
        double h_ik;
        gsl_blas_ddot(&qi.vector, q, &h_ik);
        //printf("(%i,%i), h = %g \n",i,k,h_ik);
        gsl_matrix_set(H, k, k, h_ik);
        gsl_blas_daxpy(-h_ik, &qi.vector, q);
    }
    else {
        for (int i = k - 1; i <= k; i++) {
            gsl_vector_view qi = gsl_matrix_column(Q, i);
            double h_ik;
            gsl_blas_ddot(&qi.vector, q, &h_ik);
            //printf("(%i,%i), h = %g \n",i,k,h_ik);
            gsl_matrix_set(H, i, k, h_ik);
            gsl_blas_daxpy(-h_ik, &qi.vector, q);
        }
    }
}


//q0 is destroyed, but can be found in Q
void lanczos(gsl_matrix* A, gsl_matrix* Q, gsl_matrix* H, gsl_vector* q0){

    //Checking the sizes of the matrices
    assert(A->size2 == A->size1 && A->size2 == Q->size1 && Q->size2 == H->size1);
    int m = A->size1;
    int n = Q->size2;
    gsl_matrix_set_zero(H);
    // Normalising the initial random vector given through q0
    double norm_q0 = gsl_blas_dnrm2(q0);
    gsl_vector_scale(q0,1./norm_q0); // q0 = 1/normq0*q0

    // Setting 0th coloumn of Q to Q0
    gsl_matrix_set_col(Q,0,q0);

    for(int k=0; k<n; k++){
        //q_(k+1) = A*q_k
        gsl_vector_view qk = gsl_matrix_column(Q,k);

        gsl_blas_dgemv(CblasNoTrans,1.,A,&qk.vector,0.,q0);
        //printf("k = %i:\n",k);
        //gsl_vector_fprintf(stdout,q0,"%g");

        //Orthogonalising q_0 to the others using Gram Scmidt
        GS_ortho(q0,Q,H,k);

        // Normalising the new q0 and saving the norm
        double norm= gsl_blas_dnrm2(q0);
        gsl_vector_scale(q0,1./norm);
        if(k!=(n-1)) {
            gsl_matrix_set(H, k + 1, k, norm);

            //Adding the new q0 to the set of Q's
            gsl_matrix_set_col(Q, k + 1, q0);
        }
    }
}

