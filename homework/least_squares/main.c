#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include "declarations.h"

typedef double (*f)(double x);

double f0(double x){
    return 1.;
}

double f1(double x){
    return x;
}

void lst_sqr_fit(gsl_matrix* xydy, f fun[], int func_nr, gsl_vector* ck,gsl_matrix* cov){
    // Allocating a matrix to be destroyed by QR -factorisation
    gsl_matrix* A = gsl_matrix_alloc(xydy->size1,func_nr);
    gsl_matrix* R = gsl_matrix_alloc(func_nr,func_nr);
    gsl_matrix* covcpy = gsl_matrix_alloc(cov->size1,cov->size2);
    gsl_matrix* covR = gsl_matrix_alloc(cov->size1,cov->size2);
    gsl_vector_view y = gsl_matrix_column(xydy,1);
    for(int i =0 ; i<A->size1;i++){
        double xi = gsl_matrix_get(xydy,i,0);
        for(int j=0; j<A->size2;j++){
            double Aij= fun[j](xi);
            gsl_matrix_set(A,i,j,Aij);
        }
    }
    //Decomposing the A matrix
    gs_decomp(A,R);

    //Solve Rx=Q^T b
    gs_solve(A,R,(&y.vector),ck);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,R,R,0,covcpy);
    gs_decomp(covcpy,covR);
    gs_invers(covcpy,covR,cov);

    //freeing allocated space
    gsl_matrix_free(A);
    gsl_matrix_free(R);
    gsl_matrix_free(covcpy);
    gsl_matrix_free(covR);

}
int main(){
    FILE* cks = fopen("out.ckandHalflife.txt","w");
    FILE* fit = fopen("out.fitfunc.txt","w");
    FILE* dat = fopen("out.fitdat.txt","w");

    int n = 9;
    int func_n = 2;
    double x[] = {1.,2.,3.,4.,6.,9.,10.,13.,15};
    double y[] = {117.,100., 88., 72., 53.,
                  29.5, 25.2,15.2,11.1};
    f func[]={&f0, &f1};
    gsl_matrix* xydy = gsl_matrix_alloc(n,3);
    gsl_matrix* cov = gsl_matrix_alloc(func_n,func_n);
    gsl_vector* ck = gsl_vector_alloc(func_n);
    for(int i=0; i<n; i++){
        gsl_matrix_set(xydy,i,0,x[i]);
        gsl_matrix_set(xydy,i,1,log(y[i]));
        gsl_matrix_set(xydy,i,2,(y[i]/20)/y[i]);
    }
    lst_sqr_fit(xydy,func,func_n,ck,cov);
    double c0=gsl_vector_get(ck,0);
    double c1=gsl_vector_get(ck,1);
    double sigma_c0 = sqrt(gsl_matrix_get(cov,0,0));
    double sigma_c1 = sqrt(gsl_matrix_get(cov,1,1));
    fprintf(cks,"Covariance Matix:\n");
    matrix_print(cov,cks);
    fprintf(cks,"Params and Halflife:\n");
    fprintf(cks,"c0 = %g\nc1 (lambda) = %g\ncalc T_(1/2) = (%g +/- %g) days\nreal T_(1/2) = %g days\n", c0,c1, log(2)/(-c1)
    ,sigma_c1, 3.63);

    for(int i=0; i<xydy->size1;i++){
        double xi, yi, dyi;
        xi=gsl_matrix_get(xydy,i,0);
        yi=gsl_matrix_get(xydy,i,1);
        dyi=gsl_matrix_get(xydy,i,2);
        fprintf(dat,"%g %g %g\n",xi,yi,dyi);
    }
    int Steps = 100;
    for(int i=0; i<Steps;i++){
        double xi = ((double) i)/ Steps*x[n-1];
        double fi = c0*f0(xi)+c1*f1(xi);
        double fip = (c0+sigma_c0)*f0(xi)+(c1+sigma_c1)*f1(xi);
        double fim = (c0-sigma_c0)*f0(xi)+(c1-sigma_c1)*f1(xi);
        fprintf(fit,"%g %g %g %g\n",xi,fi, fip, fim);
    }

	return 0;
}
