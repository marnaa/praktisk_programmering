#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>
#include "declarations.h"
#include <gsl/gsl_blas.h>
#include <float.h>


void rkstep12(void f(double x, gsl_vector* y, gsl_vector* dydx)
        , double x, gsl_vector* yx, double h,
              gsl_vector* yh, gsl_vector* err) {
    int n = yx->size;
    gsl_vector *k0 = gsl_vector_alloc(n);
    gsl_vector *k1 = gsl_vector_alloc(n);
    gsl_vector *yt = gsl_vector_alloc(n);
    //Calculate first order
    f(x, yx, k0);
    for (int i = 0; i < n; i++) {
        double yxi = gsl_vector_get(yx, i);
        double k0i = gsl_vector_get(k0, i);
        double yti = yxi + h / 2 * k0i;
        gsl_vector_set(yt, i, yti);
    }
    //Calculate second order
    f(x + h / 2, yt, k1);
    for (int i = 0; i < n; i++) {
        double yxi = gsl_vector_get(yx, i);
        double k1i = gsl_vector_get(k1, i);
        double yhi = yxi + h * k1i;
        gsl_vector_set(yh, i, yhi);
    }
    //Error estimate
    for (int i = 0; i < n; i++) {
        double k0i = gsl_vector_get(k0, i);
        double k1i = gsl_vector_get(k1, i);
        double erri = (k0i - k1i) * h / 2;
        gsl_vector_set(err,i,erri);
    }

    gsl_vector_free(k0);
    gsl_vector_free(k1);
    gsl_vector_free(yt);
}

int bby_driver(void f(double x, gsl_vector* y, gsl_vector* dydx)
        ,double a, gsl_vector* ya, double b, double h,
               double acc, double eps){
    /*Y is space for the solution it should be of the size n*m where
     m is assumed number of steps and n is the dimensionality of the problem*/
    /* x is an empty vector with a many entries ass there are columns
     in Y */
    int k = 0;
    double dy,normy,tol;
    int n = ya-> size;
    gsl_vector* err = gsl_vector_alloc(n);
    gsl_vector* yplaceholder = gsl_vector_alloc(n);
    double xpos = a;
    while (xpos < b ){
        if (xpos+h>b) h=b-xpos;
        rkstep12(f,xpos,ya,h, yplaceholder,err);
        dy = gsl_blas_dnrm2(err);
        normy = gsl_blas_dnrm2(yplaceholder);
        tol = (normy*eps+acc)*sqrt(h/(b-a));
        if(dy<tol){
            xpos+=h;
            gsl_vector_memcpy(ya,yplaceholder);
        }
        if(dy>0) h*=pow(tol/dy,0.25)*0.95; else h*=2;
    }
    return k;
}

// Same as bby-driver but now with the added feature of being able to store the differentials of the ODE and not
//just the result
int adult_driver(void f(double x, gsl_vector* y, gsl_vector* dydx)
        ,double a, gsl_vector* ya, double **Yal,
               int steps, double **Xal, double b, double h,
               double acc, double eps){
    /*Y is space for the solution it should be of the size n*m where
     m is assumed number of steps and n is the dimensionality of the problem*/
    /* x is an empty vector with a many entries ass there are columns
     in Y */
    int k = 0;
    double dy,normy,tol;
    int n = ya-> size;
    gsl_vector* err = gsl_vector_alloc(n);
    gsl_matrix_view Y = gsl_matrix_view_array(*Yal, steps, n);
    gsl_vector_view X = gsl_vector_view_array(*Xal,steps);
    gsl_matrix_set_row(&Y.matrix,0,ya);
    gsl_vector* yplaceholder = gsl_vector_alloc(n);
    double xpos = a;
    while (xpos < b ){
        if (xpos+h>b) h=b-xpos;
        gsl_vector_view yx = gsl_matrix_row(&Y.matrix,k);
        //printf("%i %i %i\n",(&Y.matrix)->size1,(&Y.matrix)->size2,(&yx.vector)->size);
        rkstep12(f,xpos,&yx.vector,h, yplaceholder,err);
        //printf("ny\n");
        //gsl_vector_fprintf(stdout,yplaceholder,"%g");
        dy = gsl_blas_dnrm2(err);
        normy = gsl_blas_dnrm2(yplaceholder);
        tol = (normy*eps+acc)*sqrt(h/(b-a));
        if(dy<tol){
            k++;
            if(k>=steps){
                printf("%i",k);
                return -k; //Too few steps;
                printf("Her %i\n",k);
                Yal = realloc(*Yal, sizeof(double)*(n*(k+1)));
                Xal = realloc(*Xal, sizeof(double)*(k+1));
                Y = gsl_matrix_view_array(*Yal, k+1, n);
                printf("Her %i\n",k);
                X = gsl_vector_view_array(*Xal,k+1);
                printf("Her %i\n",k);
                printf("%g",gsl_vector_get(&X.vector,0));
            }
            xpos+=h;
            //printf("xpos: %g h: %g\n",xpos, h);
            //printf("first row at k: %i\n", k);
            //printf("%li %i\n",(&X.vector)->size,k);
            //gsl_vector_fprintf(stdout,&X.vector,"%g");
            gsl_vector_set(&X.vector,k,xpos);
            gsl_matrix_set_row(&Y.matrix,k,yplaceholder);
            //gsl_vector_view fr = gsl_matrix_row(&Y.matrix,0);
            //matrix_print(&Y.matrix,stdout);
            //gsl_vector_fprintf(stdout,&fr.vector,"%g");

        }
        if(dy>0) h*=pow(tol/dy,0.25)*0.95; else h*=2;
    }
    return k+1;
}
