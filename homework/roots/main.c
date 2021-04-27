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

static double E;
static double R_max;
static double R_max_bound;
static double E_bound;

void N_rootfinder(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){
    int n = x->size;
    gsl_matrix* J = gsl_matrix_alloc(n,n);
    gsl_matrix* R = gsl_matrix_alloc(n,n);
    gsl_vector* Dx = gsl_vector_alloc(n);
    gsl_vector* xLamb = gsl_vector_alloc(n);
    gsl_vector* fx = gsl_vector_alloc(n);
    gsl_vector* Dfx = gsl_vector_alloc(n);
    f(x,fx);
    int N_MAX =0;
    while(gsl_blas_dnrm2(fx)>eps && N_MAX<10000) {
        N_MAX++;
        gsl_vector_memcpy(Dx, x);
        for (int j = 0; j < n; j++) {
            double xj = gsl_vector_get(x, j);
            gsl_vector_set(Dx, j, xj + sqrt(DBL_EPSILON));
            f(Dx, Dfx);
            for (int i = 0; i < n; i++) {
                double J_ij = (gsl_vector_get(Dfx, i) - gsl_vector_get(fx, i))/sqrt(DBL_EPSILON);
                gsl_matrix_set(J, i, j, -J_ij);
            }
        }
        double lambda = 1.;
        gs_decomp(J,R);
        gs_solve(J, R, fx, Dx);
        //gsl_linalg_HH_solve(J, fx, Dx);
        for (int i = 0; i < n; i++) {
            double xi = gsl_vector_get(x, i);
            double Dxi = gsl_vector_get(Dx, i);
            gsl_vector_set(xLamb, i, xi + lambda * Dxi);
        }
        f(xLamb, Dfx);
        while (gsl_blas_dnrm2(Dfx) > (1 - lambda / 2) * gsl_blas_dnrm2(fx) && lambda > 1. / 64) {
            lambda /= 2;
            for (int i = 0; i < n; i++) {
                double xi = gsl_vector_get(x, i);
                double Dxi = gsl_vector_get(Dx, i);
                gsl_vector_set(xLamb, i, xi + lambda * Dxi);
            }
            f(xLamb, Dfx);
        }
        gsl_vector_memcpy(x, xLamb);
        f(x, fx);
    }

    gsl_vector_free(Dx);
    gsl_vector_free(xLamb);
    gsl_vector_free(fx);
    gsl_vector_free(Dfx);
    gsl_matrix_free(J);
}

void xAnden(gsl_vector* x, gsl_vector* fx){
    gsl_vector_set(fx,0,pow(gsl_vector_get(x,0),3));
    gsl_vector_set(fx,1,pow(gsl_vector_get(x,1),3));
}

void Rosen_diff(gsl_vector* x, gsl_vector* fx){
    double X = gsl_vector_get(x,0);
    double Y = gsl_vector_get(x,1);
    gsl_vector_set(fx,0,2*(1-X)-4*X*100*(Y-X*X));
    gsl_vector_set(fx,1,2*100*(Y-X*X));
}


void sch_eq(double x, gsl_vector* y, gsl_vector* dydx){
    gsl_vector_set(dydx,0,gsl_vector_get(y,1));
    double fx = gsl_vector_get(y,0);
    double ddy = -2.*(E+1./x)*fx;
    gsl_vector_set(dydx,1,ddy);
}
void ground_state(gsl_vector* eps, gsl_vector* fx){
    double a = 0.01;
    double acc = 0.001;
    double rel_acc = 0.001;
    E=gsl_vector_get(eps,0);
    gsl_vector* ya = gsl_vector_alloc(2);
    gsl_vector_set(ya,0,(a-a*a));
    gsl_vector_set(ya,1,(1-2.*a));
    bby_driver(sch_eq,a,ya,R_max,0.1,acc,rel_acc);
    gsl_vector_set(fx,0,gsl_vector_get(ya,0));
}
void ground_state_bound(gsl_vector* eps, gsl_vector* fx){
    double a = 0.01;
    double acc = 0.001;
    double rel_acc = 0.001;
    E_bound=gsl_vector_get(eps,0);
    gsl_vector* ya = gsl_vector_alloc(2);
    gsl_vector_set(ya,0,(a-a*a));
    gsl_vector_set(ya,1,(1-2.*a));
    bby_driver(sch_eq,a,ya,R_max_bound,0.1,acc,rel_acc);
    long double Ffoer = gsl_vector_get(ya,0);
    long double Fefter =(R_max_bound*exp(-sqrt(-2*E_bound)*R_max_bound));
    //printf("%g %g %g\n",Ffoer,Fefter,E_bound);
    gsl_vector_set(fx,0,Ffoer-Fefter);
}

int main(){
    FILE* psi0 = fopen("out.psi0.txt","w");
    FILE* conv = fopen("out.conv.txt","w");
    FILE* bound = fopen("out.convBound.txt","w");
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector* y = gsl_vector_alloc(1);
    gsl_vector* y_bound = gsl_vector_alloc(1);
    gsl_vector_set(x,0,1.3);
    gsl_vector_set(x,1,1.3);
    gsl_vector_set(y,0,-3);
    gsl_vector_set(y_bound,0,-1);
    N_rootfinder(Rosen_diff,x,1e-4);
    printf("The extremum of Rosen is at (x= %g,y=%g) \n", gsl_vector_get(x,0),gsl_vector_get(x,1));

    R_max = 8.;
    R_max_bound = 0.5;
    N_rootfinder(ground_state,y,0.001);
    N_rootfinder(ground_state_bound,y_bound,0.001);
    printf("Lowest energy: E0 = %g %g\n",E, E_bound);
    double eps = 0.001;
    double abs = 0.001;
    double a=0.01, b=8.;
    gsl_vector* ya = gsl_vector_calloc(2);
    gsl_vector_set(ya,0,(a-a*a));
    gsl_vector_set(ya,1,(1-2.*a));
    int n=2, m=200;
    double* Yal = malloc(sizeof(double)*((m*n)));
    double* Xal = malloc(sizeof(double)*m);
    adult_driver(sch_eq,a,ya,&Yal,m,&Xal,b,(b-a)/10,abs,eps);
    gsl_matrix_view Y = gsl_matrix_view_array(Yal, m, n);
    gsl_vector_view X = gsl_vector_view_array(Xal, m);
    for(int i = 0; i<m; i++){
        double xi = gsl_vector_get(&X.vector,i);
        if(xi!=0){
            double th_exp = xi*exp(-xi);
            double calc = gsl_matrix_get(&Y.matrix,i,0);
            fprintf(psi0,"%g %g %g\n",xi, calc, th_exp);
        }
    }

    gsl_vector_set(y_bound,0,-1);
    gsl_vector_set(y,0,-3);
    for(int i = 1; i<50; i++){
        R_max = (8.)*(double)i/50.;
        R_max_bound = (2.3)*(double)i/50.;
        N_rootfinder(ground_state,y,0.001);
        N_rootfinder(ground_state_bound,y_bound,0.001);
        fprintf(conv,"%g %g \n", R_max,E+0.5);
        fprintf(bound,"%g %g \n", R_max_bound,E_bound+0.5);
        gsl_vector_set(y_bound,0,-1);
        gsl_vector_set(y,0,-3);
    }

    fclose(psi0);
    fclose(conv);
    fclose(bound);
    free(Yal);
    free(Xal);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(y_bound);
    gsl_vector_free(ya);
    return 0;
}
