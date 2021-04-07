#include <gsl/gsl_matrix.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>

void matrix_print(gsl_matrix* A,FILE* fil) {
    for (int i = 0; i < (A->size1); i++) {
        for (int j = 0; j < (A->size2); j++) {
            double Aij = gsl_matrix_get(A, i, j);
            fprintf(fil, "%0.2f  ", Aij);
        }
        fprintf(fil, "\n");
    }
}

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
        ,double a, gsl_vector* ya, double **Yal,
        double steps, double **Xal, double b, double h,
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
        rkstep12(f,xpos,&yx.vector,h, yplaceholder,err);
        dy = gsl_blas_dnrm2(err);
        normy = gsl_blas_dnrm2(yplaceholder);
        tol = (normy*eps+acc)*sqrt(h/(b-a));
        if(dy<tol){
            k++;
            if(k>=steps){
                //printf("Her %i\n",k);
                #warning original memory to small additional is added
                *Yal = realloc(*Yal, sizeof(double)*(n*(k+1)));
                *Xal = realloc(*Xal, sizeof(double)*(k+1));
                Y = gsl_matrix_view_array(*Yal, k+1, n);
                X = gsl_vector_view_array(*Xal,k+1);

            }
            xpos+=h;
            //printf("first row at k: %i\n", k);
            //matrix_print(&Y.matrix,stdout);
            gsl_vector_set(&X.vector,k,xpos);
            gsl_matrix_set_row(&Y.matrix,k,yplaceholder);
            //gsl_vector_view fr = gsl_matrix_row(&Y.matrix,0);
            //printf("efter:\n",k);
            //matrix_print(&Y.matrix,stdout);
            //gsl_vector_fprintf(stdout,&fr.vector,"%g");

        }
        if(dy>0) h*=pow(tol/dy,0.25)*0.95; else h*=2;
    }
    return k+1;
    }


void harmosc(double x, gsl_vector* y, gsl_vector* dydx){
    double yi = gsl_vector_get(y,0);
    double dy = gsl_vector_get(y,1);
    gsl_vector_set(dydx,1,-yi);
    gsl_vector_set(dydx,0,dy);
}
void SIR(double x, gsl_vector* y, gsl_vector* dydx){
    double N = 5.5e6;
    double contact = 2. ;
    double Tr = 9.; //days
    double Tc = Tr/contact;
    double s = gsl_vector_get(y,0);
    double i = gsl_vector_get(y,1);
    double dsdt = -i*s/(N*Tc);
    assert(dsdt<0);
    double didt = i*s/(N*Tc)-i/Tr;
    double drdt = i/Tr;
    gsl_vector_set(dydx,0,dsdt);
    gsl_vector_set(dydx,1,didt);
    gsl_vector_set(dydx,2,drdt);
}

int main(){
    FILE* harmOsc = fopen("out.harmosc.txt","w");
    FILE* disease = fopen("out.disease.txt","w");
    double eps = 0.005;
    double abs = 0.005;
    double a=0, b=M_PI;
    int k;
    gsl_vector* ya = gsl_vector_calloc(2);
    gsl_vector_set(ya,1,1.);
    int n=2, m=40;
    double* Yal = malloc(sizeof(double)*(m*n));
    double* Xal = malloc(sizeof(double)*(m));
    k = bby_driver(harmosc,a,ya,&Yal,m,&Xal,b,(b-a)/30,abs,eps);
    gsl_matrix_view Y = gsl_matrix_view_array(Yal, k, n);
    gsl_vector_view X = gsl_vector_view_array(Xal, k);


    //Disease part
    int dk;
    double da=0., db=100.;
    gsl_vector* dya = gsl_vector_calloc(3);
    double I = 661., hadDisease=234318.;
    double vaccinated = 409.187, R=hadDisease+vaccinated;
    double S = 5.5e6-I-R;
    gsl_vector_set(dya,0,S);
    gsl_vector_set(dya,1,I);
    gsl_vector_set(dya,2,R);
    int dn=3, dm=6;
    double* Dal = malloc(sizeof(double)*(dm*dn));
    double* Tal = malloc(sizeof(double)*(dm));
    dk = bby_driver(SIR,da,dya,&Dal,dm,&Tal,db,(db-da)/dm,abs,eps);
    gsl_matrix_view D = gsl_matrix_view_array(Dal,dk,dn);
    gsl_vector_view T = gsl_vector_view_array(Tal, dk);
    //printf("%i",dk);
    //printf("her");
    //gsl_vector_fprintf(stdout, dya,"%g");
    //printf("Her\n");
    //matrix_print(&D.matrix,stdout);
    for(int i =0; i<dk;i++){
        double ti = gsl_vector_get(&T.vector,i);
        fprintf(disease,"%g %g %g %g\n",ti,gsl_matrix_get(&D.matrix,i,0),gsl_matrix_get(&D.matrix,i,1),gsl_matrix_get(&D.matrix,i,2));
        //printf("%g %g %g %g \n",ti,gsl_matrix_get(&D.matrix,i,0),gsl_matrix_get(&D.matrix,i,1),gsl_matrix_get(&D.matrix,i,2));

    }
    for(int i =0; i<k;i++){
        double xi = gsl_vector_get(&X.vector,i);
        fprintf(harmOsc,"%g %g %g\n",xi,gsl_matrix_get(&Y.matrix,i,0),sin(xi));
        //printf("%g %g %g\n",xi,gsl_matrix_get(&Y.matrix,0,i),sin(xi));
    }

    free(Xal);
    free(Yal);
    free(Tal);
    free(Dal);
	return 0;
}
