#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
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
                return -k ;
                printf("Her %i\n",k);
                Yal = realloc(*Yal, sizeof(double)*(n*(k+1)));
                Xal = realloc(*Xal, sizeof(double)*(k+1));
                Y = gsl_matrix_view_array(*Yal, k+1, n);
                X = gsl_vector_view_array(*Xal,k+1);
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

void SIR3(double x, gsl_vector* y, gsl_vector* dydx){
    double N = 5.5e6;
    double contact = 3. ;
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
void threebody( double x, gsl_vector* y, gsl_vector* dydx){
    double G = 1.;
    double m1 = 1.;
    double m2 = 1.;
    double m3 = 1.;
    double r1x = gsl_vector_get(y,0);
    double r1y = gsl_vector_get(y,1);
    double r2x = gsl_vector_get(y,2);
    double r2y = gsl_vector_get(y,3);
    double r3x = gsl_vector_get(y,4);
    double r3y = gsl_vector_get(y,5);
    double dr1x = gsl_vector_get(y,6);
    double dr1y = gsl_vector_get(y,7);
    double dr2x = gsl_vector_get(y,8);
    double dr2y = gsl_vector_get(y,9);
    double dr3x = gsl_vector_get(y,10);
    double dr3y = gsl_vector_get(y,11);
    double r1r2 = sqrt(pow((r1x-r2x),2)+pow((r1y-r2y),2));
    double r1r3 = sqrt(pow((r1x-r3x),2)+pow((r1y-r3y),2));
    double r2r3 = sqrt(pow((r2x-r3x),2)+pow((r2y-r3y),2));
    assert(r1r2>0);
    assert(r2r3>0);
    assert(r1r3>0);

    double ddr1x = -G*m2/pow(r1r2,3)*(r1x-r2x)-G*m3/pow(r1r3,3)*(r1x-r3x);
    double ddr1y = -G*m2/pow(r1r2,3)*(r1y-r2y)-G*m3/pow(r1r3,3)*(r1y-r3y);
    double ddr2x = -G*m3/pow(r2r3,3)*(r2x-r3x)-G*m1/pow(r1r2,3)*(r2x-r1x);
    double ddr2y = -G*m3/pow(r2r3,3)*(r2y-r3y)-G*m1/pow(r1r2,3)*(r2y-r1y);
    double ddr3x = -G*m1/pow(r1r3,3)*(r3x-r1x)-G*m2/pow(r2r3,3)*(r3x-r2x);
    double ddr3y = -G*m1/pow(r1r3,3)*(r3y-r1y)-G*m2/pow(r2r3,3)*(r3y-r2y);
    gsl_vector_set(dydx,0,dr1x);
    gsl_vector_set(dydx,1,dr1y);
    gsl_vector_set(dydx,2,dr2x);
    gsl_vector_set(dydx,3,dr2y);
    gsl_vector_set(dydx,4,dr3x);
    gsl_vector_set(dydx,5,dr3y);
    gsl_vector_set(dydx,6,ddr1x);
    gsl_vector_set(dydx,7,ddr1y);
    gsl_vector_set(dydx,8,ddr2x);
    gsl_vector_set(dydx,9,ddr2y);
    gsl_vector_set(dydx,10,ddr3x);
    gsl_vector_set(dydx,11,ddr3y);

}

int main(){
    FILE* harmOsc = fopen("out.harmosc.txt","w");
    FILE* disease = fopen("out.disease.txt","w");
    FILE* disease3 = fopen("out.disease3.txt","w");
    FILE* body3 = fopen("out.3body.txt","w");
    double eps = 0.005;
    double abs = 0.005;
    double a=0, b=M_PI;
    int k;
    gsl_vector* ya = gsl_vector_calloc(2);
    gsl_vector_set(ya,1,1.);
    int n=2, m=100;
    double* Yal = malloc(sizeof(double)*(m*n));
    double* Xal = malloc(sizeof(double)*(m));
    k = bby_driver(harmosc,a,ya,&Yal,m,&Xal,b,(b-a)/30,abs,eps);
    gsl_matrix_view Y = gsl_matrix_view_array(Yal, k, n);
    gsl_vector_view X = gsl_vector_view_array(Xal, k);


    //Disease part
    int dk, dk3;
    double da=0., db=100.;
    gsl_vector* dya = gsl_vector_calloc(3);
    double I = 661., hadDisease=234318.;
    double vaccinated = 409.187, R=hadDisease+vaccinated;
    double S = 5.5e6-I-R;
    gsl_vector_set(dya,0,S);
    gsl_vector_set(dya,1,I);
    gsl_vector_set(dya,2,R);
    int dn=3, dm=300;
    double* Dal = malloc(sizeof(double)*(dm*dn));
    double* Tal = malloc(sizeof(double)*(dm));
    double* Dal3 = malloc(sizeof(double)*(dm*dn));
    double* Tal3 = malloc(sizeof(double)*(dm));
    dk = bby_driver(SIR,da,dya,&Dal,dm,&Tal,db,(db-da)/dm,abs,eps);
    dk3 = bby_driver(SIR3,da,dya,&Dal3,dm,&Tal3,db,(db-da)/dm,abs,eps);
    gsl_matrix_view D = gsl_matrix_view_array(Dal,dk,dn);
    gsl_vector_view T = gsl_vector_view_array(Tal, dk);
    gsl_matrix_view D3 = gsl_matrix_view_array(Dal3,dk3,dn);
    gsl_vector_view T3 = gsl_vector_view_array(Tal3, dk3);
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
    for(int i =0; i<dk3;i++){
        double ti3 = gsl_vector_get(&T3.vector,i);
        fprintf(disease3,"%g %g %g %g\n",ti3,gsl_matrix_get(&D3.matrix,i,0),gsl_matrix_get(&D3.matrix,i,1),gsl_matrix_get(&D3.matrix,i,2));
        //printf("%g %g %g %g \n",ti,gsl_matrix_get(&D.matrix,i,0),gsl_matrix_get(&D.matrix,i,1),gsl_matrix_get(&D.matrix,i,2));
    }

    for(int i =0; i<k;i++){
        double xi = gsl_vector_get(&X.vector,i);
        fprintf(harmOsc,"%g %g %g\n",xi,gsl_matrix_get(&Y.matrix,i,0),sin(xi));
        //printf("%g %g %g\n",xi,gsl_matrix_get(&Y.matrix,0,i),sin(xi));
    }


    //3-body system

    double c=0., d=9.;
    gsl_vector* vec0 = gsl_vector_calloc(12);
    gsl_vector* vecplace = gsl_vector_calloc(12);
    gsl_vector* err = gsl_vector_calloc(12);
    double r1x,r1y,r2x,r2y,r3x,r3y;
    double dr1x,dr1y,dr2x,dr2y,dr3x,dr3y;
    r1x=0.97000436;
    r1y=-0.24308753;
    r2x=-0.97000436;
    r2y=0.24308753;
    r3x=0.;
    r3y=0.;
    dr3x=-0.93240737;
    dr3y=-0.86473146;
    dr1x=-dr3x/2;
    dr1y=-dr3y/2;
    dr2x=-dr3x/2;
    dr2y=-dr3y/2;

    gsl_vector_set(vec0,0,r1x);
    gsl_vector_set(vec0,1,r1y);
    gsl_vector_set(vec0,2,r2x);
    gsl_vector_set(vec0,3,r2y);
    gsl_vector_set(vec0,4,r3x);
    gsl_vector_set(vec0,5,r3y);
    gsl_vector_set(vec0,6,dr1x);
    gsl_vector_set(vec0,7,dr1y);
    gsl_vector_set(vec0,8,dr2x);
    gsl_vector_set(vec0,9,dr2y);
    gsl_vector_set(vec0,10,dr3x);
    gsl_vector_set(vec0,11,dr3y);
    double N=12;
    double M=200;
    double ABS =0.01;
    double EPS = 0.01;
    double* Ral = malloc(sizeof(double)*(M*N));
    double* Zal = malloc(sizeof(double)*(M));
    int j = bby_driver(threebody,c,vec0,&Ral,M,&Zal,d,0.001 ,ABS,EPS);
    gsl_matrix_view r = gsl_matrix_view_array(Ral, M, N);
    gsl_vector_view Z = gsl_vector_view_array(Zal, M);
    for(int i =0; i<M;i++){
        double zi = gsl_vector_get(&Z.vector,i);
        fprintf(body3,"%g %g %g %g %g %g %g %g %g %g %g %g %g\n",zi,gsl_matrix_get(&r.matrix,i,0),gsl_matrix_get(&r.matrix,i,1),
                gsl_matrix_get(&r.matrix,i,2),gsl_matrix_get(&r.matrix,i,3),
                gsl_matrix_get(&r.matrix,i,4),gsl_matrix_get(&r.matrix,i,5),
                gsl_matrix_get(&r.matrix,i,6),gsl_matrix_get(&r.matrix,i,7),
                gsl_matrix_get(&r.matrix,i,8),gsl_matrix_get(&r.matrix,i,9),
                gsl_matrix_get(&r.matrix,i,10),gsl_matrix_get(&r.matrix,i,11));
    }

    gsl_vector_free(vec0);
    free(Xal);
    free(Yal);
    free(Tal);
    free(Dal);
    free(Tal3);
    free(Dal3);
    gsl_vector_free(ya);
    gsl_vector_free(dya);
    gsl_vector_free(vecplace);
    gsl_vector_free(err);
    fclose(disease);
    fclose(disease3);
    fclose(harmOsc);
    fclose(body3);
	return 0;
}
