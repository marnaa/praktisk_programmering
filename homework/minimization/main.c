#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <float.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

static int N_MAX;

static double E[30]={101,103,105,107,109,111,113,115,117,
                     119,121,123,125,127,129,131,133,135,137,139,
                     141,143,145,147,149,151,153,155,157,159};

static double sig[30] = {-0.25,-0.30,-0.15,-1.71,0.81,
                  0.65,-0.91,0.91,0.96,-2.52,
                  -1.01,2.01,4.83,4.58,1.26,
                  1.01,-1.26,0.45,0.15,-0.91,
                  -0.81,-1.41,1.36,0.50,-0.45,
                  1.61,-2.21,-1.86, 1.76,-0.50};

static double dsig[30]={2.0,2.0,1.9,1.9,1.9,1.9,1.9,1.9,1.6,
                 1.6,1.6,1.6,1.6,1.6,1.3,1.3,1.3,1.3,
                 1.3,1.3,1.1,1.1,1.1,1.1,1.1,1.1,1.1,0.9,0.9,0.9};

void gradient(double f(gsl_vector* x),gsl_vector* x, gsl_vector* gradient,
              gsl_vector* Dx){
    int n = x -> size;
    gsl_vector_memcpy(Dx,x);
    for (int j = 0; j < n; j++) {
        double xj = gsl_vector_get(x, j);
        gsl_vector_set(Dx, j, xj+sqrt(DBL_EPSILON));
        double J_j = (f(Dx)-f(x))/sqrt(DBL_EPSILON);
        gsl_vector_set(gradient, j, J_j);
        gsl_vector_memcpy(Dx,x);
    }
}
void qnewton( double f(gsl_vector* x),
              gsl_vector* x, double eps){
    int n = x->size;
    double alpha = 1e-4;
    double prikprod;
    double prikprodsy;
    gsl_matrix* B = gsl_matrix_alloc(n,n);
    gsl_vector* grad = gsl_vector_alloc(n);
    gsl_vector* xs = gsl_vector_alloc(n);
    gsl_vector* grad_vec = gsl_vector_alloc(n);
    gsl_vector* s = gsl_vector_alloc(n);
    gsl_vector* grad_xs=gsl_vector_alloc(n);
    gsl_vector* y=gsl_vector_alloc(n);
    gsl_vector* u=gsl_vector_alloc(n);
    N_MAX = 0;
    gsl_matrix_set_identity(B);
    gradient(f,x,grad,grad_vec);
    while(gsl_blas_dnrm2(grad)>eps && N_MAX<10000){
        N_MAX++;
        double lambda = 1;
        gsl_blas_dgemv(CblasNoTrans,-1.,B,grad,0.,s);
        for(int i=0; i<n; i++){
            gsl_vector_set(xs,i, gsl_vector_get(s,i)+gsl_vector_get(x,i));
        }
        gsl_vector_memcpy(xs,x);
        gsl_vector_add(xs,s);
        gsl_blas_ddot(s,grad,&prikprod);
        double fx = f(x);
        while(f(xs)>fx+alpha*prikprod){
            if(lambda<sqrt(DBL_EPSILON)){
                gsl_matrix_set_identity(B);
                break;
            }
            lambda/=2;
            //printf("%g %g\n",f(xs), f(x)+alpha*prikprod);
            gsl_vector_scale(s,0.5);
            gsl_vector_memcpy(xs,x);
            gsl_vector_add(xs,s);
            gsl_blas_ddot(s,grad,&prikprod);
        }
        gradient(f,xs,grad_xs,grad_vec);
        gsl_vector_memcpy(y,grad_xs);
        gsl_blas_daxpy(-1,grad,y); //y = grad(x+s)-grad(x)
        gsl_vector_memcpy(u,s);
        gsl_blas_dgemv(CblasNoTrans,-1.,B,y,1.,u);
        gsl_blas_ddot(s,y,&prikprodsy);
        if(fabs(prikprodsy)>1e-12){
            gsl_blas_dger(1./prikprodsy,u,s,B);
        }
        gsl_vector_memcpy(x,xs);
        gsl_vector_memcpy(grad ,grad_xs);
    }
    gsl_vector_free(grad);
    gsl_matrix_free(B);
    gsl_vector_free(grad_vec);
    gsl_vector_free(xs);
    gsl_vector_free(y);
    gsl_vector_free(grad_xs);
    gsl_vector_free(u);
    gsl_vector_free(s);
}



double xiAnden(gsl_vector* xs){
    double x = gsl_vector_get(xs,0);
    double y = gsl_vector_get(xs,1);
    double z = gsl_vector_get(xs,2);
    return sqrt(x*x+y*y+z*z);
}
double rosenbrock(gsl_vector* x){
    double xs = gsl_vector_get(x,0);
    double ys = gsl_vector_get(x,1);
    return pow((1-xs),2)+100.*pow((ys-xs*xs),2);
}
double himmelblau(gsl_vector* xs){
    double x = gsl_vector_get(xs,0);
    double y = gsl_vector_get(xs,1);
    return pow((x*x+y-11),2)+pow((x+y*y-7),2);
}
double breitWigner(gsl_vector* xs,double E){
    double m = gsl_vector_get(xs,0);
    double Gamma = gsl_vector_get(xs,1);
    double A = gsl_vector_get(xs,2);
    return A/(pow(E-m,2)+Gamma*Gamma/4);
}
double deviation_BW(gsl_vector* xs){
    double val=0;
    for(int i =0; i<30; i++){
        val+= pow(breitWigner(xs,E[i])-sig[i],2)/(dsig[i]*dsig[i]);
    }
    return val;
}
int main(){
    FILE* exB = fopen("out.exB.txt","w");
    //vector for norm function
    gsl_vector* x= gsl_vector_alloc(3);
    gsl_vector_set(x,0,3.);
    gsl_vector_set(x,1,3.);
    gsl_vector_set(x,2,3.);

    qnewton(xiAnden,x,0.01);
    printf("min euler norm: (%g, %g, %g)\n", gsl_vector_get(x,0), gsl_vector_get(x,1), gsl_vector_get(x,2));

    //vector for himmelblau and rosenbrock
    gsl_vector* y= gsl_vector_alloc(2);
    gsl_vector_set(y,0,0);
    gsl_vector_set(y,1,0);

    qnewton(rosenbrock,y,0.001);
    printf("rosenbrock min: (%g, %g) used steps = %i\n", gsl_vector_get(y,0), gsl_vector_get(y,1),N_MAX);
    gsl_vector_set(y,0,0);
    gsl_vector_set(y,1,0);
    qnewton(himmelblau,y,0.001);
    printf("himmelblau min: (%g, %g) used steps = %i\n", gsl_vector_get(y,0), gsl_vector_get(y,1),N_MAX);

    gsl_vector_set(x,0,100.);
    gsl_vector_set(x,1,10.);
    gsl_vector_set(x,2,10.);
    qnewton(deviation_BW,x,0.001);
    for(int i =0; i<30;i++){
        fprintf(exB,"%g %g %g %g\n",E[i],sig[i],dsig[i], breitWigner(x,E[i]));
    }

    gsl_vector_free(x);
    gsl_vector_free(y);
    fclose(exB);
	return 0;
}
