#include "min_declarations.h"
#include "integration_declarations.h"
#include <gsl/gsl_integration.h>
//Using gaussian wavelet as activation fuction
double gaussWave(double x){
    return x*exp(-x*x);
}

double gaussWave_diff(double x) {
    return exp(-x * x) - 2 * x * x * exp(-x * x);
}
double gaussWave_diffdiff(double x){
    return 2*x*exp(-x*x)*(2*x*x-3);
}

double gaussWave_int(double x){
    return -0.5*exp(-x*x);
}

//Allocating the memory to the ann structure
//Note to handle all the requirements we both integral and derivative, first and second, of the
//training function
ann*   ann_alloc   (int n,double(*f)(double),double(*f_diff)(double),double(*f_diffdiff)(double),double(* f_int)(double )){
    ann* network = malloc(sizeof(ann));
    gsl_vector* params = gsl_vector_alloc(3*n);
    network->params = params;
    network->f =f;
    network->f_diff = f_diff;
    network->f_diffdiff = f_diffdiff;
    network -> f_int = f_int;
    return network;
}

void   ann_free    (ann* network){
    gsl_vector_free(network->params);
}


double ann_response(ann* network,double x){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        //Getting the parameters assuming an (a0,b0,w0),(a1,b1,w1)...
        // structure
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=((network->f)((x-ai)/bi))*wi;
    }
    return Fp;
}

//From a trained netwotk find the first derivative of the function
double ann_diff(ann* network, double x){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        //Getting the parameters assuming an (a0,b0,w0),(a1,b1,w1)...
        // structure
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=((network->f_diff)((x-ai)/bi))*wi/bi;
    }
    return Fp;
}

//From a trained netwotk find the second derivative of the function
double ann_diffdiff(ann* network, double x){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        //Getting the parameters assuming an (a0,b0,w0),(a1,b1,w1)...
        // structure
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=((network->f_diffdiff)((x-ai)/bi))*wi/(bi*bi);
    }
    return Fp;
}

//From a trained netwotk find the integral of the function from a to b
double ann_int(ann* network, double a, double b){
    int n = (network->params)->size;
    double Fp = 0;
    for(int i=0; i<n;i+=3){
        //Getting the parameters assuming an (a0,b0,w0),(a1,b1,w1)...
        // structure
        double ai = gsl_vector_get(network->params,i);
        double bi = gsl_vector_get(network->params,i+1);
        double wi = gsl_vector_get(network->params,i+2);
        Fp+=(((network->f_int)((b-ai)/bi))*wi*bi)-(((network->f_int)((a-ai)/bi))*wi*bi);
    }
    return Fp;
}


//cost function for the trained interpolator network
double ann_cost(ann* network, gsl_vector* xs, gsl_vector* ys){
    int N = xs->size;
    double cost = 0;
    for(int i=0; i<N; i++){
        double xi = gsl_vector_get(xs,i);
        double yi = gsl_vector_get(ys,i);
        //printf("response:\n");
        //printf("%g\n",ann_response(network,xi));
        cost += pow((ann_response(network,xi)-yi),2);
    }
    return cost/N;
}

//The cost function for the untrained neural network, therefore wild.
double annWild_cost(ann* network, double diffeq_pow2(double responseofX, ann* network),
                    double a, double b,double boundary_x,double boundary_y
                    ,double boundary_ydot){
    double abserr;
    int eta;
    double result = CC_integrater(diffeq_pow2,a,b,network,1e-7,1e-7,&eta, &abserr);
    double y_err = pow(ann_response(network,boundary_x)-boundary_y,2)*(b-a);
    double ydot_err = pow(ann_diff(network,boundary_x)-boundary_ydot,2)*(b-a);
    double val = result+y_err+ydot_err;
    /*
    printf("result: %g\n",result);
    printf("y_err: %g\n",y_err);
    printf("ydot_err: %g\n",ydot_err);
    */
    return val;
}

//train functions for the trained and the wild neural networks, respectively.
void   ann_train   (ann* network,gsl_vector* xs,gsl_vector* ys){
    ann_amoeba(ann_cost,network,xs,ys,1e-7);
}
void   annWild_train  (ann* network, double diffeq_pow2(double responseofX, ann* network),
                       double a, double b,double boundary_x,double boundary_y
                       ,double boundary_ydot){
    annWild_amoeba(annWild_cost,network,diffeq_pow2,a, b, boundary_x, boundary_y,boundary_ydot,1e-7);
}