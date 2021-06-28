#include "ann_declarations.h"
#include <gsl/gsl_specfunc.h>

//Defining the training functions
double gauss(double x){
    return exp(-x*x);
}
double gauss_int(double a, double b){
    return 0.5*sqrt(M_PI)*(erf(b)-erf(a));
}
double gauss_diff(double x) {
    return -2 * x * exp(-x * x);
}

double diffeq_pow2(double x, ann* network){
    double y = ann_response(network,x);
    double yd = ann_diff(network,x);
    double ydd = ann_diffdiff(network,x);
    return pow(ydd+y,2);
}

int main(){
    FILE* compare = fopen("out.compare.txt","w");
    FILE* diffeq = fopen("out.diffeq.txt","w");

    gsl_vector* xs = gsl_vector_alloc(20);
    gsl_vector* ys = gsl_vector_alloc(20);
    for(int i =-10; i<10;i++){
        double xi = (double) i/2;
        double yi = gauss(xi);
        gsl_vector_set(xs,i+10,xi);
        gsl_vector_set(ys,i+10,yi);
    }
    ann* interp = ann_alloc(3,gaussWave, gaussWave_diff,gaussWave_diffdiff, gaussWave_int);
    ann* wild = ann_alloc(4,gaussWave, gaussWave_diff,gaussWave_diffdiff, gaussWave_int);


    //PART A AND B
    double a_interp = -5;
    double b_interp = 5;
    int n_neuron_interp = ((interp->params)->size)/3;

    for(int i=0; i<n_neuron_interp;i++){
        double ai = -5+(b_interp-a_interp)*(i)/(n_neuron_interp-1);
        double bi = 1;
        double wi = 1;

        gsl_vector_set((interp->params),3*i,1);
        gsl_vector_set((interp->params),3*i+1,bi);
        gsl_vector_set((interp->params),3*i+2,wi);
    }


    //TRAINING THE INTERPOLATER FUNCTION
    ann_train(interp,xs,ys);

    for(int i=0; i<100; i++){
        double xi = -5. +10.* (double) i/100;
        double interp_i = ann_response(interp,xi);
        double interp_int_i = ann_int(interp,-5,xi);
        double interp_diff_i = ann_diff(interp,xi);
        double exact_diff_i = gauss_diff(xi);
        double exact_int_i = gauss_int(-5,xi);
        double exact_i = gauss(xi);
        fprintf(compare,"%g %g %g %g %g %g %g\n",xi,interp_i,exact_i,interp_diff_i,exact_diff_i,
                interp_int_i,exact_int_i);
    }

    //PART C
    double a = 0;
    double b = 2*M_PI;
    int n_neuron = ((wild->params)->size)/3;
    for(int i=0; i<n_neuron;i++){
        double ai = a+(b-a)*i/(n_neuron-1);
        double bi = 1;
        double wi = 1;

        gsl_vector_set((wild->params),3*i,ai);
        gsl_vector_set((wild->params),3*i+1,bi);
        gsl_vector_set((wild->params),3*i+2,wi);
    }
    annWild_train(wild,diffeq_pow2,a,b,0,0,1);


    for(int i=0; i<100; i++){
        double xi = 2*M_PI* (double) i/100;
        double wild_i = ann_response(wild,xi);
        double exact_i = sin(xi);
        fprintf(diffeq,"%g %g %g\n",xi,wild_i,exact_i);
    }
    fclose(compare);
    fclose(diffeq);
	return 0;
}
