//
// Created by Martin on 21-05-2021.
//
#include "min_declarations.h"
#include <assert.h>

double intrun(double f(double x, ann* network), double f2, double f3, double a, double b, ann* network,
              double acc, double eps,double nrec, int *eta, double* error){
    //finding the higher order value
    assert(nrec<10000);
    double x1 = a+1./6*(b-a);
    double x4 = a+5./6*(b-a);
    double f1=f(x1,network) , f4 = f(x4,network);
    double Q = (b-a)/(6.)*(2*f1+f2+f3+2*f4);
    double q = (b-a)/(4.)*(f1+f2+f3+f4);
    double tol = acc + eps*fabs(Q);
    double err = fabs(Q-q);
    if(err< tol){
        if(*eta<nrec){
            *eta=nrec;
        }
        *error += err*err;
        return Q;
    }
    else{
        double Q1 = intrun(f,f1,f2,a,(a+b)/2,network, acc/sqrt(2), eps,nrec+1,eta,error);
        double Q2 =intrun(f,f3,f4,(a+b)/2,b,network, acc/sqrt(2), eps,nrec+1,eta,error);
        return Q1+Q2;
    }
}
double integrater(double f(double x, ann* network), double a,
                  double b, ann* network, double acc, double eps, int* eta,double* error){
    double f2 = f(a+2.*(b-a)/6.,network), f3 = f(a+3.*(b-a)/6,network);
    int nrec = 0;
    *error = 0.;
    double val = intrun(f,f2,f3,a,b,network,acc,eps,nrec,eta,error);
    *error = sqrt(*error);
    return val;

}