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
double trapezrun(double h(double f(double x, ann* network), double A, double B, double x, ann* network),
                 double f(double, ann* network), double h2, double h3, double a, double b,double A, double B, ann* network,
                 double acc, double eps,double nrec, int* eta, double* error) {
    //finding the higher order value
    assert(nrec < 10000);
    double x1 = a + 1. / 6 * (b - a);
    double x4 = a + 5. / 6 * (b - a);
    double h1 = h(f, A, B, x1, network), h4 = h(f, A, B, x4, network);
    double Q = (b - a) / (6.) * (2 * h1 + h2 + h3 + 2 * h4);
    double q = (b - a) / (4.) * (h1 + h2 + h3 + h4);
    double tol = acc + eps * fabs(Q);
    double err = fabs(Q - q);
    if (err < tol) {
        if(*eta<nrec){*eta = nrec;
        }
        *error+=err*err;
        return Q;
    } else {
        double Q1 = trapezrun(h, f,h1, h2, a, (a + b) / 2,A,B, network, acc / sqrt(2), eps, nrec+1,eta,error);
        double Q2 = trapezrun(h, f,h3, h4, (a + b) / 2,b, A,B,network, acc / sqrt(2), eps, nrec+1,eta,error);
        return Q1 + Q2;
    }
}
//x -> ((b-a)x+(b+a))/2
double h(double f (double x, ann* network), double a, double b,double x,ann* network){
    double gcosx=((b-a)/2*cos(x)+(b+a)/2);
    double hx = f(gcosx, network)*sin(x)*(b+a)/2;
    return hx;
}

double infinf(double f (double x, ann* network),double a, double b, double x, ann* network) {
    double gcosx = (cos(x)/(1-cos(x)*cos(x)));
    double hx = f(gcosx, network) * sin(x) * (1+cos(x)*cos(x))/pow((1-cos(x)*cos(x)),2);
    return hx;
}
double bIsInf(double f (double x, ann* network), double a, double b,double x,ann* network) {
    double gcosx = a+(cos(x)+1)/(1-cos(x)) ;
    double hx = f(gcosx, network) * sin(x) * 2/pow(cos(x)-1,2);
    return hx;
}
double aIsInf(double f (double x, ann* network), double a, double b,double x,ann* network) {
    double gcosx = b-(1-cos(x))/(1+cos(x)) ;
    double hx = f(gcosx, network) * sin(x) * 2/pow(cos(x)+1,2);
    return hx;
}
double CC_integrater(double f(double x, ann* network), double a,
                     double b, ann* network, double acc, double eps,int* eval, double* error){
    double Q, x2, x3, f2, f3;
    int nrec = 0;
    *error = 0.;
    if(isinf(a)==-1 && isinf(b)==1){
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = infinf(f,0.,0.,x2,network);
        f3 = infinf(f,0.,0.,x3,network);
        Q = trapezrun(infinf, f, f2, f3, 0., M_PI, a, b, network, acc, eps, nrec, eval, error);

    }
    else if(isinf(a)==0 && isinf(b)==1){
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = bIsInf(f,a,0.,x2,network);
        f3 = bIsInf(f,a,0.,x3,network);
        Q = trapezrun(bIsInf, f, f2, f3, 0., M_PI, a, b, network, acc, eps, nrec, eval, error);
    }
    else if(isinf(a)==-1 && isinf(b)==0){
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = aIsInf(f,0.,b,x2,network);
        f3 = aIsInf(f,0.,b,x3,network);
        Q = trapezrun(aIsInf, f, f2, f3, 0., M_PI, a, b, network, acc, eps, nrec, eval, error);
    }
    else {
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = h(f, a, b, x2, network);
        f3 = h(f, a, b, x3, network);
        Q = trapezrun(h, f, f2, f3, 0., M_PI, a, b,network, acc, eps, nrec, eval, error);
    }
    *error = sqrt(*error);
    return Q;
}