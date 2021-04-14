#include <math.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>
#include <gsl/gsl_integration.h>


//TO DO:
/* double integrate(double f(double), double a, double b, double δ, double ε)
{
double Q = higher_order_rule, q = lower_order_rule, err=|Q-q|;
if (err < δ+ε|Q|) return Q;
else return integrate(f,a,(a+b)/2,δ/√2,ε)+
            integrate(f,(a+b)/2,b,δ/√2,ε);
}
 But open
 */
//Open integrater using polynomials
double intrun(double f(double), double f2, double f3, double a, double b, double acc, double eps,double nrec){
    //finding the higher order value
    assert(nrec<10000);
    double x1 = a+1./6*(b-a);
    double x4 = a+5./6*(b-a);
    double f1=f(x1) , f4 = f(x4);
    double Q = (b-a)/(6.)*(2*f1+f2+f3+2*f4);
    double q = (b-a)/(4.)*(f1+f2+f3+f4);
    double tol = acc + eps*fabs(Q);
    double err = fabs(Q-q);
    if(err< tol){
        return Q;
    }
    else{
        double Q1 = intrun(f,f1,f2,a,(a+b)/2,acc/sqrt(2),eps,nrec+1);
        double Q2 =intrun(f,f3,f4,(a+b)/2,b,acc/sqrt(2),eps,nrec+1);
        return Q1+Q2;
    }
}
double integrater(double f(double), double a,
                  double b, double acc, double eps){
    double f2 = f(a+2.*(b-a)/6.), f3 = f(a+3*(b-a)/6);
    printf("%g %g\n",f2,f3);
    int nrec = 0;
    return intrun(f,f2,f3,a,b,acc,eps,nrec);
}

double trapezrun(double h(double f(double), double A, double B, double x),
              double f(double), double h2, double h3, double a, double b,double A, double B, double acc, double eps,double nrec) {
    //finding the higher order value
    assert(nrec < 10000);
    double x1 = a + 1. / 6 * (b - a);
    double x4 = a + 5. / 6 * (b - a);
    double h1 = h(f,A,B,x1), h4 = h(f,A,B,x4);
    double Q = (b - a) / (6.) * (2 * h1 + h2 + h3 + 2 * h4);
    double q = (b - a) / (4.) * (h1 + h2 + h3 + h4);
    double tol = acc + eps * fabs(Q);
    double err = fabs(Q - q);
    if (err < tol) {
        return Q;
    } else {
        double Q1 = trapezrun(h, f,h1, h2, a, (a + b) / 2,A,B, acc / sqrt(2), eps, nrec + 1);
        double Q2 = trapezrun(h, f,h3, h4, (a + b) / 2,b, A,B, acc / sqrt(2), eps, nrec + 1);
        return Q1 + Q2;
    }
}
//x -> ((b-a)x+(b+a))/2
double h(double f (double), double a, double b,double x){
    double gcosx=((b-a)/2*cos(x)+(b+a)/2);
    double hx = f(gcosx)*sin(x)*(b+a)/2;
    return hx;
}
double CC_integrater(double f(double), double a,
                  double b, double acc, double eps){
    int nrec = 0;
    double x2 = 2./6*M_PI;
    double x3 = 4./6*M_PI;
    double f2=h(f,a,b,x2);
    double f3=h(f,a,b,x3);
    double Q = trapezrun(h,f,f2,f3,0.,M_PI,a,b,acc,eps,nrec);
    return Q;
}


double Fa1(double x){
    return sqrt(x);
}
double Fa2(double x){
    return 4.*sqrt(1-x*x);
}

double Fb1(double x){
    return 1/sqrt(x);
}
double Fb2(double x){
    return log(x)/sqrt(x);
}


int main(){
    FILE* pi_compare = fopen("out.piCompare.txt","w");
    double a = 0., b=1., acc = 0.0001, eps = 0.;
    double fa1 = integrater(Fa1,a,b,acc,eps);
    double fa2 = integrater(Fa2,a,b,acc,eps);
    double ha1 = CC_integrater(Fa1,a,b,acc,eps);
    double ha2 = CC_integrater(Fa2,a,b,acc,eps);
    double hb1 = CC_integrater(Fb1,a,b,acc,eps);
    double hb2 = CC_integrater(Fb2,a,b,acc,eps);
    printf("No var change: int(sqrt(x),0..1)= %.16f\n",fa1);
    printf("var change: int(sqrt(x),0..1)= %.16f ",ha1);
    printf("var change: int(1/sqrt(x),0..1)= %.16f int(ln(x)/sqrt(x),0..1)=%.16f\n",hb1,hb2);
    printf("acc= %g ; eps = %g",acc,eps);
    printf( "Int: int(4*sqrt(1-x^2),0..1)=%.16f\n", fa2);
    printf("CC-int: int(4*sqrt(1-x^2),0..1)=%.16f\n",ha2);
    for(int i=0; i<10; i++){
        for
    }
    return 0;
}
