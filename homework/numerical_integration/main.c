#include <math.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>
#include <gsl/gsl_integration.h>



//Open integrater using polynomials
double intrun(double f(double), double f2, double f3, double a, double b,
              double acc, double eps,double nrec, int *eta, double* error){
    //finding the higher order value
    assert(nrec<10000);
    double x1 = a+1./6*(b-a);
    double x4 = a+5./6*(b-a);
    double f1=f(x1) , f4 = f(x4);
    double Q = (b-a)/(6.)*(2*f1+f2+f3+2*f4);
    double q = (b-a)/(4.)*(f1+f2+f3+f4);
    double tol = acc + eps*fabs(Q);
    double err = fabs(Q-q);
    //Assesment of the error compared to the tolerance
    //If error okay, if not further splitting
    if(err< tol){
        if(*eta<nrec){
            *eta=nrec;
        }
        *error += err*err;
        return Q;
    }
    else{
        double Q1 = intrun(f,f1,f2,a,(a+b)/2,acc/sqrt(2),eps,nrec+1,eta,error);
        double Q2 =intrun(f,f3,f4,(a+b)/2,b,acc/sqrt(2),eps,nrec+1,eta,error);
        return Q1+Q2;
    }
}

//Driver for the integration routine
double integrater(double f(double), double a,
                  double b, double acc, double eps, int* eta,double* error){
    double f2 = f(a+2.*(b-a)/6.), f3 = f(a+3.*(b-a)/6);
    int nrec = 0;
    *error = 0.;
    double val = intrun(f,f2,f3,a,b,acc,eps,nrec,eta,error);
    *error = sqrt(*error);
    return val;

}

//Trapez run only different from intrun in llowing subtitutions through h
double trapezrun(double h(double f(double), double A, double B, double x),
              double f(double), double h2, double h3, double a, double b,double A, double B,
              double acc, double eps,double nrec, int* eta, double* error) {
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
        if(*eta<nrec){*eta = nrec;
        }
        *error+=err*err;
        return Q;
    } else {
        double Q1 = trapezrun(h, f,h1, h2, a, (a + b) / 2,A,B, acc / sqrt(2), eps, nrec+1,eta,error);
        double Q2 = trapezrun(h, f,h3, h4, (a + b) / 2,b, A,B, acc / sqrt(2), eps, nrec+1,eta,error);
        return Q1 + Q2;
    }
}


//x -> ((b-a)x+(b+a))/2
//Function for doing the clebshaw-cutis substitution
double h(double f (double), double a, double b,double x){
    double gcosx=((b-a)/2*cos(x)+(b+a)/2);
    double hx = f(gcosx)*sin(x)*(b+a)/2;
    return hx;
}


//Functions for when different limits are infinity
double infinf(double f (double),double a, double b, double x) {
    double gcosx = (cos(x)/(1-cos(x)*cos(x)));
    double hx = f(gcosx) * sin(x) * (1+cos(x)*cos(x))/pow((1-cos(x)*cos(x)),2);
    return hx;
}
double bIsInf(double f (double), double a, double b,double x) {
    double gcosx = a+(cos(x)+1)/(1-cos(x)) ;
    double hx = f(gcosx) * sin(x) * 2/pow(cos(x)-1,2);
    return hx;
}
double aIsInf(double f (double), double a, double b,double x) {
    double gcosx = b-(1-cos(x))/(1+cos(x)) ;
    double hx = f(gcosx) * sin(x) * 2/pow(cos(x)+1,2);
    return hx;
}


//Driver for clebshaw-curtis integration with possibility for things being infinite
double CC_integrater(double f(double), double a,
                  double b, double acc, double eps,int* eval, double* error){
    double Q, x2, x3, f2, f3;
    int nrec = 0;
    *error = 0.;
    if(isinf(a)==-1 && isinf(b)==1){
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = infinf(f,0.,0.,x2);
        f3 = infinf(f,0.,0.,x3);
        Q = trapezrun(infinf, f, f2, f3, 0., M_PI, a, b, acc, eps, nrec, eval, error);

    }
    else if(isinf(a)==0 && isinf(b)==1){
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = bIsInf(f,a,0.,x2);
        f3 = bIsInf(f,a,0.,x3);
        Q = trapezrun(bIsInf, f, f2, f3, 0., M_PI, a, b, acc, eps, nrec, eval, error);
    }
    else if(isinf(a)==-1 && isinf(b)==0){
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = aIsInf(f,0.,b,x2);
        f3 = aIsInf(f,0.,b,x3);
        Q = trapezrun(aIsInf, f, f2, f3, 0., M_PI, a, b, acc, eps, nrec, eval, error);
    }
    else {
        x2 = 2. / 6 * M_PI;
        x3 = 4. / 6 * M_PI;
        f2 = h(f, a, b, x2);
        f3 = h(f, a, b, x3);
        Q = trapezrun(h, f, f2, f3, 0., M_PI, a, b, acc, eps, nrec, eval, error);
    }
    *error = sqrt(*error);
    return Q;
}


//Testing functions
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

double ekspM(double x){
    return exp(-x);
}
double ekspP(double x){
    return exp(x);
}
double ekspSq(double x){
    return exp(-(x*x));
}


int main(){
    FILE* pi_compare = fopen("out.piCompare.txt","w");
    FILE* inf_compare = fopen("out.infCompare.txt","w");
    int lol, lol_CC;
    double error_CC, error_norm;
    double a = 0., b=1., acc = 0.0001, eps = 0.;


    //Integrating the test functions with the normal integrator and CC
    double fa1 = integrater(Fa1,a,b,acc,eps,&lol,&error_norm);
    double fa2 = integrater(Fa2,a,b,acc,eps,&lol,&error_norm);
    double ha1 = CC_integrater(Fa1,a,b,acc,eps,&lol_CC,&error_CC);
    double ha2 = CC_integrater(Fa2,a,b,acc,eps,&lol_CC,&error_CC);
    double hb1 = CC_integrater(Fb1,a,b,acc,eps,&lol_CC,&error_CC);
    double hb2 = CC_integrater(Fb2,a,b,acc,eps,&lol_CC,&error_CC);


    printf("No var change: int(sqrt(x),0..1)= %.16f\n",fa1);
    printf("var change: int(sqrt(x),0..1)= %.16f \n",ha1);
    printf("var change: int(1/sqrt(x),0..1)= %.16f \n and int(ln(x)/sqrt(x),0..1)=%.16f\n",hb1,hb2);
    printf("AT acc= %g ; eps = %g TRYING THE TWO DIFFERENT ROUTINES\n",acc,eps);
    printf( "No var: int(4*sqrt(1-x^2),0..1)=%.16f\n", fa2);
    printf("Var change: int(4*sqrt(1-x^2),0..1)=%.16f and error = %.16f\n",ha2,error_CC);

    //Testing against gsl_integrator etc.
    double loop_acc[]={0.0001,0.0005,0.001,0.005,0.01,0.05,0.06,0.07,0.08};
    gsl_integration_workspace * w
            = gsl_integration_workspace_alloc (10000);
    gsl_integration_workspace * v
            = gsl_integration_workspace_alloc (100000);
    gsl_integration_workspace * u
            = gsl_integration_workspace_alloc (10000);
    gsl_integration_workspace * p
            = gsl_integration_workspace_alloc (10000);


    double err, err_CC;
    int eta, eta_CC;
    int ekspM_eval, ekspP_eval, ekspSq_eval;
    double ekspM_error, ekspP_error, ekspSq_error;
    double ekspM_error_gsl, ekspP_error_gsl, ekspSq_error_gsl;
    double ekspM_gsl, ekspP_gsl, ekspSq_gsl;
    double result, err_gsl;
    gsl_function F,EKSP_SQ;
    F.function = &Fa2;
    EKSP_SQ.function = &ekspSq;
    printf("\n\nERROR BELONGING TO int(4*sqrt(1-x*x),0..1) FOR THE PICOMPARE.PNG\n\n");
    for( int i =0; i<4; i++){
        eta=0;
        eta_CC=0;
        gsl_integration_qags (&F, 0, 1, loop_acc[i], 0, 10000,
                              w, &result, &err_gsl);

        double compare = fabs(M_PI-integrater(Fa2,a,b,loop_acc[i],0.,&eta,&err));
        double CC_compare = fabs(M_PI-CC_integrater(Fa2,a,b,loop_acc[i],0.,&eta_CC,&err_CC));
        fprintf(pi_compare,"%g %g %g %g %i %i\n",loop_acc[i],compare,
                CC_compare,fabs(M_PI-result),eta, eta_CC);

        printf("Error at acc %g: _CC= %g _normal = %g  _gsl= %g\n",loop_acc[i],err_CC,err, err_gsl);


        gsl_integration_qagiu(&EKSP_SQ,0,loop_acc[i],0.,10000,u,&ekspM_gsl,&ekspM_error_gsl);
        gsl_integration_qagil(&EKSP_SQ,0,loop_acc[i],0.,10000,p,&ekspP_gsl,&ekspP_error_gsl);
        gsl_integration_qagi(&EKSP_SQ,loop_acc[i],0.,10000,v,&ekspSq_gsl,&ekspSq_error_gsl);


        double M = CC_integrater(ekspSq,0,INFINITY,loop_acc[i],0.,&ekspM_eval,&ekspM_error);
        double P = CC_integrater(ekspSq,-INFINITY,0,loop_acc[i],0.,&ekspP_eval,&ekspP_error);
        double Sq = CC_integrater(ekspSq,-INFINITY,INFINITY,loop_acc[i],0.,&ekspSq_eval,&ekspSq_error);

        fprintf(inf_compare,"%g %g %g %g %g %g %g %i %i %i %g %g %g %g %g %g \n",
                loop_acc[i], M,P,Sq,ekspM_gsl, ekspP_gsl,
                ekspSq_gsl, ekspM_eval,ekspP_eval,
                ekspSq_eval, ekspM_error,ekspP_error, ekspSq_error,ekspM_error_gsl, ekspP_error_gsl, ekspSq_error_gsl );

    }
    gsl_integration_workspace_free(w);
    gsl_integration_workspace_free(v);
    gsl_integration_workspace_free(p);
    gsl_integration_workspace_free(u);
    return 0;
}
