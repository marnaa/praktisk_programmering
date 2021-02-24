#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

typedef struct {double z; double n;} besselparam;


double logoversqrtint(double x, void* params){
	double logoversqrtint = log(x)/(sqrt(x));
	return logoversqrtint;
 }

double erfint(double x, void* params){
	double erfint = 2/sqrt(M_PI)*exp(-1*pow(x,2));
	return erfint;
}

double besselint(double x,void* params){
	besselparam * par = (besselparam*)params;
	double z = (par -> z);
	double n = (par -> n);
	double besselint = 1/M_PI*cos(x*n-z*sin(x));
	return besselint;
}

double logoversqrt(){
	gsl_function F;
	F.function = &logoversqrtint;
	int limit = 999;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	double a=0, b=1, acc=1e-6, eps=1e-6, result, error;
	gsl_integration_qag(&F, a, b, acc, eps, limit, 3, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}
double Erf(double z){
        gsl_function F;
        F.function = &erfint;
        int limit = 999;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
        double a=0, acc=1e-6, eps=1e-6, result, error;
        gsl_integration_qag(&F, 0, z, acc, eps, limit, 3, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}
double Bessel(double z, double n){
        gsl_function F;
	besselparam params = {z, n};
        F.function = &besselint;
	F.params = &params;
        int limit = 999;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
        double a=0, b=M_PI, acc=1e-6, eps=1e-6, result, error;
        gsl_integration_qag(&F, a, b, acc, eps, limit, 3, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
}

int main(){
	FILE* outA=fopen("out.valueExerciseA.txt","w");
	FILE* outErf=fopen("out.erfValues.txt","w");
	FILE* outBessel=fopen("out.besselValues.txt","w");
	double logsqrt = logoversqrt();
	fprintf(outA,"int(log(x)/sqrt(x), x=0..1) = %g\n", logsqrt);
	for(double x=-5; x<=5; x+=1.0/8){
		fprintf(outErf, "%10g %10g\n",x,Erf(x));
	}
	for(double x=-6; x<=6; x+=1.0/8){
		fprintf(outBessel, "%10g %10g %10g %10g\n",x, Bessel(x,0.0),
		Bessel(x,1.0), Bessel(x,2.0));
	}
	return 0;
}
