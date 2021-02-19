#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>

double Erf(double);
double Gamma(double);

int main(){
	double xmin=-2.1, xmax=4.2;
	double x=xmin;
	while(x<=xmax){
	printf("%10g %10g %10g %10g %10g %10g %10g\n",x, erf(x), gsl_sf_erf(x), Erf(x),tgamma(x), gsl_sf_gamma(x), Gamma(x));
	x+=1.0/8;
}

	return 0;
}

