#include <tgmath.h>
#include <stdio.h>
int main(){
	double gam=tgamma(5.0);
	double bes=j1(0.5);
	complex sqqq=csqrt(-2);
	complex eipi=cexp(I*M_PI);
	complex ei=cexp(I);
	complex ie=cpow(I,M_E);
	complex ii=cpow(I,I);
	printf("gammafunc(5)= %g\n", gam);
	printf("j1(0.5)= %g\n", bes);
	printf("sqrt-2= %g+I%g\n", creal(sqqq),cimag(sqqq));
	printf("e^ipi: %g+I%g\n", creal(eipi),cimag(eipi));
	printf("e^i: %g+I%g\n", creal(ei),cimag(ei));
	printf("i^e: %g+I%g\n", creal(ie),cimag(ie));
	printf("i^i: %g+I%g\n", creal(ii),cimag(ii));
// precision of the diferent formats
	float x_f=1.f/9;
	double x_d=1./9;
	long double x_ld=1.L/9;
	printf("float: %.25g \n",x_f);
	printf("double: %.25lg\n", x_d);
	printf("long double: %.365Lg\n",x_ld);
	return 0;
}
