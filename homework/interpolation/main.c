#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <assert.h>

int binsearch(gsl_vector* x, double z){
	/* locates the interval for z by bisection */
	int n = (*x).size;
	assert(gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,n-2));
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
		}
	return i;
	}

double linterp(gsl_vector* x, gsl_vector* y, double z){
	int i = binsearch(x, z);
	double a = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/
	(gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
	double fz = a*(z-gsl_vector_get(x,i))+gsl_vector_get(y,i);
	return fz;
	}
double linterp_integ(gsl_vector* x, gsl_vector* y, double z){
	int i = binsearch(x, z);
        for(int j=0; j<i; j++){
	
	}
	double a = (gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/
        (gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
        double x0 = gsl_vector_get(x,0);
	double intfz = 0;
	double xi = gsl_vector_get(x,i);
	double yi = gsl_vector_get(y,i);
        intfz += 0.5*a*(z*z-x0*x0)-a*xi*(z-x0)+yi*(z-x0);
	}
	return intfz;
	}

int main(){
	int N=10;
	gsl_vector* x = gsl_vector_alloc(N);
	gsl_vector* y = gsl_vector_alloc(N);
	for(int i=0; i<N; i++){
		gsl_vector_set(x,i,i);
		gsl_vector_set(y,i,i*i);
	}
	double z = 7.5;
	double interz = linterp(x,y,z);
	double integralz = linterp_integ(x, y, z);
	printf("x:\n"); 
	gsl_vector_fprintf(stdout,x,"%g");
	printf("y:\n");
	gsl_vector_fprintf(stdout,y,"%g");
	printf("interpolated z val: %g\n",interz);
	printf("integral x0 to z: %g\n", integralz);
	gsl_vector_free(x);
	gsl_vector_free(y);
	return 0;
}
