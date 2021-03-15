#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

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
	double intfz = 0.;
        for(int j=0; j<i; j++){
        double xj = gsl_vector_get(x,j), xjj = gsl_vector_get(x,j+1);
        double yj = gsl_vector_get(y,j), yjj = gsl_vector_get(y,j+1);
	double a = (yjj-yj)/(xjj-xj);
        intfz += 0.5*a*(xjj*xjj-xj*xj)-a*xj*(xjj-xj)+yj*(xjj-xj);
	}
	double xi = gsl_vector_get(x,i), xii = gsl_vector_get(x,i+1);
        double yi = gsl_vector_get(y,i), yii = gsl_vector_get(y,i+1);
        double a = (yii-yi)/(xii-xi);
	intfz += 0.5*a*(z*z-xi*xi)-a*xi*(z-xi)+yi*(z-xi);
	return intfz;
	}

int main(){
	FILE* xandy = fopen("out.xy.txt","w");
	FILE* out10 = fopen("out.10.txt","w");
	FILE* out100 = fopen("out.100.txt","w");
	FILE* out500 = fopen("out.500.txt","w");
	double xmax = 5;
	FILE* files[]={out10, out100, out500};
	double Ns[]={10, 100, 500};
	for(int n = 0; n<3; n++){
		double stepsize = xmax/(double)Ns[n];
		gsl_vector* x = gsl_vector_alloc(Ns[n]);
		gsl_vector* y = gsl_vector_alloc(Ns[n]);
		double xa[(int)Ns[n]];
		double ya[(int)Ns[n]];
		for(int i=0; i<Ns[n]; i++){
			gsl_vector_set(x,i,stepsize*i);
			gsl_vector_set(y,i,(sin(stepsize*i)));
			xa[i] = stepsize*i;
			ya[i] = sin(stepsize*i);
		}
		gsl_interp *intel = gsl_interp_alloc(gsl_interp_linear,(int)Ns[n]);
		gsl_interp_init(intel, xa, ya, (int) Ns[n]);
		gsl_interp_accel *accel = gsl_interp_accel_alloc();
		int count=0;
		while(count<100){
			count+=1;
			double z = (double)rand()/(double) RAND_MAX*
			gsl_vector_get(x,Ns[n]-2);
			double interz = linterp(x,y,z);
			double integralz = linterp_integ(x, y, z);
			double gsl_interz = gsl_interp_eval(intel,xa,ya,z,accel);
			double gsl_ing = gsl_interp_eval_integ(intel,xa,ya,xa[0]
			,z,accel);
			fprintf(files[n],"%g %g %g %g %g\n", z, 
			interz, -integralz+1, gsl_interz, -gsl_ing+1);
		}
		gsl_interp_free(intel);
		gsl_vector_free(x);
		gsl_vector_free(y);
		}
	for(int i=0; i<=200;i++){
		fprintf(xandy, "%g %g %g\n",(double)i*xmax/200., 
		sin((double)i*xmax/200.),cos((double) i*xmax/200.));
	}
	return 0;
}
