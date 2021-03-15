#include <math.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//Declaring binsearch found in main
int binsearch(gsl_vector* x, double z){
        /* locates the interval for z by bisection */
        int n = (*x).size;
        int i=0, j=n-1;
        while(j-i>1){
                int mid=(i+j)/2;
                if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
                }
	return i;
	}

//Defining the qsline structure
typedef struct{gsl_vector *x, *y, *b, *c;} qspline;

qspline* qspline_make(gsl_vector* x, gsl_vector* y){
	//Initially allocating memory for spline
	int n = x->size;
	int i;
	qspline* s = (qspline*)malloc(sizeof(qspline));
	s->b = gsl_vector_alloc(n-1);
	s->c = gsl_vector_alloc(n-1);
	s->x = gsl_vector_alloc(n);
	s->y = gsl_vector_alloc(n);
	for(i = 0; i<n; i++){
		double xi = gsl_vector_get(x,i);
		double yi = gsl_vector_get(y,i);
		gsl_vector_set(s->x,i,xi);
		gsl_vector_set(s->y,i,yi);
	}
	//Finding p_i and h_i (really, dx_i)
	gsl_vector* p = gsl_vector_alloc(n-1);
	gsl_vector* h = gsl_vector_alloc(n-1);
	for(i=0; i<n-1; i++){
		double dx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
		double dy = gsl_vector_get(y,i+1)-gsl_vector_get(y,i);
		gsl_vector_set(h,i,dx);
		gsl_vector_set(p,i,dy/dx);
	}
	//Using the one DOF to set c_0 = 0
	gsl_vector_set(s->c,0,0.0);
	//upwards recursion
	for(i=0; i<n-2; i++){
		//Regretting doing this with vectors and not arrays
		double ci = gsl_vector_get(s->c,i);
		double dp = gsl_vector_get(p,i+1)-gsl_vector_get(p,i);
		double hi = gsl_vector_get(h,i);
		double hii = gsl_vector_get(h,i+1);
		double cii = (dp-ci*hi)/hii;
		gsl_vector_set(s->c,i+1,cii);
	}

	//downwards recursion
	gsl_vector_set(s->c,n-2,gsl_vector_get(s->c,n-2)/2);
	for(i=n-3; i>=0; i--){
                double cii = gsl_vector_get(s->c,i+1);
                double dp = gsl_vector_get(p,i+1)-gsl_vector_get(p,i);
                double hi = gsl_vector_get(h,i);
                double hii = gsl_vector_get(h,i+1);
                double ci = (dp-cii*hii)/hi;
                gsl_vector_set(s->c,i,ci);
        }
	for(i=0; i<n-1; i++){
		double pi = gsl_vector_get(p,i);
		double ci = gsl_vector_get(s->c,i);
		double hi = gsl_vector_get(h,i);
		gsl_vector_set(s->b,i,pi-ci*hi);
	}
	return s;
}

//FUNCTION FOR EVALUATING THE SPLINE IN A GIVEN Z
double qeval(qspline* s, double z){
        int i = binsearch(s->x, z);
        double h = z - gsl_vector_get(s->x,i);
	double bi = gsl_vector_get(s->b,i);
	double ci = gsl_vector_get(s->c,i);
	double yi = gsl_vector_get(s->y,i);
	double val = yi+h*(bi+h*ci);
	return val;
        }
//FREEING THE ALLOCATED MEMORY
void qspline_free(qspline* s){
	gsl_vector_free(s->x);
	gsl_vector_free(s->y);
	gsl_vector_free(s->c);
	gsl_vector_free(s->b);
	free(s);
	}
int main(){
	int n = 5;
	gsl_vector* x = gsl_vector_alloc(5);
	gsl_vector* y = gsl_vector_alloc(5);
	gsl_vector* yx = gsl_vector_alloc(5);
	gsl_vector* yxx = gsl_vector_alloc(5);
	for(int i=0; i<n; i++){
	gsl_vector_set(x,i, (double)i);
	gsl_vector_set(y,i, 1.0);
	gsl_vector_set(yx,i,(double)i);
	gsl_vector_set(yxx,i,(double)i*i);
        }
	qspline* s = qspline_make(x,y);
	qspline* sx = qspline_make(x,yx);
	qspline* sxx = qspline_make(x,yxx);
	for(int i=0; i<n; i++){
		double z = (double) i+0.3;
		printf("%g %g %g\n", qeval(s,z),qeval(sx,z),qeval(sxx,z));
	}
	qspline_free(s);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(yx);
	gsl_vector_free(yxx);
	return 0;
}
