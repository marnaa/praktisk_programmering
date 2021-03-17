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
double qeval(qspline* s, double z) {
    int i = binsearch(s->x, z);
    double h = z - gsl_vector_get(s->x, i);
    double bi = gsl_vector_get(s->b, i);
    double ci = gsl_vector_get(s->c, i);
    double yi = gsl_vector_get(s->y, i);
    double val = yi + h * (bi + h * ci);
    return val;
}

double qint(qspline* s, double z){
    int j = binsearch(s->x,z);
    double intval = 0;
    for(int i=0; i<j; i++){
        double yi = gsl_vector_get(s->y,i);
        double bi = gsl_vector_get(s->b,i);
        double ci = gsl_vector_get(s->c,i);
        double xi = gsl_vector_get(s->x,i), xii = gsl_vector_get(s->x,i+1);
        double a = yi-bi*xi+ci*xi*xi;
        double b = bi-2*xi*ci;
        double stepval = ((a*(xii-xi)+0.5*b*(xii*xii-xi*xi)+1./3*ci*(xii*xii*xii-xi*xi*xi)));
        intval+=stepval;
    }
    double yj = gsl_vector_get(s->y,j);
    double bj = gsl_vector_get(s->b,j);
    double cj = gsl_vector_get(s->c,j);
    double xj = gsl_vector_get(s->x,j);
    double a = yj-bj*xj+cj*xj*xj;
    double b = bj-2*xj*cj;
    double stepval = a*(z-xj)+0.5*b*(z*z-xj*xj)+1./3*cj*(z*z*z-xj*xj*xj);
    intval+=stepval;
    return intval;
}

double qdiff(qspline* s, double z){
    int j = binsearch(s->x,z);
    double bj = gsl_vector_get(s->b,j);
    double cj = gsl_vector_get(s->c,j);
    double xj = gsl_vector_get(s->x,j);
    double b = bj-2*xj*cj;
    return b+2*cj*z;
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
	int n = 100;
	FILE* qxandy = fopen("out.qxy.txt","w");
    FILE* qout = fopen("out.qdata.txt","w");
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* cosy = gsl_vector_alloc(n);
	gsl_vector* siny = gsl_vector_alloc(n);
	for(int i=0; i<n; i++){
	    double xi = ((double) i)/n*2*M_PI;
	    gsl_vector_set(x,i,xi);
	    gsl_vector_set(siny,i, sin(xi));
	    gsl_vector_set(cosy,i, cos(xi));
	}
	qspline* s = qspline_make(x,siny);
	int count=0;
	printf("her");
	while(count<100){
	    count+=1;
        double z = (double)rand()/(double) RAND_MAX*gsl_vector_get(x,n-2);
		fprintf(qout,"%g %g %g %g\n", z, qeval(s,z),-qint(s,z)+1,qdiff(s,z));
	}
    for(int i=0; i<n;i++) {
        fprintf(qxandy, "%g %g %g\n", gsl_vector_get(x, i),
                gsl_vector_get(siny, i), gsl_vector_get(cosy, i));
    }
	qspline_free(s);
	gsl_vector_free(x);
	gsl_vector_free(cosy);
	gsl_vector_free(siny);
	return 0;
}