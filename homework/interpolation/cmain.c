#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
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
typedef struct{gsl_vector *x, *y, *b, *c, *d;} cspline;

cspline* cspline_make(gsl_vector* x, gsl_vector* y){
	//Initially allocating memory for spline
	int n = x->size;
	int i;
	cspline* s = (cspline*)malloc(sizeof(cspline));
	s->b = gsl_vector_calloc(n);
	s->c = gsl_vector_calloc(n-1);
	s->d = gsl_vector_calloc(n-1);
	s->x = gsl_vector_calloc(n);
	s->y = gsl_vector_calloc(n);
	for(i = 0; i<n; i++){
		double xi = gsl_vector_get(x,i);
		double yi = gsl_vector_get(y,i);
		gsl_vector_set(s->x,i,xi);
		gsl_vector_set(s->y,i,yi);
	}
	gsl_vector* p = gsl_vector_calloc(n-1);
	gsl_vector* h = gsl_vector_calloc(n-1);
	gsl_vector* D = gsl_vector_calloc(n);
	gsl_vector* Q = gsl_vector_calloc(n-1);
	gsl_vector* B = gsl_vector_calloc(n);
	//Finding p_i and h_i (really, dx_i)
	for(i=0; i<n-1; i++){
		double dx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
		double dy = gsl_vector_get(y,i+1)-gsl_vector_get(y,i);
		gsl_vector_set(h,i,dx);
		gsl_vector_set(p,i,dy/dx);
	}
	//Finding values for D, Q and B
	gsl_vector_set(D,0,2.);
	for(i=0; i<n-2; i++){
		double hi = gsl_vector_get(h,i);
		double hii = gsl_vector_get(h,i+1);
		gsl_vector_set(D,i+1,2*hi/hii+2);
	}
	gsl_vector_set(D,n-1,2.);
	gsl_vector_set(Q,0,1.);
        for(i=0; i<n-2; i++){
                double hi = gsl_vector_get(h,i);
                double hii = gsl_vector_get(h,i+1);
                gsl_vector_set(Q,i+1,hi/hii);
        }
	for(i=0; i<n-2; i++){
		double hi = gsl_vector_get(h,i);
		double hii = gsl_vector_get(h,i+1);
		double pi = gsl_vector_get(p,i);
		double pii = gsl_vector_get(p,i+1);
		gsl_vector_set(B,i+1,3*(pi+pii*hi/hii));
	}
	gsl_vector_set(B,0,3.*gsl_vector_get(p,0));
	gsl_vector_set(B,n-1,3.*gsl_vector_get(p,n-2));
	for(i=1; i<n; i++){
		double di = gsl_vector_get(D,i);
		double di_i = gsl_vector_get(D,i-1);
		double qi_i = gsl_vector_get(Q,i-1);
		gsl_vector_set(D,i,di-qi_i/di_i);
		double bi = gsl_vector_get(B,i);
                double bi_i = gsl_vector_get(B,i-1);
                gsl_vector_set(B,i,bi-bi_i/di_i);
	}
	gsl_vector_set(s->b,n-1,gsl_vector_get(B,n-1)/(gsl_vector_get(D,n-1)));
	for(i=n-2; i>=0; i--){
		double Bi = gsl_vector_get(B,i);
		double Qi = gsl_vector_get(Q,i);
		double Di = gsl_vector_get(D,i);
		double bii = gsl_vector_get(s->b,i+1);
		gsl_vector_set(s->b,i,(Bi-Qi*bii)/Di);
	}
	for(i=0; i<n-1; i++){
		double hi = gsl_vector_get(h,i);
                double pi = gsl_vector_get(p,i);
		double bi = gsl_vector_get(s->b,i);
		double bii = gsl_vector_get(s->b,i+1);
		gsl_vector_set(s->c,i,(-2*bi-bii+3*pi)/hi);
		gsl_vector_set(s->d,i,(bi+bii-2*pi)/hi/hi);
	}
	return s;
}



//FUNCTION FOR EVALUATING THE SPLINE IN A GIVEN Z
double ceval(cspline* s, double z) {
    int i = binsearch(s->x, z);
    double h = z - gsl_vector_get(s->x, i);
    double bi = gsl_vector_get(s->b, i);
    double ci = gsl_vector_get(s->c, i);
    double di = gsl_vector_get(s->d, i);
    double val = gsl_vector_get(s->y, i) + h * (bi + h * (ci + h * di));
    return val;
}

double cint(cspline* s, double z){
    int j = binsearch(s->x,z);
    double intval = 0;
    for(int i=0; i<j; i++){
        double yi = gsl_vector_get(s->y,i);
        double bi = gsl_vector_get(s->b,i);
        double ci = gsl_vector_get(s->c,i);
        double di = gsl_vector_get(s->d,i);
        double xi = gsl_vector_get(s->x,i), xii = gsl_vector_get(s->x,i+1);
        double a = yi-bi*xi+ci*xi*xi-di*pow(xi,3);
        double b = bi-2*xi*ci+3*di*xi*xi;
        double c = ci-3*xi*di;
        double stepval = a*(xii-xi)+0.5*b*(xii*xii-xi*xi)+1./3*c*(xii*xii*xii-xi*xi*xi)+1./4*di*(pow(xii,4)-pow(xi,4));
        intval+=stepval;
    }
    double yj = gsl_vector_get(s->y,j);
    double bj = gsl_vector_get(s->b,j);
    double cj = gsl_vector_get(s->c,j);
    double dj = gsl_vector_get(s->d,j);
    double xj = gsl_vector_get(s->x,j);
    double a = yj-bj*xj+cj*xj*xj-dj*pow(xj,3);
    double b = bj-2*xj*cj+3*dj*xj*xj;
    double c = cj-3*xj*dj;
    double stepval = a*(z-xj)+0.5*b*(z*z-xj*xj)+1./3*c*(z*z*z-xj*xj*xj)+1./4*dj*(pow(z,4)-pow(xj,4));
    intval+=stepval;
    return intval;
}

double cdiff(cspline* s, double z){
    int j = binsearch(s->x,z);
    double bj = gsl_vector_get(s->b,j);
    double cj = gsl_vector_get(s->c,j);
    double dj = gsl_vector_get(s->d,j);
    double xj = gsl_vector_get(s->x,j);
    double b = bj-2*xj*cj+3*dj*xj*xj;
    double c = cj-3*xj*dj;
    return b+2*c*z+3*dj*z*z;
}


//FREEING THE ALLOCATED MEMORY
void cspline_free(cspline* s){
	gsl_vector_free(s->x);
	gsl_vector_free(s->y);
	gsl_vector_free(s->c);
	gsl_vector_free(s->b);
	gsl_vector_free(s->d);
	free(s);
	}
int main(){
    FILE* cxandy = fopen("out.cxy.txt","w");
    FILE* cout = fopen("out.cdata.txt","w");
	int n = 100;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* siny = gsl_vector_alloc(n);
	gsl_vector* cosy = gsl_vector_alloc(n);
	double xa[n];
	double ya[n];
	for(int i=0; i<n; i++){
	    double xi = ((double) i)/n*2*M_PI;
	    gsl_vector_set(x,i,xi);
	    gsl_vector_set(siny,i,sin(xi));
	    gsl_vector_set(cosy,i,cos(xi));
	    xa[i]=xi;
	    ya[i]=sin(xi);
	}
	cspline* s = cspline_make(x,siny);
    gsl_interp* c = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp_init(c,xa,ya,n);
	int count=0;
	while(count<100){
	    count++;
        double z = (double)rand()/(double) RAND_MAX*gsl_vector_get(x,n-2);
		double za=gsl_interp_eval(c,xa,ya,z,NULL);
		double zdiff = gsl_interp_eval_deriv(c,xa,ya,z,NULL);
		double zint = gsl_interp_eval_integ(c,xa,ya,gsl_vector_get(x,0),z,NULL);
		fprintf(cout,"%g %g %g %g %g %g %g\n",z,ceval(s,z), -cint(s,z)+1,cdiff(s,z),za,-zint+1,zdiff);
	}
    for(int i=0; i<n;i++) {
        fprintf(cxandy, "%g %g %g\n", gsl_vector_get(x, i),
                gsl_vector_get(siny, i), gsl_vector_get(cosy, i));
    }
	cspline_free(s);
	gsl_vector_free(x);
	gsl_vector_free(cosy);
	gsl_vector_free(siny);
	gsl_interp_free(c);
	return 0;
}
