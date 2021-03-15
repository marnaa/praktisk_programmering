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
double ceval(cspline* s, double z){
        int i = binsearch(s->x, z);
        double h = z - gsl_vector_get(s->x,i);
	double bi = gsl_vector_get(s->b,i);
	double ci = gsl_vector_get(s->c,i);
	double di = gsl_vector_get(s->d,i);
	double val = gsl_vector_get(s->y,i)+h*(bi+h*(ci+h*di));
        return val;
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
	int n = 5;
	gsl_vector* x = gsl_vector_alloc(5);
	gsl_vector* y = gsl_vector_alloc(5);
	gsl_vector* yx = gsl_vector_alloc(5);
	gsl_vector* yxx = gsl_vector_alloc(5);
	gsl_vector* y3x = gsl_vector_alloc(5);
	double xa[n];
	double yk[n];
	double ya[n];
	double yaa[n];
	double y3a[n];
	for(int i=0; i<n; i++){
	gsl_vector_set(x,i,(double)i);
	gsl_vector_set(y,i,1.0);
	gsl_vector_set(yx,i,(double)i);
	gsl_vector_set(yxx,i,(double)(i*i));
	gsl_vector_set(y3x,i,(double)(i*i*i));
	xa[i]=(double) i;
	yk[i]=1.;
	ya[i]=(double) i;
	yaa[i]=(double)(i*i);
	y3a[i]=(double) (i*i*i);
	}
	cspline* s = cspline_make(x,y);
	cspline* sx = cspline_make(x,yx);
	cspline* sxx = cspline_make(x,yxx);
	cspline* s3x = cspline_make(x,y3x);
	gsl_interp* c = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp* ca = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp* caa = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp* c3a = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp_init(c,xa,yk,n);
	gsl_interp_init(ca,xa,ya,n);
	gsl_interp_init(caa,xa,yaa,n);
	gsl_interp_init(c3a,xa,y3a,n);
	for(int i=0; i<n-1; i++){
		double z = ((double)i+0.5) ;
		printf("%g\n",z);
		double zk=gsl_interp_eval(c,xa,yk,z,NULL);
		double za=gsl_interp_eval(ca,xa,ya,z,NULL);
		double zaa=gsl_interp_eval(caa,xa,yaa,z,NULL);
		double z3a=gsl_interp_eval(c3a,xa,y3a,z,NULL);
		printf("%g %g %g %g %g %g %g %g\n", 
		ceval(s,z),ceval(sx,z),ceval(sxx,z),ceval(s3x,z),
		zk, za, zaa, z3a);
	}
	printf("Spline values:\n");
	for(int i=0; i<n-1; i++){
		double a=gsl_vector_get(sxx->y,i);
		double b=gsl_vector_get(sxx->b,i);
		double c=gsl_vector_get(sxx->c,i);
		double d=gsl_vector_get(sxx->d,i);
		printf("%g %g %g %g\n",a,b,c,d);
	}
	cspline_free(s);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(yx);
	gsl_vector_free(yxx);
	gsl_vector_free(y3x);
	gsl_interp_free(c);
	gsl_interp_free(ca);
	gsl_interp_free(caa);
	gsl_interp_free(c3a);
	return 0;
}
