#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define  RND ((double) rand()/ RAND_MAX)

double corput(int n, int base){
    double q=0, bk=(double ) 1 / base ;
    while ( n>0){
        q += ( n % base)* bk ;
        n /= base;
        bk /= base;
    }
    return q;
}
void halton(int n, int dim, double x[],double xerr[]){
    int base[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    int base2[]={5,7,11,13,17,19,23,29,31,37,41,43,47,53,59};
    int maxdim = sizeof(base)/sizeof(int);
    assert(dim<=maxdim);
    for(int i=0; i<dim; i++){
        x[i]=corput(n,base[i]);
        xerr[i] = corput(n,base2[i]);
    }

}

void quasiMC_int(int dim, int N, double a[],double b[],
                double f(int dim,double x[]), double* result, double* error){
    double mean_f = 0, mean_ff = 0, V = 1, x[dim], xerr[dim];
    for(int i = 0; i<dim; i++){
        V*=b[i]-a[i];
    }
    for(int i=0; i<N; i++){
        halton(i,dim,x,xerr);
        for(int j=0; j<dim; j++){
            x[j]=a[j]+x[j]*(b[j]-a[j]);
            xerr[j]=a[j]+xerr[j]*(b[j]-a[j]);
        }
        double fx = f(dim,x);
        double fxerr = f(dim,xerr);
        if(isinf(fxerr)==0 && isinf(fx)==0) {
            mean_f += fx;
            mean_ff += fxerr;
        }
        //printf(" means %g %g\n",mean_f,mean_ff);
    }
    mean_f*=1./N;
    mean_ff*=1./N;
    *error = fabs(mean_f-mean_ff)*V;
    *result = mean_f*V;
}

void randMC_int(int dim, int N, double a[],double b[],
                double f(int dim,double x[]), double* result, double* error){
    double mean_f = 0, mean_ff = 0, V = 1, x[dim];
    for(int i = 0; i<dim; i++){
        V*=b[i]-a[i];
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<dim; j++){
            x[j]=a[j]+RND*(b[j]-a[j]);
        }
        mean_f+=f(dim,x);
        mean_ff+=f(dim,x)*f(dim,x);
    }
    mean_f*=1./N;
    mean_ff*=1./N;
    *error = sqrt(mean_ff-mean_f*mean_f)*V/sqrt(N);
    *result = mean_f*V;
}

double strata(int dim, int N, double a[],double b[],
             double f(int dim,double x[]), double* error,
             double acc, double eps,int n_reuse, double mean_reuse) {
    double V = 1, mean = 0;
    double maxvar;
    int n = 0.2*N;
    int idiv = 0;
    int n_left[dim], n_right[dim];
    double x[dim], mean_left[dim], mean_right[dim];
    for (int i = 0; i < dim; i++) {
        V *= b[i] - a[i];
        mean_left[i] = 0;
        n_left[i] = 0;
        mean_right[i] = 0;
        n_right[i] = 0;
    }
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < dim; i++) {
            x[i] = a[i] + RND * (b[i] - a[i]);
        }
        double fx = f(dim, x);
        if(isinf(fx)==0){
            mean += fx;
            for (int i = 0; i < dim; i++) {
                if (x[i] < (a[i] + b[i]) / 2) {
                    n_left[i]++;
                    mean_left[i] += fx;
                }
                else {
                n_right[i]++;
                mean_right[i] += fx;
                }
            }
        }
    }
    mean/=n;

    for(int i = 0; i<dim; i++){
        mean_left[i]/=n_left[i];
        mean_right[i]/=n_right[i];
    }
    for(int i = 0; i<dim; i++){
        double var=fabs(mean_right[i]-mean_left[i]);
        if(var>maxvar){
            maxvar = var;
            idiv = i;
        }
    }
    double integ = (mean*n+mean_reuse*n_reuse)/(n+n_reuse)*V;
    double err = fabs(mean_reuse-mean)*V;
    double tol = acc+fabs(integ)*eps;
    if(err<tol){
        *error += err;
        return integ;
    }
    double a2[dim], b2[dim];
    for(int i=0; i<dim; i++){
        a2[i] = a[i];
        b2[i] = b[i];
    }
    a2[idiv] = (a[idiv]+b[idiv])/2;
    b2[idiv] = (a[idiv]+b[idiv])/2;
    double integ_left =
            strata(dim,N,a,b2,f,error,acc/sqrt(2),eps,n_left[idiv],mean_left[idiv]);
    double integ_right = strata(dim,N,a2,b,f,error,acc/sqrt(2),eps,n_right[idiv],mean_right[idiv]);
    return integ_left+integ_right;
}


double norm_dist_sig1(int dim,double x[]){
    double sigma = 1.;
    double val = 1./(sigma* sqrt(2*M_PI))*exp(-0.5*pow((x[0]/sigma),2));
    if(x[0]<1 && x[0]>-1){
        return val;
    }
    else{
        return 0.;
    }
}double norm_dist_sig2(int dim,double x[]){
    double sigma = 1.;
    double val = 1./(sigma* sqrt(2*M_PI))*exp(-0.5*pow((x[0]/sigma),2));
    if(x[0]<2 && x[0]>-2){
        return val;
    }
    else{
        return 0.;
    }
}
double gam(int dim, double x[]){
    for(int i = 0; i<dim; i++){
        if(x[i]>M_PI ||x[i]<0 ){
            return 0;
        }
    }
    double val = 1./((pow(M_PI,3))*(1-cos(x[0])*cos(x[1])*cos(x[2])));
    return val;
}double dim_special(int dim, double x[]){
    if(x[0]*x[0]+x[1]*x[1]<1.){
        double val = 1.;
        return val;
    }
    else return 0;
}

int main(){
    FILE* compare = fopen("out.compare.txt","w");
    double a[] = {-10};
    double b[] = {10};
    double a_gam[] = {0,0,0};
    double b_gam[] = {M_PI,M_PI,M_PI};
    double a_dim[] = {-1.,-1.};
    double b_dim[] = {1.,1.};
    double norm_sig1, err_sig1;
    double norm_sig2, err_sig2;
    double res_gam, err_gam;
    double res_gam_str, err_gam_str;
    double res_gam_quasi, err_gam_quasi;
    int dim = 1, N = 10000,dimGam=3, dimDim=2;
    for(int i = 0;i<50;i++){
        int N_gam = 10000+10000*i;
        randMC_int(dimDim,N_gam,a_dim,b_dim,dim_special,&res_gam,&err_gam);
        quasiMC_int(dimDim,N_gam,a_dim,b_dim,dim_special,&res_gam_quasi,&err_gam_quasi);
        fprintf(compare,"%i %g %g %g %g %g\n",N_gam, res_gam,res_gam_quasi,err_gam,err_gam_quasi,M_PI);
    }

    randMC_int(dim,N,a,b,norm_dist_sig1,&norm_sig1, &err_sig1);
    //randMC_int(dim,N,a,b,norm_dist_sig2,&norm_sig2, &err_sig2);
    quasiMC_int(dim,N,a,b,norm_dist_sig2,&norm_sig2,&err_sig2);
    printf("int of normal dist 1 sigma = %g error = %g N = %i\n",norm_sig1, err_sig1, N);
    printf("int of normal dist 2 sigma [quasi]= %g error = %g N = %i\n",norm_sig2, err_sig2, N);
    printf("Val of ex b) integral (exact val = 1.3932): \n");
    for(int i = 0;i<5;i++){
        int N_gam = 10000+10000*i;
        randMC_int(dimGam,N_gam,a_gam,b_gam,gam,&res_gam,&err_gam);
        quasiMC_int(dimGam,N_gam,a_gam,b_gam,gam,&res_gam_quasi,&err_gam_quasi);
        res_gam_str = strata(dimGam,N_gam,a_gam,b_gam,gam,&err_gam_str,0.005,0.01,0,0);
        printf("rand = %g (err = %g), quasi = %g (err = %g), strata = %g (err = %g) at N = %i \n", res_gam, err_gam,res_gam_quasi,err_gam_quasi,res_gam_str, err_gam_str,N_gam);
    }
    fclose(compare);
    return 0;
}
