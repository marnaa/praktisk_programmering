#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
// Part don't work on my pc with gcc 
int equal(double a, double b, double tau, double epsilon);
int main(){
/*
	int i=1;
	while(i+1>i){
		i++
	}
	printf("my max int = %i \n",i);
	printf("my INT_MAX = %i \n", INT_MAX);
*/	
	double x=1;
	while(1+x!=1){
	x/=2;
	}
	x*=2;
	printf("Calculated eps double= %g\n",x); 
	printf("From float.h eps double= %g \n", DBL_EPSILON);
	float fx=1;
        while(1+fx!=1){
        fx/=2;
        }
        fx*=2;
        printf("Calculated eps float= %f\n",fx);
        printf("From float.h eps float= %f \n", FLT_EPSILON);
	long double lx=1;
        while(1+lx!=1){
        lx/=2;
        }
        lx*=2;
        printf("Calculated eps long double= %Lg\n",lx);
        printf("From float.h eps long double= %Lg \n", LDBL_EPSILON);
	//Exercise 2)
	printf("2)\n");
	int max=INT_MAX /2;
	float sum_up_float=0;
	float sum_down_float=0;
	for (int i=1; i<=max; i++){
		sum_up_float+=1.0f/i;
		sum_down_float+=1.0f/(max-(i-1));
	}
	printf("max = %i \n",max); 
	printf("sum_up_float=%f\n ",sum_up_float);
	printf("sum_down_float= %f \n", sum_down_float);
        double sum_up_double=0;
        double sum_down_double=0;
        for (int i=1; i<=max; i++){
                sum_up_double+=1.0/i;
                sum_down_double+=1.0/(max-(i-1));
        }
        printf("sum_up_float=%f\n ",sum_up_double);
        printf("sum_down_float= %f \n", sum_down_double);
	// Exercise 3 
	printf("3)\n"); 
	double a=1.0;
	double b=2.0;
	double eps=4.0;
	double tau=2.0;	
	int result = equal(a,b,tau, eps);
	printf("Result= %d ", result);
	return 0;
}
