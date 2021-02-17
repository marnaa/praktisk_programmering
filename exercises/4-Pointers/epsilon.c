#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
int equal(double a, double b, double tau, double epsilon);
int main(){
	printf("exercise epsilon \n");
	
	printf("1)\n");
	int i=1;
	while(i+1>i){
		i++;
	}
	printf("my max int = %i \n",i);
	printf("my INT_MAX = %i \n", INT_MAX);
	
	i=1;
	while(i-1<i){
		i--;  
	}
	printf("my min int = %i \n",i);
	printf("my INT_MIN = %i \n", INT_MIN);

	double x=1;
	while(1+x!=1){
	x/=2;
	}
	x*=2;
	printf("Calculated eps double while = %g\n",x);	

	// With do loops
	x=1;
	do{  
	x/=2;
	}
	while(1+x!=1);
	x*=2;
	printf("Calculated eps double do = %g\n",x); 


	// With for loops
	double e;
	for(e=1; 1+e!=1; e/=2){  
	}
	e*=2;
	printf("Calculated eps double for = %g\n", e); 

	printf("From float.h eps double= %g \n", DBL_EPSILON);
	printf("The rest of for and do's are ignored because of runtime\
	consideratons \n");


	float fx=1;
        while(1+fx!=1){
        fx/=2;
        }
        fx*=2;
        printf("Calculated eps float= %g\n",fx);
	/*	
	// Do loop
	fx=1;
	do{fx/=2;
	}while(1+x!=1);
	printf("From eps float for = %f \n", fx);

	
	//For loop 
	float fe;
	for(fe=1; 1+fe!=1; fe/=2){}
	fe*=2;	
	printf("Calculated eps float for = %f \n", fe);
	*/
        printf("From float.h eps float= %g \n", FLT_EPSILON);
	
	long double lx=1;
        while(1+lx!=1){
        lx/=2;
        }
        lx*=2;
        printf("Calculated eps long double= %Lg\n",lx);
       
	/*
	//Do loop
	lx=1;
	do{fx/=2;
        }while(1+x!=1);
        printf("From eps long double do = %Lg \n", lx);
	*/

	/*
	//For loop
	long double le;
	for(le=1; le+1!=1; le/=2){}
	le*=2;
	printf("From eps long double for = %Lg \n", le);
	*/

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
	double b=3.0;
	double eps=1.0;
	double tau=1.0;	
	int result = equal(a,b,tau, eps);
	printf("Result= %d \n", result);
	return 0;
}
