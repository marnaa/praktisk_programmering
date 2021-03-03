#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(){
	double N = 1e8;
	double n1=N/3, n2=N/3, n3=N/3;
	double calN= n1+n2+n3;
	unsigned int seed1=0, seed2=0, seed3=0;
	double count1 = 0, count2 = 0, count3 = 0;
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for(double i=0; i<=n1; i+=1){
			double x = (double) rand_r(&seed1)/RAND_MAX;
			double y = (double) rand_r(&seed1)/RAND_MAX;
			if(x*x+y*y<=1){
				count1+=1;
                	}
                	else{
                	}
			}
		}
		#pragma omp section
		{
			for(double i=0; i<=n2; i+=1){
				double x = (double) rand_r(&seed2)/RAND_MAX;
				double y = (double) rand_r(&seed2)/RAND_MAX;
				if(x*x+y*y<=1){
					count2+=1;
                		}
                		else{
                		}
       			}
		}
		#pragma omp section
		{
			for(double i=0; i<=n3; i+=1){
			double x = (double) rand_r(&seed3)/RAND_MAX;
			double y = (double) rand_r(&seed3)/RAND_MAX;
			if(x*x+y*y<=1){
				count3+=1;
                	}
                	else{
                	}
       		}
		}


	}
	double final_count = count1 + count2 + count3;
	double est_pi = 4*final_count/(calN);
	printf("est pi = %g and math pi = %g\n", est_pi, M_PI);
	return 0;
}
