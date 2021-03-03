#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {double N; double count; unsigned int seed;} params; 

void* function(void* arg){
	params* p = (params*) arg;
	double N = (*p).N;
	double pse_count=0;
	unsigned int seed = (*p).seed;
	for(double i=0; i<=N; i+=1){
		double x = (double) rand_r(&seed)/RAND_MAX;
		double y = (double) rand_r(&seed)/RAND_MAX;
		if(x*x+y*y<=1){
			pse_count+=1;
		}
		else{
		}
	}
	p->count = pse_count;
	return NULL;
}

int main(){
	double N = 1e8;
	double n1=N/3, n2=N/3, n3=N/3;
	double calN= n1+n2+n3;
	unsigned int seed1=0, seed2=0, seed3=0;
	double count1 = 0, count2 = 0, count3 = 0;
	pthread_t t1, t2, t3;
	params p1 = {.N = n1, .count = count1, .seed = seed1};
	params p2 = {.N = n2, .count = count2, .seed = seed2};
	params p3 = {.N = n3, .count = count3, .seed = seed3};
	pthread_create(&t1, NULL, function, (void*)&p1);
	pthread_create(&t2, NULL, function, (void*)&p2);
	pthread_create(&t3, NULL, function, (void*)&p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	double final_count = p1.count + p2.count + p3.count;
	double est_pi = 4*final_count/(calN);
	printf("est pi = %g and math pi = %g\n", est_pi, M_PI);
	return 0;
}
