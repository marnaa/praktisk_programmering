#include <math.h>
#include <stdio.h>
double ex(double x){
if(x<0)return 1/ex(-x);
if(x>1./8)return pow(ex(x/2),2);
return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}
int main(){
	int N=200;
	for(int i=-N+1; i<N; i++){
	double x=5* (double) i/ ((double) N);
	printf("%g %g %g\n",x, ex(x), exp(x));
	}
	return 0;
}
