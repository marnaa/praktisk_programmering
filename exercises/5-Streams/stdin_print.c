#include <stdio.h>
#include <math.h>
int main(){
	double x;
	int items;
	while(items!=EOF){
	items=fscanf(stdin,"%lg",&x);
	printf("cos(x)=%g sin(x)=%g\n",cos(x),sin(x));
	}
	return 0;
}

