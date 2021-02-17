#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char** argv){
	if(argc<2) fprintf(stderr, "%s: there were no arguments \n", argv[0]);
	else{
		
		for (int i=1; i<argc; i++){
			double x = atof(argv[i]);
			printf("cos(%g)=%g and sin(%g)=%g\n",x,cos(x),x,sin(x));
		}
	}

	return 0;
}
