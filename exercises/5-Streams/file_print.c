#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(){
	int items;
	double x;
	FILE* infile = fopen("infile.txt","r");
	FILE* outfile = fopen("outfile.txt","w");
	while(items!=EOF){
		items=fscanf(infile, "%lg", &x);
		fprintf(outfile,"cos(%g)=%g and sin(%g)=%g\n",x,cos(x),x,sin(x));
	}
	fclose(infile);
	fclose(outfile);

} 
