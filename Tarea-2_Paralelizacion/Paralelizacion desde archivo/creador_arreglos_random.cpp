#include <stdio.h>
#include <stdlib.h>


int main()
{
	char dir[] = "/home/damorgal/Documents/Clase Met_Num/Paralelizacion/arreglos2.txt";
	
	FILE *fp;
    int n=100000000;

	fp = fopen(dir, "wb");
	if (fp==NULL) {fputs ("File error",stderr); exit (1);}
	fprintf(fp, "%d\n", n);
	for(int i=0; i<n; i++)
		fprintf(fp, "%f ", ((rand()/(double)RAND_MAX)*2)-1);
	fprintf(fp, "\n");
	for(int i=0; i<n; i++)
		fprintf(fp, "%f ", ((rand()/(double)RAND_MAX)*2)-1);
	fclose(fp);

	return 0;
}
