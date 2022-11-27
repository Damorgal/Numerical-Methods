#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

void suma(double* a, double* b, double* c, int n);

int lee_long(const char *dir);
void lee_arreglos(const char *dir, double *a, double *b, int n);
void tiempoTranscurrido(struct timeval t1, struct timeval  t2);

int main()  {
    double *a, *b, *c;
    struct timeval t1, t2;
    int n;
    // Cambiar la direccion a donde se encuentre el archivo
    char dir[] = "/home/damorgal/Documents/Clase Met_Num/Paralelizacion/arreglos2.txt";

    n = lee_long(dir);
    a = (double*)malloc(sizeof(double)*n);
	b = (double*)malloc(sizeof(double)*n);
	c = (double*)malloc(sizeof(double)*n);
	lee_arreglos(dir,a,b,n);

    gettimeofday(&t1, NULL);
    suma(a, b, c, n);
    gettimeofday(&t2, NULL);

    tiempoTranscurrido(t1, t2);

	free(c);
    free(b);
    free(a);
	return 0;
}

// Funcion que suma los vectores a y b y los guarda en c (le debes pasar la longitud n de los vectores)
void suma(double *a, double *b, double *c, int n) {
	#pragma omp parallel for
	for (int i=0; i<n; i++)
        c[i] = a[i] + b[i];
}

// Funcion que lee la longitud de los arreglos desde un archivo
int lee_long(const char *dir)   {
    FILE *fp;
    int n;

	fp = fopen(dir, "rb");
	if (fp==NULL) {fputs ("File error",stderr); exit (1);}
	fscanf(fp, "%d", &n);
	fclose(fp);
	return n;
}

// Funcion que lee de archivo los arreglos a y b dada la longitud n antes leida
void lee_arreglos(const char *dir, double *a, double *b, int n) {
    FILE *fp;

	fp = fopen(dir, "rb");
	if (fp==NULL) {fputs ("File error",stderr); exit (1);}
	int aux; fscanf(fp,"%d",&aux);  // Auxiliar para volver a leer n

	// Leer arreglos
	for(int i=0; i<n; i++)
        fscanf(fp, "%lf", &a[i]);
    for(int i=0; i<n; i++)
        fscanf(fp, "%lf", &b[i]);

	fclose(fp);
}

void tiempoTranscurrido(struct timeval t1, struct timeval  t2) {
    double elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;  // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;      // us to ms
    printf("Tiempo transcurrido %6.4f mseg. \n", elapsedTime);
}
