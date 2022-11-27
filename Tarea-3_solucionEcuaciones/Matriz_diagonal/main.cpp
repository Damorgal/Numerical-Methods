#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* resuelve_diagonal(double* a, double* b, int n, int *ind);
void calcula_error(double *a, double *b, double *x, int n);

double *lee_arreglo(char *dir, int *n);
double* crea_vector(int n);

int main(int narg, char **varg)  {
    int n,m;
    char *dirmat = varg[1];
    char *dirvec = varg[2];

    if(narg != 3) {printf("\nIntroduce la direccion de los DOS archivos binarios (matriz y sol. b)\n"); return 0;}
    // Leemos los archivos
    double *a = lee_arreglo(dirmat, &n);
    double *b = lee_arreglo(dirvec, &m);
    if(n != m) {printf("Error, la dimension del vector b y la matriz no coinciden.\n"); free(a); free(b); return 0;}

	// Resolvemos el sistema
	// Indicador nos dira si tiene infinitas soluciones
    int ind=0;
	double *x = resuelve_diagonal(a,b,n,&ind);

    // Salida del programa
    if(x == NULL) printf("La matriz no tiene solucion.\n");
    else if(ind == 1) printf("La matriz tiene infinitas soluciones.\n");
    else    {
        printf("\nEl tamagno de la matriz es: %dx%d\n",n,n);
        printf("La matriz SI tiene solucion\n");
        calcula_error(a,b,x,n);
    }
	free(x);
    free(b);
    free(a);
	return 0;
}

// Resuelve una matriz diagonal (el vector diagonal)
double* resuelve_diagonal(double* a, double* b, int n, int *ind)   {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        if(a[i] != 0)
            x[i] = b[i]/a[i];
        else if(b[i] == 0)
            *ind=1;
        else {free(x); return NULL;}
    }
    return x;
}

// Calcula el error de que un a por un vector x de la solucion del vector b
void calcula_error(double *a, double *b, double *x, int n)  {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += std::pow(a[i]*x[i] - b[i], 2);
    sum = std::sqrt(sum);
    printf("Y el error calculado de la solucion es: %lf\n", sum);
}

// Funcion que lee de archivo un vector
double *lee_arreglo(char *dir, int *n) {
    FILE *fp;
    double *vec;

	fp = fopen(dir, "rb");
	if (fp==NULL) {printf("El archivo %s no puedo ser leido.\n", dir); exit (1);}

	// Pedir memoria
	fread(n, sizeof(int), 1, fp);
    vec = crea_vector(*n);

	// Leer arreglo
	fread(vec, sizeof(double), *n, fp);

	fclose(fp);
	return(vec);
}

double* crea_vector(int n)  {
    double* aux;
    aux = (double*)malloc(sizeof(double)*n);
    return aux;
}
