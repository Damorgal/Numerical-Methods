#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double ind, double* b, int n);
void sust_foward(double* d1, double* d3, double* b, double* x, int n);
void sust_baward(double* d1, double* d3, double* b, double* x, int n);
void convierte_tridiagonal(double* d1, double* d3, double ind, double* b, double* x, int n);
void revierte_tridiagonal(double* d1, double* d3, double ind, double* b, double* x, int n);
void calcula_error(double* x, double* s, int n);

double* crea_s(int n);
double* crea_b(int n);
double* crea_diagonal(double val, int n);
double* crea_vector(int n);

int main(int narg, char **varg)
{
    int n;
    if(narg < 2) {printf("\nIntroduce el numero n que es la dimension de la matriz.\n"); return 0;}
    n = atoi(varg[1]);

    // Creamos los 3 vectores 1D de la matriz tridiagonal -1 0 1
    double *dia1 = crea_diagonal(-1,n-2);
    double *dia2 = crea_diagonal(0,n-2);
    double *dia3 = crea_diagonal(1,n-2);
    // Termino independiente
    double ti = 1;
    // Generamos los vectores s y b
    double *s = crea_s(n);
    double *b = crea_b(n);

    // Resolvemos el sistema tridiagonal
    double *fi = resuelve_tridiagonal(dia1,dia2,dia3,ti,b,n);

    // Calculamos el error
    calcula_error(fi,s,n);

    // Imprimimos s y fi para visualizar lo que hizo
    for(int i=0; i<n; i++) {
        printf("s[%d]= %lf   fi[%d]= %lf\n",i,s[i],i,fi[i]);
    }

    free(fi);
    free(s); free(b);
    free(dia1);
    free(dia2);
    free(dia3);
    return 0;
}

// Este metodo esta basado en lo descrito en el archivo solTricero.txt
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double ind, double* b, int n) {
    double *x = crea_vector(n);
    // Convertimos la matriz a en una matriz tridiagonal
    convierte_tridiagonal(d1,d3,ind,b,x,n);

    // Como la diagonal principal es cero aplicamos este metodo
    // Este metodo esta basado en lo descrito en el archivo solTricero.txt
    sust_foward(d1,d3,b,x,n);
    sust_baward(d1,d3,b,x,n);

    // Regresamos la matriz A como estaba (sus diagonales y el vector b)
    revierte_tridiagonal(d1,d3,ind,b,x,n);
    return x;
}

void sust_foward(double* d1, double* d3, double* b, double* x, int n)   {
    x[2] = d3[0]*b[1];
    for(int i=3; i<n-2; i=i+2)
        x[i+1] = (b[i] - x[i-1]*d1[i-1]) / d3[i-1];
}

void sust_baward(double* d1, double* d3, double* b, double* x, int n)   {
    x[n-3] = d1[n-3]*b[n-2];
    for(int i=n-4; i>0; i=i-2)
        x[i-1] = (b[i] - x[i+1]*d3[i-1]) / d1[i-1];
}

void convierte_tridiagonal(double* d1, double* d3, double ind, double* b, double* x, int n) {
    x[0] = b[0] / ind;
    x[n-1] = b[n-1] / ind;
    b[1] = b[1] - d1[0]*x[0];
    d1[0] = 0;
    b[n-2] = b[n-2] - d3[n-3]*x[n-1];
    d3[n-2] = 0;
}

void revierte_tridiagonal(double* d1, double* d3, double ind, double* b, double* x, int n) {
    d1[0] = -1;
    b[1] = b[1] + d1[0]*x[0];
    d3[n-3] = 1;
    b[n-2] = b[n-2] + d3[n-3]*x[n-1];
}

void calcula_error(double* x, double* s, int n)    {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += (x[i]-s[i]) * (x[i]-s[i]);
    sum = sqrt(sum);
    printf("El error calculado ||X-S|| = %lf\n\n", sum);
}

double* crea_s(int n)   {
    double *x = crea_vector(n);
    double xi;
    for(int i=0; i<n; i++)  {
        xi = i/(double)(n-1);
        x[i] = exp(xi*xi);
    }
    return x;
}

double* crea_b(int n)   {
    double *x = crea_vector(n);
    x[0] = 1;
    x[n-1] = exp(1);
    double xi;
    for(int i=1; i<n-1; i++)  {
        xi = i/(double)(n-1);
        x[i] = (4.0 * xi) * (exp(xi*xi)) / (double)(n-1);

    }
    return x;
}

double* crea_diagonal(double val, int n)    {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)
        x[i] = val;
    return x;
}

double* crea_vector(int n)  {
    double* aux;
    aux = (double*)malloc(sizeof(double)*n);
    return aux;
}
