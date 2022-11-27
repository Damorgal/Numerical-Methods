#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define k 0.5
#define a 0.0
#define b 5.0

double *calcula_d(double *x, double dx, int n);
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double* d, int n);
void calcula_error(double* fiv, double* x, int n);
double* discretizacion(double dx, int n);

void guarda_matriz_txt(double **mat, int nr, int nc, char *cfile);
double* crea_diagonal(double val, double ind, int n);
double* crea_vector(int n);
void freeMat(double** mat);
double** crea_matriz(int nr, int nc);

// fi(x) = 2 + xe^x
double fi(double val)   {
    double x = 2.0 + val*exp(val);
    return x;
}

// q(x) = -k fi''(x) y fi''(x) = 2e^x + xe^x pues fi(x) = 2 + xe^x
double q(double val)    {
    double x = -0.5 * (2.0*exp(val) + val*exp(val));
    return x;
}

int main(int narg, char **varg)
{
    int n;
    if(narg < 2) {printf("\nIntroduce el numero n (el numero de particiones del intervalo [a,b]).\n"); return 0;}
    n = atoi(varg[1]);

    // Creamos los 3 vectores 1D de la matriz tridiagonal -1 0 1
    double *dia1 = crea_diagonal(1,0,n+1);
    double *dia2 = crea_diagonal(-2,1,n+1);
    double *dia3 = crea_diagonal(1,0,n+1);

    // Generamos la discretizacion
    double dx = (b - a) / n;
    double *x = discretizacion(dx,n+1);

    double *d = calcula_d(x,dx,n+1);
    // Resolvemos el sistema tridiagonal usando la formula que debimos usar
    double *fiv = resuelve_tridiagonal(dia1,dia2,dia3,d,n+1);

    // Calculamos el error
    calcula_error(fiv,x,n+1);

    // Guarda matriz en un archivo
    double **mat = crea_matriz(3,n+1);
    mat[0] = x;
    for(int i=0; i<n+1; i++)
        mat[1][i] = fi(x[i]);
    mat[2] = fiv;
    char s[] = "DatosEcuPoisson100.txt";
    guarda_matriz_txt(mat,3,n+1,s);

    freeMat(mat);
    free(d);
    free(dia1);
    free(dia2);
    free(dia3);
    return 0;
}

// Lo resolvemos con el caso particular de la gaussiana visto en clase
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double* d, int n) {
    double *bi = crea_vector(n);
    double *ci = crea_vector(n);
    double *di = crea_vector(n);
    // Procedemos con el algoritmo para crear las entradas modificadas
    for(int i=0; i<n; i++)   {
        if(i == 0)  {
            bi[i] = d2[i];
            ci[i] = d3[i];
            di[i] = d[i];
        }
        else    {
            bi[i] = bi[i-1]*d2[i] - d1[i]*ci[i-1];
            ci[i] = bi[i-1]*d3[i];
            di[i] = bi[i-1]*d[i] - d1[i]*di[i-1];
        }
    }
    // Creamos las n variables fi de la solucion
    double *x = crea_vector(n);
    x[n-1] = d[n-1];
    for(int i=n-2; i>=0; i--)
        x[i] = (di[i] - ci[i]*x[i+1]) / bi[i];

    free(bi);
    free(ci);
    free(di);
    return x;
}

void calcula_error(double* fiv, double* x, int n)    {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += (fi(x[i]) - fiv[i]) * (fi(x[i]) - fiv[i]);
    sum = sqrt(sum);
    printf("El error calculado de la solucion analitica con la de ecuaciones es:\n||FI(x) - fi(xi)|| = %lf\n", sum);
}

double *calcula_d(double *x, double dx, int n)     {
    double *d = crea_vector(n);
    for(int i=0; i<n; i++)  {
        if(i == 0 || i == n-1) d[i] = fi(x[i]);
        else d[i] = -q(x[i]) * dx * dx / k;
    }
    return d;
}

double* discretizacion(double dx, int n)    {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)
        x[i] = a + i*dx;
    return x;
}

void guarda_matriz_txt(double **mat, int nr, int nc, char *cfile) {
    FILE   *fp=fopen(cfile, "wt");
    int     i, j;

    if(!fp) {
        printf("No se puede abrir el archivo\n");
        exit(0);
    }
    printf("Generando el archivo %s ...\n", cfile);

    for(i=0; i<nc; ++i) {
        for(j=0; j<nr; ++j)
            fprintf(fp, "%lf    ", mat[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

double* crea_diagonal(double val, double ind, int n)    {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        if(i == 0 || i == n-1)
            x[i] = ind;
        else x[i] = val;
    }
    return x;
}

double* crea_vector(int n)  {
    double* aux;
    aux = (double*)malloc(sizeof(double)*n);
    return aux;
}

double** crea_matriz(int nr, int nc) {
    int i;
    double** mat;

    mat = (double**)malloc(nr * sizeof(double*));
    if(mat==NULL) return(NULL);
    mat[0] = (double*)malloc(nr * nc * sizeof(double));
    if(mat[0]==NULL) return(NULL);

    for(i=1; i<nr; i++)
        mat[i] = mat[i-1] + nc;

    return mat;
}

// Libera memoria de matriz
void freeMat(double** mat)   {
    free(mat[0]);
    free(mat);
}
