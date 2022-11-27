// Examen parcial 1 de metodos numericos
// La libreria "funciones.hpp" la implemente con las funciones de todo el curso y tareas

// Para compilar g++ -o pro main.cpp -fopenmp
// Para ejecutar ./pro

// Diego Aaron Moreno Galvan

/* La funcion del ejercicio 4 de la norma infinito ya la habia implementado en tareas
anteriores, asi que esta adjunta en la libreria "funciones.hpp" que hice */

#include "funciones.hpp"

void camino1(double **A, int n);
void encuentra_minymax(double *x, double *mn, double *mx, int n);

void camino2(double **V, int n);
double** vandermonde(int n);
double** inversaX(double **A, int n);
void inicializa_ei(double* ei, int fi, int n);
void copy_vec2mat(double **A, double *x, int ci, int n);
void calcula_numCondicion(double **A, double **X, double n);
double* resta_vectores(double *a,double *b, int n);
double* construyeb(int n);
double* construyec(double *b, int n);

int main(int narg, char **varg)  {
    int nr,nc;

    if(narg != 2) {printf("\nIntroduce el numero positivo n.\n"); return 0;}
    int n = atoi(varg[1]);

    double **V = vandermonde(n);
    camino1(V,n);
    camino2(V,n);

    freeMat(V);
	return 0;
}

double** vandermonde(int n) {
    double **A = crea_matriz(n,n);
    double aux;
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++) {
            aux = ((double)i + 1)/(double)n;
            aux = pow(aux, j);
            A[i][j] = aux;
        }
    return A;
}

/****************************** Funciones para el camino 1 *******************************************/
void camino1(double **A, int n) {
    double **M = crea_matriz(n,n);
    double **At = traspuesta(A,n);
    producto_matricial(A,At,M,n);
    freeMat(At);

    double tol = 0.000001;
    jacobiData data = jacobi(M,1000,tol,n);
    calcula_error_jacobi(M,data.V,data.lm,n);
    double s1, s2, mn, mx;
    encuentra_minymax(data.lm,&mn,&mx,n);
    s1 = sqrt(mx);
    s2 = sqrt(mn);
    double k2 = s2 / s1;

    printf("El numero de operaciones realizadas para encontrar los valores propios fue:\n");
    printf("%d\nY el valor de la toleracia es: %lf para toleracia inicial de %lf\n",data.k,data.e,tol);
    printf("\nEl numero de condicion k2 = s1 / s2 es: %lf\n", k2);

    freeMat(M);
    freeMat(data.V);
    free(data.lm);
}

void encuentra_minymax(double *x, double *mn, double *mx, int n)   {
    for(int i=0; i<n; i++)  {
        if(i == 0)
            *mn = *mx = x[i];
        else    {
            if(*mx < x[i]) *mx = x[i];
            if(*mn > x[i]) *mn = x[i];
        }
    }
}

/****************************** Funciones para el camino 2 *******************************************/
void camino2(double **V, int n) {
    double **X = inversaX(V,n);
    calcula_numCondicion(V,X,n);

    double *b = construyeb(n);
    double *c = construyec(b,n);
    double *x1 = resuelve_ecuacion(V,b,n);
    double *x2 = resuelve_ecuacion(V,c,n);
    double *x = resta_vectores(x1,x2,n);
    printf("La norma ||x1 - x2|| es: %lf\n",norma(x,n));

    free(x);
    free(x1); free(x2);
    free(c); free(b);
    freeMat(X);
}

double* resta_vectores(double *a,double *b, int n)  {
    double *resta = crea_vector(n);
    for(int i=0; i<n; i++)
        resta[i] = a[i] - b[i];
    return resta;
}

double* construyeb(int n)   {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        double aux = ((double)i + 1.0)/(double)n;
        x[i] = sin(aux);
    }
    return x;
}

double* construyec(double *b, int n)   {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        double aux = pow(-1,i+1) * (0.001);
        x[i] = b[i] + aux;
    }
    return x;
}

double** inversaX(double **A, int n)    {
    double **X = crea_matriz(n,n);
    double *ei = crea_vector(n);
    for(int i=0; i<n; i++)  {
        inicializa_ei(ei,i,n);
        double *xi = resuelve_ecuacion(A,ei,n);
        copy_vec2mat(X,xi,i,n);
        free(xi);
    }
    free(ei);
    return X;
}

void calcula_numCondicion(double **A, double **X, double n)   {
    double k;
    k = norma_inf(A,n);
    k *= norma_inf(X,n);
    printf("El numero de condicion kinf = ||A||inf||A^-1||inf es: %lf\n", k);
}

void inicializa_ei(double* ei, int fi, int n)   {
    inicializa_vec(ei,0,n);
    ei[fi] = 1;
}

void copy_vec2mat(double **A, double *x, int ci, int n) {
    for(int i=0; i<n; i++)
        A[i][ci] = x[i];
}
