// Funciones hechas para el primer parcial de Metodos numericos
// funciones.hpp hecho para la declaracion de funciones que se utilizan

// Hecho por Diego Aaron Moreno Galvan

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>

/********** Estructuras ***********/
typedef struct datos_potencia_inversa {
    double mu; // Valor propio
    double *x; // Vector propio
    int k;     // Numero de operaciones realizadas
    double e;  // Norma ||w - px^|| < tolerancia

} eigpar_pot_inv;

typedef struct datos_jacobi {
    double *lm; // Valores propios de matriz
    double **V; // Vectores propios de matriz
    int k;     // Numero de operaciones realizadas
    double e;  // Entrada mas grande de matriz A^(k-1)

} jacobiData;

typedef struct indices_matriz {
    int i;
    int j;
    double mx; // Los indices del valor maximo

} indices;

/********** Funciones **********/
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

/*** Funciones generales para la creacion de estructuras ***/
void freeMat(double** mat);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
double *lee_arreglo(char *dir, int *n);
double **lee_matriz(char *dir, int *nr, int *nc);
void inicializa_mat(double** mat, double val, int nr, int nc);
void inicializa_vec(double* vec, double val, int n);
void inicializa_vector_random(double *vec, int n);
void inicializa_matriz_random(double **A, int n);
void guarda_matriz_txt(double **mat, int nr, int nc, char *cfile);

/*** Funciones de operaciones matriciales y vectoriales ***/
void suma(double* a, double* b, double* c, int n);
void multiplica(double* a, double* b, double* c, int n);
double producto_punto(double* a, double* b, int n);
void producto_MatVec(double** P, double* a, double* b, int n);
void producto_matricial(double** P, double** Q, double** R, int n);
void pivoteo_parcial(double** A, int fi, int fj, int nc);
void pivoteo_total(double** A, int fi, int fj, int ci, int cj, int nr, int nc);
double** matriz_pivoteo(int fi, int fj, int n);

/*** Funciones factorizacion LU y solucion de ecuaciones***/
int factoriza_LU(double** A, double** L, double** U, int n);
double* resuelve_triangularL(double** A, double* b, int n, int *ind);
double* resuelve_triangularU(double** A, double* b, int n, int *ind);
double* resuelve_diagonal(double* a, double* b, int n, int *ind);
double* resuelve_ecuacion(double **A, double *b, int n);
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double* d, int n);

/*** Funciones de errores ***/
void calcula_errorM(double **A, double *b, double *x, int n);
void calcula_errorLU(double **A, double **L, double **U, int n);
void calcula_error_jacobi(double **A, double **V, double *lm, int n);
double error_eigenpar(double *y, double lm, double *v, int n);
void calcula_error_fi(double* fiv, double* x, int n);
void calcula_error(double *a, double *b, double *x, int n);

/*** Funciones de normas ***/
double norma(double *vec, int n);
double norma_inf(double **A, int n);

/*** Funciones de los metodos: potencia, pot. inversa, jacobi, tridiagonal y valores propios ***/
jacobiData jacobi(double **A, int omax, double tol, int n);
double** traspuesta(double **A, int n);
double** rotacion_givens(double cos, double sen, int coori, int coorj, int n);
indices busca_mayor(double **A, int n);

void valores_propios(double **A, double d, int N, int M, double tol, int n);
eigpar_pot_inv potencia_inversa(double** A, double *x, double del, int omax, double tol, int n);
double** resta_matricial_diagonal(double** A, double del, int n);
double* discretizacion(double dx, int n);
double *calcula_d(double *x, double dx, int n);

double* potencia(double** A, int omax, double tol, int n);

void tiempoTranscurrido(struct timeval t1, struct timeval  t2);

#include "funciones.cpp"
