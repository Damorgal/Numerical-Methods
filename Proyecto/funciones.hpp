// Funciones hechas para el proyecto final

// Hecho por Diego Aaron Moreno Galvan

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>

using namespace std;

/*** Funciones generales para la creacion de estructuras ***/
void freeMat(double** mat);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
double *lee_arreglo(char *dir, int *n);
double **lee_matriz(char *dir, int *nr, int *nc);
void guarda_matriz_txt(double **mat, int nr, int nc, char *cfile);
void inicializa_mat(double** mat, double val, int nr, int nc);
void inicializa_vec(double* vec, double val, int n);

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

/******************** Funciones para la integracion de Romberg ************************/
double f1(double x);
double f(double x, int n);
double derf(double x, double h);
double trapecio_recursivo(double (*f)(double), double val, double a, double b, int ind);
double integracion(double (*f)(double), double a, double b, int n);
double newton_rapson(double x, double h, double tol, int M, int n);
double determina_h(double x, double tol);
void haz_spline(double x);
/* Splines cubicos */
double* coef_spline_cubico(double *px, double *py, int n);
double evalua_pol_cubico(double *px, double *py, double *coef, double val, int n, int ind);
double** evalua_datos_spline(double *px, double *py, double *coef, int n, int m);

/********************** Funciones para proyecto *****************************************/
typedef struct punto_tridimensional	{
	double x;
	double y; 
	double z;
} punto;

typedef struct espacio_de_vigilancia	{
	double ancho;
	double largo;
	punto poscam1;
	punto poscam2;
	punto apcam1;
	punto apcam2;
} place;

punto trayectoria_1(double t);

#include "funciones.cpp"
