// Declaracion de funciones para la tarea 9 y 10 de interpolacion
// Interpolacion de funciones

// Hecho por Diego Aaron Moreno Galvan

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>

using namespace std;

/********** Estructuras ***********/
typedef struct datos_gradiente {
    double *x; // Solucion aproximada de Ax=b
    int k;     // Numero de operaciones realizadas
    double e;  // Valor del error (r*r/n)^(1/2)

} gradienteData;

/*** Funciones generales para la creacion de estructuras ***/
void freeMat(double** mat);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
double *lee_arreglo(char *dir, int *n);
double **lee_matriz(char *dir, int *nr, int *nc);
void guarda_matriz_txt(double **mat, int nr, int nc, char *cfile);

/*** Funciones de operaciones matriciales y vectoriales ***/
void suma(double* a, double* b, double* c, int n);
void multiplica(double* a, double* b, double* c, int n);
double producto_punto(double* a, double* b, int n);
void producto_MatVec(double** P, double* a, double* b, int n);
void producto_matricial(double** P, double** Q, double** R, int n);
void inicializa_mat(double** mat, double val, int nr, int nc);
void inicializa_vec(double* vec, double val, int n);
double norma(double *vec, int n);

/********* Funciones para la interpolacion de nube de puntos *************/
/* Polinomio interpolador */
double** vandermonde(double *x, int n);
double* coef_pol_interpolador(double *abscisas, double *ordenadas, int n);
double evalua_pol(double *coef, double x, int n);
double** evalua_datos(double *x, double *c, int n);
/* Polinomio de Lagrange */
double evalua_lagrange(double *x, double val, int ind, int n);
double evalua_pol_lagrange(double *x, double *y, double val, int n);
double** evalua_datos_lagrage(double *x, double *y, int n);
/* Polinomio de diferencias */
double evalua_diferencias(double *x, double val, int ind);
double** matriz_diferencias(double *x, double *y, int n);
double evalua_pol_diferencias(double *x, double *y, double val, int n);
double** evalua_datos_diferencias(double *x, double *y, int n);
/* Nubes de puntos */
double evalua_fi(double *fis, double val, double *dis, int n);
double** evalua_datos_nubes(double *px, double *fis, double *x, int n, int m);
double* minimiza_diferencias(double *px, double *py, double *x, double lm, double h, int n, int m);
double Nk(double *x, double val, double h, int k, int n);

#include "interpolacion.cpp"
