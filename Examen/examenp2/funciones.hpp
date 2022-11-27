// Declaracion de funciones para la tarea 11 
// Integracion Romberg

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
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double* d, int n);

/******************** Funciones para la integgracion de Romberg ************************/
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
double evalua_pol_cubico(double *px, double *py, double *coef, double val, int n);
double** evalua_datos_spline(double *px, double *py, double *coef, int n, int m);

#include "funciones.cpp"
