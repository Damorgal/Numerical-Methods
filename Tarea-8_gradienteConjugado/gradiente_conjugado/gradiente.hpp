// Declaracion de funciones para la tarea 8
// Metodo del gradiente conjugado

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

/*** Funciones de operaciones matriciales y vectoriales ***/
void suma(double* a, double* b, double* c, int n);
void multiplica(double* a, double* b, double* c, int n);
double producto_punto(double* a, double* b, int n);
void producto_MatVec(double** P, double* a, double* b, int n);
void producto_matricial(double** P, double** Q, double** R, int n);
void inicializa_mat(double** mat, double val, int nr, int nc);
void inicializa_vec(double* vec, double val, int n);
double norma(double *vec, int n);

/********* Funciones para metodo del gradiente conjugado *************/
gradienteData gradiente_conjugado(double **A, double *b, double tol, int n);

#include "gradiente.cpp"
