// Declaracion de funciones para la tarea 12
// SOlucion de EDOS

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

/*** Funciones factorizacion LU y solucion de ecuaciones***/
int factoriza_LU(double** A, double** L, double** U, int n);
double* resuelve_triangularL(double** A, double* b, int n, int *ind);
double* resuelve_triangularU(double** A, double* b, int n, int *ind);
double* resuelve_diagonal(double* a, double* b, int n, int *ind);
double* resuelve_ecuacion(double **A, double *b, int n);
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double* d, int n);

/******************** Funciones para la solucion de EDOS ************************/
double Y(double x);
double f(double x, double y, double dy);
double* solve_EDO(double a, double b, int n);

#include "EDOS.cpp"
