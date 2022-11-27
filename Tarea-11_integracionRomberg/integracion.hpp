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

/******************** Funciones para la integgracion de Romberg ************************/
double f1(double x);
double f2(double x);
double trapecio_recursivo(double (*f)(double), double val, double a, double b, int ind);
double integracion(double (*f)(double), double a, double b, int n);

#include "integracion.cpp"
