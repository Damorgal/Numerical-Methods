#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
double** traspuesta(double **A, int n, int m);
void producto_MatVec(double** P, double* a, double* b, int nr, int nc);
void producto_matricial(double** P, double** Q, double** R, int nr, int nc);
int factoriza_LU(double** A, double** L, double** U, int n);
double* resuelve_triangularL(double** A, double* b, int n, int *ind);
double* resuelve_triangularU(double** A, double* b, int n, int *ind);
double* resuelve_ecuacion(double **A, double *b, int n);

/*** Funciones para el programa ***/
double* crea_valores_xi(int m);
double* crea_valores_yi(double* x, int n);
double* minimos_cuadrados(double* x, double* y, int m);
double calcula_error_mincuadrados(double* f, double* y, int m);

#include "cuadrados.cpp"
