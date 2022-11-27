#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

double **createMatrix(int nr, int nc);
void freeMatrix(double **mat);
void inicializacion(double **A, int n);
void suma(double **A, double **B, int n, double **C);
void tiempoTranscurrido(struct timeval t1, struct timeval  t2);

int main(int argc,char *argv[])
{
    int      n = 1000;
    double **A, **B,  **C;
    struct timeval t1, t2;

    if(argc>1) n = atoi(argv[1]);

    srand(time(NULL));

    A = createMatrix(n, n);
    B = createMatrix(n, n);
    C = createMatrix(n, n);

    inicializacion(A, n);
    inicializacion(B, n);

    gettimeofday(&t1, NULL);
    printf("------- Calculo de la suma ------------\n");
    suma(A, B, n, C);
    gettimeofday(&t2, NULL);
    tiempoTranscurrido(t1, t2);

    freeMatrix(A);
    freeMatrix(B);
    freeMatrix(C);

    return(0);
}

double **createMatrix(int nr, int nc) {
    int i;
    double **mat;

    // Reservamos memoria
    mat = (double **) malloc( (nr)*sizeof(double *));
    if(mat==NULL) return(NULL);
    mat[0] = (double *) malloc( nr*nc*sizeof(double));

    if(mat[0]==NULL) return(NULL);

    for(i=1; i<nr; ++i)
        mat[i] = mat[i-1] + nc;

    return(mat);
}

// Libera la memoria del arreglo bidimensional
void freeMatrix(double **mat) {
    free(mat[0]);
    free(mat);
}

void inicializacion(double **A, int n)  {
    int   i, j;

    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            A[i][j] = ((rand()/(double)RAND_MAX)*2)-1;
}

void suma(double **A, double **B, int n, double **C) {
    int   i, j;

    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            C[i][j] = A[i][j] + B[i][j];
}

void tiempoTranscurrido(struct timeval t1, struct timeval  t2) {
    double elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;  // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;      // us to ms
    printf("Tiempo transcurrido %6.4f mseg. \n", elapsedTime);
}

