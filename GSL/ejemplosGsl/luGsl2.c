
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
 


// Lee el valor de los parametros y los asigna a las variables:
// cfile1   - Nombre de la matriz
// cfile1   - Nombre del vector
void readParams(int argc, char **argv, char *cfile1, char *cfile2) {
    if(argc>1) strcpy(cfile1, argv[1]);
    if(argc>2) strcpy(cfile2, argv[2]);
}

// Lectura del vector en el archivo cfile.
// Devuelve NULL si no se pudo abrir el archivo.
double *readVector(char *cfile, int *nr) {
    double   *vec;
    FILE       *f1 = fopen(cfile, "rb");
    
    if(!f1) return(NULL);
    fread(nr, sizeof(int), 1, f1);
    vec = (double *) malloc( (*nr)*sizeof(double));
    if(vec==NULL) return(NULL);
    fread(vec, sizeof(double), *nr, f1);
    fclose(f1);
    return(vec);
}


// Reserva memoria para una matriz de tama~no nr x nc
double **createMatrix(int nr, int nc) {
    int      i;
    double **mat;

    // Reservamos memoria
    mat    = (double **) malloc( (nr)*sizeof(double *));
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

// Lectura de la matriz en el archivo cfile.
// Devuelve NULL si no se pudo abrir el archivo.
double **readMatrix(char *cfile, int *nr, int *nc) {
    double  **mat;
    FILE     *f1 = fopen(cfile, "rb");
    
    if(!f1) return(NULL);
    // Lectura del tama~no de la matriz
    fread(nr, sizeof(int), 1, f1);
    fread(nc, sizeof(int), 1, f1);
    // Reservamos memoria
    mat    = createMatrix(*nr, *nc);
    // Lectura de los datos
    fread(mat[0], sizeof(double), (*nr)*(*nc), f1);
    fclose(f1);
    return(mat);
}

// Imprime en la consola las entradas del vector
void printVector(double *vec, int dim) {
    int i; 
    for(i=0; i<dim; ++i) 
	    printf("% 8.4f   ", vec[i]);
	printf("\n");
}


// Imprime en la consola las entradas de la matriz
void printMatrix(double **mat, int nr, int nc) {
    int i, j; 
    for(i=0; i<nr; ++i) {
        for(j=0; j<nc; ++j) 
            printf("% 6.2f   ", mat[i][j]);
        printf("\n");
    }
}

void initGslVector(gsl_vector      *vec, double *v, int dim) {
    gsl_block *block=vec->block;

    vec->data   = v;
    vec->size   = dim;
    vec->owner  = 0;
    vec->stride = 1;
    if(block!=NULL) {
        block->data = v;
        block->size = dim;
    }
}



void solve(int nr, int nc, double **mA, double *vb, double *vx) {
    int              i, signum, nitems=nr*nc;
    gsl_permutation *perm;
    gsl_matrix      *matLU;
    gsl_vector      *vecb, *vecx;
    double          *p, *pA;

    // Reservamos memoria para los apuntadores
    vecb = (gsl_vector *) malloc(sizeof(gsl_vector));
    vecx = (gsl_vector *) malloc(sizeof(gsl_vector));

    // Hacemos que los vectores hagan referencia a las arreglos vb y vx.
    // En este caso no se duplica la informacion.
    initGslVector(vecb, vb, nr);
    initGslVector(vecx, vx, nr);

    // Creamos una matriz y copiamos la informacion del arreglo mA.
    // Se duplica la informacion porque al calcular la factorizacion
    // se destruye la matriz.
    matLU = gsl_matrix_alloc(nr, nc);
    
    p  = matLU->data;
    pA = mA[0];
    for(i=0; i<nitems; ++i, ++p, ++pA) *p = *pA;
 
    // Calcula la factorizacion LU y resuelve el sistema
    perm = gsl_permutation_alloc(nc);
    gsl_linalg_LU_decomp(matLU, perm, &signum);
    gsl_linalg_LU_solve(matLU, perm, vecb, vecx);

    gsl_permutation_free(perm);
    gsl_matrix_free(matLU); 
    free(vecb);
    free(vecx);
}


int main(int argc, char **argv)      {
    char    cfile1[50], cfile2[50];

    // Lectura de los parametros
    readParams(argc, argv, cfile1, cfile2);

    // Lectura e impresion de la matriz
    int     nr, nc, dim;
    double **mat= readMatrix(cfile1, &nr, &nc);
    double *vx, *vb = readVector(cfile2, &dim);
    if(vb==NULL) {
        printf("El archivo %s no puedo ser leido.\n", cfile2);
        return(0);
    }
    if(mat==NULL) {
        printf("El archivo %s no puedo ser leido.\n", cfile1);
        return(0);
    }

    printf("\nVector de entrada:\n" );
    printVector(vb, dim);
    printf("\nMatrix de entrada:\n" );
    printMatrix(mat, nr , nc);

    vx = (double *) malloc( dim*sizeof(double));

    solve(nr, nc, mat, vb, vx);

    printf("\nSolucion:\n" );
    printVector(vx, dim);

    free(vb);
    free(vx);
    freeMatrix(mat);
    return(0);
}
// Para compilar
// gcc -o luGsl2 luGsl2.c -lgsl -lgslcblas -lm 

// Para ejecutar
// ./luGsl2 matA1.bin vecb1.bin





