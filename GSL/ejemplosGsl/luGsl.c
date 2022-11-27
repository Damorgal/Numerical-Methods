#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
 
// Lee el valor de los parametros y los asigna a las variables:
// cmat   - Nombre del archivo de la matriz
// cvec   - Nombre del archivo del vector
void readParams(int argc, char **argv, char *cmat, char *cvec) {
    if(argc>1) strcpy(cmat, argv[1]);
    if(argc>2) strcpy(cvec, argv[2]);
}

// Lectura de la matriz en el archivo cfile.
// Devuelve NULL si no se pudo abrir el archivo.
gsl_matrix *readMatrix(char *cfile) {
    gsl_matrix *mat;
    int         nr, nc;
    FILE       *f1 = fopen(cfile, "rb");
    
    if(!f1) return(NULL);
    fread(&nr, sizeof(int), 1, f1);
    fread(&nc, sizeof(int), 1, f1);
    mat = gsl_matrix_alloc(nr, nc);
    gsl_matrix_fread(f1, mat);
    fclose(f1);
    return(mat);
}


// Lectura del vector en el archivo cfile.
// Devuelve NULL si no se pudo abrir el archivo.
gsl_vector *readVector(char *cfile) {
    gsl_vector *vec;
    int         nr;
    FILE       *f1 = fopen(cfile, "rb");
    
    if(!f1) return(NULL);
    fread(&nr, sizeof(int), 1, f1);
    vec = gsl_vector_alloc(nr);
    gsl_vector_fread(f1, vec);
    fclose(f1);
    return(vec);
}


// Devuelve un doble apuntador para poder trabajar con las
// entradas de las matriz mat con un arreglo bidimensional
double **gslmat2array2d(gsl_matrix *mat) {
    double **array2d;
    unsigned int i;
    
    array2d = (double **) malloc(mat->size1 * sizeof(double *));
    if(array2d==NULL) return(NULL);

    for(i=0; i<mat->size1; ++i)
       array2d[i] = mat->data + mat->size2*i;
    return(array2d);    
}



// Imprime en la consola las entradas del vector
void printVector(gsl_vector *vec) {
    unsigned int i; 
    double *ptr = vec->data;
    for(i=0; i<vec->size; ++i) 
	    printf("% 6.2f   ", ptr[i]);
	printf("\n");
}


// Calcula la norma del error Ax-b
double normError(gsl_matrix *matA, gsl_vector *vecx, gsl_vector *vecb) {
    unsigned int i, j; 
    double **mat = gslmat2array2d(matA), *vb=vecb->data, *vx=vecx->data;
    double derror, dval;

    derror = 0.0;
    for(i=0; i<matA->size1; ++i) {
        dval = -vb[i];
        for(j=0; j<matA->size2; ++j) 
            dval += mat[i][j]*vx[j];
        derror += dval*dval;
	}

    free(mat);
    return(sqrt(derror));
}


void solve(gsl_matrix *matA, gsl_vector *vecb, gsl_vector *vecx) {
    int              signum;
    gsl_permutation *perm;
    gsl_matrix      *matLU;
 
    matLU = gsl_matrix_alloc(matA->size1, matA->size2);
    
    gsl_matrix_memcpy(matLU, matA);
 
    // Calcula la factorizacion LU y resuelve el sistema
    perm = gsl_permutation_alloc(matA->size2);
    gsl_linalg_LU_decomp(matLU, perm, &signum);
    gsl_linalg_LU_solve(matLU, perm, vecb, vecx);

    gsl_permutation_free(perm);
    gsl_matrix_free(matLU); 
}
 
int main(int argc, char **argv)      {
    char    cfile1[50], cfile2[50];

    // Lectura de los parametros
    readParams(argc, argv, cfile1, cfile2);

    // Lectura e impresion de la matriz
    gsl_matrix *matA = readMatrix(cfile1);
    if(matA==NULL) {
        printf("El archivo %s no puedo ser leido.\n", cfile1);
        return(0);
    }

    // Lectura e impresion del vector
    gsl_vector *vecb = readVector(cfile2);
    if(vecb==NULL) {
        printf("El archivo %s no puedo ser leido.\n", cfile2);
        return(0);
    }
    printf("Tama~no de la matriz: %d x %d\n", (int)matA->size1, (int)matA->size2);

    gsl_vector *vecx;
    double     derror;

    vecx  = gsl_vector_alloc(matA->size2);
    gsl_block *b1=vecx->block, *b2=matA->block;

    solve(matA, vecb, vecx);

    // Imprime el vector solucion
    if(vecx->size<20) { 
        printf("Solucion del sistema:\n"); 
        printVector(vecx);
    }

    // Calculamos el error 
    derror = normError(matA, vecx, vecb);
    printf("Error en la solucion del sistema: %e\n", derror); 

    gsl_matrix_free(matA);
    gsl_vector_free(vecb); 
    gsl_vector_free(vecx); 

    return 0;
}
// Para compilar
// gcc -o luGsl luGsl.c  -lgsl -lgslcblas -lm

// Para ejecutar:
// ./luGsl matA1.bin vecb1.bin

