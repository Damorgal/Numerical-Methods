
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>


// Lee el valor de los parametros y los asigna a las variables:
// cdat   - Nombre del archivo de datos
// dft    - factor de escala
void readParams(int argc, char **argv, char *cdat, double *dft) {
    if(argc>1) strcpy(cdat, argv[1]);
    if(argc>2) *dft = atof(argv[2]);
    else   *dft = 1.0;
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

// Escritura de la matriz mat en el archivo cfile.
// Devuelve 0 en caso de exito y 1 si no.
int writeMatrix(gsl_matrix *mat, char *cfile) {
    FILE       *f1 = fopen(cfile, "wb");
    
    if(!f1) return(1);
    fwrite(&(mat->size1), sizeof(int), 1, f1);
    fwrite(&(mat->size2), sizeof(int), 1, f1);
    gsl_matrix_fwrite(f1, mat);
    fclose(f1);
    return(0);
}


// Imprime en la consola las entradas de la matriz
void printMatrix(gsl_matrix *mat) {
    unsigned int i, j; 
    for(i=0; i<mat->size1; ++i) {
        for(j=0; j<mat->size2; ++j) 
            printf("% 6.2f   ", gsl_matrix_get(mat, i, j));
        printf("\n");
    }
}

// Multiplica por dft a los elementos de la matriz
void scaleMatrix(gsl_matrix *mat, double dft) {
    unsigned int i, j; 
    double   dval;
    for(i=0; i<mat->size1; ++i) 
	for(j=0; j<mat->size2; ++j) {
	    dval = gsl_matrix_get(mat, i, j);
	    gsl_matrix_set(mat, i, j, dft*dval);
	}
}


int main(int argc, char **argv)      {
    char    cfile1[50], cfile2[50];
    double  dft;

    // Lectura de los parametros
    readParams(argc, argv, cfile1, &dft);

    // Lectura e impresion de la matriz
    gsl_matrix *mat = readMatrix(cfile1);
    if(mat==NULL) {
        printf("El archivo %s no puedo ser leido.\n", cfile1);
        return(0);
    }
    printf("\nMatriz en el archivo %s:\n", cfile1);
    printMatrix(mat);

    // Escalamiento de la matriz y almacenamiento
    scaleMatrix(mat, dft);

    printf("\nMatriz escalada por el factor %f:\n", dft);
    printMatrix(mat);

    sprintf(cfile2, "scaled%s", cfile1);
    printf("\nSalida almacena en %s\n", cfile2);
    writeMatrix(mat, cfile2);

    gsl_matrix_free(mat);
    return(0);
}
// Para compilar
// gcc -o matricesGsl1 matricesGsl1.c  -lgsl -lgslcblas -lm

// Para ejecutar
// ./matricesGsl1 matA1.bin  3


