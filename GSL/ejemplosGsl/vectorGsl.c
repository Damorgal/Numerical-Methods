#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>


// Lee el valor de los parametros y los asigna a las variables:
// cdat   - Nombre del archivo de datos
// dft    - factor de escala
void readParams(int argc, char **argv, char *cdat, double *dft) {
    if(argc>1) strcpy(cdat, argv[1]);
    if(argc>2) *dft = atof(argv[2]);
    else   *dft = 1.0;
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

// Escritura del vector vec en el archivo cfile.
// Devuelve 0 en caso de exito y 1 si no.
int writeVector(gsl_vector *vec, char *cfile) {
    FILE       *f1 = fopen(cfile, "wb");
    
    if(!f1) return(1);
    fwrite(&(vec->size), sizeof(int), 1, f1);
    gsl_vector_fwrite(f1, vec);
    fclose(f1);
    return(0);
}


// Imprime en la consola las entradas del vector
void printVector(gsl_vector *vec) {
    unsigned int i; 
    double *ptr = vec->data;
    for(i=0; i<vec->size; ++i) 
	    printf("% 6.2f   ", ptr[i]);
	printf("\n");
}

// Multiplica por dft a los elementos de la matriz
void scaleVector(gsl_vector *vec, double dft) {
    unsigned int i; 
    double *ptr = vec->data;
    for(i=0; i<vec->size; ++i)
	    ptr[i] *= dft;
}


int main(int argc, char **argv)      {
    char    cfile1[50], cfile2[50];
    double  dft;

    // Lectura de los parametros
    readParams(argc, argv, cfile1, &dft);

    // Lectura e impresion de la matriz
    gsl_vector *vec = readVector(cfile1);
    if(vec==NULL) {
        printf("El archivo %s no puedo ser leido.\n", cfile1);
        return(0);
    }
    printf("\nVector en el archivo %s:\n", cfile1);
    printVector(vec);

    // Escalamiento del vector y almacenamiento
    scaleVector(vec, dft);


    printf("\nVector escalado por el factor %f:\n", dft);
    printVector(vec);

    sprintf(cfile2, "scaled%s", cfile1);
    printf("\nSalida almacena en %s\n", cfile2);
    writeVector(vec, cfile2);

    gsl_vector_free(vec);
    return(0);
}
// Para compilar
// gcc -o vectorGsl vectorGsl.c  -lgsl -lgslcblas -lm

// Para ejecutar
// ./vectorGsl v1.bin  3



