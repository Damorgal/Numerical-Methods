#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* resuelve_triangularL(double** A, double* b, int n, int *ind);
void calcula_errorM(double **A, double *b, double *x, int n);

double *lee_arreglo(char *dir, int *n);
double **lee_matriz(char *dir, int *nr, int *nc);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
void freeMat(double** mat);

int main(int narg, char **varg)  {
    int n,nr,nc;
    char *dirmat = varg[1];
    char *dirvec = varg[2];

    if(narg != 3) {printf("\nIntroduce la direccion de los DOS archivos binarios (matriz y sol. b)\n"); return 0;}
    // Leemos los archivos
    double **A = lee_matriz(dirmat,&nr,&nc);
    if(nr != nc) {printf("Error, las matrices en el programa deben ser cuadradas.\n"); freeMat(A); return 0;}
    double *b = lee_arreglo(dirvec, &n);
    if(n != nc) {printf("Error, la dimension del vector b y la matriz no coinciden.\n"); freeMat(A); free(b); return 0;}

	// Resolvemos el sistema
	// Indicador que nos dira si el sistema tendra infinitas soluciones
    int ind=0;
	double *x = resuelve_triangularL(A,b,n,&ind);

    // Salida del programa
    if(x == NULL) printf("La matriz no tiene solucion.\n");
    else if(ind == 1) printf("La matriz tiene infinitas soluciones.\n");
    else    {
        printf("\nEl tamagno de la matriz es: %dx%d\n", nr, nc);
        printf("La matriz SI tiene solucion\n");
        calcula_errorM(A,b,x,n);
    }
	free(x);
    free(b);
    freeMat(A);
	return 0;
}

// Resuelve una matriz triangular inferior
double* resuelve_triangularL(double** A, double* b, int n, int *ind)   {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        double sum=0;
        if(i == 0)  {
            if(A[i][i] != 0)
                x[i] = b[i]/A[i][i];
            else if(b[i] == 0)
                *ind=1;
            else {free(x); return NULL;}
        }
        else    {
            for(int j=0; j<i; j++)
                sum += A[i][j]*x[j];
            if(A[i][i] != 0)
                x[i] = (b[i] - sum)/A[i][i];
            else if((b[i]-sum) == 0)
                *ind=1;
            else {free(x); return NULL;}
        }
    }
    return x;
}

// Calcula el error de que una matriz por un vector x de la solucion del vector b
void calcula_errorM(double **A, double *b, double *x, int n)  {
    double sum=0;
    for(int i=0; i<n; i++)  {
        double aux=0;
        for(int j=0; j<n; j++)
            aux += A[i][j]*x[j];
        sum += std::pow(aux - b[i], 2);
    }
    sum = std::sqrt(sum);
    printf("Y el error calculado de la solucion es: %lf\n", sum);
}

// Funcion que lee de archivo un vector
double *lee_arreglo(char *dir, int *n) {
    FILE *fp;
    double *vec;

	fp = fopen(dir, "rb");
	if (fp==NULL) {printf("El archivo %s no puedo ser leido.\n", dir); exit (1);}

	// Pedir memoria
	fread(n, sizeof(int), 1, fp);
    vec = crea_vector(*n);

	// Leer arreglo
	fread(vec, sizeof(double), *n, fp);

	fclose(fp);
	return(vec);
}

// Funcion que lee de archivo una matriz
double **lee_matriz(char *dir, int *nr, int *nc) {
    double **mat;

    FILE *fp = fopen(dir, "rb");
    if (fp==NULL) {printf("El archivo %s no puedo ser leido.\n", dir); exit (1);}

    // Lectura del tama~no de la matriz
    fread(nr, sizeof(int), 1, fp);
    fread(nc, sizeof(int), 1, fp);
    // Reservamos memoria
    mat = crea_matriz(*nr, *nc);
    // Lectura de los datos
    fread(mat[0], sizeof(double), (*nr)*(*nc), fp);
    fclose(fp);
    return(mat);
}

double* crea_vector(int n)  {
    double* aux;
    aux = (double*)malloc(sizeof(double)*n);
    return aux;
}

double** crea_matriz(int nr, int nc) {
    int i;
    double** mat;

    mat = (double**)malloc(nr * sizeof(double*));
    if(mat==NULL) return(NULL);
    mat[0] = (double*)malloc(nr * nc * sizeof(double));
    if(mat[0]==NULL) return(NULL);

    for(i=1; i<nr; i++)
        mat[i] = mat[i-1] + nc;

    return mat;
}

// Libera memoria de matriz
void freeMat(double** mat)   {
    free(mat[0]);
    free(mat);
}
